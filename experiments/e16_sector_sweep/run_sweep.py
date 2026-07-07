#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E16 -- sector angle-sweep: protocol-faithful replication of the
published COMSOL sector result (see README.md, preregistered 2026-07-07).

Per angle: first N_KEEP nonzero elastic modes of the free Kirchhoff disk
sector, S/A classified by the bisector mirror. Statistics pooled over the
sweep against a same-protocol pooled-Poisson baseline (R1); two-mesh spot
gate + fixed-angle matched-window references at 50 / 90 / 114.59 deg (R2);
hierarchy placement vs the E3c long ladder and the disk control (R3).
Primary reading on the Lambda = omega^2 scale; frequency scale reported
alongside, ungated (README amendment, pre-run)."""
import json
import os
import time

import numpy as np

from platefem import (assemble_plate, classify_parity_resolved, n_star,
                      sector_basis, solve_modes, split_rigid)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    nu=0.33, R=1.0,
    th_min_deg=30.0, th_max_deg=150.0, n_theta=200,
    n_keep=70,
    nrings=45, nrings_check=32,
    sigma=-10.0,
    drop_low=10,
    spacing_frac=0.1,
    probe_npts=800,
    n_surrogates=200,
    seed=16,
    ref_angles_deg=(50.0, 90.0, 114.59155902616465),
)
E3C_LADDER = (0.4294, 0.0068)   # E3c v2, 1600 modes, theta = 2.0 rad
DISK_CONTROL = 0.3905           # E3b semi-analytic class sequence


def sector_probes(basis, theta, R, npts, seed=17):
    """Interior points and their bisector-mirror images (E3c trick: the
    x-mirror is the identity, so 'ee' -> S and 'eo' -> A)."""
    rng = np.random.default_rng(seed)
    r = np.sqrt(rng.uniform(0.0, 0.94, npts)) * R
    th = rng.uniform(-0.49 * theta, 0.49 * theta, npts)
    pts = np.vstack([r * np.cos(th), r * np.sin(th)])
    P = basis.probes(pts)
    Pmy = basis.probes(np.vstack([pts[0], -pts[1]]))
    return P, P, Pmy


def solve_angle(theta, cfg, nrings):
    """First n_keep elastic modes + S/A labels at one sector angle."""
    mesh, basis = sector_basis(nrings, theta, cfg["R"])
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_keep"] + 9, cfg["sigma"])
    lam, V, n_rigid, rigid_max = split_rigid(lam, V)
    lam, V = lam[:cfg["n_keep"]], V[:, :cfg["n_keep"]]
    P, Pmx, Pmy = sector_probes(basis, theta, cfg["R"], cfg["probe_npts"])
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    sa = ["S" if l == "ee" else ("A" if l == "eo" else "x") for l in labels]
    return dict(lam=np.asarray(lam), labels=sa, n_rigid=int(n_rigid),
                rigid_max=float(rigid_max), min_q=float(np.min(qual)),
                resid=float(sinfo["max_resid"]), ndof=int(K.shape[0]))


def rt_ratios(levels, drop):
    """min/max consecutive-spacing ratios of one sorted sequence."""
    x = np.asarray(levels, float)[drop:]
    s = np.diff(x)
    s = s[s > 0]
    if len(s) < 2:
        return np.empty(0)
    return np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])


def pooled_stats(per_angle, cfg, scale):
    """Pool per-class ratios over the sweep. scale: 'lam' or 'freq'."""
    out = {}
    for cls in ("S", "A"):
        rs, lens = [], []
        for rec in per_angle:
            seq = [l for l, s in zip(rec["lam"], rec["labels"]) if s == cls]
            seq = np.sqrt(seq) if scale == "freq" else np.asarray(seq)
            r = rt_ratios(seq, cfg["drop_low"])
            rs.append(r)
            lens.append(max(0, len(seq) - cfg["drop_low"]))
        rs = np.concatenate(rs)
        out[cls] = dict(mean=float(rs.mean()), sem=float(rs.std() / np.sqrt(len(rs))),
                        n=int(len(rs)), lens=lens)
    both = np.concatenate([
        np.concatenate([rt_ratios(
            (np.sqrt if scale == "freq" else np.asarray)(
                [l for l, s in zip(rec["lam"], rec["labels"]) if s == cls]),
            cfg["drop_low"]) for rec in per_angle])
        for cls in ("S", "A")])
    out["pooled"] = dict(mean=float(both.mean()),
                         sem=float(both.std() / np.sqrt(len(both))),
                         n=int(len(both)))
    return out


def poisson_baseline(stats_lam, cfg, rng):
    """Same-protocol pooled-Poisson baseline: iid-exponential-spacing
    sequences of the SAME post-drop lengths, pooled identically."""
    lens = stats_lam["S"]["lens"] + stats_lam["A"]["lens"]
    means = []
    for _ in range(cfg["n_surrogates"]):
        rs = []
        for L in lens:
            if L < 3:
                continue
            lev = np.cumsum(rng.exponential(size=L))
            rs.append(rt_ratios(lev, 0))
        rs = np.concatenate(rs)
        means.append(rs.mean())
    return float(np.mean(means)), float(np.std(means))


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": dict(cfg), "gates": {}, "refs": {}, "sweep": []}

    # ---- R2 spot gates + fixed-angle references (both meshes) ----
    print("[refs] fixed-angle references + two-mesh spot gates")
    ref_ok = True
    for deg in cfg["ref_angles_deg"]:
        th = np.deg2rad(deg)
        t0 = time.time()
        fine = solve_angle(th, cfg, cfg["nrings"])
        coarse = solve_angle(th, cfg, cfg["nrings_check"])
        ns = n_star(fine["lam"], coarse["lam"], cfg["spacing_frac"])
        ok = ns >= cfg["n_keep"] and "x" not in fine["labels"]
        ref_ok &= ok
        rs = np.concatenate([rt_ratios(
            [l for l, s in zip(fine["lam"], fine["labels"]) if s == c],
            cfg["drop_low"]) for c in ("S", "A")])
        results["refs"][f"{deg:.2f}"] = dict(
            n_star=int(ns), ok=bool(ok), ndof=fine["ndof"],
            resid=fine["resid"], min_q=fine["min_q"],
            r_matched=float(rs.mean()),
            sem_matched=float(rs.std() / np.sqrt(len(rs))),
            n_ratios=int(len(rs)),
            lam=fine["lam"].tolist(), labels=fine["labels"])
        print(f"  [ref {deg:7.2f} deg] N* = {ns}/{cfg['n_keep']}, "
              f"r_matched = {rs.mean():.4f} ({len(rs)} ratios), "
              f"gate {'PASS' if ok else 'FAIL'} ({time.time()-t0:.1f} s)")
    results["gates"]["refs_ok"] = bool(ref_ok)

    # ---- the sweep ----
    angles = np.linspace(cfg["th_min_deg"], cfg["th_max_deg"], cfg["n_theta"])
    ambiguous = 0
    for i, deg in enumerate(angles):
        rec = solve_angle(np.deg2rad(deg), cfg, cfg["nrings"])
        ambiguous += rec["labels"].count("x")
        results["sweep"].append(dict(theta_deg=float(deg),
                                     lam=rec["lam"].tolist(),
                                     labels=rec["labels"],
                                     n_rigid=rec["n_rigid"],
                                     min_q=rec["min_q"], resid=rec["resid"]))
        if (i + 1) % 20 == 0 or i == len(angles) - 1:
            print(f"[sweep] {i+1}/{cfg['n_theta']} angles "
                  f"({time.time()-t00:.0f} s)")
            with open(os.path.join(HERE, "results_raw.json"), "w") as f:
                json.dump(results, f, default=float)
    results["gates"]["ambiguous_labels"] = int(ambiguous)

    # ---- statistics ----
    per_angle = [dict(lam=np.asarray(r["lam"]), labels=r["labels"])
                 for r in results["sweep"]]
    st_lam = pooled_stats(per_angle, cfg, "lam")
    st_frq = pooled_stats(per_angle, cfg, "freq")
    base_mean, base_std = poisson_baseline(st_lam, cfg, rng)

    pm, ps = st_lam["pooled"]["mean"], st_lam["pooled"]["sem"]
    sig = (pm - base_mean) / np.sqrt(ps**2 + base_std**2)
    refs = [results["refs"][k] for k in results["refs"]]
    ref_mean = float(np.mean([r["r_matched"] for r in refs]))
    ref_sem = float(np.sqrt(np.sum([r["sem_matched"]**2 for r in refs])) / len(refs))
    z_pool = (pm - ref_mean) / np.sqrt(ps**2 + ref_sem**2)

    results["stats"] = dict(
        pooled_lam=st_lam, pooled_freq=st_frq,
        poisson_baseline=dict(mean=base_mean, std=base_std),
        R1_sigma=float(sig),
        R2_ref_mean=ref_mean, R2_ref_sem=ref_sem, R2_z=float(z_pool),
        R3_e3c_ladder=E3C_LADDER, R3_disk=DISK_CONTROL)

    r1 = ("SUPPORTS the published direction" if sig >= 3 else
          ("direction positive, underpowered" if sig > 0 else
           "does NOT support"))
    r2 = ("pooling benign (within 2 sigma of fixed-angle matched-window mean)"
          if abs(z_pool) <= 2 else
          f"sweep protocol moves the statistic ({z_pool:+.1f} sigma)")
    md = [
        "# E16 -- sector angle-sweep (RESULTS)\n",
        f"Sweep theta in [{cfg['th_min_deg']:g}, {cfg['th_max_deg']:g}] deg, "
        f"{cfg['n_theta']} angles x first {cfg['n_keep']} elastic modes, "
        f"S/A by bisector mirror; drop-low {cfg['drop_low']} per class; "
        f"nrings {cfg['nrings']} (check {cfg['nrings_check']}).\n",
        f"- Gates: refs two-mesh {'PASS' if ref_ok else 'FAIL'} "
        f"(N* at {[results['refs'][k]['n_star'] for k in results['refs']]}); "
        f"ambiguous labels in sweep: {ambiguous}.",
        f"- Pooled <r~> (Lambda scale): {pm:.4f} +/- {ps:.4f} "
        f"({st_lam['pooled']['n']} ratios; S {st_lam['S']['mean']:.4f}, "
        f"A {st_lam['A']['mean']:.4f}).",
        f"- Pooled <r~> (frequency scale, source convention): "
        f"{st_frq['pooled']['mean']:.4f} +/- {st_frq['pooled']['sem']:.4f}.",
        f"- Same-protocol pooled-Poisson baseline: {base_mean:.4f} "
        f"+/- {base_std:.4f}.",
        f"- **R1**: excess {pm-base_mean:+.4f} = {sig:+.1f} sigma -> {r1}.",
        f"- **R2**: fixed-angle matched-window mean {ref_mean:.4f} "
        f"+/- {ref_sem:.4f} (50/90/114.59 deg); pooled-vs-fixed "
        f"{z_pool:+.1f} sigma -> {r2}.",
        f"- **R3**: hierarchy placement -- pooled sweep {pm:.4f} vs E3c "
        f"long ladder {E3C_LADDER[0]:.4f} +/- {E3C_LADDER[1]:.4f} and disk "
        f"control {DISK_CONTROL:.4f}.",
    ]
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, default=float)
    print("\n".join(md))


if __name__ == "__main__":
    main()
