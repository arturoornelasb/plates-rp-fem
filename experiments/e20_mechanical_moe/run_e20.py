#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E20 -- mechanical mixture-of-experts (RUNNER). See README.md (frozen).
One certified plate solve (cached), then modal-space sweeps."""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, rectangle_basis,
                      solve_modes, split_rigid)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33, mesh=(96, 60),
    M=100, M_alt=[80, 120],
    K_sweep=[4, 6, 8, 12, 16],
    K_r1=8,
    n_conn=3, n_posts=2,
    q_rel=[0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
    detune=0.05,
    seeds=[0, 1, 2, 3, 4, 5],
    window=0.25,
    seed_plate=5,
    R_POI=0.3863, R_GOE=0.5307,
)


def plate_modes():
    """One certified plate: eigenvalues + mode-shape probe machinery."""
    path = os.path.join(HERE, "plate_cache.npz")
    cfg = CFG
    mesh, basis = rectangle_basis(*cfg["mesh"], cfg["a"], cfg["b"],
                                  ElementTriArgyris)
    if os.path.exists(path):
        z = np.load(path)
        return z["lam"], z["V"], basis
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V, sinfo, _ = solve_modes(K, M, max(cfg["M_alt"]) + 43,
                                   resid_sanity=1e-3, sweeps_max=30)
    lam, V, n_rigid, _ = split_rigid(lam, V)
    np.savez(path, lam=lam, V=V)
    print(f"[plate] {len(lam)} elastic cached, resid "
          f"{sinfo['max_resid']:.1e}")
    return lam, V, basis


def build_system(lam0, PHI, K, M, q, rng, cfg):
    """Modal MoE: block-diag detuned plates + post-spring coupling.
    PHI: (n_pts_pool, n_modes) mode values at a pooled random point set;
    each post draws a row index. Returns (H, per) with per = M."""
    det = 1.0 + cfg["detune"] * (2 * rng.random(K) - 1.0)
    H = np.zeros((K * M, K * M))
    for k in range(K):
        H[k * M:(k + 1) * M, k * M:(k + 1) * M] = np.diag(lam0[:M] * det[k])
    # sparse random contact graph
    for k in range(K):
        partners = rng.choice([j for j in range(K) if j != k],
                              size=min(cfg["n_conn"], K - 1), replace=False)
        for j in partners:
            for _ in range(cfg["n_posts"]):
                p = rng.integers(0, PHI.shape[0])
                fa = PHI[p, :M]
                fb = PHI[rng.integers(0, PHI.shape[0]), :M]
                # spring q (w_a(x) - w_b(y))^2 -> modal blocks
                H[k * M:(k + 1) * M, k * M:(k + 1) * M] += q * np.outer(fa, fa)
                H[j * M:(j + 1) * M, j * M:(j + 1) * M] += q * np.outer(fb, fb)
                blk = -q * np.outer(fa, fb)
                H[k * M:(k + 1) * M, j * M:(j + 1) * M] += blk
                H[j * M:(j + 1) * M, k * M:(k + 1) * M] += blk.T
    return 0.5 * (H + H.T)


def rt_mid(ev, frac):
    n = len(ev)
    w0, w1 = int((0.5 - frac / 2) * n), int((0.5 + frac / 2) * n)
    s = np.diff(np.sort(ev)[w0:w1])
    s = s[s > 0]
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r


def nsec_mid(H, K, M, frac):
    ev, U = np.linalg.eigh(H)
    n = len(ev)
    w0, w1 = int((0.5 - frac / 2) * n), int((0.5 + frac / 2) * n)
    W = U[:, w0:w1].reshape(K, M, w1 - w0)
    w = np.sum(W ** 2, axis=1)
    w = w / (w.sum(axis=0, keepdims=True) + 1e-300)
    return 1.0 / np.sum(w ** 2, axis=0), ev, (w0, w1)


def lambda_eff(H, K, M):
    offs, gaps = [], []
    ev = np.linalg.eigvalsh(H)
    d = np.median(np.diff(np.sort(ev)))
    for k in range(K):
        for j in range(k + 1, K):
            blk = H[k * M:(k + 1) * M, j * M:(j + 1) * M]
            m = np.abs(blk[blk != 0])
            if len(m):
                offs.append(np.median(m))
    return float(np.median(offs) / d) if offs else 0.0


def main():
    t00 = time.time()
    cfg = CFG
    lam0, V, basis = plate_modes()
    rngp = np.random.default_rng(cfg["seed_plate"])
    pts = np.vstack([rngp.uniform(0.06 * cfg["a"], 0.94 * cfg["a"], 400),
                     rngp.uniform(0.06 * cfg["b"], 0.94 * cfg["b"], 400)])
    PHI = (basis.probes(pts) @ V[:, :max(cfg["M_alt"])])
    PHI = np.asarray(PHI)
    d0 = float(np.median(np.diff(lam0[:cfg["M"]])))
    results = {"config": {k: (list(v) if isinstance(v, (list, tuple))
                              else v) for k, v in cfg.items()
                          if k not in ("R_POI", "R_GOE")}}
    print(f"[setup] plate + probes ({time.time()-t00:.0f} s)")

    # ---------------- R1: the transition at K_r1 ----------------
    md = ["# E20 -- mechanical mixture-of-experts (RESULTS)\n"]
    md.append("## R1: the fabricated transition (K = 8)\n")
    md.append("| q/Delta | lambda_eff | pooled r-tilde |")
    md.append("|---|---|---|")
    r1_rows = []
    for qr in cfg["q_rel"]:
        q = qr * d0
        rs, les = [], []
        for sd in cfg["seeds"]:
            rng = np.random.default_rng(100 + sd)
            H = build_system(lam0, PHI, cfg["K_r1"], cfg["M"], q, rng, cfg)
            rs.extend(rt_mid(np.linalg.eigvalsh(H), cfg["window"]).tolist())
            les.append(lambda_eff(H, cfg["K_r1"], cfg["M"]))
        rs = np.array(rs)
        r1_rows.append(dict(q_rel=qr, lam_eff=float(np.mean(les)),
                            r=float(rs.mean()),
                            sem=float(rs.std() / np.sqrt(len(rs)))))
        md.append(f"| {qr:g} | {r1_rows[-1]['lam_eff']:.3f} | "
                  f"{r1_rows[-1]['r']:.4f}({r1_rows[-1]['sem']:.4f}) |")
        print(f"[R1 q={qr:g}] r = {r1_rows[-1]['r']:.4f} "
              f"({time.time()-t00:.0f} s)", flush=True)
    results["r1"] = r1_rows
    rvals = np.array([r["r"] for r in r1_rows])
    mono = all(np.diff(rvals) > -2 * np.array(
        [r["sem"] for r in r1_rows[1:]]))
    reach = rvals.max() >= cfg["R_GOE"] - 2 * r1_rows[int(np.argmax(rvals))]["sem"]
    r1_verdict = ("SUPPORTS (monotone Poisson -> GOE)" if mono and reach
                  else ("saturates below GOE" if mono else "non-monotone"))
    md.append(f"\n- R1: {r1_verdict}")

    # transition coupling q*: pooled r closest to midpoint
    mid = 0.5 * (cfg["R_POI"] + cfg["R_GOE"])
    qstar = cfg["q_rel"][int(np.argmin(np.abs(rvals - mid)))]
    md.append(f"- transition coupling q* = {qstar:g} (r closest to "
              f"midpoint {mid:.3f})\n")

    # ---------------- R2/R3: sector D2 at q* ----------------
    md.append(f"## R2/R3: sector multifractality at q* = {qstar:g}\n")
    md.append("| K | mean N_sec | N_sec/K |")
    md.append("|---|---|---|")
    q = qstar * d0
    nsec_K = {}
    for K in cfg["K_sweep"]:
        vals = []
        for sd in cfg["seeds"]:
            rng = np.random.default_rng(200 + sd)
            H = build_system(lam0, PHI, K, cfg["M"], q, rng, cfg)
            ns, _, _ = nsec_mid(H, K, cfg["M"], cfg["window"])
            vals.append(float(np.mean(ns)))
        nsec_K[K] = float(np.mean(vals))
        md.append(f"| {K} | {nsec_K[K]:.3f} | {nsec_K[K]/K:.3f} |")
        print(f"[R2 K={K}] N_sec = {nsec_K[K]:.2f} "
              f"({time.time()-t00:.0f} s)", flush=True)
    Ks = np.array(cfg["K_sweep"], float)
    D2 = float(np.polyfit(np.log(Ks),
                          np.log([nsec_K[k] for k in cfg["K_sweep"]]),
                          1)[0])
    results["nsec_K"] = nsec_K
    results["D2_mech"] = D2

    # truncation stability gate
    stab = {}
    for Ma in cfg["M_alt"]:
        vals_lo, vals_hi = [], []
        for sd in cfg["seeds"][:3]:
            rng = np.random.default_rng(200 + sd)
            H = build_system(lam0, PHI, 4, Ma, q, rng, cfg)
            vals_lo.append(float(np.mean(nsec_mid(H, 4, Ma,
                                                  cfg["window"])[0])))
            rng = np.random.default_rng(200 + sd)
            H = build_system(lam0, PHI, 16, Ma, q, rng, cfg)
            vals_hi.append(float(np.mean(nsec_mid(H, 16, Ma,
                                                  cfg["window"])[0])))
        stab[Ma] = float((np.log(np.mean(vals_hi))
                          - np.log(np.mean(vals_lo)))
                         / (np.log(16) - np.log(4)))
    results["D2_M_stability"] = stab

    # audit-style ladder contrast at K = 16 (protocol immunity)
    rng = np.random.default_rng(200)
    H = build_system(lam0, PHI, 16, cfg["M"], q, rng, cfg)
    H0 = np.zeros_like(H)
    for k in range(16):
        sl = slice(k * cfg["M"], (k + 1) * cfg["M"])
        H0[sl, sl] = H[sl, sl]
    w0_, U0 = np.linalg.eigh(H0)
    order = np.argsort(w0_)
    U0 = U0[:, order]
    Hrep = U0.T @ H @ U0
    lad = []
    for N in [128, 256, 512, 1024]:
        i0, i1 = int(0.4 * N), int(0.6 * N)
        wN, UN = np.linalg.eigh(0.5 * (Hrep[:N, :N] + Hrep[:N, :N].T))
        lad.append((N, float(np.mean(np.log(
            np.sum(UN[:, i0:i1] ** 4, axis=0))))))
    D2_lad = float(-np.polyfit(np.log([x[0] for x in lad]),
                               [x[1] for x in lad], 1)[0])
    results["D2_ladder_contrast"] = D2_lad

    md.append(f"\n- **D2_mech = {D2:.3f}** (K-sweep, canonical sector "
              f"protocol); M-stability {stab}; audit-style truncated "
              f"ladder on the same operator: {D2_lad:.3f}")
    AI_BAND = (0.56, 0.80)
    if 0.2 < D2 < 0.9:
        r2 = "GENUINE RP-multifractal phase, fabricated mechanically"
    elif D2 >= 0.9:
        r2 = "ERGODIC (reaches GOE fixed point)"
    else:
        r2 = "localized-mimicry"
    r3 = ("INSIDE the AI band [0.56, 0.80] -- the fabrication bridge "
          "closes QUANTITATIVELY on both substrates"
          if AI_BAND[0] <= D2 <= AI_BAND[1] else
          f"outside the AI band [0.56, 0.80]")
    md.append(f"- R2: {r2}")
    md.append(f"- R3: {r3}")
    verdict = f"R1 {r1_verdict}; R2 {r2}; R3 {r3}."
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-8:]))


if __name__ == "__main__":
    main()
