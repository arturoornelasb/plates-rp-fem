#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E12 -- P7 non-proportional damping -> AI-dagger (RUNNER + ANALYSIS).
Complex-spacing-ratio estimator and linearization ported from the repo's
p7_pretest.py (conventions identical). See README.md (preregistered)."""
import json
import os
import time

import numpy as np
from scipy.linalg import eig

from platefem import (ElementTriArgyris, assemble_plate, classify_parity_resolved,
                      probe_operators, rectangle_basis, solve_modes, split_rigid,
                      SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(96, 60), n_modes=1203,
    n_per_sector=295,
    patch_center=(0.2231, 0.1372), patch_radius=0.06, patch_npts=200,
    # v2: DISTRIBUTED dissipation -- 24 patch quartets (single-patch damping
    # proved only weakly non-proportional: commutator 3e-2, near-diagonal C,
    # Poisson-like resonances; the P7 dense case needs distributed contacts)
    n_patches=24,
    gammas=[2.0, 4.0, 8.0, 16.0],
    rayleigh=(2e-3, 2e-4),
    seed=7, n_baseline=1200, baseline_reps=6,
)


def complex_spacing_ratios(lam):
    """Ported from p7_pretest.py (bulk-filtered z-ratios)."""
    lam = np.asarray(lam, complex)
    if len(lam) < 4:
        return np.nan, np.nan, 0
    cx, cy = lam.real, lam.imag
    keep = ((cx > np.quantile(cx, .15)) & (cx < np.quantile(cx, .85)) &
            (cy > np.quantile(cy, .15)) & (cy < np.quantile(cy, .85)))
    zs = []
    for i in np.where(keep)[0]:
        d = np.abs(lam - lam[i]); d[i] = np.inf
        order = np.argsort(d)
        denom = lam[order[1]] - lam[i]
        if denom != 0:
            zs.append((lam[order[0]] - lam[i]) / denom)
    zs = np.array(zs)
    if len(zs) < 5:
        return np.nan, np.nan, 0
    return (float(np.abs(zs).mean()),
            float(-np.cos(np.angle(zs)).mean()), len(zs))


def baselines(n, reps, rng):
    out = {}
    accum = {"Poisson2D": [], "GinUE": [], "AI_dagger": []}
    for _ in range(reps):
        rad = np.sqrt(rng.random(n)); ang = 2 * np.pi * rng.random(n)
        accum["Poisson2D"].append(
            complex_spacing_ratios(rad * np.exp(1j * ang))[:2])
        G = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        accum["GinUE"].append(complex_spacing_ratios(np.linalg.eigvals(G))[:2])
        A_ = rng.standard_normal((n, n)); B_ = rng.standard_normal((n, n))
        H = (A_ + A_.T) / 2 + 1j * (B_ + B_.T) / 2
        accum["AI_dagger"].append(
            complex_spacing_ratios(np.linalg.eigvals(H))[:2])
    for k, v in accum.items():
        v = np.array(v)
        out[k] = dict(r=float(np.nanmean(v[:, 0])),
                      r_se=float(np.nanstd(v[:, 0]) / np.sqrt(reps)),
                      c=float(np.nanmean(v[:, 1])),
                      c_se=float(np.nanstd(v[:, 1]) / np.sqrt(reps)))
    return out


def qep_spectrum(Lam, C):
    """Upper-half spectrum of the QEP (modal M = I, K = diag(Lam));
    linearization as in p7_pretest."""
    N = len(Lam)
    K = np.diag(Lam)
    Z = np.zeros((N, N))
    A = np.block([[Z, K], [K, C]])
    B = np.block([[K, Z], [Z, -np.eye(N)]])
    s = eig(A, B, right=False)
    s = s[np.isfinite(s)]
    # conjugate pairing error
    su = s[np.abs(s.imag) > 1e-9]
    miss = 0
    for z in su:
        if np.min(np.abs(su - np.conj(z))) > 1e-6 * max(1.0, abs(z)):
            miss += 1
    err = miss / max(len(su), 1)
    return s[s.imag > 0], float(err)


def main():
    t00 = time.time()
    cfg = CFG
    a, b = cfg["a"], cfg["b"]
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}}

    # baselines at matched n
    results["baselines"] = baselines(cfg["n_baseline"], cfg["baseline_reps"],
                                     rng)
    print("[baselines]",
          {k: (round(v["r"], 4), round(v["c"], 4))
           for k, v in results["baselines"].items()})

    # certified modes + patch quartet damping
    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-4,
                                   sweeps_max=30)
    lam, V, n_rigid, _ = split_rigid(lam, V)
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    print(f"[modes] rigid {n_rigid}, resid {sinfo['max_resid']:.1e} "
          f"({time.time()-t00:.1f} s)")

    # patch quadrature points (disk quartet, mirror images -> sector-pure)
    quartet = []
    for _ in range(cfg["n_patches"]):
        x0 = rng.uniform(0.06 * a, 0.44 * a)
        y0 = rng.uniform(0.06 * b, 0.44 * b)
        npp = max(cfg["patch_npts"] // cfg["n_patches"], 12)
        rr = cfg["patch_radius"] * np.sqrt(rng.random(npp))
        th = 2 * np.pi * rng.random(npp)
        px = np.clip(x0 + rr * np.cos(th), 0.005 * a, 0.495 * a)
        py = np.clip(y0 + rr * np.sin(th), 0.005 * b, 0.495 * b)
        quartet += [np.vstack([px, py]), np.vstack([a - px, py]),
                    np.vstack([px, b - py]), np.vstack([a - px, b - py])]

    results["cases"] = {}
    for s in SECTORS:
        idx = [i for i, l in enumerate(labels) if l == s][:cfg["n_per_sector"]]
        Lam = lam[idx]
        d_gam = float(np.mean(np.diff(Lam)))          # in Lambda units
        # modal patch damping: C_ij = sum_quartet sum_pts phi_i phi_j
        Cpatch = np.zeros((len(idx), len(idx)))
        for q in quartet:
            S = basis.probes(q) @ V[:, idx]
            Cpatch += S.T @ S
        Cpatch /= cfg["patch_npts"]
        Cpoint = np.zeros_like(Cpatch)
        for q in quartet:
            Sp_ = basis.probes(q[:, :1]) @ V[:, idx]
            Cpoint += Sp_.T @ Sp_
        # commutator diagnostic (non-proportionality)
        Kd = np.diag(Lam)
        comm = float(np.linalg.norm(Cpatch @ Kd - Kd @ Cpatch)
                     / (np.linalg.norm(Cpatch) * np.linalg.norm(Kd)))
        results["cases"][s] = dict(commutator=comm, runs={})
        # frequency-plane spacing: resonances live in sqrt(Lambda) ~ omega;
        # widths ~ eig shifts of the QEP. Scale c0 so median width / mean
        # omega-spacing = gamma*.
        om = np.sqrt(Lam)
        d_om = float(np.mean(np.diff(om)))
        for g in cfg["gammas"]:
            # half-width ~ c0 * Cii / 2 in the omega plane (light damping)
            c0 = g * d_om * 2.0 / np.median(np.diag(Cpatch))
            sp, perr = qep_spectrum(Lam, c0 * Cpatch)
            r, c, nz = complex_spacing_ratios(sp)
            results["cases"][s]["runs"][f"g{g:g}"] = dict(
                c0=float(c0), pairing_err=perr, r=r, c=c, n=nz)
        c0 = 4.0 * d_om * 2.0 / np.median(np.diag(Cpoint))
        sp, perr = qep_spectrum(Lam, c0 * Cpoint)
        r, c, nz = complex_spacing_ratios(sp)
        results["cases"][s]["point_g4"] = dict(r=r, c=c, n=nz)
        # proportional control at gamma* = 1
        alpha, beta = cfg["rayleigh"]
        c0 = 1.0 * d_om * 2.0 / np.median(alpha + beta * Lam)
        sp, perr = qep_spectrum(Lam, c0 * (alpha * np.eye(len(idx))
                                           + beta * Kd))
        r, c, nz = complex_spacing_ratios(sp)
        results["cases"][s]["rayleigh"] = dict(r=r, c=c, n=nz,
                                               pairing_err=perr)
        print(f"[{s}] comm {comm:.2e}; nonprop g=4: "
              f"r={results['cases'][s]['runs']['g4']['r']:.4f} "
              f"c={results['cases'][s]['runs']['g4']['c']:.4f}; rayleigh: "
              f"r={r:.4f} c={c:.4f}")

    # ---------------- pooled verdict ----------------
    md = ["# E12 -- P7: non-proportional damping -> AI-dagger (RESULTS)\n",
          f"Certified FFFF modes ({cfg['n_per_sector']}/sector), viscous patch "
          f"quartet at {cfg['patch_center']}, r = {cfg['patch_radius']}. "
          f"Baselines at n = {cfg['n_baseline']}.\n"]
    bl = results["baselines"]
    md.append("| reference | <|z|> | -<cos theta> |")
    md.append("|---|---|---|")
    for k, v in bl.items():
        md.append(f"| {k} | {v['r']:.4f} +/- {v['r_se']:.4f} "
                  f"| {v['c']:.4f} +/- {v['c_se']:.4f} |")
    md.append("\n| case | <|z|> | -<cos theta> | n | pairing err |")
    md.append("|---|---|---|---|---|")
    pooled = {}
    for g in cfg["gammas"]:
        rs = [results["cases"][s]["runs"][f"g{g:g}"] for s in SECTORS]
        n_tot = sum(x["n"] for x in rs)
        r_m = float(np.sum([x["r"] * x["n"] for x in rs]) / n_tot)
        c_m = float(np.sum([x["c"] * x["n"] for x in rs]) / n_tot)
        pe = max(x["pairing_err"] for x in rs)
        pooled[g] = (r_m, c_m, n_tot)
        md.append(f"| non-prop gamma*={g:g} (pooled) | {r_m:.4f} | {c_m:.4f} "
                  f"| {n_tot} | {pe:.1e} |")
    rs = [results["cases"][s]["point_g4"] for s in SECTORS]
    n_tot = sum(x["n"] for x in rs)
    md.append(f"| point-dashpots gamma*=4 (pooled) | "
              f"{np.sum([x['r']*x['n'] for x in rs])/n_tot:.4f} | "
              f"{np.sum([x['c']*x['n'] for x in rs])/n_tot:.4f} | {n_tot} | - |")
    rs = [results["cases"][s]["rayleigh"] for s in SECTORS]
    n_tot = sum(x["n"] for x in rs)
    r_ray = float(np.sum([x["r"] * x["n"] for x in rs]) / n_tot)
    c_ray = float(np.sum([x["c"] * x["n"] for x in rs]) / n_tot)
    md.append(f"| Rayleigh control (pooled) | {r_ray:.4f} | {c_ray:.4f} "
              f"| {n_tot} | - |")

    r1, c1, n1 = pooled[4.0]
    se_est = 0.02 * np.sqrt(1200 / max(n1, 1))     # empirical marker sem scale
    d_ai = np.hypot(r1 - bl["AI_dagger"]["r"], c1 - bl["AI_dagger"]["c"])
    d_po = np.hypot(r1 - bl["Poisson2D"]["r"], c1 - bl["Poisson2D"]["c"])
    d_gi = np.hypot(r1 - bl["GinUE"]["r"], c1 - bl["GinUE"]["c"])
    md.append(f"\n- marker distances at gamma* = 4: to AI-dagger {d_ai:.4f}, "
              f"to Poisson2D {d_po:.4f}, to GinUE {d_gi:.4f} "
              f"(marker scale se ~ {se_est:.4f})")
    if d_ai < d_po and d_ai < 3 * se_est and d_po > 5 * se_est:
        verdict = ("SUPPORTS P7: the damped true plate lands on the AI-dagger "
                   "markers, separated from 2D Poisson and from the "
                   "proportional filament -- first AI-dagger realization on a "
                   "certified classical elastic model. (AI-dagger vs GinUE "
                   "fine call registered for the 2400-mode stretch.)")
    elif d_po < d_ai:
        verdict = "CHALLENGES (Poisson-like: no 2D repulsion at these widths)"
    else:
        verdict = "AMBIGUOUS -- see markers"
    md.append(f"\n**Reading: {verdict}**")
    results["pooled"] = {f"g{g:g}": pooled[g] for g in cfg["gammas"]}
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
