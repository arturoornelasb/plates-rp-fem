#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E12b -- AI-dagger vs Ginibre fine call (RUNNER). See README.md
(frozen). Reuses the frozen E12 machinery; modal basis solved once,
markers pooled over 4 patch-placement realizations."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "e12_damping_AIdagger"))
from run_e12 import (CFG as E12CFG, baselines, complex_spacing_ratios,
                     qep_spectrum)

from platefem import (ElementTriArgyris, assemble_plate,
                      classify_parity_resolved, probe_operators,
                      rectangle_basis, solve_modes, split_rigid, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

SEEDS = [7, 17, 27, 37]
GSTAR = [16.0, 4.0]
N_BASE, BASE_REPS = 1200, 24


def quartets_for(rng, cfg, a, b):
    out = []
    for _ in range(cfg["n_patches"]):
        x0 = rng.uniform(0.06 * a, 0.44 * a)
        y0 = rng.uniform(0.06 * b, 0.44 * b)
        npp = max(cfg["patch_npts"] // cfg["n_patches"], 12)
        rr = cfg["patch_radius"] * np.sqrt(rng.random(npp))
        th = 2 * np.pi * rng.random(npp)
        px = np.clip(x0 + rr * np.cos(th), 0.005 * a, 0.495 * a)
        py = np.clip(y0 + rr * np.sin(th), 0.005 * b, 0.495 * b)
        out += [np.vstack([px, py]), np.vstack([a - px, py]),
                np.vstack([px, b - py]), np.vstack([a - px, b - py])]
    return out


def main():
    t00 = time.time()
    cfg = dict(E12CFG)
    a, b = cfg["a"], cfg["b"]
    rng = np.random.default_rng(101)
    results = {"config": dict(seeds=SEEDS, gstar=GSTAR, n_base=N_BASE,
                              base_reps=BASE_REPS)}

    results["baselines"] = baselines(N_BASE, BASE_REPS, rng)
    print("[baselines]", {k: (round(v["r"], 4), round(v["c"], 4))
                          for k, v in results["baselines"].items()})

    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-4,
                                   sweeps_max=30)
    lam, V, n_rigid, _ = split_rigid(lam, V)
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    print(f"[modes] rigid {n_rigid}, resid {sinfo['max_resid']:.1e} "
          f"({time.time()-t00:.1f} s)")

    pool = {g: dict(r=[], c=[], n=[]) for g in GSTAR}
    for seed in SEEDS:
        srng = np.random.default_rng(seed)
        quartet = quartets_for(srng, cfg, a, b)
        for s in SECTORS:
            idx = [i for i, l in enumerate(labels)
                   if l == s][:cfg["n_per_sector"]]
            Lam = lam[idx]
            Cpatch = np.zeros((len(idx), len(idx)))
            for q in quartet:
                S = basis.probes(q) @ V[:, idx]
                Cpatch += S.T @ S
            Cpatch /= cfg["patch_npts"]
            om = np.sqrt(Lam)
            d_om = float(np.mean(np.diff(om)))
            for g in GSTAR:
                c0 = g * d_om * 2.0 / np.median(np.diag(Cpatch))
                sp, perr = qep_spectrum(Lam, c0 * Cpatch)
                r, c, nz = complex_spacing_ratios(sp)
                pool[g]["r"].append(r * nz)
                pool[g]["c"].append(c * nz)
                pool[g]["n"].append(nz)
        print(f"[seed {seed}] done ({time.time()-t00:.1f} s)")

    md = ["# E12b -- AI-dagger vs Ginibre fine call (RESULTS)\n",
          f"{len(SEEDS)} patch realizations x 4 sectors, "
          f"{cfg['n_per_sector']}/sector; baselines n = {N_BASE} x "
          f"{BASE_REPS} reps.\n",
          "| reference | <|z|> | -<cos theta> |", "|---|---|---|"]
    bl = results["baselines"]
    for k, v in bl.items():
        md.append(f"| {k} | {v['r']:.4f} +/- {v['r_se']:.4f} | "
                  f"{v['c']:.4f} +/- {v['c_se']:.4f} |")
    md.append("\n| case | <|z|> | -<cos theta> | n |")
    md.append("|---|---|---|---|")
    marks = {}
    for g in GSTAR:
        n_tot = int(np.sum(pool[g]["n"]))
        r_m = float(np.sum(pool[g]["r"]) / n_tot)
        c_m = float(np.sum(pool[g]["c"]) / n_tot)
        marks[g] = (r_m, c_m, n_tot)
        md.append(f"| gamma* = {g:g} (pooled) | {r_m:.4f} | {c_m:.4f} "
                  f"| {n_tot} |")

    r1, c1, n1 = marks[16.0]
    se_marker = 0.02 * np.sqrt(1200 / max(n1, 1))
    se = float(np.hypot(se_marker,
                        max(bl["AI_dagger"]["r_se"], bl["GinUE"]["r_se"])))
    d_ai = float(np.hypot(r1 - bl["AI_dagger"]["r"],
                          c1 - bl["AI_dagger"]["c"]))
    d_gi = float(np.hypot(r1 - bl["GinUE"]["r"], c1 - bl["GinUE"]["c"]))
    d_po = float(np.hypot(r1 - bl["Poisson2D"]["r"],
                          c1 - bl["Poisson2D"]["c"]))
    md.append(f"\n- at gamma* = 16: d(AI-dagger) = {d_ai:.4f}, d(GinUE) = "
              f"{d_gi:.4f}, d(Poisson2D) = {d_po:.4f}; combined marker "
              f"scale se = {se:.4f}")
    if d_ai <= 2 * se and d_gi >= 4 * se:
        verdict = ("AI-DAGGER CONFIRMED over Ginibre -- the damped "
                   "reciprocal plate selects the complex-symmetric class, "
                   "completing P7's fine call")
    elif d_gi <= 2 * se and d_ai >= 4 * se:
        verdict = "GINIBRE -- the reciprocity signature is NOT retained"
    else:
        verdict = (f"UNRESOLVED (d_AI/se = {d_ai/se:.1f}, d_Gin/se = "
                   f"{d_gi/se:.1f}); registered next = physical-damping "
                   f"refinement, not more statistics")
    md.append(f"\n**Reading: {verdict}**")
    results["marks"] = {f"g{g:g}": marks[g] for g in GSTAR}
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
