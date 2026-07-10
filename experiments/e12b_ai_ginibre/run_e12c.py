#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E12c -- drive the damping crossover to saturation (preregistered in
header, 2026-07-09). E12b located the gamma* = 16 / 24-quartet point
MID-FLOW. E12c densifies and strengthens the non-proportional
dissipation (96 quartets; gamma* in {32, 64}; 2 patch realizations) and
applies the FROZEN saturation-then-endpoint reading: SATURATED if the
gamma* = 32 and 64 pooled markers agree within the combined scale se;
if saturated, the E12b endpoint gates apply (<= 2 se to one of
AI-dagger/GinUE AND >= 4 se from the other). Otherwise STILL-FLOWING
(report the trajectory)."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "e12_damping_AIdagger"))
from run_e12 import (CFG as E12CFG, baselines, complex_spacing_ratios,
                     qep_spectrum)
from run_e12b import quartets_for

from platefem import (ElementTriArgyris, assemble_plate,
                      classify_parity_resolved, probe_operators,
                      rectangle_basis, solve_modes, split_rigid, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = [7, 17]
GSTAR = [32.0, 64.0]


def main():
    t00 = time.time()
    cfg = dict(E12CFG)
    cfg["n_patches"] = 96
    a, b = cfg["a"], cfg["b"]
    rng = np.random.default_rng(202)
    results = {"config": dict(seeds=SEEDS, gstar=GSTAR,
                              n_patches=cfg["n_patches"])}
    results["baselines"] = baselines(1200, 24, rng)
    bl = results["baselines"]
    print("[baselines]", {k: (round(v["r"], 4), round(v["c"], 4))
                          for k, v in bl.items()})

    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-4,
                                   sweeps_max=30)
    lam, V, n_rigid, _ = split_rigid(lam, V)
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    print(f"[modes] resid {sinfo['max_resid']:.1e} "
          f"({time.time()-t00:.0f} s)")

    pool = {g: dict(r=[], c=[], n=[]) for g in GSTAR}
    for seed in SEEDS:
        srng = np.random.default_rng(seed)
        quartet = quartets_for(srng, cfg, a, b)
        for s in SECTORS:
            idx = [i for i, l in enumerate(labels)
                   if l == s][:cfg["n_per_sector"]]
            Lam = lam[idx]
            Cp = np.zeros((len(idx), len(idx)))
            for q in quartet:
                S_ = basis.probes(q) @ V[:, idx]
                Cp += S_.T @ S_
            Cp /= cfg["patch_npts"]
            om = np.sqrt(Lam)
            d_om = float(np.mean(np.diff(om)))
            for g in GSTAR:
                c0 = g * d_om * 2.0 / np.median(np.diag(Cp))
                sp, _ = qep_spectrum(Lam, c0 * Cp)
                r, c, nz = complex_spacing_ratios(sp)
                pool[g]["r"].append(r * nz)
                pool[g]["c"].append(c * nz)
                pool[g]["n"].append(nz)
        print(f"[seed {seed}] done ({time.time()-t00:.0f} s)", flush=True)

    md = ["# E12c -- saturation attempt (RESULTS)\n",
          "| case | <|z|> | -<cos theta> | n |", "|---|---|---|---|"]
    marks = {}
    for g in GSTAR:
        n_tot = int(np.sum(pool[g]["n"]))
        r_m = float(np.sum(pool[g]["r"]) / n_tot)
        c_m = float(np.sum(pool[g]["c"]) / n_tot)
        marks[g] = (r_m, c_m, n_tot)
        md.append(f"| gamma* = {g:g}, 96 quartets | {r_m:.4f} | "
                  f"{c_m:.4f} | {n_tot} |")
    se_m = 0.02 * np.sqrt(1200 / marks[64.0][2])
    se = float(np.hypot(se_m, bl["AI_dagger"]["r_se"]))
    gap = float(np.hypot(marks[32.0][0] - marks[64.0][0],
                         marks[32.0][1] - marks[64.0][1]))
    r1, c1, _ = marks[64.0]
    d_ai = float(np.hypot(r1 - bl["AI_dagger"]["r"],
                          c1 - bl["AI_dagger"]["c"]))
    d_gi = float(np.hypot(r1 - bl["GinUE"]["r"], c1 - bl["GinUE"]["c"]))
    md.append(f"\n- gamma 32 vs 64 gap = {gap:.4f} (se {se:.4f}); at "
              f"gamma* = 64: d_AI = {d_ai:.4f}, d_Gin = {d_gi:.4f}")
    if gap <= 2 * se:
        if d_ai <= 2 * se and d_gi >= 4 * se:
            verdict = ("SATURATED on AI-DAGGER: the reciprocal damped "
                       "plate selects the complex-symmetric class -- the "
                       "P7 endpoint call, completed.")
        elif d_gi <= 2 * se and d_ai >= 4 * se:
            verdict = "SATURATED on GINIBRE."
        else:
            verdict = (f"SATURATED BETWEEN (d_AI/se = {d_ai/se:.1f}, "
                       f"d_Gin/se = {d_gi/se:.1f}) -- the stationary "
                       f"point is not a canonical ensemble; physical-"
                       f"damping model question, registered.")
    else:
        verdict = (f"STILL-FLOWING (gap {gap/se:.1f} se): stronger "
                   f"overlap needed; trajectory reported.")
    md.append(f"\n**Reading: {verdict}**")
    results["marks"] = {f"g{g:g}": marks[g] for g in GSTAR}
    results["verdict"] = verdict
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS_E12C.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e12c.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
