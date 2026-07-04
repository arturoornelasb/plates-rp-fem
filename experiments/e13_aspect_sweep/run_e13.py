#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E13 -- aspect-ratio robustness sweep (registered follow-up; RUNNER+ANALYSIS).

The whole campaign ran at a/b = 1.6189. This sweep closes that caveat:
FFFF per-sector <r> at five generic irrational aspect ratios, identical
protocol. Preregistered reading: <r> stable across ratios (the mechanism is
generic) vs ratio-dependent structure (would flag a special-geometry
artifact). Square-like ratios are avoided (exact degeneracies).
"""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, mean_r,
                      probe_operators, r_values, rectangle_basis, solve_modes,
                      split_rigid, SECTORS)
from platefem.stats import classify_parity_resolved

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    nu=0.33, n_modes=800,
    ratios=[1.2724, 1.4142135624, 1.6189043236, 1.8329988,
            2.0322198],
    ny=60,
)
SKIP = 10


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg), "runs": {}}
    rows = []
    for ratio in cfg["ratios"]:
        t0 = time.time()
        a, b = 1.0, 1.0 / ratio
        nx = int(round(cfg["ny"] * ratio))
        mesh, basis = rectangle_basis(nx, cfg["ny"], a, b, ElementTriArgyris)
        K, M = assemble_plate(mesh, basis, cfg["nu"])
        lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3)
        lam, V, n_rigid, _ = split_rigid(lam, V)
        P, Pmx, Pmy = probe_operators(basis, a, b)
        labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
        lam = lam[:cfg["n_modes"]]
        labels = labels[:cfg["n_modes"]]
        per, rv = {}, []
        for s in SECTORS:
            ev = np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
            per[s] = mean_r(ev, SKIP)
            rv.extend(r_values(ev[SKIP:]).tolist())
        rv = np.array(rv)
        pool = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))
        rows.append(dict(ratio=ratio, per_sector=per, pooled=pool,
                         mesh=[nx, cfg["ny"]], n_rigid=n_rigid,
                         min_quality=float(np.min(qual)),
                         wall_s=round(time.time() - t0, 1)))
        results["runs"][f"r{ratio:g}"] = rows[-1]
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)
        print(f"[a/b={ratio:g}] pooled <r> = {pool[0]:.4f} +/- {pool[1]:.4f} "
              f"({time.time()-t0:.1f} s)")

    md = ["# E13 -- aspect-ratio robustness sweep (RESULTS)\n",
          f"FFFF, nu = {cfg['nu']}, {cfg['n_modes']} modes per ratio, "
          f"identical protocol.\n"]
    md.append("| a/b | " + " | ".join(SECTORS) + " | pooled |")
    md.append("|---|---|---|---|---|---|")
    for r in rows:
        md.append(f"| {r['ratio']:.4f} | "
                  + " | ".join(f"{r['per_sector'][s][0]:.3f}" for s in SECTORS)
                  + f" | **{r['pooled'][0]:.4f} +/- {r['pooled'][1]:.4f}** |")
    pools = [r["pooled"] for r in rows]
    spread = max(p[0] for p in pools) - min(p[0] for p in pools)
    typ = float(np.mean([p[1] for p in pools]))
    md.append(f"\n- spread across ratios: {spread:.4f} (typical sem {typ:.4f})")
    verdict = ("ROBUST (mechanism generic across aspect ratios)"
               if spread < 3 * typ else
               "RATIO-DEPENDENT structure -- inspect per-sector table")
    md.append(f"\n**Reading: {verdict}**")
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
