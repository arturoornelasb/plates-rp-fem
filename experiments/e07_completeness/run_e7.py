#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E7 -- completeness tests: Protocol A mixed config + material independence.

(a) Protocol A config (ii) (paper, 'Boundary-controlled transition'):
    x-edges free, y-edges simply supported -- both reflection symmetries
    preserved. Realized with the validated edge-selective Winkler penalty
    (kappa = 1e10 on the y-edges only). Expectation under the campaign
    synthesis: intermediate statistics persist (the free x-edges still break
    the product structure), between SSSS and FFFF.

(b) Material independence (paper Prediction 'material'): per-sector <r> of
    the FFFF plate for nu in {0.25, 0.30, 0.33, 0.36} -- <r> should be
    insensitive to nu (the operator family stays orthogonal-class and the
    coupling structure is geometric).
"""
import json
import os
import time

import numpy as np

from skfem import FacetBasis
from platefem import (ElementTriArgyris, assemble_plate, boundary_matrix,
                      classify_parity_resolved, mean_r, probe_operators,
                      r_values, rectangle_basis, solve_modes, split_rigid,
                      SECTORS)
from platefem.kirchhoff import boundary_mass

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236,
    mesh=(96, 60), n_modes=800,
    nus=[0.25, 0.30, 0.33, 0.36],
    nu_protA=0.33, kappa=1e10,
)
SKIP = 10


def pooled_sector_r(lam, labels):
    rv = []
    per = {}
    for s in SECTORS:
        ev = np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
        per[s] = mean_r(ev, SKIP)
        rv.extend(r_values(ev[SKIP:]).tolist())
    rv = np.array(rv)
    return per, (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))


def main():
    t00 = time.time()
    cfg = CFG
    a, b = cfg["a"], cfg["b"]
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}, "runs": {}}

    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    P, Pmx, Pmy = probe_operators(basis, a, b)

    # ---------------- (a) Protocol A mixed: y-edges SS ----------------
    t0 = time.time()
    K, M = assemble_plate(mesh, basis, cfg["nu_protA"])
    fac_y = mesh.facets_satisfying(
        lambda x: (x[1] < 1e-10) | (x[1] > b - 1e-10))
    fb = FacetBasis(mesh, ElementTriArgyris(), facets=fac_y)
    By = boundary_mass.assemble(fb).tocsc()
    Kk = (K + cfg["kappa"] * By).tocsc()
    lam, V, sinfo, _ = solve_modes(Kk, M, cfg["n_modes"] + 1)
    # mixed config: 1 rigid-like mode remains? (x-translation is blocked by
    # w=0 on y-edges; no rigid modes survive) -- detect by gap as usual
    lam, V, n_rigid, _ = split_rigid(lam, V)
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, Kk, M)
    lam = lam[:cfg["n_modes"]]
    labels = labels[:cfg["n_modes"]]
    per, pool = pooled_sector_r(lam, labels)
    results["runs"]["protA_mixed"] = dict(
        n_rigid=n_rigid, per_sector=per, pooled=pool,
        counts={s: labels.count(s) for s in SECTORS + ["xx"]},
        min_quality=float(np.min(qual)), wall_s=round(time.time() - t0, 1))
    print(f"[protA] FSFS-like (free x, SS y): pooled <r> = {pool[0]:.4f} +/- "
          f"{pool[1]:.4f}, rigid {n_rigid} ({time.time()-t0:.1f} s)")

    # ---------------- (b) nu sweep, FFFF ----------------
    for nu in cfg["nus"]:
        t0 = time.time()
        K, M = assemble_plate(mesh, basis, nu)
        lam, V, sinfo, _ = solve_modes(K.tocsc(), M, cfg["n_modes"] + 3)
        lam, V, n_rigid, _ = split_rigid(lam, V)
        labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
        lam = lam[:cfg["n_modes"]]
        labels = labels[:cfg["n_modes"]]
        per, pool = pooled_sector_r(lam, labels)
        results["runs"][f"nu{nu:g}"] = dict(
            per_sector=per, pooled=pool, n_rigid=n_rigid,
            min_quality=float(np.min(qual)), wall_s=round(time.time() - t0, 1))
        print(f"[nu={nu:g}] pooled <r> = {pool[0]:.4f} +/- {pool[1]:.4f} "
              f"({time.time()-t0:.1f} s)")

    # ---------------- markdown ----------------
    md = ["# E7 -- completeness: Protocol A mixed config + material independence\n"]
    pa = results["runs"]["protA_mixed"]
    md.append(f"## (a) Protocol A config (ii): x-edges free, y-edges SS\n")
    md.append(f"- pooled <r> = **{pa['pooled'][0]:.4f} +/- {pa['pooled'][1]:.4f}** "
              f"(FFFF at this ladder: ~0.442; SSSS baseline: ~0.391); per sector "
              + ", ".join(f"{s} {pa['per_sector'][s][0]:.3f}" for s in SECTORS))
    md.append(f"- rigid modes detected: {pa['n_rigid']} (expected 0: w = 0 on "
              f"the y-edges blocks all rigid motions)")
    md.append("\n## (b) FFFF <r>(nu) -- material independence\n")
    md.append("| nu | " + " | ".join(SECTORS) + " | pooled |")
    md.append("|---|---|---|---|---|---|")
    for nu in cfg["nus"]:
        r = results["runs"][f"nu{nu:g}"]
        md.append(f"| {nu:g} | "
                  + " | ".join(f"{r['per_sector'][s][0]:.3f}" for s in SECTORS)
                  + f" | **{r['pooled'][0]:.4f} +/- {r['pooled'][1]:.4f}** |")
    pools = [results["runs"][f"nu{nu:g}"]["pooled"] for nu in cfg["nus"]]
    spread = max(p[0] for p in pools) - min(p[0] for p in pools)
    typ_se = float(np.mean([p[1] for p in pools]))
    md.append(f"\n- spread across nu: {spread:.4f} (typical sem {typ_se:.4f}) "
              f"-> {'consistent with material independence' if spread < 3 * typ_se else 'nu-DEPENDENT (check)'}")
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
