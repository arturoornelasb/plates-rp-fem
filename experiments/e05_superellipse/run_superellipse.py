#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E5 -- superellipse corner sweep (RUNNER). See README.md (preregistered).
p = 2 is reused from E3b v2 at analysis time; here p in {3, 4, 6, 10}."""
import json
import os
import time

import numpy as np

from platefem import (assemble_plate, centered_probe_operators,
                      classify_parity_resolved, n_star, solve_lowest,
                      solve_modes, split_rigid, superellipse_basis, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    nu=0.33,
    ratio=1.6189043236,
    refine=6, refine_check=5,
    n_modes=1600,
    ps=[3.0, 4.0, 6.0, 10.0],
    sigma=-10.0, spacing_frac=0.1,
    probe_npts=1500,
)


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def mesh_quality(mesh):
    p, t = mesh.p, mesh.t
    a, b, c = p[:, t[0]], p[:, t[1]], p[:, t[2]]
    ar = 0.5 * np.abs((b[0] - a[0]) * (c[1] - a[1])
                      - (c[0] - a[0]) * (b[1] - a[1]))
    return float(ar.min() / ar.mean())


def main():
    t00 = time.time()
    cfg = CFG
    a_ax = np.sqrt(cfg["ratio"])
    b_ax = 1.0 / a_ax
    results = {"config": dict(cfg), "gates": {}, "runs": {}}

    # G2 at the sharpest p (worst corner elements)
    p_worst = cfg["ps"][-1]
    t0 = time.time()
    mesh_c, basis_c = superellipse_basis(cfg["refine_check"], a_ax, b_ax, p_worst)
    Kc, Mc = assemble_plate(mesh_c, basis_c, cfg["nu"])
    lam_c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam_c, _, _, _ = split_rigid(lam_c)
    del Kc, Mc, mesh_c, basis_c
    print(f"[G2-prep] p={p_worst} refine-{cfg['refine_check']} reference "
          f"({time.time()-t0:.1f} s)")

    for p_exp in cfg["ps"]:
        t0 = time.time()
        mesh, basis = superellipse_basis(cfg["refine"], a_ax, b_ax, p_exp)
        q = mesh_quality(mesh)
        K, M = assemble_plate(mesh, basis, cfg["nu"])
        lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3, cfg["sigma"])
        lam, V, n_rigid, rigid_max = split_rigid(lam, V)
        P, Pmx, Pmy = centered_probe_operators(basis, a_ax, b_ax,
                                               cfg["probe_npts"])
        labels, qual, n_res = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
        counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
        rec = dict(p=p_exp, mesh_quality=q,
                   lam=np.asarray(lam)[:cfg["n_modes"]].tolist(),
                   labels=labels[:cfg["n_modes"]], counts=counts,
                   min_quality=float(np.min(qual)), n_resolved=n_res,
                   n_rigid=n_rigid, rigid_max=rigid_max,
                   solve=dict(sweeps=sinfo["sweeps"],
                              max_resid=sinfo["max_resid"]),
                   wall_s=round(time.time() - t0, 1))
        if p_exp == p_worst:
            ns = n_star(lam, lam_c, cfg["spacing_frac"])
            results["gates"]["G2"] = dict(p=p_exp, n_star=ns)
            print(f"[G2] p={p_exp}: internal N* = {ns}/{cfg['n_modes']}")
        results["runs"][f"p{p_exp:g}"] = rec
        save(results)
        print(f"[p={p_exp:g}] quality {q:.3f}, counts {counts}, rigid {n_rigid} "
              f"(max {rigid_max:.1e}), resolved {n_res}, min qual "
              f"{np.min(qual):.3f} ({time.time()-t0:.1f} s)")
        del K, M, V, mesh, basis

    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s -> results_raw.json; "
          f"now run analyze_superellipse.py")


if __name__ == "__main__":
    main()
