#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3c -- free disk sector (RUNNER). See README.md (preregistered)."""
import json
import os
import time

import numpy as np

from platefem import (assemble_plate, classify_parity_resolved, n_star,
                      sector_basis, solve_lowest, solve_modes, split_rigid)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    nu=0.33, theta=2.0, R=1.0,
    # v2 (long ladder): decider for the ambiguous +1.7 sigma v1
    nrings=127, nrings_check=90,
    n_modes=1600,
    sigma=-10.0,
    spacing_frac=0.1,
    probe_npts=1500,
)


def sector_probes(basis, theta, R, npts, seed=17):
    """Interior points of the sector and their bisector-mirror images.
    The x-mirror is the identity trick: Pmx = P (cx = 1 for every mode), so
    classify_parity_resolved maps 'ee' -> S and 'eo' -> A."""
    rng = np.random.default_rng(seed)
    r = np.sqrt(rng.uniform(0.0, 0.94, npts)) * R
    th = rng.uniform(-0.49 * theta, 0.49 * theta, npts)
    pts = np.vstack([r * np.cos(th), r * np.sin(th)])
    P = basis.probes(pts)
    Pmy = basis.probes(np.vstack([pts[0], -pts[1]]))
    return P, P, Pmy


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg), "gates": {}, "runs": {}}

    print(f"[setup] sector theta={cfg['theta']}, nrings={cfg['nrings']} ...")
    mesh, basis = sector_basis(cfg["nrings"], cfg["theta"], cfg["R"])
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    print(f"[setup] {K.shape[0]} dofs ({time.time()-t00:.1f} s)")

    # G2: internal two-mesh check
    t0 = time.time()
    mesh_c, basis_c = sector_basis(cfg["nrings_check"], cfg["theta"], cfg["R"])
    Kc, Mc = assemble_plate(mesh_c, basis_c, cfg["nu"])
    lam_c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam_c, _, nr_c, rigmax_c = split_rigid(lam_c)
    del Kc, Mc, mesh_c, basis_c
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3, cfg["sigma"])
    lam, V, n_rigid, rigid_max = split_rigid(lam, V)
    ns_int = n_star(lam, lam_c, cfg["spacing_frac"])
    print(f"[G2] nrings {cfg['nrings']} (rigid {n_rigid}, max {rigid_max:.1e}) "
          f"vs {cfg['nrings_check']} (rigid {nr_c}): internal N* = "
          f"{ns_int}/{cfg['n_modes']} ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=ns_int, rigid=[n_rigid, nr_c],
                                  rigid_max=[rigid_max, rigmax_c],
                                  ndof=int(K.shape[0]))

    P, Pmx, Pmy = sector_probes(basis, cfg["theta"], cfg["R"], cfg["probe_npts"])
    labels, qual, n_res = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    lab_sa = ["S" if l == "ee" else ("A" if l == "eo" else "xx") for l in labels]
    counts = {s: lab_sa.count(s) for s in ["S", "A", "xx"]}
    print(f"[sector] counts {counts}, resolved {n_res}, min quality "
          f"{np.min(qual):.3f}, solve sweeps {sinfo['sweeps']}, "
          f"max_resid {sinfo['max_resid']:.1e}")
    results["runs"]["sector"] = dict(
        lam=np.asarray(lam)[:cfg["n_modes"]].tolist(),
        labels=lab_sa[:cfg["n_modes"]], counts=counts,
        min_quality=float(np.min(qual)), n_resolved=n_res,
        n_rigid=n_rigid,
        solve=dict(sweeps=sinfo["sweeps"], max_resid=sinfo["max_resid"]))

    results["n_use"] = int(min(cfg["n_modes"], ns_int))
    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s; strict N_use = "
          f"{results['n_use']}; now run analyze_sector.py")


if __name__ == "__main__":
    main()
