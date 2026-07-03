#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3a -- free equilateral triangle (RUNNER). See README.md (preregistered)."""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, boundary_matrix,
                      classify_c3v, n_star, solve_lowest, solve_modes,
                      split_rigid, triangle_basis, triangle_probe_operators,
                      triangle_ss_exact)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    L=1.0, nu=0.33,
    refine=7, refine_check=6,
    n_modes=1000,
    kappa_ss=1e10,
    sigma=-10.0,
    spacing_frac=0.1,
    probe_npts=1500,
)
TRI_SECTORS = ["A1", "A2", "E"]


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg), "gates": {}, "runs": {}}

    print(f"[setup] assembling triangle refine={cfg['refine']} Argyris ...")
    mesh, basis, verts = triangle_basis(cfg["refine"], cfg["L"])
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    B = boundary_matrix(mesh, ElementTriArgyris)
    P, P_R, P_s = triangle_probe_operators(basis, cfg["L"], cfg["probe_npts"])
    print(f"[setup] {K.shape[0]} dofs ({time.time()-t00:.1f} s)")

    # ---------------- G1: SS anchor vs exact Lame^2 ----------------
    t0 = time.time()
    Kss = (K + cfg["kappa_ss"] * B).tocsc()
    lam_ss = solve_lowest(Kss, M, cfg["n_modes"], cfg["sigma"])
    lam_ex = triangle_ss_exact(cfg["L"], cfg["n_modes"])
    nstar_ss = n_star(lam_ss, lam_ex, cfg["spacing_frac"])
    err_ss = float(np.max(np.abs(lam_ss - lam_ex) / lam_ex))
    print(f"[G1] kappa=1e10 vs exact Lame^2: N* = {nstar_ss}/{cfg['n_modes']}, "
          f"max relerr {err_ss:.2e} ({time.time()-t0:.1f} s)")
    results["gates"]["G1"] = dict(n_star=nstar_ss, max_relerr=err_ss)
    save(results)

    # ---------------- G2: free-end internal check ----------------
    t0 = time.time()
    mesh_c, basis_c, _ = triangle_basis(cfg["refine_check"], cfg["L"])
    Kc, Mc = assemble_plate(mesh_c, basis_c, cfg["nu"])
    lam0c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam0c, _, n_rigid_c, rigid_max_c = split_rigid(lam0c)
    lam0_novec = solve_lowest(K, M, cfg["n_modes"] + 3, cfg["sigma"])
    lam0_novec, _, n_rigid0, rigid_max0 = split_rigid(lam0_novec)
    nstar_int = n_star(lam0_novec, lam0c, cfg["spacing_frac"])
    print(f"[G2] refine {cfg['refine']} (rigid {n_rigid0}, max {rigid_max0:.1e}) "
          f"vs {cfg['refine_check']} (rigid {n_rigid_c}): internal N* = "
          f"{nstar_int}/{cfg['n_modes']} ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=nstar_int, rigid=[n_rigid0, n_rigid_c],
                                  rigid_max=[rigid_max0, rigid_max_c])
    save(results)
    del Kc, Mc, mesh_c, basis_c

    n_use = int(min(cfg["n_modes"], nstar_ss, nstar_int))
    results["n_use"] = n_use
    print(f"[gates] N_use = {n_use}")

    # ---------------- runs: free (kappa=0) and SS (kappa=1e10) ----------------
    for tag, Kk, extra in [("free", K.tocsc(), 3), ("ss", Kss, 0)]:
        t0 = time.time()
        lam, V, sinfo, _ = solve_modes(Kk, M, cfg["n_modes"] + extra, cfg["sigma"])
        if extra:
            lam, V, n_rig, rig_max = split_rigid(lam, V)
        else:
            n_rig, rig_max = 0, 0.0
        labels, qual, n_res = classify_c3v(lam, V, P, P_R, P_s, Kk, M)
        counts = {s: labels.count(s) for s in TRI_SECTORS + ["xx"]}
        results["runs"][tag] = dict(
            kappa=0.0 if tag == "free" else cfg["kappa_ss"],
            lam=np.asarray(lam)[:cfg["n_modes"]].tolist(),
            labels=labels[:cfg["n_modes"]],
            min_quality=float(np.min(qual)), n_resolved=n_res,
            counts=counts, n_rigid=n_rig, rigid_max=rig_max,
            solve=dict(sweeps=sinfo["sweeps"], max_resid=sinfo["max_resid"]),
            wall_s=round(time.time() - t0, 1))
        save(results)
        print(f"[{tag}] counts {counts}, resolved {n_res}, min quality "
              f"{results['runs'][tag]['min_quality']:.3f}, rigid {n_rig} "
              f"(max {rig_max:.1e}) ({time.time()-t0:.1f} s)")

    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s -> results_raw.json; "
          f"now run analyze_triangle.py")


if __name__ == "__main__":
    main()
