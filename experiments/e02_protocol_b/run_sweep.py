#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E2 -- Protocol B: boundary-controlled Poisson -> RP transition (SWEEP RUNNER)
=============================================================================
Paper prediction ("Boundary-controlled transition"): tuning all four edges
symmetrically with a Winkler restraint V_n + kappa w = 0 (M_n = 0 throughout)
should interpolate the per-sector spacing statistics continuously from Poisson
(<r> = 0.3863, simply supported, kappa -> inf) to the intermediate/RP value
(free, kappa = 0). Symmetric variation preserves Z2 x Z2, so the four parity
sectors remain well defined along the whole path.

This script computes and stores the raw data (eigenvalues + sector labels per
kappa); analyze.py turns them into <r> curves and the verdict.

Accuracy gates (run first, recorded in results):
  G1: kappa = 1e10 vs exact SSSS -> N*_SS      (absolute anchor at the SS end)
  G2: kappa = 0 on 96x60 vs 128x80 -> N*_int   (internal two-mesh estimate at
      the free end; difference-based)
  G3: kappa = 0 sector-resolved lowest modes vs the Legendre-Ritz per-sector
      spectra -> classifier validation (labels must reproduce each sector's
      spectrum, not just the union)

Modes used by the analysis: N_use = min(n_modes, N*_SS, N*_int).
"""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, boundary_matrix,
                      classify_parity_resolved, n_star, probe_operators,
                      rectangle_basis, ritz, solve_lowest, solve_modes,
                      split_rigid, ssss_exact, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0,
    b=1.0 / 1.6189043236,
    nu=0.33,
    mesh=(96, 60),
    mesh_check=(128, 80),
    n_modes=400,                       # elastic modes stored per kappa
    kappas=[0.0, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 1e8, 1e10],
    sigma=-10.0,
    spacing_frac=0.1,
    ritz_nleg=64, ritz_nleg_conv=72, ritz_tol=1e-4,
    probe_grid=(48, 30),
    n_ctrl=30,                         # per-sector modes for the G3 control
)


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    results = {"config": {k: (list(v) if isinstance(v, (list, tuple)) else v)
                          for k, v in cfg.items()},
               "gates": {}, "sweep": {}}

    # ---------------- assembly (once per mesh) ----------------
    nx, ny = cfg["mesh"]
    print(f"[setup] assembling {nx}x{ny} Argyris ...")
    mesh, basis = rectangle_basis(nx, ny, a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, nu)
    B = boundary_matrix(mesh, ElementTriArgyris)
    P, Pmx, Pmy = probe_operators(basis, a, b, *cfg["probe_grid"])
    print(f"[setup] {K.shape[0]} dofs ({time.time()-t00:.1f} s)")

    # ---------------- G1: SS anchor ----------------
    t0 = time.time()
    lam_ss = solve_lowest(K + 1e10 * B, M, cfg["n_modes"], cfg["sigma"])
    lam_ex = ssss_exact(a, b, cfg["n_modes"])
    nstar_ss = n_star(lam_ss, lam_ex, cfg["spacing_frac"])
    err_ss = float(np.max(np.abs(lam_ss - lam_ex) / lam_ex))
    print(f"[G1] kappa=1e10 vs exact SSSS: N* = {nstar_ss}/{cfg['n_modes']}, "
          f"max relerr {err_ss:.2e} ({time.time()-t0:.1f} s)")
    results["gates"]["G1"] = dict(n_star=nstar_ss, max_relerr=err_ss)

    # ---------------- G2: free-end internal check ----------------
    t0 = time.time()
    lam0, V0, solve_info0 = solve_modes(K, M, cfg["n_modes"] + 3, cfg["sigma"])
    lam0, V0, n_rigid, rigid_max = split_rigid(lam0, V0)
    print(f"[G2] kappa=0 96x60: rigid {n_rigid} (max {rigid_max:.1e}), "
          f"solve {solve_info0} ({time.time()-t0:.1f} s)")
    t0 = time.time()
    mesh_c, basis_c = rectangle_basis(*cfg["mesh_check"], a, b, ElementTriArgyris)
    Kc, Mc = assemble_plate(mesh_c, basis_c, nu)
    lam0c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam0c, _, n_rigid_c, rigid_max_c = split_rigid(lam0c)
    nstar_int = n_star(lam0, lam0c, cfg["spacing_frac"])
    print(f"[G2] vs 128x80 (rigid {n_rigid_c}, max {rigid_max_c:.1e}): "
          f"internal N* = {nstar_int}/{cfg['n_modes']} ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=nstar_int, rigid=[n_rigid, n_rigid_c],
                                  rigid_max=[rigid_max, rigid_max_c],
                                  solve=solve_info0)
    del Kc, Mc, mesh_c, basis_c

    # ---------------- G3: classifier control at kappa=0 ----------------
    t0 = time.time()
    labels0, qual0, n_res0 = classify_parity_resolved(lam0, V0, P, Pmx, Pmy, K, M)
    spectra = ritz.converged_spectrum(a, b, nu, nleg=cfg["ritz_nleg"],
                                      nleg_conv=cfg["ritz_nleg_conv"],
                                      tol=cfg["ritz_tol"])
    g3 = {}
    for sec in SECTORS:
        lam_sec = np.array([l for l, s in zip(lam0, labels0) if s == sec])
        nc = min(cfg["n_ctrl"], len(lam_sec), spectra[sec]["n_conv"])
        rel = np.abs(lam_sec[:nc] - spectra[sec]["elastic"][:nc]) / \
            spectra[sec]["elastic"][:nc]
        g3[sec] = dict(n=int(nc), max_reldiff=float(np.max(rel)))
    n_amb = labels0.count("xx")
    print(f"[G3] classifier vs Ritz sectors: "
          + ", ".join(f"{s}: {g3[s]['max_reldiff']:.1e} ({g3[s]['n']})" for s in SECTORS)
          + f"; resolved {n_res0}, ambiguous {n_amb}/{len(labels0)} "
          f"({time.time()-t0:.1f} s)")
    results["gates"]["G3"] = dict(per_sector=g3, n_ambiguous=n_amb,
                                  n_resolved=n_res0)

    n_use = int(min(cfg["n_modes"], nstar_ss, nstar_int))
    results["n_use"] = n_use
    print(f"[gates] N_use = {n_use}")
    save(results)

    # ---------------- sweep ----------------
    for kap in cfg["kappas"]:
        t0 = time.time()
        if kap == 0.0:
            lam, labels, qual = lam0, labels0, qual0
            n_rig, rig_max, sinfo, n_res = n_rigid, rigid_max, solve_info0, n_res0
        else:
            Kk = (K + kap * B).tocsc()
            lam, V, sinfo = solve_modes(Kk, M, cfg["n_modes"], cfg["sigma"])
            n_rig, rig_max = 0, 0.0
            labels, qual, n_res = classify_parity_resolved(lam, V, P, Pmx, Pmy,
                                                           Kk, M)
            del Kk, V
        counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
        key = f"{kap:.0e}" if kap else "0"
        results["sweep"][key] = dict(
            kappa=kap, lam=np.asarray(lam)[:cfg["n_modes"]].tolist(),
            labels=labels[:cfg["n_modes"]],
            min_quality=float(np.min(qual)), n_resolved=n_res,
            counts=counts, n_rigid=n_rig, rigid_max=rig_max, solve=sinfo,
            wall_s=round(time.time() - t0, 1))
        save(results)
        print(f"[sweep] kappa {kap:.0e}: counts {counts}, resolved {n_res}, "
              f"min quality {results['sweep'][key]['min_quality']:.3f} "
              f"({time.time()-t0:.1f} s)")

    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s -> results_raw.json; "
          f"now run analyze.py")


if __name__ == "__main__":
    main()
