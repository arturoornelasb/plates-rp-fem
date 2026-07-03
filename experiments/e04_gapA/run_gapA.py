#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E4 -- Gap A on the true operator (RUNNER). See README.md (preregistered)."""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, boundary_matrix,
                      classify_parity_resolved, n_star, probe_operators,
                      rectangle_basis, solve_lowest, solve_modes, split_rigid,
                      ssss_exact, SECTORS)
from platefem import bases

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(128, 80), mesh_check=(96, 60),
    n_modes=1200,
    sigma=-10.0, spacing_frac=0.1,
    resid_vec=1e-5,                       # tightened: coefficients are consumed
    n_funcs_axis=48,                      # 1D family size per axis
    grid=(192, 120),                      # Gauss grid for projections
    ladder=[64, 100, 144, 200, 256, 324],
)


def save(obj, name):
    with open(os.path.join(HERE, name), "w") as f:
        json.dump(obj, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    results = {"config": {k: (list(v) if isinstance(v, (list, tuple)) else v)
                          for k, v in cfg.items()}, "gates": {}}

    nx, ny = cfg["mesh"]
    print(f"[setup] assembling {nx}x{ny} Argyris ...")
    mesh, basis = rectangle_basis(nx, ny, a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, nu)
    B = boundary_matrix(mesh, ElementTriArgyris)
    print(f"[setup] {K.shape[0]} dofs ({time.time()-t00:.1f} s)")

    # G1: exact anchor on the production mesh
    t0 = time.time()
    lam_ss = solve_lowest((K + 1e10 * B).tocsc(), M, cfg["n_modes"], cfg["sigma"])
    lam_ex = ssss_exact(a, b, cfg["n_modes"])
    ns_ss = n_star(lam_ss, lam_ex, cfg["spacing_frac"])
    print(f"[G1] kappa=1e10 vs exact SSSS: N* = {ns_ss}/{cfg['n_modes']}, "
          f"max relerr {np.max(np.abs(lam_ss-lam_ex)/lam_ex):.2e} "
          f"({time.time()-t0:.1f} s)")
    results["gates"]["G1"] = dict(n_star=ns_ss)

    # G2: internal
    t0 = time.time()
    mesh_c, basis_c = rectangle_basis(*cfg["mesh_check"], a, b, ElementTriArgyris)
    Kc, Mc = assemble_plate(mesh_c, basis_c, nu)
    lam_c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam_c, _, _, _ = split_rigid(lam_c)
    del Kc, Mc, mesh_c, basis_c

    # certified eigenpairs, tightened residuals
    lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3, cfg["sigma"],
                                   resid_sanity=cfg["resid_vec"], sweeps_max=40)
    lam, V, n_rigid, rigid_max = split_rigid(lam, V)
    ns_int = n_star(lam, lam_c, cfg["spacing_frac"])
    print(f"[G2] internal N* = {ns_int}/{cfg['n_modes']}; rigid {n_rigid} "
          f"(max {rigid_max:.1e}); solve sweeps {sinfo['sweeps']}, "
          f"max_resid {sinfo['max_resid']:.1e} ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=ns_int, rigid=n_rigid,
                                  solve=dict(sweeps=sinfo["sweeps"],
                                             max_resid=sinfo["max_resid"]))
    n_use = int(min(cfg["n_modes"], ns_ss, ns_int))
    results["n_use"] = n_use
    print(f"[gates] N_use = {n_use}")

    # sector labels
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, n_res = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
    print(f"[labels] {counts}, resolved {n_res}, min quality {np.min(qual):.3f}")
    results["labels"] = dict(counts=counts, min_quality=float(np.min(qual)))

    # projections onto the registered bases
    t0 = time.time()
    gx, wx, gy, wy = bases.gauss_grid(a, b, *cfg["grid"])
    Xg, Yg = np.meshgrid(gx, gy, indexing="ij")
    Pgrid = basis.probes(np.vstack([Xg.ravel(), Yg.ravel()]))
    Wgrid = Pgrid @ V[:, :n_use]
    print(f"[proj] grid evaluation ({time.time()-t0:.1f} s)")

    fam = {}
    Sx, kSx, pSx = bases.sine_1d(cfg["n_funcs_axis"], a, gx)
    Sy, kSy, pSy = bases.sine_1d(cfg["n_funcs_axis"], b, gy)
    Bx, kBx, pBx = bases.beam_1d(cfg["n_funcs_axis"], a, gx, wx)
    By, kBy, pBy = bases.beam_1d(cfg["n_funcs_axis"], b, gy, wy)
    t0 = time.time()
    C_sine = bases.project_modes(Wgrid, Sx, wx, Sy, wy, *cfg["grid"])
    C_beam = bases.project_modes(Wgrid, Bx, wx, By, wy, *cfg["grid"])
    print(f"[proj] coefficients ({time.time()-t0:.1f} s)")
    fam["sine"] = (C_sine, kSx, pSx, kSy, pSy)
    fam["beam"] = (C_beam, kBx, pBx, kBy, pBy)

    # per-sector ladders + spectra
    out = {"sectors": {}}
    lam_use = lam[:n_use]
    for s in SECTORS:
        idx = [i for i, l in enumerate(labels[:n_use]) if l == s]
        out["sectors"][s] = {"lam": lam_use[idx].tolist(), "n": len(idx)}
        for name, (C, kx, px, ky, py) in fam.items():
            pi, pj = bases.sector_product_order(kx, px, ky, py, s,
                                                max(cfg["ladder"]))
            rows = bases.ipr_ladder(C[idx], pi, pj, cfg["ladder"])
            out["sectors"][s][name] = rows
        print(f"[gapA] sector {s}: {len(idx)} modes, ladders done")
    out["goe"] = bases.goe_ipr_baseline(cfg["ladder"])
    results["gapA"] = out
    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results, "results_raw.json")
    print(f"\n[done] total {results['wall_time_s']} s -> results_raw.json; "
          f"now run analyze_gapA.py")


if __name__ == "__main__":
    main()
