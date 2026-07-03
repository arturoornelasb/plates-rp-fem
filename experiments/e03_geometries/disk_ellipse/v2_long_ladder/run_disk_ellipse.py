#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3b -- free disk (control) vs free ellipse (RUNNER). See README.md."""
import json
import os
import time

import numpy as np

from platefem import (assemble_plate, centered_probe_operators,
                      classify_parity_resolved, disk, disk_basis,
                      ellipse_basis, n_star, solve_lowest, solve_modes,
                      split_rigid, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    # v2 (registered decider for the ambiguous v1): ladder doubled to 1600 at
    # the same refine-6 mesh. Validity rests on the OPERATIONAL disk gate over
    # the full extended ladder (strict per-eigenvalue N* is geometry-limited,
    # but the error is smooth; v1 measured 0.0-sigma statistics agreement).
    nu=0.33,
    refine=6, refine_check=5,
    n_modes=1600,
    ellipse_ratio=1.6189043236,          # a/b; a*b = 1 (area pi, disk-matched)
    sigma=-10.0,
    spacing_frac=0.1,
    disk_x_max=145.0, disk_m_max=150,
    probe_npts=1500,
)


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg), "gates": {}, "runs": {}}

    # ---------------- semi-analytic disk reference ----------------
    t0 = time.time()
    x2 = disk.validate(cfg["nu"])
    roots = disk.free_disk_roots(cfg["nu"], cfg["disk_x_max"], cfg["disk_m_max"])
    lam_ref, mult = disk.full_spectrum(roots)
    lam_full = disk.expand_multiplicity(lam_ref, mult)
    cls = disk.class_levels(roots)
    lam_class = np.array([r[0] for r in cls])
    print(f"[disk] fundamental x^2 = {x2:.4f}; {len(lam_ref)} distinct levels, "
          f"{len(lam_full)} with multiplicity, {len(lam_class)} class levels "
          f"({time.time()-t0:.1f} s)")
    results["disk_semianalytic"] = dict(
        fundamental_x2=float(x2), n_distinct=len(lam_ref),
        n_class=len(lam_class), lam_class=lam_class.tolist(),
        lam_full=lam_full[:3 * cfg["n_modes"]].tolist(),
        mult=mult.tolist(), lam_distinct=lam_ref.tolist())
    save(results)

    # ---------------- G1: FEM disk vs semi-analytic ----------------
    fem_disk = {}
    for r in [cfg["refine_check"], cfg["refine"]]:
        t0 = time.time()
        mesh, basis = disk_basis(r)
        K, M = assemble_plate(mesh, basis, cfg["nu"])
        lam = solve_lowest(K, M, cfg["n_modes"] + 3, cfg["sigma"])
        lam, _, n_rigid, rigid_max = split_rigid(lam)
        ref = lam_full[:len(lam)]
        ns = n_star(lam, ref, cfg["spacing_frac"])
        rel = np.abs(lam[:len(ref)] - ref) / ref
        fem_disk[r] = lam
        print(f"[G1] FEM disk refine {r} ({K.shape[0]} dofs, rigid {n_rigid}, "
              f"max {rigid_max:.1e}): strict N* = {ns}/{cfg['n_modes']}, "
              f"max relerr {np.max(rel):.2e}, median {np.median(rel):.2e} "
              f"({time.time()-t0:.1f} s)")
        results["gates"][f"G1_refine{r}"] = dict(
            n_star=ns, max_relerr=float(np.max(rel)),
            med_relerr=float(np.median(rel)), n_rigid=n_rigid,
            rigid_max=rigid_max, ndof=int(K.shape[0]))
        results["runs"][f"disk_fem_r{r}"] = dict(lam=lam[:cfg["n_modes"]].tolist())
        save(results)
        del K, M, mesh, basis

    # ---------------- ellipse ----------------
    a_ax = np.sqrt(cfg["ellipse_ratio"])
    b_ax = 1.0 / a_ax
    t0 = time.time()
    mesh, basis = ellipse_basis(cfg["refine"], a_ax, b_ax)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    print(f"[ellipse] refine {cfg['refine']}: {K.shape[0]} dofs "
          f"({time.time()-t0:.1f} s)")

    t0 = time.time()
    mesh_c, basis_c = ellipse_basis(cfg["refine_check"], a_ax, b_ax)
    Kc, Mc = assemble_plate(mesh_c, basis_c, cfg["nu"])
    lam_c = solve_lowest(Kc, Mc, cfg["n_modes"] + 3, cfg["sigma"])
    lam_c, _, nr_c, _ = split_rigid(lam_c)
    del Kc, Mc, mesh_c, basis_c

    lam_e, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3, cfg["sigma"])
    lam_e, V, n_rigid, rigid_max = split_rigid(lam_e, V)
    ns_int = n_star(lam_e, lam_c, cfg["spacing_frac"])
    print(f"[G2] ellipse refine {cfg['refine']} vs {cfg['refine_check']}: "
          f"internal N* = {ns_int}/{cfg['n_modes']}, rigid {n_rigid} "
          f"(max {rigid_max:.1e}) ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=ns_int, rigid=[n_rigid, nr_c])

    P, Pmx, Pmy = centered_probe_operators(basis, a_ax, b_ax, cfg["probe_npts"])
    labels, qual, n_res = classify_parity_resolved(lam_e, V, P, Pmx, Pmy, K, M)
    counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
    print(f"[ellipse] counts {counts}, resolved {n_res}, min quality "
          f"{np.min(qual):.3f}, solve sweeps {sinfo['sweeps']}, "
          f"max_resid {sinfo['max_resid']:.1e}")
    results["runs"]["ellipse"] = dict(
        a=a_ax, b=b_ax, lam=np.asarray(lam_e)[:cfg["n_modes"]].tolist(),
        labels=labels[:cfg["n_modes"]], counts=counts,
        min_quality=float(np.min(qual)), n_resolved=n_res,
        solve=dict(sweeps=sinfo["sweeps"], max_resid=sinfo["max_resid"]))

    results["n_use"] = int(min(cfg["n_modes"], ns_int))
    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s; N_use = {results['n_use']}; "
          f"now run analyze_disk_ellipse.py")


if __name__ == "__main__":
    main()
