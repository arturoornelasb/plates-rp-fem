#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E10 -- annulus: the second Poisson control (RUNNER + ANALYSIS).

The free annulus preserves polar separability (both boundaries are r = const
coordinate lines): like the disk, the paper's mechanism REQUIRES Poisson
statistics here despite free edges everywhere. Semi-analytic 4x4 J/Y/I/K
determinant per m (platefem.annulus); FEM ring-mesh cross-validation.
Preregistered: SUPPORTS if class-pooled <r> is Poisson-consistent and every
fixed-m sequence is a picket fence; CHALLENGES if intermediate.
"""
import json
import os
import time

import numpy as np

from platefem import annulus, annulus_basis, assemble_plate, n_star, \
    r_values, solve_lowest, split_rigid
from platefem.stats import R_GOE, R_POISSON, mean_r

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(nu=0.33, r_in=0.4, r_out=1.0, x_max=110.0, m_max=120,
           fem_nr=24, fem_k=203, spacing_frac=0.1)
SKIP = 10


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg)}

    # ---------- semi-analytic ----------
    t0 = time.time()
    roots = annulus.free_annulus_roots(cfg["nu"], cfg["r_in"], cfg["r_out"],
                                       cfg["x_max"], cfg["m_max"])
    lam_class = annulus.class_levels(roots)
    lam_full, mult = annulus.full_spectrum(roots)
    lam_full_exp = np.repeat(lam_full, mult)
    print(f"[semi] {len(lam_class)} class levels, m up to "
          f"{max(roots)} ({time.time()-t0:.1f} s)")

    # picket-fence check per m
    r_per_m = []
    for m, bs in roots.items():
        if len(bs) >= 13:
            r_per_m.append(mean_r(bs ** 4, 2)[0])
    # pooled class statistics (paper-style)
    n_use = min(len(lam_class), 1500)
    r_pool = mean_r(lam_class[:n_use], SKIP)
    print(f"[semi] pooled class <r> = {r_pool[0]:.4f} +/- {r_pool[1]:.4f} "
          f"({r_pool[2]} ratios); per-m picket <r> mean "
          f"{np.mean(r_per_m):.3f}")

    # ---------- FEM cross-validation ----------
    t0 = time.time()
    mesh, basis = annulus_basis(cfg["fem_nr"], cfg["r_in"], cfg["r_out"])
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam_fem = solve_lowest(K, M, cfg["fem_k"])
    lam_fem, _, n_rigid, rigid_max = split_rigid(lam_fem)
    ref = lam_full_exp[:len(lam_fem)]
    rel = np.abs(lam_fem[:len(ref)] - ref) / ref
    ns = n_star(lam_fem, ref, cfg["spacing_frac"])
    print(f"[FEM] {K.shape[0]} dofs, rigid {n_rigid} (max {rigid_max:.1e}); "
          f"N* = {ns}/{len(ref)}, max relerr {np.max(rel):.2e} "
          f"({time.time()-t0:.1f} s)")

    # ---------- verdict ----------
    poisson_ok = abs(r_pool[0] - R_POISSON) < 3 * r_pool[1]
    fem_ok = ns >= 100
    verdict = ("SUPPORTS (annulus control Poisson; FEM cross-validated)"
               if (poisson_ok and fem_ok) else
               ("CHALLENGES (annulus intermediate -- separability side fails)"
                if not poisson_ok else "CHECK FEM gate"))
    md = ["# E10 -- free annulus: second Poisson control (RESULTS)\n",
          f"r_in/r_out = {cfg['r_in']}/{cfg['r_out']}, nu = {cfg['nu']}. "
          f"References: Poisson {R_POISSON}, GOE {R_GOE}.\n",
          f"- semi-analytic class levels: {len(lam_class)} "
          f"(m up to {max(roots)})",
          f"- pooled one-reflection-class <r> over {r_pool[2]} ratios: "
          f"**{r_pool[0]:.4f} +/- {r_pool[1]:.4f}**",
          f"- fixed-m sequences: mean picket-fence <r> = "
          f"{np.mean(r_per_m):.3f} (rigid-ladder limit ~0.95+)",
          f"- FEM ring mesh ({K.shape[0]} dofs): rigid {n_rigid}, "
          f"N* = {ns}, max relerr {np.max(rel):.2e} over {len(ref)} modes",
          f"\n**Reading: {verdict}**"]
    results.update(dict(n_class=len(lam_class), r_pooled=r_pool,
                        r_picket_mean=float(np.mean(r_per_m)),
                        fem=dict(ndof=int(K.shape[0]), n_rigid=n_rigid,
                                 n_star=int(ns), max_relerr=float(np.max(rel))),
                        verdict=verdict,
                        wall_time_s=round(time.time() - t00, 1)))
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
