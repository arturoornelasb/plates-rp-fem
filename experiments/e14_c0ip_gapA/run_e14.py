#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E14 -- Gap A reconciliation via C0-IP true-operator windows (RUNNER).
See README.md (preregistered)."""
import json
import os
import time

import numpy as np
from skfem import MeshTri

from platefem import (ElementTriArgyris, assemble_plate, n_star,
                      rectangle_basis, solve_lowest, solve_modes, split_rigid,
                      ssss_exact, SECTORS)
from platefem.c0ip import assemble_c0ip, boundary_dofs
from platefem.stats import classify_parity_resolved
from platefem import bases

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(144, 89), k_fem=4, sigma_factor=10.0,
    n_modes=2400,
    gate_argyris_mesh=(128, 80), gate_n=1200,
    spacing_frac=0.1,
    n_funcs_axis=64, grid=(256, 160),
    ladder=[256, 512],
    probe_npts=1500,
)


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}, "gates": {}}

    mesh = MeshTri.init_tensor(np.linspace(0, a, cfg["mesh"][0] + 1),
                               np.linspace(0, b, cfg["mesh"][1] + 1))
    t0 = time.time()
    K, M, space = assemble_c0ip(mesh, k=cfg["k_fem"], nu=nu,
                                sigma_factor=cfg["sigma_factor"])
    K = 0.5 * (K + K.T)
    print(f"[setup] C0-IP P{cfg['k_fem']} {space.N} dofs "
          f"({time.time()-t0:.1f} s)")

    # ---------------- G1: SSSS vs exact at production mesh ----------------
    t0 = time.time()
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    lam_ss = solve_lowest(K[I][:, I].tocsc(), M[I][:, I].tocsc(),
                          cfg["n_modes"])
    lam_ex = ssss_exact(a, b, cfg["n_modes"])
    ns_ss = n_star(lam_ss, lam_ex, cfg["spacing_frac"])
    rel = np.abs(lam_ss - lam_ex) / lam_ex
    print(f"[G1] SSSS vs exact: N* = {ns_ss}/{cfg['n_modes']}, max relerr "
          f"{np.max(rel):.2e}, median {np.median(rel):.2e} "
          f"({time.time()-t0:.1f} s)")
    results["gates"]["G1"] = dict(n_star=int(ns_ss),
                                  max_relerr=float(np.max(rel)),
                                  med_relerr=float(np.median(rel)))
    save(results)

    # ---------------- G2: FFFF vs certified Argyris ----------------
    t0 = time.time()
    mk, bk = rectangle_basis(*cfg["gate_argyris_mesh"], a, b, ElementTriArgyris)
    Ka, Ma = assemble_plate(mk, bk, nu)
    lam_arg = solve_lowest(Ka, Ma, cfg["gate_n"] + 3)
    lam_arg, _, _, _ = split_rigid(lam_arg)
    del Ka, Ma, mk, bk
    lam_ff = solve_lowest(K.tocsc(), M.tocsc(), cfg["gate_n"] + 3)
    lam_ff, _, nr, rig = split_rigid(lam_ff)
    n_cmp = min(len(lam_ff), len(lam_arg))
    ns_ff = n_star(lam_ff[:n_cmp], lam_arg[:n_cmp], cfg["spacing_frac"])
    rel = np.abs(lam_ff[:n_cmp] - lam_arg[:n_cmp]) / lam_arg[:n_cmp]
    print(f"[G2] FFFF vs Argyris({cfg['gate_argyris_mesh']}): rigid {nr} "
          f"(max {rig:.1e}); N* = {ns_ff}/{n_cmp}, max relerr "
          f"{np.max(rel):.2e} ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=int(ns_ff), rigid=nr,
                                  max_relerr=float(np.max(rel)))
    save(results)

    # ---------------- production: certified eigenpairs ----------------
    t0 = time.time()
    lam, V, sinfo, _ = solve_modes(K.tocsc(), M.tocsc(), cfg["n_modes"] + 3,
                                   resid_sanity=1e-4, sweeps_max=30)
    lam, V, n_rigid, rigid_max = split_rigid(lam, V)
    print(f"[solve] {len(lam)} elastic, rigid {n_rigid} (max {rigid_max:.1e}), "
          f"resid {sinfo['max_resid']:.1e} ({time.time()-t0:.1f} s)")
    results["solve"] = dict(resid=sinfo["max_resid"], rigid=n_rigid)

    # classification
    t0 = time.time()
    rng = np.random.default_rng(13)
    pts = np.vstack([rng.uniform(0.02 * a, 0.98 * a, cfg["probe_npts"]),
                     rng.uniform(0.02 * b, 0.98 * b, cfg["probe_npts"])])
    P = space.probes(pts)
    Pmx = space.probes(np.vstack([a - pts[0], pts[1]]))
    Pmy = space.probes(np.vstack([pts[0], b - pts[1]]))
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
    print(f"[labels] {counts}, min quality {np.min(qual):.3f} "
          f"({time.time()-t0:.1f} s)")
    results["labels"] = dict(counts=counts, min_quality=float(np.min(qual)))
    save(results)

    # ---------------- projections + ladders ----------------
    t0 = time.time()
    gx, wx, gy, wy = bases.gauss_grid(a, b, *cfg["grid"])
    Xg, Yg = np.meshgrid(gx, gy, indexing="ij")
    Pgrid = space.probes(np.vstack([Xg.ravel(), Yg.ravel()]))
    Wgrid = Pgrid @ V
    Sx, kSx, pSx = bases.sine_1d(cfg["n_funcs_axis"], a, gx)
    Sy, kSy, pSy = bases.sine_1d(cfg["n_funcs_axis"], b, gy)
    Bx, kBx, pBx = bases.beam_1d(cfg["n_funcs_axis"], a, gx, wx)
    By, kBy, pBy = bases.beam_1d(cfg["n_funcs_axis"], b, gy, wy)
    C_sine = bases.project_modes(Wgrid, Sx, wx, Sy, wy, *cfg["grid"])
    C_beam = bases.project_modes(Wgrid, Bx, wx, By, wy, *cfg["grid"])
    print(f"[proj] done ({time.time()-t0:.1f} s)")

    out = {"sectors": {}}
    for s in SECTORS:
        idx = [i for i, l in enumerate(labels) if l == s]
        rec = {"n": len(idx)}
        for name, (C, kx, px, ky, py) in [("sine", (C_sine, kSx, pSx, kSy, pSy)),
                                          ("beam", (C_beam, kBx, pBx, kBy, pBy))]:
            pi, pj = bases.sector_product_order(kx, px, ky, py, s,
                                                max(cfg["ladder"]))
            rec[name] = bases.ipr_ladder(C[idx], pi, pj, cfg["ladder"])
        out["sectors"][s] = rec
        print(f"[gapA] {s}: {len(idx)} modes, ladders done")
    out["goe"] = bases.goe_ipr_baseline(cfg["ladder"])
    results["gapA"] = out
    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)

    # ---------------- decisive comparison ----------------
    md = ["# E14 -- Gap A reconciliation (RESULTS)\n",
          f"C0-IP P{cfg['k_fem']} at {space.N} dofs, {cfg['n_modes']} certified "
          f"modes (~{cfg['n_modes']//4}/sector). E9 truncated-ladder reference "
          f"at N = 512: sine 0.0946, beam 0.1784 (falling); E4 true-operator "
          f"at N <= 324: sine ~0.17-0.19 (flat).\n"]
    md.append(f"- G1 (SSSS vs exact): N* = {results['gates']['G1']['n_star']}, "
              f"max relerr {results['gates']['G1']['max_relerr']:.2e}")
    md.append(f"- G2 (FFFF vs certified Argyris): N* = "
              f"{results['gates']['G2']['n_star']}, max relerr "
              f"{results['gates']['G2']['max_relerr']:.2e}\n")
    md.append("| sector | basis | " + " | ".join(f"IPR N={N}" for N in
              cfg["ladder"]) + " | E9 truncated at N=512 |")
    md.append("|" + "---|" * (len(cfg["ladder"]) + 3))
    e9ref = {"sine": 0.0946, "beam": 0.1784}
    verdict_rows = []
    for s in SECTORS:
        for name in ("sine", "beam"):
            rows = out["sectors"][s][name]
            iprs = [np.exp(r["mlnipr"]) for r in rows]
            if len(iprs) >= 2:
                verdict_rows.append((name, iprs[-1]))
            md.append(f"| {s} | {name} | "
                      + " | ".join(f"{v:.4f}" for v in iprs)
                      + f" | {e9ref[name]:.4f} |")
    sine512 = [v for nm, v in verdict_rows if nm == "sine"]
    mean512 = float(np.mean(sine512)) if sine512 else np.nan
    md.append(f"\n- true-operator sine IPR at N = 512 (sector mean): "
              f"{mean512:.4f}; E9 truncated: 0.0946; E4 flat level: ~0.18")
    if mean512 < 0.13:
        verdict = ("RECONCILED-RP: the true operator's window IPR falls to "
                   "the truncated-ladder level -- the scaling onset is "
                   "physical and E9's RP reading stands for the true plate")
    elif mean512 > 0.15:
        verdict = ("PROTOCOL ARTIFACT: true-operator IPR stays flat at the "
                   "E4 level through N = 512 -- E9's scaling is a property "
                   "of the truncation protocol, not the plate")
    else:
        verdict = "INTERMEDIATE -- between the two references; see table"
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    save(results)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
