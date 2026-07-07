#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17 -- Gap A: true-operator windows at N = 1024-2048 (RUNNER).
See README.md (preregistered). Adapted from E14 with: scaled mesh/mode
count, the decisive G3 two-mesh gate over the full range, chunked
projections, and the preregistered N = 2048 verdict block."""
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
    mesh=(216, 133), mesh_check=(184, 114), k_fem=4, sigma_factor=10.0,
    n_modes=5600,
    gate_argyris_mesh=(128, 80), gate_n=1200,
    spacing_frac=0.1,
    n_funcs_axis=96, grid=(384, 240),
    ladder=[256, 512, 1024, 2048],
    probe_npts=1500,
    proj_chunk=800,
)
E9_SINE_2048 = 0.052   # E9 truncated-ladder reference at N = 2048
FLAT_SINE = 0.174      # E4/E14 flat true-operator level (N <= 512)


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def build_c0ip(cfg, nxy):
    mesh = MeshTri.init_tensor(np.linspace(0, cfg["a"], nxy[0] + 1),
                               np.linspace(0, cfg["b"], nxy[1] + 1))
    K, M, space = assemble_c0ip(mesh, k=cfg["k_fem"], nu=cfg["nu"],
                                sigma_factor=cfg["sigma_factor"])
    return 0.5 * (K + K.T), M, space


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}, "gates": {}}

    K, M, space = build_c0ip(cfg, cfg["mesh"])
    print(f"[setup] C0-IP P{cfg['k_fem']} {space.N} dofs "
          f"({time.time()-t00:.1f} s)")

    # ---------------- G1: SSSS vs exact at production mesh ----------------
    t0 = time.time()
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    lam_ss = solve_lowest(K[I][:, I].tocsc(), M[I][:, I].tocsc(),
                          cfg["n_modes"])
    lam_ex = ssss_exact(a, b, cfg["n_modes"])
    ns_ss = n_star(lam_ss, lam_ex, cfg["spacing_frac"])
    print(f"[G1] SSSS vs exact: N* = {ns_ss}/{cfg['n_modes']} "
          f"({time.time()-t0:.1f} s)")
    results["gates"]["G1"] = dict(n_star=int(ns_ss))
    save(results)

    # ---------------- G2: FFFF vs certified Argyris (first 1200) ----------
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
    print(f"[G2] FFFF vs Argyris: N* = {ns_ff}/{n_cmp}, rigid {nr} "
          f"(max {rig:.1e}) ({time.time()-t0:.1f} s)")
    results["gates"]["G2"] = dict(n_star=int(ns_ff), rigid=int(nr))
    save(results)

    # ---------------- production: certified eigenpairs ----------------
    t0 = time.time()
    lam, V, sinfo, _ = solve_modes(K.tocsc(), M.tocsc(), cfg["n_modes"] + 3,
                                   resid_sanity=1e-4, sweeps_max=40)
    lam, V, n_rigid, rigid_max = split_rigid(lam, V)
    print(f"[solve] {len(lam)} elastic, rigid {n_rigid} (max {rigid_max:.1e}),"
          f" resid {sinfo['max_resid']:.1e} ({time.time()-t0:.1f} s)")
    results["solve"] = dict(resid=sinfo["max_resid"], rigid=int(n_rigid),
                            n_elastic=int(len(lam)))
    save(results)

    # ---------------- G3: decisive two-mesh gate over the full range ------
    t0 = time.time()
    Kc, Mc, space_c = build_c0ip(cfg, cfg["mesh_check"])
    lam_c = solve_lowest(Kc.tocsc(), Mc.tocsc(), cfg["n_modes"] + 3)
    lam_c, _, nr_c, _ = split_rigid(lam_c)
    del Kc, Mc, space_c
    n_cmp = min(len(lam), len(lam_c))
    ns_full = n_star(lam[:n_cmp], lam_c[:n_cmp], cfg["spacing_frac"])
    print(f"[G3] FFFF two-mesh {cfg['mesh']} vs {cfg['mesh_check']}: "
          f"N* = {ns_full}/{n_cmp} ({time.time()-t0:.1f} s)")
    results["gates"]["G3"] = dict(n_star=int(ns_full), rigid_check=int(nr_c))
    save(results)
    n_use = int(min(ns_full, len(lam)))

    # ---------------- classification ----------------
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

    # ---------------- projections (chunked) + ladders ----------------
    t0 = time.time()
    gx, wx, gy, wy = bases.gauss_grid(a, b, *cfg["grid"])
    Xg, Yg = np.meshgrid(gx, gy, indexing="ij")
    Pgrid = space.probes(np.vstack([Xg.ravel(), Yg.ravel()]))
    Sx, kSx, pSx = bases.sine_1d(cfg["n_funcs_axis"], a, gx)
    Sy, kSy, pSy = bases.sine_1d(cfg["n_funcs_axis"], b, gy)
    Bx, kBx, pBx = bases.beam_1d(cfg["n_funcs_axis"], a, gx, wx)
    By, kBy, pBy = bases.beam_1d(cfg["n_funcs_axis"], b, gy, wy)
    cs_parts, cb_parts = [], []
    for j0 in range(0, n_use, cfg["proj_chunk"]):
        Wg = Pgrid @ V[:, j0:j0 + cfg["proj_chunk"]]
        cs_parts.append(bases.project_modes(Wg, Sx, wx, Sy, wy, *cfg["grid"]))
        cb_parts.append(bases.project_modes(Wg, Bx, wx, By, wy, *cfg["grid"]))
        print(f"[proj] {min(j0+cfg['proj_chunk'], n_use)}/{n_use} "
              f"({time.time()-t0:.1f} s)")
    C_sine = np.concatenate(cs_parts, axis=0)
    C_beam = np.concatenate(cb_parts, axis=0)
    del cs_parts, cb_parts

    out = {"sectors": {}}
    for s in SECTORS:
        idx = [i for i, l in enumerate(labels[:n_use]) if l == s]
        rec = {"n": len(idx)}
        for name, (C, kx, px, ky, py) in [("sine", (C_sine, kSx, pSx, kSy, pSy)),
                                          ("beam", (C_beam, kBx, pBx, kBy, pBy))]:
            pi, pj = bases.sector_product_order(kx, px, ky, py, s,
                                                max(cfg["ladder"]))
            rec[name] = bases.ipr_ladder(C[idx], pi, pj, cfg["ladder"])
        out["sectors"][s] = rec
        print(f"[gapA] {s}: {len(idx)} modes (of gate-covered {n_use}), "
              f"ladders done")
    out["goe"] = bases.goe_ipr_baseline(cfg["ladder"])
    results["gapA"] = out
    save(results)

    # ---------------- preregistered verdict ----------------
    md = ["# E17 -- Gap A true-operator windows at N = 1024-2048 (RESULTS)\n",
          f"C0-IP P{cfg['k_fem']} at {space.N} dofs, {len(lam)} certified "
          f"elastic modes, gate-covered n_use = {n_use}. References: E9 "
          f"truncated sine at N = 2048: {E9_SINE_2048}; E4/E14 flat level: "
          f"~{FLAT_SINE}.\n"]
    for g in ("G1", "G2", "G3"):
        md.append(f"- {g}: {results['gates'][g]}")
    md.append("\n| sector | basis | " + " | ".join(f"IPR N={N}" for N in
              cfg["ladder"]) + " |")
    md.append("|" + "---|" * (len(cfg["ladder"]) + 2))
    per_rung = {}
    for s in SECTORS:
        for name in ("sine", "beam"):
            rows = out["sectors"][s][name]
            iprs = [float(np.exp(r["mlnipr"])) for r in rows]
            for r, v in zip(rows, iprs):
                per_rung.setdefault((name, int(r["N"])), []).append(v)
            md.append(f"| {s} | {name} | "
                      + " | ".join(f"{v:.4f}" for v in iprs)
                      + " |" * (len(cfg["ladder"]) - len(iprs) + 1))
    top = max(N for (nm, N) in per_rung if nm == "sine")
    i_top = float(np.mean(per_rung[("sine", top)]))
    lnN, lnI = [], []
    for N in cfg["ladder"]:
        if ("sine", N) in per_rung and N >= 512:
            lnN.append(np.log(N))
            lnI.append(np.log(np.mean(per_rung[("sine", N)])))
    d2_true = float(-np.polyfit(lnN, lnI, 1)[0]) if len(lnN) >= 2 else np.nan
    md.append(f"\n- sine sector-mean IPR at top covered rung N = {top}: "
              f"{i_top:.4f}; fitted D2_true (N >= 512): {d2_true:.3f}")
    if top < max(cfg["ladder"]):
        verdict = (f"COVERAGE-LIMITED: gates/sector coverage stop at N = "
                   f"{top} < 2048; largest covered rung reported; "
                   f"continuation mesh registered")
    elif i_top >= 0.14:
        verdict = ("PROTOCOL ARTIFACT THROUGH THE FULL REGISTERED RANGE: "
                   "true-operator windows stay sparse to N = 2048; "
                   "truncated-ladder RP is an artifact at every registered "
                   "rung")
    elif i_top <= 0.09:
        verdict = (f"SCALING ONSET: genuine delocalization inside the "
                   f"registered range (D2_true = {d2_true:.3f}; "
                   f"RECONCILED-RP if within E9's 0.42-0.50)")
    else:
        verdict = (f"INTERMEDIATE: between references at N = 2048 "
                   f"(D2_true = {d2_true:.3f}); see table for where "
                   f"flatness breaks")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
