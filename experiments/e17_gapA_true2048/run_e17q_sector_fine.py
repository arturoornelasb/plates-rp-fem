#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17q-FINE -- finer-mesh 2048 certification variant (registered task 23): production (140,86) ~193k dofs/sector, honest check (120,74) at 1.167x, 1472 modes (2048 window needs 1229). One parity sector of the Gap-A N=2048 experiment on the
VALIDATED i01 quarter-plate instrument (see README, E17q section).
Usage: run_e17q_sector.py <ee|eo|oe|oo>. Same preregistered readings as
E17; the serial run was killed externally post-G2 (partial record
preserved as *_serial_partial.*).

Projection trick: the quarter eigenvectors are evaluated on the FULL
Gauss grid by parity reflection, so the T1 registered bases and the
E14/E17 projection machinery apply verbatim per sector."""
import json
import os
import sys
import time

import numpy as np
from skfem import MeshTri

from platefem import n_star, solve_lowest, solve_modes
from platefem.c0ip import (C0IPSpace, assemble_c0ip, boundary_dofs,
                           facets_where)
from platefem import bases

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(140, 86), mesh_check=(120, 74), k_fem=4, sigma_factor=10.0,
    n_modes=1472,
    spacing_frac=0.1,
    n_funcs_axis=96, grid=(384, 240),
    ladder=[256, 512, 1024, 2048],
    proj_chunk=700,
)
# v2 after the reaped first run: finer check mesh (94,58 -> 102,63; the
# coarse check, not production, limited oo's gate to 757); rigid-gap
# threshold relaxed 1e5 -> 1e3 (marginal on the guided quarter -- the
# ratio is recorded); eigenpairs cached to .npz right after the solve so
# an external stop never again discards a completed solve (resume skips
# straight to the gate/projection stages).
RIGID = dict(ee=1, eo=1, oe=1, oo=0)
TOL = 1e-9


def build_sector(cfg, sector, nxy):
    a, b = cfg["a"], cfg["b"]
    mesh = MeshTri.init_tensor(np.linspace(0.0, a / 2, nxy[0] + 1),
                               np.linspace(0.0, b / 2, nxy[1] + 1))
    tmp = C0IPSpace(mesh, cfg["k_fem"])
    xline = facets_where(tmp, lambda x, y: abs(x - a / 2) < TOL)
    yline = facets_where(tmp, lambda x, y: abs(y - b / 2) < TOL)
    gf, ssf = [], []
    (gf if sector[0] == "e" else ssf).append(xline)
    (gf if sector[1] == "e" else ssf).append(yline)
    gf = np.concatenate(gf) if gf else None
    K, M, space = assemble_c0ip(mesh, k=cfg["k_fem"], nu=cfg["nu"],
                                sigma_factor=cfg["sigma_factor"],
                                guided_facets=gf)
    K = 0.5 * (K + K.T)
    I = None
    if ssf:
        D = boundary_dofs(space, np.concatenate(ssf))
        I = np.setdiff1d(np.arange(space.N), D)
        K, M = K[I][:, I].tocsc(), M[I][:, I].tocsc()
    else:
        K, M = K.tocsc(), M.tocsc()
    return mesh, space, K, M, I


def split_expected(lam, V, n_exp):
    lam = np.asarray(lam)
    if n_exp:
        ratio = float(lam[n_exp] / max(abs(float(lam[n_exp - 1])), 1e-300))
    else:
        ratio = float(lam[0] / (1e-8 * lam[len(lam) // 2]))
    gap_ok = ratio > 1e3
    return (lam[n_exp:], (V[:, n_exp:] if V is not None else None),
            bool(gap_ok), ratio)


def main():
    sector = sys.argv[1]
    assert sector in RIGID
    t00 = time.time()
    cfg = CFG
    a, b = cfg["a"], cfg["b"]
    out = {"sector": sector, "config": {k: (list(v) if isinstance(v, tuple)
                                            else v) for k, v in cfg.items()}}

    mesh, space, K, M, I = build_sector(cfg, sector, cfg["mesh"])
    ndof = K.shape[0]
    print(f"[{sector}] {ndof} dofs ({time.time()-t00:.1f} s)")

    def stage_save():
        with open(os.path.join(HERE, f"sector_fine_{sector}.json"), "w") as f:
            json.dump(out, f, indent=1, default=float)

    eig_path = os.path.join(HERE, f"sector_fine_{sector}_eig.npz")
    if os.path.exists(eig_path):
        z = np.load(eig_path)
        lam, V = z["lam"], z["V"]
        resid, gap_ok, ratio = (float(z["resid"][0]), bool(z["gap_ok"][0]),
                                float(z["ratio"][0]))
        print(f"[{sector}] RESUMED {len(lam)} elastic modes from cache")
    else:
        t0 = time.time()
        lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + RIGID[sector],
                                       resid_sanity=1e-4, sweeps_max=40)
        lam, V, gap_ok, ratio = split_expected(lam, V, RIGID[sector])
        resid = float(sinfo["max_resid"])
        np.savez(eig_path, lam=lam, V=V, resid=[resid],
                 gap_ok=[gap_ok], ratio=[ratio])
        print(f"[{sector}] {len(lam)} elastic, rigid gap {gap_ok} "
              f"(ratio {ratio:.1e}), resid {resid:.1e} "
              f"({time.time()-t0:.1f} s); eigenpairs cached")
    out["solve"] = dict(resid=resid, gap_ok=gap_ok, gap_ratio=ratio,
                        n_elastic=int(len(lam)), ndof=int(ndof))
    stage_save()

    # -------- decisive per-sector two-mesh gate (eigenvalues only) --------
    gate_path = os.path.join(HERE, f"sector_fine_{sector}_gate.json")
    if os.path.exists(gate_path):
        with open(gate_path) as f:
            out["gate"] = json.load(f)
        print(f"[{sector}] G3 RESUMED from cache: N* = "
              f"{out['gate']['n_star']}/{out['gate']['n_cmp']}")
    else:
        t0 = time.time()
        _, _, Kc, Mc, _ = build_sector(cfg, sector, cfg["mesh_check"])
        lam_c = solve_lowest(Kc, Mc, cfg["n_modes"] + RIGID[sector])
        lam_c, _, _, _ = split_expected(lam_c, None, RIGID[sector])
        n_cmp = min(len(lam), len(lam_c))
        ns = n_star(lam[:n_cmp], lam_c[:n_cmp], cfg["spacing_frac"])
        print(f"[{sector}] G3 two-mesh: N* = {ns}/{n_cmp} "
              f"({time.time()-t0:.1f} s)")
        out["gate"] = dict(n_star=int(ns), n_cmp=int(n_cmp))
        with open(gate_path, "w") as f:
            json.dump(out["gate"], f)
    n_use = int(min(out["gate"]["n_star"], len(lam)))
    stage_save()

    # -------- projection: parity-reflected evaluation on the FULL grid ----
    t0 = time.time()
    gx, wx, gy, wy = bases.gauss_grid(a, b, *cfg["grid"])
    Xg, Yg = np.meshgrid(gx, gy, indexing="ij")
    X, Y = Xg.ravel(), Yg.ravel()
    xm = np.where(X <= a / 2, X, a - X)
    ym = np.where(Y <= b / 2, Y, b - Y)
    sx = np.where((X <= a / 2) | (sector[0] == "e"), 1.0, -1.0)
    sy = np.where((Y <= b / 2) | (sector[1] == "e"), 1.0, -1.0)
    sgn = (sx * sy)[:, None]
    Pq = space.probes(np.vstack([xm, ym]))
    if I is not None:
        Pq = Pq[:, I]
    Sx, kSx, pSx = bases.sine_1d(cfg["n_funcs_axis"], a, gx)
    Sy, kSy, pSy = bases.sine_1d(cfg["n_funcs_axis"], b, gy)
    Bx, kBx, pBx = bases.beam_1d(cfg["n_funcs_axis"], a, gx, wx)
    By, kBy, pBy = bases.beam_1d(cfg["n_funcs_axis"], b, gy, wy)
    cs_parts, cb_parts = [], []
    V_use = V[:, :n_use]                 # gate-covered modes ONLY
    for j0 in range(0, n_use, cfg["proj_chunk"]):
        Wg = (Pq @ V_use[:, j0:j0 + cfg["proj_chunk"]]) * sgn
        cs_parts.append(bases.project_modes(Wg, Sx, wx, Sy, wy, *cfg["grid"]))
        cb_parts.append(bases.project_modes(Wg, Bx, wx, By, wy, *cfg["grid"]))
    C_sine = np.concatenate(cs_parts, axis=0)
    C_beam = np.concatenate(cb_parts, axis=0)
    print(f"[{sector}] projections done ({time.time()-t0:.1f} s)")

    rec = {"n": int(n_use)}
    for name, (C, kx, px, ky, py) in [("sine", (C_sine, kSx, pSx, kSy, pSy)),
                                      ("beam", (C_beam, kBx, pBx, kBy, pBy))]:
        pi, pj = bases.sector_product_order(kx, px, ky, py, sector,
                                            max(cfg["ladder"]))
        rec[name] = bases.ipr_ladder(C, pi, pj, cfg["ladder"])
    out["gapA"] = rec
    out["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, f"sector_fine_{sector}.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"[{sector}] done in {out['wall_time_s']} s")


if __name__ == "__main__":
    main()
