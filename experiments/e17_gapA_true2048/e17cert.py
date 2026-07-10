#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17 continuation -- 2048-window certification gates (registered).
Usage: e17cert.py <ee|eo|oe|oo> <check92|anchor110>

check92:   honest two-mesh -- quarter FFFF sector eigenvalues at (92, 57)
           (h ratio 1.196 vs production (110, 68); the (102, 63) check was
           1.08x = correlated-optimistic).
anchor110: absolute anchor -- quarter with SS OUTER edges + per-sector
           center-line BCs at the PRODUCTION mesh, vs the exact SSSS
           sector spectrum (measures absolute mesh accuracy per sector at
           production h; informational, E17-serial-G1 style)."""
import os
import sys
import time

import numpy as np

from platefem import solve_lowest
from platefem.c0ip import (C0IPSpace, assemble_c0ip, boundary_dofs,
                           facets_where)
from run_e17q_sector import CFG, RIGID, split_expected, build_sector

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-9


def build_anchor(cfg, sector, nxy):
    """Quarter with SS outer edges + parity center-line BCs (per i01)."""
    from skfem import MeshTri
    a, b = cfg["a"], cfg["b"]
    mesh = MeshTri.init_tensor(np.linspace(0.0, a / 2, nxy[0] + 1),
                               np.linspace(0.0, b / 2, nxy[1] + 1))
    tmp = C0IPSpace(mesh, cfg["k_fem"])
    xline = facets_where(tmp, lambda x, y: abs(x - a / 2) < TOL)
    yline = facets_where(tmp, lambda x, y: abs(y - b / 2) < TOL)
    gf, ssf = [], []
    (gf if sector[0] == "e" else ssf).append(xline)
    (gf if sector[1] == "e" else ssf).append(yline)
    ssf.append(facets_where(tmp, lambda x, y: abs(x) < TOL))
    ssf.append(facets_where(tmp, lambda x, y: abs(y) < TOL))
    gf = np.concatenate(gf) if gf else None
    K, M, space = assemble_c0ip(mesh, k=cfg["k_fem"], nu=cfg["nu"],
                                sigma_factor=cfg["sigma_factor"],
                                guided_facets=gf)
    K = 0.5 * (K + K.T)
    D = boundary_dofs(space, np.concatenate(ssf))
    I = np.setdiff1d(np.arange(space.N), D)
    return K[I][:, I].tocsc(), M[I][:, I].tocsc()


def main():
    sector, kind = sys.argv[1], sys.argv[2]
    t00 = time.time()
    if kind == "check92":
        _, _, K, M, _ = build_sector(CFG, sector, (92, 57))
        n_req = CFG["n_modes"] + RIGID[sector]
        lam = solve_lowest(K, M, n_req)
        lam, _, _, _ = split_expected(lam, None, RIGID[sector])
    else:
        K, M = build_anchor(CFG, sector, CFG["mesh"])
        lam = solve_lowest(K, M, CFG["n_modes"])
        lam = np.asarray(lam)
    np.savez(os.path.join(HERE, f"cert_{kind}_{sector}.npz"), lam=lam,
             ndof=[K.shape[0]])
    print(f"[cert {kind} {sector}] {len(lam)} eigenvalues, {K.shape[0]} "
          f"dofs ({time.time()-t00:.1f} s)")


if __name__ == "__main__":
    main()
