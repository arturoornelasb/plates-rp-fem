#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18 shared machinery: C3v-symmetric barycentric triangle mesh (matching
the E3a macro-triangle geometry exactly) and the job configurations."""
import numpy as np
from skfem import MeshTri

L_SIDE = 1.0
NU = 0.33
KFEM = 4
N_PROD, N_CHECK = 190, 166
N_MODES = dict(free=1056, ss=2000)
RIGID = dict(free=3, ss=0)
LADDER_A = [128, 256]
LADDER_E = [128, 256, 512, 1024]
TRI_SECTORS = ["A1", "A2", "E"]


def macro_verts(L=L_SIDE):
    return np.array([[-0.5 * L, -np.sqrt(3) / 6 * L],
                     [0.5 * L, -np.sqrt(3) / 6 * L],
                     [0.0, np.sqrt(3) / 3 * L]])


def tri_mesh_bary(n, L=L_SIDE):
    """Barycentric n-subdivision of the macro equilateral triangle --
    C3v-symmetric for every n (the lattice and the up/down triangulation
    map to themselves under the 120-degree rotation and the mirrors)."""
    v = macro_verts(L)
    idx = {}
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            idx[(i, j)] = len(pts)
            pts.append(v[0] + (i / n) * (v[1] - v[0])
                       + (j / n) * (v[2] - v[0]))
    tris = []
    for i in range(n):
        for j in range(n - i):
            tris.append([idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]])
            if i + j < n - 1:
                tris.append([idx[(i + 1, j)], idx[(i + 1, j + 1)],
                             idx[(i, j + 1)]])
    return MeshTri(np.array(pts).T.copy(),
                   np.array(tris, dtype=np.int64).T.copy())


def split_expected(lam, V, n_exp):
    lam = np.asarray(lam)
    if n_exp:
        ratio = float(lam[n_exp] / max(abs(float(lam[n_exp - 1])), 1e-300))
    else:
        ratio = float("inf")
    return (lam[n_exp:], (V[:, n_exp:] if V is not None else None),
            bool(ratio > 1e3), ratio)
