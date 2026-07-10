#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19 shared machinery: polar-structured superellipse mesh (fine
granularity; mirror-symmetric for nth divisible by 4) and the job
configuration. Geometry matches E5 exactly (a = sqrt(ratio), b = 1/a,
p = 10)."""
import numpy as np
from skfem import MeshTri

RATIO = 1.6189043236
A_AX = float(np.sqrt(RATIO))
B_AX = 1.0 / A_AX
P_EXP = 10.0
NU = 0.33
KFEM = 4
MESH_PROD = (64, 192)      # (n rings, n theta); nth % 4 == 0
MESH_CHECK = (56, 168)
N_MODES = dict(free=1336, ss=2120)
RIGID = dict(free=3, ss=0)
LADDER = [128, 256, 512]


def superellipse_polar_mesh(nr, nth, a=A_AX, b=B_AX, p=P_EXP):
    """Structured polar mesh radially mapped onto the superellipse:
    a center vertex, nr rings, nth angular divisions. Mirror-symmetric
    about both axes when nth % 4 == 0."""
    th = np.linspace(0.0, 2 * np.pi, nth, endpoint=False)
    rho_b = (np.abs(np.cos(th) / a) ** p
             + np.abs(np.sin(th) / b) ** p) ** (-1.0 / p)
    pts = [(0.0, 0.0)]
    for i in range(1, nr + 1):
        s = i / nr
        for j in range(nth):
            pts.append((s * rho_b[j] * np.cos(th[j]),
                        s * rho_b[j] * np.sin(th[j])))
    def vid(i, j):
        return 0 if i == 0 else 1 + (i - 1) * nth + (j % nth)
    tris = []
    for j in range(nth):                       # center fan
        tris.append([vid(0, 0), vid(1, j), vid(1, j + 1)])
    for i in range(1, nr):                     # ring quads split
        for j in range(nth):
            v00, v01 = vid(i, j), vid(i, j + 1)
            v10, v11 = vid(i + 1, j), vid(i + 1, j + 1)
            # alternate the diagonal by quadrant-symmetric parity to keep
            # the mesh mirror-symmetric
            if (j < nth // 4) or (nth // 2 <= j < 3 * nth // 4):
                tris.append([v00, v10, v11])
                tris.append([v00, v11, v01])
            else:
                tris.append([v01, v00, v10])
                tris.append([v01, v10, v11])
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
