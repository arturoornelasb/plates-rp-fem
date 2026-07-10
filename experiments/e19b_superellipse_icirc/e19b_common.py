#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19b shared machinery: the init_circle-mapped superellipse mesh (the
E5-validated, quality-measured family) replaces E19's failed polar-ring
construction. Ladders scoped to the refine-6 coverage."""
import numpy as np
from platefem import superellipse_basis

RATIO = 1.6189043236
A_AX = float(np.sqrt(RATIO))
B_AX = 1.0 / A_AX
P_EXP = 10.0
NU = 0.33
KFEM = 4
LEVEL_PROD, LEVEL_CHECK = 6, 5
N_MODES = dict(free=668, ss=1040)
RIGID = dict(free=3, ss=0)
LADDER = [128, 256]


def build_mesh(level):
    mesh, _ = superellipse_basis(level, A_AX, B_AX, P_EXP)
    return mesh


def split_expected(lam, V, n_exp):
    lam = np.asarray(lam)
    if n_exp:
        ratio = float(lam[n_exp] / max(abs(float(lam[n_exp - 1])), 1e-300))
    else:
        ratio = float("inf")
    return (lam[n_exp:], (V[:, n_exp:] if V is not None else None),
            bool(ratio > 1e3), ratio)
