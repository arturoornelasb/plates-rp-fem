#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19c shared constants: refine-7 free instrument + extended refine-6
SS basis with cross-mesh projection (see README.md, frozen)."""
import numpy as np
from platefem import superellipse_basis

RATIO = 1.6189043236
A_AX = float(np.sqrt(RATIO))
B_AX = 1.0 / A_AX
P_EXP = 10.0
NU = 0.33
KFEM = 4
LEVEL_FREE, LEVEL_SS = 7, 6
N_FREE, N_SS = 1500, 2300
RIGID_FREE = 3
LADDER = [128, 256, 512]
COLLAR = 2e-3


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


def collar_pullback(pts, a=A_AX, b=B_AX, p=P_EXP, delta=COLLAR):
    """Radially pull points with superellipse coordinate sigma > 1-delta
    back to sigma = 1-delta (chord-vs-curve sliver handling; README)."""
    sig = (np.abs(pts[0] / a) ** p + np.abs(pts[1] / b) ** p) ** (1.0 / p)
    scale = np.where(sig > 1.0 - delta, (1.0 - delta) / np.maximum(sig, 1e-300),
                     1.0)
    return pts * scale
