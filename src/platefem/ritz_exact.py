#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exactly-assembled Legendre-Ritz matrices for the Kirchhoff rectangle.

The quadrature-assembled Ritz reference (platefem.ritz, repo T1 protocol) has
a measured ~1e-4 round-off floor: eigenvalues go non-monotone in NLEG beyond
~96 because forming A2 = (V2 wg) V2^T mixes entries spanning many orders of
magnitude (documented in experiments/e01_validation/diag_ritz_nleg.py).

This module assembles the SAME matrices from closed forms using exact
integer arithmetic (python ints), converting to float64 only at the end:
with the classical expansion P_n' = sum_{k = n-1, n-3, ...} (2k+1) P_k, the
derivative operator in the orthonormal-Legendre basis is an exact sparse
integer-weighted matrix D, and

    A1 = D^T G D,   A2 = (D^2)^T G (D^2),   B = (D^2)^T G  ... etc.

where G = diag(2/(2n+1)) is the Legendre Gram. All products are carried out
in Fraction-free integer/rational arithmetic scaled to the orthonormal
convention at the end. Validated against the quadrature assembly at NLEG=48
(where the quadrature path is still clean) to machine precision.
"""
from fractions import Fraction

import numpy as np


def _deriv_matrix(n):
    """D[k, m] = coefficient of P_k in P_m' (exact ints), size n x n."""
    D = [[0] * n for _ in range(n)]
    for m in range(n):
        k = m - 1
        while k >= 0:
            D[k][m] = 2 * k + 1
            k -= 2
    return D


def _matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n) if A[i][k])
             for j in range(n)] for i in range(n)]


def exact_1d_matrices(nleg):
    """Exact (A0, A1, A2, B) in the ORTHONORMAL Legendre convention on [-1,1]:
    A0 = I; A1[i,j] = <phi_i', phi_j'>; A2 = <phi_i'', phi_j''>;
    B[i,j] = <phi_i'', phi_j> (matching platefem.ritz.mats_1d)."""
    D = _deriv_matrix(nleg)
    D2 = _matmul(D, D)
    # Gram of monic-normalized P_n: g_n = 2/(2n+1) (rational)
    g = [Fraction(2, 2 * k + 1) for k in range(nleg)]

    def gram_product(Xt, Y):
        # (X^T G Y)[i,j] = sum_k X[k][i] g_k Y[k][j]
        out = [[Fraction(0)] * nleg for _ in range(nleg)]
        for k in range(nleg):
            gk = g[k]
            Xk, Yk = Xt[k], Y[k]
            nz = [j for j in range(nleg) if Xk[j] or Yk[j]]
            for i in nz:
                if not Xk[i]:
                    continue
                xi = Xk[i] * gk
                for j in nz:
                    if Yk[j]:
                        out[i][j] += xi * Yk[j]
        return out

    A1f = gram_product(D, D)
    A2f = gram_product(D2, D2)
    # B[i,j] = <P_i'', P_j> = sum_k D2[k][i] g_k delta_kj = D2[j][i] g_j
    Bf = [[D2[j][i] * g[j] for j in range(nleg)] for i in range(nleg)]

    # orthonormal scaling: phi_n = sqrt((2n+1)/2) P_n
    s = np.sqrt((2.0 * np.arange(nleg) + 1.0) / 2.0)
    A1 = np.array([[float(A1f[i][j]) for j in range(nleg)]
                   for i in range(nleg)]) * np.outer(s, s)
    A2 = np.array([[float(A2f[i][j]) for j in range(nleg)]
                   for i in range(nleg)]) * np.outer(s, s)
    B = np.array([[float(Bf[i][j]) for j in range(nleg)]
                  for i in range(nleg)]) * np.outer(s, s)
    A0 = np.eye(nleg)
    return A0, A1, A2, B


def assemble_K_sector(idx_x, idx_y, A1, A2, B, a, b, nu):
    """Identical to platefem.ritz.assemble_K_sector (kept in sync)."""
    Jx, Jy = 2.0 / a, 2.0 / b
    pre = (a * b) / 4.0
    A1x, A2x, Bx = A1[np.ix_(idx_x, idx_x)], A2[np.ix_(idx_x, idx_x)], B[np.ix_(idx_x, idx_x)]
    A1y, A2y, By = A1[np.ix_(idx_y, idx_y)], A2[np.ix_(idx_y, idx_y)], B[np.ix_(idx_y, idx_y)]
    Ix, Iy = np.eye(len(idx_x)), np.eye(len(idx_y))
    K = (Jx ** 4) * np.kron(A2x, Iy) + (Jy ** 4) * np.kron(Ix, A2y)
    K += nu * (Jx ** 2) * (Jy ** 2) * (np.kron(Bx, By.T) + np.kron(Bx.T, By))
    K += 2.0 * (1.0 - nu) * (Jx ** 2) * (Jy ** 2) * np.kron(A1x, A1y)
    return pre * K
