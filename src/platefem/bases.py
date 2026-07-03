#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reference bases and projections for Gap A (paper Sec. 'Gap A'; T1 protocol).

The T1 audit established that IPR/D2 of a finite representation depends on
the REPRESENTATION BASIS, so Gap A must be computed in the registered bases:
  (b) 'sine'  -- simply supported products sqrt(2/L) sin(m pi x / L)  (the H0
                 basis of the Gap A recipe),
  (c) 'beam'  -- free-free beam-function products (classical Ritz basis for
                 FFFF plates; cosh/cos and sinh/sin forms).
Never the raw FEM nodal basis.

All 1D families are L2(dx)-orthonormal on [0, L]; parity is about x = L/2
(sine: odd m -> even parity; beam: even family includes the rigid constant,
odd family includes the rigid tilt). Projections use tensor Gauss-Legendre
quadrature; captured norm p < 1 is reported per mode and coefficient vectors
are renormalized before IPR (T1 convention).
"""
import numpy as np
from scipy.optimize import brentq


# ---------------- 1D families on [0, L], orthonormal in L2(dx) ----------------
def sine_1d(n_funcs, L, x):
    """Rows: sqrt(2/L) sin(m pi x / L), m = 1..n_funcs; parity(m) = 'e' for
    odd m, 'o' for even m (about x = L/2). Returns (values, k_proxy, parity)."""
    m = np.arange(1, n_funcs + 1)
    vals = np.sqrt(2.0 / L) * np.sin(np.outer(m, np.pi * x / L))
    k = m * np.pi / L
    par = np.where(m % 2 == 1, "e", "o")
    return vals, k, par


def _even_beam_roots(n):
    """tan(b) = -tanh(b) on xi in [-1, 1] (T1 convention): 2.36502, 5.49780..."""
    g = lambda b: np.sin(b) + np.cos(b) * np.tanh(b)
    out = []
    for k in range(n):
        out.append(brentq(g, k * np.pi + np.pi / 2 + 1e-10,
                          (k + 1) * np.pi - 1e-10, xtol=1e-14))
    return np.array(out)


def _odd_beam_roots(n):
    """tan(b) = +tanh(b): 3.92660, 7.06858, ..."""
    g = lambda b: np.sin(b) - np.cos(b) * np.tanh(b)
    out = []
    for k in range(1, n + 1):
        out.append(brentq(g, k * np.pi + 1e-10, k * np.pi + np.pi / 2 - 1e-10,
                          xtol=1e-14))
    return np.array(out)


def _cosh_ratio(b, xi):
    ax = np.abs(xi)
    return np.exp(b * (ax - 1.0)) * (1.0 + np.exp(-2.0 * b * ax)) / \
        (1.0 + np.exp(-2.0 * b))


def _sinh_ratio(b, xi):
    s = np.sign(xi)
    ax = np.abs(xi)
    return s * np.exp(b * (ax - 1.0)) * (1.0 - np.exp(-2.0 * b * ax)) / \
        (1.0 - np.exp(-2.0 * b))


def beam_1d(n_funcs, L, x, wq):
    """Free-free beam functions on [0, L] (xi = 2x/L - 1), numerically
    L2(dx)-normalized then Loewdin-orthonormalized w.r.t. the quadrature
    weights wq. Even family: const + cos(b xi) + cos(b) cosh-ratio; odd:
    tilt + sin(b xi) + sin(b) sinh-ratio. Returns (values, k_proxy, parity)."""
    xi = 2.0 * x / L - 1.0
    n_e = (n_funcs + 1) // 2
    n_o = n_funcs - n_e
    rows, ks, par = [np.ones_like(xi)], [0.0], ["e"]
    for b in _even_beam_roots(n_e - 1):
        rows.append(np.cos(b * xi) + np.cos(b) * _cosh_ratio(b, xi))
        ks.append(2.0 * b / L)
        par.append("e")
    rows.append(xi.copy())
    ks.append(1.0 * np.pi / L)          # tilt frequency proxy between 0 and b1
    par.append("o")
    for b in _odd_beam_roots(n_o - 1):
        rows.append(np.sin(b * xi) + np.sin(b) * _sinh_ratio(b, xi))
        ks.append(2.0 * b / L)
        par.append("o")
    W = np.array(rows)
    W /= np.sqrt(np.sum(wq * W ** 2, axis=1))[:, None]
    # Loewdin re-orthonormalization (families are near- but not exactly
    # orthogonal under quadrature)
    G = (W * wq) @ W.T
    ev, U = np.linalg.eigh(0.5 * (G + G.T))
    W = (U / np.sqrt(np.clip(ev, 1e-14, None))) @ U.T @ W
    return W, np.array(ks), np.array(par)


# ---------------- tensor projection of FEM modes ----------------
def gauss_grid(a, b, nx, ny):
    """Tensor Gauss-Legendre points/weights on [0,a]x[0,b]."""
    xg, wx = np.polynomial.legendre.leggauss(nx)
    yg, wy = np.polynomial.legendre.leggauss(ny)
    x = 0.5 * a * (xg + 1.0)
    y = 0.5 * b * (yg + 1.0)
    wx = wx * 0.5 * a
    wy = wy * 0.5 * b
    return x, wx, y, wy


def project_modes(Wgrid, Fx, wx, Fy, wy, nx, ny):
    """Coefficients of FEM modes in a tensor family. Wgrid: (nx*ny, n_modes)
    values on the grid (x-major). Fx: (nfx, nx), Fy: (nfy, ny). Returns
    C with shape (n_modes, nfx, nfy)."""
    n_modes = Wgrid.shape[1]
    W3 = Wgrid.reshape(nx, ny, n_modes)
    A = Fx * wx                              # (nfx, nx)
    B = Fy * wy                              # (nfy, ny)
    return np.einsum("ix,xym,jy->mij", A, W3, B, optimize=True)


def sector_product_order(kx, px, ky, py, sector, n_max):
    """Indices (i, j) of the lowest-n_max tensor products belonging to a
    parity sector, ordered by the flexural frequency proxy kx^2 + ky^2."""
    ix = np.where(px == sector[0])[0]
    iy = np.where(py == sector[1])[0]
    key = (kx[ix][:, None] ** 2 + ky[iy][None, :] ** 2).ravel()
    order = np.argsort(key, kind="stable")[:n_max]
    ii, jj = np.unravel_index(order, (len(ix), len(iy)))
    return ix[ii], iy[jj]


def ipr_ladder(C_sector, prods_i, prods_j, ladder, win=(0.4, 0.6)):
    """T1 ladder on TRUE modes: for each N, represent the window modes
    [win0*N, win1*N) of the sector in its lowest-N products, renormalize the
    captured coefficients, IPR = sum c^4. C_sector: (n_modes_sector, nfx, nfy)
    ordered by frequency. Returns list of dicts per N."""
    out = []
    n_sec = C_sector.shape[0]
    for N in ladder:
        i0, i1 = int(win[0] * N), int(win[1] * N)
        if i1 > n_sec:
            break
        pi, pj = prods_i[:N], prods_j[:N]
        iprs, caps = [], []
        for k in range(i0, i1):
            d = C_sector[k, pi, pj]
            p = float(d @ d)
            if p <= 0:
                continue
            dn = d / np.sqrt(p)
            iprs.append(float(np.sum(dn ** 4)))
            caps.append(p)
        out.append(dict(N=int(N), n_window=len(iprs),
                        mlnipr=float(np.mean(np.log(iprs))),
                        slnipr=float(np.std(np.log(iprs))),
                        mean_captured=float(np.mean(caps))))
    return out


def goe_ipr_baseline(ladder, win=(0.4, 0.6), n_real=8, seed=12345):
    """Empirical GOE mean-ln-IPR at each ladder size (T1 convention)."""
    out = {}
    for N in ladder:
        rng = np.random.default_rng(seed + N)
        g0, g1 = int(win[0] * N), int(win[1] * N)
        iprs = []
        for _ in range(n_real):
            H = rng.standard_normal((N, N))
            _, U = np.linalg.eigh(0.5 * (H + H.T))
            iprs.extend(np.sum(U[:, g0:g1] ** 4, axis=0).tolist())
        out[int(N)] = float(np.mean(np.log(iprs)))
    return out


def number_variance(lam, L_values, poly_deg=5):
    """Sigma^2(L) of one level sequence, unfolded by a polynomial fit of the
    staircase against sqrt(lambda) (Weyl-linear variable for plates)."""
    lam = np.sort(np.asarray(lam, float))
    x = np.sqrt(lam)
    stair = np.arange(1, len(lam) + 1, dtype=float)
    coef = np.polyfit(x, stair, poly_deg)
    u = np.polyval(coef, x)                  # unfolded positions
    out = []
    for L in L_values:
        lo, hi = u[0], u[-1] - L
        if hi <= lo:
            out.append(np.nan)
            continue
        starts = np.linspace(lo, hi, 200)
        counts = np.searchsorted(u, starts + L) - np.searchsorted(u, starts)
        out.append(float(np.var(counts)))
    return np.array(out)
