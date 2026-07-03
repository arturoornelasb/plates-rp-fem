#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectral statistics and symmetry-sector classification.

Reference values (Atas et al. 2013): <r> Poisson = 0.3863, GOE = 0.5307.
Sector labels follow the repo convention: parity about the plate midlines,
'ee' = even in x and y, ..., matching the Legendre-Ritz sectors.
"""
import numpy as np

R_POISSON = 0.3863
R_GOE = 0.5307
SECTORS = ["ee", "eo", "oe", "oo"]


def r_values(ev):
    """Consecutive spacing ratios r_n = min(s_n, s_n+1)/max(...) of a sorted set."""
    s = np.diff(np.sort(np.asarray(ev, float)))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r[np.isfinite(r)]


def mean_r(ev, skip_low=10):
    """<r> +/- sem of a level sequence, dropping the lowest skip_low levels."""
    ev = np.sort(np.asarray(ev, float))[skip_low:]
    rv = r_values(ev)
    if len(rv) < 2:
        return np.nan, np.nan, 0
    return float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))), len(rv)


def probe_operators(basis, a, b, npx=48, npy=30, margin=0.02):
    """Sparse evaluation operators on an interior grid and its two mirror images
    about the midlines x = a/2, y = b/2. Used for parity classification."""
    xs = np.linspace(margin * a, (1 - margin) * a, npx)
    ys = np.linspace(margin * b, (1 - margin) * b, npy)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.vstack([X.ravel(), Y.ravel()])
    P = basis.probes(pts)
    Pmx = basis.probes(np.vstack([a - pts[0], pts[1]]))
    Pmy = basis.probes(np.vstack([pts[0], b - pts[1]]))
    return P, Pmx, Pmy


def classify_parity(V, P, Pmx, Pmy, ambiguous_below=0.8):
    """Classify each eigenvector column of V into ee/eo/oe/oo by the correlation
    of its point values with the mirrored point values. Returns
    (labels list, cx array, cy array); label 'xx' if |corr| < ambiguous_below."""
    W = P @ V
    Wx = Pmx @ V
    Wy = Pmy @ V
    nrm = np.einsum("ij,ij->j", W, W)
    cx = np.einsum("ij,ij->j", W, Wx) / np.sqrt(nrm * np.einsum("ij,ij->j", Wx, Wx))
    cy = np.einsum("ij,ij->j", W, Wy) / np.sqrt(nrm * np.einsum("ij,ij->j", Wy, Wy))
    labels = []
    for i in range(V.shape[1]):
        if min(abs(cx[i]), abs(cy[i])) < ambiguous_below:
            labels.append("xx")
        else:
            labels.append(("e" if cx[i] > 0 else "o") + ("e" if cy[i] > 0 else "o"))
    return labels, cx, cy


def _sym_diag_parity(Gx, Gy):
    """Simultaneously diagonalize the (commuting) reflection matrices restricted
    to a near-degenerate cluster subspace. Returns (rotation U, sx, sy) where
    sx, sy are the diagonal parity values (should be ~ +/-1)."""
    Sx = 0.5 * (Gx + Gx.T)
    wx, U = np.linalg.eigh(Sx)
    # within each x-parity eigenspace, diagonalize the y-reflection block
    U2 = np.zeros_like(U)
    for sign in (-1.0, 1.0):
        idx = np.where(np.sign(wx) == sign)[0]
        if len(idx) == 0:
            continue
        Ug = U[:, idx]
        Syb = Ug.T @ (0.5 * (Gy + Gy.T)) @ Ug
        _, Q = np.linalg.eigh(0.5 * (Syb + Syb.T))
        U2[:, idx] = Ug @ Q
    sx = np.einsum("ik,ij,jk->k", U2, 0.5 * (Gx + Gx.T), U2)
    sy = np.einsum("ik,ij,jk->k", U2, 0.5 * (Gy + Gy.T), U2)
    return U2, sx, sy


def classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M, corr_thresh=0.95,
                             cluster_frac=0.35, quality_min=0.8):
    """Parity classification robust to eigenvector mixing inside near-degenerate
    clusters (the eigensolver may return arbitrary rotations there).

    First pass: plain mirror-correlation labels. Modes whose |corr| falls below
    corr_thresh are grouped with neighbors closer than cluster_frac x local
    spacing; within each cluster the reflection operators are simultaneously
    diagonalized (they commute with H), the rotated vectors are re-assigned to
    the cluster's eigenvalue slots by Rayleigh quotient, and labels are read
    off the diagonalized parities. Only if the resolved parity quality is still
    below quality_min does a mode keep the label 'xx'.

    Returns (labels, quality, n_resolved) -- quality[i] = min(|sx|,|sy|)."""
    labels, cx, cy = classify_parity(V, P, Pmx, Pmy, ambiguous_below=corr_thresh)
    quality = np.minimum(np.abs(cx), np.abs(cy))
    d = local_spacing(np.asarray(lam))
    n = len(lam)
    W, Wx, Wy = P @ V, Pmx @ V, Pmy @ V
    n_resolved = 0
    i = 0
    while i < n:
        if labels[i] != "xx":
            i += 1
            continue
        # grow the cluster around i
        lo = i
        while lo > 0 and (lam[lo] - lam[lo - 1]) < cluster_frac * d[lo]:
            lo -= 1
        hi = i
        while hi < n - 1 and (lam[hi + 1] - lam[hi]) < cluster_frac * d[hi]:
            hi += 1
        cl = np.arange(lo, hi + 1)
        if len(cl) < 2:
            i += 1
            continue
        Wc, Wxc, Wyc = W[:, cl], Wx[:, cl], Wy[:, cl]
        nrm = np.sqrt(np.einsum("ij,ij->j", Wc, Wc))
        Wc, Wxc, Wyc = Wc / nrm, Wxc / nrm, Wyc / nrm
        U, sx, sy = _sym_diag_parity(Wc.T @ Wxc, Wc.T @ Wyc)
        # re-assign rotated vectors to eigenvalue slots by Rayleigh quotient
        Vc = V[:, cl]
        Hc = Vc.T @ (K @ Vc)
        Mc = Vc.T @ (M @ Vc)
        rq = np.einsum("ik,ij,jk->k", U, Hc, U) / np.einsum("ik,ij,jk->k", U, Mc, U)
        slot = np.argsort(np.argsort(rq))       # rotated col -> eigenvalue rank
        for col, s in enumerate(slot):
            q = min(abs(sx[col]), abs(sy[col]))
            j = cl[s]
            if q >= quality_min:
                labels[j] = (("e" if sx[col] > 0 else "o")
                             + ("e" if sy[col] > 0 else "o"))
                quality[j] = q
                n_resolved += 1
        i = hi + 1
    return labels, quality, n_resolved


def local_spacing(ref, window=11):
    """Smoothed local mean spacing of a sorted reference spectrum."""
    d = np.gradient(ref)
    kern = np.ones(window) / window
    pad = np.pad(d, (window // 2, window // 2), mode="edge")
    return np.convolve(pad, kern, mode="valid")[:len(ref)]


def n_star(lam, ref, frac=0.1):
    """Largest n such that ALL modes 1..n have |lam - ref| < frac * local spacing."""
    n = min(len(lam), len(ref))
    ok = np.abs(lam[:n] - ref[:n]) < frac * local_spacing(ref)[:n]
    bad = np.where(~ok)[0]
    return int(bad[0]) if len(bad) else n
