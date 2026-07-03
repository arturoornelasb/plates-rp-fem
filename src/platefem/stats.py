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


def triangle_probe_operators(basis, L=1.0, npts=1500, seed=11):
    """Evaluation operators for C3v classification on the equilateral triangle
    of platefem.kirchhoff.triangle_basis (centroid at origin, vertex on +y).
    Random interior points X plus their images under the 120-degree rotation
    (R^-1 X, so P_R evaluates the rotated function) and the reflection about
    the y-axis. Returns (P, P_R, P_s)."""
    rng = np.random.default_rng(seed)
    verts = np.array([[-0.5 * L, -np.sqrt(3) / 6 * L],
                      [0.5 * L, -np.sqrt(3) / 6 * L],
                      [0.0, np.sqrt(3) / 3 * L]])
    bary = rng.dirichlet([1.0, 1.0, 1.0], size=4 * npts)
    bary = bary[bary.min(axis=1) > 0.02][:npts]
    pts = (bary @ verts).T                       # (2, npts)
    th = -2.0 * np.pi / 3.0                      # R^-1 = rotation by -120 deg
    Rinv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    P = basis.probes(pts)
    P_R = basis.probes(Rinv @ pts)
    P_s = basis.probes(np.vstack([-pts[0], pts[1]]))
    return P, P_R, P_s


def classify_c3v(lam, V, P, P_R, P_s, K, M, thresh=0.7, cluster_frac=0.35,
                 quality_min=0.8):
    """C3v irrep classification (equilateral triangle): labels 'A1', 'A2', 'E'
    (or 'xx' if unresolved). Per-mode characters c_R = <w, Rw>/<w,w> and
    c_s = <w, sw>/<w,w>: c_R = +1 for A-type (then c_s = +1 -> A1, -1 -> A2)
    and c_R = -1/2 for ANY vector inside an E doublet (rotation restricted to
    the 2D E space has trace -1). Ambiguous modes (accidental A/E clusters,
    solver mixing) are resolved by diagonalizing the symmetrized R-Gram inside
    the near-degenerate cluster -- eigenvalues split {~+1} (A) from {~-1/2}
    (E) -- then the reflection Gram inside the A part; rotated vectors are
    re-assigned to eigenvalue slots by Rayleigh quotient.

    Returns (labels, quality, n_resolved); quality is the character margin."""
    W, WR, Ws = P @ V, P_R @ V, P_s @ V
    nrm = np.einsum("ij,ij->j", W, W)
    cR = np.einsum("ij,ij->j", W, WR) / nrm
    cS = np.einsum("ij,ij->j", W, Ws) / nrm

    labels, quality = [], np.zeros(len(lam))
    for i in range(len(lam)):
        if cR[i] > thresh and cS[i] > thresh:
            labels.append("A1"); quality[i] = min(cR[i], cS[i])
        elif cR[i] > thresh and cS[i] < -thresh:
            labels.append("A2"); quality[i] = min(cR[i], -cS[i])
        elif abs(cR[i] + 0.5) < 0.5 * (1.0 - thresh):
            labels.append("E"); quality[i] = 1.0 - 2.0 * abs(cR[i] + 0.5)
        else:
            labels.append("xx")

    d = local_spacing(np.asarray(lam))
    n = len(lam)
    n_resolved = 0
    i = 0
    while i < n:
        if labels[i] != "xx":
            i += 1
            continue
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
        Wc = W[:, cl] / np.sqrt(nrm[cl])
        GR = Wc.T @ (WR[:, cl] / np.sqrt(nrm[cl]))
        GS = Wc.T @ (Ws[:, cl] / np.sqrt(nrm[cl]))
        GR = 0.5 * (GR + GR.T)
        wR, U = np.linalg.eigh(GR)
        U2 = U.copy()
        lab_cols = []
        a_idx = np.where(wR > 0.0)[0]            # ~+1 -> A-type
        e_idx = np.where(wR <= 0.0)[0]           # ~-1/2 -> E
        if len(a_idx):
            Ua = U[:, a_idx]
            Sb = Ua.T @ (0.5 * (GS + GS.T)) @ Ua
            _, Q = np.linalg.eigh(0.5 * (Sb + Sb.T))
            U2[:, a_idx] = Ua @ Q
        sR = np.einsum("ik,ij,jk->k", U2, GR, U2)
        sS = np.einsum("ik,ij,jk->k", U2, 0.5 * (GS + GS.T), U2)
        for col in range(len(cl)):
            if sR[col] > quality_min:
                q = min(sR[col], abs(sS[col]))
                lab_cols.append(("A1" if sS[col] > 0 else "A2", q))
            elif abs(sR[col] + 0.5) < 0.5 * (1.0 - quality_min):
                lab_cols.append(("E", 1.0 - 2.0 * abs(sR[col] + 0.5)))
            else:
                lab_cols.append(("xx", 0.0))
        Vc = V[:, cl]
        Hc = Vc.T @ (K @ Vc)
        Mc = Vc.T @ (M @ Vc)
        rq = np.einsum("ik,ij,jk->k", U2, Hc, U2) / \
            np.einsum("ik,ij,jk->k", U2, Mc, U2)
        slot = np.argsort(np.argsort(rq))
        for col, s in enumerate(slot):
            lab, q = lab_cols[col]
            j = cl[s]
            if lab != "xx" and q >= quality_min:
                labels[j] = lab
                quality[j] = q
                n_resolved += 1
        i = hi + 1
    return labels, quality, n_resolved


def dedupe_doublets(lam, labels, sector="E", rel_gap=1e-5):
    """Collapse exact doublets (consecutive same-sector levels closer than
    rel_gap x Lambda) to one level per pair (the mean). Returns (levels,
    n_pairs, n_odd, max_split) -- n_odd counts unpaired leftovers (should be
    0 or edge-of-window truncation)."""
    ev = [(l, i) for i, (l, s) in enumerate(zip(lam, labels)) if s == sector]
    out, splits = [], []
    i = 0
    n_pairs = n_odd = 0
    while i < len(ev):
        if i + 1 < len(ev) and abs(ev[i + 1][0] - ev[i][0]) < rel_gap * abs(ev[i][0]):
            out.append(0.5 * (ev[i][0] + ev[i + 1][0]))
            splits.append(abs(ev[i + 1][0] - ev[i][0]) / abs(ev[i][0]))
            n_pairs += 1
            i += 2
        else:
            out.append(ev[i][0])
            n_odd += 1
            i += 1
    return (np.array(out), n_pairs, n_odd,
            float(max(splits)) if splits else 0.0)


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
