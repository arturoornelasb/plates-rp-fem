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


def collar_pullback(pts, delta, a=A_AX, b=B_AX, p=P_EXP):
    """Radially pull points with superellipse coordinate sigma > 1-delta
    back to sigma = 1-delta (chord-vs-curve sliver handling; README)."""
    sig = (np.abs(pts[0] / a) ** p + np.abs(pts[1] / b) ** p) ** (1.0 / p)
    scale = np.where(sig > 1.0 - delta, (1.0 - delta) / np.maximum(sig, 1e-300),
                     1.0)
    return pts * scale


def sag_delta(mesh_coarse, a=A_AX, b=B_AX, p=P_EXP, safety=1.5):
    """Collar depth self-calibrated to the coarse mesh: the worst
    chord-vs-curve sag of its boundary facets (superellipse-coordinate
    deficit of chord midpoints), times a safety factor. Amends the fixed
    COLLAR = 2e-3 of the first freeze (documented: the sag is
    level-dependent; measured, not assumed)."""
    bf = mesh_coarse.boundary_facets()
    ends = mesh_coarse.facets[:, bf]
    mid = 0.5 * (mesh_coarse.p[:, ends[0]] + mesh_coarse.p[:, ends[1]])
    sig = (np.abs(mid[0] / a) ** p + np.abs(mid[1] / b) ** p) ** (1.0 / p)
    return float(safety * max(1.0 - sig.min(), 1e-12))


def dof_points(space):
    """Physical coordinates of every global dof of a C0IPSpace (nodal
    P_k lattice: vertices; k-1 equispaced points per facet ordered from
    the lower to the higher global vertex index; interior lattice nodes
    per element in node_kind order). Verified in-run by the nodal
    self-test (probes at dof points ~ identity)."""
    mesh, k = space.mesh, space.k
    nv = mesh.p.shape[1]
    ne = mesh.facets.shape[1]
    pts = np.zeros((2, space.N))
    pts[:, :nv] = mesh.p
    fa, fb = mesh.facets[0], mesh.facets[1]
    lo = np.minimum(fa, fb)
    hi = np.maximum(fa, fb)
    # edge dofs facet-major: dof = nv + f*n_edge + idx, param from lo
    for idx in range(space.n_edge):
        t = (idx + 1) / k
        pts[:, nv + idx:nv + ne * space.n_edge:space.n_edge] = \
            (1.0 - t) * mesh.p[:, lo] + t * mesh.p[:, hi]
    int_nodes = np.array([space.ref.nodes[ln]
                          for ln, kind in enumerate(space.node_kind)
                          if kind[0] == "i"])            # (n_int, 2)
    if len(int_nodes):
        base = nv + ne * space.n_edge
        # physical = A[e] @ node + b0[e], element-major blocks of n_int
        phys = np.einsum("eij,nj->eni", space.A, int_nodes) \
            + space.b0[:, None, :]                       # (nt, n_int, 2)
        pts[:, base:] = phys.reshape(-1, 2).T
    return pts


def tol_finder(space, pts, tol=1e-10):
    """Tolerant point location (skfem's element_finder rejects points
    exactly ON element boundaries -- i.e. every nodal dof point -- on
    float noise). KD-tree of element centroids, nearest-K candidates,
    accept the candidate with the best min-barycentric margin >= -tol."""
    from scipy.spatial import cKDTree
    mesh = space.mesh
    cent = mesh.p[:, mesh.t].mean(axis=1)                # (2, nt)
    tree = cKDTree(cent.T)
    npts = pts.shape[1]
    el = np.full(npts, -1, dtype=np.int64)
    margin = np.full(npts, -np.inf)
    todo = np.arange(npts)
    for K in (8, 64, 512):
        K = min(K, mesh.t.shape[1])
        _, cand = tree.query(pts[:, todo].T, k=K)
        cand = np.atleast_2d(cand)
        loc = np.einsum("nkij,nkj->nki", space.Ainv[cand],
                        pts[:, todo].T[:, None, :] - space.b0[cand])
        marg = np.minimum(np.minimum(loc[..., 0], loc[..., 1]),
                          1.0 - loc[..., 0] - loc[..., 1])
        best = np.argmax(marg, axis=1)
        r = np.arange(len(todo))
        el[todo] = cand[r, best]
        margin[todo] = marg[r, best]
        todo = todo[margin[todo] < -tol]
        if not len(todo):
            break
    if len(todo):
        raise ValueError(f"{len(todo)} points outside the mesh "
                         f"(worst margin {margin[todo].min():.2e})")
    return el


def vec_probes(space, pts):
    """Vectorized point evaluation (same math as C0IPSpace.probes, one
    ref.eval call instead of a per-point Python loop; tolerant finder)."""
    import scipy.sparse as sp
    el = tol_finder(space, pts)
    loc = np.einsum("nij,nj->ni", space.Ainv[el], (pts.T - space.b0[el]))
    v, _, _ = space.ref.eval(loc.T)                      # (nb, npts)
    npts = pts.shape[1]
    rows = np.repeat(np.arange(npts), space.ref.nb)
    cols = space.elem_dofs[el].ravel()
    vals = v.T.ravel()
    return sp.csr_matrix((vals, (rows, cols)), shape=(npts, space.N))
