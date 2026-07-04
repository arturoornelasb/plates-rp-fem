#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C0 interior-penalty discretization of the Kirchhoff plate (Gap A
reconciliation instrument, task #0).

Why this exists: the Argyris element's numerically-built basis breaks down
below h ~ 0.008 (measured, E1), and skfem's Lagrange elements do not carry
second derivatives -- so large-N true-operator eigenvector windows need a
custom discretization. This module implements P_k Lagrange (k = 3 or 4) on
affine triangles with analytic reference hessians, assembling the
Brenner-Sung C0-IP bilinear form

  a_h(u,v) = sum_T int [(1-nu) H(u):H(v) + nu (lap u)(lap v)]
           - sum_e int_e {M_nn(u)} [[dn v]] + {M_nn(v)} [[dn u]]
           + sum_e (sigma k^2 / h_e) int_e [[dn u]] [[dn v]] ,

with e running over INTERIOR edges only: free edges are the natural
boundary conditions (no boundary facet terms), simply supported edges are
the essential condition w = 0 on boundary nodes. M_nn = (1-nu) n.H.n +
nu lap. skfem's MeshTri supplies mesh topology (t2f, f2t, element_finder);
all element calculus is done here.

Validation gates are E1-style (exact SSSS; certified-Argyris FFFF) and a
penalty-robustness sweep; see experiments/e14_c0ip_gapA.
"""
import numpy as np
import scipy.sparse as sp


# ---------------- reference P_k basis with analytic derivatives ----------------
def _lattice(k):
    pts = []
    for i in range(k + 1):
        for j in range(k + 1 - i):
            pts.append((i / k, j / k))
    return np.array(pts)                      # (nb, 2)


def _monomials(k):
    return [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]


class RefBasis:
    """P_k Lagrange basis on the reference triangle with value / grad /
    hessian evaluation from exact monomial coefficients."""

    def __init__(self, k):
        self.k = k
        self.nodes = _lattice(k)
        self.mono = _monomials(k)
        nb = len(self.nodes)
        V = np.array([[x ** a * y ** b for (a, b) in self.mono]
                      for (x, y) in self.nodes])
        self.C = np.linalg.solve(V, np.eye(nb))   # C[:, i] = coeffs of phi_i
        self.nb = nb

    def _eval_mono(self, X, da, db):
        """d^da/dx d^db/dy of each monomial at points X (2, nq)."""
        out = np.zeros((len(self.mono), X.shape[1]))
        for m, (a, b) in enumerate(self.mono):
            if a < da or b < db:
                continue
            ca = np.prod(np.arange(a, a - da, -1)) if da else 1.0
            cb = np.prod(np.arange(b, b - db, -1)) if db else 1.0
            out[m] = ca * cb * X[0] ** (a - da) * X[1] ** (b - db)
        return out

    def eval(self, X):
        """Returns (val (nb,nq), grad (nb,2,nq), hess (nb,3,nq) [xx,xy,yy])."""
        val = self.C.T @ self._eval_mono(X, 0, 0)
        gx = self.C.T @ self._eval_mono(X, 1, 0)
        gy = self.C.T @ self._eval_mono(X, 0, 1)
        hxx = self.C.T @ self._eval_mono(X, 2, 0)
        hxy = self.C.T @ self._eval_mono(X, 1, 1)
        hyy = self.C.T @ self._eval_mono(X, 0, 2)
        return val, np.stack([gx, gy], 1), np.stack([hxx, hxy, hyy], 1)


# ---------------- dof management ----------------
class C0IPSpace:
    """Global P_k dof numbering on a skfem MeshTri: vertex dofs, (k-1)
    oriented dofs per edge (ordered from the lower to the higher global
    vertex index), and interior bubble dofs per element."""

    def __init__(self, mesh, k):
        self.mesh = mesh
        self.k = k
        self.ref = RefBasis(k)
        nv = mesh.p.shape[1]
        ne = mesh.facets.shape[1]
        nt = mesh.t.shape[1]
        self.n_edge = k - 1
        self.n_int = (k - 1) * (k - 2) // 2
        self.N = nv + ne * self.n_edge + nt * self.n_int
        # classify reference lattice nodes: vertex / edge / interior
        nodes = self.ref.nodes
        vid = {(0.0, 0.0): 0, (1.0, 0.0): 1, (0.0, 1.0): 2}
        # reference vertices in barycentric terms: v0=(0,0), v1=(1,0), v2=(0,1)
        # edges (skfem MeshTri t2f order): f0=(v0,v1), f1=(v1,v2), f2=(v0,v2)
        self.node_kind = []
        for (x, y) in nodes:
            key = (round(x, 12), round(y, 12))
            if key in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]:
                self.node_kind.append(("v", [(0.0, 0.0), (1.0, 0.0),
                                             (0.0, 1.0)].index(key)))
            elif abs(y) < 1e-12:
                self.node_kind.append(("e", 0, x))          # on v0-v1, param x
            elif abs(x + y - 1.0) < 1e-12:
                self.node_kind.append(("e", 1, y))          # on v1-v2, param y
            elif abs(x) < 1e-12:
                self.node_kind.append(("e", 2, y))          # on v0-v2, param y
            else:
                self.node_kind.append(("i", None))
        # build element dof tables
        t = mesh.t
        t2f = mesh.t2f
        facets = mesh.facets
        self.elem_dofs = np.zeros((nt, self.ref.nb), dtype=np.int64)
        int_count = 0
        edge_local = {0: (0, 1), 1: (1, 2), 2: (0, 2)}
        int_seen = {}
        for e in range(nt):
            verts = t[:, e]
            n_int_used = 0
            for ln, kind in enumerate(self.node_kind):
                if kind[0] == "v":
                    self.elem_dofs[e, ln] = verts[kind[1]]
                elif kind[0] == "e":
                    le, param = kind[1], kind[2]
                    f = t2f[le, e]
                    va, vb = verts[list(edge_local[le])]
                    # orientation: dof slots ordered from min to max global id
                    fa, fb = facets[:, f]
                    lo = min(fa, fb)
                    # param measured from local va; slot index along the edge
                    idx = int(round(param * k)) - 1          # 0..k-2 from va
                    if va != lo:
                        idx = (k - 2) - idx
                    self.elem_dofs[e, ln] = nv + f * self.n_edge + idx
                else:
                    self.elem_dofs[e, ln] = nv + ne * self.n_edge \
                        + e * self.n_int + n_int_used
                    n_int_used += 1
        # affine maps
        p = mesh.p
        v0 = p[:, t[0]]
        self.A = np.stack([p[:, t[1]] - v0, p[:, t[2]] - v0], axis=2)  # (2, nt, 2)
        self.A = np.transpose(self.A, (1, 0, 2))                       # (nt,2,2)
        self.b0 = v0.T                                                 # (nt,2)
        self.detA = (self.A[:, 0, 0] * self.A[:, 1, 1]
                     - self.A[:, 0, 1] * self.A[:, 1, 0])
        Ainv = np.empty_like(self.A)
        Ainv[:, 0, 0] = self.A[:, 1, 1]
        Ainv[:, 0, 1] = -self.A[:, 0, 1]
        Ainv[:, 1, 0] = -self.A[:, 1, 0]
        Ainv[:, 1, 1] = self.A[:, 0, 0]
        self.Ainv = Ainv / self.detA[:, None, None]                    # A^{-1}
        self.boundary_facets = np.nonzero(mesh.f2t[1] == -1)[0] \
            if (mesh.f2t[1] == -1).any() else \
            np.nonzero(mesh.f2t[1] == mesh.f2t[0])[0]

    # ---------- point evaluation (probes) ----------
    def probes(self, pts):
        finder = self.mesh.element_finder()
        el = finder(pts[0], pts[1])
        loc = np.einsum("nij,nj->ni", self.Ainv[el],
                        (pts.T - self.b0[el]))
        rows, cols, vals = [], [], []
        for i in range(pts.shape[1]):
            v, _, _ = self.ref.eval(loc[i][:, None])
            rows.extend([i] * self.ref.nb)
            cols.extend(self.elem_dofs[el[i]].tolist())
            vals.extend(v[:, 0].tolist())
        return sp.csr_matrix((vals, (rows, cols)),
                             shape=(pts.shape[1], self.N))


def _tri_quadrature(order):
    """Tensor Gauss rule on the reference triangle via the Duffy transform."""
    n = order // 2 + 2
    x, wx = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (x + 1.0)
    wx = 0.5 * wx
    X, Y = np.meshgrid(x, x, indexing="ij")
    W = np.outer(wx, wx) * (1.0 - X)
    pts = np.vstack([X.ravel(), (Y * (1.0 - X)).ravel()])
    return pts, W.ravel()


def assemble_c0ip(mesh, k=4, nu=0.33, sigma_factor=10.0):
    """Assemble (K, M, space) for the free Kirchhoff plate. Boundary
    conditions: none imposed here (free edges natural); SS via essential
    w = 0 on boundary dofs (see boundary_dofs)."""
    space = C0IPSpace(mesh, k)
    ref = space.ref
    nb, nt = ref.nb, mesh.t.shape[1]

    # ---------------- volume terms ----------------
    Xq, Wq = _tri_quadrature(2 * k + 2)
    val, grad, hess = ref.eval(Xq)                     # (nb,nq), (nb,2,nq), (nb,3,nq)
    nq = Xq.shape[1]
    Ai = space.Ainv                                    # (nt,2,2)
    # physical hessian: H_phys = Ai^T H_ref Ai  (per element, constant map)
    # build hessian transform for the 3-vector (xx, xy, yy)
    a11, a12 = Ai[:, 0, 0], Ai[:, 0, 1]
    a21, a22 = Ai[:, 1, 0], Ai[:, 1, 1]
    # H_phys components from H_ref = [hxx, hxy, hyy]:
    T = np.zeros((nt, 3, 3))
    T[:, 0, 0] = a11 * a11; T[:, 0, 1] = 2 * a11 * a21; T[:, 0, 2] = a21 * a21
    T[:, 1, 0] = a11 * a12; T[:, 1, 1] = a11 * a22 + a12 * a21
    T[:, 1, 2] = a21 * a22
    T[:, 2, 0] = a12 * a12; T[:, 2, 1] = 2 * a12 * a22; T[:, 2, 2] = a22 * a22
    Hp = np.einsum("tab,ibq->tiaq", T, hess)           # (nt, nb, 3, nq)
    lap = Hp[:, :, 0, :] + Hp[:, :, 2, :]
    dA = np.abs(space.detA)
    wq = Wq[None, :] * dA[:, None]                     # (nt, nq)
    # bending: (1-nu)(Hxx Hxx' + 2 Hxy Hxy' + Hyy Hyy') + nu lap lap'
    w3 = np.array([1.0, 2.0, 1.0])
    Ke = (1.0 - nu) * np.einsum("tiaq,a,tjaq,tq->tij", Hp, w3, Hp, wq,
                                optimize=True) \
        + nu * np.einsum("tiq,tjq,tq->tij", lap, lap, wq, optimize=True)
    Me = np.einsum("iq,jq,tq->tij", val, val, wq, optimize=True)

    rows = np.repeat(space.elem_dofs, nb, axis=1).ravel()
    cols = np.tile(space.elem_dofs, (1, nb)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(space.N, space.N))
    M = sp.csr_matrix((Me.ravel(), (rows, cols)), shape=(space.N, space.N))

    # ---------------- interior-facet terms ----------------
    f2t = mesh.f2t
    interior = np.nonzero(f2t[1] != -1)[0] if (f2t[1] == -1).any() else \
        np.nonzero(f2t[1] != f2t[0])[0]
    facets = mesh.facets
    p = mesh.p
    xq1, wq1 = np.polynomial.legendre.leggauss(k + 3)
    xq1 = 0.5 * (xq1 + 1.0)
    wq1 = 0.5 * wq1
    nq1 = len(xq1)
    sigma = sigma_factor * k * k

    rows_f, cols_f, vals_f = [], [], []
    for f in interior:
        va, vb = facets[:, f]
        pa, pb = p[:, va], p[:, vb]
        h_e = float(np.hypot(*(pb - pa)))
        tang = (pb - pa) / h_e
        nrm = np.array([tang[1], -tang[0]])
        pts_e = pa[:, None] + (pb - pa)[:, None] * xq1[None, :]
        t1, t2 = f2t[0, f], f2t[1, f]
        # orient normal outward from t1
        c1 = p[:, mesh.t[:, t1]].mean(axis=1)
        if np.dot(nrm, pa + 0.5 * (pb - pa) - c1) < 0:
            nrm = -nrm
        # evaluate both sides at the edge points
        side = []
        for tt in (t1, t2):
            loc = space.Ainv[tt] @ (pts_e - space.b0[tt][:, None])
            v_, g_, h_ = ref.eval(loc)
            # physical grad / hessian
            gp = np.einsum("ba,ibq->iaq", space.Ainv[tt], g_)
            hp = np.einsum("ab,ibq->iaq", T[tt], h_)
            dn = nrm[0] * gp[:, 0, :] + nrm[1] * gp[:, 1, :]
            nHn = (nrm[0] * nrm[0] * hp[:, 0, :]
                   + 2 * nrm[0] * nrm[1] * hp[:, 1, :]
                   + nrm[1] * nrm[1] * hp[:, 2, :])
            lapp = hp[:, 0, :] + hp[:, 2, :]
            Mnn = (1.0 - nu) * nHn + nu * lapp
            side.append((dn, Mnn))
        dofs = np.concatenate([space.elem_dofs[t1], space.elem_dofs[t2]])
        dn_j = np.concatenate([side[0][0], -side[1][0]])       # jump
        Mnn_a = 0.5 * np.concatenate([side[0][1], side[1][1]])  # average
        wq_e = wq1 * h_e
        blk = (- np.einsum("iq,jq,q->ij", Mnn_a, dn_j, wq_e)
               - np.einsum("iq,jq,q->ij", dn_j, Mnn_a, wq_e)
               + (sigma / h_e) * np.einsum("iq,jq,q->ij", dn_j, dn_j, wq_e))
        rows_f.append(np.repeat(dofs, 2 * nb))
        cols_f.append(np.tile(dofs, 2 * nb))
        vals_f.append(blk.ravel())
    Kf = sp.csr_matrix((np.concatenate(vals_f),
                        (np.concatenate(rows_f), np.concatenate(cols_f))),
                       shape=(space.N, space.N))
    K = (K + Kf).tocsc()
    return K, M.tocsc(), space


def boundary_dofs(space):
    """Dofs (vertex + edge) lying on the boundary -- the SS essential set."""
    mesh = space.mesh
    nv = mesh.p.shape[1]
    out = []
    for f in space.boundary_facets:
        va, vb = mesh.facets[:, f]
        out.extend([int(va), int(vb)])
        out.extend(range(nv + f * space.n_edge,
                         nv + (f + 1) * space.n_edge))
    return np.unique(np.array(out, dtype=np.int64))
