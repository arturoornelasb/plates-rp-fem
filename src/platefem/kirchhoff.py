#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Kirchhoff plate FEM machinery (scikit-fem), validated in experiments/e01_validation.

Conventions: flexural eigenproblem D lap^2 w = rho h omega^2 w with D = rho h = 1;
physical eigenvalue Lambda = omega^2. Free edges are the natural BCs of the
bending energy (no constraints); a Winkler edge spring kappa*w*v on the boundary
gives V_n + kappa w = 0 (Protocol B), recovering simply supported as kappa -> inf.

Validated facts to respect (see e01 RESULTS.md):
- ElementGlobal subclasses (Morley, Argyris) cache mesh-dependent arrays on the
  instance: ALWAYS instantiate a fresh element per mesh / per FacetBasis.
- The Argyris near-zero (rigid) cluster is polluted at fine mesh (~1e-2 absolute
  at 128x80): detect rigid modes by spectral gap, never fixed tolerance.
- Uniform Argyris meshes beyond ~128x80 are unusable (h^-5 Vandermonde).
"""
import numpy as np
from scipy.sparse.linalg import eigsh

from skfem import Basis, BilinearForm, FacetBasis, MeshTri
from skfem.element import ElementTriArgyris, ElementTriMorley  # noqa: F401 (re-export)
from skfem.helpers import dd, ddot, trace


def make_forms(nu):
    """Bending energy a(w,v) = int [(1-nu) w_,ij v_,ij + nu lap(w) lap(v)] and mass."""
    @BilinearForm
    def bending(u, v, w):
        return (1.0 - nu) * ddot(dd(u), dd(v)) + nu * trace(dd(u)) * trace(dd(v))

    @BilinearForm
    def mass(u, v, w):
        return u * v

    return bending, mass


@BilinearForm
def boundary_mass(u, v, w):
    return u * v


def rectangle_basis(nx, ny, a, b, element_cls=ElementTriArgyris):
    xs = np.linspace(0.0, a, nx + 1)
    ys = np.linspace(0.0, b, ny + 1)
    mesh = MeshTri.init_tensor(xs, ys)
    return mesh, Basis(mesh, element_cls())


def assemble_plate(mesh, basis, nu):
    bending, mass = make_forms(nu)
    return bending.assemble(basis).tocsc(), mass.assemble(basis).tocsc()


def boundary_matrix(mesh, element_cls=ElementTriArgyris):
    """Winkler boundary mass B; add kappa*B to K for the edge-spring plate."""
    fb = FacetBasis(mesh, element_cls())
    return boundary_mass.assemble(fb).tocsc()


def solve_lowest(K, M, k, sigma=-10.0):
    """Lowest k eigenVALUES of K x = Lambda M x (shift-invert Lanczos).
    The no-vector path of scipy's eigsh is reliable for these matrices and is
    the source of record for eigenvalues; for eigenpairs use solve_modes."""
    return np.sort(eigsh(K, k=k, M=M, sigma=sigma, which="LM",
                         return_eigenvectors=False))


def _ritz_project(Y, K, M, k_min):
    """Whitened Rayleigh-Ritz projection of a block Y. M is numerically
    semidefinite (cond ~ 1e14 for Argyris), so Br = Y'MY is whitened by its
    own eigendecomposition with truncation of the numerically-null subspace
    rather than passed to a generalized eigh. Returns (theta, Yritz)."""
    Y, _ = np.linalg.qr(Y)
    Ar = Y.T @ (K @ Y)
    Br = Y.T @ (M @ Y)
    Ar, Br = 0.5 * (Ar + Ar.T), 0.5 * (Br + Br.T)
    eb, Ub = np.linalg.eigh(Br)
    keep = eb > 1e-12 * eb[-1]
    if int(np.sum(keep)) < k_min:
        raise RuntimeError(f"projected mass matrix kept {int(np.sum(keep))} "
                           f"of {Y.shape[1]} directions (< {k_min})")
    T = Ub[:, keep] / np.sqrt(eb[keep])
    th, Wr = np.linalg.eigh(T.T @ Ar @ T)
    return th, Y @ (T @ Wr)


def solve_modes(K, M, k, sigma=-10.0, band=128, buffer=64, sweeps_max=25,
                seed=7, match_frac=0.5, resid_sanity=1e-2, Y0=None):
    """Certified lowest-k eigenPAIRS of K x = Lambda M x.

    Eigenvalues of record come from the no-vector ARPACK path (solve_lowest),
    which is reliable for these matrices. Eigenvectors are computed WITHOUT
    ARPACK: scipy 1.18.0's return_eigenvectors=True shift-invert path is
    broken for the ill-conditioned Argyris matrices (eigenvalues shifted by
    exactly -sigma when it converges, misconverged Ritz pairs or ARPACK
    error -8 otherwise; see experiments/e01 and the E2 work log).

    Method: BANDED shift-invert subspace iteration WITH DEFLATION. The
    certified spectrum is partitioned into index bands of ~`band` modes, cuts
    placed at locally large gaps (never through a near-degenerate cluster);
    each band gets its own factorization of (K - sigma_b M) at a gap midpoint
    near the band center. Bands are processed bottom-up and every sweep is
    M-orthogonalized against all previously converged vectors: without this,
    an interior shift's transfer function 1/(Lambda - sigma) is dominated by
    the ENTIRE spectrum below sigma (|1/(Lambda-sigma)| >= 1/sigma there), so
    the block wastes itself on already-known low modes and the band's upper
    tail never converges (observed as a stuck ~50% mismatch floor -- see the
    E2 v2 work log). A single global shift fails for spectra spanning many
    decades for the analogous reason.

    Certification per band: block Ritz values must contain the certified
    eigenvalues 1:1 within match_frac x local mean spacing, and the selected
    columns must pass the residual sanity gate on elastic modes. These vectors
    are used for LABELING (parity classification), never for spacing
    statistics; mixing inside near-degenerate clusters is handled downstream
    by cluster-resolved classification.

    Y0 is accepted for backward compatibility and ignored (bands re-seed
    deterministically). Returns (lam, X, info, None).
    """
    from scipy.sparse.linalg import splu
    from .stats import local_spacing

    lam = solve_lowest(K, M, k, sigma)
    d = np.maximum(local_spacing(lam), 1e-8 * abs(lam[-1]))
    elastic_all = np.abs(lam) > 1e-6 * abs(lam[-1])
    n = K.shape[0]
    rng = np.random.default_rng(seed)

    # band cut points at locally large gaps (avoid splitting clusters)
    cuts = [0]
    while cuts[-1] + band < k:
        target = cuts[-1] + band
        lo, hi = max(cuts[-1] + band // 2, target - 12), min(k - 1, target + 12)
        gaps = np.diff(lam[lo:hi + 1])
        cuts.append(lo + int(np.argmax(gaps)) + 1)
    cuts.append(k)

    X = np.empty((n, k))
    MXc = np.empty((n, k))                       # M @ converged vectors (deflation)
    infos = []
    for i0, i1 in zip(cuts[:-1], cuts[1:]):
        nb = i1 - i0
        # shift: gap midpoint nearest the band center (first band: below spectrum)
        if i0 == 0:
            sigma_b = sigma
        else:
            c = min((i0 + i1) // 2, k - 2)
            sigma_b = 0.5 * (lam[c] + lam[c + 1])
        lu = splu((K - sigma_b * M).tocsc())
        m = min(n, nb + 2 * buffer)
        Y = rng.standard_normal((n, m))
        history = []
        certified = False
        for sweep in range(1, sweeps_max + 1):
            Y = lu.solve(M @ Y)
            if i0 > 0:                           # deflate converged lower bands
                Y -= X[:, :i0] @ (MXc[:, :i0].T @ Y)
            th, Y = _ritz_project(Y, K, M, nb)
            if Y.shape[1] < m:
                Y = np.hstack([Y, rng.standard_normal((n, m - Y.shape[1]))])
            # align the sorted block Ritz values with the band's certified
            # eigenvalues: search the best contiguous offset, then gate per mode
            nth = len(th)
            best_s, best_err = 0, np.inf
            for s in range(0, nth - nb + 1):
                err = float(np.sum(np.abs(th[s:s + nb] - lam[i0:i1])
                                   / np.maximum(d[i0:i1], 1e-300)))
                if err < best_err:
                    best_s, best_err = s, err
            sel = slice(best_s, best_s + nb)
            ok_match = np.all(np.abs(th[sel] - lam[i0:i1]) < match_frac * d[i0:i1])
            Xb = Y[:, sel]
            KX = K @ Xb
            MX = M @ Xb
            rq = np.einsum("ij,ij->j", Xb, KX) / np.einsum("ij,ij->j", Xb, MX)
            resid = np.linalg.norm(KX - rq * MX, axis=0) / (
                np.linalg.norm(KX, axis=0) + np.abs(rq) * np.linalg.norm(MX, axis=0))
            el = elastic_all[i0:i1]
            max_res = float(np.max(resid[el])) if el.any() else 0.0
            history.append((int(np.sum(~(np.abs(th[sel] - lam[i0:i1])
                                         < match_frac * d[i0:i1]))), max_res))
            if ok_match and max_res < resid_sanity:
                X[:, i0:i1] = Xb
                MXc[:, i0:i1] = MX
                infos.append(dict(band=[int(i0), int(i1)], sweeps=sweep,
                                  max_resid=max_res))
                certified = True
                break
        if not certified:
            raise RuntimeError(f"solve_modes: band [{i0},{i1}) not certified "
                               f"after {sweeps_max} sweeps; history "
                               f"(n_mismatch, max_resid): {history}")
        del lu
    info = dict(bands=infos,
                max_resid=float(max(b["max_resid"] for b in infos)),
                sweeps=int(max(b["sweeps"] for b in infos)))
    return lam, X, info, None


def split_rigid(lam, V=None, max_check=10, gap_factor=1e3):
    """Gap-based rigid-mode split (see module docstring). Returns
    (elastic lam, elastic V or None, n_rigid, largest rigid |Lambda|)."""
    head = np.maximum(np.abs(lam[:max_check]), 1e-300)
    ratios = head[1:] / head[:-1]
    i_gap = int(np.argmax(ratios))
    n_rigid = i_gap + 1 if ratios[i_gap] > gap_factor else 0
    rigid_max = float(np.abs(lam[n_rigid - 1])) if n_rigid else 0.0
    Ve = V[:, n_rigid:] if V is not None else None
    return lam[n_rigid:], Ve, n_rigid, rigid_max


def ssss_exact(a, b, n_modes):
    """Exact SSSS rectangle spectrum Lambda_mn = ((m pi/a)^2 + (n pi/b)^2)^2."""
    mmax = int(np.ceil(np.sqrt(n_modes) * 8))
    m = np.arange(1, mmax + 1)
    kx2 = (m * np.pi / a) ** 2
    ky2 = (m * np.pi / b) ** 2
    lam = (kx2[:, None] + ky2[None, :]) ** 2
    return np.sort(lam.ravel())[:n_modes]


def triangle_basis(nrefine, L=1.0, element_cls=ElementTriArgyris):
    """Equilateral triangle of side L, centroid at the origin, one vertex on
    the +y axis, meshed by uniform midpoint refinement of the single macro
    triangle -- the mesh is therefore C3v-symmetric by construction (E-doublet
    splittings stay at round-off level). Returns (mesh, basis, verts)."""
    verts = np.array([[-0.5 * L, -np.sqrt(3) / 6 * L],
                      [0.5 * L, -np.sqrt(3) / 6 * L],
                      [0.0, np.sqrt(3) / 3 * L]])
    mesh = MeshTri(verts.T.copy(), np.array([[0, 1, 2]]).T.copy()).refined(nrefine)
    return mesh, Basis(mesh, element_cls()), verts


def disk_basis(nrefine, element_cls=ElementTriArgyris):
    """Unit disk (R = 1) via MeshTri.init_circle (boundary vertices projected
    onto the circle at every refinement -- inscribed-polygon geometry, so the
    boundary error is O(h^2) and must be GATED against the semi-analytic disk
    reference, see platefem.disk)."""
    mesh = MeshTri.init_circle(nrefine)
    return mesh, Basis(mesh, element_cls())


def ellipse_basis(nrefine, a, b, element_cls=ElementTriArgyris):
    """Ellipse with semi-axes (a, b), centroid at origin, axes on x/y:
    the init_circle mesh scaled anisotropically (boundary vertices land
    exactly on the ellipse)."""
    m0 = MeshTri.init_circle(nrefine)
    p = m0.p.copy()
    p[0] *= a
    p[1] *= b
    mesh = MeshTri(p, m0.t.copy())
    return mesh, Basis(mesh, element_cls())


def annulus_basis(n_r, r_in=0.4, r_out=1.0, element_cls=ElementTriArgyris):
    """Annulus r_in <= r <= r_out: structured polar grid (n_r radial cells,
    angular spacing ~ radial spacing), quads split into triangles; vertices
    exactly on both circles."""
    radii = np.linspace(r_in, r_out, n_r + 1)
    h = (r_out - r_in) / n_r
    n_th = int(np.ceil(2 * np.pi * 0.5 * (r_in + r_out) / h))
    th = np.linspace(0.0, 2 * np.pi, n_th, endpoint=False)
    pts = np.array([[r * np.cos(t), r * np.sin(t)]
                    for r in radii for t in th]).T
    tris = []
    for i in range(n_r):
        for j in range(n_th):
            j2 = (j + 1) % n_th
            v00 = i * n_th + j
            v01 = i * n_th + j2
            v10 = (i + 1) * n_th + j
            v11 = (i + 1) * n_th + j2
            tris.append([v00, v10, v11])
            tris.append([v00, v11, v01])
    mesh = MeshTri(pts.copy(), np.array(tris).T.copy())
    return mesh, Basis(mesh, element_cls())


def superellipse_basis(nrefine, a, b, p, element_cls=ElementTriArgyris):
    """Superellipse |x/a|^p + |y/b|^p = 1 (p = 2: ellipse; p -> inf:
    rectangle), centered, axes-aligned: the init_circle mesh mapped radially
    by rho(theta), so boundary vertices land exactly on the superellipse.
    Corner-region curvature grows with p at fixed smoothness -- the E5 knob."""
    m0 = MeshTri.init_circle(nrefine)
    x, y = m0.p
    th = np.arctan2(y, x)
    rho = (np.abs(np.cos(th) / a) ** p + np.abs(np.sin(th) / b) ** p) ** (-1.0 / p)
    pcoord = np.vstack([x * rho, y * rho])
    mesh = MeshTri(pcoord, m0.t.copy())
    return mesh, Basis(mesh, element_cls())


def sector_basis(nrings, theta=2.0, R=1.0, element_cls=ElementTriArgyris):
    """Disk sector of opening angle theta (radians), symmetric about the
    x-axis, apex at the origin. Polar-graded point rings (arc spacing ~ radial
    spacing) triangulated by Delaunay -- the domain is convex for theta < pi,
    so the triangulation fills it exactly; arc vertices lie exactly on the
    circle (inscribed-polygon arc, straight radial edges exact)."""
    from scipy.spatial import Delaunay
    assert theta < np.pi, "sector mesh assumes a convex sector (theta < pi)"
    pts = [(0.0, 0.0)]
    for k in range(1, nrings + 1):
        r = R * k / nrings
        n_th = max(3, int(np.ceil(theta * k)) + 1)
        ang = np.linspace(-0.5 * theta, 0.5 * theta, n_th)
        pts.extend(zip(r * np.cos(ang), r * np.sin(ang)))
    p = np.array(pts).T
    tri = Delaunay(p.T)
    t = tri.simplices.T
    a, b, c = p[:, t[0]], p[:, t[1]], p[:, t[2]]
    area2 = np.abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
    t = t[:, area2 > 1e-12]
    mesh = MeshTri(p.copy(), np.ascontiguousarray(t))
    return mesh, Basis(mesh, element_cls())


def triangle_ss_exact(L, n_modes):
    """Exact simply-supported equilateral-triangle PLATE spectrum. On a
    straight-edged polygon the SS (Navier) plate factorizes through the
    Dirichlet Helmholtz problem, whose equilateral-triangle eigenvalues are
    the Lame lattice lam_H = (16 pi^2 / (9 L^2)) (m^2 + m n + n^2), m,n >= 1
    (both (m,n) and (n,m) count; m != n is the doublet). Plate: Lambda = lam_H^2."""
    mmax = int(np.ceil(np.sqrt(n_modes) * 8))
    m = np.arange(1, mmax + 1)
    q = m[:, None] ** 2 + m[:, None] * m[None, :] + m[None, :] ** 2
    lam_h = (16.0 * np.pi ** 2 / (9.0 * L ** 2)) * q
    return np.sort(lam_h.ravel() ** 2)[:n_modes]
