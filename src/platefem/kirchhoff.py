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


def solve_modes(K, M, k, sigma=-10.0, block_extra=0.25, sweeps_max=15, seed=7,
                match_frac=0.5, resid_sanity=1e-2, Y0=None):
    """Certified lowest-k eigenPAIRS of K x = Lambda M x.

    Eigenvalues of record come from the no-vector ARPACK path (solve_lowest),
    which is reliable for these matrices. Eigenvectors are computed WITHOUT
    ARPACK: scipy 1.18.0's return_eigenvectors=True shift-invert path is
    broken for the ill-conditioned Argyris matrices (eigenvalues shifted by
    exactly -sigma when it converges, misconverged Ritz pairs or ARPACK
    error -8 otherwise; see experiments/e01 and the E2 work log). Instead:
    block inverse iteration on the factorized (K - sigma M) with a dense
    Rayleigh-Ritz projection each sweep -- deterministic, and near-degenerate
    clusters are resolved exactly within the projected subspace.

    Certification per sweep: Rayleigh quotients of the Ritz vectors must match
    the certified eigenvalues 1:1 within match_frac x local mean spacing, and
    relative residuals must pass resid_sanity for every elastic mode. These
    vectors are used for LABELING (parity classification), never for spacing
    statistics; residual mixing inside near-degenerate clusters is handled by
    cluster-resolved classification downstream.

    Y0: optional warm-start block (from a previous nearby solve, e.g. the
    preceding kappa in a sweep). Returns (lam, X, info, Y) with Y the full
    converged block for warm-starting the next call.
    """
    from scipy.sparse.linalg import splu
    import scipy.linalg as sla
    from .stats import local_spacing

    lam = solve_lowest(K, M, k, sigma)
    d = np.maximum(local_spacing(lam), 1e-8 * abs(lam[-1]))
    elastic = np.abs(lam) > 1e-6 * abs(lam[-1])
    n = K.shape[0]
    m = min(n, k + max(16, int(block_extra * k)))
    lu = splu((K - sigma * M).tocsc())
    rng = np.random.default_rng(seed)
    if Y0 is not None and Y0.shape[0] == n:
        Y = np.hstack([Y0, rng.standard_normal((n, max(0, m - Y0.shape[1])))])[:, :m]
    else:
        Y = rng.standard_normal((n, m))
    info = None
    for sweep in range(1, sweeps_max + 1):
        Y = lu.solve(M @ Y)                      # inverse-iteration step
        Y, _ = np.linalg.qr(Y)
        Ar = Y.T @ (K @ Y)
        Br = Y.T @ (M @ Y)
        Ar, Br = 0.5 * (Ar + Ar.T), 0.5 * (Br + Br.T)
        try:
            th, C = sla.eigh(Ar, Br)
        except np.linalg.LinAlgError:            # M is numerically semidefinite
            jit = 1e-12 * np.trace(Br) / m
            th, C = sla.eigh(Ar, Br + jit * np.eye(m))
        Y = Y @ C                                # Ritz vectors, ascending theta
        X = Y[:, :k]
        KX = K @ X
        MX = M @ X
        rq = np.einsum("ij,ij->j", X, KX) / np.einsum("ij,ij->j", X, MX)
        resid = np.linalg.norm(KX - rq * MX, axis=0) / (
            np.linalg.norm(KX, axis=0) + np.abs(rq) * np.linalg.norm(MX, axis=0))
        mismatch = np.abs(rq - lam) > match_frac * d
        max_res = float(np.max(resid[elastic])) if elastic.any() else 0.0
        info = dict(sweeps=sweep, n_mismatch=int(np.sum(mismatch)),
                    max_resid=max_res)
        if not mismatch.any() and max_res < resid_sanity:
            return lam, X, info, Y
    raise RuntimeError(f"solve_modes: not certified after {sweeps_max} sweeps; "
                       f"last diagnostics {info}")


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
