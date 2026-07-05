#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""In-plane (plane-stress) elastodynamics of a free 2D rotor -- E15 machinery.

The flat Kirchhoff plate is Coriolis-blind for rotation about its normal
(Omega zhat x wdot zhat = 0): no first-order gyroscopic term. The correct
minimal continuum for the Coriolis / GOE->GUE program is IN-PLANE elasticity of
a free 2D body spinning about z, whose velocities are in-plane so the Coriolis
term is a genuine first-order gyroscopic bilinear form.

Conventions (E = 1, nu = 0.33, rho = 1 unless overridden):
  strain eps = sym(grad u),  plane-stress  sigma = lam* tr(eps) I + 2 mu eps
  with  lam* = E nu / (1 - nu^2),  mu = E / (2 (1 + nu)).
  K(u,v) = int sigma(u):eps(v) dA        (real symmetric, 3 rigid modes free)
  M(u,v) = int rho u.v dA                (real symmetric PD)
  G0(u,v) = 2 rho int (u_x v_y - u_y v_x) dA   (real ANTIsymmetric; Omega=1 form)
  rotating pencil:  (K - omega^2 M + i omega Omega G0) phi = 0
  -> Hermitian quadratic eigenvalue problem (i Omega G0 is Hermitian), real omega
     for the stable gyroscopic system.

Centrifugal terms (spin softening -Omega^2 M, geometric prestress) are real
symmetric (class-preserving) and O((Omega/omega)^2); v1 defers them (recorded
approximation, checked a posteriori via the Omega/omega ratio).

Free edges are the natural BCs (no constraints). Reuse kirchhoff.solve_lowest /
solve_modes (operator-agnostic) for the Omega = 0 certified spectrum.
"""
import numpy as np
from skfem import Basis, BilinearForm, ElementVector, ElementTriP2, MeshTri
from skfem.helpers import sym_grad, trace, ddot, dot


# ----------------------------- forms & assembly -----------------------------
def _lame_planestress(E, nu):
    return E * nu / (1.0 - nu * nu), E / (2.0 * (1.0 + nu))   # (lam*, mu)


def make_forms(nu, E=1.0, rho=1.0):
    lam, mu = _lame_planestress(E, nu)

    @BilinearForm
    def stiffness(u, v, w):
        eu, ev = sym_grad(u), sym_grad(v)
        return lam * trace(eu) * trace(ev) + 2.0 * mu * ddot(eu, ev)

    @BilinearForm
    def mass(u, v, w):
        return rho * dot(u, v)

    @BilinearForm
    def gyro0(u, v, w):
        # 2 rho (u_x v_y - u_y v_x); the Omega factor is applied outside.
        return 2.0 * rho * (u[0] * v[1] - u[1] * v[0])

    return stiffness, mass, gyro0


def assemble_elastic(mesh, basis, nu, E=1.0, rho=1.0):
    """Return (K, M, G0) as CSC. G0 is the Omega = 1 gyroscopic matrix."""
    stiffness, mass, gyro0 = make_forms(nu, E, rho)
    K = stiffness.assemble(basis).tocsc()
    M = mass.assemble(basis).tocsc()
    G0 = gyro0.assemble(basis).tocsc()
    return K, M, G0


# ------------------------------- domains ------------------------------------
def _vector_basis(mesh):
    return Basis(mesh, ElementVector(ElementTriP2()))


def disk_basis(nrefine, R=1.0):
    """Free unit disk (in-plane), init_circle mesh, vector P2."""
    m = MeshTri.init_circle(nrefine)
    if R != 1.0:
        m = MeshTri(R * m.p.copy(), m.t.copy())
    return m, _vector_basis(m)


def star_basis(nrefine, harmonics, R=1.0):
    """Star domain r(theta) = R (1 + sum_k a_k cos(k theta + phi_k)), the
    init_circle mesh mapped radially so boundary vertices land exactly on the
    star. `harmonics` = iterable of (k, a_k, phi_k). Amplitudes must keep
    r(theta) > 0 (checked). Symmetry by construction:
      even k only, no reflection  -> C2 (R_pi about z), chiral (D_prot)
      one term / mirror-aligned    -> a reflection axis (D_mirror)
      odd + even, generic phases    -> no C2, no mirror (D_chir)
    """
    m0 = MeshTri.init_circle(nrefine)
    x, y = m0.p
    th = np.arctan2(y, x)
    rmap = np.ones_like(th)
    for k, a, phi in harmonics:
        rmap = rmap + a * np.cos(k * th + phi)
    if np.min(rmap) <= 0.05:
        raise ValueError(f"star radial map non-positive/too thin: min={np.min(rmap):.3f}")
    p = np.vstack([x * rmap * R, y * rmap * R])
    m = MeshTri(p, m0.t.copy())
    return m, _vector_basis(m)


# ----------------------- point-mass mistuning (modal) -----------------------
def point_mass_modal(basis, X, points):
    """Modal mistuning matrix M_pts,ij = sum_k phi_i(x_k) . phi_j(x_k), where
    phi_i are the certified Omega=0 modes (columns of X) and points is (2, npts).
    Returns the real-symmetric N x N modal matrix (unit total point mass /
    normalized so that ||M_pts|| is O(mean spacing) after the c_delta scaling)."""
    P = basis.probes(np.asarray(points, float))     # (npts * ncomp, ndof) sparse
    Vals = P @ X                                     # (npts*2, Nmodes)
    npts = points.shape[1]
    Vx = Vals[0::2, :]                               # x-component at each point
    Vy = Vals[1::2, :]
    # sum over points of outer product of the vector value phi_i(x_k).phi_j(x_k)
    Mpts = Vx.T @ Vx + Vy.T @ Vy                     # (Nmodes, Nmodes)
    assert Vx.shape[0] == npts
    return 0.5 * (Mpts + Mpts.T)


# --------------------------- rotating pencil solve --------------------------
def solve_rotor(Lam, G0m, Omega, Mpts=None, delta=0.0, real_tol=1e-6):
    """Real positive frequencies of the modal rotating pencil
        (diag(Lam) - omega^2 M + i omega Omega G0m) phi = 0,
    with M = I + delta * Mpts (mistuning as a mass perturbation).

    Lam    : (N,) modal eigenvalues omega_n^2 (Omega=0 certified spectrum).
    G0m    : (N,N) real antisymmetric modal gyroscopic matrix (Omega=1).
    Returns dict: omega (sorted positive real freqs), max_imag (realness gate),
    pair_err (|+omega vs -omega| symmetry of the full spectrum)."""
    N = len(Lam)
    M = np.eye(N)
    if Mpts is not None and delta != 0.0:
        M = M + delta * Mpts
    K = np.diag(Lam)
    H = Omega * G0m                                  # real antisymmetric
    # QEP (omega^2 M - i omega H - K) phi = 0 ; companion A z = omega B z
    I = np.eye(N)
    Z = np.zeros((N, N))
    A = np.block([[Z, I], [K, 1j * H]])
    B = np.block([[I, Z], [Z, M]])
    from scipy.linalg import eig
    ev = eig(A, B, right=False)
    ev = ev[np.isfinite(ev)]
    max_imag = float(np.max(np.abs(ev.imag)) / (np.max(np.abs(ev.real)) + 1e-300))
    w = ev.real
    pos = np.sort(w[w > 0])
    neg = np.sort(-w[w < 0])
    n = min(len(pos), len(neg))
    pair_err = float(np.max(np.abs(pos[:n] - neg[:n])) /
                     (pos[n - 1] if n else 1.0)) if n else np.nan
    return dict(omega=pos, max_imag=max_imag, pair_err=pair_err,
                n_pos=len(pos), n_neg=len(neg))


def modal_reduce(K, M, G0, X):
    """Project (K, M, G0) onto the certified modes X (M-orthonormalized).
    Returns (Lam diag values, G0_modal). Assumes X already M-orthonormal;
    re-orthonormalizes defensively."""
    MX = M @ X
    B = X.T @ MX
    # symmetric orthonormalization w.r.t. M
    ev, U = np.linalg.eigh(0.5 * (B + B.T))
    T = U / np.sqrt(np.clip(ev, 1e-14, None))
    Xn = X @ T
    Lam = np.diag(Xn.T @ (K @ Xn))
    G0m = Xn.T @ (G0 @ Xn)
    G0m = 0.5 * (G0m - G0m.T)                        # enforce antisymmetry
    return np.asarray(Lam), G0m, Xn
