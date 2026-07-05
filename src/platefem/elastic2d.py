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


def star_polar_basis(nrings, n_th, harmonics, R=1.0):
    """Structured polar mesh of a star r(theta) = R (1 + sum a_k cos(k th+phi_k)):
    a center node + nrings rings of n_th angular nodes. n_th MUST be even. The
    node set respects, EXACTLY as index permutations,
      sigma_v (theta -> -theta)      when r is even (cos, phi in {0, pi}),
      R_pi    (theta -> theta + pi)  when r has even k only,
    so modes carry exact symmetry parity -- required for clean per-class <r>.
    Returns (mesh, basis, meta) with meta['theta'], meta['n_th'], meta['nrings']."""
    assert n_th % 2 == 0, "n_th must be even for exact R_pi / sigma_v node maps"
    th = 2.0 * np.pi * np.arange(n_th) / n_th
    rmap = np.ones(n_th)
    for k, a, phi in harmonics:
        rmap = rmap + a * np.cos(k * th + phi)
    if rmap.min() <= 0.05:
        raise ValueError(f"star radial map too thin: min={rmap.min():.3f}")
    # ring nodes: index(i,j) = 1 + (i-1)*n_th + j, i = 1..nrings ; center = 0.
    pts = [(0.0, 0.0)]
    for i in range(1, nrings + 1):
        frac = i / nrings
        for j in range(n_th):
            r = R * frac * rmap[j]
            pts.append((r * np.cos(th[j]), r * np.sin(th[j])))

    def rn(i, j):
        return 1 + (i - 1) * n_th + (j % n_th)

    tris = []
    for j in range(n_th):                         # center fan (already symmetric)
        tris.append([0, rn(1, j), rn(1, j + 1)])
    # ring quads split from their CENTROID into 4 tris: symmetric under BOTH
    # sigma_v (theta->-theta) and R_pi (theta->theta+pi), unlike a fixed diagonal.
    cbase = len(pts)
    for i in range(1, nrings):
        for j in range(n_th):
            A, B = rn(i, j), rn(i, j + 1)
            C, D = rn(i + 1, j + 1), rn(i + 1, j)
            c = cbase + (i - 1) * n_th + j
            pts.append((0.25 * (pts[A][0] + pts[B][0] + pts[C][0] + pts[D][0]),
                        0.25 * (pts[A][1] + pts[B][1] + pts[C][1] + pts[D][1])))
            tris += [[c, A, B], [c, B, C], [c, C, D], [c, D, A]]
    p = np.array(pts).T
    mesh = MeshTri(p.copy(), np.array(tris).T.copy())
    return mesh, _vector_basis(mesh), dict(theta=th, n_th=n_th, nrings=nrings)


# ----------------------- symmetry-class separation --------------------------
def _probe_vals(basis, pts):
    """(Vx, Vy) sparse->dense probe operators returning per-point components.
    basis.probes on a vector element returns interleaved (x,y) rows."""
    pts = np.asarray(pts, float)
    npts = pts.shape[1]
    P = basis.probes(pts).tocsr()                # rows blocked: [all x; all y]
    return P[:npts, :], P[npts:, :]              # (Vx, Vy), each (npts, ndof)


def parity_classes(basis, X, kind, harmonics, R=1.0, npts=600, inset=0.92, seed=3):
    """Classify certified modes X (columns, Omega=0) by an EXACT domain symmetry.
    kind='mirror_x' (sigma_v about x): image (x,-y), vector (u_x,-u_y).
    kind='rpi'      (R_pi about z):    image (-x,-y), vector (-u_x,-u_y).
    Sample points strictly inside the star r(theta)=R(1+sum a_k cos(k th+phi))
    (radius = u * inset * r(theta)); the image of an interior star point is also
    interior for the relevant symmetry, so probes never leave the mesh.
    Returns (labels in {+1,-1,0}, c); 0 = ambiguous (|c| < 0.8)."""
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0, 2 * np.pi, npts)
    rth = np.ones(npts)
    for k, a, phi in harmonics:
        rth = rth + a * np.cos(k * ang + phi)
    rad = R * inset * rth * np.sqrt(rng.uniform(0, 1, npts))
    p = np.vstack([rad * np.cos(ang), rad * np.sin(ang)])
    if kind == "mirror_x":
        pim = np.vstack([p[0], -p[1]]); sx, sy = 1.0, -1.0
    elif kind == "rpi":
        pim = np.vstack([-p[0], -p[1]]); sx, sy = -1.0, -1.0
    else:
        raise ValueError(kind)
    Vx, Vy = _probe_vals(basis, p)
    Wx, Wy = _probe_vals(basis, pim)
    A = np.vstack([Vx @ X, Vy @ X])              # field at p        (2*npts, Nm)
    B = np.vstack([sx * (Wx @ X), sy * (Wy @ X)])  # S-transformed field
    num = np.einsum("ij,ij->j", A, B)
    c = num / np.sqrt(np.einsum("ij,ij->j", A, A) * np.einsum("ij,ij->j", B, B))
    lab = np.where(np.abs(c) < 0.8, 0, np.sign(c)).astype(int)
    return lab, c


def build_symop(basis, kind, tol=1e-7):
    """Exact discrete symmetry operator S (signed dof permutation) for a vector
    P2 basis on a mesh that respects the symmetry:
      kind='mirror_x' (sigma_v about x): loc (x,y)->(x,-y), (u_x,u_y)->(u_x,-u_y)
      kind='rpi'      (R_pi about z):    loc (x,y)->(-x,-y), (u_x,u_y)->(-u_x,-u_y)
    Returns scipy.sparse S (ndof x ndof), orthogonal involution (S@S = I) IFF the
    mesh is symmetric to `tol`. Raises if an image dof location is missing."""
    from scipy.sparse import coo_matrix
    from scipy.spatial import cKDTree
    ix, iy = basis.split_indices()               # x-comp dofs, y-comp dofs
    L = basis.doflocs[:, ix]                      # scalar-dof locations (2, m)
    m = L.shape[1]
    if kind == "mirror_x":
        tgt = np.vstack([L[0], -L[1]]); sgx, sgy = 1.0, -1.0
    elif kind == "rpi":
        tgt = np.vstack([-L[0], -L[1]]); sgx, sgy = -1.0, -1.0
    else:
        raise ValueError(kind)
    tree = cKDTree(L.T)
    dist, perm = tree.query(tgt.T, k=1)
    if dist.max() > tol:
        raise RuntimeError(f"symop {kind}: image dof missing (max match dist "
                           f"{dist.max():.1e} > tol {tol}); mesh not symmetric")
    rows = np.concatenate([ix, iy])
    cols = np.concatenate([ix[perm], iy[perm]])
    data = np.concatenate([np.full(m, sgx), np.full(m, sgy)])
    return coo_matrix((data, (rows, cols)), shape=(basis.N, basis.N)).tocsr()


def sym_residuals(S, K, M, G0):
    """T5-style check: how exactly does S act on the operators.
      S involution : ||S@S - I||
      commutes K   : ||S K S - K|| / ||K||     (want ~0: sigma_v, R_pi keep K,M)
      commutes M   : ||S M S - M|| / ||M||
      G0 sign      : ||S G0 S - G0|| and ||S G0 S + G0|| (Coriolis: sigma_v
                     ANTI-commutes -> +G0 residual ~0; R_pi COMMUTES -> -G0 ~0)."""
    import scipy.sparse as sp
    n = K.shape[0]
    I = sp.identity(n, format="csr")
    def rel(A, B):
        return float(sp.linalg.norm(A - B) / (sp.linalg.norm(B) + 1e-300))
    SKS, SMS, SGS = S @ K @ S, S @ M @ S, S @ G0 @ S
    return dict(involution=float(sp.linalg.norm(S @ S - I)),
                commute_K=rel(SKS, K), commute_M=rel(SMS, M),
                G0_anticommute=rel(SGS, -G0), G0_commute=rel(SGS, G0))


# ----------------------- point-mass mistuning (modal) -----------------------
def point_mass_modal(basis, X, points):
    """Modal mistuning matrix M_pts,ij = sum_k phi_i(x_k) . phi_j(x_k), where
    phi_i are the certified Omega=0 modes (columns of X) and points is (2, npts).
    Returns the real-symmetric N x N modal matrix (unit total point mass /
    normalized so that ||M_pts|| is O(mean spacing) after the c_delta scaling)."""
    points = np.asarray(points, float)
    npts = points.shape[1]
    P = basis.probes(points).tocsr()                 # rows blocked: [all x; all y]
    Vals = P @ X                                      # (npts*2, Nmodes)
    Vx = Vals[:npts, :]                               # x-component at each point
    Vy = Vals[npts:, :]
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


def parity_adapt_reduce(K, M, G0, X, S):
    """Reduce (K,M,G0) onto the certified span in a sigma-PARITY-ADAPTED,
    M-orthonormal basis (using the exact symmetry operator S). Because S is
    M-orthogonal (S^T M S = M) the even/odd subspaces are M-orthogonal; splitting
    the (possibly cluster-mixed) certified modes X into clean parity blocks makes
    the reduced G0m EXACTLY block-off-diagonal when S anticommutes with G0
    (sigma_v case) -- restoring the sigma_v*T orthogonal-class structure that
    cluster-mixed vectors destroy. Returns (Lam, G0m, labels, Xn)."""
    def orth(Y):
        B = 0.5 * (Y.T @ (M @ Y) + (Y.T @ (M @ Y)).T)
        ev, U = np.linalg.eigh(B)
        keep = ev > 1e-10 * ev[-1]
        return Y @ (U[:, keep] / np.sqrt(ev[keep]))
    SX = S @ X
    # diagonalize K WITHIN each parity block (K commutes with S) -> true
    # eigenvalues + clean parity eigenvectors (diag(X'KX) alone would discard the
    # off-diagonal K coupling of cluster-mixed vectors and fake Poisson).
    cols, Lam, lab = [], [], []
    for Yb, s in [(orth(0.5 * (X + SX)), 1), (orth(0.5 * (X - SX)), -1)]:
        eb, Wb = np.linalg.eigh(0.5 * (Yb.T @ (K @ Yb) + (Yb.T @ (K @ Yb)).T))
        cols.append(Yb @ Wb)
        Lam.append(eb)
        lab.append(np.full(len(eb), s))
    Xn = np.hstack(cols)
    Lam = np.concatenate(Lam)
    lab = np.concatenate(lab)
    G0m = Xn.T @ (G0 @ Xn)
    G0m = 0.5 * (G0m - G0m.T)
    return np.asarray(Lam), G0m, lab, Xn


def modal_reduce(K, M, G0, X):
    """Project (K, M, G0) onto the certified modes X (M-orthonormalized).
    Returns (Lam diag values, G0_modal). Assumes X already M-orthonormal;
    re-orthonormalizes defensively."""
    MX = M @ X
    B = X.T @ MX
    # symmetric orthonormalization w.r.t. M
    ev, U = np.linalg.eigh(0.5 * (B + B.T))
    keep = ev > 1e-10 * ev[-1]
    T = U[:, keep] / np.sqrt(ev[keep])
    Xn = X @ T
    # diagonalize the reduced K (NOT diag(): that would drop the off-diagonal
    # coupling of cluster-mixed certified vectors and fake a Poisson spectrum).
    Km = 0.5 * (Xn.T @ (K @ Xn) + (Xn.T @ (K @ Xn)).T)
    Lam, W = np.linalg.eigh(Km)
    Xn = Xn @ W
    G0m = Xn.T @ (G0 @ Xn)
    G0m = 0.5 * (G0m - G0m.T)                        # enforce antisymmetry
    return np.asarray(Lam), G0m, Xn
