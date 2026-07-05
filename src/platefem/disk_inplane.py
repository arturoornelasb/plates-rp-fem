#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Semi-analytic in-plane (plane-stress) free circular disk -- E15 exact anchor (G1).

Helmholtz potentials u = grad(phi) + curl(psi zhat), with
  phi = J_n(alpha r) cos n(theta),  psi = J_n(beta r) sin n(theta),
alpha = omega/c_p, beta = omega/c_s. Plane stress (E = rho = 1):
  c_p^2 = 1/(1-nu^2),  c_s^2 = 1/(2(1+nu))  =>  alpha^2 = omega^2 (1-nu^2),
  beta^2 = 2 omega^2 (1+nu);  mu = 1/(2(1+nu)), lam* = nu/(1-nu^2).

Traction-free edge sigma_rr = sigma_rt = 0 at r = R gives a 2x2 determinant in
the potential amplitudes (A, B). Derived from sigma_rr = 2 mu U' - lam* alpha^2 f
and sigma_rt/mu = -2n/r (f' - f/r) + 2/r g' + (beta^2 - 2n^2/r^2) g, using the
Bessel equation to remove second derivatives (verified analytically on the n=0
radial mode; the n=0 torsional branch reduces to J_2(beta R) = 0).

Eigenvalue Lambda = omega^2. For n >= 1 each root is a doublet (cos/sin);
n = 0 yields a radial (dilatational) and a torsional (shear) branch. Rigid modes
(2 translations at n=1, rotation at n=0) sit at omega = 0 and are excluded by the
scan floor. Cross-validates the elastic2d FEM disk (gate G1).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import jv


def _Jp(n, x):
    return 0.5 * (jv(n - 1, x) - jv(n + 1, x))


def _det(n, omega, nu, R=1.0):
    """2x2 traction-free determinant for azimuthal order n at frequency omega."""
    mu = 1.0 / (2.0 * (1.0 + nu))
    a = omega * np.sqrt(1.0 - nu * nu) * R          # alpha R
    b = omega * np.sqrt(2.0 * (1.0 + nu)) * R        # beta R
    al, be = a / R, b / R
    Jna, Jnb = jv(n, a), jv(n, b)
    Jpa, Jpb = _Jp(n, a), _Jp(n, b)
    w2 = omega * omega
    # sigma_rr row (x A, x B)
    M11 = -2.0 * mu * al * Jpa - (w2 - 2.0 * mu * n * n / R ** 2) * Jna
    M12 = 2.0 * mu * n * (be * Jpb - Jnb / R) / R
    # sigma_rt row / mu
    M21 = 2.0 * n * (-al * Jpa + Jna / R) / R
    M22 = 2.0 * be * Jpb / R + (be ** 2 - 2.0 * n ** 2 / R ** 2) * Jnb
    return M11 * M22 - M12 * M21


def inplane_free_disk(nu, omega_max, n_max, R=1.0, omega_min=0.3, spu=40):
    """Roots Lambda = omega^2 of the in-plane free-disk determinant.
    Returns (levels sorted, mult) where a class-list convention is used:
    n = 0 branches once each; n >= 1 as doublets (mult 2). `levels` is the sorted
    Lambda array WITH multiplicities expanded (as the FEM sees them)."""
    lam, mult = [], []
    for n in range(n_max + 1):
        xs = np.linspace(omega_min, omega_max,
                         max(200, int(spu * (omega_max - omega_min))))
        vals = np.array([_det(n, w, nu, R) for w in xs])
        sign = np.sign(vals)
        idx = np.where(sign[:-1] * sign[1:] < 0)[0]
        for i in idx:
            try:
                w = brentq(lambda w: _det(n, w, nu, R), xs[i], xs[i + 1],
                           xtol=1e-11, rtol=1e-13)
            except ValueError:
                continue
            lam.append(w * w)
            mult.append(1 if n == 0 else 2)
    order = np.argsort(lam)
    lam = np.array(lam)[order]
    mult = np.array(mult)[order]
    return lam, mult


def expand(lam, mult):
    return np.repeat(lam, mult)


def validate(nu=0.33):
    """n=0 torsional branch must satisfy J_2(beta R) = 0 (analytic identity)."""
    from scipy.special import jn_zeros
    b1 = jn_zeros(2, 1)[0]                            # first zero of J_2 ~ 5.1356
    w = b1 / np.sqrt(2.0 * (1.0 + nu))                # beta = b/R, R=1
    d = _det(0, w, nu)
    return b1, w, d
