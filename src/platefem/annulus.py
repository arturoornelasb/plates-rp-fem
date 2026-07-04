#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Semi-analytic free-free Kirchhoff annulus (second Poisson control, E10).

Radial solution W = A J_m(br) + B Y_m(br) + C I_m(br) + D K_m(br); free-edge
conditions (M_r = 0, V_r = 0; same operators as platefem.disk, validated
there against the literature fundamental) at BOTH radii give a 4x4
determinant per azimuthal m. Exponential scaling: the I column is scaled by
exp(-x_out) and the K column by exp(+x_in) (one constant per column, so
determinant roots are unchanged); the residual factors exp(-(x_out - x_in))
are always <= 1.

The annulus preserves polar separability with free edges -- like the full
disk, it must give Poisson when the independent m sequences of one
reflection class are pooled.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import ive, jv, kve, yv


def _rowvals(m, x, nu, kind):
    """(M_r, V_r) condition values for one radial function at argument x.
    kind in {'J','Y','I','K'}; I/K are the scipy exponentially scaled forms
    (ive, kve) -- the caller accounts for the column scaling."""
    mm = m * m
    if kind in ("J", "Y"):
        f = jv if kind == "J" else yv
        F = f(m, x)
        Fp = 0.5 * (f(m - 1, x) - f(m + 1, x))
        Fpp = -Fp / x + (mm / (x * x) - 1.0) * F
        lap_sign = -1.0
    elif kind == "I":
        F = ive(m, x)
        Fp = 0.5 * (ive(m - 1, x) + ive(m + 1, x))
        Fpp = -Fp / x + (mm / (x * x) + 1.0) * F
        lap_sign = +1.0
    else:
        F = kve(m, x)
        Fp = -0.5 * (kve(m - 1, x) + kve(m + 1, x))
        Fpp = -Fp / x + (mm / (x * x) + 1.0) * F
        lap_sign = +1.0
    Mrow = Fpp + nu * (Fp / x - mm * F / (x * x))
    Vrow = lap_sign * Fp - (1.0 - nu) * (mm / (x * x)) * (Fp - F / x)
    return Mrow, Vrow


def _det(m, beta, nu, r_in, r_out):
    xo, xi = beta * r_out, beta * r_in
    scale = np.exp(-(xo - xi))
    rows = []
    for x, sc_I, sc_K in [(xo, 1.0, scale), (xi, scale, 1.0)]:
        MJ, VJ = _rowvals(m, x, nu, "J")
        MY, VY = _rowvals(m, x, nu, "Y")
        MI, VI = _rowvals(m, x, nu, "I")
        MK, VK = _rowvals(m, x, nu, "K")
        rows.append([MJ, MY, MI * sc_I, MK * sc_K])
        rows.append([VJ, VY, VI * sc_I, VK * sc_K])
    return np.linalg.det(np.array(rows))


def free_annulus_roots(nu, r_in, r_out, x_max, m_max, x_min=0.3,
                       samples_per_unit=40):
    """Roots beta*r_out = x of the free-free annulus determinant.
    Returns dict m -> ascending roots (in beta, physical)."""
    out = {}
    for m in range(m_max + 1):
        lo = max(x_min, 0.5 * m)         # turning-point guard (outer radius)
        if lo >= x_max:
            continue
        xs = np.linspace(lo, x_max, max(64, int(samples_per_unit * (x_max - lo))))
        vals = np.array([_det(m, x / r_out, nu, r_in, r_out) for x in xs])
        good = np.isfinite(vals)
        xs, vals = xs[good], vals[good]
        sign = np.sign(vals)
        idx = np.where(sign[:-1] * sign[1:] < 0)[0]
        roots = []
        for i in idx:
            try:
                roots.append(brentq(
                    lambda x: _det(m, x / r_out, nu, r_in, r_out),
                    xs[i], xs[i + 1], xtol=1e-12, rtol=1e-14))
            except ValueError:
                pass
        if roots:
            out[m] = np.array(roots) / r_out          # beta values
    return out


def class_levels(roots):
    """One-reflection-class Lambda list (m taken once), sorted."""
    lam = []
    for m, bs in roots.items():
        lam.extend((bs ** 4).tolist())
    return np.sort(np.array(lam))


def full_spectrum(roots):
    lam, mult = [], []
    for m, bs in roots.items():
        for b_ in bs:
            lam.append(b_ ** 4)
            mult.append(1 if m == 0 else 2)
    order = np.argsort(lam)
    return np.array(lam)[order], np.array(mult)[order]
