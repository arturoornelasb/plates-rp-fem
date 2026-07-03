#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Semi-analytic free-edge Kirchhoff circular plate (the E3 Poisson control).

For each azimuthal number m the radial factor is W = A J_m(br) + B I_m(br)
and the free-edge conditions at r = R (with x = bR),

  M_r = 0:  W'' + nu (W'/r - m^2 W / r^2) = 0
  V_r = 0:  (lap_m W)' - (1-nu) (m^2/r^2)(W' - W/r) = 0
  (V_r = Q_r + (1/r) dM_rt/dtheta; the twisting-moment term carries the
   minus sign in this normalization),

give a 2x2 determinant whose roots x_mn yield plate eigenvalues
Lambda = (x/R)^4 (D = rho h = 1). lap_m J = -b^2 J and lap_m I = +b^2 I.
The I column is evaluated with scipy's exponentially scaled ive (a common
positive factor exp(x) per column does not move determinant roots).

This is the same construction the paper's executed disk control used
(validated there against tabulated frequency parameters to < 0.4%); the
classic fundamental for nu = 0.33 is x^2 = 5.25 at m = 2 (Leissa 1969,
Itao & Crandall 1979) and is asserted below.

Statistics convention (paper): one reflection class = all m = 0 levels plus
each m >= 1 level ONCE (the cos class); each fixed-m sequence is a
Berry-Tabor picket fence; pooling the class gives Poisson.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import ive, jv


def _cols(m, x, nu):
    """Free-edge condition rows evaluated for the J and (scaled) I columns."""
    J = jv(m, x)
    Jp = 0.5 * (jv(m - 1, x) - jv(m + 1, x))
    Jpp = -Jp / x + (m * m / (x * x) - 1.0) * J
    Ie = ive(m, x)
    Iep = 0.5 * (ive(m - 1, x) + ive(m + 1, x))
    Iepp = -Iep / x + (m * m / (x * x) + 1.0) * Ie
    # M_r row (per column)
    mJ = Jpp + nu * (Jp / x - m * m * J / (x * x))
    mI = Iepp + nu * (Iep / x - m * m * Ie / (x * x))
    # V_r row
    vJ = -Jp - (1.0 - nu) * (m * m / (x * x)) * (Jp - J / x)
    vI = Iep - (1.0 - nu) * (m * m / (x * x)) * (Iep - Ie / x)
    return mJ * vI - mI * vJ


def free_disk_roots(nu, x_max, m_max, x_min=0.35, samples_per_unit=40):
    """Roots x_mn of the free-disk determinant for m = 0..m_max, x <= x_max.
    Returns dict m -> array of roots (ascending). Rigid modes (piston m=0,
    tilt m=1 at x=0) are excluded by the scan floor."""
    out = {}
    for m in range(m_max + 1):
        lo = max(x_min, 0.6 * m)                 # below the Bessel turning point
        if lo >= x_max:                          # (x ~ m) the determinant is
            continue                             # numerical noise; roots sit above
        xs = np.linspace(lo, x_max, max(64, int(samples_per_unit * (x_max - lo))))
        vals = _cols(m, xs, nu)
        roots = []
        sign = np.sign(vals)
        idx = np.where(sign[:-1] * sign[1:] < 0)[0]
        for i in idx:
            try:
                r = brentq(lambda x: _cols(m, x, nu), xs[i], xs[i + 1],
                           xtol=1e-12, rtol=1e-14)
                roots.append(r)
            except ValueError:
                pass
        if roots:
            out[m] = np.array(roots)
    return out


def class_levels(roots, R=1.0):
    """One-reflection-class level list: (Lambda, m, n) for every root, m >= 0
    taken once. Returns levels sorted by Lambda."""
    rows = []
    for m, rr in roots.items():
        for n, x in enumerate(rr):
            rows.append(((x / R) ** 4, m, n))
    rows.sort()
    return rows


def full_spectrum(roots, R=1.0):
    """Full free-disk spectrum with multiplicities: m = 0 once, m >= 1 twice.
    Returns (sorted Lambda array, multiplicity array aligned with it)."""
    lam, mult = [], []
    for m, rr in roots.items():
        for x in rr:
            lam.append((x / R) ** 4)
            mult.append(1 if m == 0 else 2)
    order = np.argsort(lam)
    return np.array(lam)[order], np.array(mult)[order]


def expand_multiplicity(lam, mult):
    """Spectrum as the FEM sees it (doublets as two consecutive levels)."""
    return np.repeat(lam, mult)


def dedupe_by_reference(lam_fem, lam_ref, mult_ref):
    """Collapse FEM doublets guided by the reference multiplicity pattern.
    Walks both sorted lists; a mult-2 reference level consumes two FEM values
    (their mean). Returns (class-sequence array, max relative pair splitting)."""
    out, splits = [], []
    i = 0
    for L, mu in zip(lam_ref, mult_ref):
        if i + mu > len(lam_fem):
            break
        vals = lam_fem[i:i + mu]
        out.append(float(np.mean(vals)))
        if mu == 2:
            splits.append(float((vals[1] - vals[0]) / vals[0]))
        i += mu
    return np.array(out), (max(splits) if splits else 0.0)


def validate(nu=0.33):
    """Fundamental of the free disk: x^2 ~ 5.25 at m = 2 for nu = 0.33."""
    roots = free_disk_roots(nu, x_max=6.0, m_max=3)
    x2 = roots[2][0] ** 2
    assert 5.1 < x2 < 5.4, f"free-disk fundamental x^2 = {x2:.3f}, expected ~5.25"
    return x2
