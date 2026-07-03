#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ritz_reference.py -- spectral Rayleigh-Ritz reference for the FFFF Kirchhoff rectangle
======================================================================================
Assembly adapted verbatim from
zenodo_v4/code_pretests/blindaje_v4_suite/T1_basis_audit/t1_basis_audit.py
(same operator, same conventions). This module only exposes what step 1 needs:
convergence-filtered *physical* eigenvalues per parity sector.

Conventions (identical to the repo pretests):
  - rectangle with sides a (x) and b (y), reference square xi,eta in [-1,1]
  - Kirchhoff flexural eigenproblem  D lap^2 w = rho h omega^2 w  with D = rho h = 1
  - orthonormal Legendre Ritz basis; mass matrix = (a*b/4) * I, so the physical
    eigenvalue is  Lambda = omega^2 = eig(K) / (a*b/4)
  - parity sectors ee/eo/oe/oo; rigid modes: piston (ee), y-tilt (eo), x-tilt (oe)
"""
import numpy as np
import numpy.polynomial.legendre as npleg
from scipy.linalg import eigh

SECTORS = ["ee", "eo", "oe", "oo"]
N_RIGID = {"ee": 1, "eo": 1, "oe": 1, "oo": 0}


def legendre_tables(nleg, ngrid):
    xg, wg = npleg.leggauss(ngrid)
    V0 = np.zeros((nleg, ngrid)); V1 = np.zeros_like(V0); V2 = np.zeros_like(V0)
    for i in range(nleg):
        c = np.zeros(i + 1); c[i] = 1.0
        nrm = np.sqrt((2 * i + 1) / 2.0)
        V0[i] = npleg.legval(xg, c) * nrm
        V1[i] = npleg.legval(xg, npleg.legder(c, 1)) * nrm
        V2[i] = npleg.legval(xg, npleg.legder(c, 2)) * nrm
    return xg, wg, V0, V1, V2


def mats_1d(wg, V0, V1, V2):
    A0 = (V0 * wg) @ V0.T
    A1 = (V1 * wg) @ V1.T
    A2 = (V2 * wg) @ V2.T
    B = (V2 * wg) @ V0.T
    return A0, A1, A2, B


def assemble_K_sector(idx_x, idx_y, A1, A2, B, a, b, nu):
    Jx, Jy = 2.0 / a, 2.0 / b
    pre = (a * b) / 4.0
    A1x, A2x, Bx = A1[np.ix_(idx_x, idx_x)], A2[np.ix_(idx_x, idx_x)], B[np.ix_(idx_x, idx_x)]
    A1y, A2y, By = A1[np.ix_(idx_y, idx_y)], A2[np.ix_(idx_y, idx_y)], B[np.ix_(idx_y, idx_y)]
    Ix, Iy = np.eye(len(idx_x)), np.eye(len(idx_y))
    K = (Jx ** 4) * np.kron(A2x, Iy) + (Jy ** 4) * np.kron(Ix, A2y)
    K += nu * (Jx ** 2) * (Jy ** 2) * (np.kron(Bx, By.T) + np.kron(Bx.T, By))
    K += 2.0 * (1.0 - nu) * (Jx ** 2) * (Jy ** 2) * np.kron(A1x, A1y)
    return pre * K


def sector_indices(sec, nleg):
    ix = np.arange(0 if sec[0] == "e" else 1, nleg, 2)
    iy = np.arange(0 if sec[1] == "e" else 1, nleg, 2)
    return ix, iy


def sector_spectrum(sec, a, b, nu, nleg, tables_cache):
    """Physical eigenvalues Lambda = omega^2 of one parity sector (rigid included)."""
    if nleg not in tables_cache:
        xg, wg, V0, V1, V2 = legendre_tables(nleg, max(400, 2 * nleg + 50))
        tables_cache[nleg] = mats_1d(wg, V0, V1, V2)
    A0, A1, A2, B = tables_cache[nleg]
    ix, iy = sector_indices(sec, nleg)
    K = assemble_K_sector(ix, iy, A1, A2, B, a, b, nu)
    ev = eigh(0.5 * (K + K.T), eigvals_only=True) / ((a * b) / 4.0)
    return ev


def converged_spectrum(a, b, nu, nleg=64, nleg_conv=72, tol=1e-4):
    """Convergence-filtered elastic spectra per sector (repo protocol: NLEG vs
    NLEG_CONV, keep modes with relative shift < tol). Note: per-mode shifts
    plateau around ~1e-4 -- diag_ritz_nleg.py shows the eigenvalues go
    NON-monotone in NLEG beyond ~96 (round-off corrupts the high-order
    Legendre assembly), so ~1e-4 is this reference's floor, NLEG 64/72 its
    sweet spot, and a tighter tol empties the sectors. Returns dict:
      sector -> dict(elastic=Lambda array, rel=per-mode NLEG-shift array,
                     n_conv=int, n_rigid=int)
    """
    cache = {}
    out = {}
    for sec in SECTORS:
        evA = sector_spectrum(sec, a, b, nu, nleg, cache)
        evB = sector_spectrum(sec, a, b, nu, nleg_conv, cache)
        ref = abs(evA[min(5, len(evA) - 1)])
        nzA = int(np.sum(evA < 1e-6 * ref))
        nzB = int(np.sum(evB < 1e-6 * ref))
        assert nzA == N_RIGID[sec] and nzB == N_RIGID[sec], \
            f"rigid count {sec}: {nzA},{nzB} != {N_RIGID[sec]}"
        wA, wB = evA[nzA:], evB[nzB:]
        nmin = min(len(wA), len(wB))
        rel = np.abs(wA[:nmin] - wB[:nmin]) / wB[:nmin]
        bad = np.where(rel > tol)[0]
        n_conv = int(bad[0]) if len(bad) else nmin
        if n_conv == 0:
            raise RuntimeError(f"sector {sec}: no modes converged at tol {tol:g} "
                               f"(first-mode shift {rel[0]:.2e}); relax tol")
        out[sec] = dict(elastic=wA[:n_conv], rel=rel[:n_conv],
                        n_conv=n_conv, n_rigid=nzA)
    return out


def union_spectrum(spectra):
    """Sorted union of the elastic sector spectra (= full-plate elastic spectrum,
    valid up to the smallest per-sector convergence horizon). Returns
    (pooled Lambda, aligned per-mode reference uncertainty, horizon)."""
    horizon = min(s["elastic"][-1] for s in spectra.values())
    lam = np.concatenate([s["elastic"] for s in spectra.values()])
    unc = np.concatenate([s["rel"] for s in spectra.values()])
    order = np.argsort(lam)
    lam, unc = lam[order], unc[order]
    keep = lam <= horizon
    return lam[keep], unc[keep], horizon
