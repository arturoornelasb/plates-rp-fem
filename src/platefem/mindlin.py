#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reissner-Mindlin plate (E11: the paper's thick-plate prediction).

Fields (theta, w); nondimensionalization D = rho h = 1 with thickness ratio
t: shear stiffness S = kappa_s G h = 5(1-nu)/t^2 (kappa_s = 5/6), rotary
inertia t^2/12. Free edges are natural BCs (no constraints); the Kirchhoff
limit is t -> 0. Discretization: P2 rotations + P2 deflection, with the
shear term assembled under REDUCED integration (selective reduced
integration) to avoid locking in the thin limit -- the thin-limit gate
against the certified Kirchhoff spectrum is the arbiter.
"""
import numpy as np

from skfem import Basis, BilinearForm
from skfem.element import ElementComposite, ElementTriP2, ElementVector
from skfem.helpers import ddot, div, dot, grad, sym_grad


def mindlin_element():
    return ElementComposite(ElementVector(ElementTriP2()), ElementTriP2())


def assemble_mindlin(mesh, nu, t, intorder_full=4, intorder_shear=2):
    """Returns (basis, K, M) for the free Mindlin plate at thickness t."""
    S = 5.0 * (1.0 - nu) / (t * t)

    @BilinearForm
    def bend(theta, w, eta, v, _):
        return (1.0 - nu) * ddot(sym_grad(theta), sym_grad(eta)) \
            + nu * div(theta) * div(eta)

    @BilinearForm
    def shear(theta, w, eta, v, _):
        return S * dot(grad(w) - theta, grad(v) - eta)

    @BilinearForm
    def mass(theta, w, eta, v, _):
        return w * v + (t * t / 12.0) * dot(theta, eta)

    basis_f = Basis(mesh, mindlin_element(), intorder=intorder_full)
    basis_r = Basis(mesh, mindlin_element(), intorder=intorder_shear)
    K = (bend.assemble(basis_f) + shear.assemble(basis_r)).tocsc()
    M = mass.assemble(basis_f).tocsc()
    return basis_f, K, M


def w_probe_operator(basis, pts):
    """Point-evaluation operator for the w field of the composite basis:
    a sparse matrix P with (P @ x)[k] = w(pts[:, k]). Uses split_indices /
    split_bases (validated in the E11 smoke: the piston rigid mode probes to
    a constant)."""
    import scipy.sparse as sp
    i_th, i_w = basis.split_indices()
    b_th, b_w = basis.split_bases()
    Pw = b_w.probes(pts)
    emb = sp.lil_matrix((Pw.shape[0], basis.N))
    emb[:, i_w] = Pw
    return emb.tocsr()
