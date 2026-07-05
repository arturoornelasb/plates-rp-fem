#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 stretch: can the Omega=0 in-plane baseline reach GOE?

The pooled baseline was Poisson (~0.38): the in-plane spectrum superposes
dilatational (P) and shear (S) families. Test whether STRONG chaos (rough
boundary -> strong P<->S mode conversion, or a compact Sinai mass inclusion)
couples them into a single GOE (~0.53) spectrum. All domains chiral (no symmetry)
=> single class => clean <r>.
"""
import sys, numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_lowest, solve_modes, split_rigid
from platefem.stats import mean_r, R_GOE, R_POISSON

NU = 0.33
NRINGS, N_TH, NMODES = 18, 60, 340
print(f"refs: Poisson {R_POISSON}  GOE {R_GOE}\n")

# --- roughness sweep: increasingly rough chiral stars ---
rng = np.random.default_rng(7)
configs = {
    "mild  (k2-5, A~0.46)":  [(2, 0.15, 0.5), (3, 0.13, 1.7), (4, 0.10, 2.9), (5, 0.08, 0.9)],
    "rough (k2-7, A~0.66)":  [(2, 0.17, 0.5), (3, 0.15, 1.7), (4, 0.12, 2.9),
                              (5, 0.10, 0.9), (6, 0.08, 3.4), (7, 0.04, 1.2)],
    "vrough(k2-8, A~0.78)":  [(2, 0.18, 0.3), (3, 0.16, 2.1), (4, 0.14, 1.0),
                              (5, 0.12, 3.5), (6, 0.09, 0.6), (7, 0.06, 2.7), (8, 0.03, 1.9)],
}
for name, harm in configs.items():
    try:
        m, b, _ = e2.star_polar_basis(NRINGS, N_TH, harm)
    except ValueError as ex:
        print(f"  {name}: skipped ({ex})"); continue
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    lam = solve_lowest(K, M, NMODES)
    el, _, nr, _ = split_rigid(lam, None)
    r, sem, n = mean_r(np.sqrt(np.abs(el)))
    rmap = 1 + sum(a * np.cos(k * np.linspace(0, 2 * np.pi, 400) + p) for k, a, p in harm)
    print(f"  {name}: <r>_Omega0 = {r:.4f} +/- {sem:.4f}  (min r/R={rmap.min():.2f}, n={n})")

# --- Sinai inclusion: strong compact added mass on a mild chiral star ---
print("\nSinai inclusion (strong compact added mass, full perturbed solve):")
harm = configs["mild  (k2-5, A~0.46)"]
m, b, _ = e2.star_polar_basis(NRINGS, N_TH, harm)
K, M, G0 = e2.assemble_elastic(m, b, NU)
# a few compact mass lumps via a boundary-mass-like penalty at scattered points:
# add w * sum_k (phi.phi at x_k) to M, realized as a low-rank dense modal bump is
# not enough at full rank; instead load nodes nearest chosen points directly.
from scipy.spatial import cKDTree
ix, iy = b.split_indices()
L = b.doflocs[:, ix]
tree = cKDTree(L.T)
pts = np.array([[0.35, -0.15, 0.1], [0.2, 0.3, -0.25]])   # 3 lump centers
import scipy.sparse as sp
for w in [0.0, 30.0, 100.0, 300.0]:
    Madd = M.tolil()
    for c in range(pts.shape[1]):
        d, near = tree.query(pts[:, c], k=8)
        for nd in near:
            Madd[ix[nd], ix[nd]] += w
            Madd[iy[nd], iy[nd]] += w
    Madd = Madd.tocsc()
    lam = solve_lowest(K, Madd, NMODES)
    el, _, nr, _ = split_rigid(lam, None)
    # drop the (few) heavy localized lump modes pushed high; use the bulk window
    ev = np.sort(np.abs(el))
    win = ev[len(ev) // 5: 4 * len(ev) // 5]
    r, sem, n = mean_r(np.sqrt(win))
    print(f"  w={w:6.1f}: bulk <r> = {r:.4f} +/- {sem:.4f} (n={n})")
