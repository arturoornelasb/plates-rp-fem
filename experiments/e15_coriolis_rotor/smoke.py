#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 Phase-A smoke: validate elastic2d before the certified runs.

Gates checked here (the anchor-free ones):
  S1 G0 antisymmetry               ||G0+G0'|| / ||G0|| ~ 0
  S2 rigid trio at Omega=0         3 near-zero eigenvalues, then a gap
  S3 two-mesh convergence          low in-plane freqs stable refine 3 -> 4
  S4 pencil realness / pairing     solve_rotor: real omega, +-omega symmetric
  S5 star domains build            D_prot (C2), D_mirror, D_chir mesh OK
"""
import sys
import numpy as np

sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_lowest, solve_modes, split_rigid

NU = 0.33


def low_freqs(nrefine, k=40):
    mesh, basis = e2.disk_basis(nrefine)
    K, M, G0 = e2.assemble_elastic(mesh, basis, NU)
    lam = solve_lowest(K, M, k)
    return mesh, basis, K, M, G0, lam


print("=== S1/S2/S3: disk in-plane ===")
mesh3, basis3, K3, M3, G03, lam3 = low_freqs(3, 60)
asym = np.abs(G03 + G03.T).sum() / (np.abs(G03).sum() + 1e-300)
print(f"S1 G0 antisymmetry rel = {asym:.2e}  (want ~0)")
print(f"   ndof = {basis3.N}")

el, _, nr, rmax = split_rigid(lam3, None)
print(f"S2 rigid modes detected = {nr} (want 3); |Lam| rigid max = {rmax:.2e}")
print(f"   lowest 8 Lam = {np.array2string(lam3[:8], precision=4)}")
print(f"   first 5 elastic freqs omega = "
      f"{np.array2string(np.sqrt(np.abs(el[:5])), precision=4)}")

# S3 two-mesh
_, _, _, _, _, lam4 = low_freqs(4, 60)
el4, _, _, _ = split_rigid(lam4, None)
w3 = np.sqrt(np.abs(el[:12]))
w4 = np.sqrt(np.abs(el4[:12]))
relerr = np.abs(w3 - w4) / w4
print(f"S3 two-mesh (refine 3 vs 4) low-12 omega max relerr = {relerr.max():.2e}")

# S4 pencil realness / pairing on a modal reduction
print("\n=== S4: rotating pencil ===")
N = 200
lamN = solve_lowest(K3, M3, N + 5)
# certified vectors for modal reduction
lamv, X, info, _ = solve_modes(K3, M3, N + 5)
Lam, G0m, Xn = e2.modal_reduce(K3, M3, G03, X)
# drop the 3 rigid modes from the modal set
elmask = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
Lam_e, G0m_e = Lam[elmask], G0m[np.ix_(elmask, elmask)]
mean_sp = np.mean(np.diff(np.sort(Lam_e)))
for Om in [0.0, 0.3 * np.sqrt(mean_sp), 1.0 * np.sqrt(mean_sp)]:
    res = e2.solve_rotor(Lam_e, G0m_e, Om)
    print(f"   Omega={Om:8.4f}: n_pos={res['n_pos']} max_imag={res['max_imag']:.1e} "
          f"pair_err={res['pair_err']:.1e}")

# S5 star domains
print("\n=== S5: star domains build ===")
doms = {
    "D_prot  (C2, even k, no mirror)": [(2, 0.18, 0.4), (4, 0.10, 1.1), (6, 0.06, 2.3)],
    "D_mirror (one mirror axis)":      [(2, 0.18, 0.0), (4, 0.10, 0.0)],
    "D_chir  (odd+even, no C2)":       [(1, 0.12, 0.5), (2, 0.14, 1.3), (3, 0.09, 2.7)],
}
for name, harm in doms.items():
    try:
        m, b = e2.star_basis(3, harm)
        Ks, Ms, Gs = e2.assemble_elastic(m, b, NU)
        lw = solve_lowest(Ks, Ms, 20)
        _, _, nrs, _ = split_rigid(lw, None)
        print(f"   {name}: ndof={b.N} rigid={nrs} w1={np.sqrt(np.abs(lw[nrs])):.4f}")
    except Exception as ex:
        print(f"   {name}: FAILED -> {ex}")

print("\nSMOKE DONE")
