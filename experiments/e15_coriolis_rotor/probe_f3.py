#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 F3 (the P6 headline): mistuning BREAKS the sigma_v*T protection.

Take the protected mirror rotor (D_mirror, GOE at all Omega). Add point-mass
mistuning at sigma_v-ASYMMETRIC positions -> couples the parity blocks -> destroys
sigma_v*T -> the rotor crosses to GUE. This is 'chirally mistuned rotating body ->
unitary class'. At Omega=0 the mistuned mass is real symmetric -> stays GOE-side;
only rotation + broken signature reaches GUE.
"""
import sys, numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r, R_GOE

NU, R_GUE = 0.33, 0.5996
NRINGS, N_TH, NMODES = 16, 48, 300
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]

m, b, meta = e2.star_polar_basis(NRINGS, N_TH, MIRR)
K, M, G0 = e2.assemble_elastic(m, b, NU)
S = e2.build_symop(b, "mirror_x", tol=1e-6)
lam, X, info, _ = solve_modes(K, M, NMODES)
Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S)
el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
o = np.argsort(Lam[el])
Lam_e, G0m_e = Lam[el][o], G0m[np.ix_(el, el)][np.ix_(o, o)]
Xe = Xn[:, el][:, o]

# sigma_v-ASYMMETRIC mistuning points (generic, not mirror pairs): breaks sigma_v
rng = np.random.default_rng(11)
ang = np.array([0.6, 1.9, 3.3, 5.1])            # not symmetric about x-axis
rad = np.array([0.55, 0.42, 0.63, 0.48])
pts = np.vstack([rad * np.cos(ang), rad * np.sin(ang)])
Mpts = e2.point_mass_modal(b, Xe, pts)
Mpts = Mpts / (np.mean(np.abs(np.diag(Mpts))) + 1e-30)   # ~unit diagonal scale
dsp = np.mean(np.diff(np.sort(Lam_e)))
sq = np.sqrt(dsp)

print(f"refs: GOE {R_GOE}  GUE {R_GUE}\n")
print(f"D_mirror, mistuning at 4 sigma_v-asymmetric points; N={len(Lam_e)}\n")
print("  c_delta \\ c_Omega:   0.0      1.0      2.0")
for cd in [0.0, 0.5, 1.0, 2.0, 4.0]:
    delta = cd * sq                              # mistuning strength
    row = f"  c_delta={cd:4.1f} :  "
    for cO in [0.0, 1.0, 2.0]:
        res = e2.solve_rotor(Lam_e, G0m_e, cO * sq, Mpts=Mpts, delta=delta)
        r, sem, n = mean_r(res["omega"], skip_low=max(10, len(res["omega"]) // 10))
        row += f"{r:.3f}    "
    print(row)
print("\n(expect: c_delta=0 row stays GOE ~0.53 at all Omega [protected];")
print(" strong c_delta + rotation -> GUE ~0.60 [protection broken by mistuning])")
