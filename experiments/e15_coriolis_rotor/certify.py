#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 certification: G1 anchor, prestress check, disk control, F2 robustness."""
import sys, numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2, disk_inplane as di
from platefem.kirchhoff import solve_lowest, solve_modes, split_rigid
from platefem.stats import mean_r, r_values, R_GOE

NU, R_GUE = 0.33, 0.5996

# ---- G1: exact in-plane disk anchor vs FEM ----
lam_a, mult = di.inplane_free_disk(NU, omega_max=8.0, n_max=12)
ana = np.sort(di.expand(lam_a, mult))[:24]
m, b = e2.disk_basis(5)
K, M, G0 = e2.assemble_elastic(m, b, NU)
el = np.sort(split_rigid(solve_lowest(K, M, 60), None)[0])[:24]
print(f"G1 anchor (refine 5): first-24 Lambda rel-err median="
      f"{np.median(np.abs(el-ana)/ana):.1e} max={np.max(np.abs(el-ana)/ana):.1e}")

# ---- disk integrable control: class sequence <r> (dedupe +-m doublets) ----
lamA, multA = di.inplane_free_disk(NU, omega_max=20.0, n_max=40)
# one reflection class = each (n) root once
r, sem, n = mean_r(np.sqrt(lamA))
print(f"disk control (integrable): class-sequence <r>={r:.4f} +/- {sem:.4f} (n={n}); "
      f"expect Poisson ~0.39")

# ---- prestress ratio at the crossover (deferred centrifugal ~ (Omega/omega)^2) ----
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
mm, bb, _ = e2.star_polar_basis(16, 48, MIRR)
Kk, Mm, G0m0 = e2.assemble_elastic(mm, bb, NU)
S = e2.build_symop(bb, "mirror_x", tol=1e-6)
lam, X, info, _ = solve_modes(Kk, Mm, 300)
Lam, G0mm, lab, Xn = e2.parity_adapt_reduce(Kk, Mm, G0m0, X, S)
elm = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
Lam_e = np.sort(Lam[elm])
sq = np.sqrt(np.mean(np.diff(Lam_e)))
for cO in [0.5, 1.0, 2.0]:
    Om = cO * sq
    wmed = np.sqrt(np.median(Lam_e))
    print(f"prestress c_Omega={cO}: Omega={Om:.4f} omega_med={wmed:.3f} "
          f"Omega/omega_med={Om/wmed:.4f} (centrifugal ~{(Om/wmed)**2:.1e}, deferred OK if <<1)")

# ---- F2 robustness: higher resolution + second seed, mirror vs chir ----
print("\nF2 robustness (polar 22x72, N=380):")
CHIR = [(2, 0.16, 1.1), (3, 0.12, 2.4), (4, 0.10, 0.7), (5, 0.07, 3.5)]   # new phases
for name, harm, kind in [("D_mirror (sigma_v*T -> GOE)", MIRR, "mirror_x"),
                          ("D_chir   (-> GUE)", CHIR, None)]:
    m2, b2, _ = e2.star_polar_basis(22, 72, harm)
    K2, M2, G2 = e2.assemble_elastic(m2, b2, NU)
    lam2, X2, info2, _ = solve_modes(K2, M2, 380)
    if kind:
        S2 = e2.build_symop(b2, kind, tol=1e-6)
        L2, Gm2, lab2, Xn2 = e2.parity_adapt_reduce(K2, M2, G2, X2, S2)
    else:
        L2, Gm2, Xn2 = e2.modal_reduce(K2, M2, G2, X2)
    e2m = np.abs(L2) > 1e-6 * np.abs(L2).max()
    o = np.argsort(L2[e2m]); Le = L2[e2m][o]; Ge = Gm2[np.ix_(e2m, e2m)][np.ix_(o, o)]
    sq2 = np.sqrt(np.mean(np.diff(np.sort(Le))))
    line = f"  {name}: "
    for cO in [0.0, 1.0, 2.0]:
        res = e2.solve_rotor(Le, Ge, cO * sq2)
        rr, ss, nn = mean_r(res["omega"], skip_low=max(10, len(res["omega"]) // 10))
        line += f"cO={cO}:{rr:.3f}  "
    print(line)
