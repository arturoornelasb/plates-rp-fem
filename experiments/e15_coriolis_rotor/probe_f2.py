#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 F2 decisive test: the protection theorem, with the corrected protector
(sigma_v*T, NOT R_pi*T -- proven by build_symop/sym_residuals at 1e-15) and a
parity-ADAPTED modal reduction that removes the cluster-mixing artifact.

Prediction (now with the exact symmetry): under rotation
  D_mirror (sigma_v*T protected)  -> STAYS orthogonal, GOE (~0.53), NOT 0.60
  D_prot   (R_pi, NOT protected)  -> crosses to GUE (~0.60)   [PLAN mislabels it]
  D_chir   (no symmetry)          -> crosses to GUE (~0.60)
"""
import sys, numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r, R_GOE

NU, R_GUE = 0.33, 0.5996
NRINGS, N_TH, NMODES = 16, 48, 300
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
PROT = [(2, 0.18, 0.4), (4, 0.10, 1.1), (6, 0.06, 2.3)]  # even k -> C2/R_pi
CHIR = [(2, 0.15, 0.5), (3, 0.13, 1.7), (4, 0.10, 2.9), (5, 0.08, 0.9)]
cOs = [0.0, 0.5, 1.0, 2.0]


def sweep_r(Lam, G0m):
    sq = np.sqrt(np.mean(np.diff(np.sort(Lam))))
    out = []
    for cO in cOs:
        res = e2.solve_rotor(Lam, G0m, cO * sq)
        r, sem, n = mean_r(res["omega"], skip_low=max(10, len(res["omega"]) // 10))
        out.append((cO, r, sem, res["max_imag"]))
    return out


def prep(harm, kind):
    m, b, meta = e2.star_polar_basis(NRINGS, N_TH, harm)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    lam, X, info, _ = solve_modes(K, M, NMODES)
    if kind is None:
        Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        lab = None
    else:
        S = e2.build_symop(b, kind, tol=1e-6)
        Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S)
    el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
    o = np.argsort(Lam[el])
    Lam_e = Lam[el][o]
    G0m_e = G0m[np.ix_(el, el)][np.ix_(o, o)]
    lab_e = lab[el][o] if lab is not None else None
    return Lam_e, G0m_e, lab_e, info


print(f"refs: GOE {R_GOE}  GUE {R_GUE}\n")

# D_chir: no symmetry -> single class, pooled
Lam, G0m, _, info = prep(CHIR, None)
print(f"[D_chir  (no symmetry -> unitary route F4)]  N={len(Lam)} resid={info['max_resid']:.0e}")
for cO, r, sem, mi in sweep_r(Lam, G0m):
    print(f"   c_Omega={cO:4.1f}  <r>={r:.4f} +/- {sem:.4f}")
print()

# D_mirror: sigma_v anticommutes with G0 -> single class, sigma_v*T PROTECTS -> GOE
Lam, G0m, lab, info = prep(MIRR, "mirror_x")
print(f"[D_mirror (sigma_v*T PROTECTED -> GOE)]  N={len(Lam)} resid={info['max_resid']:.0e}")
for cO, r, sem, mi in sweep_r(Lam, G0m):
    print(f"   c_Omega={cO:4.1f}  <r>={r:.4f} +/- {sem:.4f}")
print()

# D_prot: R_pi COMMUTES with G0 -> block-diagonal -> read PER parity class -> GUE
Lam, G0m, lab, info = prep(PROT, "rpi")
print(f"[D_prot   (R_pi NOT protected -> GUE per class; PLAN mislabels it)]  "
      f"N={len(Lam)} resid={info['max_resid']:.0e}")
for s, nm in [(1, "even"), (-1, "odd")]:
    sel = lab == s
    Ls, Gs = Lam[sel], G0m[np.ix_(sel, sel)]
    line = f"   R_pi-{nm} (n={int(sel.sum())}): "
    for cO, r, sem, mi in sweep_r(Ls, Gs):
        line += f"cO={cO:.1f}:{r:.3f}  "
    print(line)
