#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 first physics probe: does a rotating in-plane elastic rotor cross GOE->GUE?

Focus on D_chir (odd+even harmonics, NO C2, NO mirror): a single symmetry class,
so <r> is clean without sector separation -- the cleanest demonstration of the
GEOMETRIC route (F4). Baseline at Omega=0 should be GOE-side (in-plane
elastodynamic chaos); rotation alone (no masses) should move it toward GUE since
no antiunitary signature survives (no C2, no mirror).

Empirical protection cross-check: D_prot (C2, R_pi) and D_mirror (one mirror) at
delta=0 under rotation -- but these carry 2 symmetry classes each, so the pooled
<r> is depressed by class superposition; reported here only as a coarse look
(clean protection test needs class separation -> follow-up).
"""
import sys, time
import numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes, split_rigid
from platefem.stats import mean_r, R_POISSON, R_GOE

NU = 0.33
R_GUE = 0.5996
NREF = 4
NMODES = 380

DOMS = {
    "D_chir":   [(1, 0.12, 0.5), (2, 0.14, 1.3), (3, 0.09, 2.7)],   # no C2, no mirror
    "D_prot":   [(2, 0.18, 0.4), (4, 0.10, 1.1), (6, 0.06, 2.3)],   # C2, no mirror
    "D_mirror": [(2, 0.18, 0.0), (4, 0.10, 0.0)],                   # one mirror axis
}


def certify(harm):
    m, b = e2.star_basis(NREF, harm)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    lam, X, info, _ = solve_modes(K, M, NMODES)
    Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
    el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
    Lam_e, G0m_e = Lam[el], G0m[np.ix_(el, el)]
    order = np.argsort(Lam_e)
    Lam_e, G0m_e = Lam_e[order], G0m_e[np.ix_(order, order)]
    return Lam_e, G0m_e, info["max_resid"]


print(f"config: refine {NREF}, {NMODES} certified modes, nu={NU}")
print(f"refs: Poisson {R_POISSON}  GOE {R_GOE}  GUE {R_GUE}\n")

for name, harm in DOMS.items():
    t0 = time.time()
    Lam, G0m, resid = certify(harm)
    dsp = np.mean(np.diff(np.sort(Lam)))
    sq = np.sqrt(dsp)                                  # Omega unit ~ sqrt(spacing)
    row = []
    for cO in [0.0, 0.5, 1.0, 2.0, 4.0]:
        res = e2.solve_rotor(Lam, G0m, cO * sq)
        w = res["omega"]
        r, sem, n = mean_r(w, skip_low=max(10, len(w) // 10))
        row.append((cO, r, sem, res["max_imag"], res["pair_err"]))
    print(f"[{name}]  N_el={len(Lam)}  resid={resid:.1e}  dt={time.time()-t0:.0f}s")
    for cO, r, sem, mi, pe in row:
        print(f"   c_Omega={cO:4.1f}  <r>={r:.4f} +/- {sem:.4f}  "
              f"(imag {mi:.0e}, pair {pe:.0e})")
    print()
