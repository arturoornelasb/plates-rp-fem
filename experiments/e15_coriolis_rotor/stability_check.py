#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 full-scale operational gate: the strict two-mesh N* fails at this
scale (per-eigenvalue error > 0.1 spacing beyond mode ~50), as expected for
P2 + chord-boundary O(h^2) smooth error. The statistics-grade criterion
(established on the disk in E3b) is REFINEMENT STABILITY OF THE STATISTICS:
repeat the decisive sweeps on the coarse mesh and compare <r>(c)."""
import os
import sys
import numpy as np

sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r

NU, NM = 0.33, 1203
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
CHIR = [(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)]
FINE = {  # from results_full.json (fine 32x96)
    "chir": {0.0: (0.5303, 0.0077), 1.0: (0.5985, 0.0072),
             2.0: (0.5992, 0.0071)},
    "mirror": {0.0: (0.4180, 0.0085), 1.0: (0.5278, 0.0078),
               2.0: (0.5233, 0.0076)},
}

for name, harms, kind in [("chir", CHIR, None), ("mirror", MIRR, "mirror_x")]:
    m, b, _ = e2.star_polar_basis(24, 72, harms)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    lam, X, info, _ = solve_modes(K, M, NM, resid_sanity=1e-3, sweeps_max=30)
    cap = float(np.max(lam))
    if kind is None:
        Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
    else:
        S = e2.build_symop(b, kind, tol=1e-6)
        Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S, lam_cap=cap)
    el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
    o = np.argsort(Lam[el])
    Lam_e = Lam[el][o]
    G0m_e = G0m[np.ix_(el, el)][np.ix_(o, o)]
    sq = np.sqrt(np.mean(np.diff(Lam_e)))
    print(f"[{name}] coarse 24x72: {b.N} dofs, resid {info['max_resid']:.1e}, "
          f"N = {len(Lam_e)}")
    for c in (0.0, 1.0, 2.0):
        out = e2.solve_rotor(Lam_e, G0m_e, c * sq)
        r, sem, n = mean_r(out["omega"], skip_low=max(10, len(out["omega"]) // 10))
        rf, sf = FINE[name][c]
        dsig = (r - rf) / np.hypot(sem, sf)
        print(f"  c_Om={c:g}: coarse <r> = {r:.4f} +/- {sem:.4f} vs fine "
              f"{rf:.4f} +/- {sf:.4f} -> {dsig:+.1f} sigma")
