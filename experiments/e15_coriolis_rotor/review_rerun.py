#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 REVIEW RERUN (2026-07-05): the executed probes fed the FULL
parity-adapted sequence to the statistics, which includes ~170 spurious
noise-parity directions above the certified band (measured benign, <= 0.5
sigma on verdict cells; see review notes in RESULTS). This rerun repeats the
decisive F2/F3 tables with the lam_cap fix (band-truncated reduction) and
quotes n and sem explicitly."""
import sys
import numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r, R_GOE

NU, R_GUE = 0.33, 0.5996
NRINGS, N_TH, NMODES = 16, 48, 300
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
CHIR = [(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)]


def reduced(harmonics, kind):
    m, b, meta = e2.star_polar_basis(NRINGS, N_TH, harmonics)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    lam, X, info, _ = solve_modes(K, M, NMODES)
    cap = float(np.max(lam))
    if kind is None:
        Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        lab = None
    else:
        S = e2.build_symop(b, kind, tol=1e-6)
        Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S, lam_cap=cap)
    el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
    o = np.argsort(Lam[el])
    return (Lam[el][o], G0m[np.ix_(el, el)][np.ix_(o, o)], b,
            Xn[:, el][:, o], info["max_resid"])


def rrow(Lam_e, G0m_e, Mpts=None, cds=(0.0,), cos_=(0.0, 0.5, 1.0, 2.0)):
    dsp = np.mean(np.diff(np.sort(Lam_e)))
    sq = np.sqrt(dsp)
    out = {}
    for cd in cds:
        for cO in cos_:
            res = e2.solve_rotor(Lam_e, G0m_e, cO * sq, Mpts=Mpts,
                                 delta=cd * sq if Mpts is not None else 0.0)
            r, sem, n = mean_r(res["omega"],
                               skip_low=max(10, len(res["omega"]) // 10))
            out[(cd, cO)] = (r, sem, n)
    return out


print("=== F2 with lam_cap (band-truncated), n and sem quoted ===")
Lam_c, G0_c, b_c, X_c, res_c = reduced(CHIR, None)
row = rrow(Lam_c, G0_c)
print(f"D_chir (resid {res_c:.1e}, N={len(Lam_c)}):")
for cO in (0.0, 0.5, 1.0, 2.0):
    r, sem, n = row[(0.0, cO)]
    print(f"  c_Om={cO:3.1f}: <r> = {r:.4f} +/- {sem:.4f} (n={n})")

Lam_m, G0_m, b_m, X_m, res_m = reduced(MIRR, "mirror_x")
row = rrow(Lam_m, G0_m)
print(f"D_mirror protected (resid {res_m:.1e}, N={len(Lam_m)}):")
for cO in (0.0, 0.5, 1.0, 2.0):
    r, sem, n = row[(0.0, cO)]
    print(f"  c_Om={cO:3.1f}: <r> = {r:.4f} +/- {sem:.4f} (n={n})")

print("\n=== F3 two-axis with lam_cap ===")
ang = np.array([0.6, 1.9, 3.3, 5.1])
rad = np.array([0.55, 0.42, 0.63, 0.48])
pts = np.vstack([rad * np.cos(ang), rad * np.sin(ang)])
Mpts = e2.point_mass_modal(b_m, X_m, pts)
Mpts = Mpts / (np.mean(np.abs(np.diag(Mpts))) + 1e-30)
tab = rrow(Lam_m, G0_m, Mpts=Mpts, cds=(0.0, 1.0), cos_=(0.0, 1.0, 2.0))
for cd in (0.0, 1.0):
    line = f"  c_delta={cd:3.1f}: "
    for cO in (0.0, 1.0, 2.0):
        r, sem, n = tab[(cd, cO)]
        line += f"{r:.4f}+/-{sem:.4f}  "
    print(line + f"(n={tab[(cd, 2.0)][2]})")
print(f"\nrefs: GOE {R_GOE}, GUE {R_GUE}")
