#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18 follow-up analysis: the registered Dq-flatness (RP-vs-PBRM)
separator applied to the triangle's GENUINE true-operator slope -- the
only real delocalization the campaign found. Calibration (E9, registered):
flat/RP-like D_1.5 - D_4 = 0.069-0.100; PBRM rival band 0.20-0.28.
Runs on the cached eigenpairs."""
import json
import os
import time

import numpy as np

from platefem import n_star
from platefem.stats import classify_c3v, triangle_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e18_common import (KFEM, LADDER_A, LADDER_E, L_SIDE, NU, N_PROD,
                        TRI_SECTORS, tri_mesh_bary)

HERE = os.path.dirname(os.path.abspath(__file__))
QS = [1.5, 2.0, 4.0]
WIN = (0.4, 0.6)
FLAT_BAND = (0.069, 0.100)
PBRM_BAND = (0.20, 0.28)


def main():
    t00 = time.time()
    zf = np.load(os.path.join(HERE, "eig_free_prod.npz"))
    zs = np.load(os.path.join(HERE, "eig_ss_prod.npz"))
    lam_f, V_f = zf["lam"], zf["V"]
    lam_s, V_s = zs["lam"], zs["V"]
    n_use_f, n_use_s = 1056, 1977          # E18 gates

    mesh = tri_mesh_bary(N_PROD)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    P, P_R, P_s = triangle_probe_operators(space, L_SIDE, 1500)
    lab_f, _, _ = classify_c3v(lam_f, V_f, P, P_R, P_s, K.tocsc(), M.tocsc())
    lab_s, _, _ = classify_c3v(lam_s, V_s, P[:, I], P_R[:, I], P_s[:, I],
                               K[I][:, I].tocsc(), M[I][:, I].tocsc())
    MVf = (M @ V_f[:, :n_use_f])[I]
    print(f"[setup] labels + products ({time.time()-t00:.0f} s)")

    md = ["# E18 -- Dq flatness of the GENUINE true-operator slope "
          "(RESULTS)\n",
          f"q in {QS}; calibration: flat/RP {FLAT_BAND}, PBRM "
          f"{PBRM_BAND}.\n",
          "| sector | D_1.5 | D_2 | D_4 | D_1.5 - D_4 | class |",
          "|---|---|---|---|---|---|"]
    results = {}
    for s in TRI_SECTORS:
        idx_s = [i for i, l in enumerate(lab_s[:n_use_s]) if l == s]
        idx_f = [i for i, l in enumerate(lab_f[:n_use_f]) if l == s]
        Cmat = (V_s[:, idx_s].T @ MVf[:, idx_f]).T
        ladder = LADDER_E if s == "E" else LADDER_A
        Ns, mln = [], {q: [] for q in QS}
        for N in ladder:
            i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
            if i1 > Cmat.shape[0] or N > Cmat.shape[1]:
                break
            vals = {q: [] for q in QS}
            for k in range(i0, i1):
                d = Cmat[k, :N]
                p2 = float(d @ d)
                if p2 <= 0:
                    continue
                dn2 = (d / np.sqrt(p2)) ** 2
                for q in QS:
                    vals[q].append(float(np.sum(dn2 ** q)))
            Ns.append(N)
            for q in QS:
                mln[q].append(float(np.mean(np.log(vals[q]))))
        Dq = {}
        for q in QS:
            slope = float(np.polyfit(np.log(Ns), mln[q], 1)[0])
            Dq[q] = -slope / (q - 1.0)
        spread = Dq[1.5] - Dq[4.0]
        if spread <= FLAT_BAND[1]:
            cls = "flat/RP-like"
        elif spread >= PBRM_BAND[0]:
            cls = "PBRM-like"
        else:
            cls = "between bands"
        results[s] = dict(Ns=Ns, Dq={str(q): Dq[q] for q in QS},
                          spread=spread, cls=cls)
        md.append(f"| {s} | {Dq[1.5]:.3f} | {Dq[2.0]:.3f} | {Dq[4.0]:.3f} "
                  f"| {spread:.3f} | {cls} |")
        print(f"[{s}] Dq = {[round(Dq[q], 3) for q in QS]}, spread "
              f"{spread:.3f} -> {cls}")

    spreads = [results[s]["spread"] for s in TRI_SECTORS]
    if all(sp <= FLAT_BAND[1] for sp in spreads):
        verdict = ("FLAT Dq: the triangle's weak genuine delocalization "
                   "is RP-LIKE by the registered separator (spread "
                   "within the flat band in every sector) -- the "
                   "unadapted tier's true operator sits in a weak "
                   "RP-type regime, not a PBRM one.")
    elif any(sp >= PBRM_BAND[0] for sp in spreads):
        verdict = ("PBRM-LEANING in at least one sector -- see table; "
                   "the genuine slope's multifractal family differs "
                   "from RP.")
    else:
        verdict = "BETWEEN BANDS -- see table."
    md.append(f"\n**Reading: {verdict}**")
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS_DQ.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_dq.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-4:]))


if __name__ == "__main__":
    main()
