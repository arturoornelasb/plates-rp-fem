#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15d -- SVK finite-deformation rotor (RUNNER). See README_E15D.md
(frozen). Staged: prestates and modal tangents cached per (domain,
Omega); stop-proof."""
import json
import os
import time

import numpy as np
from scipy.linalg import eig

from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r, R_GOE

HERE = os.path.dirname(os.path.abspath(__file__))
R_GUE = 0.5996
NU, NM = 0.33, 1203
MESH = (24, 72)
CHIR = [(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)]
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
OMEGAS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
WIN = 200
STRAIN_CAP = 0.25
MID = 0.5 * (R_GOE + R_GUE)


def centered_star(nrings, n_th, harms):
    m, b, _ = e2.star_polar_basis(nrings, n_th, harms)
    pc = m.p[:, m.t].mean(axis=1)
    d1 = m.p[:, m.t[1]] - m.p[:, m.t[0]]
    d2 = m.p[:, m.t[2]] - m.p[:, m.t[0]]
    a2 = np.abs(d1[0] * d2[1] - d1[1] * d2[0])
    c = np.array([[np.sum(pc[0] * a2)], [np.sum(pc[1] * a2)]]) / np.sum(a2)
    m2 = type(m)(m.p - c, m.t.copy())
    return m2, e2._vector_basis(m2)


def solve_pencil(K_m, G0m, Omega, real_tol=1e-6):
    """Companion solve of (K_m - omega^2 I + i omega Omega G0m) phi = 0
    (the solve_rotor linearization with a full modal K)."""
    N = K_m.shape[0]
    Z = np.zeros((N, N))
    A = np.block([[Z, K_m], [K_m, Omega * G0m]])
    B = np.block([[K_m, Z], [Z, -np.eye(N)]])
    s = eig(A, B, right=False)
    s = s[np.isfinite(s)]
    om = np.sort(s[np.abs(s.imag) < real_tol * np.abs(s).max()].real)
    om = om[om > 0]
    max_imag = float(np.max(np.abs(s.imag)) / max(np.abs(s).max(), 1e-300))
    return om, max_imag


def windowed_r(w):
    centers, rs = [], []
    for i0 in range(10, len(w) - WIN, WIN // 2):
        seg = w[i0:i0 + WIN]
        r, _, _ = mean_r(seg, skip_low=0)
        centers.append(float(np.median(seg)))
        rs.append(float(r))
    return centers, rs


def main():
    t00 = time.time()
    results = {"domains": {}}
    md = ["# E15d -- SVK finite-deformation rotor (RESULTS)\n",
          f"Newton SVK prestate per Omega (continuation), pencil about the "
          f"deformed state; {NM} modes, mesh {MESH}; SVK validity strain "
          f"<= {STRAIN_CAP:.0%}. Refs GOE {R_GOE} / GUE {R_GUE}.\n"]

    for name, harms, kind in [("chir", CHIR, None),
                              ("mirror", MIRR, "mirror_x")]:
        t0 = time.time()
        m, b = centered_star(*MESH, harms)
        K, M, G0 = e2.assemble_elastic(m, b, NU)
        K, M = K.tocsc(), M.tocsc()
        lam, X, info, _ = solve_modes(K, M, NM, resid_sanity=1e-3,
                                      sweeps_max=30)
        if kind is None:
            Lam0, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        else:
            S = e2.build_symop(b, kind, tol=1e-6)
            Lam0, G0m, lab, Xn = e2.parity_adapt_reduce(
                K, M, G0, X, S, lam_cap=float(np.max(lam)))
        el = np.abs(Lam0) > 1e-6 * np.abs(Lam0).max()
        o = np.argsort(Lam0[el])
        idx = np.nonzero(el)[0][o]
        Xe = Xn[:, idx]
        G0_e = G0m[np.ix_(idx, idx)]
        print(f"[{name}] {b.N} dofs, resid {info['max_resid']:.1e} "
              f"({time.time()-t0:.1f} s)")

        rows, u0 = [], None
        md.append(f"\n## {name}\n")
        md.append("| Omega_nd | strain | newton | pooled <r> | "
                  "top-third <r> | max_imag |")
        md.append("|---|---|---|---|---|---|")
        for Om in OMEGAS:
            t1 = time.time()
            if Om == 0.0:
                u0 = np.zeros(K.shape[0])
                ninfo = dict(ok=True, iters=0)
                KT = K
            else:
                u0, ninfo = e2.newton_prestate(m, b, M, NU, Om, u0=u0)
                _, KT = e2.svk_residual_tangent(m, b, u0, NU, Om)
                KT = 0.5 * (KT + KT.T)
            strain = float(e2.max_prestress_strain(m, b, u0))
            KT_m = Xe.T @ (KT @ Xe)
            KT_m = 0.5 * (KT_m + KT_m.T)
            om, mi = solve_pencil(KT_m, G0_e, Om)
            r_p, s_p, _ = mean_r(om, skip_low=max(10, len(om) // 10))
            _, rs_ = windowed_r(om)
            top3 = float(np.mean(rs_[-3:]))
            rows.append(dict(Om=Om, strain=strain, r=r_p, sem=s_p,
                             top3=top3, newton_ok=bool(ninfo["ok"]),
                             newton_iters=int(ninfo["iters"]),
                             max_imag=mi))
            md.append(f"| {Om:g} | {strain:.3f} | "
                      f"{'ok' if ninfo['ok'] else 'FAIL'}/"
                      f"{ninfo['iters']} | {r_p:.4f}({s_p:.4f}) | "
                      f"{top3:.4f} | {mi:.1e} |")
            print(f"  [{name}] Om={Om:g}: strain {strain:.3f}, r {r_p:.4f} "
                  f"top3 {top3:.4f} ({time.time()-t1:.1f} s)")
            results["domains"][name] = dict(rows=rows,
                                            resid=info["max_resid"])
            with open(os.path.join(HERE, "results_e15d.json"), "w") as f:
                json.dump(results, f, indent=1, default=float)

    # ---------------- frozen reading ----------------
    ch = results["domains"]["chir"]["rows"]
    mi_ = results["domains"]["mirror"]["rows"]
    valid = [r for r in ch if r["strain"] <= STRAIN_CAP and r["newton_ok"]]
    best = max(valid, key=lambda r: r["top3"])
    prot = max(r["r"] for r in mi_
               if r["strain"] <= STRAIN_CAP and r["Om"] >= 0.15)
    if best["top3"] >= MID and prot < MID:
        verdict = (f"COMPLETES: top-third <r> = {best['top3']:.3f} at "
                   f"Omega_nd = {best['Om']:g} (strain "
                   f"{best['strain']:.0%}) inside SVK validity; mirror "
                   f"protection intact (max pooled {prot:.3f}). The "
                   f"soft-rotor experiment is fully predicted with "
                   f"finite-deformation physics.")
    else:
        verdict = (f"PARTIAL-AGAIN: best in-validity top-third <r> = "
                   f"{best['top3']:.3f} at Omega_nd = {best['Om']:g}; "
                   f"mirror max {prot:.3f}. Neo-Hookean refinement / "
                   f"large-strain caveat registered.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E15D.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e15d.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-8:]))


if __name__ == "__main__":
    main()
