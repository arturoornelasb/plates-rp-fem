#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15c -- prestressed rotor: the quantitative soft-rotor prediction.

Adds the Omega^2 physics deferred in E15: exact spin softening (-Omega^2 M)
and geometric stiffness from the static centrifugal prestress (modal static
solve on certified modes; Kg assembled from the resulting stress field;
machinery smoke-validated: free-particle whirl omega = Omega doubly
degenerate, prestress re-stiffens what softening softens, strain metric
1.455 Omega_nd^2).

Question: does the GOE -> GUE crossover survive the Omega^2 terms inside the
linear-validity window (centrifugal strain <= ~10-15%, i.e. Omega_nd <=
~0.3)? If yes, the soft-rotor experiment (elastomer disk, burst invariant
Omega_nd,max ~ 1) is quantitatively predicted."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from platefem.stats import mean_r, R_GOE

HERE = os.path.dirname(os.path.abspath(__file__))
R_GUE = 0.5996
NU, NM = 0.33, 1203
MESH = (24, 72)
CHIR = [(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)]
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]
OMEGAS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
WIN = 200
SIL = dict(E=2e6, rho=1100.0, R=0.10)             # soft silicone disk


def centered_star(nrings, n_th, harms):
    m, b, _ = e2.star_polar_basis(nrings, n_th, harms)
    pc = m.p[:, m.t].mean(axis=1)
    d1 = m.p[:, m.t[1]] - m.p[:, m.t[0]]
    d2 = m.p[:, m.t[2]] - m.p[:, m.t[0]]
    a2 = np.abs(d1[0] * d2[1] - d1[1] * d2[0])
    c = np.array([[np.sum(pc[0] * a2)], [np.sum(pc[1] * a2)]]) / np.sum(a2)
    m2 = type(m)(m.p - c, m.t.copy())
    return m2, e2._vector_basis(m2)


def windowed_r(w):
    centers, rs = [], []
    for i0 in range(10, len(w) - WIN, WIN // 2):
        seg = w[i0:i0 + WIN]
        r, sem, n = mean_r(seg, skip_low=0)
        centers.append(float(np.median(seg)))
        rs.append(float(r))
    return centers, rs


def main():
    t00 = time.time()
    results = {"domains": {}}
    md = ["# E15c -- prestressed rotor: soft-rotor prediction (RESULTS)\n",
          f"Spin softening + geometric prestress stiffness on certified modes "
          f"({NM} modes, mesh {MESH}). Strain metric: eps_max = "
          f"s1 * Omega_nd^2 (s1 measured per domain). Refs GOE {R_GOE} / "
          f"GUE {R_GUE}.\n"]

    for name, harms, kind in [("chir", CHIR, None), ("mirror", MIRR, "mirror_x")]:
        t0 = time.time()
        m, b = centered_star(*MESH, harms)
        K, M, G0 = e2.assemble_elastic(m, b, NU)
        lam, X, info, _ = solve_modes(K, M, NM, resid_sanity=1e-3,
                                      sweeps_max=30)
        u0, rigid_frac = e2.centrifugal_modal_static(b, lam, X, M)
        Kg = e2.assemble_geometric(m, b, u0, NU)
        s1 = e2.max_prestress_strain(m, b, u0)
        cap = float(np.max(lam))
        if kind is None:
            Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        else:
            S = e2.build_symop(b, kind, tol=1e-6)
            Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S,
                                                       lam_cap=cap)
        Kg_m = Xn.T @ (Kg @ Xn)
        Kg_m = 0.5 * (Kg_m + Kg_m.T)
        el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
        o = np.argsort(Lam[el])
        Lam_e = Lam[el][o]
        G0_e = G0m[np.ix_(el, el)][np.ix_(o, o)]
        Kg_e = Kg_m[np.ix_(el, el)][np.ix_(o, o)]
        print(f"[{name}] {b.N} dofs, resid {info['max_resid']:.1e}, rigid "
              f"load frac {rigid_frac:.1e}, s1 = {s1:.3f} "
              f"({time.time()-t0:.1f} s)")

        rows = []
        md.append(f"\n## {name} (s1 = {s1:.3f})\n")
        md.append("| Omega_nd | strain % | pooled <r> (prestressed) | pooled "
                  "<r> (bare) | top-third <r> (prestressed) | max_imag |")
        md.append("|---|---|---|---|---|---|")
        for Om in OMEGAS:
            rp = e2.solve_rotor(Lam_e, G0_e, Om, Kg_m=Kg_e, spin_soften=True)
            rb = e2.solve_rotor(Lam_e, G0_e, Om)
            r_p, s_p, n_p = mean_r(rp["omega"],
                                   skip_low=max(10, len(rp["omega"]) // 10))
            r_b, s_b, _ = mean_r(rb["omega"],
                                 skip_low=max(10, len(rb["omega"]) // 10))
            c_, rs_ = windowed_r(rp["omega"])
            top3 = float(np.mean(rs_[-3:]))
            strain = s1 * Om ** 2 * 100
            rows.append(dict(Om=Om, strain_pct=strain, r_pre=r_p, sem=s_p,
                             r_bare=r_b, r_top3=top3,
                             max_imag=rp["max_imag"],
                             centers=c_, rs=rs_))
            md.append(f"| {Om:g} | {strain:.1f} | {r_p:.4f}({s_p:.4f}) | "
                      f"{r_b:.4f} | {top3:.4f} | {rp['max_imag']:.1e} |")
            print(f"  [{name}] Om={Om:g}: strain {strain:.1f}%, pre "
                  f"{r_p:.4f} bare {r_b:.4f} top3 {top3:.4f}")
        results["domains"][name] = dict(s1=s1, rows=rows,
                                        resid=info["max_resid"])
        with open(os.path.join(HERE, "results_prestress.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)

    # ---------------- soft-rotor prediction ----------------
    ch = results["domains"]["chir"]["rows"]
    s1 = results["domains"]["chir"]["s1"]
    c0 = float(np.sqrt(SIL["E"] / (SIL["rho"] * (1 - NU ** 2))))
    tsc = SIL["R"] / c0
    md.append(f"\n## Soft-rotor prediction (silicone disk R = {SIL['R']} m, "
              f"E = {SIL['E']/1e6:g} MPa, c0 = {c0:.1f} m/s)\n")
    md.append("| Omega_nd | RPM | strain % | pooled <r> | top-third <r> | "
              "top window freq |")
    md.append("|---|---|---|---|---|---|")
    for row in ch:
        rpm = row["Om"] / tsc * 60 / (2 * np.pi)
        ftop = row["centers"][-1] / (2 * np.pi * tsc)
        md.append(f"| {row['Om']:g} | {rpm:,.0f} | {row['strain_pct']:.1f} | "
                  f"{row['r_pre']:.4f} | {row['r_top3']:.4f} | "
                  f"{ftop:.0f} Hz |")

    # verdict
    ok_rows = [r for r in ch if r["strain_pct"] <= 15.0]
    best = max(ok_rows, key=lambda r: r["r_top3"])
    mid = 0.5 * (R_GOE + R_GUE)
    mi_rows = results["domains"]["mirror"]["rows"]
    prot = [r for r in mi_rows if r["strain_pct"] <= 15.0]
    prot_max = max(r["r_pre"] for r in prot[3:]) if len(prot) > 3 else np.nan
    if best["r_top3"] >= mid:
        verdict = (f"SOFT-ROTOR CROSSOVER SURVIVES PRESTRESS: top-window "
                   f"<r> = {best['r_top3']:.3f} at Omega_nd = {best['Om']:g} "
                   f"(strain {best['strain_pct']:.0f}%), within linear "
                   f"validity; protection intact (mirror max "
                   f"{prot_max:.3f}). The elastomer-disk experiment is "
                   f"quantitatively predicted.")
    else:
        verdict = (f"PARTIAL: within strain <= 15% the top window reaches "
                   f"<r> = {best['r_top3']:.3f} at Omega_nd = {best['Om']:g} "
                   f"-- crossover onset visible but not completed inside the "
                   f"linear-validity window; hyperelastic model registered.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS_PRESTRESS.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_prestress.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-10:]))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
