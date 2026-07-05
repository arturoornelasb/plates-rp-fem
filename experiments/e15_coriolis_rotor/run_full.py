#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 FULL-SCALE run (registered ladder: ~1200 certified modes/domain).

Same physics and conventions as the executed probes (reviewed 2026-07-05),
at the PLAN's registered scale, with the review lam_cap fix active and all
gates quoted: G2 two-mesh internal N*, G3 modal-truncation stability
(N = 600 vs 1200 at a strong-coupling cell), G4 realness/pairing, prestress
ratio. Decisive questions at this power (sem ~ 0.008):
  (a) D_chir GOE -> GUE crossover at ~8 sigma;
  (b) protection of the mirror rotor at all speeds;
  (c) the two-axis protected-vs-mistuned contrast;
  (d) does the small super-GUE overshoot of the probe runs (+2-3 sigma at
      n ~ 265) persist at n ~ 1100? (finite-n vs real rigidity excess).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", ".."))
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_lowest, solve_modes
from platefem.stats import mean_r, n_star, R_GOE

HERE = os.path.dirname(os.path.abspath(__file__))
R_GUE = 0.5996

CFG = dict(
    nu=0.33, n_modes=1203,
    mesh=(32, 96), mesh_check=(24, 72),          # (nrings, n_th)
    mirr=[(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)],
    chir=[(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)],
    mist_ang=[0.6, 1.9, 3.3, 5.1], mist_rad=[0.55, 0.42, 0.63, 0.48],
    c_omegas=[0.0, 0.25, 0.5, 1.0, 2.0, 4.0],
    two_axis_cd=[0.0, 0.5, 1.0, 2.0], two_axis_co=[0.0, 1.0, 2.0],
    spacing_frac=0.1,
)


def save(res):
    with open(os.path.join(HERE, "results_full.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)


def stats_of(omega):
    r, sem, n = mean_r(omega, skip_low=max(10, len(omega) // 10))
    return dict(r=r, sem=sem, n=n)


def main():
    t00 = time.time()
    cfg = CFG
    res = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                      for k, v in cfg.items()}, "gates": {}, "domains": {}}

    for name, harms, kind in [("chir", cfg["chir"], None),
                              ("mirror", cfg["mirr"], "mirror_x")]:
        t0 = time.time()
        m, b, _ = e2.star_polar_basis(*cfg["mesh"], harms)
        K, M, G0 = e2.assemble_elastic(m, b, cfg["nu"])
        mc, bc, _ = e2.star_polar_basis(*cfg["mesh_check"], harms)
        Kc, Mc, _ = e2.assemble_elastic(mc, bc, cfg["nu"])
        lam_c = solve_lowest(Kc, Mc, cfg["n_modes"])
        lam, X, info, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-3,
                                      sweeps_max=30)
        ns = n_star(np.sort(lam), np.sort(lam_c), cfg["spacing_frac"])
        print(f"[{name}] {b.N} dofs, resid {info['max_resid']:.1e}, "
              f"two-mesh N* = {ns}/{cfg['n_modes']} ({time.time()-t0:.1f} s)")
        cap = float(np.max(lam))
        if kind is None:
            Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        else:
            S = e2.build_symop(b, kind, tol=1e-6)
            sr = e2.sym_residuals(S, K, M, G0)
            print(f"[{name}] symop: K {sr['commute_K']:.1e}, "
                  f"G0 anticommute {sr['G0_anticommute']:.1e}")
            res["gates"][f"{name}_symop"] = sr
            Lam, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S,
                                                       lam_cap=cap)
        el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
        o = np.argsort(Lam[el])
        Lam_e = Lam[el][o]
        G0m_e = G0m[np.ix_(el, el)][np.ix_(o, o)]
        Xe = Xn[:, el][:, o]
        dsp = float(np.mean(np.diff(Lam_e)))
        sq = np.sqrt(dsp)
        dom = {"ndof": int(b.N), "n_levels": int(len(Lam_e)),
               "two_mesh_nstar": int(ns), "resid": info["max_resid"],
               "sweep": {}}

        # prestress ratio at the strongest speed used for the verdict (c=2)
        om_med = float(np.sqrt(np.median(Lam_e)))
        dom["prestress_ratio_c2"] = float((2.0 * sq / om_med) ** 2)

        for cO in cfg["c_omegas"]:
            out = e2.solve_rotor(Lam_e, G0m_e, cO * sq)
            st = stats_of(out["omega"])
            st.update(max_imag=out["max_imag"], pair_err=out["pair_err"])
            dom["sweep"][f"{cO:g}"] = st
            print(f"  [{name}] c_Om={cO:g}: <r> = {st['r']:.4f} +/- "
                  f"{st['sem']:.4f} (n={st['n']}, imag {st['max_imag']:.1e})")

        # G3 truncation stability at c = 1 (lowest 600 of the reduced model)
        half = len(Lam_e) // 2
        out6 = e2.solve_rotor(Lam_e[:half], G0m_e[:half, :half], 1.0 * sq)
        st6 = stats_of(out6["omega"])
        dom["trunc_gate_c1_N600"] = st6
        print(f"  [{name}] G3 (N={half}) c_Om=1: <r> = {st6['r']:.4f} +/- "
              f"{st6['sem']:.4f}")

        if kind is not None:
            ang = np.array(cfg["mist_ang"])
            rad = np.array(cfg["mist_rad"])
            pts = np.vstack([rad * np.cos(ang), rad * np.sin(ang)])
            Mpts = e2.point_mass_modal(b, Xe, pts)
            Mpts = Mpts / (np.mean(np.abs(np.diag(Mpts))) + 1e-30)
            grid = {}
            for cd in cfg["two_axis_cd"]:
                for cO in cfg["two_axis_co"]:
                    out = e2.solve_rotor(Lam_e, G0m_e, cO * sq, Mpts=Mpts,
                                         delta=cd * sq)
                    grid[f"{cd:g}_{cO:g}"] = stats_of(out["omega"])
            dom["two_axis"] = grid
            print(f"  [{name}] two-axis grid done")
        res["domains"][name] = dom
        save(res)

    # ---------------- verdict ----------------
    ch = res["domains"]["chir"]["sweep"]
    mi = res["domains"]["mirror"]["sweep"]
    ta = res["domains"]["mirror"]["two_axis"]
    cross = (ch["1"]["r"] - ch["0"]["r"]) / np.hypot(ch["1"]["sem"],
                                                     ch["0"]["sem"])
    vs_gue_1 = (ch["1"]["r"] - R_GUE) / ch["1"]["sem"]
    vs_gue_2 = (ch["2"]["r"] - R_GUE) / ch["2"]["sem"]
    prot_dev = max(abs(mi[f"{c:g}"]["r"] - R_GOE) / mi[f"{c:g}"]["sem"]
                   for c in [0.5, 1.0, 2.0, 4.0])
    contrast = (ta["1_2"]["r"] - ta["0_2"]["r"]) / np.hypot(
        ta["1_2"]["sem"], ta["0_2"]["sem"])
    md = ["# E15 -- FULL-SCALE RESULTS (registered ladder)\n",
          f"~1200 certified modes/domain, meshes {cfg['mesh']} vs "
          f"{cfg['mesh_check']}, lam_cap review fix active. Refs: GOE "
          f"{R_GOE}, GUE {R_GUE}.\n"]
    md.append("## D_chir (no symmetry): geometric GOE -> GUE\n")
    md.append("| c_Omega | <r> | sem | n |")
    md.append("|---|---|---|---|")
    for c in cfg["c_omegas"]:
        s = ch[f"{c:g}"]
        md.append(f"| {c:g} | {s['r']:.4f} | {s['sem']:.4f} | {s['n']} |")
    md.append(f"\n- crossover 0 -> 1: {cross:+.1f} sigma; vs GUE at c=1: "
              f"{vs_gue_1:+.1f} sigma; at c=2: {vs_gue_2:+.1f} sigma")
    md.append("\n## D_mirror (sigma_v*T protected)\n")
    md.append("| c_Omega | <r> | sem | n |")
    md.append("|---|---|---|---|")
    for c in cfg["c_omegas"]:
        s = mi[f"{c:g}"]
        md.append(f"| {c:g} | {s['r']:.4f} | {s['sem']:.4f} | {s['n']} |")
    md.append(f"\n- max deviation from GOE over c >= 0.5: {prot_dev:.1f} sigma")
    md.append("\n## Two-axis (mirror rotor + asymmetric point masses)\n")
    md.append("| c_delta \\ c_Omega | " + " | ".join(
        f"{c:g}" for c in cfg["two_axis_co"]) + " |")
    md.append("|---|" + "---|" * len(cfg["two_axis_co"]))
    for cd in cfg["two_axis_cd"]:
        md.append(f"| {cd:g} | " + " | ".join(
            f"{ta[f'{cd:g}_{c:g}']['r']:.4f}({ta[f'{cd:g}_{c:g}']['sem']:.4f})"
            for c in cfg["two_axis_co"]) + " |")
    md.append(f"\n- protected-vs-mistuned contrast at c_Omega=2: "
              f"{contrast:+.1f} sigma")
    md.append(f"- gates: two-mesh N* chir "
              f"{res['domains']['chir']['two_mesh_nstar']}, mirror "
              f"{res['domains']['mirror']['two_mesh_nstar']}; G3 truncation "
              f"c=1: chir {res['domains']['chir']['trunc_gate_c1_N600']['r']:.4f}"
              f" vs {ch['1']['r']:.4f}; prestress ratio (c=2) "
              f"{res['domains']['chir']['prestress_ratio_c2']:.2e}")
    ok = cross > 3 and prot_dev < 3 and contrast > 3
    verdict = ("SUPPORTS P6/P11 at full scale (crossover, protection, and "
               "two-axis contrast all resolved)" if ok else
               "CHECK -- see gates/tables")
    md.append(f"\n**Reading: {verdict}**")
    res["verdict"] = verdict
    res["wall_time_s"] = round(time.time() - t00, 1)
    save(res)
    with open(os.path.join(HERE, "RESULTS_FULL.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[done] {res['wall_time_s']} s")


if __name__ == "__main__":
    main()
