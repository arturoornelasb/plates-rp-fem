#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18 analysis stage: gates, C3v classification, true-vs-truncated
ladders per sector, frozen-dichotomy verdict. Runs after the four solve
jobs have cached their eigenpairs (cheap; re-runnable)."""
import json
import os
import time

import numpy as np

from platefem import n_star, triangle_ss_exact
from platefem.stats import classify_c3v, triangle_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e18_common import (KFEM, LADDER_A, LADDER_E, L_SIDE, NU, N_PROD,
                        TRI_SECTORS, tri_mesh_bary)

HERE = os.path.dirname(os.path.abspath(__file__))
SPACING_FRAC = 0.1
WIN = (0.4, 0.6)
SLOPE_FLAT, SLOPE_FALL = 0.15, -0.15


def load(job):
    return np.load(os.path.join(HERE, f"eig_{job}.npz"))


def window_ipr(Cmat, N, win=WIN):
    """Window-mean ln IPR of rows of Cmat (modes x basis), truncated to the
    first N basis functions, renormalized. Returns (mlnipr, se, n_win)."""
    i0, i1 = int(win[0] * N), int(win[1] * N)
    if i1 > Cmat.shape[0] or N > Cmat.shape[1]:
        return None
    iprs = []
    for kk in range(i0, i1):
        d = Cmat[kk, :N]
        p = float(d @ d)
        if p <= 0:
            continue
        dn = d / np.sqrt(p)
        iprs.append(float(np.sum(dn ** 4)))
    ln = np.log(iprs)
    return (float(np.mean(ln)), float(np.std(ln) / np.sqrt(len(ln))),
            len(iprs))


def fit_slope(Ns, ms, ses):
    """Weighted LS slope of mlnipr vs ln N, with its standard error."""
    x = np.log(Ns)
    w = 1.0 / np.maximum(np.array(ses), 1e-6) ** 2
    W = np.sum(w)
    xb = np.sum(w * x) / W
    yb = np.sum(w * np.array(ms)) / W
    sxx = np.sum(w * (x - xb) ** 2)
    slope = float(np.sum(w * (x - xb) * (np.array(ms) - yb)) / sxx)
    se = float(np.sqrt(1.0 / sxx))
    return slope, se


def main():
    t00 = time.time()
    results = {"gates": {}, "sectors": {}}

    zf, zs = load("free_prod"), load("ss_prod")
    lam_f, V_f = zf["lam"], zf["V"]
    lam_s, V_s = zs["lam"], zs["V"]
    lam_fc, lam_sc = load("free_check")["lam"], load("ss_check")["lam"]

    # ---------------- gates ----------------
    n1 = min(len(lam_f), len(lam_fc))
    ns_f2m = n_star(lam_f[:n1], lam_fc[:n1], SPACING_FRAC)
    n2 = min(len(lam_s), len(lam_sc))
    ns_s2m = n_star(lam_s[:n2], lam_sc[:n2], SPACING_FRAC)
    ex = triangle_ss_exact(L_SIDE, len(lam_s))
    ns_anchor = n_star(lam_s, ex, SPACING_FRAC)
    with open(os.path.join(
            HERE, "..", "e03_geometries", "triangle",
            "results_raw.json")) as f:
        e3a = json.load(f)
    lam_arg = np.array(e3a["runs"]["free"]["lam"])
    n3 = min(len(lam_f), len(lam_arg))
    ns_arg = n_star(lam_f[:n3], lam_arg[:n3], SPACING_FRAC)
    results["gates"] = dict(free_two_mesh=int(ns_f2m),
                            ss_two_mesh=int(ns_s2m),
                            ss_exact_anchor=int(ns_anchor),
                            free_vs_e3a_argyris=int(ns_arg), n_e3a=int(n3))
    n_use_f = int(min(ns_f2m, len(lam_f)))
    n_use_s = int(min(ns_s2m, ns_anchor, len(lam_s)))
    print(f"[gates] free 2-mesh {ns_f2m}, ss 2-mesh {ns_s2m}, ss exact "
          f"{ns_anchor}, free-vs-Argyris {ns_arg}/{n3} -> n_use free "
          f"{n_use_f}, ss {n_use_s}")

    # ---------------- classification ----------------
    t0 = time.time()
    mesh = tri_mesh_bary(N_PROD)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    P, P_R, P_s = triangle_probe_operators(space, L_SIDE, 1500)
    lab_f, qual_f, _ = classify_c3v(lam_f, V_f, P, P_R, P_s,
                                    K.tocsc(), M.tocsc())
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    lab_s, qual_s, _ = classify_c3v(lam_s, V_s, P[:, I], P_R[:, I],
                                    P_s[:, I], K[I][:, I].tocsc(),
                                    M[I][:, I].tocsc())
    cf = {s: lab_f[:n_use_f].count(s) for s in TRI_SECTORS + ["xx"]}
    cs = {s: lab_s[:n_use_s].count(s) for s in TRI_SECTORS + ["xx"]}
    results["labels"] = dict(free=cf, ss=cs,
                             min_q_free=float(np.min(qual_f)),
                             min_q_ss=float(np.min(qual_s)))
    print(f"[labels] free {cf}, ss {cs} ({time.time()-t0:.1f} s)")

    # ---------------- per-sector ladders ----------------
    t0 = time.time()
    # M-inner products of free modes with zero-extended SS modes: compute
    # the FULL M product first (boundary-to-interior coupling matters for
    # free modes), then slice to interior rows.
    MVf = (M @ V_f[:, :n_use_f])[I]              # (n_interior, n_use_f)
    md = ["# E18 -- unadapted-tier protocol scrutiny: the triangle "
          "(RESULTS)\n",
          f"C0-IP P4, production n = {N_PROD} (~{space.N} dofs), gates: "
          f"{results['gates']}\n"]
    verdicts = {}
    for s in TRI_SECTORS:
        idx_s = [i for i, l in enumerate(lab_s[:n_use_s]) if l == s]
        idx_f = [i for i, l in enumerate(lab_f[:n_use_f]) if l == s]
        n_bas, n_fr = len(idx_s), len(idx_f)
        ladder = LADDER_E if s == "E" else LADDER_A
        # coefficient matrix: free modes (rows) x ss basis (cols), M-inner
        Vs_sec = V_s[:, idx_s]                       # interior dofs
        Cmat = (Vs_sec.T @ MVf[:, idx_f]).T          # (n_fr, n_bas)
        # truncated operator in the sector basis (zero-extended = slice K)
        T = Vs_sec.T @ (K[I][:, I] @ Vs_sec)
        T = 0.5 * (T + T.T)
        rows_t, rows_u = [], []
        for N in ladder:
            r_true = window_ipr(Cmat, N)
            if N <= n_bas:
                w, U = np.linalg.eigh(T[:N, :N])
                i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
                ln = np.log(np.sum(U[:, i0:i1] ** 4, axis=0))
                r_trunc = (float(np.mean(ln)),
                           float(np.std(ln) / np.sqrt(len(ln))), i1 - i0)
            else:
                r_trunc = None
            rows_t.append((N, r_true))
            rows_u.append((N, r_trunc))
        covered = [(N, rt, ru) for (N, rt), (_, ru) in zip(rows_t, rows_u)
                   if rt is not None and ru is not None]
        rec = dict(n_basis=n_bas, n_free=n_fr,
                   ladder=[dict(N=N, true_mln=rt[0], true_se=rt[1],
                                trunc_mln=ru[0], trunc_se=ru[1])
                           for N, rt, ru in covered])
        if len(covered) >= 2:
            Ns = [c[0] for c in covered]
            st, se_t = fit_slope(Ns, [c[1][0] for c in covered],
                                 [c[1][1] for c in covered])
            su, se_u = fit_slope(Ns, [c[2][0] for c in covered],
                                 [c[2][1] for c in covered])
            comb = float(np.sqrt(se_t ** 2 + se_u ** 2))
            if abs(st) < SLOPE_FLAT and su < SLOPE_FALL:
                v = "PROTOCOL ARTIFACT"
            elif su < SLOPE_FALL and abs(st - su) <= 2 * comb:
                v = "GENUINE SCALING"
            else:
                v = "INTERMEDIATE"
            rec.update(slope_true=st, se_true=se_t, slope_trunc=su,
                       se_trunc=se_u, verdict=v,
                       d2_true=-st, d2_trunc=-su)
            verdicts[s] = v
        else:
            rec["verdict"] = "COVERAGE-LIMITED"
            verdicts[s] = "COVERAGE-LIMITED"
        results["sectors"][s] = rec
        md.append(f"\n## sector {s} (basis {n_bas}, free {n_fr})\n")
        md.append("| N | true IPR | trunc IPR |")
        md.append("|---|---|---|")
        for N, rt, ru in covered:
            md.append(f"| {N} | {np.exp(rt[0]):.4f} | "
                      f"{np.exp(ru[0]):.4f} |")
        if "slope_true" in rec:
            md.append(f"\n- slope_true = {rec['slope_true']:+.3f} "
                      f"({rec['se_true']:.3f}); slope_trunc = "
                      f"{rec['slope_trunc']:+.3f} ({rec['se_trunc']:.3f}) "
                      f"-> D2_true = {rec['d2_true']:.3f}, D2_trunc = "
                      f"{rec['d2_trunc']:.3f}")
        md.append(f"- **{rec['verdict']}**")
        print(f"[{s}] {rec.get('verdict')} "
              f"(true {rec.get('slope_true', float('nan')):+.3f}, trunc "
              f"{rec.get('slope_trunc', float('nan')):+.3f})")

    arts = sum(v == "PROTOCOL ARTIFACT" for v in verdicts.values())
    gens = sum(v == "GENUINE SCALING" for v in verdicts.values())
    if gens == 0 and arts >= 2:
        overall = ("PROTOCOL ARTIFACT on the unadapted tier: combined with "
                   "E14/E17, Gap A closes negatively for ALL tiers at "
                   "accessible N; the superellipse contingency does NOT "
                   "trigger")
    elif gens >= 1:
        overall = (f"GENUINE SCALING in {gens} sector(s) -- the unadapted "
                   f"tier carries a real RP-candidate eigenvector phase; "
                   f"superellipse contingency TRIGGERS")
    else:
        overall = "MIXED/INTERMEDIATE -- see per-sector table"
    md.append(f"\n**Reading: {overall}**")
    results["verdict"] = dict(per_sector=verdicts, overall=overall)
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time (analysis): {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
