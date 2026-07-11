#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19c analysis: gates (E5-cross primary on free; sigma-20 penalty
order on the SS basis), Z2xZ2 classification, CROSS-MESH projection
(refine-6 SS functions interpolated to refine-7 dof points with the
collar pullback, per-sector nested Cholesky re-orthonormalization --
README, frozen), 3-rung true-operator ladders, frozen readings."""
import json
import os
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular

from platefem import n_star
from platefem.stats import classify_parity_resolved, centered_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e19c_common import (A_AX, B_AX, KFEM, LADDER, LEVEL_FREE, LEVEL_SS,
                         NU, build_mesh, collar_pullback)

HERE = os.path.dirname(os.path.abspath(__file__))
SECT = ["ee", "eo", "oe", "oo"]
WIN = (0.4, 0.6)
QS = [1.5, 2.0, 4.0]
TRI_REF, RECT_REF = 0.138, 0.011
BLK = 256


def wls(Ns, ms, ses):
    x = np.log(Ns)
    w = 1.0 / np.maximum(np.array(ses), 1e-6) ** 2
    W = np.sum(w)
    xb = np.sum(w * x) / W
    yb = np.sum(w * np.array(ms)) / W
    sxx = np.sum(w * (x - xb) ** 2)
    return (float(np.sum(w * (x - xb) * (np.array(ms) - yb)) / sxx),
            float(np.sqrt(1.0 / sxx)))


def main():
    t00 = time.time()
    zf = np.load(os.path.join(HERE, "eig_free7.npz"))
    lam_f, V_f = zf["lam"], zf["V"]
    zs = np.load(os.path.join(HERE, "eig_ss_ext.npz"))
    lam_s, V_s = zs["lam"], zs["V"]
    lam_sig = np.load(os.path.join(HERE, "eig_sig20_ext.npz"))["lam"]

    # ---- gates ----
    with open(os.path.join(HERE, "..", "e05_superellipse",
                           "results_raw.json")) as f:
        e5 = json.load(f)
    lam_arg = np.array(e5["runs"]["p10"]["lam"])
    n2 = min(len(lam_f), len(lam_arg))
    ns_cross = n_star(lam_f[:n2], lam_arg[:n2], 0.1)
    n_use_f = int(ns_cross if ns_cross < n2 else len(lam_f))
    lam_f6 = np.load(os.path.join(HERE, "..", "e19b_superellipse_icirc",
                                  "eig_free_prod.npz"))["lam"]
    n1 = min(len(lam_f), len(lam_f6))
    ns_f2m = n_star(lam_f[:n1], lam_f6[:n1], 0.1)   # informative only
    n4 = min(len(lam_s), len(lam_sig))
    ns_pen = n_star(lam_s[:n4], lam_sig[:n4], 1.0)
    n_use_s = int(min(ns_pen, len(lam_s)))
    gates = dict(free_vs_e5=int(ns_cross), n_e5=int(n2),
                 free_two_mesh_r6_informative=int(ns_f2m),
                 ss_penalty_order=int(ns_pen),
                 n_use_f=n_use_f, n_use_s=n_use_s)
    print(f"[gates] {gates} ({time.time()-t00:.0f} s)", flush=True)

    # ---- spaces and classification ----
    mesh7 = build_mesh(LEVEL_FREE)
    K7, M7, space7 = assemble_c0ip(mesh7, k=KFEM, nu=NU)
    K7 = 0.5 * (K7 + K7.T)
    K7, M7 = K7.tocsc(), M7.tocsc()
    P7, Pmx7, Pmy7 = centered_probe_operators(space7, A_AX, B_AX, 3000)
    lab_f, qf, _ = classify_parity_resolved(lam_f, V_f, P7, Pmx7, Pmy7,
                                            K7, M7)
    print(f"[labels free] "
          f"{dict((s, lab_f[:n_use_f].count(s)) for s in SECT)} "
          f"({time.time()-t00:.0f} s)", flush=True)

    mesh6 = build_mesh(LEVEL_SS)
    K6, M6, space6 = assemble_c0ip(mesh6, k=KFEM, nu=NU)
    K6 = 0.5 * (K6 + K6.T)
    D6 = boundary_dofs(space6)
    I6 = np.setdiff1d(np.arange(space6.N), D6)
    P6, Pmx6, Pmy6 = centered_probe_operators(space6, A_AX, B_AX, 3000)
    lab_s, qs_, _ = classify_parity_resolved(
        lam_s, V_s, P6[:, I6], Pmx6[:, I6], Pmy6[:, I6],
        K6[I6][:, I6].tocsc(), M6[I6][:, I6].tocsc())
    print(f"[labels ss] "
          f"{dict((s, lab_s[:n_use_s].count(s)) for s in SECT)} "
          f"({time.time()-t00:.0f} s)", flush=True)

    # ---- cross-mesh projection operator ----
    pts7 = collar_pullback(np.array(space7.doflocs))
    P67 = space6.probes(pts7)                       # (N7, N6) sparse
    Vfull6 = np.zeros((space6.N, n_use_s))
    Vfull6[I6] = V_s[:, :n_use_s]
    W7 = M7 @ V_f[:, :n_use_f]
    del V_f
    print(f"[projection op] built ({time.time()-t00:.0f} s)", flush=True)

    md = ["# E19c -- superellipse p = 10, refine-7 crisp standalone "
          "(RESULTS)\n", f"Gates: {gates}\n",
          "| sector | N | true IPR |", "|---|---|---|"]
    results = dict(gates=gates, sectors={})
    slopes, ses = [], []
    for s in SECT:
        idx_s = [i for i, l in enumerate(lab_s[:n_use_s]) if l == s]
        idx_f = [i for i, l in enumerate(lab_f[:n_use_f]) if l == s]
        if not idx_s or not idx_f:
            results["sectors"][s] = dict(coverage_limited=True)
            md.append(f"| {s} | - | COVERAGE-LIMITED (0 in gates) |")
            continue
        U = np.empty((space7.N, len(idx_s)))
        for b0 in range(0, len(idx_s), BLK):
            blk = idx_s[b0:b0 + BLK]
            U[:, b0:b0 + len(blk)] = P67 @ Vfull6[:, blk]
        MU = M7 @ U
        G = U.T @ MU
        G = 0.5 * (G + G.T)
        C = U.T @ W7[:, idx_f]
        del U, MU
        try:
            L = cholesky(G, lower=True)
        except np.linalg.LinAlgError:
            results["sectors"][s] = dict(gram_fail=True)
            md.append(f"| {s} | - | GRAM-DISQUALIFIED |")
            continue
        cond = float((np.max(np.diag(L)) / np.min(np.diag(L))) ** 2)
        coeff = solve_triangular(L, C, lower=True)
        Ns, mlns, sess = [], [], []
        dq_acc = {q: [] for q in QS}
        for N in LADDER:
            i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
            if i1 > len(idx_f) or N > len(idx_s):
                break
            vals = {q: [] for q in QS}
            for k in range(i0, i1):
                d = coeff[:N, k]
                p2 = float(d @ d)
                if p2 <= 0:
                    continue
                dn2 = (d / np.sqrt(p2)) ** 2
                for q in QS:
                    vals[q].append(float(np.sum(dn2 ** q)))
            Ns.append(N)
            ln2 = np.log(vals[2.0])
            mlns.append(float(np.mean(ln2)))
            sess.append(float(np.std(ln2) / np.sqrt(len(ln2))))
            for q in QS:
                dq_acc[q].append(float(np.mean(np.log(vals[q]))))
            md.append(f"| {s} | {N} | {np.exp(mlns[-1]):.4f} |")
        if len(Ns) < 2:
            results["sectors"][s] = dict(Ns=Ns, coverage_limited=True,
                                         n_free=len(idx_f),
                                         n_basis=len(idx_s))
            md.append(f"| {s} | - | COVERAGE-LIMITED ({len(idx_f)} free, "
                      f"{len(idx_s)} basis in gates) |")
            continue
        sl, se = wls(Ns, mlns, sess)
        Dq = {q: float(-np.polyfit(np.log(Ns), dq_acc[q], 1)[0] / (q - 1))
              for q in QS}
        results["sectors"][s] = dict(
            Ns=Ns, mln=mlns, slope=sl, se=se, D2=-sl,
            Dq={str(q): Dq[q] for q in QS},
            dq_spread=Dq[1.5] - Dq[4.0], gram_cond=cond,
            n_free=len(idx_f), n_basis=len(idx_s))
        slopes.append(sl)
        ses.append(se)
        md.append(f"| {s} | slope | {sl:+.3f} ({se:.3f}) |")
        print(f"[{s}] slope {sl:+.3f} ({se:.3f}), gram cond {cond:.1e} "
              f"({time.time()-t00:.0f} s)", flush=True)

    if not slopes:
        verdict = "COVERAGE-LIMITED: no sector covers two rungs."
    else:
        w = 1.0 / np.array(ses) ** 2
        sl_pool = float(np.sum(w * np.array(slopes)) / np.sum(w))
        se_pool = float(np.sqrt(1.0 / np.sum(w)))
        sig = abs(sl_pool) / se_pool
        d2 = -sl_pool
        sameside = all(s < 0 for s in slopes)
        md.append(f"\n- pooled slope = {sl_pool:+.3f} ({se_pool:.3f}), "
                  f"{sig:.1f} sigma from flat -> D2_true(p10) = {d2:.3f}; "
                  f"hierarchy refs: rectangle ~{RECT_REF}, triangle "
                  f"{TRI_REF} +/- 0.029; E19b standalone: 0.205 +/- 0.121")
        md.append("- Dq spreads (D_1.5 - D_4): "
                  + ", ".join(f"{s}: {results['sectors'][s]['dq_spread']:.3f}"
                              for s in SECT
                              if "dq_spread" in results["sectors"][s]))
        if sameside and sig >= 3:
            cmpr = ("ABOVE the triangle (hierarchy prediction holds)"
                    if d2 >= TRI_REF else
                    "BELOW the triangle (hierarchy prediction FAILS in "
                    "order)")
            verdict = (f"THIRD-POINT-CONFIRMED: genuine true-operator "
                       f"delocalization, D2_true = {d2:.3f} at {sig:.1f} "
                       f"sigma -- {cmpr}.")
        elif abs(sl_pool) < 0.05:
            verdict = ("FLAT: the strongest-coupling geometry shows no "
                       "genuine true-operator delocalization -- the "
                       "E5/E8 hierarchy prediction FAILS on the true "
                       "operator.")
        else:
            verdict = f"INTERMEDIATE ({sig:.1f} sigma; see table)."
        results["pooled"] = dict(slope=sl_pool, se=se_pool, D2=d2,
                                 sigma=sig)
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    md.append(f"\nWall (analysis): {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
