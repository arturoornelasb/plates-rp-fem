#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19 analysis stage: gates (two-mesh + E5-Argyris cross-instrument +
SS order-stability), Z2xZ2 classification, true-operator ladders,
slopes + Dq, frozen readings. Cheap; runs on the .npz caches."""
import json
import os
import time

import numpy as np

from platefem import n_star
from platefem.stats import classify_parity_resolved, centered_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e19b_common import (A_AX, B_AX, KFEM, LADDER, LEVEL_PROD, NU,
                        build_mesh)

HERE = os.path.dirname(os.path.abspath(__file__))
SECT = ["ee", "eo", "oe", "oo"]
WIN = (0.4, 0.6)
QS = [1.5, 2.0, 4.0]
TRI_REF, RECT_REF = 0.138, 0.011


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
    zf = np.load(os.path.join(HERE, "eig_free_prod.npz"))
    zs = np.load(os.path.join(HERE, "eig_ss_prod.npz"))
    lam_f, V_f = zf["lam"], zf["V"]
    lam_s, V_s = zs["lam"], zs["V"]
    lam_fc = np.load(os.path.join(HERE, "eig_free_check.npz"))["lam"]
    lam_sc = np.load(os.path.join(HERE, "eig_ss_check.npz"))["lam"]

    n1 = min(len(lam_f), len(lam_fc))
    ns_f2m = n_star(lam_f[:n1], lam_fc[:n1], 0.1)
    with open(os.path.join(HERE, "..", "e05_superellipse",
                           "results_raw.json")) as f:
        e5 = json.load(f)
    lam_arg = np.array(e5["runs"]["p10"]["lam"])
    n2 = min(len(lam_f), len(lam_arg))
    ns_cross = n_star(lam_f[:n2], lam_arg[:n2], 0.1)
    n3 = min(len(lam_s), len(lam_sc))
    ns_s_strict = n_star(lam_s[:n3], lam_sc[:n3], 0.1)
    ns_s_order = n_star(lam_s[:n3], lam_sc[:n3], 1.0)
    n_use_f = int(min(ns_f2m, ns_cross if ns_cross < n2 else len(lam_f)))
    n_use_s = int(min(ns_s_order, len(lam_s)))
    gates = dict(free_two_mesh=int(ns_f2m), free_vs_e5=int(ns_cross),
                 n_e5=int(n2), ss_strict=int(ns_s_strict),
                 ss_order=int(ns_s_order),
                 n_use_f=n_use_f, n_use_s=n_use_s)
    print(f"[gates] {gates}")

    mesh = build_mesh(LEVEL_PROD)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    P, Pmx, Pmy = centered_probe_operators(space, A_AX, B_AX, 1500)
    lab_f, qf, _ = classify_parity_resolved(lam_f, V_f, P, Pmx, Pmy,
                                            K.tocsc(), M.tocsc())
    lab_s, qs_, _ = classify_parity_resolved(
        lam_s, V_s, P[:, I], Pmx[:, I], Pmy[:, I],
        K[I][:, I].tocsc(), M[I][:, I].tocsc())
    MVf = (M @ V_f[:, :n_use_f])[I]
    print(f"[labels] free {dict((s, lab_f[:n_use_f].count(s)) for s in SECT)} "
          f"({time.time()-t00:.0f} s)")

    md = ["# E19b -- superellipse p = 10 true-operator windows, init_circle mesh (RESULTS)\n",
          f"Gates: {gates}\n",
          "| sector | N | true IPR |", "|---|---|---|"]
    results = dict(gates=gates, sectors={})
    slopes, ses = [], []
    for s in SECT:
        idx_s = [i for i, l in enumerate(lab_s[:n_use_s]) if l == s]
        idx_f = [i for i, l in enumerate(lab_f[:n_use_f]) if l == s]
        Cmat = (V_s[:, idx_s].T @ MVf[:, idx_f]).T
        Ns, mlns, sess = [], [], []
        dq_acc = {q: [] for q in QS}
        for N in LADDER:
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
            ln2 = np.log(vals[2.0])
            mlns.append(float(np.mean(ln2)))
            sess.append(float(np.std(ln2) / np.sqrt(len(ln2))))
            for q in QS:
                dq_acc[q].append(float(np.mean(np.log(vals[q]))))
            md.append(f"| {s} | {N} | {np.exp(mlns[-1]):.4f} |")
        if len(Ns) < 2:
            results["sectors"][s] = dict(Ns=Ns, coverage_limited=True)
            md.append(f"| {s} | - | COVERAGE-LIMITED ({len(idx_f)} free, "
                      f"{len(idx_s)} basis in gates) |")
            continue
        sl, se = wls(Ns, mlns, sess)
        Dq = {q: float(-np.polyfit(np.log(Ns), dq_acc[q], 1)[0] / (q - 1))
              for q in QS}
        results["sectors"][s] = dict(Ns=Ns, mln=mlns, slope=sl, se=se,
                                     D2=-sl, Dq={str(q): Dq[q] for q in QS},
                                     dq_spread=Dq[1.5] - Dq[4.0])
        slopes.append(sl)
        ses.append(se)
        md.append(f"| {s} | slope | {sl:+.3f} ({se:.3f}) |")
    if not slopes:
        verdict = ("COVERAGE-LIMITED: no sector covers even the first "
                   "rung -- the instrument fails its gates (see gates "
                   "line); instrument still insufficient at refine 6.")
        results["verdict"] = verdict
        md.append("\n**Reading: " + verdict + "**")
        with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(md) + "\n")
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)
        print(verdict)
        return
    w = 1.0 / np.array(ses) ** 2
    sl_pool = float(np.sum(w * np.array(slopes)) / np.sum(w))
    se_pool = float(np.sqrt(1.0 / np.sum(w)))
    sig = abs(sl_pool) / se_pool
    d2 = -sl_pool
    sameside = all(s < 0 for s in slopes)
    md.append(f"\n- pooled slope = {sl_pool:+.3f} ({se_pool:.3f}), "
              f"{sig:.1f} sigma from flat -> D2_true(p10) = {d2:.3f}; "
              f"hierarchy refs: rectangle ~{RECT_REF}, triangle "
              f"{TRI_REF} +/- 0.029")
    dqs = [results["sectors"][s]["dq_spread"] for s in SECT
           if "dq_spread" in results["sectors"][s]]
    md.append(f"- Dq spreads (D_1.5 - D_4): "
              + ", ".join(f"{s}: {results['sectors'][s]['dq_spread']:.3f}"
                          for s in SECT))
    if sameside and sig >= 3:
        cmpr = ("ABOVE the triangle (hierarchy prediction holds)"
                if d2 >= TRI_REF else
                "BELOW the triangle (hierarchy prediction FAILS in order)")
        verdict = (f"THIRD-POINT-CONFIRMED: genuine true-operator "
                   f"delocalization, D2_true = {d2:.3f} at {sig:.1f} "
                   f"sigma -- {cmpr}.")
    elif abs(sl_pool) < 0.05:
        verdict = ("FLAT: the strongest-coupling geometry shows no "
                   "genuine true-operator delocalization -- the E5/E8 "
                   "hierarchy prediction FAILS on the true operator.")
    else:
        verdict = f"INTERMEDIATE ({sig:.1f} sigma; see table)."
    md.append(f"\n**Reading: {verdict}**")
    results["pooled"] = dict(slope=sl_pool, se=se_pool, D2=d2, sigma=sig)
    results["verdict"] = verdict
    md.append(f"\nWall (analysis): {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
