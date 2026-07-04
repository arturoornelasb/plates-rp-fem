#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E11 -- Mindlin thick plates: thickness sweep (RUNNER + ANALYSIS).

Paper Prediction 'thick plates': for thicker plates additional coupling
channels (shear) modify the effective coupling; SUPPORTS = systematic
variation of per-sector <r> with thickness (potentially toward GOE);
CHALLENGES = statistics insensitive to thickness.

Discretization: P2/P2 with selective reduced integration (platefem.mindlin),
validated: 3 rigid modes; rigid trio probes to exact planes; thin limit
converges to the certified Kirchhoff spectrum (gate G1 below quantifies on
the production mesh). FFFF throughout (free edges = natural BCs).
"""
import json
import os
import time

import numpy as np

from skfem import MeshTri
from platefem import (ElementTriArgyris, assemble_plate, mean_r,
                      r_values, rectangle_basis, solve_lowest, solve_modes,
                      split_rigid, SECTORS)
from platefem.mindlin import assemble_mindlin, w_probe_operator
from platefem.stats import classify_parity_resolved

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(96, 60),
    t_gate=0.005, n_gate=203,
    ts=[0.02, 0.05, 0.10, 0.15],
    n_modes=400,
    probe_npts=1500,
)
SKIP = 10


def probes_trio(basis, a, b, npts, seed=13):
    rng = np.random.default_rng(seed)
    pts = np.vstack([rng.uniform(0.02 * a, 0.98 * a, npts),
                     rng.uniform(0.02 * b, 0.98 * b, npts)])
    P = w_probe_operator(basis, pts)
    Pmx = w_probe_operator(basis, np.vstack([a - pts[0], pts[1]]))
    Pmy = w_probe_operator(basis, np.vstack([pts[0], b - pts[1]]))
    return P, Pmx, Pmy


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}, "runs": {}}
    mesh = MeshTri.init_tensor(np.linspace(0, a, cfg["mesh"][0] + 1),
                               np.linspace(0, b, cfg["mesh"][1] + 1))

    # ---------------- G1: thin-limit gate vs certified Kirchhoff ----------------
    t0 = time.time()
    mk, bk = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    Kk, Mk = assemble_plate(mk, bk, nu)
    lam_k = solve_lowest(Kk, Mk, cfg["n_gate"])
    lam_k, _, _, _ = split_rigid(lam_k)
    del Kk, Mk, mk, bk
    basis_f, K, M = assemble_mindlin(mesh, nu, cfg["t_gate"])
    lam_t = solve_lowest(K, M, cfg["n_gate"])
    lam_t, _, nr, _ = split_rigid(lam_t)
    n_cmp = min(len(lam_t), len(lam_k), 200)
    rel = np.abs(lam_t[:n_cmp] - lam_k[:n_cmp]) / lam_k[:n_cmp]
    print(f"[G1] t={cfg['t_gate']}: rigid {nr}; vs Kirchhoff over {n_cmp} "
          f"modes: max relerr {np.max(rel):.2e}, median {np.median(rel):.2e} "
          f"({time.time()-t0:.1f} s)")
    results["gates"] = dict(G1=dict(max_relerr=float(np.max(rel)),
                                    med_relerr=float(np.median(rel)),
                                    rigid=nr, ndof=int(K.shape[0])))
    del K, M

    # ---------------- thickness sweep ----------------
    rows = []
    for t in cfg["ts"]:
        t0 = time.time()
        basis_f, K, M = assemble_mindlin(mesh, nu, t)
        lam, V, sinfo, _ = solve_modes(K, M, cfg["n_modes"] + 3)
        lam, V, n_rigid, _ = split_rigid(lam, V)
        P, Pmx, Pmy = probes_trio(basis_f, a, b, cfg["probe_npts"])
        labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
        lam = lam[:cfg["n_modes"]]
        labels = labels[:cfg["n_modes"]]
        per, rv = {}, []
        for s in SECTORS:
            ev = np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
            per[s] = mean_r(ev, SKIP)
            rv.extend(r_values(ev[SKIP:]).tolist())
        rv = np.array(rv)
        pool = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))
        counts = {s: labels.count(s) for s in SECTORS + ["xx"]}
        rows.append(dict(t=t, per_sector=per, pooled=pool, counts=counts,
                         n_rigid=n_rigid, min_quality=float(np.min(qual)),
                         solve_resid=sinfo["max_resid"],
                         wall_s=round(time.time() - t0, 1)))
        results["runs"][f"t{t:g}"] = rows[-1]
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)
        print(f"[t={t:g}] pooled <r> = {pool[0]:.4f} +/- {pool[1]:.4f}, "
              f"counts {counts}, rigid {n_rigid} ({time.time()-t0:.1f} s)")
        del K, M, V

    # ---------------- report ----------------
    md = ["# E11 -- Mindlin thick plates: thickness sweep (RESULTS)\n",
          f"FFFF rectangle (campaign geometry), P2/P2 SRI, "
          f"{results['gates']['G1']['ndof']} dofs, {cfg['n_modes']} modes/t. "
          f"Kirchhoff FFFF reference at this ladder: ~0.442.\n"]
    md.append(f"- G1 thin-limit gate (t = {cfg['t_gate']}): max relerr "
              f"{results['gates']['G1']['max_relerr']:.2e} (median "
              f"{results['gates']['G1']['med_relerr']:.2e}) vs certified "
              f"Kirchhoff over 200 modes")
    md.append("\n| t | " + " | ".join(SECTORS) + " | pooled |")
    md.append("|---|---|---|---|---|---|")
    for r in rows:
        md.append(f"| {r['t']:g} | "
                  + " | ".join(f"{r['per_sector'][s][0]:.3f}" for s in SECTORS)
                  + f" | **{r['pooled'][0]:.4f} +/- {r['pooled'][1]:.4f}** |")
    pools = [r["pooled"] for r in rows]
    drift = pools[-1][0] - pools[0][0]
    sig = drift / np.sqrt(pools[-1][1] ** 2 + pools[0][1] ** 2)
    md.append(f"\n- end-to-end drift t = {cfg['ts'][0]:g} -> {cfg['ts'][-1]:g}: "
              f"{drift:+.4f} ({sig:+.1f} sigma)")
    if abs(sig) >= 3:
        verdict = (f"SUPPORTS the thick-plate prediction (systematic "
                   f"{'increase' if drift > 0 else 'decrease'} of <r> with "
                   f"thickness)")
    elif abs(sig) < 1.5:
        verdict = "CHALLENGES (statistics insensitive to thickness at this range)"
    else:
        verdict = "AMBIGUOUS (drift below 3 sigma)"
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
