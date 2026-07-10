#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E20b -- grid extension + the coupler-density dial (amendment,
2026-07-10; frozen READINGS unchanged). E20's q-grid missed the Poisson
end (weakest point lambda_eff = 0.24, r = 0.46), so q* landed on the
grid edge and the M-stability gate failed there. E20b: (i) extends the
grid downward to resolve the full S-curve and re-selects q* INSIDE it;
(ii) re-measures the K-sweep D2 at the resolved q* with the M-stability
gate; (iii) the DENSITY DIAL -- if the sparse coupler explains
D2_mech < D2_AI, denser contact graphs should RAISE D2 toward the AI
band: n_(conn, posts) in {(3,2), (5,4), (7,8)} at matched
mid-transition coupling."""
import json
import os
import time

import numpy as np

from run_e20 import (CFG, build_system, lambda_eff, nsec_mid, plate_modes,
                     rt_mid)

HERE = os.path.dirname(os.path.abspath(__file__))
Q_EXT = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
DENSITIES = [(3, 2), (5, 4), (7, 8)]


def r_at(q, K, cfg, PHI, lam0, dens=None):
    c = dict(cfg)
    if dens:
        c["n_conn"], c["n_posts"] = dens
    rs = []
    for sd in cfg["seeds"]:
        rng = np.random.default_rng(100 + sd)
        H = build_system(lam0, PHI, K, cfg["M"], q, rng, c)
        rs.extend(rt_mid(np.linalg.eigvalsh(H), cfg["window"]).tolist())
    rs = np.array(rs)
    return float(rs.mean()), float(rs.std() / np.sqrt(len(rs)))


def d2_at(q, cfg, PHI, lam0, dens=None):
    c = dict(cfg)
    if dens:
        c["n_conn"], c["n_posts"] = dens
    ns = {}
    for K in cfg["K_sweep"]:
        vals = []
        for sd in cfg["seeds"]:
            rng = np.random.default_rng(200 + sd)
            H = build_system(lam0, PHI, K, cfg["M"], q, rng, c)
            vals.append(float(np.mean(nsec_mid(H, K, cfg["M"],
                                               cfg["window"])[0])))
        ns[K] = float(np.mean(vals))
    Ks = np.array(cfg["K_sweep"], float)
    D2 = float(np.polyfit(np.log(Ks),
                          np.log([ns[k] for k in cfg["K_sweep"]]), 1)[0])
    return D2, ns


def main():
    t00 = time.time()
    cfg = CFG
    lam0, V, basis = plate_modes()
    rngp = np.random.default_rng(cfg["seed_plate"])
    pts = np.vstack([rngp.uniform(0.06 * cfg["a"], 0.94 * cfg["a"], 400),
                     rngp.uniform(0.06 * cfg["b"], 0.94 * cfg["b"], 400)])
    PHI = np.asarray(basis.probes(pts) @ V[:, :max(cfg["M_alt"])])
    d0 = float(np.median(np.diff(lam0[:cfg["M"]])))
    results = {}

    md = ["# E20b -- grid extension + density dial (RESULTS)\n",
          "## Extended transition curve (K = 8)\n",
          "| q/Delta | pooled r-tilde |", "|---|---|"]
    curve = []
    for qr in Q_EXT:
        r, se = r_at(qr * d0, cfg["K_r1"], cfg, PHI, lam0)
        curve.append(dict(q_rel=qr, r=r, sem=se))
        md.append(f"| {qr:g} | {r:.4f}({se:.4f}) |")
        print(f"[curve q={qr:g}] r = {r:.4f} ({time.time()-t00:.0f} s)",
              flush=True)
    results["curve"] = curve
    mid = 0.5 * (cfg["R_POI"] + cfg["R_GOE"])
    rs = np.array([c["r"] for c in curve])
    qstar = Q_EXT[int(np.argmin(np.abs(rs - mid)))]
    interior = Q_EXT[0] < qstar < Q_EXT[-1]
    md.append(f"\n- resolved q* = {qstar:g} "
              f"({'interior' if interior else 'STILL AT EDGE'}); Poisson "
              f"end reached: r(min q) = {curve[0]['r']:.4f} vs "
              f"{cfg['R_POI']}")

    # ---------------- D2 at the resolved q*, with M gate ----------------
    D2, ns = d2_at(qstar * d0, cfg, PHI, lam0)
    stab = {}
    for Ma in cfg["M_alt"]:
        c = dict(cfg)
        c["M"] = Ma
        vals = {}
        for K in (4, 16):
            vv = []
            for sd in cfg["seeds"][:3]:
                rng = np.random.default_rng(200 + sd)
                H = build_system(lam0, PHI, K, Ma, qstar * d0, rng, cfg)
                vv.append(float(np.mean(nsec_mid(H, K, Ma,
                                                 cfg["window"])[0])))
            vals[K] = float(np.mean(vv))
        stab[Ma] = float((np.log(vals[16]) - np.log(vals[4]))
                         / (np.log(16) - np.log(4)))
    md.append(f"\n## D2 at resolved q* = {qstar:g}\n")
    md.append(f"- **D2_mech = {D2:.3f}**; N_sec: "
              + ", ".join(f"K{k}: {ns[k]:.2f}" for k in cfg["K_sweep"]))
    md.append(f"- M-stability (2-pt): {stab} "
              f"(gate: spread <= 0.1 -> "
              f"{'PASS' if abs(stab[80]-stab[120]) <= 0.1 else 'FAIL'})")
    results["qstar"] = qstar
    results["D2_resolved"] = D2
    results["M_stability"] = stab

    # ---------------- the density dial ----------------
    md.append("\n## The coupler-density dial (at each config's own "
              "mid-transition q)\n")
    md.append("| (n_conn, n_posts) | q_mid | r(q_mid) | D2 |")
    md.append("|---|---|---|---|")
    dial = []
    for dens in DENSITIES:
        best, bq = None, None
        for qr in Q_EXT:
            r, se = r_at(qr * d0, cfg["K_r1"], cfg, PHI, lam0, dens)
            if best is None or abs(r - mid) < abs(best - mid):
                best, bq = r, qr
        D2d, _ = d2_at(bq * d0, cfg, PHI, lam0, dens)
        dial.append(dict(dens=list(dens), q_mid=bq, r=best, D2=D2d))
        md.append(f"| {dens} | {bq:g} | {best:.4f} | {D2d:.3f} |")
        print(f"[dial {dens}] D2 = {D2d:.3f} ({time.time()-t00:.0f} s)",
              flush=True)
    results["dial"] = dial
    d2s = [d["D2"] for d in dial]
    tunable = all(np.diff(d2s) > 0)
    top = d2s[-1]
    md.append(f"\n- monotone in density: {tunable}; densest D2 = "
              f"{top:.3f} (AI band [0.56, 0.80])")
    if tunable and top >= 0.56:
        concl = ("D2 IS AN ENGINEERING DIAL: coupler density tunes the "
                 "fabricated dimension INTO the AI band -- the "
                 "fabrication bridge closes quantitatively with density "
                 "as the matching knob.")
    elif tunable:
        concl = ("D2 rises with coupler density (dial confirmed) but the "
                 "tested densities stop short of the AI band -- extend "
                 "the dial (registered).")
    else:
        concl = ("Density does not tune D2 monotonically -- the "
                 "mech-vs-AI dimension gap is not a simple density "
                 "effect (report).")
    md.append(f"\n**Reading: {concl}**")
    results["verdict"] = concl
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E20B.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e20b.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
