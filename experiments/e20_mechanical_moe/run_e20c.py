#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E20c -- the MATCHED-POSITION density dial (preregistered in header,
2026-07-10). E20b's dial was confounded: configs sat at different
positions on their own transition curves. E20c bisects q per config to
a COMMON pooled r-tilde target (0.460 +/- 0.005, mid-transition), then
measures the K-sweep D2 there. Configs: (n_conn, n_posts) in
{(3,2), (5,4), (7,8), (3,8)} -- the last isolates post-RANK from graph
density. FROZEN reading: DIAL-CONFIRMED-INTO-BAND if D2 rises
monotonically with (conn x posts) among the first three AND the densest
>= 0.56; DIAL-CONFIRMED-SHORT if monotone but < 0.56; NOT-MONOTONE
otherwise. The rank config is reported ungated (rank vs density
attribution)."""
import json
import os
import time

import numpy as np

from run_e20 import CFG, build_system, nsec_mid, plate_modes, rt_mid

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIGS = [(3, 2), (5, 4), (7, 8), (3, 8)]
R_TARGET, R_TOL = 0.460, 0.005


def pooled_r(q, cfg, PHI, lam0, dens):
    c = dict(cfg)
    c["n_conn"], c["n_posts"] = dens
    rs = []
    for sd in cfg["seeds"]:
        rng = np.random.default_rng(100 + sd)
        H = build_system(lam0, PHI, cfg["K_r1"], cfg["M"], q, rng, c)
        rs.extend(rt_mid(np.linalg.eigvalsh(H), cfg["window"]).tolist())
    return float(np.mean(rs))


def bisect_q(cfg, PHI, lam0, dens, d0):
    lo, hi = 1e-4 * d0, 3.0 * d0
    r_lo, r_hi = (pooled_r(lo, cfg, PHI, lam0, dens),
                  pooled_r(hi, cfg, PHI, lam0, dens))
    for _ in range(14):
        mid = np.sqrt(lo * hi)
        r_m = pooled_r(mid, cfg, PHI, lam0, dens)
        if abs(r_m - R_TARGET) <= R_TOL:
            return mid, r_m
        if r_m < R_TARGET:
            lo, r_lo = mid, r_m
        else:
            hi, r_hi = mid, r_m
    return mid, r_m


def main():
    t00 = time.time()
    cfg = CFG
    lam0, V, basis = plate_modes()
    rngp = np.random.default_rng(cfg["seed_plate"])
    pts = np.vstack([rngp.uniform(0.06 * cfg["a"], 0.94 * cfg["a"], 400),
                     rngp.uniform(0.06 * cfg["b"], 0.94 * cfg["b"], 400)])
    PHI = np.asarray(basis.probes(pts) @ V[:, :max(cfg["M_alt"])])
    d0 = float(np.median(np.diff(lam0[:cfg["M"]])))

    md = ["# E20c -- matched-position density dial (RESULTS)\n",
          f"All configs bisected to pooled r-tilde = {R_TARGET} +/- "
          f"{R_TOL}.\n",
          "| (conn, posts) | q_eq/Delta | r achieved | D2 |",
          "|---|---|---|---|"]
    results = {"configs": []}
    for dens in CONFIGS:
        q_eq, r_ach = bisect_q(cfg, PHI, lam0, dens, d0)
        ns = {}
        for K in cfg["K_sweep"]:
            vals = []
            for sd in cfg["seeds"]:
                rng = np.random.default_rng(200 + sd)
                c = dict(cfg)
                c["n_conn"], c["n_posts"] = dens
                H = build_system(lam0, PHI, K, cfg["M"], q_eq, rng, c)
                vals.append(float(np.mean(nsec_mid(H, K, cfg["M"],
                                                   cfg["window"])[0])))
            ns[K] = float(np.mean(vals))
        Ks = np.array(cfg["K_sweep"], float)
        D2 = float(np.polyfit(np.log(Ks),
                              np.log([ns[k] for k in cfg["K_sweep"]]),
                              1)[0])
        results["configs"].append(dict(dens=list(dens),
                                       q_eq_rel=float(q_eq / d0),
                                       r=r_ach, D2=D2, nsec=ns))
        md.append(f"| {dens} | {q_eq/d0:.4f} | {r_ach:.4f} | {D2:.3f} |")
        print(f"[{dens}] q_eq/D = {q_eq/d0:.4f}, r = {r_ach:.4f}, "
              f"D2 = {D2:.3f} ({time.time()-t00:.0f} s)", flush=True)

    d2_density = [c["D2"] for c in results["configs"][:3]]
    mono = all(np.diff(d2_density) > 0)
    top = d2_density[-1]
    rank_d2 = results["configs"][3]["D2"]
    if mono and top >= 0.56:
        verdict = (f"DIAL-CONFIRMED-INTO-BAND: matched-position D2 rises "
                   f"{d2_density[0]:.2f} -> {top:.2f} with coupler "
                   f"density, entering the AI band [0.56, 0.80] -- the "
                   f"fabricated dimension is an engineering dial and the "
                   f"bridge closes quantitatively.")
    elif mono:
        verdict = (f"DIAL-CONFIRMED-SHORT: monotone "
                   f"{d2_density[0]:.2f} -> {top:.2f}, below the AI band "
                   f"-- extend the dial (denser graphs) to test entry.")
    else:
        verdict = (f"NOT-MONOTONE at matched position "
                   f"({[round(x, 2) for x in d2_density]}): the dimension "
                   f"gap is not a density effect.")
    md.append(f"\n- rank config (3, 8): D2 = {rank_d2:.3f} (vs (3, 2) "
              f"{d2_density[0]:.3f} -- the RANK contribution)")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E20C.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e20c.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
