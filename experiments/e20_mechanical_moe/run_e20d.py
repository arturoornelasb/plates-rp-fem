#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E20d -- the complete-graph (router-like) dial endpoint (preregistered
in header, 2026-07-10). E20c: at matched position D2 is MONOTONE in
graph degree (0.26 -> 0.38) and posts-per-pair is inert ((3,8) ~ (3,2))
-- the dial variable is connectivity. The AI MoE's router couples every
expert to every other; the mechanical analog is the COMPLETE graph
(n_conn = 99 -> min(n_conn, K-1) = K-1 at every K), which also removes
E20c's K-varying-degree caveat on dense configs. Configs: (99, 2) and
(99, 8), both bisected to pooled r-tilde = 0.460 +/- 0.005; K-sweep D2
as in E20/E20c; M-stability spot-check (M = 80 vs 120, K in {4, 16}
slope) on the denser config, reported as the systematic. FROZEN reading:
DIAL-ENTERS-BAND if either config's D2 lands in [0.56, 0.80] (the
fabrication bridge closes QUANTITATIVELY at router-like coupling);
DIAL-SATURATES-SHORT if both < 0.56 (the mechanical substrate cannot
reach the AI dimension by connectivity alone -- report the saturation
value); ABOVE-BAND if > 0.80."""
import json
import os
import time

import numpy as np

from run_e20 import CFG, build_system, nsec_mid, plate_modes
from run_e20c import bisect_q

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIGS = [(99, 2), (99, 8)]


def d2_sweep(cfg, PHI, lam0, dens, q, Ks, M, seeds):
    ns = {}
    for K in Ks:
        vals = []
        for sd in seeds:
            rng = np.random.default_rng(200 + sd)
            c = dict(cfg)
            c["n_conn"], c["n_posts"] = dens
            H = build_system(lam0, PHI, K, M, q, rng, c)
            vals.append(float(np.mean(nsec_mid(H, K, M,
                                               cfg["window"])[0])))
        ns[K] = float(np.mean(vals))
    D2 = float(np.polyfit(np.log(np.array(Ks, float)),
                          np.log([ns[k] for k in Ks]), 1)[0])
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

    md = ["# E20d -- complete-graph (router-like) dial endpoint "
          "(RESULTS)\n",
          "Complete coupling graph (degree = K-1 at every K); matched "
          "position r-tilde = 0.460 +/- 0.005.\n",
          "| (conn, posts) | q_eq/Delta | r achieved | D2 |",
          "|---|---|---|---|"]
    results = {"configs": []}
    for dens in CONFIGS:
        q_eq, r_ach = bisect_q(cfg, PHI, lam0, dens, d0)
        D2, ns = d2_sweep(cfg, PHI, lam0, dens, q_eq, cfg["K_sweep"],
                          cfg["M"], cfg["seeds"])
        results["configs"].append(dict(dens=list(dens),
                                       q_eq_rel=float(q_eq / d0),
                                       r=r_ach, D2=D2, nsec=ns))
        md.append(f"| {dens} | {q_eq/d0:.4f} | {r_ach:.4f} | {D2:.3f} |")
        print(f"[{dens}] q_eq/D = {q_eq/d0:.4f}, r = {r_ach:.4f}, "
              f"D2 = {D2:.3f} ({time.time()-t00:.0f} s)", flush=True)

    # M-stability spot-check on the denser config
    dens = CONFIGS[1]
    q_eq = results["configs"][1]["q_eq_rel"] * d0
    stab = {}
    for Ma in cfg["M_alt"]:
        D2m, _ = d2_sweep(cfg, PHI, lam0, dens, q_eq, [4, 16], Ma,
                          cfg["seeds"][:3])
        stab[Ma] = D2m
    results["M_stability"] = stab
    md.append(f"\n- M-stability (denser config, 2-point slope): {stab}")

    d2s = [c["D2"] for c in results["configs"]]
    top = max(d2s)
    if any(0.56 <= d <= 0.80 for d in d2s):
        verdict = (f"DIAL-ENTERS-BAND: complete-graph D2 = "
                   f"{[round(d, 2) for d in d2s]} reaches the AI band "
                   f"[0.56, 0.80] at matched position -- the fabrication "
                   f"bridge closes QUANTITATIVELY at router-like "
                   f"coupling.")
    elif top < 0.56:
        verdict = (f"DIAL-SATURATES-SHORT: complete graph gives D2 = "
                   f"{[round(d, 2) for d in d2s]} < 0.56 -- connectivity "
                   f"alone cannot lift the mechanical substrate to the "
                   f"AI dimension (saturation ~ {top:.2f}).")
    else:
        verdict = (f"ABOVE-BAND: D2 = {[round(d, 2) for d in d2s]} "
                   f"overshoots [0.56, 0.80].")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E20D.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e20d.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
