#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E16b -- anatomy of the theta ~ 75-deg dip (preregistered in header,
2026-07-09). E16's angle-resolved curve shows a pronounced minimum near
75 deg (binned mean ~0.35, at the Poisson baseline) between two maxima.
FROZEN reading: (i) REPRODUCES-AND-LOCALIZES if the fine sweep (60
angles in [65, 85] deg, same instrument and windows) shows a minimum
deeper than baseline+2 sigma_bin with a localized center; report the
center angle and whether it matches a candidate special angle
(rational-pi fraction / near-separable geometry); (ii) WASHES-OUT if the
fine curve is flat within noise -- the E16 dip was a fluctuation.
Per-class (S/A) resolution reported to see if one symmetry class drives
it. Reuses the frozen E16 machinery unchanged."""
import json
import os
import time

import numpy as np

from run_sweep import CFG as E16CFG, solve_angle, rt_ratios

HERE = os.path.dirname(os.path.abspath(__file__))
TH = np.linspace(65.0, 85.0, 60)


def main():
    t00 = time.time()
    cfg = dict(E16CFG)
    rows = []
    for i, deg in enumerate(TH):
        rec = solve_angle(np.deg2rad(deg), cfg, cfg["nrings"])
        out = {"theta_deg": float(deg)}
        for c in ("S", "A"):
            seq = [l for l, s in zip(rec["lam"], rec["labels"]) if s == c]
            r = rt_ratios(seq, cfg["drop_low"])
            out[c] = float(r.mean()) if len(r) else float("nan")
        both = np.concatenate([rt_ratios(
            [l for l, s in zip(rec["lam"], rec["labels"]) if s == c],
            cfg["drop_low"]) for c in ("S", "A")])
        out["pooled"] = float(both.mean())
        out["sem"] = float(both.std() / np.sqrt(len(both)))
        rows.append(out)
        if (i + 1) % 10 == 0:
            print(f"[dip] {i+1}/{len(TH)} ({time.time()-t00:.0f} s)",
                  flush=True)
            with open(os.path.join(HERE, "results_e16b.json"), "w") as f:
                json.dump(rows, f, indent=1)

    th = np.array([r["theta_deg"] for r in rows])
    rp = np.array([r["pooled"] for r in rows])
    # 5-angle running mean
    kern = np.ones(5) / 5
    sm = np.convolve(rp, kern, mode="valid")
    thm = th[2:-2]
    imin = int(np.argmin(sm))
    base = 0.3865
    sig_bin = float(np.std(rp) / np.sqrt(5))
    depth_ref = float(np.median(sm) - sm[imin])
    md = ["# E16b -- dip anatomy (RESULTS)\n",
          f"60 angles in [65, 85] deg, first 70 modes, same E16 "
          f"instrument.\n",
          f"- smoothed minimum: theta = {thm[imin]:.2f} deg, pooled "
          f"<r~> = {sm[imin]:.4f} (median of window {np.median(sm):.4f}; "
          f"depth {depth_ref:.4f}; pooled-Poisson baseline ~{base})",
          f"- per-class at the minimum: S = "
          f"{rows[imin+2]['S']:.4f}, A = {rows[imin+2]['A']:.4f}"]
    # candidate special angles in-window
    cands = {"72 = 2pi/5": 72.0, "75 = 5pi/12": 75.0, "80 = 4pi/9": 80.0,
             "67.5 = 3pi/8": 67.5}
    near = min(cands.items(), key=lambda kv: abs(kv[1] - thm[imin]))
    md.append(f"- nearest candidate special angle: {near[0]} "
              f"(|delta| = {abs(near[1]-thm[imin]):.2f} deg)")
    if depth_ref > 2 * sig_bin and sm[imin] < base + 0.02:
        verdict = (f"REPRODUCES-AND-LOCALIZES: a real minimum at "
                   f"{thm[imin]:.1f} deg reaching the Poisson baseline "
                   f"(depth {depth_ref:.3f} > 2 sigma_bin "
                   f"{2*sig_bin:.3f}); nearest special angle {near[0]}.")
    elif depth_ref > 2 * sig_bin:
        verdict = (f"REPRODUCES (shallow): minimum at {thm[imin]:.1f} "
                   f"deg, depth {depth_ref:.3f}, above baseline.")
    else:
        verdict = "WASHES-OUT: no minimum beyond 2 sigma_bin."
    md.append(f"\n**Reading: {verdict}**")
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS_E16B.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e16b.json"), "w") as f:
        json.dump(dict(rows=rows, verdict=verdict), f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
