#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E16c -- does the sub-Poisson dip persist DEEP in the spectrum?
(preregistered in header, 2026-07-10). E16b localized a sub-Poisson
minimum at theta ~ 74.5 deg (first 70 modes, nearest candidate 5pi/12).
If the dip marks a quasi-separable/commensurate angle, the level
clustering should PERSIST (or sharpen) deep in the spectrum; if it is a
low-mode boundary effect, it should relax to Poisson. Instrument:
n_keep = 400 modes at nrings = 60 with an nrings = 45 cross-mesh ladder
per angle; certified depth = first index where the two ladders disagree
by > 1e-3 (minus 5 buffer) -- windows are capped there. Angles: the E16b
smoothed dip center (recomputed from its frozen results), 75.0 (5pi/12),
and flanks 73.0 / 77.0. Class-resolved ladders (S/A by bisector mirror);
unclassified modes ('x') break BOTH class sequences at that energy
(segment-wise ratio pooling -- no punctured-ladder spacings). Windows on
the class index: LOW = [10, 35) (the E16/E16b window), DEEP = [50,
min(195, certified cap)). FROZEN reading (pooled r-tilde, Poisson
0.3863): PERSISTS if the dip angle's DEEP window is < Poisson - 2 sigma
while both flanks are not; LOW-MODE-ONLY if its LOW window is below
baseline + 2 sigma but the DEEP window is at/above Poisson - 2 sigma;
NOT-REPRODUCED if the LOW window shows no dip. INSTRUMENT-LIMITED tag if
pooled certified depth < 250 modes. Deeper-of {dip, 75.0} in the DEEP
window reported ungated (angle pinning)."""
import json
import os
import time

import numpy as np

from run_sweep import CFG as E16CFG, solve_angle, rt_ratios

HERE = os.path.dirname(os.path.abspath(__file__))
R_POI = 0.3863
CERT_TOL = 1e-3


def dip_center():
    with open(os.path.join(HERE, "results_e16b.json")) as f:
        rows = json.load(f)["rows"]
    th = np.array([r["theta_deg"] for r in rows])
    rp = np.array([r["pooled"] for r in rows])
    sm = np.convolve(rp, np.ones(5) / 5, mode="valid")
    return float(th[2:-2][int(np.argmin(sm))])


def class_window_ratios(lam, labels, cls, w0, w1):
    """Segment-wise spacing ratios of one class, class index in [w0, w1);
    an 'x' mode breaks the sequence at its energy."""
    segs, cur, ci = [], [], 0
    for l, s in zip(lam, labels):
        if s == cls:
            if w0 <= ci < w1:
                cur.append(l)
            elif ci >= w1:
                break
            ci += 1
        elif s == "x" and cur:
            segs.append(cur)
            cur = []
    if cur:
        segs.append(cur)
    return np.concatenate([rt_ratios(s, 0) for s in segs]
                          + [np.empty(0)])


def main():
    t00 = time.time()
    cfg = dict(E16CFG)
    cfg["n_keep"] = 400
    cfg["probe_npts"] = 2000
    th_dip = dip_center()
    angles = sorted({round(th_dip, 2), 75.0, 73.0, 77.0})
    md = ["# E16c -- deep-spectrum dip persistence (RESULTS)\n",
          f"E16b smoothed dip center: {th_dip:.2f} deg. 400 modes, "
          f"nrings 60 (cross 45), certified-depth-capped windows.\n",
          "| theta (deg) | cert. depth | frac x | LOW r~ (sem) | "
          "DEEP r~ (sem) |", "|---|---|---|---|---|"]
    out = {"th_dip": th_dip, "angles": {}}
    for deg in angles:
        t0 = time.time()
        fine = solve_angle(np.deg2rad(deg), cfg, 60)
        coarse = solve_angle(np.deg2rad(deg), cfg, 45)
        nc = min(len(fine["lam"]), len(coarse["lam"]))
        rel = np.abs(fine["lam"][:nc] - coarse["lam"][:nc]) \
            / np.abs(fine["lam"][:nc])
        bad = np.nonzero(rel > CERT_TOL)[0]
        cert = int((bad[0] if len(bad) else nc) - 5)
        cap = max(0, cert // 2 - 5)
        fx = float(np.mean([s == "x" for s in fine["labels"]]))
        rec = {"cert_depth": cert, "frac_x": fx,
               "resid": fine["resid"], "ndof": fine["ndof"]}
        for wname, (w0, w1) in [("LOW", (10, min(35, cap))),
                                ("DEEP", (50, min(195, cap)))]:
            if w1 - w0 < 15:
                rec[wname] = None
                continue
            r = np.concatenate([class_window_ratios(
                fine["lam"], fine["labels"], c, w0, w1)
                for c in ("S", "A")])
            rec[wname] = dict(r=float(r.mean()),
                              sem=float(r.std() / np.sqrt(len(r))),
                              n=int(len(r)), win=[w0, w1])
        out["angles"][f"{deg:g}"] = rec
        fmt = lambda w: (f"{w['r']:.4f} ({w['sem']:.4f})"
                         if w else "insufficient")
        md.append(f"| {deg:g} | {cert} | {fx:.3f} | {fmt(rec['LOW'])} | "
                  f"{fmt(rec['DEEP'])} |")
        print(f"[{deg:g} deg] cert {cert}, LOW {fmt(rec['LOW'])}, DEEP "
              f"{fmt(rec['DEEP'])} ({time.time()-t0:.0f} s)", flush=True)
        with open(os.path.join(HERE, "results_e16c.json"), "w") as f:
            json.dump(out, f, indent=1, default=float)

    dip = out["angles"][f"{round(th_dip, 2):g}"]
    flanks = [out["angles"]["73"], out["angles"]["77"]]
    lim = min(a["cert_depth"] for a in out["angles"].values()) < 250
    if lim:
        md.append("\n- INSTRUMENT-LIMITED: certified depth < 250 at "
                  "some angle; DEEP windows capped accordingly.")
    verdict = "INSTRUMENT-LIMITED (no certified DEEP window at the dip)"
    if dip["DEEP"] and all(f["DEEP"] for f in flanks):
        sub_deep = dip["DEEP"]["r"] < R_POI - 2 * dip["DEEP"]["sem"]
        fl_ok = all(f["DEEP"]["r"] >= R_POI - 2 * f["DEEP"]["sem"]
                    for f in flanks)
        low_dip = (dip["LOW"] and
                   dip["LOW"]["r"] < R_POI + 2 * dip["LOW"]["sem"])
        if sub_deep and fl_ok:
            verdict = (f"PERSISTS: the dip angle stays sub-Poisson deep "
                       f"in the spectrum (DEEP r~ = {dip['DEEP']['r']:.4f}"
                       f" < {R_POI}) while the flanks do not -- a genuine "
                       f"special-angle clustering, not a low-mode effect.")
        elif low_dip:
            verdict = (f"LOW-MODE-ONLY: the dip reproduces in the LOW "
                       f"window (r~ = {dip['LOW']['r']:.4f}) but the DEEP "
                       f"window relaxes to r~ = {dip['DEEP']['r']:.4f} -- "
                       f"a boundary/low-mode effect.")
        else:
            verdict = (f"NOT-REPRODUCED: LOW window at the dip angle "
                       f"shows no sub-Poisson depression "
                       f"(r~ = {dip['LOW']['r']:.4f}).")
        d75 = out["angles"]["75"]
        if d75["DEEP"]:
            deeper = (f"{th_dip:.2f}" if dip["DEEP"]["r"]
                      < d75["DEEP"]["r"] else "75.00 (5pi/12)")
            md.append(f"\n- angle pinning (DEEP, ungated): deeper at "
                      f"{deeper} ({dip['DEEP']['r']:.4f} vs "
                      f"{d75['DEEP']['r']:.4f})")
    md.append(f"\n**Reading: {verdict}**")
    out["verdict"] = verdict
    out["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {out['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E16C.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e16c.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
