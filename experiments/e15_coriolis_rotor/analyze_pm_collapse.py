#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 post-hoc analysis (labeled as such): do the three rotor sweeps
(bare Coriolis, linear prestress, SVK finite deformation) collapse onto
the Pandey-Mehta GOE->GUE master curve <r>(alpha) with ONE fitted scale
alpha = c * Omega each? Reading: UNIVERSAL-COLLAPSE if the rms residual
of each dataset is <= 2x its sem (~0.008); report fitted scales."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N, REPS = 1200, 3
# dense in the PHYSICAL crossover variable alpha*sqrt(N) in [0, ~3.5]
ALPHAS = (np.array([0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.7,
                    2.2, 3.0, 4.5]) / np.sqrt(1200)).tolist()


def rtilde(ev):
    s = np.diff(np.sort(ev))
    s = s[s > 0]
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    n = len(r)
    return float(np.mean(r[n // 10: 9 * n // 10]))


def master_curve(rng):
    out = []
    for a in ALPHAS:
        vals = []
        for _ in range(REPS):
            S = rng.standard_normal((N, N))
            S = (S + S.T) / np.sqrt(2)
            A = rng.standard_normal((N, N))
            H = (S + 1j * a * (A - A.T) / np.sqrt(2)) / np.sqrt(1 + a * a)
            vals.append(rtilde(np.linalg.eigvalsh(H)))
        out.append(float(np.mean(vals)))
    return np.array(ALPHAS), np.array(out)


def fit_scale(oms, rs, al, rm):
    """One-parameter fit alpha = c*Omega on the INCREMENT r - r(0)
    (removes the small endpoint offset between the rotor's window and
    the synthetic ensemble; c remains the only fitted parameter)."""
    dr = rs - rs[0]
    drm = rm - rm[0]
    cs = np.geomspace(0.005, 2.0, 240)
    best = (None, 1e9)
    for c in cs:
        pred = np.interp(np.clip(c * oms, al[0], al[-1]), al, drm)
        rms = float(np.sqrt(np.mean((pred - dr) ** 2)))
        if rms < best[1]:
            best = (float(c), rms)
    return best


def main():
    rng = np.random.default_rng(9)
    al, rm = master_curve(rng)
    print("[master]", dict(zip(al.tolist(), np.round(rm, 4).tolist())))

    with open(os.path.join(HERE, "results_prestress.json")) as f:
        pc = json.load(f)["domains"]["chir"]["rows"]
    with open(os.path.join(HERE, "results_e15d.json")) as f:
        sd = json.load(f)["domains"]["chir"]["rows"]
    data = {
        "bare (E15c)": (np.array([r["Om"] for r in pc]),
                        np.array([r["r_bare"] for r in pc])),
        "linear prestress (E15c)": (np.array([r["Om"] for r in pc]),
                                    np.array([r["r_pre"] for r in pc])),
        "SVK (E15d)": (np.array([r["Om"] for r in sd]),
                       np.array([r["r"] for r in sd])),
    }
    md = ["# E15 -- Pandey-Mehta collapse of the rotor sweeps (post-hoc "
          "analysis)\n",
          f"Master curve: PM interpolation at N = {N} ({REPS} reps), "
          f"mid-window r-tilde. One fitted scale alpha = c Omega per "
          f"dataset.\n",
          "| dataset | c | rms | verdict (<= 0.016) |", "|---|---|---|---|"]
    results = {"master": dict(alpha=al.tolist(), r=rm.tolist()),
               "fits": {}}
    all_ok = True
    for name, (oms, rs) in data.items():
        c, rms = fit_scale(oms, rs, al, rm)
        ok = rms <= 0.016
        all_ok &= ok
        results["fits"][name] = dict(c=c, rms=rms, ok=bool(ok))
        md.append(f"| {name} | {c:.3f} | {rms:.4f} | "
                  f"{'collapse' if ok else 'DEVIATES'} |")
        print(f"[{name}] c = {c:.3f}, rms = {rms:.4f}")
    verdict = ("UNIVERSAL-COLLAPSE: all three sweeps sit on the "
               "Pandey-Mehta one-parameter family within twice their "
               "sem -- the mechanical crossover follows the universal "
               "GOE->GUE law, not merely its endpoints."
               if all_ok else
               "PARTIAL: see rms column -- at least one sweep deviates "
               "beyond 2x sem from the one-parameter family.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    with open(os.path.join(HERE, "RESULTS_PM_COLLAPSE.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_pm_collapse.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md[-4:]))


if __name__ == "__main__":
    main()
