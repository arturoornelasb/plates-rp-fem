#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17q finalize -- aggregate the four sector runs and apply the
PREREGISTERED E17 verdict block (identical references and thresholds)."""
import json
import os

import numpy as np

from platefem import SECTORS
from platefem import bases

HERE = os.path.dirname(os.path.abspath(__file__))
E9_SINE_2048 = 0.052
FLAT_SINE = 0.174
LADDER = [256, 512, 1024, 2048]


def main():
    secs = {}
    for s in SECTORS:
        with open(os.path.join(HERE, f"sector_{s}.json")) as f:
            secs[s] = json.load(f)

    md = ["# E17q -- Gap A true-operator windows at N = 1024-2048 "
          "(RESULTS; i01 quarter-plate instrument)\n",
          "Same preregistered readings as E17 (serial run killed "
          "externally post-G2; gates G1 = 1836/5600 informational and "
          "G2 = 1200/1200 preserved in *_serial_partial.json). Four "
          "parity sectors solved independently on the validated i01 "
          "instrument.\n",
          "| sector | dofs | elastic | G3 two-mesh N* | resid | wall s |",
          "|---|---|---|---|---|---|"]
    for s in SECTORS:
        d = secs[s]
        md.append(f"| {s} | {d['solve']['ndof']} | "
                  f"{d['solve']['n_elastic']} | {d['gate']['n_star']}/"
                  f"{d['gate']['n_cmp']} | {d['solve']['resid']:.1e} | "
                  f"{d['wall_time_s']:.0f} |")

    md.append("\n| sector | basis | "
              + " | ".join(f"IPR N={N}" for N in LADDER) + " |")
    md.append("|" + "---|" * (len(LADDER) + 2))
    per_rung = {}
    for s in SECTORS:
        for name in ("sine", "beam"):
            rows = secs[s]["gapA"][name]
            iprs = [float(np.exp(r["mlnipr"])) for r in rows]
            for r, v in zip(rows, iprs):
                per_rung.setdefault((name, int(r["N"])), []).append(v)
            md.append(f"| {s} | {name} | "
                      + " | ".join(f"{v:.4f}" for v in iprs)
                      + " |" * (len(LADDER) - len(iprs) + 1))
    md.append("\nGOE baselines (mean ln IPR): "
              + ", ".join(f"N={N}: {np.exp(v):.4f}" for N, v in
                          bases.goe_ipr_baseline(LADDER).items()))

    top = max(N for (nm, N) in per_rung if nm == "sine")
    i_top = float(np.mean(per_rung[("sine", top)]))
    lnN, lnI = [], []
    for N in LADDER:
        if ("sine", N) in per_rung and N >= 512:
            lnN.append(np.log(N))
            lnI.append(np.log(np.mean(per_rung[("sine", N)])))
    d2_true = float(-np.polyfit(lnN, lnI, 1)[0]) if len(lnN) >= 2 else np.nan
    md.append(f"\n- sine sector-mean IPR at top covered rung N = {top}: "
              f"{i_top:.4f}; fitted D2_true (N >= 512): {d2_true:.3f}; "
              f"references: E9 truncated {E9_SINE_2048} (N = 2048), "
              f"flat E4/E14 ~{FLAT_SINE}")
    if top < max(LADDER):
        verdict = (f"COVERAGE-LIMITED: sector coverage stops at N = {top} "
                   f"< 2048; largest covered rung reported")
    elif i_top >= 0.14:
        verdict = ("PROTOCOL ARTIFACT THROUGH THE FULL REGISTERED RANGE: "
                   "true-operator windows stay sparse to N = 2048; "
                   "truncated-ladder RP is an artifact at every registered "
                   "rung")
    elif i_top <= 0.09:
        verdict = (f"SCALING ONSET: genuine delocalization inside the "
                   f"registered range (D2_true = {d2_true:.3f}; "
                   f"RECONCILED-RP if within E9's 0.42-0.50)")
    else:
        verdict = (f"INTERMEDIATE: between references at N = 2048 "
                   f"(D2_true = {d2_true:.3f})")
    md.append(f"\n**Reading: {verdict}**")

    results = dict(sectors=secs, verdict=verdict,
                   top_rung=int(top), i_top=i_top, d2_true=d2_true)
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
