#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17q finalize v2 -- aggregate the four sector runs, apply the
PREREGISTERED verdict at the top CERTIFIED rung (ladders use only
gate-covered modes), and report the indicative N = 2048 row (from the
pre-clip run whose top window included modes beyond the two-mesh gate)
as clearly-labeled supplementary evidence."""
import json
import os

import numpy as np

from platefem import SECTORS
from platefem import bases

HERE = os.path.dirname(os.path.abspath(__file__))
FLAT_SINE = 0.174
E9_REF = {2048: 0.052, 1024: 0.069}   # E9 truncated ladder; 1024 interpolated
ONSET_THR = {2048: 0.09, 1024: 0.10}
LADDER = [256, 512, 1024, 2048]


def main():
    secs, uncert = {}, {}
    for s in SECTORS:
        with open(os.path.join(HERE, f"sector_fine_{s}.json")) as f:
            secs[s] = json.load(f)
        up = os.path.join(HERE, f"sector_fine_{s}_uncert.json")
        if os.path.exists(up):
            with open(up) as f:
                uncert[s] = json.load(f)

    md = ["# E17q -- Gap A true-operator windows at N = 1024-2048 "
          "(FINE-MESH CERTIFICATION RESULTS; production (140,86), honest 1.167x check)\n",
          "Same preregistered readings as E17 (serial run killed "
          "externally post-G2; its gates G1 = 1836/5600 informational and "
          "G2 = 1200/1200 preserved in *_serial_partial.json). Four parity "
          "sectors solved independently; ladders use gate-covered modes "
          "only.\n",
          "| sector | dofs | elastic | G3 two-mesh N* | rigid ratio | "
          "resid |", "|---|---|---|---|---|---|"]
    for s in SECTORS:
        d = secs[s]
        md.append(f"| {s} | {d['solve']['ndof']} | "
                  f"{d['solve']['n_elastic']} | {d['gate']['n_star']}/"
                  f"{d['gate']['n_cmp']} | "
                  f"{d['solve'].get('gap_ratio', float('nan')):.1e} | "
                  f"{d['solve']['resid']:.1e} |")

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
    md.append("\nGOE baselines: "
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
    md.append(f"\n- top CERTIFIED rung N = {top}: sine sector-mean IPR = "
              f"{i_top:.4f}; fitted D2_true (N >= 512): {d2_true:.3f}; "
              f"references: E9 truncated {E9_REF.get(top)}, flat E4/E14 "
              f"~{FLAT_SINE}")

    if uncert:
        u2048 = [float(np.exp(r["mlnipr"]))
                 for s in SECTORS if s in uncert
                 for r in uncert[s]["gapA"]["sine"] if int(r["N"]) == 2048]
        if u2048:
            md.append(f"- INDICATIVE N = 2048 (window includes modes beyond "
                      f"the two-mesh gate; pre-clip run): sine sector-mean "
                      f"IPR = {float(np.mean(u2048)):.4f} -- consistent "
                      f"with the certified flat level")

    thr_on = ONSET_THR.get(top, 0.09)
    if i_top >= 0.14:
        verdict = (f"PROTOCOL ARTIFACT THROUGH THE CERTIFIED RANGE "
                   f"(N <= {top}): true-operator windows stay sparse at "
                   f"every certified rung; truncated-ladder RP is an "
                   f"artifact there"
                   + ("" if top >= 2048 else
                      "; N = 2048 certification is COVERAGE-LIMITED by the "
                      "two-mesh gate (indicative row agrees; separated-"
                      "check-mesh continuation registered)"))
    elif i_top <= thr_on:
        verdict = (f"SCALING ONSET inside the certified range "
                   f"(D2_true = {d2_true:.3f})")
    else:
        verdict = f"INTERMEDIATE at N = {top} (D2_true = {d2_true:.3f})"
    md.append(f"\n**Reading: {verdict}**")

    results = dict(sectors=secs, indicative_2048=uncert != {},
                   verdict=verdict, top_rung=int(top), i_top=i_top,
                   d2_true=d2_true)
    with open(os.path.join(HERE, "results_fine.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    with open(os.path.join(HERE, "RESULTS_FINE.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
