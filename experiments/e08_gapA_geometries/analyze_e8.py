#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E8 -- unified Gap A + Dq (ANALYSIS). Preregistered readings in README.md."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
QS = ["1.5", "2", "4"]


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    se = float(np.sqrt(np.sum(res ** 2) / max(len(x) - 2, 1)
                       / np.sum((x - x.mean()) ** 2)))
    return float(coef[0]), se


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    goe = res["goe"]
    md = ["# E8 -- unified Gap A across geometries + Dq (RESULTS)\n",
          "H0 basis = the simply supported eigenbasis of the same domain "
          "(exact M-inner-product coefficients). Ladder D_q from "
          "P_q ~ N^{-(q-1) D_q}, window [0.4N, 0.6N), GOE baselines at "
          "identical N. Registered separator: D_1.5 - D_4 = 0.05--0.07 (RP) "
          "vs 0.20--0.28 (PBRM).\n"]
    summary = {}
    for gname, g in res["geoms"].items():
        md.append(f"\n## {gname}\n")
        md.append("| sector | n_free | leakage | mean IPR by N | D2 ladder "
                  "| +/- | D_1.5 | D_4 | D_1.5 - D_4 | captured |")
        md.append("|" + "---|" * 10)
        summary[gname] = {}
        for s, sec in g["sectors"].items():
            rows = sec["ladder"]
            if len(rows) < 3:
                continue
            lnN = np.log([r["N"] for r in rows])
            d = {}
            for q in QS:
                y = np.array([r[f"mlnP{q}"] for r in rows])
                sl, se = linfit(lnN, y)
                qv = float(q)
                d[q] = (-sl / (qv - 1.0), se / abs(qv - 1.0))
            iprs = " ".join(f"{np.exp(r['mlnP2']):.3f}" for r in rows)
            dq_gap = d["1.5"][0] - d["4"][0]
            summary[gname][s] = dict(D2=d["2"][0], D2_se=d["2"][1],
                                     dq_gap=float(dq_gap),
                                     ipr_last=float(np.exp(rows[-1]["mlnP2"])))
            md.append(f"| {s} | {sec['n_free']} | {sec['leakage']:.3f} | "
                      f"{iprs} | {d['2'][0]:.3f} | {d['2'][1]:.3f} | "
                      f"{d['1.5'][0]:.3f} | {d['4'][0]:.3f} | "
                      f"**{dq_gap:.3f}** | "
                      f"{np.mean([r['mean_captured'] for r in rows]):.3f} |")
        goe_iprs = " ".join(f"{np.exp(goe[str(r['N'])]['2']):.4f}"
                            for r in list(g["sectors"].values())[0]["ladder"])
        md.append(f"| GOE | - | - | {goe_iprs} | ~1 | - | - | - | - | - |")

    # ---------------- readings ----------------
    md.append("\n## Verdict (preregistered readings)\n")
    for gname, ss in summary.items():
        d2s = [v["D2"] for v in ss.values()]
        iprs = [v["ipr_last"] for v in ss.values()]
        gaps = [v["dq_gap"] for v in ss.values()]
        if all(i > 0.5 for i in iprs):
            state = "SPARSE/quasi-diagonal (modes ~ few SS modes)"
        elif all(d < 0.2 for d in d2s):
            state = "sparse regime (flat IPR, D2 ~ 0)"
        elif any(0.2 < d < 0.9 for d in d2s):
            flat = np.mean(gaps)
            state = (f"INTERMEDIATE fractality (D2 in (0.2,0.9)); "
                     f"mean D_1.5-D_4 = {flat:.3f} -> "
                     + ("RP-like (flat Dq)" if flat < 0.15 else
                        "PBRM-like (multifractal)"))
        else:
            state = "mixed -- see table"
        md.append(f"- {gname}: " + ", ".join(
            f"{s} D2={v['D2']:.2f}" for s, v in ss.items()) + f" -> {state}")

    # figure
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = {"rectangle": "tab:gray", "triangle": "tab:red",
              "ellipse": "tab:blue", "superellipse10": "tab:orange"}
    for gname, g in res["geoms"].items():
        sec0 = list(g["sectors"].values())[0]
        Ns = [r["N"] for r in sec0["ladder"]]
        allsec = list(g["sectors"].values())
        y = [np.mean([np.exp(s_["ladder"][k]["mlnP2"]) for s_ in allsec
                      if k < len(s_["ladder"])]) for k in range(len(Ns))]
        ax.loglog(Ns, y, "o-", color=colors.get(gname, "k"), label=gname)
    Ns = [int(k) for k in sorted(goe, key=int)]
    ax.loglog(Ns, [np.exp(goe[str(N)]["2"]) for N in Ns], "k--", label="GOE")
    ax.set_xlabel("N (SS-mode representation truncation)")
    ax.set_ylabel("sector-mean window IPR")
    ax.set_title("E8: Gap A in the same-domain SS eigenbasis, all geometries")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "ipr_geometries.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / ipr_geometries.png in {HERE}")


if __name__ == "__main__":
    main()
