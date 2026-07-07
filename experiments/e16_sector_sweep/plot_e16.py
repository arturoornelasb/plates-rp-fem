#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E16 figure: per-angle pooled r-tilde vs sector angle, with the pooled
sweep value, the same-protocol pooled-Poisson baseline, and the
fixed-angle references."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_sweep import rt_ratios

HERE = os.path.dirname(os.path.abspath(__file__))
R_GOE = 0.5307


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    cfg = res["config"]
    st = res["stats"]

    th, rbar, rsem = [], [], []
    for rec in res["sweep"]:
        rs = np.concatenate([rt_ratios(
            [l for l, s in zip(rec["lam"], rec["labels"]) if s == c],
            cfg["drop_low"]) for c in ("S", "A")])
        th.append(rec["theta_deg"])
        rbar.append(rs.mean())
        rsem.append(rs.std() / np.sqrt(len(rs)))
    th, rbar = np.array(th), np.array(rbar)

    nb = 20
    edges = np.linspace(th[0], th[-1], nb + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bmean = [rbar[(th >= a) & (th < b)].mean() for a, b in zip(edges, edges[1:])]
    bsem = [rbar[(th >= a) & (th < b)].std()
            / np.sqrt(max(1, ((th >= a) & (th < b)).sum()))
            for a, b in zip(edges, edges[1:])]

    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(th, rbar, ".", color="tab:blue", ms=3, alpha=0.35,
            label="per-angle (first 70 modes)")
    ax.errorbar(centers, bmean, yerr=bsem, color="tab:blue", marker="D",
                ms=5, lw=1.8, label="binned mean")
    pm, ps = st["pooled_lam"]["pooled"]["mean"], st["pooled_lam"]["pooled"]["sem"]
    ax.axhspan(pm - ps, pm + ps, color="tab:blue", alpha=0.15, lw=0)
    ax.text(31, pm + 0.004, f"pooled sweep {pm:.3f}", fontsize=8,
            color="tab:blue")
    bp, bs = st["poisson_baseline"]["mean"], st["poisson_baseline"]["std"]
    ax.axhline(bp, color="gray", ls="--", lw=1)
    ax.text(148, bp + 0.003, "pooled-Poisson baseline", fontsize=8,
            color="gray", ha="right")
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.text(148, R_GOE + 0.003, "GOE", fontsize=8, color="gray", ha="right")
    for k, r in res["refs"].items():
        ax.plot(float(k), r["r_matched"], marker="s", color="tab:red", ms=7,
                mfc="none", mew=1.6)
    ax.plot([], [], marker="s", color="tab:red", ms=7, mfc="none", mew=1.6,
            ls="none", label="fixed-angle refs (two-mesh gated)")
    ax.set_xlabel("sector angle theta (deg)")
    ax.set_ylabel("pooled r-tilde (modes 11-70)")
    ax.set_title("E16: sector angle sweep, faithful protocol "
                 "(200 angles x 70 modes)")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_vs_theta.png"), dpi=140)
    print("wrote r_vs_theta.png")


if __name__ == "__main__":
    main()
