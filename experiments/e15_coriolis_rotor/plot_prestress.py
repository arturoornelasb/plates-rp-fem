#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15c figure: prestressed vs bare-Coriolis <r> for the chiral and mirror
rotors, with the linear-validity boundary and GOE/GUE references."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R_GOE, R_GUE = 0.5307, 0.5996
STRAIN_CAP = 0.15


def main():
    with open(os.path.join(HERE, "results_prestress.json")) as f:
        res = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, name, title in [(axes[0], "chir", "chiral rotor"),
                            (axes[1], "mirror", "mirror rotor ($\\sigma_v T$-protected)")]:
        dom = res["domains"][name]
        rows = dom["rows"]
        om = np.array([r["Om"] for r in rows])
        pre = np.array([r["r_pre"] for r in rows])
        sem = np.array([r["sem"] for r in rows])
        bare = np.array([r["r_bare"] for r in rows])
        ax.errorbar(om, pre, yerr=sem, color="tab:blue", marker="D", ms=5,
                    lw=1.8, label="prestressed ($-\\Omega^2 M + \\Omega^2 K_g$)")
        ax.plot(om, bare, color="tab:red", marker="o", ms=4, lw=1.2, ls="--",
                label="bare Coriolis (E15)")
        if name == "chir":
            top3 = np.array([r["r_top3"] for r in rows])
            ax.plot(om, top3, color="tab:blue", marker=".", ms=4, lw=1,
                    ls=":", alpha=0.75, label="top-third window (prestressed)")
        om_cap = np.sqrt(STRAIN_CAP / dom["s1"])
        if om_cap < om[-1]:
            ax.axvspan(om_cap, om[-1] + 0.01, color="gray", alpha=0.12, lw=0)
            ax.text(om_cap + 0.004, 0.592, "strain > 15%", fontsize=8,
                    color="gray", rotation=90, va="top")
        ax.axhline(R_GOE, color="gray", ls="--", lw=1)
        ax.axhline(R_GUE, color="gray", ls=":", lw=1)
        ax.text(0.355, R_GOE + 0.002, "GOE", fontsize=8, color="gray",
                ha="right")
        ax.text(0.355, R_GUE + 0.002, "GUE", fontsize=8, color="gray",
                ha="right")
        ax.set_xlim(-0.012, 0.362)
        ax.set_xlabel("$\\Omega_{\\rm nd} = \\Omega R / c_0$")
        ax.set_title(f"E15c: {title}  (s1 = {dom['s1']:.2f})", fontsize=10)
        ax.legend(fontsize=8,
                  loc="lower right" if name == "chir" else "upper left")
    axes[0].set_ylabel("pooled $\\langle r \\rangle$")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "rotor_prestress.png"), dpi=140)
    print("wrote rotor_prestress.png")


if __name__ == "__main__":
    main()
