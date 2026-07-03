#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E4 -- Gap A on the true operator (ANALYSIS). Preregistered readings in
README.md: D2 ladders per sector and basis vs GOE, Sigma^2(L) positioning."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import SECTORS
from platefem.bases import number_variance

HERE = os.path.dirname(os.path.abspath(__file__))
BASES = ["sine", "beam"]


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    n = len(x)
    s2 = np.sum((y - yhat) ** 2) / max(n - 2, 1)
    se = float(np.sqrt(s2 / np.sum((x - np.mean(x)) ** 2)))
    return float(coef[0]), se


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    cfg = res["config"]
    g = res["gates"]
    goe = {int(k): v for k, v in res["gapA"]["goe"].items()}
    md = ["# E4 -- Gap A on the true operator (RESULTS)\n",
          f"Rectangle a/b = {cfg['a']/cfg['b']:.6f}, nu = {cfg['nu']}, Argyris "
          f"{cfg['mesh'][0]}x{cfg['mesh'][1]}, N_use = {res['n_use']} certified "
          f"modes (G1 N* = {g['G1']['n_star']}, G2 N* = {g['G2']['n_star']}, "
          f"vector residuals <= {g['G2']['solve']['max_resid']:.1e}). Ladder "
          f"window [0.4N, 0.6N); GOE baseline 8 realizations/size.\n"]
    md.append(f"- labels: {res['labels']['counts']}, min quality "
              f"{res['labels']['min_quality']:.3f}\n")

    # ---------------- IPR tables and D2 fits ----------------
    D2 = {}
    for name in BASES:
        md.append(f"\n## Basis: {name}\n")
        md.append("| sector | " + " | ".join(f"N={r['N']}" for r in
                  res["gapA"]["sectors"]["ee"][name]) + " | D2 ladder | +/- se "
                  "| D2 fixed-N (max N) | mean captured p |")
        md.append("|" + "---|" * (len(res["gapA"]["sectors"]["ee"][name]) + 5))
        for s in SECTORS:
            rows = res["gapA"]["sectors"][s][name]
            lnN = np.log([r["N"] for r in rows])
            y = np.array([r["mlnipr"] for r in rows])
            sl, se = linfit(lnN, y)
            Nmax = rows[-1]["N"]
            d2_fixed = 1.0 - (rows[-1]["mlnipr"] - goe[Nmax]) / np.log(Nmax)
            D2.setdefault(name, {})[s] = dict(ladder=-sl, se=se,
                                              fixedN=float(d2_fixed))
            md.append(f"| {s} | "
                      + " | ".join(f"{np.exp(r['mlnipr']):.4f}" for r in rows)
                      + f" | {-sl:.3f} | {se:.3f} | {d2_fixed:.3f} "
                      f"| {np.mean([r['mean_captured'] for r in rows]):.3f} |")
        goe_row = " | ".join(f"{np.exp(goe[r['N']]):.4f}"
                             for r in res["gapA"]["sectors"]["ee"][name])
        md.append(f"| GOE | {goe_row} | ~1 | - | 1.000 (def.) | - |")

    # ---------------- Sigma^2(L) ----------------
    Lv = np.arange(1, 21)
    md.append("\n## Long-range statistics Sigma^2(L) per sector\n")
    md.append("| L | " + " | ".join(SECTORS) + " | Poisson | GOE |")
    md.append("|" + "---|" * (len(SECTORS) + 3))
    sig = {s: number_variance(res["gapA"]["sectors"][s]["lam"], Lv)
           for s in SECTORS}
    goe_sig = (2.0 / np.pi ** 2) * (np.log(2 * np.pi * Lv) + 0.5772 + 1.0
                                    - np.pi ** 2 / 8.0)
    for i, L in enumerate([1, 2, 5, 10, 15, 20]):
        j = L - 1
        md.append(f"| {L} | " + " | ".join(f"{sig[s][j]:.2f}" for s in SECTORS)
                  + f" | {L:.1f} | {goe_sig[j]:.2f} |")

    # ---------------- verdict ----------------
    lads = [D2[n][s]["ladder"] for n in BASES for s in SECTORS]
    fixs = [D2[n][s]["fixedN"] for n in BASES for s in SECTORS]
    frac_sig20 = np.mean([sig[s][19] / 20.0 for s in SECTORS])
    nonergodic = all(0.05 < d < 0.95 for d in fixs)
    extended = all(d > 0.05 for d in lads)
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- ladder D2 range across bases/sectors: "
              f"[{min(lads):.3f}, {max(lads):.3f}]; fixed-N D2 range: "
              f"[{min(fixs):.3f}, {max(fixs):.3f}] (P12 reference 0.76 +/- 0.15; "
              f"T1 Ritz-based ladder 0.20--0.43 in these bases)")
    md.append(f"- Sigma^2(20)/20 (Poisson = 1, GOE ~ 0.1): mean over sectors "
              f"{frac_sig20:.2f}")
    md.append(f"- eigenvectors extended (D2 ladder > 0 in all cells): {extended}")
    md.append(f"- non-ergodic (fixed-N D2 strictly inside (0,1)): {nonergodic}")
    if extended and nonergodic and 0.1 < frac_sig20 < 1.0:
        verdict = ("RP-COMPATIBLE: extended non-ergodic eigenvectors with "
                   "intermediate long-range statistics on the true operator")
    elif not extended:
        verdict = "LOCALIZED/BANDED-LEANING (IPR ~ N-independent)"
    else:
        verdict = "MIXED -- see tables; basis spread and window caveats apply"
    md.append(f"\n**Reading: {verdict}**")
    md.append("\nCaveats: D2 from a representational ladder at these N carries "
              "the T1 basis-dependence systematic (both registered bases "
              "reported); ladders reach N = 324 (paper registers 256--2048; "
              "larger rungs need quarter-plate reduction -- future work).")

    # ---------------- figures ----------------
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, name in zip(axes, BASES):
        for s, mk in zip(SECTORS, "osv^"):
            rows = res["gapA"]["sectors"][s][name]
            ax.loglog([r["N"] for r in rows],
                      [np.exp(r["mlnipr"]) for r in rows], mk + "-", ms=4,
                      label=f"{s} (D2={D2[name][s]['ladder']:.2f})")
        Ns = [r["N"] for r in res["gapA"]["sectors"]["ee"][name]]
        ax.loglog(Ns, [np.exp(goe[N]) for N in Ns], "k--", label="GOE")
        ax.set_xlabel("N (representation truncation)")
        ax.set_ylabel("mean IPR (window)")
        ax.set_title(f"E4 Gap A -- {name} basis")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "ipr_ladders.png"), dpi=140)

    fig2, ax2 = plt.subplots(figsize=(6, 4.2))
    for s, mk in zip(SECTORS, "osv^"):
        ax2.plot(Lv, sig[s], mk + "-", ms=3, label=s)
    ax2.plot(Lv, Lv, "k--", lw=1, label="Poisson")
    ax2.plot(Lv, goe_sig, "k:", lw=1, label="GOE")
    ax2.set_xlabel("L")
    ax2.set_ylabel("Sigma^2(L)")
    ax2.set_title("E4: number variance per sector (true operator)")
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(os.path.join(HERE, "sigma2.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(D2=D2, sigma2={s: sig[s].tolist() for s in SECTORS},
                       verdict=verdict), f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / figures in {HERE}")


if __name__ == "__main__":
    main()
