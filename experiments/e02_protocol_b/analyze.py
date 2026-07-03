#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E2 -- Protocol B (ANALYSIS): <r> per sector vs kappa, verdict per the paper.

Reads results_raw.json (from run_sweep.py). For each kappa: per-sector level
sequences (only the first N_use modes; lowest 10 per sector dropped), spacing
ratios, <r> +/- sem per sector and pooled; also low/high spectral windows.
Preregistered readings (paper, Prediction 'Boundary-controlled transition'):
  SUPPORTS   -- monotonic transition from ~0.3863 (Poisson) to the RP value.
  CHALLENGES -- abrupt transition, non-monotonic behavior, or intermediate
                statistics already at nearly-simply-supported (kappa >= 1e8).
"""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import R_GOE, R_POISSON, SECTORS, mean_r, r_values

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_LOW = 10


def sector_series(rec, n_use):
    lam = np.array(rec["lam"])[:n_use]
    labels = rec["labels"][:n_use]
    return {s: np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
            for s in SECTORS}


def pooled_r(series, lo_frac=0.0, hi_frac=1.0):
    """r-values pooled over sectors, using the [lo_frac, hi_frac) index window
    of each sector's elastic sequence (after the skip-low cut)."""
    rv = []
    for s in SECTORS:
        ev = series[s][SKIP_LOW:]
        i0, i1 = int(lo_frac * len(ev)), int(hi_frac * len(ev))
        if i1 - i0 >= 3:
            rv.extend(r_values(ev[i0:i1]).tolist())
    rv = np.array(rv)
    return float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))), len(rv)


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    n_use = res["n_use"]
    cfg = res["config"]
    md = ["# E2 -- Protocol B: boundary-controlled transition (RESULTS)\n",
          f"Rectangle a/b = {cfg['a']/cfg['b']:.6f}, nu = {cfg['nu']}, Argyris "
          f"{cfg['mesh'][0]}x{cfg['mesh'][1]}, {cfg['n_modes']} modes computed, "
          f"N_use = {n_use} (accuracy gates), lowest {SKIP_LOW}/sector dropped. "
          f"References: Poisson {R_POISSON}, GOE {R_GOE}.\n"]

    g = res["gates"]
    md.append("## Accuracy gates\n")
    md.append(f"- G1 (kappa=1e10 vs exact SSSS): N* = {g['G1']['n_star']}, "
              f"max relerr {g['G1']['max_relerr']:.2e}")
    md.append(f"- G2 (kappa=0, 96x60 vs 128x80): internal N* = {g['G2']['n_star']}, "
              f"rigid {g['G2']['rigid']}")
    md.append(f"- G3 (classifier vs Ritz sectors, lowest modes): "
              + ", ".join(f"{s} {g['G3']['per_sector'][s]['max_reldiff']:.1e}"
                          for s in SECTORS)
              + f"; ambiguous {g['G3']['n_ambiguous']}")

    kaps, rows = [], []
    for key, rec in sorted(res["sweep"].items(), key=lambda kv: kv[1]["kappa"]):
        series = sector_series(rec, n_use)
        per_sec = {s: mean_r(series[s], SKIP_LOW) for s in SECTORS}
        pool = pooled_r(series)
        low = pooled_r(series, 0.0, 0.5)
        high = pooled_r(series, 0.5, 1.0)
        kaps.append(rec["kappa"])
        rows.append(dict(kappa=rec["kappa"], per_sector=per_sec, pooled=pool,
                         low=low, high=high, counts=rec["counts"],
                         min_quality=rec["min_quality"]))

    md.append("\n## <r> vs kappa (pooled over sectors; windows are index halves "
              "of each sector ladder)\n")
    md.append("| kappa | " + " | ".join(SECTORS)
              + " | pooled | low half | high half | n_r | min quality |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    for row in rows:
        ps = row["per_sector"]
        md.append(f"| {row['kappa']:.0e} | "
                  + " | ".join(f"{ps[s][0]:.3f}+/-{ps[s][1]:.3f}" for s in SECTORS)
                  + f" | **{row['pooled'][0]:.4f} +/- {row['pooled'][1]:.4f}**"
                  f" | {row['low'][0]:.3f} | {row['high'][0]:.3f}"
                  f" | {row['pooled'][2]} | {row['min_quality']:.2f} |")

    # ---------------- verdict ----------------
    r_free = rows[0]["pooled"][0]
    r_ss = rows[-1]["pooled"][0]
    se_free = rows[0]["pooled"][1]
    se_ss = rows[-1]["pooled"][1]
    r_seq = np.array([r["pooled"][0] for r in rows])
    se_seq = np.array([r["pooled"][1] for r in rows])
    # monotone within noise: fit means decreasing kappa->0 ... sequence sorted by
    # kappa ascending; expect r to DEcrease with kappa (free end high, SS end low)
    viol = np.sum(np.diff(r_seq) > 3 * np.sqrt(se_seq[1:] ** 2 + se_seq[:-1] ** 2))
    near_ss_intermediate = rows[-2]["pooled"][0] > R_POISSON + 5 * rows[-2]["pooled"][1] \
        and rows[-2]["kappa"] >= 1e8
    ok_endpoints = (abs(r_ss - R_POISSON) < 4 * se_ss) and \
        (r_free > R_POISSON + 5 * se_free)
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- endpoints: <r>(kappa=0) = {r_free:.4f} +/- {se_free:.4f} "
              f"(intermediate expected), <r>(kappa=1e10) = {r_ss:.4f} +/- {se_ss:.4f} "
              f"(Poisson {R_POISSON} expected): {'OK' if ok_endpoints else 'NOT OK'}")
    md.append(f"- monotonicity (3-sigma violations along the kappa grid): {viol}")
    md.append(f"- intermediate already at kappa >= 1e8 (challenge condition): "
              f"{near_ss_intermediate}")
    supports = ok_endpoints and viol == 0 and not near_ss_intermediate
    md.append(f"\n**Reading: {'SUPPORTS the hypothesis' if supports else 'CHECK -- see flags above'}**")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(7, 4.6))
    kplot = [max(k, 3e1) for k in kaps]          # kappa=0 plotted at left edge
    for s, mk in zip(SECTORS, "osv^"):
        y = [r["per_sector"][s][0] for r in rows]
        e = [r["per_sector"][s][1] for r in rows]
        ax.errorbar(kplot, y, yerr=e, marker=mk, ms=4, lw=1, alpha=0.6, label=s)
    y = [r["pooled"][0] for r in rows]
    e = [r["pooled"][1] for r in rows]
    ax.errorbar(kplot, y, yerr=e, color="k", marker="D", ms=5, lw=2, label="pooled")
    ax.axhline(R_POISSON, color="gray", ls="--", lw=1)
    ax.text(kplot[-1], R_POISSON + 0.002, "Poisson", fontsize=8, color="gray")
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.text(kplot[-1], R_GOE + 0.002, "GOE", fontsize=8, color="gray")
    ax.set_xscale("log")
    ax.set_xlabel("kappa  (leftmost point = free plate, kappa = 0)")
    ax.set_ylabel("<r> per sector")
    ax.set_title("E2: boundary-controlled transition (Winkler spring, all edges)")
    ax.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_vs_kappa.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(n_use=n_use, rows=rows, supports=bool(supports)), f,
                  indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / r_vs_kappa.png in {HERE}")


if __name__ == "__main__":
    main()
