#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3a -- free equilateral triangle (ANALYSIS). Preregistered readings in
README.md: free-triangle per-sector <r> vs the same-protocol SS baseline."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import R_GOE, R_POISSON, dedupe_doublets, mean_r, r_values

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_LOW = 10
TRI_SECTORS = ["A1", "A2", "E"]


def sector_series(rec, n_use):
    lam = np.array(rec["lam"])[:n_use]
    labels = rec["labels"][:n_use]
    out, dedup_info = {}, {}
    for s in ["A1", "A2"]:
        out[s] = np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
    ev, n_pairs, n_odd, max_split = dedupe_doublets(lam, labels, "E")
    out["E"] = ev
    dedup_info = dict(n_pairs=n_pairs, n_odd=n_odd, max_split=max_split)
    return out, dedup_info


def pooled_r(series, lo_frac=0.0, hi_frac=1.0):
    rv = []
    for s in TRI_SECTORS:
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
    g = res["gates"]
    md = ["# E3a -- free equilateral triangle (RESULTS)\n",
          f"Side L = {cfg['L']}, nu = {cfg['nu']}, Argyris refine-{cfg['refine']} "
          f"(C3v-symmetric mesh), {cfg['n_modes']} modes computed, N_use = {n_use}, "
          f"lowest {SKIP_LOW}/sector dropped after E-dedup. References: Poisson "
          f"{R_POISSON}, GOE {R_GOE}.\n"]
    md.append("## Accuracy gates\n")
    md.append(f"- G1 (kappa=1e10 vs exact Lame^2): N* = {g['G1']['n_star']}, "
              f"max relerr {g['G1']['max_relerr']:.2e}")
    md.append(f"- G2 (free, refine {cfg['refine']} vs {cfg['refine_check']}): "
              f"internal N* = {g['G2']['n_star']}, rigid {g['G2']['rigid']}")

    rows = {}
    for tag in ["free", "ss"]:
        rec = res["runs"][tag]
        series, dd = sector_series(rec, n_use)
        per_sec = {s: mean_r(series[s], SKIP_LOW) for s in TRI_SECTORS}
        rows[tag] = dict(per_sector=per_sec, pooled=pooled_r(series),
                         low=pooled_r(series, 0.0, 0.5),
                         high=pooled_r(series, 0.5, 1.0),
                         dedup=dd, counts=rec["counts"],
                         min_quality=rec["min_quality"])

    md.append("\n## Classifier / degeneracy diagnostics\n")
    for tag in ["free", "ss"]:
        rec, row = res["runs"][tag], rows[tag]
        md.append(f"- {tag}: counts {row['counts']}, min quality "
                  f"{row['min_quality']:.3f}, E-doublets {row['dedup']['n_pairs']} "
                  f"pairs + {row['dedup']['n_odd']} unpaired, max pair splitting "
                  f"{row['dedup']['max_split']:.1e}, rigid {rec['n_rigid']}")

    md.append("\n## <r> per sector (E deduplicated)\n")
    md.append("| case | A1 | A2 | E | pooled | low half | high half | n_r |")
    md.append("|---|---|---|---|---|---|---|---|")
    for tag in ["free", "ss"]:
        row = rows[tag]
        ps = row["per_sector"]
        md.append(f"| {tag} | "
                  + " | ".join(f"{ps[s][0]:.3f}+/-{ps[s][1]:.3f}" for s in TRI_SECTORS)
                  + f" | **{row['pooled'][0]:.4f} +/- {row['pooled'][1]:.4f}** "
                  f"| {row['low'][0]:.3f} | {row['high'][0]:.3f} | {row['pooled'][2]} |")

    # ---------------- verdict ----------------
    rf, sf = rows["free"]["pooled"][0], rows["free"]["pooled"][1]
    rs, ss_ = rows["ss"]["pooled"][0], rows["ss"]["pooled"][1]
    amp = (rf - rs) / np.sqrt(sf ** 2 + ss_ ** 2)
    per_sector_sep = {s: (rows["free"]["per_sector"][s][0] - rows["ss"]["per_sector"][s][0])
                      / np.sqrt(rows["free"]["per_sector"][s][1] ** 2
                                + rows["ss"]["per_sector"][s][1] ** 2)
                      for s in TRI_SECTORS}
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- free pooled <r> = {rf:.4f} +/- {sf:.4f}; SS baseline (same "
              f"protocol) = {rs:.4f} +/- {ss_:.4f}")
    md.append(f"- separation free vs SS baseline: {amp:.1f} sigma (preregistered "
              f"threshold: 3)")
    md.append(f"- per-sector separations: "
              + ", ".join(f"{s} {per_sector_sep[s]:.1f} sigma" for s in TRI_SECTORS))
    all_positive = all(v > 0 for v in per_sector_sep.values())
    if amp >= 3 and all_positive:
        verdict = "SUPPORTS the geometry prediction (intermediate statistics " \
                  "in the non-separable free triangle)"
    elif amp < 1:
        verdict = "CHALLENGES (free triangle consistent with the separable baseline)"
    else:
        verdict = "AMBIGUOUS (separation below preregistered threshold; " \
                  "the paper itself flags this geometry as speculative)"
    md.append(f"\n**Reading: {verdict}**")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    xpos = np.arange(len(TRI_SECTORS) + 1)
    for tag, off, color in [("free", -0.12, "tab:red"), ("ss", 0.12, "tab:blue")]:
        row = rows[tag]
        y = [row["per_sector"][s][0] for s in TRI_SECTORS] + [row["pooled"][0]]
        e = [row["per_sector"][s][1] for s in TRI_SECTORS] + [row["pooled"][1]]
        ax.errorbar(xpos + off, y, yerr=e, fmt="o", ms=6, color=color,
                    label="free (kappa=0)" if tag == "free" else "SS (kappa=1e10)")
    ax.axhline(R_POISSON, color="gray", ls="--", lw=1)
    ax.text(xpos[-1] + 0.2, R_POISSON, "Poisson", fontsize=8, color="gray")
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.text(xpos[-1] + 0.2, R_GOE, "GOE", fontsize=8, color="gray")
    ax.set_xticks(xpos, TRI_SECTORS + ["pooled"])
    ax.set_ylabel("<r>")
    ax.set_title("E3a: free vs simply-supported equilateral triangle, per C3v sector")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_triangle.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(n_use=n_use, rows=rows, amplitude_sigma=float(amp),
                       per_sector_sigma={k: float(v) for k, v in per_sector_sep.items()},
                       verdict=verdict), f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / r_triangle.png in {HERE}")


if __name__ == "__main__":
    main()
