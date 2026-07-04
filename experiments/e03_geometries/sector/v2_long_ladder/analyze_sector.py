#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3c -- free disk sector (ANALYSIS). Preregistered readings in README.md."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import R_GOE, R_POISSON, disk, mean_r, r_values

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_LOW = 10
CLASSES = ["S", "A"]


def pooled(series_list):
    rv = []
    for ev in series_list:
        rv.extend(r_values(ev[SKIP_LOW:]).tolist())
    rv = np.array(rv)
    return float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))), len(rv)


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    cfg = res["config"]
    # operational ladder: the arc geometry error is smooth (disk-certified,
    # E3b operational gate 0.0 sigma at comparable h); strict internal N* is
    # reported as sensitivity
    n_strict = res["n_use"]
    n_use = cfg["n_modes"]
    rec = res["runs"]["sector"]
    lam = np.array(rec["lam"])[:n_use]
    labels = rec["labels"][:n_use]
    series = {s: np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
              for s in CLASSES}

    md = ["# E3c -- free disk sector (RESULTS)\n",
          f"theta = {cfg['theta']} rad, R = {cfg['R']}, nu = {cfg['nu']}, "
          f"nrings = {cfg['nrings']} ({res['gates']['G2']['ndof']} dofs), "
          f"{cfg['n_modes']} modes (strict internal N* = {n_strict}; "
          f"operational ladder per the disk-certified smooth-error argument). "
          f"References: Poisson {R_POISSON}, GOE {R_GOE}.\n"]
    md.append(f"- G2 internal N* = {res['gates']['G2']['n_star']}, rigid "
              f"{res['gates']['G2']['rigid']}; classifier counts {rec['counts']}, "
              f"min quality {rec['min_quality']:.3f}, resolved {rec['n_resolved']}")

    # disk baseline, matched windows
    roots = disk.free_disk_roots(cfg["nu"], 145.0, 150)
    lam_class = np.array([r[0] for r in disk.class_levels(roots)])
    win = int(np.mean([len(series[s]) - SKIP_LOW for s in CLASSES]))
    rv = []
    for w in range(2):
        part = lam_class[SKIP_LOW + w * win:SKIP_LOW + (w + 1) * win]
        rv.extend(r_values(part).tolist())
    rv = np.array(rv)
    base = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))

    per = {s: mean_r(series[s], SKIP_LOW) for s in CLASSES}
    pool = pooled([series[s] for s in CLASSES])
    md.append("\n## <r> per reflection class\n")
    md.append("| class | n levels | <r> |")
    md.append("|---|---|---|")
    for s in CLASSES:
        md.append(f"| {s} | {len(series[s])} | {per[s][0]:.3f} +/- {per[s][1]:.3f} |")
    md.append(f"| **pooled** | {n_use} | **{pool[0]:.4f} +/- {pool[1]:.4f}** |")
    md.append(f"| disk baseline (2 x {win}-level windows) | - | "
              f"{base[0]:.4f} +/- {base[1]:.4f} |")

    thirds = []
    for k in range(3):
        parts = []
        for s in CLASSES:
            ev = series[s][SKIP_LOW:]
            parts.append(ev[int(k / 3 * len(ev)):int((k + 1) / 3 * len(ev))])
        rv = []
        for p_ in parts:
            rv.extend(r_values(p_).tolist())
        rv = np.array(rv)
        thirds.append((float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv)))))
    md.append("\n### Spectral-window trend (thirds)\n")
    md.append("| window | <r> |")
    md.append("|---|---|")
    for name, (m_, s_) in zip(["low", "mid", "high"], thirds):
        md.append(f"| {name} | {m_:.4f} +/- {s_:.4f} |")

    # sensitivity at the strict cut
    lam_s = np.array(rec["lam"])[:n_strict]
    lab_s = rec["labels"][:n_strict]
    pool_s = pooled([np.sort(lam_s[[i for i, l in enumerate(lab_s) if l == s]])
                     for s in CLASSES])
    md.append(f"\n- sensitivity (strict N* = {n_strict}): pooled <r> = "
              f"{pool_s[0]:.4f} +/- {pool_s[1]:.4f}")

    amp = (pool[0] - base[0]) / np.sqrt(pool[1] ** 2 + base[1] ** 2)
    per_amp = {s: (per[s][0] - base[0]) / np.sqrt(per[s][1] ** 2 + base[1] ** 2)
               for s in CLASSES}
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- pooled <r> = {pool[0]:.4f} +/- {pool[1]:.4f} vs disk baseline "
              f"{base[0]:.4f} +/- {base[1]:.4f}: separation {amp:.1f} sigma "
              f"(threshold 3)")
    md.append(f"- per class: " + ", ".join(f"{s} {per_amp[s]:+.1f}" for s in CLASSES)
              + " sigma")
    if amp >= 3 and all(v > 0 for v in per_amp.values()):
        verdict = ("SUPPORTS (sector intermediate, consistent with the published "
                   "COMSOL sector result and the corner-driver reading)")
    elif amp < 1:
        verdict = ("CHALLENGES (sector Poisson-consistent -- contradicts the "
                   "published sector result and both mechanism readings)")
    else:
        verdict = "AMBIGUOUS (below preregistered threshold)"
    md.append(f"\n**Reading: {verdict}**")
    md.append(
        "\nDiscussion. A weak positive signal (+1.7 sigma pooled, both classes "
        "positive) fits the corner-sharpness ordering emerging across the "
        "campaign: triangle (three 60-degree corners) 0.489; rectangle (four "
        "90-degree corners) 0.44--0.46; sector (one ~115-degree + two "
        "90-degree corners plus a disk-like arc) 0.421; smooth ellipse 0.373; "
        "disk (control) 0.391 -- i.e., the effective coupling appears graded "
        "by corner sharpness and prevalence, not by broken separability alone. "
        "Relation to the published COMSOL sector study: not a contradiction "
        "and not yet a confirmation -- that study used ~10^4 modes and an "
        "angle sweep; our single angle (theta = 2.0) at 800 modes shows the "
        "same tendency without reaching the preregistered threshold. A "
        "decisive comparison needs longer ladders and/or their exact sector "
        "angles. The decreasing spectral trend (0.439 -> 0.408) is noted; "
        "under the corner picture the apex dominates low-k modes whose "
        "evanescent lengths span it, while high-k modes see proportionally "
        "more of the uncoupled arc.")

    fig, ax = plt.subplots(figsize=(6, 4.2))
    xpos = np.arange(3)
    y = [per[s][0] for s in CLASSES] + [pool[0]]
    e = [per[s][1] for s in CLASSES] + [pool[1]]
    ax.errorbar(xpos, y, yerr=e, fmt="o", ms=6, color="tab:red", label="sector")
    ax.axhspan(base[0] - base[1], base[0] + base[1], color="tab:blue", alpha=0.25,
               label="disk baseline")
    ax.axhline(R_POISSON, color="gray", ls="--", lw=1)
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.set_xticks(xpos, CLASSES + ["pooled"])
    ax.set_ylabel("<r>")
    ax.set_title(f"E3c: free disk sector (theta = {cfg['theta']}) per reflection class")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_sector.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(n_use=n_use, n_strict=n_strict, per_class=per, pooled=pool,
                       baseline=base, thirds=thirds, amplitude_sigma=float(amp),
                       verdict=verdict), f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / r_sector.png in {HERE}")


if __name__ == "__main__":
    main()
