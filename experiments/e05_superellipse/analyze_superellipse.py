#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E5 -- superellipse corner sweep (ANALYSIS). p = 2 point reused from the
E3b v2 ellipse (identical protocol)."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import R_GOE, R_POISSON, SECTORS, disk, r_values

HERE = os.path.dirname(os.path.abspath(__file__))
E3B_V2 = os.path.normpath(os.path.join(
    HERE, "..", "e03_geometries", "disk_ellipse", "v2_long_ladder",
    "results_raw.json"))
SKIP_LOW = 10


def stats_of(lam, labels, n_use):
    lam = np.array(lam)[:n_use]
    labels = labels[:n_use]
    series = {s: np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
              for s in SECTORS}
    rv = []
    for s in SECTORS:
        rv.extend(r_values(series[s][SKIP_LOW:]).tolist())
    rv = np.array(rv)
    pool = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))
    thirds = []
    for k in range(3):
        rw = []
        for s in SECTORS:
            ev = series[s][SKIP_LOW:]
            rw.extend(r_values(ev[int(k / 3 * len(ev)):
                                  int((k + 1) / 3 * len(ev))]).tolist())
        rw = np.array(rw)
        thirds.append((float(np.mean(rw)), float(np.std(rw) / np.sqrt(len(rw)))))
    return pool, thirds


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    cfg = res["config"]
    n_use = cfg["n_modes"]
    with open(E3B_V2) as f:
        e3b = json.load(f)

    rows = [dict(p=2.0, **dict(zip(("pooled", "thirds"),
            stats_of(e3b["runs"]["ellipse"]["lam"],
                     e3b["runs"]["ellipse"]["labels"], n_use))))]
    for key in sorted(res["runs"], key=lambda k: res["runs"][k]["p"]):
        rec = res["runs"][key]
        pool, thirds = stats_of(rec["lam"], rec["labels"], n_use)
        rows.append(dict(p=rec["p"], pooled=pool, thirds=thirds,
                         quality=rec["mesh_quality"]))

    # disk baseline (matched windows, ~390 levels x 4)
    roots = disk.free_disk_roots(cfg["nu"], 145.0, 150)
    lam_class = np.array([r[0] for r in disk.class_levels(roots)])
    win = (n_use // 4) - SKIP_LOW
    rv = []
    for w in range(4):
        rv.extend(r_values(lam_class[SKIP_LOW + w * win:
                                     SKIP_LOW + (w + 1) * win]).tolist())
    rv = np.array(rv)
    base = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))

    md = ["# E5 -- superellipse corner sweep (RESULTS)\n",
          f"|x/a|^p + |y/b|^p = 1, a/b = {cfg['ratio']:.6f}, nu = {cfg['nu']}, "
          f"refine-{cfg['refine']}, {n_use} modes per p (p = 2 reused from "
          f"E3b v2, identical protocol). Disk baseline {base[0]:.4f} +/- "
          f"{base[1]:.4f}. References: Poisson {R_POISSON}, GOE {R_GOE}.\n"]
    if "G2" in res["gates"]:
        md.append(f"- G2 (p = {res['gates']['G2']['p']}, worst corners): "
                  f"internal N* = {res['gates']['G2']['n_star']}")
    md.append("\n## <r>(p)\n")
    md.append("| p | pooled <r> | low third | mid | high third | vs baseline |")
    md.append("|---|---|---|---|---|---|")
    for r in rows:
        sep = (r["pooled"][0] - base[0]) / np.sqrt(r["pooled"][1] ** 2
                                                   + base[1] ** 2)
        md.append(f"| {r['p']:g} | **{r['pooled'][0]:.4f} +/- {r['pooled'][1]:.4f}** | "
                  + " | ".join(f"{t[0]:.4f}" for t in r["thirds"])
                  + f" | {sep:+.1f} sigma |")

    # ---------------- verdict ----------------
    pool_seq = [r["pooled"][0] for r in rows]
    se_seq = [r["pooled"][1] for r in rows]
    sep10 = (pool_seq[-1] - base[0]) / np.sqrt(se_seq[-1] ** 2 + base[1] ** 2)
    diffs = np.diff(pool_seq)
    mono_up = np.all(diffs > -2 * np.sqrt(np.array(se_seq[1:]) ** 2
                                          + np.array(se_seq[:-1]) ** 2))
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- <r>(p) sequence: "
              + ", ".join(f"p={r['p']:g}: {r['pooled'][0]:.4f}" for r in rows))
    md.append(f"- p = 10 vs disk baseline: {sep10:+.1f} sigma")
    md.append(f"- monotone non-decreasing within noise: {mono_up}")
    if sep10 >= 3 and mono_up:
        verdict = ("CURVATURE-GRADED COUPLING: <r> rises with corner sharpness; "
                   "concentrated boundary curvature drives the evanescent "
                   "coupling (mechanism refinement)")
    elif abs(sep10) < 2 and all(abs((r["pooled"][0] - base[0]))
                                < 3 * np.sqrt(r["pooled"][1] ** 2 + base[1] ** 2)
                                for r in rows):
        verdict = ("TRUE-CORNER THRESHOLD: Poisson-consistent through p = 10; "
                   "only genuine curvature singularities generate effective "
                   "coupling at these mode numbers (mechanism refinement)")
    else:
        verdict = "MIXED/PARTIAL -- see table; check per-p gates before reading"
    md.append(f"\n**Reading: {verdict}**")

    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    ps = [r["p"] for r in rows]
    ax.errorbar(ps, pool_seq, yerr=se_seq, fmt="o-", color="tab:red",
                label="pooled <r>(p)")
    hi = [r["thirds"][2][0] for r in rows]
    ax.plot(ps, hi, "s--", color="tab:orange", ms=4, label="high third")
    ax.axhspan(base[0] - base[1], base[0] + base[1], color="tab:blue",
               alpha=0.25, label="disk baseline")
    ax.axhline(R_POISSON, color="gray", ls="--", lw=1)
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.set_xlabel("superellipse exponent p (corner sharpness)")
    ax.set_ylabel("<r>")
    ax.set_title("E5: corner sweep at fixed broken separability")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_vs_p.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(rows=rows, baseline=base, sep10_sigma=float(sep10),
                       verdict=verdict), f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / r_vs_p.png in {HERE}")


if __name__ == "__main__":
    main()
