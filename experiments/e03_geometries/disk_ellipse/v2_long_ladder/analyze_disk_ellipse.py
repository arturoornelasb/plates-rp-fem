#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3b -- free disk (control) vs free ellipse (ANALYSIS). README readings."""
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import R_GOE, R_POISSON, SECTORS, disk, mean_r, r_values

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_LOW = 10


def windows_r(seq, win_len, n_win):
    """<r> pooled over n_win disjoint windows of length win_len (after the
    skip-low cut) of one long sequence -- the matched-protocol baseline."""
    seq = np.asarray(seq)[SKIP_LOW:]
    rv = []
    for w in range(n_win):
        part = seq[w * win_len:(w + 1) * win_len]
        if len(part) >= 3:
            rv.extend(r_values(part).tolist())
    rv = np.array(rv)
    return float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))), len(rv)


def main():
    with open(os.path.join(HERE, "results_raw.json")) as f:
        res = json.load(f)
    cfg = res["config"]
    n_use = res["n_use"]
    md = ["# E3b -- free disk (Poisson control) vs free ellipse (RESULTS)\n",
          f"nu = {cfg['nu']}, Argyris refine-{cfg['refine']} (check "
          f"refine-{cfg['refine_check']}), {cfg['n_modes']} modes, N_use = {n_use}, "
          f"ellipse a/b = {cfg['ellipse_ratio']:.6f} (area pi). References: "
          f"Poisson {R_POISSON}, GOE {R_GOE}.\n"]

    md.append("## Gates\n")
    for r in [cfg["refine_check"], cfg["refine"]]:
        g = res["gates"][f"G1_refine{r}"]
        md.append(f"- G1 FEM disk refine {r} ({g['ndof']} dofs): strict N* = "
                  f"{g['n_star']}, max relerr {g['max_relerr']:.2e} "
                  f"(median {g['med_relerr']:.2e}), rigid {g['n_rigid']}")
    md.append(f"- G2 ellipse internal: N* = {res['gates']['G2']['n_star']}, "
              f"rigid {res['gates']['G2']['rigid']}")

    # ---------------- disk statistics: semi-analytic and FEM ----------------
    lam_class = np.array(res["disk_semianalytic"]["lam_class"])
    lam_distinct = np.array(res["disk_semianalytic"]["lam_distinct"])
    mult = np.array(res["disk_semianalytic"]["mult"])
    # paper-style full-control value on the class sequence
    n_paper = min(len(lam_class), 1460 + SKIP_LOW)
    r_ctrl = mean_r(lam_class[:n_paper], SKIP_LOW)
    md.append(f"\n## Disk control (semi-analytic)\n")
    md.append(f"- pooled one-reflection-class <r> over {r_ctrl[2]} ratios: "
              f"**{r_ctrl[0]:.4f} +/- {r_ctrl[1]:.4f}** (paper's executed "
              f"control: 0.386 +/- 0.007 on 1460 levels)")

    # FEM disk: dedupe by reference multiplicities, compare statistics
    md.append("\n## Operational gate G1s: FEM-disk statistics vs semi-analytic\n")
    g1s_pass = {}
    for r in [cfg["refine_check"], cfg["refine"]]:
        lam_fem = np.array(res["runs"][f"disk_fem_r{r}"]["lam"])
        seq_fem, max_split = disk.dedupe_by_reference(lam_fem, lam_distinct, mult)
        n_cmp = min(len(seq_fem), len(lam_class))
        r_fem = mean_r(seq_fem[:n_cmp], SKIP_LOW)
        r_sa = mean_r(lam_class[:n_cmp], SKIP_LOW)
        sep = abs(r_fem[0] - r_sa[0]) / np.sqrt(r_fem[1] ** 2 + r_sa[1] ** 2)
        # v2: the UPPER half of the extended ladder is where the verdict lives,
        # so it gets its own statistics gate (discretization error grows with
        # mode index and is not smooth -- must be shown harmless there).
        half = n_cmp // 2
        rf_hi = mean_r(seq_fem[half:n_cmp], 0)
        rs_hi = mean_r(lam_class[half:n_cmp], 0)
        sep_hi = abs(rf_hi[0] - rs_hi[0]) / np.sqrt(rf_hi[1] ** 2 + rs_hi[1] ** 2)
        g1s_pass[r] = sep < 2.0 and sep_hi < 2.0
        md.append(f"- refine {r}: FEM class <r> = {r_fem[0]:.4f} +/- {r_fem[1]:.4f} "
                  f"vs semi-analytic {r_sa[0]:.4f} +/- {r_sa[1]:.4f} on the same "
                  f"{n_cmp} levels ({sep:.1f} sigma apart; upper half "
                  f"{rf_hi[0]:.4f} vs {rs_hi[0]:.4f}, {sep_hi:.1f} sigma; doublet "
                  f"max splitting {max_split:.1e}): "
                  f"{'PASS' if g1s_pass[r] else 'FAIL'}")

    # ---------------- ellipse sectors ----------------
    # N_use: per the preregistered protocol the OPERATIONAL gate governs --
    # geometry error is smooth and G1s certifies the statistics on the control
    # over the full ladder, so the full n_modes are used when G1s passes; the
    # strict internal N* is reported as a sensitivity cut below.
    n_strict = n_use
    if g1s_pass[cfg["refine"]]:
        n_use = cfg["n_modes"]
    rec = res["runs"]["ellipse"]
    lam_e = np.array(rec["lam"])[:n_use]
    labels = rec["labels"][:n_use]
    series = {s: np.sort(lam_e[[i for i, l in enumerate(labels) if l == s]])
              for s in SECTORS}
    per_sec = {s: mean_r(series[s], SKIP_LOW) for s in SECTORS}
    rv_all = []
    for s in SECTORS:
        rv_all.extend(r_values(series[s][SKIP_LOW:]).tolist())
    rv_all = np.array(rv_all)
    pool = (float(np.mean(rv_all)), float(np.std(rv_all) / np.sqrt(len(rv_all))),
            len(rv_all))

    # matched-window disk baseline
    win_len = int(np.mean([len(series[s]) - SKIP_LOW for s in SECTORS]))
    base = windows_r(lam_class, win_len, 4)

    md.append("\n## Ellipse sectors vs matched disk baseline\n")
    md.append("| sector | n levels | <r> |")
    md.append("|---|---|---|")
    for s in SECTORS:
        md.append(f"| {s} | {len(series[s])} | {per_sec[s][0]:.3f} +/- {per_sec[s][1]:.3f} |")
    md.append(f"| **pooled** | {n_use} | **{pool[0]:.4f} +/- {pool[1]:.4f}** |")
    md.append(f"| disk baseline (4 x {win_len}-level windows) | - "
              f"| {base[0]:.4f} +/- {base[1]:.4f} |")
    md.append(f"\n- classifier: counts {rec['counts']}, min quality "
              f"{rec['min_quality']:.3f}, resolved {rec['n_resolved']}")
    md.append(f"- N_use = {n_use} (operational, G1s-certified); strict internal "
              f"cut {n_strict} as sensitivity: pooled <r> there = ")

    # window trend (thirds of each sector ladder) + strict-cut sensitivity
    def pool_window(nu_, lo, hi):
        rv = []
        for s in SECTORS:
            ev = np.sort(np.array(rec["lam"])[:nu_][
                [i for i, l in enumerate(rec["labels"][:nu_]) if l == s]])[SKIP_LOW:]
            i0, i1 = int(lo * len(ev)), int(hi * len(ev))
            if i1 - i0 >= 3:
                rv.extend(r_values(ev[i0:i1]).tolist())
        rv = np.array(rv)
        return float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv)))
    strict_pool = pool_window(n_strict, 0.0, 1.0)
    md[-1] += f"{strict_pool[0]:.4f} +/- {strict_pool[1]:.4f}"
    thirds = [pool_window(n_use, k / 3.0, (k + 1) / 3.0) for k in range(3)]
    md.append("\n### Spectral-window trend (thirds of each sector ladder)\n")
    md.append("| window | <r> |")
    md.append("|---|---|")
    for name, (m_, s_) in zip(["low", "mid", "high"], thirds):
        md.append(f"| {name} | {m_:.4f} +/- {s_:.4f} |")
    trend = (thirds[2][0] - thirds[0][0]) / np.sqrt(thirds[2][1] ** 2
                                                    + thirds[0][1] ** 2)
    md.append(f"\n- low -> high trend: {thirds[2][0] - thirds[0][0]:+.4f} "
              f"({trend:+.1f} sigma)")

    # ---------------- verdict ----------------
    amp = (pool[0] - base[0]) / np.sqrt(pool[1] ** 2 + base[1] ** 2)
    disk_poisson = abs(r_ctrl[0] - R_POISSON) < 3 * r_ctrl[1]
    per_pos = {s: (per_sec[s][0] - base[0]) / np.sqrt(per_sec[s][1] ** 2 + base[1] ** 2)
               for s in SECTORS}
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- disk control Poisson-consistent (semi-analytic): {disk_poisson} "
              f"({r_ctrl[0]:.4f} vs {R_POISSON})")
    md.append(f"- FEM-disk statistics gate at production refinement: "
              f"{'PASS' if g1s_pass[cfg['refine']] else 'FAIL'}")
    md.append(f"- ellipse pooled <r> = {pool[0]:.4f} +/- {pool[1]:.4f} vs disk "
              f"baseline {base[0]:.4f} +/- {base[1]:.4f}: separation {amp:.1f} sigma "
              f"(threshold 3)")
    md.append(f"- per-sector separations: "
              + ", ".join(f"{s} {per_pos[s]:+.1f}" for s in SECTORS) + " sigma")
    md.append(f"- spectral trend low -> high: {trend:+.1f} sigma; low window "
              f"sub-Poisson: {thirds[0][0] < R_POISSON}")
    if not g1s_pass[cfg["refine"]]:
        verdict = "VOID (geometry error corrupts statistics at this refinement)"
    elif not disk_poisson:
        verdict = "CHALLENGES (disk control not Poisson -- separability side fails)"
    elif amp >= 3 and all(v > 0 for v in per_pos.values()):
        verdict = ("SUPPORTS the geometry dichotomy (disk Poisson, ellipse "
                   "intermediate)")
    elif trend >= 2.0 and thirds[0][0] < R_POISSON:
        verdict = ("AMBIGUOUS -- weak-coupling regime: sub-Poisson at low modes "
                   "with a significant upward spectral trend. This matches the "
                   "paper's preregistered caveat that at weak coupling level "
                   "repulsion 'has not developed and can even reverse sign, so a "
                   "single weak-coupling run must not be read as refutation'. "
                   "The registered decider is a longer ladder (v2).")
    elif amp < 1 and abs(trend) < 2.0:
        verdict = "CHALLENGES (ellipse consistent with the separable control)"
    else:
        verdict = "AMBIGUOUS (separation below preregistered threshold)"
    md.append(f"\n**Reading: {verdict}**")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    xpos = np.arange(len(SECTORS) + 1)
    y = [per_sec[s][0] for s in SECTORS] + [pool[0]]
    e = [per_sec[s][1] for s in SECTORS] + [pool[1]]
    ax.errorbar(xpos, y, yerr=e, fmt="o", ms=6, color="tab:red",
                label="free ellipse")
    ax.axhspan(base[0] - base[1], base[0] + base[1], color="tab:blue", alpha=0.25,
               label="disk baseline (matched windows)")
    ax.axhline(R_POISSON, color="gray", ls="--", lw=1)
    ax.text(xpos[-1] + 0.2, R_POISSON, "Poisson", fontsize=8, color="gray")
    ax.axhline(R_GOE, color="gray", ls=":", lw=1)
    ax.text(xpos[-1] + 0.2, R_GOE, "GOE", fontsize=8, color="gray")
    ax.set_xticks(xpos, SECTORS + ["pooled"])
    ax.set_ylabel("<r>")
    ax.set_title("E3b: free ellipse per parity sector vs free-disk control")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "r_disk_ellipse.png"), dpi=140)

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(n_use=n_use, disk_control=r_ctrl, disk_baseline=base,
                       ellipse_per_sector=per_sec, ellipse_pooled=pool,
                       amplitude_sigma=float(amp),
                       per_sector_sigma={k: float(v) for k, v in per_pos.items()},
                       g1s_pass={str(k): bool(v) for k, v in g1s_pass.items()},
                       verdict=verdict), f, indent=1, default=float)
    print("\n".join(md))
    print(f"\nWrote RESULTS.md / results.json / r_disk_ellipse.png in {HERE}")


if __name__ == "__main__":
    main()
