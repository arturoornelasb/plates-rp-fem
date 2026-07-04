#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E6 -- designed contact coupling + inverse design (RUNNER + ANALYSIS, v3).

Two identical free plates coupled by posts (paper's contact experiment).
EXACT structure exploited: with V = [[C, -C], [-C, C]], the symmetric
combinations (A+B)/sqrt2 lie in ker(V) (levels unshifted at H0), and the
antisymmetric ones feel H0 + 2*lambda*C. The physical doubled spectrum is
the union of both blocks; only the N x N antisymmetric block needs
diagonalizing.

v3 upgrades over v2 (git history): all four parity sectors pooled (mirror
quartets couple only within a sector: the quartet image signs square out for
same-sector pairs and cancel across sectors), and the inverse-design control
law includes the SECOND-ORDER repulsion term -- still closed-form from the
nodal maps:

    Delta(lambda) = 2 lam q phi_i^2
                  + sum_{j != i} (2 lam q phi_i phi_j)^2 / (E_i - E_j).

Preregistered thresholds unchanged: dense -> GOE-ward vs banded Poisson-side
at > 3 sigma; inverse design within 10% (now for the 2nd-order law; the
1st-order-only error is reported alongside to document why it is needed).
"""
import json
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import (ElementTriArgyris, assemble_plate, classify_parity_resolved,
                      probe_operators, r_values, rectangle_basis, solve_modes,
                      split_rigid, R_GOE, R_POISSON, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(96, 60), n_modes=403,
    n_per_sector=95,
    n_quartets=120,
    strip=(0.08, 0.14),
    c_ratios=np.logspace(-2.0, 2.0, 17).tolist(),
    seed=23,
)
Q = 4.0


def main():
    t00 = time.time()
    cfg = CFG
    a, b = cfg["a"], cfg["b"]
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}}

    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V_, sinfo, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-4,
                                    sweeps_max=30)
    lam, V_, n_rigid, _ = split_rigid(lam, V_)
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, _ = classify_parity_resolved(lam, V_, P, Pmx, Pmy, K, M)
    sec_idx = {s: [i for i, l in enumerate(labels) if l == s][:cfg["n_per_sector"]]
               for s in SECTORS}
    print(f"[modes] per-sector counts "
          f"{ {s: len(v) for s, v in sec_idx.items()} }, resid "
          f"{sinfo['max_resid']:.1e} ({time.time()-t00:.1f} s)")

    def shapes_at(pts, idx):
        return basis.probes(pts) @ V_[:, idx]

    def quartet_pts(n, x_range=(0.05, 0.45)):
        return np.vstack([rng.uniform(x_range[0] * a, x_range[1] * a, n),
                          rng.uniform(0.05 * b, 0.45 * b, n)])

    # ---------- pattern sweeps, pooled over sectors ----------
    def sweep_pattern(pts):
        rows = []
        secdata = {}
        for s in SECTORS:
            idx = sec_idx[s]
            S = shapes_at(pts, idx)
            C = Q * S.T @ S
            H0s = lam[idx]
            secdata[s] = (H0s, C,
                          float(np.mean(np.diff(H0s))),
                          float(np.mean(np.abs(C[np.triu_indices_from(C, 1)]))))
        for c in cfg["c_ratios"]:
            rv, iprs = [], []
            for s in SECTORS:
                H0s, C, d_s, m_s = secdata[s]
                lam_c = c * d_s / m_s
                w_anti, U = np.linalg.eigh(np.diag(H0s) + 2.0 * lam_c * C)
                spectrum = np.sort(np.concatenate([H0s, w_anti]))
                rv.extend(r_values(spectrum[10:]).tolist())
                n2 = len(w_anti)
                iprs.extend(np.sum(U ** 4, axis=0)[n2 // 4: 3 * n2 // 4].tolist())
            rv = np.array(rv)
            rows.append(dict(c=float(c), r_mean=float(np.mean(rv)),
                             r_sem=float(np.std(rv) / np.sqrt(len(rv))),
                             ipr_anti=float(np.mean(iprs))))
        ranks = {s: int(np.sum(np.linalg.svd(secdata[s][1], compute_uv=False)
                               > 0.01 * np.linalg.svd(secdata[s][1],
                                                      compute_uv=False)[0]))
                 for s in SECTORS}
        return rows, ranks

    results["sweeps"] = {}
    for name, pts in [("dense", quartet_pts(cfg["n_quartets"])),
                      ("banded", quartet_pts(cfg["n_quartets"],
                                             x_range=cfg["strip"])),
                      ("central", np.array([[a / 2], [b / 2]]))]:
        rows, ranks = sweep_pattern(pts)
        results["sweeps"][name] = rows
        results["sweeps"][name + "_ranks"] = ranks
        print(f"[{name}] ranks {ranks}, sweep done")

    # ---------- inverse design (2nd-order control law, ee sector) ----------
    idx = sec_idx["ee"]
    H0s = lam[idx]
    i_t = 50
    d_loc = float(0.5 * (H0s[i_t + 1] - H0s[i_t - 1]))
    gxs = np.linspace(0.03 * a, 0.47 * a, 60)
    gys = np.linspace(0.03 * b, 0.47 * b, 40)
    GX, GY = np.meshgrid(gxs, gys, indexing="ij")
    Sg = shapes_at(np.vstack([GX.ravel(), GY.ravel()]), idx)
    near = [j for j in range(len(idx)) if j != i_t and abs(j - i_t) <= 10]
    sel = Sg[:, i_t] ** 2 / (np.sum(Sg[:, near] ** 2, axis=1) + 1e-30)
    k_best = int(np.argmax(sel))
    x_star = (float(GX.ravel()[k_best]), float(GY.ravel()[k_best]))
    phi = Sg[k_best]                            # phi_j(x*) for all ee modes
    rows = []
    for fac in [0.25, 0.5, 1.0]:
        Delta = fac * d_loc
        lam_c = Delta / (2.0 * Q * phi[i_t] ** 2)      # 1st-order inversion
        C1 = Q * np.outer(phi, phi)
        w_anti, U = np.linalg.eigh(np.diag(H0s) + 2.0 * lam_c * C1)
        k_i = int(np.argmax(U[i_t] ** 2))
        meas = float(w_anti[k_i] - H0s[i_t])           # symmetric partner exact
        pred1 = float(2.0 * lam_c * Q * phi[i_t] ** 2)
        pt2 = sum((2.0 * lam_c * Q * phi[i_t] * phi[j]) ** 2
                  / (H0s[i_t] - H0s[j]) for j in range(len(H0s)) if j != i_t)
        pred2 = float(pred1 + pt2)
        rows.append(dict(
            factor=fac, lam=float(lam_c), measured=meas,
            pred_1st=pred1, err_1st=float(abs(meas - pred1) / abs(meas)),
            pred_2nd=pred2, err_2nd=float(abs(meas - pred2) / abs(meas)),
            content=float(U[i_t, k_i] ** 2)))
        print(f"[inverse] Delta = {fac:g} d_loc: measured {meas:.4f}; 1st-order "
              f"pred {pred1:.4f} (err {rows[-1]['err_1st']:.1%}); 2nd-order "
              f"pred {pred2:.4f} (err {rows[-1]['err_2nd']:.1%}); content "
              f"{rows[-1]['content']:.3f}")
    results["inverse_design"] = dict(mode=i_t, x_star=x_star, d_local=d_loc,
                                     rows=rows)

    # ---------- verdict + report ----------
    md = ["# E6 -- designed contact coupling and inverse design (RESULTS, v3)\n",
          f"Two identical free plates, all four parity sectors pooled "
          f"(~{cfg['n_per_sector']}/sector), doubled spectra = symmetric "
          f"kernel block (unshifted) + antisymmetric block H0 + 2 lambda C. "
          f"References: Poisson {R_POISSON}, GOE {R_GOE}.\n"]
    md.append("| c | dense <r> | dense IPR(anti) | banded <r> | banded IPR | central <r> |")
    md.append("|---|---|---|---|---|---|")
    for k in range(0, len(cfg["c_ratios"]), 2):
        rd = results["sweeps"]["dense"][k]
        rb = results["sweeps"]["banded"][k]
        rc = results["sweeps"]["central"][k]
        md.append(f"| {rd['c']:.2f} | {rd['r_mean']:.3f} | {rd['ipr_anti']:.3f} "
                  f"| {rb['r_mean']:.3f} | {rb['ipr_anti']:.3f} "
                  f"| {rc['r_mean']:.3f} |")
    md.append(f"\n- effective ranks/sector: dense "
              f"{results['sweeps']['dense_ranks']}, banded "
              f"{results['sweeps']['banded_ranks']}, central "
              f"{results['sweeps']['central_ranks']}")

    rd_max = max(results["sweeps"]["dense"], key=lambda r: r["r_mean"])
    rb_at = min(results["sweeps"]["banded"],
                key=lambda r: abs(r["c"] - rd_max["c"]))
    sep = (rd_max["r_mean"] - rb_at["r_mean"]) / np.sqrt(
        rd_max["r_sem"] ** 2 + rb_at["r_sem"] ** 2)
    ok_dich = rd_max["r_mean"] > 0.48 and sep > 3
    inv = results["inverse_design"]["rows"]
    ok_inv = all(r["err_2nd"] < 0.10 for r in inv)
    md.append("\n## Verdict elements\n")
    md.append(f"- dense reaches <r> = {rd_max['r_mean']:.3f} +/- "
              f"{rd_max['r_sem']:.3f} at c = {rd_max['c']:.1f}; banded there "
              f"{rb_at['r_mean']:.3f} +/- {rb_at['r_sem']:.3f} "
              f"(separation {sep:.1f} sigma): {'OK' if ok_dich else 'NOT OK'}")
    md.append(f"- control law (2nd order): max rel err "
              f"{max(r['err_2nd'] for r in inv):.1%} "
              f"(1st order alone: {max(r['err_1st'] for r in inv):.1%}): "
              f"{'OK' if ok_inv else 'NOT OK'}")
    verdict = ("SUPPORTS the contact prediction and demonstrates the "
               "second-order control law" if (ok_dich and ok_inv)
               else "PARTIAL -- see elements above")
    md.append(f"\n**Reading: {verdict}**")
    md.append("\nControl law (P10): to split the doubled level i by Delta, "
              "place the post quartet at the nodal-map point x* maximizing "
              "local selectivity and set lambda = Delta / (2 q phi_i(x*)^2), "
              "then correct with the closed-form second-order term "
              "sum_j (2 lam q phi_i phi_j)^2 / (E_i - E_j) -- all quantities "
              "readable from the mode shapes, none requiring the solver.")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for name, color in [("dense", "tab:red"), ("banded", "tab:blue"),
                        ("central", "tab:gray")]:
        rows_ = results["sweeps"][name]
        axes[0].errorbar([r["c"] for r in rows_], [r["r_mean"] for r in rows_],
                         yerr=[r["r_sem"] for r in rows_], marker="o", ms=3,
                         color=color, label=name)
        axes[1].loglog([r["c"] for r in rows_], [r["ipr_anti"] for r in rows_],
                       "o-", ms=3, color=color, label=name)
    axes[0].set_xscale("log")
    axes[0].axhline(R_POISSON, color="gray", ls="--", lw=1)
    axes[0].axhline(R_GOE, color="gray", ls=":", lw=1)
    axes[0].set_xlabel("c = lambda <|C|> / d (per sector)")
    axes[0].set_ylabel("<r> pooled over sectors")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("c")
    axes[1].set_ylabel("antisym-block mid IPR")
    axes[1].legend(fontsize=8)
    fig.suptitle("E6 v3: designed contact coupling, two coupled plates")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "designed_coupling.png"), dpi=140)

    results["wall_time_s"] = round(time.time() - t00, 1)
    results["verdict"] = verdict
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-9:]))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
