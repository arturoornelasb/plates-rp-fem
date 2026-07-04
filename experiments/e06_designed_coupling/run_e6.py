#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E6 -- designed contact coupling + inverse design (RUNNER + ANALYSIS).
Single script: one FEM solve, then pure numpy. See README.md (preregistered)."""
import json
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from platefem import (ElementTriArgyris, assemble_plate, classify_parity_resolved,
                      mean_r, probe_operators, r_values, rectangle_basis,
                      solve_modes, split_rigid, R_GOE, R_POISSON)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    mesh=(96, 60), n_modes=403,
    sector="ee", N=100,
    n_quartets_dense=24,
    strip=(0.08, 0.14),                 # banded: quartets inside this x-strip
    lambdas=np.logspace(-2.5, 1.5, 17).tolist(),
    seed=23,
)


def quartet_points(rng, a, b, n, x_range=(0.05, 0.45), y_range=(0.05, 0.45)):
    """Base points in the lower-left quadrant; mirror quartets are implied."""
    xs = rng.uniform(x_range[0] * a, x_range[1] * a, n)
    ys = rng.uniform(y_range[0] * b, y_range[1] * b, n)
    return np.vstack([xs, ys])


def quartet_V(phi_at_pts):
    """V from mirror quartets for same-sector (ee) modes: each base point
    contributes 4 phi_i phi_j (all four images equal for ee)."""
    return 4.0 * phi_at_pts.T @ phi_at_pts


def sweep(H0diag, V, lambdas):
    out = []
    d = np.mean(np.diff(H0diag))
    for lam in lambdas:
        w, U = np.linalg.eigh(np.diag(H0diag) + lam * V)
        rv = r_values(w)
        ipr = np.mean(np.sum(U ** 4, axis=0)[len(w) // 4: 3 * len(w) // 4])
        out.append(dict(lam=float(lam),
                        coupling_ratio=float(lam * np.mean(np.abs(
                            V[np.triu_indices_from(V, 1)])) / d),
                        r_mean=float(np.mean(rv)),
                        r_sem=float(np.std(rv) / np.sqrt(len(rv))),
                        ipr_mid=float(ipr)))
    return out


def main():
    t00 = time.time()
    cfg = CFG
    a, b = cfg["a"], cfg["b"]
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items() if k != "lambdas"}}
    results["config"]["lambdas"] = cfg["lambdas"]

    # ---------- certified modes ----------
    mesh, basis = rectangle_basis(*cfg["mesh"], a, b, ElementTriArgyris)
    K, M = assemble_plate(mesh, basis, cfg["nu"])
    lam, V_, sinfo, _ = solve_modes(K, M, cfg["n_modes"], resid_sanity=1e-4,
                                    sweeps_max=30)
    lam, V_, n_rigid, _ = split_rigid(lam, V_)
    P, Pmx, Pmy = probe_operators(basis, a, b)
    labels, qual, _ = classify_parity_resolved(lam, V_, P, Pmx, Pmy, K, M)
    idx = [i for i, l in enumerate(labels) if l == cfg["sector"]][:cfg["N"]]
    H0 = lam[idx]
    d_mean = float(np.mean(np.diff(H0)))
    print(f"[modes] {len(idx)} {cfg['sector']} modes, solve resid "
          f"{sinfo['max_resid']:.1e} ({time.time()-t00:.1f} s)")

    def shapes_at(pts):
        return (basis.probes(pts) @ V_[:, idx])       # (n_pts, N)

    # normalize mode shapes to unit L2 (they are M-orthonormal; the dof
    # normalization is inherited -- consistent across patterns, which is all
    # the designed-V comparison needs)

    # ---------- patterns ----------
    patterns = {}
    pts_d = quartet_points(rng, a, b, cfg["n_quartets_dense"])
    patterns["dense"] = quartet_V(shapes_at(pts_d))
    pts_b = quartet_points(rng, a, b, cfg["n_quartets_dense"],
                           x_range=cfg["strip"])
    patterns["banded"] = quartet_V(shapes_at(pts_b))
    patterns["central"] = quartet_V(shapes_at(
        np.array([[a / 2], [b / 2]])))    # single central point (rank 1 in ee)

    results["sweeps"] = {}
    for name, Vmat in patterns.items():
        results["sweeps"][name] = sweep(H0, Vmat, cfg["lambdas"])
        rk = int(np.linalg.matrix_rank(Vmat, tol=1e-8 * np.linalg.norm(Vmat)))
        results["sweeps"][name + "_rank"] = rk
        print(f"[{name}] rank {rk}, sweep done")

    # ---------- inverse design ----------
    # target: the closest same-sector pair in the mid-spectrum
    gaps = np.diff(H0)
    j = int(np.argmin(gaps[30:80])) + 30
    i_t, j_t = j, j + 1
    # scan a grid for x* maximizing |phi_i phi_j| / mean|phi_k phi_l| selectivity
    gxs = np.linspace(0.03 * a, 0.47 * a, 60)
    gys = np.linspace(0.03 * b, 0.47 * b, 40)
    GX, GY = np.meshgrid(gxs, gys, indexing="ij")
    shapes_grid = shapes_at(np.vstack([GX.ravel(), GY.ravel()]))
    tgt = np.abs(shapes_grid[:, i_t] * shapes_grid[:, j_t])
    others = np.mean(np.abs(shapes_grid) ** 2, axis=1)
    sel = tgt / (others + 1e-30)
    k_best = int(np.argmax(sel))
    x_star = (float(GX.ravel()[k_best]), float(GY.ravel()[k_best]))
    phi_i = float(shapes_grid[k_best, i_t])
    phi_j = float(shapes_grid[k_best, j_t])
    Vt = quartet_V(shapes_at(np.array([[x_star[0]], [x_star[1]]])))
    # control law: lambda for target splitting Delta = 5 x natural gap
    gap0 = float(H0[j_t] - H0[i_t])
    Delta_target = 5.0 * gap0
    # 2x2 degenerate-perturbation prediction: splitting(lambda) =
    # sqrt(gap0^2 + ... ) with off-diagonal v = lambda*4*phi_i*phi_j and
    # diagonal shifts lambda*4*phi^2
    lam_ctrl = []
    for fac in [0.5, 1.0, 2.0, 4.0]:
        Delta = fac * Delta_target
        v_needed = 0.5 * np.sqrt(max(Delta ** 2 - gap0 ** 2, 0.0))
        lam_c = v_needed / abs(4.0 * phi_i * phi_j)
        w, U = np.linalg.eigh(np.diag(H0) + lam_c * Vt)
        # measured splitting of the targeted pair (track by overlap)
        ov_i = np.argmax(np.abs(U[i_t]))
        ov_j = np.argmax(np.abs(U[j_t]))
        d11 = lam_c * Vt[i_t, i_t]
        d22 = lam_c * Vt[j_t, j_t]
        v12 = lam_c * Vt[i_t, j_t]
        pred = float(np.sqrt((gap0 + d22 - d11) ** 2 + 4 * v12 ** 2))
        meas = float(abs(w[ov_j] - w[ov_i]))
        mix = float(U[i_t, ov_j] ** 2 + U[j_t, ov_j] ** 2)
        lam_ctrl.append(dict(factor=fac, lam=float(lam_c), predicted=pred,
                             measured=meas,
                             rel_err=float(abs(meas - pred) / pred),
                             pair_content_of_mixed_state=mix))
        print(f"[inverse] Delta x{fac:g}: lam = {lam_c:.3e}, predicted "
              f"{pred:.4f}, measured {meas:.4f} "
              f"(rel err {abs(meas-pred)/pred:.1%})")
    results["inverse_design"] = dict(pair=[int(i_t), int(j_t)], gap0=gap0,
                                     x_star=x_star, phi_i=phi_i, phi_j=phi_j,
                                     selectivity=float(sel[k_best]),
                                     rows=lam_ctrl)

    # ---------- markdown + figure ----------
    md = ["# E6 -- designed contact coupling and inverse design (RESULTS)\n",
          f"FFFF rectangle (campaign geometry), {cfg['sector']} sector, "
          f"N = {cfg['N']}, mean spacing {d_mean:.1f}; contact patterns as "
          f"preregistered. References: Poisson {R_POISSON}, GOE {R_GOE}.\n"]
    md.append("## <r> and mid-spectrum IPR vs coupling (per pattern)\n")
    md.append("| lambda*<|V|>/d | dense <r> | dense IPR | banded <r> | "
              "banded IPR | central <r> |")
    md.append("|---|---|---|---|---|---|")
    for k in range(0, len(cfg["lambdas"]), 2):
        rd = results["sweeps"]["dense"][k]
        rb = results["sweeps"]["banded"][k]
        rc = results["sweeps"]["central"][k]
        md.append(f"| {rd['coupling_ratio']:.2e} | {rd['r_mean']:.3f} | "
                  f"{rd['ipr_mid']:.3f} | {rb['r_mean']:.3f} | "
                  f"{rb['ipr_mid']:.3f} | {rc['r_mean']:.3f} |")
    md.append(f"\n- ranks: dense {results['sweeps']['dense_rank']}, banded "
              f"{results['sweeps']['banded_rank']}, central "
              f"{results['sweeps']['central_rank']}")

    rd_max = max(results["sweeps"]["dense"], key=lambda r: r["r_mean"])
    rb_at = min(results["sweeps"]["banded"],
                key=lambda r: abs(r["lam"] - rd_max["lam"]))
    md.append(f"\n## Verdict elements\n")
    md.append(f"- dense pattern reaches <r> = {rd_max['r_mean']:.3f} +/- "
              f"{rd_max['r_sem']:.3f} (GOE {R_GOE}) with mid IPR "
              f"{rd_max['ipr_mid']:.3f}")
    md.append(f"- banded at the same lambda: <r> = {rb_at['r_mean']:.3f} +/- "
              f"{rb_at['r_sem']:.3f}, IPR {rb_at['ipr_mid']:.3f}")
    ok_dichotomy = (rd_max["r_mean"] > 0.48
                    and rb_at["r_mean"] < rd_max["r_mean"] - 3 * np.sqrt(
                        rd_max["r_sem"] ** 2 + rb_at["r_sem"] ** 2))
    inv = results["inverse_design"]["rows"]
    ok_inverse = all(r["rel_err"] < 0.10 for r in inv)
    md.append(f"- designed-structure dichotomy (dense -> GOE-ward, banded "
              f"stays Poisson-side, > 3 sigma apart): "
              f"{'OK' if ok_dichotomy else 'NOT OK'}")
    md.append(f"- inverse design: pair ({inv[0]['lam']:.2e}-scale posts), "
              f"max rel err of predicted vs measured splitting = "
              f"{max(r['rel_err'] for r in inv):.1%} (threshold 10%): "
              f"{'OK' if ok_inverse else 'NOT OK'}")
    verdict = ("SUPPORTS the contact prediction and demonstrates the control "
               "law" if (ok_dichotomy and ok_inverse) else
               "PARTIAL -- see elements above")
    md.append(f"\n**Reading: {verdict}**")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for name, color in [("dense", "tab:red"), ("banded", "tab:blue"),
                        ("central", "tab:gray")]:
        rows = results["sweeps"][name]
        axes[0].errorbar([r["coupling_ratio"] for r in rows],
                         [r["r_mean"] for r in rows],
                         yerr=[r["r_sem"] for r in rows], marker="o", ms=3,
                         color=color, label=name)
        axes[1].semilogx([r["coupling_ratio"] for r in rows],
                         [r["ipr_mid"] for r in rows], "o-", ms=3,
                         color=color, label=name)
    axes[0].set_xscale("log")
    axes[0].axhline(R_POISSON, color="gray", ls="--", lw=1)
    axes[0].axhline(R_GOE, color="gray", ls=":", lw=1)
    axes[0].set_xlabel("lambda <|V_ij|> / mean spacing")
    axes[0].set_ylabel("<r>")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("lambda <|V_ij|> / mean spacing")
    axes[1].set_ylabel("mid-spectrum IPR (H0 basis)")
    axes[1].legend(fontsize=8)
    fig.suptitle("E6: designed contact coupling (ee sector, N = 100)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "designed_coupling.png"), dpi=140)

    results["wall_time_s"] = round(time.time() - t00, 1)
    results["verdict"] = verdict
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-8:]))
    print(f"\n[done] {results['wall_time_s']} s -> RESULTS.md / results.json")


if __name__ == "__main__":
    main()
