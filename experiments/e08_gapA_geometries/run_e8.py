#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E8 -- unified Gap A across geometries + Dq (RUNNER). See README.md."""
import json
import os
import time

import numpy as np

from platefem import (ElementTriArgyris, assemble_plate, boundary_matrix,
                      centered_probe_operators, classify_c3v,
                      classify_parity_resolved, probe_operators, solve_modes,
                      split_rigid, rectangle_basis, superellipse_basis,
                      triangle_basis, triangle_probe_operators, SECTORS)

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    n_modes=803, kappa_ss=1e10, resid=1e-4,
    ladder=[64, 100, 144, 200],
    qs=[1.5, 2.0, 4.0],
    win=(0.4, 0.6),
    n_goe_real=8, seed=12345,
)
RECT = dict(a=1.0, b=1.0 / 1.6189043236, nu=0.33, mesh=(96, 60))
RATIO = 1.6189043236


def pq_ladder(C_sec, ladder, qs, win):
    """Generalized-IPR ladder on sector coefficients C_sec (n_free, n_ss),
    both frequency-ordered. Representation = lowest-N SS modes."""
    out = []
    n_free, n_ss = C_sec.shape
    for N in ladder:
        i0, i1 = int(win[0] * N), int(win[1] * N)
        if i1 > n_free or N > n_ss:
            break
        rows = dict(N=int(N), n_window=i1 - i0)
        caps, pqs = [], {q: [] for q in qs}
        for k in range(i0, i1):
            d = C_sec[k, :N]
            p = float(d @ d)
            if p <= 0:
                continue
            dn = np.abs(d) / np.sqrt(p)
            caps.append(p)
            for q in qs:
                pqs[q].append(float(np.sum(dn ** (2 * q))))
        rows["mean_captured"] = float(np.mean(caps))
        for q in qs:
            rows[f"mlnP{q:g}"] = float(np.mean(np.log(pqs[q])))
        out.append(rows)
    return out


def goe_pq(ladder, qs, win, n_real, seed):
    out = {}
    for N in ladder:
        rng = np.random.default_rng(seed + N)
        g0, g1 = int(win[0] * N), int(win[1] * N)
        acc = {q: [] for q in qs}
        for _ in range(n_real):
            H = rng.standard_normal((N, N))
            _, U = np.linalg.eigh(0.5 * (H + H.T))
            A = np.abs(U[:, g0:g1])
            for q in qs:
                acc[q].extend(np.sum(A ** (2 * q), axis=0).tolist())
        out[int(N)] = {f"{q:g}": float(np.mean(np.log(acc[q]))) for q in qs}
    return out


def certified_pair(K, M, B, cfg, extra_rigid):
    """(free eigenpairs, SS eigenpairs) with vectors, certified."""
    lamF, VF, sF, _ = solve_modes(K.tocsc(), M, cfg["n_modes"] + extra_rigid,
                                  resid_sanity=cfg["resid"], sweeps_max=30)
    lamF, VF, nrF, _ = split_rigid(lamF, VF)
    Kss = (K + cfg["kappa_ss"] * B).tocsc()
    lamS, VS, sS, _ = solve_modes(Kss, M, cfg["n_modes"],
                                  resid_sanity=cfg["resid"], sweeps_max=30)
    return (lamF, VF, sF), (lamS, VS, sS), Kss


def save(results):
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)


def main():
    t00 = time.time()
    cfg = CFG
    results = {"config": dict(cfg), "geoms": {},
               "goe": goe_pq(cfg["ladder"], cfg["qs"], cfg["win"],
                             cfg["n_goe_real"], cfg["seed"])}

    geoms = {}
    m, bas = rectangle_basis(*RECT["mesh"], RECT["a"], RECT["b"],
                             ElementTriArgyris)
    geoms["rectangle"] = (m, bas, RECT["nu"], "z2",
                          lambda b_: probe_operators(b_, RECT["a"], RECT["b"]))
    m, bas, _ = triangle_basis(7, 1.0)
    geoms["triangle"] = (m, bas, 0.33, "c3v",
                         lambda b_: triangle_probe_operators(b_, 1.0))
    a_ax = np.sqrt(RATIO)
    b_ax = 1.0 / a_ax
    m, bas = __import__("platefem").ellipse_basis(6, a_ax, b_ax)
    geoms["ellipse"] = (m, bas, 0.33, "z2",
                        lambda b_: centered_probe_operators(b_, a_ax, b_ax))
    m, bas = superellipse_basis(6, a_ax, b_ax, 10.0)
    geoms["superellipse10"] = (m, bas, 0.33, "z2",
                               lambda b_: centered_probe_operators(b_, a_ax, b_ax))

    for gname, (mesh, basis, nu, sym, probe_fn) in geoms.items():
        t0 = time.time()
        K, M = assemble_plate(mesh, basis, nu)
        B = boundary_matrix(mesh, ElementTriArgyris)
        (lamF, VF, sF), (lamS, VS, sS), Kss = certified_pair(K, M, B, cfg, 3)
        P, Pmx, Pmy = probe_fn(basis)
        if sym == "z2":
            labF, qF, _ = classify_parity_resolved(lamF, VF, P, Pmx, Pmy, K, M)
            labS, qS, _ = classify_parity_resolved(lamS, VS, P, Pmx, Pmy, Kss, M)
            secs = SECTORS
        else:
            labF, qF, _ = classify_c3v(lamF, VF, P, Pmx, Pmy, K, M)
            labS, qS, _ = classify_c3v(lamS, VS, P, Pmx, Pmy, Kss, M)
            secs = ["A1", "A2", "E"]
        MVS = M @ VS
        g = {"solve_resid": [sF["max_resid"], sS["max_resid"]],
             "counts_free": {s: labF.count(s) for s in secs + ["xx"]},
             "counts_ss": {s: labS.count(s) for s in secs + ["xx"]},
             "sectors": {}}
        for s in secs:
            iF = [i for i, l in enumerate(labF) if l == s]
            iS = [i for i, l in enumerate(labS) if l == s]
            C = VF[:, iF].T @ MVS[:, iS]           # exact M-inner products
            # cross-sector leakage control
            leak = 1.0 - float(np.mean(np.sum(C[:len(iF) // 2] ** 2, axis=1)))
            g["sectors"][s] = dict(
                n_free=len(iF), n_ss=len(iS), leakage=leak,
                lam_free=lamF[iF][:400].tolist(),
                ladder=pq_ladder(np.abs(C), cfg["ladder"], cfg["qs"],
                                 cfg["win"]))
        results["geoms"][gname] = g
        save(results)
        print(f"[{gname}] done: resid {sF['max_resid']:.1e}/{sS['max_resid']:.1e}, "
              f"free counts {g['counts_free']} ({time.time()-t0:.1f} s)")

    results["wall_time_s"] = round(time.time() - t00, 1)
    save(results)
    print(f"\n[done] total {results['wall_time_s']} s -> results_raw.json; "
          f"now run analyze_e8.py")


if __name__ == "__main__":
    main()
