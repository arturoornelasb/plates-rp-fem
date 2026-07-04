#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E9 -- Gap A at the registered N = 2048 rung (RUNNER + ANALYSIS).
T1 truncated-operator ladder with exactly-assembled Ritz matrices."""
import json
import os
import time

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from platefem import ritz_exact

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0, b=1.0 / 1.6189043236, nu=0.33,
    nleg=160, ngrid=1200, mbasis=96,
    ladder=[256, 512, 1024, 2048],
    qs=[1.5, 2.0, 4.0],
    win=(0.4, 0.6),
    sectors=["ee", "oo"],
    n_goe_real=6, seed=12345,
)


# ---------- representation families (ported from the repo T1 audit) ----------
def sine_family(m_vals, xg, L):
    x = 0.5 * L * xg
    S = np.array([np.sqrt(2.0 / L) * np.sin(m * np.pi * (x + L / 2) / L)
                  for m in m_vals])
    return S, np.array([m * np.pi / L for m in m_vals])


def beam_roots(n_roots, even=True):
    g = (lambda bb: np.sin(bb) + np.cos(bb) * np.tanh(bb)) if even else \
        (lambda bb: np.sin(bb) - np.cos(bb) * np.tanh(bb))
    roots = []
    k = 0
    while len(roots) < n_roots:
        lo = k * np.pi + (np.pi / 2 if even else 0.0) + 1e-10
        hi = lo + np.pi / 2 - 2e-10
        try:
            roots.append(brentq(g, lo, hi, xtol=1e-14))
        except ValueError:
            pass
        k += 1
    return np.array(roots)


def cosh_ratio(bb, x):
    ax = np.abs(x)
    return np.exp(bb * (ax - 1.0)) * (1.0 + np.exp(-2.0 * bb * ax)) / \
        (1.0 + np.exp(-2.0 * bb))


def sinh_ratio(bb, x):
    s = np.sign(x)
    ax = np.abs(x)
    return s * np.exp(bb * (ax - 1.0)) * (1.0 - np.exp(-2.0 * bb * ax)) / \
        (1.0 - np.exp(-2.0 * bb))


def beam_family(n_funcs, xg, L, even=True):
    Ws, ks = [], []
    if even:
        Ws.append(np.ones_like(xg)); ks.append(0.0)
        for bb in beam_roots(n_funcs - 1, True):
            Ws.append(np.cos(bb * xg) + np.cos(bb) * cosh_ratio(bb, xg))
            ks.append(2.0 * bb / L)
    else:
        Ws.append(xg.copy()); ks.append(np.pi / L)
        for bb in beam_roots(n_funcs - 1, False):
            Ws.append(np.sin(bb * xg) + np.sin(bb) * sinh_ratio(bb, xg))
            ks.append(2.0 * bb / L)
    return np.array(Ws), np.array(ks)


def loewdin(Fam, wg, L):
    G = (L / 2.0) * (Fam * wg) @ Fam.T
    evg, Ug = np.linalg.eigh(0.5 * (G + G.T))
    return ((Ug / np.sqrt(np.clip(evg, 1e-14, None))) @ Ug.T) @ Fam


def pq_of(dn, qs):
    return {f"{q:g}": float(np.sum(np.abs(dn) ** (2 * q))) for q in qs}


def main():
    t00 = time.time()
    cfg = CFG
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    nleg = cfg["nleg"]
    results = {"config": {k: (list(v) if isinstance(v, (tuple, list)) else v)
                          for k, v in cfg.items()}, "sectors": {}}

    t0 = time.time()
    _, A1, A2, B = ritz_exact.exact_1d_matrices(nleg)
    print(f"[exact] 1D matrices NLEG={nleg} ({time.time()-t0:.1f} s)")

    # quadrature grid for representation overlaps
    xg, wg = np.polynomial.legendre.leggauss(cfg["ngrid"])
    import numpy.polynomial.legendre as npleg
    V0 = np.zeros((nleg, cfg["ngrid"]))
    for i in range(nleg):
        c = np.zeros(i + 1); c[i] = 1.0
        V0[i] = npleg.legval(xg, c) * np.sqrt((2 * i + 1) / 2.0)

    for sec in cfg["sectors"]:
        t0 = time.time()
        par = 0 if sec[0] == "e" else 1
        idx_x = np.arange(par, nleg, 2)
        idx_y = np.arange(0 if sec[1] == "e" else 1, nleg, 2)
        K = ritz_exact.assemble_K_sector(idx_x, idx_y, A1, A2, B, a, b, nu)
        K = 0.5 * (K + K.T)
        nx, ny = len(idx_x), len(idx_y)
        kx = 2.0 * idx_x / a
        ky = 2.0 * idx_y / b
        keyL = (kx[:, None] ** 2 + ky[None, :] ** 2).ravel()
        ordL = np.argsort(keyL, kind="stable")
        print(f"[{sec}] K {K.shape} assembled ({time.time()-t0:.1f} s)")

        # representation families for this sector's parities
        M = cfg["mbasis"]
        even_x = sec[0] == "e"
        even_y = sec[1] == "e"
        m_vals_x = np.arange(1 if even_x else 2, 2 * M + 2, 2)[:M]
        m_vals_y = np.arange(1 if even_y else 2, 2 * M + 2, 2)[:M]
        Sx, ksx = sine_family(m_vals_x, xg, a)
        Sy, ksy = sine_family(m_vals_y, xg, b)
        Bx, kbx = beam_family(M, xg, a, even_x)
        By, kby = beam_family(M, xg, b, even_y)
        Sx, Sy = loewdin(Sx, wg, a), loewdin(Sy, wg, b)
        Bx, By = loewdin(Bx, wg, a), loewdin(By, wg, b)
        Lx = np.sqrt(2.0 / a) * V0[idx_x]
        Ly = np.sqrt(2.0 / b) * V0[idx_y]
        Tsx = (a / 2.0) * (Sx * wg) @ Lx.T
        Tsy = (b / 2.0) * (Sy * wg) @ Ly.T
        Tbx = (a / 2.0) * (Bx * wg) @ Lx.T
        Tby = (b / 2.0) * (By * wg) @ Ly.T
        ordS = np.argsort((ksx[:, None] ** 2 + ksy[None, :] ** 2).ravel(),
                          kind="stable")[:max(cfg["ladder"])]
        spi, spj = np.unravel_index(ordS, (M, M))
        ordB = np.argsort((kbx[:, None] ** 2 + kby[None, :] ** 2).ravel(),
                          kind="stable")[:max(cfg["ladder"])]
        bpi, bpj = np.unravel_index(ordB, (M, M))

        rows = []
        for N in cfg["ladder"]:
            t1 = time.time()
            sel = ordL[:N]
            ev, U = eigh(K[np.ix_(sel, sel)])
            ref = ev[min(6, N - 1)]
            nz = int(np.sum(ev < 1e-6 * abs(ref)))
            evE, UE = ev[nz:], U[:, nz:]
            i0, i1 = int(cfg["win"][0] * len(evE)), int(cfg["win"][1] * len(evE))
            pp, qq = sel // ny, sel % ny
            acc = {bs: {f"{q:g}": [] for q in cfg["qs"]} for bs in ("sin", "beam")}
            caps = {"sin": [], "beam": []}
            for k in range(i0, i1):
                C = np.zeros((nx, ny)); C[pp, qq] = UE[:, k]
                Ds = (Tsx @ C @ Tsy.T)[spi[:N], spj[:N]]
                p = float(Ds @ Ds); caps["sin"].append(p)
                for q_, v in pq_of(Ds / np.sqrt(p), cfg["qs"]).items():
                    acc["sin"][q_].append(v)
                Db = (Tbx @ C @ Tby.T)[bpi[:N], bpj[:N]]
                p = float(Db @ Db); caps["beam"].append(p)
                for q_, v in pq_of(Db / np.sqrt(p), cfg["qs"]).items():
                    acc["beam"][q_].append(v)
            row = dict(N=int(N), n_rigid=nz, n_window=i1 - i0)
            for bs in ("sin", "beam"):
                row[bs] = {q_: float(np.mean(np.log(vv)))
                           for q_, vv in acc[bs].items()}
                row[bs + "_captured"] = float(np.mean(caps[bs]))
            rows.append(row)
            print(f"[{sec}] N={N}: window {i1-i0}, IPR sin "
                  f"{np.exp(row['sin']['2']):.4f} beam "
                  f"{np.exp(row['beam']['2']):.4f} ({time.time()-t1:.1f} s)")
        results["sectors"][sec] = rows
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)

    # GOE baselines
    goe = {}
    for N in cfg["ladder"]:
        t1 = time.time()
        rng = np.random.default_rng(cfg["seed"] + N)
        g0, g1 = int(cfg["win"][0] * N), int(cfg["win"][1] * N)
        acc = {f"{q:g}": [] for q in cfg["qs"]}
        for _ in range(cfg["n_goe_real"]):
            H = rng.standard_normal((N, N))
            _, Ug = np.linalg.eigh(0.5 * (H + H.T))
            A_ = np.abs(Ug[:, g0:g1])
            for q in cfg["qs"]:
                acc[f"{q:g}"].extend(np.sum(A_ ** (2 * q), axis=0).tolist())
        goe[int(N)] = {q_: float(np.mean(np.log(v))) for q_, v in acc.items()}
        print(f"[GOE] N={N} ({time.time()-t1:.1f} s)")
    results["goe"] = goe

    # ---------- report ----------
    def linfit(x, y):
        A_ = np.vstack([x, np.ones_like(x)]).T
        c, *_ = np.linalg.lstsq(A_, y, rcond=None)
        r = y - A_ @ c
        se = float(np.sqrt(np.sum(r ** 2) / max(len(x) - 2, 1)
                           / np.sum((x - x.mean()) ** 2)))
        return float(c[0]), se

    md = ["# E9 -- Gap A at the registered N = 2048 rung (RESULTS)\n",
          f"T1 truncated-operator ladder, exact assembly, NLEG = {cfg['nleg']} "
          f"(sector dim {(cfg['nleg']//2)**2}), window [0.4N, 0.6N). "
          f"P12 reference D2 = 0.76 +/- 0.15.\n"]
    summary = {}
    for sec, rows in results["sectors"].items():
        md.append(f"\n## sector {sec}\n")
        md.append("| basis | " + " | ".join(f"IPR N={r['N']}" for r in rows)
                  + " | D2 ladder | +/- | D_1.5 - D_4 | captured |")
        md.append("|" + "---|" * (len(rows) + 5))
        lnN = np.log([r["N"] for r in rows])
        for bs in ("sin", "beam"):
            d = {}
            for q in ("1.5", "2", "4"):
                sl, se = linfit(lnN, np.array([r[bs][q] for r in rows]))
                d[q] = (-sl / (float(q) - 1.0), se / abs(float(q) - 1.0))
            gap = d["1.5"][0] - d["4"][0]
            summary[f"{sec}_{bs}"] = dict(D2=d["2"][0], D2_se=d["2"][1],
                                          dq_gap=float(gap))
            md.append(f"| {bs} | "
                      + " | ".join(f"{np.exp(r[bs]['2']):.4f}" for r in rows)
                      + f" | {d['2'][0]:.3f} | {d['2'][1]:.3f} | {gap:.3f} | "
                      f"{np.mean([r[bs + '_captured'] for r in rows]):.3f} |")
        md.append("| GOE | "
                  + " | ".join(f"{np.exp(goe[r['N']]['2']):.5f}" for r in rows)
                  + " | ~1 | - | - | - |")

    d2s = [v["D2"] for v in summary.values()]
    gaps = [v["dq_gap"] for v in summary.values()]
    md.append("\n## Verdict (preregistered readings)\n")
    md.append(f"- ladder D2 range: [{min(d2s):.3f}, {max(d2s):.3f}] "
              f"(P12 0.76 +/- 0.15); D_1.5 - D_4 range: "
              f"[{min(gaps):.3f}, {max(gaps):.3f}] (RP < 0.15 < PBRM)")
    if max(d2s) < 0.2:
        verdict = ("SPARSE REGIME PERSISTS through N = 2048: IPR essentially "
                   "flat -- the banded/sparse reading of E4 extends to the "
                   "full registered ladder for the rectangle")
    elif min(d2s) > 0.2 and max(gaps) < 0.15:
        verdict = "RP PHASE DEVELOPS (fractal, flat Dq) at large N"
    elif min(d2s) > 0.2:
        verdict = "SCALING DEVELOPS but multifractal (PBRM-leaning)"
    else:
        verdict = "MIXED across sectors/bases -- see table"
    md.append(f"\n**Reading: {verdict}**")
    results["summary"] = summary
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md))
    print(f"\n[done] {results['wall_time_s']} s")


if __name__ == "__main__":
    main()
