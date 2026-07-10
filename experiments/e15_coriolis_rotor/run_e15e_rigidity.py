#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15e -- the eigenvector face of the mechanical GOE->GUE crossover:
phase rigidity (preregistered in header, 2026-07-09).

Phase rigidity of a complex eigenvector phi: R = |phi^T phi| / (phi^+ phi)
-- R = 1 for real (orthogonal-class) vectors, and in the GUE limit R has
the known broad distribution with <R> -> 0 as the class saturates
(mesoscopic literature). FROZEN reading: <R> falls MONOTONICALLY with
Omega for the chiral rotor, tracking the spacing crossover (the
eigenvector face of E15), while the sigma_v T-protected mirror stays at
R ~ 1 at every speed. Reported: full P(R) at the endpoint vs the GOE/GUE
reference ensembles at matched N. Uses the E15d SVK tangents (cached) --
so this is the finite-deformation crossover's eigenvector face.

Analysis experiment on frozen caches + one modal solve per domain."""
import json
import os
import time

import numpy as np
from scipy.linalg import eig

from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from run_e15d import centered_star, CHIR, MIRR, MESH, NM, NU

HERE = os.path.dirname(os.path.abspath(__file__))
OMEGAS = [0.0, 0.1, 0.2, 0.3, 0.35, 0.5]
WIN = (0.4, 0.6)


def pencil_vectors(K_m, G0m, Omega):
    """Right eigenvectors of (K_m - om^2 I + i om Omega G0m) phi = 0 via
    the Hermitian-companion route in the K_m eigenbasis."""
    w_e, Q = np.linalg.eigh(K_m)
    Gp = Q.T @ G0m @ Q
    Gp = 0.5 * (Gp - Gp.T)
    N = len(w_e)
    # linearization: z = (om phi, phi); A z = om z with
    # A = [[-i Omega Gp, diag(w_e)], [I, 0]]  (om^2 I + i om Omega Gp' ...)
    A = np.block([[-1j * Omega * Gp, np.diag(w_e)],
                  [np.eye(N), np.zeros((N, N))]])
    vals, vecs = eig(A)
    keep = vals.real > 1e-9
    om = vals[keep].real
    phi = vecs[N:, keep]                       # phi block, K_m eigenbasis
    order = np.argsort(om)
    return om[order], phi[:, order]


def rigidity(phi):
    num = np.abs(np.einsum("ij,ij->j", phi, phi))
    den = np.einsum("ij,ij->j", np.conj(phi), phi).real
    return num / den


def reference_PR(N, n_real, rng, kind):
    out = []
    for _ in range(n_real):
        S = rng.standard_normal((N, N))
        S = (S + S.T) / 2
        if kind == "GUE":
            Aa = rng.standard_normal((N, N))
            H = S + 1j * (Aa - Aa.T) / 2
        else:
            H = S.astype(complex)
        _, V = np.linalg.eig(H) if kind == "GUE" else np.linalg.eigh(H)
        out.extend(rigidity(np.asarray(V, complex)).tolist())
    return np.array(out)


def main():
    t00 = time.time()
    results = {"domains": {}}
    md = ["# E15e -- phase rigidity across the mechanical GOE->GUE "
          "crossover (RESULTS)\n",
          f"SVK tangents from the E15d cache; window {WIN} of the "
          f"positive spectrum; R = |phi^T phi| / (phi^+ phi).\n"]
    for name, harms in [("chir", CHIR), ("mirror", MIRR)]:
        t0 = time.time()
        m, b = centered_star(*MESH, harms)
        K, M, G0 = e2.assemble_elastic(m, b, NU)
        K, M = K.tocsc(), M.tocsc()
        lam, X, info, _ = solve_modes(K, M, NM, resid_sanity=1e-3,
                                      sweeps_max=30)
        if name == "chir":
            Lam0, G0m, Xn = e2.modal_reduce(K, M, G0, X)
        else:
            S = e2.build_symop(b, "mirror_x", tol=1e-6)
            Lam0, G0m, lab, Xn = e2.parity_adapt_reduce(
                K, M, G0, X, S, lam_cap=float(np.max(lam)))
        el = np.abs(Lam0) > 1e-6 * np.abs(Lam0).max()
        o = np.argsort(Lam0[el])
        idx = np.nonzero(el)[0][o]
        G0_e = G0m[np.ix_(idx, idx)]
        print(f"[{name}] basis ready ({time.time()-t0:.0f} s)")
        rows = []
        md.append(f"\n## {name}\n")
        md.append("| Omega_nd | <R> (window) | median R | frac R > 0.9 |")
        md.append("|---|---|---|---|")
        for Om in OMEGAS:
            cpath = os.path.join(HERE, "e15d_cache", f"{name}_{Om:g}.npz")
            KT_m = np.load(cpath)["KT_m"]
            KT_e = KT_m[np.ix_(idx, idx)] if KT_m.shape[0] != len(idx) \
                else KT_m
            om, phi = pencil_vectors(0.5 * (KT_e + KT_e.T), G0_e, Om)
            i0, i1 = int(WIN[0] * len(om)), int(WIN[1] * len(om))
            R = rigidity(phi[:, i0:i1])
            rows.append(dict(Om=Om, mean_R=float(R.mean()),
                             med_R=float(np.median(R)),
                             frac_hi=float(np.mean(R > 0.9)),
                             n=int(len(R))))
            md.append(f"| {Om:g} | {R.mean():.4f} | {np.median(R):.4f} | "
                      f"{np.mean(R > 0.9):.3f} |")
            print(f"  [{name}] Om={Om:g}: <R> = {R.mean():.4f} "
                  f"({time.time()-t0:.0f} s)")
        results["domains"][name] = rows
        with open(os.path.join(HERE, "results_e15e.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)

    rng = np.random.default_rng(3)
    for kind in ("GOE", "GUE"):
        PR = reference_PR(300, 4, rng, kind)
        i0, i1 = int(WIN[0] * len(PR)), int(WIN[1] * len(PR))
        results[f"ref_{kind}"] = dict(mean=float(PR.mean()),
                                      med=float(np.median(PR)))
    md.append(f"\nReferences (N = 300): GOE <R> = "
              f"{results['ref_GOE']['mean']:.4f}, GUE <R> = "
              f"{results['ref_GUE']['mean']:.4f}.")

    ch = results["domains"]["chir"]
    mono = all(ch[i + 1]["mean_R"] <= ch[i]["mean_R"] + 0.02
               for i in range(len(ch) - 1))
    drop = ch[0]["mean_R"] - min(r["mean_R"] for r in ch)
    prot = min(r["mean_R"] for r in results["domains"]["mirror"])
    if mono and drop > 0.3 and prot > 0.9:
        verdict = (f"SUPPORTS: chiral <R> falls monotonically by "
                   f"{drop:.2f} toward the GUE reference while the "
                   f"protected mirror stays at R >= {prot:.2f} -- the "
                   f"eigenvector face of the mechanical GOE->GUE "
                   f"crossover, first measurement.")
    else:
        verdict = (f"MIXED: monotone={mono}, chiral drop {drop:.2f}, "
                   f"mirror min <R> {prot:.2f} -- see table.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E15E.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e15e.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
