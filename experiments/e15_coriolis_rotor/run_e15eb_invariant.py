#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15e-b -- the CORRECT protected observable for the sigma_v T rotor
(preregistered in header, 2026-07-10). E15e showed the naive phase
rigidity R = |phi^T phi|/(phi^+ phi) is NOT pinned to 1 on the mirror
rotor: the protecting antiunitary is sigma_v T, not T, so the invariant
is the S-weighted rigidity

    R_S = |phi^T S phi| / (phi^+ phi),   S = parity metric (diag +/-1),

which is eigenvector-phase invariant and equals 1 EXACTLY for any
nondegenerate sigma_v T-symmetric state (phi = e^{i a} S conj(phi)).
Two arms, one prediction each:
  (i)  PROTECTED: mirror rotor, all cached Omega -- median R_S >= 0.999
       at every speed (machine-precision protection at the eigenvector
       level, where the naive R fell);
  (ii) BROKEN (falsification): add a real symmetric mistuning V coupling
       ONLY opposite-parity pairs (RMS element = window spacing d) at
       Omega = 0.3 -- median R_S < 0.99.
FROZEN reading: PROTECTION-EXACT-AND-FALSIFIABLE if both hold;
PROTECTION-EXACT if only (i); NOT-EXACT otherwise. Structure gates
(sanity, must pass): in the parity basis KT is parity-block-diagonal
(opp-parity fraction < 1e-6) and G0 is purely opposite-parity
(same-parity fraction < 1e-6) -- the algebra behind the protection; the
1e-6 floor is the Newton-prestate mirror-symmetry tolerance, not the
FEM's.
Analysis experiment on the frozen E15d caches + one modal solve."""
import json
import os
import time

import numpy as np
from scipy.linalg import eig

from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes
from run_e15d import centered_star, MIRR, MESH, NM, NU

HERE = os.path.dirname(os.path.abspath(__file__))
OMEGAS = [0.0, 0.1, 0.2, 0.3, 0.35, 0.5]
WIN = (0.4, 0.6)


def pencil_vectors_lab(K_m, G0m, Omega):
    """Same validated companion route as E15e, but vectors returned in
    the INPUT (parity-adapted) basis so S = diag(lab) applies."""
    w_e, Q = np.linalg.eigh(K_m)
    Gp = Q.T @ G0m @ Q
    Gp = 0.5 * (Gp - Gp.T)
    N = len(w_e)
    A = np.block([[-1j * Omega * Gp, np.diag(w_e)],
                  [np.eye(N), np.zeros((N, N))]])
    vals, vecs = eig(A)
    keep = vals.real > 1e-9
    om = vals[keep].real
    phi = vecs[N:, keep]
    order = np.argsort(om)
    return om[order], Q @ phi[:, order]


def rig_S(phi, s):
    num = np.abs(np.einsum("ij,i,ij->j", phi, s, phi))
    den = np.einsum("ij,ij->j", np.conj(phi), phi).real
    return num / den


def rig_plain(phi):
    num = np.abs(np.einsum("ij,ij->j", phi, phi))
    den = np.einsum("ij,ij->j", np.conj(phi), phi).real
    return num / den


def main():
    t00 = time.time()
    m, b = centered_star(*MESH, MIRR)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    K, M = K.tocsc(), M.tocsc()
    lam, X, info, _ = solve_modes(K, M, NM, resid_sanity=1e-3,
                                  sweeps_max=30)
    S = e2.build_symop(b, "mirror_x", tol=1e-6)
    Lam0, G0m, lab, Xn = e2.parity_adapt_reduce(K, M, G0, X, S,
                                                lam_cap=float(np.max(lam)))
    el = np.abs(Lam0) > 1e-6 * np.abs(Lam0).max()
    o = np.argsort(Lam0[el])
    idx = np.nonzero(el)[0][o]
    G0_e = G0m[np.ix_(idx, idx)]
    lab_e = np.asarray(lab)[idx].astype(float)
    n = len(idx)
    print(f"[mirror] basis ready, n_e = {n} ({time.time()-t00:.0f} s)",
          flush=True)

    results = {"n_e": int(n)}
    md = ["# E15e-b -- the sigma_v T invariant R_S (RESULTS)\n",
          "R_S = |phi^T S phi| / (phi^+ phi), S = parity metric; window "
          f"{WIN}. Naive R shown for contrast (E15e).\n"]

    # ---- structure gates ----
    opp = np.outer(lab_e, lab_e) < 0
    kt0 = np.load(os.path.join(HERE, "e15d_cache", "mirror_0.3.npz"))["KT_m"]
    KT03 = kt0[np.ix_(idx, idx)] if kt0.shape[0] != n else kt0
    KT03 = 0.5 * (KT03 + KT03.T)
    fK = float(np.linalg.norm(KT03[opp]) / np.linalg.norm(KT03))
    fG = float(np.linalg.norm(G0_e[~opp]) / (np.linalg.norm(G0_e) + 1e-300))
    results["gate_KT_oppfrac"] = fK
    results["gate_G0_samefrac"] = fG
    gates_ok = fK < 1e-6 and fG < 1e-6
    md.append(f"- structure gates: KT opp-parity fraction = {fK:.2e}, "
              f"G0 same-parity fraction = {fG:.2e} -> "
              f"{'PASS' if gates_ok else 'FAIL'}\n")
    md.append("| Omega_nd | median R_S | p5 R_S | frac R_S < 0.999 | "
              "median R (naive) |")
    md.append("|---|---|---|---|---|")

    # ---- arm (i): protected ----
    rows = []
    for Om in OMEGAS:
        z = np.load(os.path.join(HERE, "e15d_cache",
                                 f"mirror_{Om:g}.npz"))["KT_m"]
        KT_e = z[np.ix_(idx, idx)] if z.shape[0] != n else z
        om, phi = pencil_vectors_lab(0.5 * (KT_e + KT_e.T), G0_e, Om)
        i0, i1 = int(WIN[0] * len(om)), int(WIN[1] * len(om))
        RS = rig_S(phi[:, i0:i1], lab_e)
        RN = rig_plain(phi[:, i0:i1])
        rows.append(dict(Om=Om, med_RS=float(np.median(RS)),
                         p5_RS=float(np.percentile(RS, 5)),
                         frac_low=float(np.mean(RS < 0.999)),
                         med_R=float(np.median(RN))))
        md.append(f"| {Om:g} | {rows[-1]['med_RS']:.6f} | "
                  f"{rows[-1]['p5_RS']:.6f} | {rows[-1]['frac_low']:.3f} "
                  f"| {rows[-1]['med_R']:.4f} |")
        print(f"  Om={Om:g}: med R_S = {rows[-1]['med_RS']:.6f}, "
              f"naive R = {rows[-1]['med_R']:.4f} "
              f"({time.time()-t00:.0f} s)", flush=True)
    results["protected"] = rows

    # ---- arm (ii): broken (falsification) ----
    Om_b = 0.3
    w_all = np.linalg.eigvalsh(KT03)
    j0, j1 = int(WIN[0] * n), int(WIN[1] * n)
    d = float(np.median(np.diff(w_all[j0:j1])))
    rng = np.random.default_rng(11)
    B = rng.standard_normal((n, n))
    B = 0.5 * (B + B.T)
    V = np.where(opp, B, 0.0)
    V *= d / np.sqrt(np.mean(V[opp] ** 2))
    om, phi = pencil_vectors_lab(KT03 + V, G0_e, Om_b)
    i0, i1 = int(WIN[0] * len(om)), int(WIN[1] * len(om))
    RSb = rig_S(phi[:, i0:i1], lab_e)
    results["broken"] = dict(Om=Om_b, rms_over_d=1.0,
                             med_RS=float(np.median(RSb)),
                             p5_RS=float(np.percentile(RSb, 5)))
    md.append(f"\n- BROKEN arm (odd-parity mistuning, RMS element = d, "
              f"Omega = {Om_b:g}): median R_S = "
              f"{results['broken']['med_RS']:.4f} "
              f"(p5 = {results['broken']['p5_RS']:.4f})")

    prot_ok = all(r["med_RS"] >= 0.999 for r in rows)
    broke_ok = results["broken"]["med_RS"] < 0.99
    if prot_ok and broke_ok and gates_ok:
        verdict = ("PROTECTION-EXACT-AND-FALSIFIABLE: the sigma_v T "
                   "invariant R_S stays at 1 (median >= 0.999) at every "
                   "speed while the naive R falls, and an odd-parity "
                   "mistuning collapses it to "
                   f"{results['broken']['med_RS']:.3f} -- the correct "
                   "eigenvector-level protection observable, measured "
                   "and falsified in one experiment.")
    elif prot_ok and gates_ok:
        verdict = (f"PROTECTION-EXACT (weak falsification arm: broken "
                   f"median R_S = {results['broken']['med_RS']:.4f}).")
    else:
        verdict = (f"NOT-EXACT or gate failure: protected min median "
                   f"{min(r['med_RS'] for r in rows):.6f}, gates "
                   f"{'PASS' if gates_ok else 'FAIL'} -- see table.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    results["wall_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall: {results['wall_s']} s.")
    with open(os.path.join(HERE, "RESULTS_E15EB.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_e15eb.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-6:]))


if __name__ == "__main__":
    main()
