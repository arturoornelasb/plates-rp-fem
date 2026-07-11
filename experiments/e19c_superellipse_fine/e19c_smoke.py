#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19c SMOKE -- validate the cross-mesh projection machinery on cached
E19b data before the big solves. Project a level-5 SS basis (fresh
small solve, 600 modes) onto the level-6 mesh (collar pullback +
per-sector nested Cholesky) and compute the rung-128 window IPR of the
E19b level-6 free modes in it; compare against the NATIVE level-6 SS
basis (E19b cache) through the SAME code path. PASS gates: no NaN;
Gram cond < 1e3 per sector; |delta mean-ln-IPR| <= 0.05 per covered
sector (level-5 vs native functions approximate the same continuum
family; O(h^4) eigenfunction error). PASS -> launch the production
solves; FAIL -> fix before spending."""
import os
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular

from platefem import solve_modes
from platefem.stats import classify_parity_resolved, centered_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e19c_common import (A_AX, B_AX, KFEM, NU, build_mesh, collar_pullback)

HERE = os.path.dirname(os.path.abspath(__file__))
E19B = os.path.join(HERE, "..", "e19b_superellipse_icirc")
SECT = ["ee", "eo", "oe", "oo"]
WIN = (0.4, 0.6)
N_RUNG = 128


def rung_mln(coeff, idx_f_count, N):
    i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
    if i1 > idx_f_count or N > coeff.shape[0]:
        return None
    out = []
    for k in range(i0, i1):
        d = coeff[:N, k]
        p2 = float(d @ d)
        if p2 <= 0:
            continue
        dn2 = (d / np.sqrt(p2)) ** 2
        out.append(float(np.log(np.sum(dn2 ** 2))))
    return float(np.mean(out))


def sector_coeff(U_or_V, M6, W6, idx_s, idx_f):
    U = U_or_V[:, idx_s]
    MU = M6 @ U
    G = U.T @ MU
    G = 0.5 * (G + G.T)
    L = cholesky(G, lower=True)
    cond = float((np.max(np.diag(L)) / np.min(np.diag(L))) ** 2)
    C = U.T @ W6[:, idx_f]
    return solve_triangular(L, C, lower=True), cond


def main():
    t00 = time.time()
    # small level-5 ss solve (cached, stop-proof)
    p5 = os.path.join(HERE, "smoke_ss5.npz")
    mesh5 = build_mesh(5)
    K5, M5, space5 = assemble_c0ip(mesh5, k=KFEM, nu=NU)
    K5 = 0.5 * (K5 + K5.T)
    D5 = boundary_dofs(space5)
    I5 = np.setdiff1d(np.arange(space5.N), D5)
    if os.path.exists(p5):
        z = np.load(p5)
        lam5, V5 = z["lam"], z["V"]
    else:
        lam5, V5, sinfo, _ = solve_modes(
            K5[I5][:, I5].tocsc(), M5[I5][:, I5].tocsc(), 600,
            resid_sanity=1e-2, sweeps_max=40)
        np.savez(p5, lam=lam5, V=V5)
        print(f"[ss5] 600 modes, resid {sinfo['max_resid']:.1e} "
              f"({time.time()-t00:.0f} s)", flush=True)

    # level-6 target space + E19b caches
    mesh6 = build_mesh(6)
    K6, M6, space6 = assemble_c0ip(mesh6, k=KFEM, nu=NU)
    K6 = 0.5 * (K6 + K6.T)
    K6, M6 = K6.tocsc(), M6.tocsc()
    zf = np.load(os.path.join(E19B, "eig_free_prod.npz"))
    lam_f, V_f = zf["lam"], zf["V"]
    zs = np.load(os.path.join(E19B, "eig_ss_prod.npz"))
    lam_s6, V_s6 = zs["lam"], zs["V"]
    D6 = boundary_dofs(space6)
    I6 = np.setdiff1d(np.arange(space6.N), D6)
    print(f"[spaces] N5 = {space5.N}, N6 = {space6.N} "
          f"({time.time()-t00:.0f} s)", flush=True)

    P6, Pmx6, Pmy6 = centered_probe_operators(space6, A_AX, B_AX, 2000)
    lab_f, _, _ = classify_parity_resolved(lam_f, V_f, P6, Pmx6, Pmy6,
                                           K6, M6)
    lab_s6, _, _ = classify_parity_resolved(
        lam_s6, V_s6, P6[:, I6], Pmx6[:, I6], Pmy6[:, I6],
        K6[I6][:, I6].tocsc(), M6[I6][:, I6].tocsc())
    P5, Pmx5, Pmy5 = centered_probe_operators(space5, A_AX, B_AX, 2000)
    lab_s5, _, _ = classify_parity_resolved(
        lam5, V5, P5[:, I5], Pmx5[:, I5], Pmy5[:, I5],
        K5[I5][:, I5].tocsc(), M5[I5][:, I5].tocsc())
    print(f"[labels] done ({time.time()-t00:.0f} s)", flush=True)

    # cross-mesh operator 5 -> 6 (the machinery under test)
    t1 = time.time()
    pts6 = collar_pullback(np.array(space6.doflocs))
    P56 = space5.probes(pts6)
    t_probes = time.time() - t1
    Vfull5 = np.zeros((space5.N, V5.shape[1]))
    Vfull5[I5] = V5
    U6 = P56 @ Vfull5
    nan_frac = float(np.mean(~np.isfinite(U6)))
    W6 = M6 @ V_f
    Vfull6 = np.zeros((space6.N, V_s6.shape[1]))
    Vfull6[I6] = V_s6
    print(f"[projection] probes {t_probes:.0f} s for {space6.N} pts, "
          f"nan frac {nan_frac:.2e} ({time.time()-t00:.0f} s)",
          flush=True)

    lines = [f"probes wall for {space6.N} pts: {t_probes:.0f} s; "
             f"nan frac {nan_frac:.2e}"]
    ok = nan_frac == 0.0
    for s in SECT:
        idx_f = [i for i, l in enumerate(lab_f) if l == s]
        i5 = [i for i, l in enumerate(lab_s5) if l == s]
        i6 = [i for i, l in enumerate(lab_s6) if l == s]
        if min(len(i5), len(i6)) < N_RUNG or len(idx_f) < 0.6 * N_RUNG:
            lines.append(f"{s}: skipped (coverage {len(i5)}/{len(i6)}/"
                         f"{len(idx_f)})")
            continue
        cx, cond_x = sector_coeff(U6, M6, W6, i5, idx_f)
        cn, cond_n = sector_coeff(Vfull6, M6, W6, i6, idx_f)
        mx = rung_mln(cx, len(idx_f), N_RUNG)
        mn = rung_mln(cn, len(idx_f), N_RUNG)
        d = abs(mx - mn)
        ok = ok and cond_x < 1e3 and d <= 0.05
        lines.append(f"{s}: cross mln {mx:+.4f} vs native {mn:+.4f} "
                     f"(|d| = {d:.4f}); gram cond cross {cond_x:.1e} / "
                     f"native {cond_n:.1e}")
    verdict = "SMOKE PASS" if ok else "SMOKE FAIL"
    lines.append(verdict)
    with open(os.path.join(HERE, "SMOKE.md"), "w", encoding="utf-8") as f:
        f.write("# E19c smoke -- cross-mesh projection validation\n\n"
                + "\n".join(f"- {l}" for l in lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
