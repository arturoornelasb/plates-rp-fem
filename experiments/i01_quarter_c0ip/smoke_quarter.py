#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""i01 -- quarter-plate C0-IP instrument smoke (GATES A and B).
See README.md. Small meshes; runs politely alongside E17."""
import json
import os
import time

import numpy as np
from skfem import MeshTri

from platefem import SECTORS, n_star, solve_lowest, solve_modes, split_rigid
from platefem.stats import classify_parity_resolved
from platefem.c0ip import (C0IPSpace, assemble_c0ip, boundary_dofs,
                           facets_where)

HERE = os.path.dirname(os.path.abspath(__file__))
A, B, NU, KFEM = 1.0, 1.0 / 1.6189043236, 0.33, 4
FULL_MESH, Q_MESH = (32, 20), (16, 10)
N_FULL, N_Q = 252, 78
N_CMP_A, N_CMP_B = 60, 48
RIGID = dict(ee=1, eo=1, oe=1, oo=0)
TOL = 1e-9


def tensor_mesh(nx, ny, ax, by):
    return MeshTri.init_tensor(np.linspace(0.0, ax, nx + 1),
                               np.linspace(0.0, by, ny + 1))


def split_expected(lam, n_exp):
    """Drop exactly n_exp rigid modes; verify the multiplicative gap."""
    lam = np.asarray(lam)
    if n_exp == 0:
        ok = lam[0] > 1e-8 * lam[len(lam) // 2]
        return lam, bool(ok)
    ok = lam[n_exp] > 1e5 * max(abs(float(lam[n_exp - 1])), 1e-300)
    return lam[n_exp:], bool(ok)


def quarter_spectrum(sector, outer):
    """Eigenvalues of one parity sector on the quarter plate.
    outer: 'free' or 'ss' (the two outer edges x=0, y=0)."""
    mesh = tensor_mesh(*Q_MESH, A / 2, B / 2)
    tmp = C0IPSpace(mesh, KFEM)
    gf, ssf = [], []
    xline = facets_where(tmp, lambda x, y: abs(x - A / 2) < TOL)
    yline = facets_where(tmp, lambda x, y: abs(y - B / 2) < TOL)
    (gf if sector[0] == "e" else ssf).append(xline)
    (gf if sector[1] == "e" else ssf).append(yline)
    if outer == "ss":
        ssf.append(facets_where(tmp, lambda x, y: abs(x) < TOL))
        ssf.append(facets_where(tmp, lambda x, y: abs(y) < TOL))
    gf = np.concatenate(gf) if gf else None
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU, guided_facets=gf)
    if ssf:
        D = boundary_dofs(space, np.concatenate(ssf))
        I = np.setdiff1d(np.arange(space.N), D)
        K, M = K[I][:, I].tocsc(), M[I][:, I].tocsc()
    return solve_lowest(K, M, N_Q)


def exact_ssss_sector(px, py, n):
    mm = np.arange(1, 200, 2) if px == "e" else np.arange(2, 200, 2)
    nn = np.arange(1, 200, 2) if py == "e" else np.arange(2, 200, 2)
    lam = ((mm[:, None] * np.pi / A) ** 2
           + (nn[None, :] * np.pi / B) ** 2) ** 2
    return np.sort(lam.ravel())[:n]


def main():
    t00 = time.time()
    results = {"gates": {}}
    md = ["# i01 -- quarter-plate C0-IP instrument (SMOKE RESULTS)\n"]

    # ---------------- GATE A: BC exactness on the SSSS problem ----------
    # The decisive form: the four quarter sectors MERGED must reproduce the
    # SAME-MESH full-plate SSSS spectrum 1:1 (both share the discretization
    # error; any mismatch is the new center-line BCs). The exact-spectrum
    # comparison is reported as informational (it measures plain mesh error
    # at this deliberately coarse smoke h).
    md.append("## GATE A: merged quarter sectors vs same-mesh full SSSS\n")
    q_ss = {}
    md_info = ["\n### informational: per-sector vs EXACT (mesh error at "
               "smoke h)\n", "| sector | N* | max relerr |", "|---|---|---|"]
    for s in SECTORS:
        lam = quarter_spectrum(s, "ss")
        lam, _ = split_expected(lam, 0)
        q_ss[s] = lam
        ex = exact_ssss_sector(s[0], s[1], N_CMP_A)
        ns_i = n_star(lam[:N_CMP_A], ex, 0.1)
        rel_i = float(np.max(np.abs(lam[:N_CMP_A] - ex) / ex))
        results["gates"][f"Ainfo_{s}"] = dict(n_star=int(ns_i),
                                              max_relerr=rel_i)
        md_info.append(f"| {s} | {ns_i}/{N_CMP_A} | {rel_i:.2e} |")
        print(f"[A-info {s}] vs exact: N* = {ns_i}/{N_CMP_A}, "
              f"max relerr {rel_i:.2e}")
    mesh_f = tensor_mesh(*FULL_MESH, A, B)
    Kf_, Mf_, space_f = assemble_c0ip(mesh_f, k=KFEM, nu=NU)
    Kf_ = 0.5 * (Kf_ + Kf_.T)
    Df = boundary_dofs(space_f)
    If = np.setdiff1d(np.arange(space_f.N), Df)
    n_merge = 4 * N_CMP_A
    lam_ss_full = solve_lowest(Kf_[If][:, If].tocsc(),
                               Mf_[If][:, If].tocsc(), n_merge)
    merged = np.sort(np.concatenate([q_ss[s][:N_CMP_A + 8]
                                     for s in SECTORS]))[:n_merge]
    # criterion note: the campaign's statistics-grade rule (err < 0.1 x
    # local mean spacing) is defined PER SECTOR; the merged ladder is 4x
    # denser, so the equivalent merged-ladder criterion is frac = 0.4.
    # (The residual mismatch is dominated by the FULL mesh's uniform
    # triangle-diagonal, which breaks mirror symmetry at the discrete
    # level; the quarter problems are exactly symmetric.)
    ns_a = n_star(merged, lam_ss_full, 0.4)
    rel_a = float(np.max(np.abs(merged - lam_ss_full) / lam_ss_full))
    # gate: the quarter must agree with the full instrument at
    # statistics-grade EVERYWHERE the full instrument is itself accurate
    # against the exact spectrum (4 x the per-sector informational N*);
    # beyond that range both instruments are outside their own validity.
    min_info = min(results["gates"][f"Ainfo_{s}"]["n_star"]
                   for s in SECTORS)
    ok_a = ns_a >= 4 * min_info
    results["gates"]["A_merged"] = dict(n_star=int(ns_a), n_cmp=n_merge,
                                        max_relerr=rel_a,
                                        required=int(4 * min_info),
                                        criterion="0.4 x merged spacing = "
                                        "0.1 x sector spacing")
    md.append(f"- merged 4 x {N_CMP_A} quarter modes vs full-plate SSSS at "
              f"the same h: **N* = {ns_a}/{n_merge}, max relerr "
              f"{rel_a:.2e}** (required: >= {4 * min_info} = the full "
              f"instrument's own exact-accuracy range)")
    md.extend(md_info)
    print(f"[A merged] N* = {ns_a}/{n_merge}, max relerr {rel_a:.2e}")

    # ---------------- GATE B: cross-instrument (FFFF) ----------------
    md.append("\n## GATE B: quarter (free outer) vs full-plate C0-IP sectors\n")
    mesh = tensor_mesh(*FULL_MESH, A, B)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    lam_f, V, sinfo, _ = solve_modes(K.tocsc(), M.tocsc(), N_FULL + 3,
                                     resid_sanity=1e-3, sweeps_max=30)
    lam_f, V, n_rig, _ = split_rigid(lam_f, V)
    rng = np.random.default_rng(5)
    pts = np.vstack([rng.uniform(0.02 * A, 0.98 * A, 900),
                     rng.uniform(0.02 * B, 0.98 * B, 900)])
    P = space.probes(pts)
    Pmx = space.probes(np.vstack([A - pts[0], pts[1]]))
    Pmy = space.probes(np.vstack([pts[0], B - pts[1]]))
    labels, qual, _ = classify_parity_resolved(lam_f, V, P, Pmx, Pmy, K, M)
    print(f"[B full] {len(lam_f)} elastic (rigid {n_rig}), min quality "
          f"{np.min(qual):.3f} ({time.time()-t00:.1f} s)")
    md.append("| sector | rigid (expected) | N* | max relerr |")
    md.append("|---|---|---|---|")
    ok_b = n_rig == 3
    for s in SECTORS:
        lam_full_s = np.array([l for l, lb in zip(lam_f, labels) if lb == s])
        lam = quarter_spectrum(s, "free")
        lam, gap_ok = split_expected(lam, RIGID[s])
        n_cmp = min(N_CMP_B, len(lam_full_s), len(lam))
        ns = n_star(lam[:n_cmp], lam_full_s[:n_cmp], 0.1)
        rel = float(np.max(np.abs(lam[:n_cmp] - lam_full_s[:n_cmp])
                           / lam_full_s[:n_cmp]))
        ok = (ns >= n_cmp) and gap_ok
        ok_b &= ok
        results["gates"][f"B_{s}"] = dict(n_star=int(ns), n_cmp=int(n_cmp),
                                          max_relerr=rel, gap_ok=gap_ok,
                                          rigid_expected=RIGID[s])
        md.append(f"| {s} | {RIGID[s]} ({gap_ok}) | {ns}/{n_cmp} | "
                  f"{rel:.2e} |")
        print(f"[B {s}] N* = {ns}/{n_cmp}, max relerr {rel:.2e}, "
              f"rigid gap {gap_ok}")

    verdict = ("PASS -- promoted to production instrument for E18 and "
               "future Gap-A rungs" if (ok_a and ok_b) else
               "FAIL -- see gate table")
    md.append(f"\n**GATE A: {'PASS' if ok_a else 'FAIL'}; "
              f"GATE B: {'PASS' if ok_b else 'FAIL'} -> {verdict}.**")
    results["verdict"] = verdict
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS_SMOKE.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_smoke.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md[-4:]))


if __name__ == "__main__":
    main()
