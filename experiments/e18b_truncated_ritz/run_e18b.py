#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18b -- faithful E9-analog truncated protocol (RUNNER).
See README.md (frozen). Reuses E18's cached eigenpairs."""
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import eval_jacobi

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "e18_triangle_protocol"))
from e18_common import KFEM, L_SIDE, NU, N_PROD, macro_verts, tri_mesh_bary

from platefem import n_star
from platefem.stats import classify_c3v, triangle_probe_operators
from platefem.c0ip import assemble_c0ip, boundary_dofs

HERE = os.path.dirname(os.path.abspath(__file__))
E18 = os.path.join(HERE, "..", "e18_triangle_protocol")

P_MAX = 58
LADDER = dict(A1=[128, 256], A2=[128, 256], E=[128, 256, 512, 1024])
WIN = (0.4, 0.6)
SLOPE_FLAT, SLOPE_FALL = 0.15, -0.15
TRUE_SLOPES = dict(A1=(-0.250, 0.116), A2=(-0.147, 0.183),
                   E=(-0.138, 0.029))          # E18 valid true ladders


# ---------------- Dubiner basis on the macro triangle ----------------
def dubiner_ref(x, y, p_max):
    """Values of the graded Dubiner family at reference points (x, y)
    (reference triangle (0,0)-(1,0)-(0,1)). Returns (N, npts) and the
    degree of each function."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    onemy = np.maximum(1.0 - y, 1e-14)
    eta1 = 2.0 * x / onemy - 1.0
    eta2 = 2.0 * y - 1.0
    vals, degs = [], []
    for d in range(p_max + 1):
        for m in range(d + 1):
            n = d - m
            v = (eval_jacobi(m, 0, 0, eta1) * onemy ** m
                 * eval_jacobi(n, 2 * m + 1, 0, eta2))
            vals.append(v)
            degs.append(d)
    return np.array(vals), np.array(degs)


def to_ref(pts_phys):
    """Physical (macro triangle) -> reference coordinates."""
    v = macro_verts(L_SIDE)
    A = np.column_stack([v[1] - v[0], v[2] - v[0]])
    return np.linalg.solve(A, (pts_phys.T - v[0]).T)


def dof_coords(space):
    """Physical coordinates of every C0IP dof (P4 Lagrange nodes)."""
    mesh, k = space.mesh, space.k
    coords = np.zeros((space.N, 2))
    ref_nodes = space.ref.nodes
    for e in range(mesh.t.shape[1]):
        phys = (space.A[e] @ ref_nodes.T).T + space.b0[e]
        coords[space.elem_dofs[e]] = phys
    return coords


def window_ipr_rows(C, N):
    i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
    iprs = []
    for kk in range(i0, i1):
        d = C[kk, :N]
        pcap = float(d @ d)
        if pcap <= 0:
            continue
        dn = d / np.sqrt(pcap)
        iprs.append(float(np.sum(dn ** 4)))
    ln = np.log(iprs)
    return float(np.mean(ln)), float(np.std(ln) / np.sqrt(len(ln))), len(iprs)


def fit_slope(Ns, ms, ses):
    x = np.log(Ns)
    w = 1.0 / np.maximum(np.array(ses), 1e-6) ** 2
    W = np.sum(w)
    xb = np.sum(w * x) / W
    yb = np.sum(w * np.array(ms)) / W
    sxx = np.sum(w * (x - xb) ** 2)
    return (float(np.sum(w * (x - xb) * (np.array(ms) - yb)) / sxx),
            float(np.sqrt(1.0 / sxx)))


def main():
    t00 = time.time()
    results = {"gates": {}, "sectors": {}}

    # ---------------- FEM operators + cached eigenpairs ----------------
    mesh = tri_mesh_bary(N_PROD)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    D = boundary_dofs(space)
    I = np.setdiff1d(np.arange(space.N), D)
    zs = np.load(os.path.join(E18, "eig_ss_prod.npz"))
    zf = np.load(os.path.join(E18, "eig_free_prod.npz"))
    lam_s, V_s = zs["lam"], zs["V"]
    lam_f, V_f = zf["lam"], zf["V"]
    n_use_s = 1977                                # E18 gates
    print(f"[setup] {space.N} dofs; caches loaded ({time.time()-t00:.1f} s)")

    # classification (SS for the representation basis; free for Ritz sanity)
    P, P_R, P_s_op = triangle_probe_operators(space, L_SIDE, 1500)
    lab_s, _, _ = classify_c3v(lam_s, V_s, P[:, I], P_R[:, I], P_s_op[:, I],
                               K[I][:, I].tocsc(), M[I][:, I].tocsc())
    lab_f, _, _ = classify_c3v(lam_f, V_f, P, P_R, P_s_op,
                               K.tocsc(), M.tocsc())
    print(f"[labels] done ({time.time()-t00:.1f} s)")

    # ---------------- Dubiner interpolants ----------------
    t0 = time.time()
    xy = dof_coords(space)
    ref = to_ref(xy.T)
    Pall, degs = dubiner_ref(ref[0], ref[1], P_MAX)
    Pall = Pall.T                                  # (Ndof, Npoly)
    # validation: degree-4 polys are exactly representable
    rng = np.random.default_rng(2)
    pts = rng.uniform(-0.2, 0.2, (2, 50))
    pref = to_ref(pts)
    inside = (pref[0] > 0.02) & (pref[1] > 0.02) & (pref[0] + pref[1] < 0.98)
    pts = pts[:, inside]
    Pv = space.probes(pts)
    test_ref, _ = dubiner_ref(*to_ref(pts), 4)
    err_interp = float(np.max(np.abs(Pv @ Pall[:, :15] - test_ref.T)))
    assert err_interp < 1e-9, err_interp
    results["gates"]["interp_deg4"] = err_interp
    print(f"[dubiner] {Pall.shape[1]} functions, deg<=4 reproduction "
          f"{err_interp:.1e} ({time.time()-t0:.1f} s)")

    # ---------------- C3v representation + sector projectors ----------
    t0 = time.time()
    nl = 2 * P_MAX
    ii, jj = np.meshgrid(np.arange(nl + 1), np.arange(nl + 1), indexing="ij")
    keep = (ii + jj <= nl) & (ii + jj >= 1)
    Q = np.vstack([(ii[keep] / nl), (jj[keep] / nl)]) * 0.96 + 0.01
    keep2 = Q.sum(axis=0) < 0.99
    Q = Q[:, keep2]
    v = macro_verts(L_SIDE)
    Aff = np.column_stack([v[1] - v[0], v[2] - v[0]])
    Qphys = (Aff @ Q).T + v[0]                     # (nq, 2)
    Vq, _ = dubiner_ref(Q[0], Q[1], P_MAX)
    Vq = Vq.T                                      # (nq, Npoly)
    th = 2 * np.pi / 3
    Rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    Mir = np.diag([-1.0, 1.0])
    Ds = {}
    for name, g in [("r", Rot), ("r2", Rot @ Rot), ("s", Mir),
                    ("sr", Mir @ Rot), ("sr2", Mir @ Rot @ Rot)]:
        Qg = (g.T @ Qphys.T)                       # g^{-1} x  (orthogonal)
        rg = to_ref(Qg)
        Vg, _ = dubiner_ref(rg[0], rg[1], P_MAX)
        Ds[name], *_ = np.linalg.lstsq(Vq, Vg.T, rcond=None)
    Iden = np.eye(Pall.shape[1])
    PA1 = (Iden + Ds["r"] + Ds["r2"] + Ds["s"] + Ds["sr"] + Ds["sr2"]) / 6.0
    PA2 = (Iden + Ds["r"] + Ds["r2"] - Ds["s"] - Ds["sr"] - Ds["sr2"]) / 6.0
    Cs = {}
    for sname, Proj in [("A1", PA1), ("A2", PA2),
                        ("E", Iden - PA1 - PA2)]:
        cols = []
        for d in range(P_MAX + 1):
            idx = np.nonzero(degs == d)[0]
            B = 0.5 * (Proj[np.ix_(idx, idx)] + Proj[np.ix_(idx, idx)].T)
            w, U = np.linalg.eigh(B)
            sel = U[:, w > 0.5]
            if sel.shape[1]:
                block = np.zeros((Pall.shape[1], sel.shape[1]))
                block[idx] = sel
                cols.append(block)
        Cs[sname] = np.concatenate(cols, axis=1)
        print(f"[proj] {sname}: {Cs[sname].shape[1]} functions "
              f"({time.time()-t0:.1f} s)")
    results["gates"]["sector_counts"] = {s: int(Cs[s].shape[1]) for s in Cs}

    # ---------------- per-sector nested Ritz ----------------
    MVf_int = None
    md = ["# E18b -- faithful truncated protocol (RESULTS)\n",
          f"Dubiner p <= {P_MAX} ({Pall.shape[1]} functions), sector counts "
          f"{results['gates']['sector_counts']}; deg<=4 reproduction "
          f"{err_interp:.1e}.\n"]
    verdicts = {}
    for s in ("A1", "A2", "E"):
        t0 = time.time()
        Ps = Pall @ Cs[s]
        n_avail = Ps.shape[1]
        G = Ps.T @ (M @ Ps)
        G = 0.5 * (G + G.T)
        Lch = cholesky(G, lower=True)
        Ks = Ps.T @ (K @ Ps)
        Ks = 0.5 * (Ks + Ks.T)
        Kor = solve_triangular(Lch, solve_triangular(
            Lch, Ks, lower=True).T, lower=True)
        Kor = 0.5 * (Kor + Kor.T)
        # SS sector basis (first n per eigen-order within the sector)
        idx_ss = [i for i, l in enumerate(lab_s[:n_use_s]) if l == s]
        Phi = V_s[:, idx_ss]                       # interior dofs
        # Ritz sanity vs certified free sector eigenvalues
        idx_f = [i for i, l in enumerate(lab_f) if l == s][:20]
        wfull = np.linalg.eigvalsh(Kor)
        rel = np.abs(wfull[:20] - lam_f[idx_f]) / lam_f[idx_f]
        results["gates"][f"ritz_sanity_{s}"] = float(np.median(rel))
        rows = []
        for N in LADDER[s]:
            if N > n_avail or N > len(idx_ss):
                break
            w, U = np.linalg.eigh(Kor[:N, :N])
            i0, i1 = int(WIN[0] * N), int(WIN[1] * N)
            # window eigenvectors -> FEM vectors -> SS representation
            Uw = solve_triangular(Lch[:N, :N], U[:, i0:i1],
                                  lower=True, trans="T")
            Wfem = Ps[:, :N] @ Uw                   # (Ndof, nwin)
            Crep = (Phi.T @ (M @ Wfem)[I]).T        # (nwin, n_ss_sector)
            iprs = []
            for r_ in range(Crep.shape[0]):
                d_ = Crep[r_, :N]
                pcap = float(d_ @ d_)
                if pcap <= 0:
                    continue
                dn = d_ / np.sqrt(pcap)
                iprs.append(float(np.sum(dn ** 4)))
            ln = np.log(iprs)
            rows.append((N, float(np.mean(ln)),
                         float(np.std(ln) / np.sqrt(len(ln)))))
            print(f"[{s}] N={N}: trunc IPR {np.exp(rows[-1][1]):.4f} "
                  f"({time.time()-t0:.1f} s)")
        st_true, se_true = TRUE_SLOPES[s]
        su, se_u = fit_slope([r[0] for r in rows], [r[1] for r in rows],
                             [r[2] for r in rows])
        comb = float(np.sqrt(se_true ** 2 + se_u ** 2))
        if abs(st_true) < SLOPE_FLAT and su < SLOPE_FALL:
            v = "PROTOCOL ARTIFACT"
        elif su < SLOPE_FALL and abs(st_true - su) <= 2 * comb:
            v = "GENUINE SCALING"
        else:
            v = "INTERMEDIATE"
        verdicts[s] = v
        results["sectors"][s] = dict(
            rows=[dict(N=N, mln=m, se=se) for N, m, se in rows],
            slope_trunc=su, se_trunc=se_u,
            slope_true=st_true, se_true=se_true,
            ratio=float(su / st_true) if st_true else None,
            ritz_sanity=results["gates"][f"ritz_sanity_{s}"], verdict=v)
        md.append(f"\n## sector {s} (polys {n_avail}; Ritz sanity "
                  f"{results['gates'][f'ritz_sanity_{s}']:.1e})\n")
        md.append("| N | trunc IPR | (true IPR, E18) |")
        md.append("|---|---|---|")
        e18r = {r["N"]: r for r in json.load(open(os.path.join(
            E18, "results_raw.json")))["sectors"][s]["ladder"]}
        for N, m, se in rows:
            tt = f"{np.exp(e18r[N]['true_mln']):.4f}" if N in e18r else "-"
            md.append(f"| {N} | {np.exp(m):.4f} | {tt} |")
        md.append(f"\n- slope_trunc = {su:+.3f} ({se_u:.3f}); slope_true "
                  f"(E18) = {st_true:+.3f} ({se_true:.3f}); ratio "
                  f"trunc/true = {su/st_true:+.2f}")
        md.append(f"- **{v}**")
        del Ps, G, Ks, Kor
    arts = sum(v == "PROTOCOL ARTIFACT" for v in verdicts.values())
    gens = sum(v == "GENUINE SCALING" for v in verdicts.values())
    if gens >= 1:
        overall = (f"GENUINE SCALING in {gens} sector(s): the unadapted "
                   f"tier carries a real RP-candidate eigenvector phase; "
                   f"superellipse contingency TRIGGERS")
    elif arts >= 2:
        overall = ("PROTOCOL ARTIFACT on the unadapted tier; combined "
                   "with E14/E17, Gap A closes negatively for ALL tiers "
                   "at accessible N (the triangle's weak true slope "
                   "notwithstanding -- see ratio)")
    else:
        overall = "MIXED/INTERMEDIATE -- see per-sector ratios"
    md.append(f"\n**Reading: {overall}**")
    results["verdict"] = dict(per_sector=verdicts, overall=overall)
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md[-8:]))


if __name__ == "__main__":
    main()
