#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step1_validate.py -- STEP 1: validate an open-source FEM stack for the Kirchhoff plate
======================================================================================
Goal: establish that scikit-fem on this machine reproduces the plate eigenproblem
to statistics-grade accuracy before using it for the paper's COMSOL-targeted tests
(Protocol B, non-rectangular geometries, Gap A on the true operator).

Geometry throughout = the repo pretest conventions:
a = 1, b = 1/1.6189043236, nu = 0.33, D = rho*h = 1.

Four parts:

  (A) Morley SSSS vs the EXACT spectrum Lambda_mn = ((m pi/a)^2 + (n pi/b)^2)^2.
      Machinery check incl. the essential-BC path (w = 0 at boundary vertices).
  (B) Morley FFFF vs the repo's spectral Rayleigh-Ritz reference (T1 protocol:
      NLEG 64 vs 72, tol 1e-4; union of the four parity sectors). Free edges are
      natural BCs: no constraints. Expect exactly 3 rigid modes.
  (C) Argyris FFFF vs the same Ritz reference. PRODUCTION configuration: the
      quintic C1-conforming Argyris element converges eigenvalues at high order,
      which is what the high mode ladders of Gap A need. No essential BCs.
  (D) Argyris + Winkler edge spring (kappa * w * v on the boundary) at large
      kappa vs the exact SSSS spectrum. This validates the Protocol-B boundary
      term (step 2 uses it) and anchors Argyris against an exact reference:
      free edge (M_n = 0, V_n = 0) + spring -> V_n + kappa w = 0; kappa -> inf
      recovers simply supported (w = 0, M_n = 0) with O(1/kappa) error.

Pass criterion (statistics-grade, per the spectral program): eigenvalue error
must be small against the LOCAL MEAN LEVEL SPACING -- err_i < 0.1 * spacing_i
for all compared modes -- plus tight absolute parity for the production element.
"""
import json
import os
import time

import numpy as np
from scipy.sparse.linalg import eigsh

from skfem import Basis, BilinearForm, FacetBasis, MeshTri
from skfem.element import ElementTriArgyris, ElementTriMorley
from skfem.helpers import dd, ddot, trace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ritz_reference as ritz

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    a=1.0,
    b=1.0 / 1.6189043236,
    nu=0.33,
    meshes_morley=[(64, 40), (128, 80), (256, 160)],
    meshes_argyris=[(32, 20), (64, 40), (128, 80)],
    n_modes=200,                                 # elastic modes to compare (A,B,C)
    n_modes_winkler=100,                         # modes for the kappa check (D)
    winkler_mesh=(64, 40),
    winkler_kappas=[1e6, 1e8, 1e10],
    ritz_nleg=64, ritz_nleg_conv=72, ritz_tol=1e-4,   # repo T1 protocol; see note
    sigma=-10.0,                                 # shift-invert target (below spectrum)
    rigid_tol=1e-6,
    spacing_frac=0.1,                            # statistics-grade: err < frac * spacing
)


# ------------------------------- forms -------------------------------
def make_forms(nu):
    @BilinearForm
    def bending(u, v, w):
        return (1.0 - nu) * ddot(dd(u), dd(v)) + nu * trace(dd(u)) * trace(dd(v))

    @BilinearForm
    def mass(u, v, w):
        return u * v

    return bending, mass


@BilinearForm
def boundary_mass(u, v, w):
    return u * v


def assemble(nx, ny, a, b, nu, element_cls):
    # NOTE: ElementGlobal subclasses (Morley, Argyris) cache mesh-dependent
    # arrays on the instance -- always instantiate fresh per mesh.
    xs = np.linspace(0.0, a, nx + 1)
    ys = np.linspace(0.0, b, ny + 1)
    mesh = MeshTri.init_tensor(xs, ys)
    basis = Basis(mesh, element_cls())
    bending, mass = make_forms(nu)
    K = bending.assemble(basis)
    M = mass.assemble(basis)
    return mesh, basis, K, M


def solve_lowest(K, M, k, sigma):
    vals = eigsh(K.tocsc(), k=k, M=M.tocsc(), sigma=sigma, which="LM",
                 return_eigenvectors=False)
    return np.sort(vals)


def split_rigid(lam, n_expected):
    """Gap-based rigid-mode detection. The numerically built Argyris basis
    (ElementGlobal) pollutes the near-zero cluster at fine mesh (up to ~1e-2
    absolute at 128x80), so a fixed tolerance misclassifies; the multiplicative
    gap rigid -> first elastic stays > 4 decades and is robust. Returns the
    elastic spectrum, the rigid count, and the largest rigid eigenvalue
    (pollution quality metric; clean would be ~1e-10)."""
    head = np.maximum(np.abs(lam[:10]), 1e-300)
    ratios = head[1:] / head[:-1]
    i_gap = int(np.argmax(ratios))
    n_rigid = i_gap + 1 if ratios[i_gap] > 1e3 else 0
    rigid_max = float(np.abs(lam[n_rigid - 1])) if n_rigid else 0.0
    return lam[n_rigid:], n_rigid, rigid_max


# --------------------------- references / metrics ---------------------------
def ssss_exact(a, b, n_modes):
    mmax = int(np.ceil(np.sqrt(n_modes) * 8))
    m = np.arange(1, mmax + 1)
    kx2 = (m * np.pi / a) ** 2
    ky2 = (m * np.pi / b) ** 2
    lam = (kx2[:, None] + ky2[None, :]) ** 2
    return np.sort(lam.ravel())[:n_modes]


def relerr(x, ref):
    n = min(len(x), len(ref))
    return np.abs(x[:n] - ref[:n]) / np.abs(ref[:n])


def local_spacing(ref):
    """Smoothed local mean spacing of the reference spectrum."""
    d = np.gradient(ref)
    w = 11
    kern = np.ones(w) / w
    pad = np.pad(d, (w // 2, w // 2), mode="edge")
    return np.convolve(pad, kern, mode="valid")[:len(ref)]


def n_star(lam, ref, frac):
    """Largest n such that ALL modes 1..n have |err| < frac * local spacing."""
    n = min(len(lam), len(ref))
    ok = np.abs(lam[:n] - ref[:n]) < frac * local_spacing(ref)[:n]
    bad = np.where(~ok)[0]
    return int(bad[0]) if len(bad) else n


def richardson(lams):
    n = min(len(x) for x in lams)
    c, m, f = (x[:n] for x in lams[-3:])
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.log2(np.abs((c - m) / (m - f)))
    p = np.where(np.isfinite(p), p, 2.0)
    p_use = np.clip(np.median(p), 0.5, 8.0)
    ext = f + (f - m) / (2.0 ** p_use - 1.0)
    return ext, p, float(p_use)


# ------------------------------- runners -------------------------------
def run_mesh_family(tag, element_cls, case, meshes, lam_ref, cfg, md, results):
    """Run one element/case over nested meshes vs a reference spectrum."""
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    n_cmp = min(cfg["n_modes"], len(lam_ref))
    lams, rows, rigid_info = [], [], []
    print(f"[{tag}] {case}, {element_cls.__name__}")
    for nx, ny in meshes:
        t0 = time.time()
        mesh, basis, K, M = assemble(nx, ny, a, b, nu, element_cls)
        if case == "SSSS":
            D = np.unique(basis.get_dofs().nodal["u"])
            I = np.setdiff1d(np.arange(K.shape[0]), D)
            K, M = K[I][:, I], M[I][:, I]
            lam = solve_lowest(K, M, n_cmp, cfg["sigma"])
            n_rigid, rigid_max = 0, 0.0
        else:
            lam = solve_lowest(K, M, n_cmp + 3, cfg["sigma"])
            lam, n_rigid, rigid_max = split_rigid(lam, 3)
        dt = time.time() - t0
        lams.append(lam)
        rigid_info.append(n_rigid)
        e = relerr(lam, lam_ref)
        ns = n_star(lam, lam_ref, cfg["spacing_frac"])
        rows.append(dict(nx=nx, ny=ny, ndof=int(K.shape[0]), wall_s=round(dt, 1),
                         max_relerr=float(np.max(e[:n_cmp])),
                         med_relerr=float(np.median(e[:n_cmp])),
                         n_star=ns, n_rigid=n_rigid, rigid_max=rigid_max))
        print(f"    {nx}x{ny} ({K.shape[0]} dofs, {dt:.1f} s): "
              f"max relerr {np.max(e[:n_cmp]):.2e}, median {np.median(e[:n_cmp]):.2e}, "
              f"N* = {ns}/{n_cmp}, rigid {n_rigid} (max {rigid_max:.1e})")

    ext, p_arr, p_use = richardson(lams) if len(lams) >= 3 else (lams[-1], None, None)
    e_ext = relerr(ext, lam_ref)
    ns_ext = n_star(ext, lam_ref, cfg["spacing_frac"])

    md.append(f"\n## ({tag}) {case} -- {element_cls.__name__}\n")
    md.append(f"| mesh | dofs | wall (s) | max relerr ({n_cmp}) | median | N* (err < "
              f"{cfg['spacing_frac']} x spacing) |")
    md.append("|---|---|---|---|---|---|")
    for r in rows:
        md.append(f"| {r['nx']}x{r['ny']} | {r['ndof']} | {r['wall_s']} "
                  f"| {r['max_relerr']:.2e} | {r['med_relerr']:.2e} | {r['n_star']}/{n_cmp} |")
    if p_use is not None:
        md.append(f"| Richardson (p = {p_use:.2f}) | - | - | {np.max(e_ext[:n_cmp]):.2e} "
                  f"| {np.median(e_ext[:n_cmp]):.2e} | {ns_ext}/{n_cmp} |")
        md.append(f"\n- observed convergence order: median {np.median(p_arr):.2f} "
                  f"(IQR {np.percentile(p_arr, 25):.2f}--{np.percentile(p_arr, 75):.2f})")
    if case == "FFFF":
        md.append(f"- rigid modes per mesh (expected 3): {rigid_info}; largest rigid "
                  f"eigenvalue (pollution metric, clean ~1e-10): "
                  + ", ".join(f"{r['rigid_max']:.1e}" for r in rows))

    results[f"{tag}_{case}_{element_cls.__name__}"] = dict(
        meshes=rows, p_used=p_use,
        p_observed_median=(float(np.median(p_arr)) if p_arr is not None else None),
        relerr_extrap=dict(max=float(np.max(e_ext[:n_cmp])),
                           median=float(np.median(e_ext[:n_cmp]))),
        n_star_extrap=ns_ext, n_star_finest=rows[-1]["n_star"],
        lam_finest=lams[-1][:n_cmp].tolist(), lam_extrap=ext[:n_cmp].tolist(),
    )
    return lams, ext, rows


def run_winkler(cfg, md, results, lam_exact):
    """(D) Argyris + Winkler spring at large kappa vs exact SSSS."""
    a, b, nu = cfg["a"], cfg["b"], cfg["nu"]
    nx, ny = cfg["winkler_mesh"]
    n_cmp = cfg["n_modes_winkler"]
    print(f"[D] Winkler spring kappa sweep, Argyris, mesh {nx}x{ny}")
    mesh, basis, K, M = assemble(nx, ny, a, b, nu, ElementTriArgyris)
    fb = FacetBasis(mesh, ElementTriArgyris())
    Bmat = boundary_mass.assemble(fb)

    md.append(f"\n## (D) Protocol-B boundary term -- Argyris + Winkler spring, "
              f"mesh {nx}x{ny}, vs exact SSSS\n")
    md.append("| kappa | max relerr vs SSSS (100) | median | N* |")
    md.append("|---|---|---|---|")
    out = []
    for kap in cfg["winkler_kappas"]:
        t0 = time.time()
        lam = solve_lowest(K + kap * Bmat, M, n_cmp, cfg["sigma"])
        e = relerr(lam, lam_exact[:n_cmp])
        ns = n_star(lam, lam_exact[:n_cmp], cfg["spacing_frac"])
        out.append(dict(kappa=kap, max_relerr=float(np.max(e)),
                        med_relerr=float(np.median(e)), n_star=ns,
                        wall_s=round(time.time() - t0, 1)))
        md.append(f"| {kap:.0e} | {np.max(e):.2e} | {np.median(e):.2e} | {ns}/{n_cmp} |")
        print(f"    kappa {kap:.0e}: max relerr {np.max(e):.2e}, "
              f"median {np.median(e):.2e}, N* {ns}/{n_cmp} ({out[-1]['wall_s']} s)")
    md.append("\n- expected: error falls ~ 1/kappa toward the exact simply "
              "supported limit (spring: V_n + kappa w = 0, M_n = 0 throughout)")
    results["D_winkler"] = out
    return out


# ------------------------------------ main ------------------------------------
def main():
    t00 = time.time()
    cfg = CFG
    md = ["# Step 1 -- open-source FEM validation for the Kirchhoff rectangle\n",
          f"Geometry a = {cfg['a']}, b = {cfg['b']:.10f} (a/b = {cfg['a']/cfg['b']:.6f}), "
          f"nu = {cfg['nu']}, D = rho*h = 1. {cfg['n_modes']} elastic modes compared; "
          f"statistics-grade criterion: err < {cfg['spacing_frac']} x local mean spacing.\n"]
    results = {"config": {k: (list(v) if isinstance(v, (list, tuple)) else v)
                          for k, v in cfg.items()}}

    lam_exact = ssss_exact(cfg["a"], cfg["b"], cfg["n_modes"])

    t0 = time.time()
    spectra = ritz.converged_spectrum(cfg["a"], cfg["b"], cfg["nu"],
                                      nleg=cfg["ritz_nleg"],
                                      nleg_conv=cfg["ritz_nleg_conv"],
                                      tol=cfg["ritz_tol"])
    pooled, pooled_unc, _ = ritz.union_spectrum(spectra)
    lam_ritz = pooled[:cfg["n_modes"]]
    unc_ritz = pooled_unc[:cfg["n_modes"]]
    md.append(f"Ritz reference (T1 protocol, NLEG {cfg['ritz_nleg']} vs "
              f"{cfg['ritz_nleg_conv']}): converged modes/sector "
              + ", ".join(f"{s}: {spectra[s]['n_conv']}" for s in ritz.SECTORS)
              + f"; union {len(pooled)} modes ({time.time()-t0:.1f} s). "
              f"Reference own uncertainty over the {len(lam_ritz)} compared modes: "
              f"median {np.median(unc_ritz):.1e}, max {np.max(unc_ritz):.1e}. "
              f"NOTE (diag_ritz_nleg.py): the Ritz eigenvalues are NON-monotone in "
              f"NLEG beyond ~96 (lambda_1 ee: 450.898 at 48 -> 451.571 at 120), "
              f"impossible for a nested variational basis -- round-off corrupts the "
              f"high-order Legendre assembly, so ~1e-4 is this reference's floor and "
              f"NLEG 64/72 (the repo choice) is its sweet spot. FEM-vs-Ritz "
              f"differences at or below ~1e-4 therefore mean agreement to within "
              f"reference accuracy; the FEM's true accuracy is anchored separately "
              f"by the exact SSSS references (parts A and D).")
    results["ritz_converged"] = {s: int(spectra[s]["n_conv"]) for s in ritz.SECTORS}
    results["ritz_uncertainty"] = dict(median=float(np.median(unc_ritz)),
                                       max=float(np.max(unc_ritz)),
                                       per_mode=unc_ritz.tolist())

    # (A) Morley SSSS vs exact
    lams_a, ext_a, _ = run_mesh_family("A", ElementTriMorley, "SSSS",
                                       cfg["meshes_morley"], lam_exact, cfg, md, results)
    # (B) Morley FFFF vs Ritz
    lams_b, ext_b, _ = run_mesh_family("B", ElementTriMorley, "FFFF",
                                       cfg["meshes_morley"], lam_ritz, cfg, md, results)
    # (C) Argyris FFFF vs Ritz  (production configuration)
    lams_c, ext_c, rows_c = run_mesh_family("C", ElementTriArgyris, "FFFF",
                                            cfg["meshes_argyris"], lam_ritz, cfg, md, results)
    # (D) Argyris + Winkler vs exact SSSS
    wink = run_winkler(cfg, md, results, lam_exact)

    # ---------- side-by-side table for the production element ----------
    md.append("\n### First 10 FFFF elastic modes, Argyris finest mesh vs Ritz "
              "(omega = sqrt(Lambda))\n")
    md.append("| mode | omega Ritz | omega FEM | reldiff |")
    md.append("|---|---|---|---|")
    lam_c = np.array(results["C_FFFF_ElementTriArgyris"]["lam_finest"])
    for i in range(10):
        md.append(f"| {i+1} | {np.sqrt(lam_ritz[i]):.6f} | {np.sqrt(lam_c[i]):.6f} "
                  f"| {abs(lam_c[i]-lam_ritz[i])/lam_ritz[i]:.2e} |")

    # ---------- figure ----------
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    n = cfg["n_modes"]
    ax.semilogy(1 + np.arange(n), relerr(lams_b[-1], lam_ritz), ".", ms=3,
                label="Morley FFFF finest vs Ritz")
    ax.semilogy(1 + np.arange(n), relerr(ext_b, lam_ritz), ".", ms=3,
                label="Morley FFFF extrap vs Ritz")
    for lam, (nx, ny) in zip(lams_c, cfg["meshes_argyris"]):
        ax.semilogy(1 + np.arange(min(n, len(lam))), relerr(lam, lam_ritz), ".", ms=3,
                    label=f"Argyris FFFF {nx}x{ny} vs Ritz")
    sp = local_spacing(lam_ritz) * cfg["spacing_frac"] / lam_ritz
    ax.semilogy(1 + np.arange(n), sp, "k--", lw=1,
                label=f"{cfg['spacing_frac']} x local spacing (criterion)")
    ax.semilogy(1 + np.arange(n), unc_ritz, "-", color="gray", lw=1, alpha=0.7,
                label="Ritz reference own uncertainty")
    ax.set_xlabel("elastic mode index")
    ax.set_ylabel("relative eigenvalue error")
    ax.set_title("Step 1: FEM vs spectral Ritz reference (FFFF)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "relerr_vs_mode.png"), dpi=140)

    # ---------- verdict ----------
    argy = results["C_FFFF_ElementTriArgyris"]
    ok_machinery = (results["A_SSSS_ElementTriMorley"]["relerr_extrap"]["median"] < 5e-3
                    and all(r["n_rigid"] == 3 for r in results["B_FFFF_ElementTriMorley"]["meshes"])
                    and all(r["n_rigid"] == 3 for r in argy["meshes"]))
    # production: some Argyris mesh reaches full N* with max relerr < 2e-4
    # (2e-4 = 5x margin below the 0.1-spacing criterion at mode 200)
    ok_production = any(r["n_star"] >= cfg["n_modes"] and r["max_relerr"] < 2e-4
                        for r in argy["meshes"])
    ok_winkler = (wink[-1]["max_relerr"] < wink[0]["max_relerr"]
                  and wink[-1]["n_star"] >= cfg["n_modes_winkler"])
    verdict = "PASS" if (ok_machinery and ok_production and ok_winkler) else "FAIL"
    md.append(f"\n## Verdict: {verdict}\n")
    md.append(f"- machinery (Morley SSSS median extrap < 5e-3; 3 rigid modes everywhere): "
              f"{ok_machinery}")
    best = min(argy["meshes"], key=lambda r: r["max_relerr"])
    md.append(f"- production (Argyris FFFF: full N* and max relerr < 2e-4 on some mesh; "
              f"best: {best['nx']}x{best['ny']} with N* = {best['n_star']}/{cfg['n_modes']}, "
              f"max relerr {best['max_relerr']:.2e}): {ok_production}")
    md.append(f"- Protocol-B term (Winkler error falls with kappa, N* full at "
              f"kappa = {cfg['winkler_kappas'][-1]:.0e}): {ok_winkler}")
    md.append(f"\nTotal wall time: {time.time()-t00:.1f} s")
    results["verdict"] = verdict

    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nVerdict: {verdict}")
    print(f"Wrote RESULTS.md / results.json / relerr_vs_mode.png in {HERE}")


if __name__ == "__main__":
    main()
