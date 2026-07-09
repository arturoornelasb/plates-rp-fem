#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18 solve stage (RUNNER). Usage: run_e18_solve.py <free|ss>_<prod|check>.
Staged and stop-proof: eigenpairs cached to .npz the moment the solve
finishes. prod jobs solve WITH vectors; check jobs eigenvalues-only."""
import json
import os
import sys
import time

import numpy as np

from platefem import solve_lowest, solve_modes
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e18_common import (KFEM, NU, N_CHECK, N_MODES, N_PROD, RIGID,
                        split_expected, tri_mesh_bary)

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    job = sys.argv[1]                      # e.g. "ss_prod"
    bc, stage = job.split("_")
    t00 = time.time()
    n_sub = N_PROD if stage == "prod" else N_CHECK
    mesh = tri_mesh_bary(n_sub)
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    I = None
    if bc == "ss":
        D = boundary_dofs(space)
        I = np.setdiff1d(np.arange(space.N), D)
        K, M = K[I][:, I].tocsc(), M[I][:, I].tocsc()
    else:
        K, M = K.tocsc(), M.tocsc()
    print(f"[{job}] {K.shape[0]} dofs ({time.time()-t00:.1f} s)")

    n_req = N_MODES[bc] + RIGID[bc]
    if stage == "prod":
        # free: the rigid cluster's residual floor at ~290k dofs is ~4e-4
        # (assembly roundoff; convergence history clean) -- E15c/E3c-grade
        # threshold; eigenvalue accuracy is gated independently (two-mesh
        # + exact Lame^2 anchor + Argyris cross-check). ss (no rigid
        # modes) keeps the strict E14/E17 threshold.
        sanity = 1e-3 if bc == "free" else 1e-4
        lam, V, sinfo, _ = solve_modes(K, M, n_req, resid_sanity=sanity,
                                       sweeps_max=40)
        lam, V, gap_ok, ratio = split_expected(lam, V, RIGID[bc])
        np.savez(os.path.join(HERE, f"eig_{job}.npz"), lam=lam, V=V,
                 resid=[float(sinfo["max_resid"])], gap_ok=[gap_ok],
                 ratio=[ratio], ndof=[K.shape[0]])
        print(f"[{job}] {len(lam)} elastic cached, resid "
              f"{sinfo['max_resid']:.1e}, rigid ratio {ratio:.1e} "
              f"({time.time()-t00:.1f} s)")
    else:
        lam = solve_lowest(K, M, n_req)
        lam, _, gap_ok, ratio = split_expected(lam, None, RIGID[bc])
        np.savez(os.path.join(HERE, f"eig_{job}.npz"), lam=lam,
                 ratio=[ratio], ndof=[K.shape[0]])
        print(f"[{job}] {len(lam)} eigenvalues cached "
              f"({time.time()-t00:.1f} s)")
    with open(os.path.join(HERE, f"stage_{job}.json"), "w") as f:
        json.dump(dict(job=job, done=True,
                       wall_s=round(time.time() - t00, 1)), f)


if __name__ == "__main__":
    main()
