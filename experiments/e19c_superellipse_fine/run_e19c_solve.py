#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19c solve stage. Usage: run_e19c_solve.py <free7|ss_ext|sig20_ext>.
Staged/stop-proof (eigenpairs to .npz right after the solve)."""
import json
import os
import sys
import time

import numpy as np

from platefem import solve_lowest, solve_modes
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e19c_common import (KFEM, LEVEL_FREE, LEVEL_SS, N_FREE, N_SS, NU,
                         RIGID_FREE, build_mesh, split_expected)

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    job = sys.argv[1]
    t00 = time.time()
    if job == "free7":
        mesh = build_mesh(LEVEL_FREE)
        K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
        K = 0.5 * (K + K.T)
        K, M = K.tocsc(), M.tocsc()
        print(f"[{job}] {K.shape[0]} dofs ({time.time()-t00:.1f} s)",
              flush=True)
        lam, V, sinfo, _ = solve_modes(K, M, N_FREE + RIGID_FREE,
                                       resid_sanity=3e-2, sweeps_max=40)
        lam, V, gap_ok, ratio = split_expected(lam, V, RIGID_FREE)
        np.savez(os.path.join(HERE, "eig_free7.npz"), lam=lam, V=V,
                 resid=[float(sinfo["max_resid"])], ratio=[ratio],
                 ndof=[K.shape[0]])
        print(f"[{job}] {len(lam)} elastic cached, resid "
              f"{sinfo['max_resid']:.1e}, rigid ratio {ratio:.1e} "
              f"({time.time()-t00:.1f} s)")
    else:
        sigma_factor = 20.0 if job == "sig20_ext" else 10.0
        mesh = build_mesh(LEVEL_SS)
        K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU,
                                    sigma_factor=sigma_factor)
        K = 0.5 * (K + K.T)
        D = boundary_dofs(space)
        I = np.setdiff1d(np.arange(space.N), D)
        K, M = K[I][:, I].tocsc(), M[I][:, I].tocsc()
        print(f"[{job}] {K.shape[0]} dofs ({time.time()-t00:.1f} s)",
              flush=True)
        if job == "ss_ext":
            lam, V, sinfo, _ = solve_modes(K, M, N_SS, resid_sanity=1e-2,
                                           sweeps_max=40)
            np.savez(os.path.join(HERE, "eig_ss_ext.npz"), lam=lam, V=V,
                     resid=[float(sinfo["max_resid"])],
                     ndof=[K.shape[0]])
            print(f"[{job}] {len(lam)} cached, resid "
                  f"{sinfo['max_resid']:.1e} ({time.time()-t00:.1f} s)")
        else:
            lam = solve_lowest(K, M, N_SS)
            np.savez(os.path.join(HERE, "eig_sig20_ext.npz"), lam=lam,
                     ndof=[K.shape[0]])
            print(f"[{job}] {len(lam)} eigenvalues cached "
                  f"({time.time()-t00:.1f} s)")
    with open(os.path.join(HERE, f"stage_{job}.json"), "w") as f:
        json.dump(dict(job=job, done=True,
                       wall_s=round(time.time() - t00, 1)), f)


if __name__ == "__main__":
    main()
