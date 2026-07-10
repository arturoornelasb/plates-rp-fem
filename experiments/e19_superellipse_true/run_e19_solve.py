#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19 solve stage. Usage: run_e19_solve.py <free|ss>_<prod|check>.
Staged/stop-proof (eigenpairs to .npz right after the solve)."""
import json
import os
import sys
import time

import numpy as np

from platefem import solve_lowest, solve_modes
from platefem.c0ip import assemble_c0ip, boundary_dofs

from e19_common import (KFEM, MESH_CHECK, MESH_PROD, N_MODES, NU, RIGID,
                        split_expected, superellipse_polar_mesh)

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    job = sys.argv[1]
    bc, stage = job.split("_")
    t00 = time.time()
    mesh = superellipse_polar_mesh(
        *(MESH_PROD if stage == "prod" else MESH_CHECK))
    K, M, space = assemble_c0ip(mesh, k=KFEM, nu=NU)
    K = 0.5 * (K + K.T)
    if bc == "ss":
        D = boundary_dofs(space)
        I = np.setdiff1d(np.arange(space.N), D)
        K, M = K[I][:, I].tocsc(), M[I][:, I].tocsc()
    else:
        K, M = K.tocsc(), M.tocsc()
    print(f"[{job}] {K.shape[0]} dofs ({time.time()-t00:.1f} s)", flush=True)

    n_req = N_MODES[bc] + RIGID[bc]
    if stage == "prod":
        sanity = 1e-3 if bc == "free" else 1e-4
        lam, V, sinfo, _ = solve_modes(K, M, n_req, resid_sanity=sanity,
                                       sweeps_max=40)
        lam, V, gap_ok, ratio = split_expected(lam, V, RIGID[bc])
        np.savez(os.path.join(HERE, f"eig_{job}.npz"), lam=lam, V=V,
                 resid=[float(sinfo["max_resid"])], ratio=[ratio],
                 ndof=[K.shape[0]])
        print(f"[{job}] {len(lam)} elastic cached, resid "
              f"{sinfo['max_resid']:.1e}, rigid ratio {ratio:.1e} "
              f"({time.time()-t00:.1f} s)")
    else:
        lam = solve_lowest(K, M, n_req)
        lam, _, _, _ = split_expected(lam, None, RIGID[bc])
        np.savez(os.path.join(HERE, f"eig_{job}.npz"), lam=lam,
                 ndof=[K.shape[0]])
        print(f"[{job}] {len(lam)} eigenvalues cached "
              f"({time.time()-t00:.1f} s)")
    with open(os.path.join(HERE, f"stage_{job}.json"), "w") as f:
        json.dump(dict(job=job, done=True,
                       wall_s=round(time.time() - t00, 1)), f)


if __name__ == "__main__":
    main()
