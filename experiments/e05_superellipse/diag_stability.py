"""E5 validity: are the p=3 and p=10 STATISTICS refinement-stable?
Compute the full pooled-<r> protocol on refine-5 meshes and compare with the
production refine-6 values (0.4944 and 0.5159)."""
import json
import time
import numpy as np
from platefem import (assemble_plate, centered_probe_operators,
                      classify_parity_resolved, r_values, solve_modes,
                      split_rigid, superellipse_basis, SECTORS)

ratio = 1.6189043236
a_ax = np.sqrt(ratio); b_ax = 1.0 / a_ax
nu, n_modes, SKIP = 0.33, 1600, 10

with open(r"C:\Github\plates-rp-fem\experiments\e05_superellipse\results_raw.json") as f:
    res = json.load(f)

for p_exp in [3.0, 10.0]:
    t0 = time.time()
    mesh, basis = superellipse_basis(5, a_ax, b_ax, p_exp)
    K, M = assemble_plate(mesh, basis, nu)
    lam, V, info, _ = solve_modes(K, M, n_modes + 3)
    lam, V, n_rigid, _ = split_rigid(lam, V)
    P, Pmx, Pmy = centered_probe_operators(basis, a_ax, b_ax)
    labels, qual, _ = classify_parity_resolved(lam, V, P, Pmx, Pmy, K, M)
    lam = lam[:n_modes]; labels = labels[:n_modes]
    rv = []
    for s in SECTORS:
        ev = np.sort(lam[[i for i, l in enumerate(labels) if l == s]])
        rv.extend(r_values(ev[SKIP:]).tolist())
    rv = np.array(rv)
    r5 = (float(np.mean(rv)), float(np.std(rv) / np.sqrt(len(rv))))
    rec6 = res["runs"][f"p{p_exp:g}"]
    lam6 = np.array(rec6["lam"]); lab6 = rec6["labels"]
    rv6 = []
    for s in SECTORS:
        ev = np.sort(lam6[[i for i, l in enumerate(lab6) if l == s]])
        rv6.extend(r_values(ev[SKIP:]).tolist())
    rv6 = np.array(rv6)
    r6 = (float(np.mean(rv6)), float(np.std(rv6) / np.sqrt(len(rv6))))
    dsig = abs(r6[0] - r5[0]) / np.sqrt(r6[1] ** 2 + r5[1] ** 2)
    print(f"p={p_exp:g}: refine-5 <r> = {r5[0]:.4f} +/- {r5[1]:.4f} "
          f"({K.shape[0]} dofs) vs refine-6 {r6[0]:.4f} +/- {r6[1]:.4f} "
          f"-> {dsig:.1f} sigma apart ({time.time()-t0:.0f} s)")
