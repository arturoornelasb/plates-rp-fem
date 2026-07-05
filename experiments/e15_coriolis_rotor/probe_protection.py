#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15 decisive probe: the PROTECTION theorem on structured symmetric meshes.

Physics (the corrected mechanism, T5): sigma_v ANTI-commutes with the Coriolis
form G0 (rotation maps sigma_v -> -sigma_v), so sigma_v*T survives cranking and
pins a MIRROR-symmetric rotor to the ORTHOGONAL class (GOE) at all Omega. R_pi
COMMUTES with G0, so R_pi*T does NOT survive -> a chiral rotor reaches the
UNITARY class (GUE) under rotation. Prediction:
  D_chir (no symmetry)         : GOE at Omega=0 -> GUE under rotation.
  D_mirror (exact sigma_v mesh) : stays GOE (~0.53) at all Omega (protected).
The first probe saw D_mirror reach GUE because its init_circle mesh was NOT
mirror-symmetric (sigma_v not discrete-exact). Here the mesh is structured and
symmetric, so sigma_v is exact.
"""
import sys, time
import numpy as np
sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes, split_rigid
from platefem.stats import mean_r, R_POISSON, R_GOE

NU = 0.33
R_GUE = 0.5996
NRINGS, N_TH = 20, 60
NMODES = 340

# strong, rich deformation for a chaotic Omega=0 baseline
CHIR = [(2, 0.15, 0.5), (3, 0.13, 1.7), (4, 0.10, 2.9), (5, 0.08, 0.9)]  # no sym
MIRR = [(2, 0.15, 0.0), (3, 0.13, 0.0), (4, 0.10, 0.0), (5, 0.08, 0.0)]  # sigma_v(x)


def build(harm):
    mesh, basis, meta = e2.star_polar_basis(NRINGS, N_TH, harm)
    K, M, G0 = e2.assemble_elastic(mesh, basis, NU)
    lam, X, info, _ = solve_modes(K, M, NMODES)
    Lam, G0m, Xn = e2.modal_reduce(K, M, G0, X)
    el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
    o = np.argsort(Lam[el])
    Lam_e = Lam[el][o]
    G0m_e = G0m[np.ix_(el, el)][np.ix_(o, o)]
    return mesh, basis, K, M, G0, Lam, X, Xn, Lam_e, G0m_e, info


def sweep(Lam_e, G0m_e, cOs):
    dsp = np.mean(np.diff(np.sort(Lam_e)))
    sq = np.sqrt(dsp)
    out = []
    for cO in cOs:
        res = e2.solve_rotor(Lam_e, G0m_e, cO * sq)
        w = res["omega"]
        r, sem, n = mean_r(w, skip_low=max(10, len(w) // 10))
        out.append((cO, r, sem, res["max_imag"]))
    return out


cOs = [0.0, 0.5, 1.0, 2.0, 4.0]
print(f"config: polar mesh {NRINGS}x{N_TH}, {NMODES} modes, nu={NU}")
print(f"refs: Poisson {R_POISSON}  GOE {R_GOE}  GUE {R_GUE}\n")

# --- D_chir ---
t0 = time.time()
*_, Lam_c, G0_c, info_c = build(CHIR)
print(f"[D_chir] N_el={len(Lam_c)} resid={info_c['max_resid']:.1e} dt={time.time()-t0:.0f}s"
      f"  (no symmetry -> single class)")
for cO, r, sem, mi in sweep(Lam_c, G0_c, cOs):
    print(f"   c_Omega={cO:4.1f}  <r>={r:.4f} +/- {sem:.4f}")

# --- D_mirror: check sigma_v is exact, then per-class + pooled sweep ---
t0 = time.time()
mesh, basis, K, M, G0, Lam, X, Xn, Lam_m, G0_m, info_m = build(MIRR)
labm, cm = e2.parity_classes(basis, Xn[:, :], "mirror_x", MIRR)
# restrict labels to the elastic, frequency-ordered set
el = np.abs(Lam) > 1e-6 * np.abs(Lam).max()
o = np.argsort(Lam[el])
lab_e = labm[el][o]
clean = float(np.mean(np.abs(cm[el]) > 0.9))
print(f"\n[D_mirror] N_el={len(Lam_m)} resid={info_m['max_resid']:.1e} dt={time.time()-t0:.0f}s")
print(f"   sigma_v classification: {int(np.sum(lab_e>0))} even / {int(np.sum(lab_e<0))} odd"
      f" / {int(np.sum(lab_e==0))} ambiguous ; |c|>0.9 frac = {clean:.2f}")

# per-class <r> at Omega=0 (each should be GOE if chaotic)
for s, nm in [(1, "even"), (-1, "odd")]:
    ev = np.sort(Lam_m[lab_e == s])
    r, sem, n = mean_r(np.sqrt(np.abs(ev)))
    print(f"   Omega=0  sigma_v-{nm}: <r>={r:.4f} +/- {sem:.4f}  (n={n})")

# pooled sweep: protection => stays GOE, does NOT reach GUE
print("   pooled <r> vs Omega (protection => plateau at GOE ~0.53, NOT 0.60):")
for cO, r, sem, mi in sweep(Lam_m, G0_m, cOs):
    print(f"     c_Omega={cO:4.1f}  <r>={r:.4f} +/- {sem:.4f}")
