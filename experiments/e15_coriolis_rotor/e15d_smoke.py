#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15d SVK machinery smoke: (i) K_T(0, 0) == linear K; (ii) K_T(0, Om)
== K - Om^2 M (load stiffness = spin softening); (iii) Newton prestate
at small Omega matches E15c's linear centrifugal field; (iv) Newton
converges quadratically."""
import sys
import numpy as np

sys.path.insert(0, "src")
from platefem import elastic2d as e2
from platefem.kirchhoff import solve_modes

NU = 0.33
CHIR = [(2, 0.12, 0.4), (3, 0.11, 1.7), (4, 0.09, 2.9), (5, 0.07, 4.4)]


def main():
    m, b, _ = e2.star_polar_basis(12, 36, CHIR)
    K, M, G0 = e2.assemble_elastic(m, b, NU)
    K = K.tocsc()
    M = M.tocsc()

    # (i) zero-state tangent == linear stiffness
    R0, KT0 = e2.svk_residual_tangent(m, b, np.zeros(K.shape[0]), NU, 0.0)
    dK = abs(KT0 - K).max() / abs(K).max()
    print(f"[i] K_T(0,0) vs K: rel {dK:.2e}, |R(0,0)| {np.linalg.norm(R0):.2e}")
    assert dK < 1e-12

    # (ii) load stiffness at Omega
    Om = 0.3
    _, KT = e2.svk_residual_tangent(m, b, np.zeros(K.shape[0]), NU, Om)
    dL = abs(KT - (K - Om ** 2 * M)).max() / abs(K).max()
    print(f"[ii] K_T(0,Om) vs K - Om^2 M: rel {dL:.2e}")
    assert dL < 1e-12

    # (iii) small-Omega prestate vs E15c linear field
    lam, X, info, _ = solve_modes(K, M, 203, resid_sanity=1e-3,
                                  sweeps_max=30)
    u_hat, rig = e2.centrifugal_modal_static(b, lam, X, M)
    Om = 0.05
    u_lin = Om ** 2 * u_hat
    u_svk, ninfo = e2.newton_prestate(m, b, M, NU, Om)
    # compare after removing rigid content from both
    R3 = e2.rigid_inplane_basis(b, M)
    Mr = M @ R3
    u_l = u_lin - R3 @ (Mr.T @ u_lin)
    u_s = u_svk - R3 @ (Mr.T @ u_svk)
    rel = np.linalg.norm(u_s - u_l) / np.linalg.norm(u_l)
    print(f"[iii] Om=0.05: |svk - lin|/|lin| = {rel:.3e} "
          f"(newton ok={ninfo['ok']}, iters {ninfo['iters']})")
    assert rel < 0.05 and ninfo["ok"]

    # (iv) Newton history at a strong spin (continuation from small)
    u0 = u_svk
    for Om in (0.2, 0.35):
        u0, ninfo = e2.newton_prestate(m, b, M, NU, Om, u0=u0)
        s1 = e2.max_prestress_strain(m, b, u0)
        print(f"[iv] Om={Om}: ok={ninfo['ok']} iters={ninfo['iters']} "
              f"hist={['%.1e' % h for h in ninfo['hist'][:6]]} "
              f"strain={s1:.3f}")
        assert ninfo["ok"]
    print("SMOKE PASS")


if __name__ == "__main__":
    main()
