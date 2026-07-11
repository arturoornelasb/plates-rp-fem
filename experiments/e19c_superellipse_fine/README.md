# E19c -- superellipse p = 10: the crisp standalone point (preregistered)

Frozen 2026-07-10, before any run. The registered E19b follow-up
("refine-7 + cross-mesh projection machinery"). E19b's INTERMEDIATE
(pooled slope -0.205 +/- 0.121, 1.7 sigma) was limited by exactly two
things: (i) only 668 free modes solved (2-rung ladders; the E5-cross
gate held at FULL coverage and never bound), and (ii) the oo sector's
SS basis (250 < 256) killed its second rung. E19c removes both.

## Instrument

- FREE side (where the windows live): refine-7 init_circle-mapped
  superellipse (~520k P4 C0-IP dofs), 1500 elastic modes with vectors.
  Certification: the E5 certified-Argyris cross-instrument gate
  (independent discretization, 1600 modes available) is PRIMARY, as
  declared in the E19b gate amendment; the refine-6 two-mesh N* (vs the
  cached E19b 668-mode ladder) is reported, informative only.
- SS basis (a representation device): the refine-6 solve EXTENDED to
  2300 modes with vectors (from E19b's 1040). Order gate: penalty
  robustness (sigma_factor 20 vs 10) at spacing_frac 1.0, now via an
  EIGENVALUES-ONLY sigma-20 solve (the gate consumes only eigenvalue
  order; declared here, pre-run).
- CROSS-MESH PROJECTION (the new machinery): refine-6 SS functions are
  interpolated onto the refine-7 P4 dof points (Lagrange point
  evaluation via basis.probes). Refine-7 points in the chord-vs-curve
  sliver (superellipse coordinate sigma > 1 - delta, delta = 2e-3) are
  radially pulled back to sigma = 1 - delta before evaluation: SS
  functions are Dirichlet-zero at the true boundary, so the collar
  perturbation is O(delta * grad w), uniform over the basis.
  [AMENDMENT pre-production, 2026-07-10, at smoke stage: delta is
  SELF-CALIBRATED to the coarse mesh -- 1.5x its worst measured
  boundary chord sag -- instead of the fixed 2e-3, which is
  level-dependent and wrong for the level-5 smoke; and dof locations /
  point evaluation use the nodal-lattice reconstruction + vectorized
  probes of e19c_common (C0IPSpace has no doflocs; per-point probes too
  slow at 520k), verified in-run by the nodal self-test
  (probes-at-dof-points = identity).] The
  interpolated sector family is then re-orthonormalized NESTED (per
  sector, Cholesky of the refine-7 M-Gram in SS eigenvalue order --
  the E18b pattern), so the first-N spans are exact and nested and the
  frozen statistic is unchanged. Gram conditioning reported; a Cholesky
  failure disqualifies the sector (gate).

## Ladders and readings (frozen, inherited verbatim from E19)

Per sector (ee/eo/oe/oo), true-operator windows [0.4N, 0.6N) at
N in {128, 256, 512} (coverage permitting), coefficients in the first-N
same-domain SS sector functions, renormalized; window-mean ln IPR;
weighted slope over covered rungs; Dq spread (q = 1.5, 2, 4) alongside.

- THIRD-POINT-CONFIRMED: sector-consistent negative slopes with pooled
  |slope|/se >= 3; report D2_true(p10) and its hierarchy place
  (prediction: >= triangle's 0.138 +/- 0.029).
- FLAT: pooled |slope| < 0.05.
- INTERMEDIATE otherwise.

Coverage arithmetic (declared): 1500 free / 2300 basis with E19b's
sector shares (oo lowest at ~24%) give every sector >= 308 free and
>= 512 basis -- full 3-rung ladders expected in all four sectors; a
sector missing rung 512 contributes its covered rungs (weighted
pooling), as in E19. Crispness target (declared, not gated): se_pool
~ 0.05-0.07, i.e. a ~2x sharper standalone point; this run SUPERSEDES
E19b's INTERMEDIATE as the standalone p10 number either way.

## Budget

Three detached, staged, stop-proof solves (parallel): free7 ~4-12 h
(the long pole), ss_ext ~4-8 h, sig20_ext (eigenvalues-only) ~1-3 h;
analysis on caches ~minutes. Residual floors: free 3e-2, ss 1e-2 (the
E19 convention: floors confirm convergence-to-floor; accuracy is
arbited by the independent gates above).
