# E18b -- the faithful E9-analog truncated protocol on the triangle (preregistered)

Frozen 2026-07-09, before any run. E18's truncated comparator was invalid
by construction (SS eigenvectors diagonalize the interior operator); the
E9 protocol truncates an INCOMPLETE nested polynomial trial space. E18b
implements exactly that and compares against E18's VALID true-operator
ladders (E-sector slope -0.138 +/- 0.029; A-sectors -0.25/-0.15).

## Protocol

- Trial family: Dubiner/Koornwinder polynomials on the macro triangle,
  ordered by TOTAL DEGREE (the frequency proxy) up to p = 58
  (N_full = 1770), realized as their P4 interpolants on the production
  mesh (dof-coordinate evaluation; interpolation error of the highest
  degrees ~0.7% perturbs Ritz eigenvalues at second order -- negligible
  for eigenvector-structure statistics). Built-in validation: exact
  reproduction of degree <= 4 polynomials through the probe operator
  (1e-10 gate).
- Sector adaptation: C3v group representation matrices computed
  numerically on the polynomial space (least squares on a fine lattice;
  the linear group action preserves total degree, so projectors are
  applied DEGREE-BLOCKWISE, keeping the family graded/nested per
  sector). Sector ladders: A1/A2 rungs {128, 256}; E rungs
  {128, 256, 512, 1024}.
- Per rung N (per sector): FEM-M-orthonormalize the first N sector
  polynomials (Cholesky of the Gram -- lower-triangular, so leading
  submatrices span nested graded subspaces), diagonalize the projected
  free-plate operator (the E9 rung), take window [0.4N, 0.6N)
  eigenvectors, represent them in the first N SS sector eigenmodes
  (M-inner products, full-M then interior slice), renormalize captured,
  window-mean ln IPR -- identical windows and representation as E18's
  true ladders.
- Instrument sanity (reported, loosely gated at 1e-3 relative): the
  lowest ~20 rung-eigenvalues per sector must approximate the certified
  free-plate sector eigenvalues from above (Ritz).

## Preregistered readings (the E18 frozen dichotomy, unchanged, now with
a valid comparator; per sector over shared rungs)

- PROTOCOL ARTIFACT: |slope_true| < 0.15 AND slope_trunc < -0.15.
- GENUINE SCALING: slope_trunc < -0.15 AND
  |slope_true - slope_trunc| <= 2 sigma_combined.
- INTERMEDIATE otherwise.
- Additionally reported (ungated, the attribution number): the ratio
  slope_trunc / slope_true per sector -- E18 found the triangle's true
  operator delocalizes slowly (4.8 sigma from flat in E), so beyond the
  binary verdict the quantitative statement is how much the truncated
  protocol EXAGGERATES a weak real effect.

## Inputs and budget

Reuses E18's cached eigenpairs (eig_free_prod.npz, eig_ss_prod.npz) and
gates (ss n_use = 1977); no new large eigensolves. Runtime ~15-30 min
(Dubiner evaluation, Grams, rung eighs, SS representations), ~8 GB peak.
