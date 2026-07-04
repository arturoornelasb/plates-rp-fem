# E12 -- P7: non-proportional damping and the AI-dagger class (preregistered)

Paper Prediction 7 (non-Hermitian extension): a freely vibrating plate with
NON-proportional damping is a quadratic eigenproblem with real symmetric
(M, C, K); its linearization is COMPLEX-SYMMETRIC, i.e. transpose-symmetric
non-Hermitian -- the AI-dagger universality class -- and the complex
resonances should show 2D level repulsion with the AI-dagger markers of the
complex spacing ratio (Sa-Ribeiro-Prosen), NOT Ginibre and NOT 2D Poisson.
The repo's p7_pretest validated the pipeline on a toy (diagonal K, synthetic
dense C); this experiment realizes it on the TRUE certified plate operator
with a PHYSICAL damping mechanism.

## Design

- Certified eigenpairs of the campaign FFFF rectangle (96x60, 1200 modes,
  banded-deflated vectors); modal K = diag(Lambda), M = I per parity sector.
- Physical non-proportional damping: a localized viscous patch quartet
  (mirror images preserve the Z2xZ2 classification), C_ij = c0 sum_patch
  int phi_i phi_j -- spatially localized, hence dense and non-commuting with
  K in the modal basis (non-proportionality is physical, not synthetic).
- Strength scan: c0 set so the median resonance half-width over the mean
  spacing gamma* in {0.5, 1, 2} (the overlapping-resonance regime).
- Controls: Rayleigh proportional damping (C = alpha M + beta K: resonances
  on a 1D filament -- no 2D statistics), conjugate-pairing integrity,
  synthetic baselines (AI-dagger, GinUE, 2D Poisson) recomputed at matched n
  with the pretest's estimator.

## Preregistered readings

- SUPPORTS P7: pooled (4 sectors) complex-ratio markers <|z|>, -<cos theta>
  of the non-proportional case land on the AI-dagger baseline, separated
  from 2D Poisson AND from the proportional-control filament by >= 3
  combined sigma at n ~ 1200 resonances.
- CHALLENGES: markers Poisson-like (no 2D repulsion) or
  filament-indistinguishable.
- STRETCH (registered): AI-dagger vs GinUE fine discrimination needs
  ~1.5-2k bulk resonances -> v2 with 2400 modes at 128x80 if v1 supports.

First realization target: the AI-dagger class in a classical elastic wave
system built from a certified true-operator model.
