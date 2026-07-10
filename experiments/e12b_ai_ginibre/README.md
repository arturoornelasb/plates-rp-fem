# E12b -- the AI-dagger vs Ginibre fine call (preregistered)

Frozen 2026-07-09, before any run. E12 v3 established the 2D-Poisson ->
AI-dagger crossover (distributed dissipation at gamma* ~ 16 within ~3
sigma of the AI-dagger markers, 2D-Poisson excluded at ~5 sigma); the
registered stretch is the FINE discrimination between AI-dagger (complex
symmetric, the P7 prediction for reciprocal damped media) and GinUE
(no symmetry) at ~4x the resonance statistics.

## Design

- Same certified instrument as E12 (Argyris (96, 60), 1203 modes,
  ~295/sector; classifier + modal patch-quartet machinery imported from
  the frozen E12 runner).
- The modal basis is damping-independent: solve ONCE, then pool the
  complex-spacing-ratio markers over 4 independent patch-placement
  realizations (seeds 7, 17, 27, 37) x 4 sectors at the verdict point
  gamma* = 16 (distributed 24-quartet dissipation), ~4600 bulk z-ratios
  (marker scale se ~ 0.010 by the E12 empirical calibration
  0.02 x sqrt(1200/n)). gamma* = 4 pooled alongside for continuity.
- Baselines (Poisson2D / GinUE / AI-dagger) at n = 1200 with 24
  repetitions (reference sems ~3x sharper than E12).

## Preregistered reading (at gamma* = 16, pooled)

Let d_AI and d_Gin be the Euclidean marker distances
(<|z|>, -<cos theta>) to the AI-dagger and GinUE baseline centroids, and
se the combined marker scale (empirical scale + baseline sem in
quadrature).

- AI-DAGGER CONFIRMED: d_AI <= 2 se AND d_Gin >= 4 se.
- GINIBRE: d_Gin <= 2 se AND d_AI >= 4 se.
- UNRESOLVED: otherwise (report both distances; the registered next step
  would be the physical-damping model refinement, not more statistics).

## Budget

~1 h (one 1203-mode certified solve + 16 modal QEPs + baselines).
