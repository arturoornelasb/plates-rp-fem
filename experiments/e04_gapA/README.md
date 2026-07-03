# E4 -- Gap A on the true operator (preregistered before execution)

The paper's decisive test (Sec. 'Gap A', observable 1): compute the
eigenmodes of a free rectangular plate "numerically (COMSOL or finite
elements)", restrict to a single symmetry sector, and measure the
eigenvectors' inverse participation ratios and fractal dimension D2 as a
function of sector size N. Rosenzweig-Porter predicts non-ergodic FRACTAL
eigenvectors (D2 strictly between 0 and 1 with characteristic N-scaling);
a banded coupling matrix instead localizes them. "The verdict rests on
observable 1."

This is the first execution on the TRUE operator: the paper's own in-silico
pretest (and the T1 audit) used Rayleigh-Ritz modes, whose mid-spectrum
window at truncation N is not converged ("inherent to the ladder protocol").
Here the modes are certified FEM eigenmodes of the free plate; the ladder is
purely representational (lowest-N reference-basis products, T1 window
[0.4N, 0.6N)), which removes that caveat.

## Constraints honored (from T1)

- D2 of a finite representation depends on the basis (T1 spread up to 0.78 vs
  0.38): we measure in BOTH registered bases -- SS sine products (the H0 of
  the Gap A recipe) and free-free beam products -- and report both, never the
  FEM nodal basis. Captured norms p are reported; coefficients renormalized.
- Empirical GOE baseline at identical N (8 seeded realizations), giving both
  the ladder D2 (slope) and the fixed-N effective D2, exactly as in T1.
- Eigenvector quality: the banded-deflated solver is run to a tightened
  residual gate (1e-5) since coefficients, not just labels, are consumed.

## Design

Rectangle a = 1, b = 1/1.6189043236, nu = 0.33 (campaign geometry), Argyris
128x80, n_modes = 1200 -> ~300 per parity sector. Gates: G1 = Winkler
kappa = 1e10 vs exact SSSS on the same mesh (all 1200 modes); G2 = internal
vs 96x60. Ladder N in {64, 100, 144, 200, 256, 324}; long-range Sigma^2(L)
to L = 20 per sector vs Poisson and GOE.

## Preregistered readings (paper Sec. gapA + P12)

- RP-SUPPORTING: sector eigenvectors extended but non-ergodic -- ladder D2
  strictly between 0 and 1 in both bases, IPR well above GOE at fixed N but
  falling with N; Sigma^2(L) between Poisson and GOE at L <= 20.
- BANDED/LOCALIZED: IPR ~ N-independent (D2 -> 0) -- would favor the
  banded-matrix rival and, per the paper, change the framework.
- Reference points: P12 preregistered D2 = 0.76 +/- 0.15 (fraction-of-
  ergodic, ladder), T1 Ritz-based values (sine 0.28-0.39, beam 0.20-0.43
  ladder; ~0.35-0.54 fixed-N) with the basis-spread caveat.
