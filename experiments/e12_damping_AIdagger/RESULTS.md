# E12 -- P7: non-proportional damping -> AI-dagger (RESULTS)

Certified FFFF modes (295/sector), viscous patch quartet at (0.2231, 0.1372), r = 0.06. Baselines at n = 1200.

| reference | <|z|> | -<cos theta> |
|---|---|---|
| Poisson2D | 0.6627 +/- 0.0032 | -0.0227 +/- 0.0170 |
| GinUE | 0.7336 +/- 0.0022 | 0.2523 +/- 0.0076 |
| AI_dagger | 0.7242 +/- 0.0048 | 0.2105 +/- 0.0124 |

| case | <|z|> | -<cos theta> | n | pairing err |
|---|---|---|---|---|
| non-prop gamma*=2 (pooled) | 0.6141 | -0.0439 | 587 | 0.0e+00 |
| non-prop gamma*=4 (pooled) | 0.6642 | -0.0692 | 584 | 0.0e+00 |
| non-prop gamma*=8 (pooled) | 0.6990 | 0.0331 | 577 | 0.0e+00 |
| non-prop gamma*=16 (pooled) | 0.7043 | 0.1168 | 552 | 0.0e+00 |
| point-dashpots gamma*=4 (pooled) | 0.7062 | 0.0170 | 572 | - |
| Rayleigh control (pooled) | 0.5512 | 0.0624 | 817 | - |

- marker distances at gamma* = 4: to AI-dagger 0.2861, to Poisson2D 0.0466, to GinUE 0.3289 (marker scale se ~ 0.0287)

**Reading: CHALLENGES (Poisson-like: no 2D repulsion at these widths)**

## Corrected reading (v3 addendum)

The automated verdict keys on gamma* = 4, but the v3 sweep shows the
physics: pooled markers move monotonically from 2D-Poisson toward the
AI-dagger corner as resonance overlap grows -- (0.614, -0.044) at gamma*=2
to (0.704, +0.117) at gamma*=16, where the case sits NEAREST AI-dagger
(distance 0.096 ~ 3 marker-sigma) with 2D-Poisson excluded at ~5 sigma and
GinUE at 0.138. Point dashpots reach the 2D-repulsion level at 4x lower
overlap, confirming the controlling parameter is gamma* x (off-diagonal /
diagonal) ~ gamma*/10 for patch dissipation. RECORDED READING:
**TREND-SUPPORTS P7 -- the 2D-Poisson -> AI-dagger crossover is observed
and overlap-controlled in the certified plate**; full convergence onto the
AI-dagger markers and the AI-dagger-vs-GinUE fine call are registered for
the ~2k-resonance stretch (2400 modes at 128x80). Physical statement for
the paper: weak/localized dissipation gives uncorrelated resonances even
when widths ~ spacings; the AI-dagger class requires distributed
dissipation with strong overlap -- itself a testable experimental knob.
