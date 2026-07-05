# E12 -- P7: non-proportional damping -> AI-dagger (RESULTS)

Certified FFFF modes (295/sector), viscous patch quartet at (0.2231, 0.1372), r = 0.06. Baselines at n = 1200.

| reference | <|z|> | -<cos theta> |
|---|---|---|
| Poisson2D | 0.6627 +/- 0.0032 | -0.0227 +/- 0.0170 |
| GinUE | 0.7336 +/- 0.0022 | 0.2523 +/- 0.0076 |
| AI_dagger | 0.7242 +/- 0.0048 | 0.2105 +/- 0.0124 |

| case | <|z|> | -<cos theta> | n | pairing err |
|---|---|---|---|---|
| non-prop gamma*=0.5 (pooled) | 0.5820 | 0.0386 | 585 | 0.0e+00 |
| non-prop gamma*=1 (pooled) | 0.6200 | -0.0923 | 595 | 0.0e+00 |
| non-prop gamma*=2 (pooled) | 0.6770 | -0.0608 | 587 | 0.0e+00 |
| Rayleigh control (pooled) | 0.5512 | 0.0624 | 817 | - |

- marker distances at gamma* = 1: to AI-dagger 0.3202, to Poisson2D 0.0816, to GinUE 0.3628 (marker scale se ~ 0.0284)

**Reading: CHALLENGES (Poisson-like: no 2D repulsion at these widths)**
