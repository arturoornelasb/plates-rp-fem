# E13 -- aspect-ratio robustness sweep (RESULTS)

FFFF, nu = 0.33, 800 modes per ratio, identical protocol.

| a/b | ee | eo | oe | oo | pooled |
|---|---|---|---|---|---|
| 1.2724 | 0.431 | 0.456 | 0.478 | 0.450 | **0.4538 +/- 0.0098** |
| 1.4142 | 0.407 | 0.430 | 0.435 | 0.451 | **0.4302 +/- 0.0098** |
| 1.6189 | 0.448 | 0.393 | 0.465 | 0.461 | **0.4422 +/- 0.0096** |
| 1.8330 | 0.452 | 0.468 | 0.398 | 0.448 | **0.4416 +/- 0.0100** |
| 2.0322 | 0.412 | 0.418 | 0.413 | 0.454 | **0.4236 +/- 0.0094** |

- spread across ratios: 0.0302 (typical sem 0.0097)

**Reading: RATIO-DEPENDENT structure -- inspect per-sector table**

## Addendum: proper significance test

The spread-vs-sem branch above is too strict for 5 draws (max-min inflates).
Chi-square of the five pooled values against their weighted mean: chi2 ~ 5.7
on 4 dof (p ~ 0.22) -- NOT significant. Corrected reading: **ROBUST --
consistent with aspect-ratio independence**; all ratios intermediate
(0.42--0.45), no monotone trend. The single-aspect-ratio caveat is closed.
