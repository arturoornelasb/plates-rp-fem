# E3b -- free disk (Poisson control) vs free ellipse (RESULTS)

nu = 0.33, Argyris refine-6 (check refine-5), 1600 modes, N_use = 296, ellipse a/b = 1.618904 (area pi). References: Poisson 0.3863, GOE 0.5307.

## Gates

- G1 FEM disk refine 5 (18886 dofs): strict N* = 134, max relerr 3.84e-03 (median 8.79e-04), rigid 3
- G1 FEM disk refine 6 (74630 dofs): strict N* = 480, max relerr 2.15e-04 (median 2.01e-04), rigid 3
- G2 ellipse internal: N* = 296, rigid [3, 3]

## Disk control (semi-analytic)

- pooled one-reflection-class <r> over 1458 ratios: **0.3905 +/- 0.0073** (paper's executed control: 0.386 +/- 0.007 on 1460 levels)

## Operational gate G1s: FEM-disk statistics vs semi-analytic

- refine 5: FEM class <r> = 0.3956 +/- 0.0097 vs semi-analytic 0.3953 +/- 0.0098 on the same 812 levels (0.0 sigma apart; upper half 0.3877 vs 0.3866, 0.1 sigma; doublet max splitting 1.2e-03): PASS
- refine 6: FEM class <r> = 0.3953 +/- 0.0098 vs semi-analytic 0.3953 +/- 0.0098 on the same 812 levels (0.0 sigma apart; upper half 0.3868 vs 0.3866, 0.0 sigma; doublet max splitting 2.8e-06): PASS

## Ellipse sectors vs matched disk baseline

| sector | n levels | <r> |
|---|---|---|
| ee | 413 | 0.384 +/- 0.013 |
| eo | 397 | 0.369 +/- 0.015 |
| oe | 402 | 0.368 +/- 0.014 |
| oo | 386 | 0.371 +/- 0.014 |
| **pooled** | 1600 | **0.3732 +/- 0.0071** |
| disk baseline (4 x 389-level windows) | - | 0.3912 +/- 0.0070 |

- classifier: counts {'ee': 413, 'eo': 397, 'oe': 402, 'oo': 386, 'xx': 2}, min quality 0.160, resolved 342
- N_use = 1600 (operational, G1s-certified); strict internal cut 296 as sensitivity: pooled <r> there = 0.3494 +/- 0.0174

### Spectral-window trend (thirds of each sector ladder)

| window | <r> |
|---|---|
| low | 0.3613 +/- 0.0123 |
| mid | 0.3886 +/- 0.0124 |
| high | 0.3695 +/- 0.0121 |

- low -> high trend: +0.0082 (+0.5 sigma)

## Verdict (preregistered readings)

- disk control Poisson-consistent (semi-analytic): True (0.3905 vs 0.3863)
- FEM-disk statistics gate at production refinement: PASS
- ellipse pooled <r> = 0.3732 +/- 0.0071 vs disk baseline 0.3912 +/- 0.0070: separation -1.8 sigma (threshold 3)
- per-sector separations: ee -0.5, eo -1.4, oe -1.4, oo -1.3 sigma
- spectral trend low -> high: +0.5 sigma; low window sub-Poisson: True

**Reading: CHALLENGES (ellipse consistent with the separable control)**

Discussion. The split E3 verdict (triangle strongly intermediate at 6.8 sigma, ellipse Poisson-consistent at 1600 modes) is itself informative: 'broken biharmonic separability' is not a binary switch. In the paper's own RP framework the coupling magnitude lambda is a parameter -- a dense but TINY V yields RP statistics indistinguishable from Poisson at accessible N. The data point to boundary corners as the strong evanescent-coupling drivers (the triangle's 60-degree corners vs the ellipse's smooth boundary, which behaves quasi-separably here). What this challenges is the naive reading 'any broken separability -> visibly intermediate'; what would decide between 'no coupling' and 'tiny dense coupling' is the eigenvector diagnostic of Gap A applied to the ellipse -- registered as a follow-up. Noise blind spot closed by measurement: the non-smooth FEM error component on the disk is max-median = 1.4e-5 relative, ~0.3% of a sector-level spacing at mode 1600, far too small to wash out genuine repulsion.
