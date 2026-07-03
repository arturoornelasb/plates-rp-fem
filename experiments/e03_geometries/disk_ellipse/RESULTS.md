# E3b -- free disk (Poisson control) vs free ellipse (RESULTS)

nu = 0.33, Argyris refine-6 (check refine-5), 800 modes, N_use = 296, ellipse a/b = 1.618904 (area pi). References: Poisson 0.3863, GOE 0.5307.

## Gates

- G1 FEM disk refine 5 (18886 dofs): strict N* = 134, max relerr 1.03e-03 (median 8.05e-04), rigid 3
- G1 FEM disk refine 6 (74630 dofs): strict N* = 480, max relerr 2.01e-04 (median 2.01e-04), rigid 3
- G2 ellipse internal: N* = 296, rigid [3, 3]

## Disk control (semi-analytic)

- pooled one-reflection-class <r> over 1458 ratios: **0.3905 +/- 0.0073** (paper's executed control: 0.386 +/- 0.007 on 1460 levels)

## Operational gate G1s: FEM-disk statistics vs semi-analytic

- refine 5: FEM class <r> = 0.4035 +/- 0.0141 vs semi-analytic 0.4040 +/- 0.0141 on the same 408 levels (0.0 sigma apart; doublet max splitting 6.6e-06): PASS
- refine 6: FEM class <r> = 0.4040 +/- 0.0141 vs semi-analytic 0.4040 +/- 0.0141 on the same 408 levels (0.0 sigma apart; doublet max splitting 2.8e-06): PASS

## Ellipse sectors vs matched disk baseline

| sector | n levels | <r> |
|---|---|---|
| ee | 211 | 0.385 +/- 0.019 |
| eo | 198 | 0.355 +/- 0.021 |
| oe | 201 | 0.376 +/- 0.021 |
| oo | 190 | 0.397 +/- 0.021 |
| **pooled** | 800 | **0.3781 +/- 0.0102** |
| disk baseline (4 x 190-level windows) | - | 0.3985 +/- 0.0101 |

- classifier: counts {'ee': 211, 'eo': 198, 'oe': 201, 'oo': 190, 'xx': 0}, min quality 0.941, resolved 113
- N_use = 800 (operational, G1s-certified); strict internal cut 296 as sensitivity: pooled <r> there = 0.3494 +/- 0.0174

### Spectral-window trend (thirds of each sector ladder)

| window | <r> |
|---|---|
| low | 0.3537 +/- 0.0176 |
| mid | 0.3762 +/- 0.0177 |
| high | 0.4128 +/- 0.0183 |

- low -> high trend: +0.0591 (+2.3 sigma)

## Verdict (preregistered readings)

- disk control Poisson-consistent (semi-analytic): True (0.3905 vs 0.3863)
- FEM-disk statistics gate at production refinement: PASS
- ellipse pooled <r> = 0.3781 +/- 0.0102 vs disk baseline 0.3985 +/- 0.0101: separation -1.4 sigma (threshold 3)
- per-sector separations: ee -0.6, eo -1.9, oe -1.0, oo -0.1 sigma
- spectral trend low -> high: +2.3 sigma; low window sub-Poisson: True

**Reading: AMBIGUOUS -- weak-coupling regime: sub-Poisson at low modes with a significant upward spectral trend. This matches the paper's preregistered caveat that at weak coupling level repulsion 'has not developed and can even reverse sign, so a single weak-coupling run must not be read as refutation'. The registered decider is a longer ladder (v2).**
