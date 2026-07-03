# Step 1 -- open-source FEM validation for the Kirchhoff rectangle

Geometry a = 1.0, b = 0.6177017291 (a/b = 1.618904), nu = 0.33, D = rho*h = 1. 200 elastic modes compared; statistics-grade criterion: err < 0.1 x local mean spacing.

Ritz reference (T1 protocol, NLEG 64 vs 72): converged modes/sector ee: 171, eo: 162, oe: 174, oo: 164; union 652 modes (0.9 s). Reference own uncertainty over the 200 compared modes: median 8.4e-10, max 8.4e-05. NOTE (diag_ritz_nleg.py): the Ritz eigenvalues are NON-monotone in NLEG beyond ~96 (lambda_1 ee: 450.898 at 48 -> 451.571 at 120), impossible for a nested variational basis -- round-off corrupts the high-order Legendre assembly, so ~1e-4 is this reference's floor and NLEG 64/72 (the repo choice) is its sweet spot. FEM-vs-Ritz differences at or below ~1e-4 therefore mean agreement to within reference accuracy; the FEM's true accuracy is anchored separately by the exact SSSS references (parts A and D).

## (A) SSSS -- ElementTriMorley

| mesh | dofs | wall (s) | max relerr (200) | median | N* (err < 0.1 x spacing) |
|---|---|---|---|---|---|
| 64x40 | 10241 | 2.1 | 2.58e-01 | 1.49e-01 | 7/200 |
| 128x80 | 40961 | 14.4 | 9.04e-02 | 4.33e-02 | 20/200 |
| 256x160 | 163841 | 77.2 | 2.46e-02 | 1.13e-02 | 38/200 |
| Richardson (p = 1.76) | - | - | 6.25e-03 | 1.88e-03 | 81/200 |

- observed convergence order: median 1.76 (IQR 1.60--1.87)

## (B) FFFF -- ElementTriMorley

| mesh | dofs | wall (s) | max relerr (200) | median | N* (err < 0.1 x spacing) |
|---|---|---|---|---|---|
| 64x40 | 10449 | 2.4 | 2.14e-01 | 1.05e-01 | 17/200 |
| 128x80 | 41377 | 14.5 | 6.76e-02 | 2.96e-02 | 30/200 |
| 256x160 | 164673 | 77.4 | 1.80e-02 | 7.65e-03 | 46/200 |
| Richardson (p = 1.82) | - | - | 4.53e-03 | 8.84e-04 | 131/200 |

- observed convergence order: median 1.82 (IQR 1.69--1.92)
- rigid modes per mesh (expected 3): [3, 3, 3]; largest rigid eigenvalue (pollution metric, clean ~1e-10): 2.2e-09, 3.7e-08, 1.3e-06

## (C) FFFF -- ElementTriArgyris

| mesh | dofs | wall (s) | max relerr (200) | median | N* (err < 0.1 x spacing) |
|---|---|---|---|---|---|
| 32x20 | 6130 | 2.8 | 1.07e-04 | 3.47e-06 | 200/200 |
| 64x40 | 23774 | 18.4 | 9.91e-05 | 9.60e-08 | 200/200 |
| 128x80 | 93622 | 106.0 | 7.48e-05 | 3.75e-06 | 200/200 |
| Richardson (p = 0.50) | - | - | 1.13e-04 | 1.30e-05 | 200/200 |

- observed convergence order: median 0.00 (IQR -3.22--2.43)
- rigid modes per mesh (expected 3): [3, 3, 3]; largest rigid eigenvalue (pollution metric, clean ~1e-10): 2.1e-09, 3.2e-06, 1.4e-02

## (D) Protocol-B boundary term -- Argyris + Winkler spring, mesh 64x40, vs exact SSSS

| kappa | max relerr vs SSSS (100) | median | N* |
|---|---|---|---|
| 1e+06 | 2.00e-02 | 9.00e-03 | 26/100 |
| 1e+08 | 2.44e-04 | 9.87e-05 | 100/100 |
| 1e+10 | 3.05e-06 | 1.40e-06 | 100/100 |

- expected: error falls ~ 1/kappa toward the exact simply supported limit (spring: V_n + kappa w = 0, M_n = 0 throughout)

### First 10 FFFF elastic modes, Argyris finest mesh vs Ritz (omega = sqrt(Lambda))

| mode | omega Ritz | omega FEM | reldiff |
|---|---|---|---|
| 1 | 21.205140 | 21.205493 | 3.32e-05 |
| 2 | 21.235403 | 21.234609 | 7.48e-05 |
| 3 | 48.379721 | 48.379925 | 8.42e-06 |
| 4 | 57.878227 | 57.878075 | 5.25e-06 |
| 5 | 58.296553 | 58.296493 | 2.05e-06 |
| 6 | 73.802224 | 73.802220 | 1.18e-07 |
| 7 | 88.277997 | 88.278603 | 1.37e-05 |
| 8 | 106.347678 | 106.347850 | 3.25e-06 |
| 9 | 120.494828 | 120.495248 | 6.97e-06 |
| 10 | 144.024007 | 144.024507 | 6.94e-06 |

## Verdict: PASS

- machinery (Morley SSSS median extrap < 5e-3; 3 rigid modes everywhere): True
- production (Argyris FFFF: full N* and max relerr < 2e-4 on some mesh; best: 128x80 with N* = 200/200, max relerr 7.48e-05): True
- Protocol-B term (Winkler error falls with kappa, N* full at kappa = 1e+10): True

Total wall time: 335.9 s
