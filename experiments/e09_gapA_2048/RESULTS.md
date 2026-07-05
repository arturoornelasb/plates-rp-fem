# E9 -- Gap A at the registered N = 2048 rung (RESULTS)

T1 truncated-operator ladder, exact assembly, NLEG = 160 (sector dim 6400), window [0.4N, 0.6N). P12 reference D2 = 0.76 +/- 0.15.


## sector ee

| basis | IPR N=256 | IPR N=512 | IPR N=1024 | IPR N=2048 | D2 ladder | +/- | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|
| sin | 0.1233 | 0.0946 | 0.0667 | 0.0522 | 0.422 | 0.019 | 0.069 | 0.883 |
| beam | 0.2432 | 0.1784 | 0.1241 | 0.0947 | 0.461 | 0.017 | 0.087 | 0.995 |
| GOE | 0.01145 | 0.00582 | 0.00292 | 0.00146 | ~1 | - | - | - |

## sector oo

| basis | IPR N=256 | IPR N=512 | IPR N=1024 | IPR N=2048 | D2 ladder | +/- | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|
| sin | 0.1161 | 0.0884 | 0.0597 | 0.0458 | 0.460 | 0.025 | 0.069 | 0.877 |
| beam | 0.2313 | 0.1767 | 0.1141 | 0.0843 | 0.500 | 0.032 | 0.100 | 0.985 |
| GOE | 0.01145 | 0.00582 | 0.00292 | 0.00146 | ~1 | - | - | - |

## Verdict (preregistered readings)

- ladder D2 range: [0.422, 0.500] (P12 0.76 +/- 0.15); D_1.5 - D_4 range: [0.069, 0.100] (RP < 0.15 < PBRM)

**Reading: RP PHASE DEVELOPS (fractal, flat Dq) at large N**

## Discussion (post-run addendum)

The registered rung decides in favor of the RP reading, with precision:
scaling is genuine (three octaves, both sectors, both registered bases,
sine captured 0.88 / beam 0.99), and the registered rival separator is
FLAT (0.069--0.100, vs PBRM's 0.20--0.28 calibration) -- the power-law
banded alternative is excluded at these N. Two qualifications. (1) D2 =
0.46 +/- 0.03 sits below the P12 preregistered 0.76 +/- 0.15: the
non-ergodic fractal PHASE is confirmed; the preregistered dimension value
is not. (2) The truncated-operator ladder (this experiment; the standard
RP-numerics protocol, and the paper's registered one) shows scaling where
the true-operator windows at N <= 324 (E4/E8) showed flat IPR -- exactly
the protocol difference T1 flagged. Physically consistent: the effective
coupling grows with mode index, and the truncation ladder samples deeper
windows as N grows. The reconciliation measurement -- true-operator
eigenvectors at N ~ 1024-2048 -- requires the C0-interior-penalty
discretization (Argyris ceiling measured at h ~ 0.008; skfem Lagrange
elements lack second derivatives, so this needs custom reference-hessian
assembly) and remains the registered follow-up.

## Reinterpretation (E14, 2026-07-04)

The E14 reconciliation measurement decides the protocol question: on TRUE
operator eigenvectors (C0-IP P4, 206k dofs, gates G2 = 1200/1200 vs
certified Argyris at 1.2e-4), the window IPR stays FLAT at the E4 level
(sine 0.174 sector-mean at N = 512; beam 0.63-0.72, i.e. modes remain ~1-2
beam products) while this experiment's truncated ladder falls to 0.095 at
the same N. The scaling reported above is therefore a property of the
TRUNCATION PROTOCOL, not of the plate: diagonalizing frequency-ordered
principal submatrices progressively mixes the representation in a way the
physical operator does not. Consequence for the program: any Gap A executed
by truncated-ladder numerics (including COMSOL pipelines following the
registered P12 recipe) would report an RP-like phase as an artifact; Gap A
must be stated on true-operator windows, where the plate is sparse /
quasi-separable through at least N = 512 per sector.
