# E9 -- Gap A at the registered N = 2048 rung (preregistered)

The paper registers ladder sizes N = 256--2048 ("COMSOL-realistic mode
counts"); E4/E8 reached N <= 324 on true-operator modes and found the sparse
regime for the rectangle, with the explicit caveat that the RP fractal phase
may develop at larger N (effective coupling grows with mode index). This
experiment executes the full registered ladder in the T1 sense: the ladder
rung IS the truncated operator (principal submatrix of the sector Ritz
matrix in the frequency-proxy order) -- "the truncated operator IS the
system at size N, as in RP numerics" (T1).

## Why this is now possible

Two measured obstructions are removed. (i) The ~1e-4 floor of the
quadrature-assembled Ritz matrices: platefem.ritz_exact assembles the SAME
matrices from closed-form integer arithmetic (validated to 3e-11 against
quadrature at NLEG = 48). (ii) The residual high-NLEG corruption is now
attributed by measurement to LAPACK backward error against ||K|| ~ n^7 of
the FULL matrix -- which the ladder protocol never touches: rung
submatrices have physical norms (||K_2048|| ~ 3e10, eigh error ~ 3e-6, ten
orders below the local spacing).

## Protocol

Campaign rectangle, sectors ee and oo, NLEG = 160/axis (sector dim 6400);
ladder N in {256, 512, 1024, 2048}; window [0.4N, 0.6N); representation
bases = T1's SS sine products and FF beam products (Loewdin-orthonormalized,
overlaps by high-order quadrature); IPR and D_q (q = 1.5, 2, 4); empirical
GOE baselines at identical N.

## Preregistered readings

- RP PHASE DEVELOPS: IPR resumes falling with N beyond ~324 (ladder D2
  rising toward the P12 reference 0.76 +/- 0.15) with flat D_q
  (D_1.5 - D_4 < 0.15).
- SPARSE REGIME PERSISTS: IPR stays flat through N = 2048 (D2 ~ 0) -- the
  banded/sparse reading of E4 extends to the full registered range.
- PBRM-LIKE: D2 intermediate but D_1.5 - D_4 > 0.15.
Caveat (T1's): mid-spectrum modes of a truncated Ritz matrix are not
converged plate modes; this is inherent to the registered ladder protocol
and is exactly what E4/E8 complement with true-operator windows.
