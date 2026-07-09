# E18 -- unadapted-tier protocol scrutiny: the triangle (RESULTS)

C0-IP P4, production n = 190 (~289941 dofs), gates: {'free_two_mesh': 1056, 'ss_two_mesh': 1977, 'ss_exact_anchor': 1999, 'free_vs_e3a_argyris': 1000, 'n_e3a': 1000}


## sector A1 (basis 346, free 188)

| N | true IPR | trunc IPR |
|---|---|---|
| 128 | 0.1262 | 1.0000 |
| 256 | 0.1061 | 1.0000 |

- slope_true = -0.250 (0.116); slope_trunc = -0.000 (0.000) -> D2_true = 0.250, D2_trunc = 0.000
- **INTERMEDIATE**

## sector A2 (basis 313, free 166)

| N | true IPR | trunc IPR |
|---|---|---|
| 128 | 0.1479 | 1.0000 |
| 256 | 0.1335 | 1.0000 |

- slope_true = -0.147 (0.183); slope_trunc = -0.000 (0.000) -> D2_true = 0.147, D2_trunc = 0.000
- **INTERMEDIATE**

## sector E (basis 1318, free 702)

| N | true IPR | trunc IPR |
|---|---|---|
| 128 | 0.0805 | 0.8968 |
| 256 | 0.0788 | 0.9256 |
| 512 | 0.0698 | 0.8530 |
| 1024 | 0.0631 | 0.8415 |

- slope_true = -0.138 (0.029); slope_trunc = -0.047 (0.014) -> D2_true = 0.138, D2_trunc = 0.047
- **INTERMEDIATE**

**Reading: MIXED/INTERMEDIATE -- see per-sector table**

Wall time (analysis): 92.2 s.

## Post-run correction (2026-07-09): the truncated column is INVALID by construction; the true ladders stand

- **Truncated column degenerate:** projecting the free stiffness onto the
  zero-extended SS eigenvectors is trivially diagonal in a conforming
  discretization (they are eigenvectors of the SAME interior operator;
  the free/SS difference lives in the boundary dofs, not the matrix) --
  hence IPR = 1.0000 in A1/A2; the E-sector 0.84-0.93 is arbitrary
  rotation inside exact doublets. This is NOT the E9 protocol (which
  diagonalizes nested rungs of an INCOMPLETE polynomial trial space and
  represents them in the reference basis). Constructor error, mine;
  the INTERMEDIATE verdicts above compare against a broken comparator
  and carry no weight.
- **What stands:** the gates (free two-mesh 1056/1056 FULL coverage;
  SS vs exact Lame^2 1999/2000; free vs certified Argyris 1000/1000)
  and the TRUE-operator ladders -- the first real departure from the
  rectangle: E-sector window IPR falls 0.081 -> 0.063 over
  N = 128 -> 1024 (slope -0.138 +/- 0.029, i.e. 4.8 sigma from flat;
  D2_true ~ 0.14), A-sectors consistent (slopes -0.25/-0.15, noisier).
  The rectangle at the same protocol was flat (-0.01). The unadapted
  tier's true operator DOES delocalize slowly.
- **Registered: E18b**, the faithful E9 analog -- free-plate Ritz in
  nested Koornwinder (total-degree-ordered) polynomial spaces via their
  P4 interpolants (M-orthonormalized through the FEM mass; C0-IP form on
  smooth interpolants = the bending form to interpolation error), rung
  eigenvectors sector-classified with the validated C3v machinery, then
  represented in the SS eigenbasis; same windows, same frozen dichotomy
  against the (valid) true ladders above.
