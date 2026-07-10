# E18b -- faithful truncated protocol (RESULTS)

Dubiner p <= 58 (1770 functions), sector counts {'A1': 310, 'A2': 280, 'E': 1180}; deg<=4 reproduction 1.5e-14.


## sector A1 (polys 310; Ritz sanity 2.2e-01)

| N | trunc IPR | (true IPR, E18) |
|---|---|---|
| 128 | 0.1081 | 0.1262 |
| 256 | 0.0553 | 0.1061 |

- slope_trunc = -0.967 (0.176); slope_true (E18) = -0.250 (0.116); ratio trunc/true = +3.87
- **INTERMEDIATE**

## sector A2 (polys 280; Ritz sanity 1.7e-02)

| N | trunc IPR | (true IPR, E18) |
|---|---|---|
| 128 | 0.0902 | 0.1479 |
| 256 | 0.0422 | 0.1335 |

- slope_trunc = -1.096 (0.181); slope_true (E18) = -0.147 (0.183); ratio trunc/true = +7.46
- **PROTOCOL ARTIFACT**

## sector E (polys 1180; Ritz sanity 3.7e-01)

| N | trunc IPR | (true IPR, E18) |
|---|---|---|
| 128 | 0.1401 | 0.0805 |
| 256 | 0.0826 | 0.0788 |
| 512 | 0.0483 | 0.0698 |
| 1024 | 0.0210 | 0.0631 |

- slope_trunc = -0.961 (0.048); slope_true (E18) = -0.138 (0.029); ratio trunc/true = +6.96
- **PROTOCOL ARTIFACT**

**Reading: PROTOCOL ARTIFACT on the unadapted tier; combined with E14/E17, Gap A closes negatively for ALL tiers at accessible N (the triangle's weak true slope notwithstanding -- see ratio)**

Wall time: 254.9 s.

## Post-run notes (2026-07-09)

- **Headline numbers:** the faithful truncated protocol produces
  near-unit slopes in EVERY sector (D2_trunc ~ 0.96-1.10 -- near-ergodic
  scaling), 4-7x the true operator's genuine slopes (ratios 3.9 / 7.5 /
  7.0). E-sector is the cleanest: trunc starts ABOVE the true level at
  N = 128 (0.140 vs 0.081) and crashes far below by N = 1024 (0.021 vs
  0.063) -- the truncation-scaling signature on a second geometry.
- **Ritz-sanity check bug (does not affect the ladders):** the lowest-20
  comparison is rigid-shifted in A1/E -- the polynomial space contains
  the constant (A1) and the linear doublet (E) that the free solve's
  rigid split excluded; the window rungs [0.4N, 0.6N) sit far above.
  A2 (no rigid content) shows ~1.7e-2: genuine slow polynomial
  convergence to CORNER-SINGULAR free-plate modes -- which is precisely
  the artifact's mechanism (smooth truncated trial spaces misrepresent
  corner/boundary structure and spread modes across many reference
  functions).
- **Refined Gap-A statement for the unadapted tier:** a weak, genuine
  true-operator delocalization exists (E18: E slope -0.138 +/- 0.029,
  4.8 sigma) and truncated-ladder protocols exaggerate it ~7-fold to
  near-ergodic scaling. A1's INTERMEDIATE is the flat-bar (0.15) sitting
  below its noisy true slope, not a divergent story (ratio 3.9x).
- Superellipse contingency: NOT triggered (no GENUINE SCALING sector).
