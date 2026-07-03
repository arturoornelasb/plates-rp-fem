# E3a -- free equilateral triangle (RESULTS)

Side L = 1.0, nu = 0.33, Argyris refine-7 (C3v-symmetric mesh), 1000 modes computed, N_use = 1000, lowest 10/sector dropped after E-dedup. References: Poisson 0.3863, GOE 0.5307.

## Accuracy gates

- G1 (kappa=1e10 vs exact Lame^2): N* = 1000, max relerr 2.97e-05
- G2 (free, refine 7 vs 6): internal N* = 1000, rigid [3, 3]

## Classifier / degeneracy diagnostics

- free: counts {'A1': 178, 'A2': 156, 'E': 666, 'xx': 0}, min quality 0.790, E-doublets 333 pairs + 0 unpaired, max pair splitting 7.5e-06, rigid 3
- ss: counts {'A1': 179, 'A2': 155, 'E': 666, 'xx': 0}, min quality 0.707, E-doublets 333 pairs + 0 unpaired, max pair splitting 5.7e-07, rigid 0

## <r> per sector (E deduplicated)

| case | A1 | A2 | E | pooled | low half | high half | n_r |
|---|---|---|---|---|---|---|---|
| free | 0.490+/-0.017 | 0.437+/-0.022 | 0.512+/-0.014 | **0.4889 +/- 0.0099** | 0.522 | 0.458 | 631 |
| ss | 0.381+/-0.024 | 0.375+/-0.025 | 0.373+/-0.020 | **0.3753 +/- 0.0134** | 0.421 | 0.333 | 631 |

## Verdict (preregistered readings)

- free pooled <r> = 0.4889 +/- 0.0099; SS baseline (same protocol) = 0.3753 +/- 0.0134
- separation free vs SS baseline: 6.8 sigma (preregistered threshold: 3)
- per-sector separations: A1 3.8 sigma, A2 1.9 sigma, E 5.6 sigma

**Reading: SUPPORTS the geometry prediction (intermediate statistics in the non-separable free triangle)**
