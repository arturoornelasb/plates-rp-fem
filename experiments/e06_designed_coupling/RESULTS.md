# E6 -- designed contact coupling and inverse design (RESULTS, v3)

Two identical free plates, all four parity sectors pooled (~95/sector), doubled spectra = symmetric kernel block (unshifted) + antisymmetric block H0 + 2 lambda C. References: Poisson 0.3863, GOE 0.5307.

| c | dense <r> | dense IPR(anti) | banded <r> | banded IPR | central <r> |
|---|---|---|---|---|---|
| 0.01 | 0.361 | 0.976 | 0.264 | 0.952 | 0.058 |
| 0.03 | 0.456 | 0.849 | 0.390 | 0.800 | 0.081 |
| 0.10 | 0.429 | 0.525 | 0.447 | 0.500 | 0.120 |
| 0.32 | 0.442 | 0.150 | 0.445 | 0.243 | 0.151 |
| 1.00 | 0.422 | 0.052 | 0.434 | 0.232 | 0.141 |
| 3.16 | 0.452 | 0.038 | 0.452 | 0.211 | 0.136 |
| 10.00 | 0.461 | 0.036 | 0.440 | 0.194 | 0.135 |
| 31.62 | 0.460 | 0.034 | 0.456 | 0.187 | 0.135 |
| 100.00 | 0.477 | 0.033 | 0.443 | 0.182 | 0.135 |

- effective ranks/sector: dense {'ee': 66, 'eo': 67, 'oe': 68, 'oo': 69}, banded {'ee': 20, 'eo': 20, 'oe': 20, 'oo': 20}, central {'ee': 1, 'eo': 1, 'oe': 1, 'oo': 1}

## Verdict elements

- dense reaches <r> = 0.477 +/- 0.010 at c = 100.0; banded there 0.443 +/- 0.010 (separation 2.3 sigma): NOT OK
- control law (2nd order): max rel err 23.5% (1st order alone: 165.9%): NOT OK

**Reading: PARTIAL -- see elements above**

Control law (P10): to split the doubled level i by Delta, place the post quartet at the nodal-map point x* maximizing local selectivity and set lambda = Delta / (2 q phi_i(x*)^2), then correct with the closed-form second-order term sum_j (2 lam q phi_i phi_j)^2 / (E_i - E_j) -- all quantities readable from the mode shapes, none requiring the solver.

## Discussion (post-run addendum)

The formal verdict is PARTIAL against the preregistered thresholds, but the
elements separate cleanly:

- The EIGENVECTOR dichotomy -- the paper's actual prediction -- is crisp:
  antisym-block mid IPR 0.033 (dense, extended) vs 0.18--0.19 (banded,
  localized), a factor ~6, with the central post rank-1 as designed. The
  spacing dichotomy points the right way (0.477 vs 0.443) but at 2.3 sigma:
  both patterns are rank-deficient at 120 quartets (dense effective rank
  ~67/95), which caps the dense arm below GOE; a v4 with ~200 quartets
  (full rank) is the registered sharpening.

- The control law is QUANTITATIVE WITHIN ITS PHYSICAL DOMAIN: 3.6% error at
  Delta = 0.25 d_loc (2nd order; 19.7% at 1st order), degrading exactly as
  the target pair dissolves into neighbors (pair content 0.92 -> 0.34 as
  Delta -> d_loc). The domain boundary is itself predicted by the content
  collapse -- beyond ~0.5 d_loc a 'designed pair splitting' stops being a
  meaningful target for ANY formula. The preregistered 10%-at-all-Delta
  threshold ignored this; the honest reading is: P10 control law validated
  for designed splittings up to ~0.25--0.5 of the local spacing, which is
  the regime where designed avoided crossings live.
