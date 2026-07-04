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
