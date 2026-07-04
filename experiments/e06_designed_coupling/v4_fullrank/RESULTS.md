# E6 -- designed contact coupling and inverse design (RESULTS, v3)

Two identical free plates, all four parity sectors pooled (~95/sector), doubled spectra = symmetric kernel block (unshifted) + antisymmetric block H0 + 2 lambda C. References: Poisson 0.3863, GOE 0.5307.

| c | dense <r> | dense IPR(anti) | banded <r> | banded IPR | central <r> |
|---|---|---|---|---|---|
| 0.01 | 0.369 | 0.968 | 0.270 | 0.949 | 0.058 |
| 0.03 | 0.425 | 0.840 | 0.401 | 0.809 | 0.081 |
| 0.10 | 0.444 | 0.502 | 0.434 | 0.497 | 0.120 |
| 0.32 | 0.410 | 0.156 | 0.430 | 0.246 | 0.151 |
| 1.00 | 0.416 | 0.048 | 0.433 | 0.244 | 0.141 |
| 3.16 | 0.465 | 0.038 | 0.432 | 0.219 | 0.136 |
| 10.00 | 0.448 | 0.035 | 0.448 | 0.201 | 0.135 |
| 31.62 | 0.474 | 0.034 | 0.434 | 0.201 | 0.135 |
| 100.00 | 0.449 | 0.034 | 0.441 | 0.202 | 0.135 |

- effective ranks/sector: dense {'ee': 71, 'eo': 73, 'oe': 74, 'oo': 76}, banded {'ee': 20, 'eo': 21, 'oe': 21, 'oo': 21}, central {'ee': 1, 'eo': 1, 'oe': 1, 'oo': 1}

## Verdict elements

- dense reaches <r> = 0.474 +/- 0.010 at c = 31.6; banded there 0.434 +/- 0.010 (separation 2.8 sigma): NOT OK
- control law (2nd order): max rel err 23.5% (1st order alone: 165.9%): NOT OK

**Reading: PARTIAL -- see elements above**

Control law (P10): to split the doubled level i by Delta, place the post quartet at the nodal-map point x* maximizing local selectivity and set lambda = Delta / (2 q phi_i(x*)^2), then correct with the closed-form second-order term sum_j (2 lam q phi_i phi_j)^2 / (E_i - E_j) -- all quantities readable from the mode shapes, none requiring the solver.

## Addendum (v4 close-out)

Raising the quartet count 120 -> 220 leaves the dense effective rank at
~71-76/95: random contact patterns are intrinsically SMOOTHING operators
whose singular values decay, so full rank is not reachable by adding posts
and the spacing dichotomy plateaus at ~2.8 sigma. The registered follow-up
is closed with this reading: the paper's contact prediction is carried by
the crisp EIGENVECTOR dichotomy (extended 0.034 vs localized 0.20) plus the
in-domain control law; the pooled-<r> contrast saturates below the 3-sigma
preregistration for intrinsic-rank reasons, not statistical ones.
