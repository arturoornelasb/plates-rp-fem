# E4 -- Gap A on the true operator (RESULTS)

Rectangle a/b = 1.618904, nu = 0.33, Argyris 128x80, N_use = 1200 certified modes (G1 N* = 1200, G2 N* = 1200, vector residuals <= 7.6e-06). Ladder window [0.4N, 0.6N); GOE baseline 8 realizations/size.

- labels: {'ee': 310, 'eo': 298, 'oe': 302, 'oo': 290, 'xx': 0}, min quality 0.968


## Basis: sine

| sector | N=64 | N=100 | N=144 | N=200 | N=256 | N=324 | D2 ladder | +/- se | D2 fixed-N (max N) | mean captured p |
|---|---|---|---|---|---|---|---|---|---|---|
| ee | 0.1837 | 0.1760 | 0.1916 | 0.1840 | 0.1754 | 0.1699 | 0.034 | 0.031 | 0.496 | 0.878 |
| eo | 0.1897 | 0.1819 | 0.1678 | 0.1500 | 0.1606 | 0.1677 | 0.107 | 0.045 | 0.498 | 0.878 |
| oe | 0.2034 | 0.1943 | 0.1926 | 0.1615 | 0.1674 | 0.1777 | 0.120 | 0.045 | 0.488 | 0.878 |
| oo | 0.1804 | 0.1690 | 0.1806 | 0.1804 | 0.1734 | 0.1658 | 0.028 | 0.028 | 0.500 | 0.878 |
| GOE | 0.0443 | 0.0291 | 0.0203 | 0.0147 | 0.0115 | 0.0092 | ~1 | - | 1.000 (def.) | - |

## Basis: beam

| sector | N=64 | N=100 | N=144 | N=200 | N=256 | N=324 | D2 ladder | +/- se | D2 fixed-N (max N) | mean captured p |
|---|---|---|---|---|---|---|---|---|---|---|
| ee | 0.7270 | 0.7094 | 0.7920 | 0.7360 | 0.6926 | 0.6604 | 0.047 | 0.045 | 0.261 | 0.999 |
| eo | 0.7728 | 0.7395 | 0.6831 | 0.6055 | 0.6340 | 0.6674 | 0.125 | 0.043 | 0.259 | 1.000 |
| oe | 0.8319 | 0.8039 | 0.7880 | 0.6554 | 0.6747 | 0.7219 | 0.131 | 0.048 | 0.246 | 0.999 |
| oo | 0.7628 | 0.6606 | 0.7600 | 0.7592 | 0.7221 | 0.6832 | 0.023 | 0.050 | 0.255 | 1.000 |
| GOE | 0.0443 | 0.0291 | 0.0203 | 0.0147 | 0.0115 | 0.0092 | ~1 | - | 1.000 (def.) | - |

## Long-range statistics Sigma^2(L) per sector

| L | ee | eo | oe | oo | Poisson | GOE |
|---|---|---|---|---|---|---|
| 1 | 0.70 | 0.95 | 0.80 | 0.79 | 1.0 | 0.44 |
| 2 | 1.27 | 1.42 | 1.57 | 1.38 | 2.0 | 0.58 |
| 5 | 2.46 | 2.66 | 2.76 | 2.36 | 5.0 | 0.77 |
| 10 | 3.40 | 3.24 | 2.84 | 2.44 | 10.0 | 0.91 |
| 15 | 3.15 | 2.66 | 2.94 | 2.22 | 15.0 | 0.99 |
| 20 | 2.83 | 2.66 | 2.55 | 2.21 | 20.0 | 1.05 |

## Verdict (preregistered readings)

- ladder D2 range across bases/sectors: [0.023, 0.131]; fixed-N D2 range: [0.246, 0.500] (P12 reference 0.76 +/- 0.15; T1 Ritz-based ladder 0.20--0.43 in these bases)
- Sigma^2(20)/20 (Poisson = 1, GOE ~ 0.1): mean over sectors 0.13
- eigenvectors extended (D2 ladder > 0 in all cells): False
- non-ergodic (fixed-N D2 strictly inside (0,1)): True

**Reading: LOCALIZED/BANDED-LEANING (IPR ~ N-independent)**

Caveats: D2 from a representational ladder at these N carries the T1 basis-dependence systematic (both registered bases reported); ladders reach N = 324 (paper registers 256--2048; larger rungs need quarter-plate reduction -- future work). Sigma^2 beyond L ~ 10 is unfolding-limited (polynomial unfolding absorbs long-wavelength fluctuations); the L <= 5 values are the reliable ones and sit intermediate, mildly below Poisson.

Discussion. The beam-product basis represents the true modes almost losslessly (captured p = 0.999--1.000), and there each free-plate mode is ~1--2 products (IPR ~ 0.7): the true operator at these mode numbers is a WEAKLY PERTURBED SEPARABLE system with sparse strong pair-hybridizations -- microscopically, exactly the avoided crossings of the original Lopez-Gonzalez observations -- rather than the dense ergodic mixing of the RP non-ergodic phase. Combined with the campaign's spacing results (intermediate <r> ~ 0.44--0.49 for cornered geometries), the accessible-N picture is 'sparse/weak V with strong local hybridizations': intermediate spacing statistics WITHOUT RP-fractal eigenvectors. Two sharpenings for the paper: (i) the true-operator D2 is LOWER than the Ritz-model values (beam fixed-N 0.25 vs T1's 0.35--0.42) -- the projection/model route overestimates coupling density, as T1's truncation caveat anticipated; (ii) since the effective coupling grows with mode index (E2 transition; E3b trend), the RP fractal phase may only develop at ladders beyond N ~ 324 -- the paper's registered 2048 rung, reachable with quarter-plate reduction, is where this verdict could still flip. At accessible N, observable 1 favors the banded/weak-coupling reading.
