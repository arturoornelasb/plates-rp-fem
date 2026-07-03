# E2 -- Protocol B: boundary-controlled transition (RESULTS)

Rectangle a/b = 1.618904, nu = 0.33, Argyris 128x80, 800 modes computed, N_use = 800 (accuracy gates), lowest 10/sector dropped. References: Poisson 0.3863, GOE 0.5307.

## Accuracy gates

- G1 (kappa=1e10 vs exact SSSS): N* = 800, max relerr 3.30e-04
- G2 (kappa=0, 96x60 vs 128x80): internal N* = 800, rigid [3, 3]
- G3 (classifier vs Ritz sectors, lowest modes): ee 7.5e-05, eo 8.0e-06, oe 8.6e-06, oo 3.3e-05; ambiguous 0

## Finite-size baselines (identical protocol, ~200 levels/sector)

- exact SSSS (separable -> Poisson limit): pooled <r> = **0.3908 +/- 0.0104** (asymptotic 0.3863; the excess is pure finite-size, measured on the EXACT spectrum)
- Legendre-Ritz FFFF (162/sector): pooled <r> = **0.4434 +/- 0.0107** (independent method, free-plate reference)

## <r> vs kappa (pooled over sectors; windows are index halves of each sector ladder)

| kappa | ee | eo | oe | oo | pooled | low half | high half | n_r | min quality |
|---|---|---|---|---|---|---|---|---|---|
| 0e+00 | 0.448+/-0.018 | 0.393+/-0.019 | 0.465+/-0.018 | 0.461+/-0.020 | **0.4422 +/- 0.0096** | 0.455 | 0.430 | 752 | 0.96 |
| 1e+04 | 0.441+/-0.018 | 0.398+/-0.019 | 0.458+/-0.019 | 0.460+/-0.020 | **0.4392 +/- 0.0096** | 0.454 | 0.424 | 752 | 0.97 |
| 3e+05 | 0.463+/-0.019 | 0.398+/-0.019 | 0.440+/-0.020 | 0.430+/-0.021 | **0.4334 +/- 0.0099** | 0.435 | 0.435 | 752 | 0.96 |
| 1e+06 | 0.414+/-0.019 | 0.428+/-0.019 | 0.402+/-0.022 | 0.386+/-0.021 | **0.4080 +/- 0.0101** | 0.391 | 0.426 | 752 | 0.98 |
| 3e+06 | 0.410+/-0.020 | 0.399+/-0.020 | 0.409+/-0.020 | 0.401+/-0.022 | **0.4052 +/- 0.0102** | 0.429 | 0.378 | 752 | 0.99 |
| 1e+08 | 0.411+/-0.020 | 0.354+/-0.021 | 0.416+/-0.020 | 0.369+/-0.022 | **0.3881 +/- 0.0104** | 0.406 | 0.368 | 752 | 0.98 |
| 1e+10 | 0.418+/-0.020 | 0.354+/-0.021 | 0.419+/-0.020 | 0.369+/-0.022 | **0.3908 +/- 0.0104** | 0.407 | 0.372 | 752 | 0.96 |

## Verdict (preregistered readings, finite-size-corrected)

- SS end: <r>(kappa=1e10) = 0.3908 +/- 0.0104 vs exact-SSSS baseline 0.3908 +/- 0.0104: consistent
- free end: <r>(kappa=0) = 0.4422 +/- 0.0096 vs Ritz baseline 0.4434 +/- 0.0107: consistent
- transition amplitude: free end sits 3.6 sigma above the separable (exact SSSS) baseline
- monotonicity (3-sigma violations along the kappa grid): 0
- intermediate already at kappa >= 1e8 (challenge condition): False

**Reading: SUPPORTS the hypothesis**

Notes: the finite-size excess of the separable baseline over asymptotic Poisson (0.391 vs 0.3863) is measured on the exact spectrum, not assumed. v2 doubles v1's ladder to ~200 levels/sector: the separable baseline dropped toward Poisson exactly as expected (0.417 -> 0.391) and the transition amplitude cleared the preregistered 3 sigma.
