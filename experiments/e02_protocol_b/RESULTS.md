# E2 -- Protocol B: boundary-controlled transition (RESULTS)

Rectangle a/b = 1.618904, nu = 0.33, Argyris 96x60, 400 modes computed, N_use = 400 (accuracy gates), lowest 10/sector dropped. References: Poisson 0.3863, GOE 0.5307.

## Accuracy gates

- G1 (kappa=1e10 vs exact SSSS): N* = 400, max relerr 2.68e-05
- G2 (kappa=0, 96x60 vs 128x80): internal N* = 400, rigid [3, 3]
- G3 (classifier vs Ritz sectors, lowest modes): ee 9.9e-05, eo 7.9e-06, oe 8.0e-06, oo 2.6e-06; ambiguous 0

## Finite-size baselines (identical protocol, ~100 levels/sector)

- exact SSSS (separable -> Poisson limit): pooled <r> = **0.4167 +/- 0.0155** (asymptotic 0.3863; the excess is pure finite-size, measured on the EXACT spectrum)
- Legendre-Ritz FFFF (106/sector): pooled <r> = **0.4535 +/- 0.0134** (independent method, free-plate reference)

## <r> vs kappa (pooled over sectors; windows are index halves of each sector ladder)

| kappa | ee | eo | oe | oo | pooled | low half | high half | n_r | min quality |
|---|---|---|---|---|---|---|---|---|---|
| 0e+00 | 0.467+/-0.027 | 0.446+/-0.028 | 0.476+/-0.027 | 0.448+/-0.027 | **0.4595 +/- 0.0136** | 0.423 | 0.492 | 352 | 1.00 |
| 1e+02 | 0.468+/-0.027 | 0.446+/-0.028 | 0.479+/-0.027 | 0.448+/-0.027 | **0.4605 +/- 0.0136** | 0.428 | 0.493 | 352 | 0.98 |
| 3e+02 | 0.468+/-0.027 | 0.446+/-0.028 | 0.479+/-0.027 | 0.447+/-0.027 | **0.4604 +/- 0.0136** | 0.428 | 0.493 | 352 | 1.00 |
| 1e+03 | 0.466+/-0.027 | 0.446+/-0.028 | 0.478+/-0.027 | 0.447+/-0.027 | **0.4599 +/- 0.0136** | 0.427 | 0.493 | 352 | 1.00 |
| 3e+03 | 0.464+/-0.027 | 0.447+/-0.028 | 0.478+/-0.027 | 0.447+/-0.027 | **0.4592 +/- 0.0137** | 0.425 | 0.493 | 352 | 1.00 |
| 1e+04 | 0.458+/-0.026 | 0.448+/-0.028 | 0.471+/-0.028 | 0.449+/-0.028 | **0.4569 +/- 0.0137** | 0.420 | 0.493 | 352 | 1.00 |
| 3e+04 | 0.448+/-0.025 | 0.436+/-0.029 | 0.478+/-0.029 | 0.473+/-0.029 | **0.4586 +/- 0.0140** | 0.425 | 0.491 | 352 | 1.00 |
| 1e+05 | 0.411+/-0.026 | 0.438+/-0.029 | 0.454+/-0.030 | 0.421+/-0.028 | **0.4309 +/- 0.0141** | 0.407 | 0.455 | 352 | 0.99 |
| 3e+05 | 0.457+/-0.027 | 0.421+/-0.028 | 0.409+/-0.029 | 0.451+/-0.030 | **0.4348 +/- 0.0144** | 0.422 | 0.443 | 352 | 0.99 |
| 1e+06 | 0.429+/-0.028 | 0.374+/-0.029 | 0.362+/-0.033 | 0.416+/-0.032 | **0.3957 +/- 0.0153** | 0.391 | 0.406 | 352 | 0.97 |
| 3e+06 | 0.477+/-0.028 | 0.401+/-0.030 | 0.438+/-0.028 | 0.427+/-0.033 | **0.4366 +/- 0.0150** | 0.404 | 0.476 | 352 | 1.00 |
| 1e+07 | 0.455+/-0.030 | 0.380+/-0.030 | 0.433+/-0.029 | 0.399+/-0.034 | **0.4179 +/- 0.0154** | 0.407 | 0.434 | 352 | 1.00 |
| 1e+08 | 0.454+/-0.029 | 0.374+/-0.031 | 0.422+/-0.029 | 0.409+/-0.034 | **0.4159 +/- 0.0154** | 0.408 | 0.429 | 352 | 0.98 |
| 1e+10 | 0.458+/-0.030 | 0.374+/-0.031 | 0.419+/-0.029 | 0.411+/-0.034 | **0.4167 +/- 0.0155** | 0.408 | 0.430 | 352 | 0.99 |

## Verdict (preregistered readings, finite-size-corrected)

- SS end: <r>(kappa=1e10) = 0.4167 +/- 0.0155 vs exact-SSSS baseline 0.4167 +/- 0.0155: consistent
- free end: <r>(kappa=0) = 0.4595 +/- 0.0136 vs Ritz baseline 0.4535 +/- 0.0134: consistent
- transition amplitude: free end sits 2.1 sigma above the separable (exact SSSS) baseline
- monotonicity (3-sigma violations along the kappa grid): 0
- intermediate already at kappa >= 1e8 (challenge condition): False

**Reading: CHECK -- see flags above**

Caveats: ~100 levels/sector is a short ladder -- the finite-size excess of the separable baseline over asymptotic Poisson (0.417 vs 0.3863) is measured, not assumed; a longer-ladder rerun (800+ modes) is the natural sharpening step.
