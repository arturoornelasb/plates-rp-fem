# E7 -- completeness: Protocol A mixed config + material independence

## (a) Protocol A config (ii): x-edges free, y-edges SS

- pooled <r> = **0.3973 +/- 0.0103** (FFFF at this ladder: ~0.442; SSSS baseline: ~0.391); per sector ee 0.376, eo 0.399, oe 0.382, oo 0.434
- rigid modes detected: 0 (expected 0: w = 0 on the y-edges blocks all rigid motions)

## (b) FFFF <r>(nu) -- material independence

| nu | ee | eo | oe | oo | pooled |
|---|---|---|---|---|---|
| 0.25 | 0.425 | 0.378 | 0.457 | 0.421 | **0.4208 +/- 0.0094** |
| 0.3 | 0.439 | 0.389 | 0.464 | 0.451 | **0.4361 +/- 0.0095** |
| 0.33 | 0.448 | 0.393 | 0.465 | 0.461 | **0.4422 +/- 0.0096** |
| 0.36 | 0.458 | 0.401 | 0.472 | 0.466 | **0.4497 +/- 0.0096** |

- spread across nu: 0.0290 (typical sem 0.0095) -> nu-DEPENDENT (check)

## Discussion (post-run addendum)

(a) The preregistered expectation ("between SSSS and FFFF") was wrong in an
instructive way: the mixed configuration (free x-edges, SS y-edges) is the
classical LEVY-SEPARABLE plate -- w = X(x) sin(n pi y / b) solves the free
x-edge conditions exactly, so the spectrum is a union of independent 1D
quartic pencils per n. Under the campaign's adaptedness hierarchy this
configuration belongs to the exactly-reduced tier and must be Poisson; the
measured 0.3973 +/- 0.0103 (vs the SSSS finite-size baseline 0.391) confirms
it. One free edge PAIR does not suffice: intermediate statistics require
free edges that break every available product structure.

(b) The nu sweep reproduces the T1 Ritz bonus value-by-value (T1: 0.4209,
0.4376, 0.4412, 0.4470; FEM: 0.4208, 0.4361, 0.4422, 0.4497 for
nu = 0.25..0.36) -- two independent methods agreeing to ~0.003. The
monotone +0.029 end-to-end shift (~3 sigma pooled here, vs T1's 1.7 sigma)
establishes that the effective free-edge coupling is MATERIAL-TUNABLE
through the nu-dependence of the boundary conditions, supporting the
Prediction-5 recast (lambda(nu) as a control knob) rather than strict
material independence of <r>.
