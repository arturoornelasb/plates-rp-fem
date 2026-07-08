# i01 -- quarter-plate C0-IP instrument (SMOKE RESULTS)

## GATE A: merged quarter sectors vs same-mesh full SSSS

- merged 4 x 60 quarter modes vs full-plate SSSS at the same h: **N* = 239/240, max relerr 1.81e-03** (required: >= 188 = the full instrument's own exact-accuracy range)

### informational: per-sector vs EXACT (mesh error at smoke h)

| sector | N* | max relerr |
|---|---|---|
| ee | 51/60 | 5.29e-03 |
| eo | 50/60 | 6.05e-03 |
| oe | 50/60 | 6.32e-03 |
| oo | 47/60 | 8.03e-03 |

## GATE B: quarter (free outer) vs full-plate C0-IP sectors

| sector | rigid (expected) | N* | max relerr |
|---|---|---|---|
| ee | 1 (True) | 48/48 | 1.68e-04 |
| eo | 1 (True) | 48/48 | 1.93e-04 |
| oe | 1 (True) | 48/48 | 1.73e-04 |
| oo | 0 (True) | 48/48 | 7.47e-04 |

**GATE A: PASS; GATE B: PASS -> PASS -- promoted to production instrument for E18 and future Gap-A rungs.**

Wall time: 53.0 s.
