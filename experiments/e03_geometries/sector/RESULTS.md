# E3c -- free disk sector (RESULTS)

theta = 2.0 rad, R = 1.0, nu = 0.33, nrings = 90 (74166 dofs), 800 modes (strict internal N* = 800; operational ladder per the disk-certified smooth-error argument). References: Poisson 0.3863, GOE 0.5307.

- G2 internal N* = 800, rigid [3, 3]; classifier counts {'S': 410, 'A': 390, 'xx': 0}, min quality 0.983, resolved 0

## <r> per reflection class

| class | n levels | <r> |
|---|---|---|
| S | 410 | 0.417 +/- 0.013 |
| A | 390 | 0.426 +/- 0.014 |
| **pooled** | 800 | **0.4213 +/- 0.0095** |
| disk baseline (2 x 390-level windows) | - | 0.3972 +/- 0.0100 |

### Spectral-window trend (thirds)

| window | <r> |
|---|---|
| low | 0.4386 +/- 0.0165 |
| mid | 0.4166 +/- 0.0164 |
| high | 0.4078 +/- 0.0168 |

- sensitivity (strict N* = 800): pooled <r> = 0.4213 +/- 0.0095

## Verdict (preregistered readings)

- pooled <r> = 0.4213 +/- 0.0095 vs disk baseline 0.3972 +/- 0.0100: separation 1.7 sigma (threshold 3)
- per class: S +1.2, A +1.7 sigma

**Reading: AMBIGUOUS (below preregistered threshold)**

Discussion. A weak positive signal (+1.7 sigma pooled, both classes positive) fits the corner-sharpness ordering emerging across the campaign: triangle (three 60-degree corners) 0.489; rectangle (four 90-degree corners) 0.44--0.46; sector (one ~115-degree + two 90-degree corners plus a disk-like arc) 0.421; smooth ellipse 0.373; disk (control) 0.391 -- i.e., the effective coupling appears graded by corner sharpness and prevalence, not by broken separability alone. Relation to the published COMSOL sector study: not a contradiction and not yet a confirmation -- that study used ~10^4 modes and an angle sweep; our single angle (theta = 2.0) at 800 modes shows the same tendency without reaching the preregistered threshold. A decisive comparison needs longer ladders and/or their exact sector angles. The decreasing spectral trend (0.439 -> 0.408) is noted; under the corner picture the apex dominates low-k modes whose evanescent lengths span it, while high-k modes see proportionally more of the uncoupled arc.
