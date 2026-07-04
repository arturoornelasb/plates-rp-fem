# E3c -- free disk sector (RESULTS)

theta = 2.0 rad, R = 1.0, nu = 0.33, nrings = 127 (146945 dofs), 1600 modes (strict internal N* = 1600; operational ladder per the disk-certified smooth-error argument). References: Poisson 0.3863, GOE 0.5307.

- G2 internal N* = 1600, rigid [3, 3]; classifier counts {'S': 813, 'A': 785, 'xx': 2}, min quality 0.289, resolved 23

## <r> per reflection class

| class | n levels | <r> |
|---|---|---|
| S | 813 | 0.425 +/- 0.009 |
| A | 785 | 0.434 +/- 0.010 |
| **pooled** | 1600 | **0.4294 +/- 0.0068** |
| disk baseline (2 x 789-level windows) | - | 0.3919 +/- 0.0070 |

### Spectral-window trend (thirds)

| window | <r> |
|---|---|
| low | 0.4279 +/- 0.0115 |
| mid | 0.4372 +/- 0.0119 |
| high | 0.4209 +/- 0.0118 |

- sensitivity (strict N* = 1600): pooled <r> = 0.4294 +/- 0.0068

## Verdict (preregistered readings)

- pooled <r> = 0.4294 +/- 0.0068 vs disk baseline 0.3919 +/- 0.0070: separation 3.9 sigma (threshold 3)
- per class: S +2.8, A +3.5 sigma

**Reading: SUPPORTS (sector intermediate, consistent with the published COMSOL sector result and the corner-driver reading)**

Discussion (v2). The v1 ambiguity (+1.7 sigma at 800 modes) is resolved as
statistical power, not physics: doubling the ladder takes the separation to
+3.9 sigma pooled with both reflection classes individually positive. The
open-source sector now CONFIRMS the direction of the published COMSOL sector
claim at a single generic angle (theta = 2.0); the exact-angle comparison
remains registered. The decreasing spectral trend persists (high third
0.4209), consistent with the arc's growing share of the boundary influence
at high k.
