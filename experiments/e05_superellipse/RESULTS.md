# E5 -- superellipse corner sweep (RESULTS)

|x/a|^p + |y/b|^p = 1, a/b = 1.618904, nu = 0.33, refine-6, 1600 modes per p (p = 2 reused from E3b v2, identical protocol). Disk baseline 0.3915 +/- 0.0070. References: Poisson 0.3863, GOE 0.5307.

- G2 (p = 10.0, worst corners): internal N* = 61

## <r>(p)

| p | pooled <r> | low third | mid | high third | vs baseline |
|---|---|---|---|---|---|
| 2 | **0.3732 +/- 0.0071** | 0.3613 | 0.3886 | 0.3695 | -1.8 sigma |
| 3 | **0.4944 +/- 0.0069** | 0.5117 | 0.4946 | 0.4794 | +10.4 sigma |
| 4 | **0.4893 +/- 0.0069** | 0.5012 | 0.5069 | 0.4591 | +9.9 sigma |
| 6 | **0.4942 +/- 0.0067** | 0.4971 | 0.4941 | 0.4929 | +10.5 sigma |
| 10 | **0.5159 +/- 0.0066** | 0.5396 | 0.5117 | 0.5004 | +12.9 sigma |

## Verdict (preregistered readings)

- <r>(p) sequence: p=2: 0.3732, p=3: 0.4944, p=4: 0.4893, p=6: 0.4942, p=10: 0.5159
- p = 10 vs disk baseline: +12.9 sigma
- monotone non-decreasing within noise: True
- shape: jump p=2 -> 3 = +12.2 sigma; plateau spread p=3..6 = 0.7 sigma-units; p=6 -> 10 rise = +2.3 sigma

**Reading: SEPARABLE-BOUNDARY STEP: the coupling switches on fully as soon as the boundary ceases to be a coordinate line of a separable system (already at p = 3, with soft smooth corners), then plateaus with only a mild further rise toward sharper corners. NEITHER preregistered shape fits: not curvature-graded, not a true-corner threshold. The operative variable is the boundary's compatibility with separable coordinates -- corners are secondary.**

Validity: statistics are refinement-stable at 0.1 sigma for both the headline point (p = 3: refine-5 0.4957 vs refine-6 0.4944) and the flagged worst-corner point (p = 10: 0.5166 vs 0.5159) -- diag_stability.py; the strict per-eigenvalue internal N* fails at p = 10 (61/1600) because coarse-mesh corner-region geometry error is large, but the spacing statistics do not inherit it (same smooth-error mechanism measured on the disk in E3b).

Discussion. This result REINTERPRETS the E3b ellipse verdict and supersedes the corner-sharpness reading of E3c: the ellipse is not 'a smooth boundary with weak corners' but the unique quasi-separable member of this family (confocal boundary = elliptic-coordinate line, giving the two biharmonic parts a near-diagonal angular structure). The campaign hierarchy consistent with all seven geometries: (i) exact reduction (disk) -> Poisson; (ii) coordinate-adapted / quasi-separable boundaries (ellipse 0.373, sector 0.421, rectangle 0.44 -- whose true modes E4 showed are ~single beam products) -> weak sparse coupling; (iii) unadapted boundaries (triangle 0.489, superellipse p != 2 at 0.49--0.52) -> full coupling. This SHARPENS the paper's separability mechanism: what matters is adaptedness of the boundary to ANY separable structure for the fourth-order operator, and the free-edge coupling magnitude is set by how badly that fails -- not by corners per se.
