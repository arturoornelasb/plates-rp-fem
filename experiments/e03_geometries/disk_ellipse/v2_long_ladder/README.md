# E3b v2 -- long ladder (registered decider)

v1 (parent folder) confirmed the disk Poisson control and found the ellipse
AMBIGUOUS: sub-Poisson at low modes (0.354) with a significant upward
spectral trend (+2.3 sigma; high third at 0.4128) -- the paper's
preregistered weak-coupling caveat ("a single weak-coupling run must not be
read as refutation"). Physically plausible: the ellipse's smooth boundary
(no corners) couples more weakly at low k than the triangle's 60-degree
corners.

v2 doubles the ladder to 1600 modes at the same refine-6 mesh. Validity
rests on the operational disk gate, now applied ALSO to the upper half of
the extended ladder (discretization error grows with mode index and is not
smooth, so it must be shown harmless exactly where the verdict lives).

Preregistered readings for v2:
- SUPPORTS: upper-half (or high-third) ellipse <r> >= 3 sigma above the
  matched disk baseline, continuing the v1 trend.
- CHALLENGES: trend absent in the extended ladder (high windows fall back
  to the separable baseline).
- VOID: the upper-half disk statistics gate fails at refine 6 (then E3b
  needs refine 7 or curved elements before any reading).
