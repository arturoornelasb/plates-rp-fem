# E2 v2 -- long ladder

v1 (parent folder) established the full kappa map with all accuracy gates
green, endpoints matching their finite-size baselines exactly (SS end = exact
SSSS to 4 decimals; free end = independent Ritz), but a transition amplitude
of only 2.1 sigma at ~100 levels/sector -- underpowered against the
preregistered 3 sigma.

v2 doubles the ladder (800 modes, ~200 levels/sector) on the finer production
mesh (Argyris 128x80) with a reduced kappa grid {0, 1e4, 3e5, 1e6, 3e6, 1e8,
1e10}: endpoints plus the transition core. Expected: the separable baseline
falls toward Poisson with ladder length while sem halves, projecting the
amplitude to ~4-5 sigma if the free-plate value holds.

Same scripts as v1 (frozen there), config block adjusted; analyze.py computes
the same finite-size baselines at the new ladder length.
