# E5 -- superellipse corner sweep (preregistered; novel experiment)

Motivated by the E3 split verdict (free triangle strongly intermediate at
6.8 sigma; free ellipse Poisson-consistent at 1600 modes): what actually
drives the evanescent coupling at free edges? The family
|x/a|^p + |y/b|^p = 1 separates the candidate mechanisms cleanly, because
biharmonic separability is equally broken for EVERY p in this family
(including p = 2), while the boundary's corner-region curvature grows
continuously with p at fixed C-infinity smoothness.

## Design

p in {2, 3, 4, 6, 10}; p = 2 is the E3b v2 ellipse (data reused, identical
protocol); p in {3, 4, 6, 10} run here. Same campaign parameters throughout:
a/b = 1.6189043236, a*b = 1, nu = 0.33, Argyris refine-6 (radially mapped
init_circle; boundary vertices exactly on the superellipse), 1600 certified
modes, Z2xZ2 parity sectors via the validated cluster-resolved classifier,
lowest 10 per sector dropped, matched-window disk baseline (E3b machinery).
Gates: mesh quality per p; internal two-refinement N* at p = 10 (worst
corner-region elements -- if the sharpest case passes, milder ones do);
the smooth-geometry-error statistics argument is inherited from the E3b
operational gate at the same refinement and h.

## Preregistered readings

- CURVATURE-GRADED COUPLING: <r>(p) rises monotonically with p, with p = 10
  significantly (>= 3 sigma) above the disk baseline. Refines the mechanism:
  coupling scales with concentrated boundary curvature.
- TRUE-CORNER THRESHOLD: <r>(p) stays Poisson-consistent through p = 10.
  Refines the mechanism differently: only genuine corners (curvature
  singularities) generate effective evanescent coupling at these mode
  numbers.
- Either reading SHARPENS the paper's hypothesis (both identify the
  microscopic driver of V); an erratic non-monotone pattern instead flags
  artifacts and voids the sweep.
