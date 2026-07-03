# E3b -- free disk (Poisson control) vs free ellipse (preregistered)

Paper prediction ("geometry"): the full free disk is the program's
ZERO-COUPLING CONTROL -- free edges without broken separability must give
Poisson (each fixed-m radial sequence a Berry-Tabor picket fence; pooling the
independent m-sectors within one reflection class gives Poisson). The ellipse
breaks biharmonic separability at free edges -> intermediate statistics
expected. This is the falsifiable dichotomy: free edges ALONE must not
produce intermediate statistics.

## Design

- Disk control: semi-analytic (per-m J/I free-edge determinant, nu = 0.33,
  fundamental validated at x^2 = 5.262 vs literature ~5.25) -- the same
  construction the paper's executed control used. Gives an exact-grade
  reference spectrum AND the Poisson baseline sequences.
- FEM on curved boundaries: inscribed-polygon meshes cap eigenvalue accuracy
  at O(h^2) geometry error (measured: 4.0x reduction per refinement). The
  disk is therefore ALSO the geometry-error gate: (i) strict N* of FEM vs the
  semi-analytic spectrum; (ii) operationally, the FEM-disk class-sequence
  <r> must reproduce the semi-analytic Poisson value -- geometry error is
  smooth, so spacing statistics can remain valid beyond the strict N*; this
  is measured on the control, then transferred to the ellipse (same mesher,
  same refinement).
- Ellipse: semi-axes a/b = 1.6189043236 (campaign aspect), a*b = 1
  (area = pi, disk-matched). Z2xZ2 parity sectors via the validated
  cluster-resolved classifier (centered mirrors). FEM doublet handling on the
  disk: reference-multiplicity-guided dedup.

## Gates

- G1: FEM free disk vs semi-analytic at two refinements -> strict N* and the
  O(h^2) trend.
- G1s (operational): FEM-disk pooled class <r> vs semi-analytic same-window
  value -- must agree within 2 sigma for the production refinement.
- G2: ellipse internal two-refinement N*.
- G3: classifier diagnostics on the ellipse (ambiguous count, quality).

## Preregistered readings

- SUPPORTS: disk control Poisson-consistent (both semi-analytic and FEM) AND
  ellipse per-sector <r> >= 3 sigma (pooled) above the disk baseline measured
  with matched window lengths.
- CHALLENGES: ellipse Poisson-consistent (separability broken but no
  intermediate statistics), OR disk control NOT Poisson (free edges alone
  suffice -- would refute the separability side of the mechanism).
- The FEM-disk statistics gate failing (G1s) voids the ellipse reading at
  this refinement (geometry error corrupts statistics) rather than counting
  for or against the hypothesis.

## Protocol

n_modes = 800 certified eigenvalues; ellipse labels from the banded-deflated
eigenvector solve; lowest 10 per sector dropped; disk baseline = pooled <r>
over disjoint windows of the semi-analytic class sequence with lengths
matched to the ellipse sector ladders.
