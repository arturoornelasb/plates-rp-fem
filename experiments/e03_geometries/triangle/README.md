# E3a -- free equilateral triangle (preregistered before execution)

Paper prediction ("geometry"): geometries integrable for the second-order
Helmholtz billiard but NOT separable for the fourth-order biharmonic plate
should show intermediate statistics at free edges; the paper names the
equilateral triangle and flags it as a "stronger and more speculative test,
since the integrable base is not assured there".

## Design

- Equilateral triangle, side L = 1, D = rho*h = 1, nu = 0.33, Argyris on the
  C3v-symmetric mesh (uniform refinement of the single macro triangle).
- Symmetry: C3v -> irreps A1, A2 (singlets) and E (exact doublets). Classifier
  uses the rotation character (<w,Rw> = +1 for A-type, -1/2 for any vector of
  an E doublet) then the reflection character for A1/A2; E doublets are
  deduplicated to one level per pair before any spacing statistics.
- Exact anchor: on a straight-edged polygon the simply supported plate
  factorizes through the Dirichlet Helmholtz problem; for the equilateral
  triangle that is the Lame lattice, so Lambda_SS = lam_H^2 EXACTLY. The SS
  limit is reached with the validated Winkler spring at kappa = 1e10.

## Gates

- G1: kappa = 1e10 spectrum vs exact Lame^2 -> N*_SS (primary accuracy gate).
- G2: kappa = 0 on refine-7 vs refine-6 -> internal N*.
- G3: classifier sanity at both kappa: character clustering at {+1, -1/2},
  E-doublet splittings at round-off, rigid modes = piston (A1) + tilt pair (E),
  sector fractions ~ 1:1:2 (A1:A2:E individuals).

## Preregistered readings

- SUPPORTS: per-sector <r> of the free triangle (A1, A2, deduped E) sits
  significantly (>= 3 sigma pooled) above the same-protocol SS-triangle
  baseline measured at kappa = 1e10 with identical windows.
- CHALLENGES: free-triangle sectors Poisson-consistent (no separation from
  the SS baseline).
- Per the paper's own caveat, an intermediate-but-weak or sector-inconsistent
  outcome is reported as ambiguous rather than forced into either reading.

## Protocol

n_modes = 1000 certified eigenvalues per kappa in {0, 1e10}; labels from the
banded-deflated eigenvector solve; lowest 10 levels per sector dropped after
E-dedup; <r> per sector and pooled, low/high spectral halves reported.
N_use = min(n_modes, N*_SS, N*_int).
