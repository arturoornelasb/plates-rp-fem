# plates_fem — open-source FEM campaign for the v4 plates paper

Goal: execute the COMSOL-targeted numerical tests of
`zenodo_v4/paper/why_rosenzweig_porter_plates.tex` on this machine with an
open-source stack, so the results are reproducible and independent of COMSOL.
The paper itself sanctions this ("COMSOL **or finite elements**", Gap A,
observable 1).

Stack: **scikit-fem (Morley element) + scipy shift-invert Lanczos**, conda env
`plates-fem` (Python 3.12, scikit-fem 12.0.2). Rationale: the Kirchhoff operator
is fourth order (C1 problem); Morley is the standard nonconforming element for
it; free edges (M_n = 0, V_n = 0) are the *natural* BCs of the bending energy, so
the FFFF plate needs no constraints at all; simply supported is `w = 0` at
boundary deflection dofs; Protocol B's Winkler spring is one boundary term
`kappa * w * v ds`. FEniCSx was considered and set aside: its native Windows
build lacks PETSc/SLEPc, and biharmonic there needs an interior-penalty DG
formulation whose free-edge terms are delicate.

## Plan (agreed 2026-07-03)

1. **step1_validation/** — validate the stack: SSSS rectangle vs exact analytic
   spectrum; FFFF rectangle vs the repo's spectral Rayleigh–Ritz reference
   (T1 conventions: a/b = 1.6189043236, nu = 0.33, D = rho*h = 1).
2. **Protocol B** — Winkler-spring kappa sweep on all four edges
   (paper, Prediction "boundary-controlled transition"): per-sector <r> from
   Poisson (~0.386) toward the RP value as kappa: inf -> 0.
3. **Geometry tests** — ellipse, triangle (broken biharmonic separability ->
   expect intermediate), annulus and full disk (Poisson controls), disk sector.
4. **Gap A on the true operator** — eigenvector statistics (IPR/D2, Dq flatness,
   Sigma^2(L)) of FEM eigenmodes per symmetry sector. NOTE the T1 audit result:
   D2 depends strongly on the representation basis; FEM nodal coefficients are
   yet another basis, so modes must be projected onto the registered reference
   basis (SS sines / FF beams) before any IPR is computed.

Caveats tracked from step 1 onward:

- Statistics-grade accuracy criterion: eigenvalue error < 0.1 x local mean level
  spacing (relative spacing ~ 2/n at mode n), not a blanket relative tolerance.
- Morley converges at O(h^2): fine as a machinery cross-check, too slow for
  high mode ladders. Production element is Argyris (quintic C1).
- Argyris in scikit-fem is an ElementGlobal (basis built numerically from local
  Vandermonde systems mixing value/derivative dofs, condition ~ h^-5). Empirical
  pollution of the near-zero rigid cluster: Rayleigh quotient of interpolated
  w=1 is 2.4e-6 at mesh 64x40 but 1.0e-2 at 128x80 (elastic modes still clean;
  Jacobi equilibration does NOT help -- it is assembly error, not solver error).
  Consequences: (i) rigid modes are split off by the multiplicative GAP in the
  low spectrum (>4 decades), never by fixed tolerance or index; (ii) uniform
  meshes beyond ~128x80 are unusable with this element -- the multi-thousand-mode
  Gap A ladders will need quarter-plate symmetry reduction and/or a C0
  interior-penalty formulation. Ceiling documented in step1 RESULTS.md.
