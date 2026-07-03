# Conventions and accuracy policy

## Physical conventions

- Kirchhoff-Love flexural eigenproblem: `D lap^2 w = rho h omega^2 w` with
  `D = rho h = 1`. Physical eigenvalue `Lambda = omega^2`.
- Canonical rectangle (from the paper's pretests): `a = 1`,
  `b = 1 / 1.6189043236` (generic irrational aspect ratio ~1.618904, avoids
  accidental degeneracies), Poisson ratio `nu = 0.33`.
- Boundary conditions. Free: `M_n = 0, V_n = 0` -- the natural BCs of the
  bending energy, no constraints imposed. Simply supported: `w = 0` essential,
  `M_n = 0` natural. Winkler edge spring (Protocol B): add `kappa * u * v` on
  the boundary -> `V_n + kappa w = 0` with `M_n = 0` retained; `kappa = 0` is
  free, `kappa -> inf` is simply supported (error ~ 1/kappa; 3e-6 at 1e10,
  measured in E1).
- Symmetry sectors: parity about the midlines `x = a/2`, `y = b/2`; labels
  ee/eo/oe/oo (even/odd in x, then y), matching the Legendre-Ritz sectors.
  Rigid modes: piston (ee), y-tilt (eo), x-tilt (oe) -> 3 for the free plate.

## Spectral statistics

- Spacing ratios `r_n = min(s_n, s_n+1) / max(s_n, s_n+1)`, no unfolding.
  References (Atas et al. 2013): Poisson 0.3863, GOE 0.5307.
- The lowest 10 levels per sector are dropped (non-universal).

## Accuracy policy (statistics-grade)

- Criterion: eigenvalue error < 0.1 x local mean level spacing for every mode
  used. Relative spacing ~ 2/n at mode n, so this tightens with the ladder.
- Production element: Argyris (quintic C1). Morley is O(h^2) -- cross-checks
  only. Fresh element instance per mesh/FacetBasis (skfem ElementGlobal caches
  mesh arrays on the instance).
- Rigid modes are split by the multiplicative gap in the low spectrum
  (> 3 decades), never by fixed tolerance or index: the Argyris rigid cluster
  is polluted at fine mesh (~1e-2 absolute at 128x80, h^-5 Vandermonde
  round-off) while elastic modes stay clean. Uniform meshes beyond ~128x80
  are UNUSABLE with this element.
- Reference floors. Exact SSSS: machine precision. Legendre-Ritz FFFF
  (src/platefem/ritz.py, T1 protocol NLEG 64 vs 72, tol 1e-4): floor ~1e-4;
  its eigenvalues go non-monotone in NLEG beyond ~96 (round-off; evidenced in
  experiments/e01_validation/diag_ritz_nleg.py), so do not tighten the filter.

## Environment

conda env `plates-fem` (Python 3.12, scikit-fem 12.0.2, scipy 1.18, numpy 2.5)
-- see environment.yml. On Windows run scripts with the env's python.exe
directly; `conda run` swallows stdout. `pip install -e .` exposes `platefem`.
