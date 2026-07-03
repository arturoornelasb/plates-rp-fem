# Research plan and experiment registry

Goal: execute, with an open-source FEM stack, the numerical tests that
"Why Rosenzweig-Porter?" (the v4 plates paper, Zenodo) assigns to COMSOL,
so every result is independently reproducible. The paper itself sanctions
finite elements for its decisive test ("COMSOL or finite elements", Gap A,
observable 1).

Each experiment lives in `experiments/eNN_*/` with its own code, RESULTS.md,
results.json and figures. Executed experiments are FROZEN (never edited after
their run); shared machinery is factored into `src/platefem` for later
experiments. Accuracy policy and conventions: see CONVENTIONS.md.

## E1 -- Stack validation [DONE 2026-07-03 -- PASS]

Morley (machinery check) and Argyris (production element) against two
references: exact SSSS spectrum and the paper's own Legendre-Ritz FFFF
reference (T1 protocol). Statistics-grade criterion: eigenvalue error
< 0.1 x local mean level spacing.

Result: Argyris FFFF passes on all meshes for 200 modes (median parity down
to ~1e-7); the Winkler edge-spring term reproduces exact SSSS to 3e-6 at
kappa = 1e10 (clean O(1/kappa) over four decades). Three documented pitfalls:
skfem ElementGlobal instance caching; Argyris rigid-cluster pollution at fine
mesh (gap-based detection required, uniform-mesh ceiling ~128x80); the
Legendre-Ritz reference has a ~1e-4 round-off floor (non-monotone in NLEG
beyond ~96 -- diag_ritz_nleg.py), so NLEG 64/72 is its sweet spot.

## E2 -- Protocol B: boundary-controlled transition [RUNNING]

Paper prediction ("Boundary-controlled transition"): per-sector <r> should
interpolate continuously from Poisson (0.3863) at simply supported to the
intermediate/RP value at free, as the symmetric Winkler restraint
V_n + kappa w = 0 (with M_n = 0 throughout) is tuned from kappa = inf to 0.

Method: full rectangle (E1 geometry), Argyris 96x60, ~400 elastic modes with
eigenvectors, kappa grid {0, 1e2 ... 1e8 half-decade, 1e10}; modes classified
into ee/eo/oe/oo by mirror-correlation of probed mode shapes; <r> per sector
(drop lowest 10). Accuracy gates before the sweep: kappa = 1e10 vs exact SSSS
(N*), kappa = 0 vs a finer 128x80 mesh (internal N*), and kappa = 0
sector-resolved check against the Ritz per-sector spectra (classifier
validation).

Preregistered readings (from the paper): SUPPORTS = monotonic transition
0.3863 -> intermediate. CHALLENGES = abrupt transition, non-monotonic
behavior, or intermediate statistics already at nearly-simply-supported
(kappa >= 1e8).

## E3 -- Geometry tests [PLANNED]

Free-edge ellipse and triangle (biharmonic separability broken -> expect
intermediate), full disk and annulus (separable -> Poisson controls), disk
sector (paper's RP case). Curved boundaries need mesh generation beyond
init_tensor and a convergence study on curved-element geometry error.

## E4 -- Gap A on the true operator [PLANNED]

Eigenvector statistics (IPR / D2 vs sector size N, Dq flatness, Sigma^2(L) to
L >= 20) of FEM eigenmodes per symmetry sector. Two known constraints:
(i) T1 basis-dependence -- FEM modes must be projected onto the registered
reference basis (SS sines / FF beam products) before any IPR is computed,
never the raw nodal basis; (ii) the Argyris uniform-mesh ceiling caps full-
plate ladders, so large-N ladders need quarter-plate symmetry reduction
(sector BCs on the midlines) and possibly a C0 interior-penalty formulation.
