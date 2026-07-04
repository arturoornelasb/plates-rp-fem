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

## E2 -- Protocol B: boundary-controlled transition [DONE 2026-07-03 -- SUPPORTS]

Executed twice: v1 (400 modes, full 14-point kappa map) and v2 (800 modes,
7-point grid). All accuracy gates green in both. v1 established the map and
the finite-size methodology (endpoints match their baselines exactly; exact
SSSS gives <r> = 0.417 at ~100 levels/sector, NOT the asymptotic 0.3863) but
was underpowered (2.1 sigma). v2 (~200 levels/sector): SS end 0.3908 = exact-
SSSS baseline 0.3908; free end 0.4422 = Ritz baseline 0.4434; transition
monotone within errors; amplitude 3.6 sigma > preregistered 3 sigma; no
intermediate statistics at nearly-simply-supported. READING: SUPPORTS the
boundary-controlled-transition prediction.

Paper prediction ("Boundary-controlled transition"): per-sector <r> should
interpolate continuously from Poisson at simply supported to the
intermediate/RP value at free, as the symmetric Winkler restraint
V_n + kappa w = 0 (with M_n = 0 throughout) is tuned from kappa = inf to 0.
Method and preregistered readings: see the experiment folder. Endpoints are
judged against finite-size baselines measured with the identical protocol
(exact SSSS spectrum; independent Legendre-Ritz), never against asymptotic
constants.

## E3 -- Geometry tests [E3a triangle DONE 2026-07-03 -- SUPPORTS; rest planned]

E3a result (all gates green: N* = 1000/1000 against the EXACT Lame^2
reference at 3.0e-5; zero ambiguous labels; 333 exact E-doublets, max
splitting 7.5e-6): free triangle pooled <r> = 0.4889 +/- 0.0099
(intermediate) vs same-protocol SS baseline 0.3753 +/- 0.0134
(Poisson-consistent) -- separation 6.8 sigma (preregistered threshold 3).
Per sector: A1 +3.8 sigma, E +5.6 sigma, A2 +1.9 sigma (all positive; A2 is
the shortest sequence). READING: SUPPORTS the geometry prediction --
breaking biharmonic separability alone (Helmholtz-integrable domain, free
edges) produces intermediate statistics, now in a symmetry group (C3v)
disjoint from the rectangle's Z2xZ2.

Free-edge ellipse and triangle (biharmonic separability broken -> expect
intermediate), full disk and annulus (separable -> Poisson controls), disk
sector (paper's RP case).

E3a (equilateral triangle, preregistered in its README before execution):
straight edges keep the Argyris machinery exact; C3v symmetry needs the new
irrep classifier (rotation character +1 for A-type, exactly -1/2 inside an E
doublet; E doublets deduplicated before statistics; C3v-symmetric mesh by
uniform refinement of the macro triangle). Exact anchor: SS polygon plate =
(Dirichlet Helmholtz)^2 = Lame lattice squared, reached through the validated
kappa = 1e10 Winkler limit -- smoke-validated to 1.7e-6 including doublets.

E3b (disk control + ellipse; v1 + v2 decider executed 2026-07-03): the
curved-boundary strategy is measurement-based -- the semi-analytic free disk
(per-m J/I determinant, fundamental 5.262 vs lit. ~5.25) gates the FEM both
strictly (flat 2e-4 = smooth O(h^2) geometry error through mode 1600) and
operationally (FEM-disk class <r> = semi-analytic value to 0.0 sigma on both
halves of the extended ladder; non-smooth error component max-median =
1.4e-5, ~0.3% of a sector spacing -- cannot wash out repulsion).

DISK CONTROL: Poisson CONFIRMED, 0.3905 +/- 0.0073 (paper: 0.386 +/- 0.007).
ELLIPSE (v2 decider, 1600 modes): pooled <r> = 0.3732 +/- 0.0071,
Poisson-consistent; the v1 upward trend (+2.3 sigma at 800 modes) did NOT
continue (thirds 0.361/0.389/0.370, +0.5 sigma). READING per the
preregistered v2 criteria: CHALLENGES the naive geometry reading ("any
broken separability -> visibly intermediate"). Split E3 verdict (triangle
6.8 sigma intermediate, ellipse Poisson-like) points to boundary CORNERS as
the strong evanescent-coupling drivers; within the RP framework a dense but
tiny V (small lambda) is indistinguishable from Poisson at accessible N.
Registered follow-up: Gap A eigenvector diagnostics on the ellipse to
separate 'no coupling' from 'tiny dense coupling'.

E3c (disk sector, theta = 2.0, executed 2026-07-03): all gates green
(internal N* = 800/800, zero ambiguous S/A labels). Pooled <r> = 0.4213 +/-
0.0095 vs disk baseline 0.3972 +/- 0.0100: +1.7 sigma, both classes positive
-- AMBIGUOUS per preregistration (weak intermediate). Fits the emerging
corner-sharpness ordering (triangle 60deg: 0.489 > rectangle 90deg:
0.44-0.46 > sector 90-115deg + arc: 0.421 > smooth ellipse: 0.373 ~ disk
control: 0.391). Not a contradiction of the published COMSOL sector study
(10^4 modes, angle sweep) and not yet a confirmation; a decisive comparison
needs longer ladders or their exact angles.

## E4 -- Gap A on the true operator [DONE 2026-07-03 -- LOCALIZED-LEANING at N <= 324]

First execution of the paper's decisive test on TRUE-operator modes (1200
certified eigenpairs, vector residuals <= 8e-6, all gates 1200/1200, zero
ambiguous labels). Result: IPR is N-INDEPENDENT in both registered bases
(ladder D2 ~ 0.02--0.13 vs GOE ~1); in the near-lossless beam basis
(captured p = 0.999--1.000) each free-plate mode is ~1--2 beam products
(IPR ~ 0.7). Per the preregistered dichotomy this is the localized/banded
side at accessible N: sparse strong pair-hybridizations (microscopically,
the avoided crossings of the original observations) producing intermediate
<r> WITHOUT RP-fractal eigenvectors. Sharpenings: the true-operator D2 is
LOWER than the T1 Ritz-model values (the model route overestimates coupling
density); and since effective coupling grows with mode index, the verdict
could flip at the paper's registered N ~ 2048 rung -- quarter-plate
reduction is the registered instrument for that follow-up. Sigma^2(L <= 5)
intermediate, mildly below Poisson; larger L unfolding-limited.

## E5 -- Superellipse sweep [DONE 2026-07-03 -- SEPARABLE-BOUNDARY STEP]

Novel discriminator (|x/a|^p + |y/b|^p = 1, p in {2,3,4,6,10}, 1600 modes
each, identical protocol; p = 2 = the E3b ellipse). Result: <r>(p) = 0.373,
0.494, 0.489, 0.494, 0.516 -- a +12.2 sigma STEP between p = 2 and p = 3
(soft smooth corners!), a flat plateau p = 3..6 (0.7 sigma spread), and a
mild +2.3 sigma rise at p = 10. Neither preregistered shape fits; the
operative variable is the boundary's compatibility with separable
coordinates (the ellipse is the family's unique quasi-separable member),
corners are secondary. Statistics refinement-stable at 0.1 sigma for p = 3
and p = 10 (diag_stability.py).

## Campaign synthesis (2026-07-03)

Seven geometries, one hierarchy: (i) exact reduction -> Poisson (disk
0.391); (ii) coordinate-adapted / quasi-separable boundaries -> weak sparse
coupling (ellipse 0.373, sector 0.421, rectangle 0.442--0.46, whose true
modes are ~single beam products per E4); (iii) unadapted boundaries -> full
coupling (triangle 0.489, superellipse p != 2: 0.49--0.52). The
boundary-controlled transition (E2, 3.6 sigma) confirms the coupling lives
in the free-edge condition. Gap A on the true operator (E4) shows the
accessible-N mixing is SPARSE (avoided-crossing hybridization), not
RP-fractal -- the model-based (Ritz) route overestimated coupling density.
Open registered follow-ups: (a) the N ~ 2048 Gap A rung via quarter-plate
reduction (where the RP phase could still appear); (b) longer sector
ladders / the published study's exact angles; (c) Dq flatness and form-
factor rival exclusion on the strongest-coupling geometries (triangle,
superellipse p = 10).
