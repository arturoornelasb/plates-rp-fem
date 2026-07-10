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

## E6 -- Designed contact coupling + control law [DONE 2026-07-04 -- PARTIAL/IN-DOMAIN]

Two identical free plates coupled by post quartets (paper's contact test,
built directly from certified mode shapes; exact antisym-block reduction).
Eigenvector dichotomy CRISP (dense IPR 0.033 extended vs banded 0.18-0.19
localized; central rank-1 control); spacing dichotomy directionally right at
2.3 sigma (rank-deficiency capped; v4 with ~200 quartets registered). P10
CONTROL LAW validated within its physical domain: lambda =
Delta/(2 q phi_i(x*)^2) + closed-form 2nd-order correction gives 3.6% at
Delta = 0.25 d_loc; the domain boundary (Delta ~ 0.5 d_loc) is predicted by
pair-content collapse.

## E7 -- Completeness [DONE 2026-07-04]

(a) Protocol-A mixed config (free-x / SS-y) = the Levy-separable plate ->
lands ON the Poisson baseline (0.3973 +/- 0.0103): the adaptedness hierarchy
required this; one free edge pair does not suffice. (b) nu sweep: FEM
reproduces the T1 Ritz bonus value-by-value (~0.003 agreement); the
effective coupling is MATERIAL-TUNABLE (+0.029 over nu = 0.25..0.36,
~3 sigma) -- Prediction-5 recast supported over strict material
independence.

## E8 -- Unified Gap A across geometries + Dq [DONE 2026-07-04]

Same-domain SS-eigenbasis coefficients (exact M-inner products). The
eigenvector hierarchy MIRRORS the E5 spacing hierarchy: adapted boundaries
have flat-IPR sparse vectors (rectangle 0.18, cross-validating E4;
ellipse 0.40 -- the registered decider lands on TINY SPARSE coupling, ~2.5
SS modes per free mode, N-independent); unadapted boundaries have genuinely
SCALING vectors (triangle D2 = 0.27/0.36; superellipse p=10 D2 =
0.22-0.50) -- the RP-candidate regime sits exactly where the spacing
statistics are strongest. RP-vs-PBRM adjudication split at these N
(superellipse Dq-flat/RP-leaning; triangle multifractal-leaning); decision
registered for the N ~ 2048 rung (C0-IP build).

## E10/E11 [running 2026-07-04]: annulus control; Mindlin thickness sweep.

## E9 -- Gap A at the registered N = 2048 rung [DONE 2026-07-04 -- RP PHASE DEVELOPS]

T1 truncated-operator ladder with exact integer-arithmetic assembly
(platefem/ritz_exact.py; quadrature floor removed, eigh applied only to
physically-normed rung submatrices). IPR falls cleanly over three octaves
(N = 256..2048) in both sectors and both registered bases: ladder D2 =
0.42-0.50 (SE ~ 0.02-0.03) with FLAT Dq (D_1.5 - D_4 = 0.069-0.100 -- the
registered PBRM rival band 0.20-0.28 excluded). The non-ergodic fractal
(RP) phase develops at the registered ladder; the P12 preregistered D2
value (0.76 +/- 0.15) is NOT confirmed (measured 0.46). Protocol
reconciliation with the true-operator windows (flat at N <= 324, E4/E8)
is the registered follow-up (C0-IP discretization; skfem Lagrange lacks
second derivatives -> custom reference-hessian assembly needed).

## E10 -- Annulus control [DONE 2026-07-04 -- SUPPORTS]

Second exactly-reduced Poisson control confirmed: semi-analytic 4x4
J/Y/I/K determinant gives pooled class <r> = 0.3734 +/- 0.0076 (1333
ratios), per-m picket fences; FEM ring mesh cross-validates at N* =
200/200, max relerr 4.3e-4.

## E11 -- Mindlin thick plates [DONE 2026-07-04 -- SUPPORTS, +5.4 sigma]

P2/P2 SRI Mindlin (validated: rigid trio, plane-probes, thin-limit
convergence): pooled <r> rises monotonically 0.4540 -> 0.5568 over
t = 0.02 -> 0.15, continuous with Kirchhoff at the thin end and reaching
GOE by t ~ 0.1 -- the paper's thick-plate prediction (shear channels
increase net coupling toward GOE) confirmed.

## E14 -- Gap A reconciliation [DONE 2026-07-04 -- PROTOCOL ARTIFACT]

The C0-IP instrument (custom analytic-hessian P4 Lagrange + Brenner-Sung
interior penalty; no ElementGlobal ceiling) resolves the E4/E8-vs-E9
contradiction: true-operator window IPR is FLAT through N = 512/sector
(sine 0.174 sector-mean; beam 0.63-0.72 = modes remain ~1-2 beam products)
where the registered truncated ladder falls to 0.095. E9's RP-like scaling
is a property of the truncation protocol. GAP A, stated cleanly: at
accessible mode numbers the true free-plate operator is in the
sparse/quasi-separable regime; truncated-ladder numerics (incl. the P12
recipe on COMSOL data) would report RP as an artifact. Gates: G2 FFFF vs
certified Argyris N* = 1200/1200 at 1.2e-4; G1 SSSS N* = 1099/2400
(ladders use modes within gate range). Registered next: true-operator
windows at N ~ 1024-2048 (larger C0-IP meshes) to bound where/if genuine
scaling begins.

## E3c v2 / E13 / E6 v4 [DONE 2026-07-04]

Sector at 1600 modes: SUPPORTS +3.9 sigma (S +2.8, A +3.5) -- the v1
ambiguity was statistical power; confirms the published sector direction at
a generic angle. Aspect sweep: chi2 p ~ 0.22 across five ratios -> robust.
E6 v4: dense contact rank saturates intrinsically (~75/95); spacing
dichotomy plateaus at 2.8 sigma; carried result = eigenvector dichotomy +
in-domain control law.

## E12 -- P7 damping [v1 CHALLENGES-as-run; v2 distributed running]

Single-patch damping is only weakly non-proportional (commutator 3e-2,
near-diagonal modal C): resonances shift quasi-independently -> 2D-Poisson
markers. Physical lesson recorded; v2 with 24 distributed quartets (the
dense-contact analog) is the fair realization of the P7 dense case.

## E12b -- AI-dagger vs Ginibre fine call [DONE 2026-07-09 -- UNRESOLVED per frozen gates; mid-flow located]

4 patch realizations x 4 sectors pooled (~2214 bulk ratios at gamma* =
16), baselines sharpened 3x (AI-dagger/GinUE separation 0.046 vs pooled
scale se = 0.015 -- the call is resolvable IF saturated). Verdict:
UNRESOLVED -- the pooled point (0.713, 0.117) sits MID-FLOW on the
Poisson -> AI-dagger crossover, closer to AI-dagger (6.0 se) than to
GinUE (9.0 se) or Poisson (9.2 se), but not saturated. The v3 "within
~3 sigma of AI-dagger" reflected its coarser marker scale. Registered
next (per the frozen reading): drive the crossover to saturation via a
damping-model refinement (denser/stronger non-proportional
dissipation); statistics are no longer the limit.

## E15 -- Coriolis rotor -> GUE [PLANNED; plan preregistered 2026-07-04]

Full execution plan in experiments/e15_coriolis_rotor/PLAN.md. Core design:
in-plane (plane-stress) elastodynamics of a free 2D rotor (the flat
Kirchhoff plate is Coriolis-blind about its normal); Hermitian gyroscopic
pencil with real frequencies; two-axis (Omega, mistuning) sweep on a
C2-symmetric chiral star realizing the CORRECTED mechanism from the P11
deep-research sweep and toy pretest: rotation alone is signature-protected
(R_pi x T) and must stay GOE; breaking the signature (point-mass mistuning
or fully asymmetric geometry) reaches GUE. Controls: protected star, mirror
star, fully asymmetric star, disk. Anchor: classical in-plane free-disk
determinant. Budget: ~1 session (build 2-3 h, runs 2-3 h, analysis 1 h).

## E15 -- Coriolis rotor -> GUE [DONE 2026-07-04/05 -- SUPPORTS AT FULL SCALE]

Executed by a second agent, independently reviewed (spurious-band lam_cap
fix), then run at the registered scale (~1080 spacings/point, 48k dofs,
refinement-stability gate PASSED at 1.4 sigma worst cell). Physics
correction to the plan, verified two ways: the surviving antiunitary is
sigma_v * T (G0 anticommutes with mirrors, commutes with R_pi; measured at
1e-15). Results: chiral rotor GOE (0.5303) -> exactly GUE (0.5992 +/-
0.0071, -0.1 sigma; +6.5 sigma crossover; probe-scale overshoot resolved as
finite-n); mirror rotor pinned to GOE at every speed (GUE excluded at 9-10
sigma); two-axis mistuning theorem: both ingredients necessary, together
sufficient (mistuned cell GUE +0.7 sigma; protected-vs-mistuned +7.8
sigma). First in-silico mechanical GOE -> GUE with its protection theorem
on the same operator. Registered next: prestressed rotor model;
bladed-disk FE re-analysis (experimental variant).

## E15b -- Industrial-speed feasibility [DONE 2026-07-05]

Reduced-model sweep (N = 1200 certified chiral modes, windowed <r>) across
machine speeds and materials. Rotation speed IS the controlling knob (the
crossover tracks Omega_nd in every run), but the burst limit makes the
achievable dimensionless speed a radius- and RPM-free MATERIAL invariant,
Omega_nd,max ~ v_max / c0 = sqrt(sigma (1 - nu^2) / E): steel 0.067,
Al 7075 0.080, CFRP 0.122 -- all 2-4x below the bulk crossover ~0.293, so a
survivable stiff rotor keeps its low spectrum GOE at any speed. The
observable is frequency-resolved (f*(Omega) falls with mode density: low
spectrum GOE, high spectrum GUE, windowed <r>). Two quantified routes:
(1) deep ultrasonic windows in metals at the rim limit (steel R = 0.5 m:
f* ~ 294 kHz at mode ~23,000, modal-resolution Q ~ 5e4 -- marginal, needs
high-Q RUS in vacuum); (2) SOFT ROTORS (elastomers, Omega_nd,max ~ 1-1.5):
crossover in the first ~100 modes at fully safe speeds. Caveat: at
Omega_nd ~ 1 the centrifugal prestress deferred in E15 is NOT small for
soft materials -- E15c is the registered prerequisite for a quantitative
soft-rotor prediction.

## E15c -- Prestressed rotor: the soft-rotor prediction [machinery ready; RUNNING 2026-07-07]

Adds the Omega^2 physics deferred in E15: exact spin softening
(-Omega^2 M) and the geometric stiffness Kg of the static centrifugal
prestress (unit-Omega^2 modal static solve on the certified modes; Kg
assembled from the resulting stress field). Machinery smoke-validated:
free-particle whirl omega = Omega doubly degenerate; prestress re-stiffens
what softening softens; strain metric eps_max = s1 * Omega_nd^2. Question:
does the GOE -> GUE crossover survive the Omega^2 terms inside the
linear-validity window (centrifugal strain <= ~10-15%)? Preregistered
reading (verdict logic frozen in run_prestress.py before execution):
SURVIVES if the prestressed chiral rotor's top-third-window <r> reaches
the GOE/GUE midpoint (0.565) within strain <= 15% with the mirror rotor's
sigma_v T protection intact; otherwise PARTIAL (onset visible but not
completed inside linear validity) and the hyperelastic/large-deformation
model is registered as the follow-up. Setup: chiral + mirror stars, polar
mesh (24, 72), 1203 certified modes, same protocol as E15 full scale;
bare-Coriolis pencil computed alongside as the contrast.

## E15c -- Prestressed rotor [DONE 2026-07-07 -- PARTIAL; prestress ACCELERATES the crossover]

Executed at the preregistered setup (resid <= 9.2e-4; rigid load frac
~1e-14; Om = 0 prestressed == bare exactly -- built-in null). Measured
strain scales: s1 = 1.719 (chiral) / 0.984 (mirror). Findings: the
prestressed pooled <r> exceeds the bare pencil at EVERY Om > 0 (up to
+0.030) -- the Omega^2 physics (spin softening + centrifugal Kg) ADVANCES
the GOE->GUE crossover rather than washing it out; the sigma_v T
protection survives prestress (mirror 0.43-0.49 at all speeds, GUE
excluded). Frozen-criteria verdict PARTIAL: within strain <= 15% the
top-third window reaches <r> = 0.547 at Om = 0.25 (10.7%); the GOE/GUE
midpoint is crossed at Om = 0.30 (15.5% -- just outside). Soft-rotor
prediction delivered: silicone disk R = 0.1 m, E = 2 MPa -> crossover
onset at ~1,100-1,300 RPM, measured window ~2.4 kHz -- laboratory-trivial
speeds. REGISTERED NEXT (E15d): hyperelastic / finite-deformation rotor
(neo-Hookean prestate) to carry the quantitative prediction through the
15-25% strain regime where completion occurs.

## E15d -- SVK finite-deformation rotor [DONE 2026-07-09 -- COMPLETES]

Instrument (validated by e15d_smoke: tangent consistency 3e-16, exact
spin softening, 0.16% vs linear at small Omega, quadratic Newton to
machine floor): SVK plane-stress prestate by Newton continuation
(translations projected = deformed-CoM pinning; rotation torque-free),
linearized rotating-frame pencil about u0 solved via the validated E15
companion in the tangent eigenbasis (run 1's hand-rolled companion
mistreated the gyroscopic block -> NaN; fixed, documented). RESULT:
**the crossover COMPLETES inside SVK validity** -- chiral top-third
<r> = 0.5671 at Omega_nd = 0.35 (strain 17.6% <= 25%), midpoint 0.565
crossed; mirror protection intact under finite deformation (pooled
<= 0.501 at all speeds to Omega_nd = 0.5); pencil real to ~1e-14
everywhere. Soft-rotor experiment fully predicted: silicone disk
R = 10 cm crosses at ~1,500 RPM. Beyond validity (strain 27-32%) the
window statistic declines -- reported, unclaimed. Neo-Hookean
refinement NOT needed for the verdict.

## E16 -- Sector angle-sweep (protocol-faithful COMSOL replication) [PREREGISTERED 2026-07-07]

The published COMSOL-2024 sector result turns out to be a pooled
ANGLE-SWEEP statistic, not a fixed-angle ladder: ~70 out-of-plane modes
per angle, angle swept (range unpublished; 50 deg shown), S/A classes,
14,000 modes, pooled P(r) envelope ~ RP (source PDF archived in the
dossier repo, 05_avoided_crossings_comsol/). This replaces the earlier
"sector at the published angles" follow-up, which is ill-posed as stated.
E16 replicates the protocol faithfully: theta in [30, 150] deg x 200
values, first 70 elastic modes per angle, S/A by the bisector mirror,
pooled statistics in both ratio conventions. Preregistered readings in
experiments/e16_sector_sweep/README.md: R1 replication vs same-protocol
pooled-Poisson baseline (>= 3 sigma SUPPORTS); R2 protocol scrutiny (the
E14 lesson applied to the source: pooled sweep vs matched-window
fixed-angle references at 50 / 90 / 114.59 deg); R3 hierarchy consistency
(sector = coordinate-adapted tier; compare the E3c ladder 0.4294 and disk
0.391). Run after E15c frees the CPU.

**[DONE 2026-07-07 -- SUPPORTS +9.5 sigma; pooling BENIGN;
hierarchy-consistent.]** All gates green (two-mesh N* = 70/70 at all
three reference angles; zero ambiguous labels over 200 x 70 modes).
R1: pooled <r~> = 0.4266 +/- 0.0027 (9200 ratios; S 0.4336 / A 0.4179)
vs same-protocol pooled-Poisson 0.3865 +/- 0.0032 -> +9.5 sigma: the
published direction REPLICATES under the faithful protocol. R2: pooled
vs matched-window fixed-angle mean +0.1 sigma -> the sweep/pooling
protocol is BENIGN (the E14 worry does NOT apply to the published sector
methodology). R3: 0.4266 vs E3c ladder 0.4294 +/- 0.0068 and disk
0.3905 -- the coordinate-adapted tier is protocol- and window-consistent.
Notable: strong real angle gradient in the low-mode window (50 deg 0.347
Poisson-consistent; 90 deg 0.443; 114.6 deg 0.483) averaged by the pooled
statistic; Lambda- vs frequency-scale conventions agree in practice
(0.4266 vs 0.4261).

## E17 -- Gap A: true-operator windows at N = 1024-2048 [PREREGISTERED 2026-07-07; runner ready]

The registered continuation of E14 over the full ladder range of the
paper's registry (N = 256..2048 per sector). Instrument: the same C0-IP
P4, scaled to mesh (216, 133) ~ 461k dofs, 5600 certified modes
(~1400/sector; the N = 2048 window [819, 1229) needs >= 1229 in the
smallest sector); T1 projection bases at n_funcs_axis = 96 (ratio-form
beam functions, stable), Gauss grid (384, 240), chunked projections.
Gates: G1 SSSS-vs-exact (informational), G2 FFFF vs certified Argyris
first 1200, and the DECISIVE G3 two-mesh eigenvalue gate (216,133) vs
(184,114) over all 5600 -- ladders use only gate-covered modes, and a
coverage shortfall is reported as COVERAGE-LIMITED with the continuation
mesh registered. Preregistered readings (frozen in
experiments/e17_gapA_true2048/README.md + verdict block in run_e17.py):
references E9-truncated sine IPR 0.052 at N = 2048 vs the E4/E14 flat
level ~0.174; I2048 >= 0.14 -> PROTOCOL ARTIFACT THROUGH THE FULL
REGISTERED RANGE; I2048 <= 0.09 -> SCALING ONSET (report fitted D2_true
over N >= 512; RECONCILED-RP if in E9's 0.42-0.50 band); else
INTERMEDIATE. Budget: overnight (~4-10 h solve, ~25 GB peak RAM); run
serially after E15c and E16.

## E17/E17q -- Gap A at N = 1024-2048 [DONE 2026-07-09 -- PROTOCOL ARTIFACT THROUGH THE CERTIFIED RANGE (N <= 1024); indicative 2048 agrees]

Executed on the validated i01 quarter-plate instrument after the serial
run was killed externally post-G2 (its gates preserved: G1 = 1836/5600
informational, G2 = 1200/1200). Four sectors x 1408 elastic modes at
~120k dofs, resid <= 1e-4; decisive per-sector two-mesh gate
((110,68) vs (102,63)): N* = 951/986/983/952 -> ladders certified
through N = 1024 (window [409,614] fully covered; the 2048 window needs
1229). RESULT: the true operator's sine-basis window IPR is FLAT --
0.16-0.19 across N = 256..1024 in all four sectors independently
(fitted D2_true = -0.011); beam basis 0.63-0.72 (each mode ~1-2 beam
products) -- vs the registered truncated ladder's 0.069 at N = 1024
(interpolated) and 0.052 at 2048 (D2 = 0.42-0.50). The INDICATIVE
N = 2048 row (pre-clip windows including modes beyond the gate) sits at
0.1824 -- identical to the certified flat level. Gap A on the rectangle
closes at every registered scale: sparse avoided-crossing hybridization
throughout; truncated-ladder numerics (incl. the P12 recipe) manufacture
the RP phase as a protocol artifact at every rung. CERTIFICATION
HARDENED (2026-07-09, CERT_RESULTS.md): two independent gates agree --
honest 1.196x two-mesh N* = 680-766, SSSS-exact absolute anchor
N* = 742-777 (max relerr ~1e-3 at mode 1408) -- the N = 1024 window is
certified by BOTH gates in all sectors with >= 66 margin; the original
1.08x check (951-986) confirmed correlated-optimistic (methods lesson);
the 2048 window is honestly NOT certifiable at this mesh (needs
1229/sector; finer production ~(140,86)/quarter registered, low
priority given the uniform flatness of the indicative row).
**FINE-MESH RUN EXECUTED (2026-07-10, E17q-FINE, RESULTS_FINE.md):**
production (140,86) ~193k dofs/sector, honest 1.167x check -> N* =
1025/1093/1064/1029 (N = 1024 window re-certified at a SECOND
independent production mesh) and the certified-rung IPRs agree with the
coarse mesh to +/-0.001 -- TOTAL mesh convergence of the physics. The
2048 window (needs 1229) remains formally indicative; the certification
question is DECLARED CLOSED: two production meshes agree to 3-4 decimals
on every certified window and on the (coarse) indicative 2048 row
(0.182 = the flat level), so a third, ~250k-dof/sector mesh would stamp
already-double-converged physics. Threshold note: the certification
floor scales with dofs (7.7e-5 @120k, ~1.4e-4 @193k on rigid-containing
bands) -- resid_sanity 1e-3 for free quarter solves at >150k dofs, with
accuracy arbited by the two-mesh gate (the E18/E15c-grade convention). Notes: oo's "rigid
ratio 14" is a diagnostic-formula artifact (the oo twist mode is
legitimately low; ratio matches (mid/first)^2; n_exp = 0, nothing
dropped). Ops lessons recorded: session-harness reaps long background
tasks -> detached OS-owned launches + stage caching are the pattern.

## i01 -- Quarter-plate C0-IP instrument [VALIDATED 2026-07-08 -- PASS, promoted]

Sector reduction for the rectangle's Gap-A runs: per-edge BCs in
platefem/c0ip.py (free / ss per edge / NEW guided = one-sided
Brenner-Sung Nitsche block enforcing w_n = 0 on symmetry lines). The four
parity sectors become four independent quarter problems (~1/4 dofs,
~1/4 modes each; orthogonalization cost ~ modes^2 drops ~16x;
parallelizable across processes -- an E17-class run drops from ~46 h
serial to ~2-3 h). Gates (experiments/i01_quarter_c0ip): A = merged
quarter sectors reproduce same-mesh full SSSS at the per-sector
statistics-grade criterion through 239/240 (required >= 188, the full
instrument's own exact-accuracy range; residual top-edge mismatch is the
FULL mesh's diagonal asymmetry -- the quarter is exactly symmetric);
B = per-sector FFFF quarter-vs-full 48/48 in all four sectors (relerr
~2e-4), rigid counts exact (ee/eo/oe 1, oo 0). E18 and future rectangle
rungs use this instrument.

## E18 -- Unadapted-tier protocol scrutiny (triangle first) [REGISTERED 2026-07-07, design level]

The E8 verdicts that place genuinely SCALING eigenvectors on the
unadapted tier (triangle D2 = 0.27/0.36; superellipse p=10 D2 =
0.22-0.50) were measured in a same-domain representation basis -- the
protocol family E14 showed manufactures RP-like scaling on the rectangle.
Until an E14-style comparison runs on an unadapted geometry, "genuine RP
phase in free plates" remains unverified everywhere. Design: C0-IP P4 on
the free equilateral triangle (straight edges, so no curved-boundary
complication; SSSS gate = the exact Lame^2 anchor of E3a; FFFF gate =
the certified Argyris E3a instrument), then the matched-window
comparison: (i) true-operator windows in the same-domain SS eigenbasis at
increasing N vs (ii) the truncated-operator ladder built in the SAME
basis (the E9 analog: project the free operator onto the first N SS
eigenfunctions, diagonalize rung submatrices). Superellipse p = 10
follows the same design if the triangle verdict warrants it. Runner and
frozen readings deferred until E17 executes (serial discipline);
readings will follow the E14 dichotomy.

**[EXECUTED 2026-07-09 -- true ladders VALID and NEW; truncated column
INVALID by construction; verdict deferred to E18b.]** Gates spectacular
(free two-mesh 1056/1056 FULL; SS vs exact Lame^2 1999/2000; free vs
Argyris 1000/1000). Finding: the triangle's TRUE-operator eigenvectors
delocalize slowly -- E-sector window IPR 0.081 -> 0.063 over
N = 128..1024, slope -0.138 +/- 0.029 (4.8 sigma from flat, D2_true ~
0.14); A-sectors consistent -- the first real departure from the
rectangle (flat -0.01). The truncated comparator collapsed to the
diagonal identity (SS modes are eigenvectors of the same interior
operator; constructor error, documented in RESULTS). E18b registered:
the faithful E9 analog via nested Koornwinder polynomial Ritz (P4
interpolants, FEM-M orthonormalized, sector-classified rung
eigenvectors, SS-basis representation), same frozen dichotomy vs the
valid true ladders.

## E18b -- Faithful truncated protocol [DONE 2026-07-09 -- PROTOCOL ARTIFACT; Gap A closes negatively for ALL tiers]

Dubiner nested-Ritz (p <= 58, 1770 functions; deg<=4 reproduction
1.5e-14; degree-blockwise C3v projection: A1 310 / A2 280 / E 1180).
RESULT: the truncated protocol manufactures near-ergodic scaling in
every sector -- slopes -0.97 / -1.10 / -0.96 (D2_trunc ~ 1.0) vs the
true operator's -0.25 / -0.15 / -0.138: EXAGGERATION RATIOS 3.9x / 7.5x
/ 7.0x. Frozen gates: A2 + E = PROTOCOL ARTIFACT, A1 = INTERMEDIATE
(flat-bar below its noisy true slope; same ratio story). E-sector
textbook: trunc starts ABOVE true at N = 128 (0.140 vs 0.081), crashes
below by N = 1024 (0.021 vs 0.063). Superellipse contingency NOT
triggered. CAMPAIGN STATEMENT: Gap A closes protocol-negative across
the whole tier hierarchy -- rectangle flat (E17q), triangle weakly
delocalizing (4.8 sigma) but protocol-exaggerated ~7x; the artifact's
mechanism is visible in the Ritz quality itself (corner-singular free
modes are poorly captured by smooth truncated trial spaces -- A2 Ritz
error ~2e-2 at p = 58). Note: the Ritz-sanity check as coded is
rigid-shifted in A1/E (constant/linear content; windows unaffected).
