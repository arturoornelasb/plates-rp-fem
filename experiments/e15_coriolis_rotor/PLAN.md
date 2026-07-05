# E15 -- Coriolis rotor -> GUE (P6/P11): execution plan (preregistered design)

Goal: the first GOE -> GUE crossover in a mechanical elastic system, with the
CORRECTED mechanism from the repo's deep-research sweep and toy pretest
(p11_pretest.py, executed and favorable): uniform rotation alone breaks the
naive T but CONSERVES the antiunitary signature symmetry R_pi * T, so clean
rotors stay ORTHOGONAL for all Omega; reaching the unitary class requires
ALSO breaking the signature -- mistuning or geometric chirality. The toy
confirmed both faces on an exact block model (protected: <r> = 0.52 at all
Omega; mistuned: -> 0.59). E15 realizes the full two-axis phase diagram
(Omega, delta) on a true elastic continuum.

## 1. Physical system (why in-plane)

The flat Kirchhoff plate is Coriolis-BLIND for rotation about its normal
(Omega z-hat x w-dot z-hat = 0): no first-order gyroscopic term. The correct
minimal continuum is IN-PLANE (plane-stress) elastodynamics of a free 2D
rotor spinning about z: velocities are in-plane, so the Coriolis term is a
genuine first-order gyroscopic bilinear form

    G(u, v) = 2 rho Omega int (u_x v_y - u_y v_x) dA      (real antisymmetric),

and the rotating eigenproblem (K - omega^2 M + i omega G) phi = 0 is a
HERMITIAN quadratic pencil with real frequencies (stable gyroscopic
system) -- ordinary <r> statistics apply (GOE 0.5307, GUE 0.5996; separation
~ 0.07, so ~1200 bulk modes give ~8 sigma discrimination per grid point).
Centrifugal terms (spin softening -Omega^2, prestress geometric stiffness)
are REAL SYMMETRIC (class-preserving) and of order (Omega/omega_n)^2; the
statistics crossover happens at Omega ~ mean spacing / coupling, which at
mode n ~ 1000 is << omega_1 -- so v1 justifiably defers prestress (recorded
approximation, checked a posteriori via the Omega/omega ratio).

## 2. Symmetry design (the heart of the experiment)

Domains via the radial-map star machinery (existing): r(theta) = 1 +
sum_k a_k cos(k theta + phi_k).

- D_prot ("protected"): C2-symmetric chiral star -- EVEN harmonics only
  (k = 2, 4, 6, generic phases, no mirrors). Possesses R_pi, no reflections.
  Signature-protected: rotation alone must KEEP GOE at every Omega.
- Mistuning knob delta: point masses at asymmetric positions (breaks C2),
  projected in modal space via probes (same machinery as E6/E12):
  delta * M_pts, M_pts,ij = sum_k phi_i(x_k) . phi_j(x_k).
- D_chir ("fully asymmetric"): odd + even harmonics, no C2, no mirrors --
  geometric route to GUE without added masses (cross-check).
- D_mirror: one mirror axis retained -> sigma_v * T protection -> GOE at all
  Omega (second protection control).
- Disk: integrable control (m good; per-m splitting; pooled Poisson-ish).

## 3. Numerical architecture (all proven components)

- platefem/elastic2d.py: P2 vector Lagrange (first derivatives only -- fully
  skfem-native, no hessian issue), plane-stress K (E = 1, nu = 0.33,
  rho = 1), M, gyroscopic form G, star_basis radial map. Free edges =
  traction-free = natural; 3 rigid modes at Omega = 0 (gap detection).
- Anchor (G1): classical in-plane free-disk determinant (Love/Onoe: J_m 2x2
  in the longitudinal/shear wavenumbers) -> platefem/disk_inplane.py;
  validates elements + forms exactly, like platefem.disk did for bending.
- Omega = 0 certified eigenpairs per domain (solve_lowest / solve_modes are
  operator-agnostic), k ~ 1203 at two meshes (internal gate G2).
- Modal reduction: project G and M_pts onto N = 1200 certified modes; the
  (Omega, delta) sweep is then dense linear algebra: Hermitian pencil
  linearization, 2N = 2400 complex eig ~ seconds per grid point.
  Truncation validated by N = 600 vs 1200 stability (gate G3), plus
  spectrum realness / +-omega pairing (gate G4).
- Sweep: dimensionless c_Omega = Omega <|G_ij|> / mean spacing and
  c_delta = delta <|M_pts,ij|> / mean spacing, grid {0, 0.5, 1, 2, 4} x
  {0, 0.5, 1, 2}; <r> on the bulk (drop lowest 10%).

## 4. Preregistered readings (the toy's faces, on the continuum)

- F1 baseline: D_prot and D_chir at Omega = 0 are GOE-side (in-plane
  elastodynamic chaos). If the baseline lands intermediate instead, verdict
  logic switches to SHIFTS relative to the own baseline (recorded fallback).
- F2 protection (the corrected claim, decisive and falsifiable): D_prot
  with delta = 0 stays at its baseline class within 2 sigma for ALL Omega.
  Failure here would CHALLENGE the deep-research correction itself --
  informative either way.
- F3 crossover (headline): D_prot with delta > 0 and Omega > 0 moves
  monotonically toward GUE, reaching >= 3 sigma above GOE (toward 0.5996)
  when both couplings ~ spacing; delta > 0 at Omega = 0 stays GOE-side
  (real symmetric).
- F4 geometric route: D_chir + rotation (no masses) also -> GUE-ward;
  D_mirror + rotation stays GOE-side (mirror protection).
- Disk control: pooled Poisson-ish at all Omega (integrability preserved).
SUPPORTS P6/P11 = F2 + F3 (+F4 consistent). The two-axis contrast
(protected flat vs mistuned crossing) is the publishable core: first
mechanical GOE -> GUE, with the protection theorem demonstrated on the same
system.

## 5. Phases, budget, risks

- Phase A (build, ~2-3 h): elastic2d.py + disk_inplane.py + smokes
  (disk anchor low modes; rigid trio; G antisymmetry; Hermiticity).
- Phase B (runs, ~2-3 h background): certified Omega = 0 solves (3 star
  domains x 2 meshes + disk), modal projections.
- Phase C (~1 h): sweeps (numpy), statistics, figures, RESULTS + verdict,
  draft section, commit.
Risks: (i) Omega = 0 star baseline not GOE (mitigation: richer harmonics;
fallback reading F1); (ii) modal truncation at strong coupling (gate G3);
(iii) disk-anchor root conventions (validated against tabulated
fundamental modes); (iv) whirl behavior of the low modes at high c_Omega
(dropped with the low cut, monitored).

## 6. Relation to the program

P11's framing: the second independent plate <-> AI bridge, on the Dyson-class
axis (the naive version fails in both systems for the same structural
reason: a surviving antiunitary symmetry; a specific ingredient -- mistuning
/ causal mask -- breaks it). The bladed-disk FE re-analysis remains the
"medio-solo" experimental variant for the collaboration; E15 is the
self-contained in-silico realization.
