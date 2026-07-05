# E15 Phase A -- machinery built + validated, first crossover probe (checkpoint)

Status: **Phase A DONE (machinery validated); first physics probe run.** Phases B/C
(certified long ladders, full two-axis sweep, per-class F2/F3/F4 verdict) remain.
This note is a progress checkpoint, not the frozen RESULTS.

## What was built

`src/platefem/elastic2d.py` -- in-plane (plane-stress) elastodynamics:
- plane-stress K, mass M, gyroscopic Coriolis form `G0(u,v)=2 rho int(u_x v_y - u_y v_x)`
  (real antisymmetric, Omega=1); rotating Hermitian pencil `(K - w^2 M + i w Omega G0)`.
- `star_basis` radial-map domains `r(theta)=R(1+sum a_k cos(k theta+phi_k))`:
  D_prot (even k only -> C2, no mirror), D_mirror (mirror-aligned), D_chir (odd+even).
- modal reduction + `solve_rotor` (2N companion, scipy.linalg.eig) + point-mass
  mistuning `M_pts` in modal space.
- Reuses `kirchhoff.solve_lowest / solve_modes / split_rigid` (operator-agnostic).

## Validation gates passed (anchor-free), `smoke.py`

| gate | result |
|---|---|
| S1 G0 antisymmetry | `||G0+G0'||/||G0|| = 0.0` (exact) |
| S2 rigid trio (Omega=0) | 3 rigid modes, `|Lam|_rigid = 1.8e-14`, clean gap |
| S3 two-mesh convergence | low-12 elastic `Lam` refine 3 vs 4 agree to **~0.5%** (O(h^2)) |
| S4 pencil realness / pairing | `max_imag ~ 1e-14`, `+-omega pair_err ~ 1e-14` at all Omega |
| S5 star domains | D_prot / D_mirror / D_chir all mesh, 3 rigid each |

(The disk in-plane modes come out as clean `+-m` doublets, as they must.)

## First crossover probe, `probe_crossover.py` (380 modes, refine 4, pooled per domain)

`<r>` on the positive-frequency bulk (refs: Poisson 0.386, GOE 0.531, GUE 0.600):

| domain | c_Om=0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| D_chir  | 0.379 | 0.588 | 0.589 | 0.599 | 0.590 |
| D_prot  | 0.375 | 0.594 | 0.595 | 0.604 | 0.590 |
| D_mirror| 0.401 | 0.606 | 0.587 | 0.583 | 0.599 |

## Reading (honest, with the two subtleties the probe exposed)

1. **The machinery is right and the rotation-induced repulsion is genuine.** All
   pencils are real and `+-omega`-paired to machine precision. Rotation moves `<r>`
   UP to the GUE level (~0.60); since class superposition can only *depress* `<r>`,
   an upward move to GUE cannot be a superposition artifact -- it is genuine
   unitary-class level repulsion induced by the Coriolis coupling. Encouraging
   for P6/P11.

2. **The Omega=0 baseline is Poisson (~0.38), not GOE.** The mildly-deformed stars
   (a_k ~ 0.1-0.18) are **near-integrable at rest** -- the F1 recorded fallback. So
   what this probe shows is a Poisson->GUE-level move, not yet a clean GOE->GUE
   (the So-Anlage-Ott-Oerter analog). To realize a clean GOE->GUE the Omega=0
   baseline must first be made genuinely chaotic (stronger / richer deformation, or
   a Sinai-like inclusion as T5 found necessary for ergodicity).

3. **The protection test (F2) CANNOT be read from these pooled numbers.** D_prot
   (C2) and D_mirror carry TWO symmetry classes each; pooling them corrupts `<r>`.
   That D_prot and D_mirror also reach ~0.60 at delta=0 is therefore NOT evidence
   that "protection fails" -- it is the class-pooling artifact. The clean protection
   test requires **separating the R_pi (and mirror) symmetry classes** and reading
   `<r>` per class. This is exactly the delicate symmetry bookkeeping the whole
   P6 prediction hinges on, and the reason the PLAN's `R_pi*T` protector label
   must be checked against the actual operator (T5 found `R_z(pi)*T` does NOT protect
   a cranked operator; only `sigma_v*T` does -- see the review flag).

## Remaining work (sharpened by the probe)

- **B1**: tune the deformation so the Omega=0 baseline is genuinely GOE (chaotic),
  or adopt F1's shift-vs-own-baseline logic.
- **B2**: implement per-class separation (R_pi parity for D_prot; mirror parity for
  D_mirror) so F2/F3/F4 are read within a single class -- the load-bearing step.
- **B3**: certified long ladders (~1200 modes, two meshes, gate G2) + the exact
  in-plane free-disk (Love/Onoe) anchor `disk_inplane.py` (gate G1, still to build).
- **C**: full two-axis (Omega, delta) sweep, per-class `<r>`, F2/F3/F4 verdict.

Compute so far: ~10 s/domain at 380 modes (refine 4) on the desktop; the full
ladder + sweep is the Phase-B/C budget in the PLAN.
