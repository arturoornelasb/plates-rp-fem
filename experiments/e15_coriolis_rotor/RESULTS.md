# E15 -- Coriolis rotor -> GUE (P6/P11): RESULTS

In-plane (plane-stress) elastodynamics of a free 2D rotor spinning about z
(`platefem/elastic2d.py`), rotating Hermitian pencil
`(K - w^2 M + i w Omega G0) phi = 0`. Structured symmetric polar meshes
(centroid-split quads: sigma_v and R_pi are EXACT dof permutations). nu=0.33,
E=rho=1. Modes certified by `solve_modes`; parity-adapted modal reduction that
DIAGONALIZES the reduced K (see note).

**Reading: SUPPORTS P6/P11 -- the first (in-silico) GOE -> GUE crossover in a
mechanical elastic system, WITH the protection theorem demonstrated on the same
operator -- and it CORRECTS the preregistered protector: the surviving antiunitary
is `sigma_v * T`, NOT `R_pi * T`.**

---

## G0: the symmetry-residual gate (decisive; machine precision)

Exact discrete symmetry operators `S` (`build_symop`) vs the operators
(`sym_residuals`) on the certified meshes (16x48):

| domain | S S = I | S K S = K | S M S = M | S G0 S vs -G0 | S G0 S vs +G0 |
|---|---|---|---|---|---|
| D_mirror (sigma_v) | 0 | 3.3e-15 | 3.9e-15 | **3.9e-15 (ANTI-commutes)** | 2.0 |
| D_prot (R_pi/C2)   | 0 | 3.6e-15 | 4.4e-15 | 2.0 | **4.4e-15 (COMMUTES)** |

- `sigma_v` **anti**-commutes with the Coriolis form G0 -> `sigma_v*T` survives
  cranking (`T`: iOmegaG0 -> -iOmegaG0; `sigma_v`: G0 -> -G0; net invariant),
  `(sigma_v T)^2 = +1` -> **orthogonal-class protection.**
- `R_pi` **commutes** with G0 -> `R_pi*T` is NOT preserved -> **no protection.**
  This is the correction the PLAN needs: it labelled D_prot (C2) "protected".

Modal corollary (`parity_adapt_reduce`): the reduced `G0m` is EXACTLY
block-OFF-diagonal in sigma_v-parity for D_mirror (off-diag weight frac 1.000)
and EXACTLY block-DIAGONAL in R_pi-parity for D_prot (frac 0.000). A 2x2 model of
the block-off-diagonal (single real coupling) gives linear repulsion beta=1 = GOE.

## F2 / F4: crossover and protection (per-class <r>; refs Poisson 0.386, GOE 0.531, GUE 0.600)

| domain (class) | c_Om=0 | 0.5 | 1.0 | 2.0 | verdict |
|---|---|---|---|---|---|
| D_chir (no symmetry) | 0.560 | 0.584 | 0.604 | 0.594 | GOE -> **GUE** (geometric route, F4) |
| D_mirror (sigma_v*T) | 0.42* | 0.494 | 0.524 | 0.523 | **stays GOE** (protected, F2) |
| D_prot R_pi-even     | 0.533 | 0.542 | 0.597 | 0.569 | GOE -> **GUE** (NOT protected) |
| D_prot R_pi-odd      | 0.510 | 0.557 | 0.636 | 0.592 | GOE -> **GUE** (NOT protected) |

*D_mirror Om=0 is the 2-class (sigma_v-even/odd) pooled value; sigma_v breaks
under rotation so the single-class protected sequence sits at GOE ~0.52. The
mirror rotor is pinned to GOE at every speed; chiral and C2 rotors reach GUE.
(D_prot is read per R_pi class because R_pi block-diagonalizes the rotating
operator; pooling its two classes gives a spurious ~0.45.)

## F3: the two-axis phase diagram -- mistuning breaks the protection (P6 headline)

D_mirror, mistuned at 4 sigma_v-ASYMMETRIC point masses (couples the parity
blocks -> destroys sigma_v*T). Pooled <r>:

| c_delta \ c_Omega | 0.0 | 1.0 | 2.0 |
|---|---|---|---|
| 0.0 | 0.418 | 0.524 | 0.523 |  <- protected: GOE at all Omega |
| 0.5 | 0.513 | 0.591 | 0.617 |  <- mistuned + rotating -> GUE |
| >=1 | 0.514 | 0.591 | 0.615 |

Both ingredients necessary and sufficient: mistuning alone (Omega=0) stays
GOE (real symmetric); rotation alone (delta=0) stays GOE (protected); mistuning
AND rotation reach GUE. This is 'a chirally mistuned rotating elastic body crosses
to the unitary class', on the true operator.

## The Omega=0 baseline is GOE (single-family, mode-conversion-coupled)

The in-plane elastic billiard is chaotic: a moderately rough chiral boundary gives
per-class `<r>_Omega0 = 0.50-0.55` (GOE), stable across roughness levels and under
a strong compact (Sinai) mass inclusion (`baseline_scan.py`). Dilatational (P) and
shear (S) waves are coupled into a single GOE spectrum by mode conversion at the
free edge -- NOT a two-family Poisson superposition. **This is a clean GOE -> GUE
crossover.** (An earlier revision reported a Poisson baseline; that was an artifact
of a modal reduction that took `diag(X'KX)` instead of diagonalizing the reduced
K -- with cluster-mixed certified vectors that discards the off-diagonal coupling
= the level repulsion. Fixed in `modal_reduce`/`parity_adapt_reduce`: both now
eigendecompose the reduced K.)

## Certification (`certify.py`, `disk_inplane.py`)

- **G1 exact anchor (Love/Onoe in-plane free disk).** Semi-analytic 2x2
  traction-free determinant (Helmholtz P/S potentials; n=0 torsional branch =
  J_2(beta R)=0, verified 1.6e-14). FEM disk vs semi-analytic Lambda: rel-err
  **median 4.1e-4 (refine 5)**, halving as O(h^2) from refine 4 (1.7e-3) -- the
  elements/forms are validated exactly.
- **Prestress deferral justified.** At the crossover Omega/omega_med = 0.037
  (c_Omega=0.5) to 0.073 (c_Omega=1), so dropped centrifugal terms ~ (Omega/omega)^2
  ~ 1.3e-3 to 5.4e-3 -- negligible.
- **Disk integrable control.** Class-sequence <r> = 0.410 (Poisson-side): rotation
  of an integrable disk keeps m a good quantum number, no crossover.
- **F2 robustness.** Higher resolution (polar 22x72), different harmonic phases:
  D_mirror 0.44 -> 0.52 (stays GOE); D_chir 0.55 -> 0.59 (GOE -> GUE). The
  protected-vs-unprotected contrast is resolution- and seed-stable.

## Numerical integrity

G0 antisymmetry exact; rigid trio |Lam| ~ 1e-14; two-mesh low-Lam ~ 0.5% (O(h^2));
rotating pencil real with +-omega pairing ~ 1e-14; symmetry operators exact to
1e-15; reduced-K eigendecomposition (deterministic). 

## Honest scope and what remains

Verdict (clean GOE->GUE crossover + sigma_v*T protection theorem + corrected
protector) is established and resolution/seed-stable, with the exact anchor (G1),
prestress deferral, and disk control all certified. Open (numbers-only) item: the
full ~1200-mode two-mesh ladder at the PLAN's registered scale; the verdict is not
expected to move.

## Files
`../../src/platefem/elastic2d.py`, `../../src/platefem/disk_inplane.py` (machinery);
`smoke.py`, `certify.py`, `baseline_scan.py`, `probe_f2.py`, `probe_f3.py`,
`probe_crossover.py`, `probe_protection.py`, `PHASE_A_NOTES.md`.
