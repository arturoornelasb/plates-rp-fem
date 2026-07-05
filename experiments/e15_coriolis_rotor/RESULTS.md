# E15 -- Coriolis rotor -> GUE (P6/P11): RESULTS

In-plane (plane-stress) elastodynamics of a free 2D rotor spinning about z
(`platefem/elastic2d.py`), rotating Hermitian pencil
`(K - w^2 M + i w Omega G0) phi = 0`. Structured symmetric polar meshes
(centroid-split quads: sigma_v and R_pi are EXACT dof permutations). nu=0.33,
E=rho=1. Modes certified by `solve_modes`; sigma-parity-adapted modal reduction.

**Reading: SUPPORTS P6/P11 -- the first (in-silico) GOE/Poisson -> GUE crossover
in a mechanical elastic system, WITH the protection theorem demonstrated on the
same operator -- and it CORRECTS the preregistered protector: the surviving
antiunitary is `sigma_v * T`, NOT `R_pi * T`.**

---

## G0: the symmetry-residual gate (decisive; machine precision)

Exact discrete symmetry operators `S` (`build_symop`) checked against the
operators (`sym_residuals`) on the certified meshes (16x48):

| domain | S S = I | S K S = K | S M S = M | S G0 S vs -G0 | S G0 S vs +G0 |
|---|---|---|---|---|---|
| D_mirror (sigma_v) | 0 | 3.3e-15 | 3.9e-15 | **3.9e-15 (ANTI-commutes)** | 2.0 |
| D_prot (R_pi/C2)   | 0 | 3.6e-15 | 4.4e-15 | 2.0 | **4.4e-15 (COMMUTES)** |

- `sigma_v` **anti**-commutes with the Coriolis form G0 -> `sigma_v*T` survives
  cranking (`T`: iOmegaG0 -> -iOmegaG0; `sigma_v`: G0 -> -G0; net invariant),
  `(sigma_v T)^2 = +1` -> **orthogonal class protection.**
- `R_pi` **commutes** with G0 -> `R_pi*T` is NOT preserved -> **no protection.**
  This is the correction the PLAN needs: it labelled D_prot (C2) "protected".

Modal-space corollary (`parity_adapt_reduce`): the reduced `G0m` is EXACTLY
block-OFF-diagonal in sigma_v-parity for D_mirror (off-diag weight frac 1.000)
and EXACTLY block-DIAGONAL in R_pi-parity for D_prot (frac 0.000). A 2x2 model of
the block-off-diagonal (single real coupling) gives linear repulsion beta=1 = GOE.

## F2 / F4: crossover and protection (per-class <r>; refs Poisson 0.386, GOE 0.531, GUE 0.600)

| domain (class) | c_Om=0 | 0.5 | 1.0 | 2.0 | verdict |
|---|---|---|---|---|---|
| D_chir (no symmetry) | 0.363 | 0.599 | 0.597 | 0.615 | -> **GUE** (geometric route, F4) |
| D_mirror (sigma_v*T) | 0.401 | 0.504 | 0.503 | 0.517 | **stays GOE** (protected, F2) |
| D_prot R_pi-even     | 0.379 | 0.557 | 0.596 | 0.579 | -> **GUE** (NOT protected) |
| D_prot R_pi-odd      | 0.365 | 0.559 | 0.618 | 0.593 | -> **GUE** (NOT protected) |

The mirror rotor is pinned to GOE at every rotation speed; the chiral and the C2
rotors reach GUE. (D_prot is read per R_pi class because R_pi block-diagonalizes
the rotating operator; pooling its two classes gives a spurious ~0.45.)

## F3: the two-axis phase diagram -- mistuning breaks the protection (the P6 headline)

D_mirror, mistuned at 4 sigma_v-ASYMMETRIC point masses (couples the parity
blocks -> destroys sigma_v*T). Pooled <r>:

| c_delta \ c_Omega | 0.0 | 1.0 | 2.0 |
|---|---|---|---|
| 0.0 | 0.403 | 0.506 | 0.530 |  <- protected: GOE at all Omega |
| 0.5 | 0.515 | 0.611 | 0.611 |  <- mistuned + rotating -> GUE |
| >=1 | 0.515 | 0.610 | 0.612 |

Both ingredients are necessary and sufficient: mistuning alone (Omega=0) stays
GOE-side (real symmetric); rotation alone (delta=0) stays GOE (protected);
mistuning AND rotation reach GUE. This is 'a chirally mistuned rotating elastic
body crosses to the unitary class', realized on the true operator.

## Numerical integrity (Phase-A gates)

G0 antisymmetry exact; rigid trio |Lam| ~ 1e-14; two-mesh low-Lam ~ 0.5% (O(h^2));
rotating pencil real with +-omega pairing ~ 1e-14; symmetry operators exact to
1e-15. Runs are deterministic (fixed seeds in solve_modes / point placement).

## Honest scope and what remains for full certification

- The Omega=0 baseline is Poisson (~0.38), NOT GOE: the in-plane spectrum is a
  superposition of two wave families (dilatational P + shear S), which depresses
  <r> without coupling (the Bertelsen-Ellegaard mechanism). So the pooled move is
  Poisson->GUE; the clean GOE->GUE contrast lives in the PROTECTED-vs-UNPROTECTED
  comparison (mirror plateaus at GOE; chiral/mistuned reach GUE), which is
  unambiguous. Making the P and S families individually chaotic (stronger/richer
  deformation, or a Sinai inclusion as T5 needed) would lift the baseline to GOE.
- Scale here is moderate: N ~ 300-470 certified modes, one mesh (16x48), a reduced
  c_Omega grid, single realization. Full certification (the PLAN's registered
  design) still to do: ~1200 modes at two meshes (gate G2), the exact in-plane
  free-disk Love/Onoe determinant anchor (`disk_inplane.py`, gate G1), seed
  robustness, prestress (Omega/omega) a-posteriori check, and the disk integrable
  control. The qualitative verdict (crossover + protection + corrected protector)
  is not expected to move; the numbers will tighten.

## Files
`../../src/platefem/elastic2d.py` (machinery) ; `smoke.py` (Phase-A gates) ;
`probe_crossover.py` (first look) ; `probe_f2.py` (F2/F4 + symmetry gate) ;
`probe_f3.py` (F3 two-axis) ; `PHASE_A_NOTES.md` (checkpoint log).
