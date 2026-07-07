# E17 -- Gap A: true-operator windows at N = 1024-2048 (preregistered)

Preregistered 2026-07-07, before any run. The registered continuation of
E14: bound where, if anywhere, genuine eigenvector scaling begins on the
TRUE operator, over the full ladder range the paper registers
(N = 256..2048 per sector). E14 established flatness through N = 512; E9's
truncated-operator ladder falls to sine IPR ~ 0.052 by N = 2048
(D2 = 0.42-0.50). This experiment closes the remaining range.

## Instrument

Same C0-IP P4 instrument as E14 (platefem/c0ip.py), scaled: mesh
(216, 133) -> ~461k dofs, n_modes = 5600 certified (~1400/sector; the
N = 2048 window [819, 1229) per sector requires >= 1229 modes in the
smallest sector). Projection bases per the T1 registered convention:
sine + free-free beam products, n_funcs_axis = 96 (>= 2048 products per
parity sector), Gauss grid (384, 240) (4x oversampling per axis, as E14).
Beam functions are ratio-form + Loewdin-orthonormalized (bases.beam_1d),
stable at this order. Projection chunked over modes to bound memory.

## Gates

- G1: SSSS (essential w = 0) vs the exact spectrum at the production mesh
  (report N*; informational, as E14).
- G2: FFFF vs certified Argyris 128x80 over the first 1200 modes
  (instrument cross-check, as E14; threshold: N* = 1200/1200 at
  spacing_frac 0.1).
- G3 (decisive accuracy gate): FFFF two-mesh eigenvalue check --
  production (216, 133) vs check mesh (184, 114), eigenvalues-only on the
  check mesh, N* at spacing_frac 0.1 over all 5600. Ladders use only modes
  within min(G3 N*, sector coverage); if that cuts below the N = 2048
  window, the run reports the largest covered rung and the continuation
  mesh is registered.

## Preregistered readings

References: E9 truncated ladder at N = 2048: sine IPR ~ 0.052 (falling,
D2 = 0.42-0.50, flat Dq). Flat true-operator level (E4/E14): sine ~ 0.17
(E14 sector-mean 0.174 at N = 512); beam 0.63-0.72 (modes ~ 1-2 beam
products). GOE baselines computed at each rung (T1 convention).

Let I2048 = sector-mean sine window IPR at N = 2048 on the true operator:

- I2048 >= 0.14: PROTOCOL ARTIFACT THROUGH THE FULL REGISTERED RANGE --
  the plate's true eigenvectors are sparse at every registered rung;
  E9-style ladders (including the P12 recipe on simulation data)
  manufacture the RP phase everywhere in the registry range.
- I2048 <= 0.09: SCALING ONSET -- genuine delocalization begins inside
  the registered range; report the fitted D2_true (slope of mean-ln-IPR
  vs ln N over N = 512..2048) and the onset rung; RECONCILED-RP if
  D2_true falls in E9's 0.42-0.50 band.
- Otherwise INTERMEDIATE: report D2_true, the rung where flatness breaks,
  and the beam-basis picture (does the ~1-2-beam-product structure
  persist?).

Secondary (ungated, reported): beam-basis ladders + capture fractions at
each rung; per-sector spread; the Dq flatness of any scaling segment.

## Budget

Overnight run: certified solve of 5600 modes at 461k dofs (est. 4-10 h) +
check-mesh eigenvalues + chunked projections. Peak RAM ~ 25 GB (V is
~21 GB). Run serially AFTER E15c and E16; nothing else on the machine.
