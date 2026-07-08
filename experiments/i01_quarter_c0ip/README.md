# i01 -- quarter-plate C0-IP instrument (validation gates)

Instrument upgrade, built 2026-07-08 while E17 runs. Purpose: the
rectangle's four parity sectors decouple exactly, so the full-plate
eigenproblem (461k dofs x 5600 modes, ~46 h serial in E17) can be solved
as four independent quarter-plate problems (~115k dofs x ~1400 modes
each, parallelizable across processes; orthogonalization cost ~ modes^2
drops ~16x). E4 registered quarter-plate reduction as the intended
instrument for high Gap-A rungs; this implements it for the C0-IP
discretization.

## Method

Per-edge boundary conditions in platefem/c0ip.py:
- free: natural (existing);
- ss (w = 0): essential dofs, now selectable per edge
  (`boundary_dofs(space, facet_subset)`);
- **guided** (symmetry line: w_n = 0, shear natural): NEW -- one-sided
  Brenner-Sung/Nitsche facet block on the listed boundary facets (jump ->
  one-sided trace, average -> one-sided value, same sigma k^2/h_e
  penalty).

Sector map on the quarter [0, a/2] x [0, b/2] (parities about the two
center lines): even -> guided on that center line; odd -> ss. Rigid-mode
counts per sector (FFFF outer): ee 1 (piston), eo 1, oe 1 (tilts), oo 0
-- asserted exactly by the smoke.

## Gates (smoke_quarter.py; small meshes, matched h)

- GATE A (BC exactness): the four quarter sectors (SS outer edges),
  MERGED, must reproduce the same-mesh full-plate SSSS spectrum at the
  campaign's statistics-grade criterion expressed on the merged ladder
  (err < 0.1 x sector mean spacing = 0.4 x merged spacing; the merged
  ladder is 4x denser), over AT LEAST the range where the full
  instrument is itself accurate against the exact spectrum (4x the
  per-sector informational N* at the smoke h) -- beyond that range both
  instruments are outside their own validity and the comparison is
  unconstrained. The per-sector comparison against the EXACT spectra is
  reported as informational (plain mesh error at the deliberately coarse
  smoke h). Note: the full tensor mesh's uniform triangle-diagonal
  breaks mirror symmetry at the discrete level, so quarter-vs-full
  agreement at the top of the range is bounded by the FULL instrument's
  asymmetry; the quarter problems are exactly symmetric by construction.
- GATE B (cross-instrument): quarter with FREE outer edges vs the full
  free plate solved with the SAME C0-IP instrument at matched h,
  per-sector via the mirror classifier. Requires N* = full coverage over
  the first 48 elastic modes/sector + exact rigid-count match.

PASS on both gates promotes this to the production instrument for E18
and future Gap-A rungs (per-sector processes in parallel).
