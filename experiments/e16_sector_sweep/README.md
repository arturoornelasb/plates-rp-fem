# E16 -- sector angle-sweep: protocol-faithful replication of the published COMSOL sector result

Preregistered 2026-07-07 (before any run). Runner: `run_sweep.py` (to be
frozen with this README once executed).

## Source protocol (extracted from the paper, archived in the dossier repo)

Lopez-Gonzalez, "Discovery of avoided crossings in plate vibrations using
COMSOL", COMSOL Conference 2024 Boston (PDF:
`lopez-gonzalez-research-analysis/05_avoided_crossings_comsol/`): the
sector result is a pooled ANGLE-SWEEP statistic -- ~70 out-of-plane modes
per angle, the angle swept (range unpublished; the one concrete angle shown
is 50 deg), modes classified S/A by the bisector mirror, 14,000 modes
total, and the pooled spacing-ratio histogram envelope compared to the
Rosenzweig-Porter model. "Sector at the published angles" (the earlier
registry phrasing) is therefore ill-posed: the published object is a pooled
sweep of short low-frequency sequences, NOT a fixed-angle long ladder
(which is what E3c measured: theta = 2.0 rad, 1600 modes, +3.9 sigma).

## Design

- Free Kirchhoff disk sector, R = 1, nu = 0.33 (campaign conventions),
  Argyris instrument with the E3b operational curved-boundary gate.
- Angle sweep: theta in [30, 150] deg, 200 equally spaced values
  (preregistered choice bracketing the published 50 deg; the source does
  not publish its range).
- Per angle: first 70 nonzero elastic out-of-plane modes (rigid modes
  dropped by the multiplicative gap rule); S/A classification by the
  bisector mirror; zero-ambiguity requirement.
- Pooled statistics per class over the sweep, in BOTH ratio conventions:
  P(r) with r = s_n/s_{n+1} (the paper's Fig. 5b convention, for the
  envelope comparison) and r-tilde = min/max (the campaign statistic).

## Preregistered readings

- R1 (replication of the published direction): pooled per-class <r-tilde>
  against the SAME-PROTOCOL pooled Poisson baseline (identical pooling of
  per-angle independent Poisson sequences of matched lengths, 200
  surrogate draws). SUPPORTS if the excess is >= 3 sigma; the RP-envelope
  comparison of P(r) is reported as a fit quality, not gated.
- R2 (protocol scrutiny, the E14 lesson applied to the source): fixed-angle
  references at theta = 50, 90, 114.59 deg (= 2.0 rad, the E3c anchor)
  measured in the MATCHED window (modes 11-70). If the pooled-sweep
  statistic agrees with the fixed-angle matched-window mean within 2 sigma,
  pooling is benign; if not, the sweep protocol itself moves the statistic
  (quantify the bias and its sign).
- R3 (hierarchy consistency): the sector is coordinate-adapted (polar
  boundary lines), so the low-mode sweep should sit in the weak-coupling
  tier: report its position against the E3c long-ladder value (0.4294 +/-
  0.0068) and the disk control (0.391).

## Gates

- Operational curved-boundary gate (established E3b): the FEM disk's class
  statistics reproduce the semi-analytic Poisson control at the matched
  window.
- Spot two-mesh refinement check at three angles (50, 90, 114.59 deg):
  pooled-window statistic stable within 1 sigma.
- Rigid-mode separation by gap (> 4 decades) at every angle; any ambiguous
  S/A label fails the run.

## Budget

70 low modes x 200 angles: small per-angle solves; estimated well under an
hour on the 8-core desktop. Run AFTER E15c frees the CPU.
