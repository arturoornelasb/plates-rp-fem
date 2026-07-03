# E3c -- free disk sector (preregistered before execution)

The disk sector is the geometry of Lopez-Gonzalez's published COMSOL study
(COMSOL Conference 2024): straight radial edges break the azimuthal
separation of the disk, and the reported statistics approach RP within each
reflection class. This experiment reproduces that setting open-source.

It is also the sharpest available test of the E3 split-verdict
interpretation (corners drive the evanescent coupling): the sector has three
corners (apex plus two arc-radial corners) AND broken separability, so both
readings predict intermediate statistics here -- but a Poisson outcome would
contradict the published COMSOL result and both mechanisms at once.

## Design

- Sector: R = 1, opening angle theta = 2.0 rad (generic, convex, symmetric
  about the x-axis), free edges everywhere. Polar-graded Delaunay mesh: arc
  vertices exactly on the circle, radial edges exactly straight.
- Symmetry: single reflection about the bisector -> classes S (symmetric)
  and A (antisymmetric), classified by the validated mirror-correlation
  machinery (y-mirror only). Rigid modes: piston + x-tilt in S, y-tilt in A.
- Curved-boundary validity: only the arc is approximate; the smooth O(h^2)
  geometry error was shown on the disk (E3b) to leave spacing statistics
  intact (operational gate 0.0 sigma over 1600 modes at comparable h). Gates
  here: internal two-mesh N*, rigid-gap integrity, classifier ambiguity.
- Baseline: matched-window sequences of the semi-analytic free-disk class
  spectrum (E3b machinery).

## Preregistered readings

- SUPPORTS: both classes' <r> >= 3 sigma (pooled) above the matched disk
  baseline -- consistent with the published COMSOL claim and with the
  corner-driver reading.
- CHALLENGES: Poisson-consistent classes (would contradict the published
  sector result and both mechanism readings).
- Windowed trend and per-class values reported regardless.

## Protocol

n_modes = 800 certified eigenvalues; labels from the banded-deflated
eigenvector solve (identity x-mirror trick: 'ee' -> S, 'eo' -> A);
lowest 10 per class dropped; meshes nrings = 90 (production) and 63 (check).
