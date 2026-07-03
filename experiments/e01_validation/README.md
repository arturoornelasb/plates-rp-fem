# E1 -- Stack validation [FROZEN, executed 2026-07-03, PASS]

Validates Morley (machinery) and Argyris (production) elements for the
Kirchhoff rectangle against (A) the exact SSSS spectrum, (B/C) the paper's
Legendre-Ritz FFFF reference, and (D) the exact SSSS limit of the Winkler
edge spring at large kappa. See RESULTS.md for the verdict and
docs/CONVENTIONS.md for the accuracy policy distilled from this experiment.

Files are kept exactly as executed (self-contained: `ritz_reference.py` is
the local copy that `step1_validate.py` imports; `src/platefem/ritz.py` is
the packaged version for later experiments). `diag_ritz_nleg.py` is the
diagnostic that established the ~1e-4 round-off floor of the Ritz reference.

Environment: conda env plates-fem (scikit-fem 12.0.2), Ryzen 7 5800XT,
total wall time ~336 s.
