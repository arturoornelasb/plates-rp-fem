"""platefem -- FEM machinery for spectral statistics of Kirchhoff plates.

Validated in experiments/e01_validation (2026-07-03). See docs/CONVENTIONS.md.
"""
from .kirchhoff import (ElementTriArgyris, ElementTriMorley, assemble_plate,
                        boundary_matrix, make_forms, rectangle_basis,
                        solve_lowest, solve_modes, split_rigid, ssss_exact)
from .stats import (R_GOE, R_POISSON, SECTORS, classify_parity,
                    classify_parity_resolved, local_spacing, mean_r, n_star,
                    probe_operators, r_values)
from . import ritz

__version__ = "0.1.0"
