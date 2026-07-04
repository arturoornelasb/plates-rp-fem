# E8 -- unified Gap A across geometries + Dq rival exclusion (preregistered)

Two registered follow-ups executed together:
(1) the E3b decider between 'no coupling' and 'tiny dense coupling' for the
ellipse, and the E4 question of whether the STRONG-coupling geometries
(triangle, superellipse p = 10) have RP-fractal eigenvectors where the
rectangle's are sparse;
(2) the paper's Dq-flatness separator against the power-law-banded rival
(registered calibration: D_{1.5} - D_4 = 0.05--0.07 for RP vs 0.20--0.28 for
PBRM at comparable D_2).

## Method

For each geometry the H0 reference basis is the SIMPLY SUPPORTED eigenbasis
of the SAME domain (reached by the validated kappa = 1e10 Winkler limit),
computed on the same mesh: coefficients are exact M-inner products in the
FEM space, C = Phi_free^T M Phi_SS, restricted per symmetry sector -- no
probes, no quadrature projection, no basis-construction ambiguity. For the
rectangle this basis IS the sine basis, giving a direct cross-check against
E4 (registered E4 sine values: IPR ~ 0.17--0.19 flat, ladder D2 ~ 0.03--0.12).

T1 ladder protocol (window [0.4N, 0.6N), lowest-N SS modes of the sector as
the representation), generalized IPRs P_q = sum |c|^{2q} for
q in {1.5, 2, 4}; ladder D_q from P_q ~ N^{-(q-1) D_q}; empirical GOE
baselines at identical N; captured norms reported, coefficients renormalized.

Geometries: rectangle 96x60 (cross-check), triangle refine-7 (C3v sectors,
E-doublet handling), ellipse refine-6, superellipse p = 10 refine-6
(Z2xZ2). 800+ certified modes each, both solves (free and SS) vector-
certified at 1e-4.

## Preregistered readings

- Ellipse decider: mean off-diagonal weight of C distinguishes 'no coupling'
  (C ~ identity up to reordering; IPR_free-in-SS ~ 1) from 'tiny dense
  coupling' (IPR < 1 with small but nonzero spread growing with N).
- Fractality: any geometry with ladder D2 clearly inside (0.2, 0.9) in this
  basis is RP-candidate; flat-IPR (D2 ~ 0) geometries are in the sparse
  regime at accessible N (as the rectangle in E4).
- Rival exclusion: where D2 is intermediate, D_{1.5} - D_4 < 0.15 favors RP
  (fractal, flat Dq); > 0.15 favors the multifractal PBRM rival.
