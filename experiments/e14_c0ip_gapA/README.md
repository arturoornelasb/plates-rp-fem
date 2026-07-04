# E14 -- Gap A reconciliation: true-operator windows at large N (preregistered)

THE open computation (task #0): E4/E8 measured flat IPR (sparse regime) on
TRUE-operator eigenvectors at N <= 324/sector, while E9's registered
truncated-operator ladder shows the RP fractal phase developing
(D2 = 0.42--0.50) through N = 2048. Until the protocols are compared at
matching N, Gap A cannot be claimed cleanly.

## Instrument

C0 interior-penalty P4 discretization (platefem/c0ip.py): custom analytic
reference hessians (skfem Lagrange carries none), Brenner-Sung interior
facet terms, free edges natural. Smoke-validated: 3 rigid modes; FFFF
elastic values match certified Argyris to 7e-7 at equal h; eigenvalues
penalty-insensitive across sigma = 5..20; SSSS at 5e-5 on a coarse mesh.
No ElementGlobal Vandermonde -> no h ~ 0.008 ceiling.

## Protocol

Full campaign rectangle, P4 C0-IP at h ~ 0.007 (144x89, ~206k dofs),
n_modes = 2400 certified -> ~600/sector. Gates: G1 = SSSS (essential w = 0)
vs exact spectrum at the production mesh; G2 = FFFF vs certified Argyris
128x80 over the first 1200 modes. Then E4's registered-basis projection
machinery (sine + beam), ladder N in {256, 512} per sector, window
[0.4N, 0.6N), GOE baselines.

## The decisive comparison (preregistered)

At N = 512, window [205, 307): E9's truncated-operator ladder gives window
IPR ~ 0.095 (sine) / 0.178 (beam), falling. E4's true-operator windows at
N <= 324 gave ~ 0.17--0.19 (sine), flat.

- RECONCILED-RP: true-operator window IPR at N = 512 drops toward the
  truncated values (scaling onset is physical; E9's RP reading stands for
  the true plate).
- PROTOCOL ARTIFACT: true-operator IPR stays ~ flat at 0.17--0.19 through
  N = 512 (the truncated ladder's scaling reflects the truncation itself,
  not the plate; E9's reading must be reinterpreted as a property of the
  registered protocol only).
- Either outcome resolves the contradiction and fixes how Gap A is stated
  in the draft.
