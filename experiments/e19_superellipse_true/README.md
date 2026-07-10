# E19 -- superellipse p = 10: the true-operator genuine-D2 point (preregistered)

Frozen 2026-07-10, before any run. Purpose: the THIRD point of the
genuine-delocalization hierarchy. Measured so far on true operators with
the E18 protocol: rectangle slope -0.01 (flat; E17q), triangle
-0.138 +/- 0.029 (weak genuine delocalization, RP-leaning Dq; E18).
The superellipse p = 10 is the E5/E8 hierarchy's strongest-coupling
member; the hierarchy PREDICTION is |slope(p10)| >= |slope(triangle)|.
This run measures the TRUE slope only (the truncated-artifact side was
established on two geometries; not repeated here).

## Instrument

C0-IP P4 on a polar-structured superellipse mesh (fine granularity;
mirror-symmetric; boundary vertices exactly on the superellipse;
geometry matches E5: a = sqrt(1.6189), b = 1/a, p = 10). Production
(64 rings x 192) ~196k dofs; check (56 x 168) ~148k (h ratio ~1.14).
Same-domain SS eigenbasis (essential w = 0), sectors Z2 x Z2 via the
centered mirror classifier. Free-solve certification threshold 1e-3
(the documented rigid-band floor at >150k dofs); ss 1e-4.

## Gates (declared)

- FREE modes (where the windows live): spacing-grade
  n_use_f = min(two-mesh N* at 0.1, cross-instrument N* vs the E5
  certified Argyris p = 10 spectrum over its 1600 modes).
- SS basis (a representation device, not a spacings source): the
  ORDER-STABILITY gate -- two-mesh N* at spacing_frac 1.0 (eigenvalue
  agreement to within one local spacing preserves the nested ordering);
  reported alongside the strict N* for transparency. This relaxation is
  declared HERE, pre-run, with the rationale: representation fidelity
  needs stable identity/order of basis functions, not spacing-grade
  eigenvalues (the E8 precedent for curved same-domain bases).

## Ladders and readings (frozen)

Per sector (ee/eo/oe/oo), true-operator windows [0.4N, 0.6N) at
N in {128, 256, 512} (coverage permitting), coefficients in the first-N
same-domain SS sector functions, renormalized; window-mean ln IPR;
weighted slope over covered rungs; Dq spread (q = 1.5, 2, 4) alongside.

- THIRD-POINT-CONFIRMED: sector-consistent negative slopes with pooled
  |slope|/se >= 3; report D2_true(p10) and its place in the hierarchy
  (prediction: >= triangle's 0.14).
- FLAT: pooled |slope| < 0.05 -- the hierarchy prediction FAILS on the
  true operator (a finding against the E5/E8 ordering).
- INTERMEDIATE otherwise.

## Budget

Overnight-class: ss solve ~2120 modes @196k dofs (the long pole,
~6-9 h detached), free ~1336 (~3-4 h), checks eigenvalues-only.
Staged/stop-proof as E17q/E18.
