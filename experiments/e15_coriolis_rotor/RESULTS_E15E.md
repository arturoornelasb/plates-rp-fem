# E15e -- phase rigidity across the mechanical GOE->GUE crossover (RESULTS)

SVK tangents from the E15d cache; window (0.4, 0.6) of the positive spectrum; R = |phi^T phi| / (phi^+ phi).


## chir

| Omega_nd | <R> (window) | median R | frac R > 0.9 |
|---|---|---|---|
| 0 | 1.0000 | 1.0000 | 1.000 |
| 0.1 | 0.9816 | 0.9922 | 0.967 |
| 0.2 | 0.9499 | 0.9734 | 0.904 |
| 0.3 | 0.9236 | 0.9513 | 0.775 |
| 0.35 | 0.8830 | 0.9401 | 0.683 |
| 0.5 | 0.8607 | 0.9237 | 0.621 |

## mirror

| Omega_nd | <R> (window) | median R | frac R > 0.9 |
|---|---|---|---|
| 0 | 1.0000 | 1.0000 | 1.000 |
| 0.1 | 0.9380 | 0.9849 | 0.829 |
| 0.2 | 0.8657 | 0.9484 | 0.667 |
| 0.3 | 0.7900 | 0.8970 | 0.483 |
| 0.35 | 0.7465 | 0.8451 | 0.388 |
| 0.5 | 0.6386 | 0.7431 | 0.192 |

References (N = 300): GOE <R> = 1.0000, GUE <R> = 0.0728.

**Reading: MIXED: monotone=True, chiral drop 0.14, mirror min <R> 0.64 -- see table.**

Wall: 880.0 s.

## Post-run correction (2026-07-09): the frozen mirror expectation was wrong physics

- What stands: the CHIRAL rigidity falls monotonically (1.0 -> 0.86)
  -- time-reversal breaking is visible on the eigenvector side, the
  qualitative eigenvector face of the crossover exists.
- What was wrong: sigma_v T-invariance forces phi = S conj(phi), i.e.
  real part sigma_v-EVEN and imaginary part sigma_v-ODD -- phi^T phi is
  then REAL (= |a|^2 - |b|^2) but the naive rigidity
  |phi^T phi|/|phi|^2 = ||a|^2-|b|^2|/(|a|^2+|b|^2) is NOT pinned to 1.
  The mirror's faster fall is this parity-split complexity, not a
  protection violation (its SPACING statistics stay GOE, as E15
  established at 9-10 sigma).
- The correct protected eigenvector observable is the ANTIUNITARY
  INVARIANCE DEFECT, d = ||phi - S conj(phi)|| / ||phi|| (exactly 0
  under protection at every speed; no analog constraint for the chiral
  rotor). REGISTERED as E15e-b: compute d(Omega) for both rotors from
  the same caches (S available from build_symop) -- expected mirror
  d ~ 0, chiral d = O(1); that pair, not naive rigidity, is the
  class-discriminant eigenvector face.
- The chiral <R> not reaching the linear-GUE reference (0.07) at the
  crossover is expected: the pencil is quadratic and the window
  parametrization differs; the discriminating content is the monotone
  fall plus the (registered) invariance-defect contrast.
