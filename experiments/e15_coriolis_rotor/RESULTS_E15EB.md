# E15e-b -- the sigma_v T invariant R_S (RESULTS)

R_S = |phi^T S phi| / (phi^+ phi), S = parity metric; window (0.4, 0.6). Naive R shown for contrast (E15e).

- structure gates: KT opp-parity fraction = 1.76e-12, G0 same-parity fraction = 2.17e-15 -> PASS

| Omega_nd | median R_S | p5 R_S | frac R_S < 0.999 | median R (naive) |
|---|---|---|---|---|
| 0 | 1.000000 | 0.999705 | 0.025 | 1.0000 |
| 0.1 | 1.000000 | 1.000000 | 0.000 | 0.9849 |
| 0.2 | 1.000000 | 1.000000 | 0.000 | 0.9474 |
| 0.3 | 1.000000 | 1.000000 | 0.000 | 0.8921 |
| 0.35 | 1.000000 | 1.000000 | 0.000 | 0.8715 |
| 0.5 | 1.000000 | 1.000000 | 0.000 | 0.7257 |

- BROKEN arm (odd-parity mistuning, RMS element = d, Omega = 0.3): median R_S = 0.2620 (p5 = 0.0682)

**Reading: PROTECTION-EXACT-AND-FALSIFIABLE: the sigma_v T invariant R_S stays at 1 (median >= 0.999) at every speed while the naive R falls, and an odd-parity mistuning collapses it to 0.262 -- the correct eigenvector-level protection observable, measured and falsified in one experiment.**

Wall: 622.3 s.
