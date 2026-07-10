# E20 -- mechanical mixture-of-experts (RESULTS)

## R1: the fabricated transition (K = 8)

| q/Delta | lambda_eff | pooled r-tilde |
|---|---|---|
| 0.03 | 0.241 | 0.4614(0.0080) |
| 0.1 | 0.634 | 0.5175(0.0076) |
| 0.3 | 1.618 | 0.5228(0.0077) |
| 1 | 5.069 | 0.5302(0.0075) |
| 3 | 14.765 | 0.5487(0.0074) |
| 10 | 49.383 | 0.5566(0.0074) |
| 30 | 148.407 | 0.5548(0.0074) |

- R1: SUPPORTS (monotone Poisson -> GOE)
- transition coupling q* = 0.03 (r closest to midpoint 0.458)

## R2/R3: sector multifractality at q* = 0.03

| K | mean N_sec | N_sec/K |
|---|---|---|
| 4 | 1.456 | 0.364 |
| 6 | 1.722 | 0.287 |
| 8 | 1.923 | 0.240 |
| 12 | 2.024 | 0.169 |
| 16 | 2.134 | 0.133 |

- **D2_mech = 0.268** (K-sweep, canonical sector protocol); M-stability {80: 0.39376973239724467, 120: 0.24552214039284415}; audit-style truncated ladder on the same operator: -0.645
- R2: GENUINE RP-multifractal phase, fabricated mechanically
- R3: outside the AI band [0.56, 0.80]

**Reading: R1 SUPPORTS (monotone Poisson -> GOE); R2 GENUINE RP-multifractal phase, fabricated mechanically; R3 outside the AI band [0.56, 0.80].**

Wall: 134.6 s.
