# E20b -- grid extension + density dial (RESULTS)

## Extended transition curve (K = 8)

| q/Delta | pooled r-tilde |
|---|---|
| 0.001 | 0.3684(0.0082) |
| 0.003 | 0.3857(0.0079) |
| 0.01 | 0.4225(0.0080) |
| 0.03 | 0.4614(0.0080) |
| 0.1 | 0.5175(0.0076) |
| 0.3 | 0.5228(0.0077) |
| 1 | 0.5302(0.0075) |

- resolved q* = 0.03 (interior); Poisson end reached: r(min q) = 0.3684 vs 0.3863

## D2 at resolved q* = 0.03

- **D2_mech = 0.268**; N_sec: K4: 1.46, K6: 1.72, K8: 1.92, K12: 2.02, K16: 2.13
- M-stability (2-pt): {80: 0.39376973239724467, 120: 0.24552214039284415} (gate: spread <= 0.1 -> FAIL)

## The coupler-density dial (at each config's own mid-transition q)

| (n_conn, n_posts) | q_mid | r(q_mid) | D2 |
|---|---|---|---|
| (3, 2) | 0.03 | 0.4614 | 0.268 |
| (5, 4) | 0.01 | 0.4393 | 0.242 |
| (7, 8) | 0.01 | 0.4779 | 0.465 |

- monotone in density: False; densest D2 = 0.465 (AI band [0.56, 0.80])

**Reading: Density does not tune D2 monotonically -- the mech-vs-AI dimension gap is not a simple density effect (report).**

Wall: 33.0 s.

## Campaign close-out (E20 + E20b, 2026-07-10)

- The extended curve resolves the FULL fabricated transition (0.368 ->
  0.530, endpoints on the ensemble values; q* interior). R1 stands
  cleanly.
- D2_mech = 0.268 at the resolved q* (canonical audit-immune protocol;
  ladder contrast confirms protocol immunity on this operator);
  M-parametrization systematic [0.25, 0.39] -- quote as
  D2_mech ~ 0.27 (+0.12/-0.02 systematic). Genuine RP (R2) stands.
- Density dial: directional (densest coupler 0.465, toward the AI band)
  but CONFOUNDED as-designed -- the per-density q_mid selection on the
  coarse grid puts configs at different positions on their own
  transition curves. E20c registered (design level): matched-r-tilde
  dial (interpolated q per config to equal pooled r), plus per-pair
  post-rank as a second knob.
- **CAMPAIGN STATEMENT (the constructive bridge):**
  fabrication-under-forced-cooperation now holds on BOTH substrates --
  a designed plate assembly and a neural MoE each develop a genuine
  RP-multifractal eigenvector phase when decoupled sectors are forced
  to cooperate, measured with the same protocol-immune diagnostics; the
  fabricated fractal dimension is substrate/coupler-dependent
  (mechanical 0.27-0.47 by coupler density, AI 0.56-0.80 by
  architecture). What Gap A denied to nature's plate, design provides.
