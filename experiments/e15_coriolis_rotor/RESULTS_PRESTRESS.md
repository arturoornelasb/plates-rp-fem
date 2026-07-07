# E15c -- prestressed rotor: soft-rotor prediction (RESULTS)

Spin softening + geometric prestress stiffness on certified modes (1203 modes, mesh (24, 72)). Strain metric: eps_max = s1 * Omega_nd^2 (s1 measured per domain). Refs GOE 0.5307 / GUE 0.5996.


## chir (s1 = 1.719)

| Omega_nd | strain % | pooled <r> (prestressed) | pooled <r> (bare) | top-third <r> (prestressed) | max_imag |
|---|---|---|---|---|---|
| 0 | 0.0 | 0.5147(0.0078) | 0.5147 | 0.5230 | 0.0e+00 |
| 0.05 | 0.4 | 0.5232(0.0079) | 0.5157 | 0.5399 | 8.0e-14 |
| 0.1 | 1.7 | 0.5247(0.0077) | 0.5181 | 0.5171 | 6.7e-14 |
| 0.15 | 3.9 | 0.5282(0.0077) | 0.5215 | 0.5076 | 5.2e-14 |
| 0.2 | 6.9 | 0.5466(0.0075) | 0.5259 | 0.5334 | 3.7e-14 |
| 0.25 | 10.7 | 0.5479(0.0077) | 0.5308 | 0.5471 | 2.6e-14 |
| 0.3 | 15.5 | 0.5664(0.0072) | 0.5364 | 0.5669 | 2.2e-14 |
| 0.35 | 21.1 | 0.5608(0.0074) | 0.5422 | 0.5589 | 1.9e-14 |

## mirror (s1 = 0.984)

| Omega_nd | strain % | pooled <r> (prestressed) | pooled <r> (bare) | top-third <r> (prestressed) | max_imag |
|---|---|---|---|---|---|
| 0 | 0.0 | 0.4303(0.0084) | 0.4303 | 0.4277 | 0.0e+00 |
| 0.05 | 0.2 | 0.4421(0.0085) | 0.4339 | 0.4342 | 9.4e-15 |
| 0.1 | 1.0 | 0.4488(0.0085) | 0.4423 | 0.4345 | 1.3e-14 |
| 0.15 | 2.2 | 0.4419(0.0083) | 0.4517 | 0.4326 | 8.9e-15 |
| 0.2 | 3.9 | 0.4471(0.0080) | 0.4607 | 0.4605 | 9.9e-15 |
| 0.25 | 6.2 | 0.4566(0.0080) | 0.4688 | 0.4593 | 1.2e-14 |
| 0.3 | 8.9 | 0.4726(0.0079) | 0.4763 | 0.4714 | 9.1e-15 |
| 0.35 | 12.1 | 0.4856(0.0081) | 0.4819 | 0.4729 | 1.0e-14 |

## Soft-rotor prediction (silicone disk R = 0.1 m, E = 2 MPa, c0 = 45.2 m/s)

| Omega_nd | RPM | strain % | pooled <r> | top-third <r> | top window freq |
|---|---|---|---|---|---|
| 0 | 0 | 0.0 | 0.5147 | 0.5230 | 2403 Hz |
| 0.05 | 216 | 0.4 | 0.5232 | 0.5399 | 2405 Hz |
| 0.1 | 431 | 1.7 | 0.5247 | 0.5171 | 2410 Hz |
| 0.15 | 647 | 3.9 | 0.5282 | 0.5076 | 2419 Hz |
| 0.2 | 863 | 6.9 | 0.5466 | 0.5334 | 2431 Hz |
| 0.25 | 1,078 | 10.7 | 0.5479 | 0.5471 | 2447 Hz |
| 0.3 | 1,294 | 15.5 | 0.5664 | 0.5669 | 2465 Hz |
| 0.35 | 1,510 | 21.1 | 0.5608 | 0.5589 | 2486 Hz |

**Reading: PARTIAL: within strain <= 15% the top window reaches <r> = 0.547 at Omega_nd = 0.25 -- crossover onset visible but not completed inside the linear-validity window; hyperelastic model registered.**
