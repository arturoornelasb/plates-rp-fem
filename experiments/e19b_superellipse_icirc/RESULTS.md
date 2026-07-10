# E19b -- superellipse p = 10 true-operator windows, init_circle mesh (RESULTS)

Gates: {'free_two_mesh': 61, 'free_vs_e5': 668, 'n_e5': 668, 'ss_strict': 26, 'ss_order': 396, 'n_use_f': 668, 'n_use_s': 1040, 'ss_penalty_order': 1040}

| sector | N | true IPR |
|---|---|---|
| ee | 128 | 0.1706 |
| ee | 256 | 0.1311 |
| ee | slope | -0.380 (0.233) |
| eo | 128 | 0.1353 |
| eo | 256 | 0.1175 |
| eo | slope | -0.204 (0.212) |
| oe | 128 | 0.1307 |
| oe | 256 | 0.1230 |
| oe | slope | -0.088 (0.190) |
| oo | 128 | 0.1426 |
| oo | - | COVERAGE-LIMITED (160 free, 250 basis in gates) |

- pooled slope = -0.205 (0.121), 1.7 sigma from flat -> D2_true(p10) = 0.205; hierarchy refs: rectangle ~0.011, triangle 0.138 +/- 0.029
- Dq spreads (D_1.5 - D_4): ee: 0.132, eo: 0.160, oe: 0.182

**Reading: INTERMEDIATE (1.7 sigma; see table).**

Wall (analysis): 52 s.
