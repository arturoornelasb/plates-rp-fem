# E19c -- superellipse p = 10, refine-7 crisp standalone (RESULTS)

Gates: {'free_vs_e5_order': 1500, 'free_vs_e5_spacing_reported': 224, 'n_e5': 1500, 'free_two_mesh_r6_informative': 668, 'uniform_drift_vs_e5': 0.0003222053585456395, 'ss_penalty_order': 1256, 'n_use_f': 1500, 'n_use_s': 2300, 'ss_sector_order': {'ee': 417, 'eo': 467, 'oe': 473, 'oo': 561}, 'collar_delta': 0.0012692979846251617}

| sector | N | true IPR |
|---|---|---|
| ee | 128 | 0.1700 |
| ee | 256 | 0.1305 |
| ee | slope | -0.382 (0.232) |
| eo | 128 | 0.1350 |
| eo | 256 | 0.1173 |
| eo | slope | -0.203 (0.211) |
| oe | 128 | 0.1312 |
| oe | 256 | 0.1225 |
| oe | slope | -0.099 (0.190) |
| oo | 128 | 0.1426 |
| oo | 256 | 0.1221 |
| oo | 512 | 0.0990 |
| oo | slope | -0.271 (0.090) |

- pooled slope = -0.249 (0.072), 3.5 sigma from flat -> D2_true(p10) = 0.249; hierarchy refs: rectangle ~0.011, triangle 0.138 +/- 0.029; E19b standalone: 0.205 +/- 0.121
- Dq spreads (D_1.5 - D_4): ee: 0.133, eo: 0.158, oe: 0.178, oo: 0.100

**Reading: THIRD-POINT-CONFIRMED: genuine true-operator delocalization, D2_true = 0.249 at 3.5 sigma -- ABOVE the triangle (hierarchy prediction holds).**

Wall (analysis): 306 s.
