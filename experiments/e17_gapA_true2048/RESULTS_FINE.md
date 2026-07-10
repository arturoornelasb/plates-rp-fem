# E17q -- Gap A true-operator windows at N = 1024-2048 (FINE-MESH CERTIFICATION RESULTS; production (140,86), honest 1.167x check)

Same preregistered readings as E17 (serial run killed externally post-G2; its gates G1 = 1836/5600 informational and G2 = 1200/1200 preserved in *_serial_partial.json). Four parity sectors solved independently; ladders use gate-covered modes only.

| sector | dofs | elastic | G3 two-mesh N* | rigid ratio | resid |
|---|---|---|---|---|---|
| ee | 193545 | 1472 | 1025/1472 | 1.6e+03 | 9.7e-04 |
| eo | 192984 | 1472 | 1093/1472 | 8.6e+03 | 9.9e-05 |
| oe | 193200 | 1472 | 1064/1472 | 1.0e+04 | 9.8e-05 |
| oo | 192640 | 1472 | 1029/1472 | 1.3e+01 | 8.2e-05 |

| sector | basis | IPR N=256 | IPR N=512 | IPR N=1024 | IPR N=2048 |
|---|---|---|---|---|---|
| ee | sine | 0.1755 | 0.1852 | 0.1807 | |
| ee | beam | 0.6926 | 0.7007 | 0.6649 | |
| eo | sine | 0.1606 | 0.1706 | 0.1741 | |
| eo | beam | 0.6339 | 0.6709 | 0.6623 | |
| oe | sine | 0.1674 | 0.1708 | 0.1755 | |
| oe | beam | 0.6747 | 0.6660 | 0.6600 | |
| oo | sine | 0.1734 | 0.1702 | 0.1717 | |
| oo | beam | 0.7221 | 0.6866 | 0.6694 | |

GOE baselines: N=256: 0.0115, N=512: 0.0058, N=1024: 0.0029, N=2048: 0.0015

- top CERTIFIED rung N = 1024: sine sector-mean IPR = 0.1755; fitted D2_true (N >= 512): -0.011; references: E9 truncated 0.069, flat E4/E14 ~0.174

**Reading: PROTOCOL ARTIFACT THROUGH THE CERTIFIED RANGE (N <= 1024): true-operator windows stay sparse at every certified rung; truncated-ladder RP is an artifact there; N = 2048 certification is COVERAGE-LIMITED by the two-mesh gate (indicative row agrees; separated-check-mesh continuation registered)**
