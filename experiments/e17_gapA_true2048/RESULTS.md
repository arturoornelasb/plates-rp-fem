# E17q -- Gap A true-operator windows at N = 1024-2048 (RESULTS; i01 quarter-plate instrument)

Same preregistered readings as E17 (serial run killed externally post-G2; its gates G1 = 1836/5600 informational and G2 = 1200/1200 preserved in *_serial_partial.json). Four parity sectors solved independently; ladders use gate-covered modes only.

| sector | dofs | elastic | G3 two-mesh N* | rigid ratio | resid |
|---|---|---|---|---|---|
| ee | 120393 | 1408 | 951/1408 | 1.3e+04 | 7.7e-05 |
| eo | 119952 | 1408 | 986/1408 | 6.3e+04 | 9.6e-05 |
| oe | 120120 | 1408 | 983/1408 | 1.5e+05 | 9.3e-05 |
| oo | 119680 | 1408 | 952/1408 | 1.4e+01 | 8.2e-05 |

| sector | basis | IPR N=256 | IPR N=512 | IPR N=1024 | IPR N=2048 |
|---|---|---|---|---|---|
| ee | sine | 0.1755 | 0.1852 | 0.1806 | |
| ee | beam | 0.6926 | 0.7008 | 0.6645 | |
| eo | sine | 0.1606 | 0.1706 | 0.1743 | |
| eo | beam | 0.6339 | 0.6708 | 0.6629 | |
| oe | sine | 0.1674 | 0.1708 | 0.1754 | |
| oe | beam | 0.6747 | 0.6661 | 0.6594 | |
| oo | sine | 0.1734 | 0.1702 | 0.1717 | |
| oo | beam | 0.7221 | 0.6865 | 0.6689 | |

GOE baselines: N=256: 0.0115, N=512: 0.0058, N=1024: 0.0029, N=2048: 0.0015

- top CERTIFIED rung N = 1024: sine sector-mean IPR = 0.1755; fitted D2_true (N >= 512): -0.011; references: E9 truncated 0.069, flat E4/E14 ~0.174
- INDICATIVE N = 2048 (window includes modes beyond the two-mesh gate; pre-clip run): sine sector-mean IPR = 0.1824 -- consistent with the certified flat level

**Reading: PROTOCOL ARTIFACT THROUGH THE CERTIFIED RANGE (N <= 1024): true-operator windows stay sparse at every certified rung; truncated-ladder RP is an artifact there; N = 2048 certification is COVERAGE-LIMITED by the two-mesh gate (indicative row agrees; separated-check-mesh continuation registered)**
