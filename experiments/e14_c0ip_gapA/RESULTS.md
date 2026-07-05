# E14 -- Gap A reconciliation (RESULTS)

C0-IP P4 at 205989 dofs, 2400 certified modes (~600/sector). E9 truncated-ladder reference at N = 512: sine 0.0946, beam 0.1784 (falling); E4 true-operator at N <= 324: sine ~0.17-0.19 (flat).

- G1 (SSSS vs exact): N* = 1099, max relerr 1.06e-03
- G2 (FFFF vs certified Argyris): N* = 1200, max relerr 1.23e-04

| sector | basis | IPR N=256 | IPR N=512 | E9 truncated at N=512 |
|---|---|---|---|---|
| ee | sine | 0.1755 | 0.1855 | 0.0946 |
| ee | beam | 0.6926 | 0.7018 | 0.1784 |
| eo | sine | 0.1606 | 0.1705 | 0.0946 |
| eo | beam | 0.6339 | 0.6705 | 0.1784 |
| oe | sine | 0.1674 | 0.1712 | 0.0946 |
| oe | beam | 0.6746 | 0.6679 | 0.1784 |
| oo | sine | 0.1735 | 0.1701 | 0.0946 |
| oo | beam | 0.7223 | 0.6866 | 0.1784 |

- true-operator sine IPR at N = 512 (sector mean): 0.1743; E9 truncated: 0.0946; E4 flat level: ~0.18

**Reading: PROTOCOL ARTIFACT: true-operator IPR stays flat at the E4 level through N = 512 -- E9's scaling is a property of the truncation protocol, not the plate**
