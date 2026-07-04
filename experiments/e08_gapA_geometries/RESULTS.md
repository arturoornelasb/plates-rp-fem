# E8 -- unified Gap A across geometries + Dq (RESULTS)

H0 basis = the simply supported eigenbasis of the same domain (exact M-inner-product coefficients). Ladder D_q from P_q ~ N^{-(q-1) D_q}, window [0.4N, 0.6N), GOE baselines at identical N. Registered separator: D_1.5 - D_4 = 0.05--0.07 (RP) vs 0.20--0.28 (PBRM).


## rectangle

| sector | n_free | leakage | mean IPR by N | D2 ladder | +/- | D_1.5 | D_4 | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|---|
| ee | 210 | 0.092 | 0.184 0.176 0.192 0.184 | -0.023 | 0.047 | 0.016 | -0.048 | **0.065** | 0.862 |
| eo | 196 | 0.097 | 0.190 0.182 0.168 | 0.150 | 0.036 | 0.191 | 0.110 | **0.081** | 0.851 |
| oe | 203 | 0.094 | 0.203 0.194 0.193 0.161 | 0.180 | 0.073 | 0.219 | 0.140 | **0.079** | 0.862 |
| oo | 194 | 0.098 | 0.180 0.169 0.181 | 0.005 | 0.093 | 0.033 | -0.010 | **0.044** | 0.852 |
| GOE | - | - | 0.0443 0.0291 0.0203 0.0147 | ~1 | - | - | - | - | - |

## triangle

| sector | n_free | leakage | mean IPR by N | D2 ladder | +/- | D_1.5 | D_4 | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|---|
| A1 | 143 | 0.101 | 0.141 0.123 0.113 | 0.272 | 0.022 | 0.337 | 0.194 | **0.142** | 0.866 |
| E | 537 | 0.104 | 0.135 0.115 0.098 0.091 | 0.361 | 0.029 | 0.464 | 0.209 | **0.255** | 0.778 |
| GOE | - | - | 0.0443 0.0291 0.0203 | ~1 | - | - | - | - | - |

## ellipse

| sector | n_free | leakage | mean IPR by N | D2 ladder | +/- | D_1.5 | D_4 | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|---|
| ee | 211 | 0.084 | 0.398 0.400 0.380 0.396 | 0.019 | 0.032 | 0.042 | 0.011 | **0.032** | 0.877 |
| eo | 199 | 0.087 | 0.405 0.394 0.390 0.390 | 0.034 | 0.010 | 0.058 | 0.022 | **0.036** | 0.874 |
| oe | 202 | 0.086 | 0.376 0.389 0.372 0.369 | 0.025 | 0.028 | 0.037 | 0.023 | **0.013** | 0.878 |
| oo | 191 | 0.090 | 0.411 0.415 0.399 | 0.034 | 0.037 | 0.061 | 0.021 | **0.040** | 0.871 |
| GOE | - | - | 0.0443 0.0291 0.0203 0.0147 | ~1 | - | - | - | - | - |

## superellipse10

| sector | n_free | leakage | mean IPR by N | D2 ladder | +/- | D_1.5 | D_4 | D_1.5 - D_4 | captured |
|---|---|---|---|---|---|---|---|---|---|
| ee | 209 | 0.089 | 0.204 0.162 0.128 0.130 | 0.429 | 0.097 | 0.474 | 0.375 | **0.099** | 0.868 |
| eo | 198 | 0.094 | 0.175 0.120 0.118 | 0.501 | 0.228 | 0.542 | 0.450 | **0.092** | 0.857 |
| oe | 203 | 0.091 | 0.174 0.142 0.133 0.136 | 0.223 | 0.081 | 0.288 | 0.137 | **0.151** | 0.867 |
| oo | 193 | 0.094 | 0.149 0.163 0.122 | 0.229 | 0.281 | 0.288 | 0.152 | **0.136** | 0.857 |
| GOE | - | - | 0.0443 0.0291 0.0203 0.0147 | ~1 | - | - | - | - | - |

## Verdict (preregistered readings)

- rectangle: ee D2=-0.02, eo D2=0.15, oe D2=0.18, oo D2=0.00 -> sparse regime (flat IPR, D2 ~ 0)
- triangle: A1 D2=0.27, E D2=0.36 -> INTERMEDIATE fractality (D2 in (0.2,0.9)); mean D_1.5-D_4 = 0.199 -> PBRM-like (multifractal)
- ellipse: ee D2=0.02, eo D2=0.03, oe D2=0.03, oo D2=0.03 -> sparse regime (flat IPR, D2 ~ 0)
- superellipse10: ee D2=0.43, eo D2=0.50, oe D2=0.22, oo D2=0.23 -> INTERMEDIATE fractality (D2 in (0.2,0.9)); mean D_1.5-D_4 = 0.120 -> RP-like (flat Dq)

## Discussion (post-run addendum)

The eigenvector hierarchy mirrors the E5 spacing hierarchy exactly, which is
the strongest internal consistency check the campaign has produced:

- ADAPTED boundaries: rectangle IPR ~ 0.18 flat (D2 ~ 0; cross-validates the
  E4 registered-sine-basis result method-to-method) and ellipse IPR ~ 0.40
  flat (D2 ~ 0.03). The ELLIPSE DECIDER lands between its preregistered
  poles: not 'no coupling' (IPR would be ~1) and not 'tiny DENSE coupling'
  (spread would grow with N) but TINY SPARSE coupling -- each free-ellipse
  mode hybridizes ~2.5 SS modes independently of N, even fewer than the
  rectangle's ~5.5. Weak sparse hybridization with Poisson-like spacings.

- UNADAPTED boundaries: triangle D2 = 0.27 (A1) / 0.36 (E) and superellipse
  p=10 D2 = 0.22--0.50 -- genuine eigenvector SCALING, the RP-candidate
  regime, precisely where the spacing statistics are strongest. Spacing
  ratios alone could never separate the rectangle (intermediate <r>, sparse
  vectors) from the triangle (intermediate <r>, scaling vectors): this is
  the paper's Gap A doing exactly the job it was designed for.

- RP vs PBRM at these N: split verdict. Superellipse Dq gap 0.09--0.15
  (flat-ish, RP-leaning); triangle 0.14--0.26 (multifractal-leaning). Both
  D2 values sit well below the P12 preregistered 0.76 +/- 0.15, and the
  registered separator calibration was performed at D2 ~ 0.7, so the
  adjudication needs larger ladders (the N ~ 2048 rung) before either label
  sticks. Caveats: triangle A2 lacks ladder rungs (few SS-A2 modes at these
  energies); captured norms ~ 0.78--0.88 (SS-basis tails), renormalized per
  the T1 convention.
