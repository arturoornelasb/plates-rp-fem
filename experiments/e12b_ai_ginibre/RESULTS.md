# E12b -- AI-dagger vs Ginibre fine call (RESULTS)

4 patch realizations x 4 sectors, 295/sector; baselines n = 1200 x 24 reps.

| reference | <|z|> | -<cos theta> |
|---|---|---|
| Poisson2D | 0.6665 +/- 0.0020 | -0.0101 +/- 0.0071 |
| GinUE | 0.7403 +/- 0.0019 | 0.2484 +/- 0.0056 |
| AI_dagger | 0.7242 +/- 0.0019 | 0.2057 +/- 0.0044 |

| case | <|z|> | -<cos theta> | n |
|---|---|---|---|
| gamma* = 16 (pooled) | 0.7134 | 0.1174 | 2214 |
| gamma* = 4 (pooled) | 0.6589 | -0.0394 | 2339 |

- at gamma* = 16: d(AI-dagger) = 0.0889, d(GinUE) = 0.1336, d(Poisson2D) = 0.1359; combined marker scale se = 0.0148

**Reading: UNRESOLVED (d_AI/se = 6.0, d_Gin/se = 9.0); registered next = physical-damping refinement, not more statistics**

## Post-run note (2026-07-09)

The sharpened statistics relocate the gamma* = 16 point MID-FLOW on the
Poisson -> AI-dagger crossover rather than on a fixed point: it is
closer to AI-dagger than to GinUE or Poisson2D (ordering preserved from
E12 v3), but 6 se short of saturation -- the v3 "within ~3 sigma"
reflected its coarser marker scale, not saturation. The gamma* = 4 point
(0.659, -0.039) confirms the flow direction. The AI-dagger-vs-Ginibre
ENDPOINT call therefore requires driving the crossover to saturation
(denser/stronger non-proportional dissipation -- a damping-model
refinement, registered), after which the present baselines (separation
0.046 vs se 0.015) make the call trivially resolvable.
