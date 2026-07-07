# E18 -- Unadapted-tier protocol scrutiny: the free triangle (preregistered)

README frozen 2026-07-07, BEFORE E17's verdict and before any E18 run
(the runner will be committed before its execution; these readings are
binding). Motivation: E8 placed genuinely SCALING eigenvectors on the
unadapted tier (triangle D2 = 0.27/0.36, superellipse p=10 D2 =
0.22-0.50) -- measured in a same-domain representation basis, i.e. the
protocol family E14 proved manufactures RP-like scaling on the
rectangle. Until an unadapted geometry gets the E14 treatment, "genuine
RP eigenvector phase in free plates" is unverified everywhere.

## Design

- Instrument: C0-IP P4 on the free equilateral triangle (straight edges,
  so no curved-boundary complication); C3v-symmetric mesh by uniform
  refinement of the macro triangle; A1/A2/E sectors via the E3a
  classifier (exact E doublets deduplicated before statistics).
- Gates: G1 = SSSS vs the exact Lame^2 anchor (E3a route); G2 = FFFF vs
  the certified Argyris E3a instrument over shared modes; two-mesh spot
  check at the production mesh.
- The matched-window comparison, per sector:
  (i) TRUE-operator windows: certified C0-IP free modes projected onto
  the same-domain SS (Dirichlet-Helmholtz squared / Lame) eigenbasis at
  increasing truncation N; T1 window [0.4N, 0.6N) IPR ladder.
  (ii) TRUNCATED-operator ladder in the SAME basis (the E9 analog):
  free-operator matrix elements in the first N SS eigenfunctions via
  their FEM interpolants (K_ij = phi_i^T K phi_j, M_ij likewise), rung
  submatrices diagonalized, same window IPR.
- Ladder N in {128, 256, 512} per sector to start (~2400 certified modes
  total); extension rung registered if the verdict needs it.

## Preregistered readings (the E14 dichotomy, applied to the tier where RP could still live)

- PROTOCOL ARTIFACT: the true-operator window IPR stays flat at the E8
  fixed-N level while the truncated ladder falls -> the unadapted tier's
  "scaling" is a representation-protocol artifact; combined with
  E14/E17, Gap A closes negatively for ALL tiers at accessible N.
- GENUINE SCALING: the true-operator IPR falls with N, consistent with
  the truncated ladder's slope within 2 sigma -> the unadapted tier
  carries a real RP-candidate eigenvector phase; report D2_true and the
  registered Dq-flatness (RP vs PBRM) separator.
- INTERMEDIATE otherwise: report both curves and where they part.
- Contingency (registered): the superellipse p = 10 gets the same
  treatment ONLY if the triangle reads GENUINE SCALING.

## Budget

Comparable to E14 (~1-3 h). Runs AFTER E17 completes (serial); the
runner freezes before execution.
