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

## Runner concretization (pre-run, 2026-07-09; readings above UNCHANGED)

- Instrument: C0-IP P4 on a C3v-symmetric barycentric n-subdivision of
  the E3a macro triangle (side 1, centroid origin, vertex on +y);
  production n = 190 (~290k dofs), check n = 166 (h ratio 1.145 -- a
  properly separated check, per the E17 lesson).
- The same-domain SS eigenbasis is realized as the FEM SS eigenmodes
  (essential w = 0; the discrete Dirichlet-Helmholtz-squared problem):
  M-orthonormal, exact discrete inner products -- the E8 protocol. The
  Lame^2 EXACT spectrum gates the SS solve (the anchor the rectangle
  never had); the free solve is cross-gated against the certified E3a
  Argyris values and its own two-mesh check.
- Truncated ladder (the E9 analog): T = Phi_ss^T K_free Phi_ss per
  sector (zero-extended SS eigenvectors; M_trunc = I), rung submatrices
  eigh'd, window [0.4N, 0.6N) IPR in the basis. True-operator windows:
  free-mode coefficients Phi_ss^T M w_free per sector, first-N truncated,
  renormalized, same windows.
- Sectors via the validated E3a C3v classifier (probe operators through
  the C0-IP space's .probes()). Sector fractions ~1/6 (A1), ~1/6 (A2),
  ~2/3 (E) set the coverage: with 2000 SS modes and 1056 free modes the
  ladders are A-sectors {128, 256} and E {128, 256, 512, 1024}
  (gate-covered ranges permitting); rungs beyond coverage report per the
  COVERAGE rules, extension registered.
- Stop-proof staging (E17q pattern): four detached solve jobs
  (free/ss x prod/check), eigenpairs cached to .npz immediately;
  analysis is a separate cheap stage.
- Slope-based operationalization of the frozen dichotomy, per sector
  over shared rungs: slope_true and slope_trunc = d ln IPR / d ln N
  (window-mean, with window-spread standard errors). PROTOCOL ARTIFACT:
  |slope_true| < 0.15 AND slope_trunc < -0.15. GENUINE SCALING:
  slope_trunc < -0.15 AND |slope_true - slope_trunc| <= 2 sigma.
  INTERMEDIATE otherwise. (D2 = -slope reported per sector alongside the
  E8 fixed-N reference levels.)
