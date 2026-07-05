# plates-rp-fem

Open-source finite-element execution of the numerical program of the paper
*"Why Rosenzweig-Porter? Evanescent Coupling and a Random Matrix Hypothesis
for Intermediate Statistics in Freely Vibrating Plates"* (Ornelas Brand,
Zenodo, v4 series). The paper identifies the calculations that decide its
hypothesis and assigns them to COMSOL "or finite elements"; this repository
performs them with a fully open stack (scikit-fem + scipy) so every number is
reproducible on commodity hardware.

## Scientific context

Freely vibrating thin plates show intermediate (Rosenzweig-Porter) spectral
statistics within each symmetry sector, while simply supported plates are
Poissonian (Lopez-Gonzalez et al.). The paper's hypothesis: evanescent-wave
coupling at free edges builds an effectively dense GOE-like perturbation V on
an integrable base H0. The tests here probe that mechanism: the
boundary-controlled Poisson-to-RP transition (E2), geometry
dichotomies (E3), and the eigenvector statistics that decide dense-vs-banded
V (E4, "Gap A").

## Repository structure

```
src/platefem/          shared, validated FEM + statistics machinery
  kirchhoff.py         forms, assembly, eigensolver, rigid-mode handling
  ritz.py              spectral Legendre-Ritz reference (paper T1 protocol)
  stats.py             spacing ratios, sector classification, accuracy metrics
experiments/
  e01_validation/      E1: stack validation vs exact + Ritz references [PASS]
  e02_protocol_b/      E2: Winkler kappa-sweep, Poisson -> RP transition
docs/
  RESEARCH_PLAN.md     experiment registry with preregistered pass criteria
  CONVENTIONS.md       physical conventions and the accuracy policy
```

Executed experiments are frozen with their RESULTS.md, results.json and
figures; each is self-contained and re-runnable.

## Quickstart

```
conda env create -f environment.yml
conda activate plates-fem
python experiments/e01_validation/step1_validate.py
```

(On Windows, invoke the env's `python.exe` directly rather than `conda run`.)

## Status

| experiment | question | status |
|---|---|---|
| E1 | is the stack statistics-grade accurate? | PASS (2026-07-03) |
| E2 | boundary-controlled Poisson -> RP transition | SUPPORTS, 3.6 sigma (2026-07-03) |
| E3a | free equilateral triangle (C3v, non-separable) | SUPPORTS, 6.8 sigma (2026-07-03) |
| E3b | disk control (Poisson) + free ellipse | disk CONFIRMED Poisson; ellipse CHALLENGES the naive reading (Poisson-like at 1600 modes; corners, not curvature, drive coupling) |
| E3c | disk sector (paper's RP case) | AMBIGUOUS +1.7 sigma (weak intermediate; fits corner-sharpness ordering) |
| E4 | Gap A: eigenvector statistics of the true operator | LOCALIZED-LEANING at N<=324: modes ~1-2 beam products; sparse hybridization, not RP-fractal (could flip at N~2048) |
| E5 | superellipse sweep (novel discriminator) | SEPARABLE-BOUNDARY STEP: +12.2 sigma jump at p=2->3, plateau after; coupling set by boundary adaptedness to separable coordinates, corners secondary |
| E6 | designed contact coupling + P10 control law | eigenvector dichotomy crisp (IPR 0.033 vs 0.19); control law 3.6% in-domain; spacing dichotomy 2.3 sigma (v4 full-rank registered) |
| E7 | Protocol A mixed + nu sweep | mixed config = Levy-separable -> ON Poisson (hierarchy confirmed); coupling material-tunable, matches T1 value-by-value |
| E8 | unified Gap A across geometries + Dq | eigenvector hierarchy mirrors E5: adapted -> flat IPR (ellipse = tiny SPARSE coupling); unadapted -> scaling D2 0.22-0.50; RP-vs-PBRM split, decision at N~2048 |
| E10 | annulus (second Poisson control) | SUPPORTS: Poisson 0.3734+/-0.0076; 4x4 J/Y/I/K determinant FEM-cross-validated at 4.3e-4 |
| E9 | Gap A at the registered N=2048 rung | truncated ladder shows D2 = 0.42-0.50 flat-Dq -- REINTERPRETED by E14 as a protocol artifact |
| E12 | P7 damping -> AI-dagger | v1 single-patch: Poisson-like (weak non-proportionality); v2 distributed running |
| E13 | aspect-ratio sweep | ROBUST (chi2 p ~ 0.22); caveat closed |
| E3c v2 | sector long ladder | SUPPORTS +3.9 sigma -- confirms the published sector direction |
| E14 | Gap A reconciliation (C0-IP true operator) | PROTOCOL ARTIFACT: true-operator IPR flat through N = 512/sector; Gap A must be stated on true-operator windows |
| E11 | Mindlin thick plates (thickness sweep) | SUPPORTS +5.4 sigma: <r> 0.454 -> 0.557 over t = 0.02 -> 0.15, GOE-ward |

## License

Code: MIT. The paper and its Zenodo record are licensed separately (CC-BY-4.0).
