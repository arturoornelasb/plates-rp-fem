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
| E5 | superellipse corner sweep (novel discriminator) | queued |

## License

Code: MIT. The paper and its Zenodo record are licensed separately (CC-BY-4.0).
