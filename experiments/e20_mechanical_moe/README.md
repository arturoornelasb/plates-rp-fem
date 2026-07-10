# E20 -- the mechanical mixture-of-experts: fabricating the RP phase by design (preregistered)

Frozen 2026-07-10, before any run. The campaign's constructive inverse
of Gap A. Nature's free plate does not carry an RP eigenvector phase
(E14/E17/E18/E19: flat-to-weak true-operator spreads, ladder scaling an
artifact). The AI program's deepest result is that the Poisson -> GOE
transition IS fabricable when decoupled sectors are forced to cooperate
(MoE; and the fabricated transition is GENUINELY RP at the eigenvector
level: sector D2 = 0.76, audit-hardened, architecture band
[0.56, 0.80]). E20 builds the mechanical mirror: K detuned plates
coupled by a DESIGNED contact-spring network (the E6-validated modal
contact machinery) -- a mechanical MoE -- and asks the same two
questions with the same observables.

## Design

- Sectors: K identical-planform free plates (certified Argyris modal
  data, one solve, (96, 60), first M = 100 elastic modes per plate),
  material-detuned (eigenvalue scale factors 1 + 0.05 u_k, mode shapes
  shared) -- the mechanical analog of per-expert random init; the
  uncoupled spectrum is a superposition of K detuned sequences (born
  Poisson).
- Coupling: sparse random contact graph -- each plate connects to
  n_c = 3 random partners with n_p = 2 posts each at random interior
  points; post spring energy q (w_a(x) - w_b(x))^2 assembled modally via
  probes (E6/E12 machinery). Posts at generic points break all
  symmetries. Per-plate coupling structure is HELD FIXED as K grows
  (the confound-free axis, mirroring d2_ai_sector_clean).
- Sweeps: coupling q over a log grid spanning the transition
  (reported as the measured lambda_eff = median inter-block element /
  median local spacing); K in {4, 6, 8, 12, 16} at fixed M for the D2
  measurement; 6 disorder realizations (detuning + graph + posts).
- Observables (the AI paper's, verbatim): pooled r-tilde of the
  mid-window (25%) spectrum; sector participation
  N_sec = 1 / sum_s w_s^2 per mid-window eigenvector (canonical modal
  basis -- the audit-proven protocol; no truncation, no reference-basis
  ladder); D2 from ln mean N_sec vs ln K at the transition coupling.

## Frozen readings

- R1 (the fabricated transition): pooled r-tilde(q) rises monotonically
  from the Poisson baseline to within 2 sigma of the GOE value at
  strong coupling, at every K -- the mechanical Poisson -> GOE under
  forced cooperation. CHALLENGES if it saturates below GOE - 4 sigma.
- R2 (the eigenvector face, Gap-A grade): at the transition coupling
  q* (pooled r-tilde closest to the Poisson/GOE midpoint 0.458), the
  K-sweep exponent D2_mech in (0.2, 0.9) -> a GENUINE RP-multifractal
  phase, fabricated mechanically; >= 0.9 ergodic; <= 0.2
  localized-mimicry.
- R3 (the bridge number): D2_mech inside the AI architecture band
  [0.56, 0.80] -> the fabrication bridge closes QUANTITATIVELY on both
  substrates; outside -> the two fabrications differ in dimension
  (report; still a two-substrate fabrication result if R2 holds).
- Instrument gates (all reported): modal-truncation stability
  (M = 80 vs 120: D2 within 0.1); seed spread; and the audit-style
  ladder contrast at K = 16 (the plate-style truncated protocol on the
  same operator SHOULD manufacture inflated scaling -- protocol-immunity
  demonstration, symmetry of scrutiny).

## Budget

One certified plate solve (~20 min) + modal sweeps (eighs <= 1920^2;
minutes). ~45 min total, detached.
