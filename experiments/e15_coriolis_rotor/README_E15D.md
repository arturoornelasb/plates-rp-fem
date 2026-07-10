# E15d -- the finite-deformation (SVK) rotor: completing the soft-rotor prediction (preregistered)

Frozen 2026-07-09, before any run. E15c (linear prestress) read PARTIAL:
the crossover onset sits inside the <= 15% linear-validity window and the
GOE/GUE midpoint is crossed just outside it (15.5%). E15d carries the
prediction through the finite-strain regime with a St. Venant-Kirchhoff
plane-stress prestate (exact polynomial tangents; neo-Hookean registered
as a refinement ONLY if the verdict proves model-sensitive at the
validity boundary).

## Instrument (validated by e15d_smoke.py, committed)

- Static finite-deformation centrifugal prestate by Newton in the
  rotating frame (total Lagrangian SVK; centrifugal follower load on the
  deformed positions -- its load stiffness is the EXACT -Omega^2 M spin
  softening). Translations projected (pins the deformed center of mass
  to the axis = zero net force identically); rotation needs no
  projection (torque-free load; -Omega^2 load-stiffness eigenvalue).
- Smoke gates passed: K_T(0,0) = K to 3e-16; K_T(0, Om) = K - Om^2 M to
  4e-16; small-Omega prestate matches E15c's linear field to 0.16%;
  quadratic convergence to machine floor through 15% strain.
- Dynamics: linearized rotating-frame pencil about u0(Omega) --
  M w'' + 2 Omega G w' + K_T(u0, Omega) w = 0, modal-reduced on the
  certified Omega = 0 basis, same Hermitian companion solve as E15/E15c.

## Design

Chiral + mirror stars, mesh (24, 72), 1203 certified modes (as
E15c). Sweep Omega_nd = 0 -> 0.50 in steps of 0.05 with Newton
continuation. Strain metric: max principal strain of u0 (engineering
measure, as E15c). SVK validity cap: strain <= 25% (moderate-strain
regime); points beyond reported, flagged.

## Preregistered readings

- COMPLETES: the chiral top-third-window <r> reaches the GOE/GUE
  midpoint (0.565) at some Omega with strain <= 25%, with the mirror
  rotor's protection intact (its pooled <r> stays below the midpoint at
  every speed). The soft-rotor experiment is then fully predicted with
  finite-deformation physics.
- PARTIAL-AGAIN: onset visible but the midpoint is not reached within
  25% strain -> the neo-Hookean refinement (and/or the physical
  large-strain caveat) becomes the registered next step.
- Reported alongside (ungated): the SVK-vs-linear-prestress drift at
  matched Omega (E15c comparison -- how much the finite-deformation
  terms move <r> beyond the linear Kg picture).
