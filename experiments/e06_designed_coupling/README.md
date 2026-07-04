# E6 -- designed contact coupling and inverse design (preregistered)

The paper's contact-coupling test (i), executed in silico with the certified
true-operator mode shapes: build the coupling matrix DIRECTLY from mode-shape
values at contact points,

    V_ij = sum_k phi_i(x_k) phi_j(x_k)        (unit contact stiffness),

assemble H(lambda) = H0 + lambda V within one parity sector (contacts placed
in mirror quartets so the Z2xZ2 classification survives), and sweep lambda
and the contact pattern. This avoids the projection-route infidelity the
paper flags: no basis projection anywhere -- phi are the certified FEM modes
evaluated at the physical points, and H0 is their certified spectrum.

This experiment is also the P10 'programmable coupling' thread: the same
formula, inverted, is a CONTROL LAW -- to hybridize a chosen mode pair
(i, j) with target splitting Delta, place the post quartet at the point x*
that maximizes |phi_i phi_j| on the nodal maps and set
lambda = Delta / (2 q phi_i(x*) phi_j(x*)) with q the quartet multiplicity.

## Patterns (paper's three)

- DENSE: n_c mirror quartets at generic positions -> V dense in the
  frequency-ordered basis.
- BANDED: quartets confined to a narrow strip -> V effectively banded /
  rank-deficient.
- CENTRAL POST: a single center contact -- within the ee sector a rank-1 V
  (symmetry-selection control).

## Preregistered readings (paper, Prediction 'contact')

- SUPPORTS: statistics track the designed structure -- dense toward GOE with
  extended eigenvectors (in the H0 mode basis, which is unambiguous here);
  banded staying on the Poisson side (possibly below, via unshifted
  near-degeneracies) with localized eigenvectors; rank-1 central post blind
  to most of the sector.
- CHALLENGES: dense and banded produce the same statistics, or neither limit
  matches its prediction at any coupling.
- INVERSE DESIGN (P10 demo, beyond the paper): the designed avoided
  crossing's measured splitting must match the closed-form prediction within
  the perturbative regime (target: < 10% relative error), and the targeted
  pair must hybridize ~50/50 while untargeted couplings stay parametrically
  small.

## Protocol

FFFF rectangle (campaign geometry), 96x60 Argyris, 403 certified eigenpairs;
ee sector, N = 100 modes; lambda swept so the mean coupling matrix element
crosses the mean level spacing; IPR measured in the H0 eigenbasis.
