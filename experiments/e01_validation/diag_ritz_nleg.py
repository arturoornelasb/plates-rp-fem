"""Is the Ritz ee-sector converging monotonically in NLEG, and how fast?
Nested even-Legendre bases -> lambda_k(NLEG) must be non-increasing.
Monotone + slowly decaying increments = corner-singularity algebraic convergence
(reference-limited); non-monotone = numerical corruption at high order."""
import sys
sys.path.insert(0, r"<STEP1_PROTOTYPE_DIR>")  # early cross-check vs the step-1 prototype (path neutralized)
import numpy as np
import ritz_reference as ritz

a, b, nu = 1.0, 1.0 / 1.6189043236, 0.33
cache = {}
NLEGS = [48, 64, 72, 96, 108, 120]
vals = {}
for n in NLEGS:
    ev = ritz.sector_spectrum("ee", a, b, nu, n, cache)
    ev = ev[ev > 1e-6 * ev[5]]          # drop piston
    vals[n] = ev

print("lambda_k(NLEG), ee sector, elastic modes 1..5:")
print("NLEG  " + "  ".join(f"{'mode '+str(k+1):>16s}" for k in range(5)))
for n in NLEGS:
    print(f"{n:4d}  " + "  ".join(f"{vals[n][k]:16.8f}" for k in range(5)))

print("\nconsecutive relative shifts (prev -> next), mode 1:")
for i in range(1, len(NLEGS)):
    p, q = NLEGS[i - 1], NLEGS[i]
    d = (vals[p][0] - vals[q][0]) / vals[q][0]
    print(f"  {p:3d} -> {q:3d}: {d:+.3e}   (positive = monotone decrease, variational OK)")

print("\nmode-1 distance to NLEG=120 value:")
for n in NLEGS[:-1]:
    d = (vals[n][0] - vals[NLEGS[-1]][0]) / vals[NLEGS[-1]][0]
    print(f"  NLEG {n:3d}: {d:+.3e}")
