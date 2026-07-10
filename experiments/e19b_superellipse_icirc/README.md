# E19b -- superellipse p = 10 on the E5-validated mesh (preregistered)

Frozen 2026-07-10, before its run. E19's polar-ring instrument failed
its own gates (E5-cross N* = 224; corner-element conditioning). E19b
reruns the SAME frozen protocol and readings (see
../e19_superellipse_true/README.md, unchanged) on the init_circle-mapped
superellipse mesh -- the family E5 quality-measured and certified to
1600 Argyris modes -- at refine 6 (~131k P4 dofs; check at refine 5,
h ratio 2, honest). Ladders scoped to the expected refine-6 coverage:
{128, 256} per sector (slopes from 2 rungs; pooled over 4 sectors).
Solve budgets: free 668 modes (~1 h), ss 1040 (~2 h), checks fast.
All other design elements (gates incl. E5-cross + SS order-stability,
windows, readings THIRD-POINT-CONFIRMED / FLAT / INTERMEDIATE)
inherited verbatim.
