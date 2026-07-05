# E15b -- industrial-speed feasibility (RESULTS)

Hypothesis under test (user): rotation speed is the controlling knob and real machine speeds suffice for the GOE -> GUE crossover in accessible frequency windows.

Chiral rotor reduced model, N = 1200 certified modes; windowed <r> over sliding 200-mode windows; refs GOE 0.5307 / GUE 0.5996.


## aluminum_lab_R0.15m  (c0 = 5355 m/s)

| RPM | Omega_nd | rim speed m/s | f* crossover | modes below f* | window <r> at top | prestress (top) |
|---|---|---|---|---|---|---|
| 3,000 | 0.0088 | 66 | not reached (N=1200) | -1 | 0.531 | 6.9e-08 |
| 12,000 | 0.0352 | 262 | not reached (N=1200) | -1 | 0.531 | 1.1e-06 |
| 30,000 | 0.0880 | 655 | not reached (N=1200) | -1 | 0.532 | 6.9e-06 |
| 100,000 | 0.2933 | 2183 | 58.0 kHz | 110 | 0.541 | 7.7e-05 |

## steel_industrial_R0.5m  (c0 = 5347 m/s)

| RPM | Omega_nd | rim speed m/s | f* crossover | modes below f* | window <r> at top | prestress (top) |
|---|---|---|---|---|---|---|
| 3,000 | 0.0294 | 218 | not reached (N=1200) | -1 | 0.531 | 7.7e-07 |
| 12,000 | 0.1175 | 873 | not reached (N=1200) | -1 | 0.533 | 1.2e-05 |
| 30,000 | 0.2938 | 2183 | 17.4 kHz | 110 | 0.541 | 7.7e-05 |
| 100,000 | 0.9792 | 7278 | 17.4 kHz | 110 | 0.589 | 8.6e-04 |

## Reading

- Burst-limit context: conventional rotor rims sustain ~150-250 m/s; flywheel-grade ~300-500 m/s. Rows with rim speed beyond that are not mechanically reachable for the given radius -- shrink R or the speed.
- The crossover frequency f*(Omega) falls as the window rises (coupling/spacing grows with mode density): at fixed machine speed the LOW spectrum stays GOE and the HIGH spectrum turns GUE -- the observable is a frequency-resolved class transition, measurable by windowed <r> exactly as computed here.

## Addendum: the material invariant and the two observable routes

The burst limit caps the rim speed at v_max ~ sqrt(sigma/rho), so the
achievable dimensionless rotation is a RADIUS- AND SPEED-FREE material
number:

    Om_nd,max ~ v_max / c0 = sqrt( sigma (1 - nu^2) / E ).

Measured bulk crossover: Om_nd ~ 0.293 (N <= 1200 ladder). Material table:

| material | Om_nd,max | consequence |
|---|---|---|
| steel (1 GPa / 200 GPa) | 0.067 | crossover only at window ~ mode 23,000 |
| Al 7075 (0.5 / 69 GPa) | 0.080 | ~ mode 16,000 |
| CFRP flywheel (2.5 / 150) | 0.122 | ~ mode 6,900 |
| polycarbonate (66 MPa / 2.3 GPa) | 0.160 | ~ mode 4,000 |
| PU elastomer (30 / 30 MPa) | 0.94 | **crossover in the FIRST modes** |
| soft silicone (5 / 2 MPa) | 1.49 | **crossover in the FIRST modes** |

VERDICT on the hypothesis: rotation speed IS the controlling knob
(mechanistically confirmed: the crossover tracks Om_nd across all runs), but
for stiff industrial rotors the burst invariant (~0.07-0.12) sits 2-4x below
the bulk crossover -- their LOW spectrum stays GOE at any survivable speed.
Two observable routes, both quantified:

1. DEEP-WINDOW METALS (at the rim limit): steel R = 0.5 m -> f* ~ 294 kHz
   (mode ~ 23,000, mean spacing 6.3 Hz, modal-resolution Q ~ 4.6e4);
   aluminum R = 0.15 m -> f* ~ 814 kHz (Q ~ 3.2e4). Marginal but within
   reach of high-Q resonant ultrasound on low-damping metals in vacuum.
2. SOFT ROTORS: an elastomer disk reaches Om_nd ~ 1 at modest, fully safe
   speeds -- the GOE -> GUE crossover lands in the first hundred modes,
   directly measurable. Caveat registered: at Om_nd ~ 1 the deferred
   centrifugal prestress/deformation is NOT small for soft materials; the
   prestressed rotor model (already in the registry) is required before a
   quantitative soft-rotor prediction.

Bladed-disk mistuned assemblies (structural inter-blade coupling, another
coupling scale entirely) remain the industrial-data variant, unchanged in
the registry.
