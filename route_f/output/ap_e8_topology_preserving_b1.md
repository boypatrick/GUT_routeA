# AP-E8 topology-preserving finite-grid B=1 card

Status: `mechanical_pass_finite_grid_blocker_changed_continuum_fail_closed`

## Main unanchored stationary cases

| N | L | a | E | gradient density | min dot | B rows | R_rms/a |
|---:|---:|---:|---:|---:|---:|:---:|---:|
| 15 | 6.0 | 0.428571 | 75.090097547 | 1.651e-07 | 0.091845 | [1, 1, 1] | 3.0704 |
| 17 | 6.0 | 0.375000 | 72.444298275 | 3.421e-07 | 0.090442 | [1, 1, 1] | 3.3504 |
| 19 | 6.0 | 0.333333 | 70.831336466 | 2.028e-07 | 0.083178 | [1, 1, 1] | 3.6277 |
| 21 | 6.0 | 0.300000 | 69.844122881 | 7.704e-07 | 0.077786 | [1, 1, 1] | 3.8980 |
| 16 | 6.0 | 0.400000 | 73.627393329 | 2.150e-07 | 0.090722 | [1, 1, 1] | 3.2634 |
| 21 | 8.0 | 0.400000 | 73.346860060 | 1.052e-07 | 0.095418 | [1, 1, 1] | 3.5780 |
| 26 | 10.0 | 0.400000 | 73.289478682 | 1.228e-07 | 0.095975 | [1, 1, 1] | 3.8613 |

## Theorem-level result

For each fixed grid and each nonempty admissible B=1 component, the infinite barrier makes bounded sublevels compactly contained in that component. With the coercive breathing-field mass term, the action attains an interior unanchored minimum.

This is a grid-indexed finite-a stationary family, not a proof of Gamma convergence, regulator independence, or a continuum soliton.

Checks: **13/13**.
Physics promotion: `False`.
Portal start: `False`.

## Checks

- `[PASS]` P0_provenance: the canonical AP-E6 continuum seed has the frozen checksum — checksum=81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59
- `[PASS]` P1_action: the subtracted logarithmic barrier is nonnegative, vacuum-flat, convex, and divergent at the floor — b(1)=0.000e+00; b'(1)=0.000e+00; min b''=1.020e+00
- `[PASS]` P1_action: pairwise admissibility gives a uniform no-zero bound for normalized affine interpolation — sqrt((1+3 epsilon)/4)=0.507444578255
- `[PASS]` P1_calculus: the full AP-E8 chart gradient agrees with a centred directional derivative — best relative residual=9.377e-11
- `[PASS]` P1_calculus: the optimizer-only invalid-domain continuation returns a derivative consistent with its energy — best relative residual=2.809e-11
- `[PASS]` P1_calculus: the intrinsic barrier second variation agrees with a centred geodesic second difference — best relative residual=8.159e-09
- `[PASS]` P2_stationarity: every main final field is unanchored, admissible, three-target B=1, and stationary by the direct tangent test — max gradient=7.704e-07; min pair margin=6.779e-02; B rows=[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]
- `[PASS]` P2_sequence: the fixed-box sequence contains four independently re-relaxed spacings — (N,a)=[(15, 0.4285714285714284), (17, 0.375), (19, 0.3333333333333335), (21, 0.2999999999999998)]
- `[PASS]` P2_sequence: the fixed-spacing L=6,8,10 sequence is stationary with controlled finite-volume energy spread — relative energy spread=4.589523e-03
- `[PASS]` P3_regulator_controls: nearby gamma and epsilon choices retain admissible unanchored B=1 stationary representatives — controls=[(0.5, 0.01, 1.6142235821369239e-07, 0.03476695498081362), (2.0, 0.01, 1.7944472821244817e-07, 0.21121993191627725), (1.0, 0.005, 1.4873487957096304e-07, 0.08051904818413147), (1.0, 0.02, 2.7901910638483354e-07, 0.08439859599539702)]
- `[PASS]` P3_basin_control: a deterministic tangent-plus-scalar perturbation returns to a numerically indistinguishable same-energy admissible stationary B=1 representative — seed=20260719; relative energy difference=4.920518612508454e-14
- `[PASS]` P3_smooth_field: the barrier has a stable approximately quadratic power on a compact-supported boundary-compatible profile — slope=2.159102; R2=0.999141; local slopes=[2.335595393897472, 2.1276499577950196, 2.0694530874120183, 2.0439932562383656]
- `[PASS]` P4_controls: the exact vacuum has zero barrier and B=0 — barrier=0.000e+00; B=[0, 0, 0]
