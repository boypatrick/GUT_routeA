# P54PQ-v2 P1 stationary point, Hessian, and spectrum

Status: **14/14 checks passed**; benchmark is `semidefinite_radial_global_vacuum_with_one_tuned_doublet`.

The machine Hessian covers all `328 = 54 + 252 + 20 + 2` real scalar coordinates and includes the new `chi7` invariant.

- full-gradient residual: `3.799e-16`
- gauge orbit: rank `33`; gauge plus physical PQ: rank `34`
- scalar spectrum: `38` zero, `0` negative modes at tolerance `8.500e-07`
- Hessian wall time: `15.51 s`

## Scalar eigenvalue groups

| m^2 | multiplicity | spread |
|---:|---:|---:|
| 6.098395359e-16 | 38 | 3.83e-15 |
| 3.828030653e-02 | 4 | 4.93e-16 |
| 8.888081554e-02 | 1 | 0.00e+00 |
| 1.121882733e-01 | 12 | 3.61e-15 |
| 1.159441884e-01 | 1 | 0.00e+00 |
| 1.866355208e-01 | 6 | 1.03e-15 |
| 2.116296538e-01 | 6 | 1.67e-15 |
| 2.600000000e-01 | 8 | 1.33e-15 |
| 2.722422827e-01 | 12 | 1.28e-15 |
| 3.939174318e-01 | 4 | 1.17e-15 |
| 4.006848071e-01 | 6 | 1.61e-15 |
| 5.400000000e-01 | 3 | 6.66e-16 |
| 5.522408476e-01 | 6 | 8.88e-16 |
| 1.421000000e+00 | 36 | 2.00e-15 |
| 1.470000000e+00 | 50 | 5.11e-15 |
| 1.518802622e+00 | 44 | 3.77e-15 |
| 1.521099814e+00 | 4 | 8.88e-16 |
| 1.572197378e+00 | 44 | 4.44e-15 |
| 1.757669498e+00 | 6 | 2.22e-15 |
| 1.826678073e+00 | 6 | 8.88e-16 |
| 2.058007717e+00 | 12 | 4.44e-15 |
| 2.058009152e+00 | 6 | 4.00e-15 |
| 2.158061727e+00 | 12 | 3.55e-15 |
| 2.833424996e+00 | 1 | 0.00e+00 |

## Complete SM-Casimir spectrum

`3-pair`, `6-pair`, and similar labels denote a representation plus its conjugate in the real Hessian.

| m^2 | SU(3) | SU(2) dim | |Y| | real mult. |
|---:|---|---:|---:|---:|
| 6.098395359e-16 | 1 | 1 | 2.9802e-09 | 2 |
| 6.098395359e-16 | 1 | 2 | 0.5 | 4 |
| 6.098395359e-16 | 3-pair | 2 | 0.166667 | 12 |
| 6.098395359e-16 | 3-pair | 1 | 0.666667 | 6 |
| 6.098395359e-16 | 1 | 1 | 1 | 2 |
| 6.098395359e-16 | 3-pair | 2 | 0.833333 | 12 |
| 3.828030653e-02 | 1 | 2 | 0.5 | 4 |
| 8.888081554e-02 | 1 | 1 | 0 | 1 |
| 1.121882733e-01 | 3-pair | 2 | 0.166667 | 12 |
| 1.159441884e-01 | 1 | 1 | 0 | 1 |
| 1.866355208e-01 | 3-pair | 1 | 0.333333 | 6 |
| 2.116296538e-01 | 3-pair | 1 | 0.333333 | 6 |
| 2.600000000e-01 | 8 | 1 | 1.14045e-28 | 8 |
| 2.722422827e-01 | 6-pair | 1 | 0.666667 | 12 |
| 3.939174318e-01 | 1 | 2 | 0.5 | 4 |
| 4.006848071e-01 | 3-pair | 1 | 0.333333 | 6 |
| 5.400000000e-01 | 1 | 3 | 6.88941e-28 | 3 |
| 5.522408476e-01 | 1 | 3 | 1 | 6 |
| 1.421000000e+00 | 6-pair | 1 | 0.333333 | 12 |
| 1.421000000e+00 | 3-pair | 3 | 0.333333 | 18 |
| 1.421000000e+00 | 3-pair | 1 | 1.33333 | 6 |
| 1.470000000e+00 | 6-pair | 3 | 0.333333 | 36 |
| 1.470000000e+00 | 6-pair | 1 | 1.33333 | 12 |
| 1.470000000e+00 | 1 | 1 | 2 | 2 |
| 1.518802622e+00 | 8 | 2 | 0.5 | 32 |
| 1.518802622e+00 | 3-pair | 2 | 1.16667 | 12 |
| 1.521099814e+00 | 1 | 2 | 0.5 | 4 |
| 1.572197378e+00 | 8 | 2 | 0.5 | 32 |
| 1.572197378e+00 | 3-pair | 2 | 1.16667 | 12 |
| 1.757669498e+00 | 3-pair | 1 | 0.333333 | 6 |
| 1.826678073e+00 | 3-pair | 1 | 0.333333 | 6 |
| 2.058007717e+00 | 6-pair | 1 | 0.666667 | 12 |
| 2.058009152e+00 | 1 | 3 | 1 | 6 |
| 2.158061727e+00 | 3-pair | 2 | 0.166667 | 12 |
| 2.833424996e+00 | 1 | 1 | 0 | 1 |

## Checks

| Group | Check | Result |
|---|---|---|
| stationarity | the declared 328-field background is stationary | PASS |
| stationarity | radial solver recovers the declared nonzero branch | PASS |
| stationarity | the declared orientation is the lowest enumerated radial branch | PASS |
| hessian | the full 328x328 Hessian is symmetric | PASS |
| goldstone | the gauge orbit has exactly 33 broken directions | PASS |
| goldstone | all broken-generator vectors are Hessian zero modes | PASS |
| goldstone | one independent physical PQ Goldstone remains | PASS |
| goldstone | the projected PQ direction is a Hessian zero mode | PASS |
| spectrum | the scalar eigenvalue census is 34 symmetry zeros plus one complex doublet | PASS |
| doublet | the four non-symmetry zero modes form one (1,2,+/-1/2) real multiplet | PASS |
| spectrum | the SM-Casimir irrep ledger covers all 328 real scalar coordinates | PASS |
| vector | vector spectrum has 12 massless and 33 massive generators | PASS |
| hessian | field census is complete | PASS |
| hessian | the newly required chi7 invariant enters the computed Hessian | PASS |
