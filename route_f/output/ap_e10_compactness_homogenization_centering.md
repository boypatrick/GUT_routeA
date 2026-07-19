# AP-E10 compactness, cell formula, and translation quotient card

Status: `production_exact_stencil_no_go_fail_closed`

## Exact compensated-compactness verdict

singular concentration is excluded by the L^(4/3) bound, but compensated compactness fails through a nonzero diffuse oscillatory defect.

| cells | a | sup distance to vacuum | forward current | predicted limit | L4/3 ratio |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.125000 | 4.327e-02 | -0.00295055 | -0.00390625 | 0.758227 |
| 12 | 0.083333 | 2.886e-02 | -0.00344925 | -0.00390625 | 0.803097 |
| 16 | 0.062500 | 2.165e-02 | -0.00364242 | -0.00390625 | 0.819115 |
| 24 | 0.041667 | 1.443e-02 | -0.00378678 | -0.00390625 | 0.830656 |
| 32 | 0.031250 | 1.082e-02 | -0.00383861 | -0.00390625 | 0.834714 |
| 48 | 0.020833 | 7.217e-03 | -0.00387605 | -0.00390625 | 0.837620 |

## Finite-R periodic cell problem

the finite-R barrier density is computed, but triangulation dependence survives homogenization and therefore blocks regulator universality.

The exact six/five ratios are `1.33333333333` on rank-one x and `3` on equal x-y.

## Dynamic-centering and quotient cases

| id | N | L | a | mesh | offset/a | E | grad | B | center/a | rms | pass |
|:--|---:|---:|---:|:--|:--|---:|---:|:--|---:|---:|:--:|
| joint_N25_L6 | 25 | 6.000 | 0.250000 | uniform_ppp | [0.0, 0.0, 0.0] | 76.67205795 | 3.04e-07 | [1, 1, 1] | 5.399e-01 | 1.18660 | True |
| joint_N29_L6.5 | 29 | 6.500 | 0.232143 | uniform_ppp | [0.0, 0.0, 0.0] | 77.67954538 | 1.64e-06 | [1, 1, 1] | 6.314e-01 | 1.22405 | True |
| joint_N33_L7 | 33 | 7.000 | 0.218750 | uniform_ppp | [0.0, 0.0, 0.0] | 78.55223095 | 5.34e-07 | [1, 1, 1] | 5.784e-01 | 1.26582 | True |
| joint_N37_L7.5 | 37 | 7.500 | 0.208333 | uniform_ppp | [0.0, 0.0, 0.0] | 78.98784681 | 1.70e-06 | [1, 1, 1] | 1.234e-01 | 1.35422 | True |
| translation_0_N29 | 29 | 6.500 | 0.232143 | uniform_ppp | [0.37, -0.23, 0.41] | 77.67805475 | 2.97e-07 | [1, 1, 1] | 4.824e-01 | 1.22479 | True |
| translation_1_N29 | 29 | 6.500 | 0.232143 | uniform_ppp | [0.5, 0.5, 0.5] | 77.67954538 | 3.72e-07 | [1, 1, 1] | 6.314e-01 | 1.22405 | True |
| translation_2_N29 | 29 | 6.500 | 0.232143 | uniform_ppp | [-0.41, 0.19, 0.33] | 77.67805475 | 3.77e-07 | [1, 1, 1] | 4.824e-01 | 1.22479 | True |
| mesh_five_tet_phase0_center | 29 | 6.500 | 0.232143 | five_tet_phase0 | [0.0, 0.0, 0.0] | 73.73161128 | 1.01e-06 | [1, 1, 1] | 7.396e-01 | 1.09990 | True |
| mesh_five_tet_phase0_translated | 29 | 6.500 | 0.232143 | five_tet_phase0 | [0.37, -0.23, 0.41] | 73.73161128 | 1.23e-06 | [1, 1, 1] | 7.396e-01 | 1.09990 | True |
| mesh_five_tet_phase1_center | 29 | 6.500 | 0.232143 | five_tet_phase1 | [0.0, 0.0, 0.0] | 73.79708494 | 6.26e-07 | [1, 1, 1] | 5.725e-01 | 1.08543 | True |
| mesh_five_tet_phase1_translated | 29 | 6.500 | 0.232143 | five_tet_phase1 | [0.37, -0.23, 0.41] | 73.68438195 | 5.14e-07 | [1, 1, 1] | 6.208e-01 | 1.09266 | True |

## Fail-closed gate

Numerical quotient gate: `False`.
Discrete compensated compactness: `False`.
Finite-R cell density computed: `True`.
Triangulation-independent finite-R density: `False`.
Regulator-stable continuum background: `False`.
Same-action Hessian allowed: `False`.
Regulated determinant variation allowed: `False`.

Checks: **10/10**.
