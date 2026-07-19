# AP-E9 gamma-scaling and triangulation card

Status: `production_scan_continuum_and_hessian_fail_closed`

## Exact scaling diagnosis

For fixed epsilon, the uniform-edge-margin but smooth-vanishing window is `1 <= p < 2`. The matched scan uses `p=1`, so `w=gamma*a=0.4` while `R=gamma*a^2/d^2 -> 0`; this window alone does not prove degree-sector compactness.

| p | fitted smooth slope | expected | R2 |
|---:|---:|---:|---:|
| 0 | 2.159102 | 2.0 | 0.999141 |
| 1 | 1.159102 | 1.0 | 0.997025 |
| 2 | 0.159102 | 0.0 | 0.863298 |

## Re-relaxed cases

| id | roles | N | L | a | p | mesh | offset/a | E | grad | min margin | B | barycentre | centered RMS | pass |
|:--|:--|---:|---:|---:|---:|:--|:--|---:|---:|---:|:--|---:|---:|:--:|
| p1_fixed_box_N15 | matched_p1_fixed_box | 15 | 6.0 | 0.428571 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.79519804 | 2.36e-07 | 7.515e-02 | [1, 1, 1] | 1.383e-01 | 1.29840 | True |
| p1_fixed_box_N17 | matched_p1_fixed_box | 17 | 6.0 | 0.375000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.58692061 | 1.94e-07 | 8.879e-02 | [1, 1, 1] | 1.096e-01 | 1.26128 | True |
| p1_fixed_box_N19 | matched_p1_fixed_box | 19 | 6.0 | 0.333333 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.94370250 | 1.08e-06 | 9.775e-02 | [1, 1, 1] | 1.054e-01 | 1.23291 | True |
| p1_fixed_box_N21 | matched_p1_fixed_box,triangulation_fine_centered | 21 | 6.0 | 0.300000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 74.65624552 | 2.28e-07 | 1.094e-01 | [1, 1, 1] | 1.139e-01 | 1.21145 | True |
| p1_fixed_box_N23 | matched_p1_fixed_box | 23 | 6.0 | 0.272727 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 75.59504464 | 3.21e-07 | 1.237e-01 | [1, 1, 1] | 1.253e-01 | 1.19604 | True |
| p1_fixed_box_N25 | matched_p1_fixed_box | 25 | 6.0 | 0.250000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 76.67205795 | 5.16e-07 | 1.409e-01 | [1, 1, 1] | 1.350e-01 | 1.18660 | True |
| p1_fixed_a_N16_L6 | matched_p1_fixed_spacing,triangulation_centered | 16 | 6.0 | 0.400000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.62739333 | 1.22e-07 | 8.072e-02 | [1, 1, 1] | 2.824e-01 | 1.27446 | True |
| p1_fixed_a_N21_L8 | matched_p1_fixed_spacing | 21 | 8.0 | 0.400000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.34686006 | 1.08e-07 | 8.542e-02 | [1, 1, 1] | 1.251e-01 | 1.42571 | True |
| p1_fixed_a_N26_L10 | matched_p1_fixed_spacing | 26 | 10.0 | 0.400000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.28947868 | 8.65e-08 | 8.597e-02 | [1, 1, 1] | 3.234e-01 | 1.51029 | True |
| p1_fixed_a_N31_L12 | matched_p1_fixed_spacing | 31 | 12.0 | 0.400000 | 1 | uniform_ppp | [0.0, 0.0, 0.0] | 73.27343710 | 1.88e-07 | 8.623e-02 | [1, 1, 1] | 1.286e-01 | 1.55850 | True |
| tri_translated_N16_uniform_ppp | triangulation_translated_start | 16 | 6.0 | 0.400000 | 1 | uniform_ppp | [0.37, -0.23, 0.41] | 73.61492043 | 1.44e-07 | 8.241e-02 | [1, 1, 1] | 2.342e-01 | 1.27648 | True |
| tri_center_N16_five_tet_phase0 | triangulation_centered | 16 | 6.0 | 0.400000 | 1 | five_tet_phase0 | [0.0, 0.0, 0.0] | 69.63022481 | 3.41e-07 | 8.597e-02 | [1, 1, 1] | 3.565e-01 | 1.21021 | True |
| tri_translated_N16_five_tet_phase0 | triangulation_translated_start | 16 | 6.0 | 0.400000 | 1 | five_tet_phase0 | [0.37, -0.23, 0.41] | 69.63022481 | 5.42e-07 | 8.597e-02 | [1, 1, 1] | 3.565e-01 | 1.21021 | True |
| tri_center_N16_five_tet_phase1 | triangulation_centered | 16 | 6.0 | 0.400000 | 1 | five_tet_phase1 | [0.0, 0.0, 0.0] | 69.62087833 | 8.69e-07 | 8.617e-02 | [1, 1, 1] | 3.015e-01 | 1.21277 | True |
| tri_translated_N16_five_tet_phase1 | triangulation_translated_start | 16 | 6.0 | 0.400000 | 1 | five_tet_phase1 | [0.37, -0.23, 0.41] | 69.62087833 | 8.64e-07 | 8.617e-02 | [1, 1, 1] | 3.015e-01 | 1.21277 | True |
| tri_translated_N21_uniform_ppp | triangulation_fine_translated_start | 21 | 6.0 | 0.300000 | 1 | uniform_ppp | [0.37, -0.23, 0.41] | 74.65624552 | 3.80e-07 | 1.094e-01 | [1, 1, 1] | 1.139e-01 | 1.21145 | True |
| tri_center_N21_five_tet_phase0 | triangulation_fine_centered | 21 | 6.0 | 0.300000 | 1 | five_tet_phase0 | [0.0, 0.0, 0.0] | 70.65934985 | 3.47e-07 | 1.009e-01 | [1, 1, 1] | 2.352e-01 | 1.12813 | True |
| tri_translated_N21_five_tet_phase0 | triangulation_fine_translated_start | 21 | 6.0 | 0.300000 | 1 | five_tet_phase0 | [0.37, -0.23, 0.41] | 70.65934985 | 2.56e-07 | 1.009e-01 | [1, 1, 1] | 2.352e-01 | 1.12813 | True |
| tri_center_N21_five_tet_phase1 | triangulation_fine_centered | 21 | 6.0 | 0.300000 | 1 | five_tet_phase1 | [0.0, 0.0, 0.0] | 70.65816864 | 1.72e-06 | 1.010e-01 | [1, 1, 1] | 1.136e-01 | 1.12828 | True |
| tri_translated_N21_five_tet_phase1 | triangulation_fine_translated_start | 21 | 6.0 | 0.300000 | 1 | five_tet_phase1 | [0.37, -0.23, 0.41] | 70.65816864 | 1.56e-06 | 1.010e-01 | [1, 1, 1] | 1.136e-01 | 1.12828 | True |
| uniform_half_cell_start | translated_start_control | 16 | 6.0 | 0.400000 | 1 | uniform_ppp | [0.5, 0.5, 0.5] | 73.62739333 | 3.69e-07 | 8.072e-02 | [1, 1, 1] | 2.824e-01 | 1.27446 | True |
| p0_control_N21 | p0_control | 21 | 6.0 | 0.300000 | 0 | uniform_ppp | [0.0, 0.0, 0.0] | 69.84412288 | 7.70e-07 | 6.779e-02 | [1, 1, 1] | 1.098e-01 | 1.16424 | True |
| p2_control_N21 | p2_control | 21 | 6.0 | 0.300000 | 2 | uniform_ppp | [0.0, 0.0, 0.0] | 79.68835932 | 8.86e-07 | 2.061e-01 | [1, 1, 1] | 1.102e-01 | 1.26864 | True |

## Fail-closed gate

Numerical regulator thresholds: `False`.
Fixed-box strong-L2 equicoercivity: `True`.
Barrier R->0 Gamma limit on the declared dual-cell embedding: `True`.
Positive-R barrier degree compactness: `True`.
Finite-R homogenized barrier density computed: `False`.
Degree-sector/topological-current compactness: `False`.
Full Gamma limit: `False`.
Regulator-stable continuum background: `False`.
Same-action Riemann Hessian allowed: `False`.
Regulated determinant variation allowed: `False`.

The base Dirichlet term closes fixed-box L2 equicoercivity. The p=1 argument additionally proves a uniform edge-gap consequence and smooth-profile recovery scaling. Neither result yet proves degree-sector closure or the full Skyrme liminf/recovery theorem.

Checks: **8/8**.
