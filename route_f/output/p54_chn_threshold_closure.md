# CHN finite vertex and sequential operator closure

Checks: 20/20.

**This is not a complete finite seesaw matcher. Existing diagnostic trajectories are not modified.**

The no-half CHN interaction gives M(rho)=M-2 C rho. A direct Yukawa/CHN bubble has finite vertex 2 Y diag[M(1+log(mu^2/M^2))] C/(16 pi^2). Its all-removed linear-CHN C5 and finite Higgs-mass subsets pass independent RG, integral, determinant and covariance checks at leading qH=0.

At a partially removed Takagi block the exact tree Schur map generates Y_eff(rho)=Y_r+D6 rho, D6=2 Y_h M_h^-1 C_hr. The operator -(D6 LHN_r)(HdaggerH)+h.c. has dimension six, but retained Majorana masses allow it to mix down into C5. The explicit mixed-block beta residual equals this missing down-mixing. No hierarchy suppression is silently assumed.

| Actual diagnostic case | Threshold | D6 norm | Mixed CHN norm |
|---|---:|---:|---:|
| 0 | 0 | 0.00011505006 | 1.1635675e-05 |
| 0 | 1 | 7.9314103e-05 | 7.0659491e-06 |
| 0 | 2 | 0 | 0 |
| 1 | 0 | 0.0001189891 | 1.8081431e-05 |
| 1 | 1 | 0.00015040141 | 1.122946e-05 |
| 1 | 2 | 0 | 0 |

The scalar-tree CHN can align with MR at the input, but actual running generates mixed Takagi entries. These numbers use the pre-existing synthetic families, not fitted data.

An additional controlled expansion in retained/removed Majorana masses can consistently place the D6 feedback beyond leading dimension-five accuracy. Keeping retained-mass effects without that expansion instead needs the generated operator or an explicit remainder bound; no fitted-trajectory bound is claimed here.

## Missing

- full dimension-six operator matching/running and its mass-suppressed feedback
- nonzero-qH terms and finite correlated MR,lambda and field maps
- remaining pre-existing-C5/mixed-sterile/box diagrams and lower scalar/gauge matching

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| massive_massless_vertex_integral_0 | 3.225e-16 | True |
| massive_massless_vertex_integral_1 | 2.448e-16 | True |
| massive_massless_vertex_integral_2 | 1.110e-16 | True |
| all_removed_CHN_C5_scale_derivative_beta_difference | 5.881e-17 | True |
| terminal_CHN_C5_family_covariance | 9.905e-20 | True |
| terminal_CHN_Higgs_tadpole_family_invariance | 3.970e-23 | True |
| CHN_finite_Higgs_mass_direct_determinant | 2.404e-18 | True |
| CHN_finite_Higgs_mass_scale_derivative | 1.661e-19 | True |
| degenerate_CHN_vertex_O3_covariance | 6.550e-22 | True |
| partial_CHN_matching_RG_closes_with_D6_1 | 1.882e-16 | True |
| tree_LNH_equals_exact_field_Schur_derivative_1 | 9.326e-12 | True |
| dimension5_only_mixed_block_residual_detected_1 | 0.000e+00 | True |
| partial_CHN_matching_RG_closes_with_D6_2 | 1.963e-16 | True |
| tree_LNH_equals_exact_field_Schur_derivative_2 | 8.516e-12 | True |
| dimension5_only_mixed_block_residual_detected_2 | 0.000e+00 | True |
| partial_CHN_matching_RG_closes_with_D6_3 | 2.484e-16 | True |
| actual_trajectory_three_events_0 | 0.000e+00 | True |
| actual_running_mixed_CHN_generates_D6_0 | 0.000e+00 | True |
| actual_trajectory_three_events_1 | 0.000e+00 | True |
| actual_running_mixed_CHN_generates_D6_1 | 0.000e+00 | True |
