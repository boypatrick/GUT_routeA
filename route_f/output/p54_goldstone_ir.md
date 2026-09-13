# P54 same-action Goldstone / IR background audit

**The one-loop scalar IR and momentum subgate is computed; full gauged background stability is NOT decided.**

Checks: 65/65.

All 34 Goldstone hard mass shifts are canceled by the existing common fixed-VEV counterterms. Adding a positive Goldstone mass by hand is not a same-scheme repair. Actual soft tadpoles vanish as the regulator tends to zero, but the radial potential Hessian has a nonzero positive-Gram coefficient multiplying log(epsilon/mu^2). It is not a finite physical mass matrix.

The full scalar one-loop momentum dependence on the radial slice is computed from the actual cubic/quartic action tensors. External Euclidean momentum replaces the divergent log by log(pE^2/mu^2)-2 for the massless bubble. The hard scalar momentum difference is positive semidefinite, but this is not the full gauged pole kernel.

| pE^2 / omega^2 | Lowest generalized eigenvalue of diagnostic kernel |
|---:|---:|
| 0.001 | -0.758431138 |
| 0.01 | -0.724978561 |
| 0.05 | -0.584213748 |
| 0.1 | -0.413688269 |
| 0.5 | 0.596231807 |
| 1 | 1.09858795 |

The diagnostic includes the old zero-momentum vector hard term, not its dynamical gauge/ghost completion. Its negative values cannot be promoted to a physical tachyon; positive values at high Euclidean momentum cannot establish a stable vacuum either. No default scalar parameters or old matching trajectories were changed.

## Background control and a ruled-out repair ray

The hard scalar kinetic correction has generalized eigenvalues [0.009008482863914856, 0.27686404964260186, 2.610319889177881]. The largest exceeds one in the declared MS convention. A finite canonical field redefinition can absorb this metric, but does not establish a controlled loop remainder or change the inertia of a fixed curvature matrix.

For the exact existing-parameter ray V0 -> z V0 at fixed radii/gauge/scale, H_hard(z)=z H_tree+z^2(S+log(z) L_S)+V and deltaZ_h(z)=z deltaZ_h(1). The sigma coordinate gives an analytic concave-quadratic upper bound for every 0<z<=1:

`H_sigma,sigma(z) <= -0.00388757181 omega^2 < 0`.

This excludes that entire BOSONIC repair ray, not arbitrary scalar parameters or fermionic corrections at new points. The separate scalar_reselection report tests a few less-hierarchical engineering points; none overwrites the old inputs.

## Missing

- same-prescription vector/Goldstone/ghost mixed and fermionic momentum maps with Nielsen consistency
- controlled scalar parameter feedback or an explicitly reorganized perturbative counting
- complete Wilson/box, lower gauge/ghost/Yukawa matching and finite CHN seesaw

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| frozen_tree_stationarity | 2.099e-16 | True |
| 38_tree_soft_modes | 2.584e-15 | True |
| positive_hard_gap | 0.000e+00 | True |
| 33_gauge_1_PQ_4_light_dimensions | 0.000e+00 | True |
| actual_soft_projector_decomposition | 3.898e-14 | True |
| goldstone_tree_zero_modes | 4.152e-15 | True |
| hard_Ward_mass_matrix_symmetric | 9.756e-17 | True |
| common_CT_cancels_all_34_hard_Ward_masses | 7.812e-15 | True |
| actual_cubic_Goldstone_Ward_0 | 6.115e-15 | True |
| actual_cubic_Goldstone_Ward_1 | 2.193e-14 | True |
| actual_cubic_Goldstone_Ward_2 | 3.253e-14 | True |
| unbroken_vectors_stay_massless_along_radial_0 | 1.187e-15 | True |
| unbroken_vectors_stay_massless_along_radial_1 | 3.327e-15 | True |
| unbroken_vectors_stay_massless_along_radial_2 | 0.000e+00 | True |
| gauge_IR_Gram_positive | 2.410e-20 | True |
| PQ_IR_Gram_positive | 0.000e+00 | True |
| Goldstone_IR_Gram_positive | 0.000e+00 | True |
| Higgs_IR_Gram_positive | 2.448e-23 | True |
| all_soft_IR_Gram_positive | 0.000e+00 | True |
| Goldstone_Higgs_radial_mixing_vanishes | 1.710e-15 | True |
| all_soft_IR_Gram_equals_Goldstone_plus_Higgs | 1.133e-16 | True |
| soft_tadpole_vanishes_with_regulator | 0.000e+00 | True |
| IR_log_coefficient_nonzero_on_old_negative_witness | 0.000e+00 | True |
| regulated_Hessian_log_slope_1e-06 | 3.182e-11 | True |
| regulated_Hessian_log_slope_1e-08 | 3.010e-13 | True |
| soft_regulated_Hessian_independent_spectral_difference | 4.438e-07 | True |
| all_scalar_zero_momentum_hard_part_reproduces_old_Frechet | 1.381e-16 | True |
| grouped_soft_bubble_equals_actual_IR_Gram | 1.117e-16 | True |
| bubble_independent_integral_0 | 2.039e-16 | True |
| bubble_independent_integral_1 | 0.000e+00 | True |
| bubble_independent_integral_2 | 7.666e-16 | True |
| bubble_independent_integral_3 | 2.390e-15 | True |
| momentum_IR_cancellation_identity_0.0001 | 0.000e+00 | True |
| momentum_IR_cancellation_identity_1e-07 | 1.540e-16 | True |
| momentum_IR_cancellation_identity_1e-10 | 1.540e-16 | True |
| scalar_hard_momentum_difference_positive_Gram_0.001 | 0.000e+00 | True |
| scalar_hard_momentum_difference_positive_Gram_0.01 | 0.000e+00 | True |
| scalar_hard_momentum_difference_positive_Gram_0.05 | 0.000e+00 | True |
| scalar_hard_momentum_difference_positive_Gram_0.1 | 0.000e+00 | True |
| scalar_hard_momentum_difference_positive_Gram_0.5 | 0.000e+00 | True |
| scalar_hard_momentum_difference_positive_Gram_1.0 | 0.000e+00 | True |
| scalar_momentum_quadrature_64_vs_256 | 9.327e-15 | True |
| hard_scalar_kinetic_matches_momentum_derivative_1e-05 | 4.505e-06 | True |
| hard_scalar_kinetic_matches_momentum_derivative_5e-06 | 2.253e-06 | True |
| hard_scalar_kinetic_positive | 0.000e+00 | True |
| scalar_vector_CT_split_reconstructs_old_radial | 1.172e-16 | True |
| uniform_ray_bound_sign_assumptions | 0.000e+00 | True |
| entire_bosonic_weakening_ray_has_negative_sigma_witness | 0.000e+00 | True |
| uniform_scalar_rescaling_stationarity_0.01 | 7.630e-18 | True |
| uniform_scalar_rescaling_exact_action_radial_0.01_1.0 | 3.469e-18 | True |
| uniform_scalar_rescaling_exact_action_radial_0.01_0.93 | 1.735e-18 | True |
| uniform_scalar_rescaling_Frechet_identity_0.01 | 4.338e-19 | True |
| uniform_scalar_rescaling_stationarity_0.03 | 2.063e-18 | True |
| uniform_scalar_rescaling_exact_action_radial_0.03_1.0 | 3.469e-18 | True |
| uniform_scalar_rescaling_exact_action_radial_0.03_0.93 | 1.041e-17 | True |
| uniform_scalar_rescaling_Frechet_identity_0.03 | 4.338e-18 | True |
| uniform_scalar_rescaling_stationarity_0.1 | 9.000e-17 | True |
| uniform_scalar_rescaling_exact_action_radial_0.1_1.0 | 1.388e-17 | True |
| uniform_scalar_rescaling_exact_action_radial_0.1_0.93 | 1.388e-17 | True |
| uniform_scalar_rescaling_Frechet_identity_0.1 | 2.946e-17 | True |
| uniform_scalar_rescaling_stationarity_0.3 | 6.290e-17 | True |
| uniform_scalar_rescaling_exact_action_radial_0.3_1.0 | 5.551e-17 | True |
| uniform_scalar_rescaling_exact_action_radial_0.3_0.93 | 0.000e+00 | True |
| uniform_scalar_rescaling_Frechet_identity_0.3 | 1.469e-16 | True |
| soft_bubble_negative_semidefinite_below_exp2_mu2 | 0.000e+00 | True |
