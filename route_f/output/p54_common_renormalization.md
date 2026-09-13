# Common bosonic renormalization and tadpole consistency

**Same-bare-parameter algebra is checked; full matching and physical control remain open.**

Checks: 46/46.

A finite fixed-VEV counterterm cannot be deleted at unchanged renormalized input. The equivalent MS input is p_MS=p_fixed+delta_p. Tree backgrounds, mass kernels and Wilson coefficients must be transported as well. The full first-order response is invariant; its division into tree and loop pieces is not.

| Parameter | Fixed-VEV input | Finite shift | Same-bare MS input |
|---|---:|---:|---:|
| mu2 | 1.40952534 | -0.0563307096 | 1.35319463 |
| nu2 | 0.381902025 | -5.46974409 | -5.08784206 |
| mus2 | 0.17192027 | -0.0342377144 | 0.137682556 |

## Broken radial bosonic block

Tree generalized squared masses: [0.013520470032803957, 0.09971569018567213, 2.8291658647815243]

Corrected bosonic hard+CT squared masses: [-0.7584479503161328, 0.09411893766041655, 3.0931051448460263]

Relative insertions: [-56.77588837144839, -0.05334352856629382, 0.09633903804792816]

These are fixed-gauge hard-potential curvatures, not pole masses or a full-model stability certificate. They exclude the fermionic and soft-IR completion.

## No Majorana-family rescue within the declared hard truncation

At fixed f_M, the finite fixed-VEV Majorana correction is only in the sigma,sigma entry: -4 sum_i M_i^4 log(M_i^2/mu^2)/(16 pi^2 sigma^2). The exact scalar inequality -x^2 log(x/mu^2) <= mu^4/(2e) bounds all three families without fitting or sampling them.

Their maximal positive curvature is 0.0875112373 omega^2. Even the matrix with this maximal correction has lowest generalized eigenvalue -0.71469528 omega^2. Thus no choice of the three Majorana masses rescues this frozen point in the declared one-loop hard-potential truncation. This is NOT an all-orders/full-model exclusion; a consistently reorganized scalar/IR treatment would be a different calculation.

The formal MS tree VEV response is [0.3066418847042487, -25.42887383994685, 1.6604721076627635]; the explicit loop response cancels it at first order. Tiny formal loop-parameter continuation tests this identity, not a new phenomenological vacuum scan.

## Still missing

- Full fermionic/common parameter feedback and controlled vacuum expansion.
- Complete Wilson/box graphs and lower gauge/ghost/Yukawa matching.
- Actual CHN/pre-existing-C5/mixed-sterile finite matching.

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| same_parameter_dictionary_at_both_endpoints | 0.000e+00 | True |
| canonical_radial_metric | 1.915e-16 | True |
| historical_full_bosonic_tadpole_reproduced | 1.110e-16 | True |
| fixed_VEV_radial_counterterm_cancels_once | 1.110e-16 | True |
| counterterm_radial_map_equals_full_field_contraction | 0.000e+00 | True |
| light_doublet_retuning_cannot_change_radial_block | 0.000e+00 | True |
| spectral_tadpole_finite_difference_2e-05 | 4.470e-09 | True |
| spectral_radial_Hessian_finite_difference_2e-05 | 4.543e-08 | True |
| spectral_tadpole_finite_difference_1e-05 | 1.117e-09 | True |
| spectral_radial_Hessian_finite_difference_1e-05 | 9.767e-08 | True |
| finite_scheme_radial_potential_identity_0 | 1.110e-16 | True |
| finite_scheme_radial_gradient_identity_0 | 5.554e-17 | True |
| finite_scheme_radial_potential_identity_1 | 0.000e+00 | True |
| finite_scheme_radial_gradient_identity_1 | 2.776e-16 | True |
| finite_scheme_radial_potential_identity_2 | 1.110e-16 | True |
| finite_scheme_radial_gradient_identity_2 | 2.221e-16 | True |
| broken_MS_tree_and_explicit_loop_VEV_shifts_cancel | 3.580e-15 | True |
| formal_stationary_root_1e-05 | 7.409e-16 | True |
| formal_stationary_root_5e-06 | 3.124e-16 | True |
| formal_stationary_root_2.5e-06 | 2.725e-17 | True |
| formal_stationary_conversion_remainder_is_second_order | 4.031e-03 | True |
| upper_common_valley_shift_equals_scheme_plus_loop | 1.735e-18 | True |
| upper_mass_scheme_transport_identity | 9.849e-17 | True |
| Wilson_parameter_scheme_transport_0 | 1.735e-18 | True |
| Wilson_scheme_tree_derivative_0 | 2.568e-10 | True |
| Wilson_parameter_scheme_transport_1 | 0.000e+00 | True |
| Wilson_scheme_tree_derivative_1 | 6.482e-10 | True |
| general_UV_and_EFT_finite_scheme_matching_chain_rule | 1.098e-16 | True |
| mass_cluster_decomposition_reconstructs_delta_nu2 | 1.624e-16 | True |
| negative_radial_witness_is_gauge_horizontal | 0.000e+00 | True |
| negative_radial_witness_is_not_PQ_phase | 0.000e+00 | True |
| one_massive_spinor_direction_per_family_at_broken_vacuum | 0.000e+00 | True |
| no_h_raw_mass_on_radial_background | 0.000e+00 | True |
| f_raw_radial_mass_only_depends_on_sigma | 0.000e+00 | True |
| Majorana_cap_saturated_by_three_equal_masses | 4.163e-17 | True |
| universal_radial_instability_witness_metric_normalization | 2.220e-16 | True |
| universal_radial_instability_witness_Rayleigh_identity | 1.110e-16 | True |
| negative_direction_survives_maximal_three_Majorana_curvature | 0.000e+00 | True |
| Majorana_fixed_VEV_curvature_direct_potential_0 | 2.821e-10 | True |
| Majorana_cap_test_0 | 0.000e+00 | True |
| Majorana_fixed_VEV_curvature_direct_potential_1 | 3.841e-08 | True |
| Majorana_cap_test_1 | 0.000e+00 | True |
| Majorana_fixed_VEV_curvature_direct_potential_2 | 3.529e-08 | True |
| Majorana_cap_test_2 | 0.000e+00 | True |
| Majorana_fixed_VEV_curvature_direct_potential_3 | 2.151e-08 | True |
| Majorana_cap_test_3 | 0.000e+00 | True |
