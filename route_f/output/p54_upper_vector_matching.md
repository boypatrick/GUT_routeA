# P54 upper heavy-vector finite matching subset

The existing upper PS saddle and generators are used. **This is not full (K,J,C) matching and not a physical fit.**

DR, MSbar, background-field Landau; fixed-VEV prescription. One-loop vertices and kinetic terms only.

Checks: 38/38. M_V^2 = 0.317035178968; mu = 0.563058770439; g = 0.562602889985.

Kpsi eigenvalue range: [-0.009019795749311221, -0.009019795749311221]. Kphi eigenvalue range: [-0.030065985831037402, -0.0008214493460420365].

Eliminated scalar masses are kept exactly in the vector-scalar bubble; retained scalar masses are expanded at leading hard order. The largest |m_retained^2|/M_V^2 is 0.328605; no precision error bound is claimed.

The -2, -1/2 and -3/2 finite terms in the vertex, scalar and fermion formulas are derived before taking d=4. The two spurion tensors satisfy the coset Ward identity and the complete gauge scale derivative; a scale check alone cannot fix finite constants.

## Finite linear coefficients in raw family units

| Spurion | Parent | Direct | Mirror |
|---|---|---:|---:|
| h_raw | H | 0.01771675106+0j | -7.838271578e-21+0.0001280086893j |
| h_raw | F | 0+0j | 0+0j |
| h_raw | L | 0+0j | 0+0j |
| h_raw | R | 0+0j | 0+0j |
| f_raw | H | 0+0j | 0+0j |
| f_raw | F | 0.08885260368-2.808647836e-36j | 2.066343517e-19-0.0003536675766j |
| f_raw | L | 0.06803155989+0j | 0+0j |
| f_raw | R | 0.06803155989+0j | 0+0j |

These are action-derived group contractions, not fitted observables or additional UV family parameters.

## Remaining requirements

- upper scalar-cubic kinetic, potential/tadpole and Wilson matching
- upper full source/contact hard-minus-EFT matching and operator mixing
- controlled retained-mass corrections for the stated accuracy
- lower covariant complete finite matching
- finite sequential seesaw matching
- fermionic light eigenstate feedback with matched families

## Regressions

| Check | Residual | Pass |
|---|---:|:---:|
| upper_coset_24_degenerate_masses | 1.598e-16 | True |
| coset_spinor_Casimir_3 | 1.480e-16 | True |
| all_upper_scalar_masses_positive | 0.000e+00 | True |
| active_heavy_tree_mass_block_zero | 0.000e+00 | True |
| Landau_Goldstone_Yukawa_h_zero | 0.000e+00 | True |
| Landau_Goldstone_Yukawa_f_zero | 0.000e+00 | True |
| active_linear_vector_mass_vertex_zero | 0.000e+00 | True |
| active_vector_current_no_eaten_scalar | 0.000e+00 | True |
| actual_parent_coset_Casimirs | 4.707e-16 | True |
| B0_quadrature_0.37_0.0 | 4.441e-16 | True |
| B0_quadrature_0.37_0.37 | 0.000e+00 | True |
| B0_quadrature_0.37_2.4 | 0.000e+00 | True |
| B0_quadrature_2.4_0.37 | 0.000e+00 | True |
| d_dim_vertex_finite_minus2 | 1.462e-05 | True |
| d_dim_scalar_finite_minus_half | 3.636e-06 | True |
| d_dim_fermion_finite_minus_three_halves | 5.643e-06 | True |
| Abelian_neutral_mass_Rxi_cancellation_0.0 | 0.000e+00 | True |
| Abelian_neutral_mass_Rxi_cancellation_0.3 | 0.000e+00 | True |
| Abelian_neutral_mass_Rxi_cancellation_1.0 | 0.000e+00 | True |
| Abelian_neutral_mass_Rxi_cancellation_2.0 | 0.000e+00 | True |
| Kpsi_Hermitian | 0.000e+00 | True |
| Kphi_real_symmetric | 0.000e+00 | True |
| positive_fermion_kinetic_metric | 0.000e+00 | True |
| positive_scalar_kinetic_metric | 0.000e+00 | True |
| PS_covariance_of_scalar_kinetic | 4.206e-17 | True |
| h_tree_intertwiner_dictionary | 0.000e+00 | True |
| h_coset_Yukawa_Ward | 2.093e-16 | True |
| h_six_invariant_closure | 1.388e-17 | True |
| h_full_gauge_matching_scale_derivative | 1.339e-14 | True |
| h_omitting_scalar_legs_is_detectable | 0.000e+00 | True |
| f_tree_intertwiner_dictionary | 0.000e+00 | True |
| f_coset_Yukawa_Ward | 4.358e-15 | True |
| f_six_invariant_closure | 3.894e-15 | True |
| f_full_gauge_matching_scale_derivative | 1.143e-13 | True |
| f_omitting_scalar_legs_is_detectable | 0.000e+00 | True |
| fit_rejects_this_partial_package | 0.000e+00 | True |
| fit_rejects_scheme_switch | 0.000e+00 | True |
| old_boolean_flags_do_not_bypass_KJC | 0.000e+00 | True |
