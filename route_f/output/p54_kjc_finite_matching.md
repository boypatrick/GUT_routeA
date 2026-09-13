# P54 one-loop (K,J,C) transport checkpoint

**Complete algebra, partial diagram content. No physical fit.**

108/108 checks pass. The actual stationary upper PS chart has 304 physical real scalars; 248 active parents plus PQ remain after eliminating 55 upper scalars. The final 129-dimensional chart retains whole H/F parents and PQ. Further elimination is an off-shell algebra regression, not a claim that L/R can physically be decoupled at the upper scale.

The current slots retain seven P-SPEC1 pairs and add two action-defined different-colour diquark pairs, with the unknown two UV families factored out. The earlier same-colour pairs do not see the upper triplets: antisymmetric colour contractions vanish. The new probes detect their finite induced contact response. One-loop contributions are the computed active vector vertex, universal vector fermion external legs (including heavy-source endpoints), and active vector-scalar kinetic term. Missing self-energies, heavy-source 1PI vertices, boxes, operator mixing and lower matching remain unknown.

The initial delta C is zero **for this selected additive set of graphs only**; the induced delta C is nonzero. This is not a zero boundary condition for the complete finite contact matching.

| z=p_E^2 | Induced upper delta C norm | Drop final delta C response error | Double-count scalar leg error |
|---|---:|---:|---:|
| (0.013+0j) | 0.2483693 | 0.9998932 | 0.2710652 |
| (0.27+0j) | 0.1336376 | 0.961138 | 0.3845887 |
| (-0.11+0.07j) | 0.4106111 | 0.6078394 | 0.355261 |
| (0.03+0.12j) | 0.2166991 | 0.9651812 | 0.3883261 |

Errors use ||A-B||/max(1,||A||,||B||); they are not physical cross-section errors.

The Frechet derivative agrees with direct finite differences and both elimination orders. Canonicalizing the field metric also transforms the mass matrix and sources; normalizing J while retaining the unnormalized K double counts external scalar legs.

## Fit gate

physical fit map incomplete: upper_finite_gauge, lower_finite_gauge, upper_finite_Yukawa, lower_finite_Yukawa, all_active_PS_flow, lower_EFT_flow, sequential_Weinberg_matching, same_action_scalar_feedback, upper_full_KJC, lower_covariant_full_KJC, all_six_PS_Yukawas_run, sequential_finite_seesaw, fermion_full_light_feedback, retained_mass_error_control

## Regression ledger

| Check | Residual | Pass |
|---|---:|:---:|
| 304_physical_upper_chart_orthonormal | 2.191e-16 | True |
| chart_excludes_24_eaten_upper_directions | 0.000e+00 | True |
| PQ_source_zero | 0.000e+00 | True |
| PQ_loop_source_zero_in_computed_subset | 0.000e+00 | True |
| different_colour_probes_detect_upper_scalar_exchange | 0.000e+00 | True |
| z=0.013_context_checked_sum | 0.000e+00 | True |
| z=0.013_upper_induced_delta_C_nonzero | 0.000e+00 | True |
| z=0.013_order0_K_staged | 0.000e+00 | True |
| z=0.013_order0_K_reverse | 0.000e+00 | True |
| z=0.013_order0_J_staged | 0.000e+00 | True |
| z=0.013_order0_J_reverse | 0.000e+00 | True |
| z=0.013_order0_C_staged | 5.733e-17 | True |
| z=0.013_order0_C_reverse | 5.733e-17 | True |
| z=0.013_order0_response | 1.050e-20 | True |
| z=0.013_order1_K_staged | 0.000e+00 | True |
| z=0.013_order1_K_reverse | 0.000e+00 | True |
| z=0.013_order1_J_staged | 0.000e+00 | True |
| z=0.013_order1_J_reverse | 0.000e+00 | True |
| z=0.013_order1_C_staged | 2.216e-17 | True |
| z=0.013_order1_C_reverse | 2.216e-17 | True |
| z=0.013_order1_response | 2.152e-21 | True |
| z=0.013_K_finite_difference | 3.212e-13 | True |
| z=0.013_J_finite_difference | 1.916e-12 | True |
| z=0.013_C_finite_difference | 1.694e-10 | True |
| z=0.013_canonical_response | 7.351e-17 | True |
| z=0.013_canonical_mass_insertion | 2.971e-18 | True |
| z=0.013_double_scalar_leg_detected | 0.000e+00 | True |
| z=0.013_dropping_induced_loop_C_detected | 0.000e+00 | True |
| z=0.27_context_checked_sum | 0.000e+00 | True |
| z=0.27_upper_induced_delta_C_nonzero | 0.000e+00 | True |
| z=0.27_order0_K_staged | 0.000e+00 | True |
| z=0.27_order0_K_reverse | 0.000e+00 | True |
| z=0.27_order0_J_staged | 0.000e+00 | True |
| z=0.27_order0_J_reverse | 0.000e+00 | True |
| z=0.27_order0_C_staged | 1.909e-18 | True |
| z=0.27_order0_C_reverse | 2.898e-18 | True |
| z=0.27_order0_response | 8.875e-19 | True |
| z=0.27_order1_K_staged | 0.000e+00 | True |
| z=0.27_order1_K_reverse | 0.000e+00 | True |
| z=0.27_order1_J_staged | 0.000e+00 | True |
| z=0.27_order1_J_reverse | 0.000e+00 | True |
| z=0.27_order1_C_staged | 4.246e-19 | True |
| z=0.27_order1_C_reverse | 4.540e-19 | True |
| z=0.27_order1_response | 3.874e-19 | True |
| z=0.27_K_finite_difference | 1.350e-12 | True |
| z=0.27_J_finite_difference | 1.916e-12 | True |
| z=0.27_C_finite_difference | 4.268e-11 | True |
| z=0.27_canonical_response | 7.504e-17 | True |
| z=0.27_canonical_mass_insertion | 1.155e-17 | True |
| z=0.27_double_scalar_leg_detected | 0.000e+00 | True |
| z=0.27_dropping_induced_loop_C_detected | 0.000e+00 | True |
| z=(-0.11+0.07j)_context_checked_sum | 0.000e+00 | True |
| z=(-0.11+0.07j)_upper_induced_delta_C_nonzero | 0.000e+00 | True |
| z=(-0.11+0.07j)_order0_K_staged | 0.000e+00 | True |
| z=(-0.11+0.07j)_order0_K_reverse | 0.000e+00 | True |
| z=(-0.11+0.07j)_order0_J_staged | 0.000e+00 | True |
| z=(-0.11+0.07j)_order0_J_reverse | 0.000e+00 | True |
| z=(-0.11+0.07j)_order0_C_staged | 2.416e-17 | True |
| z=(-0.11+0.07j)_order0_C_reverse | 2.421e-17 | True |
| z=(-0.11+0.07j)_order0_response | 0.000e+00 | True |
| z=(-0.11+0.07j)_order1_K_staged | 0.000e+00 | True |
| z=(-0.11+0.07j)_order1_K_reverse | 0.000e+00 | True |
| z=(-0.11+0.07j)_order1_J_staged | 0.000e+00 | True |
| z=(-0.11+0.07j)_order1_J_reverse | 0.000e+00 | True |
| z=(-0.11+0.07j)_order1_C_staged | 1.297e-18 | True |
| z=(-0.11+0.07j)_order1_C_reverse | 1.327e-17 | True |
| z=(-0.11+0.07j)_order1_response | 0.000e+00 | True |
| z=(-0.11+0.07j)_K_finite_difference | 3.066e-13 | True |
| z=(-0.11+0.07j)_J_finite_difference | 1.916e-12 | True |
| z=(-0.11+0.07j)_C_finite_difference | 3.706e-11 | True |
| z=(-0.11+0.07j)_canonical_response | 1.707e-16 | True |
| z=(-0.11+0.07j)_canonical_mass_insertion | 4.093e-18 | True |
| z=(-0.11+0.07j)_double_scalar_leg_detected | 0.000e+00 | True |
| z=(-0.11+0.07j)_dropping_induced_loop_C_detected | 0.000e+00 | True |
| z=(0.03+0.12j)_context_checked_sum | 0.000e+00 | True |
| z=(0.03+0.12j)_upper_induced_delta_C_nonzero | 0.000e+00 | True |
| z=(0.03+0.12j)_order0_K_staged | 0.000e+00 | True |
| z=(0.03+0.12j)_order0_K_reverse | 0.000e+00 | True |
| z=(0.03+0.12j)_order0_J_staged | 0.000e+00 | True |
| z=(0.03+0.12j)_order0_J_reverse | 0.000e+00 | True |
| z=(0.03+0.12j)_order0_C_staged | 4.544e-19 | True |
| z=(0.03+0.12j)_order0_C_reverse | 3.508e-17 | True |
| z=(0.03+0.12j)_order0_response | 0.000e+00 | True |
| z=(0.03+0.12j)_order1_K_staged | 0.000e+00 | True |
| z=(0.03+0.12j)_order1_K_reverse | 0.000e+00 | True |
| z=(0.03+0.12j)_order1_J_staged | 0.000e+00 | True |
| z=(0.03+0.12j)_order1_J_reverse | 0.000e+00 | True |
| z=(0.03+0.12j)_order1_C_staged | 2.779e-17 | True |
| z=(0.03+0.12j)_order1_C_reverse | 2.269e-17 | True |
| z=(0.03+0.12j)_order1_response | 0.000e+00 | True |
| z=(0.03+0.12j)_K_finite_difference | 4.822e-13 | True |
| z=(0.03+0.12j)_J_finite_difference | 1.916e-12 | True |
| z=(0.03+0.12j)_C_finite_difference | 6.560e-11 | True |
| z=(0.03+0.12j)_canonical_response | 8.782e-17 | True |
| z=(0.03+0.12j)_canonical_mass_insertion | 3.058e-18 | True |
| z=(0.03+0.12j)_double_scalar_leg_detected | 0.000e+00 | True |
| z=(0.03+0.12j)_dropping_induced_loop_C_detected | 0.000e+00 | True |
| reject_merge_different_background | 0.000e+00 | True |
| reject_merge_different_scale | 0.000e+00 | True |
| reject_merge_different_scheme | 0.000e+00 | True |
| reject_merge_different_current_basis | 0.000e+00 | True |
| reject_double_counted_diagram | 0.000e+00 | True |
| synthetic_dense_deltaB_deltaD_K | 1.983e-11 | True |
| synthetic_dense_deltaB_deltaD_J | 1.303e-11 | True |
| synthetic_dense_deltaB_deltaD_C | 3.733e-11 | True |
| synthetic_dense_response | 3.126e-16 | True |
| singular_elimination_rejected | 0.000e+00 | True |
| partial_KJC_cannot_enable_fit | 0.000e+00 | True |
