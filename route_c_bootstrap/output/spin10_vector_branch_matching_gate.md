# Route-C P14-V Spin(10) Vector-Branch Matching Gate

P14-V returns to Branch V and builds a current-current vector matching gate before proton limits.  It does not reuse the P7 scalar/source chiral-pair poles, and it does not compute physical proton lifetimes.

## Matching Formula

- current: `J_X^mu = sum_a c_a psi_target(a)^dagger bar_sigma^mu psi_source(a)`
- tree-level gate: `L_eff[X] = - g_X^2 J_X^mu J_X,mu^dagger / M_X^2`

## Summary

| quantity                                            | value                                                                                                                                                                                                                                                                        |
| --------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| broken_vector_maps                                  | 32                                                                                                                                                                                                                                                                           |
| broken_vector_maps_by_family                        | `{"Pati-Salam SU(4)_C leptoquark": 6, "Spin(10)/Pati-Salam off-face generator": 24, "broken SU(2)_R charged current": 2}`                                                                                                                                                    |
| current_transition_entries                          | 80                                                                                                                                                                                                                                                                           |
| current_transition_entries_by_family                | `{"Pati-Salam SU(4)_C leptoquark": 24, "Spin(10)/Pati-Salam off-face generator": 48, "broken SU(2)_R charged current": 8}`                                                                                                                                                   |
| current_pair_rows                                   | 72                                                                                                                                                                                                                                                                           |
| current_pair_rows_by_family                         | `{"Pati-Salam SU(4)_C leptoquark": 36, "Spin(10)/Pati-Salam off-face generator": 24, "broken SU(2)_R charged current": 12}`                                                                                                                                                  |
| matching_gate_counts                                | `{"off_face_diquark_like_pair_crossing_required": 12, "off_face_mixed_pair_crossing_required": 12, "ps_leptoquark_mixed_pair_crossing_required": 24, "ps_leptoquark_same-sector_pair_noncanonical_until_fierz": 12, "right_current_B_and_L_conserving_not_proton_seed": 12}` |
| p10_mass_matrix_layer_passes                        | True                                                                                                                                                                                                                                                                         |
| p11_p7_rows_remaining_scalar_or_source_mass_debt    | 18                                                                                                                                                                                                                                                                           |
| p11_p7_rows_with_candidate_adjoint_mass_blocks      | 0                                                                                                                                                                                                                                                                            |
| vector_branch_does_not_reuse_p7_scalar_source_poles | True                                                                                                                                                                                                                                                                         |
| physical_proton_bounds_evaluable_now                | 0                                                                                                                                                                                                                                                                            |
| branch_v_matching_gate_passes                       | True                                                                                                                                                                                                                                                                         |

## Gate Counts

| gate                                                      | count |
| --------------------------------------------------------- | ----- |
| `off_face_diquark_like_pair_crossing_required`            | 12    |
| `off_face_mixed_pair_crossing_required`                   | 12    |
| `ps_leptoquark_mixed_pair_crossing_required`              | 24    |
| `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` | 12    |
| `right_current_B_and_L_conserving_not_proton_seed`        | 12    |

## Sample Current-Pair Rows

| family                        | mass   | transition pair        | gate                                                      |
| ----------------------------- | ------ | ---------------------- | --------------------------------------------------------- |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; Q->L`          | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; nu^c->u^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; e^c->d^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; nu^c->u^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; e^c->d^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `nu^c->u^c ; e^c->d^c` | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `L->Q ; L->Q`          | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `L->Q ; u^c->nu^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `L->Q ; d^c->e^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `L->Q ; u^c->nu^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `L->Q ; d^c->e^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `u^c->nu^c ; d^c->e^c` | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; Q->L`          | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; nu^c->u^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; e^c->d^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; nu^c->u^c`     | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `Q->L ; e^c->d^c`      | `ps_leptoquark_mixed_pair_crossing_required`              |
| Pati-Salam SU(4)_C leptoquark | `M_LQ` | `nu^c->u^c ; e^c->d^c` | `ps_leptoquark_same-sector_pair_noncanonical_until_fierz` |

## Proton-Limit Blockers

- current-current to all-left chiral Wilson basis recoupling
- physical flavor rotations
- numerical vector masses and couplings
- RG running from matching scale to hadronic scale
- lattice or chiral hadronic matrix elements
- experimental channel selection and lifetime limit

## Parallel Branch Status

- scalar/source branch: P13-S retained as a separate scalar/source Wilson-basis branch
- vector branch: Branch V current-current matching gate is now explicit

## Next Stage

`P15_V_current_to_all_left_fierz_and_flavor_gate`

Choose the convention-fixed current-current to all-left operator map for the vector branch, or keep P13-S scalar/source recoupling as the active proton-matching path.  Do not insert proton limits until one branch has flavor, RG, and hadronic inputs.
