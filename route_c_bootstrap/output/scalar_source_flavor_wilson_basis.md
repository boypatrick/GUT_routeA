# Route-C P13-S Scalar/Source Flavor Tensors and Wilson Basis

P13-S chooses a symbolic flavor-tensor and chiral/Fierz Wilson basis for the P12-S BNV rows.  It does not choose numerical couplings, physical flavor rotations, RG factors, hadronic matrix elements, or proton lifetime limits.  Branch V remains open.

## Summary

| quantity | value |
| --- | ---: |
| p13_branch | scalar_source_flavor_tensor_and_wilson_basis |
| branch_v_still_required | True |
| bnv_or_sterile_bnv_rows_processed | 6 |
| canonical_wilson_basis_count | 4 |
| canonical_wilson_basis_ids | ["O_QQQL_S", "O_QQQL_T", "O_UDDN", "O_UUDE"] |
| direct_source_channel_rows | 4 |
| fierz_or_color_recoupling_required_rows | 2 |
| canonical_family_projection_required_rows | 1 |
| flavor_tensors_introduced | 11 |
| independent_complex_vertex_components | 78 |
| physical_proton_bounds_evaluable_now | 0 |

## Canonical Chiral Wilson Basis

| basis | operator | independent components | symmetry |
| --- | --- | ---: | --- |
| `O_QQQL_S` | `epsilon_abc epsilon_ij epsilon_kl (Q_p^{a i} Q_r^{b j})(Q_s^{c k} L_t^l)` | 54 | C_{prst}=C_{rpst} for the first two Q-family indices |
| `O_QQQL_T` | `epsilon_abc (tau^A epsilon)_ij (tau^A epsilon)_kl (Q_p^{a i} Q_r^{b j})(Q_s^{c k} L_t^l)` | 27 | C_{prst}=-C_{rpst} for the first two Q-family indices |
| `O_UDDN` | `epsilon_abc (u^c_p^a d^c_r^b)(d^c_s^c nu^c_t)` | 27 | C_{prst}=-C_{psrt} for the two d^c-family indices after canonical projection |
| `O_UUDE` | `epsilon_abc (u^c_p^a u^c_r^b)(d^c_s^c e^c_t)` | 27 | C_{prst}=-C_{rpst} for the two u^c-family indices |

## Flavor Tensors

| tensor | pair | source | symmetry | independent components |
| --- | --- | --- | --- | ---: |
| `Lambda_QQ__S_3_1_Ym1o3_BLm2o3` | `Q Q` | `S_3_1_Ym1o3_BLm2o3` | symmetric_family_pair | 6 |
| `Lambda_QL__S_3_1_Ym1o3_BLm2o3` | `Q L` | `S_3_1_Ym1o3_BLm2o3` | general_3x3_family_matrix | 9 |
| `Lambda_QQ__S_3_3_Ym1o3_BLm2o3` | `Q Q` | `S_3_3_Ym1o3_BLm2o3` | antisymmetric_family_pair | 3 |
| `Lambda_QL__S_3_3_Ym1o3_BLm2o3` | `Q L` | `S_3_3_Ym1o3_BLm2o3` | general_3x3_family_matrix | 9 |
| `Lambda_ucuc__S_b3_1_Y4o3_BL2o3` | `u^c u^c` | `S_b3_1_Y4o3_BL2o3` | antisymmetric_family_pair | 3 |
| `Lambda_dcec__S_b3_1_Y4o3_BL2o3` | `d^c e^c` | `S_b3_1_Y4o3_BL2o3` | general_3x3_family_matrix | 9 |
| `Lambda_ucdc__S_b3_1_Y1o3_BL2o3` | `u^c d^c` | `S_b3_1_Y1o3_BL2o3` | general_3x3_family_matrix | 9 |
| `Lambda_ucec__S_b3_1_Y1o3_BL2o3` | `u^c e^c` | `S_b3_1_Y1o3_BL2o3` | general_3x3_family_matrix | 9 |
| `Lambda_dcnuc__S_b3_1_Y1o3_BL2o3` | `d^c nu^c` | `S_b3_1_Y1o3_BL2o3` | general_3x3_family_matrix | 9 |
| `Lambda_ucnuc__S_3_1_Y2o3_BLm2o3` | `u^c nu^c` | `S_3_1_Y2o3_BLm2o3` | general_3x3_family_matrix | 9 |
| `Lambda_dcdc__S_3_1_Y2o3_BLm2o3` | `d^c d^c` | `S_3_1_Y2o3_BLm2o3` | antisymmetric_family_pair | 3 |

## BNV Wilson Matching Rows

| operator | basis | source | recoupling | projection | coefficient |
| --- | --- | --- | --- | --- | --- |
| `QQQL` | `O_QQQL_S` | `S_3_1_Ym1o3_BLm2o3` | direct_source_channel | encoded_by_vertex_pair_symmetry | `C[O_QQQL_S]_prst += Lambda_QQ__S_3_1_Ym1o3_BLm2o3_ij Lambda_QL__S_3_1_Ym1o3_BLm2o3_kl^*/M_S_3_1_Ym1o3_BLm2o3^2` |
| `QQQL` | `O_QQQL_T` | `S_3_3_Ym1o3_BLm2o3` | direct_source_channel | encoded_by_vertex_pair_symmetry | `C[O_QQQL_T]_prst += Lambda_QQ__S_3_3_Ym1o3_BLm2o3_ij Lambda_QL__S_3_3_Ym1o3_BLm2o3_kl^*/M_S_3_3_Ym1o3_BLm2o3^2` |
| `u^c u^c d^c e^c` | `O_UUDE` | `S_b3_1_Y4o3_BL2o3` | direct_source_channel | encoded_by_vertex_pair_symmetry | `C[O_UUDE]_prst += Lambda_ucuc__S_b3_1_Y4o3_BL2o3_ij Lambda_dcec__S_b3_1_Y4o3_BL2o3_kl^*/M_S_b3_1_Y4o3_BL2o3^2` |
| `u^c u^c d^c e^c` | `O_UUDE` | `S_b3_1_Y1o3_BL2o3` | requires_chiral_fierz_color_recoupling | canonical_projection_included_in_recoupling_matrix | `C[O_UUDE]_prst += rho_O_UUDE__S_b3_1_Y1o3_BL2o3 Lambda_ucdc__S_b3_1_Y1o3_BL2o3_ij Lambda_ucec__S_b3_1_Y1o3_BL2o3_kl^*/M_S_b3_1_Y1o3_BL2o3^2` |
| `u^c d^c d^c nu^c` | `O_UDDN` | `S_b3_1_Y1o3_BL2o3` | direct_source_channel | requires_d_pair_antisymmetrization_projection | `C[O_UDDN]_prst += Pi_d_asym_rs[Lambda_ucdc__S_b3_1_Y1o3_BL2o3_pr Lambda_dcnuc__S_b3_1_Y1o3_BL2o3_st^*/M_S_b3_1_Y1o3_BL2o3^2]` |
| `u^c d^c d^c nu^c` | `O_UDDN` | `S_3_1_Y2o3_BLm2o3` | requires_chiral_fierz_color_recoupling | canonical_projection_included_in_recoupling_matrix | `C[O_UDDN]_prst += rho_O_UDDN__S_3_1_Y2o3_BLm2o3 Lambda_ucnuc__S_3_1_Y2o3_BLm2o3_ij Lambda_dcdc__S_3_1_Y2o3_BLm2o3_kl^*/M_S_3_1_Y2o3_BLm2o3^2` |

## Verification

| check | value |
| --- | ---: |
| all_p12_bnv_rows_mapped | True |
| every_row_has_basis_id | True |
| every_row_has_left_and_right_flavor_tensor | True |
| branch_v_retained | True |
| physical_proton_bounds_evaluable_now | 0 |
| p13_scalar_source_basis_layer_passes | True |

## Next Stage

`P14_scalar_source_recoupling_matrix_or_branch_v`

Either compute the chiral/color Fierz recoupling matrices for the two non-direct scalar/source rows and then add flavor rotations, or return to the completed current-current vector Branch V before proton limits.
