# Route-C P12-S Scalar/Source Chiral-Pair Branch

P12-S assigns scalar/source mass denominators to the P7 chiral-pair poles.  It does not compute physical proton lifetimes, and it does not close the separate current-current vector branch.

## Branch Decision

- branches to try: `completed_current_current_vector_branch, scalar_source_chiral_pair_branch`
- activated first: `scalar_source_chiral_pair_branch`
- reason: P11 showed that all 18 P7 chiral-pair pole rows remain scalar/source mass debt rather than adjoint-vector mass matches.

## Ansatz

```text
L ⊃ lambda_{ab,R} psi_a psi_b S_R + h.c. + M_{S_R}^2 |S_R|^2
C6(12|34;R)=lambda_{12,R} lambda^*_{34,R}/M_{S_R}^2
M_{S_R}^2 > 0 for every activated source sector
```

## Summary

| quantity | value |
| --- | ---: |
| p12_branch_choice | scalar_source_chiral_pair_branch |
| other_required_future_branch | completed_current_current_vector_branch |
| p11_p7_rows_remaining_scalar_source_debt | 18 |
| p7_rows_converted_to_scalar_source_rows | 18 |
| unique_scalar_source_sectors | 12 |
| operator_status_counts | {"BNV_with_sterile_neutrino": 2, "B_and_L_conserving": 12, "standard_BNV_seed": 4} |
| standard_bnv_seed_rows | 4 |
| sterile_bnv_rows | 2 |
| baryon_violating_or_sterile_bnv_rows | 6 |
| physical_proton_bounds_evaluable_now | 0 |
| all_scalar_source_masses_symbolic | True |

## Scalar/Source Sectors

| source | mass | rows | operators |
| --- | --- | ---: | --- |
| `S_1_1_Y1_BL2` | `M_S_1_1_Y1_BL2` | 1 | `LL e^c nu^c` |
| `S_1_2_Y1o2_BL0` | `M_S_1_2_Y1o2_BL0` | 3 | `LL e^c nu^c`, `QL u^c e^c`, `QQ u^c d^c` |
| `S_1_2_Ym1o2_BL0` | `M_S_1_2_Ym1o2_BL0` | 1 | `QL d^c nu^c` |
| `S_3_1_Y2o3_BLm2o3` | `M_S_3_1_Y2o3_BLm2o3` | 1 | `u^c d^c d^c nu^c` |
| `S_3_1_Ym1o3_BLm2o3` | `M_S_3_1_Ym1o3_BLm2o3` | 2 | `QQ u^c d^c`, `QQQL` |
| `S_3_3_Ym1o3_BLm2o3` | `M_S_3_3_Ym1o3_BLm2o3` | 1 | `QQQL` |
| `S_8_2_Y1o2_BL0` | `M_S_8_2_Y1o2_BL0` | 1 | `QQ u^c d^c` |
| `S_b3_1_Y1o3_BL2o3` | `M_S_b3_1_Y1o3_BL2o3` | 4 | `QL d^c nu^c`, `QL u^c e^c`, `u^c d^c d^c nu^c`, `u^c u^c d^c e^c` |
| `S_b3_1_Y4o3_BL2o3` | `M_S_b3_1_Y4o3_BL2o3` | 1 | `u^c u^c d^c e^c` |
| `S_b3_2_Ym1o6_BLm4o3` | `M_S_b3_2_Ym1o6_BLm4o3` | 1 | `QL d^c nu^c` |
| `S_b3_2_Ym7o6_BLm4o3` | `M_S_b3_2_Ym7o6_BLm4o3` | 1 | `QL u^c e^c` |
| `S_b6_1_Ym1o3_BLm2o3` | `M_S_b6_1_Ym1o3_BLm2o3` | 1 | `QQ u^c d^c` |

## BNV Scalar/Source Rows

| operator | source | coefficient |
| --- | --- | --- |
| `QQQL` | `S_3_1_Ym1o3_BLm2o3` | `C6_S[QQQL; X[(3,1);Y=-1/3;B-L=-2/3]] = lambda_QQ_3_1_Ym1o3_BLm2o3 lambda_QL_3_1_Ym1o3_BLm2o3^*/M_S_3_1_Ym1o3_BLm2o3^2` |
| `QQQL` | `S_3_3_Ym1o3_BLm2o3` | `C6_S[QQQL; X[(3,3);Y=-1/3;B-L=-2/3]] = lambda_QQ_3_3_Ym1o3_BLm2o3 lambda_QL_3_3_Ym1o3_BLm2o3^*/M_S_3_3_Ym1o3_BLm2o3^2` |
| `u^c u^c d^c e^c` | `S_b3_1_Y4o3_BL2o3` | `C6_S[u^c u^c d^c e^c; X[(bar3,1);Y=4/3;B-L=2/3]] = lambda_ucuc_b3_1_Y4o3_BL2o3 lambda_dcec_b3_1_Y4o3_BL2o3^*/M_S_b3_1_Y4o3_BL2o3^2` |
| `u^c u^c d^c e^c` | `S_b3_1_Y1o3_BL2o3` | `C6_S[u^c u^c d^c e^c; X[(bar3,1);Y=1/3;B-L=2/3]] = lambda_ucdc_b3_1_Y1o3_BL2o3 lambda_ucec_b3_1_Y1o3_BL2o3^*/M_S_b3_1_Y1o3_BL2o3^2` |
| `u^c d^c d^c nu^c` | `S_b3_1_Y1o3_BL2o3` | `C6_S[u^c d^c d^c nu^c; X[(bar3,1);Y=1/3;B-L=2/3]] = lambda_ucdc_b3_1_Y1o3_BL2o3 lambda_dcnuc_b3_1_Y1o3_BL2o3^*/M_S_b3_1_Y1o3_BL2o3^2` |
| `u^c d^c d^c nu^c` | `S_3_1_Y2o3_BLm2o3` | `C6_S[u^c d^c d^c nu^c; X[(3,1);Y=2/3;B-L=-2/3]] = lambda_ucnuc_3_1_Y2o3_BLm2o3 lambda_dcdc_3_1_Y2o3_BLm2o3^*/M_S_3_1_Y2o3_BLm2o3^2` |

## Verification

| check | value |
| --- | ---: |
| all_p11_scalar_source_debt_rows_accounted_for | True |
| every_scalar_row_has_mass_symbol | True |
| every_source_sector_has_positive_mass_condition | True |
| physical_proton_bounds_evaluable_now | 0 |
| p12_scalar_source_layer_passes | True |

## Next Stage

`P13_scalar_source_flavor_tensor_and_operator_basis`

Choose scalar/source flavor tensors and a chiral/Fierz Wilson basis for the BNV rows, or return to the parallel current-vector branch before inserting proton limits.
