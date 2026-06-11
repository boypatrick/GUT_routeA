# Route-C P16 Physical Hermitian Vector Multiplet Assembly

P16 assembles charge-conjugate root maps into symbolic physical Hermitian vector multiplets and audits whether canonical vector BNV skeletons can appear.  It does not fix numerical Clebsch phases, flavor rotations, RG factors, hadronic matrix elements, or proton lifetime limits.

## Physical Vector Rule

- root pair: `E_q + E_-q with symbolic cross-root phase xi`
- product: `J_q^mu J_-q,mu`
- prefactor: `-2 g_X^2 xi c_q c_-q / M_X^2`

## Summary

| quantity                                          | value                                                                                                                                     |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| p15_same_map_rows_were_canonical_bnv              | 0                                                                                                                                         |
| physical_cross_root_rows                          | 448                                                                                                                                       |
| physical_cross_root_rows_by_family                | `{"Pati-Salam SU(4)_C leptoquark": 144, "Spin(10)/Pati-Salam off-face generator": 288, "broken SU(2)_R charged current": 16}`             |
| physical_pair_status_counts                       | `{"B_and_L_conserving_physical_pair": 160, "BminusL_preserving_BNV_noncanonical_pair": 144, "canonical_vector_BNV_basis_candidate": 144}` |
| canonical_vector_bnv_basis_counts                 | `{"O_V_QQ_DbarNbar": 72, "O_V_QQ_UbarEbar": 72, "none": 304}`                                                                             |
| bnv_candidate_rows_before_flavor                  | 288                                                                                                                                       |
| canonical_vector_bnv_candidate_rows_before_flavor | 144                                                                                                                                       |
| physical_proton_bounds_evaluable_now              | 0                                                                                                                                         |
| all_rows_have_cross_root_phase                    | True                                                                                                                                      |
| all_rows_have_charge_audit                        | True                                                                                                                                      |
| branch_v_physical_multiplet_gate_passes           | True                                                                                                                                      |

## Sample Canonical Vector BNV Candidates

| basis             | mass        | fields                   | phase                                        |
| ----------------- | ----------- | ------------------------ | -------------------------------------------- |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_rg_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_rg_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_rb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_rb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_gb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_u_minus__6_2_2_C_gb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_rg_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_rg_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_rb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_rb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_gb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_rell_d_minus__6_2_2_C_gb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_rg_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_rg_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_rb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_rb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_gb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_u_minus__6_2_2_C_gb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_rg_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_rg_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_rb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_rb_d_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_gb_u_plus` |
| `O_V_QQ_DbarNbar` | `M_(6,2,2)` | `Q Q bar(nu^c) bar(d^c)` | `xi_6_2_2_C_gell_d_minus__6_2_2_C_gb_d_plus` |

## Interpretation

- From P15: P15 same-map J J^dagger rows were all B/L conserving.  P16 allows charge-conjugate root-map products J_q J_-q, which can produce vector-mediated BNV skeletons.
- Proton limits: The gate still has symbolic cross-root phases, no physical flavor rotations, no RG evolution, and no hadronic/channel inputs.
- Branch S: The scalar/source Branch S remains available and is not discarded by the vector result.

## Next Stage

`P17_vector_clebsch_projection_and_flavor_interface`

Fix a physical Hermitian generator normalization and cross-root Clebsch phase convention for the canonical vector BNV rows, then export a flavor-rotation interface.  Do not insert proton limits until RG and hadronic inputs are supplied.
