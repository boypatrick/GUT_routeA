# Route-C P15-V Vector Current to All-Left Fierz Gate

P15-V fixes a two-component Fierz/crossing convention for Branch V and audits the crossed all-left skeletons.  It does not compute physical proton lifetimes, and it does not override the scalar/source Branch S.

## Fixed Convention

- input: `(psi_t1^dagger bar_sigma^mu psi_s1)(psi_s2^dagger bar_sigma_mu psi_t2)`
- Fierz: `(chi^dagger bar_sigma^mu psi)(eta^dagger bar_sigma_mu xi) = 2 (chi^dagger eta^dagger)(psi xi)`
- crossed skeleton: `(psi_s1 psi_t2)(bar(psi_t1) bar(psi_s2))`
- coefficient prefactor: `-2 g_X^2 c1 c2^*/M_X^2`

## Summary

| quantity                             | value                                                                                                                                                                                                                                                                        |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| p14_current_pair_rows_consumed       | 72                                                                                                                                                                                                                                                                           |
| p15_rows_generated                   | 72                                                                                                                                                                                                                                                                           |
| all_left_status_counts               | `{"B_and_L_conserving_crossed_skeleton": 72}`                                                                                                                                                                                                                                |
| canonical_bnv_basis_counts           | `{"none": 72}`                                                                                                                                                                                                                                                               |
| p14_gate_counts_replayed             | `{"off_face_diquark_like_pair_crossing_required": 12, "off_face_mixed_pair_crossing_required": 12, "ps_leptoquark_mixed_pair_crossing_required": 24, "ps_leptoquark_same-sector_pair_noncanonical_until_fierz": 12, "right_current_B_and_L_conserving_not_proton_seed": 12}` |
| canonical_bnv_rows_evaluable_now     | 0                                                                                                                                                                                                                                                                            |
| physical_proton_bounds_evaluable_now | 0                                                                                                                                                                                                                                                                            |
| all_rows_have_crossed_fields         | True                                                                                                                                                                                                                                                                         |
| all_rows_have_charge_audit           | True                                                                                                                                                                                                                                                                         |
| branch_v_all_left_gate_passes        | True                                                                                                                                                                                                                                                                         |

## Sample Crossed Rows

| row | mass   | fields                       | status                                |
| --- | ------ | ---------------------------- | ------------------------------------- |
| 0   | `M_LQ` | `Q L bar(L) bar(Q)`          | `B_and_L_conserving_crossed_skeleton` |
| 1   | `M_LQ` | `Q u^c bar(L) bar(nu^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 2   | `M_LQ` | `Q d^c bar(L) bar(e^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 3   | `M_LQ` | `Q u^c bar(L) bar(nu^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 4   | `M_LQ` | `Q d^c bar(L) bar(e^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 5   | `M_LQ` | `nu^c d^c bar(u^c) bar(e^c)` | `B_and_L_conserving_crossed_skeleton` |
| 6   | `M_LQ` | `L Q bar(Q) bar(L)`          | `B_and_L_conserving_crossed_skeleton` |
| 7   | `M_LQ` | `L nu^c bar(Q) bar(u^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 8   | `M_LQ` | `L e^c bar(Q) bar(d^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 9   | `M_LQ` | `L nu^c bar(Q) bar(u^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 10  | `M_LQ` | `L e^c bar(Q) bar(d^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 11  | `M_LQ` | `u^c e^c bar(nu^c) bar(d^c)` | `B_and_L_conserving_crossed_skeleton` |
| 12  | `M_LQ` | `Q L bar(L) bar(Q)`          | `B_and_L_conserving_crossed_skeleton` |
| 13  | `M_LQ` | `Q u^c bar(L) bar(nu^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 14  | `M_LQ` | `Q d^c bar(L) bar(e^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 15  | `M_LQ` | `Q u^c bar(L) bar(nu^c)`     | `B_and_L_conserving_crossed_skeleton` |
| 16  | `M_LQ` | `Q d^c bar(L) bar(e^c)`      | `B_and_L_conserving_crossed_skeleton` |
| 17  | `M_LQ` | `nu^c d^c bar(u^c) bar(e^c)` | `B_and_L_conserving_crossed_skeleton` |

## Interpretation

- Branch V: Branch V is now in an all-left crossed-skeleton bookkeeping basis.  Under the current sparse-map pairing, no row is yet a canonical BNV proton operator.
- Branch S: P13-S remains a separate scalar/source Wilson-basis branch and must not be discarded.
- Proton limits: The gate fixes the Fierz convention and global-charge audit, but physical proton limits require a canonical external-state crossing convention, physical flavor rotations, RG factors, hadronic matrix elements, and channel limits.

## Next Stage

`P16_choose_vector_physical_multiplet_or_return_to_scalar_source`

Either assemble physical hermitian vector multiplets and their cross-root Clebsch phases to search for canonical BNV vector operators, or return to the P13-S scalar/source recoupling path. Do not insert proton limits yet.
