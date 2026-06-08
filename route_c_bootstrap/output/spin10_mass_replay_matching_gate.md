# Route-C P11 Staged Spin(10) Mass Replay

P11 attaches P10 staged mass blocks to the matching gates.  It does not convert chiral-pair symbolic poles into adjoint-vector exchange unless the mediator quantum numbers pass the adjoint-current gate, and it does not compute numerical proton lifetimes.

## Mass Replacement Dictionary

| symbolic mass | P10 replacement |
| --- | --- |
| `M_(6,2,2)` | `sqrt(kappa_PS) g_10 v_PS` |
| `M_LQ` | `sqrt(2/3) g_4 v_4` |
| `M_WR` | `g_R v_R` |
| `M_Zprime` | `0.5 v_Y sqrt(g_R^2+g_BL^2)` |

## Summary

| quantity | value |
| --- | ---: |
| current_vector_mass_blocks_available | 4 |
| p7_rows_replayed | 18 |
| p7_operator_status_counts | `{"BNV_with_sterile_neutrino": 2, "B_and_L_conserving": 12, "standard_BNV_seed": 4}` |
| p7_adjoint_mass_gate_counts | `{"not_an_adjoint_current_mass_match": 18}` |
| p7_rows_with_candidate_adjoint_mass_blocks | 0 |
| p7_rows_remaining_scalar_or_source_mass_debt | 18 |
| baryon_violating_or_sterile_bnv_rows | 6 |
| physical_proton_bounds_evaluable_now | 0 |

## Current-Vector Replay

| mass block | maps | status |
| --- | ---: | --- |
| `M_(6,2,2)` | 24 | mass_denominator_supplied; Goldstone/Higgs amplitude and spin numerator still required for a full high-energy check |
| `M_LQ` | 6 | mass_denominator_supplied; Goldstone/Higgs amplitude and spin numerator still required for a full high-energy check |
| `M_WR` | 2 | mass_denominator_supplied; Goldstone/Higgs amplitude and spin numerator still required for a full high-energy check |
| `M_Zprime` | 1 | neutral mass denominator supplied; physical charges and normalization must be chosen before low-energy matching |

## P7 Chiral-Pair Mass Gate

| operator | mediator | gate | P11 denominator status |
| --- | --- | --- | --- |
| `QQ u^c d^c` | `X[(bar6,1);Y=-1/3;B-L=-2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QQQL` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QQ u^c d^c` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QQQL` | `X[(3,3);Y=-1/3;B-L=-2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL u^c e^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL d^c nu^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QQ u^c d^c` | `X[(1,2);Y=1/2;B-L=0]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL u^c e^c` | `X[(1,2);Y=1/2;B-L=0]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QQ u^c d^c` | `X[(8,2);Y=1/2;B-L=0]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL d^c nu^c` | `X[(1,2);Y=-1/2;B-L=0]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL d^c nu^c` | `X[(bar3,2);Y=-1/6;B-L=-4/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `QL u^c e^c` | `X[(bar3,2);Y=-7/6;B-L=-4/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `LL e^c nu^c` | `X[(1,1);Y=1;B-L=2]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `LL e^c nu^c` | `X[(1,2);Y=1/2;B-L=0]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `u^c u^c d^c e^c` | `X[(bar3,1);Y=4/3;B-L=2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `u^c u^c d^c e^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `u^c d^c d^c nu^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |
| `u^c d^c d^c nu^c` | `X[(3,1);Y=2/3;B-L=-2/3]` | `not_an_adjoint_current_mass_match` | scalar/source mass debt |

## Interpretation

- improved: The adjoint current-vector sectors now carry explicit staged Spin(10) mass denominators M_(6,2,2), M_LQ, M_WR, and M_Zprime.
- not improved: The P7 chiral-pair pole rows do not automatically become adjoint current-vector exchanges.  Rows failing the adjoint gate remain scalar/source mass debt.
- next: Either rewrite the relevant BNV amplitudes in a completed current-current vector basis, or choose scalar/source masses for the chiral-pair pole rows and then supply flavor/RG/hadronic inputs.

## Next Stage

`P12_choose_vector_current_or_scalar_source_proton_branch`

Decide whether proton matching proceeds through completed Spin(10) current-vector exchange or through scalar/source chiral-pair poles; the two branches require different Wilson basis maps.
