# Route-C P10 Spin(10) Breaking Mass-Matrix Ledger

P10 chooses a staged orthogonal source branch and supplies a diagonalizable broken-vector mass matrix.  It does not specify a complete scalar potential, threshold spectrum, flavor fit, or proton lifetime.

## Branch

- name: `staged_orthogonal_source_branch`
- reason: It is the minimal conservative branch that turns the P9 symbolic mass inputs into positive diagonal blocks while preserving the hypercharge zero mode.

Breaking chain:

- Spin(10) -> SU(4)_C x SU(2)_L x SU(2)_R
- SU(4)_C -> SU(3)_C x U(1)_{B-L}
- SU(2)_R -> U(1)_{T3R}
- U(1)_{T3R} x U(1)_{B-L} -> U(1)_Y

## Symbolic Mass Blocks

| block | dimension | mass-squared |
| --- | ---: | --- |
| (6,2,2) | 24 | `M_(6,2,2)^2 = kappa_PS g_10^2 v_PS^2` |
| SU(4)_C leptoquark roots | 6 | `M_LQ^2 = (2/3) g_4^2 v_4^2` |
| SU(2)_R charged roots | 2 | `M_WR^2 = g_R^2 v_R^2` |
| neutral U(1) mixing | 2 | `v_Y^2/4 [[g_R^2, -g_R g_BL],[-g_R g_BL, g_BL^2]]` |

## Demo Diagonalization

| quantity | value |
| --- | ---: |
| active matrix dimension, including hypercharge zero mode | 34 |
| positive broken eigenvalue count | 33 |
| zero eigenvalue count | 1 |
| M_(6,2,2)^2 demo | 1.0 |
| M_LQ^2 demo | 0.666666666666667 |
| M_WR^2 demo | 1.0 |
| neutral eigenvalues demo | [0.0, 0.5] |

The massless neutral eigenvector is the hypercharge gauge boson:

```text
B_Y proportional to g_BL W_R^3 + g_R B_{B-L}
g_Y = g_R g_BL / sqrt(g_R^2 + g_BL^2)
```

## Verification

| check | value |
| --- | ---: |
| p9_generator_layer_passes | True |
| all_p9_generator_maps_assigned_to_mass_blocks | True |
| positive_broken_count_matches_spin10_to_sm | True |
| one_hypercharge_zero_mode | True |
| neutral_block_determinant_demo | 0.0 |
| neutral_block_trace_demo | 0.5 |
| p10_mass_matrix_layer_passes | True |

## P6/P7 Replay Interface

| symbolic mass | replacement |
| --- | --- |
| `M_(6,2,2)` | `sqrt(kappa_PS) g_10 v_PS` |
| `M_LQ` | `sqrt(2/3) g_4 v_4` |
| `M_WR` | `g_R v_R` |
| `M_Zprime` | `0.5 v_Y sqrt(g_R^2+g_BL^2)` |

Still required:

- choose numerical scales/couplings or scan ranges
- choose scalar/source branch for the 8 non-adjoint P4 pairs if activated
- provide physical flavor rotations
- insert RG factors and hadronic matrix elements
- insert current experimental proton lifetime limits

## Next Stage

`P11_replay_P6_P7_with_staged_spin10_masses`

Choose symbolic or numerical ranges for v_PS, v_4, v_R, v_Y and couplings, then rerun the P6/P7 matching gates with these mass blocks.
