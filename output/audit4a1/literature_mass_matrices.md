# Audit 4a.1 CMSGUT Literature Mass-Matrix Import

This artifact transcribes the Aulakh--Girdhar doublet/triplet mass
matrices and runs the first numerical gates needed before Audit 2.

## Digest

- card sha256: `8c5af27aeb06d7d094c8f7d1d22d9032c11bfe1bca8e0b2bbe41cfea90d5e2f2`
- upstream mass-export schema sha256: `94f14ed3cb2eac9cfa582a59dd938e9aeca9b1d2d728f453294e9faa4e4ebea1`

## Imported Blocks

| block | shape | source lines | status |
| --- | --- | --- | --- |
| `H_doublet_4x4` | `[4, 4]` | `2601-2624` | transcribed_from_aulakh_girdhar_appendix |
| `T_triplet_5x5` | `[5, 5]` | `2629-2648` | transcribed_from_aulakh_girdhar_appendix |
| `G_mixed_neutral_6x6` | `[6, 6]` | `2654-2669` | transcribed_from_aulakh_girdhar_appendix |
| `E_mixed_4x4` | `[4, 4]` | `2671-2685` | transcribed_from_aulakh_girdhar_appendix |
| `F_mixed_3x3` | `[3, 3]` | `2688-2699` | transcribed_from_aulakh_girdhar_appendix |
| `J_mixed_4x4` | `[4, 4]` | `2701-2714` | transcribed_from_aulakh_girdhar_appendix |
| `X_mixed_3x3` | `[3, 3]` | `2718-2731` | transcribed_from_aulakh_girdhar_appendix |

## Numeric Gates

- Mixed chiral/gauge Goldstone gates all pass: `True`.
- Triplet inverse gate pass: `True`.
- Triplet inverse residual norm: `2.243e-15`.
- Triplet condition number: `4.037e+01`.

| mixed block | chiral null ratio | full determinant | pass |
| --- | --- | --- | --- |
| `G_mixed_neutral_6x6` | `5.460e-33` | `{'re': 19.126896386195, 'im': 0.0}` | `True` |
| `E_mixed_4x4` | `7.951e-17` | `{'re': -9.664268911298, 'im': 0.0}` | `True` |
| `F_mixed_3x3` | `0.000e+00` | `{'re': 0.85034843469, 'im': 0.0}` | `True` |
| `J_mixed_4x4` | `2.065e-17` | `{'re': -4.848831190909, 'im': 0.0}` | `True` |
| `X_mixed_3x3` | `1.869e-17` | `{'re': -3.541563628851, 'im': 0.0}` | `True` |

Audit-2-required numeric inverse entries on the smoke sample:

| entry | value |
| --- | --- |
| `S_1^1` | `1.239746430932 + 0.047498160261 i` |
| `S_1^2` | `-0.06736333461 + -0.017947220088 i` |
| `S_2^1` | `0.251612633761 + 0.067035685527 i` |
| `S_2^2` | `-0.012961357573 + -0.006788367186 i` |
| `S_1^4` | `0.162173968105 + 0.375375135201 i` |
| `S_2^4` | `0.015848221077 + -0.326618671218 i` |

## Boundary

- The BMSV local note in this workspace gives only qualitative
  proton-decay/triplet-mixing guidance; it is recorded as a
  cross-reference, not as a second matrix-entry table.
- The triplet inverse exported here is a numerical smoke inversion.
  The symbolic rational inverse is specified by an adjugate/cofactor
  contract and remains pending expansion.
- The Goldstone gate verifies all G,E,F,J,X mixed chiral/gauge blocks;
  scalar Hessian Goldstone directions remain pending.
