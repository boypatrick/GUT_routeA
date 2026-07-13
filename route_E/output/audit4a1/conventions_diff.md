# Audit 4a.1 CMSGUT Conventions-Diff

This file is intentionally unresolved.  It records the sign and
normalization items that must be closed before the CMSGUT-like branch can
claim a literature-compatible vacuum cubic or heavy spectrum.

## Pati-Salam Singlet VEV Dictionary

| symbol | field | PS block | role |
| --- | --- | --- | --- |
| `p` | `210_H` | `(1,1,1)` | Pati-Salam singlet vev |
| `a` | `210_H` | `(15,1,1)` | B-L-like adjoint vev |
| `omega` | `210_H` | `(15,1,3)` | right-triplet-aligned adjoint vev |
| `sigma` | `126_H` | `(overline{10},1,3)` | B-L breaking source direction |
| `bar_sigma` | `overline{126}_H` | `(10,1,3)` | post-B-L P_nu^c source and D-flat conjugate |

## D-Flatness Convention

- Required condition: `|sigma| = |bar_sigma|`.
- Working gauge convention allowed for stage 1: `sigma = bar_sigma` after
  fixing the conjugate source phase.
- This is only a convention statement; it is not a solved F/D-flat vacuum.

## Unresolved Convention Items

- `normalization_of_p_a_omega`: unresolved; must be mapped to the chosen BMSV/Aulakh-style table (F-term cubic coefficients and mass matrices are convention-sensitive.)
- `sigma_bar_sigma_phase`: stage-1 convention: record |sigma|=|bar_sigma| and allow sigma=bar_sigma after gauge fixing (D-flatness fixes magnitudes, while the relative phase is a gauge/source convention.)
- `branch_variable_x_and_xi`: fixed to Aulakh-Girdhar 2005: x=-lambda*omega/m, xi=lambda*M/(eta*m) (The literature cubic P3(x;xi)=0 has multiple equivalent definitions.)
- `cubic_polynomial_coefficients`: fixed to 8*x^3 - 15*x^2 + 14*x - 3 = -xi*(1-x)^2 (Acceptance depends on reproducing enhanced-symmetry limits, not on a custom cubic.)
- `special_point_values`: source-named SU(5), flipped SU(5), and G_LR points validated; no full-Pati-Salam point claimed (SU(5), flipped SU(5), and left-right/Pati-Salam-related limits are the stage-1 validation anchors.)
- `BMSV_vs_Aulakh_symbol_map`: conventions-diff file created; unresolved until dual-source comparison is filled (Historical spectrum tables differ by signs and normalizations.)

## Literature Special-Point Acceptance

Stage 1 reproduces the source-named enhanced-symmetry special points
of the chosen Aulakh-Girdhar convention.  The source labels `xi=3` as
`G_LR`; this scaffold therefore does not claim a full Pati-Salam
special point.

| point | source name | x | xi | residual | pass |
| --- | --- | --- | --- | --- | --- |
| `SU5_x_half` | SU(5) | 1/2 | -5 | 0 | True |
| `SU5_x_minus_one` | SU(5) | -1 | 10 | 0 | True |
| `GLR_x_zero` | G_LR | 0 | 3 | 0 | True |
| `flipped_SU5_x_third` | flipped SU(5) x U(1) | 1/3 | -2/3 | 0 | True |
| `full_Pati_Salam_requested_anchor` | full Pati-Salam | None | None | None | None |

## Boundary

No threshold, proton-decay, or full heavy-spectrum conclusion may use this
stage-1 scaffold until the cubic coefficients, special-point values,
Goldstone count, and triplet inverse block are exported by a later audit.
