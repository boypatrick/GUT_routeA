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
| `sigma` | `126_H` | `(10,1,3)` | B-L breaking source direction |
| `bar_sigma` | `overline{126}_H` | `(overline{10},1,3)` | post-B-L P_nu^c source and D-flat conjugate |

## D-Flatness Convention

- Required condition: `|sigma| = |bar_sigma|`.
- Working gauge convention allowed for stage 1: `sigma = bar_sigma` after
  fixing the conjugate source phase.
- This is only a convention statement; it is not a solved F/D-flat vacuum.

## Unresolved Convention Items

- `normalization_of_p_a_omega`: unresolved; must be mapped to the chosen BMSV/Aulakh-style table (F-term cubic coefficients and mass matrices are convention-sensitive.)
- `sigma_bar_sigma_phase`: stage-1 convention: record |sigma|=|bar_sigma| and allow sigma=bar_sigma after gauge fixing (D-flatness fixes magnitudes, while the relative phase is a gauge/source convention.)
- `branch_variable_x_and_xi`: placeholder; exact definitions must be imported before any special-point pass (The literature cubic P3(x;xi)=0 has multiple equivalent definitions.)
- `cubic_polynomial_coefficients`: not fixed in this scaffold (Acceptance depends on reproducing enhanced-symmetry limits, not on a custom cubic.)
- `special_point_values`: required from literature table; no pass claimed here (SU(5), flipped SU(5), and Pati-Salam limits are the stage-1 validation anchors.)
- `BMSV_vs_Aulakh_symbol_map`: conventions-diff file created; unresolved until dual-source comparison is filled (Historical spectrum tables differ by signs and normalizations.)

## Literature Special-Point Acceptance

Stage 1 must reproduce the enhanced-symmetry special points from the
chosen literature convention before any CMSGUT spectrum is accepted.

| point | required status |
| --- | --- |
| `SU5_unflipped` | not evaluated; literature value required |
| `flipped_SU5` | not evaluated; literature value required |
| `Pati_Salam` | not evaluated; literature value required |

## Boundary

No threshold, proton-decay, or full heavy-spectrum conclusion may use this
stage-1 scaffold until the cubic coefficients, special-point values,
Goldstone count, and triplet inverse block are exported by a later audit.
