# Audit 4a.1 CMSGUT Vacuum-Branch Stage-1 Scaffold

This artifact fixes the Pati-Salam singlet-vev dictionary and the
literature-special-point acceptance contract for the CMSGUT-like branch.
It does not export a heavy spectrum.

## Digest

- card sha256: `fc6bfdddb5de973a469f0e1b00d8f5c82f9dae2c9f8335f35c4c76843ed4bb78`
- upstream Audit 4a sha256: `be735ff9f65a62f53a7cd1cbf504b4113d0f65dcaec553a4f6bb55e410471fbc`
- upstream Audit 0 sha256: `a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc`

## Registered Singlet VEVs

| symbol | field | Pati-Salam block | role |
| --- | --- | --- | --- |
| `p` | `210_H` | `(1,1,1)` | Pati-Salam singlet vev |
| `a` | `210_H` | `(15,1,1)` | B-L-like adjoint vev |
| `omega` | `210_H` | `(15,1,3)` | right-triplet-aligned adjoint vev |
| `sigma` | `126_H` | `(10,1,3)` | B-L breaking source direction |
| `bar_sigma` | `overline{126}_H` | `(overline{10},1,3)` | post-B-L P_nu^c source and D-flat conjugate |

## Stage Gates

| gate | value |
| --- | --- |
| `ps_singlet_vev_registry_complete` | `True` |
| `d_flatness_condition_recorded` | `True` |
| `cubic_coefficients_fixed_from_literature` | `False` |
| `special_points_reproduced` | `False` |
| `symbolic_mass_matrices_exported` | `False` |
| `goldstone_count_33_exported` | `False` |
| `doublet_mass_matrix_exported` | `False` |
| `triplet_inverse_block_exported` | `False` |
| `heavy_spectrum_json_nonplaceholder` | `False` |

## Downstream Boundary

- Audit 3 remains blocked until a non-placeholder heavy spectrum exports
  masses and beta vectors.
- Audit 2 remains blocked until physical flavor rotations and the colored
  triplet inverse block are exported.
