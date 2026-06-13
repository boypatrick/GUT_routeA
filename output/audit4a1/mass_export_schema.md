# Audit 4a.1 CMSGUT Mass-Export Stage-2 Schema

This artifact starts the mass-export layer after the vacuum-cubic
validator.  It proves only the group-theoretic Goldstone count and
exports the doublet/triplet block schema; it does not import the full
CMSGUT mass matrices.

## Digest

- card sha256: `94f14ed3cb2eac9cfa582a59dd938e9aeca9b1d2d728f453294e9faa4e4ebea1`
- upstream vacuum sha256: `8aaa5979c5dc0832ee7edecfd7531ded34a0dadd9c6118406487502880a240ff`
- upstream Audit 4a sha256: `be735ff9f65a62f53a7cd1cbf504b4113d0f65dcaec553a4f6bb55e410471fbc`

## Goldstone Count

- Parent dimension: `45`.
- Unbroken dimension: `12`.
- Broken dimension: `33`.
- Sector sum: `33`.
- Pass: `True`.

| sector | SM reps | count |
| --- | --- | --- |
| `SU4C_over_SU3C_U1BL_offdiagonal` | (3,1,2/3), (bar3,1,-2/3) | 6 |
| `SU2R_charged` | (1,1,1), (1,1,-1) | 2 |
| `orthogonal_U1R_BL` | (1,1,0)_orthogonal_to_Y | 1 |
| `PS_broken_bidoublet` | (3,2,1/6), (bar3,2,-1/6), (3,2,-5/6), (bar3,2,5/6) | 24 |

## Mass-Block Export Status

| block | shape | status |
| --- | --- | --- |
| `doublet_mass_matrix` | `[4, 4]` | schema_placeholder_pending_literature_import |
| `triplet_mass_matrix` | `[5, 5]` | schema_placeholder_pending_literature_import |
| `triplet_inverse_block` | `[5, 5]` | blocked_until_triplet_mass_entries_imported_and_inverted |

## Boundary

Audit 4a.1 now exports the Goldstone-count gate and the mass-block schema.
It still does not export CMSGUT mass entries, heavy beta vectors, a
non-placeholder heavy spectrum, or a dimension-five Wilson input.
