# Audit 4a.1 Triplet Symbolic Inverse Entries

This artifact expands the Aulakh--Girdhar 5x5 proton-decay triplet
inverse contract into explicit Leibniz cofactor entries required by
Audit 2.  It is hand-audited symbolic data, not a CAS-simplified
polynomial.

## Digest

- card sha256: `177df81d1f6a9fed2cb284aed989f28401261831fdba0d97f9e2d77252bd6801`
- source literature card sha256: `8c5af27aeb06d7d094c8f7d1d22d9032c11bfe1bca8e0b2bbe41cfea90d5e2f2`

## Numeric Gate

- determinant term count: `120`
- det expansion absolute error: `7.994e-15`
- inverse-entry max absolute error: `1.010e-15`
- pass: `True`

## Audit-2 Required Entries

| entry | minor | terms | abs error | pass |
| --- | --- | --- | --- | --- |
| `S_1^1` | remove row `1`, col `1` | `24` | `2.259e-16` | `True` |
| `S_1^2` | remove row `2`, col `1` | `24` | `4.163e-17` | `True` |
| `S_2^1` | remove row `1`, col `2` | `24` | `2.861e-16` | `True` |
| `S_2^2` | remove row `2`, col `2` | `24` | `1.637e-17` | `True` |
| `S_1^4` | remove row `4`, col `1` | `24` | `1.010e-15` | `True` |
| `S_2^4` | remove row `4`, col `2` | `24` | `7.228e-16` | `True` |

## Boundary

- The formulas are exact finite Leibniz expansions in the source-anchored
  triplet matrix entries `T_ij`.
- The artifact does not simplify or factor the polynomials; that remains
  optional CAS polish rather than an Audit-2 blocker.
- Scalar-Hessian Goldstone directions and the non-placeholder heavy
  spectrum remain separate Audit 4a gates.
