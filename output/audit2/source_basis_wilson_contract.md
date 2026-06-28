# Audit 2 Source-Basis Wilson Contract

This artifact consumes the Audit 4a.1 symbolic triplet inverse entries
and exports formal source-basis `C5L` and `C5R` Wilson-tensor contracts.
It is not yet a mass-basis proton-decay calculation.

## Digest

- card sha256: `ed7f07bd2fdfdec51982cbfa9461f007cb2f809d8f8e33cf248fb999c83703e1`
- symbolic inverse card: `177df81d1f6a9fed2cb284aed989f28401261831fdba0d97f9e2d77252bd6801`

## Formal Contract

- `C5L`: `C5L_{ab,cd}=sum_{(i,j)} S_i^j (Y_QQ^i)_{ab} (Y_QL^j)_{cd}`
- `C5R`: `C5R_{ab,cd}=sum_{(i,j)} S_i^j (Y_UE^i)_{ab} (Y_UD^j)_{cd}`

## Source Entries

| entry | left column source | right row source | source error |
| --- | --- | --- | --- |
| `S_1^1` | `t1 = H_4` | `bar_t1 = H^4` | `2.259e-16` |
| `S_1^2` | `t1 = H_4` | `bar_t2 = overlineSigma_{(a)}^4` | `4.163e-17` |
| `S_2^1` | `t2 = overlineSigma_{(a)4}` | `bar_t1 = H^4` | `2.861e-16` |
| `S_2^2` | `t2 = overlineSigma_{(a)4}` | `bar_t2 = overlineSigma_{(a)}^4` | `1.637e-17` |
| `S_1^4` | `t1 = H_4` | `bar_t4 = Sigma^4_{R0}` | `1.010e-15` |
| `S_2^4` | `t2 = overlineSigma_{(a)4}` | `bar_t4 = Sigma^4_{R0}` | `7.228e-16` |

## Numeric Gate

- random seed: `20260616`
- `C5L` symbolic/direct max error: `0.000e+00`
- `C5R` symbolic/direct max error: `0.000e+00`
- pass: `True`

## Boundary

- This is a source-basis Wilson contract only.
- Physical flavor rotations, triplet Yukawa maps, SUSY dressing, lattice
  matrix elements, and channel widths remain deferred.
