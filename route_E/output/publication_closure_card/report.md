# Publication closure candidate card

No web lookup was used.

## Selected branch

Row: `all10_radius_0p30_ckmheavy_1`.
Mode: `all10_radius_0p30_ckmheavy`.

## Recomputed local gates

| gate | pass |
|---|---:|
| `strict_ckm` | `True` |
| `loose_mass` | `True` |
| `seesaw_residual` | `True` |
| `future_1e35_d5_margin` | `True` |
| `post_spin10_source_symmetry` | `True` |
| `dterm_unitary_lock` | `True` |
| `hidden_radial_lock` | `True` |
| `endpoint_vectorlike_safe` | `True` |

## Numerical observables

CKM score: `7.634719e-04`.
Mass score: `1.490597e-01`.
Seesaw residual: `5.173996e-12`.
Future d=5 margin: `1.668117e+10`.

## Reproducibility diffs

| quantity | absolute diff |
|---|---:|
| `ckm_score_abs_diff` | 0.000000e+00 |
| `mass_score_abs_diff` | 0.000000e+00 |
| `seesaw_residual_abs_diff` | 0.000000e+00 |
| `d5_margin_abs_diff` | 0.000000e+00 |

## Verdict

Local source-consistent candidate passes: `True`.
Publication-level complete: `False`.

The source-consistent crossed-120 row closes the local strict CKM, mass, seesaw, and 1e35 yr d=5 stress gates when paired with the post-Spin(10) source projector and unitary/radial/endpoint completion audits.  This supersedes the older full_flavor_d5_pipeline bottleneck for the local branch, but it is still not a publication-final proof because the full channel table and final unified phenomenology input manifest are not yet rebuilt around this card.

Machine-readable outputs:

- `output/publication_closure_card/summary.json`
- `output/publication_closure_card/publication_closure_card.json`
- `output/publication_closure_card/input_manifest.csv`
