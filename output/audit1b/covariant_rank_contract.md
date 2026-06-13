# Audit 1b Covariant-Rank Contract

This is a design/rank check, not a flavor fit.

- card sha256: `1bd7d80dd8cd4896bbfa51960c6d3301d09acf581226726e58ab7786bac3154b`
- all rank checks pass: `True`

## Rank Rows

| k | ambient C dim | rank | expected | constraints C | expected | fiber C |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 5 | 5 | 0 | 0 | 2 |
| 2 | 10 | 7 | 7 | 3 | 3 | 2 |
| 3 | 15 | 9 | 9 | 6 | 6 | 2 |
| 4 | 20 | 11 | 11 | 9 | 9 | 2 |
| 5 | 25 | 13 | 13 | 12 | 12 | 2 |

## Design Implications

- Do not run the literature `(h,f)` two-tensor test as a pass/fail audit when the family basis is unknown; its net real constraints are negative after scanning a full `U(3)` orientation.
- Use the internal companion test with at least three aligned family tensors after the `K_tr` contact direction fixes the basis up to a residual `O(3)`.
- Report tensor-by-tensor residual increments; a pooled residual can hide whether the expected 6 real constraints per added tensor are doing work.
