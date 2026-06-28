# Route-A Lean Reproducibility Manifest

This manifest covers only the lean Route-A note. It does not claim to
reproduce the deferred full flavor, proton-decay, threshold, or source-sector
audits.

## Paper

- Source: `paper/gut_framework.tex`
- Bibliography: `paper/refs.bib`
- Compiled PDF: `paper/gut_framework.pdf`

Build from the `paper/` directory:

```sh
pdflatex gut_framework.tex
bibtex gut_framework
pdflatex gut_framework.tex
pdflatex gut_framework.tex
```

## Numerical Checks Kept In The Lean Note

1. D5 half-spin branching and hypercharge normalization:

```sh
python3 code/verify_d5_half_spin_hypercube.py
```

Expected output artifacts:

- `output/d5_half_spin/d5_half_spin_summary.json`
- `output/d5_half_spin/d5_half_spin_weights.csv`
- `output/d5_half_spin/d5_half_spin_report.md`

Core checks:

- 16 half-spin weights.
- Pati-Salam counts `(4,2,1): 8`, `(bar4,1,2): 8`.
- Field multiplicities `Q:6, L:2, u^c:3, d^c:3, nu^c:1, e^c:1`.
- `Tr Y^2 = 10/3`, `Tr T3L^2 = 2`, `k_Y = 5/3`.

2. Route-B algebraic hidden-zeta check:

```sh
python3 code/verify_hidden_zeta_origin.py
```

Expected output artifact:

- `output/hidden_zeta_origin/routeB_hidden_zeta_verification.json`

This verifies only the algebraic Schur-complement matching data for
`zeta K_tr`; it is not a full phenomenology audit.

## Deferred Audit Scaffold Cards

These files are reproducible scaffolds for future companion audits. They are
not publication-grade flavor, proton-decay, threshold, or source-sector
phenomenology results.

```sh
python3 code/audit4a1_cmsgut_vacuum_branches.py
python3 code/audit4a1_cmsgut_mass_export.py
python3 code/audit4a1_cmsgut_literature_mass_import.py
```

Expected output artifacts:

- `output/audit4a1/vacuum_branches.json`
- `output/audit4a1/vacuum_branch_report.md`
- `output/audit4a1/conventions_diff.md`
- `output/audit4a1/mass_export_schema.json`
- `output/audit4a1/mass_export_schema.md`
- `output/audit4a1/literature_mass_matrices.json`
- `output/audit4a1/literature_mass_matrices.md`
- `output/audit4a1/triplet_symbolic_inverse.json`
- `output/audit4a1/triplet_symbolic_inverse.md`
- `output/audit2/source_basis_wilson_contract.json`
- `output/audit2/source_basis_wilson_contract.md`

Core scaffold checks:

- Aulakh-Girdhar vacuum cubic convention is fixed.
- Source-named SU(5), flipped-SU(5), and `G_LR` special points pass.
- The generic `Spin(10) -> G_SM` broken-generator count is `45-(8+3+1)=33`.
- Aulakh-Girdhar doublet/triplet mass entries are transcribed with source-line
  anchors.  The local BMSV note is recorded as qualitative cross-reference
  rather than a second matrix-entry table.
- All mixed chiral/gauge Goldstone smoke gates `G,E,F,J,X` and a numerical
  triplet inverse gate pass on a generic F-flat sample.
- The Audit-2-required triplet inverse entries
  `S_1^1`, `S_1^2`, `S_2^1`, `S_2^2`, `S_1^4`, and `S_2^4` are exported as
  hand-audited Leibniz cofactor expansions.  The determinant expansion has
  `120` terms, each required minor has `24` terms, and the numeric gate gives
  det error `7.994039446406413e-15` with max inverse-entry error
  `1.0103182026100664e-15`.
- Audit 2 now has a source-basis `C5L/C5R` Wilson-tensor contract consuming
  those six symbolic inverse entries.  The random-tensor gate gives
  `C5L` and `C5R` symbolic/direct max errors `0.0`; physical flavor rotations,
  dressing, lattice matrix elements, and channel widths remain deferred.
- Scalar-Hessian Goldstone directions and the non-placeholder heavy spectrum
  remain pending.

## Status Log

- Roadmap and task ledger: `roadmap.md`
