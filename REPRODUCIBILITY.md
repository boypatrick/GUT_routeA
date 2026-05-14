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

## Status Log

- Roadmap and task ledger: `roadmap.md`
