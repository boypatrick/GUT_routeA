# AP-E12 graph-density and regular-minimizer audit

Status: `production_endpoint_density_regular_minimizer_fail_closed`.

## Exact results

- degree continuous in the complete-minor graph norm: `True`.
- maximum singular-value identity residual: `3.553e-14`.
- maximum rank-one symbol residual: `2.842e-13`.
- maximum periodic strong-quasiconvex identity residual: `3.553e-15`.
- power-law Malý obstruction overlap: `False`.
- logarithmic-borderline overlap: `False`.

## Numerical evidence

- all re-relaxed backgrounds and perturbations pass: `True`.
- finite-grid isolation proxy: `True`.
- perturbed-start energy spread: `3.379369e-11`.
- maximum quotient profile distance: `5.702611e-07`.
- maximum translation/target-$SO(3)$ quotient field RMS: `1.307658e-05`.

## Fail-closed theorem boundary

- endpoint graph-L2 density proved: `False`.
- endpoint density disproved: `False`.
- local recovery for relaxed minimizer: `False`.
- partial regularity: `False`.
- classicality: `False`.
- continuum isolation: `False`.
- Hessian gate: `False`.

Checks: `9/9`.
