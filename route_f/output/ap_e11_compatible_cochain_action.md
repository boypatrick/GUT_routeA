# AP-E11 compatible tetrahedral cochain action

Status: `production_compatible_action_density_and_quotient_audit`.

## Exact algebra

- d^2 residual: `0.000e+00`.
- Leibniz residual: `0.000e+00`.
- third-cup determinant residual: `8.327e-16`.
- minimum Whitney-Hodge eigenvalue: `2.485453e-02`.

## Cell formula

`Q_R,K(A)=1/2 |A|^2 + R/2 |wedge^2 A|^2 + K/2 |wedge^3 A|^2`; the periodic corrector is exactly zero on all tested conforming meshes.

## Quotient scan

- cases: `13`; all pass: `True`.
- production coverage: `True`.
- threshold pass: `True`.
- quotient continuum gate: `True`.

## Authorization

- relaxed fixed-degree Gamma theorem: `True`.
- classical unrelaxed density theorem: `False`.
- relaxed regulator stable: `True`.
- regulator stable: `False`.
- same-action Hessian gate: `False`.
- determinant-variation gate: `False`.

Checks: `11/11`.
