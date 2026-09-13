# Lower finite scalar determinant subtraction

**Actual finite constant-background scalar result; full lower covariant matching remains open.**

Checks: 22/22.

The 303-dimensional upper-vector quotient has 290 positive scalar poles and 13 retained zeros. Its exact quadratic Schur determinant factorizes into 55 upper scalar poles and the nonlocal 248-dimensional PS operator. The local two-derivative PS operator is not interchangeable with it inside a loop trace.

```json
{
  "V_full": -0.04307001295455113,
  "V_C": -0.005964439308191641,
  "V_PS_nonlocal": -0.03710557364635949,
  "V_PS_local": -0.03705904934774003,
  "finite_local_truncation_difference": -4.652429861946167e-05
}
```

Values are in omega^4 reference units. The nonlocal-minus-local difference is a finite matching contribution in the declared scalar determinant, not a particle mass prediction. A constant vacuum term alone does not provide its source derivatives or the required Wilson coefficients.

| p_E^2 | Exact logdet defect | Nonlocal minus local logdet |
|---|---:|---:|
| 0.0001 | 3.979e-13 | -4.2138765e-08 |
| 0.003 | -7.958e-13 | -3.39588e-05 |
| 0.03 | -2.274e-13 | -0.0019212543 |
| 0.3 | 3.695e-13 | -0.031932074 |
| 3.0 | -2.842e-13 | -0.090942221 |
| 30.0 | 9.095e-13 | -0.10773673 |

## Missing components

- background/source derivatives of the finite functional
- transverse gauge and ghost/Goldstone package
- upper loop Wilson insertion and lower EFT subtraction
- finite kinetic/Yukawa/box operator matching

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| upper_heavy_basis_horizontal_at_broken_background | 0.000e+00 | True |
| full_kinetic_block_structure | 1.688e-16 | True |
| quotient_full_290_positive_scalar_poles | 7.606e-16 | True |
| full_303_chart_has_13_zero_modes | 0.000e+00 | True |
| local_248_chart_has_13_zero_modes | 0.000e+00 | True |
| local_metric_includes_heavy_response | 1.056e-16 | True |
| p2=0.0001_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=0.0001_exact_nonlocal_logdet_subtraction | 3.979e-13 | True |
| p2=0.003_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=0.003_exact_nonlocal_logdet_subtraction | 7.958e-13 | True |
| p2=0.03_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=0.03_exact_nonlocal_logdet_subtraction | 2.274e-13 | True |
| p2=0.3_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=0.3_exact_nonlocal_logdet_subtraction | 3.695e-13 | True |
| p2=3.0_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=3.0_exact_nonlocal_logdet_subtraction | 2.842e-13 | True |
| p2=30.0_positive_Euclidean_kernels | 0.000e+00 | True |
| p2=30.0_exact_nonlocal_logdet_subtraction | 9.095e-13 | True |
| full_generalized_spectrum_logdet | 2.021e-15 | True |
| finite_subtraction_scale_derivative | 2.788e-14 | True |
| local_determinant_is_not_exact | 0.000e+00 | True |
| finite_scalar_subtraction_broken_PS_covariance | 4.163e-17 | True |
