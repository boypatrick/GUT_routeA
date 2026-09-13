# Upper scalar hard matching and common tadpole anchor

**Actual scalar kinetic/mass subset, not full upper Wilson or lower finite matching.**

Checks: 14/14.

The eight real PS-invariant scalar bilinears are derived from the existing intertwiners. Eight polynomial-action probes reconstruct the full 248x248 metric and mass correction; independent complex/mixed-parent directions test the reconstruction.

Scalar-cubic kinetic eigenvalue range: [9.388830851244331e-05, 0.00450338488471728]. Hard scalar mass eigenvalue range: [-0.14082030326785688, -0.011778343742355847].

Heavy-heavy and mixed-heavy-soft propagators are included at leading retained-mass hard order. Pure-soft graphs are not included in this hard coefficient. The gauge-fixed upper Goldstone directions are kept in the scalar state sum, not mistaken for physical retained PS scalars.

## Common fixed-VEV prescription

The pre-existing invariant mass counterterms are anchored at the broken vacuum and cannot be independently adjusted at the upper PS saddle. The following diagnostic uses the exact counterterm reference scale, not the nearby upper matching scale:

```json
{
  "mu": 0.5626028899853225,
  "scalar_hard": [
    -0.022006161825818504,
    -0.00210401894327551
  ],
  "vector_hard": [
    -0.01518933975705058,
    -0.0
  ],
  "counterterm": [
    0.135303251276017,
    0.00870588716765
  ],
  "residual": [
    0.09810774969314792,
    0.006601868224374491
  ],
  "radial_hessian": [
    [
      6.798242799348019,
      0.061076097553767035
    ],
    [
      0.061076097553767035,
      0.10345143557471
    ]
  ],
  "hard_plus_CT_linear_displacement": [
    -0.013931904939631713,
    -0.05559093314870954
  ],
  "full_upper_quantum_stationary_point": false,
  "reason": "retained EFT loops and fermionic fixed-VEV counterterms remain absent"
}
```

This displacement is only the hard-plus-counterterm contribution. Retained EFT-loop tadpoles are still necessary for a full upper quantum stationary background. It must not be substituted as a completed background into a fit.

## Missing items

- source-dependent three/four-point upper Wilson loops
- retained-mass expansion control
- upper EFT soft-loop tadpoles in common global prescription
- lower gauge-covariant finite matching
- finite matching with nonzero CHN and partially retained sterile fields

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| eight_PS_invariant_real_symmetric_bilinears | 5.914e-16 | True |
| full_heavy_soft_state_sum | 5.920e-16 | True |
| probe_matrix_invertible | 0.000e+00 | True |
| scalar_kinetic_positive_Gram | 0.000e+00 | True |
| matched_scalar_kinetic_metric_positive | 0.000e+00 | True |
| independent_kinetic_mass_invariant_reconstruction_0 | 2.109e-15 | True |
| independent_kinetic_mass_invariant_reconstruction_1 | 1.915e-15 | True |
| quartic_action_two_step_jet_identity | 1.129e-13 | True |
| full_PS_Ward_for_scalar_kinetic_and_mass | 3.866e-18 | True |
| upper_hard_tadpole_linear_shift_equation | 1.388e-17 | True |
| combined_bosonic_kinetic_metric_positive | 0.000e+00 | True |
| canonical_mass_is_real_symmetric | 5.816e-20 | True |
| h_scalar_leg_six_invariant_closure | 2.453e-18 | True |
| f_scalar_leg_six_invariant_closure | 7.078e-16 | True |
