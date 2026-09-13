# Upper factorizable scalar-exchange Wilson matching

**One-loop scalar-exchange subset, not a complete four-fermion matching.**

Checks: 16/16.

All 24 Yukawa-coupled real upper six modes are retained, including the ten real PS-invariant mass mixings between their four copies. Actual scalar/vector mass diagrams, the broken-anchored bosonic counterterms and the hard heavy-valley response enter the same propagator. Pure-Yukawa/vector proper vertices and external fermion legs enter J. Heavy scalar external-leg normalization is not inserted twice.

C0=J^T D^-1 J; delta C=delta J^T D^-1 J+J^T D^-1 delta J-J^T D^-1 delta D D^-1 J.

| Synthetic one-family case | ||C0|| | ||delta C|| | Ratio |
|---|---:|---:|---:|
| 0 | 0.43315694 | 0.019121882 | 0.044145389 |
| 1 | 0.24196642 | 0.087235788 | 0.36052849 |

These are diagnostic source kernels, not proton lifetimes or fitted physical amplitudes. The ratio alone is not a full perturbativity assessment because important diagrams and parameter feedback remain absent.

## Non-small mass insertion: do not promote to a fit

The spectral radius of D^(-1/2) deltaD D^(-1/2) is 6.7559262. Its Neumann expansion does not converge at unit loop parameter for this computed subset, even though the resummed subset mass is positive. Smaller corrections in selected source projections cannot certify the whole propagator. Missing diagrams, common parameter feedback and the declared truncation must be resolved before judging the complete model. The algebra checks below do not assert perturbative control.

Contributions on the dominant normalized mass direction (not separate eigenvalues):

```json
{
  "scalar": -0.5389464501617361,
  "vector": -0.016986813080004683,
  "common_bosonic_CT": 7.346373994963831,
  "hard_plus_CT_valley_response": -0.03451453761034232
}
```

## Missing contributions

- nonfactorizable scalar/vector/fermion boxes and direct vector exchange
- complete Lorentz/gauge/family operator basis and mixing
- fermionic fixed-VEV counterterms and final light retuning
- upper scalar three/four-point potential Wilson vertices
- controlled retained-mass corrections

## Checks

| Check | Residual | Pass |
|---|---:|:---:|
| 24_Yukawa_coupled_upper_real_modes | 0.000e+00 | True |
| ten_real_symmetric_six_copy_invariants | 3.386e-16 | True |
| all_six_copy_intertwiners_PS_covariant | 3.695e-15 | True |
| independent_full_six_mass_reconstruction | 7.425e-15 | True |
| six_linear_VV_mass_vertices_vanish | 1.112e-14 | True |
| same_action_tree_six_mass_block | 7.556e-16 | True |
| case0_tree_C_positive_semidefinite | 0.000e+00 | True |
| case0_finite_C_transpose_symmetry | 9.858e-18 | True |
| case0_exact_resolvent_finite_difference | 2.012e-08 | True |
| case0_heavy_kinetic_redefinition_cancels_in_C | 1.743e-17 | True |
| case0_self_energy_insertion_is_not_zero | 0.000e+00 | True |
| case1_tree_C_positive_semidefinite | 0.000e+00 | True |
| case1_finite_C_transpose_symmetry | 1.458e-17 | True |
| case1_exact_resolvent_finite_difference | 5.137e-08 | True |
| case1_heavy_kinetic_redefinition_cancels_in_C | 3.729e-17 | True |
| case1_self_energy_insertion_is_not_zero | 0.000e+00 | True |
