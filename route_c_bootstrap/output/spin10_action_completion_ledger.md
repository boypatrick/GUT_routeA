# Route-C P9 Spin(10)-Locked Action-Completion Ledger

P9 starts action-level completion for the P8-leading Spin(10) branch.  It supplies generator and interface data, but it does not yet choose a full Higgs potential, mediator spectrum, flavor fit, or numerical proton-bound evaluation.

## Summary

| quantity | value |
| --- | ---: |
| basis dimension | 16 |
| broken generator sparse maps | 32 |
| SU(4) leptoquark root maps | 6 |
| SU(2)_R charged root maps | 2 |
| Spin(10)/Pati-Salam off-face Clebsch maps | 24 |
| P4 broken transition pairs | 30 |
| Spin(10) adjoint vector transition pairs | 22 |
| non-adjoint P4 completion pairs | 8 |
| generated pairs are subset of P4 broken pairs | True |
| P9 generator layer passes | True |

## Generator Families

| family | maps | role |
| --- | ---: | --- |
| Pati-Salam SU(4)_C leptoquark | 6 | completes Q/L and conjugate-sector leptoquark transitions |
| broken SU(2)_R charged current | 2 | completes u^c-d^c and nu^c-e^c transitions |
| Spin(10)/Pati-Salam off-face (6,2,2) | 24 | connects (4,2,1) to (bar4,1,2) via antisymmetric SU(4) Clebsch maps |

P9 intentionally does not force every P4 charge-bookkeeping transition to be
a Spin(10) adjoint gauge vector.  The non-adjoint P4 pairs remain
scalar/source/auxiliary completion candidates if a later branch activates
them.

## Goldstone/Higgs Completion Templates

| sector | broken generators | status |
| --- | --- | --- |
| Spin(10) -> Pati-Salam | (6,2,2) | template_only_order_parameter_not_chosen |
| SU(4)_C -> SU(3)_C x U(1)_{B-L} | SU(4) leptoquark roots E_{ell a}, E_{a ell} | template_only_order_parameter_not_chosen |
| SU(2)_R -> U(1)_{T3R} | T_R^+, T_R^- | template_only_order_parameter_not_chosen |

The Ward target for every broken vector remains

```text
p_mu M^mu(V_X) = M_X M(phi_X)
```

## Required Symbolic Inputs

| input | status |
| --- | --- |
| `g_10` | symbolic_required |
| `M_{(6,2,2)}` | symbolic_required |
| `M_{LQ}` | symbolic_required |
| `M_{W_R}` | symbolic_required |
| `M_{T}, Y_T` | not_supplied_in_p9 |

## Physical Flavor and Proton-Bound Boundary

- flavor status: `template_only_no_numeric_flavor_fit`
- proton-bound status: `P9 supplies the interface and Spin(10) generator data.  It does not insert numerical current experimental limits or compute lifetimes.`

P9 supplies the Spin(10) action-completion interface.  It does not claim
that proton lifetimes or full flavor observables are now evaluable.

## Next Stage

`P10_choose_spin10_breaking_branch_and_mass_spectrum`

Choose a concrete Spin(10)-breaking order parameter/source branch, then diagonalize broken-vector masses and replay P6/P7 with explicit masses, couplings, and flavor rotations.
