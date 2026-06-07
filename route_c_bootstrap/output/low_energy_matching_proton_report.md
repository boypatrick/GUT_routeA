# Route-C P7 Low-Energy Matching and Proton-Bound Readiness Report

P7 consumes the P2 symbolic pole ledger, P3 factorization ledger, and P6
high-energy-growth audit.  It builds conditional dimension-six matching
entries and marks every physical proton bound as conditional until the
missing completion/flavor/hadronic inputs are supplied.

## Summary

| quantity | value |
| --- | ---: |
| allowed symbolic poles processed | 18 |
| baryon-violating symbolic poles | 6 |
| BNV vector poles blocked by P6 completion | 6 |
| physical proton bounds evaluable now | 0 |
| all physical bounds conditional | True |

Operator status counts:

| status | count |
| --- | ---: |
| B_and_L_conserving | 12 |
| standard_BNV_seed | 4 |
| BNV_with_sterile_neutrino | 2 |

## Operator Groups

| operator | status | B | L | B-L | symbolic poles | bound evaluable |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `QQ u^c d^c` | B_and_L_conserving | 0 | 0 | 0 | 4 | False |
| `QQQL` | standard_BNV_seed | 1 | 1 | 0 | 2 | False |
| `QL u^c e^c` | B_and_L_conserving | 0 | 0 | 0 | 3 | False |
| `QL d^c nu^c` | B_and_L_conserving | 0 | 0 | 0 | 3 | False |
| `LL e^c nu^c` | B_and_L_conserving | 0 | 0 | 0 | 2 | False |
| `u^c u^c d^c e^c` | standard_BNV_seed | -1 | -1 | 0 | 2 | False |
| `u^c d^c d^c nu^c` | BNV_with_sterile_neutrino | -1 | -1 | 0 | 2 | False |

## Conditional Matching Rows

| external multiset | mediator | operator status | P6 gate | symbolic coefficient |
| --- | --- | --- | --- | --- |
| `Q Q d^c u^c` | `X[(bar6,1);Y=-1/3;B-L=-2/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QQ u^c d^c; X[(bar6,1);Y=-1/3;B-L=-2/3]] = g_{QQX} g_{u^cd^cX}^*/M_X[(bar6,1);Y=-1/3;B-L=-2/3]^2` |
| `L Q Q Q` | `X[(3,1);Y=-1/3;B-L=-2/3]` | standard_BNV_seed | blocked_as_vector_until_completion | `C6[QQQL; X[(3,1);Y=-1/3;B-L=-2/3]] = g_{QQX} g_{QLX}^*/M_X[(3,1);Y=-1/3;B-L=-2/3]^2` |
| `Q Q d^c u^c` | `X[(3,1);Y=-1/3;B-L=-2/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QQ u^c d^c; X[(3,1);Y=-1/3;B-L=-2/3]] = g_{QQX} g_{u^cd^cX}^*/M_X[(3,1);Y=-1/3;B-L=-2/3]^2` |
| `L Q Q Q` | `X[(3,3);Y=-1/3;B-L=-2/3]` | standard_BNV_seed | blocked_as_vector_until_completion | `C6[QQQL; X[(3,3);Y=-1/3;B-L=-2/3]] = g_{QQX} g_{QLX}^*/M_X[(3,3);Y=-1/3;B-L=-2/3]^2` |
| `L Q e^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QL u^c e^c; X[(bar3,1);Y=1/3;B-L=2/3]] = g_{QLX} g_{u^ce^cX}^*/M_X[(bar3,1);Y=1/3;B-L=2/3]^2` |
| `L Q d^c nu^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QL d^c nu^c; X[(bar3,1);Y=1/3;B-L=2/3]] = g_{QLX} g_{d^cnu^cX}^*/M_X[(bar3,1);Y=1/3;B-L=2/3]^2` |
| `Q Q d^c u^c` | `X[(1,2);Y=1/2;B-L=0]` | B_and_L_conserving | not_a_proton_seed | `C6[QQ u^c d^c; X[(1,2);Y=1/2;B-L=0]] = g_{Qu^cX} g_{Qd^cX}^*/M_X[(1,2);Y=1/2;B-L=0]^2` |
| `L Q e^c u^c` | `X[(1,2);Y=1/2;B-L=0]` | B_and_L_conserving | not_a_proton_seed | `C6[QL u^c e^c; X[(1,2);Y=1/2;B-L=0]] = g_{Qu^cX} g_{Le^cX}^*/M_X[(1,2);Y=1/2;B-L=0]^2` |
| `Q Q d^c u^c` | `X[(8,2);Y=1/2;B-L=0]` | B_and_L_conserving | not_a_proton_seed | `C6[QQ u^c d^c; X[(8,2);Y=1/2;B-L=0]] = g_{Qu^cX} g_{Qd^cX}^*/M_X[(8,2);Y=1/2;B-L=0]^2` |
| `L Q d^c nu^c` | `X[(1,2);Y=-1/2;B-L=0]` | B_and_L_conserving | not_a_proton_seed | `C6[QL d^c nu^c; X[(1,2);Y=-1/2;B-L=0]] = g_{Qd^cX} g_{Lnu^cX}^*/M_X[(1,2);Y=-1/2;B-L=0]^2` |
| `L Q d^c nu^c` | `X[(bar3,2);Y=-1/6;B-L=-4/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QL d^c nu^c; X[(bar3,2);Y=-1/6;B-L=-4/3]] = g_{Qnu^cX} g_{Ld^cX}^*/M_X[(bar3,2);Y=-1/6;B-L=-4/3]^2` |
| `L Q e^c u^c` | `X[(bar3,2);Y=-7/6;B-L=-4/3]` | B_and_L_conserving | not_a_proton_seed | `C6[QL u^c e^c; X[(bar3,2);Y=-7/6;B-L=-4/3]] = g_{Qe^cX} g_{Lu^cX}^*/M_X[(bar3,2);Y=-7/6;B-L=-4/3]^2` |
| `L L e^c nu^c` | `X[(1,1);Y=1;B-L=2]` | B_and_L_conserving | not_a_proton_seed | `C6[LL e^c nu^c; X[(1,1);Y=1;B-L=2]] = g_{LLX} g_{nu^ce^cX}^*/M_X[(1,1);Y=1;B-L=2]^2` |
| `L L e^c nu^c` | `X[(1,2);Y=1/2;B-L=0]` | B_and_L_conserving | not_a_proton_seed | `C6[LL e^c nu^c; X[(1,2);Y=1/2;B-L=0]] = g_{Lnu^cX} g_{Le^cX}^*/M_X[(1,2);Y=1/2;B-L=0]^2` |
| `d^c e^c u^c u^c` | `X[(bar3,1);Y=4/3;B-L=2/3]` | standard_BNV_seed | blocked_as_vector_until_completion | `C6[u^c u^c d^c e^c; X[(bar3,1);Y=4/3;B-L=2/3]] = g_{u^cu^cX} g_{d^ce^cX}^*/M_X[(bar3,1);Y=4/3;B-L=2/3]^2` |
| `d^c e^c u^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | standard_BNV_seed | blocked_as_vector_until_completion | `C6[u^c u^c d^c e^c; X[(bar3,1);Y=1/3;B-L=2/3]] = g_{u^cd^cX} g_{u^ce^cX}^*/M_X[(bar3,1);Y=1/3;B-L=2/3]^2` |
| `d^c d^c nu^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | BNV_with_sterile_neutrino | blocked_as_vector_until_completion | `C6[u^c d^c d^c nu^c; X[(bar3,1);Y=1/3;B-L=2/3]] = g_{u^cd^cX} g_{d^cnu^cX}^*/M_X[(bar3,1);Y=1/3;B-L=2/3]^2` |
| `d^c d^c nu^c u^c` | `X[(3,1);Y=2/3;B-L=-2/3]` | BNV_with_sterile_neutrino | blocked_as_vector_until_completion | `C6[u^c d^c d^c nu^c; X[(3,1);Y=2/3;B-L=-2/3]] = g_{u^cnu^cX} g_{d^cd^cX}^*/M_X[(3,1);Y=2/3;B-L=-2/3]^2` |

## Proton-Bound Interface

The physical width template is

```text
Gamma(p -> channel a) = sum_ij C_i^phys C_j^{phys *} H_ij^(a)
```

The corresponding bound template is

```text
tau_a^exp > 1/Gamma_a implies C^phys dagger H^(a) C^phys < 1/tau_a^exp
```

Required inputs before numerical proton bounds can be quoted:

- completed high-energy branch passing P6
- mediator masses and couplings
- flavor rotations to physical fermion basis
- operator basis and chiral contractions
- RG factors from matching scale to hadronic scale
- lattice/chiral matrix elements
- current experimental lower limit for the selected proton channel

## P7 Boundary

P7 is a conditional matching ledger.  It does not compute physical proton lifetimes because P6 completion, flavor rotations, mediator masses/couplings, RG factors, and hadronic matrix elements are not yet supplied.
