# Route-C P4 Ward-Identity Bookkeeping Ledger

This is a bookkeeping layer for Ward identities.  It does not yet
compute complete amplitudes.  It identifies which current sectors are
ready for unbroken SM Ward checks and which broken/off-face sectors need
Goldstone/Higgs completion before a real Ward test can be claimed.

## Ward Targets

For an unbroken gauge boson, the target check is

```text
epsilon_mu -> p_mu,    p_mu M^mu = 0.
```

For a broken massive vector, the target is instead the generalized
Ward/Slavnov-Taylor relation

```text
p_mu M^mu(V_X) = M_X M(phi_X),
```

so the Goldstone/Higgs/source sector must be supplied before Route C can
test the amplitude.

## Summary

| quantity | value |
| --- | ---: |
| unbroken SM current sectors | 3 |
| current transitions checked | 36 |
| massless diagonal transitions ready | 6 |
| broken/off-face transitions needing completion | 30 |
| mediator sectors inspected | 26 |
| mediator sectors needing completion if vector | 26 |

## Unbroken SM Current Readiness

| current sector | acts on | bookkeeping pass | full component Ward proved | status |
| --- | --- | --- | --- | --- |
| SU(3)_C current | `Q, u^c, d^c` | True | False | bookkeeping-ready: preserves multiplet label and acts inside color components; full check needs component-level SU(3) generators |
| SU(2)_L current | `Q, L` | True | False | bookkeeping-ready: preserves multiplet label and acts inside weak doublets; full check needs component-level SU(2) generators |
| U(1)_Y current | `Q, L, u^c, d^c, nu^c, e^c` | True | False | bookkeeping-ready: diagonal hypercharge current with exact P1 charge conservation |

## Current Transition Classification

| transition | Y carried | B-L carried | classification | completion needed |
| --- | ---: | ---: | --- | --- |
| `(Q)^dagger -> Q` | 0 | 0 | unbroken diagonal SM-face current | False |
| `(Q)^dagger -> L` | -2/3 | -4/3 | broken Pati-Salam leptoquark current | True |
| `(Q)^dagger -> u^c` | -5/6 | -2/3 | broken Spin(10) off-face current | True |
| `(Q)^dagger -> d^c` | 1/6 | -2/3 | broken Spin(10) off-face current | True |
| `(Q)^dagger -> nu^c` | -1/6 | 2/3 | broken Spin(10) off-face current | True |
| `(Q)^dagger -> e^c` | 5/6 | 2/3 | broken Spin(10) off-face current | True |
| `(L)^dagger -> Q` | 2/3 | 4/3 | broken Pati-Salam leptoquark current | True |
| `(L)^dagger -> L` | 0 | 0 | unbroken diagonal SM-face current | False |
| `(L)^dagger -> u^c` | -1/6 | 2/3 | broken Spin(10) off-face current | True |
| `(L)^dagger -> d^c` | 5/6 | 2/3 | broken Spin(10) off-face current | True |
| `(L)^dagger -> nu^c` | 1/2 | 2 | broken Spin(10) off-face current | True |
| `(L)^dagger -> e^c` | 3/2 | 2 | broken Spin(10) off-face current | True |
| `(u^c)^dagger -> Q` | 5/6 | 2/3 | broken Spin(10) off-face current | True |
| `(u^c)^dagger -> L` | 1/6 | -2/3 | broken Spin(10) off-face current | True |
| `(u^c)^dagger -> u^c` | 0 | 0 | unbroken diagonal SM-face current | False |
| `(u^c)^dagger -> d^c` | 1 | 0 | broken SU(2)_R charged current | True |
| `(u^c)^dagger -> nu^c` | 2/3 | 4/3 | broken Pati-Salam leptoquark current | True |
| `(u^c)^dagger -> e^c` | 5/3 | 4/3 | broken Spin(10) off-face current | True |
| `(d^c)^dagger -> Q` | -1/6 | 2/3 | broken Spin(10) off-face current | True |
| `(d^c)^dagger -> L` | -5/6 | -2/3 | broken Spin(10) off-face current | True |
| `(d^c)^dagger -> u^c` | -1 | 0 | broken SU(2)_R charged current | True |
| `(d^c)^dagger -> d^c` | 0 | 0 | unbroken diagonal SM-face current | False |
| `(d^c)^dagger -> nu^c` | -1/3 | 4/3 | broken Spin(10) off-face current | True |
| `(d^c)^dagger -> e^c` | 2/3 | 4/3 | broken Pati-Salam leptoquark current | True |
| `(nu^c)^dagger -> Q` | 1/6 | -2/3 | broken Spin(10) off-face current | True |
| `(nu^c)^dagger -> L` | -1/2 | -2 | broken Spin(10) off-face current | True |
| `(nu^c)^dagger -> u^c` | -2/3 | -4/3 | broken Pati-Salam leptoquark current | True |
| `(nu^c)^dagger -> d^c` | 1/3 | -4/3 | broken Spin(10) off-face current | True |
| `(nu^c)^dagger -> nu^c` | 0 | 0 | unbroken diagonal SM-face current | False |
| `(nu^c)^dagger -> e^c` | 1 | 0 | broken SU(2)_R charged current | True |
| `(e^c)^dagger -> Q` | -5/6 | -2/3 | broken Spin(10) off-face current | True |
| `(e^c)^dagger -> L` | -3/2 | -2 | broken Spin(10) off-face current | True |
| `(e^c)^dagger -> u^c` | -5/3 | -4/3 | broken Spin(10) off-face current | True |
| `(e^c)^dagger -> d^c` | -2/3 | -4/3 | broken Pati-Salam leptoquark current | True |
| `(e^c)^dagger -> nu^c` | -1 | 0 | broken SU(2)_R charged current | True |
| `(e^c)^dagger -> e^c` | 0 | 0 | unbroken diagonal SM-face current | False |

## P3 Mediator Ward Requirements

| mediator | channels | Ward target | completion needed if vector |
| --- | ---: | --- | --- |
| `X[(1,1);Y=-1;B-L=-2]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,1);Y=-2;B-L=-2]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,1);Y=0;B-L=-2]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,1);Y=1;B-L=2]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,2);Y=-1/2;B-L=0]` | 2 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,2);Y=1/2;B-L=0]` | 2 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(1,3);Y=1;B-L=2]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,1);Y=-1/3;B-L=-2/3]` | 3 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,1);Y=-4/3;B-L=-2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,1);Y=2/3;B-L=-2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,2);Y=1/6;B-L=4/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,2);Y=7/6;B-L=4/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(3,3);Y=-1/3;B-L=-2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(6,1);Y=-2/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(6,1);Y=1/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(6,1);Y=4/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(8,2);Y=-1/2;B-L=0]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(8,2);Y=1/2;B-L=0]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,1);Y=-2/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,1);Y=1/3;B-L=2/3]` | 2 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,1);Y=4/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,2);Y=-1/6;B-L=-4/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,2);Y=-7/6;B-L=-4/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar3,3);Y=1/3;B-L=2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar6,1);Y=-1/3;B-L=-2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |
| `X[(bar6,3);Y=-1/3;B-L=-2/3]` | 1 | `p_mu M^mu(V_X)=M_X M(phi_X)` | True |

## P4 Boundary

P4 is not a complete Ward-identity proof.  It verifies that the
unbroken SM current sectors are bookkeeping-ready and reports all
broken/off-face sectors that need explicit generator, mass,
Goldstone, Higgs, or source data.  The actual cancellations belong
to later P4 refinements and P6 high-energy-growth checks.
