# Route-C P2 Symbolic Four-Point Pole Ansatz

This is the first charge-filtered pole ledger.  It is not yet a
unitarity, positivity, or Ward-identity proof.  A pole listed here is
only a symbolic term not forbidden by the P1 SM-face charge table.

## Ansatz

For all-incoming left-handed Weyl bilinears, the seed ansatz is

```text
A_4(12|34) = sum_X R_s(X)/(s-M_X^2) + C6(12|34),
R_s(X) = g_{12X} g_{34Xbar}.
```

The script declares the quantum numbers of `X` from the pair product
and keeps every residue symbolic.

## Summary

| quantity | value |
| --- | ---: |
| pair types | 21 |
| pair-product components | 31 |
| tested pair-component pairings | 496 |
| allowed symbolic pole components | 18 |
| forbidden by Abelian charge | 455 |
| forbidden by non-Abelian conjugacy after Abelian pass | 23 |

## Representative Allowed Poles

| left pair | right pair | mediator | external multiset |
| --- | --- | --- | --- |
| `Q+Q` in `(6,1)` | `u^c+d^c` in `(bar6,1)` | `X[(bar6,1);Y=-1/3;B-L=-2/3]` | `Q Q d^c u^c` |
| `Q+Q` in `(bar3,1)` | `Q+L` in `(3,1)` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `L Q Q Q` |
| `Q+Q` in `(bar3,1)` | `u^c+d^c` in `(3,1)` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `Q Q d^c u^c` |
| `Q+Q` in `(bar3,3)` | `Q+L` in `(3,3)` | `X[(3,3);Y=-1/3;B-L=-2/3]` | `L Q Q Q` |
| `Q+L` in `(3,1)` | `u^c+e^c` in `(bar3,1)` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `L Q e^c u^c` |
| `Q+L` in `(3,1)` | `d^c+nu^c` in `(bar3,1)` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `L Q d^c nu^c` |
| `Q+u^c` in `(1,2)` | `Q+d^c` in `(1,2)` | `X[(1,2);Y=1/2;B-L=0]` | `Q Q d^c u^c` |
| `Q+u^c` in `(1,2)` | `L+e^c` in `(1,2)` | `X[(1,2);Y=1/2;B-L=0]` | `L Q e^c u^c` |
| `Q+u^c` in `(8,2)` | `Q+d^c` in `(8,2)` | `X[(8,2);Y=1/2;B-L=0]` | `Q Q d^c u^c` |
| `Q+d^c` in `(1,2)` | `L+nu^c` in `(1,2)` | `X[(1,2);Y=-1/2;B-L=0]` | `L Q d^c nu^c` |
| `Q+nu^c` in `(3,2)` | `L+d^c` in `(bar3,2)` | `X[(bar3,2);Y=-1/6;B-L=-4/3]` | `L Q d^c nu^c` |
| `Q+e^c` in `(3,2)` | `L+u^c` in `(bar3,2)` | `X[(bar3,2);Y=-7/6;B-L=-4/3]` | `L Q e^c u^c` |
| `L+L` in `(1,1)` | `nu^c+e^c` in `(1,1)` | `X[(1,1);Y=1;B-L=2]` | `L L e^c nu^c` |
| `L+nu^c` in `(1,2)` | `L+e^c` in `(1,2)` | `X[(1,2);Y=1/2;B-L=0]` | `L L e^c nu^c` |
| `u^c+u^c` in `(3,1)` | `d^c+e^c` in `(bar3,1)` | `X[(bar3,1);Y=4/3;B-L=2/3]` | `d^c e^c u^c u^c` |
| `u^c+d^c` in `(3,1)` | `u^c+e^c` in `(bar3,1)` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `d^c e^c u^c u^c` |
| `u^c+d^c` in `(3,1)` | `d^c+nu^c` in `(bar3,1)` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `d^c d^c nu^c u^c` |
| `u^c+nu^c` in `(bar3,1)` | `d^c+d^c` in `(3,1)` | `X[(3,1);Y=2/3;B-L=-2/3]` | `d^c d^c nu^c u^c` |

## Proton-Operator Seeds

These entries are useful later for low-energy matching.  Their
presence here means only that a charge-allowed symbolic pole exists;
P3--P7 must still check residues, Ward identities, high-energy
behavior, and bounds.

| external multiset | allowed pole components |
| --- | ---: |
| `Q Q Q L` | 2 |
| `u^c u^c d^c e^c` | 2 |
| `Q Q u^c e^c` | 0 |
| `Q L u^c d^c` | 0 |

## Current-Transition Charge Ledger

The current ledger records additive charges for schematic
`psi_a^dagger psi_b` transitions.  It is included for later
vector-exchange Ward checks, but it does not yet provide explicit
GUT generator matrices.

| transition | Y current charge | B-L current charge | classification |
| --- | ---: | ---: | --- |
| `(Q)^dagger -> Q` | 0 | 0 | diagonal SM-face or Cartan current |
| `(Q)^dagger -> L` | -2/3 | -4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(Q)^dagger -> u^c` | -5/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(Q)^dagger -> d^c` | 1/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(Q)^dagger -> nu^c` | -1/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(Q)^dagger -> e^c` | 5/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(L)^dagger -> Q` | 2/3 | 4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(L)^dagger -> L` | 0 | 0 | diagonal SM-face or Cartan current |
| `(L)^dagger -> u^c` | -1/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(L)^dagger -> d^c` | 5/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(L)^dagger -> nu^c` | 1/2 | 2 | off-diagonal candidate current; requires declared GUT generator |
| `(L)^dagger -> e^c` | 3/2 | 2 | off-diagonal candidate current; requires declared GUT generator |
| `(u^c)^dagger -> Q` | 5/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(u^c)^dagger -> L` | 1/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(u^c)^dagger -> u^c` | 0 | 0 | diagonal SM-face or Cartan current |
| `(u^c)^dagger -> d^c` | 1 | 0 | off-diagonal candidate current; requires declared GUT generator |
| `(u^c)^dagger -> nu^c` | 2/3 | 4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(u^c)^dagger -> e^c` | 5/3 | 4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(d^c)^dagger -> Q` | -1/6 | 2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(d^c)^dagger -> L` | -5/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(d^c)^dagger -> u^c` | -1 | 0 | off-diagonal candidate current; requires declared GUT generator |
| `(d^c)^dagger -> d^c` | 0 | 0 | diagonal SM-face or Cartan current |
| `(d^c)^dagger -> nu^c` | -1/3 | 4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(d^c)^dagger -> e^c` | 2/3 | 4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(nu^c)^dagger -> Q` | 1/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(nu^c)^dagger -> L` | -1/2 | -2 | off-diagonal candidate current; requires declared GUT generator |
| `(nu^c)^dagger -> u^c` | -2/3 | -4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(nu^c)^dagger -> d^c` | 1/3 | -4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(nu^c)^dagger -> nu^c` | 0 | 0 | diagonal SM-face or Cartan current |
| `(nu^c)^dagger -> e^c` | 1 | 0 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> Q` | -5/6 | -2/3 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> L` | -3/2 | -2 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> u^c` | -5/3 | -4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> d^c` | -2/3 | -4/3 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> nu^c` | -1 | 0 | off-diagonal candidate current; requires declared GUT generator |
| `(e^c)^dagger -> e^c` | 0 | 0 | diagonal SM-face or Cartan current |

## P2 Boundary

Allowed here means allowed by SM-face charge conservation and
non-Abelian conjugacy of the pair products.  It does not mean the
mediator exists in a candidate action, that the residue has positive
norm, or that Ward identities cancel.  Those are P3--P6 tasks.
