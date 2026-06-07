# Route-C P3 Residue Factorization and Positivity Ledger

This ledger turns the P2 charge-allowed symbolic poles into
two-particle channel blocks.  It verifies the conditional Gram-matrix
positivity statement for positive-norm mediator exchange.

## Mathematical Statement

For a declared mediator sector `X`, let `alpha` label all two-particle
channels that can couple to `X`.  Near the pole, unitarity requires

```text
Res_X A_{alpha beta} = g_{alpha X} g^*_{beta X}.
```

This is a Gram matrix.  With a positive-norm mediator metric it is
positive semidefinite.  With a wrong-sign metric, the same unit-coupling
test has a negative eigenvalue and fails.

## Summary

| quantity | value |
| --- | ---: |
| mediator sectors | 26 |
| total vertex channels | 31 |
| P2 allowed poles checked | 18 |
| factorization checks passed | 18 |
| PSD blocks under positive metric | 26 |
| largest channel-block size | 3 |

## Mediator Blocks

| mediator | channel count | unit-coupling positive spectrum | wrong-sign stress spectrum |
| --- | ---: | --- | --- |
| `X[(1,1);Y=-1;B-L=-2]` | 1 | `[1]` | `[-1]` |
| `X[(1,1);Y=-2;B-L=-2]` | 1 | `[1]` | `[-1]` |
| `X[(1,1);Y=0;B-L=-2]` | 1 | `[1]` | `[-1]` |
| `X[(1,1);Y=1;B-L=2]` | 1 | `[1]` | `[-1]` |
| `X[(1,2);Y=-1/2;B-L=0]` | 2 | `[2, 0]` | `[-2, 0]` |
| `X[(1,2);Y=1/2;B-L=0]` | 2 | `[2, 0]` | `[-2, 0]` |
| `X[(1,3);Y=1;B-L=2]` | 1 | `[1]` | `[-1]` |
| `X[(3,1);Y=-1/3;B-L=-2/3]` | 3 | `[3, 0, 0]` | `[-3, 0, 0]` |
| `X[(3,1);Y=-4/3;B-L=-2/3]` | 1 | `[1]` | `[-1]` |
| `X[(3,1);Y=2/3;B-L=-2/3]` | 1 | `[1]` | `[-1]` |
| `X[(3,2);Y=1/6;B-L=4/3]` | 1 | `[1]` | `[-1]` |
| `X[(3,2);Y=7/6;B-L=4/3]` | 1 | `[1]` | `[-1]` |
| `X[(3,3);Y=-1/3;B-L=-2/3]` | 1 | `[1]` | `[-1]` |
| `X[(6,1);Y=-2/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(6,1);Y=1/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(6,1);Y=4/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(8,2);Y=-1/2;B-L=0]` | 1 | `[1]` | `[-1]` |
| `X[(8,2);Y=1/2;B-L=0]` | 1 | `[1]` | `[-1]` |
| `X[(bar3,1);Y=-2/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar3,1);Y=1/3;B-L=2/3]` | 2 | `[2, 0]` | `[-2, 0]` |
| `X[(bar3,1);Y=4/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar3,2);Y=-1/6;B-L=-4/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar3,2);Y=-7/6;B-L=-4/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar3,3);Y=1/3;B-L=2/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar6,1);Y=-1/3;B-L=-2/3]` | 1 | `[1]` | `[-1]` |
| `X[(bar6,3);Y=-1/3;B-L=-2/3]` | 1 | `[1]` | `[-1]` |

## Factorization Checks

| external multiset | mediator | left channel | right channel | passes |
| --- | --- | --- | --- | --- |
| `Q Q d^c u^c` | `X[(bar6,1);Y=-1/3;B-L=-2/3]` | `Q+Q_(6,1;Y=1/3;B-L=2/3)` | `u^c+d^c_(bar6,1;Y=-1/3;B-L=-2/3)` | True |
| `L Q Q Q` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `Q+Q_(bar3,1;Y=1/3;B-L=2/3)` | `Q+L_(3,1;Y=-1/3;B-L=-2/3)` | True |
| `Q Q d^c u^c` | `X[(3,1);Y=-1/3;B-L=-2/3]` | `Q+Q_(bar3,1;Y=1/3;B-L=2/3)` | `u^c+d^c_(3,1;Y=-1/3;B-L=-2/3)` | True |
| `L Q Q Q` | `X[(3,3);Y=-1/3;B-L=-2/3]` | `Q+Q_(bar3,3;Y=1/3;B-L=2/3)` | `Q+L_(3,3;Y=-1/3;B-L=-2/3)` | True |
| `L Q e^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `Q+L_(3,1;Y=-1/3;B-L=-2/3)` | `u^c+e^c_(bar3,1;Y=1/3;B-L=2/3)` | True |
| `L Q d^c nu^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `Q+L_(3,1;Y=-1/3;B-L=-2/3)` | `d^c+nu^c_(bar3,1;Y=1/3;B-L=2/3)` | True |
| `Q Q d^c u^c` | `X[(1,2);Y=1/2;B-L=0]` | `Q+u^c_(1,2;Y=-1/2;B-L=0)` | `Q+d^c_(1,2;Y=1/2;B-L=0)` | True |
| `L Q e^c u^c` | `X[(1,2);Y=1/2;B-L=0]` | `Q+u^c_(1,2;Y=-1/2;B-L=0)` | `L+e^c_(1,2;Y=1/2;B-L=0)` | True |
| `Q Q d^c u^c` | `X[(8,2);Y=1/2;B-L=0]` | `Q+u^c_(8,2;Y=-1/2;B-L=0)` | `Q+d^c_(8,2;Y=1/2;B-L=0)` | True |
| `L Q d^c nu^c` | `X[(1,2);Y=-1/2;B-L=0]` | `Q+d^c_(1,2;Y=1/2;B-L=0)` | `L+nu^c_(1,2;Y=-1/2;B-L=0)` | True |
| `L Q d^c nu^c` | `X[(bar3,2);Y=-1/6;B-L=-4/3]` | `Q+nu^c_(3,2;Y=1/6;B-L=4/3)` | `L+d^c_(bar3,2;Y=-1/6;B-L=-4/3)` | True |
| `L Q e^c u^c` | `X[(bar3,2);Y=-7/6;B-L=-4/3]` | `Q+e^c_(3,2;Y=7/6;B-L=4/3)` | `L+u^c_(bar3,2;Y=-7/6;B-L=-4/3)` | True |
| `L L e^c nu^c` | `X[(1,1);Y=1;B-L=2]` | `L+L_(1,1;Y=-1;B-L=-2)` | `nu^c+e^c_(1,1;Y=1;B-L=2)` | True |
| `L L e^c nu^c` | `X[(1,2);Y=1/2;B-L=0]` | `L+nu^c_(1,2;Y=-1/2;B-L=0)` | `L+e^c_(1,2;Y=1/2;B-L=0)` | True |
| `d^c e^c u^c u^c` | `X[(bar3,1);Y=4/3;B-L=2/3]` | `u^c+u^c_(3,1;Y=-4/3;B-L=-2/3)` | `d^c+e^c_(bar3,1;Y=4/3;B-L=2/3)` | True |
| `d^c e^c u^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `u^c+d^c_(3,1;Y=-1/3;B-L=-2/3)` | `u^c+e^c_(bar3,1;Y=1/3;B-L=2/3)` | True |
| `d^c d^c nu^c u^c` | `X[(bar3,1);Y=1/3;B-L=2/3]` | `u^c+d^c_(3,1;Y=-1/3;B-L=-2/3)` | `d^c+nu^c_(bar3,1;Y=1/3;B-L=2/3)` | True |
| `d^c d^c nu^c u^c` | `X[(3,1);Y=2/3;B-L=-2/3]` | `u^c+nu^c_(bar3,1;Y=-2/3;B-L=2/3)` | `d^c+d^c_(3,1;Y=2/3;B-L=-2/3)` | True |

## P3 Boundary

The PSD statement is conditional.  It assumes that the mediator sector
exists in the candidate action with positive norm and canonical
unitarity normalization.  P3 does not yet check explicit generator
algebra, spin-statistics projections, Ward identities, high-energy
growth, or proton bounds.
