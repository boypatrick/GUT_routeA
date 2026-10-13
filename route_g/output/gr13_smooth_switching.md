# G-R13 one smooth pulse on the fixed 21-state mediator

70/70 checks pass. One physical pulse, no optimization.

Fix gamma=.03, Delta=80 and g_peak=.1. The relay envelope is f=sin²(pi*s/T);
the compensation is f²*Vc0, not f*Vc0. Peak coupling and effective preparation
area are unchanged. Duration increases from 6.25 to 16.66666667, a factor 8/3.
The continuous control parameter is assumed, not derived time or a rounded update count.

| quantity | compensated square | one smooth pulse |
|---|---:|---:|
| unconditional preparation fidelity | 0.9963111171 | 0.9981207221 |
| endpoint auxiliary probability | 0.003335852134 | 5.020910281e-13 |
| endpoint auxiliary energy | 0.2802115792 | 4.217564636e-11 |
| net controller work on full model | 0.2635368684 | -0.003314373769 |
| positive input, separate control ports | 0.2635855303 | 0.1712514185 |
| returned work, separate control ports | 4.866189573e-05 | 0.1745657923 |
| absolute work, separate control ports | 0.2636341922 | 0.3458172108 |
| duration | 6.25 | 16.66666667 |

## Complete energy ledger

Integrated relay work = -0.0381709643229; compensation work = 0.0348565905537.
Their sum matches the endpoint energy change, with the first law also checked
throughout the pulse. Compensation is not free or silently cancelled from the ledger.
Negative net work means energy is returned from this prepared input; it does not
mean zero input work, free initial-state preparation, zero dissipation or a free controller.
The sampled transient auxiliary probability reaches 0.00184626,
so the mediator remains active even though the endpoint residue is suppressed.

## What improved and what did not

For this single fixed case, fidelity, endpoint residue and ideal positive input work
all improve. Duration and absolute exchanged work increase. Actual hardware energy
consumption has not been supplied. If both protocols share return efficiency eta and
constant overhead power p0, define cost=W_plus-eta*W_minus+p0*duration using separate ports.
The smooth pulse costs less only if p0 < 0.00886407473 + 0.0167536445*eta
in model power units. No eta, p0 or physical scale is fitted or inferred.

## Verification and stopping boundary

Adaptive DOP853 tolerance refinement, independent unitary CF4 step refinement,
and power quadrature refinement test the same waveform. The tiny endpoint residue
is nonzero and numerically resolved, not a universal scaling law or a noise tolerance.
The unchanged G-R11 reader keeps every outcome and its 57/82 versus 25/82 disturbance.
This one-pulse question is bounded-done. Do not begin pulse optimization. A real
controller's losses/overhead or an explicit autonomous trigger are separate physical
choices; no time-arrow, Ed or Route-F claim is promoted.
