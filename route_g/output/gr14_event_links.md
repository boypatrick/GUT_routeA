# G-R14 Event links from one binary register rule

223/223 checks pass. Fixed fixtures; no physical fit.

U(Q)=exp(-i a Q tensor |1><1|), a=2 pi/3; fresh registers have active probability p.
Every response comes from F_p(x)=1-p+p exp(-i x); shared events use the SUM of probe charges.

| p | mean progress N | phase rate A / bare A | phase rate B / bare B | loop phase (rad) | loop visibility |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.00 | -0.000000000 | -0.000000000 | 0.000000000 | 1.000000000 |
| 0.25 | 0.25 | 0.245641993 | 0.231631437 | -1.000419517 | 0.289379050 |
| 0.50 | 0.50 | 0.500000000 | 0.500000000 | -3.141592654 | 0.125000000 |
| 0.75 | 0.75 | 0.754358007 | 0.768368563 | 1.000419517 | 0.289379050 |
| 1.00 | 1.00 | 1.000000000 | 1.000000000 | 0.000000000 | 1.000000000 |

## Exact constraints and the universality boundary

With z=F_p(2 pi/3), r=|z| and W=z^3=r^3 exp(i Phi),
r^2=1-3p(1-p) and 1-3r^2+2r^3 cos(Phi)=0.
Clock coherence is F_p(a q)^n, not an independently chosen damping function.
A common mean progress p does not generally give a common fractional phase rate.
For all sufficiently small q it is exact only at p=0, 1/2 or 1; at p=1/2 it holds
on |a q|<pi, with visibility |cos(a q/2)|^n. Phase is undefined at zeros.

## The memory counterprotocol and joint prediction

Reusing one record on three links gives F_p(3a)=1, instead of F_p(a)^3.
At p=1/2, fresh records give phase pi and visibility 1/8; the reused record gives phase 0 and visibility 1.
The marginals of one event are identical. The inter-event memory rule changes, not a fitted phase.
A simultaneous clock/loop coherence is F_p(a(q+1))^3, not F_p(a q)^3 F_p(a)^3.
Equal-frequency clocks sharing the same records retain unit relative coherence even when their individual dials fade.

## Interpretation and stopping boundary

This is a bounded quantum collision model inspired by event comparison, not a proof of the source theory.
The three links are not spatial dimensions. The reduced loop phase is operational interference, not a Maxwell field.
Update ordering, register supply, universal coupling, probe charges, quantum probabilities and readout remain inputs.
Tracing records is reversible loss of access; no irreversible decoherence or autonomous event generation is proved.
For the commuting observables tested, classical correlated binary noise reproduces the reduced predictions.
Next constrain memory correlations or require cross-probe agreement without independent phase/noise fits.
G-R13 controller calibration remains an implementation side branch, not a prerequisite.
