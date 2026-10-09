# G-R9b paired channels and two ways to one outlet

586/586 checks pass. All inputs are dimensionless toy assumptions; no fit.

## Result

G-R9's extra outlet is not an accident of theta=pi/2. Flipping the sign of every
link maps L to 4-L; on a three-node ring this shifts the loop phase by pi, so
spec L(pi-theta) = 4 - spec L(theta). The lower gap of one 3-ring equals the upper
gap of its mirror partner and vice versa. Hence for ANY two flux 3-rings, resonance
of |1,1> with |0,2> forces theta_B = pi - theta_A and makes |2,0> degenerate too.
Checked on 49 theta_A values in (0,pi); partner found by bisection.

Every real density-density contact couples |1,1> to the two outlets with equal
magnitude, and the resonant 3x3 block is circulant with loop phase Phi=3 arg z.
Symmetric contacts give Phi=0.

## Exact bound for the G-R9 vertex contact

The vertex contact is rank one. Two combinations of the degenerate triple are exact
dark eigenstates, so P(0,2) = |G_bb(s) - exp(-i E0 s)|^2 / 9 with |G_bb| <= 1.
Therefore P(0,2) <= 4/9 EXACTLY, for every gamma and every update, for identical
rings and for every mirror pair. This sharpens G-R9's exact 1/2 and promotes its
weak-contact 4/9 envelope to an exact bound. General rank-one form: P <= 4 w_i w_t.
Both outlets have IDENTICAL amplitudes at every update, also for mirror pairs,
which have no swap symmetry: the vertex contact cannot prefer either outlet.

| theta_A | gamma | max sampled P(0,2) | max sampled P(2,0) |
|---:|---:|---:|---:|
| 1.570796 | 0.09 | 0.444418019 | 0.444418019 |
| 1.570796 | 0.03 | 0.444431999 | 0.444431999 |
| 1.570796 | 0.01 | 0.444444088 | 0.444444088 |
| 1.570796 | 0.5 | 0.439608352 | 0.439608352 |
| 1.570796 | 2.0 | 0.437157653 | 0.437157653 |
| 1.570796 | 10.0 | 0.417446514 | 0.417446514 |
| 1.047198 | 0.09 | 0.444298572 | 0.444298572 |
| 1.047198 | 0.03 | 0.444425910 | 0.444425910 |
| 1.047198 | 0.01 | 0.444443243 | 0.444443243 |
| 1.047198 | 0.5 | 0.439938292 | 0.439938292 |
| 1.047198 | 2.0 | 0.386793114 | 0.386793114 |
| 1.047198 | 10.0 | 0.394901401 | 0.394901401 |
| 0.600000 | 0.09 | 0.444101439 | 0.444101439 |
| 0.600000 | 0.03 | 0.444409724 | 0.444409724 |
| 0.600000 | 0.01 | 0.444439018 | 0.444439018 |
| 0.600000 | 0.5 | 0.430555584 | 0.430555584 |
| 0.600000 | 2.0 | 0.391353447 | 0.391353447 |
| 0.600000 | 10.0 | 0.412216115 | 0.412216115 |

## Exit A: distinguishable spectra (3-ring + 4-ring)

A: 3-ring at theta_A=pi/2 (gaps sqrt3, sqrt3). B: 4-ring with
theta_B = 4 arcsin(sqrt3/4) = 1.791329587716; its spectrum is
{2 -+ sqrt13/2, 2 -+ sqrt3/2}, so B's middle gap is sqrt3.
Only |1,1> and |0,2> share the initial energy; the reverse outlet is detuned by
(3 sqrt3 - sqrt13)/2 = 0.795300574. Over the whole resonant domain
theta_A in (0.275643, pi) the reverse detuning stays >= 0.3185.
With the SAME G-R9 vertex contact, P(0,2) = |G_bb - exp(-i E0 s)|^2 / 4 can reach one.

| gamma | weak-peak update | best n | P(0,2) | P(2,0) | (1-P)/gamma^2 |
|---:|---:|---:|---:|---:|---:|
| 0.09 | 837.8 | 856 | 0.999543095 | 9.7e-05 | 0.0564 |
| 0.03 | 2513.3 | 2535 | 0.999936858 | 4.33e-07 | 0.0702 |
| 0.01 | 7539.8 | 7558 | 0.999994927 | 3.08e-07 | 0.0507 |

## Exit B: oriented contact on identical rings

Contact gamma sum h((j-l) mod 3) P_j x P_l with h=(2,0,1) keeps both copies
identical but orients the link. It commutes with the joint translation, so the
resonant triple is an exactly closed sector (no leakage at any gamma). Its loop
phase is -pi/2: a chiral triangle, which transfers |1,1> -> |0,2> with unit
probability at t* = 2 pi/(3 gamma). Reversing the orientation selects |2,0>.
The symmetric profile h=(1,0,0) has loop phase 0 and stays at 4/9.

| gamma | t* | update n | P(0,2) at n | P(2,0) at n |
|---:|---:|---:|---:|---:|
| 0.09 | 23.2711 | 93 | 0.9999976057 | 1.2e-06 |
| 0.03 | 69.8132 | 279 | 0.9999976057 | 1.2e-06 |
| 0.01 | 209.4395 | 838 | 0.9999997561 | 1.22e-07 |

Integer-update values differ from one only through the sampling offset of t*.

## What changed physically

One-for-one exchange needs the two copies to be distinguishable: by their spectra
(Exit A) or by an orientation of the contact (Exit B). Identical 3-rings, or a
mirror pair, joined by an unoriented local contact split the quantum between two
outlets: exactly for identical copies (swap, <=1/2) and for the vertex contact
(<=4/9), and in the weak-contact block otherwise.
Both exits still contain one supplied choice (a resonant flux; an oriented
contact ratio). Neither is derived from a dynamical principle.

## Next bounded decision

G-R10 built the all-channel readout for the vertex contact and found a conserved
2/3 dark weight of |1,1>. Lemma 3 is the same structure: |1,1> has bright weight 1/3.
Exit B does not lift G-R10's preparation obstruction: every translation-covariant
contact conserves the |T> weight 1/3 of |1,1>, and D+ is orthogonal to |T>, so
the fidelity to D+ never exceeds 2/3 (checked for oriented, symmetric and random
profiles). Exits A and B are declared alternatives for one-outlet EXCHANGE only.
G-R10's own next choice (a declared preparation/record instrument, or the
unfiltered stationary mixture) is unchanged by this side branch.
No physical seconds, records, metric or particle spectrum are claimed.
