# G-R3: density candidate, interaction and relative clocks

Status: bounded-done; 274/274 checks.
This is a declared model, not a derivation of Chen's undefined Ed.

## Definition and observable

Ed = epsilon_D N_D / nu(D), with nu(D)=1 explicitly addressed matter cell.
Use ground-subtracted excitation energy, not assumed spatial volume,
record entropy, a timestamp, or the entire information in the joint state.
The one-cell test cannot distinguish density from total excitation energy;
a multi-cell claim must fix the averaging measure and coupling profile.
E0 is epsilon_D per cell. V=kappa H_B Ed/E0 is one fixed dispersive
interaction. Reference A is shielded from D; its calibration stays fixed.

A single stationary 15x3x2 model gives rate r_B/A=1+kappa d for d=0,1.
The same reference spectrum and readout serve kappa=-1/4,0,+1/4:
excited-source rates are 3/4,1,5/4. These are synthetic checks, not fitted
physics. Common energy rescaling leaves the ratio invariant; equal
fractional coupling of both clocks produces no relative-rate effect.

## Exact finite-clock predictions

At reference phase pi/2, for the +1/4 card, B's three orthogonal dial
probabilities are [1,0,0] with d=0 and
[1/9,(1+sqrt(3))^2/9,(1-sqrt(3))^2/9] with d=1.
For an unobserved equal source mixture or coherent superposition, P0=5/9.
Substituting mean Ed into the Hamiltonian incorrectly predicts
(3+2sqrt(2))/9=0.647603014.
Density fluctuations produce multiple conditional rates and reduced
clock coherence, not necessarily one slower clock. At reference phase
pi, adjacent coherence vanishes but B is not fully dephased (purity5/9).
This does not by itself distinguish quantum source coherence from a
classical mixture: B alone sees source populations.

## Microscopic sign audit

Jaynes-Cummings exchange, in its controlled dispersive regime, yields
chi=g^2/Delta and calibrated r=(omega+chi)/(omega-chi).
Detuning +0.5 gives r=1.003606492;
detuning -0.5 gives r=0.996406468.
Exact 2x2 spectra verify both signs and the second-order remainder bound.
Higher-order clock gaps are not equal: a universal rigid time dilation
does not follow from this mechanism. These are model units, not a
transmon calibration or an exact realization of the illustrative kappa=1/4.

## Closure boundary and next decision

The candidate/interaction/two-clock algebra is closed within its assumptions.
Finite phase likelihoods require repeated calibrated comparisons and
event/epoch labels, not an exact time from one outcome.
No external time is used in the stationary conditional construction.
The source/reference state, commensurate spectrum and reference
shielding are declared inputs. The G-R2 history constraint and the new
stationary constraint have NOT been conflated into one energy-complete
recording apparatus.

Next choose a physically grounded source/clock realization and determine
coupling sign, or test universality across independent clock transitions.
Do not fit a density slowdown law or promote gravity/Route-F conclusions.

Run: python3 route_g/code/verify_gr3_density_clocks.py.
Every check, fixed input and prior-file hash is stored in the JSON.
