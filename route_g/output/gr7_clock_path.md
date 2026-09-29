# G-R7: frequency, phase and clock-path correlations

Status: conservative kernel bounded-done; full apparatus and empirical gates OPEN.
66/66 checks pass. No fit or new time law.

## One evolution, three linked observables

Same G-R4/G-R6 gaps: theta_ir=2*pi*nu_i*integral[(1+b_i(T_r))*d_tau_r].
The phase derivative gives frequency. The overlap
kappa=Tr(rho_clock U_B^dagger U_A) fixes the internal contribution
to fringe phase and ideal visibility. It is not an independently fitted
cosine. The full port phase also contains scalar motional/ground/guide
and laser contributions, which are not supplied by differential polarizabilities.

For (g+e_i)/sqrt(2), kappa=exp(-i*delta_i/2)*cos(delta_i/2).
The ideal visibility is |cos(delta_i/2)|; the phase is undefined at zeros.
Pure joint-state negativity is |sin(delta_i/2)|/2.
For one shared-ground qutrit use (1+exp(-i*delta_2)+exp(-i*delta_3))/3,
not the product of two independent qubit overlaps.

## Conditional hold calculation, not an apparatus forecast

Hold 0.01 s, illustrative 1 m separation and 300/310 K baths.
Other trajectory segments and control phases are not invented.
Numbers below are for the ideal separate two-level preparations.

| Hold case | Clock | Delta frequency (Hz) | Delta clock phase (rad) | 1 - ideal visibility |
|---|---|---:|---:|---:|
| gravity_only | E2 | +0.075109 | +0.0047193 | 2.784e-06 |
| gravity_only | E3 | +0.070064 | +0.0044023 | 2.422e-06 |
| thermal_only | E2 | -0.050506 | -0.0031734 | 1.259e-06 |
| thermal_only | E3 | -0.0064999 | -0.0004084 | 2.085e-08 |
| gravity_and_thermal | E2 | +0.024604 | +0.0015459 | 2.987e-07 |
| gravity_and_thermal | E3 | +0.063564 | +0.0039939 | 1.994e-06 |

The static-response phase combination
sum_i w_i*delta_i/(2*pi*nu_i) = Delta tau,
with the unchanged G-R5 weights, rejects the same leading thermal
integral even with path-dependent proper-time weighting. This requires
matched trajectories/readout conventions and calibrated relative phases.
It is not a linear combination of visibilities, nor evidence for a
free-standing time law. Internal scalar energy-zero conventions cannot
change the total predicted fringe when the scalar action is transformed too.

## Visibility does not certify entanglement

For an algebra fixture Ub=Z, Ua=I, the pure input |+> gives a
maximally entangled joint state with negativity 1/2. The energy-dephased
input I/2 gives a separable joint state with negativity zero. Both have
the same completely incoherent reduced path state and zero path visibility.
The joint witness (I-X_path*Z_clock-Z_path*X_clock)/2 is -1/2 for the
pure fixture and zero for the mixed fixture; the separable bound is zero.
This fixture is not a realized pi-phase Yb+ experiment.

A common final clock unitary, or any unconditioned trace-preserving local
clock channel, cannot change the reduced path state. Clock decoherence
is therefore not automatically an extra multiplicative path-visibility
factor. A path-controlled inverse can undo the ideal correlation, but
requires specified controls; postselection must retain its success rates.

## Physical limits and next useful action

At 1 m the ideal gravity-only first visibility zeros need about
6.66 s (E2) and 7.14 s (E3).
Published E2 excited-state lifetimes are on the tens-of-ms scale.
Do not extrapolate that unperturbed E2 hold to seconds. The lifetime
references in CLOCK_PATH_CARD.json are cautions, not a pooled decay rate
or complete environmental model. E3 longevity alone does not certify
charged-ion path coherence, guide closure, or laser stability.

The conservative real BBR shift does not determine dissipative field
correlations, photon which-path information or spontaneous decay.
Only under additional factorization assumptions may a separate motional/
environment overlap multiply kappa. No physical open-system visibility is
assigned. A conventional freely falling interferometer can also cancel the
putative uniform-field proper-time contrast; derive complete paths and
pulse phases before substituting g*height*time.

Next prioritize an explicit closed-path/guide/readout sequence and
state-resolved differential phase (linear in delta), with coherent vs
dephased controls and joint correlations when claiming entanglement.
Visibility loss is quadratic for small delta and is a consistency
observable, not the sole or first discovery channel.
No quantum-gravity, new Ed law, experiment, or Route-F promotion is claimed.
Derivation: tex/route_g_clock_path.tex; sources: SOURCES.md.
