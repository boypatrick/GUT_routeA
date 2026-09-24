# G-R1: a finite relational clock and an energy-density-only counterexample

Verification: **109/109 checks passed**.

The relational-time line is now Route G's priority. Older KK/TOF files are retained unchanged as comparison models. No dimensional energy, measured particle mass, external elapsed time or spatial geometry is fitted.

## Exact finite construction

For n=0,...,d-1, HC|n>=epsilon(d-1-n)|n>, HS|n>=epsilon n|n>. Both energies are nonnegative. The state |Psi>=sum a_n |n,n> satisfies (HC+HS-(d-1)epsilon)|Psi>=0.
Clock |phi>=sum exp(i n phi)|n>/sqrt(d); its d grid readings are orthonormal, while the continuous phase dial is a normalized POVM, not uncountably many orthogonal states.
Conditioning gives |psi(phi)>=sum a_n exp(-i n phi)|n>, so i epsilon d_phi |psi>=HS|psi>. Only after calibration tau=hbar phi/epsilon does this become the usual Schrödinger equation. No geometric rotation is introduced.

## Three-level example and exact counterexample

Uniform amplitudes, epsilon=1 in model units. Compare the coherent |Psi><Psi| with sum |n,n><n,n|/3. Both have identical local marginals, identical joint energy statistics, total energy2 and variance0. Measure |+>=(|0>+|1>+|2>)/sqrt(3) conditionally:

| clock phase / pi | coherent P(+|phi) | dephased P(+|phi) |
|---:|---:|---:|
| 0 | 1 | 0.333333333 |
| 0.3333333 | 0.444444444 | 0.333333333 |
| 0.6666667 | -4.62592927e-18 | 0.333333333 |
| 1 | 0.111111111 | 0.333333333 |
| 0.2269549 | 0.701577071 | 0.333333333 |

Exact laws: P_coherent=(1+2 cos(phi))²/9; P_dephased=1/3. At phi=0 the probabilities differ by2/3 despite identical energy data. The missing information is joint coherence, not a decimal correction to the energy spectrum.
Even after supplying a fixed volume V0, mean energy/V0 would be identical. This rules out that density-only identification, not a richer Ei/Ed that explicitly includes relational coherences and has its own dynamics.

## Guardrails that passed

- Constraint, positivity, clock-grid completeness, continuous POVM and history-state identities.
- Direct projection versus unitary conditional state, finite-difference conditional equation, arbitrary normalized nonuniform amplitudes.
- Constraint-invariant joint effects obtained by exact Fourier twirling preserve the same conditional probabilities.
- Finite clock commutator has zero trace; no impossible exact canonical time operator is assumed.
- Common energy scaling changes calibration but not phase-conditioned correlations; a local energy-zero shift with matching total shift leaves the constraint unchanged.
- A naive clock/system product generally fails the constraint. An arbitrary clock interaction does not preserve the old stationary state.
- All comparison code and selected earlier result hashes are unchanged.

## Not closed / next step

The model assumes quantum mechanics, the tensor split, Hamiltonians, a state and a clock POVM. It is an exact conditional construction, not an unconditional derivation of time from nothing. The clock is cyclic and single-event conditioning is not a sequential measurement history.
G-R2 should add a minimal explicit record/measurement model before interpreting the original T as accumulated event count. G-R3 can then define an Ei/Ed candidate and one shared clock-matter constraint, and compare two physical clocks. Do not insert a fitted slowdown function or relabel the phase as an extra spatial coordinate.
Full derivation and dictionary: tex/route_g_relational_clock.tex and RELATIONAL_TIME_CARD.md. The earlier Higgs/xenon rejection is neither undone nor a bound on this unspecified new realization.
