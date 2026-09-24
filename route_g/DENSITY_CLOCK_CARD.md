# G-R3 density / interaction / clock convention card

Date: 2026-09-24. Status: bounded candidate test, not Chen's completed Ed law.

## Definitions fixed before calculating

| Symbol | Definition | What is not implied |
|---|---|---|
| D | One addressed two-level matter subsystem | A spatial voxel or G-R2's memory |
| nu(D)=1 | Counting measure: one specified operational cell | An emergent physical volume |
| N_D | Excitation projector, eigenvalues 0 and 1 | Entropy, elapsed time or total event history |
| epsilon_D | Energy gap above the same source ground state | A fitted particle mass |
| Ed operator | epsilon_D N_D / nu(D) | Complete Ei; universal energy density |
| E0 | epsilon_D / nu(D) | A freely retuned per-clock normalization |
| H_B | epsilon diag(0,1,2) | A physical frequency calibration in hertz |
| V | kappa H_B Ed/E0 | Gravity or a preselected slowdown law |
| Reference A | Fixed complementary spectrum, no D coupling | A universal ideal reference for arbitrary interactions |

A candidate full information specification contains the joint state,
event algebra, instruments and reference measure. Ed is only one
observable of that specification. It need not determine the history.
With one cell, total excitation energy and density cannot be empirically
distinguished. A many-cell proposal must justify its addressed set and
coupling profile; changing the averaging measure is not a physical effect.

## Fixed exact card

- Model energy unit epsilon=1; epsilon_D=1; delta=1/4.
- A has 15 levels, energies (14-j)delta; E*=14delta.
- B has three levels; D has two.
- Main kappa=1/4; audits kappa=0 and -1/4 use the same reference and POVMs.
- K=H_A+H_B+H_D+V-E*. The matched state is an input.
- Its support uses j_nd=4n(1+kappa d)+4d, an integer on these cards.
- A phase basis has exp(+ij phi); B phase basis has exp(-in theta).
- Calibrations are tau_A=hbar phi/delta and tau_B=hbar theta/epsilon.
- Rate r_B/A=1+kappa d is the slope of the conditional likelihood peak,
  on a declared unwrapped phase branch. A single read is not an exact rate.

The rational card supports an exact finite stationary construction;
no arbitrary-spectrum completion is inferred. Positive clock gaps
are required (kappa>-1). No negative physical clock energy is needed.

## Main findings and denial conditions

1. A relative-rate effect requires differential coupling/calibration.
   Equal multiplicative response of the two clocks gives ratio one.
2. Sign is dynamical: positive coupling speeds B relative to A;
   negative coupling slows it on the declared positive-gap domain.
3. Source fluctuations are not equivalent to substituting mean Ed.
   The exact unread-source P0=5/9 differs from (3+2sqrt(2))/9.
4. B alone cannot distinguish source coherence from a classical mixture
   with the same populations under this nondemolition coupling.
5. Ground-energy offsets and common energy-unit changes do not change
   the relational prediction.
6. The separate Jaynes-Cummings check derives the dispersive term and
   both signs from detuning. It is not an exact implementation of the
   rational card or a transmon design; dressed operators matter.
7. Nonuniform shifts of different clock transitions refute a rigid
   universal-rate reading of this mechanism, even when a frequency
   shift itself is real.

## Relation to saved records

G-R2 established a finite initial-record/final-readout protocol and
the need for epoch/event information at recurrence. It did not supply
an energy-complete apparatus. G-R3 does not silently identify its
stationary energy constraint with the G-R2 history constraint.
Event stamps, phase calibration and measurement backreaction remain
explicit; a joint autonomous implementation is still absent.

## Deliverables and next decision

[Derivation](tex/route_g_density_clocks.tex),
[PDF](output/pdf/route_g_density_clocks.pdf),
[executable](code/verify_gr3_density_clocks.py),
[report](output/gr3_density_clocks.md).

## Subsequent G-R4 physical control (2026-09-24)

The source/sign and cross-transition audit is now completed separately:
[pinned BBR inputs](data/BBR_CLOCK_INPUTS.json),
[derivation](tex/route_g_bbr_universality.tex),
[report](output/gr4_bbr_universality.md).
One isotropic thermal bath shifts two Yb+ transitions through -d.E.
Measured differential polarizabilities independently fix both signs
as negative, while their leading fractional responses differ by about7.25.
The selected thermal response is not solely a common clock rescaling.

u_gamma is energy per physical volume in conventional QED, NOT a
redefinition of this card's Ed or a derivation of emergent space.
An additional exactly common multiplicative factor cancels from the
frequency ratio, so this comparison cannot exclude it. Raw ratio data
and source calibration are absent; the result is a literature-pinned
prediction, not a new measurement. Next use such data or specify an
independently responding reference for the common-factor question.
Do not infer gravity, fit an Ed slowdown function, revive detector
scans, or promote Route-F gates.
