# Route G Roadmap — Relational-time mainline and comparison models

Created: 2026-09-22
G1 verified and PDF visually reviewed: 2026-09-23
Current priority (2026-09-24): relational/emergent time; G-R4 physical-source audit bounded-done.
Earlier G0–G2-T: preserved comparison/regression models, not the current mainline.
Status labels: open, in-progress, bounded-done, rejected-for-this-model.

## Authority and research discipline

The user explicitly chose Scheme A and requested a new Route G.
This is not the historical Route A and does not replace Route F.
Keep all input assumptions and claimed interpretations visible.
Do not import arbitrary particle-specific radii, projectors or mass shifts.
Do not turn every unfinished UV detail into another prerequisite paper.
No existing P54 parameters, fits or regression results are changed.

## New governing decision — G-R relational time (2026-09-23)

The user selected the approach closer to the source document: one
subsystem's events label another subsystem's changes, without assuming
a common physical external time. The old KK model explicitly assumed
Minkowski spacetime; G2-T calculated travel within that spacetime. Those
are not derivations of the source's time concept.

The original sections1.2 and1.8 use Ei accumulation, Ed and a time-event
quantity T. These are not silently identified with an energy eigenvalue,
ordinary energy density, flight time or a KK circle. The source's
negative-energy information and dimensional projection are likewise not
assigned to the clock Hamiltonians or Hilbert dimension.

This is a new conditional quantum construction, **not an edit of the
old KK action**. Preserve all earlier outputs and conditional exclusions.
No new source, detector scan, extra spatial dimension or Route-F bridge
is a prerequisite. Use [the interpretation card](RELATIONAL_TIME_CARD.md).

### G-R1 — Finite Page-Wootters clock and density-only counterexample

Status: bounded-done. Assumed inputs are ordinary quantum mechanics,
a clock/system tensor factorization, two complementary finite energy
spectra, a fixed-total-energy constraint, a matched state and a clock POVM.
They are not derived from Chen's text or from gravity.

For n=0,...,d−1, HC=epsilon(d−1−n), HS=epsilon n and E*=(d−1)epsilon,
the normalized state sum a_n |n,n> obeys K|Psi>=0 for
K=HC+HS−E*. Both local Hamiltonians are nonnegative.
Condition on |phi>=sum exp(i n phi)|n>/sqrt(d):
|psi(phi)>=sum a_n exp(−i n phi)|n> and
i epsilon d_phi |psi>=HS|psi>.
Only after calibration tau=hbar phi/epsilon does this become the
usual time-unit Schrödinger equation. No absolute physical gap is fitted.

The d clock-grid events are orthogonal; the continuous dial is a POVM,
not uncountably many orthogonal readings. The finite clock is cyclic
and cannot obey an exact canonical time/energy commutator. An averaged
joint effect commutes with K and reproduces the conditional probabilities.
One-event statistics are certified, not sequential measurement records.

New diagnostic result: a uniform three-level coherent state and its
energy-dephased mixture have identical local reduced density matrices,
complete joint energy statistics and total energy2epsilon with zero
variance. Yet for a fixed equal-superposition system projector,
P_coh(+|phi)=(1+2 cos(phi))²/9 and P_diag(+|phi)=1/3.
At phi=0 these are1 versus1/3. Thus an ordinary energy density alone
does not specify the relational history. This does not refute a richer
Ei/Ed containing additional state/coherence data, nor prove density
causes slowdown. Candidate Ei must retain relational information if
it is to encode this kind of conditional dynamics.

Common energy rescaling changes the conversion from phase to calibrated
time but not phase correlations; it is not relative clock dilation
without a reference. An energy-zero shift with a matching E* shift
also preserves K in this nongravitational model.

A general clock interaction yields an integral kernel in clock labels.
It does not automatically give a local Ed-dependent clock speed; the
joint state must be solved again and normalized conditioning rechecked.
The interaction projection identity is checked, but no interacting
stationary completion is claimed.

Artifacts: [interpretation card](RELATIONAL_TIME_CARD.md),
[report](output/gr1_relational_clock.md),
[code](code/verify_gr1_relational_clock.py),
[TeX](tex/route_g_relational_clock.tex),
[PDF](output/pdf/route_g_relational_clock.pdf).
New checks109/109; all4115 previous checks remain passing, total4224/4224.
The five-page PDF was compiled twice and every page visually reviewed;
no final layout or reference warnings remain. These checks certify
finite algebra and implementation, not physical realization.

### G-R2 — One explicit record and one final readout

Status: **bounded-done (2026-09-24)**. No G-R3 interaction is activated.

One initial unitary record and one terminal unitary readout act on a
three-level system and two blank-or-binary pointers. The same finite
model gives p(a,b)=tr[N_b U M_a rho M_a† U† N_b†]; it does not multiply
unperturbed one-reading probabilities. A fixed three-level example gives
the joint table [[1,2],[1,5]]/9 and terminal P(+)=2/9, versus zero when
the initial interaction is omitted. A different initial instrument
with exactly the same POVM effects gives P(+)=4/9. Thus both the event
effects AND their recording backreaction must be specified.

**Explicit change of assumptions:** an open history-state propagation
constraint K_hist replaces the noninteracting G-R1 constraint. It has
a unique zero mode on 162 dimensions, with ordered gates, an input
boundary and one epoch bit. G-R1's free U is retained, but its old
fixed-total-local-energy proof is not asserted for this new interaction.
The final readout does not commute with bare HS; apparatus/work
energetics and microscopic gate durations are not supplied.

The cyclic labels k=0,1,2,0,1,2 coexist with saved-event counts
N_evt=0,1,1,2,2,2. The two pointers carry fixed scheduled stamps
(epoch0,k1) and (epoch1,k0), not an arbitrary unknown timestamp.
Once written, their outcomes persist under the declared future links.
Conditioning only on recurrent k=0 mixes blank and completed histories;
including the epoch or completed-record marker resolves this ambiguity.
Full-state recurrence with the same blank input cannot also retain a
fresh nonblank record: the two states are orthogonal (closure residual
sqrt(2)). This does not rule out cyclic clocks with reset, longer cycles
or additional memory. No irreversible arrow, universal event counter
or empirical memory lifetime is derived.

Artifacts: [report](output/gr2_two_event_records.md),
[code](code/verify_gr2_two_event_records.py),
[TeX](tex/route_g_two_event_records.tex),
[PDF](output/pdf/route_g_two_event_records.pdf).
New checks139/139; complete Route-G regression4363/4363.
The five-page PDF is compiled twice and visually reviewed page by page.
These checks certify the declared finite protocol, not Chen's Ed law
or a physically energy-complete autonomous clock.

This stage is closed within its stated scope. Do not expand it into
an unlimited measurement-theory prerequisite. The subsequently authorized
G-R3 candidate/interaction/two-clock test is recorded below.

### G-R3 — An Ei/Ed candidate and an observable clock comparison

Status: **bounded-done (2026-09-24), for one declared candidate**.
No universal density-time law or original-paper identification is proved.

Define Ed=epsilon_D N_D/nu(D), with nu(D)=1 addressed matter cell and
ground-subtracted excitation energy. Units are energy per specified
cell, not assumed spatial volume. Full information still includes
the joint state, event algebra, instruments and reference measure;
Ed is only one observable. The one-cell example cannot experimentally
distinguish density from total excitation energy.

Choose V=kappa H_B Ed/E0, with E0=epsilon_D per cell and a shielded
reference A. The common stationary constraint K=H_A+H_B+H_D+V-E*
has an exact 15x3x2 matched state, positive energies, normalized phase
POVMs and constraint-invariant joint effects. The reference spectrum
is fixed for the synthetic kappa=-1/4,0,+1/4 audit cards. Commensurate
energies, preparation and phase calibration are assumptions, not fitted
data or a solution for arbitrary interacting clocks.

The conditional likelihood peak gives r_B/A=1+kappa d under fixed
calibration. At excited source occupation d=1 the cards yield 3/4,1,5/4.
Common energy rescaling leaves the ratio unchanged; equal fractional
response of two clocks cancels. A single finite-clock outcome does not
give an exact rate or an unambiguous cycle count.

**New physical insight:** density does not fix the direction of the
rate change, and density fluctuations do not generally define one
clock at the mean density. For the +1/4 card an unread equal source
mixture gives a specified outcome P0=5/9; replacing Ed by its mean
incorrectly gives (3+2sqrt(2))/9. Reduced visibility and multiple rate
branches must be retained. B-only data do not distinguish source
coherence from a classical mixture with the same populations.

A separate Jaynes-Cummings exchange calculation derives chi=g²/Delta
and kappa_eff=2chi/(omega-chi). Exact excitation-block spectra verify
that either detuning sign is possible and control the dispersive
energy remainder. This does not implement the rational kappa=1/4 card.
Higher-order gaps are transition dependent; dressed density/readout
operators matter. The result is an environment-dependent clock shift,
not universal time dilation, gravity or a fitted slowdown law.

Artifacts: [density/convention card](DENSITY_CLOCK_CARD.md),
[report](output/gr3_density_clocks.md),
[code](code/verify_gr3_density_clocks.py),
[TeX](tex/route_g_density_clocks.tex),
[PDF](output/pdf/route_g_density_clocks.pdf).
New checks273/273; complete Route-G regression4636/4636.
The six-page PDF is compiled twice and visually checked on every page.

The G-R2 history constraint and this stationary energy constraint are
not silently merged: an autonomous energy-complete recording apparatus
remains absent. That is a limitation, not a reason to reopen every
measurement-theory detail before reporting the bounded result.

### G-R4 — Independently fixed physical source and universality criterion

Status: **bounded-done (2026-09-24), source/response audit only**.
No new clock experiment, physical apparatus or universal-time effect.

Choose one isotropic blackbody radiation bath and two transitions of
the same 171Yb+ ion: E2 and E3 from the common S ground level.
H_int=-d.E fixes the mechanism. Independently measured differential
scalar polarizabilities, converted to alpha_excited-alpha_ground, are
6.9(1.4)e-40 and 0.888(0.016)e-40 J m^2/V^2. Both positive means both
thermal electric-dipole shifts are negative. The E2 paper originally
uses the opposite sign convention. These are pinned historical inputs,
not a claimed newest calibration or a fit to a temperature-cycle signal.

**Explicit boundary:** conventional radiation density u_gamma=a_rad T^4
has units J/m^3 and assumes physical space. It is NOT silently equated
to G-R3 energy per operational cell or Chen Ed. This is a physical
control on the interpretation, not a derived realization of emergent time.

From second-order perturbation theory derive
delta_nu_i=-(2h epsilon0)^(-1) integral Delta_alpha_i(omega) u_omega d_omega.
The leading fractional coefficient is
beta_i=-Delta_alpha_i(0)/(2h epsilon0 nu_i). At 300 K, the static scalar
electric-dipole prediction gives fractions -5.235e-16 and -7.223e-17:
their magnitudes differ by about7.25, not by numerical rounding.
The illustrative ratio modulation R(310 K)/R(300 K)-1 is +6.325e-17,
where R=nu_E3/nu_E2. This is a prediction, not an observed signal.

For a connected graph of levels, equal fractional gap shifts require
delta_H_secular=f H0+C I. The independently fixed thermal coefficients
fail that necessary condition. Thus this selected thermal response
cannot be explained solely as one universal multiplicative time factor.
An additional exactly common factor F(u) remains **unidentifiable**:
nu_i=F(u) nu_i^0(1+s_i) cancels F from the ratio exactly. Neither
nonconstant nor constant same-bath ratios alone settle universal time.

The leading T^4 coefficient mismatch is separate from full dynamic
corrections. The published E3 correction is recorded; complete E2
spectral response, covariance and apparatus controls remain absent.
A +/-3 published-sigma coefficient box plus assumed +/-10% independent
dynamic factors preserves the mismatch at all16 corners. This is a
conditional stress bound, not statistical significance or a derived
bound on omitted physics. Do not reuse BBR-corrected outputs as raw
temperature-response observations.

Artifacts: [pinned inputs](data/BBR_CLOCK_INPUTS.json),
[report](output/gr4_bbr_universality.md),
[code](code/verify_gr4_bbr_universality.py),
[TeX](tex/route_g_bbr_universality.tex),
[PDF](output/pdf/route_g_bbr_universality.pdf).
New checks91/91; complete Route-G regression4728/4728.
The old G-R3 script now counts274 checks because its automatic file-
preservation inventory also includes G-R4; its original273-check stage
snapshot above is historical. No earlier verifier code was changed.
The five-page PDF is compiled twice and visually checked page by page.
These validate the calculation, not an empirical discovery.

### Next physical decision — requires observations or a new observable

Priority A: actual common-bath radiometry plus interleaved raw E2/E3
readings at two setpoints, retaining the intentional BBR modulation
while bounding temperature-correlated non-BBR systematics. No data
are fabricated and no additional temperature/Gaussian scan is needed.

Priority B, if the target is a universal factor: specify an independent
reference outside the modified bath or a predicted non-clock observable,
including its response and comparison protocol. Assuming an unaffected
reference by definition cannot establish a universal time effect.
Optional C: equal total radiation density with different spectral content,
to test scalar-density sufficiency via the same spectral response kernel.

Use interaction-derived spectral response as the current explanation;
do not fit an arbitrary Ed slowdown law or add per-species free factors.
No conventional-spacetime control is promoted to an emergent-space proof.

### Beyond the bounded clock program

Space, Lorentzian causal structure, universal proper time, gravity,
mass generation, original negative-mass/time-jump claims, and a physical
implementation remain unproved. They are not all made mandatory gates
for G-R2, and no existing Route-F result is promoted.

## Preserved comparison program — G0–G2-T

The historical priorities and next steps below apply only within their
declared conventional-spacetime branch. They do not supersede G-R or
silently become prerequisites for relational time. G3/G4 at the end
likewise remain optional continuations of the old KK prototype.

## G0 — Action and interpretation card

Status: bounded-done.

- One fixed M4 × S1 background, one complex scalar, positive kinetic term.
- Shared R, M5 and alpha; alpha modulo integers is a flat holonomy.
- One pair of source couplings at 0 and pi, with a specified path and
  parallel transport; no per-mode detector tuning.
- Four numerical engineering cards, chosen without experimental masses.
- Scope: free-field spectrum and linear source response only. No dynamical
  gravity, gauge bosons, radion potential or stabilization is assumed.
- The Chen final document's common-energy-level idea motivates the
  candidate; its other assertions are not treated as established axioms.

## G1 — Spectrum, visibility and covariance

Status: bounded-done; implementation independently rerun, 331/331 checks.

Derive and verify:

1. The KK mass tower and positive-energy Hamiltonian.
2. The fixed probes' full 2 × 2 residue matrix, including cross-response.
3. Passive basis covariance, small-gauge covariance and integer
   large-gauge relabeling.
4. A finite-difference spectrum independent of the analytic formula,
   and a bounded KK response truncation.
5. The two physically meaningful restrictions:
   - signed-branch squared masses have second difference 2/R²;
   - the two probe residues have a nonzero fixed sum for every signed
     Fourier branch, not every superposition in a degenerate subspace.

Do not misread adjacent SORTED energy levels as adjacent signed KK labels.
At degeneracy, compare spectral projectors/residue sums, not arbitrary
eigenvectors. Relative probe weights are not branching fractions.

Completion requires one executable report with spectra, response and
failed-claim flags. Numerical checks certify the implementation, not the
existence of an extra dimension in nature.

Delivered:

- One concise action-to-response derivation, TeX and typeset PDF.
- Four predeclared cards, full complex two-channel residues and responses.
- Independent 32/64/128-site diagonalization. The last spacing-halving
  error ratios are 3.9846--3.9882, as required by the O(h²) symbol.
- Gauge-transformed response error <=1.4e-15 in the recorded run;
  passive-basis response error <=2.8e-16.
- Infinite image-sum reference and rigorous per-entry truncation bound,
  tested at 16/32/64/128/256 modes on each side and three momenta.
- A physically important correction: at alpha=0 the n=±1 coupling map
  has rank 1, so a sine standing wave is dark to BOTH probes. At alpha=1/2
  the n=-1,0 map has rank 2. Full pole residues, not an arbitrary choice
  of degenerate eigenvectors, are basis-invariant.

This milestone is closed within its declared free-field scope. Do not
reopen it merely to gain more decimal precision or fit particle masses.

## G2 — One genuine transition, not a change of coordinates

Status: bounded-done (tree-level EFT); 262/262 independently rerun
implementation checks. G1 regression: 331/331.

Chosen action: add a neutral complex dynamical field X with positive
kinetic term, mass MD and potential kappa5 |Phi|² |X|², kappa5>=0.
This is local in five dimensions. Localization belongs to the detector
quantum state, not an externally pinned coupling or infinite-mass sink.
The classical vacuum remains Phi=X=0; G1 is the tree quadratic baseline.

The exact circle overlap gives g=kappa5/(2 pi R) and

\[
\phi_n+X_l\longrightarrow\phi_m+X_j,\quad j=l+n-m,\quad
\mathcal M_{\rm tree}=-g.
\]

Every component obeys four-momentum and compact-momentum conservation.
The open-channel cross section is g² p_f/(16 pi s p_i); no identical
particle factor is needed. Inverse detailed balance uses the SAME s.

The earlier G1 candidate remains rejected kinematically: for the
bulk vertex chi Phi-dagger Phi with a massless neutral five-dimensional
chi, KK momentum conservation requires q=n-m and m_chi=|n-m|/R.
The triangle inequality gives

\[
m_n < m_m+|n-m|/R \quad(M_5>0,\ n\ne m).
\]

Thus the proposed one-to-two decay is closed. At M5=0 it is forbidden
or threshold-only, not an open two-body decay. This excludes that process,
not all interactions or scattering. No stable interacting completion
is inferred from a cubic vertex alone.

The G2 repair is physical: for l=0 the detector's internal recoil costs
sqrt(MD²+(n-m)²/R²)-MD rather than |n-m|/R. With
Delta=mn-mm>0 the rest-limit channel is exothermic iff

\[
M_D>\frac{(n-m)^2/R^2-\Delta^2}{2\Delta}.
\]

The engineering card R=1,M5=.5,alpha=.25,MD=5,g=.05,p_i=.2 gives
Phi_1+X_0 -> Phi_0+X_1 with p_f=1.02175458 and
sigma=6.27164081e-6. Its rest-limit critical MD is only 0.24146563:
this opening is not a near-threshold decimal adjustment.
The released Phi rest energy .78727421 pays .09901951 internal detector
recoil plus .68825469 additional ordinary kinetic energy.

Preparation: a finite normalized detector packet with l=-4,...,4,
c_l proportional to exp[-l²/(4*1.2²)] exp[-i l*.7]. It freely spreads.
The exact first-order packet amplitude and traced-recoil probability
are derived in the new TeX. Different l give orthogonal outgoing j
for fixed n,m. The inclusive rate cannot measure the initial packet
phase or center against a single incident KK plane wave.

Numerical scope: component cross sections, all open final modes,
conservation laws, inverse detailed balance, localization diagnostics
and a common-overlap averaged sigma*v. Do NOT label the latter an
absolute finite-encounter probability; four-dimensional envelopes and
encounter geometry are still needed. The five-dimensional EFT has no
certified loop/cutoff uncertainty or radiative naturalness result.
No interacting all-orders G1 identity or physical particle fit is promoted.

Results: the common-overlap packet rate is 1.11933296e-6 for target m=0
and 1.85661445e-6 summed over all scattering channels, including elastic.
The 0.60288929 ratio is a prepared-state rate fraction, not an absolute
encounter probability. Conditional incoming/outgoing detector mean
momenta differ by exactly +1; weighting both by the selected event rate
avoids confusing preparation bias with recoil.
Free packet localization decreases from .91640536 to .48041872 in
the first circular moment between t=0 and 5, with energy and norm
unchanged. This is not a pinned defect or a derived localized ground state.

Artifacts: [derivation](tex/route_g_local_conversion.tex),
[PDF](output/pdf/route_g_local_conversion.pdf),
[numerical report](output/g2_local_conversion.md),
[executable](code/verify_g2_local_conversion.py).

### G2-R — Recoil-resolved continuation

Status: bounded-done at tree level (2026-09-23). The reproducible verifier
passes 357/357 checks, independently rerun; the four-page derivation PDF
has been compiled and visually reviewed. This closes the joint prediction
under the declared preparation, not an experimental detector model.

Keep the G2 action, g, all masses, alpha, common incoming momentum and
finite detector packet unchanged. For every open channel, compute
W(m,j)=w_l K(l,m), l=m+j-n. Normalize separately over all scattering
events and over mode-changing events m!=n. Elastic scattering is not
the identity/no-encounter probability.

Deliver full signed joint matrices, mass-only folded matrices, marginals,
Bayes conditionals, event-conditioned incoming l, ordinary recoil
energies/momenta and the isotropic tree angular law. All 25 open channels
are retained. The all-scattering output weights for m=(-2,-1,0,1) are
(.00785985,.27119405,.60288929,.11805680); m=1 is elastic.
Mode-changing events comprise .88194320 of scattering events under
the common-overlap prescription, not of arbitrary incident packets.

New physical insight: the X mass alone determines |j| and loses its
sign. At this nondegenerate Phi card, a known incoming p_i and the
outgoing ordinary p_f can reconstruct the sign for m!=n. With d=m-n,

\[
l^2=R^2[(E_{\rm out}-E_n)^2-M_D^2-p_i^2],
\qquad j=(l^2-d^2-|j|^2)/(2d).
\]

This is MODEL-DEPENDENT reconstruction. Do not infer j from this formula
and then claim the assumed conservation law has been independently
verified. Independent observables are the predicted joint line positions,
fractions and angular dependence. The m=0, |j|=1 pair has ordinary
momenta 1.02175458 and 1.36191317 (gap .34015860). A deterministic
momentum error below half the gap is a sufficient ideal resolution
requirement, not an assigned experimental uncertainty.
Elastic events cannot resolve the sign this way. At a degenerate Phi
holonomy, even the initial output-mode assignment needs a new audit.

No packet phase or localization is witnessed by this diagonal joint
distribution: the dephased preparation gives the same prediction.
No actual detector response, efficiency, loop-error model or absolute
4D encounter probability has been calculated. No parameters are fitted.

Artifacts: [joint report](output/g2_recoil_joint.md),
[TeX](tex/route_g_recoil_joint.tex),
[PDF](output/pdf/route_g_recoil_joint.pdf),
[executable](code/verify_g2_recoil_joint.py).

The finite-readout continuation and activated drive are now recorded
below. They do not overwrite this undriven joint prediction.

### G2-S — Finite-resolution recoil discrimination

Status: bounded-done (2026-09-23); 272/272 implementation checks.

Use a classical incoherent incoming COM momentum ensemble with
parent-Gaussian location p0=.2, widths sigma_i=0,.02,.1,.2,.4, truncated
to p>=0. sigma_i is not the truncated distribution's actual deviation.
Use a Gaussian real-valued reconstructed momentum estimator with
sigma_d=.02,.05,.1,.2,.4; a negative fitted estimator is allowed.
Mass labels m=0 and |j|=1 are assumed exact, with equal acceptance,
common COM frame and the same external-overlap prescription.

Integrate the unnormalized rates w_l f(p) K_lmq(p) R(x|p_f) first.
Then normalize the selected pair, derive the posterior and compute the
optimal equal-cost sign error integral min(h_plus,h_minus)/Z.
The engineering decision criterion is <=5%, not a discovery threshold.
Compare to the event-selected majority-only baseline (~22.5%).
The 50-card scan covers both undriven and driven cases, with ideal
independent q-tagging as a mathematical information bound only.

| sigma_i | sigma_d | Undriven error | Driven q-unknown error | Ideal q-tag error |
|---:|---:|---:|---:|---:|
| .1 | .1 | 3.77384% | 5.94485% | 3.77589% |
| .2 | .1 | 4.84993% | 6.83823% | 4.85210% |
| .4 | .1 | 11.0294% | 12.2710% | 11.0331% |

At (.1,.1), the undriven card passes the declared 5% criterion, but
the driven card without q information fails. At (.1,.05) the driven
untagged error drops to 2.20225%. Extra detector Gaussian noise cannot
improve optimal classification (L1 contraction); beam broadening also
changes selected priors, so do not assert the same general monotonicity.

Numerical controls: exact sharp-Gaussian Bayes comparison, G2-R hash/rate
regressions, truncated-normal moments, nonlinear vs narrow-width
derivatives, quadrature-order/boundary-grid comparison, and analytic
incident/response tail bounds. Independent adaptive integration confirms
the main card. Root searches are convergence-tested, not certified
exhaustive for arbitrary mixture parameters. Gaussian tails are bounded
under the formal tree kernel, not evidence for infinite UV validity.

Artifacts: [report](output/g2_resolution.md),
[code](code/verify_g2_resolution.py), [full data](output/g2_resolution.json).

### G2-D — Activated time-dependent energy reservoir

Status: bounded-done at first Born/long-time order (2026-09-23);
565/565 checks. User activation supersedes the former backup-only status.

Choose the minimal explicit local modulation
g(t)=g[1+epsilon cos(Omega t)], epsilon=.5, Omega=.2.
It is spatially/circle homogeneous, with pump rest frame equal to the
incoming COM. The instantaneous quartic remains positive; free masses
are unchanged. The new external reservoir is a hypothesis, not a
closed UV completion or dynamically selected geometry.

At one vertex g0=g, g(+/-1)=g epsilon/2, Eout=Ein+q Omega, and
K_lmq=|gq|^2 p_f/(16 pi E_n E_l Eout). Spatial and compact momenta
remain conserved. Every event records pump work q Omega; matter alone
is not energetically closed. The epsilon=0 rate regresses to G2.

All 66 open rows are retained (14/25/27 in q=-1/0/+1).
Total K=2.10563491267e-6; mode-conversion share=.87381701.
Mean absorbed pump work/scattering=.00900097215 and work coefficient
1.89527612162e-8 are not a luminosity-normalized power.
Mode-preserving m=n,q!=0 events are inelastic; only m=n,q=0 is elastic.

The selected two passive lines become six components. Opposite-sign
sidebands at p=1.20257369 and 1.18924247 have gap .01333122.
No exact degeneracy occurs at the chosen card; perfect sharp-beam
readout still distinguishes the full catalogue. The passive inversion
must include unknown q, and finite resolution can worsen classification.
Omega=.1,.2,.4 are sensitivity cards, not a frequency fit.
The analytic exact-crossing frequency .19243951642 is a diagnostic,
not used to tune the background.

The code checks first-order Fourier coefficients, independent radial
phase space, full enumeration, gauge/compact/spatial conservation,
on-shell pump work, positive quartic and finite top-hat Fourier factors.
The latter are NOT finite-collision probabilities. Fixed Omega>0 and
long-time rates are assumed; Omega=0 requires coherent harmonic addition.
Coherent finite beams/pulses need amplitude-level interference and cannot
be replaced by this classical momentum histogram.

Artifacts: [drive report](output/g2_driven.md),
[code](code/verify_g2_driven.py), [full data](output/g2_driven.json),
[shared derivation](tex/route_g_resolution_drive.tex),
[PDF](output/pdf/route_g_resolution_drive.pdf).

The readout choice is made in G2-M below. Finite 4D packets remain
unnecessary for these conditional long-time rates; they are needed
for absolute encounters, finite pulses or coherent interference.

Checkpoint 6bbd6fa saved and pushed G2-R before this continuation,
as requested. G2-S/D does not promote a Route-F fit or identify SM fields.

### G2-M — Momentum-only readout selected; noisy tag comparison

Status: bounded-done as an engineering contract (2026-09-23).
Frozen-readout checks: 154/154. Noisy-tag checks: 125/125.
Both verifiers were independently rerun; the four-page TeX/PDF was
compiled and visually reviewed.
Physical apparatus realization remains open, not claimed by a Gaussian
response or the name "detector" attached to the recoil field X.

Main contract: unchanged driven action, sigma_i=.1, nominal sigma_d=.05,
same ideal m=0, |j|=1 selection. Freeze the nominal rule
j_hat=+1 iff x<1.25892910. Nominal total error is 2.20225%, with
class-conditional errors (1.0315%,6.2279%) for true signs (+1,-1).
Aggregate <=5% does NOT establish >=95% recall for each sign.
Do not silently change the loss function or priors to hide this.

Calibration response: x=(1+a)p_f+b+N(0,sigma^2).
The continuous stress box is |a|<=.01, |b|<=.02, .04<=sigma<=.06,
in model units. The classifier is NOT retrained at each unknown
calibration. The worst of 27 sampled cards is 3.8648%, not a proof
of a global supremum. The unsplit CDF bound is loose (5.3902%).
One predetermined 2x2x2 shared-nuisance cover reduces the upper envelope
to about 4.46%, sufficient for the overall 5% criterion.
Kink-aware integration is convergence-tested with a separate analytic
tail bound, not interval-certified. No repeated refinement or fit.

Oracle retraining is displayed separately. With true a,b known,
an affine transformation makes oracle error depend only on sigma/(1+a).
This is not deployed robustness. At a=b=0 the nominal 5% noise limits
are about .0907 with retraining and .0891 with the frozen rule;
they are not guarantees for arbitrary calibration.

Alternative: T_eta(z|q)=(1-eta)delta_zq+eta/3, true mistag=2eta/3.
At sigma_d=.1 the 5% crossing lies in eta=[.1758,.1768].
This requires the FULL specified matrix, not merely 88.2% tag accuracy:
constant z=0 already scores 88.98% but adds no information and leaves
5.94485% recoil-sign error. Tests preserve rates/priors and verify
stochastic degradation and ideal/no-tag limits. These tag comparisons
assume a known calibrated matrix and its optimized classifier,
not robustness against unknown tag drift.
At primary sigma_d=.05 no extra tag is needed for the overall criterion.

A tag made only from the same measured x cannot add information;
a drive clock records phase, not event work. Extra sensor data need
one joint likelihood over shared latent momentum, not a product of
separately marginalized likelihoods. No real sensor coupling, physical
energy unit, mass-confusion/acceptance calibration or quantum pump
counter has been derived.

Artifacts: [frozen-readout report](output/g2_readout_contract.md),
[code](code/verify_g2_readout_contract.py),
[noisy-tag report](output/g2_noisy_tag.md),
[tag code](code/verify_g2_noisy_tag.py),
[derivation](tex/route_g_readout_contract.tex),
[PDF](output/pdf/route_g_readout_contract.pdf).

Next physical handoff: choose a concrete energy scale/realization and
sensor interaction, or provide an actual apparatus's beam/response
calibration. Decide explicitly whether the target is overall error or
per-sign recall. Do not introduce another mathematical prerequisite
or finer scan in place of that missing physical choice.
Only an absolute-yield/finite-pulse/coherence question activates a
four-dimensional packet and pump-envelope specification.

### G2-P — Physical energy anchor and visible-probe candidate

Status: bounded-done as a scale/coupling audit (2026-09-23);
physical realization and measured calibration remain open.
Scale checks: 877/877; candidate probe/matching checks: 113/113.
Both were independently rerun. Existing G1/G2/G2-R/G2-S/G2-D/G2-M/tag
regressions also pass (2066/2066 combined). The four-page TeX/PDF was
compiled and every page visually reviewed. These are implementation
checks, not empirical validation or calibrated hardware qualification.

The current dimensionless forecast has a scale degeneracy. With
Lambda=1/R, energies and resolutions scale as Lambda, K as Lambda^-2,
mean work/event as Lambda, and normalized classification error is
unchanged. The 877/877 synthetic checks rescale every one of the 66
driven open channels and recover the same 2.20225% error. No eV/GeV
choice or measured calibration is inferred from this exercise.

Three externally calibrated consecutive signed-mode mass² values
determine A=Lambda², B=2 alpha Lambda², C=M_eff²+alpha² Lambda².
The fourth must satisfy Delta³(m_n²)=0 within specified covariance and
theory error. One mass only anchors Lambda if all dimensionless ratios
and its mode identification are already assumed. Sorted mass order is
not signed-mode order. A common positive mass² shift preserves this
closure test and cannot repair a wrong spacing pattern.

Candidate, not yet adopted: an SM-singlet tower with a uniform
`-H†H(lambda_Phi sum|phi_n|² + lambda_X sum|x_j|²)` interaction.
It is a KK-diagonal, charge-preserving 4D EFT. Writing one 4D Higgs
uniformly around the circle is not a completed local 5D SM embedding.
A brane-localized alternative is off diagonal, mixes masses and lets
the apparatus absorb compact momentum; old full-system selection rules
cannot be retained without accounting for that change.

The candidate audit derives common EWSB shifts, the finite-t free-nucleon
kernel, its NR limit, and the complete kinematically open complex-tower
Higgs width. Their shared coupling yields a conditional upper bound on
scattering from an applicable invisible-width budget. No cosmological
abundance, experimental limit value or instrument performance is imported.
At fixed bare masses the bound is implicit because thresholds change.

Two corrections prevent false promotion:

1. Keeping physical MD²=25 Lambda² with the retained action (no X
   self-quartic) requires lambda_X<=50 Lambda²/v². A negative underlying
   X mass runs away along H=Phi=0; adding a stabilizer changes the model.
2. If both portals are nonzero, Higgs exchange adds an elastic Phi-X
   amplitude that interferes with -g. Old inclusive fractions then need
   recomputation. An X-only coupling at a stated matching scale is a
   minimal candidate; a zero Phi portal is not radiatively protected.

The collision and readout verifiers are unchanged. These normalization
checks are synthetic, not evidence of physical realizability. A nucleon
cross section is not an instrument response: one recoil deposit retains
an unknown angle; nuclear response, efficiencies, geometry, flux and
mass-label confusion remain physical inputs. The drive is still a
prescribed reservoir, not an event-resolved work sensor.

Deliverables: [physical input card](PHYSICAL_INPUT_CARD.md),
[scale code/report](output/g2_scale_anchor.md),
[probe code/report](output/g2_probe_coupling.md),
[TeX](tex/route_g_physical_anchor.tex),
[PDF](output/pdf/route_g_physical_anchor.pdf).

Next: choose a neutral-sector physical scale and target to compare a
required interaction strength with stability and experimental constraints;
alternatively supply actual calibration data, or explicitly select an
analogue-mode platform with a more limited interpretation. Do not replace
this physical choice by a finer resolution scan or an unrequested UV
construction. If the required coupling is incompatible with the declared
constraints, reject that candidate. No Route-F promotion follows.

### G2-X — External scale/target choice and actual response fold

Status: bounded-done (2026-09-23); a conditional feasibility calculation,
not a completed production/detection apparatus.

Authorized choice: X-only circle-uniform Higgs coupling at a stated tree
matching prescription, physical intercept MD=5 Lambda, Lambda=1–30 GeV,
natural liquid xenon. External motivation is the Higgs mass threshold
Lambda_h=mh/10=12.508 GeV, not a fit of a separate radius to each particle.
Input mh125.08 GeV, v246 GeV, GammaSM4.10 MeV comes from the pinned PDG
review provenance; use ATLAS2023 observed B_inv<.107 at95%CL with SM
production/escape/drive-off assumptions. Sum all open complex X modes.

Actual response is prioritized: use the official XENON1T S2-only2019/2020
release, commit5a364bc8709f2561e5a013ddea6993a5a7c8e313. Keep raw CSVs,
hashes, license metadata and usage instructions. Restrict true NR to
.7–50 keV and full S2 bins29–198 (150.027–3000 PE). The exact raw search
exposure is356770 kg day, and the matrix already includes selection
losses. No extra fiducial fraction, fake calibration, halo density,
halo-WIMP exclusion or new Gaussian-resolution scan is introduced.

The coherent scalar-current/Helm kernel keeps finite incident kinematics
and natural-isotope weights. It is only used at admitted low momentum
transfer; unmodeled low/high recoil contributions are not assigned zero
physical response or included through uncontrolled total-cross-section
fractions. Source weights remain the rate-conditioned G2 driven mixture.

Map three distinct items: conditional Higgs exclusion; failure of the
retained action when lambda_X>50 Lambda²/v²; and the at-detector selected
particle flux required for three expected accepted-window recoils.
No actual flux is known. The third item is not discovery sensitivity,
95% exclusion, or an achievable beam. Response, form-factor, matching
and source uncertainties are not replaced by numerical quadrature error.

New readout result: p_X,min=.8091444644 Lambda gives a recoil endpoint
above the complete response window throughout the chosen range. Thus
d_sigma/d_ER=(lambda_X²/p_X²)H_A(ER); every common energy-only response
preserves identical normalized S2 shapes for j=+1 and j=-1. Accepted
priors are reweighted by1/p_X², but S2 provides no further sign information
within this kernel/window. Accepted priors are 85.943%/14.057%; the best
S2-only error is **14.057%**, merely the minority prior. The minimum
endpoint is 9.543 MeV. Monotonicity in Lambda bounds the whole interval,
not just numerical sample points. Reject transferring the old2.20% momentum
classifier claim to this sensor; it is not repaired by finer Gaussian fits.

The folded cross section is 7.99989e-39 cm² times
lambda_X²/(Lambda/GeV)². At Lambda=1, 5, 30 GeV the three-event
requirements at the conditional coupling ceiling are 5.07e9, 2.15e10,
4.32e6 cm⁻²s⁻¹. The minimum here optimizes only the declared coupling
domain, not production or full experimental sensitivity. The new
response/kinematics/tower verifier passes 573/573; all prior 3056
checks pass unchanged, total **3629/3629**. The figure uses analytic
formula rendering, not additional Gaussian scans. Data, formula,
implementation and plot hashes are preserved for reproduction.
The six-page derivation PDF and standalone figure were rendered and
visually reviewed; the final TeX build has no layout/reference warnings.

Physical source issue exposed by setting units: the activated drive
Omega=.2–6 GeV corresponds to frequencies about4.84e22–1.45e24 Hz.
A prescribed pump is not a realized supply or an ordinary oscillator.
The drive remains in the mathematical branch, but no hardware promotion
is made. X-only is not loop-protected; this is a tree feasibility test,
not a claim of completed quantum matching or certified UV control.

Next decision is now between two physical aims, not another resolution
scan: (1) event-count search, requiring production/at-detector flux and
background inference; (2) sign identification, requiring an independent
timing/momentum channel, direction or endpoint-sensitive target/readout.
If a real source cannot meet the required flux, reject that realization;
if sign separation is mandatory, this retained S2-only channel fails it.

Artifacts: [figure](output/figures/g2_xenon_feasibility.png),
[report](output/g2_xenon_feasibility.md),
[TeX](tex/route_g_xenon_feasibility.tex),
[PDF](output/pdf/route_g_xenon_feasibility.pdf),
[raw-response provenance](data/xenon1t_s2only/PROVENANCE.md).
Route-F physics and fits remain independent and unpromoted.

### G2-H — Physical source budget and production/readout compatibility

Status: bounded-done (2026-09-23). The on-shell-Higgs/single-pass-xenon
retained elastic NR realization is **rejected-for-this-model** at the
published Run-2 source budget. This is not a theory-wide no-go.

Source selection: ordinary pp→h+anything, h→X_j anti-X_j, no classical
drive and no assumption of the old selected Phi-X collision ensemble.
Use one actual13 TeV interaction-point dataset,139 fb^-1, and the
SM production input55.6 pb from ATLAS2023. The measured visible-channel
extraction assumes SM branching ratios; it is not treated as a separate
model-independent normalization after changing invisible decays.
The retained conditional B_inv<=.107 allowance gives
N_h=7.7284e6 and N_X+anti-X<=1.6538776e6. No actual dark flux is observed.

Source weights now follow B_j=Gamma_j/(GammaSM+sum Gamma), not the
old rate-conditioned collision priors. The source-to-fluence expression
includes the actual Higgs boost distribution, geometry/transport and
time profile as physical inputs. No full correlated distribution is
fabricated. The integrated budget must not be multiplied by the old
kg-day exposure again.

The decisive result bypasses unknown beam details without assuming them:
for every momentum, the retained-window scalar cross section obeys
sigma_acc<=lambda_X² fN² mN²/(4pi mh⁴)
sum_A xi_A A² mA²/(mX+mA)². This follows by bounding only the admitted
kernel, not by extrapolating coherent physics to all recoil energies.
Give all produced particles the maximum chord of the published xenon
cylinder and readout efficiency one. Use mX>=5Lambda, the existing bare
mass condition and monotonic Lambda⁴/(5Lambda+mA)². Across1–30 GeV,
**N_acc<3.96e-6**; an explicit generous input stress gives5.35e-6.
This is an upper envelope, not a simultaneously attainable threshold
point, an experimental confidence limit or a measured event rate.
Above12.508 GeV the on-shell source is closed; the abs(j)=1 source
already closes at12.2651 GeV.

New information gate: for direct scalar Higgs production and the
diagonal X-only probe, opposite signs have equal priors and identical
full four-dimensional kinematic distributions. Timing, direction and
independent momentum cannot distinguish them: I(j;z)=0 and optimal
error50%. This is source/channel-specific, not an all-orders symmetry
claim for the Phi-holonomy action. Do not transfer the old asymmetric
source's14.057% S2 prior error or2.20% ideal momentum error.

Implementation: new source/bound/endpoint/response/Lorentz checks389/389;
all prior3629 checks remain passing, total4018/4018. Actual NR response
folds are tested for slow particles whose endpoints cross bin edges;
the principal arbitrary-position bound correctly uses efficiency<=1.
No Gaussian optimization, beam Monte Carlo or added coupling is used.
The four-page PDF was compiled twice and every page visually reviewed;
the final build has no layout or unresolved-reference warnings.

Stop: do not refine this rejected chain's beam geometry or sensor noise.
Next options are (1) production/missing-momentum observables that avoid
a second rare interaction; (2) a physical sign-sensitive preparation or
coupling if signed conversion is required. On-shell invisible-Higgs
production measures only an inclusive width, not the tower or signs.
Off-shell sources, astrophysical populations, other targets, electronic/
inelastic and outside-window recoils, and recirculation are untested,
not quietly ruled out or appended to the same model.

Artifacts: [source card](data/HIGGS_SOURCE_CARD.md),
[report](output/g2_higgs_source.md),
[TeX](tex/route_g_higgs_source.tex),
[PDF](output/pdf/route_g_higgs_source.pdf).
Route F and its gates remain independent.

### G2-T — TOF/momentum metrology and signed-mode audit

Status: bounded-done (2026-09-23), **mathematical measurement conditions
only**. Actual source timing, independent X momentum/energy sensing,
directional calibration and signed preparation remain unspecified.
The action, source strength and prior response are unchanged.

The user selected two separate lanes:

1. Different masses: derive m²=p²[(cT/L)²−1], its full logarithmic
   covariance and exact bounded-error intervals. TOF alone fixes m/p;
   two timed hits still do not determine mass without another constraint.
   Unknown t0 and source-time/momentum correlations must enter the joint
   likelihood. No Gaussian or independent-source-time assumption is made.
2. Equal mass, opposite sign: retain the source-specific indistinguishability
   theorem. Direct Higgs production plus the diagonal portal remains
   sign-blind even with perfect timing, direction and momentum.
   The old asymmetric collision is a different preparation.

New practical obstruction: X is neutral; a charged-track momentum
resolution cannot be imported. A missing-momentum pair sum is not an
individual-particle momentum. Existing S2 response is not X TOF, and
electron drift time is not the source-to-interaction flight time.

A same-action alternative is derived, not promoted to hardware.
For known target mass M at rest, q²=ER(ER+2M) and
D=beta q cos(theta)−ER,
E=M ER/D, p=beta E, m=E sqrt(1−beta²).
It requires recoil direction, incident direction, beta, energy and a
target-isotope treatment. At Lambda=5 GeV, E=62.54 GeV, A=132, ER=10 keV,
the angular-only 1% mass-error budget is about0.896 arcsec.
This conditioning result, with other inputs exact, is not a measured
detector resolution or a recommendation to optimize xenon directionality.

The at-rest,10 m illustration gives a139.214 ps j=0/abs(j)=1 time gap.
With bounded dp/p<=.001 and dL<=1 mm, the exact reconstructed-mass boxes
first touch at a timing halfwidth48.672 ps;20 ps boxes are disjoint.
These are illustrative requirements, not sigmas or collider confidence
intervals. The +1/−1 boxes are identical.

Revalidated, not rediscovered: G2-R's undriven source correlation
j=n+l−m, or its full COM kinematic inversion for unobserved l and
m−n nonzero. A real Phi source/tag is still missing; driven events need
the work exchange in the inference. It is not an independent test of
the conservation law used to infer j.

New-interaction guardrail:
|Me+j Mo|²−|Me−j Mo|²=4j Re(Me* Mo).
A sign-odd amplitude by itself, a relative phase pi/2, or orthogonal
final states do not supply the required interference. Any new coupling
must exhibit an observable signed response and its physical reference.
No new interaction or external field is activated in this milestone.

The joint source/response measure includes arbitrary t0–p correlations.
Additional gates lie in[0,1], so the previous arbitrary-momentum,
efficiency-one retained-elastic-event ceiling3.9596863e-6 is unchanged.
No revived LHC-to-xenon experiment, new flux or Gaussian scan is claimed.

Stop/next: mass lane needs a real source tag plus independent momentum/
energy mechanism and a surviving event budget before response forecasts.
Sign lane first needs a physically accessible source correlation, or
an explicitly proposed interaction for separate audit. Do not refine
the rejected Higgs/xenon chain; collider production observables avoiding
a second tiny scattering probability remain a distinct useful lead.
Neither lane is a new Route-F prerequisite.

Artifacts: [input card](data/TOF_READOUT_CARD.md),
[report](output/g2_tof_mode_audit.md),
[TeX](tex/route_g_tof_mode_audit.tex),
[PDF](output/pdf/route_g_tof_mode_audit.pdf).
Verification:97/97 new checks and all4018 prior checks pass, total4115.
The five-page PDF was compiled twice, has no final layout/reference
warnings, and every page was visually reviewed. Algebra tests are not
empirical validation, and none of the synthetic requirements is calibration.

## G3 — Background selection

Status: open, not a prerequisite for the conditional G1 calculation.

Specify what fixes R and alpha, then check a stable stationary background
and whether the field content actually favors that holonomy. A potential
chosen solely to force a desired particle mass is not a derivation.
Do not assume a freely adjustable nonzero alpha is a selected vacuum.

## G4 — Particle-content and Route-F interface

Status: open; no identification presently authorized.

Before comparing with elementary-particle masses, derive spin, charges,
chirality and multiplicities from a single field/boundary content.
The present scalar KK tower cannot be relabeled as the SM fermions.
Any link to Route-F h_D, f_D, Higgs overlaps and Majorana scale must come
from shared normalized modes and one action, not independent fits.
Full string embedding is optional and separate, not a new mandatory gate.

## Explicit rejection conditions

| Test | What fails |
|---|---|
| A passive basis/gauge change alters poles after all sources are transformed | The implementation or proposed interpretation. |
| A free-circle tower requires mode-dependent R, M5 or alpha | This common-action ansatz; do not repair with a separate geometry per particle. |
| Claimed neighboring signed branches violate the common second-difference law | The free-circle mass assignment at its declared accuracy. |
| A pole disappears from both fixed probes despite their nonzero sum rule | The G1 result or the proposed disappearance claim. |
| An interaction needs energy creation or forbidden charge/KK-momentum changes | That transition model. |
| Norm positivity or bounded energy is lost | That proposed completion; negative mode number is not a negative-energy solution. |

G1's failure or success is not a proof or disproof of every possible
energy-spectrum model. A more general internal operator is a new
explicit hypothesis and must earn its extra parameters with predictions.
