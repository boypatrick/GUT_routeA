# Route G Roadmap — Scheme A

Created: 2026-09-22
G1 verified and PDF visually reviewed: 2026-09-23
Current scope: a separate, bounded common-spectrum prototype.
Status labels: open, in-progress, bounded-done, rejected-for-this-model.

## Authority and research discipline

The user explicitly chose Scheme A and requested a new Route G.
This is not the historical Route A and does not replace Route F.
Keep all input assumptions and claimed interpretations visible.
Do not import arbitrary particle-specific radii, projectors or mass shifts.
Do not turn every unfinished UV detail into another prerequisite paper.
No existing P54 parameters, fits or regression results are changed.

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
