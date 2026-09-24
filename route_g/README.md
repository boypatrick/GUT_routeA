# Route G — Relational time and the common-spectrum comparison

Created: 2026-09-22. Separate research branch, explicitly requested by the user.

Route G 回到「同一內部結構，產生不同能譜與可見性」的構想。它不把 Spin(10) 的表示維度
當成額外時空，也不把轉換座標誤認為轉換粒子質量。
Route F 的四維 P54 模型與先前計算保留不動；Route G 尚未導出或取代它。

## Current mainline — G-R4 physical source and universality (2026-09-24)

G-R4 is **bounded-done as a physical-source audit**. Select one isotropic
blackbody bath and the E2/E3 clock transitions of the same 171Yb+ ion.
The interaction -d.E and independently measured polarizabilities fix
both shift signs; no desired slowdown is fitted. Conventional spatial
radiation density is an explicit control, not an identification of Chen
Ed or G-R3's energy-per-cell candidate.

At 300 K, leading static electric-dipole shifts are -0.3604 Hz and
-0.04638 Hz, or fractions -5.235e-16 and -7.223e-17. Both go downward,
but the fractions differ by about7.25. This thermal channel is therefore
not solely a universal clock-rate rescaling. **An additional exactly
common factor is invisible in the same-bath frequency ratio** and is
not excluded. Equal fractions would be necessary, not sufficient, for
the sole-universal-factor interpretation.

The two-setpoint 300 -> 310 K ratio modulation is predicted to be
+6.325e-17 in the retained static approximation; no such new observations
are analyzed. Complete dynamic response, actual radiometry and correlated
apparatus systematics are not claimed. The numerical stress box is
conditional, not an empirical confidence statement.

Read the [five-page PDF](output/pdf/route_g_bbr_universality.pdf),
[full TeX](tex/route_g_bbr_universality.tex),
[pinned inputs](data/BBR_CLOCK_INPUTS.json),
[numerical report](output/gr4_bbr_universality.md) and [roadmap](ROADMAP.md).
Run: python3 route_g/code/verify_gr4_bbr_universality.py.
New checks91/91; complete Route-G regression4728/4728.
PDF compiled twice and visually reviewed page by page.

Next useful work is raw interleaved E2/E3 data plus local radiation
calibration, not more time-law fitting. To study a common factor itself,
first specify an independent reference or non-clock observable; ratio
data alone cannot identify it. No Route-F gate is changed.

## G-R3 foundation — density and relative clocks (2026-09-24)

G-R3 is **bounded-done for one explicit candidate**. Define
Ed=epsilon_D N_D per addressed matter cell, with ground-subtracted
excitation energy; this is not an assumed spatial energy density or
the full information in a state. Couple B through V=kappa H_B Ed/E0,
while reference A is shielded. A common finite stationary constraint
gives normalized joint phase likelihoods and r_B/A=1+kappa d.
Calibration stays fixed; no slowdown function is fitted.

Two useful findings:

- The interaction sets the sign. A microscopic dispersive exchange
  model permits both faster and slower clocks, depending on detuning.
  It does not establish that higher density universally slows time.
- Unread density fluctuations produce a mixture of clock phases,
  not generally the phase at mean Ed. The exact test gives P0=5/9,
  whereas mean-density substitution gives (3+2sqrt(2))/9.

This is an environment-dependent clock-frequency shift, not yet
universal proper time. The single-cell definition does not distinguish
density from total excitation energy. Reference preparation, finite
commensurate spectra and measurement calibration remain assumptions.
The G-R2 record apparatus is not silently promoted to an energy-complete
implementation of the new stationary constraint.

Read the [density/convention card](DENSITY_CLOCK_CARD.md),
[six-page PDF](output/pdf/route_g_density_clocks.pdf),
[TeX](tex/route_g_density_clocks.tex) and
[numerical report](output/gr3_density_clocks.md).
Run: python3 route_g/code/verify_gr3_density_clocks.py.
New checks273/273; complete Route-G regression4636/4636.
The PDF is compiled twice and each page visually reviewed.
The subsequently authorized G-R4 physical-source/cross-transition audit
is recorded above; it does not silently replace this candidate's Ed.

## G-R2 foundation — saved event records (2026-09-24)

G-R2 is **bounded-done**: one initial record and one terminal readout,
with exact joint probabilities and explicit measurement backreaction.
The cyclic dial 0,1,2,0,1,2 is distinct from the stored-event count
0,1,1,2,2,2. Fixed event stamps and a one-bit epoch label distinguish
the completed episode from the initially blank one.

The exact terminal P(+) is 0 without the initial record and 2/9 with it.
Another instrument with the same initial effects gives 4/9: an event
effect alone does not specify what a record does to the system.
Pointer probabilities survive all declared links after writing.
Full-state recurrence cannot simultaneously return to blank memory
and retain a fresh record; this is not a no-go theorem for cyclic clocks.

The construction explicitly adds a finite open-history constraint,
ordered gates and an input boundary. It is NOT the unchanged G-R1
noninteracting energy constraint. Stamps are fixed by this scheduled
protocol, not measured arbitrary timestamps. No physical apparatus
energy budget, irreversible arrow, Ed interaction or clock-rate law
is claimed; the bare-energy-noncommuting terminal gate is flagged.

Read the [five-page derivation PDF](output/pdf/route_g_two_event_records.pdf),
[TeX](tex/route_g_two_event_records.tex),
[numerical report](output/gr2_two_event_records.md) and
[roadmap](ROADMAP.md).
Run: python3 route_g/code/verify_gr2_two_event_records.py.
New checks139/139; all4363 Route-G checks pass.
The PDF was compiled twice and every page visually inspected.
The subsequent G-R3 candidate and relative-clock test is now recorded
above; G-R2 itself did not define or fit a density-rate law.

## G-R1 foundation — relational quantum clock (2026-09-23)

The user has chosen **relational/emergent time**, closer to the original
document's Ei/Ed and time-event discussion. The previous G0–G2-T
KK/TOF work is retained as a comparison program, not a derivation of that
time concept and not a prerequisite for the new line.

The first bounded stage is complete. A finite Page-Wootters construction
has positive local clock/system energies and one stationary joint state.
Conditioning on the clock's event gives exact subsystem Schrödinger
evolution. Clock phase is a measurement label, not a geometric rotation
or a new spatial dimension; conversion to seconds still needs calibration.
Quantum mechanics, the subsystem split and spectra are assumed, not
derived from the original paper.

**Key finding:** identical energy statistics do not determine a relational
history. In the explicit three-level pair, coherent and energy-dephased
states have the same local marginals and joint energies, but
P(+|phi)=(1+2cos(phi))²/9 versus1/3. Ordinary energy density alone therefore
cannot be the entire time-defining information. This does not establish
or universally refute the source's still-undefined Ed law.

G-R1 also fixes the finite-clock/POVM normalization and constraint-invariant
observables. It does not yet construct sequential measurements, accumulated
event records, an arrow of time, emergent space, mass or gravity.
No new detector calculation or freely chosen slowdown law is introduced.

The subsequent **G-R2** record task is now completed above within an
explicit finite protocol extension. G-R3 above is a separate density model.

Read [the interpretation card](RELATIONAL_TIME_CARD.md),
[derivation PDF](output/pdf/route_g_relational_clock.pdf),
[TeX](tex/route_g_relational_clock.tex) and
[numerical report](output/gr1_relational_clock.md).
Run: python3 route_g/code/verify_gr1_relational_clock.py.
New algebra checks109/109 pass; all4115 comparison checks also pass,
total4224/4224. The five-page PDF was compiled twice and each page
visually checked using the PDF workflow, with no final layout/reference
warnings. No empirical validation is claimed.

## Preserved comparison results

The sections below record the previous program. Their conditional
conclusions remain valid, but their next-step suggestions are scoped
to that program and do not override the relational-time priority.

## G2-T — Separate mass metrology from signed-mode readout (2026-09-23)

Bounded-done: the mathematical measurement audit is complete; a physical
instrument is not. The action, earlier source budget and response remain
unchanged. No new Gaussian scan or sign-sensitive coupling was activated.

- **Different masses:** TOF plus independently constrained momentum gives
  m²=p²[(cT/L)²−1]. The full (p,T,L) covariance, source start time and
  exact finite-error boxes are derived. Neutral X has no ordinary
  magnetic-curvature momentum readout; the formula alone does not
  supply the missing measurement.
- **A same-action alternative:** TOF plus recoil energy AND direction
  algebraically determines mass for a known target at rest. The existing
  S2-only response supplies neither the recoil vector nor tagged X TOF.
  At one illustrative 10 keV xenon recoil, an angular-only 1% mass budget
  would require about0.896 arcsec, with all other inputs exact.
  This is a severe conditioning requirement, not achieved hardware.
- **Same mass, opposite signs:** retain G2-R's conditional Phi-source
  correlation, not an assumed preparation. A purely sign-odd amplitude
  is still rate-blind: an even/odd interference term in the same channel,
  or another calibrated signed response, is needed. The direct Higgs
  source remains at50% optimal sign error.

The at-rest illustration Lambda=5 GeV, E=62.54 GeV, L=10 m gives
139.214 ps separation between j=0 and abs(j)=1, not between +1 and −1.
This is not a real LHC beam or a timing calibration.
Any additional gate can only reduce the previous3.96e-6 count ceiling.
Do not optimize the rejected on-shell-Higgs/xenon chain.

Next, require an actual source timestamp and momentum/energy measurement
mechanism before an instrument forecast. If signs matter, prioritize a
physically tagged Phi-source correlation; a new interaction is a separate
explicit proposal. Production observables without rare rescattering
remain the better lead for the rejected source.

Artifacts: [measurement input card](data/TOF_READOUT_CARD.md),
[numerical report](output/g2_tof_mode_audit.md),
[TeX](tex/route_g_tof_mode_audit.tex),
[PDF](output/pdf/route_g_tof_mode_audit.pdf).
New checks97/97; complete Route-G regression4115/4115.
The five-page PDF was compiled twice and every page visually checked
under the PDF workflow, with no final layout or reference warnings.

## G2-H — Actual accelerator source budget: reject this rescattering chain (2026-09-23)

Bounded-done. Use ordinary on-shell LHC Higgs production as a **separate,
drive-off source** for the same X-only portal, not the original prepared
Phi-X mode-conversion collision. One published Run-2 interaction-point
dataset has 139 fb^-1; the SM production input55.6 pb gives7.7284 million
Higgs bosons. Under the retained invisible allowance, at most1.6539 million
X plus anti-X can be produced. This is an integrated source budget,
not a measured dark flux or a built LHC-to-XENON beam.

The same coupling determines production AND nuclear scattering. Even
assigning every particle the maximum xenon chord, arbitrary favorable
Higgs boosts, and efficiency one, a continuous bound over Lambda=1–30 GeV
gives **fewer than3.96e-6 events** in the retained elastic NR window.
A generous input stress gives5.35e-6. Actual geometry, attenuation,
detector live time and response losses cannot improve that single-pass,
nonmultiplying bound. No detailed beam Monte Carlo is justified here.
The high-scale part of the earlier map is not a rescue: on-shell Higgs
production itself closes at Lambda=12.508 GeV (12.2651 for abs(j)=1).

The source also changes the mode-inference question. Direct scalar Higgs
decay makes j=+1 and -1 equally likely with identical full four-dimensional
momentum distributions. The diagonal probe is sign-blind. **Timing,
direction and independent momentum readout all remain at50% optimal
sign error for this source**, not the old asymmetric source's14.057%.
Different abs(j) masses and opposite signs of one degenerate mass must
not be conflated. No statement of all-orders full-theory symmetry is made.

Reject this source/readout realization, not Route G as a whole. Next
consider production/missing-momentum inference without a second tiny
scattering probability; if signed conversion is essential, first require
a genuinely sign-sensitive preparation/coupling. Off-shell production,
other response channels and new interactions are not calculated or
silently ruled out. No additional Gaussian scan or GeV pump was added.

Artifacts: [source card](data/HIGGS_SOURCE_CARD.md),
[numerical report](output/g2_higgs_source.md),
[TeX](tex/route_g_higgs_source.tex),
[PDF](output/pdf/route_g_higgs_source.pdf).
New checks389/389; complete Route-G regression4018/4018.

## G2-X — Real xenon response and conditional physical feasibility (2026-09-23)

The selected conditional branch is an X-only, KK-diagonal Higgs portal,
Lambda=1–30 GeV, and natural liquid xenon. The range crosses the measured
Higgs mass's pair threshold Lambda=mh/10=12.508 GeV; it does not determine
the radius or identify a known particle with X. The old collision kernels
remain unchanged after an explicitly declared tree mass matching.

We use the official **XENON1T S2-only 2019 analysis/2020 response release**,
not a fabricated efficiency or the latest halo-WIMP exclusion curve.
The pinned NR/S2 matrices include all selection losses; the exposure is
356770 kg day. True NR energies are limited to0.7–50 keV, and complete
observed bins to150.027–3000 PE. We do not extrapolate response outside
that interval or multiply a second fiducial efficiency.

The map combines the conditional ATLAS B_inv<.107 bound on the full
open complex X tower, the retained-action vacuum condition, and the
required at-detector flux for **three expected accepted recoils**.
That count target is not a discovery significance or a confidence limit.
No physical incoming flux has been supplied, so no actual detectable
region or xenon coupling exclusion is invented.

**New physical finding:** below every recoil endpoint, the scalar kernel
factorizes as d_sigma/d_ER=(lambda_X²/p_X²) H(ER). Folding the real
energy-only S2 response preserves identical normalized shapes for the
two recoil signs. It changes their event priors but supplies no sign
information in this retained window. Accepted priors are 85.943%/14.057%,
and the best S2-only error is **14.057%**, merely a majority guess.
The old2.20% direct-momentum error
cannot be transferred to this readout by improving an S2 resolution.

The response-weighted cross section is 7.99989e-39 cm² times
lambda_X²/(Lambda/GeV)². At Lambda=1, 5, 30 GeV the conditional maximum
couplings require respectively 5.07e9, 2.15e10, 4.32e6 cm⁻²s⁻¹ incident
flux for three accepted events. These are requirements, not achievable
beams. The fN-only flux variation is about −11%/+13%, not a total error.
The minimum recoil endpoint is 9.543 MeV, proving the entire 0.7–50 keV
component is below every endpoint across the continuous scale range.

There is also no realized source: the retained drive now corresponds to
Omega=.2–6 GeV, approximately4.84e22–1.45e24 Hz. Its prescribed reservoir
is not a specified laboratory oscillator. If the goal is event counting,
next supply production/transport and actual at-target flux; if the goal
is sign identification, change the observable to independently measured
momentum/timing, direction or a reachable recoil endpoint.

Artifacts: [feasibility figure](output/figures/g2_xenon_feasibility.png),
[numerical report](output/g2_xenon_feasibility.md),
[TeX](tex/route_g_xenon_feasibility.tex),
[PDF](output/pdf/route_g_xenon_feasibility.pdf), and
[data provenance](data/xenon1t_s2only/PROVENANCE.md).
The new verifier passes 573/573 checks; the complete Route-G regression
passes **3629/3629**. This milestone is bounded-done as a conditional
physical feasibility test, not a validated particle or detector design.

## G2-P — Physical scale and visible-probe audit (2026-09-23)

The dimensionless readout forecast does not determine an energy unit.
With Lambda=1/R, all energies and resolutions co-scale, K scales as
Lambda^-2, and the conditional error stays 2.20225%. Synthetic rescalings
(.1,1,10) pass 877/877 checks; these are not proposed GeV scales.
Three calibrated, consecutive **signed** mode masses determine the
quadratic spectrum. A fourth must satisfy the out-of-fit test
`m_2² - 3 m_1² + 3 m_0² - m_-1² = 0` at declared tree accuracy.
A universal mass² shift cannot repair wrong mode spacings.

We derive, but do not activate, a candidate circle-uniform Higgs-density
coupling to the neutral scalar towers. It is a four-dimensional EFT,
not a completed five-dimensional SM embedding or a Route-F portal.
It preserves KK labels but shifts the common masses. Coupling BOTH
towers also adds Higgs-mediated elastic Phi-X scattering, changing the
old inclusive fractions. A detector-facing X-only coupling is the
minimal matching-scale option, not an all-orders protected zero for Phi.

This produces a real physical tradeoff: the same coupling governs
nucleon scattering and the summed open-tower Higgs decay width.
Keeping the old X mass with no added X self-quartic also requires
`lambda_X <= 50 Lambda²/v²` to avoid a negative underlying X mass and
runaway. Increasing a nominal probe strength is therefore not a free
way to make the toy recoil partner measurable.

No actual energy anchor, sensor calibration or target has been supplied.
A single recoil-energy deposit is not a momentum measurement; the
scattering-angle and target response must be included. See the
[physical input card](PHYSICAL_INPUT_CARD.md),
[derivation](tex/route_g_physical_anchor.tex),
[PDF](output/pdf/route_g_physical_anchor.pdf),
[scale tests](output/g2_scale_anchor.md) and
[probe-normalization tests](output/g2_probe_coupling.md).
Next choose a physical sector/scale for a coupling-versus-constraint
feasibility test, or supply real apparatus response data. This bounded
audit is not hardware realization, a mass fit, or further Gaussian tuning.
The scale and candidate-probe verifiers pass 877/877 and 113/113 checks;
all 2066 existing checks also pass. The four-page PDF was visually reviewed.

## G2-M — Selected readout contract (2026-09-23)

Choose momentum-only readout, keeping the drive on, sigma_i=.1 and
nominal sigma_d=.05. Freeze the rule: predict j=+1 for reconstructed
x<1.25892910, otherwise j=-1. Nominal overall error is 2.20225%.
The class-conditional errors are 1.0315% and 6.2279%: overall <=5%
does NOT imply at least 95% recall for each sign.

For true x=(1+a)p_f+b+Gaussian(0,sigma^2), stress
|a|<=.01, |b|<=.02 and .04<=sigma<=.06 without retraining.
The 27-point sampled maximum is 3.8648%; an eight-cell continuous-box
CDF envelope is about 4.46%. Its inequality covers the nuisance box
with convergence-tested integration, not interval-certified arithmetic.
The frozen-rule verifier passes 154/154 checks.

The noisy-tag alternative passes 125/125 checks. At sigma_d=.1 the
declared symmetric ternary T_eta=(1-eta)I+eta U needs eta about .1765
or less. This is NOT a generic "88% tag accuracy" requirement:
always reporting q=0 is 88.98% accurate here but adds zero information,
leaving recoil-sign error at 5.94485%. The full confusion matrix matters.

This closes the engineering readout calculation, not hardware design.
X is a recoil field, not a completed instrument. The physical energy
scale, sensor coupling, mass-label response and acceptance are missing.
Next choose a physical realization or provide calibration data,
not indefinitely finer Gaussian scans.

## G2-S/D — Resolution and an activated powered background (2026-09-23)

The user has now activated the time-dependent option. The undriven
G2/G2-R results remain intact as the reference; the separate driven
branch explicitly changes the quartic to g(t)=g[1+epsilon cos(Omega t)].
The engineering pump card is epsilon=.5, Omega=.2. It exchanges energy
q Omega, not compact momentum, and leaves the free mass tower unchanged.
This is a prescribed external reservoir, not a derived closed-system pump.

G2-S adds an incoherent, p>=0 truncated-normal incoming distribution
around p0=.2 and a Gaussian reconstructed-momentum readout. All rates
are integrated before normalizing the selected m=0, |j|=1 sample.
The 50-card scan uses a declared 5% optimal sign-error criterion,
not an experimental confidence level or an assumed instrument.

| Parent beam width sigma_i | Readout width sigma_d | No drive | Drive, q unknown | Drive, ideal q tag |
|---:|---:|---:|---:|---:|
| .1 | .1 | 3.77384% | 5.94485% | 3.77589% |
| .2 | .1 | 4.84993% | 6.83823% | 4.85210% |
| .4 | .1 | 11.0294% | 12.2710% | 11.0331% |

The drive creates six selected recoil lines. Two opposite-sign
components lie at 1.20257369 and 1.18924247, only .01333122 apart.
Thus adding power does not automatically improve identifiability.
These lines are distinct at zero noise; finite-resolution inference
uses their full weights, not a smallest-gap rule alone.
An ideal q tag is an information benchmark, not a constructed detector.
At sigma_i=.1, reducing sigma_d to .05 gives 2.20225% driven untagged error.

The sharp-beam driven ledger retains all 66 open (l,m,q) channels.
Total K=2.10563491e-6 and mean matter-absorbed work per scattering is
.00900097215; without luminosity neither is an absolute probability
or a pump power. G2-S passes 272/272 checks and G2-D passes 565/565.
These certify the implementation within its tree/long-time assumptions.
No four-dimensional collision packet, finite-pulse prediction or
experimental mass fit has been supplied.

## G2-R — Joint output and recoil (2026-09-23)

Bounded-done at tree level: 357/357 reproducible checks pass, including
regression against the unchanged G2 rates. The four-page derivation PDF
has been compiled and visually reviewed.

The same action, parameter card and detector preparation now give a
joint conditional event distribution for the output Phi mode m and
detector recoil j. No interaction or fitted parameter is added:

\[
P(m,j\mid{\rm scattering})
=\frac{w_{m+j-n}K_{m+j-n,m}}{\sum_{l,m'}w_lK_{lm'}}.
\]

The report separates all scattering events (including elastic collisions)
from mode-changing events. It includes full joint matrices, marginals,
conditional distributions, recoil momenta/energies and the ideal
isotropic two-body angular law. It is not a probability of an arbitrary
incoming packet having an encounter.

The physically important distinction is the readout: a detector mass
measurement gives |j|, not the sign of j. On the alpha=.25 card, adding
outgoing COM momentum distinguishes the two signs in nonelastic
channels, provided the incoming momentum and model parameters are known.
For m=0, |j|=1, the two predicted momentum lines are 1.02175458 and
1.36191317. Elastic events retain the sign ambiguity.

This is model-dependent kinematic reconstruction, not an independent
measurement of compact momentum or an entanglement/localization witness.
Finite beam width, detector resolution and efficiencies are not assigned
fictitious values. A calibrated readout can be applied next; absolute
encounter probabilities still require four-dimensional packets.

## G2 — Local dynamical conversion (2026-09-23)

G2 adds one neutral complex detector field X and the positive local
potential kappa5 |Phi|² |X|². The detector is localized by a quantum
wavepacket, not a pinned external defect. Its recoil makes the process

\[
\phi_n+X_l\longrightarrow\phi_m+X_{l+n-m}
\]

possible, with one tree amplitude M=-g, g=kappa5/(2 pi R).
Total four-dimensional energy/momentum, compact momentum and both
charges are conserved. The outgoing excitation occupies another existing
mass branch; no coordinate transformation changes a mass eigenvalue.
This action is an explicit minimal candidate, not uniquely derived from
the motivating paper.

On the unchanged G1 quarter-holonomy card (R=1, M5=0.5, alpha=0.25),
add MD=5, g=0.05 and incoming COM momentum 0.2. The component
Phi_1 + X_0 -> Phi_0 + X_1 has outgoing momentum 1.02175458 and
tree cross section 6.27164081e-6 in inverse-square common energy units.
The Phi rest-energy decrease 0.78727421 pays for detector internal
recoil 0.09901951; the remaining 0.68825469 becomes ordinary kinetic
energy. No external power source is present.

The finite detector preparation is localized around theta=0.7 with
coefficients c_l proportional to exp[-l²/(4*1.2²)] exp[-i l theta0],
l=-4,...,4. Different components have different collision energies.
Tracing over recoil removes cross terms for a fixed incoming/outgoing
Phi mode. A common-overlap averaged sigma*v is a rate coefficient,
not an absolute encounter probability; the latter still needs specified
four-dimensional envelopes and geometry.

The result is **tree-level EFT only**. The classical vacuum is stable
and the G1 quadratic tower is unchanged, but interacting loop masses,
counterterms and higher operators have not been calculated. G1's exact
free response is not promoted to an all-orders interacting identity.
No measured particle data were fitted.
The G2 verifier passes **262/262 implementation checks**; the G1
regression still passes 331/331. The tree conversion milestone is
complete within this scope, not a closure of the full physical theory.

## First bounded milestone: G1

The initial prototype is a positive-energy complex scalar on fixed
four-dimensional Minkowski spacetime times one circle. One shared
compactification radius R, bulk mass M5 and background holonomy alpha
determine the entire tower:

\[
m_n^2=M_5^2+(n+\alpha)^2/R^2,\qquad n\in\mathbb Z.
\]

Two probes at opposite points of the same circle have fixed interactions,
including the parallel transporter between them. They see the SAME poles
with different residues. Their normalized relative weights are

\[
p_{\pm,n}=\frac{1\pm\cos[\pi(n+\alpha)]}{2},\qquad
p_{+,n}+p_{-,n}=1.
\]

These are probe weights, not decay probabilities. At a degenerate mass,
residues are summed over the complete eigenspace.
In particular, at alpha=0 a degenerate sine standing wave can be dark
to BOTH probes: the nonzero branch sum does not ensure visibility of
every superposition. The mass pole remains visible through the other
combination.

**The useful distinction:** changing a basis leaves masses and the
correctly transformed response invariant; changing a physical holonomy
or radius changes the tower. A mode dark to one probe still exists.
No measured particle masses are used to select the example parameters.

## Files and reproduction

- [Research plan and rejection gates](ROADMAP.md)
- [G2-X feasibility calculation](output/g2_xenon_feasibility.md), [code](code/verify_g2_xenon_feasibility.py), [figure](output/figures/g2_xenon_feasibility.png), [TeX](tex/route_g_xenon_feasibility.tex) and [PDF](output/pdf/route_g_xenon_feasibility.pdf)
- [G2-P physical input card](PHYSICAL_INPUT_CARD.md), [TeX](tex/route_g_physical_anchor.tex) and [PDF](output/pdf/route_g_physical_anchor.pdf)
- [Scale-anchor report](output/g2_scale_anchor.md) and [code](code/verify_g2_scale_anchor.py)
- [Candidate probe report](output/g2_probe_coupling.md) and [code](code/verify_g2_probe_coupling.py)
- [G2-M frozen rule and calibration report](output/g2_readout_contract.md), [code](code/verify_g2_readout_contract.py) and [data](output/g2_readout_contract.json)
- [Noisy-tag comparison](output/g2_noisy_tag.md), [code](code/verify_g2_noisy_tag.py) and [data](output/g2_noisy_tag.json)
- [Readout-contract derivation](tex/route_g_readout_contract.tex) and [PDF](output/pdf/route_g_readout_contract.pdf)
- [G2-S finite-resolution report](output/g2_resolution.md), [code](code/verify_g2_resolution.py) and [data](output/g2_resolution.json)
- [G2-D powered-channel report](output/g2_driven.md), [code](code/verify_g2_driven.py) and [data](output/g2_driven.json)
- [G2-S/D derivation](tex/route_g_resolution_drive.tex) and [PDF](output/pdf/route_g_resolution_drive.pdf)
- [G2-R joint prediction report](output/g2_recoil_joint.md)
- [G2-R derivation](tex/route_g_recoil_joint.tex) and [PDF](output/pdf/route_g_recoil_joint.pdf)
- [G2-R executable](code/verify_g2_recoil_joint.py) and [full joint data](output/g2_recoil_joint.json)
- [G2 local action, recoil and packet derivation](tex/route_g_local_conversion.tex)
- [G2 typeset derivation](output/pdf/route_g_local_conversion.pdf)
- [G2 numerical report](output/g2_local_conversion.md)
- [G2 executable](code/verify_g2_local_conversion.py) and [saved tests](output/g2_local_conversion.json)
- [Action and full short derivation](tex/route_g_scheme_a.tex)
- [Typeset derivation](output/pdf/route_g_scheme_a.pdf)
- [Reproducible numerical checks](code/verify_g1_circle_spectrum.py)
- [Numerical report](output/g1_circle_spectrum.md)
- [Saved spectra, residues, tests and source hash](output/g1_circle_spectrum.json)
- [Sources and interpretation boundary](SOURCES.md)

~~~sh
python3 route_g/code/verify_g1_circle_spectrum.py
python3 route_g/code/verify_g2_local_conversion.py
python3 route_g/code/verify_g2_recoil_joint.py
python3 route_g/code/verify_g2_resolution.py
python3 route_g/code/verify_g2_driven.py
python3 route_g/code/verify_g2_readout_contract.py
python3 route_g/code/verify_g2_noisy_tag.py
python3 route_g/code/verify_g2_scale_anchor.py
python3 route_g/code/verify_g2_probe_coupling.py
python3 route_g/code/verify_g2_xenon_feasibility.py
python3 route_g/code/verify_g2_higgs_source.py
python3 route_g/code/plot_g2_xenon_feasibility.py
~~~

G1 uses NumPy; G2 uses NumPy and SciPy; the G2-X figure uses Matplotlib.
From the repository root, after generating the figure, build its note:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=route_g/output/pdf route_g/tex/route_g_xenon_feasibility.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=route_g/output/pdf route_g/tex/route_g_xenon_feasibility.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=route_g/output/pdf route_g/tex/route_g_higgs_source.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=route_g/output/pdf route_g/tex/route_g_higgs_source.tex
~~~

The four engineering cards are R=1, M5=0.5 with alpha=0, 0.25, 0.5,
and R=2, M5=0.5, alpha=0.25, in common arbitrary energy units.
One spectral computation is checked independently against a link-covariant
finite-difference operator and explicit gauge/basis transformations.
This is not a new AP-E soliton lattice scan.
The independently rerun verifier passes **331/331 implementation checks**.
Halving the grid spacing gives low-spectrum error ratios 3.98--3.99,
consistent with the derived second-order convergence; the truncated
responses satisfy their analytic tail bounds. These are code checks,
not experimental evidence for an extra dimension.

## What is and is not being claimed

| Claim | Status |
|---|---|
| One action produces a correlated tower and probe couplings | Analytically derived; numerical verification in G1. |
| Coordinates alone change the particle masses | Rejected by basis covariance. |
| Darkness in one channel means the mode ceases to exist | Rejected for the two-probe prototype. |
| A physical change of alpha or R can change the spectrum | Demonstrated between the prescribed static cards; symmetries can still relate equal spectra. |
| Every state must be visible in at least one of these two probes | False in a degenerate eigenspace: a joint dark standing wave exists at alpha=0. |
| The background dynamically selects alpha and R | Not yet modeled. |
| A collision physically changes the occupied Phi mode | G2 gives an allowed tree-level process with dynamical detector recoil. G1 alone does not. |
| A complete finite-encounter probability has been numerically computed | Not claimed: cross sections and a scoped common-overlap rate are supplied. |
| An actual nuclear-recoil response has been folded | G2-X uses the pinned XENON1T release; source flux and a new signal likelihood are still missing. |
| S2-only readout realizes the old momentum sign classifier | Rejected in the retained scalar-current model/window: normalized shapes coincide; best error is 14.057%. |
| A Run-2 on-shell Higgs source supplies useful xenon rescattering counts | Rejected for the retained elastic single-pass channel: even perfect interception gives fewer than3.96e-6 events. |
| Timing or momentum resolves the direct Higgs source's opposite signs | Rejected at the declared tree level: full four-dimensional likelihoods coincide; best sign error is50%. |
| Electron, muon, quarks or three families are identified with this tower | Not claimed; the present field is scalar, not a chiral fermion. |
| The original paper's Ei/Ep/En have been derived as these fields | Not claimed. |
| Negative energy, emergent time, superluminal travel, black holes or galaxy curves are explained | Not included. |
| This is a complete superstring compactification or a Spin(10) embedding | Not claimed. |

R here is a physical circle radius, not an identification with the symbol
R or the dimensionality discussion in the source paper. Alpha is a
modeling hypothesis motivated by internal phases, not a measured Ei phase.
The ordinary background spacetime is assumed in G1, not derived as emergent.

## Connection to Route F

The previous Route-F low-order flavor/seesaw tension remains unchanged.
Only a later explicit map from shared geometry and normalized modes to
the P54 fields, couplings and thresholds could connect the two programs.
That map would be a new tested hypothesis, not an inference from similar
pictures. No full matching or physical fit is promoted by G1.

G2 implements the dynamical-detector option. Its heavy finite recoil
bypasses the previously closed one-body emission process by supplying
an incoming collision partner, not by violating the old inequality.
G2-R calculates recoil-resolved predictions and mass-only coarse
graining. G2-S/D adds finite-resolution forecasts and a powered branch.
G2-M selects the momentum-only frozen readout and quantifies calibration
margin; noisy sideband tagging remains a conditional alternative.
G2-P adds the scale-identifiability proof and a conditional visible-sector
coupling audit. No energy scale or measured response is fabricated.
G2-X now makes a conditional Higgs/xenon choice, folds a real public
response and computes the required flux under the declared constraints.
G2-H now tests an actual LHC production budget and rejects its xenon
rescattering realization with a boost-independent upper bound.
Next examine production-only/missing-momentum observables, or require
a physical sign-sensitive preparation if signed conversion is essential.
Specified finite-encounter
packets and a pump envelope are needed for absolute probabilities or
finite-pulse interference, not to reopen the completed conditional
rate calculation. Route F remains independent.
