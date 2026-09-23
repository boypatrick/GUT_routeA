# Route G — Scheme A: common spectrum and internal modes

Created: 2026-09-22. Separate research branch, explicitly requested by the user.

Route G 回到「同一內部結構，產生不同能譜與可見性」的構想。它不把 Spin(10) 的表示維度
當成額外時空，也不把轉換座標誤認為轉換粒子質量。
Route F 的四維 P54 模型與先前計算保留不動；Route G 尚未導出或取代它。

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
~~~

G1 uses NumPy; G2 uses NumPy and SciPy.
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
Next choose the physical sector/scale and test its coupling against
constraints, or provide real calibration. Specified finite-encounter
packets and a pump envelope are needed for absolute probabilities or
finite-pulse interference, not to reopen the completed conditional
rate calculation. Route F remains independent.
