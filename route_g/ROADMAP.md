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

Status: open; next bounded physical option after the completed G2 result.

Prioritize one joint observable: outgoing Phi mass branch and detector
recoil, with j-l=n-m and the predicted p_f. This distinguishes actual
transfer from a detector-only visibility change. If an absolute
encounter probability is requested, specify normalized four-dimensional
wavepackets and integrate the existing first-order amplitude.
Do not treat this optional next measurement as reopening the G2 tree
existence calculation. No automatic large matching program is required.

Backup: an explicitly energy-accounted time-dependent background.
It is not needed to obtain the present allowed collision and is not
implemented. G3 and G4 remain independent further hypotheses.

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
