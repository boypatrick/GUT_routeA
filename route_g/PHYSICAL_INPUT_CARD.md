# Route-G comparison-branch physical input card

Date: 2026-09-23. Updated by **G2-T: separate mass and sign measurements**.

Priority notice: the user subsequently selected the relational-time
mainline G-R. Its assumptions and missing inputs are in
[RELATIONAL_TIME_CARD.md](RELATIONAL_TIME_CARD.md).
This file remains the conditional KK/Higgs/xenon/TOF comparison card.
Its source/readout bounds are retained, not transferred to a different
relational-clock realization or made a prerequisite for it.

## Current measurement decision — G2-T

Different masses and opposite signs are separate tasks. TOF plus an
independent single-particle momentum/energy constraint can determine mass;
neutral X currently has no specified sensor for that constraint.
No charged-track resolution, pair missing momentum or S2 drift time is
substituted for it. A source timestamp/path and joint correlations are
also required. See the [measurement input card](data/TOF_READOUT_CARD.md).

TOF plus recoil energy and direction admits an exact same-action elastic
inversion for a known target at rest, but the actual S2-only response
does not supply those observables. The illustrative angular-only budget
of0.896 arcsec for1% mass uncertainty is a requirement, not calibration.

Signed-mode inference still needs a measured Phi-source correlation or
a genuinely signed response. No new interaction was activated. Purely
odd amplitudes remain rate-blind without interference or another signed
observable. Direct Higgs production retains its equal-sign likelihoods.
Any added time gate reduces or preserves signal count; the source
decision below is unchanged. The mathematical audit is bounded-done.

## Current source decision — G2-H

The first concrete source is on-shell Higgs production at one LHC
Run-2 interaction point: actual luminosity139 fb^-1, SM production55.6 pb,
and the same conditional invisible allowance. See the
[external source card](data/HIGGS_SOURCE_CARD.md). It supplies an upper
budget of1.6539 million hypothetical X plus anti-X, not an observed
dark beam. This drive-off decay preparation does not implement the old
Phi-X collision or its external pump.

The integrated source is insufficient even under perfect interception:
the maximum-chord, efficiency-one bound is fewer than3.96e-6 retained
elastic NR events across the whole declared scale domain. Do not reuse
the XENON kg-day exposure to count time twice, apply the volume-average
matrix to an arbitrary narrow beam, or invent a live-time flux.

Direct Higgs production has equal j=+-1 priors and identical time,
direction and momentum likelihoods at the declared tree level.
Alternative four-dimensional kinematic readout alone cannot label
this exact sign degeneracy; its optimal error is50%. This is a new
source-specific result, not a revision of the old collision priors.
The [source/readout test](output/g2_higgs_source.md) is bounded-done
and this realization is rejected, while off-shell/other-source and
other-response hypotheses remain untested.

## Current selected branch — G2-X

The user authorized choosing a scale range and target. We now select
the **X-only uniform Higgs-portal candidate**, as a conditional tree EFT,
and **natural liquid xenon** with the official XENON1T S2-only NR response.
The range Lambda=1–30 GeV spans the externally known Higgs pair threshold
Lambda=mh/10=12.508 GeV; it is a motivated search interval, not a measured
radius or a fit to known particle masses. Inputs mh=125.08 GeV, v=246 GeV
and SM Higgs width4.10 MeV come from the cited PDG review. The observed
ATLAS2023 B_inv<.107 bound is used only under its stated SM-production
and invisible-final-state assumptions.

Real response data are at [the pinned release](data/xenon1t_s2only/PROVENANCE.md):
raw exposure356770 kg day, true recoil0.7–50 keV, complete S2 bins
150.027–3000 PE, with selection losses already included. This is the
2019 analysis/2020 response release, not claimed to be the newest result.
It calibrates target nuclear recoil, **not** the old direct momentum
estimator. The response is folded rather than replaced by a new Gaussian.

The resulting [feasibility calculation](output/g2_xenon_feasibility.md)
reports required **at-detector selected-particle flux for three expected
accepted recoils**, not a measured flux or discovery sensitivity. Source
production, transport, backgrounds and a realized GeV-frequency drive
are still missing. Within the retained window and scalar nuclear kernel,
the S2-only signal shapes for j=+1 and j=-1 coincide, so this readout
does not implement the old sign-classification contract. The best
S2-only error is 14.057%, a majority guess, not the old 2.20% momentum
readout. At Lambda=1, 5, 30 GeV the required flux at the conditional
coupling ceiling is 5.07e9, 2.15e10, 4.32e6 cm⁻²s⁻¹ respectively.
Those ceilings and event requirements are not evidence of an attainable
source. See the [feasibility figure](output/figures/g2_xenon_feasibility.png).

The following G2-P intake sections are retained to document the original
requirements. Items above are now provisionally specified; the new EFT
choice is not proof of a five-dimensional SM embedding or hardware.
The synthetic tests in `output/g2_scale_anchor.*` and
`output/g2_probe_coupling.*` must never be loaded as measured calibration.

## 1. Select the intended physical claim

- **New neutral particle sector:** a candidate SM Higgs-density coupling
  is derived in `tex/route_g_physical_anchor.tex`; it is not yet adopted.
  X is a recoil partner, not an instrument. The existing external U(1)
  has not been identified with electromagnetism.
- **Existing experimental system:** name the particles or excitations,
  apparatus, interaction and independently measured observables.
- **Analogue mode-conversion system:** give the platform's Hamiltonian
  and unit map. This tests a spectral mechanism, not extra dimensions
  or elementary-particle mass conversion.

## 2. Supply a physical energy anchor

Required: observable definition, numerical value, unit, uncertainty,
source/provenance, and its relation to a model mode or pump frequency.
Do not identify a peak by mass proximity alone.

- With the frozen dimensionless mass ratios, a single identified mass
  sets Lambda=m_n/mhat_n, but does not validate those ratios.
- Without freezing those ratios, three consecutive **signed** mode
  masses determine the quadratic spectrum; a fourth is an out-of-fit
  check. Supply the full mass covariance, label uncertainty and theory
  uncertainty. Consecutive sorted masses are not sufficient.
- If an actual pump frequency is used, explicitly impose and justify
  Omega/Lambda=.2; its mere numerical appearance in the old card is
  not a measured energy scale.

Once anchored: sigma_i=.1 Lambda, sigma_d=.05 Lambda, decision threshold
c=1.25892910 Lambda. The old 2.20225% result is a conditional forecast,
not measured performance. Its minority-sign error is 6.23%.

## 3. Supply a coupling and target, not just an error width

For the optional uniform Higgs-density candidate:

- Declare lambda_X and lambda_Phi, matching scale, EFT cutoff and
  post-EWSB physical masses. Prefer studying lambda_X alone first at
  a declared matching scale; lambda_Phi=0 is not loop-protected.
- Include common mass shifts lambda_S v²/2. If both couplings are
  nonzero, recompute the extra Higgs-mediated elastic Phi-X amplitude.
- Keeping the old physical MD²=25 Lambda² with no X self-quartic requires
  lambda_X<=50 Lambda²/v². Otherwise the pre-EWSB X mass is negative
  and the retained action runs away along H=Phi=0. A stabilizing new
  self-quartic is a model change, not an unnoticed mass subtraction.
- Give a target species, column density/geometry and scattering kernel;
  a free-nucleon formula is not a nuclear or condensed-matter response.
- An applicable Higgs width budget and the full tower of accessible
  modes give a conditional upper bound on scattering strength. Check
  escape/decay/acceptance assumptions before calling a channel invisible.
- No Galactic density or dark-matter abundance is assumed for this beam.

A brane-localized coupling is a different extension: it mixes KK modes
and breaks the apparatus's compact-translation symmetry. Do not preserve
the old full-system selection rule by fiat.

## 4. Actual response/calibration fields

Provide the following, with data provenance and statistical/systematic
covariance, if available:

| Input | Needed definition | Current status |
|---|---|---|
| Incident ensemble | Momentum/energy distribution, frame, angular spread, species preparation | Only synthetic COM beam |
| Output labels | Mass identification response, confusion and acceptance for m and abs(j) | Idealized only |
| Readout | Measured variable (track, recoil deposit, timing, etc.), units, conditional response to true state | G2-X uses real XENON1T NR-to-S2 response; old momentum Gaussian not calibrated |
| Calibration | Gain, offset, non-Gaussian tails, drift, correlation with labels and efficiency | Stress tolerances, not data |
| Count normalization | Exposure, incident flux, target density/geometry or collision luminosity | XENON1T exposure supplied; selected-particle flux missing |
| Pump | Calibrated frequency, modulation amplitude, actual energy coupling | Prescribed reservoir only |
| Optional work tag | Full response T(z|q,p,...) and correlations with recoil readout | Hypothetical ternary matrix only |
| Decision goal | Overall error or per-class recall, backgrounds and cost of rejection | Old overall 5% engineering convention |

A deposited recoil energy is not incident momentum: integrate over the
unknown scattering angle and target response. Do not replace the
required response by a relabeled Gaussian width. Do not infer a work tag
from the same recoil datum and count it as new independent information.

## Bounded next action

Do not refine transport or resolution for the rejected on-shell
Higgs-to-xenon rescattering chain. Prefer a collider production/missing-
momentum feasibility test that does not pay a second small interaction
probability. Such an inclusive signal alone does not resolve a KK tower.
If signed-mode conversion is indispensable, specify an actual
sign-sensitive preparation or interaction before building a timing or
momentum instrument. No new coupling or source is authorized by a
change of observable alone.
