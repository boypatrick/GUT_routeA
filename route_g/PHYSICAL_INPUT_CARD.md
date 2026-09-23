# G2-P physical input card — not calibration data

Date: 2026-09-23. Status: **unfilled physical inputs**.
No eV/GeV scale, observed particle identification, measured sensor response,
target material or experimental performance has been assigned.
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
| Readout | Measured variable (track, recoil deposit, timing, etc.), units, conditional response to true state | Gaussian ansatz only |
| Calibration | Gain, offset, non-Gaussian tails, drift, correlation with labels and efficiency | Stress tolerances, not data |
| Count normalization | Exposure, incident flux, target density/geometry or collision luminosity | Missing |
| Pump | Calibrated frequency, modulation amplitude, actual energy coupling | Prescribed reservoir only |
| Optional work tag | Full response T(z|q,p,...) and correlations with recoil readout | Hypothetical ternary matrix only |
| Decision goal | Overall error or per-class recall, backgrounds and cost of rejection | Old overall 5% engineering convention |

A deposited recoil energy is not incident momentum: integrate over the
unknown scattering angle and target response. Do not replace the
required response by a relabeled Gaussian width. Do not infer a work tag
from the same recoil datum and count it as new independent information.

## Bounded next action

Prefer a coupling/constraint feasibility check over more resolution scans.
Choose one physical sector and an externally justified scale, or provide
real response data. Evaluate whether the allowed interaction strength and
target can supply any useful events before claiming the old sign classifier
is realizable. Only then translate the conditional Gaussian requirement
into an apparatus-specific reconstruction requirement.
