# G2 xenon feasibility: counts per incident flux and a failed S2 sign readout

Verification: **573/573 checks passed**. The matrix is a published experimental response; the source and portal remain a hypothetical candidate. No actual beam flux or signal confidence is claimed.

The chosen physical card is Lambda=1..30GeV, mX=√26Lambda, mh=125.08GeV, v=246GeV, fN=.308±.018. The portal couples X only at matching. Physical MD=5Lambda is held fixed by bare matching; no new X stabilizing quartic is introduced.

## The physical obstruction

The source minimum momentum is 0.809144Lambda. Every exact recoil endpoint in the scan exceeds 9543.0keV, above the entire0.7–50keV window. Consequently dσ/dER=lambda_X² H_A(ER)/pX² throughout the window. The two signs have **identical normalized S2 spectra** after the real response is folded.

Source sign priors: +77.470%, −22.530%. Accepted-event priors: +85.943%, −14.057%. The best S2-only error is **14.057%**, simply guessing the accepted majority sign; no S2-bin classifier improves it in this model/window. The previous2.20% ideal momentum result does not describe a xenon energy-deposit measurement.

## Response and rate normalization

Use complete S2 bins from 150.027 to 3000PE. Nuclear-recoil response rows are piecewise constant over supplied energy-bin edges, clipped to0.7–50keV; the theoretical kernel is integrated within each bin. All full-volume averaged analysis selections are already in the matrix. Do not multiply another fiducial efficiency.

N = Phi_pair ×86400×356770kg-days×(1000 N_A/131.293)×<sigma_acc_cm²>. Phi_pair is the stationary uniform selected-particle flux at the detector in cm⁻²s⁻¹, not source production or a halo prediction. No additional velocity factor is used. Three expected events require Phi_3=3/(N/Phi_pair); this is not a discovery or exclusion threshold.

## Conditional coupling budget and required flux

lambda_max=min(1,50Lambda²/v²,lambda_Higgs when a tower mode is open). The entire signed complex-X tower is counted in the width. The ATLAS B_inv=.107 allowance implies Gamma_X<=B_inv/(1-B_inv)×.00410GeV only under its invisible-channel and SM-width premises. At Lambda>=12.508GeV this particular Higgs-decay constraint disappears; stability and the chosen scan cap remain.

| Lambda [GeV] | mX [GeV] | Conditional lambda_max | Active bound | Minimum flux for3 expected events [cm⁻²s⁻¹] |
|---:|---:|---:|---|---:|
| 1 | 5.099 | 0.00072299 | conditional_invisible_width | 5.0741e+09 |
| 2 | 10.198 | 0.0010333 | conditional_invisible_width | 9.9368e+09 |
| 5 | 25.495 | 0.0017555 | conditional_invisible_width | 2.1516e+10 |
| 10 | 50.990 | 0.0038194 | conditional_invisible_width | 1.8182e+10 |
| 12 | 61.188 | 0.0085654 | conditional_invisible_width | 5.2058e+09 |
| 12.508 | 63.779 | 0.12926 | bare_mass_stability | 2.4834e+07 |
| 15 | 76.485 | 0.1859 | bare_mass_stability | 1.7268e+07 |
| 30 | 152.971 | 0.7436 | bare_mass_stability | 4.317e+06 |

At lambda_X=1,Lambda=1GeV the accepted effective cross section is 7.99989e-39cm²; scale it by lambda_X²/Lambda_GeV². This normalization point itself need not satisfy the physical coupling constraints.

The fN-only yield factors are 0.8865 and 1.1203; reciprocal flux factors are 1.1280 and 0.8926. These are not total uncertainties.

## What is and is not established

The calculation checks coherent scalar-current normalization, the low-speed limit, exact finite-velocity endpoints, isotope-number weighting, Helm F(0)=1, response clipping, independent energy integration, source regression and full signed-tower thresholds. Natural-xenon isotope masses use A×.93149410242GeV; binding/electron corrections and nuclear theory uncertainties are not promoted into a precision model.

Useful event counts still require an independently attainable detector-incident flux, a controlled UV/EFT and an applicable invisible-width interpretation. Public detector calibration does not provide that source. An energy-only ROI below all endpoints is not a momentum spectrometer: recovering sign information needs an independent observable or physically accessible endpoint information, not further fitting of a Gaussian error width.
