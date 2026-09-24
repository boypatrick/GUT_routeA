# G2-H: an actual accelerator source budget rejects the xenon rescattering chain

Verification: **389/389 checks passed**.

Use one 13 TeV LHC interaction point with the published Run-2 luminosity of 139 fb^-1 and SM Higgs production 55.6 pb. The luminosity is actual; the production cross section is a theory input supported by measurements; X production is hypothetical.

N_h = 7728400; N_X <= 2 N_h B_inv = **1653878** for B_inv<=.107. This is a particle budget across the full dataset, not a fabricated cm^-2 s^-1 flux.

A 2 tonne xenon cylinder with radius47.9cm and height97cm gives maximum chord 136.333cm and column 1.78874e+24 nuclei/cm². Assigning EVERY produced particle that column and 100% readout efficiency gives the continuous-scale ceiling **N_NR < 3.95969e-06**. The stress input (60pb, fN=.326, column2e24) gives 5.35247e-06.

This excludes the chosen source/readout chain as a useful event-count experiment even before geometry, attenuation, timing overlap, backgrounds or calibration losses. It is not an exclusion of the particle theory, off-shell sources, or unmodeled recoil channels.

## Same-coupling source and detection

Gamma_j=lambda_X² v² beta_j/(16 pi mh); B_j=Gamma_j/(Gamma_SM+sum Gamma_j). Source and scattering strengths are not independent: event yield scales as lambda_X^4/(Gamma_SM+lambda_X² Gamma_unit) at fixed masses.

| Lambda [GeV] | lambda_X ceiling | Produced X+anti-X | Produced abs(j)=1 | Arbitrary-boost event upper bound |
|---:|---:|---:|---:|---:|
| 1 | 0.00072299 | 1.6539e+06 | 33763 | 1.9518e-10 |
| 2 | 0.0010333 | 1.6539e+06 | 68268 | 3.8782e-10 |
| 5 | 0.0017555 | 1.6539e+06 | 1.8237e+05 | 9.9821e-10 |
| 10 | 0.0038194 | 1.6539e+06 | 5.474e+05 | 3.8393e-09 |
| 12 | 0.0085654 | 1.6539e+06 | 9.8323e+05 | 1.7737e-08 |
| 12.26 | 0.014125 | 1.6539e+06 | 3.7293e+05 | 4.7777e-08 |
| 12.4 | 0.019729 | 1.6539e+06 | 0 | 9.2781e-08 |
| 12.508 | 0.12926 | 0 | 0 | 0 |
| 15 | 0.1859 | 0 | 0 | 0 |
| 30 | 0.7436 | 0 | 0 | 0 |

The j=+-1 on-shell source closes at Lambda=12.2651GeV, earlier than the j=0 threshold12.508GeV. Above the latter, all on-shell Higgs production is closed even where the old coupling/flux map looked less restrictive.

## Mode information is source-dependent

Direct scalar Higgs decay gives equal populations and identical four-dimensional distributions for j=+1 and -1. The diagonal probe is likewise identical. Thus timing, direction AND independent momentum measurements have optimal sign error50% for this source; they cannot lift an exact sign degeneracy. This differs from the old asymmetric collision source's14.057% S2 prior error.

## Stop and redirect

Do not run detailed beam transport or narrower-Gaussian scans for this rejected chain. A collider production/missing-momentum search avoids the second tiny scattering probability, but on-shell missing-Higgs production alone measures an inclusive width, not signed KK modes. A threshold/recoil-mass strategy or a genuinely mode-sensitive interaction is a separate next hypothesis, not an already feasible instrument.

Full source, transport functional, bound, sign theorem and limitations: tex/route_g_higgs_source.tex. No physical Phi preparation or GeV pump has been built.
