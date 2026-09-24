# G2-T: mass measurement and signed-mode information are separate

Verification: **97/97 checks passed**.

No action, coupling, physical source, response matrix or old result is changed. No Gaussian scan or actual hardware performance is claimed.

## Different masses: TOF plus a genuinely independent constraint

m²=p²[(cT/L)²-1]. TOF alone fixes m/p, not m. Neutral X has no ordinary magnetic-curvature momentum readout. Two time hits still only determine velocity.

The logarithmic Jacobian for (p,T,L) is (1,gamma²,-gamma²); the full covariance must be used. Negative noisy m² estimates are not silently clipped.

At-rest Higgs illustration only: Lambda=5 GeV, E=62.54 GeV, L=10 m. This is not a real LHC source distribution.

| j | mass [GeV] | p [GeV] | TOF [ns] | Exact reconstructed-mass range [GeV] |
|---:|---:|---:|---:|---:|
| 0 | 25.000000 | 57.325837 | 36.390395 | [24.873308, 25.126568] |
| 1 | 25.495098 | 57.107369 | 36.529609 | [25.370216, 25.619875] |
| -1 | 25.495098 | 57.107369 | 36.529609 | [25.370216, 25.619875] |
| 5 | 35.355339 | 51.587320 | 40.438423 | [35.254238, 35.456502] |

Intervals assume bounded |delta p|/p<=0.1%, |delta T|<=20 ps and |delta L|<=1 mm. These are uncalibrated requirements, not Gaussian sigmas, confidence intervals or event-by-event guarantees for the collider.
The j=0 and abs(j)=1 nominal TOFs differ by 139.213887 ps. Their exact error boxes first touch at a timing halfwidth of 48.672148 ps with the other two bounds fixed.

## Same action, alternative constraint: directional elastic recoil

For a known target mass M at rest, q=sqrt(ER(ER+2M)), D=beta q cos(theta)-ER, E=M ER/D, p=beta E and m=E sqrt(1-beta²). Require 0<beta<1, D>0 and positive outgoing energy. The reconstructed p is correlated with beta and ER, not an additional independent datum.
At A=132, ER=10 keV and the abs(j)=1 example, theta=89.96247023 degrees. The local angular budget alone for 1% mass error is 0.895561 arcsec, assuming exact beta and ER. This severe conditioning is a diagnostic, not detector performance.
Natural xenon also has an unobserved isotope mixture. The existing S2-only response has neither recoil-vector information nor a tagged source-to-hit X flight time. No directional or timing calibration is imported.

## Same mass, opposite sign

Existing undriven Phi-X source correlations can label j=n+l-m if the needed incoming/output labels are measured. With unobserved l, G2-R's full COM kinematic inversion works only for d=m-n nonzero and known incoming conditions. This is revalidated, not a constructed source or independent conservation-law test. The driven branch needs the work exchange too.
A new sign-odd amplitude alone is insufficient: |Me+j Mo|²-|Me-j Mo|²=4j Re(Me* Mo). There must be coherent even/odd contributions in the same observed channel, a calibrated signed external response, or an informative preparation. No such new interaction is activated.
Direct tree Higgs production plus the diagonal portal retains equal signed distributions and 50% optimal sign error.

## Signal budget and next gate

The previous arbitrary-momentum, efficiency-one ceiling stays 3.9596863e-06. Any added timing/momentum/direction gate <=1 can only lower it. No timing improvement revives the rejected on-shell-Higgs/single-pass-xenon chain.
Next require an actual production tag and independent momentum/energy mechanism before instrument forecasts; for sign choose source correlation or specify a new sign-sensitive interaction for separate authorization/audit. Production observables avoiding rescattering remain preferable for the rejected source. Full response, background and coherent finite-pulse likelihoods are not supplied.
