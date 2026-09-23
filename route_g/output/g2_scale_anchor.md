# G2 scale anchor: what would make the energy unit physical?

Verification: **877/877 checks passed**. No measured masses, detector data or physical energy-unit assignment is used.

## Scale is not determined by dimensionless event fractions

Let Lambda=1/R_phys, with the existing dimensionless R=1. The frozen spectrum is m_n^2=Lambda^2[0.25+(n+0.25)^2], mu_j^2=Lambda^2[25+j^2]. M5, MD, p, Omega and all momentum resolutions scale with Lambda; R and kappa5=2 pi Rg scale with Lambda^-1. The dimensionless g, alpha and epsilon are unchanged.

The local kernel Kq has dimension energy^-2 and scales with Lambda^-2. The work per sideband event q Omega scales with Lambda, while the work-weighted rate coefficient sum w Kq q Omega scales with Lambda^-1. Conditional fractions cancel the common rate scale. Under x=Lambda xhat, rho_Lambda(x)=rho_hat(x/Lambda)/Lambda, so integral min(rho_+,rho_-) dx is invariant if detector resolution and the classifier threshold also scale.

Consequently dimensionless data alone have dP/dlog Lambda=0 and zero multinomial Fisher information for log Lambda. A calibrated drive frequency or momentum, a known mass, or an absolute rate with controlled coupling/luminosity is an external dimensional anchor; none is present in this test.

| Synthetic rescaling | Open channels | Total K | Mean pump work/event | Conditional Bayes error |
|---:|---:|---:|---:|---:|
| 0.1 | 66 | 0.00021056349 | 0.00090009722 | 2.20225% |
| 1 | 66 | 2.1056349e-06 | 0.0090009722 | 2.20225% |
| 10 | 66 | 2.1056349e-08 | 0.090009722 | 2.20225% |

The numbers 0.1, 1 and 10 are synthetic unit-rescalings, not GeV values or recommended energy scales. Holding an actual detector resolution fixed would not leave the classifier error invariant.

## Minimal and stronger anchor contracts

A single independently measured and identified mass can set Lambda=m_calibrated/m_hat **only after** freezing the dimensionless ratios. This calibrates a free scale; it does not predict the anchor mass. Likewise, inserting a known drive frequency into the assumed Omega*R=.2 card only defines a candidate radius; that arbitrary card is not a derived radius or a measured spectral spacing.

A stronger test uses three assigned consecutive SIGNED mode masses, not three masses sorted by size. For y_n=m_n^2 and n=-1,0,+1:

A=(y_1-2y_0+y_-1)/2, B=(y_1-y_-1)/2, C=y_0; R=1/sqrt(A), alpha=B/(2A), M5^2=C-B^2/(4A).

Require A>0 and M5^2>=0 for the declared positive-bulk-mass card. A fourth signed n=2 anchor must satisfy y_2-3y_1+3y_0-y_-1=0. Thus three points determine parameters, while a fourth provides a falsifiable closure test.

The synthetic baseline reconstructs A=Lambda^2, alpha=.25 and M5^2=.25 Lambda^2. Deliberately increasing only the fourth squared mass by .07 Lambda^2 violates closure and is rejected. No experimental mass has been tested or excluded.

## What a common mass shift can and cannot hide

A uniform positive delta M^2 changes C and M5^2 but leaves A, B, alpha and the second/third finite differences unchanged. The synthetic shift .3 Lambda^2 verifies this algebra. A portal producing a truly mode-independent shift therefore cannot repair a failed spectral-spacing test merely by adjusting the common mass. This script does not derive a portal or its dynamics.

Integer relabeling alpha->alpha+k,n->n-k leaves masses invariant, so alpha is only defined modulo integers until mode conventions are fixed. Uncalibrated orientation adds a sign ambiguity. The neutral detector still has mu_j=mu_-j; assigning an energy unit does not distinguish those recoil signs.

Next physical input: specify an externally calibrated mass, drive frequency or momentum reference with a justified mode identification, or provide actual detector calibration data. Until then retain all energy requirements as dimensionless ratios. Any physical beam/EFT cutoff and detector coupling require independent justification, not a chosen numerical unit.
