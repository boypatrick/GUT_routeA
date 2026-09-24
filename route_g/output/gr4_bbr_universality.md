# G-R4: physical thermal source and cross-transition test

Status: bounded-source-audit-done. Checks: 91/91.
No fitted parameters or new experimental data.

## Selected source and independently fixed sign

Same-ion 171Yb+ E2 and E3 clock transitions in one ideal isotropic
blackbody environment. Microscopic coupling is -d.E. Published static
polarizabilities, converted to alpha_excited-alpha_ground, are
6.9(1.4)e-40 and 0.888(0.016)e-40 J m^2/V^2.
Both are positive, so the thermal electric-dipole shifts are negative.
The 2005 E2 paper uses the opposite polarizability sign convention.
Source files, independent calibration methods and caveats are pinned
in data/BBR_CLOCK_INPUTS.json; these are not claimed to be the newest
or most precise possible constants.

## Static leading-order prediction at 300 K

Thermal electromagnetic density u=6.12824e-06 J/m^3;
E_rms=831.943 V/m. This assumes conventional spatial QED and
does not identify u with Chen Ed or the previous energy-per-cell candidate.

| Transition | Shift / Hz | Fractional shift |
|---|---:|---:|
| E2 | -0.360371 | -5.23523e-16 |
| E3 | -0.0463782 | -7.22266e-17 |

The fractional responses differ by a factor 7.248.
R=(nu_E3/nu_E2) gives R(300)/R(0)-1=4.51296e-16.
An illustrative 300 -> 310 K modulation gives
R(310)/R(300)-1=6.32488e-17, not an observed signal.

Static leading-order quantities are deliberately separated from full
finite-temperature predictions. The published E3 eta(300)=-0.0015
changes its 300 K shift to -0.0463087 Hz. A complete E2 dynamic
response and response covariance are not supplied.
The unequal leading T^4 coefficients already prevent exact equality
over a low-temperature interval; T^6 terms cannot change that coefficient.

## Robustness, not a significance fit

Across the full 3-published-sigma polarizability box and independent
assumed +/-10% dynamic factors, the leading ratio offset at 300 K
stays above 1.00627e-16. This is a conditional deterministic
stress statement, NOT a confidence level or a derived bound on every
neglected atomic correction. Unknown cross-paper covariance is retained.

## Essential interpretation limit

The ordinary thermal Stark channel cannot be represented solely as
H -> f(u) H + c(u) I on these three clock levels.
It does NOT exclude an additional exactly common time factor:
nu_i(u)=F(u) nu_i^0[1+s_i(u)] makes F cancel from nu_E3/nu_E2.
Subtracting published BBR corrections from both outputs and then finding
a constant ratio is not an independent test of those corrections.
No raw temperature-modulated observations have been analyzed.

## Next useful work

Before an empirical test, specify actual common-bath radiometry and
interleaved raw E2/E3 readings, retain the deliberately modulated BBR
term, and control temperature-correlated non-BBR shifts.
A claimed common time factor needs an independently specified reference
outside the changed environment or another non-clock observable.
Alternatively test equal total radiation density with different spectra;
the full response is spectral, not a scalar-density law.

Complete formulas/proofs: tex/route_g_bbr_universality.tex.
Run: python3 route_g/code/verify_gr4_bbr_universality.py.
