# G-R6 equivalence-principle baseline

Status: theory baseline bounded-done; empirical validation OPEN.
49/49 implementation checks pass.
No new measured redshift, anomaly fit or confidence bound.

## Fixed-input example, not experimental observations

Same pinned G-R4 static E2/E3 polarizabilities and G-R5 external reference.
Local temperatures 300/310 K; height separation 1 m; conventional g scale
9.80665 m/s^2 is not laboratory geodesy. Every entry is q_cell - q_00.

| Setting | E2 log-frequency change | E3 log-frequency change |
|---|---:|---:|
| T0_Phi0 | +0.00000e+00 | +0.00000e+00 |
| T1_Phi0 | -7.33714e-17 | -1.01225e-17 |
| T0_Phi1 | +1.09114e-16 | +1.09114e-16 |
| T1_Phi1 | +3.57423e-17 | +9.89912e-17 |

Temperature contrast: E2 -7.33714e-17, E3 -1.01225e-17.
Potential contrast: both +1.09114e-16.
Log interaction: zero analytically for the separable static model;
the JSON retains harmless floating-point residuals, not fitted zeros.
E3/E2 thermal contrast +6.32488e-17; its gravity contrast is zero.

## Independently measured external potential anchor

Grotti et al. (Physical Review Applied 21, L061001, 2024;
https://arxiv.org/abs/2309.14953v3) report an independent geodetic potential
difference magnitude of 3915.88(0.30) m^2/s^2. We do not use their
clock-inferred potential as calibration. The EEP prediction is
|Delta log N| = 4.3570041e-14, with propagated potential uncertainty
3.3380e-18; the same fractional prediction applies to E2 and E3
if compared across that calibrated potential with other shifts controlled.
The actual published clocks were Sr, not our E2/E3 four-setting experiment.
This is a source-backed gravity scale, not borrowed apparatus performance
or a new EEP measurement, significance test or limit.

## What becomes separable

In a static metric, local spectral response and the lapse ratio multiply.
Their logarithms add. A two-temperature/two-potential design distinguishes
their columns with independently supplied thermal and geodetic calibration.
The G-R5 weighted composite remains a leading thermal-rejection diagnostic,
but matched-temperature height comparisons do not require its cancellation.
Even an uncertain fixed BBR shift cancels at matched local radiation states.

A +0.1 K mismatch at the upper site at both setpoints leaks
-7.3447e-19 (E2) and -1.0133e-19 (E3) into the height contrast.
Thus radiation matching, transport-induced systematics and reference/link
drift, not an arbitrary new time law, are the relevant controls.
These are conditional sensitivities, not an achieved error budget.

## Boundaries and empirical handoff

Use independently measured DeltaPhi; extracting it from these clocks and
then verifying their redshift would be circular. For a 1e-18 absolute
redshift uncertainty, the potential contribution alone requires sigmaPhi
about 0.08988 m^2/s^2 (about
0.916 cm at this g), before other uncertainties.

An unconstrained height-correlated common bias is exactly degenerate with
a common anomalous redshift. A nonzero log interaction is a failure of the
controlled separability model, not by itself an EEP violation. Two chosen
transitions do not verify all EEP components or quantum off-diagonal tests.

The optional symmetric eight-plateau order cancels ideal linear drift,
not quadratic drift or state-locked bias. Real sensitivity windows,
settling, reversal and repeated blocks remain required. Four cell means
alone do not establish a noise distribution or experimental significance.

Independent thermostats are not one global equilibrium bath. The latter
instead obeys T_local*N=constant. Its illustrative 300 K/1 m temperature
change is -3.2734e-14 K; we do not demand its direct measurement
or confuse it with the imposed 10 K change.

The required acquisition contract is data/EEP_BASELINE_CARD.json.
The G-R5 missing-record gate stays open; no laboratory was contacted.
Derivation: tex/route_g_equivalence_baseline.tex.
