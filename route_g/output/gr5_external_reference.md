# G-R5: data acquisition and external-reference observation model

Status: reference-model-bounded-done; empirical-thermal-test-open. 46/46 checks pass.
No new thermal measurement, no fitted universal-time parameter.

## Acquired evidence and empirical gate

The PTB 2021 archive contains 11 published ratio points. The full ROCIT
2025 release contains 34 selected same-package E2/E3 files with 2,324,201
rows, all marked valid, zero nonfinite values and zero duplicate timestamps.
The selected fields do not include paired bath-state labels, independent
radiometry or a per-run BBR correction ledger. Data integrity is not
eligibility for a thermal-response test. No such fit is run.
Source attribution, version/hash and original-vs-processed distinctions:
data/clock_comparison_audit/PROVENANCE.md.

## Selected reference and physical model

Use a second independent Yb+ E3 ion clock outside the modulated target
enclosure, connected by a comb and stabilized optical paths. Its own
radiation response, drift, common-factor response and path/gravity terms
remain in the model. The measurable candidate is Delta log(F_S/F_R),
not an absolute F_S. Two target-to-reference channels add common-mode
information; their internally derived ratio is not a third independent datum.

At leading static order a composite with E2/E3 weights
(-0.160043, 1.160043) cancels the target's T^4 response
and retains common target-reference response with coefficient one.
It still retains reference drift and common path/gravity biases.
An unconstrained such bias is exactly degenerate with the desired signal.

For the illustrative 300 -> 310 K target change, uncertainties of the
pinned polarizabilities alone give a composite coefficient standard-error
range [2.171e-18, 2.594e-18] over unknown
correlation. This is not an achieved total uncertainty or a confidence
interval. These old inputs do not support a sub-1e-18 common-mode claim.
A 0.1 K reference-bath change contributes -9.635e-20; a 1 cm relative
height change gives about 1.091e-18 in the weak-field model.

## Concrete next step

Obtain the missing clock/thermal/correction records described in
data/EXTERNAL_REFERENCE_CARD.json. Use matched effective time weights
and low-high-high-low blocks to suppress linear drift; nonlinear and
heater-correlated bias do not cancel automatically. Keep radiometry
independent of the frequency law under test. No lab contact has been sent.
Without those records, the empirical test stays open; do not replace it
with another fit, dataset-wide Gaussian scan, or a universal-time claim.

All numerical fixtures test algebra, not fabricated observations.
Full derivation: tex/route_g_external_reference.tex.
