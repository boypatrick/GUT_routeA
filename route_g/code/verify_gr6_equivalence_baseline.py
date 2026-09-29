#!/usr/bin/env python3
"""G-R6: fixed-input EEP/thermal factorial baseline. No experimental fit.

Dimensionless large-amplitude fixtures test identities, never stand in for
observations. Physical residuals are scaled before comparison so tolerances
cannot silently accept a completely wrong ~1e-16 prediction.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr6_equivalence_baseline"


def run():
    card = json.loads((ROOT / "data/EEP_BASELINE_CARD.json").read_text())
    bbr = json.loads((ROOT / "data/BBR_CLOCK_INPUTS.json").read_text())
    const = bbr["constants"]
    h, kb, c, eps0 = [const[k] for k in
                      ("h_J_s", "kB_J_K", "c_m_s", "epsilon0_F_m")]
    nu = np.array([v["frequency_Hz"] for v in bbr["transitions"]])
    alpha = np.array([v["delta_alpha_SI"] for v in bbr["transitions"]])
    arad = 8 * np.pi**5 * kb**4 / (15 * h**3 * c**3)
    beta = -alpha / (2 * h * eps0 * nu)
    inp = card["illustrative_inputs"]
    T = np.array(inp["local_radiation_temperatures_K"])
    heights = np.array(inp["target_height_settings_m"])
    delta_phi = inp["g_scale_m_s2"] * np.diff(heights)[0]
    gamma = delta_phi / c**2
    anchor = card["published_independent_potential_anchor"]
    anchor_gamma = anchor["potential_difference_magnitude_m2_s2"] / c**2
    anchor_sigma = anchor["published_uncertainty_m2_s2"] / c**2
    shifts = beta[:, None] * arad * T[None, :]**4
    B = np.log1p(shifts)
    theta = B[:, 1] - B[:, 0]
    # q_ab - q_00. Static spectral shift is retained in local proper units.
    q = np.column_stack([np.zeros(2), theta,
                         np.full(2, gamma), theta + gamma])
    W = np.array(list(card["contrasts"].values()))
    expected = np.column_stack([theta, np.full(2, gamma), np.zeros(2)])
    weights = np.array([beta[1], -beta[0]]) / (beta[1] - beta[0])
    d = np.array([-1., 1.])
    checks = []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, scale=1.0):
        check(name, np.allclose(np.asarray(a) / scale, np.asarray(b) / scale,
                               rtol=2e-12, atol=2e-12))

    for rel, digest in card["preserved_sha256"].items():
        check("preserved:" + rel,
              hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest)
    check("scope:no_observational_fit", not card["empirical_data_ready"]
          and not card["new_time_law"])
    close("factorial:three_contrasts", q @ W.T, expected, 1e-16)
    close("factorial:constant_offset_cancels", W @ np.ones(4), np.zeros(3))
    close("factorial:contrast_orthogonality",
          W @ W.T, np.diag([1., 1., 4.]))
    design = np.array([[1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1.]])
    check("factorial:additive_rank_three", np.linalg.matrix_rank(design) == 3)
    close("factorial:interaction_annihilates_additive_model",
          W[2] @ design, np.zeros(3))
    saturated = np.column_stack([design, [0, 0, 0, 1]])
    check("factorial:interaction_fit_saturates_four_cells",
          np.linalg.matrix_rank(saturated) == 4)
    close("ratio:gravity_cancels", (d @ q) @ W[1], 0, 1e-16)
    close("ratio:retains_differential_thermal", (d @ q) @ W[0],
          theta[1] - theta[0], 1e-16)
    close("composite:retains_gravity", (weights @ q) @ W[1], gamma, 1e-16)
    close("composite:leading_thermal_null", weights @ (beta / beta[0]), 0)
    check("sign:raised_target_faster", gamma > 0)
    check("sign:both_heated_transitions_lower", np.all(theta < 0))
    check("anchor:independent_geodetic_not_clock_inferred",
          anchor["potential_difference_magnitude_m2_s2"] == 3915.88
          and anchor["published_uncertainty_m2_s2"] == .30)
    close("anchor:universal_redshift_for_both_transitions",
          np.full(2, anchor_gamma) @ d, 0, 1e-14)
    close("anchor:potential_uncertainty_propagation",
          anchor_sigma * c**2, .30)
    close("regression:GR5_thermal_shifts", shifts[:, 1] - shifts[:, 0],
          json.loads((ROOT / "output/gr5_external_reference.json").read_text())
          ["static_target_300_to_310_fractional_shifts"], 1e-16)
    close("regression:GR5_composite_weights", weights,
          json.loads((ROOT / "output/gr5_external_reference.json").read_text())
          ["composite_weights_E2_E3"])
    # Exact finite-factor fixtures: log separability vs a raw product.
    thermal, grav = -.07, .03
    R = np.array([1, np.exp(thermal), np.exp(grav), np.exp(thermal + grav)])
    close("finite:log_interaction_zero", W[2] @ np.log(R), 0)
    close("finite:raw_cross_is_product", W[2] @ (R - 1),
          np.expm1(thermal) * np.expm1(grav))
    check("finite:raw_cross_not_zero", abs(W[2] @ (R - 1)) > 1e-4)
    close("metric:time_coordinate_rescaling_cancels",
          np.log((1.03 * 7) / (1.01 * 7)), np.log(1.03 / 1.01))
    # Thermal uncertainties need not set a floor on the matched-temperature
    # height contrast. Scaled fixtures cover arbitrary two thermal shifts.
    for level_pair in ([2., 5.], [-9., 4.], [0., 0.]):
        vals = np.array([level_pair[0], level_pair[1],
                         level_pair[0] + .2, level_pair[1] + .2])
        close("matched_bath:height_independent_of_B:" + str(level_pair),
              W[1] @ vals, .2)
    # A temperature mismatch at the raised site does leak and must be measured.
    delta_T = .1
    raised_T = T + delta_T
    raised_B = np.log1p(beta[:, None] * arad * raised_T[None, :]**4)
    mismatch_q = np.column_stack([B[:, 0], B[:, 1],
                                 raised_B[:, 0] + gamma, raised_B[:, 1] + gamma])
    mismatch = .5 * np.sum(raised_B - B, axis=1)
    close("mismatched_bath:height_leakage", mismatch_q @ W[1] - gamma,
          mismatch, 1e-16)
    check("mismatched_bath:not_universal", abs(mismatch[0]) > abs(mismatch[1]) > 0)
    dBdT = 4 * shifts / (T[None, :] * (1 + shifts))
    small = .001
    fd = (np.log1p(beta[:, None] * arad * (T[None, :] + small)**4)
          - np.log1p(beta[:, None] * arad * (T[None, :] - small)**4)) / (2 * small)
    check("radiometry:analytic_derivative",
          np.allclose(fd / 1e-18, dBdT / 1e-18, rtol=1e-8, atol=1e-8))
    # A height-locked shared bias has the same design column as an anomaly.
    H = np.array([0., 0., 1., 1.])
    G = np.vstack([H, H]).reshape(-1)
    ident = np.column_stack([G, G])
    check("identifiability:common_gravity_and_bias_rank_one",
          np.linalg.matrix_rank(ident) == 1)
    close("identifiability:exact_null_direction", ident @ [1., -1.], np.zeros(8))
    fixture_alpha = np.array([.02, -.03])
    fixture_q = np.outer(1 + fixture_alpha, H) + .4 * H
    close("anomaly:differential_rejects_shared_bias",
          (d @ fixture_q) @ W[1], fixture_alpha[1] - fixture_alpha[0])
    close("anomaly:composite_retains_shared_bias",
          (weights @ fixture_q) @ W[1] - 1, weights @ fixture_alpha + .4)
    # Pair-symmetric acquisition cancels an ideal linear temporal drift.
    order = np.array(card["protocol"]["suggested_cell_order"])
    t = np.arange(-7., 8., 2.)
    average = np.array([(order == j) / 2 for j in range(4)])
    close("protocol:each_cell_visited_twice", average.sum(axis=1), np.ones(4))
    close("protocol:symmetric_time_centers", average @ t, np.zeros(4))
    close("protocol:constant_and_linear_drift_cancel",
          W @ average @ (3 + 2 * t), np.zeros(3))
    close("protocol:quadratic_drift_survives",
          W @ average @ t**2, [-8., -32., 32.])
    close("protocol:height_locked_bias_survives",
          W @ average @ H[order], [0., 1., 0.])
    # Covariance fixture, ordered [E2: four cells, E3: four cells].
    # Reference varies by cell but is shared by contemporaneous transitions.
    Kref = 9 * np.eye(4)
    Vref = np.kron(np.ones((2, 2)), Kref)
    a_diff, a_common = np.kron(d, W[1]), np.kron(weights, W[1])
    close("covariance:matched_reference_removed_in_difference",
          a_diff @ Vref @ a_diff, 0)
    close("covariance:changing_reference_remains_in_common",
          a_common @ Vref @ a_common, 9)
    close("covariance:constant_reference_offset_cancels",
          a_common @ np.ones((8, 8)) @ a_common, 0)
    Wmismatch = W[1] + np.array([.1, -.1, 0., 0.])
    amismatch = np.concatenate([-W[1], Wmismatch])
    check("covariance:unmatched_windows_leak_reference",
          amismatch @ Vref @ amismatch > 0)
    # Tolman equilibrium is an alternative boundary condition, not imposed
    # simultaneously with independent local thermostats.
    lapse_fixture = np.exp(np.array([0., .02]))
    temp_fixture = 300 / lapse_fixture
    close("equilibrium:T_times_lapse_constant",
          temp_fixture * lapse_fixture, [300., 300.])
    check("equilibrium:equal_local_T_not_global_equilibrium",
          not np.allclose(300 * lapse_fixture, [300., 300.]))
    tolman_delta_T = T[0] * np.expm1(-gamma)
    potential_sigma_for_1e18 = c**2 * 1e-18
    height_sigma_for_1e18 = potential_sigma_for_1e18 / inp["g_scale_m_s2"]
    close("units:one_cm_redshift", gamma * .01, 1.0911369672198218e-18, 1e-18)
    close("units:potential_uncertainty", potential_sigma_for_1e18 / c**2, 1e-18, 1e-18)
    check("scope:G-R5_empirical_gate_remains_open",
          not json.loads((ROOT / "data/EXTERNAL_REFERENCE_CARD.json").read_text())
          ["empirical_data_ready"])
    failed = [v["name"] for v in checks if not v["passed"]]
    result = dict(
        stage="G-R6", status="theory-baseline-bounded-done; empirical-test-open",
        measurement_claim=False, fitted_parameters=False,
        input_card_sha256=hashlib.sha256((ROOT / "data/EEP_BASELINE_CARD.json").read_bytes()).hexdigest(),
        delta_phi_illustrative_m2_s2=float(delta_phi), gamma_illustrative=float(gamma),
        published_geodetic_anchor_redshift_magnitude=float(anchor_gamma),
        published_geodetic_anchor_propagated_uncertainty=float(anchor_sigma),
        log_thermal_contrast_E2_E3=theta.tolist(),
        relative_log_four_cells_E2_E3=q.tolist(),
        contrast_order=["temperature", "potential", "interaction"],
        contrasts_E2_E3=(q @ W.T).tolist(),
        E3_over_E2_temperature_contrast=float(d @ theta),
        common_weights_E2_E3=weights.tolist(),
        raw_fractional_interaction_E2_E3=(np.expm1(gamma) * np.expm1(theta)).tolist(),
        higher_site_plus_0p1K_BBR_leakage_E2_E3=mismatch.tolist(),
        local_temperature_derivatives_per_K=dBdT.tolist(),
        tolman_300K_deltaT_over_1m_K=float(tolman_delta_T),
        potential_sigma_for_1e18_m2_s2=float(potential_sigma_for_1e18),
        height_sigma_for_1e18_m=float(height_sigma_for_1e18),
        summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
        checks=checks)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    rows = "\n".join(
        f"| {name} | {q[0,j]:+.5e} | {q[1,j]:+.5e} |"
        for j, name in enumerate(card["factorial_order"]))
    OUT.with_suffix(".md").write_text(f"""# G-R6 equivalence-principle baseline

Status: theory baseline bounded-done; empirical validation OPEN.
{len(checks)-len(failed)}/{len(checks)} implementation checks pass.
No new measured redshift, anomaly fit or confidence bound.

## Fixed-input example, not experimental observations

Same pinned G-R4 static E2/E3 polarizabilities and G-R5 external reference.
Local temperatures 300/310 K; height separation 1 m; conventional g scale
9.80665 m/s^2 is not laboratory geodesy. Every entry is q_cell - q_00.

| Setting | E2 log-frequency change | E3 log-frequency change |
|---|---:|---:|
{rows}

Temperature contrast: E2 {theta[0]:+.5e}, E3 {theta[1]:+.5e}.
Potential contrast: both {gamma:+.5e}.
Log interaction: zero analytically for the separable static model;
the JSON retains harmless floating-point residuals, not fitted zeros.
E3/E2 thermal contrast {d@theta:+.5e}; its gravity contrast is zero.

## Independently measured external potential anchor

Grotti et al. (Physical Review Applied 21, L061001, 2024;
https://arxiv.org/abs/2309.14953v3) report an independent geodetic potential
difference magnitude of 3915.88(0.30) m^2/s^2. We do not use their
clock-inferred potential as calibration. The EEP prediction is
|Delta log N| = {anchor_gamma:.7e}, with propagated potential uncertainty
{anchor_sigma:.4e}; the same fractional prediction applies to E2 and E3
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
{mismatch[0]:+.4e} (E2) and {mismatch[1]:+.4e} (E3) into the height contrast.
Thus radiation matching, transport-induced systematics and reference/link
drift, not an arbitrary new time law, are the relevant controls.
These are conditional sensitivities, not an achieved error budget.

## Boundaries and empirical handoff

Use independently measured DeltaPhi; extracting it from these clocks and
then verifying their redshift would be circular. For a 1e-18 absolute
redshift uncertainty, the potential contribution alone requires sigmaPhi
about {potential_sigma_for_1e18:.4g} m^2/s^2 (about
{100*height_sigma_for_1e18:.3g} cm at this g), before other uncertainties.

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
change is {tolman_delta_T:.4e} K; we do not demand its direct measurement
or confuse it with the imposed 10 K change.

The required acquisition contract is data/EEP_BASELINE_CARD.json.
The G-R5 missing-record gate stays open; no laboratory was contacted.
Derivation: tex/route_g_equivalence_baseline.tex.
""")
    print(json.dumps(result["summary"]))
    print("gamma:", gamma, "thermal:", theta, "mismatch:", mismatch)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
