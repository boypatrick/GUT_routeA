#!/usr/bin/env python3
"""Read-only audit of historical P2 beta/threshold compatibility.

Reconstruct product-group beta coefficients from the exported parent Casimirs,
then test matching-scale covariance in the unbroken SU(2)L factor.  This writes
only a new audit ledger; it does not rerun or modify the historical P2 fit.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "route_f/output"
CENSUS = OUT / "p54_p2_two_site_matching.json"
P2 = OUT / "p54_p2_running_thresholds_cosmology.json"
OUTPUT = OUT / "p54_ps_active_census.json"


def rational(value: float | int) -> Q:
    return Q(float(value)).limit_denominator(10000)


def beta_coefficients(scalars: list[dict]) -> tuple[list[Q], list[list[Q]]]:
    """Weyl fermions and complex scalar convention, ordering (4,L,R)."""
    adjoint_c = [Q(4), Q(2), Q(2)]
    # Three generations. T includes spectator multiplicities and generations.
    fermions = [
        ([Q(15, 8), Q(3, 4), Q(0)], [Q(3), Q(6), Q(0)]),
        ([Q(15, 8), Q(0), Q(3, 4)], [Q(3), Q(0), Q(6)]),
    ]
    scalar_ct = [(row["C"], row["T_complex"]) for row in scalars]
    a = [
        -Q(11, 3) * adjoint_c[i]
        + Q(2, 3) * sum(t[i] for _, t in fermions)
        + Q(1, 3) * sum(t[i] for _, t in scalar_ct)
        for i in range(3)
    ]
    b = [
        [
            (-Q(34, 3) * adjoint_c[i] ** 2 if i == j else Q(0))
            + sum(
                (2 * c[j] + (Q(10, 3) * adjoint_c[i] if i == j else 0)) * t[i]
                for c, t in fermions
            )
            + sum(
                (4 * c[j] + (Q(2, 3) * adjoint_c[i] if i == j else 0)) * t[i]
                for c, t in scalar_ct
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    return a, b


def encode(value):
    if isinstance(value, Q):
        return {"exact": str(value), "value": float(value)}
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [encode(v) for v in value]
    return value


def run() -> dict:
    parents = json.loads(CENSUS.read_text())
    historical = json.loads(P2.read_text())
    rows = {row["label"]: row for row in parents["parent_census"]}
    labels = [
        "phi10:(1,2,2)",
        "Sigma126:(15,2,2)",
        "Sigma126:(10-pair,3,1)",
        "Sigma126:(10-pair,1,3)",
    ]

    def scalar(label):
        row = rows[label]
        casimir = [rational(row[k]) for k in ("C4", "C2L", "C2R")]
        # All selected parents belong to complex 10 or complex 126.
        dimension = Q(row["dimension"], 2)
        indices = [dimension * c / d for c, d in zip(casimir, (15, 3, 3))]
        return {"label": label, "dimension_complex": dimension, "C": casimir, "T_complex": indices}

    active = [scalar(label) for label in labels]
    a, b = beta_coefficients(active)
    old_a = [rational(x) for x in historical["rge_coefficients"]["a_PS"]]
    old_b = [[rational(x) for x in r] for r in historical["rge_coefficients"]["b_PS"]]
    da = [old_a[i] - a[i] for i in range(3)]
    db = [[old_b[i][j] - b[i][j] for j in range(3)] for i in range(3)]
    with_six_a, with_six_b = beta_coefficients(active + [scalar("Sigma126:(6,1,1)")])

    a_sm_l = rational(historical["rge_coefficients"]["a_SM"][1])
    a_ps_l = old_a[1]
    required_lambda_slope = -6 * (a_ps_l - a_sm_l)
    delta_r = rows[parents["intermediate_breaking_parent"]]
    # C_L P_DeltaR=0 exactly. Intermediate vectors are also SU(2)L singlets.
    actual_lambda_slope = Q(0)
    residual_exact_numerator = a_ps_l - a_sm_l + actual_lambda_slope / 6
    residual = float(residual_exact_numerator) / (2 * math.pi)

    def match_residual(t: float, slope: float) -> float:
        alpha_sm = 40.0 - float(a_sm_l) * t / (2 * math.pi)
        alpha_ps = 40.0 - float(a_ps_l) * t / (2 * math.pi)
        return alpha_sm - alpha_ps + slope * t / (12 * math.pi)

    step = 0.1
    fd_bad = (match_residual(step, 0) - match_residual(-step, 0)) / (2 * step)
    fd_repaired = (
        match_residual(step, float(required_lambda_slope))
        - match_residual(-step, float(required_lambda_slope))
    ) / (2 * step)
    # Independent scalar census at MI: bidoublets and DeltaL, minus one SM
    # complex doublet. DeltaR has no SU(2)L charge.
    delta_a_l_from_census = sum(row["T_complex"][1] for row in active) / 3 - Q(1, 6)
    # Minimal stage-consistent repair. Hypercharge is GUT normalized:
    # alpha_1^-1=(2/5)alpha_4^-1+(3/5)alpha_R^-1.
    def project_sm(v):
        return [v[0], v[1], Q(2, 5) * v[0] + Q(3, 5) * v[2]]

    active_real_indices = [2 * sum(row["T_complex"][i] for row in active) for i in range(3)]
    full_real_indices = [
        sum(Q(row["dimension"]) * rational(row[k]) / d for row in rows.values())
        for k, d in zip(("C4", "C2L", "C2R"), (15, 3, 3))
    ]
    light_sm_real_indices = [Q(0), Q(1), Q(3, 5)]
    vector_i = [Q(1), Q(0), Q(14, 5)]
    vector_u = [Q(4), Q(6), Q(6)]
    all_heavy_i = [x - y for x, y in zip(project_sm(active_real_indices), light_sm_real_indices)]
    physical_scalar_i = [x - y for x, y in zip(all_heavy_i, vector_i)]
    physical_scalar_u = [x - y - z for x, y, z in zip(full_real_indices, active_real_indices, vector_u)]
    repaired_slope_i = [21 * v - s for v, s in zip(vector_i, physical_scalar_i)]
    repaired_slope_u = [21 * v - s for v, s in zip(vector_u, physical_scalar_u)]
    a_sm = [rational(x) for x in historical["rge_coefficients"]["a_SM"]]
    a_so10 = -Q(11, 3) * 8 + Q(2, 3) * 6 + Q(1, 6) * 84
    required_slope_i = [-6 * (x - y) for x, y in zip(project_sm(a), a_sm)]
    required_slope_u = [-6 * (a_so10 - x) for x in a]
    checks = [
        ("parent dimensions and Casimirs give (2/3,26/3,26/3)", a == [Q(2, 3), Q(26, 3), Q(26, 3)]),
        ("two-loop discrepancy is only b44=38/3", db == [[Q(38, 3), Q(0), Q(0)], [Q(0)] * 3, [Q(0)] * 3]),
        ("one-loop discrepancy is only a4=1/3", da == [Q(1, 3), Q(0), Q(0)]),
        ("one extra complex six exactly reproduces historical a and b", with_six_a == old_a and with_six_b == old_b),
        ("selected intermediate parent is SU2L singlet", rational(delta_r["C2L"]) == 0),
        ("exported scalar and vector MI SU2L thresholds vanish", abs(parents["thresholds"]["lambda_scalar_MI_3_2_1"][1]) < 1e-20 and abs(parents["vector"]["lambda_vector_MI_3_2_1"][1]) < 1e-20),
        ("required SU2L threshold slope is -71", required_lambda_slope == -71),
        ("independent scalar census gives same beta discontinuity", delta_a_l_from_census == a_ps_l - a_sm_l),
        ("finite difference reproduces nonzero historical matching-scale defect", abs(fd_bad - residual) < 1e-12),
        ("required logarithmic coefficient cancels scale derivative", abs(fd_repaired) < 1e-12),
        ("complete scalar real trace index is 84 in all PS factors", full_real_indices == [Q(84)] * 3),
        ("minimal repaired MI slope cancels all three gauge directions", repaired_slope_i == required_slope_i),
        ("minimal repaired MU slope cancels all three gauge directions", repaired_slope_u == required_slope_u),
        ("repaired slopes are (-46,-71,-41/5) and (72,120,120)", repaired_slope_i == [Q(-46), Q(-71), Q(-41, 5)] and repaired_slope_u == [Q(72), Q(120), Q(120)]),
    ]
    return encode({
        "schema": "route-f-p54-ps-active-census-v1",
        "date": "2026-09-05",
        "scope": "algebraic consistency audit; no refit and no change to historical P2 ledgers",
        "active_scalar_census": active,
        "reconstructed_a_PS": a,
        "reconstructed_b_PS_gauge_and_matter_only": b,
        "historical_minus_reconstructed_a": da,
        "historical_minus_reconstructed_b": db,
        "extra_six_interpretation": "historical coefficients require one additional complex (6,1,1); its 10/126 origin is not determined by gauge indices",
        "SU2L_matching_scale_audit": {
            "matching_identity": "alpha_SM^-1=alpha_PS^-1-lambda_I/(12*pi)",
            "beta_discontinuity": a_ps_l - a_sm_l,
            "required_dlambda_dlogmu": required_lambda_slope,
            "actual_dlambda_dlogmu_from_DeltaR_only_projector": actual_lambda_slope,
            "residual_derivative_exact": "71/(12*pi)",
            "residual_derivative": residual,
            "finite_difference_derivative": fd_bad,
            "derivative_after_required_log_coefficient": fd_repaired,
        },
        "minimal_stage_consistent_repair": {
            "status": "all group-theoretic logarithmic sum rules closed; finite threshold matrices and scale refit not performed",
            "active_parent_labels": labels,
            "full_scalar_real_indices_PS": full_real_indices,
            "active_scalar_real_indices_PS": active_real_indices,
            "light_scalar_real_indices_SM": light_sm_real_indices,
            "heavy_vector_real_indices_MI_SM": vector_i,
            "heavy_vector_real_indices_MU_PS": vector_u,
            "Goldstone_prescription": "exclude one real Goldstone per massive vector from scalar trace; use the same site-resolved orbit projector as the vector",
            "physical_heavy_scalar_real_indices_MI_SM": physical_scalar_i,
            "physical_heavy_scalar_real_indices_MU_PS": physical_scalar_u,
            "scalar_plus_vector_lambda_slope_MI": repaired_slope_i,
            "required_lambda_slope_MI": required_slope_i,
            "scalar_plus_vector_lambda_slope_MU": repaired_slope_u,
            "required_lambda_slope_MU": required_slope_u,
            "a_SO10_full_model": a_so10,
            "finite_background_warning": "at simultaneous nonzero vevs, parent and mass/orbit projectors need not commute; the finite matching must be derived in a common staged scheme rather than inferred from these index identities alone",
        },
        "physical_gates": {
            "historical_beta_table_matches_stated_four_parent_active_census": False,
            "historical_MI_matching_is_scale_covariant_at_one_loop": False,
            "whole_P54_model_excluded": False,
            "required_next_step": "declare one PS active census, regenerate its betas and complementary U/I thresholds, then solve scales again",
        },
        "checks": [{"name": name, "pass": ok} for name, ok in checks],
        "summary": {"passed": sum(ok for _, ok in checks), "total": len(checks), "all_pass": all(ok for _, ok in checks)},
        "sources": [
            {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in (CENSUS, P2, Path(__file__).resolve())
        ],
    })


if __name__ == "__main__":
    report = run()
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"PS active-census audit: {report['summary']['passed']}/{report['summary']['total']} checks")
    print(f"Historical MI matching scale derivative: {report['SU2L_matching_scale_audit']['residual_derivative']:.12g}")
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)
