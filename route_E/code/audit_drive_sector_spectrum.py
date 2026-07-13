#!/usr/bin/env python3
"""Drive-sector representation and threshold audit.

No web lookup is used.  This upgrades the Route-B link-driving construction by
assigning explicit Spin(10) representations to the driving fields:

* three 54-link/54-driver pairs for AA, AB, AC alignment;
* one 210-link/210-driver pair for AA alignment;
* two full 54 relative-source drivers Y_AB, Y_AC;
* gauge singlet drivers for the link norms and singlet/54 ratios.

The script checks whether the orthogonal alignment modes create non-singlet
thresholds.  The key point is that a full 54 or 210 driver pairs every
orthogonal component with the corresponding link fluctuation, so the threshold
is a complete Spin(10) multiplet threshold.  It has zero projection onto the
two-dimensional gauge-unification mismatch plane, although a large universal
shift can still change alpha_G and therefore dimension-6 proton decay.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "drive_sector_spectrum"


SPIN10_REPS = {
    "54": {
        "dimension": 54,
        "dynkin_index": 12.0,
        "aligned_sm_singlets": 1,
        "ps_decomposition": [
            "(1,1,1)",
            "(20',1,1)",
            "(6,2,2)",
            "(1,3,3)",
        ],
        "su5_decomposition": ["24", "15", "15bar"],
    },
    "210": {
        "dimension": 210,
        "dynkin_index": 56.0,
        "aligned_sm_singlets": 1,
        "ps_decomposition": [
            "(1,1,1)",
            "(15,1,1)",
            "(15,3,1)",
            "(15,1,3)",
            "(6,2,2)",
            "(10,2,2)",
            "(10bar,2,2)",
        ],
        "su5_decomposition": [
            "1",
            "5",
            "5bar",
            "10",
            "10bar",
            "24",
            "40",
            "40bar",
            "75",
        ],
    },
}


def full_pair_beta(rep: str) -> np.ndarray:
    """One link field plus one driver field in the same real Spin(10) rep."""

    index = float(SPIN10_REPS[rep]["dynkin_index"])
    return np.array([2.0 * index, 2.0 * index, 2.0 * index], dtype=float)


def full_field_beta(rep: str) -> np.ndarray:
    """One full chiral field in a Spin(10) representation."""

    index = float(SPIN10_REPS[rep]["dynkin_index"])
    return np.array([index, index, index], dtype=float)


def drive_sector_fields() -> list[dict[str, Any]]:
    fields: list[dict[str, Any]] = []
    for endpoint in ["AA", "AB", "AC"]:
        rep = SPIN10_REPS["54"]
        fields.append(
            {
                "sector": f"{endpoint}_54_alignment",
                "link_field": f"L_{endpoint}^54",
                "driver_field": f"D_{endpoint}^54",
                "spin10_representation": "54 + 54_D",
                "ps_decomposition": rep["ps_decomposition"],
                "aligned_direction": "(1,1,1) PS singlet",
                "orthogonal_complex_components_per_field": rep["dimension"]
                - rep["aligned_sm_singlets"],
                "mass_pair": "kappa_54 M_G",
                "beta_vector_for_link_driver_pair": full_pair_beta("54").tolist(),
                "projected_beta_l2": float(np.linalg.norm(base.PROJECTOR @ full_pair_beta("54"))),
            }
        )

    rep = SPIN10_REPS["210"]
    fields.append(
        {
            "sector": "AA_210_alignment",
            "link_field": "L_AA^210",
            "driver_field": "D_AA^210",
            "spin10_representation": "210 + 210_D",
            "ps_decomposition": rep["ps_decomposition"],
            "aligned_direction": "weak-volume (1,1,1) PS singlet",
            "orthogonal_complex_components_per_field": rep["dimension"]
            - rep["aligned_sm_singlets"],
            "mass_pair": "kappa_210 M_G",
            "beta_vector_for_link_driver_pair": full_pair_beta("210").tolist(),
            "projected_beta_l2": float(np.linalg.norm(base.PROJECTOR @ full_pair_beta("210"))),
        }
    )

    fields.append(
        {
            "sector": "singlet_ratio_and_norm_drivers",
            "link_field": "L_AA^1, L_AB^1, L_AC^1, L_BC^1 and parallel link coordinates",
            "driver_field": "X_AA^1, X_AA^54, X_AA^210, X_AB^n, X_AB^r, X_AC^n, X_AC^r, X_BC^1",
            "spin10_representation": "Spin(10) singlets",
            "ps_decomposition": ["(1,1,1) only"],
            "aligned_direction": "gauge singlet",
            "orthogonal_complex_components_per_field": 0,
            "mass_pair": "singlet F-term masses",
            "beta_vector_for_link_driver_pair": [0.0, 0.0, 0.0],
            "projected_beta_l2": 0.0,
        }
    )
    fields.append(
        {
            "sector": "relative_source_alignment",
            "link_field": "S_AA^54, S_AB^54, S_AC^54 plus L_e^54",
            "driver_field": "D_AA^54, D_AB^54, D_AC^54, Y_AB^54, Y_AC^54",
            "spin10_representation": "11 complete 54 fields in W_src+W_rel",
            "ps_decomposition": SPIN10_REPS["54"]["ps_decomposition"],
            "aligned_direction": "S_AB parallel S_AA and S_AC parallel S_AA",
            "orthogonal_complex_components_per_field": SPIN10_REPS["54"]["dimension"]
            - SPIN10_REPS["54"]["aligned_sm_singlets"],
            "mass_pair": "eigenvalues of the 11x11 source/link/driver block",
            "beta_vector_for_link_driver_pair": (11.0 * full_field_beta("54")).tolist(),
            "projected_beta_l2": float(
                np.linalg.norm(base.PROJECTOR @ (11.0 * full_field_beta("54")))
            ),
        }
    )
    return fields


def source_link_driver_matrix(source_mass: float, ell: float, kappa: float = 1.0) -> np.ndarray:
    """Basis [S, L, D] for one complete representation component."""

    return np.array(
        [
            [source_mass, 0.0, -kappa * ell],
            [0.0, 0.0, kappa],
            [-kappa * ell, kappa, 0.0],
        ],
        dtype=float,
    )


def relative_source_alignment_matrix(
    source_mass: float = 1.0,
    kappa_drive: float = 1.0,
    kappa_rel: float = 1.0,
    v_ab: float = 1.0,
    v_ac: float = 1.0,
) -> np.ndarray:
    """Combined 54 source/link/driver/relative-driver Hessian.

    Basis:
        S_AA, S_AB, S_AC, L_AA, L_AB, L_AC, D_AA, D_AB, D_AC, Y_AB, Y_AC.

    The same matrix is repeated for every component of a full 54 if the
    coefficients are Spin(10)-scalar.
    """

    card = comp.load_card()
    pars = card["vev_parameters_and_couplings"]
    ells = [float(pars["a54"]), float(pars["eta"]), float(pars["eta"])]
    matrix = np.zeros((11, 11), dtype=float)
    for endpoint, ell in enumerate(ells):
        s_pos = endpoint
        l_pos = 3 + endpoint
        d_pos = 6 + endpoint
        matrix[s_pos, s_pos] = source_mass
        matrix[d_pos, l_pos] = kappa_drive
        matrix[l_pos, d_pos] = kappa_drive
        matrix[d_pos, s_pos] = -kappa_drive * ell
        matrix[s_pos, d_pos] = -kappa_drive * ell

    y_ab = 9
    y_ac = 10
    matrix[y_ab, 1] = matrix[1, y_ab] = kappa_rel
    matrix[y_ab, 0] = matrix[0, y_ab] = -kappa_rel * v_ab
    matrix[y_ac, 2] = matrix[2, y_ac] = kappa_rel
    matrix[y_ac, 0] = matrix[0, y_ac] = -kappa_rel * v_ac
    return 0.5 * (matrix + matrix.T)


def log_sum_abs_eigenvalues(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    eig = np.linalg.eigvalsh(0.5 * (matrix + matrix.T))
    log_sum = float(np.sum(np.log(1.0 / np.maximum(np.abs(eig), 1.0e-14))))
    return eig, log_sum


def relative_source_spectrum() -> dict[str, Any]:
    """Spectrum and threshold of the Y_AB/Y_AC relative-alignment sector."""

    rel_matrix = relative_source_alignment_matrix()
    rel_eig, rel_log = log_sum_abs_eigenvalues(rel_matrix)

    card = comp.load_card()
    ell210 = float(card["vev_parameters_and_couplings"]["b210"])
    matrix210 = source_link_driver_matrix(1.0, ell210)
    eig210, log210 = log_sum_abs_eigenvalues(matrix210)

    delta54 = full_field_beta("54") * rel_log / (2.0 * math.pi)
    delta210 = full_field_beta("210") * log210 / (2.0 * math.pi)
    total = delta54 + delta210
    return {
        "basis_54": [
            "S_AA",
            "S_AB",
            "S_AC",
            "L_AA",
            "L_AB",
            "L_AC",
            "D_AA",
            "D_AB",
            "D_AC",
            "Y_AB",
            "Y_AC",
        ],
        "matrix_54_over_MG": rel_matrix.tolist(),
        "eigenvalues_54_over_MG": rel_eig.tolist(),
        "log_sum_54": rel_log,
        "delta_54": delta54.tolist(),
        "matrix_210_over_MG": matrix210.tolist(),
        "eigenvalues_210_over_MG": eig210.tolist(),
        "log_sum_210": log210,
        "delta_210": delta210.tolist(),
        "delta_total": total.tolist(),
        "universal_delta": float(np.mean(total)),
        "projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ total)),
        "all_eigenvalues_lifted": bool(np.min(np.abs(rel_eig)) > 1.0e-6),
    }


def drive_threshold_vector(kappa54: float, kappa210: float) -> tuple[np.ndarray, float]:
    """Return the one-loop threshold vector and its universal scalar.

    For each complete pair,
        Delta_i = b_i log(M_G/M_pair)/(2 pi).
    """

    log54 = math.log(1.0 / kappa54)
    log210 = math.log(1.0 / kappa210)
    beta = 3.0 * full_pair_beta("54") * log54 + full_pair_beta("210") * log210
    delta = beta / (2.0 * math.pi)
    universal = float(np.mean(delta))
    return delta, universal


def exact_link_metrics() -> dict[str, Any]:
    card = comp.load_card()
    ops = route_a.baseline_operators(card)
    baseline_raw = route_a.spectrum_metrics(ops["H0"], ops)
    ops["baseline_raw_heavy_delta"] = baseline_raw["heavy_threshold_delta"]
    ops["card_heavy_delta"] = card["mediator_heavy_threshold"]["delta_full"]
    metrics = route_a.spectrum_metrics(ops["H0"], ops)
    raw = np.array(metrics["heavy_threshold_delta"], dtype=float)
    calibrated = (
        np.array(ops["card_heavy_delta"], dtype=float)
        + raw
        - np.array(ops["baseline_raw_heavy_delta"], dtype=float)
    )
    metrics["heavy_threshold_delta"] = calibrated.tolist()
    metrics["heavy_threshold_projected_l2"] = float(np.linalg.norm(base.PROJECTOR @ calibrated))
    return metrics


def fixed_spectrum_rge_scan_with_universal(
    metrics: dict[str, Any],
    grid: dict[str, np.ndarray],
    prefactor: float,
    universal_delta: float,
) -> dict[str, Any]:
    """Replay the Route-A fixed-spectrum scan including a universal threshold."""

    alpha_inv = grid["alpha_inv"]
    delta_med = np.array(metrics["heavy_threshold_delta"], dtype=float)
    unwanted = np.array(metrics["unwanted_threshold_delta"], dtype=float)
    fixed_chiral = (
        route_a.BETAS["Sigma_L"] * float(metrics["log_sigma3_eff"])
        + route_a.BETAS["Sigma_8_octet"] * float(metrics["log_sigma8_eff"])
    ) / (2.0 * math.pi)
    fixed = delta_med + unwanted + fixed_chiral + universal_delta * np.ones(3)

    target = (alpha_inv - fixed) @ base.PROJECTOR.T
    mat = base.PROJECTOR @ route_a.BETAS["H_C+H_Cbar"] / (2.0 * math.pi)
    denom = float(mat @ mat)
    log_hc = target @ mat / denom
    hc_threshold = np.outer(log_hc, route_a.BETAS["H_C+H_Cbar"] / (2.0 * math.pi))
    residual_vec = (alpha_inv - fixed - hc_threshold) @ base.PROJECTOR.T
    residual = np.linalg.norm(residual_vec, axis=1)
    alpha_g_inv = np.mean(alpha_inv - fixed - hc_threshold, axis=1)

    lambda_t = np.exp(-log_hc)
    k3 = max(float(metrics["kappa3_eff"]), 1.0e-300)
    k8 = max(float(metrics["kappa8_eff"]), 1.0e-300)
    ratio = k3 / k8
    chi = (ratio - 1.0) / (2.0 + 3.0 * ratio)
    lambda_s = k3 / (1.0 + 2.0 * chi)
    cancellation = min(abs(1.0 + 2.0 * chi), abs(1.0 - 3.0 * chi))
    perturbative_fixed = (
        max(k3, k8, lambda_s, float(np.min(lambda_t))) <= 3.0
        and min(k3, k8, lambda_s, float(np.max(lambda_t))) >= 0.03
        and abs(chi) < 0.49
    )
    alpha_perturbative = alpha_g_inv > 5.0

    mx_min = np.array(
        [
            base.mx_required(float(a), prefactor) if a > 0.0 else float("inf")
            for a in alpha_g_inv
        ]
    )
    c6_gauge = np.array(
        [
            (4.0 * math.pi / a) / (mg * mg) if a > 0.0 else float("inf")
            for a, mg in zip(alpha_g_inv, grid["MG_GeV"])
        ]
    )
    tau_d6 = np.array(
        [
            base.lifetime_from_c6(float(c), prefactor) if math.isfinite(c) else 0.0
            for c in c6_gauge
        ]
    )
    m_t = lambda_t * grid["MG_GeV"]
    loop = (1.0 / 25.0) / (4.0 * math.pi)
    c6_limit = base.c6_max_for_target(prefactor)
    denom_d5 = (grid["yt"] * grid["yb"] / m_t) * loop * base.MWINO_GEV / (
        grid["MSUSY_GeV"] * grid["MSUSY_GeV"]
    )
    filter_required = c6_limit / np.maximum(denom_d5, 1.0e-300)
    c6_d5 = base.TRIPLET_FILTER_TARGET * grid["yt"] * grid["yb"] / m_t * loop * base.MWINO_GEV / (
        grid["MSUSY_GeV"] * grid["MSUSY_GeV"]
    )
    tau_d5 = np.array([base.lifetime_from_c6(float(c), prefactor) for c in c6_d5])

    proton_safe = grid["MG_GeV"] >= mx_min
    d5_safe = base.TRIPLET_FILTER_TARGET <= filter_required
    residual_ok = residual < 1.0e-3
    all_safe = proton_safe & d5_safe & residual_ok & perturbative_fixed & alpha_perturbative
    score = (
        residual * 1000.0
        + np.maximum(0.0, np.log(mx_min / grid["MG_GeV"])) ** 2 * 20.0
        + np.maximum(0.0, np.log(base.TRIPLET_FILTER_TARGET / filter_required)) ** 2 * 20.0
        + np.abs(log_hc)
        + 0.1 * abs(chi)
        + np.maximum(0.0, 5.0 - alpha_g_inv) ** 2
    )
    safe_count = int(np.sum(all_safe))
    pick = int(np.argmin(np.where(all_safe, score, np.inf))) if safe_count else int(np.argmin(score))
    return {
        "safe_count_residual_lt_1e_3": safe_count,
        "picked_index": pick,
        "best_residual_l2": float(residual[pick]),
        "best_MG_GeV": float(grid["MG_GeV"][pick]),
        "best_alphaG_inv": float(alpha_g_inv[pick]),
        "best_log_HC": float(log_hc[pick]),
        "best_M_HC_GeV": float(lambda_t[pick] * grid["MG_GeV"][pick]),
        "best_tau_dim6_years": float(tau_d6[pick]),
        "best_tau_dim5_ST_1e_5_years": float(tau_d5[pick]),
        "best_filter_required": float(filter_required[pick]),
        "proton_safe_at_pick": bool(proton_safe[pick]),
        "d5_safe_at_pick": bool(d5_safe[pick]),
        "residual_ok_at_pick": bool(residual_ok[pick]),
        "alphaG_inv_gt_5_at_pick": bool(alpha_perturbative[pick]),
        "perturbative_fixed_spectrum": bool(perturbative_fixed),
        "chi_from_fixed_kappas": float(chi),
        "lambda_S_from_fixed_kappas": float(lambda_s),
        "cancellation_index": float(cancellation),
    }


def replay_scenarios() -> list[dict[str, Any]]:
    metrics = exact_link_metrics()
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    scenarios = [
        ("canonical_at_MG", 1.0, 1.0),
        ("mild_light_complete_pairs", 0.5, 0.5),
        ("mild_heavy_complete_pairs", 2.0, 2.0),
        ("split_complete_pairs", 0.7, 1.3),
        ("stress_0p3_complete_pairs", 0.3, 0.3),
    ]
    rows: list[dict[str, Any]] = []
    for name, k54, k210 in scenarios:
        delta, universal = drive_threshold_vector(k54, k210)
        rge = fixed_spectrum_rge_scan_with_universal(metrics, grid, prefactor, universal)
        rows.append(
            {
                "scenario": name,
                "kappa54": k54,
                "kappa210": k210,
                "threshold_delta": delta.tolist(),
                "universal_delta": universal,
                "projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ delta)),
                "rge_replay": rge,
                "nonuniversal_replay_needed": bool(np.linalg.norm(base.PROJECTOR @ delta) > 1.0e-12),
                "safe": bool(rge["safe_count_residual_lt_1e_3"] > 0),
            }
        )
    return rows


def common_kappa_bounds(alpha_g_inv_ref: float, tau_ref: float) -> dict[str, float]:
    coeff = 3.0 * float(full_pair_beta("54")[0]) + float(full_pair_beta("210")[0])

    def kappa_for_alpha(alpha_min: float) -> float:
        return math.exp(-max(alpha_g_inv_ref - alpha_min, 0.0) * 2.0 * math.pi / coeff)

    alpha_for_tau_1e34 = alpha_g_inv_ref * math.sqrt(1.0e34 / tau_ref)
    alpha_for_tau_2p4e34 = alpha_g_inv_ref * math.sqrt(2.4e34 / tau_ref)
    return {
        "total_complete_pair_beta_coefficient": coeff,
        "kappa_min_for_alphaG_inv_gt_5": kappa_for_alpha(5.0),
        "kappa_min_for_tau_d6_gt_1e34": kappa_for_alpha(alpha_for_tau_1e34),
        "kappa_min_for_tau_d6_gt_2p4e34": kappa_for_alpha(alpha_for_tau_2p4e34),
        "alphaG_inv_needed_for_tau_d6_gt_1e34": alpha_for_tau_1e34,
        "alphaG_inv_needed_for_tau_d6_gt_2p4e34": alpha_for_tau_2p4e34,
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "scenario",
        "kappa54",
        "kappa210",
        "universal_delta",
        "projected_delta_l2",
        "safe_count",
        "best_alphaG_inv",
        "best_residual_l2",
        "tau_d6",
        "tau_d5",
        "safe",
    ]
    with (OUT / "drive_sector_threshold_replay.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            rge = row["rge_replay"]
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "kappa54": row["kappa54"],
                    "kappa210": row["kappa210"],
                    "universal_delta": row["universal_delta"],
                    "projected_delta_l2": row["projected_delta_l2"],
                    "safe_count": rge["safe_count_residual_lt_1e_3"],
                    "best_alphaG_inv": rge["best_alphaG_inv"],
                    "best_residual_l2": rge["best_residual_l2"],
                    "tau_d6": rge["best_tau_dim6_years"],
                    "tau_d5": rge["best_tau_dim5_ST_1e_5_years"],
                    "safe": row["safe"],
                }
            )


def write_relative_csv(relative: dict[str, Any], replay: dict[str, Any]) -> None:
    with (OUT / "drive_source_relative_eigenvalues.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["sector", "index", "eigenvalue_over_MG"])
        writer.writeheader()
        for idx, value in enumerate(relative["eigenvalues_54_over_MG"]):
            writer.writerow(
                {
                    "sector": "54_source_relative",
                    "index": idx,
                    "eigenvalue_over_MG": value,
                }
            )
        for idx, value in enumerate(relative["eigenvalues_210_over_MG"]):
            writer.writerow(
                {
                    "sector": "210_source_alignment",
                    "index": idx,
                    "eigenvalue_over_MG": value,
                }
            )

    with (OUT / "drive_source_relative_replay.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scenario",
                "universal_delta",
                "projected_delta_l2",
                "safe_count",
                "best_alphaG_inv",
                "best_residual_l2",
                "tau_d6",
                "tau_d5",
                "safe",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "scenario": "canonical_with_W_rel_source_block",
                "universal_delta": relative["universal_delta"],
                "projected_delta_l2": relative["projected_delta_l2"],
                "safe_count": replay["safe_count_residual_lt_1e_3"],
                "best_alphaG_inv": replay["best_alphaG_inv"],
                "best_residual_l2": replay["best_residual_l2"],
                "tau_d6": replay["best_tau_dim6_years"],
                "tau_d5": replay["best_tau_dim5_ST_1e_5_years"],
                "safe": replay["safe_count_residual_lt_1e_3"] > 0,
            }
        )


def write_report(payload: dict[str, Any]) -> None:
    fields = payload["drive_sector_fields"]
    rows = payload["threshold_replay"]
    bounds = payload["common_kappa_bounds"]
    relative = payload["relative_source_alignment_spectrum"]
    relative_replay = payload["relative_source_alignment_replay"]
    lines: list[str] = []
    lines.append("# Drive-sector representation and threshold audit")
    lines.append("")
    lines.append("No web lookup was used.  The link-driving card is promoted to explicit")
    lines.append("Spin(10) representation content and its orthogonal alignment modes are")
    lines.append("checked as heavy chiral thresholds.")
    lines.append("")
    lines.append("## Representation content")
    lines.append("")
    lines.append("| sector | Spin(10) content | aligned direction | orthogonal components | pair beta vector |")
    lines.append("|---|---|---|---:|---:|")
    for item in fields:
        lines.append(
            f"| `{item['sector']}` | `{item['spin10_representation']}` | "
            f"{item['aligned_direction']} | {item['orthogonal_complex_components_per_field']} | "
            f"`{item['beta_vector_for_link_driver_pair']}` |"
        )
    lines.append("")
    lines.append("The relevant branchings are")
    lines.append("")
    lines.append("```text")
    lines.append("54  -> " + " + ".join(SPIN10_REPS["54"]["ps_decomposition"]))
    lines.append("210 -> " + " + ".join(SPIN10_REPS["210"]["ps_decomposition"]))
    lines.append("```")
    lines.append("")
    lines.append("For every endpoint, the orthogonal Hessian has the local form")
    lines.append("")
    lines.append("```text")
    lines.append("W_perp = kappa M_G <D_perp, delta L_perp>,")
    lines.append("M_perp/M_G = [[0, kappa], [kappa, 0]],")
    lines.append("physical masses = |kappa| M_G.")
    lines.append("```")
    lines.append("")
    lines.append("Thus the orthogonal modes are not extra light states unless the drive")
    lines.append("couplings kappa are deliberately small.")
    lines.append("")
    lines.append("## Relative source-alignment drivers")
    lines.append("")
    lines.append("The endpoint source copies require two additional full `54` drivers,")
    lines.append("`Y_AB^54` and `Y_AC^54`, with")
    lines.append("")
    lines.append("```text")
    lines.append("W_rel = <Y_AB, S_AB - v_AB S_AA> + <Y_AC, S_AC - v_AC S_AA>.")
    lines.append("```")
    lines.append("")
    lines.append("For each non-singlet component of `54`, the combined source/link/driver")
    lines.append("basis is")
    lines.append("")
    lines.append("```text")
    lines.append(", ".join(relative["basis_54"]))
    lines.append("```")
    lines.append("")
    lines.append("and the same 11x11 Hessian is repeated across every `54` fragment.")
    lines.append("Numerically,")
    lines.append("")
    lines.append("```text")
    lines.append(
        f"min |lambda_54|      = {min(abs(x) for x in relative['eigenvalues_54_over_MG']):.3e}"
    )
    lines.append(f"log sum 54           = {relative['log_sum_54']:.6e}")
    lines.append(f"log sum 210          = {relative['log_sum_210']:.6e}")
    lines.append(f"P Delta_source+rel   = {relative['projected_delta_l2']:.3e}")
    lines.append(f"universal Delta      = {relative['universal_delta']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("Thus the new `Y_AB,Y_AC` non-singlet modes are lifted and form complete")
    lines.append("54 multiplets.  They can shift only the universal value of `alpha_G^{-1}`.")
    lines.append("The cached RGE/proton replay with this universal source block gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"safe points = {relative_replay['safe_count_residual_lt_1e_3']}")
    lines.append(f"alphaG^-1   = {relative_replay['best_alphaG_inv']:.6g}")
    lines.append(f"residual    = {relative_replay['best_residual_l2']:.3e}")
    lines.append(f"tau_d6      = {relative_replay['best_tau_dim6_years']:.3e} yr")
    lines.append(f"tau_d5      = {relative_replay['best_tau_dim5_ST_1e_5_years']:.3e} yr")
    lines.append("```")
    lines.append("")
    lines.append("## Threshold theorem")
    lines.append("")
    lines.append("A complete 54 link-driver pair contributes `(24,24,24)` and a complete")
    lines.append("210 link-driver pair contributes `(112,112,112)` to one-loop chiral beta")
    lines.append("coefficients.  With three 54 pairs and one 210 pair,")
    lines.append("")
    lines.append("```text")
    lines.append("Delta_drive = [72 log(1/kappa54) + 112 log(1/kappa210)]/(2 pi) * (1,1,1).")
    lines.append("P Delta_drive = 0.")
    lines.append("```")
    lines.append("")
    lines.append("Therefore there is no non-universal threshold to add to the matching scan.")
    lines.append("Only the universal alpha_G shift needs a proton-decay sanity check if the")
    lines.append("complete drive pairs are far below M_G.")
    lines.append("")
    lines.append("## RGE/proton replay")
    lines.append("")
    lines.append("| scenario | k54 | k210 | universal Delta | projected l2 | safe points | alphaG^-1 | residual | tau_d6 [yr] | tau_d5 [yr] |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        rge = row["rge_replay"]
        lines.append(
            f"| `{row['scenario']}` | {row['kappa54']:.3g} | {row['kappa210']:.3g} | "
            f"{row['universal_delta']:.6e} | {row['projected_delta_l2']:.3e} | "
            f"{rge['safe_count_residual_lt_1e_3']} | {rge['best_alphaG_inv']:.6g} | "
            f"{rge['best_residual_l2']:.3e} | {rge['best_tau_dim6_years']:.3e} | "
            f"{rge['best_tau_dim5_ST_1e_5_years']:.3e} |"
        )
    lines.append("")
    lines.append("For a common drive-pair mass `kappa M_G`, the total complete-pair beta")
    lines.append(f"coefficient is `{bounds['total_complete_pair_beta_coefficient']:.0f}`.")
    lines.append("The approximate lower bounds are")
    lines.append("")
    lines.append("```text")
    lines.append(f"kappa > {bounds['kappa_min_for_alphaG_inv_gt_5']:.3f} for alphaG^-1 > 5")
    lines.append(f"kappa > {bounds['kappa_min_for_tau_d6_gt_1e34']:.3f} for tau_d6 > 1e34 yr")
    lines.append(f"kappa > {bounds['kappa_min_for_tau_d6_gt_2p4e34']:.3f} for tau_d6 > 2.4e34 yr")
    lines.append("```")
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append("At the canonical mass assignment `kappa54=kappa210=1`, the drive sector")
    lines.append("adds no non-universal threshold at all.  With `W_rel`, the relative")
    lines.append("source block contributes only a complete-multiplet universal shift.")
    lines.append("Even if the complete drive pairs are moved,")
    lines.append("their non-singlet content is complete and hence projected-threshold silent.")
    lines.append("The only remaining restriction is to avoid a very light complete drive")
    lines.append("sector, roughly `kappa < 0.33`, which would lower alphaG enough to stress")
    lines.append("the dimension-6 proton bound.")
    (OUT / "drive_sector_spectrum_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fields = drive_sector_fields()
    rows = replay_scenarios()
    relative = relative_source_spectrum()
    metrics = exact_link_metrics()
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    relative_replay = fixed_spectrum_rge_scan_with_universal(
        metrics, grid, prefactor, float(relative["universal_delta"])
    )
    canonical = next(row for row in rows if row["scenario"] == "canonical_at_MG")
    bounds = common_kappa_bounds(
        canonical["rge_replay"]["best_alphaG_inv"],
        canonical["rge_replay"]["best_tau_dim6_years"],
    )
    payload = {
        "note": "No web lookup used. Explicit Spin(10) representation audit of the link-driving sector.",
        "superpotential_lift": {
            "54_alignment": "sum_{e=AA,AB,AC} kappa_e^54 M_G <D_e^54, L_e^54 - ell_e S_54>",
            "210_alignment": "kappa_AA^210 M_G <D_AA^210, L_AA^210 - ell_AA^210 Omega_210>",
            "singlet_drivers": "X fields enforce norms and singlet/54 ratios.",
        },
        "drive_sector_fields": fields,
        "threshold_formula": {
            "one_54_link_driver_pair_beta": full_pair_beta("54").tolist(),
            "one_210_link_driver_pair_beta": full_pair_beta("210").tolist(),
            "three_54_plus_one_210_common_beta": (
                3.0 * full_pair_beta("54") + full_pair_beta("210")
            ).tolist(),
            "projected_beta_l2": float(
                np.linalg.norm(base.PROJECTOR @ (3.0 * full_pair_beta("54") + full_pair_beta("210")))
            ),
            "formula": "Delta = [72 log(1/kappa54)+112 log(1/kappa210)]/(2pi) * (1,1,1)",
        },
        "threshold_replay": rows,
        "relative_source_alignment_spectrum": relative,
        "relative_source_alignment_replay": relative_replay,
        "common_kappa_bounds": bounds,
        "passes": {
            "all_drive_pair_betas_project_universal": all(
                item["projected_beta_l2"] < 1.0e-12 for item in fields
            ),
            "canonical_projected_threshold_zero": canonical["projected_delta_l2"] < 1.0e-12,
            "canonical_safe": bool(canonical["safe"]),
            "mild_light_safe": bool(
                next(row for row in rows if row["scenario"] == "mild_light_complete_pairs")["safe"]
            ),
            "relative_source_projected_threshold_zero": relative["projected_delta_l2"] < 1.0e-12,
            "relative_source_eigenvalues_lifted": relative["all_eigenvalues_lifted"],
            "relative_source_replay_safe": relative_replay["safe_count_residual_lt_1e_3"] > 0,
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "drive_sector_spectrum_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_relative_csv(relative, relative_replay)
    write_report(payload)

    print("Drive-sector representation audit")
    print(f"  canonical projected threshold: {canonical['projected_delta_l2']:.3e}")
    print(f"  canonical safe points: {canonical['rge_replay']['safe_count_residual_lt_1e_3']}")
    print(f"  canonical tau_d6: {canonical['rge_replay']['best_tau_dim6_years']:.3e} yr")
    print(
        "  common kappa min for tau_d6>2.4e34: "
        f"{bounds['kappa_min_for_tau_d6_gt_2p4e34']:.3f}"
    )
    print(f"  W_rel source projected threshold: {relative['projected_delta_l2']:.3e}")
    print(f"  W_rel source safe points: {relative_replay['safe_count_residual_lt_1e_3']}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
