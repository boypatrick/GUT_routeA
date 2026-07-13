#!/usr/bin/env python3
"""Dynamic S_54/Omega_210 source mixing audit.

No web lookup is used.  This script puts the source fields S_54 and Omega_210
back into the same local superpotential as the link/driver sector and checks
whether their mixing spoils the complete-multiplet threshold structure.

For each non-singlet component of a representation R, the local quadratic
superpotential is

    W_R = 1/2 m_R S_R^2 + sum_e kappa_e D_e (L_e - ell_e S_R).

If m_R is a Spin(10)-scalar mass, every component of R has the same finite
matrix, so each eigenvalue is a complete Spin(10) multiplet and the projected
gauge-threshold vector vanishes.  The script then adds controlled PS-fragment
splittings to m_R to quantify how strongly a non-scalar source Hessian would
damage the unification/proton branch.
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

import audit_drive_sector_spectrum as drive  # noqa: E402
import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "dynamic_source_mixing"


FRAGMENTS_54 = {
    "(1,1,1)": {"beta": np.array([0.0, 0.0, 0.0]), "split": 0.0},
    "(20',1,1)": {"beta": np.array([16.0 / 5.0, 0.0, 8.0]), "split": 1.0},
    "(6,2,2)": {"beta": np.array([26.0 / 5.0, 6.0, 4.0]), "split": -1.0},
    "(1,3,3)": {"beta": np.array([18.0 / 5.0, 6.0, 0.0]), "split": 0.5},
}

FRAGMENTS_210 = {
    "(1,1,1)": {"beta": np.array([0.0, 0.0, 0.0]), "split": 0.0},
    "(15,1,1)": {"beta": np.array([8.0 / 5.0, 0.0, 4.0]), "split": 1.0},
    "(15,3,1)": {"beta": np.array([24.0 / 5.0, 30.0, 12.0]), "split": -0.5},
    "(15,1,3)": {"beta": np.array([114.0 / 5.0, 0.0, 12.0]), "split": 0.7},
    "(6,2,2)": {"beta": np.array([26.0 / 5.0, 6.0, 4.0]), "split": -1.0},
    "(10,2,2)+(10bar,2,2)": {
        "beta": np.array([108.0 / 5.0, 20.0, 24.0]),
        "split": 0.3,
    },
}


def check_fragment_tables() -> dict[str, Any]:
    beta54 = sum((item["beta"] for item in FRAGMENTS_54.values()), np.zeros(3))
    beta210 = sum((item["beta"] for item in FRAGMENTS_210.values()), np.zeros(3))
    return {
        "sum_54_fragment_betas": beta54.tolist(),
        "sum_210_fragment_betas": beta210.tolist(),
        "expected_54_beta": [12.0, 12.0, 12.0],
        "expected_210_beta": [56.0, 56.0, 56.0],
        "54_projected_sum_l2": float(np.linalg.norm(base.PROJECTOR @ beta54)),
        "210_projected_sum_l2": float(np.linalg.norm(base.PROJECTOR @ beta210)),
    }


def source_link_driver_matrix(source_mass: float, ells: list[float], kappas: list[float]) -> np.ndarray:
    """Mass matrix for [S, L_1..L_n, D_1..D_n]."""

    n = len(ells)
    m = np.zeros((1 + 2 * n, 1 + 2 * n), dtype=float)
    m[0, 0] = source_mass
    for idx, (ell, kappa) in enumerate(zip(ells, kappas)):
        l_pos = 1 + idx
        d_pos = 1 + n + idx
        m[d_pos, l_pos] = kappa
        m[l_pos, d_pos] = kappa
        m[d_pos, 0] = -kappa * ell
        m[0, d_pos] = -kappa * ell
    return 0.5 * (m + m.T)


def eig_abs_logs(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    vals = np.linalg.eigvalsh(0.5 * (matrix + matrix.T))
    abs_vals = np.maximum(np.abs(vals), 1.0e-14)
    return vals, float(np.sum(np.log(1.0 / abs_vals)))


def card_ells() -> dict[str, Any]:
    card = comp.load_card()
    pars = card["vev_parameters_and_couplings"]
    return {
        "source_card": str(comp.CARD),
        "R": float(card["R"]),
        "ells_54": [
            float(pars["a54"]),
            float(pars["eta"]),
            float(pars["eta"]),
        ],
        "ells_54_labels": ["AA_54/a54", "AB_54/eta", "AC_54/eta"],
        "ells_210": [float(pars["b210"])],
        "ells_210_labels": ["AA_210/b210"],
        "kappas_54": [1.0, 1.0, 1.0],
        "kappas_210": [1.0],
    }


def threshold_from_fragments(
    fragments: dict[str, dict[str, Any]],
    source_mass: float,
    epsilon: float,
    ells: list[float],
    kappas: list[float],
) -> tuple[np.ndarray, dict[str, Any]]:
    delta = np.zeros(3, dtype=float)
    rows: dict[str, Any] = {}
    for name, item in fragments.items():
        mass = source_mass * (1.0 + epsilon * float(item["split"]))
        matrix = source_link_driver_matrix(mass, ells, kappas)
        evals, log_sum = eig_abs_logs(matrix)
        beta = np.array(item["beta"], dtype=float)
        delta += beta * log_sum / (2.0 * math.pi)
        rows[name] = {
            "source_mass": mass,
            "split_coefficient": float(item["split"]),
            "eigenvalues_over_MG": evals.tolist(),
            "log_sum": log_sum,
            "beta": beta.tolist(),
        }
    return delta, rows


def source_mixing_threshold(epsilon54: float, epsilon210: float) -> dict[str, Any]:
    pars = card_ells()
    delta54, rows54 = threshold_from_fragments(
        FRAGMENTS_54,
        source_mass=1.0,
        epsilon=epsilon54,
        ells=pars["ells_54"],
        kappas=pars["kappas_54"],
    )
    delta210, rows210 = threshold_from_fragments(
        FRAGMENTS_210,
        source_mass=1.0,
        epsilon=epsilon210,
        ells=pars["ells_210"],
        kappas=pars["kappas_210"],
    )
    total = delta54 + delta210

    # Degenerate-source theorem diagnostic: all fragment eigenvalue arrays must
    # match inside a given Spin(10) representation when epsilon=0.
    eig54 = [np.array(row["eigenvalues_over_MG"]) for row in rows54.values()]
    eig210 = [np.array(row["eigenvalues_over_MG"]) for row in rows210.values()]
    spread54 = max(float(np.max(np.abs(x - eig54[0]))) for x in eig54)
    spread210 = max(float(np.max(np.abs(x - eig210[0]))) for x in eig210)
    return {
        "epsilon54": epsilon54,
        "epsilon210": epsilon210,
        "delta54": delta54.tolist(),
        "delta210": delta210.tolist(),
        "delta_total": total.tolist(),
        "projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ total)),
        "universal_delta": float(np.mean(total)),
        "fragment_eigenvalue_spread_54": spread54,
        "fragment_eigenvalue_spread_210": spread210,
        "fragments54": rows54,
        "fragments210": rows210,
    }


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


def fixed_spectrum_rge_scan_with_extra_delta(
    metrics: dict[str, Any],
    grid: dict[str, np.ndarray],
    prefactor: float,
    extra_delta: np.ndarray,
) -> dict[str, Any]:
    alpha_inv = grid["alpha_inv"]
    delta_med = np.array(metrics["heavy_threshold_delta"], dtype=float)
    unwanted = np.array(metrics["unwanted_threshold_delta"], dtype=float)
    fixed_chiral = (
        route_a.BETAS["Sigma_L"] * float(metrics["log_sigma3_eff"])
        + route_a.BETAS["Sigma_8_octet"] * float(metrics["log_sigma8_eff"])
    ) / (2.0 * math.pi)
    fixed = delta_med + unwanted + fixed_chiral + extra_delta

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
    }


def replay_scenarios() -> list[dict[str, Any]]:
    metrics = exact_link_metrics()
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    scenarios = [
        ("degenerate_dynamic_sources", 0.0, 0.0),
        ("tiny_fragment_split_1e-4", 1.0e-4, 1.0e-4),
        ("small_fragment_split_1e-3", 1.0e-3, 1.0e-3),
        ("stress_fragment_split_1e-2", 1.0e-2, 1.0e-2),
        ("asymmetric_54_only_1e-3", 1.0e-3, 0.0),
        ("asymmetric_210_only_1e-3", 0.0, 1.0e-3),
    ]
    rows: list[dict[str, Any]] = []
    for name, eps54, eps210 in scenarios:
        threshold = source_mixing_threshold(eps54, eps210)
        rge = fixed_spectrum_rge_scan_with_extra_delta(
            metrics,
            grid,
            prefactor,
            np.array(threshold["delta_total"], dtype=float),
        )
        rows.append(
            {
                "scenario": name,
                **threshold,
                "rge_replay": rge,
                "safe": bool(rge["safe_count_residual_lt_1e_3"] > 0),
            }
        )
    return rows


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "scenario",
        "epsilon54",
        "epsilon210",
        "projected_delta_l2",
        "universal_delta",
        "fragment_eigenvalue_spread_54",
        "fragment_eigenvalue_spread_210",
        "safe_count",
        "best_alphaG_inv",
        "best_residual_l2",
        "tau_d6",
        "tau_d5",
        "safe",
    ]
    with (OUT / "dynamic_source_mixing_replay.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            rge = row["rge_replay"]
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "epsilon54": row["epsilon54"],
                    "epsilon210": row["epsilon210"],
                    "projected_delta_l2": row["projected_delta_l2"],
                    "universal_delta": row["universal_delta"],
                    "fragment_eigenvalue_spread_54": row["fragment_eigenvalue_spread_54"],
                    "fragment_eigenvalue_spread_210": row["fragment_eigenvalue_spread_210"],
                    "safe_count": rge["safe_count_residual_lt_1e_3"],
                    "best_alphaG_inv": rge["best_alphaG_inv"],
                    "best_residual_l2": rge["best_residual_l2"],
                    "tau_d6": rge["best_tau_dim6_years"],
                    "tau_d5": rge["best_tau_dim5_ST_1e_5_years"],
                    "safe": row["safe"],
                }
            )


def write_report(payload: dict[str, Any]) -> None:
    rows = payload["threshold_replay"]
    deg = next(row for row in rows if row["scenario"] == "degenerate_dynamic_sources")
    lines: list[str] = []
    lines.append("# Dynamic S54/Omega210 source-mixing audit")
    lines.append("")
    lines.append("No web lookup was used.  This pass puts `S_54` and `Omega_210`")
    lines.append("back into the same local superpotential as the link/driver sector.")
    lines.append("")
    lines.append("## Quadratic source-link-driver block")
    lines.append("")
    lines.append("For one Spin(10) component of a representation `R`,")
    lines.append("")
    lines.append("```text")
    lines.append("W_R = 1/2 m_R S_R^2 + sum_e kappa_e D_e (L_e - ell_e S_R).")
    lines.append("```")
    lines.append("")
    lines.append("For `54`, `e=(AA,AB,AC)` and")
    lines.append(f"`ell={payload['link_parameters']['ells_54']}`.")
    lines.append("For `210`, `e=AA` and")
    lines.append(f"`ell={payload['link_parameters']['ells_210']}`.")
    lines.append("")
    lines.append("If `m_R` is a Spin(10)-scalar, the same finite matrix is repeated for")
    lines.append("every component.  Each eigenvalue is then a complete Spin(10) multiplet.")
    lines.append("")
    lines.append("## Fragment beta checks")
    lines.append("")
    table = payload["fragment_table_check"]
    lines.append("```text")
    lines.append(f"sum beta(54 fragments)  = {table['sum_54_fragment_betas']}")
    lines.append(f"sum beta(210 fragments) = {table['sum_210_fragment_betas']}")
    lines.append(f"P sum beta(54) l2       = {table['54_projected_sum_l2']:.3e}")
    lines.append(f"P sum beta(210) l2      = {table['210_projected_sum_l2']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Degenerate-source theorem check")
    lines.append("")
    lines.append("```text")
    lines.append(f"54 fragment eigenvalue spread  = {deg['fragment_eigenvalue_spread_54']:.3e}")
    lines.append(f"210 fragment eigenvalue spread = {deg['fragment_eigenvalue_spread_210']:.3e}")
    lines.append(f"projected threshold l2         = {deg['projected_delta_l2']:.3e}")
    lines.append(f"universal Delta                = {deg['universal_delta']:.6e}")
    lines.append(f"safe points                    = {deg['rge_replay']['safe_count_residual_lt_1e_3']}")
    lines.append(f"alphaG_inv                     = {deg['rge_replay']['best_alphaG_inv']:.6g}")
    lines.append(f"tau_d6                         = {deg['rge_replay']['best_tau_dim6_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("## Fragment-split stress scan")
    lines.append("")
    lines.append("| scenario | eps54 | eps210 | projected l2 | safe points | alphaG^-1 | residual | tau_d6 [yr] |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        rge = row["rge_replay"]
        lines.append(
            f"| `{row['scenario']}` | {row['epsilon54']:.1e} | {row['epsilon210']:.1e} | "
            f"{row['projected_delta_l2']:.3e} | {rge['safe_count_residual_lt_1e_3']} | "
            f"{rge['best_alphaG_inv']:.6g} | {rge['best_residual_l2']:.3e} | "
            f"{rge['best_tau_dim6_years']:.3e} |"
        )
    lines.append("")
    lines.append("The stress split is a diagnostic deformation of the source Hessian, not a")
    lines.append("new fitted model.  It shows that scalar Spin(10)-complete source masses")
    lines.append("preserve the complete-pair structure exactly; non-scalar PS-fragment masses")
    lines.append("generate a non-universal threshold at first order in the splitting.")
    (OUT / "dynamic_source_mixing_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = replay_scenarios()
    payload = {
        "note": "No web lookup used. Dynamic S_54/Omega_210 source-link-driver mixing audit.",
        "superpotential": {
            "source_link_driver_block": "1/2 m_R S_R^2 + sum_e kappa_e D_e(L_e-ell_e S_R)",
            "54_block": "R=54, e=(AA,AB,AC)",
            "210_block": "R=210, e=AA",
        },
        "link_parameters": card_ells(),
        "fragment_table_check": check_fragment_tables(),
        "threshold_replay": rows,
        "passes": {
            "fragment_tables_sum_to_complete_multiplets": (
                check_fragment_tables()["54_projected_sum_l2"] < 1.0e-12
                and check_fragment_tables()["210_projected_sum_l2"] < 1.0e-12
            ),
            "degenerate_54_eigenvalues_identical": next(
                row for row in rows if row["scenario"] == "degenerate_dynamic_sources"
            )["fragment_eigenvalue_spread_54"]
            < 1.0e-12,
            "degenerate_210_eigenvalues_identical": next(
                row for row in rows if row["scenario"] == "degenerate_dynamic_sources"
            )["fragment_eigenvalue_spread_210"]
            < 1.0e-12,
            "degenerate_projected_threshold_zero": next(
                row for row in rows if row["scenario"] == "degenerate_dynamic_sources"
            )["projected_delta_l2"]
            < 1.0e-12,
            "degenerate_branch_safe": next(
                row for row in rows if row["scenario"] == "degenerate_dynamic_sources"
            )["safe"],
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "dynamic_source_mixing_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(payload)

    deg = next(row for row in rows if row["scenario"] == "degenerate_dynamic_sources")
    print("Dynamic S54/Omega210 source mixing audit")
    print(f"  degenerate projected threshold: {deg['projected_delta_l2']:.3e}")
    print(f"  degenerate safe points: {deg['rge_replay']['safe_count_residual_lt_1e_3']}")
    print(f"  degenerate alphaG_inv: {deg['rge_replay']['best_alphaG_inv']:.6g}")
    print(f"  degenerate tau_d6: {deg['rge_replay']['best_tau_dim6_years']:.3e} yr")
    stress = next(row for row in rows if row["scenario"] == "small_fragment_split_1e-3")
    print(f"  eps=1e-3 projected threshold: {stress['projected_delta_l2']:.3e}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
