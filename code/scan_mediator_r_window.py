#!/usr/bin/env python3
"""Scan the finite-mediator R window against the two-loop RGE/proton fit.

No web lookup is used.  This script promotes the fixed R=50 mediator benchmark
to a one-parameter superpotential-spectrum scan.  For each mediator ratio R it
solves the exact 3-by-3 mediator mass matrices, computes the finite heavy
threshold vector, and reruns the cached two-loop gauge/Yukawa/proton matching.
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

import construct_uv_projector_mediator as uv  # noqa: E402
import scan_mediator_threshold_rge as med  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "mediator_r_window"
R_GRID = [5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0]


def load_targets() -> tuple[dict[str, float], dict[str, float]]:
    source = json.loads(uv.IN.read_text(encoding="utf-8"))
    target = source["scenarios"]["target"]
    initial = source["scenarios"]["projector_lifted_solution"]["coefficients"]
    targets = {
        "Sigma_L": float(target["kappa_Sigma3"]),
        "Sigma_R": 1.0,
        "Sigma8_block": float(target["kappa_Sigma8"]),
        "X_622": 1.0,
    }
    return targets, initial


def scan_cached_rge(delta_med: np.ndarray) -> dict[str, Any]:
    prefactor = base.proton_prefactor()
    rows: list[dict[str, float | bool]] = []
    best: dict[str, float | bool] | None = None
    with med.BASE_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for source_row in reader:
            row = med.analyze_cached_row_with_mediator(source_row, prefactor, delta_med)
            rows.append(row)
            if best is None or float(row["score"]) < float(best["score"]):
                best = row
    assert best is not None

    safe_rows = [
        row
        for row in rows
        if bool(row["triplet_filter_safe"])
        and float(row["rho_X"]) <= 1.000001
        and bool(row["perturbative_superpotential"])
    ]
    factor2_safe_rows = [row for row in safe_rows if bool(row["single_scale_factor2"])]
    best_safe = min(safe_rows, key=lambda row: float(row["score"]), default=None)
    return {
        "total_points": len(rows),
        "safe_points": len(safe_rows),
        "safe_single_scale_factor2_points": len(factor2_safe_rows),
        "best": best,
        "best_safe": best_safe,
    }


def row_from_solution(r_med: float, solution: dict[str, Any], threshold: dict[str, Any], rge: dict[str, Any]) -> dict[str, Any]:
    best = rge["best_safe"] if rge["best_safe"] is not None else rge["best"]
    delta = np.array(threshold["delta_full"], dtype=float)
    projected = np.array(threshold["projected_delta"], dtype=float)
    return {
        "R": float(r_med),
        "eta": float(solution["parameters"]["eta"]),
        "eta_over_sqrt_R": float(solution["parameters"]["eta"] / math.sqrt(r_med)),
        "schur_delta_X": float(solution["parameters"]["schur_delta_X"]),
        "light_fit_residual_l2": float(solution["residual_l2"]),
        "delta_med_1": float(delta[0]),
        "delta_med_2": float(delta[1]),
        "delta_med_3": float(delta[2]),
        "projected_delta_1": float(projected[0]),
        "projected_delta_2": float(projected[1]),
        "projected_delta_3": float(projected[2]),
        "projected_l2": float(threshold["projected_l2"]),
        "R_times_projected_l2": float(r_med * threshold["projected_l2"]),
        "total_points": int(rge["total_points"]),
        "safe_points": int(rge["safe_points"]),
        "safe_single_scale_factor2_points": int(rge["safe_single_scale_factor2_points"]),
        "best_is_safe": bool(rge["best_safe"] is not None),
        "tan_beta": float(best["tan_beta"]),
        "MSUSY_GeV": float(best["MSUSY_GeV"]),
        "MG_GeV": float(best["MG_GeV"]),
        "alphaG_inv": float(best["alphaG_inv"]),
        "lambda_T": float(best["lambda_T"]),
        "lambda_S": float(best["lambda_S"]),
        "chi": float(best["chi"]),
        "kappa_3": float(best["kappa_3"]),
        "kappa_8": float(best["kappa_8"]),
        "cancellation_index": float(best["cancellation_index"]),
        "M_HC_GeV": float(best["M_HC_GeV"]),
        "M_Sigma3_GeV": float(best["M_Sigma3_GeV"]),
        "M_Sigma8_GeV": float(best["M_Sigma8_GeV"]),
        "log_HC": float(best["log_HC"]),
        "log_Sigma3": float(best["log_Sigma3"]),
        "log_Sigma8": float(best["log_Sigma8"]),
        "residual_l2_after_mediator": float(best["residual_l2_after_mediator"]),
        "tau_dim6_years": float(best["tau_dim6_years"]),
        "tau_dim5_target_filter_years": float(best["tau_dim5_target_filter_years"]),
    }


def scan() -> dict[str, Any]:
    targets, initial = load_targets()
    rows = []
    detailed = {}
    for r_med in R_GRID:
        solution = uv.solve_finite_r(targets, r_med, initial)
        threshold = uv.heavy_threshold_residual(solution)
        rge = scan_cached_rge(np.array(threshold["delta_full"], dtype=float))
        row = row_from_solution(r_med, solution, threshold, rge)
        rows.append(row)
        detailed[str(r_med)] = {
            "solution": solution,
            "heavy_threshold": threshold,
            "rge_diagnostics": {
                "total_points": rge["total_points"],
                "safe_points": rge["safe_points"],
                "safe_single_scale_factor2_points": rge["safe_single_scale_factor2_points"],
                "best": rge["best"],
                "best_safe": rge["best_safe"],
            },
        }

    tail = [row for row in rows if row["R"] >= 20.0]
    log_r = np.log([row["R"] for row in tail])
    log_res = np.log([row["projected_l2"] for row in tail])
    slope, intercept = np.polyfit(log_r, log_res, 1)
    best_safe_rows = [row for row in rows if row["best_is_safe"]]
    best = min(best_safe_rows, key=lambda row: (row["projected_l2"], row["cancellation_index"]))

    return {
        "note": "No web lookup used.  This scans the finite mediator R window and reuses the cached two-loop RGE/Yukawa grid.",
        "R_grid": R_GRID,
        "targets": targets,
        "rows": rows,
        "detailed": detailed,
        "decoupling_fit_R_ge_20": {
            "slope_log_projected_l2_vs_log_R": float(slope),
            "intercept": float(intercept),
            "expected_slope": -1.0,
        },
        "best_by_smallest_nonuniversal_threshold": best,
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    with (OUT / "mediator_r_window_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    rows = summary["rows"]
    fit = summary["decoupling_fit_R_ge_20"]
    best = summary["best_by_smallest_nonuniversal_threshold"]
    lines: list[str] = []
    lines.append("# Mediator R-window scan")
    lines.append("")
    lines.append("No web lookup was used.  For each `R=M_med/MG`, the exact mediator")
    lines.append("mass matrices are solved, the finite heavy threshold vector is computed,")
    lines.append("and the cached two-loop RGE/Yukawa/proton scan is replayed.")
    lines.append("")
    lines.append("## Schur-fixed decoupling test")
    lines.append("")
    lines.append("The scan keeps the light projector spectrum fixed by solving the exact")
    lines.append("3-by-3 mediator mass matrix at each `R`.  Numerically this realizes")
    lines.append("approximately fixed `50 eta^2/(9R)`, hence `eta ~ sqrt(R)`.  The")
    lines.append("non-universal heavy threshold should decay as `1/R` because the leading")
    lines.append("heavy states are almost complete `45` blocks and their common threshold is")
    lines.append("absorbed into `alpha_G^-1`.")
    lines.append("")
    lines.append("```text")
    lines.append(
        "fit over R >= 20: log ||P Delta_med|| = "
        f"{fit['intercept']:+.6f} {fit['slope_log_projected_l2_vs_log_R']:+.6f} log R"
    )
    lines.append("expected slope: -1")
    lines.append("```")
    lines.append("")
    lines.append("| R | eta/sqrt(R) | ||P Delta_med|| | R ||P Delta_med|| | safe points | alphaG^-1 | tau_d6 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['R']:.1f} | {row['eta_over_sqrt_R']:.9f} | "
            f"{row['projected_l2']:.6e} | {row['R_times_projected_l2']:.6e} | "
            f"{row['safe_points']} | {row['alphaG_inv']:.6f} | {row['tau_dim6_years']:.6e} |"
        )
    lines.append("")
    lines.append("## Best decoupled benchmark")
    lines.append("")
    lines.append("Choosing the smallest non-universal mediator threshold in the displayed")
    lines.append("window gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"R = {best['R']:.1f}")
    lines.append(f"eta = {best['eta']:.9f}")
    lines.append(f"||P Delta_med||_2 = {best['projected_l2']:.6e}")
    lines.append(f"safe cached two-loop points = {best['safe_points']}")
    lines.append(f"MG = {best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {best['alphaG_inv']:.6f}")
    lines.append(f"M_Sigma3 = {best['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"M_Sigma8 = {best['M_Sigma8_GeV']:.6e} GeV")
    lines.append(f"tau_dim6 = {best['tau_dim6_years']:.6e} yr")
    lines.append(f"tau_dim5(S_T=1e-5) = {best['tau_dim5_target_filter_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("No factor-two single-scale point appears anywhere in this R-window.  The")
    lines.append("viable branch remains the intermediate-Sigma3 branch.")
    (OUT / "mediator_r_window_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = scan()
    write_csv(summary["rows"])
    (OUT / "mediator_r_window_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    best = summary["best_by_smallest_nonuniversal_threshold"]
    fit = summary["decoupling_fit_R_ge_20"]
    print("Mediator R-window scan")
    print(f"  R grid: {R_GRID}")
    print(f"  decoupling slope R>=20: {fit['slope_log_projected_l2_vs_log_R']:.6f}")
    print(f"  best R by non-universal threshold: {best['R']:.1f}")
    print(f"  ||P Delta_med||: {best['projected_l2']:.6e}")
    print(f"  safe points: {best['safe_points']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
