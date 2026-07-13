#!/usr/bin/env python3
"""Goldstone-locking and dimension-5 filter scan.

No web lookup is used.  This script tests how tightly the colored Goldstone
pair in (15,1,1) must remain locked to MG once finite mediator thresholds and
dimension-5 proton-decay filters are scanned together.
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

import scan_mediator_r_window as rwin  # noqa: E402
import scan_mediator_threshold_rge as med  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "goldstone_locking"
R_GRID = [20.0, 50.0, 100.0, 200.0]
S_T_GRID = [1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6]
EPS_MULTIPLIERS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
SIGNS = [-1.0, 1.0]


def goldstone_projected_l2(epsilon_g: float) -> float:
    unit = base.PROJECTOR @ med.SIGMA8_COLORED_PAIR_B / (2.0 * math.pi)
    return float(abs(epsilon_g) * np.linalg.norm(unit))


def mediator_threshold_by_r() -> dict[float, dict[str, Any]]:
    targets, initial = rwin.load_targets()
    out: dict[float, dict[str, Any]] = {}
    for r_med in R_GRID:
        solution = rwin.uv.solve_finite_r(targets, r_med, initial)
        threshold = rwin.uv.heavy_threshold_residual(solution)
        out[r_med] = {
            "solution": solution,
            "threshold": threshold,
            "projected_l2": float(threshold["projected_l2"]),
        }
    return out


def analyze_cached_rows(delta_total: np.ndarray) -> list[dict[str, float | bool]]:
    prefactor = base.proton_prefactor()
    rows: list[dict[str, float | bool]] = []
    with med.BASE_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for source in reader:
            mg = float(source["MG_GeV"])
            msusy = float(source["MSUSY_GeV"])
            alpha_inv = np.array(
                [
                    float(source["alpha1_inv_MG"]),
                    float(source["alpha2_inv_MG"]),
                    float(source["alpha3_inv_MG"]),
                ],
                dtype=float,
            )
            yt_mg = float(source["yt_MG"])
            yb_mg = float(source["yb_MG"])
            logs, _threshold, residual, alpha_g_inv = med.exact_threshold_logs_with_mediator(alpha_inv, delta_total)
            params = base.superpotential_from_logs(logs)
            mx_min = base.mx_required(alpha_g_inv, prefactor)
            m_t = params["lambda_T"] * mg
            loop = (1.0 / 25.0) / (4.0 * math.pi)
            c6_limit = base.c6_max_for_target(prefactor)
            denominator = (yt_mg * yb_mg / m_t) * loop * base.MWINO_GEV / (msusy * msusy)
            filter_required = c6_limit / max(denominator, 1.0e-300)
            c6_gauge = (4.0 * math.pi / alpha_g_inv) / (mg * mg)
            d5_c6_per_unit_st = denominator
            rows.append(
                {
                    "tan_beta": float(source["tan_beta"]),
                    "MSUSY_GeV": msusy,
                    "MG_GeV": mg,
                    "alphaG_inv": alpha_g_inv,
                    "lambda_T": float(params["lambda_T"]),
                    "lambda_S": float(params["lambda_S"]),
                    "chi": float(params["chi"]),
                    "kappa_3": float(params["kappa_3"]),
                    "kappa_8": float(params["kappa_8"]),
                    "cancellation_index": float(params["cancellation_index"]),
                    "perturbative_superpotential": bool(params["perturbative"]),
                    "single_scale_factor2": bool(params["single_scale_factor2"]),
                    "log_HC": float(logs[0]),
                    "log_Sigma3": float(logs[1]),
                    "log_Sigma8": float(logs[2]),
                    "M_HC_GeV": float(params["lambda_T"] * mg),
                    "M_Sigma3_GeV": float(params["kappa_3"] * mg),
                    "M_Sigma8_GeV": float(params["kappa_8"] * mg),
                    "M_X_min_GeV": float(mx_min),
                    "rho_X": float(max(1.0, mx_min / mg)),
                    "residual_l2_after_thresholds": float(residual),
                    "filter_required": float(filter_required),
                    "tau_dim6_years": float(base.lifetime_from_c6(c6_gauge, prefactor)),
                    "d5_c6_per_unit_ST": float(d5_c6_per_unit_st),
                    "max_abs_log": float(np.max(np.abs(logs))),
                }
            )
    return rows


def summarize_rows(rows: list[dict[str, float | bool]], s_t: float) -> dict[str, Any]:
    prefactor = base.proton_prefactor()
    safe = [
        row
        for row in rows
        if bool(row["perturbative_superpotential"])
        and float(row["rho_X"]) <= 1.000001
        and float(row["filter_required"]) >= s_t
    ]
    factor2 = [row for row in safe if bool(row["single_scale_factor2"])]
    best = min(
        safe,
        key=lambda row: (
            float(row["max_abs_log"]),
            abs(float(row["chi"])),
            -float(row["tau_dim6_years"]),
        ),
        default=None,
    )
    if best is None:
        return {"safe_points": 0, "safe_single_scale_factor2_points": 0, "best": None}
    tau_d5 = base.lifetime_from_c6(s_t * float(best["d5_c6_per_unit_ST"]), prefactor)
    best_out = dict(best)
    best_out["S_T"] = s_t
    best_out["tau_dim5_years"] = tau_d5
    best_out["triplet_filter_safe"] = bool(s_t <= float(best["filter_required"]))
    return {
        "safe_points": len(safe),
        "safe_single_scale_factor2_points": len(factor2),
        "best": best_out,
    }


def scan() -> dict[str, Any]:
    med_by_r = mediator_threshold_by_r()
    unit_l2 = goldstone_projected_l2(1.0)
    grid_rows: list[dict[str, Any]] = []
    best_by_r_st: dict[str, Any] = {}

    for r_med in R_GRID:
        med_data = med_by_r[r_med]
        med_l2 = med_data["projected_l2"]
        epsilon_lock_max = med_l2 / unit_l2
        eps_values = [m * epsilon_lock_max for m in EPS_MULTIPLIERS]
        for eps_abs in eps_values:
            for sign in SIGNS if eps_abs > 0 else [1.0]:
                eps_signed = sign * eps_abs
                delta_gold = med.SIGMA8_COLORED_PAIR_B * eps_signed / (2.0 * math.pi)
                delta_total = np.array(med_data["threshold"]["delta_full"], dtype=float) + delta_gold
                analyzed = analyze_cached_rows(delta_total)
                gold_l2 = goldstone_projected_l2(eps_abs)
                locked = gold_l2 <= med_l2 * (1.0 + 1.0e-12)
                for s_t in S_T_GRID:
                    summary = summarize_rows(analyzed, s_t)
                    best = summary["best"]
                    row = {
                        "R": r_med,
                        "S_T": s_t,
                        "epsilon_abs": eps_abs,
                        "epsilon_signed": eps_signed,
                        "epsilon_over_lock_max": eps_abs / epsilon_lock_max if epsilon_lock_max > 0 else math.inf,
                        "goldstone_projected_l2": gold_l2,
                        "mediator_projected_l2": med_l2,
                        "goldstone_locked": locked,
                        "safe_points": summary["safe_points"],
                        "safe_single_scale_factor2_points": summary["safe_single_scale_factor2_points"],
                        "has_allowed_point": bool(locked and summary["safe_points"] > 0),
                    }
                    if best is not None:
                        for key in [
                            "MG_GeV",
                            "alphaG_inv",
                            "M_Sigma3_GeV",
                            "M_Sigma8_GeV",
                            "tau_dim6_years",
                            "tau_dim5_years",
                            "filter_required",
                            "log_HC",
                            "log_Sigma3",
                            "log_Sigma8",
                            "chi",
                            "cancellation_index",
                        ]:
                            row[key] = best[key]
                    grid_rows.append(row)

        for s_t in S_T_GRID:
            allowed = [
                row
                for row in grid_rows
                if row["R"] == r_med and row["S_T"] == s_t and row["has_allowed_point"]
            ]
            key = f"R={r_med:g},S_T={s_t:.1e}"
            if allowed:
                best_allowed = max(allowed, key=lambda row: row["epsilon_abs"])
                best_by_r_st[key] = best_allowed
            else:
                best_by_r_st[key] = None

    return {
        "note": "No web lookup used.  Goldstone-locking scan adds the colored-pair threshold and varies S_T.",
        "R_grid": R_GRID,
        "S_T_grid": S_T_GRID,
        "epsilon_multipliers": EPS_MULTIPLIERS,
        "goldstone_pair_b_vector": med.SIGMA8_COLORED_PAIR_B.tolist(),
        "goldstone_projected_l2_per_unit_epsilon": unit_l2,
        "rows": grid_rows,
        "max_allowed_by_R_and_ST": best_by_r_st,
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with (OUT / "goldstone_locking_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Goldstone-locking and dimension-5 scan")
    lines.append("")
    lines.append("No web lookup was used.  The scan adds a possible colored-Goldstone")
    lines.append("threshold")
    lines.append("")
    lines.append("```text")
    lines.append("Delta_G = (8/5, 0, 1) log(MG/M_Goldstone)/(2 pi)")
    lines.append("```")
    lines.append("")
    lines.append("before solving the H_C/Sigma3/Sigma8 matching logs.  The locking criterion")
    lines.append("is that this accidental colored-pair threshold is no larger than the")
    lines.append("finite mediator non-universal threshold:")
    lines.append("")
    lines.append("```text")
    lines.append("||P Delta_G||_2 <= ||P Delta_med(R)||_2")
    lines.append("```")
    lines.append("")
    lines.append(
        f"Unit conversion: `||P Delta_G||_2 = {summary['goldstone_projected_l2_per_unit_epsilon']:.9e} |epsilon_G|`."
    )
    lines.append("")
    lines.append("## Maximum allowed Goldstone displacement")
    lines.append("")
    lines.append("| R | S_T | max abs(epsilon_G) | safe points at max | tau_d6 [yr] | tau_d5 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for key, row in summary["max_allowed_by_R_and_ST"].items():
        if row is None:
            r_part, st_part = key.split(",")
            lines.append(f"| {r_part.split('=')[1]} | {st_part.split('=')[1]} | none | 0 | -- | -- |")
            continue
        lines.append(
            f"| {row['R']:.0f} | {row['S_T']:.1e} | {row['epsilon_abs']:.6e} | "
            f"{row['safe_points']} | {row.get('tau_dim6_years', float('nan')):.6e} | "
            f"{row.get('tau_dim5_years', float('nan')):.6e} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("The allowed `epsilon_G` is set by the threshold-locking condition and")
    lines.append("shrinks with the mediator residual.  Stronger decoupling in `R` therefore")
    lines.append("forces the colored Goldstone pair closer to `MG`.  For the displayed")
    lines.append("S_T <= 1e-5 branch, dimension-5 proton decay remains safe throughout the")
    lines.append("locked region; for larger `S_T`, the allowed region can disappear because")
    lines.append("the triplet filter is too weak.")
    (OUT / "goldstone_locking_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = scan()
    write_csv(summary["rows"])
    (OUT / "goldstone_locking_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    print("Goldstone-locking scan")
    print(f"  unit ||P Delta_G|| per epsilon: {summary['goldstone_projected_l2_per_unit_epsilon']:.9e}")
    for key, row in summary["max_allowed_by_R_and_ST"].items():
        if row is not None and row["S_T"] == 1.0e-5:
            print(
                f"  {key}: eps_max={row['epsilon_abs']:.6e}, safe={row['safe_points']}, "
                f"tau_d5={row.get('tau_dim5_years', float('nan')):.3e}"
            )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
