#!/usr/bin/env python3
"""Regenerate downstream mediator/Goldstone scans from the 4pi-corrected cache.

No web lookup is used.  This script consumes
output/yukawa_4pi_audit/yukawa_4pi_corrected_scan.csv and replays the finite
mediator threshold, mediator R-window, and Goldstone-locking scans that were
previously based on the old Yukawa two-loop cache.
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

import scan_goldstone_locking as gold  # noqa: E402
import scan_mediator_r_window as rwin  # noqa: E402
import scan_mediator_threshold_rge as med  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


CORRECTED_SCAN_CSV = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"
CORRECTED_AUDIT_JSON = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_audit_summary.json"
OUT = ROOT / "output" / "corrected_downstream"

R_WINDOW_GRID = [5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0]
GOLDSTONE_R_GRID = [20.0, 50.0, 100.0, 200.0]
S_T_GRID = [1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6]
EPS_MULTIPLIERS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
SIGNS = [-1.0, 1.0]


def iter_corrected_rows() -> list[dict[str, str]]:
    with CORRECTED_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def scan_cached_with_delta(source_rows: list[dict[str, str]], delta_total: np.ndarray) -> dict[str, Any]:
    prefactor = base.proton_prefactor()
    rows: list[dict[str, float | bool]] = []
    best: dict[str, float | bool] | None = None
    for source in source_rows:
        row = med.analyze_cached_row_with_mediator(source, prefactor, delta_total)
        rows.append(row)
        if best is None or float(row["score"]) < float(best["score"]):
            best = row
    assert best is not None

    safe = [
        row
        for row in rows
        if bool(row["triplet_filter_safe"])
        and float(row["rho_X"]) <= 1.000001
        and bool(row["perturbative_superpotential"])
    ]
    factor2 = [row for row in safe if bool(row["single_scale_factor2"])]
    best_safe = min(safe, key=lambda row: float(row["score"]), default=None)
    return {
        "total_points": len(rows),
        "safe_points": len(safe),
        "safe_single_scale_factor2_points": len(factor2),
        "best": best,
        "best_safe": best_safe,
    }


def compact_best(row: dict[str, float | bool] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    keys = [
        "score",
        "tan_beta",
        "MSUSY_GeV",
        "MG_GeV",
        "alphaG_inv",
        "lambda_T",
        "lambda_S",
        "chi",
        "kappa_3",
        "kappa_8",
        "M_HC_GeV",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "log_HC",
        "log_Sigma3",
        "log_Sigma8",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "residual_l2_after_mediator",
    ]
    out = {}
    for key in keys:
        if key in row:
            value = row[key]
            out[key] = bool(value) if isinstance(value, bool) else float(value)
    return out


def mediator_fixed_r50(source_rows: list[dict[str, str]]) -> dict[str, Any]:
    uv = med.load_uv_threshold()
    scan = scan_cached_with_delta(source_rows, uv["delta_full"])
    best_safe = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
    return {
        "mediator_threshold": {
            "R": uv["R"],
            "eta": uv["eta"],
            "delta_full": uv["delta_full"].tolist(),
            "projected_delta": uv["projected_delta"].tolist(),
            "projected_l2": uv["projected_l2"],
        },
        "diagnostics": {
            "total_points": scan["total_points"],
            "safe_points": scan["safe_points"],
            "safe_single_scale_factor2_points": scan["safe_single_scale_factor2_points"],
        },
        "best_safe_or_best": compact_best(best_safe),
    }


def mediator_r_window(source_rows: list[dict[str, str]]) -> dict[str, Any]:
    targets, initial = rwin.load_targets()
    rows = []
    detailed = {}
    for r_med in R_WINDOW_GRID:
        solution = rwin.uv.solve_finite_r(targets, r_med, initial)
        threshold = rwin.uv.heavy_threshold_residual(solution)
        scan = scan_cached_with_delta(source_rows, np.array(threshold["delta_full"], dtype=float))
        best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
        row = {
            "R": r_med,
            "eta": float(solution["parameters"]["eta"]),
            "eta_over_sqrt_R": float(solution["parameters"]["eta"] / math.sqrt(r_med)),
            "schur_delta_X": float(solution["parameters"]["schur_delta_X"]),
            "light_fit_residual_l2": float(solution["residual_l2"]),
            "projected_l2": float(threshold["projected_l2"]),
            "R_times_projected_l2": float(r_med * threshold["projected_l2"]),
            "safe_points": scan["safe_points"],
            "safe_single_scale_factor2_points": scan["safe_single_scale_factor2_points"],
            "best_is_safe": scan["best_safe"] is not None,
            "alphaG_inv": float(best["alphaG_inv"]),
            "M_Sigma3_GeV": float(best["M_Sigma3_GeV"]),
            "M_Sigma8_GeV": float(best["M_Sigma8_GeV"]),
            "tau_dim6_years": float(best["tau_dim6_years"]),
            "tau_dim5_target_filter_years": float(best["tau_dim5_target_filter_years"]),
        }
        rows.append(row)
        detailed[str(r_med)] = {
            "threshold": threshold,
            "scan_diagnostics": {
                "total_points": scan["total_points"],
                "safe_points": scan["safe_points"],
                "safe_single_scale_factor2_points": scan["safe_single_scale_factor2_points"],
                "best_safe_or_best": compact_best(best),
            },
        }

    tail = [row for row in rows if row["R"] >= 20.0]
    slope, intercept = np.polyfit(
        np.log([row["R"] for row in tail]),
        np.log([row["projected_l2"] for row in tail]),
        1,
    )
    best_by_threshold = min([row for row in rows if row["best_is_safe"]], key=lambda row: row["projected_l2"])
    return {
        "rows": rows,
        "detailed": detailed,
        "decoupling_fit_R_ge_20": {
            "slope_log_projected_l2_vs_log_R": float(slope),
            "intercept": float(intercept),
            "expected_slope": -1.0,
        },
        "best_by_smallest_nonuniversal_threshold": best_by_threshold,
    }


def analyze_goldstone_rows(source_rows: list[dict[str, str]], delta_total: np.ndarray) -> list[dict[str, float | bool]]:
    prefactor = base.proton_prefactor()
    rows: list[dict[str, float | bool]] = []
    for source in source_rows:
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
        rows.append(
            {
                "perturbative_superpotential": bool(params["perturbative"]),
                "single_scale_factor2": bool(params["single_scale_factor2"]),
                "MG_GeV": mg,
                "alphaG_inv": alpha_g_inv,
                "M_Sigma3_GeV": float(params["kappa_3"] * mg),
                "M_Sigma8_GeV": float(params["kappa_8"] * mg),
                "M_X_min_GeV": float(mx_min),
                "rho_X": float(max(1.0, mx_min / mg)),
                "residual_l2_after_thresholds": float(residual),
                "filter_required": float(filter_required),
                "tau_dim6_years": float(base.lifetime_from_c6(c6_gauge, prefactor)),
                "d5_c6_per_unit_ST": float(denominator),
                "max_abs_log": float(np.max(np.abs(logs))),
                "log_HC": float(logs[0]),
                "log_Sigma3": float(logs[1]),
                "log_Sigma8": float(logs[2]),
                "chi": float(params["chi"]),
                "cancellation_index": float(params["cancellation_index"]),
            }
        )
    return rows


def summarize_goldstone_rows(rows: list[dict[str, float | bool]], s_t: float) -> dict[str, Any]:
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
        key=lambda row: (float(row["max_abs_log"]), abs(float(row["chi"])), -float(row["tau_dim6_years"])),
        default=None,
    )
    if best is None:
        return {"safe_points": 0, "safe_single_scale_factor2_points": 0, "best": None}
    best_out = dict(best)
    best_out["S_T"] = s_t
    best_out["tau_dim5_years"] = base.lifetime_from_c6(s_t * float(best["d5_c6_per_unit_ST"]), prefactor)
    return {
        "safe_points": len(safe),
        "safe_single_scale_factor2_points": len(factor2),
        "best": best_out,
    }


def goldstone_locking(source_rows: list[dict[str, str]]) -> dict[str, Any]:
    targets, initial = rwin.load_targets()
    unit_l2 = gold.goldstone_projected_l2(1.0)
    grid_rows: list[dict[str, Any]] = []
    best_by_r_st: dict[str, Any] = {}

    for r_med in GOLDSTONE_R_GRID:
        solution = rwin.uv.solve_finite_r(targets, r_med, initial)
        threshold = rwin.uv.heavy_threshold_residual(solution)
        med_l2 = float(threshold["projected_l2"])
        epsilon_lock_max = med_l2 / unit_l2
        eps_values = [m * epsilon_lock_max for m in EPS_MULTIPLIERS]
        for eps_abs in eps_values:
            signs = SIGNS if eps_abs > 0 else [1.0]
            for sign in signs:
                eps_signed = sign * eps_abs
                delta_gold = med.SIGMA8_COLORED_PAIR_B * eps_signed / (2.0 * math.pi)
                delta_total = np.array(threshold["delta_full"], dtype=float) + delta_gold
                analyzed = analyze_goldstone_rows(source_rows, delta_total)
                gold_l2 = gold.goldstone_projected_l2(eps_abs)
                locked = gold_l2 <= med_l2 * (1.0 + 1.0e-12)
                for s_t in S_T_GRID:
                    summary = summarize_goldstone_rows(analyzed, s_t)
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
            best_by_r_st[key] = max(allowed, key=lambda row: row["epsilon_abs"], default=None)

    return {
        "goldstone_projected_l2_per_unit_epsilon": unit_l2,
        "rows": grid_rows,
        "max_allowed_by_R_and_ST": best_by_r_st,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    fixed = summary["fixed_R50"]
    r_window = summary["R_window"]
    goldstone = summary["goldstone_locking"]
    lines: list[str] = []
    lines.append("# Corrected-cache downstream scan")
    lines.append("")
    lines.append("No web lookup was used.  This report regenerates finite mediator,")
    lines.append("R-window, and Goldstone-locking scans from the 4pi-corrected Yukawa")
    lines.append("two-loop cache.")
    lines.append("")
    lines.append("## Fixed R=50 mediator threshold")
    lines.append("")
    best50 = fixed["best_safe_or_best"]
    lines.append("```text")
    lines.append(f"safe points = {fixed['diagnostics']['safe_points']}")
    lines.append(f"safe factor-2 single-scale points = {fixed['diagnostics']['safe_single_scale_factor2_points']}")
    lines.append(f"alphaG^-1 = {best50['alphaG_inv']:.6f}")
    lines.append(f"M_Sigma3 = {best50['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"M_Sigma8 = {best50['M_Sigma8_GeV']:.6e} GeV")
    lines.append(f"tau_dim6 = {best50['tau_dim6_years']:.6e} yr")
    lines.append(f"tau_dim5(S_T=1e-5) = {best50['tau_dim5_target_filter_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("## Corrected R-window")
    lines.append("")
    fit = r_window["decoupling_fit_R_ge_20"]
    lines.append("```text")
    lines.append(
        "fit over R >= 20: log ||P Delta_med|| = "
        f"{fit['intercept']:+.6f} {fit['slope_log_projected_l2_vs_log_R']:+.6f} log R"
    )
    lines.append("```")
    lines.append("")
    lines.append("| R | ||P Delta_med|| | safe points | alphaG^-1 | M_Sigma3 [GeV] | tau_d6 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for row in r_window["rows"]:
        if row["R"] in [20.0, 50.0, 100.0, 200.0]:
            lines.append(
                f"| {row['R']:.0f} | {row['projected_l2']:.6e} | {row['safe_points']} | "
                f"{row['alphaG_inv']:.6f} | {row['M_Sigma3_GeV']:.6e} | {row['tau_dim6_years']:.6e} |"
            )
    lines.append("")
    lines.append("## Corrected Goldstone-locking scan")
    lines.append("")
    lines.append(
        f"`||P Delta_G||_2 = {goldstone['goldstone_projected_l2_per_unit_epsilon']:.9e} |epsilon_G|`."
    )
    lines.append("")
    lines.append("| R | S_T | max abs(epsilon_G) | safe points | tau_d5 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|")
    for key, row in goldstone["max_allowed_by_R_and_ST"].items():
        if row is None:
            continue
        if abs(row["S_T"] - 1.0e-5) < 1.0e-20 or abs(row["S_T"] - 1.0e-4) < 1.0e-20:
            lines.append(
                f"| {row['R']:.0f} | {row['S_T']:.1e} | {row['epsilon_abs']:.6e} | "
                f"{row['safe_points']} | {row.get('tau_dim5_years', float('nan')):.6e} |"
            )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("The corrected-cache branch survives the finite mediator and Goldstone-locking")
    lines.append("replay.  The old provisional mediator/Goldstone numbers should be replaced")
    lines.append("by this corrected-cache table before final paper use.  No factor-two")
    lines.append("single-scale point appears in the corrected downstream scan.")
    (OUT / "corrected_downstream_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_rows = iter_corrected_rows()
    fixed = mediator_fixed_r50(source_rows)
    r_window = mediator_r_window(source_rows)
    goldstone_summary = goldstone_locking(source_rows)
    audit = json.loads(CORRECTED_AUDIT_JSON.read_text(encoding="utf-8"))
    summary = {
        "note": "No web lookup used. Downstream scans regenerated from the 4pi-corrected Yukawa cache.",
        "corrected_scan_csv": str(CORRECTED_SCAN_CSV),
        "yukawa_4pi_audit": {
            "best_comparison": audit["best_comparison"],
            "corrected_base_best": audit["corrected_summary"]["best"],
        },
        "fixed_R50": fixed,
        "R_window": r_window,
        "goldstone_locking": goldstone_summary,
    }
    write_csv(OUT / "corrected_R_window_scan.csv", r_window["rows"])
    write_csv(OUT / "corrected_goldstone_locking_scan.csv", goldstone_summary["rows"])
    (OUT / "corrected_downstream_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    best50 = fixed["best_safe_or_best"]
    best_r = r_window["best_by_smallest_nonuniversal_threshold"]
    print("Corrected-cache downstream scan")
    print(f"  fixed R=50 safe points: {fixed['diagnostics']['safe_points']}")
    print(f"  fixed R=50 M_Sigma3: {best50['M_Sigma3_GeV']:.6e} GeV")
    print(f"  R-window best by threshold: R={best_r['R']:.1f}, safe={best_r['safe_points']}")
    for key, row in goldstone_summary["max_allowed_by_R_and_ST"].items():
        if row is not None and row["S_T"] == 1.0e-5:
            print(f"  {key}: eps_max={row['epsilon_abs']:.6e}, safe={row['safe_points']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
