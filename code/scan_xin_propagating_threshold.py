#!/usr/bin/env python3
"""Fallback scan for a propagating conormal Xi_N normal-pair threshold.

If the conormal multiplier is auxiliary/composite, it has no threshold.
If it is promoted to a propagating chiral normal-pair sector, the local Hessian
audit implies at least

    b_XiN = 2 b_54,phys = (68/5, 12, 16).

This script adds that threshold at M_XiN=kappa_XiN M_G and replays the
corrected RGE/proton scan.
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

import scan_corrected_downstream as down  # noqa: E402
import scan_kappa54_global_fixedline as k54  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "xin_propagating_threshold"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
CORRECTED_SCAN_CSV = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"

B_XIN = np.array([68.0 / 5.0, 12.0, 16.0], dtype=float)
KAPPA_XIN_GRID = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0]
R_GRID = [50.0, 200.0]


def iter_corrected_rows() -> list[dict[str, str]]:
    with CORRECTED_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_vacuum() -> dict[str, Any]:
    return json.loads(VACUUM.read_text(encoding="utf-8"))


def mediator_delta_for_r(payload: dict[str, Any], r_value: float) -> np.ndarray:
    card = payload["benchmark_cards"][f"{r_value:.1f}"]
    return np.array(card["mediator_heavy_threshold"]["delta_full"], dtype=float)


def delta_xin(kappa_xin: float) -> np.ndarray:
    return B_XIN * math.log(1.0 / kappa_xin) / (2.0 * math.pi)


def compact_best(row: dict[str, float | bool]) -> dict[str, Any]:
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
        "log_HC",
        "log_Sigma3",
        "log_Sigma8",
        "M_HC_GeV",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "residual_l2_after_mediator",
    ]
    return {key: (bool(row[key]) if isinstance(row[key], bool) else float(row[key])) for key in keys if key in row}


def scan_rows() -> list[dict[str, Any]]:
    vacuum = load_vacuum()
    source_rows = iter_corrected_rows()
    rows = []
    for r_value in R_GRID:
        med_delta = mediator_delta_for_r(vacuum, r_value)
        for kappa in KAPPA_XIN_GRID:
            dx = delta_xin(kappa)
            scan = down.scan_cached_with_delta(source_rows, med_delta + dx)
            best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
            rows.append(
                {
                    "R": r_value,
                    "kappa_XiN": kappa,
                    "M_XiN_over_MG": kappa,
                    "delta_XiN": dx.tolist(),
                    "projected_l2_XiN": float(np.linalg.norm(base.PROJECTOR @ dx)),
                    "total_projected_l2": float(np.linalg.norm(base.PROJECTOR @ (med_delta + dx))),
                    "safe_points": int(scan["safe_points"]),
                    "safe_single_scale_factor2_points": int(scan["safe_single_scale_factor2_points"]),
                    "best_is_safe": scan["best_safe"] is not None,
                    "best": compact_best(best),
                }
            )
    return rows


def allowed_window(projected_limit: float) -> dict[str, float]:
    coeff = float(np.linalg.norm(base.PROJECTOR @ (B_XIN / (2.0 * math.pi))))
    width = projected_limit / coeff
    return {
        "projected_limit": projected_limit,
        "coefficient_norm": coeff,
        "log_width": width,
        "kappa_min": math.exp(-width),
        "kappa_max": math.exp(width),
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "R",
        "kappa_XiN",
        "projected_l2_XiN",
        "total_projected_l2",
        "safe_points",
        "safe_single_scale_factor2_points",
        "best_is_safe",
        "alphaG_inv",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "score",
    ]
    with (OUT / "xin_propagating_threshold_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            best = row["best"]
            writer.writerow(
                {
                    "R": row["R"],
                    "kappa_XiN": row["kappa_XiN"],
                    "projected_l2_XiN": row["projected_l2_XiN"],
                    "total_projected_l2": row["total_projected_l2"],
                    "safe_points": row["safe_points"],
                    "safe_single_scale_factor2_points": row["safe_single_scale_factor2_points"],
                    "best_is_safe": row["best_is_safe"],
                    "alphaG_inv": best["alphaG_inv"],
                    "M_Sigma3_GeV": best["M_Sigma3_GeV"],
                    "M_Sigma8_GeV": best["M_Sigma8_GeV"],
                    "tau_dim6_years": best["tau_dim6_years"],
                    "tau_dim5_target_filter_years": best["tau_dim5_target_filter_years"],
                    "triplet_filter_required": best["triplet_filter_required"],
                    "score": best["score"],
                }
            )


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Propagating Xi_N threshold fallback scan")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Threshold")
    lines.append("")
    lines.append("`b_XiN = (68/5, 12, 16)` if the normal-bundle multiplier propagates as a chiral normal pair.")
    lines.append("")
    lines.append("## Clean-threshold windows")
    lines.append("")
    for key, window in payload["allowed_windows"].items():
        lines.append(
            f"- `{key}`: {window['kappa_min']:.6f} < kappa_XiN < {window['kappa_max']:.6f} "
            f"(coefficient {window['coefficient_norm']:.9f})"
        )
    lines.append("")
    lines.append("## Selected scan rows")
    lines.append("")
    lines.append("| R | kappa_XiN | ||P Delta_XiN|| | safe points | alphaG^-1 | M_Sigma3 [GeV] | tau_d6 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    selected = [row for row in payload["rows"] if row["kappa_XiN"] in [0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0]]
    for row in selected:
        best = row["best"]
        lines.append(
            f"| {row['R']:.0f} | {row['kappa_XiN']:.2f} | {row['projected_l2_XiN']:.6e} | "
            f"{row['safe_points']} | {best['alphaG_inv']:.6f} | "
            f"{best['M_Sigma3_GeV']:.6e} | {best['tau_dim6_years']:.6e} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "xin_propagating_threshold_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = scan_rows()
    windows = {
        "PDelta_lt_1e_minus_2": allowed_window(1.0e-2),
        "PDelta_lt_1e_minus_3": allowed_window(1.0e-3),
    }
    best_clean = min(
        [row for row in rows if row["best_is_safe"] and row["projected_l2_XiN"] < 1.0e-2],
        key=lambda row: (float(row["best"]["score"]), row["total_projected_l2"]),
    )
    best_any = min([row for row in rows if row["best_is_safe"]], key=lambda row: float(row["best"]["score"]))
    payload = {
        "note": "No web lookup used. Fallback scan for a propagating conormal normal-pair threshold.",
        "b_XiN": B_XIN.tolist(),
        "rows": rows,
        "allowed_windows": windows,
        "best_clean": best_clean,
        "best_any": best_any,
        "verdict": (
            "A propagating Xi_N normal pair is viable only if its mass is very near M_G "
            "or if its threshold is explicitly included as a new scan parameter. The "
            "auxiliary/composite Xi_N route remains cleaner because it has Delta b=0."
        ),
    }
    write_csv(rows)
    (OUT / "xin_propagating_threshold_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Propagating Xi_N threshold scan complete")
    print(f"  b_XiN={B_XIN.tolist()}")
    print(
        "  clean window 1e-2="
        f"{windows['PDelta_lt_1e_minus_2']['kappa_min']:.6f}..{windows['PDelta_lt_1e_minus_2']['kappa_max']:.6f}"
    )
    print(
        f"  best clean R={best_clean['R']:.0f} kappa={best_clean['kappa_XiN']:.3f} "
        f"safe={best_clean['safe_points']}"
    )
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
