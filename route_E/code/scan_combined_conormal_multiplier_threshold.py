#!/usr/bin/env python3
"""Fallback threshold scan for propagating combined 54/210 conormal multipliers.

The aligned conormal audit showed that the clean constrained branch needs
normal-bundle multipliers for the rank-240 constraint system.  If those
multipliers are auxiliary/composite, they add no threshold.  If they propagate
as chiral normal-pair fields, the light threshold is incomplete.

The minimal non-universal beta vector is obtained by removing the shared
SO(10)/PS orbit (6,2,2) from one 54 and one 210 and then doubling for the
normal variable plus multiplier pair:

    b_orbit(6,2,2) = (26/5, 6, 4)
    b_54_normal    = (12,12,12) - b_orbit = (34/5, 6, 8)
    b_210_normal   = (56,56,56) - b_orbit = (254/5, 50, 52)
    b_pair         = 2 (b_54_normal + b_210_normal)
                   = (576/5, 112, 120).

This script scans that threshold at M_Xi=kappa M_G and replays the corrected
RGE/proton cache.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_corrected_downstream as down  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "combined_conormal_multiplier_threshold"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
CORRECTED_SCAN_CSV = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"

B_ORBIT_622 = np.array([26.0 / 5.0, 6.0, 4.0], dtype=float)
B_54_FULL = np.array([12.0, 12.0, 12.0], dtype=float)
B_210_FULL = np.array([56.0, 56.0, 56.0], dtype=float)
B_54_NORMAL = B_54_FULL - B_ORBIT_622
B_210_NORMAL = B_210_FULL - B_ORBIT_622
B_COMBINED_NORMAL_PAIR = 2.0 * (B_54_NORMAL + B_210_NORMAL)

KAPPA_GRID = [
    0.5,
    0.8,
    0.9,
    0.95,
    0.98,
    0.99,
    0.995,
    1.0,
    1.005,
    1.01,
    1.02,
    1.05,
    1.1,
    1.25,
    2.0,
]
R_GRID = [50.0, 200.0]


def iter_corrected_rows() -> List[Dict[str, str]]:
    with CORRECTED_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_vacuum() -> Dict[str, Any]:
    return json.loads(VACUUM.read_text(encoding="utf-8"))


def mediator_delta_for_r(payload: Dict[str, Any], r_value: float) -> np.ndarray:
    card = payload["benchmark_cards"][f"{r_value:.1f}"]
    return np.array(card["mediator_heavy_threshold"]["delta_full"], dtype=float)


def threshold_delta(beta: np.ndarray, kappa: float) -> np.ndarray:
    return beta * math.log(1.0 / kappa) / (2.0 * math.pi)


def projected_l2(vec: np.ndarray) -> float:
    return float(np.linalg.norm(base.PROJECTOR @ vec))


def allowed_window(beta: np.ndarray, projected_limit: float) -> Dict[str, float]:
    coeff = projected_l2(beta / (2.0 * math.pi))
    width = projected_limit / coeff
    return {
        "projected_limit": projected_limit,
        "coefficient_norm": coeff,
        "max_abs_log_kappa": width,
        "kappa_min": math.exp(-width),
        "kappa_max": math.exp(width),
    }


def compact_best(row: Dict[str, float | bool]) -> Dict[str, Any]:
    keys = [
        "score",
        "tan_beta",
        "MSUSY_GeV",
        "MG_GeV",
        "alphaG_inv",
        "M_HC_GeV",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "residual_l2_after_mediator",
    ]
    return {key: (bool(row[key]) if isinstance(row[key], bool) else float(row[key])) for key in keys if key in row}


def scan_rows() -> List[Dict[str, Any]]:
    vacuum = load_vacuum()
    source_rows = iter_corrected_rows()
    rows: List[Dict[str, Any]] = []
    for r_value in R_GRID:
        med = mediator_delta_for_r(vacuum, r_value)
        for kappa in KAPPA_GRID:
            dxi = threshold_delta(B_COMBINED_NORMAL_PAIR, kappa)
            scan = down.scan_cached_with_delta(source_rows, med + dxi)
            best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
            rows.append(
                {
                    "R": r_value,
                    "kappa_Xi_combined": kappa,
                    "delta_Xi": dxi.tolist(),
                    "projected_l2_Xi": projected_l2(dxi),
                    "total_projected_l2": projected_l2(med + dxi),
                    "safe_points": int(scan["safe_points"]),
                    "safe_single_scale_factor2_points": int(scan["safe_single_scale_factor2_points"]),
                    "best_is_safe": scan["best_safe"] is not None,
                    "best": compact_best(best),
                }
            )
    return rows


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    fields = [
        "R",
        "kappa_Xi_combined",
        "projected_l2_Xi",
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            best = row["best"]
            writer.writerow(
                {
                    "R": row["R"],
                    "kappa_Xi_combined": row["kappa_Xi_combined"],
                    "projected_l2_Xi": row["projected_l2_Xi"],
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


def make_report(payload: Dict[str, Any]) -> str:
    lines = [
        "# Combined conormal multiplier threshold fallback",
        "",
        "No web lookup was used.",
        "",
        "## Beta vector",
        "",
        f"`b_orbit_622 = {payload['b_orbit_622']}`",
        f"`b_54_normal = {payload['b_54_normal']}`",
        f"`b_210_normal = {payload['b_210_normal']}`",
        f"`b_combined_normal_pair = {payload['b_combined_normal_pair']}`",
        "",
        "## Clean windows",
        "",
    ]
    for key, window in payload["allowed_windows"].items():
        lines.append(
            f"- `{key}`: {window['kappa_min']:.6f} < kappa_Xi < {window['kappa_max']:.6f} "
            f"(coefficient {window['coefficient_norm']:.9f})"
        )
    lines.extend(
        [
            "",
            "## Selected rows",
            "",
            "| R | kappa | ||P Delta_Xi|| | safe points | alphaG^-1 | M_Sigma3 [GeV] | tau_d6 [yr] |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    selected_kappas = {0.9, 0.99, 0.995, 1.0, 1.005, 1.01, 1.1}
    for row in payload["rows"]:
        if row["kappa_Xi_combined"] not in selected_kappas:
            continue
        best = row["best"]
        lines.append(
            f"| {row['R']:.0f} | {row['kappa_Xi_combined']:.3f} | "
            f"{row['projected_l2_Xi']:.6e} | {row['safe_points']} | "
            f"{best['alphaG_inv']:.6f} | {best['M_Sigma3_GeV']:.6e} | "
            f"{best['tau_dim6_years']:.6e} |"
        )
    lines.extend(["", "## Verdict", "", payload["verdict"], ""])
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = scan_rows()
    windows = {
        "PDelta_lt_1e_minus_2": allowed_window(B_COMBINED_NORMAL_PAIR, 1.0e-2),
        "PDelta_lt_1e_minus_3": allowed_window(B_COMBINED_NORMAL_PAIR, 1.0e-3),
    }
    clean = [
        row
        for row in rows
        if row["best_is_safe"] and row["projected_l2_Xi"] < 1.0e-2
    ]
    best_clean = min(clean, key=lambda row: (float(row["best"]["score"]), row["total_projected_l2"]))
    payload: Dict[str, Any] = {
        "note": "No web lookup used. Propagating combined 54/210 conormal normal-pair threshold fallback.",
        "b_orbit_622": B_ORBIT_622.tolist(),
        "b_54_normal": B_54_NORMAL.tolist(),
        "b_210_normal": B_210_NORMAL.tolist(),
        "b_combined_normal_pair": B_COMBINED_NORMAL_PAIR.tolist(),
        "allowed_windows": windows,
        "rows": rows,
        "best_clean": best_clean,
        "verdict": (
            "The combined normal-bundle multiplier may propagate only if its mass "
            "is locked extremely close to M_G: roughly a one-percent window for "
            "||P Delta||<1e-2 and a per-mille window for ||P Delta||<1e-3.  This "
            "is tighter than the 54-only Xi_N fallback.  Therefore the preferred "
            "paper branch should state that the conormal multipliers are auxiliary "
            "or composite; a propagating normal-pair sector is only a tuned EFT "
            "fallback whose threshold must be scanned explicitly."
        ),
    }
    write_csv(OUT / "combined_conormal_multiplier_threshold_scan.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "report.md").write_text(make_report(payload), encoding="utf-8")
    print("Combined conormal multiplier threshold scan complete")
    print(f"b_combined={B_COMBINED_NORMAL_PAIR.tolist()}")
    print(
        "window_1e-2="
        f"{windows['PDelta_lt_1e_minus_2']['kappa_min']:.6f}.."
        f"{windows['PDelta_lt_1e_minus_2']['kappa_max']:.6f}"
    )
    print(
        "window_1e-3="
        f"{windows['PDelta_lt_1e_minus_3']['kappa_min']:.6f}.."
        f"{windows['PDelta_lt_1e_minus_3']['kappa_max']:.6f}"
    )
    print(
        f"best_clean_R={best_clean['R']:.0f} "
        f"kappa={best_clean['kappa_Xi_combined']:.3f} "
        f"safe_points={best_clean['safe_points']}"
    )
    print(f"outputs={OUT}")


if __name__ == "__main__":
    main()
