#!/usr/bin/env python3
"""Audit representation and symmetry status of the unitary triplet link.

The threshold-locked dilation proves a useful algebraic fact: a target triplet
inverse propagator W can be embedded as a visible subblock of the inverse of an
exactly degenerate 8x8 heavy mass matrix.  This script asks the next question:
can that 8x8 unitary mass block be interpreted as an ordinary perturbative
Spin(10) superpotential sector?

The answer is intentionally conservative:

* a single species of real Spin(10) 10s needs a symmetric copy-space mass
  matrix, which the audited unitary blocks do not satisfy;
* two species of complete 10s allow an arbitrary 8x8 mass matrix, but the
  extra Dynkin index loses the R=200 perturbative reach;
* a post-Spin(10) SU(5)-complete or constrained/composite unitary-link branch
  can keep the projected threshold silent, but remains conditional.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "unitary_link_spin10_embedding"

LOCKED = ROOT / "output" / "threshold_locked_triplet_lift" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"
XI = ROOT / "output" / "combined_conormal_multiplier_threshold" / "summary.json"
NLSM_SCORE = ROOT / "output" / "nlsm_uv_completion_scorecard" / "summary.json"

GAUGE_CONTRIBUTION = 24.0
T_10 = 1.0
B_TRIPLET_PAIR = np.array([2.0 / 5.0, 0.0, 1.0], dtype=float)
B_FIVE_PAIR = np.array([1.0, 1.0, 1.0], dtype=float)
PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def projected_l2(vec: np.ndarray) -> float:
    return float(np.linalg.norm(PROJECTOR @ vec))


def alpha_inv_at_ratio(alpha_inv: float, b10: float, ratio: float) -> float:
    return alpha_inv - b10 * math.log(ratio) / (2.0 * math.pi)


def landau_ratio(alpha_inv: float, b10: float) -> float:
    if b10 <= 0.0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def load_benchmark() -> dict[str, float]:
    ledger = read_json(LEDGER)
    xi = read_json(XI)
    nlsm = read_json(NLSM_SCORE)
    baseline = next(row for row in nlsm["route_summaries"] if row["route"] == "D_cutoff_shared_U_NLSM")
    return {
        "alphaG_inv_MG": float(ledger["r200_benchmark"]["alphaG_inv"]),
        "R_target": float(ledger["r200_benchmark"]["R"]),
        "MG_GeV": float(xi["best_clean"]["best"]["MG_GeV"]),
        "M_lock_GeV": float(read_json(LOCKED)["verdict"]["M_lock_GeV"]),
        "baseline_sum_T": float(baseline["sum_T"]),
        "baseline_b10": float(baseline["b10"]),
        "baseline_alpha_inv_at_R200": float(baseline["R200_alpha_inv_at_RMG"]),
        "r200_threshold_tol": float(ledger["r200_benchmark"]["total_projected_l2"]),
    }


def matrix_obstruction_rows() -> list[dict[str, Any]]:
    locked = read_json(LOCKED)
    rows = []
    for row in locked["rows"]:
        u = cmat(row["unitary_dilation_matrix"])
        n = u.shape[0]
        eye = np.eye(n, dtype=complex)
        rows.append(
            {
                "label": row["label"],
                "dimension": n,
                "unitary_residual_2norm": float(np.linalg.norm(u.conj().T @ u - eye, ord=2)),
                "symmetric_mass_residual_fro": float(np.linalg.norm(u - u.T, ord="fro")),
                "symmetric_mass_residual_relative": float(np.linalg.norm(u - u.T, ord="fro") / np.linalg.norm(u, ord="fro")),
                "complex_orthogonal_residual_2norm": float(np.linalg.norm(u.T @ u - eye, ord=2)),
                "holomorphic_F_orthogonality_possible": bool(np.linalg.norm(u.T @ u - eye, ord=2) < 1.0e-10),
                "single_species_10_mass_possible": bool(np.linalg.norm(u - u.T, ord="fro") < 1.0e-10),
                "worst_future_margin_1e35": float(row["worst_future_margin_1e35"]),
            }
        )
    return rows


def scenario_rows(bench: dict[str, float]) -> list[dict[str, Any]]:
    log_mg_over_mlock = math.log(bench["MG_GeV"] / bench["M_lock_GeV"])
    triplet_only_delta = 8.0 * B_TRIPLET_PAIR * log_mg_over_mlock / (2.0 * math.pi)
    su5_complete_delta = 8.0 * B_FIVE_PAIR * log_mg_over_mlock / (2.0 * math.pi)
    scenarios = [
        {
            "scenario": "single_species_8_complete_10s",
            "spin10_status": "renormalizable but copy-space mass must be symmetric",
            "extra_T": 8.0 * T_10,
            "projected_threshold_l2": 0.0,
            "universal_threshold_shift": float(np.mean(su5_complete_delta)),
            "matrix_condition": "requires U=U^T; fails for audited targets",
            "microscopic_level": "elementary Spin(10), but matrix form fails",
        },
        {
            "scenario": "two_species_16_complete_10s",
            "spin10_status": "renormalizable arbitrary 10_L M 10_R mass matrix",
            "extra_T": 16.0 * T_10,
            "projected_threshold_l2": 0.0,
            "universal_threshold_shift": float(np.mean(su5_complete_delta)),
            "matrix_condition": "arbitrary unitary block allowed by species charge",
            "microscopic_level": "elementary Spin(10), but UV reach fails",
        },
        {
            "scenario": "post_spin10_8_su5_complete_5pairs",
            "spin10_status": "post-breaking EFT, not a standalone Spin(10) invariant",
            "extra_T": 0.0,
            "projected_threshold_l2": projected_l2(su5_complete_delta),
            "universal_threshold_shift": float(np.mean(su5_complete_delta)),
            "matrix_condition": "arbitrary unitary block allowed after breaking",
            "microscopic_level": "conditional EFT",
        },
        {
            "scenario": "constrained_composite_unitary_link",
            "spin10_status": "unitary block treated as constrained/composite source",
            "extra_T": 0.0,
            "projected_threshold_l2": 0.0,
            "universal_threshold_shift": 0.0,
            "matrix_condition": "unitarity enforced by NLSM/Kahler/D-term data, not holomorphic W alone",
            "microscopic_level": "conditional constrained EFT",
        },
        {
            "scenario": "triplet_only_8_pairs",
            "spin10_status": "incomplete post-breaking remnant",
            "extra_T": 0.0,
            "projected_threshold_l2": projected_l2(triplet_only_delta),
            "universal_threshold_shift": float(np.mean(triplet_only_delta)),
            "matrix_condition": "arbitrary unitary triplet block allowed",
            "microscopic_level": "fails threshold unless locked at MG",
        },
    ]
    out = []
    for item in scenarios:
        sum_t = bench["baseline_sum_T"] + float(item["extra_T"])
        b10 = sum_t - GAUGE_CONTRIBUTION
        alpha_r = alpha_inv_at_ratio(bench["alphaG_inv_MG"], b10, bench["R_target"])
        lp = landau_ratio(bench["alphaG_inv_MG"], b10)
        item = dict(item)
        item.update(
            {
                "sum_T_above_MG": sum_t,
                "b10_above_MG": b10,
                "alpha_inv_at_R200": alpha_r,
                "landau_ratio_over_MG": lp,
                "reaches_R200_perturbatively": alpha_r > 0.0,
                "passes_projected_threshold": item["projected_threshold_l2"] <= bench["r200_threshold_tol"],
            }
        )
        if item["scenario"] == "single_species_8_complete_10s":
            item["status"] = "FAIL_MATRIX_AND_R200"
        elif item["scenario"] == "two_species_16_complete_10s":
            item["status"] = "FAIL_R200"
        elif item["scenario"] == "post_spin10_8_su5_complete_5pairs":
            item["status"] = "PASS_CONDITIONAL_POST_GUT_EFT"
        elif item["scenario"] == "constrained_composite_unitary_link":
            item["status"] = "PASS_CONDITIONAL_CONSTRAINED"
        else:
            item["status"] = "FAIL_THRESHOLD"
        out.append(item)
    return out


def build() -> dict[str, Any]:
    bench = load_benchmark()
    matrix_rows = matrix_obstruction_rows()
    scenarios = scenario_rows(bench)
    all_single_species_fail = all(not row["single_species_10_mass_possible"] for row in matrix_rows)
    all_holo_orth_fail = all(not row["holomorphic_F_orthogonality_possible"] for row in matrix_rows)
    conditional = [row for row in scenarios if row["status"].startswith("PASS_CONDITIONAL")]
    verdict = {
        "single_species_10_mass_fails_all_targets": all_single_species_fail,
        "holomorphic_orthogonality_fails_all_targets": all_holo_orth_fail,
        "conditional_scenarios": [row["scenario"] for row in conditional],
        "elementary_spin10_two_species_alpha_inv_R200": next(
            row["alpha_inv_at_R200"] for row in scenarios if row["scenario"] == "two_species_16_complete_10s"
        ),
        "triplet_only_projected_threshold_l2": next(
            row["projected_threshold_l2"] for row in scenarios if row["scenario"] == "triplet_only_8_pairs"
        ),
        "interpretation": (
            "The unitary dilation is a valid threshold-locked EFT block, but an "
            "ordinary elementary Spin(10) completion is not yet obtained.  A "
            "single species of complete 10s fails the symmetric-mass test; two "
            "species allow the matrix but lose the R=200 perturbative reach; "
            "triplet-only remnants fail projected thresholds.  The viable local "
            "branches are therefore post-Spin(10) SU(5)-complete EFT or a "
            "constrained/composite unitary-link sector."
        ),
    }
    return {
        "note": "No web lookup used. Representation and symmetry audit for the unitary triplet-link dilation.",
        "benchmarks": bench,
        "matrix_obstruction_rows": matrix_rows,
        "scenario_rows": scenarios,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for row in rows for key in row.keys()})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Unitary-link Spin(10) embedding audit",
        "",
        "No web lookup was used.",
        "",
        "## Matrix Obstruction",
        "",
        "| target | symmetric residual rel | U^T U-I | single 10 species | holomorphic orthogonal |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary["matrix_obstruction_rows"]:
        lines.append(
            "| `{label}` | {symmetric_mass_residual_relative:.6e} | "
            "{complex_orthogonal_residual_2norm:.6e} | {single_species_10_mass_possible} | "
            "{holomorphic_F_orthogonality_possible} |".format(**row)
        )
    lines += [
        "",
        "## Representation Scenarios",
        "",
        "| scenario | status | extra T | alpha^-1(200MG) | PDelta | comment |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in summary["scenario_rows"]:
        lines.append(
            "| `{scenario}` | `{status}` | {extra_T:.1f} | {alpha_inv_at_R200:.6e} | "
            "{projected_threshold_l2:.6e} | {matrix_condition} |".format(**row)
        )
    v = summary["verdict"]
    lines += [
        "",
        "## Verdict",
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(OUT / "matrix_obstructions.csv", summary["matrix_obstruction_rows"])
    write_csv(OUT / "scenario_scorecard.csv", summary["scenario_rows"])
    write_report(summary)
    v = summary["verdict"]
    print("Unitary-link Spin(10) embedding audit")
    print(f"  single-species 10 mass fails all targets: {v['single_species_10_mass_fails_all_targets']}")
    print(f"  holomorphic orthogonality fails all targets: {v['holomorphic_orthogonality_fails_all_targets']}")
    print(f"  conditional scenarios: {', '.join(v['conditional_scenarios'])}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
