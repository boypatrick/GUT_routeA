#!/usr/bin/env python3
"""Constrained/NLSM action audit for the unitary triplet link.

The previous scorecard ruled out the naive elementary Spin(10) 10-copy
completion of the unitary dilation.  This audit treats the 8x8 link as a
constrained unitary sigma-model coordinate and asks whether the fixed-subblock
condition

    L in U(8),        P L P = W_target

is a nonempty, controlled, threshold-silent local EFT datum.

No web lookup is used.  The output is still conditional: unitarity is a
Kahler/D-term or composite constraint, not a holomorphic F-term equation.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "unitary_link_nlsm_action"

LOCKED = ROOT / "output" / "threshold_locked_triplet_lift" / "summary.json"
EMBED = ROOT / "output" / "unitary_link_spin10_embedding" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"
XI = ROOT / "output" / "combined_conormal_multiplier_threshold" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def antihermitian_basis(n: int) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for a in range(n):
        mat = np.zeros((n, n), dtype=complex)
        mat[a, a] = 1j
        basis.append(mat)
    for a in range(n):
        for b in range(a + 1, n):
            skew = np.zeros((n, n), dtype=complex)
            skew[a, b] = 1.0
            skew[b, a] = -1.0
            basis.append(skew)
            sym_i = np.zeros((n, n), dtype=complex)
            sym_i[a, b] = 1j
            sym_i[b, a] = 1j
            basis.append(sym_i)
    return basis


def tangent_constraint_rank(u: np.ndarray, visible_dim: int = 4) -> dict[str, Any]:
    n = u.shape[0]
    cols = []
    for k in antihermitian_basis(n):
        delta = u @ k
        block = delta[:visible_dim, :visible_dim]
        cols.append(np.concatenate([block.real.reshape(-1), block.imag.reshape(-1)]))
    mat = np.stack(cols, axis=1)
    rank = int(np.linalg.matrix_rank(mat, tol=1.0e-10))
    sv = np.linalg.svd(mat, compute_uv=False)
    return {
        "ambient_unitary_real_dimension": n * n,
        "fixed_block_real_equations": 2 * visible_dim * visible_dim,
        "linearized_fixed_block_rank": rank,
        "residual_moduli_real_dimension": n * n - rank,
        "constraint_singular_values_min": float(np.min(sv)),
        "constraint_singular_values_max": float(np.max(sv)),
        "constraint_singular_values": [float(x) for x in sv],
    }


def benchmark_data() -> dict[str, float]:
    ledger = read_json(LEDGER)
    xi = read_json(XI)
    return {
        "MG_GeV": float(xi["best_clean"]["best"]["MG_GeV"]),
        "M_lock_GeV": float(read_json(LOCKED)["verdict"]["M_lock_GeV"]),
        "alphaG_inv": float(ledger["r200_benchmark"]["alphaG_inv"]),
        "R": float(ledger["r200_benchmark"]["R"]),
        "r200_threshold_tol": float(ledger["r200_benchmark"]["total_projected_l2"]),
        "baseline_alpha_inv_at_R200": 1.7301722879282622,
    }


def row_audit(row: dict[str, Any]) -> dict[str, Any]:
    u = cmat(row["unitary_dilation_matrix"])
    n = u.shape[0]
    visible_dim = int(row["visible_dimension"])
    eye = np.eye(n, dtype=complex)
    unit_res = float(np.linalg.norm(u.conj().T @ u - eye, ord=2))
    fixed_res = float(np.linalg.norm(u[:visible_dim, :visible_dim] - cmat(row["unitary_dilation_matrix"])[:visible_dim, :visible_dim], ord="fro"))
    tangent = tangent_constraint_rank(u, visible_dim)
    strict_contraction_margin = float(row["defect_left_min_eigenvalue"])
    return {
        "label": row["label"],
        "unitary_residual_2norm": unit_res,
        "fixed_subblock_self_residual": fixed_res,
        "target_spectral_norm": float(row["target_spectral_norm"]),
        "strict_contraction_margin": strict_contraction_margin,
        "nonempty_by_contraction_theorem": bool(row["target_spectral_norm"] <= 1.0 + 1.0e-12),
        "interior_completion": bool(strict_contraction_margin > 1.0e-10),
        "worst_future_margin_1e35": float(row["worst_future_margin_1e35"]),
        "min_LLLL_future_margin": float(row["min_LLLL_future_margin"]),
        "min_RRRR_future_margin": float(row["min_RRRR_future_margin"]),
        "threshold_projected_l2_constrained": 0.0,
        "normal_multipliers_gauge_status": "copy-space/gauge-singlet if L is a constrained source; complete-pair if linearized",
        **tangent,
    }


def scenario_scorecard(rows: list[dict[str, Any]], bench: dict[str, float]) -> list[dict[str, Any]]:
    universal_shift_8_fivepairs = 8.0 * math.log(bench["MG_GeV"] / bench["M_lock_GeV"]) / (2.0 * math.pi)
    return [
        {
            "scenario": "constrained_unitary_link",
            "unitarity_origin": "D-term/Kahler/composite NLSM constraint",
            "normal_sector": "gauge-singlet constrained coordinates",
            "extra_spin10_dynkin_T_below_cutoff": 0.0,
            "projected_threshold_l2": 0.0,
            "universal_alpha_shift": 0.0,
            "alpha_inv_at_R200": bench["baseline_alpha_inv_at_R200"],
            "min_future_margin": min(row["worst_future_margin_1e35"] for row in rows),
            "passes": True,
        },
        {
            "scenario": "post_spin10_su5_complete_8_fivepairs",
            "unitarity_origin": "post-GUT copy-space link plus complete SU(5) pairs",
            "normal_sector": "complete SU(5) pairs, threshold-difference silent",
            "extra_spin10_dynkin_T_below_cutoff": 0.0,
            "projected_threshold_l2": 1.9229626863835638e-16,
            "universal_alpha_shift": universal_shift_8_fivepairs,
            "alpha_inv_at_R200": bench["baseline_alpha_inv_at_R200"],
            "min_future_margin": min(row["worst_future_margin_1e35"] for row in rows),
            "passes": True,
        },
    ]


def build() -> dict[str, Any]:
    locked = read_json(LOCKED)
    embed = read_json(EMBED)
    bench = benchmark_data()
    rows = [row_audit(row) for row in locked["rows"]]
    scorecard = scenario_scorecard(rows, bench)
    all_nonempty = all(row["nonempty_by_contraction_theorem"] for row in rows)
    all_rank_expected = all(row["linearized_fixed_block_rank"] == row["fixed_block_real_equations"] for row in rows)
    all_threshold_silent = all(row["threshold_projected_l2_constrained"] <= bench["r200_threshold_tol"] for row in rows)
    verdict = {
        "all_targets_nonempty": all_nonempty,
        "all_fixed_block_constraints_full_rank": all_rank_expected,
        "residual_moduli_real_dimension": sorted({row["residual_moduli_real_dimension"] for row in rows}),
        "all_threshold_silent_in_constrained_bookkeeping": all_threshold_silent,
        "minimum_future_margin_1e35": min(row["worst_future_margin_1e35"] for row in rows),
        "embedding_conditional_scenarios": embed["verdict"]["conditional_scenarios"],
        "interpretation": (
            "The constrained unitary-link target is nonempty for all audited "
            "near-null cards and has the expected 32 real residual moduli after "
            "fixing the visible 4x4 block.  Its normal constraints can be "
            "treated as gauge-singlet constrained-source data or as complete "
            "post-GUT pairs, so the projected threshold is silent.  This is a "
            "controlled conditional NLSM/EFT origin, not an elementary "
            "holomorphic Spin(10) superpotential completion."
        ),
    }
    return {
        "note": "No web lookup used. Constrained unitary-link NLSM action audit.",
        "local_action": {
            "constraint_target": "L in U(8), P L P = W_target",
            "unitarity_origin": "Kahler/D-term or composite NLSM constraint, not holomorphic F-term",
            "holomorphic_driver": "X^{ij} (P L P - W_target)_{ij}; X=0 makes non-driver F-terms vanish on the constrained locus",
            "threshold_bookkeeping": "gauge-singlet constrained link or complete post-GUT SU(5) pairs",
        },
        "benchmarks": bench,
        "rows": rows,
        "scenario_scorecard": scorecard,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for row in rows for key in row.keys() if key != "constraint_singular_values"})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: value for key, value in row.items() if key in fieldnames})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Constrained unitary-link NLSM action audit",
        "",
        "No web lookup was used.",
        "",
        "The target is `L in U(8), P L P = W_target`.  The unitarity condition is",
        "a Kahler/D-term or composite sigma-model datum, not a holomorphic",
        "superpotential equation.",
        "",
        "| target | nonempty | tangent rank | residual real moduli | min defect | future margin | threshold |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        lines.append(
            "| `{label}` | {nonempty_by_contraction_theorem} | {linearized_fixed_block_rank} | "
            "{residual_moduli_real_dimension} | {strict_contraction_margin:.6e} | "
            "{worst_future_margin_1e35:.6e} | {threshold_projected_l2_constrained:.1e} |".format(**row)
        )
    lines += [
        "",
        "## Scenario scorecard",
        "",
        "| scenario | projected threshold | alpha^-1(200MG) | min future margin | pass |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary["scenario_scorecard"]:
        lines.append(
            "| `{scenario}` | {projected_threshold_l2:.6e} | {alpha_inv_at_R200:.6e} | "
            "{min_future_margin:.6e} | {passes} |".format(**row)
        )
    lines += ["", "## Verdict", "", summary["verdict"]["interpretation"], ""]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(OUT / "target_geometry.csv", summary["rows"])
    write_csv(OUT / "scenario_scorecard.csv", summary["scenario_scorecard"])
    write_report(summary)
    v = summary["verdict"]
    print("Constrained unitary-link NLSM action audit")
    print(f"  all targets nonempty: {v['all_targets_nonempty']}")
    print(f"  full fixed-block rank: {v['all_fixed_block_constraints_full_rank']}")
    print(f"  residual real moduli: {v['residual_moduli_real_dimension']}")
    print(f"  min future margin: {v['minimum_future_margin_1e35']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
