#!/usr/bin/env python3
"""Kappa_54 global scan and gauge-Yukawa fixed-line diagnostic.

No web lookup is used.  This script promotes kappa_54 to an explicit scan axis
and keeps the locked sub-branch

    kappa_54 = 1  <=>  lambda/g = sqrt(2)

as a named benchmark.  It also checks whether the same relation can be obtained
as an ordinary one-loop N=1 SO(10) gauge-Yukawa fixed ratio for

    W54 = lambda/3 Tr(S^3).

The answer is negative in the ordinary field content: the one-loop fixed ratios
are larger than sqrt(2), and sqrt(2) would require b_10=-26.4, more
asymptotically free than pure SO(10), hence impossible with ordinary chiral
matter.  This makes the locked branch a boundary condition unless an extended
symmetry/fixed-line mechanism beyond minimal N=1 is supplied.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_corrected_downstream as down  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "kappa54_global_scan"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
CORRECTED_SCAN_CSV = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"

N = 10
KAPPA_GRID = [
    0.5,
    0.75,
    0.9,
    0.95,
    0.99,
    1.0,
    1.01,
    1.05,
    1.1,
    1.25,
    1.5,
    2.0,
    3.0,
    5.0,
]
R_GRID = [50.0, 200.0]

B_NON_GOLDSTONE = np.array([34.0 / 5.0, 6.0, 8.0], dtype=float)
T_INDEX = {
    "10": 1.0,
    "16": 2.0,
    "45": 8.0,
    "54": 12.0,
    "120": 28.0,
    "126bar": 35.0,
    "210": 56.0,
}
GAUGE_CONTRIBUTION = 24.0
C2_54 = 10.0


def orthonormal_54_basis() -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(N):
        for j in range(i + 1, N):
            mat = np.zeros((N, N), dtype=float)
            mat[i, j] = 1.0 / math.sqrt(2.0)
            mat[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(mat)

    diagonal_vectors: list[np.ndarray] = []
    for i in range(N - 1):
        vec = np.zeros(N, dtype=float)
        vec[i] = 1.0
        vec[N - 1] = -1.0
        diagonal_vectors.append(vec)

    orth: list[np.ndarray] = []
    for vec in diagonal_vectors:
        work = vec.copy()
        for prev in orth:
            work -= np.dot(prev, work) * prev
        work /= np.linalg.norm(work)
        orth.append(work)
    for vec in orth:
        basis.append(np.diag(vec))

    gram = np.array([[np.trace(a @ b) for b in basis] for a in basis])
    if np.max(np.abs(gram - np.eye(54))) > 1.0e-12:
        raise RuntimeError("54 basis is not orthonormal")
    return basis


def cubic_invariant_contraction() -> dict[str, Any]:
    basis = orthonormal_54_basis()
    dim = len(basis)
    d_tensor = np.zeros((dim, dim, dim), dtype=float)
    for a, mat_a in enumerate(basis):
        for b, mat_b in enumerate(basis):
            ab = mat_a @ mat_b
            for c, mat_c in enumerate(basis):
                d_tensor[a, b, c] = np.trace(ab @ mat_c)
    contraction = np.einsum("apq,bpq->ab", d_tensor, d_tensor)
    eigvals = np.linalg.eigvalsh(contraction)
    a_d = float(np.trace(contraction) / dim)
    return {
        "basis_dimension": dim,
        "d_apq_d_bpq_coefficient": a_d,
        "max_off_diagonal_error": float(np.max(np.abs(contraction - a_d * np.eye(dim)))),
        "eigenvalue_counts": {str(key): val for key, val in Counter(round(x, 12) for x in eigvals).items()},
        "passes": abs(a_d - 14.0 / 5.0) < 1.0e-12
        and float(np.max(np.abs(contraction - a_d * np.eye(dim)))) < 1.0e-11,
    }


def sum_t(counts: dict[str, int]) -> float:
    return float(sum(T_INDEX[rep] * count for rep, count in counts.items()))


def fixedline_scenarios(a_d: float) -> list[dict[str, Any]]:
    scenarios = [
        {
            "name": "isolated_single_54",
            "counts": {"54": 1},
            "interpretation": "Only the single 54_H is active above M_G.",
        },
        {
            "name": "single_54_plus_three_45_projector",
            "counts": {"54": 1, "45": 3},
            "interpretation": "The low-index 54_H plus the three 45 mediator fields.",
        },
        {
            "name": "single_54_plus_minimal_flavor_higgs",
            "counts": {"54": 1, "16": 3, "10": 1, "120": 1, "126bar": 1},
            "interpretation": "Single 54_H plus minimal matter/Yukawa Higgs sector.",
        },
        {
            "name": "single_54_plus_projector_plus_flavor_higgs",
            "counts": {"54": 1, "45": 3, "16": 3, "10": 1, "120": 1, "126bar": 1},
            "interpretation": "Low-index projector sector plus minimal flavor/Higgs fields.",
        },
        {
            "name": "large_54_210_alignment_sector",
            "counts": {"45": 3, "54": 6, "210": 2},
            "interpretation": "Previously rejected large propagating source/driver tower.",
        },
    ]
    rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        total_t = sum_t(scenario["counts"])
        b10 = total_t - GAUGE_CONTRIBUTION
        # W=lambda/3 Tr S^3 = (1/6)Y_abc phi_a phi_b phi_c, Y=2 lambda d.
        # gamma = [(1/2)Y_apq Y_bpq - 2g^2 C2(54)]/(16pi^2)
        #       = [(2 A_d lambda^2) - 20 g^2]/(16pi^2).
        # beta_lambda/lambda = 3 gamma.
        yukawa_coeff = 6.0 * a_d
        gauge_coeff = 6.0 * C2_54
        fixed_lambda_over_g_sq = (gauge_coeff + b10) / yukawa_coeff
        fixed_lambda_over_g = math.sqrt(fixed_lambda_over_g_sq) if fixed_lambda_over_g_sq > 0.0 else math.nan
        rows.append(
            {
                "name": scenario["name"],
                "counts": scenario["counts"],
                "sum_T": total_t,
                "b10": b10,
                "beta_lambda_over_lambda": f"((6 A_d) lambda^2 - {gauge_coeff:.0f} g^2)/(16 pi^2)",
                "fixed_lambda_over_g_sq": fixed_lambda_over_g_sq,
                "fixed_lambda_over_g": fixed_lambda_over_g,
                "fixed_kappa54": fixed_lambda_over_g / math.sqrt(2.0)
                if math.isfinite(fixed_lambda_over_g)
                else math.nan,
                "hits_locked_sqrt2": abs(fixed_lambda_over_g - math.sqrt(2.0)) < 1.0e-9
                if math.isfinite(fixed_lambda_over_g)
                else False,
                "interpretation": scenario["interpretation"],
            }
        )

    b_needed = 2.0 * (6.0 * a_d) - 6.0 * C2_54
    rows.append(
        {
            "name": "required_b10_for_locked_sqrt2",
            "counts": {},
            "sum_T": b_needed + GAUGE_CONTRIBUTION,
            "b10": b_needed,
            "fixed_lambda_over_g_sq": 2.0,
            "fixed_lambda_over_g": math.sqrt(2.0),
            "fixed_kappa54": 1.0,
            "hits_locked_sqrt2": True,
            "interpretation": (
                "This b10 would be required for ordinary one-loop N=1 flow to fix "
                "lambda/g=sqrt(2). It is below the pure-gauge value b10=-24 and "
                "therefore impossible with ordinary positive-index chiral matter."
            ),
        }
    )
    return rows


def load_vacuum() -> dict[str, Any]:
    return json.loads(VACUUM.read_text(encoding="utf-8"))


def iter_corrected_rows() -> list[dict[str, str]]:
    with CORRECTED_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def mediator_delta_for_r(payload: dict[str, Any], r_value: float) -> np.ndarray:
    card = payload["benchmark_cards"][f"{r_value:.1f}"]
    return np.array(card["mediator_heavy_threshold"]["delta_full"], dtype=float)


def delta_54(kappa_54: float) -> np.ndarray:
    return B_NON_GOLDSTONE * math.log(1.0 / kappa_54) / (2.0 * math.pi)


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


def global_scan_rows() -> list[dict[str, Any]]:
    payload = load_vacuum()
    source_rows = iter_corrected_rows()
    rows: list[dict[str, Any]] = []
    for r_value in R_GRID:
        med_delta = mediator_delta_for_r(payload, r_value)
        for kappa in KAPPA_GRID:
            d54 = delta_54(kappa)
            scan = down.scan_cached_with_delta(source_rows, med_delta + d54)
            best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
            alpha_inv = float(best["alphaG_inv"])
            g = math.sqrt(4.0 * math.pi / alpha_inv)
            lambda_over_g = math.sqrt(2.0) * kappa
            rows.append(
                {
                    "R": r_value,
                    "kappa_54": kappa,
                    "locked_subbranch": abs(kappa - 1.0) < 1.0e-12,
                    "lambda_over_g": lambda_over_g,
                    "lambda_value": lambda_over_g * g,
                    "lambda54_perturbative_lt_sqrt4pi": lambda_over_g * g < math.sqrt(4.0 * math.pi),
                    "g_G": g,
                    "delta54": d54.tolist(),
                    "projected_l2_delta54": float(np.linalg.norm(base.PROJECTOR @ d54)),
                    "total_projected_l2": float(np.linalg.norm(base.PROJECTOR @ (med_delta + d54))),
                    "clean_threshold_1e_minus_2": float(np.linalg.norm(base.PROJECTOR @ d54)) < 1.0e-2,
                    "safe_points": int(scan["safe_points"]),
                    "safe_single_scale_factor2_points": int(scan["safe_single_scale_factor2_points"]),
                    "best_is_safe": scan["best_safe"] is not None,
                    "best": compact_best(best),
                }
            )
    return rows


def make_benchmark_cards(rows: list[dict[str, Any]]) -> dict[str, Any]:
    cards: dict[str, Any] = {}
    for r_value in R_GRID:
        locked = next(row for row in rows if row["R"] == r_value and row["locked_subbranch"])
        cards[f"R{int(r_value)}_kappa54_locked"] = {
            "R": r_value,
            "kappa_54": 1.0,
            "lambda_over_g": math.sqrt(2.0),
            "lambda_value": locked["lambda_value"],
            "g_G": locked["g_G"],
            "Delta54": [0.0, 0.0, 0.0],
            "interpretation": "locked sub-branch lambda=sqrt(2)g, no single-54 chiral threshold",
            "best": locked["best"],
            "safe_points": locked["safe_points"],
            "total_projected_l2": locked["total_projected_l2"],
            "projected_l2_delta54": locked["projected_l2_delta54"],
        }
    best_any = min(
        [row for row in rows if row["best_is_safe"]],
        key=lambda row: (float(row["best"]["score"]), row["total_projected_l2"]),
    )
    cards["best_safe_over_kappa54_grid"] = {
        "R": best_any["R"],
        "kappa_54": best_any["kappa_54"],
        "lambda_over_g": best_any["lambda_over_g"],
        "lambda_value": best_any["lambda_value"],
        "g_G": best_any["g_G"],
        "best": best_any["best"],
        "safe_points": best_any["safe_points"],
        "total_projected_l2": best_any["total_projected_l2"],
        "projected_l2_delta54": best_any["projected_l2_delta54"],
        "note": "This may prefer shifted kappa_54 because the score rewards smaller fitted logs; it is not the locked clean branch.",
    }
    best_pert = min(
        [row for row in rows if row["best_is_safe"] and row["lambda54_perturbative_lt_sqrt4pi"]],
        key=lambda row: (float(row["best"]["score"]), row["total_projected_l2"]),
    )
    cards["best_perturbative_lambda54_over_kappa54_grid"] = {
        "R": best_pert["R"],
        "kappa_54": best_pert["kappa_54"],
        "lambda_over_g": best_pert["lambda_over_g"],
        "lambda_value": best_pert["lambda_value"],
        "g_G": best_pert["g_G"],
        "best": best_pert["best"],
        "safe_points": best_pert["safe_points"],
        "total_projected_l2": best_pert["total_projected_l2"],
        "projected_l2_delta54": best_pert["projected_l2_delta54"],
        "note": "Best score subject to lambda54<sqrt(4pi); still not necessarily a clean-threshold branch.",
    }
    best_clean = min(
        [
            row
            for row in rows
            if row["best_is_safe"]
            and row["lambda54_perturbative_lt_sqrt4pi"]
            and row["clean_threshold_1e_minus_2"]
        ],
        key=lambda row: (float(row["best"]["score"]), row["total_projected_l2"]),
    )
    cards["best_clean_threshold_1e_minus_2_grid"] = {
        "R": best_clean["R"],
        "kappa_54": best_clean["kappa_54"],
        "lambda_over_g": best_clean["lambda_over_g"],
        "lambda_value": best_clean["lambda_value"],
        "g_G": best_clean["g_G"],
        "best": best_clean["best"],
        "safe_points": best_clean["safe_points"],
        "total_projected_l2": best_clean["total_projected_l2"],
        "projected_l2_delta54": best_clean["projected_l2_delta54"],
        "note": "Best score subject to lambda54<sqrt(4pi) and ||P Delta54||<1e-2.",
    }
    return cards


def write_scan_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "R",
        "kappa_54",
        "locked_subbranch",
        "lambda_over_g",
        "lambda_value",
        "lambda54_perturbative_lt_sqrt4pi",
        "g_G",
        "projected_l2_delta54",
        "total_projected_l2",
        "clean_threshold_1e_minus_2",
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
    with (OUT / "kappa54_global_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            best = row["best"]
            writer.writerow(
                {
                    "R": row["R"],
                    "kappa_54": row["kappa_54"],
                    "locked_subbranch": row["locked_subbranch"],
                    "lambda_over_g": row["lambda_over_g"],
                    "lambda_value": row["lambda_value"],
                    "lambda54_perturbative_lt_sqrt4pi": row["lambda54_perturbative_lt_sqrt4pi"],
                    "g_G": row["g_G"],
                    "projected_l2_delta54": row["projected_l2_delta54"],
                    "total_projected_l2": row["total_projected_l2"],
                    "clean_threshold_1e_minus_2": row["clean_threshold_1e_minus_2"],
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


def write_fixedline_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "name",
        "sum_T",
        "b10",
        "fixed_lambda_over_g_sq",
        "fixed_lambda_over_g",
        "fixed_kappa54",
        "hits_locked_sqrt2",
        "interpretation",
    ]
    with (OUT / "gauge_yukawa_fixedline.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def build_payload() -> dict[str, Any]:
    cubic = cubic_invariant_contraction()
    fixed = fixedline_scenarios(float(cubic["d_apq_d_bpq_coefficient"]))
    scan = global_scan_rows()
    cards = make_benchmark_cards(scan)
    locked200 = cards["R200_kappa54_locked"]
    ordinary_hits = [row for row in fixed if row["hits_locked_sqrt2"] and row["name"] != "required_b10_for_locked_sqrt2"]
    return {
        "note": "No web lookup used. Explicit kappa54 scan and gauge-Yukawa fixed-line diagnostic.",
        "cubic_invariant": cubic,
        "one_loop_formula": {
            "d_contraction": "d_apq d_bpq = (14/5) delta_ab",
            "gamma_54": "[(28/5) lambda^2 - 20 g^2]/(16 pi^2)",
            "beta_lambda_over_lambda": "[(84/5) lambda^2 - 60 g^2]/(16 pi^2)",
            "beta_g_over_g": "b10 g^2/(16 pi^2)",
            "fixed_ratio": "(lambda/g)^2 = (60+b10)/(84/5)",
            "locked_target": "(lambda/g)^2 = 2",
        },
        "fixedline_rows": fixed,
        "scan_rows": scan,
        "benchmark_cards": cards,
        "verdict": {
            "locked_R200": locked200,
            "ordinary_N1_fixedline_hits_sqrt2": bool(ordinary_hits),
            "b10_required_for_sqrt2": next(row for row in fixed if row["name"] == "required_b10_for_locked_sqrt2")["b10"],
            "fixedline_no_go": (
                "Ordinary one-loop N=1 SO(10) gauge-Yukawa flow cannot explain "
                "lambda/g=sqrt(2) with positive-index chiral matter. It would require "
                "b10=-26.4, below the pure-gauge value -24."
            ),
            "recommended_interpretation": (
                "Keep kappa54 as an explicit scan/card parameter and mark kappa54=1 "
                "as the locked lambda=sqrt(2)g sub-branch. A stronger claim requires an "
                "extended symmetry or UV fixed-line mechanism beyond ordinary N=1."
            ),
        },
    }


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Kappa54 global scan and fixed-line diagnostic")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Locked benchmark cards")
    lines.append("")
    lines.append("| branch | kappa54 | lambda/g | lambda | safe points | M_Sigma3 [GeV] | tau_d6 [yr] |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for key in ["R50_kappa54_locked", "R200_kappa54_locked"]:
        card = payload["benchmark_cards"][key]
        best = card["best"]
        lines.append(
            f"| {key} | {card['kappa_54']:.1f} | {card['lambda_over_g']:.9f} | "
            f"{card['lambda_value']:.9f} | {card['safe_points']} | "
            f"{best['M_Sigma3_GeV']:.6e} | {best['tau_dim6_years']:.6e} |"
        )
    lines.append("")
    lines.append("## Non-locked scan diagnostics")
    lines.append("")
    lines.append("| card | R | kappa54 | lambda/g | lambda | ||P Delta54|| | safe points | score |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for key in [
        "best_safe_over_kappa54_grid",
        "best_perturbative_lambda54_over_kappa54_grid",
        "best_clean_threshold_1e_minus_2_grid",
    ]:
        card = payload["benchmark_cards"][key]
        best = card["best"]
        projected = card.get("projected_l2_delta54", float("nan"))
        lines.append(
            f"| {key} | {card['R']:.0f} | {card['kappa_54']:.2f} | "
            f"{card['lambda_over_g']:.6f} | {card['lambda_value']:.6f} | "
            f"{projected:.6e} | {card['safe_points']} | {best['score']:.6f} |"
        )
    lines.append("")
    lines.append(
        "The score-only optimum can move to large kappa54, but it is not the clean branch. "
        "The clean-threshold diagnostic keeps ||P Delta54||<1e-2 and remains near kappa54=1."
    )
    lines.append("")
    lines.append("## One-loop fixed-line diagnostic")
    lines.append("")
    lines.append("The cubic invariant audit gives `d_apq d_bpq = 14/5 delta_ab`, hence")
    lines.append("")
    lines.append("```text")
    lines.append("beta_lambda/lambda = [(84/5) lambda^2 - 60 g^2]/(16 pi^2)")
    lines.append("(lambda/g)_fixed^2 = (60+b10)/(84/5)")
    lines.append("```")
    lines.append("")
    lines.append("| scenario | b10 | fixed lambda/g | fixed kappa54 |")
    lines.append("|---|---:|---:|---:|")
    for row in payload["fixedline_rows"]:
        if row["name"] != "required_b10_for_locked_sqrt2":
            lines.append(
                f"| {row['name']} | {row['b10']:.1f} | "
                f"{row['fixed_lambda_over_g']:.6f} | {row['fixed_kappa54']:.6f} |"
            )
    need = next(row for row in payload["fixedline_rows"] if row["name"] == "required_b10_for_locked_sqrt2")
    lines.append("")
    lines.append(f"To hit `lambda/g=sqrt(2)` with the ordinary one-loop formula would require `b10={need['b10']:.1f}`.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"]["fixedline_no_go"])
    lines.append(payload["verdict"]["recommended_interpretation"])
    lines.append("")
    (OUT / "kappa54_global_fixedline_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    write_scan_csv(payload["scan_rows"])
    write_fixedline_csv(payload["fixedline_rows"])
    (OUT / "kappa54_global_scan_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    # Write locked cards separately for easy inclusion in reproducibility bundle.
    (OUT / "kappa54_benchmark_cards.json").write_text(
        json.dumps(payload["benchmark_cards"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Kappa54 global scan and fixed-line diagnostic complete")
    print(f"  cubic A_d={payload['cubic_invariant']['d_apq_d_bpq_coefficient']:.9f}")
    print(f"  locked R200 lambda={payload['benchmark_cards']['R200_kappa54_locked']['lambda_value']:.9f}")
    print(f"  ordinary N=1 fixed line hits sqrt2={payload['verdict']['ordinary_N1_fixedline_hits_sqrt2']}")
    print(f"  b10 required for sqrt2={payload['verdict']['b10_required_for_sqrt2']:.3f}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
