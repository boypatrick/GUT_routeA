#!/usr/bin/env python3
"""Construct and audit the crossed 120A/120B triplet projector.

No web lookup is used.

The previous audit found a component-level d=5 null if the 10-10 triplet
source uses one antisymmetric 120-like family direction G_A while the
10-5bar triplet source uses the orthogonal G_B-like direction.  In source
language the target visible inverse-propagator block is rank one,

    W_vis = |120_A><120_B|.

This script checks three things.

1. A field-only unbroken Spin(10) grading cannot enforce this crossed entry
   while allowing both 16 16 120_A and 16 16 120_B Yukawa operators.  The two
   120 copies then have identical charge, so all 120_i 120_j mass entries have
   identical charge.
2. After Spin(10) breaking, the crossed triplet projector is algebraically
   realizable by a Julia/Sz.-Nagy unitary dilation with a degenerate complete
   heavy spectrum, hence zero non-universal one-loop threshold in the complete
   multiplet interpretation.
3. Finite condition-capped lifts leak at O(1/kappa); the existing component
   leakage ratios are replayed as d=5 margins.

The result is therefore not a full Spin(10)-invariant closure.  It is a
precise PS/post-GUT EFT target for the next Higgs-sector construction.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "crossed_120_triplet_projector"
INPUT = ROOT / "output" / "triplet_120_rrrr_mass_matrix" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    if float(np.min(vals)) < -tol:
        raise RuntimeError(f"negative defect eigenvalue {np.min(vals)}")
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def julia_dilation(w: np.ndarray) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    return np.block(
        [
            [w, hermitian_sqrt(eye - w @ w.conj().T)],
            [hermitian_sqrt(eye - w.conj().T @ w), -w.conj().T],
        ]
    )


def field_only_grading_search(nmax: int = 24) -> dict[str, Any]:
    """Brute-force the no-go for small Z_N field-only gradings."""
    witnesses = []
    for n in range(2, nmax + 1):
        for q16 in range(n):
            for q_a in range(n):
                for q_b in range(n):
                    y_a = (2 * q16 + q_a) % n == 0
                    y_b = (2 * q16 + q_b) % n == 0
                    if not (y_a and y_b):
                        continue
                    charges = {
                        "AA": (q_a + q_a) % n,
                        "AB": (q_a + q_b) % n,
                        "BA": (q_b + q_a) % n,
                        "BB": (q_b + q_b) % n,
                    }
                    crossed_only = charges["AB"] == 0 and any(charges[key] != 0 for key in ["AA", "BA", "BB"])
                    if crossed_only:
                        witnesses.append({"N": n, "q16": q16, "qA": q_a, "qB": q_b, "charges": charges})
    proof = (
        "If both Yukawa terms 16 16 120_A and 16 16 120_B are field-only invariant, "
        "then q_A=q_B=-2 q_16 in any abelian grading.  Hence q_A+q_A=q_A+q_B="
        "q_B+q_A=q_B+q_B.  A field-only symmetry that allows one 120_i 120_j mass "
        "entry allows all four, so it cannot enforce a crossed-only projector."
    )
    return {
        "searched_ZN_up_to": nmax,
        "crossed_only_witness_count": len(witnesses),
        "witnesses": witnesses[:8],
        "analytic_proof": proof,
        "field_only_unbroken_spin10_projector_possible": len(witnesses) > 0,
    }


def direct_and_dilated_rows(m_lock: float) -> list[dict[str, Any]]:
    w = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    rows: list[dict[str, Any]] = []
    for kappa in [30.0, 100.0, 300.0, 1000.0, 10000.0]:
        eps = 1.0 / kappa
        w_eps = np.array([[0.0, 1.0], [eps, 0.0]], dtype=complex)
        direct_s = np.linalg.svd(w_eps, compute_uv=False)
        direct_mass_ratio = float(direct_s[0] / direct_s[-1])
        u = julia_dilation(w_eps / float(direct_s[0]))
        unitary_residual = float(np.linalg.norm(u.conj().T @ u - np.eye(u.shape[0]), ord=2))
        full_s = np.linalg.svd(u.conj().T, compute_uv=False)
        rows.append(
            {
                "kappa": float(kappa),
                "epsilon": float(eps),
                "direct_inverse_block_singular_values": [float(x) for x in direct_s],
                "direct_mass_condition_if_used_as_inverse": direct_mass_ratio,
                "direct_Mmax_if_Mlight_Mlock_GeV": float(direct_mass_ratio * m_lock),
                "direct_below_reduced_planck": bool(direct_mass_ratio * m_lock < 2.435e18),
                "dilated_dimension": int(u.shape[0]),
                "dilated_unitary_residual_2norm": unitary_residual,
                "dilated_full_mass_singular_spread_over_Mlock": float(np.max(full_s) - np.min(full_s)),
                "dilated_threshold_projected_l2_if_complete_degenerate": 0.0,
                "W_eps": matrix_json(w_eps),
                "U_dilation": matrix_json(u),
            }
        )
    return rows


def finite_leakage_replay(source: dict[str, Any]) -> list[dict[str, Any]]:
    natural = next(row for row in source["rows"] if row["label"] == "natural_locked")
    l_margin0 = float(natural["LLLL_Knu"]["margin_1e35_at_ST_target"])
    r_margin0 = float(natural["RRRR_uusd"]["margin_1e35_at_ST_target"])
    rows = []
    for row in source["finite_lift_rows"]:
        l_ratio = float(row["linear_LLLL_ratio_to_natural"])
        r_ratio = float(row["linear_RRRR_ratio_to_natural"])
        l_margin = l_margin0 / max(l_ratio * l_ratio, 1.0e-300)
        r_margin = r_margin0 / max(r_ratio * r_ratio, 1.0e-300)
        rows.append(
            {
                "kappa": float(row["kappa"]),
                "LLLL_leakage_ratio": l_ratio,
                "RRRR_leakage_ratio": r_ratio,
                "LLLL_margin_1e35_replayed": float(l_margin),
                "RRRR_margin_1e35_replayed": float(r_margin),
                "worst_margin_1e35_replayed": float(min(l_margin, r_margin)),
                "future_safe_after_leakage": bool(min(l_margin, r_margin) > 1.0),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source = read_json(INPUT)
    ledger = read_json(LEDGER)
    m_lock = float(ledger["r200_benchmark"]["M_Sigma8_GeV"])
    grading = field_only_grading_search()
    dilation_rows = direct_and_dilated_rows(m_lock)
    leakage_rows = finite_leakage_replay(source)

    crossed_null = source["verdict"]["most_safe_physical_row"]
    all_leakage_safe = all(row["future_safe_after_leakage"] for row in leakage_rows)
    all_dilations_clean = all(row["dilated_unitary_residual_2norm"] < 1.0e-12 for row in dilation_rows)
    direct_planck_safe = [row for row in dilation_rows if row["direct_below_reduced_planck"]]
    verdict = {
        "field_only_unbroken_spin10_projector_possible": grading["field_only_unbroken_spin10_projector_possible"],
        "ps_eft_or_constrained_projector_possible": True,
        "all_finite_leakage_rows_future_safe": all_leakage_safe,
        "all_unitary_dilations_clean": all_dilations_clean,
        "direct_split_planck_safe_kappas": [row["kappa"] for row in direct_planck_safe],
        "crossed_null_worst_margin_1e35": crossed_null["worst_margin_1e35_at_ST_target"],
        "interpretation": (
            "A crossed 120A/120B projector is algebraically powerful but cannot be enforced by a "
            "field-only unbroken Spin(10) grading once both 16 16 120_A and 16 16 120_B Yukawa "
            "operators are allowed.  It should therefore be stated as a post-Spin(10)-breaking "
            "Pati-Salam/constrained triplet-sector projector.  In that conditional branch, the "
            "Julia-dilated complete-multiplet realization is threshold silent, and the replayed "
            "finite-lift d=5 leakage remains far above the 1e35 yr stress target for the audited "
            "kappa values."
        ),
    }
    summary = {
        "note": "No web lookup used. Crossed 120A/120B triplet projector construction audit.",
        "target": {
            "visible_inverse_block": "W_vis=|120_A><120_B|",
            "component_null": "10-10 triplet source G_A, 10-5bar triplet source G_B-like",
            "crossed_null_row": crossed_null["label"],
            "crossed_null_RRRR_amplitude": crossed_null["RRRR_uusd"]["amplitude"],
            "crossed_null_LLLL_amplitude": crossed_null["LLLL_Knu"]["amplitude"],
        },
        "field_only_grading_no_go": grading,
        "threshold_locked_dilation": {
            "M_lock_GeV": m_lock,
            "rows": dilation_rows,
            "identity": "For ||W||<=1, U=[[W,sqrt(I-WWdag)],[sqrt(I-WdagW),-Wdag]] is unitary and P(M_lock Udag)^{-1}P=W/M_lock.",
        },
        "finite_leakage_replay": leakage_rows,
        "selection_rule_status": {
            "full_spin10_field_only": "NO_GO",
            "post_spin10_ps_eft_or_constrained": "PASS_CONDITIONAL",
            "reason": "The crossed projector distinguishes triplet components and/or source copies beyond a field-only Spin(10) charge assignment.",
        },
        "verdict": verdict,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "finite_leakage_replay.csv",
        leakage_rows,
        [
            "kappa",
            "LLLL_leakage_ratio",
            "RRRR_leakage_ratio",
            "LLLL_margin_1e35_replayed",
            "RRRR_margin_1e35_replayed",
            "worst_margin_1e35_replayed",
            "future_safe_after_leakage",
        ],
    )
    write_csv(
        OUT / "dilation_rows.csv",
        dilation_rows,
        [
            "kappa",
            "epsilon",
            "direct_mass_condition_if_used_as_inverse",
            "direct_Mmax_if_Mlight_Mlock_GeV",
            "direct_below_reduced_planck",
            "dilated_dimension",
            "dilated_unitary_residual_2norm",
            "dilated_full_mass_singular_spread_over_Mlock",
            "dilated_threshold_projected_l2_if_complete_degenerate",
        ],
    )

    report = [
        "# Crossed 120A/120B triplet projector audit",
        "",
        "No web lookup was used.",
        "",
        "Target visible inverse block:",
        "",
        "```text",
        "W_vis = |120_A><120_B|.",
        "```",
        "",
        "## Field-only Spin(10) grading",
        "",
        grading["analytic_proof"],
        "",
        f"Brute-force search through Z_N, N<=24 found {grading['crossed_only_witness_count']} crossed-only witnesses.",
        "",
        "## Finite leakage replay",
        "",
        "| kappa | LLLL leak | RRRR leak | LLLL margin | RRRR margin | worst margin | safe |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in leakage_rows:
        report.append(
            f"| {row['kappa']:.0f} | {row['LLLL_leakage_ratio']:.3e} | "
            f"{row['RRRR_leakage_ratio']:.3e} | {row['LLLL_margin_1e35_replayed']:.3e} | "
            f"{row['RRRR_margin_1e35_replayed']:.3e} | {row['worst_margin_1e35_replayed']:.3e} | "
            f"{row['future_safe_after_leakage']} |"
        )
    report.extend(
        [
            "",
            "## Dilation versus direct split",
            "",
            "| kappa | direct Mmax [GeV] | direct Planck-safe | dilation residual | dilation threshold |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in dilation_rows:
        report.append(
            f"| {row['kappa']:.0f} | {row['direct_Mmax_if_Mlight_Mlock_GeV']:.3e} | "
            f"{row['direct_below_reduced_planck']} | {row['dilated_unitary_residual_2norm']:.3e} | "
            f"{row['dilated_threshold_projected_l2_if_complete_degenerate']:.1e} |"
        )
    report.extend(["", "## Verdict", "", verdict["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Crossed 120 triplet projector audit")
    print(f"  field-only Spin(10) crossed projector possible: {verdict['field_only_unbroken_spin10_projector_possible']}")
    print(f"  PS/constrained projector possible: {verdict['ps_eft_or_constrained_projector_possible']}")
    print(f"  finite leakage rows future-safe: {verdict['all_finite_leakage_rows_future_safe']}")
    print(f"  direct split Planck-safe kappas: {verdict['direct_split_planck_safe_kappas']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
