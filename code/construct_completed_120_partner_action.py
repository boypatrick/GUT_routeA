#!/usr/bin/env python3
"""Audit an inert 5+5bar completion of the crossed-120 source package.

No web lookup is used.

The PS crossed-120 source action showed a sharp dichotomy:

* constrained/nonpropagating source: threshold silent;
* literal triplet-only propagating mediator: nonuniversal threshold no-go.

This script checks the smallest propagating repair: add inert doublet partners
to the two triplet mediator pairs so that the finite package is two complete
5+5bar pairs at the same mass M_lock.  The audit asks:

1. how accurately the doublet partners must be mass-locked to the triplets;
2. whether an inert grading can forbid their Yukawa and physical-doublet
   mixing operators;
3. what universal alpha_G shift and SO(10) one-loop UV beta pressure this
   package creates if embedded as two 10_G copies.

The result is still conditional: it repairs one-loop nonuniversal thresholds,
but does not by itself derive the crossed projector from a microscopic
Spin(10)-invariant superpotential.
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

from audit_gtr_mediator_spectrum import B_DOUBLE_PAIR, B_FIVE_PAIR, B_TRIPLET_PAIR, projected_l2  # noqa: E402


OUT = ROOT / "output" / "completed_120_partner_action"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"

XI_GRID = [0.99, 0.995, 0.998, 0.999, 1.0, 1.001, 1.002, 1.005, 1.01]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def threshold_delta(mg: float, mt: float, md: float) -> np.ndarray:
    """Two triplet pairs at mt plus two doublet pairs at md."""

    return (
        2.0 * B_TRIPLET_PAIR * math.log(mg / mt)
        + 2.0 * B_DOUBLE_PAIR * math.log(mg / md)
    ) / (2.0 * math.pi)


def split_rows(mg: float, m_lock: float, tol: float) -> list[dict[str, Any]]:
    rows = []
    for xi in XI_GRID:
        md = xi * m_lock
        delta = threshold_delta(mg, m_lock, md)
        p_l2 = projected_l2(delta)
        rows.append(
            {
                "xi_Mdoublet_over_Mtriplet": float(xi),
                "M_triplet_GeV": float(m_lock),
                "M_doublet_GeV": float(md),
                "delta_vector": [float(x) for x in delta],
                "projected_l2": float(p_l2),
                "within_R200_budget": bool(p_l2 <= tol),
                "universal_at_exact_lock": bool(abs(xi - 1.0) < 1.0e-15 and p_l2 < 1.0e-12),
            }
        )
    return rows


def locking_window(tol: float) -> dict[str, float]:
    p_beta = projected_l2(2.0 * B_DOUBLE_PAIR)
    width = 2.0 * math.pi * tol / p_beta
    return {
        "projected_beta_l2_for_two_doublet_pairs": float(p_beta),
        "tolerance": float(tol),
        "max_abs_log_xi": float(width),
        "xi_min": float(math.exp(-width)),
        "xi_max": float(math.exp(width)),
        "percent_window": float(100.0 * (math.exp(width) - 1.0)),
    }


def grading_audit() -> dict[str, Any]:
    """Minimal post-GUT inert grading for doublet partners."""

    fields = {
        "16_i": {"Z2_inert": 0, "role": "matter"},
        "H_phys_doublets_from_120A_120B": {"Z2_inert": 0, "role": "physical doublet Yukawa sector"},
        "Phi_T_crossed_source": {"Z2_inert": 0, "role": "already-audited crossed triplet source block"},
        "L_inert_1_2": {"Z2_inert": 1, "role": "inert doublet partners completing the triplet pairs"},
        "Lbar_inert_1_2": {"Z2_inert": 1, "role": "vector partners of inert doublets"},
    }
    operators = [
        {
            "operator": "L_inert Lbar_inert",
            "charge": (fields["L_inert_1_2"]["Z2_inert"] + fields["Lbar_inert_1_2"]["Z2_inert"]) % 2,
            "desired": "allowed",
        },
        {
            "operator": "16_i 16_j L_inert",
            "charge": (2 * fields["16_i"]["Z2_inert"] + fields["L_inert_1_2"]["Z2_inert"]) % 2,
            "desired": "forbidden",
        },
        {
            "operator": "H_phys_doublet Lbar_inert",
            "charge": (
                fields["H_phys_doublets_from_120A_120B"]["Z2_inert"]
                + fields["Lbar_inert_1_2"]["Z2_inert"]
            )
            % 2,
            "desired": "forbidden",
        },
        {
            "operator": "Phi_T_crossed_source mass block",
            "charge": 0,
            "desired": "allowed",
        },
    ]
    for op in operators:
        op["actual"] = "allowed" if op["charge"] == 0 else "forbidden"
        op["passes"] = op["actual"] == op["desired"]
    return {
        "group": "Z2_inert after PS breaking",
        "fields": fields,
        "operators": operators,
        "all_selection_rules_pass": all(op["passes"] for op in operators),
        "interpretation": (
            "The inert grading can decouple the added doublet partners from matter Yukawa "
            "and physical doublet mixing while allowing their vectorlike masses.  It does "
            "not derive the crossed triplet block; it only repairs the threshold bookkeeping."
        ),
    }


def landau_ratio(alpha_inv: float, sum_t: float) -> float:
    b10 = sum_t - 24.0
    if b10 <= 0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def uv_rows(alpha_inv: float) -> list[dict[str, Any]]:
    scenarios = [
        {
            "scenario": "two_10G_completed_partner_package_only",
            "sum_T": 2.0,
            "interpretation": "two SO(10) 10_G copies carrying the two 5+5bar partner packages",
        },
        {
            "scenario": "minimal_yukawa_sector_plus_two_10G",
            "sum_T": 72.0,
            "interpretation": "three 16 matter fields plus 10_H,120_H,126bar_H and two 10_G copies",
        },
    ]
    rows = []
    for row in scenarios:
        b10 = row["sum_T"] - 24.0
        ratio = landau_ratio(alpha_inv, row["sum_T"])
        rows.append(
            {
                "scenario": row["scenario"],
                "sum_T": row["sum_T"],
                "b10": b10,
                "alphaG_inv_MG": float(alpha_inv),
                "landau_ratio_from_MG": float(ratio) if math.isfinite(ratio) else math.inf,
                "passes_R200_if_applicable": bool(ratio > 200.0 if math.isfinite(ratio) else True),
                "interpretation": row["interpretation"],
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    ledger = read_json(LEDGER)
    mg = float(ps["benchmarks"]["MG_GeV"])
    m_lock = float(ps["benchmarks"]["M_lock_GeV"])
    tol = float(ps["benchmarks"]["R200_projected_threshold_budget"])
    alpha_inv = float(ledger["r200_benchmark"]["alphaG_inv"])
    rows = split_rows(mg, m_lock, tol)
    window = locking_window(tol)
    grading = grading_audit()
    uv = uv_rows(alpha_inv)
    locked = next(row for row in rows if abs(row["xi_Mdoublet_over_Mtriplet"] - 1.0) < 1.0e-15)
    near = [row for row in rows if row["within_R200_budget"]]
    universal_shift = 2.0 * math.log(mg / m_lock) / (2.0 * math.pi)
    verdict = {
        "exact_completed_pair_threshold_silent": bool(locked["projected_l2"] < 1.0e-12),
        "inert_grading_decouples_doublets": grading["all_selection_rules_pass"],
        "mass_locking_required": True,
        "audited_split_rows_within_budget": len(near),
        "max_abs_log_xi_for_R200_budget": window["max_abs_log_xi"],
        "percent_mass_lock_window": window["percent_window"],
        "two_10G_package_uv_safe_by_itself": uv[0]["passes_R200_if_applicable"],
        "minimal_yukawa_plus_two_10G_R200_safe": uv[1]["passes_R200_if_applicable"],
        "interpretation": (
            "Adding inert doublet partners repairs the crossed-120 literal mediator threshold "
            "only if the triplet and doublet partner masses are locked at the per-mille level "
            "for the R=200 benchmark.  A Z2_inert can decouple the added doublets from matter "
            "Yukawa operators.  Embedded as two 10_G copies, the package itself is UV-benign, "
            "but this remains a completed post-GUT mediator repair rather than a microscopic "
            "derivation of the crossed projector."
        ),
    }
    return {
        "note": "No web lookup used. Inert completed-partner audit for the crossed 120 source package.",
        "input_files": {
            "ps_crossed_120_source_action": str(PS_SOURCE),
            "conditional_theorem_ledger": str(LEDGER),
        },
        "superpotential_completion": {
            "formula": "W_comp = W_T[crossed source] + M_L sum_a L_a Lbar_a",
            "exact_lock_condition": "M_L = M_T = M_lock for each inert 5+5bar package",
            "purpose": "complete two triplet mediator pairs into two threshold-silent 5+5bar pairs",
        },
        "benchmarks": {
            "MG_GeV": mg,
            "M_lock_GeV": m_lock,
            "alphaG_inv_R200": alpha_inv,
            "R200_projected_threshold_budget": tol,
            "universal_alpha_inverse_shift_at_exact_lock": float(universal_shift),
        },
        "split_scan": rows,
        "mass_locking_window": window,
        "inert_grading": grading,
        "uv_beta_scan": uv,
        "verdict": verdict,
    }


def write_report(payload: dict[str, Any]) -> None:
    lines = [
        "# Completed inert 5+5bar partner audit",
        "",
        "No web lookup was used.",
        "",
        "## Result",
        "",
        payload["verdict"]["interpretation"],
        "",
        "## Mass-locking window",
        "",
        f"R=200 budget: {payload['benchmarks']['R200_projected_threshold_budget']:.6e}",
        f"max |log xi|: {payload['mass_locking_window']['max_abs_log_xi']:.6e}",
        f"xi window: [{payload['mass_locking_window']['xi_min']:.9f}, {payload['mass_locking_window']['xi_max']:.9f}]",
        "",
        "| xi=M_D/M_T | projected l2 | within budget |",
        "|---:|---:|---:|",
    ]
    for row in payload["split_scan"]:
        lines.append(
            f"| {row['xi_Mdoublet_over_Mtriplet']:.6f} | {row['projected_l2']:.6e} | "
            f"{row['within_R200_budget']} |"
        )
    lines.extend(["", "## Inert grading", "", "| operator | actual | desired | pass |", "|---|---|---|---|"])
    for row in payload["inert_grading"]["operators"]:
        lines.append(f"| `{row['operator']}` | {row['actual']} | {row['desired']} | {row['passes']} |")
    lines.extend(["", "## UV beta scan", "", "| scenario | b10 | Landau ratio | R=200 safe |", "|---|---:|---:|---:|"])
    for row in payload["uv_beta_scan"]:
        ratio = row["landau_ratio_from_MG"]
        ratio_text = "inf" if math.isinf(ratio) else f"{ratio:.6e}"
        lines.append(f"| {row['scenario']} | {row['b10']:.1f} | {ratio_text} | {row['passes_R200_if_applicable']} |")
    lines.append("")
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "split_threshold_scan.csv",
        payload["split_scan"],
        [
            "xi_Mdoublet_over_Mtriplet",
            "M_triplet_GeV",
            "M_doublet_GeV",
            "projected_l2",
            "within_R200_budget",
            "universal_at_exact_lock",
            "delta_vector",
        ],
    )
    write_csv(
        OUT / "uv_beta_scan.csv",
        payload["uv_beta_scan"],
        [
            "scenario",
            "sum_T",
            "b10",
            "alphaG_inv_MG",
            "landau_ratio_from_MG",
            "passes_R200_if_applicable",
            "interpretation",
        ],
    )
    write_report(payload)
    verdict = payload["verdict"]
    print("Completed inert 5+5bar partner audit")
    print(f"  exact completed pair threshold silent: {verdict['exact_completed_pair_threshold_silent']}")
    print(f"  inert grading decouples doublets: {verdict['inert_grading_decouples_doublets']}")
    print(f"  percent mass-lock window: {verdict['percent_mass_lock_window']:.4f}%")
    print(f"  minimal Yukawa + two 10G R200 safe: {verdict['minimal_yukawa_plus_two_10G_R200_safe']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
