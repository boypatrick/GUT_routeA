#!/usr/bin/env python3
"""Vectorlike endpoint-gauge completion for the 8x8 clockwork meson lock.

No web lookup is used.

The previous endpoint completion was the older 4x4 crossed-link bookkeeping.
The current clockwork quotient uses

    U(8)_L x U(8)_H x U(8)_R

endpoint/copy groups with bifundamentals Q and Qtilde.  If these copy groups
are made propagating hidden gauge groups, Q and Qtilde alone are chiral under
the endpoint factors.  This audit checks the minimal vectorlike completion,
its one-loop hidden beta coefficients, and the direct visible threshold vector.
All endpoint fields are visible Spin(10) singlets by construction.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clockwork_endpoint_vectorlike_completion"
HIDDEN_MESON = ROOT / "output" / "clockwork_hidden_endpoint_meson" / "summary.json"

N = 8
TARGET_R = 200.0
ALPHA_INV_GRID = [1.0, 5.0, 10.0, 25.0]
GROUPS = ["SU8_L", "SU8_H", "SU8_R"]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def field_content(include_conjugates: bool) -> list[dict[str, Any]]:
    fields = [
        {
            "field": "Q",
            "SU8_L": "fund",
            "SU8_H": "anti",
            "SU8_R": "singlet",
            "visible_Spin10": "singlet",
            "role": "left endpoint/radial link",
        },
        {
            "field": "Qtilde",
            "SU8_L": "singlet",
            "SU8_H": "fund",
            "SU8_R": "anti",
            "visible_Spin10": "singlet",
            "role": "right endpoint/radial link",
        },
    ]
    if include_conjugates:
        fields.extend(
            [
                {
                    "field": "Qc",
                    "SU8_L": "anti",
                    "SU8_H": "fund",
                    "SU8_R": "singlet",
                    "visible_Spin10": "singlet",
                    "role": "vectorlike partner of Q",
                },
                {
                    "field": "Qtilde_c",
                    "SU8_L": "singlet",
                    "SU8_H": "anti",
                    "SU8_R": "fund",
                    "visible_Spin10": "singlet",
                    "role": "vectorlike partner of Qtilde",
                },
            ]
        )
    return fields


def rep_sign(rep: str) -> int:
    if rep == "fund":
        return 1
    if rep == "anti":
        return -1
    return 0


def multiplicity(field: dict[str, Any], group: str) -> int:
    mult = 1
    for other in GROUPS:
        if other == group:
            continue
        if field[other] in {"fund", "anti"}:
            mult *= N
    return mult


def anomaly_rows(fields: list[dict[str, Any]], label: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group in GROUPS:
        cubic = 0
        abs_chiral_index = 0
        for field in fields:
            sign = rep_sign(field[group])
            mult = multiplicity(field, group)
            cubic += sign * mult
            abs_chiral_index += abs(sign) * mult
        rows.append(
            {
                "scenario": label,
                "group": group,
                "cubic_anomaly_proxy": cubic,
                "absolute_chiral_index": abs_chiral_index,
                "anomaly_cancels": cubic == 0,
            }
        )
    return rows


def dynkin_sum(field: dict[str, Any], group: str) -> float:
    if field[group] not in {"fund", "anti"}:
        return 0.0
    return 0.5 * multiplicity(field, group)


def landau_ratio(alpha_inv: float, b: float) -> float | None:
    # d alpha^{-1}/d log mu = -b/(2 pi).
    if b <= 0.0:
        return None
    return math.exp(2.0 * math.pi * alpha_inv / b)


def beta_rows(fields: list[dict[str, Any]], label: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group in GROUPS:
        sum_t = sum(dynkin_sum(field, group) for field in fields)
        b = sum_t - 3.0 * N
        equivalent_nf_pairs = sum_t
        if equivalent_nf_pairs < 1.5 * N:
            phase = "below_conformal_window_or_quantum_deformed_edge"
        elif equivalent_nf_pairs < 3.0 * N:
            phase = "inside_conformal_window"
        elif equivalent_nf_pairs == 3.0 * N:
            phase = "one_loop_marginal_edge"
        else:
            phase = "IR_free_UV_Landau_possible"
        for alpha_inv in ALPHA_INV_GRID:
            ratio = landau_ratio(alpha_inv, b)
            rows.append(
                {
                    "scenario": label,
                    "group": group,
                    "sum_T": sum_t,
                    "equivalent_Nf_pairs": equivalent_nf_pairs,
                    "phase_label": phase,
                    "b_hidden": b,
                    "alpha_inv_at_matching": alpha_inv,
                    "uv_landau_ratio_over_matching": ratio,
                    "safe_to_R200": bool(ratio is None or ratio > TARGET_R),
                }
            )
    return rows


def visible_threshold_rows(fields: list[dict[str, Any]], label: str) -> list[dict[str, Any]]:
    return [
        {
            "scenario": label,
            "field": field["field"],
            "visible_Spin10": field["visible_Spin10"],
            "delta_b1": 0.0,
            "delta_b2": 0.0,
            "delta_b3": 0.0,
        }
        for field in fields
    ]


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def build() -> dict[str, Any]:
    hidden = read_json(HIDDEN_MESON)
    raw = field_content(include_conjugates=False)
    completed = field_content(include_conjugates=True)
    anomaly = anomaly_rows(raw, "raw_Q_Qtilde") + anomaly_rows(completed, "vectorlike_completed")
    beta = beta_rows(raw, "raw_Q_Qtilde") + beta_rows(completed, "vectorlike_completed")
    threshold = visible_threshold_rows(raw, "raw_Q_Qtilde") + visible_threshold_rows(completed, "vectorlike_completed")

    raw_free = all(row["anomaly_cancels"] for row in anomaly if row["scenario"] == "raw_Q_Qtilde")
    completed_free = all(row["anomaly_cancels"] for row in anomaly if row["scenario"] == "vectorlike_completed")
    completed_beta_safe = all(row["safe_to_R200"] for row in beta if row["scenario"] == "vectorlike_completed")
    completed_threshold = [
        sum(row[f"delta_b{i}"] for row in threshold if row["scenario"] == "vectorlike_completed")
        for i in [1, 2, 3]
    ]
    beta_alpha10 = [
        row for row in beta if row["scenario"] == "vectorlike_completed" and row["alpha_inv_at_matching"] == 10.0
    ]
    h_row = next(row for row in beta_alpha10 if row["group"] == "SU8_H")

    verdict = {
        "raw_endpoint_gauging_is_anomalous": not raw_free,
        "vectorlike_completion_cancels_all_hidden_cubic_anomalies": completed_free,
        "vectorlike_completion_beta_safe_to_R200": completed_beta_safe,
        "vectorlike_completion_visible_threshold_vector": completed_threshold,
        "hidden_mesonic_radial_lock_retained": hidden["verdict"]["finite_samples_Dflat_unitary_quantum_and_locked"],
        "residual_unitary_completion_moduli": hidden["rank_audit"]["fixed_block"]["residual_unitary_completion_moduli"],
        "SU8H_vectorlike_phase_at_alpha10": h_row["phase_label"],
        "interpretation": (
            "The 8x8 endpoint copy gauging is chiral if only Q and Qtilde are "
            "kept: SU(8)_L and SU(8)_R have nonzero cubic anomaly proxies.  "
            "Adding Qc and Qtilde_c makes every hidden copy factor vectorlike.  "
            "At one loop the vectorlike-completed beta coefficients are "
            "b_L=-16, b_H=-8, b_R=-16, so there is no UV Landau pole before "
            "R=200 for the tested matching couplings.  The SU(8)_H factor has "
            "equivalent Nf=16, inside the SQCD conformal-window range in this "
            "bookkeeping sense, while the endpoint SU(8)_L/R factors sit at "
            "the Nf=Nc edge.  All fields are visible Spin(10) singlets, so the "
            "direct visible threshold vector remains zero."
        ),
    }
    return {
        "note": "No web lookup used. 8x8 endpoint vectorlike/anomaly completion audit.",
        "field_content": {
            "raw_Q_Qtilde": raw,
            "vectorlike_completed": completed,
        },
        "anomaly_audit": anomaly,
        "beta_audit": beta,
        "visible_threshold_audit": threshold,
        "hidden_meson_input": {
            "source": "output/clockwork_hidden_endpoint_meson/summary.json",
            "residual_unitary_completion_moduli": hidden["rank_audit"]["fixed_block"]["residual_unitary_completion_moduli"],
            "max_D_residual_2norm": hidden["verdict"]["max_D_residual_2norm"],
            "max_abs_log_singular": hidden["verdict"]["max_abs_log_singular"],
        },
        "verdict": verdict,
    }


def write_report(payload: dict[str, Any]) -> None:
    v = payload["verdict"]
    beta_completed = [
        row
        for row in payload["beta_audit"]
        if row["scenario"] == "vectorlike_completed" and row["alpha_inv_at_matching"] == 10.0
    ]
    lines = [
        "# 8x8 endpoint vectorlike completion audit",
        "",
        "No web lookup was used.",
        "",
        "## Anomaly result",
        "",
        f"raw endpoint gauging anomalous: `{v['raw_endpoint_gauging_is_anomalous']}`",
        f"vectorlike completion cancels anomalies: `{v['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}`",
        "",
        "## One-loop hidden beta coefficients",
        "",
        "| group | sum T | b | phase | safe to R=200 |",
        "|---|---:|---:|---|---:|",
    ]
    for row in beta_completed:
        lines.append(
            f"| {row['group']} | {row['sum_T']:.1f} | {row['b_hidden']:.1f} | "
            f"{row['phase_label']} | {row['safe_to_R200']} |"
        )
    lines.extend(
        [
            "",
            "## Visible threshold",
            "",
            f"visible threshold vector: `{v['vectorlike_completion_visible_threshold_vector']}`",
            "",
            "## Verdict",
            "",
            v["interpretation"],
            "",
        ]
    )
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    all_fields = payload["field_content"]["raw_Q_Qtilde"] + payload["field_content"]["vectorlike_completed"]
    write_csv(
        OUT / "field_content.csv",
        all_fields,
        ["field", "SU8_L", "SU8_H", "SU8_R", "visible_Spin10", "role"],
    )
    write_csv(
        OUT / "anomaly_audit.csv",
        payload["anomaly_audit"],
        ["scenario", "group", "cubic_anomaly_proxy", "absolute_chiral_index", "anomaly_cancels"],
    )
    write_csv(
        OUT / "beta_audit.csv",
        payload["beta_audit"],
        [
            "scenario",
            "group",
            "sum_T",
            "equivalent_Nf_pairs",
            "phase_label",
            "b_hidden",
            "alpha_inv_at_matching",
            "uv_landau_ratio_over_matching",
            "safe_to_R200",
        ],
    )
    write_csv(
        OUT / "visible_threshold_audit.csv",
        payload["visible_threshold_audit"],
        ["scenario", "field", "visible_Spin10", "delta_b1", "delta_b2", "delta_b3"],
    )
    write_report(payload)
    v = payload["verdict"]
    print("8x8 endpoint vectorlike completion audit")
    print(f"  raw endpoint gauging anomalous: {v['raw_endpoint_gauging_is_anomalous']}")
    print(f"  vectorlike anomalies cancel: {v['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}")
    print(f"  beta safe to R=200: {v['vectorlike_completion_beta_safe_to_R200']}")
    print(f"  visible threshold vector: {v['vectorlike_completion_visible_threshold_vector']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
