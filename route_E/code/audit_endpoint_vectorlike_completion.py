#!/usr/bin/env python3
"""Endpoint vectorlike/anomaly completion audit for the hidden unitary link.

No web lookup is used.

The hidden radial-lock sector uses endpoint/source copy groups

    U(4)_L x U(4)_H x U(4)_R

with bifundamentals Q and Qtilde.  If the endpoint copy groups are promoted
from global/source symmetries to propagating hidden gauge groups, Q and Qtilde
alone are chiral under U(4)_L and U(4)_R.  This script checks the minimal
vectorlike completion

    Qc       ~ (anti-L, fund-H, 1),
    Qtilde_c ~ (1, anti-H, fund-R),

and computes cubic anomaly proxies, one-loop hidden beta coefficients, and the
visible threshold vector.  All fields are visible Spin(10) singlets in this
branch, so a successful endpoint completion remains threshold-silent for the
GUT matching problem.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "endpoint_vectorlike_completion"
RADIAL = ROOT / "output" / "hidden_radial_lock_sector" / "summary.json"

N = 4
TARGET_R = 200.0
ALPHA_INV_GRID = [1.0, 5.0, 10.0, 25.0]
GROUPS = ["SU4_L", "SU4_H", "SU4_R"]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def field_content(include_conjugates: bool) -> list[dict[str, Any]]:
    fields = [
        {
            "field": "Q",
            "SU4_L": "fund",
            "SU4_H": "anti",
            "SU4_R": "singlet",
            "visible_Spin10": "singlet",
            "role": "left radial link",
        },
        {
            "field": "Qtilde",
            "SU4_L": "singlet",
            "SU4_H": "fund",
            "SU4_R": "anti",
            "visible_Spin10": "singlet",
            "role": "right radial link",
        },
    ]
    if include_conjugates:
        fields.extend(
            [
                {
                    "field": "Qc",
                    "SU4_L": "anti",
                    "SU4_H": "fund",
                    "SU4_R": "singlet",
                    "visible_Spin10": "singlet",
                    "role": "vectorlike partner of Q",
                },
                {
                    "field": "Qtilde_c",
                    "SU4_L": "singlet",
                    "SU4_H": "anti",
                    "SU4_R": "fund",
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
    # d alpha^{-1}/d log mu = -b/(2 pi)
    if b <= 0.0:
        return None
    return math.exp(2.0 * math.pi * alpha_inv / b)


def beta_rows(fields: list[dict[str, Any]], label: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group in GROUPS:
        sum_t = sum(dynkin_sum(field, group) for field in fields)
        b = sum_t - 3.0 * N
        for alpha_inv in ALPHA_INV_GRID:
            ratio = landau_ratio(alpha_inv, b)
            rows.append(
                {
                    "scenario": label,
                    "group": group,
                    "sum_T": sum_t,
                    "b_hidden": b,
                    "alpha_inv_at_matching": alpha_inv,
                    "uv_landau_ratio_over_matching": ratio,
                    "safe_to_R200": True if ratio is None else ratio > TARGET_R,
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


def build() -> dict[str, Any]:
    radial = read_json(RADIAL)
    raw = field_content(include_conjugates=False)
    completed = field_content(include_conjugates=True)
    anomaly = anomaly_rows(raw, "raw_Q_Qtilde") + anomaly_rows(completed, "vectorlike_completed")
    beta = beta_rows(raw, "raw_Q_Qtilde") + beta_rows(completed, "vectorlike_completed")
    threshold = visible_threshold_rows(raw, "raw_Q_Qtilde") + visible_threshold_rows(
        completed, "vectorlike_completed"
    )

    raw_anomaly_free = all(row["anomaly_cancels"] for row in anomaly if row["scenario"] == "raw_Q_Qtilde")
    completed_anomaly_free = all(
        row["anomaly_cancels"] for row in anomaly if row["scenario"] == "vectorlike_completed"
    )
    completed_beta_safe = all(
        row["safe_to_R200"] for row in beta if row["scenario"] == "vectorlike_completed"
    )
    completed_threshold = [
        sum(row[f"delta_b{i}"] for row in threshold if row["scenario"] == "vectorlike_completed")
        for i in [1, 2, 3]
    ]

    verdict = {
        "raw_endpoint_gauging_is_anomalous": not raw_anomaly_free,
        "vectorlike_completion_cancels_all_hidden_cubic_anomalies": completed_anomaly_free,
        "vectorlike_completion_beta_safe_to_R200": completed_beta_safe,
        "vectorlike_completion_visible_threshold_vector": completed_threshold,
        "radial_lock_geometry_retained": radial["verdict"]["all_finite_samples_Dflat_unitary_and_locked"],
        "residual_unitary_completion_moduli": radial["verdict"]["residual_after_fixed_block"],
        "interpretation": (
            "Gauging the endpoint copy groups with only Q and Qtilde is chiral: "
            "SU(4)_L and SU(4)_R have nonzero cubic anomaly proxies.  Adding "
            "the conjugate hidden bifundamentals Qc and Qtilde_c makes every "
            "hidden copy gauge factor vectorlike.  The one-loop coefficients "
            "are b_L=-8, b_H=-4, b_R=-8, so there is no UV Landau obstruction "
            "to R=200 at one loop.  Since every endpoint-completion field is a "
            "visible Spin(10) singlet, the direct visible threshold vector is "
            "zero.  This closes the endpoint anomaly/beta bookkeeping "
            "conditionally, while preserving the previously audited radial "
            "unitary-link geometry."
        ),
    }

    return {
        "note": "No web lookup used. Endpoint vectorlike/anomaly completion audit.",
        "field_content": {
            "raw_Q_Qtilde": raw,
            "vectorlike_completed": completed,
        },
        "anomaly_audit": anomaly,
        "beta_audit": beta,
        "visible_threshold_audit": threshold,
        "radial_lock_input": {
            "source": "output/hidden_radial_lock_sector/summary.json",
            "residual_after_fixed_block": radial["verdict"]["residual_after_fixed_block"],
            "max_D_residual_2norm": radial["verdict"]["max_D_residual_2norm"],
            "max_abs_log_singular": radial["verdict"]["max_abs_log_singular"],
        },
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    v = payload["verdict"]
    beta_completed = [
        row
        for row in payload["beta_audit"]
        if row["scenario"] == "vectorlike_completed" and row["alpha_inv_at_matching"] == 10.0
    ]
    lines = [
        "# Endpoint vectorlike completion audit",
        "",
        "No web lookup was used.",
        "",
        "## Anomaly result",
        "",
        f"raw endpoint gauging anomalous: {v['raw_endpoint_gauging_is_anomalous']}",
        f"vectorlike completion cancels anomalies: {v['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}",
        "",
        "## One-loop hidden beta coefficients",
        "",
        "| group | sum T | b | safe to R=200 |",
        "|---|---:|---:|---:|",
    ]
    for row in beta_completed:
        lines.append(
            f"| {row['group']} | {row['sum_T']:.1f} | {row['b_hidden']:.1f} | {row['safe_to_R200']} |"
        )
    lines.extend(
        [
            "",
            "## Visible threshold",
            "",
            f"visible threshold vector: {v['vectorlike_completion_visible_threshold_vector']}",
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
        ["field", "SU4_L", "SU4_H", "SU4_R", "visible_Spin10", "role"],
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
    print("Endpoint vectorlike completion audit")
    print(f"  raw endpoint gauging anomalous: {v['raw_endpoint_gauging_is_anomalous']}")
    print(f"  vectorlike anomalies cancel: {v['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}")
    print(f"  beta safe to R=200: {v['vectorlike_completion_beta_safe_to_R200']}")
    print(f"  visible threshold vector: {v['vectorlike_completion_visible_threshold_vector']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
