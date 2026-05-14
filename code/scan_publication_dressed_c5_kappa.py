#!/usr/bin/env python3
"""Scan the crossed-triplet finite-block kappa in the dressed C5 replay.

No web lookup is used.

The publication eigenstate card stores two tensor endpoints:

    C5_rank_one  = C5[diag(1,0,0,0)],
    C5_finite30  = C5[diag(1,1/30,1/30,1/30)].

Because the contraction is linear in the inverse triplet block, the tensor for
an arbitrary common orthogonal suppression kappa is

    C5(kappa) = C5_rank_one + (30/kappa) (C5_finite30 - C5_rank_one).

This scan reuses the dressed-channel machinery and asks how close to rank-one
the source block must be before the reference filter S_T=1e-5 survives the
full soft/chiral stress grid.
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

import audit_eigenstate_d5_dressing as eig  # noqa: E402
import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import audit_mssm_mixing_d5_dressing as mssm  # noqa: E402
import audit_publication_dressed_c5_from_eigenstate_card as dressed  # noqa: E402


OUT = ROOT / "output" / "publication_dressed_c5_kappa_scan"
KAPPA_GRID = [30.0, 50.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0]
REFERENCE_KAPPA = 30.0


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def tensor4_json(tensor: np.ndarray) -> list[list[list[list[dict[str, float]]]]]:
    return [
        [
            [
                [cjson(tensor[a, b, c, d]) for d in range(3)]
                for c in range(3)
            ]
            for b in range(3)
        ]
        for a in range(3)
    ]


def tensor_for_kappa(rank: np.ndarray, finite30: np.ndarray, kappa: float) -> np.ndarray:
    return rank + (REFERENCE_KAPPA / kappa) * (finite30 - rank)


def card_for_kappa(base: dict[str, Any], kappa: float) -> dict[str, Any]:
    card = dict(base)
    for op in ["L", "R"]:
        rank = dressed.tensor4(base[f"C5{op}_rank_one"])
        finite30 = dressed.tensor4(base[f"C5{op}_finite"])
        card[f"C5{op}_finite"] = tensor4_json(tensor_for_kappa(rank, finite30, kappa))
    card["kappa"] = float(kappa)
    return card


def rows_for_kappa(base_card: dict[str, Any], kappa: float) -> list[dict[str, Any]]:
    card = card_for_kappa(base_card, kappa)
    vac = mi.vacuum_inputs()
    soft = mi.insertion_basis()
    rows: list[dict[str, Any]] = []
    block = "finite"
    for scenario, base_delta in soft["matrices"].items():
        for eps in eig.SOFT_EPS_GRID:
            spec = eig.spectrum_for_scenario(base_delta, eps)
            if not spec["positive_definite"]:
                continue
            for point in mi.SPECTRUM_POINTS:
                ewkino = mssm.electroweakino_spectrum(point)
                for row in dressed.dressed_channel_rows_for_point(
                    card, block, spec, point, vac, ewkino, scenario, eps
                ):
                    row["kappa"] = float(kappa)
                    rows.append(row)
    return rows


def summarize_kappa(rows: list[dict[str, Any]], kappa: float) -> dict[str, Any]:
    out: dict[str, Any] = {"kappa": float(kappa)}
    for case in dressed.WIDTH_CASES:
        case_rows = [row for row in rows if row["normalization_case"] == case]
        worst_1e35 = min(case_rows, key=lambda row: row["margin_1e35_at_ST_display"])
        worst_present = min(case_rows, key=lambda row: row["margin_present_at_ST_display"])
        out[f"{case}_rows"] = len(case_rows)
        out[f"{case}_unsafe_1e35_rows"] = sum(not row["passes_1e35_at_ST_display"] for row in case_rows)
        out[f"{case}_unsafe_present_rows"] = sum(not row["passes_present_at_ST_display"] for row in case_rows)
        out[f"{case}_global_STmax_1e35"] = min(row["S_T_max_1e35"] for row in case_rows)
        out[f"{case}_worst_margin_1e35"] = worst_1e35["margin_1e35_at_ST_display"]
        out[f"{case}_worst_channel"] = worst_1e35["channel"]
        out[f"{case}_worst_operator"] = worst_1e35["operator"]
        out[f"{case}_worst_scenario"] = worst_1e35["scenario"]
        out[f"{case}_worst_epsilon"] = worst_1e35["epsilon"]
        out[f"{case}_worst_spectrum"] = worst_1e35["spectrum_name"]
        out[f"{case}_worst_pair"] = worst_1e35["pair"]
        out[f"{case}_worst_present_margin"] = worst_present["margin_present_at_ST_display"]
    out["passes_central_1e35_at_ST_1e_minus_5"] = out["central_unsafe_1e35_rows"] == 0
    out["passes_max_width_1e35_at_ST_1e_minus_5"] = out["max_width_unsafe_1e35_rows"] == 0
    return out


def audit() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base = dressed.read_json(dressed.CARD)
    scan_rows: list[dict[str, Any]] = []
    for kappa in KAPPA_GRID:
        rows = rows_for_kappa(base, kappa)
        scan_rows.append(summarize_kappa(rows, kappa))

    central_pass = [row for row in scan_rows if row["passes_central_1e35_at_ST_1e_minus_5"]]
    max_pass = [row for row in scan_rows if row["passes_max_width_1e35_at_ST_1e_minus_5"]]
    summary = {
        "note": "No web lookup used. Kappa scan for dressed C5 replay using linear interpolation between rank-one and finite-kappa=30 tensors.",
        "reference_filter": dressed.DISPLAY_ST,
        "reference_kappa": REFERENCE_KAPPA,
        "kappa_grid": KAPPA_GRID,
        "scan_rows": scan_rows,
        "verdict": {
            "central_min_kappa_pass_1e35": min((row["kappa"] for row in central_pass), default=None),
            "max_width_min_kappa_pass_1e35": min((row["kappa"] for row in max_pass), default=None),
            "kappa_30_passes_central": scan_rows[0]["passes_central_1e35_at_ST_1e_minus_5"],
            "kappa_30_passes_max_width": scan_rows[0]["passes_max_width_1e35_at_ST_1e_minus_5"],
            "interpretation": (
                "If a sufficiently large kappa passes, the d=5 failure can be recast as a rank-one-locking "
                "requirement on the crossed triplet inverse block.  If no scanned kappa passes max-width, "
                "a soft-alignment condition or tighter triplet filter remains necessary."
            ),
        },
    }
    return scan_rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Publication dressed C5 kappa scan",
        "",
        "No web lookup was used.",
        "",
        "The scan uses the linear identity",
        "",
        "`C5(kappa) = C5_rank_one + (30/kappa) (C5_finite30 - C5_rank_one)`.",
        "",
        "| kappa | central unsafe 1e35 | central S_T max | max-width unsafe 1e35 | max-width S_T max |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in summary["scan_rows"]:
        lines.append(
            f"| {row['kappa']:.0f} | {row['central_unsafe_1e35_rows']} | "
            f"{row['central_global_STmax_1e35']:.6e} | {row['max_width_unsafe_1e35_rows']} | "
            f"{row['max_width_global_STmax_1e35']:.6e} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        summary["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = audit()
    write_csv(OUT / "kappa_scan.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    print("Publication dressed C5 kappa scan")
    print(f"  central min kappa pass: {summary['verdict']['central_min_kappa_pass_1e35']}")
    print(f"  max-width min kappa pass: {summary['verdict']['max_width_min_kappa_pass_1e35']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
