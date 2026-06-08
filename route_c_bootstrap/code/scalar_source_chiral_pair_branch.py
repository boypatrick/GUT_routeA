#!/usr/bin/env python3
"""Route-C P12-S scalar/source chiral-pair branch.

P12 records both possible post-P11 Wilson-basis directions, then activates the
scalar/source chiral-pair branch first.  This branch takes the P7 all-incoming
chiral-pair pole rows literally as scalar/source exchanges:

  L ⊃ lambda_{ab,R} psi_a psi_b S_R + h.c. + M_{S_R}^2 |S_R|^2

Integrating out S_R gives

  C6(12|34;R) = lambda_{12,R} lambda^*_{34,R}/M_{S_R}^2.

This is not an adjoint-vector branch and does not compute physical proton
lifetimes.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P7_JSON = OUT / "low_energy_matching_proton_report.json"
P11_JSON = OUT / "spin10_mass_replay_matching_gate.json"
OUT_JSON = OUT / "scalar_source_chiral_pair_branch.json"
OUT_MD = OUT / "scalar_source_chiral_pair_branch.md"


MED_RE = re.compile(
    r"X\[\((?P<su3>[^,]+),(?P<su2>[^\)]+)\);Y=(?P<Y>[^;]+);B-L=(?P<BL>[^\]]+)\]"
)


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_mediator(label: str) -> dict[str, str]:
    match = MED_RE.fullmatch(label)
    if not match:
        raise ValueError(f"cannot parse mediator label: {label}")
    return match.groupdict()


def sanitize(text: str) -> str:
    return (
        text.replace("^", "")
        .replace("bar", "b")
        .replace("/", "o")
        .replace("-", "m")
        .replace("+", "p")
        .replace("=", "")
        .replace(";", "_")
        .replace(",", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
        .replace(" ", "")
    )


def mass_symbol(mediator: str) -> str:
    q = parse_mediator(mediator)
    return f"M_S_{sanitize(q['su3'])}_{sanitize(q['su2'])}_Y{sanitize(q['Y'])}_BL{sanitize(q['BL'])}"


def coupling_symbol(pair: list[str], mediator: str) -> str:
    pair_name = "".join(pair).replace("^", "")
    return f"lambda_{pair_name}_{mass_symbol(mediator).removeprefix('M_S_')}"


def source_sector_rows(p7_rows: list[dict]) -> list[dict]:
    grouped: dict[str, dict] = {}
    for row in p7_rows:
        med = row["mediator"]
        q = parse_mediator(med)
        entry = grouped.setdefault(
            med,
            {
                "mediator": med,
                "source_field": "S_" + mass_symbol(med).removeprefix("M_S_"),
                "mass_symbol": mass_symbol(med),
                "quantum_numbers": q,
                "operators_supported": set(),
                "row_count": 0,
                "operator_status_counts": Counter(),
            },
        )
        entry["operators_supported"].add(row["canonical_operator"])
        entry["row_count"] += 1
        entry["operator_status_counts"][row["operator_status"]] += 1

    out = []
    for row in grouped.values():
        out.append(
            {
                "mediator": row["mediator"],
                "source_field": row["source_field"],
                "mass_symbol": row["mass_symbol"],
                "quantum_numbers": row["quantum_numbers"],
                "operators_supported": sorted(row["operators_supported"]),
                "row_count": row["row_count"],
                "operator_status_counts": dict(row["operator_status_counts"]),
                "positivity_condition": f"{row['mass_symbol']}^2 > 0",
                "branch_status": "activated_symbolically_no_numeric_mass",
            }
        )
    return sorted(out, key=lambda x: x["mass_symbol"])


def scalar_row(row: dict) -> dict:
    med = row["mediator"]
    left = row["left_pair"]
    right = row["right_pair"]
    m = mass_symbol(med)
    lam_left = coupling_symbol(left, med)
    lam_right = coupling_symbol(right, med)
    coeff = f"C6_S[{row['canonical_operator']}; {med}] = {lam_left} {lam_right}^*/{m}^2"
    return {
        "canonical_operator": row["canonical_operator"],
        "operator_status": row["operator_status"],
        "mediator": med,
        "source_field": "S_" + m.removeprefix("M_S_"),
        "mass_symbol": m,
        "left_pair": left,
        "right_pair": right,
        "scalar_source_matching_coefficient": coeff,
        "p6_longitudinal_vector_problem": False,
        "p6_scalar_auxiliary_status": "longitudinal_vector_growth_absent",
        "physical_bound_evaluable_now": False,
        "missing_inputs_before_physical_bound": [
            "numeric scalar/source masses or scan ranges",
            "coupling flavor tensors lambda_{ab,R}",
            "physical flavor rotations",
            "operator basis and Fierz/chiral contractions",
            "RG evolution",
            "hadronic matrix elements",
            "experimental proton lifetime limit for selected channel",
        ],
    }


def summarize(rows: list[dict], sectors: list[dict], p11: dict) -> dict:
    status_counts = Counter(row["operator_status"] for row in rows)
    bnv = [row for row in rows if row["operator_status"] != "B_and_L_conserving"]
    standard_bnv = [row for row in rows if row["operator_status"] == "standard_BNV_seed"]
    sterile_bnv = [row for row in rows if row["operator_status"] == "BNV_with_sterile_neutrino"]
    return {
        "p12_branch_choice": "scalar_source_chiral_pair_branch",
        "other_required_future_branch": "completed_current_current_vector_branch",
        "p11_p7_rows_remaining_scalar_source_debt": p11["summary"][
            "p7_rows_remaining_scalar_or_source_mass_debt"
        ],
        "p7_rows_converted_to_scalar_source_rows": len(rows),
        "unique_scalar_source_sectors": len(sectors),
        "operator_status_counts": dict(status_counts),
        "standard_bnv_seed_rows": len(standard_bnv),
        "sterile_bnv_rows": len(sterile_bnv),
        "baryon_violating_or_sterile_bnv_rows": len(bnv),
        "physical_proton_bounds_evaluable_now": 0,
        "all_scalar_source_masses_symbolic": True,
    }


def build_payload() -> dict:
    p7 = load_json(P7_JSON)
    p11 = load_json(P11_JSON)
    scalar_rows = [scalar_row(row) for row in p7["matching_rows"]]
    sectors = source_sector_rows(p7["matching_rows"])
    return {
        "description": "Route-C P12-S scalar/source chiral-pair branch",
        "branch_decision": {
            "both_branches_to_try": [
                "completed_current_current_vector_branch",
                "scalar_source_chiral_pair_branch",
            ],
            "activated_first": "scalar_source_chiral_pair_branch",
            "reason": (
                "P11 showed that all 18 P7 chiral-pair pole rows remain "
                "scalar/source mass debt rather than adjoint-vector mass matches."
            ),
        },
        "boundary": (
            "P12-S assigns scalar/source mass denominators to the P7 chiral-pair "
            "poles.  It does not compute physical proton lifetimes, and it does "
            "not close the separate current-current vector branch."
        ),
        "ansatz": {
            "lagrangian_template": "L ⊃ lambda_{ab,R} psi_a psi_b S_R + h.c. + M_{S_R}^2 |S_R|^2",
            "integrated_out_template": "C6(12|34;R)=lambda_{12,R} lambda^*_{34,R}/M_{S_R}^2",
            "positivity_condition": "M_{S_R}^2 > 0 for every activated source sector",
            "ward_boundary": "No broken-vector longitudinal Ward problem is introduced by scalar/source exchange.",
        },
        "scalar_source_sectors": sectors,
        "scalar_source_matching_rows": scalar_rows,
        "summary": summarize(scalar_rows, sectors, p11),
        "verification": {
            "all_p11_scalar_source_debt_rows_accounted_for": len(scalar_rows)
            == p11["summary"]["p7_rows_remaining_scalar_or_source_mass_debt"],
            "every_scalar_row_has_mass_symbol": all(row["mass_symbol"] for row in scalar_rows),
            "every_source_sector_has_positive_mass_condition": all(
                row["positivity_condition"].endswith("> 0") for row in sectors
            ),
            "physical_proton_bounds_evaluable_now": 0,
            "p12_scalar_source_layer_passes": (
                len(scalar_rows)
                == p11["summary"]["p7_rows_remaining_scalar_or_source_mass_debt"]
                and all(row["mass_symbol"] for row in scalar_rows)
            ),
        },
        "next_stage": {
            "recommended": "P13_scalar_source_flavor_tensor_and_operator_basis",
            "required_decision": (
                "Choose scalar/source flavor tensors and a chiral/Fierz Wilson "
                "basis for the BNV rows, or return to the parallel current-vector "
                "branch before inserting proton limits."
            ),
        },
    }


def write_markdown(payload: dict) -> str:
    lines = [
        "# Route-C P12-S Scalar/Source Chiral-Pair Branch",
        "",
        payload["boundary"],
        "",
        "## Branch Decision",
        "",
        f"- branches to try: `{', '.join(payload['branch_decision']['both_branches_to_try'])}`",
        f"- activated first: `{payload['branch_decision']['activated_first']}`",
        f"- reason: {payload['branch_decision']['reason']}",
        "",
        "## Ansatz",
        "",
        "```text",
        payload["ansatz"]["lagrangian_template"],
        payload["ansatz"]["integrated_out_template"],
        payload["ansatz"]["positivity_condition"],
        "```",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
    ]
    for key, value in payload["summary"].items():
        rendered = json.dumps(value, sort_keys=True) if isinstance(value, dict) else value
        lines.append(f"| {key} | {rendered} |")
    lines.extend(
        [
            "",
            "## Scalar/Source Sectors",
            "",
            "| source | mass | rows | operators |",
            "| --- | --- | ---: | --- |",
        ]
    )
    for row in payload["scalar_source_sectors"]:
        ops = ", ".join(f"`{op}`" for op in row["operators_supported"])
        lines.append(
            f"| `{row['source_field']}` | `{row['mass_symbol']}` | {row['row_count']} | {ops} |"
        )
    lines.extend(
        [
            "",
            "## BNV Scalar/Source Rows",
            "",
            "| operator | source | coefficient |",
            "| --- | --- | --- |",
        ]
    )
    for row in payload["scalar_source_matching_rows"]:
        if row["operator_status"] == "B_and_L_conserving":
            continue
        lines.append(
            f"| `{row['canonical_operator']}` | `{row['source_field']}` | "
            f"`{row['scalar_source_matching_coefficient']}` |"
        )
    lines.extend(
        [
            "",
            "## Verification",
            "",
            "| check | value |",
            "| --- | ---: |",
        ]
    )
    for key, value in payload["verification"].items():
        lines.append(f"| {key} | {value} |")
    lines.extend(
        [
            "",
            "## Next Stage",
            "",
            f"`{payload['next_stage']['recommended']}`",
            "",
            payload["next_stage"]["required_decision"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "activated_first": payload["branch_decision"]["activated_first"],
                "p7_rows_converted_to_scalar_source_rows": payload["summary"][
                    "p7_rows_converted_to_scalar_source_rows"
                ],
                "unique_scalar_source_sectors": payload["summary"][
                    "unique_scalar_source_sectors"
                ],
                "baryon_violating_or_sterile_bnv_rows": payload["summary"][
                    "baryon_violating_or_sterile_bnv_rows"
                ],
                "physical_proton_bounds_evaluable_now": payload["summary"][
                    "physical_proton_bounds_evaluable_now"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
