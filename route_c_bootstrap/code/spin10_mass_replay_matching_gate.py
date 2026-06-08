#!/usr/bin/env python3
"""Route-C P11 replay of P6/P7 with staged Spin(10) mass blocks.

P11 consumes the P10 mass-matrix ledger and asks what actually changes in the
P6/P7 matching gates.  The answer is deliberately conservative:

* the P9/P10 adjoint current-vector sectors now have explicit mass blocks;
* the P7 chiral-pair pole rows are tested against those adjoint-vector blocks;
* rows that do not match an adjoint current vector remain scalar/source mass
  debt and are not promoted to proton-bound calculations.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P7_JSON = OUT / "low_energy_matching_proton_report.json"
P9_JSON = OUT / "spin10_action_completion_ledger.json"
P10_JSON = OUT / "spin10_breaking_mass_matrix.json"
OUT_JSON = OUT / "spin10_mass_replay_matching_gate.json"
OUT_MD = OUT / "spin10_mass_replay_matching_gate.md"


MED_RE = re.compile(
    r"X\[\((?P<su3>[^,]+),(?P<su2>[^\)]+)\);Y=(?P<Y>[^;]+);B-L=(?P<BL>[^\]]+)\]"
)


MASS_BY_FAMILY = {
    "Pati-Salam SU(4)_C leptoquark": "M_LQ",
    "broken SU(2)_R charged current": "M_WR",
    "Spin(10)/Pati-Salam off-face generator": "M_(6,2,2)",
}


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_mediator(label: str) -> dict[str, str]:
    match = MED_RE.fullmatch(label)
    if not match:
        raise ValueError(f"cannot parse mediator label: {label}")
    return match.groupdict()


def adjoint_charge_signatures(p9: dict) -> dict[tuple[str, str], set[str]]:
    signatures: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in p9["broken_generator_maps"]["maps"]:
        mass = MASS_BY_FAMILY[row["generator_family"]]
        for entry in row["sparse_transition_entries"]:
            key = (entry["delta_Y"]["fraction"], entry["delta_B_minus_L"]["fraction"])
            signatures[key].add(mass)
    return signatures


def current_vector_replay_rows(p9: dict) -> list[dict]:
    grouped: dict[str, dict] = {}
    for row in p9["broken_generator_maps"]["maps"]:
        mass = MASS_BY_FAMILY[row["generator_family"]]
        entry = grouped.setdefault(
            mass,
            {
                "mass_block": mass,
                "generator_family": row["generator_family"],
                "map_count": 0,
                "transition_pair_signatures": set(),
                "matching_template": None,
            },
        )
        entry["map_count"] += 1
        for tr in row["sparse_transition_entries"]:
            entry["transition_pair_signatures"].add(
                (
                    tr["source_multiplet"],
                    tr["target_multiplet"],
                    tr["delta_Y"]["fraction"],
                    tr["delta_B_minus_L"]["fraction"],
                )
            )

    rows = []
    for mass, row in sorted(grouped.items()):
        rows.append(
            {
                "mass_block": mass,
                "generator_family": row["generator_family"],
                "map_count": row["map_count"],
                "transition_signature_count": len(row["transition_pair_signatures"]),
                "matching_template": (
                    f"C_current[{mass}] = g_10^2 T_X T_X^*/{mass}^2 "
                    "after Goldstone/Higgs Ward completion"
                ),
                "p6_status_after_p10": (
                    "mass_denominator_supplied; Goldstone/Higgs amplitude and "
                    "spin numerator still required for a full high-energy check"
                ),
            }
        )
    rows.append(
        {
            "mass_block": "M_Zprime",
            "generator_family": "neutral U(1)_{T3R} x U(1)_{B-L} broken vector",
            "map_count": 1,
            "transition_signature_count": 1,
            "matching_template": (
                "C_current[M_Zprime] = g_Zprime^2 Q_Zprime Q_Zprime^*/M_Zprime^2"
            ),
            "p6_status_after_p10": (
                "neutral mass denominator supplied; physical charges and "
                "normalization must be chosen before low-energy matching"
            ),
        }
    )
    return rows


def p7_row_replay(row: dict, signatures: dict[tuple[str, str], set[str]], p10: dict) -> dict:
    mediator = row["mediator"]
    parsed = parse_mediator(mediator)
    key = (parsed["Y"], parsed["BL"])
    direct_blocks = sorted(signatures.get(key, set()))
    conjugate_key = (neg_fraction_string(parsed["Y"]), neg_fraction_string(parsed["BL"]))
    conjugate_blocks = sorted(signatures.get(conjugate_key, set()))
    if direct_blocks or conjugate_blocks:
        gate = "charge_signature_matches_adjoint_current_but_channel_spin_must_be_checked"
        mass_blocks = direct_blocks or conjugate_blocks
        replacement = [
            p10["p6_p7_replay_interface"]["replace_symbolic_masses"][m]
            for m in mass_blocks
            if m in p10["p6_p7_replay_interface"]["replace_symbolic_masses"]
        ]
        replay_coefficient = row["symbolic_matching_coefficient"]
        for block, repl in zip(mass_blocks, replacement):
            replay_coefficient = replay_coefficient.replace(
                f"M_{mediator}", f"({repl})"
            )
    else:
        gate = "not_an_adjoint_current_mass_match"
        mass_blocks = []
        replacement = []
        replay_coefficient = (
            row["symbolic_matching_coefficient"]
            .replace(f"M_{mediator}", f"M_source[{mediator}]")
        )

    return {
        "canonical_operator": row["canonical_operator"],
        "operator_status": row["operator_status"],
        "mediator": mediator,
        "mediator_quantum_numbers": parsed,
        "symbolic_matching_coefficient": row["symbolic_matching_coefficient"],
        "p10_adjoint_mass_gate": gate,
        "candidate_adjoint_mass_blocks": mass_blocks,
        "candidate_mass_replacements": replacement,
        "p11_matching_coefficient": replay_coefficient,
        "physical_proton_bound_status_after_p11": "not_evaluable_yet",
        "reason_physical_bound_not_evaluable": (
            "P11 supplies staged vector mass gates, but physical bounds still "
            "require a completed spin/channel assignment, flavor rotations, RG "
            "factors, hadronic matrix elements, and experimental channel limits."
        ),
    }


def neg_fraction_string(value: str) -> str:
    value = value.strip()
    if value == "0":
        return "0"
    if value.startswith("-"):
        return value[1:]
    return "-" + value


def summarize(rows: list[dict], current_rows: list[dict]) -> dict:
    status = Counter(row["operator_status"] for row in rows)
    gates = Counter(row["p10_adjoint_mass_gate"] for row in rows)
    bnv = [row for row in rows if row["operator_status"] != "B_and_L_conserving"]
    assigned = [
        row
        for row in rows
        if row["p10_adjoint_mass_gate"]
        != "not_an_adjoint_current_mass_match"
    ]
    return {
        "current_vector_mass_blocks_available": len(current_rows),
        "p7_rows_replayed": len(rows),
        "p7_operator_status_counts": dict(status),
        "p7_adjoint_mass_gate_counts": dict(gates),
        "p7_rows_with_candidate_adjoint_mass_blocks": len(assigned),
        "p7_rows_remaining_scalar_or_source_mass_debt": len(rows) - len(assigned),
        "baryon_violating_or_sterile_bnv_rows": len(bnv),
        "physical_proton_bounds_evaluable_now": 0,
    }


def build_payload() -> dict:
    p7 = load_json(P7_JSON)
    p9 = load_json(P9_JSON)
    p10 = load_json(P10_JSON)
    signatures = adjoint_charge_signatures(p9)
    current_rows = current_vector_replay_rows(p9)
    p7_rows = [p7_row_replay(row, signatures, p10) for row in p7["matching_rows"]]
    return {
        "description": "Route-C P11 P6/P7 replay with staged Spin(10) masses",
        "boundary": (
            "P11 attaches P10 staged mass blocks to the matching gates.  It does "
            "not convert chiral-pair symbolic poles into adjoint-vector exchange "
            "unless the mediator quantum numbers pass the adjoint-current gate, "
            "and it does not compute numerical proton lifetimes."
        ),
        "mass_replacement_dictionary": p10["p6_p7_replay_interface"][
            "replace_symbolic_masses"
        ],
        "current_vector_replay": current_rows,
        "adjoint_current_charge_signatures": {
            f"Y={key[0]},B-L={key[1]}": sorted(values)
            for key, values in sorted(signatures.items())
        },
        "p7_mass_gate_rows": p7_rows,
        "summary": summarize(p7_rows, current_rows),
        "p11_interpretation": {
            "what_improved": (
                "The adjoint current-vector sectors now carry explicit staged "
                "Spin(10) mass denominators M_(6,2,2), M_LQ, M_WR, and M_Zprime."
            ),
            "what_did_not_improve": (
                "The P7 chiral-pair pole rows do not automatically become adjoint "
                "current-vector exchanges.  Rows failing the adjoint gate remain "
                "scalar/source mass debt."
            ),
            "next_required_step": (
                "Either rewrite the relevant BNV amplitudes in a completed "
                "current-current vector basis, or choose scalar/source masses "
                "for the chiral-pair pole rows and then supply flavor/RG/hadronic "
                "inputs."
            ),
        },
        "next_stage": {
            "recommended": "P12_choose_vector_current_or_scalar_source_proton_branch",
            "required_decision": (
                "Decide whether proton matching proceeds through completed "
                "Spin(10) current-vector exchange or through scalar/source "
                "chiral-pair poles; the two branches require different Wilson "
                "basis maps."
            ),
        },
    }


def write_markdown(payload: dict) -> str:
    def md_value(value: object) -> str:
        if isinstance(value, (dict, list)):
            return "`" + json.dumps(value, sort_keys=True) + "`"
        return str(value)

    lines = [
        "# Route-C P11 Staged Spin(10) Mass Replay",
        "",
        payload["boundary"],
        "",
        "## Mass Replacement Dictionary",
        "",
        "| symbolic mass | P10 replacement |",
        "| --- | --- |",
    ]
    for key, value in payload["mass_replacement_dictionary"].items():
        lines.append(f"| `{key}` | `{value}` |")
    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| quantity | value |",
            "| --- | ---: |",
        ]
    )
    for key, value in payload["summary"].items():
        lines.append(f"| {key} | {md_value(value)} |")
    lines.extend(
        [
            "",
            "## Current-Vector Replay",
            "",
            "| mass block | maps | status |",
            "| --- | ---: | --- |",
        ]
    )
    for row in payload["current_vector_replay"]:
        lines.append(
            f"| `{row['mass_block']}` | {row['map_count']} | {row['p6_status_after_p10']} |"
        )
    lines.extend(
        [
            "",
            "## P7 Chiral-Pair Mass Gate",
            "",
            "| operator | mediator | gate | P11 denominator status |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in payload["p7_mass_gate_rows"]:
        status = (
            ", ".join(row["candidate_adjoint_mass_blocks"])
            if row["candidate_adjoint_mass_blocks"]
            else "scalar/source mass debt"
        )
        lines.append(
            f"| `{row['canonical_operator']}` | `{row['mediator']}` | "
            f"`{row['p10_adjoint_mass_gate']}` | {status} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- improved: {payload['p11_interpretation']['what_improved']}",
            f"- not improved: {payload['p11_interpretation']['what_did_not_improve']}",
            f"- next: {payload['p11_interpretation']['next_required_step']}",
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
                "current_vector_mass_blocks_available": payload["summary"][
                    "current_vector_mass_blocks_available"
                ],
                "p7_rows_replayed": payload["summary"]["p7_rows_replayed"],
                "p7_rows_with_candidate_adjoint_mass_blocks": payload["summary"][
                    "p7_rows_with_candidate_adjoint_mass_blocks"
                ],
                "p7_rows_remaining_scalar_or_source_mass_debt": payload["summary"][
                    "p7_rows_remaining_scalar_or_source_mass_debt"
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
