#!/usr/bin/env python3
"""Route-C P14-V Spin(10) vector-branch matching gate.

P14-V returns to the current-current vector branch before proton limits are
inserted.  It consumes the P9 broken-generator maps, the P10 staged mass
blocks, and the P11 replay gate.  The output is a vector-current matching
interface:

  L_eff[X] = - g_X^2 J_X^mu J_{X,mu}^\dagger / M_X^2,

with one row for every current-pair structure inside each broken generator
map.  It deliberately does not convert these current-current rows into
physical proton lifetimes.  That still requires a convention-fixed all-left
Fierz map, physical flavor rotations, RG factors, hadronic matrix elements,
and channel limits.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P9_JSON = OUT / "spin10_action_completion_ledger.json"
P10_JSON = OUT / "spin10_breaking_mass_matrix.json"
P11_JSON = OUT / "spin10_mass_replay_matching_gate.json"
OUT_JSON = OUT / "spin10_vector_branch_matching_gate.json"
OUT_MD = OUT / "spin10_vector_branch_matching_gate.md"


MASS_BY_FAMILY = {
    "Pati-Salam SU(4)_C leptoquark": "M_LQ",
    "broken SU(2)_R charged current": "M_WR",
    "Spin(10)/Pati-Salam off-face generator": "M_(6,2,2)",
}


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def transition_signature(entry: dict) -> str:
    return f"{entry['source_multiplet']}->{entry['target_multiplet']}"


def bilinear_template(entry: dict) -> str:
    return (
        f"({entry['target_label']}^dagger bar_sigma_mu "
        f"{entry['source_label']})"
    )


def current_pair_gate(family: str, e1: dict, e2: dict) -> str:
    sigs = {transition_signature(e1), transition_signature(e2)}
    unordered = {
        tuple(sorted((e1["source_multiplet"], e1["target_multiplet"]))),
        tuple(sorted((e2["source_multiplet"], e2["target_multiplet"]))),
    }

    if family == "Pati-Salam SU(4)_C leptoquark":
        mixed_lq = (
            ("L", "Q") in unordered
            and (("nu^c", "u^c") in unordered or ("d^c", "e^c") in unordered)
        )
        if mixed_lq:
            return "ps_leptoquark_mixed_pair_crossing_required"
        return "ps_leptoquark_same-sector_pair_noncanonical_until_fierz"

    if family == "broken SU(2)_R charged current":
        return "right_current_B_and_L_conserving_not_proton_seed"

    if family == "Spin(10)/Pati-Salam off-face generator":
        if sigs == {"Q->u^c"} or sigs == {"Q->d^c"}:
            return "off_face_diquark_like_pair_crossing_required"
        if sigs in ({"L->u^c", "Q->nu^c"}, {"L->d^c", "Q->e^c"}):
            return "off_face_mixed_pair_crossing_required"
        return "off_face_current_pair_manual_fierz_required"

    return "unknown_vector_current_pair"


def build_current_rows(p9: dict, p10: dict) -> list[dict]:
    mass_replacements = p10["p6_p7_replay_interface"]["replace_symbolic_masses"]
    rows: list[dict] = []
    row_id = 0
    for vmap in p9["broken_generator_maps"]["maps"]:
        family = vmap["generator_family"]
        mass_block = MASS_BY_FAMILY[family]
        mass_expression = mass_replacements[mass_block]
        entries = vmap["sparse_transition_entries"]
        for e1, e2 in itertools.combinations(entries, 2):
            gate = current_pair_gate(family, e1, e2)
            rows.append(
                {
                    "row_id": row_id,
                    "vector_label": vmap["label"],
                    "generator_family": family,
                    "mass_block": mass_block,
                    "mass_expression": mass_expression,
                    "transition_1": {
                        "source": e1["source_label"],
                        "target": e1["target_label"],
                        "signature": transition_signature(e1),
                        "coefficient": e1["coefficient"],
                        "delta_Y": e1["delta_Y"]["fraction"],
                        "delta_B_minus_L": e1["delta_B_minus_L"]["fraction"],
                    },
                    "transition_2": {
                        "source": e2["source_label"],
                        "target": e2["target_label"],
                        "signature": transition_signature(e2),
                        "coefficient": e2["coefficient"],
                        "delta_Y": e2["delta_Y"]["fraction"],
                        "delta_B_minus_L": e2["delta_B_minus_L"]["fraction"],
                    },
                    "current_current_template": (
                        f"- g_X^2/{mass_block}^2 "
                        f"{bilinear_template(e1)} {bilinear_template(e2)}^dagger"
                    ),
                    "matching_gate": gate,
                    "canonical_proton_operator_status": (
                        "not_yet_in_all_left_basis"
                    ),
                    "required_before_proton_limit": [
                        "fix vector-current to all-left Fierz and crossing convention",
                        "apply physical flavor rotations",
                        "insert numerical masses and couplings",
                        "run to hadronic scale",
                        "contract with lattice or chiral hadronic matrix elements",
                        "choose experimental proton channel limit",
                    ],
                }
            )
            row_id += 1
    return rows


def current_map_summary(p9: dict) -> dict:
    family_counts = Counter()
    entry_counts = Counter()
    for vmap in p9["broken_generator_maps"]["maps"]:
        family = vmap["generator_family"]
        family_counts[family] += 1
        entry_counts[family] += len(vmap["sparse_transition_entries"])
    return {
        "broken_vector_maps": sum(family_counts.values()),
        "broken_vector_maps_by_family": dict(family_counts),
        "current_transition_entries": sum(entry_counts.values()),
        "current_transition_entries_by_family": dict(entry_counts),
    }


def summarize(rows: list[dict], p9: dict, p10: dict, p11: dict) -> dict:
    gates = Counter(row["matching_gate"] for row in rows)
    families = Counter(row["generator_family"] for row in rows)
    p11_summary = p11["summary"]
    return {
        **current_map_summary(p9),
        "current_pair_rows": len(rows),
        "current_pair_rows_by_family": dict(families),
        "matching_gate_counts": dict(gates),
        "p10_mass_matrix_layer_passes": p10["verification"][
            "p10_mass_matrix_layer_passes"
        ],
        "p11_p7_rows_remaining_scalar_or_source_mass_debt": p11_summary[
            "p7_rows_remaining_scalar_or_source_mass_debt"
        ],
        "p11_p7_rows_with_candidate_adjoint_mass_blocks": p11_summary[
            "p7_rows_with_candidate_adjoint_mass_blocks"
        ],
        "vector_branch_does_not_reuse_p7_scalar_source_poles": True,
        "physical_proton_bounds_evaluable_now": 0,
        "branch_v_matching_gate_passes": (
            p10["verification"]["p10_mass_matrix_layer_passes"]
            and p11_summary["p7_rows_with_candidate_adjoint_mass_blocks"] == 0
            and len(rows) > 0
        ),
    }


def build_payload() -> dict:
    p9 = load_json(P9_JSON)
    p10 = load_json(P10_JSON)
    p11 = load_json(P11_JSON)
    rows = build_current_rows(p9, p10)
    summary = summarize(rows, p9, p10, p11)
    return {
        "description": "Route-C P14-V Spin(10) vector-branch matching gate",
        "boundary": (
            "P14-V returns to Branch V and builds a current-current vector "
            "matching gate before proton limits.  It does not reuse the P7 "
            "scalar/source chiral-pair poles, and it does not compute physical "
            "proton lifetimes."
        ),
        "vector_matching_formula": {
            "current": "J_X^mu = sum_a c_a psi_target(a)^dagger bar_sigma^mu psi_source(a)",
            "tree_level_effective_lagrangian": (
                "L_eff[X] = - g_X^2 J_X^mu J_X,mu^dagger / M_X^2"
            ),
            "mass_blocks": p10["p6_p7_replay_interface"][
                "replace_symbolic_masses"
            ],
        },
        "summary": summary,
        "current_pair_rows": rows,
        "gate_interpretation": {
            "ps_leptoquark_mixed_pair_crossing_required": (
                "Pati-Salam leptoquark vector current pair with one Q/L "
                "transition and one conjugate quark/lepton transition.  It is "
                "a vector-branch candidate, but still needs an all-left Fierz "
                "map and flavor rotations before proton limits."
            ),
            "off_face_diquark_like_pair_crossing_required": (
                "Spin(10)/Pati-Salam off-face vector current pair containing "
                "two Q-to-u^c or Q-to-d^c transitions.  It is not a P7 scalar "
                "pole; it needs current-current Fierz matching."
            ),
            "off_face_mixed_pair_crossing_required": (
                "Spin(10)/Pati-Salam off-face mixed quark/lepton current pair. "
                "It remains a vector-current matching candidate, not yet a "
                "physical proton-width row."
            ),
            "right_current_B_and_L_conserving_not_proton_seed": (
                "Broken SU(2)_R charged current exchange is tracked for "
                "completeness but is not by itself a proton-decay seed."
            ),
        },
        "proton_limit_blockers": [
            "current-current to all-left chiral Wilson basis recoupling",
            "physical flavor rotations",
            "numerical vector masses and couplings",
            "RG running from matching scale to hadronic scale",
            "lattice or chiral hadronic matrix elements",
            "experimental channel selection and lifetime limit",
        ],
        "parallel_branch_status": {
            "scalar_source_branch": "P13-S retained as a separate scalar/source Wilson-basis branch",
            "vector_branch": "Branch V current-current matching gate is now explicit",
        },
        "next_stage": {
            "recommended": "P15_V_current_to_all_left_fierz_and_flavor_gate",
            "required_decision": (
                "Choose the convention-fixed current-current to all-left "
                "operator map for the vector branch, or keep P13-S scalar/source "
                "recoupling as the active proton-matching path.  Do not insert "
                "proton limits until one branch has flavor, RG, and hadronic "
                "inputs."
            ),
        },
    }


def markdown_table(rows: list[list[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(str(row[i])) for row in rows) for i in range(len(rows[0]))]
    out = []
    for idx, row in enumerate(rows):
        out.append(
            "| "
            + " | ".join(str(cell).ljust(widths[i]) for i, cell in enumerate(row))
            + " |"
        )
        if idx == 0:
            out.append("| " + " | ".join("-" * w for w in widths) + " |")
    return out


def write_markdown(payload: dict) -> str:
    summary_rows = [["quantity", "value"]]
    for key, value in payload["summary"].items():
        if isinstance(value, (dict, list)):
            display = "`" + json.dumps(value, sort_keys=True) + "`"
        else:
            display = str(value)
        summary_rows.append([key, display])

    gate_rows = [["gate", "count"]]
    for gate, count in sorted(payload["summary"]["matching_gate_counts"].items()):
        gate_rows.append([f"`{gate}`", str(count)])

    sample_rows = [
        ["family", "mass", "transition pair", "gate"],
    ]
    for row in payload["current_pair_rows"][:18]:
        pair = (
            f"{row['transition_1']['signature']} ; "
            f"{row['transition_2']['signature']}"
        )
        sample_rows.append(
            [
                row["generator_family"],
                f"`{row['mass_block']}`",
                f"`{pair}`",
                f"`{row['matching_gate']}`",
            ]
        )

    lines = [
        "# Route-C P14-V Spin(10) Vector-Branch Matching Gate",
        "",
        payload["boundary"],
        "",
        "## Matching Formula",
        "",
        f"- current: `{payload['vector_matching_formula']['current']}`",
        "- tree-level gate: "
        f"`{payload['vector_matching_formula']['tree_level_effective_lagrangian']}`",
        "",
        "## Summary",
        "",
        *markdown_table(summary_rows),
        "",
        "## Gate Counts",
        "",
        *markdown_table(gate_rows),
        "",
        "## Sample Current-Pair Rows",
        "",
        *markdown_table(sample_rows),
        "",
        "## Proton-Limit Blockers",
        "",
    ]
    for item in payload["proton_limit_blockers"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## Parallel Branch Status",
            "",
            f"- scalar/source branch: {payload['parallel_branch_status']['scalar_source_branch']}",
            f"- vector branch: {payload['parallel_branch_status']['vector_branch']}",
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
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "current_pair_rows": payload["summary"]["current_pair_rows"],
                "matching_gate_counts": payload["summary"]["matching_gate_counts"],
                "p11_p7_rows_with_candidate_adjoint_mass_blocks": payload[
                    "summary"
                ]["p11_p7_rows_with_candidate_adjoint_mass_blocks"],
                "physical_proton_bounds_evaluable_now": payload["summary"][
                    "physical_proton_bounds_evaluable_now"
                ],
                "branch_v_matching_gate_passes": payload["summary"][
                    "branch_v_matching_gate_passes"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
