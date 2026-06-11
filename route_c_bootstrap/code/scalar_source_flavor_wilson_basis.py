#!/usr/bin/env python3
"""Route-C P13-S scalar/source flavor tensors and Wilson basis.

P13-S chooses a conservative symbolic flavor-tensor basis for the scalar/source
branch activated in P12-S.  It maps every BNV or sterile-BNV scalar/source row
to a chiral all-left Weyl operator basis and records which rows still need
Fierz/color recoupling before numerical proton limits can be inserted.

This is not a proton-lifetime calculation.
"""

from __future__ import annotations

import json
import re
from collections import OrderedDict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P12_JSON = OUT / "scalar_source_chiral_pair_branch.json"
OUT_JSON = OUT / "scalar_source_flavor_wilson_basis.json"
OUT_MD = OUT / "scalar_source_flavor_wilson_basis.md"


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def sanitize(text: str) -> str:
    return (
        text.replace("^c", "c")
        .replace("^", "")
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


def pair_key(pair: list[str]) -> str:
    return "".join(sanitize(x) for x in pair)


def mediator_su2(mediator: str) -> str:
    match = re.match(r"X\[\([^,]+,([^\)]+)\);", mediator)
    if not match:
        return "unknown"
    return match.group(1)


BASIS = {
    "O_QQQL_S": {
        "canonical_operator": "QQQL",
        "basis_expression": (
            "epsilon_abc epsilon_ij epsilon_kl "
            "(Q_p^{a i} Q_r^{b j})(Q_s^{c k} L_t^l)"
        ),
        "description": "QQ weak-singlet scalar/source contraction.",
        "family_symmetry": "C_{prst}=C_{rpst} for the first two Q-family indices",
        "independent_family_components": 54,
        "proton_channels_seeded": ["p -> e+ pi0", "p -> anti-nu pi+"],
    },
    "O_QQQL_T": {
        "canonical_operator": "QQQL",
        "basis_expression": (
            "epsilon_abc (tau^A epsilon)_ij (tau^A epsilon)_kl "
            "(Q_p^{a i} Q_r^{b j})(Q_s^{c k} L_t^l)"
        ),
        "description": "QQ weak-triplet scalar/source contraction.",
        "family_symmetry": "C_{prst}=-C_{rpst} for the first two Q-family indices",
        "independent_family_components": 27,
        "proton_channels_seeded": ["p -> e+ pi0", "p -> anti-nu pi+"],
    },
    "O_UUDE": {
        "canonical_operator": "u^c u^c d^c e^c",
        "basis_expression": (
            "epsilon_abc (u^c_p^a u^c_r^b)(d^c_s^c e^c_t)"
        ),
        "description": "Right-handed-conjugate all-left BNV operator.",
        "family_symmetry": "C_{prst}=-C_{rpst} for the two u^c-family indices",
        "independent_family_components": 27,
        "proton_channels_seeded": ["p -> e+ meson"],
    },
    "O_UDDN": {
        "canonical_operator": "u^c d^c d^c nu^c",
        "basis_expression": (
            "epsilon_abc (u^c_p^a d^c_r^b)(d^c_s^c nu^c_t)"
        ),
        "description": "Sterile-neutrino BNV operator; proton interpretation needs nu^c treatment.",
        "family_symmetry": "C_{prst}=-C_{psrt} for the two d^c-family indices after canonical projection",
        "independent_family_components": 27,
        "proton_channels_seeded": ["conditional: requires nu^c mass/seesaw or light-sterile treatment"],
    },
}


def choose_basis(row: dict) -> tuple[str, str]:
    op = row["canonical_operator"]
    med = row["mediator"]
    left = row["left_pair"]
    if op == "QQQL":
        if mediator_su2(med) == "3":
            return "O_QQQL_T", "direct_source_channel"
        return "O_QQQL_S", "direct_source_channel"
    if op == "u^c u^c d^c e^c":
        if left == ["u^c", "u^c"]:
            return "O_UUDE", "direct_source_channel"
        return "O_UUDE", "requires_chiral_fierz_color_recoupling"
    if op == "u^c d^c d^c nu^c":
        if left == ["u^c", "d^c"]:
            return "O_UDDN", "direct_source_channel"
        return "O_UDDN", "requires_chiral_fierz_color_recoupling"
    raise ValueError(f"no P13-S BNV basis for operator {op}")


def tensor_symmetry(pair: list[str], row: dict) -> tuple[str, int]:
    med = row["mediator"]
    if pair == ["Q", "Q"]:
        if mediator_su2(med) == "3":
            return "antisymmetric_family_pair", 3
        return "symmetric_family_pair", 6
    if pair == ["u^c", "u^c"]:
        return "antisymmetric_family_pair", 3
    if pair == ["d^c", "d^c"]:
        return "antisymmetric_family_pair", 3
    if pair == ["L", "L"]:
        return "model_dependent_identical_pair_not_in_BNV_basis", 0
    return "general_3x3_family_matrix", 9


def tensor_symbol(pair: list[str], row: dict) -> str:
    return "Lambda_" + pair_key(pair) + "__" + row["source_field"]


def tensor_entry(pair: list[str], row: dict) -> dict:
    symmetry, components = tensor_symmetry(pair, row)
    return {
        "tensor": tensor_symbol(pair, row),
        "source_field": row["source_field"],
        "pair": pair,
        "mediator": row["mediator"],
        "family_indices": ["i", "j"],
        "symmetry": symmetry,
        "independent_complex_components": components,
    }


def coefficient_template(row: dict, basis_id: str, recoupling_status: str) -> str:
    left_tensor = tensor_symbol(row["left_pair"], row)
    right_tensor = tensor_symbol(row["right_pair"], row)
    mass = row["mass_symbol"]
    if recoupling_status == "requires_chiral_fierz_color_recoupling":
        rho = "rho_" + basis_id + "__" + row["source_field"]
        return (
            f"C[{basis_id}]_prst += {rho} "
            f"{left_tensor}_ij {right_tensor}_kl^*/{mass}^2"
        )
    if basis_id == "O_UDDN":
        return (
            f"C[{basis_id}]_prst += Pi_d_asym_rs["
            f"{left_tensor}_pr {right_tensor}_st^*/{mass}^2]"
        )
    return f"C[{basis_id}]_prst += {left_tensor}_ij {right_tensor}_kl^*/{mass}^2"


def projection_status(basis_id: str, recoupling_status: str) -> str:
    if recoupling_status == "requires_chiral_fierz_color_recoupling":
        return "canonical_projection_included_in_recoupling_matrix"
    if basis_id == "O_UDDN":
        return "requires_d_pair_antisymmetrization_projection"
    return "encoded_by_vertex_pair_symmetry"


def build_payload() -> dict:
    p12 = load_json(P12_JSON)
    bnv_rows = [
        row
        for row in p12["scalar_source_matching_rows"]
        if row["operator_status"] != "B_and_L_conserving"
    ]

    tensors: OrderedDict[str, dict] = OrderedDict()
    wilson_rows = []
    for row in bnv_rows:
        basis_id, recoupling_status = choose_basis(row)
        for pair in (row["left_pair"], row["right_pair"]):
            entry = tensor_entry(pair, row)
            tensors.setdefault(entry["tensor"], entry)
        wilson_rows.append(
            {
                "canonical_operator": row["canonical_operator"],
                "operator_status": row["operator_status"],
                "basis_id": basis_id,
                "source_field": row["source_field"],
                "mediator": row["mediator"],
                "mass_symbol": row["mass_symbol"],
                "left_pair": row["left_pair"],
                "right_pair": row["right_pair"],
                "left_tensor": tensor_symbol(row["left_pair"], row),
                "right_tensor": tensor_symbol(row["right_pair"], row),
                "coefficient_template": coefficient_template(row, basis_id, recoupling_status),
                "recoupling_status": recoupling_status,
                "canonical_projection_status": projection_status(basis_id, recoupling_status),
                "physical_bound_evaluable_now": False,
            }
        )

    used_basis = sorted({row["basis_id"] for row in wilson_rows})
    recoupling_rows = [
        row for row in wilson_rows if row["recoupling_status"] != "direct_source_channel"
    ]
    direct_rows = [
        row for row in wilson_rows if row["recoupling_status"] == "direct_source_channel"
    ]
    summary = {
        "p13_branch": "scalar_source_flavor_tensor_and_wilson_basis",
        "branch_v_still_required": True,
        "bnv_or_sterile_bnv_rows_processed": len(wilson_rows),
        "canonical_wilson_basis_count": len(used_basis),
        "canonical_wilson_basis_ids": used_basis,
        "direct_source_channel_rows": len(direct_rows),
        "fierz_or_color_recoupling_required_rows": len(recoupling_rows),
        "canonical_family_projection_required_rows": sum(
            1
            for row in wilson_rows
            if row["canonical_projection_status"]
            == "requires_d_pair_antisymmetrization_projection"
        ),
        "flavor_tensors_introduced": len(tensors),
        "independent_complex_vertex_components": sum(
            entry["independent_complex_components"] for entry in tensors.values()
        ),
        "physical_proton_bounds_evaluable_now": 0,
    }
    verification = {
        "all_p12_bnv_rows_mapped": len(wilson_rows) == 6,
        "every_row_has_basis_id": all(row["basis_id"] in BASIS for row in wilson_rows),
        "every_row_has_left_and_right_flavor_tensor": all(
            row["left_tensor"] and row["right_tensor"] for row in wilson_rows
        ),
        "branch_v_retained": True,
        "physical_proton_bounds_evaluable_now": 0,
        "p13_scalar_source_basis_layer_passes": (
            len(wilson_rows) == 6
            and all(row["basis_id"] in BASIS for row in wilson_rows)
            and all(row["left_tensor"] and row["right_tensor"] for row in wilson_rows)
        ),
    }
    return {
        "description": "Route-C P13-S scalar/source flavor tensors and chiral Wilson basis",
        "boundary": (
            "P13-S chooses a symbolic flavor-tensor and chiral/Fierz Wilson basis "
            "for the P12-S BNV rows.  It does not choose numerical couplings, "
            "physical flavor rotations, RG factors, hadronic matrix elements, "
            "or proton lifetime limits.  Branch V remains open."
        ),
        "basis_conventions": {
            "spinor_contraction": "(psi chi)=psi^alpha chi_alpha for left-handed Weyl fields",
            "family_index_range": "p,r,s,t = 1,2,3",
            "color_index_range": "a,b,c = 1,2,3",
            "weak_index_range": "i,j,k,l = 1,2",
            "recoupling_identity_boundary": (
                "Rows marked requires_chiral_fierz_color_recoupling need an explicit "
                "two-component spinor Fierz and color recoupling matrix before numerical bounds."
            ),
            "canonical_projection_boundary": (
                "Rows marked requires_d_pair_antisymmetrization_projection need the "
                "explicit antisymmetric projection on the two d^c family indices."
            ),
        },
        "canonical_wilson_basis": {key: BASIS[key] for key in used_basis},
        "flavor_tensors": list(tensors.values()),
        "wilson_matching_rows": wilson_rows,
        "summary": summary,
        "verification": verification,
        "next_stage": {
            "recommended": "P14_scalar_source_recoupling_matrix_or_branch_v",
            "required_decision": (
                "Either compute the chiral/color Fierz recoupling matrices for the two "
                "non-direct scalar/source rows and then add flavor rotations, or return "
                "to the completed current-current vector Branch V before proton limits."
            ),
        },
    }


def write_markdown(payload: dict) -> str:
    lines = [
        "# Route-C P13-S Scalar/Source Flavor Tensors and Wilson Basis",
        "",
        payload["boundary"],
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
    ]
    for key, value in payload["summary"].items():
        rendered = json.dumps(value, sort_keys=True) if isinstance(value, list) else value
        lines.append(f"| {key} | {rendered} |")

    lines.extend(
        [
            "",
            "## Canonical Chiral Wilson Basis",
            "",
            "| basis | operator | independent components | symmetry |",
            "| --- | --- | ---: | --- |",
        ]
    )
    for key, row in payload["canonical_wilson_basis"].items():
        lines.append(
            f"| `{key}` | `{row['basis_expression']}` | "
            f"{row['independent_family_components']} | {row['family_symmetry']} |"
        )

    lines.extend(
        [
            "",
            "## Flavor Tensors",
            "",
            "| tensor | pair | source | symmetry | independent components |",
            "| --- | --- | --- | --- | ---: |",
        ]
    )
    for row in payload["flavor_tensors"]:
        pair = " ".join(row["pair"])
        lines.append(
            f"| `{row['tensor']}` | `{pair}` | `{row['source_field']}` | "
            f"{row['symmetry']} | {row['independent_complex_components']} |"
        )

    lines.extend(
        [
            "",
            "## BNV Wilson Matching Rows",
            "",
            "| operator | basis | source | recoupling | projection | coefficient |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload["wilson_matching_rows"]:
        lines.append(
            f"| `{row['canonical_operator']}` | `{row['basis_id']}` | "
            f"`{row['source_field']}` | {row['recoupling_status']} | "
            f"{row['canonical_projection_status']} | "
            f"`{row['coefficient_template']}` |"
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
                "bnv_or_sterile_bnv_rows_processed": payload["summary"][
                    "bnv_or_sterile_bnv_rows_processed"
                ],
                "canonical_wilson_basis_count": payload["summary"][
                    "canonical_wilson_basis_count"
                ],
                "direct_source_channel_rows": payload["summary"][
                    "direct_source_channel_rows"
                ],
                "fierz_or_color_recoupling_required_rows": payload["summary"][
                    "fierz_or_color_recoupling_required_rows"
                ],
                "canonical_family_projection_required_rows": payload["summary"][
                    "canonical_family_projection_required_rows"
                ],
                "flavor_tensors_introduced": payload["summary"][
                    "flavor_tensors_introduced"
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
