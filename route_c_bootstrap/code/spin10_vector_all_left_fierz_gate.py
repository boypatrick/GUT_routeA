#!/usr/bin/env python3
"""Route-C P15-V all-left Fierz/flavor gate for Branch V.

P15-V consumes the P14-V current-current vector rows and fixes a two-component
Weyl Fierz/crossing convention.  It maps each row

  (psi_t1^dagger bar_sigma psi_s1)(psi_s2^dagger bar_sigma psi_t2)

to the crossed all-left skeleton

  -2 g_X^2/M_X^2 (psi_s1 psi_t2)(bar(psi_t1) bar(psi_s2)).

This is a convention gate, not a proton lifetime calculation.  Rows that do
not land in a canonical BNV all-left basis remain blocked until a physical
vector multiplet assembly, flavor rotations, RG factors, hadronic matrix
elements, and channel limits are supplied.
"""

from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P14_JSON = OUT / "spin10_vector_branch_matching_gate.json"
OUT_JSON = OUT / "spin10_vector_all_left_fierz_gate.json"
OUT_MD = OUT / "spin10_vector_all_left_fierz_gate.md"


GLOBAL_CHARGES = {
    "Q": {"B": Fraction(1, 3), "L": Fraction(0)},
    "u^c": {"B": Fraction(-1, 3), "L": Fraction(0)},
    "d^c": {"B": Fraction(-1, 3), "L": Fraction(0)},
    "L": {"B": Fraction(0), "L": Fraction(1)},
    "e^c": {"B": Fraction(0), "L": Fraction(-1)},
    "nu^c": {"B": Fraction(0), "L": Fraction(-1)},
}

CANONICAL_BNV_BASIS = {
    ("Q", "Q", "Q", "L"): "O_QQQL",
    ("u^c", "u^c", "d^c", "e^c"): "O_UUDE",
    ("u^c", "d^c", "d^c", "nu^c"): "O_UDDN",
}


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run P14-V first")
    return json.loads(path.read_text(encoding="utf-8"))


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def parse_signature(signature: str) -> tuple[str, str]:
    source, target = signature.split("->", 1)
    return source, target


def anti(field: str) -> str:
    return f"bar({field})"


def charge(field: str, conjugate: bool = False) -> dict[str, Fraction]:
    base = GLOBAL_CHARGES[field]
    sign = -1 if conjugate else 1
    return {"B": sign * base["B"], "L": sign * base["L"]}


def charge_sum(fields: list[dict]) -> dict[str, Fraction]:
    b = Fraction(0)
    l = Fraction(0)
    for field in fields:
        q = charge(field["multiplet"], field["conjugate"])
        b += q["B"]
        l += q["L"]
    return {"B": b, "L": l, "B_minus_L": b - l}


def canonical_basis_id(fields: list[dict]) -> str | None:
    if any(field["conjugate"] for field in fields):
        return None
    ordered = tuple(field["multiplet"] for field in fields)
    if ordered in CANONICAL_BNV_BASIS:
        return CANONICAL_BNV_BASIS[ordered]
    sorted_key = tuple(sorted(ordered))
    for key, basis in CANONICAL_BNV_BASIS.items():
        if tuple(sorted(key)) == sorted_key:
            return basis
    return None


def classify(fields: list[dict]) -> str:
    q = charge_sum(fields)
    basis = canonical_basis_id(fields)
    if basis is not None:
        return "canonical_bnv_all_left_basis_candidate"
    if q["B"] == 0 and q["L"] == 0:
        return "B_and_L_conserving_crossed_skeleton"
    if q["B_minus_L"] == 0 and abs(q["B"]) == 1:
        return "BNV_BminusL_preserving_but_noncanonical_crossed_skeleton"
    if q["B_minus_L"] == 0:
        return "BminusL_preserving_noncanonical_crossed_skeleton"
    return "nonstandard_BminusL_crossed_skeleton"


def crossed_fields(row: dict) -> list[dict]:
    src1, tgt1 = parse_signature(row["transition_1"]["signature"])
    src2, tgt2 = parse_signature(row["transition_2"]["signature"])
    return [
        {"multiplet": src1, "conjugate": False, "slot": "psi_s1"},
        {"multiplet": tgt2, "conjugate": False, "slot": "psi_t2"},
        {"multiplet": tgt1, "conjugate": True, "slot": "bar(psi_t1)"},
        {"multiplet": src2, "conjugate": True, "slot": "bar(psi_s2)"},
    ]


def display_fields(fields: list[dict]) -> list[str]:
    return [
        anti(field["multiplet"]) if field["conjugate"] else field["multiplet"]
        for field in fields
    ]


def build_rows(p14: dict) -> list[dict]:
    out = []
    for row in p14["current_pair_rows"]:
        fields = crossed_fields(row)
        q = charge_sum(fields)
        status = classify(fields)
        basis = canonical_basis_id(fields)
        out.append(
            {
                "row_id": row["row_id"],
                "vector_label": row["vector_label"],
                "generator_family": row["generator_family"],
                "mass_block": row["mass_block"],
                "p14_matching_gate": row["matching_gate"],
                "crossed_all_left_fields": display_fields(fields),
                "fierz_convention": (
                    "(chi^dagger bar_sigma^mu psi)(eta^dagger bar_sigma_mu xi)"
                    " = 2 (chi^dagger eta^dagger)(psi xi)"
                ),
                "coefficient_template": (
                    f"C_crossed[row_{row['row_id']}] = "
                    f"-2 g_X^2 c1 c2^*/{row['mass_block']}^2"
                ),
                "global_charges": {
                    "B": frac(q["B"]),
                    "L": frac(q["L"]),
                    "B_minus_L": frac(q["B_minus_L"]),
                },
                "all_left_status": status,
                "canonical_bnv_basis_id": basis,
                "flavor_gate_status": "symbolic_family_indices_unrotated",
                "physical_proton_bound_status": "not_evaluable_yet",
                "remaining_requirements": [
                    "assemble physical hermitian vector multiplets and Clebsch phases",
                    "map crossed conjugate fields to the chosen external-state convention",
                    "project to a canonical all-left Wilson basis if applicable",
                    "apply physical flavor rotations",
                    "insert numerical masses and couplings",
                    "run Wilson coefficients to the hadronic scale",
                    "contract with hadronic matrix elements and channel limits",
                ],
            }
        )
    return out


def summarize(rows: list[dict], p14: dict) -> dict:
    status_counts = Counter(row["all_left_status"] for row in rows)
    basis_counts = Counter(
        row["canonical_bnv_basis_id"] or "none" for row in rows
    )
    gate_counts = Counter(row["p14_matching_gate"] for row in rows)
    canonical = sum(
        1 for row in rows if row["canonical_bnv_basis_id"] is not None
    )
    return {
        "p14_current_pair_rows_consumed": p14["summary"]["current_pair_rows"],
        "p15_rows_generated": len(rows),
        "all_left_status_counts": dict(status_counts),
        "canonical_bnv_basis_counts": dict(basis_counts),
        "p14_gate_counts_replayed": dict(gate_counts),
        "canonical_bnv_rows_evaluable_now": canonical,
        "physical_proton_bounds_evaluable_now": 0,
        "all_rows_have_crossed_fields": all(
            len(row["crossed_all_left_fields"]) == 4 for row in rows
        ),
        "all_rows_have_charge_audit": all(
            set(row["global_charges"]) == {"B", "L", "B_minus_L"}
            for row in rows
        ),
        "branch_v_all_left_gate_passes": len(rows)
        == p14["summary"]["current_pair_rows"],
    }


def build_payload() -> dict:
    p14 = load_json(P14_JSON)
    rows = build_rows(p14)
    summary = summarize(rows, p14)
    return {
        "description": "Route-C P15-V vector-current to all-left Fierz/flavor gate",
        "boundary": (
            "P15-V fixes a two-component Fierz/crossing convention for Branch V "
            "and audits the crossed all-left skeletons.  It does not compute "
            "physical proton lifetimes, and it does not override the scalar/"
            "source Branch S."
        ),
        "fixed_convention": {
            "input_current_product": (
                "(psi_t1^dagger bar_sigma^mu psi_s1)"
                "(psi_s2^dagger bar_sigma_mu psi_t2)"
            ),
            "fierz_identity": (
                "(chi^dagger bar_sigma^mu psi)(eta^dagger bar_sigma_mu xi)"
                " = 2 (chi^dagger eta^dagger)(psi xi)"
            ),
            "crossed_all_left_skeleton": (
                "(psi_s1 psi_t2)(bar(psi_t1) bar(psi_s2))"
            ),
            "coefficient_prefactor": "-2 g_X^2 c1 c2^*/M_X^2",
        },
        "summary": summary,
        "all_left_rows": rows,
        "interpretation": {
            "why_no_proton_limit_yet": (
                "The gate fixes the Fierz convention and global-charge audit, "
                "but physical proton limits require a canonical external-state "
                "crossing convention, physical flavor rotations, RG factors, "
                "hadronic matrix elements, and channel limits."
            ),
            "branch_v_status": (
                "Branch V is now in an all-left crossed-skeleton bookkeeping "
                "basis.  Under the current sparse-map pairing, no row is yet a "
                "canonical BNV proton operator."
            ),
            "branch_s_status": (
                "P13-S remains a separate scalar/source Wilson-basis branch and "
                "must not be discarded."
            ),
        },
        "next_stage": {
            "recommended": "P16_choose_vector_physical_multiplet_or_return_to_scalar_source",
            "required_decision": (
                "Either assemble physical hermitian vector multiplets and their "
                "cross-root Clebsch phases to search for canonical BNV vector "
                "operators, or return to the P13-S scalar/source recoupling path. "
                "Do not insert proton limits yet."
            ),
        },
    }


def markdown_table(rows: list[list[str]]) -> list[str]:
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

    sample_rows = [["row", "mass", "fields", "status"]]
    for row in payload["all_left_rows"][:18]:
        sample_rows.append(
            [
                str(row["row_id"]),
                f"`{row['mass_block']}`",
                "`" + " ".join(row["crossed_all_left_fields"]) + "`",
                f"`{row['all_left_status']}`",
            ]
        )

    lines = [
        "# Route-C P15-V Vector Current to All-Left Fierz Gate",
        "",
        payload["boundary"],
        "",
        "## Fixed Convention",
        "",
        f"- input: `{payload['fixed_convention']['input_current_product']}`",
        f"- Fierz: `{payload['fixed_convention']['fierz_identity']}`",
        f"- crossed skeleton: `{payload['fixed_convention']['crossed_all_left_skeleton']}`",
        f"- coefficient prefactor: `{payload['fixed_convention']['coefficient_prefactor']}`",
        "",
        "## Summary",
        "",
        *markdown_table(summary_rows),
        "",
        "## Sample Crossed Rows",
        "",
        *markdown_table(sample_rows),
        "",
        "## Interpretation",
        "",
        f"- Branch V: {payload['interpretation']['branch_v_status']}",
        f"- Branch S: {payload['interpretation']['branch_s_status']}",
        f"- Proton limits: {payload['interpretation']['why_no_proton_limit_yet']}",
        "",
        "## Next Stage",
        "",
        f"`{payload['next_stage']['recommended']}`",
        "",
        payload["next_stage"]["required_decision"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "p15_rows_generated": payload["summary"]["p15_rows_generated"],
                "all_left_status_counts": payload["summary"][
                    "all_left_status_counts"
                ],
                "canonical_bnv_rows_evaluable_now": payload["summary"][
                    "canonical_bnv_rows_evaluable_now"
                ],
                "physical_proton_bounds_evaluable_now": payload["summary"][
                    "physical_proton_bounds_evaluable_now"
                ],
                "branch_v_all_left_gate_passes": payload["summary"][
                    "branch_v_all_left_gate_passes"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
