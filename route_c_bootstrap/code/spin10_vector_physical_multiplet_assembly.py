#!/usr/bin/env python3
"""Route-C P16 physical Hermitian vector multiplet assembly.

P16 goes beyond the P15 same-map J_X J_X^\dagger bookkeeping gate.  A physical
massive Hermitian vector is assembled from charge-conjugate root maps.  At the
symbolic level this permits cross-root products

  J_q^\mu J_{-q,\mu},

with an undetermined Clebsch/phase convention.  After the same two-component
Fierz identity used in P15, such cross-root rows can generate the familiar
vector-mediated BNV skeletons

  Q Q \bar{u^c}\bar{e^c},   Q Q \bar{d^c}\bar{\nu^c}

and their conjugates.

This is still a matching gate, not a proton-lifetime calculation.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P9_JSON = OUT / "spin10_action_completion_ledger.json"
P10_JSON = OUT / "spin10_breaking_mass_matrix.json"
P15_JSON = OUT / "spin10_vector_all_left_fierz_gate.json"
OUT_JSON = OUT / "spin10_vector_physical_multiplet_assembly.json"
OUT_MD = OUT / "spin10_vector_physical_multiplet_assembly.md"


MASS_BY_FAMILY = {
    "Pati-Salam SU(4)_C leptoquark": "M_LQ",
    "broken SU(2)_R charged current": "M_WR",
    "Spin(10)/Pati-Salam off-face generator": "M_(6,2,2)",
}

GLOBAL_CHARGES = {
    "Q": {"B": Fraction(1, 3), "L": Fraction(0)},
    "u^c": {"B": Fraction(-1, 3), "L": Fraction(0)},
    "d^c": {"B": Fraction(-1, 3), "L": Fraction(0)},
    "L": {"B": Fraction(0), "L": Fraction(1)},
    "e^c": {"B": Fraction(0), "L": Fraction(-1)},
    "nu^c": {"B": Fraction(0), "L": Fraction(-1)},
}

VECTOR_BNV_BASIS = {
    ("Q", "Q", "bar(e^c)", "bar(u^c)"): "O_V_QQ_UbarEbar",
    ("Q", "Q", "bar(d^c)", "bar(nu^c)"): "O_V_QQ_DbarNbar",
    ("e^c", "u^c", "bar(Q)", "bar(Q)"): "O_V_QQ_UbarEbar_dagger",
    ("d^c", "nu^c", "bar(Q)", "bar(Q)"): "O_V_QQ_DbarNbar_dagger",
}


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def qfrac(text: str) -> Fraction:
    return Fraction(text)


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def neg_key(key: tuple[str, str]) -> tuple[str, str]:
    return (frac(-qfrac(key[0])), frac(-qfrac(key[1])))


def map_charge(vmap: dict) -> tuple[str, str]:
    charges = {
        (entry["delta_Y"]["fraction"], entry["delta_B_minus_L"]["fraction"])
        for entry in vmap["sparse_transition_entries"]
    }
    if len(charges) != 1:
        raise ValueError(f"map {vmap['label']} has nonuniform vector charges: {charges}")
    return next(iter(charges))


def display_field(field: str, conjugate: bool = False) -> str:
    return f"bar({field})" if conjugate else field


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


def field_key(fields: list[dict]) -> tuple[str, ...]:
    return tuple(sorted(display_field(f["multiplet"], f["conjugate"]) for f in fields))


def basis_id(fields: list[dict]) -> str | None:
    return VECTOR_BNV_BASIS.get(field_key(fields))


def classify(fields: list[dict]) -> str:
    bid = basis_id(fields)
    q = charge_sum(fields)
    if bid and not bid.endswith("_dagger"):
        return "canonical_vector_BNV_basis_candidate"
    if bid and bid.endswith("_dagger"):
        return "canonical_vector_BNV_conjugate_candidate"
    if q["B"] == 0 and q["L"] == 0:
        return "B_and_L_conserving_physical_pair"
    if q["B_minus_L"] == 0 and abs(q["B"]) == 1:
        return "BminusL_preserving_BNV_noncanonical_pair"
    if q["B_minus_L"] == 0:
        return "BminusL_preserving_noncanonical_pair"
    return "nonstandard_BminusL_pair"


def crossed_fields(e_q: dict, e_mq: dict) -> list[dict]:
    # Product J_q J_{-q}; after Fierz:
    # (target_q^\dagger source_q)(target_-q^\dagger source_-q)
    # -> (source_q source_-q)(bar(target_q) bar(target_-q)).
    return [
        {
            "multiplet": e_q["source_multiplet"],
            "label": e_q["source_label"],
            "conjugate": False,
        },
        {
            "multiplet": e_mq["source_multiplet"],
            "label": e_mq["source_label"],
            "conjugate": False,
        },
        {
            "multiplet": e_q["target_multiplet"],
            "label": e_q["target_label"],
            "conjugate": True,
        },
        {
            "multiplet": e_mq["target_multiplet"],
            "label": e_mq["target_label"],
            "conjugate": True,
        },
    ]


def phase_symbol(label_a: str, label_b: str) -> str:
    safe = (
        label_a.replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
        .replace(",", "_")
        .replace("^", "")
        .replace("+", "p")
        .replace("-", "m")
    )
    safe_b = (
        label_b.replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
        .replace(",", "_")
        .replace("^", "")
        .replace("+", "p")
        .replace("-", "m")
    )
    return f"xi_{safe}__{safe_b}"


def build_physical_pairs(p9: dict, p10: dict) -> list[dict]:
    maps = p9["broken_generator_maps"]["maps"]
    by_charge: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for vmap in maps:
        charge_key = map_charge(vmap)
        mass_block = MASS_BY_FAMILY[vmap["generator_family"]]
        by_charge[(vmap["generator_family"], mass_block, *charge_key)].append(vmap)

    mass_replacements = p10["p6_p7_replay_interface"]["replace_symbolic_masses"]
    rows: list[dict] = []
    seen_pairs: set[tuple[str, str]] = set()
    row_id = 0
    for (family, mass_block, dy, dbl), lhs_maps in sorted(by_charge.items()):
        opposite = (family, mass_block, *neg_key((dy, dbl)))
        rhs_maps = by_charge.get(opposite, [])
        for map_q in lhs_maps:
            for map_mq in rhs_maps:
                pair_key = tuple(sorted((map_q["label"], map_mq["label"])))
                # Keep one orientation per Hermitian pair.  The opposite row is
                # the complex conjugate and is represented by the basis status.
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                phase = phase_symbol(map_q["label"], map_mq["label"])
                for e_q, e_mq in itertools.product(
                    map_q["sparse_transition_entries"],
                    map_mq["sparse_transition_entries"],
                ):
                    fields = crossed_fields(e_q, e_mq)
                    q = charge_sum(fields)
                    bid = basis_id(fields)
                    rows.append(
                        {
                            "row_id": row_id,
                            "vector_map_q": map_q["label"],
                            "vector_map_minus_q": map_mq["label"],
                            "generator_family": family,
                            "mass_block": mass_block,
                            "mass_expression": mass_replacements[mass_block],
                            "vector_charge_q": {
                                "delta_Y": dy,
                                "delta_B_minus_L": dbl,
                            },
                            "cross_root_phase": phase,
                            "transition_q": {
                                "source": e_q["source_label"],
                                "target": e_q["target_label"],
                                "signature": f"{e_q['source_multiplet']}->{e_q['target_multiplet']}",
                                "coefficient": e_q["coefficient"],
                            },
                            "transition_minus_q": {
                                "source": e_mq["source_label"],
                                "target": e_mq["target_label"],
                                "signature": f"{e_mq['source_multiplet']}->{e_mq['target_multiplet']}",
                                "coefficient": e_mq["coefficient"],
                            },
                            "crossed_all_left_fields": [
                                display_field(f["multiplet"], f["conjugate"])
                                for f in fields
                            ],
                            "coefficient_template": (
                                f"C_phys[row_{row_id}] = -2 g_X^2 "
                                f"{phase} c_q c_-q/{mass_block}^2"
                            ),
                            "global_charges": {
                                "B": frac(q["B"]),
                                "L": frac(q["L"]),
                                "B_minus_L": frac(q["B_minus_L"]),
                            },
                            "canonical_vector_bnv_basis_id": bid,
                            "physical_pair_status": classify(fields),
                            "proton_bound_status": "not_evaluable_yet",
                            "remaining_requirements": [
                                "fix physical Hermitian generator normalization",
                                "fix cross-root Clebsch phases",
                                "project color and weak indices to canonical vector BNV basis",
                                "apply physical flavor rotations",
                                "insert numerical masses and couplings",
                                "run Wilson coefficients to hadronic scale",
                                "contract with hadronic matrix elements and experimental limits",
                            ],
                        }
                    )
                    row_id += 1
    return rows


def summarize(rows: list[dict], p15: dict) -> dict:
    status_counts = Counter(row["physical_pair_status"] for row in rows)
    basis_counts = Counter(
        row["canonical_vector_bnv_basis_id"] or "none" for row in rows
    )
    family_counts = Counter(row["generator_family"] for row in rows)
    bnv = [
        row
        for row in rows
        if row["physical_pair_status"]
        in {
            "canonical_vector_BNV_basis_candidate",
            "canonical_vector_BNV_conjugate_candidate",
            "BminusL_preserving_BNV_noncanonical_pair",
        }
    ]
    canonical = [
        row
        for row in rows
        if row["physical_pair_status"]
        in {
            "canonical_vector_BNV_basis_candidate",
            "canonical_vector_BNV_conjugate_candidate",
        }
    ]
    return {
        "p15_same_map_rows_were_canonical_bnv": p15["summary"][
            "canonical_bnv_rows_evaluable_now"
        ],
        "physical_cross_root_rows": len(rows),
        "physical_cross_root_rows_by_family": dict(family_counts),
        "physical_pair_status_counts": dict(status_counts),
        "canonical_vector_bnv_basis_counts": dict(basis_counts),
        "bnv_candidate_rows_before_flavor": len(bnv),
        "canonical_vector_bnv_candidate_rows_before_flavor": len(canonical),
        "physical_proton_bounds_evaluable_now": 0,
        "all_rows_have_cross_root_phase": all(
            bool(row["cross_root_phase"]) for row in rows
        ),
        "all_rows_have_charge_audit": all(
            set(row["global_charges"]) == {"B", "L", "B_minus_L"}
            for row in rows
        ),
        "branch_v_physical_multiplet_gate_passes": len(rows) > 0
        and len(canonical) > 0
        and p15["summary"]["branch_v_all_left_gate_passes"],
    }


def build_payload() -> dict:
    p9 = load_json(P9_JSON)
    p10 = load_json(P10_JSON)
    p15 = load_json(P15_JSON)
    rows = build_physical_pairs(p9, p10)
    summary = summarize(rows, p15)
    return {
        "description": "Route-C P16 physical Hermitian vector multiplet assembly",
        "boundary": (
            "P16 assembles charge-conjugate root maps into symbolic physical "
            "Hermitian vector multiplets and audits whether canonical vector "
            "BNV skeletons can appear.  It does not fix numerical Clebsch "
            "phases, flavor rotations, RG factors, hadronic matrix elements, "
            "or proton lifetime limits."
        ),
        "physical_vector_rule": {
            "root_pair": "E_q + E_-q with symbolic cross-root phase xi",
            "effective_product": "J_q^mu J_-q,mu",
            "fierz_prefactor": "-2 g_X^2 xi c_q c_-q / M_X^2",
            "canonical_vector_bnv_basis": sorted(set(VECTOR_BNV_BASIS.values())),
        },
        "summary": summary,
        "physical_pair_rows": rows,
        "interpretation": {
            "what_changed_from_p15": (
                "P15 same-map J J^dagger rows were all B/L conserving.  P16 "
                "allows charge-conjugate root-map products J_q J_-q, which can "
                "produce vector-mediated BNV skeletons."
            ),
            "why_no_proton_limit_yet": (
                "The gate still has symbolic cross-root phases, no physical "
                "flavor rotations, no RG evolution, and no hadronic/channel "
                "inputs."
            ),
            "branch_s_status": (
                "The scalar/source Branch S remains available and is not "
                "discarded by the vector result."
            ),
        },
        "next_stage": {
            "recommended": "P17_vector_clebsch_projection_and_flavor_interface",
            "required_decision": (
                "Fix a physical Hermitian generator normalization and cross-root "
                "Clebsch phase convention for the canonical vector BNV rows, then "
                "export a flavor-rotation interface.  Do not insert proton limits "
                "until RG and hadronic inputs are supplied."
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

    candidate_rows = []
    seen_candidates: set[tuple[str, str, tuple[str, ...], str]] = set()
    for row in payload["physical_pair_rows"]:
        if row["canonical_vector_bnv_basis_id"] is None:
            continue
        key = (
            row["canonical_vector_bnv_basis_id"],
            row["mass_block"],
            tuple(row["crossed_all_left_fields"]),
            row["cross_root_phase"],
        )
        if key in seen_candidates:
            continue
        seen_candidates.add(key)
        candidate_rows.append(row)
        if len(candidate_rows) == 24:
            break
    cand_table = [["basis", "mass", "fields", "phase"]]
    for row in candidate_rows:
        cand_table.append(
            [
                f"`{row['canonical_vector_bnv_basis_id']}`",
                f"`{row['mass_block']}`",
                "`" + " ".join(row["crossed_all_left_fields"]) + "`",
                f"`{row['cross_root_phase']}`",
            ]
        )

    lines = [
        "# Route-C P16 Physical Hermitian Vector Multiplet Assembly",
        "",
        payload["boundary"],
        "",
        "## Physical Vector Rule",
        "",
        f"- root pair: `{payload['physical_vector_rule']['root_pair']}`",
        f"- product: `{payload['physical_vector_rule']['effective_product']}`",
        f"- prefactor: `{payload['physical_vector_rule']['fierz_prefactor']}`",
        "",
        "## Summary",
        "",
        *markdown_table(summary_rows),
        "",
        "## Sample Canonical Vector BNV Candidates",
        "",
        *markdown_table(cand_table),
        "",
        "## Interpretation",
        "",
        f"- From P15: {payload['interpretation']['what_changed_from_p15']}",
        f"- Proton limits: {payload['interpretation']['why_no_proton_limit_yet']}",
        f"- Branch S: {payload['interpretation']['branch_s_status']}",
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
                "physical_cross_root_rows": payload["summary"][
                    "physical_cross_root_rows"
                ],
                "physical_pair_status_counts": payload["summary"][
                    "physical_pair_status_counts"
                ],
                "canonical_vector_bnv_candidate_rows_before_flavor": payload[
                    "summary"
                ]["canonical_vector_bnv_candidate_rows_before_flavor"],
                "physical_proton_bounds_evaluable_now": payload["summary"][
                    "physical_proton_bounds_evaluable_now"
                ],
                "branch_v_physical_multiplet_gate_passes": payload["summary"][
                    "branch_v_physical_multiplet_gate_passes"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
