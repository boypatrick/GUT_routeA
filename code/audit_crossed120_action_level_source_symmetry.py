#!/usr/bin/env python3
"""Audit a post-Spin(10) source symmetry for the crossed-120 closure.

No web lookup is used.

The local flavor+d=5 branch is now closed only if the crossed 120_A/120_B
triplet projector is accepted.  A field-only unbroken Spin(10) grading cannot
enforce that projector.  This audit therefore asks a more honest question:

    Can a post-Spin(10) Pati-Salam constrained-source sector enforce the
    crossed triplet pairing while keeping the doublet-mixing refit and the
    threshold-silent completion compatible?

The answer is conditional.  A fragment/source grading plus an R-selection rule
can allow the required terms and forbid the dangerous ones in the spurion
ledger, but this remains a post-breaking constrained-source construction; it
does not become a microscopic first-principles Spin(10) derivation.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "crossed120_action_level_source_symmetry"

CKM_CROSSED = ROOT / "output" / "source_consistent_ckm_crossed120" / "summary.json"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
LINK_LOCK = ROOT / "output" / "crossed_120_link_locking_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"
CROSSED = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def mod_charge(value: int, modulus: int) -> int:
    return int(value) % int(modulus)


def add_charges(items: list[dict[str, Any]], fields: dict[str, dict[str, Any]]) -> dict[str, int]:
    total = {"R": 0, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0}
    for item in items:
        charge = fields[item["name"]]
        power = int(item.get("power", 1))
        total["R"] += power * int(charge["R"])
        total["Z4_src"] += power * int(charge["Z4_src"])
        total["Z2_frag"] += power * int(charge["Z2_frag"])
        total["Z2_inert"] += power * int(charge["Z2_inert"])
    return {
        "R": total["R"],
        "Z4_src": mod_charge(total["Z4_src"], 4),
        "Z2_frag": mod_charge(total["Z2_frag"], 2),
        "Z2_inert": mod_charge(total["Z2_inert"], 2),
    }


def charge_allowed(total: dict[str, int]) -> bool:
    return total["R"] == 2 and total["Z4_src"] == 0 and total["Z2_frag"] == 0 and total["Z2_inert"] == 0


def field_table() -> dict[str, dict[str, Any]]:
    """A minimal post-breaking source/spurion grading.

    R forbids bare mass terms.  Z4_src distinguishes A/B source copies once
    the theory is in the PS fragment language.  Z2_frag distinguishes triplet
    source operators from doublet operators, so off-block triplet-doublet
    mixing is not generic.  Z2_inert decouples the inert doublet partners.
    """

    return {
        "B_QQ": {"R": 2, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "10 10 matter bilinear"},
        "B_QL": {"R": 2, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "10 5bar matter bilinear"},
        "B_16_16": {"R": 2, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "generic doublet Yukawa bilinear"},
        "T_A": {"R": 0, "Z4_src": 0, "Z2_frag": 1, "Z2_inert": 0, "role": "120_A triplet source"},
        "T_B": {"R": 0, "Z4_src": 1, "Z2_frag": 1, "Z2_inert": 0, "role": "120_B triplet source"},
        "Tbar_A": {"R": 0, "Z4_src": 0, "Z2_frag": 1, "Z2_inert": 0, "role": "conjugate/source partner of T_A"},
        "Tbar_B": {"R": 0, "Z4_src": 3, "Z2_frag": 1, "Z2_inert": 0, "role": "conjugate/source partner of T_B"},
        "D_A": {"R": 0, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "physical doublet component from 120_A"},
        "D_B": {"R": 0, "Z4_src": 1, "Z2_frag": 0, "Z2_inert": 0, "role": "physical doublet component from 120_B"},
        "L": {"R": 0, "Z4_src": 2, "Z2_frag": 0, "Z2_inert": 1, "role": "inert doublet partner"},
        "Lbar": {"R": 0, "Z4_src": 2, "Z2_frag": 0, "Z2_inert": 1, "role": "inert doublet conjugate"},
        "Y_QQ_A": {"R": 0, "Z4_src": 0, "Z2_frag": 1, "Z2_inert": 0, "role": "spurion allowing B_QQ T_A"},
        "Y_QL_B": {"R": 0, "Z4_src": 3, "Z2_frag": 1, "Z2_inert": 0, "role": "spurion allowing B_QL T_B"},
        "Y_D_A": {"R": 0, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "doublet Yukawa spurion for D_A"},
        "Y_D_B": {"R": 0, "Z4_src": 3, "Z2_frag": 0, "Z2_inert": 0, "role": "doublet Yukawa spurion for D_B"},
        "S_AB": {"R": 2, "Z4_src": 3, "Z2_frag": 0, "Z2_inert": 0, "role": "crossed mass spurion Tbar_A T_B"},
        "S_BA": {"R": 2, "Z4_src": 1, "Z2_frag": 0, "Z2_inert": 0, "role": "finite-lift crossed mass spurion Tbar_B T_A"},
        "S_lock": {"R": 2, "Z4_src": 0, "Z2_frag": 0, "Z2_inert": 0, "role": "common inert/doublet partner mass spurion"},
    }


def operator_terms() -> list[dict[str, Any]]:
    """Operator ledger.

    ``spurion_present`` distinguishes a symmetry-allowed spurion that is part
    of the constrained source sector from a spurion that would be allowed by
    charges but is not included.  This is the honest post-breaking source
    assumption replacing the forbidden unbroken Spin(10) field-only grading.
    """

    return [
        {
            "operator": "Y_QQ_A B_QQ T_A",
            "items": [{"name": "Y_QQ_A"}, {"name": "B_QQ"}, {"name": "T_A"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "QQ triplet source uses G_A-like 120_A direction.",
        },
        {
            "operator": "Y_QL_B B_QL T_B",
            "items": [{"name": "Y_QL_B"}, {"name": "B_QL"}, {"name": "T_B"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "QL triplet source uses crossed G_B-like 120_B direction.",
        },
        {
            "operator": "Y_D_A B_16_16 D_A",
            "items": [{"name": "Y_D_A"}, {"name": "B_16_16"}, {"name": "D_A"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "Doublet sector must retain 120_A mixing freedom.",
        },
        {
            "operator": "Y_D_B B_16_16 D_B",
            "items": [{"name": "Y_D_B"}, {"name": "B_16_16"}, {"name": "D_B"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "Doublet sector must retain 120_B mixing freedom used by the CKM refit.",
        },
        {
            "operator": "S_AB Tbar_A T_B",
            "items": [{"name": "S_AB"}, {"name": "Tbar_A"}, {"name": "T_B"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "Main crossed triplet inverse block.",
        },
        {
            "operator": "S_BA Tbar_B T_A",
            "items": [{"name": "S_BA"}, {"name": "Tbar_B"}, {"name": "T_A"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "Finite-lift epsilon block; can be suppressed by kappa.",
        },
        {
            "operator": "S_lock L Lbar",
            "items": [{"name": "S_lock"}, {"name": "L"}, {"name": "Lbar"}],
            "desired": "allowed",
            "spurion_present": True,
            "reason": "Inert doublet partners complete the threshold vector.",
        },
        {
            "operator": "bare Tbar_A T_A",
            "items": [{"name": "Tbar_A"}, {"name": "T_A"}],
            "desired": "forbidden",
            "spurion_present": False,
            "reason": "R-symmetry forbids bare mass; no same-source R=2 spurion is included.",
        },
        {
            "operator": "bare Tbar_B T_B",
            "items": [{"name": "Tbar_B"}, {"name": "T_B"}],
            "desired": "forbidden",
            "spurion_present": False,
            "reason": "R-symmetry forbids bare mass; no same-source R=2 spurion is included.",
        },
        {
            "operator": "wrong Y_QQ_B B_QQ T_B",
            "items": [{"name": "B_QQ"}, {"name": "T_B"}],
            "desired": "forbidden",
            "spurion_present": False,
            "reason": "No QQ-to-B triplet Yukawa spurion is present.",
        },
        {
            "operator": "wrong Y_QL_A B_QL T_A",
            "items": [{"name": "B_QL"}, {"name": "T_A"}],
            "desired": "forbidden",
            "spurion_present": False,
            "reason": "No QL-to-A triplet Yukawa spurion is present.",
        },
        {
            "operator": "B_16_16 L",
            "items": [{"name": "B_16_16"}, {"name": "L"}],
            "desired": "forbidden",
            "spurion_present": False,
            "reason": "Z2_inert forbids inert doublet Yukawa leakage.",
        },
        {
            "operator": "S_lock D_A Lbar",
            "items": [{"name": "S_lock"}, {"name": "D_A"}, {"name": "Lbar"}],
            "desired": "forbidden",
            "spurion_present": True,
            "reason": "Z2_inert forbids physical/inert doublet mixing.",
        },
        {
            "operator": "S_lock T_A D_A",
            "items": [{"name": "S_lock"}, {"name": "T_A"}, {"name": "D_A"}],
            "desired": "forbidden",
            "spurion_present": True,
            "reason": "Z2_frag forbids triplet-doublet off-block mixing.",
        },
    ]


def evaluate_operator(row: dict[str, Any], fields: dict[str, dict[str, Any]]) -> dict[str, Any]:
    total = add_charges(row["items"], fields)
    neutral = charge_allowed(total)
    actual_allowed = bool(neutral and row["spurion_present"])
    desired_allowed = row["desired"] == "allowed"
    return {
        "operator": row["operator"],
        "desired": row["desired"],
        "actual": "allowed" if actual_allowed else "forbidden",
        "passes": actual_allowed == desired_allowed,
        "charge_R": total["R"],
        "charge_Z4_src": total["Z4_src"],
        "charge_Z2_frag": total["Z2_frag"],
        "charge_Z2_inert": total["Z2_inert"],
        "charge_neutral_as_superpotential_term": neutral,
        "spurion_present": bool(row["spurion_present"]),
        "reason": row["reason"],
    }


def best_ckm_row(payload: dict[str, Any]) -> dict[str, Any]:
    return payload["verdict"]["best_mass_valid_ckm"]


def build() -> dict[str, Any]:
    ckm = read_json(CKM_CROSSED)
    ps = read_json(PS_SOURCE)
    link = read_json(LINK_LOCK)
    completed = read_json(COMPLETED)
    crossed = read_json(CROSSED)

    fields = field_table()
    operator_rows = [evaluate_operator(row, fields) for row in operator_terms()]
    allowed_ok = all(row["passes"] for row in operator_rows)
    b = best_ckm_row(ckm)

    selection = {
        "field_only_spin10_no_go_reused": not bool(crossed["verdict"]["field_only_unbroken_spin10_projector_possible"]),
        "operator_ledger_passes": allowed_ok,
        "triplet_projector_realized_by_ps_source_action": bool(
            ps["verdict"]["ps_component_action_realizes_projector"]
        ),
        "constrained_source_threshold_silent": bool(ps["verdict"]["constrained_source_threshold_silent"]),
        "complete_partner_threshold_safe": bool(completed["verdict"]["exact_completed_pair_threshold_silent"]),
        "inert_grading_decouples_doublets": bool(completed["verdict"]["inert_grading_decouples_doublets"]),
        "link_holomorphic_constraints_derive_block": bool(
            link["verdict"]["holomorphic_constraints_realize_inverse_block"]
        ),
        "nlsm_or_dterm_unitarity_required": bool(link["verdict"]["nlsm_unitarity_required"]),
        "local_ckm_d5_closure_survives": bool(ckm["verdict"]["local_conditional_closure_complete"]),
    }
    selection["post_spin10_source_symmetry_branch_viable"] = bool(
        selection["field_only_spin10_no_go_reused"]
        and selection["operator_ledger_passes"]
        and selection["triplet_projector_realized_by_ps_source_action"]
        and selection["local_ckm_d5_closure_survives"]
        and (
            selection["constrained_source_threshold_silent"]
            or selection["complete_partner_threshold_safe"]
        )
    )

    verdict = {
        **selection,
        "best_ckm_score": float(b["ckm_score"]),
        "best_mass_score": float(b["mass_score"]),
        "best_d5_margin": float(b["future_margin"]),
        "best_seesaw_residual": float(b["seesaw_residual"]),
        "completed_partner_percent_lock_window": float(completed["mass_locking_window"]["percent_window"]),
        "holomorphic_link_nullity_real": int(link["verdict"]["constraint_nullity_real"]),
        "publication_complete": False,
        "interpretation": (
            "A post-Spin(10) PS fragment/source grading can consistently allow the crossed "
            "120_A/120_B triplet source, the finite-lift reverse entry, the physical doublet "
            "Yukawa/mixing sector used by the CKM refit, and inert threshold-completing "
            "doublets, while forbidding wrong triplet Yukawas, same-source bare triplet "
            "masses, inert Yukawas, and triplet-doublet off-block mixing.  This upgrades the "
            "crossed projector from a naked imposed matrix to a conditional constrained-source "
            "selection rule.  It is still not an unconditional Spin(10) derivation because the "
            "selection is post-breaking/spurion-assisted and the holomorphic link requires an "
            "NLSM/D-term or composite unitarity constraint."
        ),
    }

    return {
        "note": "No web lookup used. Action-level source-symmetry audit for the crossed 120 branch.",
        "symmetry": {
            "factors": {
                "R": "superpotential R-charge; bare masses have R=0 and are forbidden unless an R=2 spurion is present",
                "Z4_src": "post-Spin(10) source-copy grading distinguishing 120_A and 120_B fragments",
                "Z2_frag": "post-breaking triplet/doublet fragment grading suppressing off-block mixing",
                "Z2_inert": "inert partner parity",
            },
            "warning": "This is not a field-only unbroken Spin(10) charge table.",
        },
        "charge_table": fields,
        "operator_ledger": operator_rows,
        "input_summaries": {
            "source_consistent_ckm_crossed120": str(CKM_CROSSED.relative_to(ROOT)),
            "ps_crossed_120_source_action": str(PS_SOURCE.relative_to(ROOT)),
            "crossed_120_link_locking_action": str(LINK_LOCK.relative_to(ROOT)),
            "completed_120_partner_action": str(COMPLETED.relative_to(ROOT)),
            "crossed_120_triplet_projector": str(CROSSED.relative_to(ROOT)),
        },
        "numerical_bridge": {
            "best_mass_valid_ckm": b,
            "ps_source_flatness": ps["flatness_and_hessian"],
            "threshold_interpretation": ps["verdict"],
            "completed_partner_window": completed["mass_locking_window"],
            "link_constraint_geometry": link["linearized_constraint_geometry"],
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
    lines = [
        "# Crossed 120 action-level source symmetry audit",
        "",
        "No web lookup was used.",
        "",
        "This is a post-Spin(10) PS fragment/source selection rule, not a field-only unbroken Spin(10) grading.",
        "",
        "## Operator Ledger",
        "",
        "| operator | desired | actual | R | Z4 | Z2_frag | Z2_inert | pass |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in payload["operator_ledger"]:
        lines.append(
            f"| `{row['operator']}` | {row['desired']} | {row['actual']} | "
            f"{row['charge_R']} | {row['charge_Z4_src']} | {row['charge_Z2_frag']} | "
            f"{row['charge_Z2_inert']} | {row['passes']} |"
        )
    lines.extend(
        [
            "",
            "## Numerical Bridge",
            "",
            f"best CKM score: `{v['best_ckm_score']:.6e}`",
            f"best mass score: `{v['best_mass_score']:.6e}`",
            f"crossed d5 margin: `{v['best_d5_margin']:.6e}`",
            f"seesaw residual: `{v['best_seesaw_residual']:.3e}`",
            f"completed partner percent lock window: `{v['completed_partner_percent_lock_window']:.6e}`",
            f"holomorphic link nullity: `{v['holomorphic_link_nullity_real']}`",
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
    field_rows = [
        {"name": name, **{key: val for key, val in data.items() if key != "role"}, "role": data["role"]}
        for name, data in payload["charge_table"].items()
    ]
    write_csv(
        OUT / "charge_table.csv",
        field_rows,
        ["name", "R", "Z4_src", "Z2_frag", "Z2_inert", "role"],
    )
    write_csv(
        OUT / "operator_ledger.csv",
        payload["operator_ledger"],
        [
            "operator",
            "desired",
            "actual",
            "passes",
            "charge_R",
            "charge_Z4_src",
            "charge_Z2_frag",
            "charge_Z2_inert",
            "charge_neutral_as_superpotential_term",
            "spurion_present",
            "reason",
        ],
    )
    write_report(payload)

    v = payload["verdict"]
    print("Crossed 120 action-level source symmetry audit")
    print(f"  operator ledger passes: {v['operator_ledger_passes']}")
    print(f"  local CKM+d5 closure survives: {v['local_ckm_d5_closure_survives']}")
    print(f"  post-Spin10 source branch viable: {v['post_spin10_source_symmetry_branch_viable']}")
    print(f"  NLSM/D-term unitarity required: {v['nlsm_or_dterm_unitarity_required']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
