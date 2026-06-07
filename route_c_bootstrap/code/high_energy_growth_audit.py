#!/usr/bin/env python3
"""Route-C P6 controlled high-energy-growth audit.

P6 consumes the P3 residue ledger and P4 Ward-identity ledger.  It does not
compute complete amplitudes.  It classifies which symbolic branches have an
uncontrolled high-energy-growth risk if interpreted as isolated broken vector
exchange, which branches can remain as scalar/auxiliary exchanges at this
bookkeeping level, and when a tower-like UV completion should be treated as a
serious comparison point.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P2_JSON = OUT / "four_point_pole_ansatz.json"
P3_JSON = OUT / "residue_positivity.json"
P4_JSON = OUT / "ward_identity_ledger.json"
P5_JSON = OUT / "anomaly_chiral_consistency.json"


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def p4_by_mediator(p4: dict) -> dict[str, dict]:
    return {row["mediator"]: row for row in p4["mediator_ward_requirements"]}


def p3_by_mediator(p3: dict) -> dict[str, dict]:
    return {row["mediator"]["label"]: row for row in p3["mediator_blocks"]}


def diagnostic_seed_counts(p2: dict) -> dict[str, int]:
    seeds = p2.get("diagnostic_seeds", {})
    return {name: len(rows) for name, rows in seeds.items()}


def vector_growth_row(block: dict, ward: dict) -> dict:
    label = block["mediator"]["label"]
    requires_completion = ward["requires_goldstone_higgs_completion_if_vector"]
    if requires_completion:
        verdict = "finite_vector_incomplete"
        naive_growth = "uncancelled longitudinal piece can scale as s/M_X^2 or worse"
        required_completion = [
            "explicit broken generator and current normalization",
            "broken-vector mass relation",
            "Goldstone amplitude phi_X",
            "Higgs/source sector enforcing p_mu M^mu(V_X)=M_X M(phi_X)",
            "spin numerator and contact terms fixed by the same symmetry",
        ]
        tower_hint = (
            "not forced by P6 alone; becomes relevant if no finite "
            "Goldstone/Higgs/source completion cancels the growth"
        )
    else:
        verdict = "unbroken_vector_bookkeeping_ready"
        naive_growth = "ordinary massless Ward identity should remove p_mu terms"
        required_completion = [
            "component-level unbroken generator matrices",
            "ordinary gauge-invariant numerator check",
        ]
        tower_hint = "not indicated"

    return {
        "mediator": label,
        "channel_count": block["channel_count"],
        "vector_interpretation": ward["vector_interpretation"],
        "requires_completion_if_vector": requires_completion,
        "ward_target": ward["ward_identity_target"],
        "naive_high_energy_growth_risk": naive_growth,
        "required_finite_completion_data": required_completion,
        "p6_verdict": verdict,
        "tower_like_uv_hint": tower_hint,
    }


def scalar_auxiliary_row(block: dict) -> dict:
    return {
        "mediator": block["mediator"]["label"],
        "channel_count": block["channel_count"],
        "longitudinal_vector_growth_risk": False,
        "p6_verdict": "scalar_or_auxiliary_safe_at_longitudinal_level",
        "remaining_checks": [
            "spin-statistics projection",
            "candidate action existence",
            "low-energy matching",
            "perturbative unitarity of any induced higher-dimension contact operator",
        ],
        "boundary": (
            "P6 does not declare scalar exchange safe phenomenologically; it only "
            "says the broken-vector longitudinal Ward problem is absent."
        ),
    }


def contact_growth_summary(p2: dict) -> dict:
    allowed = p2["summary"]["allowed_symbolic_pole_components"]
    return {
        "allowed_symbolic_pole_components": allowed,
        "dimension_six_contact_scaling": "A_contact ~ s/Lambda^2 in a four-fermion EFT normalization",
        "verdict": "contact_only_is_not_a_uv_completion",
        "required_resolution": (
            "A contact term may parameterize low-energy matching, but high-energy "
            "behavior requires explicit pole, symmetry cancellation, or tower-like UV data."
        ),
    }


def transition_growth_summary(p4: dict) -> dict:
    transitions = p4["transition_checks"]
    broken = [row for row in transitions if row["requires_broken_sector_completion"]]
    unbroken = [row for row in transitions if not row["requires_broken_sector_completion"]]
    classes: dict[str, int] = {}
    for row in broken:
        classes[row["classification"]] = classes.get(row["classification"], 0) + 1
    return {
        "current_transitions_checked": len(transitions),
        "unbroken_transitions": len(unbroken),
        "broken_transitions_needing_completion": len(broken),
        "broken_transition_classes": classes,
        "interpretation": (
            "Every broken/off-face current needs a generalized Ward completion before "
            "longitudinal high-energy behavior can be judged."
        ),
    }


def candidate_implications(p5: dict) -> list[dict]:
    rows = []
    for candidate in p5["candidates"]:
        name = candidate["candidate"]
        if name == "SU(5)":
            implication = (
                "Finite field-theory completion is possible in principle, but X/Y "
                "broken-vector exchange must come with an SU(5)-breaking Higgs/Goldstone "
                "sector and proton-channel matching."
            )
            tower = "not indicated before finite SU(5) completion is tested"
        elif name == "Pati-Salam":
            implication = (
                "Leptoquark and SU(2)_R broken currents require product-group breaking "
                "data; high-energy closure cannot be judged from isolated poles."
            )
            tower = "not indicated before finite product-group completion is tested"
        elif name == "Spin(10)":
            implication = (
                "The standard Route-C minimal single-object survivor still needs the "
                "full broken D5 generator, mass, Goldstone, and Higgs/source sector for "
                "longitudinal-vector cancellations."
            )
            tower = "not indicated if the finite Spin(10) action-level completion closes"
        elif name == "E6":
            implication = (
                "Extra 10+1 states may participate in high-energy cancellations, but "
                "then exotics lifting and threshold audits become part of the same branch."
            )
            tower = "possible comparison only after finite E6 plus lifting audit"
        else:
            implication = (
                "A string-derived branch can supply an infinite tower and stronger "
                "Regge-like softness, but this requires compactification-specific spectrum "
                "and consistency data."
            )
            tower = (
                "serious comparison point if finite-pole GUT branches cannot satisfy "
                "controlled high-energy behavior"
            )
        rows.append(
            {
                "candidate": name,
                "p5_status": candidate["route_c_status"],
                "p6_implication": implication,
                "tower_like_uv_role": tower,
            }
        )
    return rows


def build_payload() -> dict:
    p2 = load_json(P2_JSON)
    p3 = load_json(P3_JSON)
    p4 = load_json(P4_JSON)
    p5 = load_json(P5_JSON)

    ward_map = p4_by_mediator(p4)
    blocks = p3["mediator_blocks"]
    vector_rows = []
    scalar_rows = []
    for block in blocks:
        label = block["mediator"]["label"]
        if label not in ward_map:
            raise SystemExit(f"P4 has no Ward row for mediator {label}")
        vector_rows.append(vector_growth_row(block, ward_map[label]))
        scalar_rows.append(scalar_auxiliary_row(block))

    requiring_completion = [
        row for row in vector_rows if row["requires_completion_if_vector"]
    ]
    unbroken_ready = [
        row for row in vector_rows if not row["requires_completion_if_vector"]
    ]

    payload = {
        "description": "Route-C P6 controlled high-energy-growth audit",
        "boundary": (
            "P6 is a growth-risk and completion-requirement ledger.  It does not "
            "prove full high-energy softness, and it does not by itself force a "
            "string/tower UV completion."
        ),
        "high_energy_targets": {
            "finite_vector_target": "no uncancelled longitudinal growth such as A ~ s/M_X^2 or A ~ s^2/M_X^4",
            "broken_vector_identity": "p_mu M^mu(V_X)=M_X M(phi_X)",
            "contact_boundary": "four-fermion contact terms are EFT data, not UV softness proofs",
            "tower_condition": (
                "tower-like UV completion becomes motivated only if finite "
                "Goldstone/Higgs/source completions fail or if stronger Regge-like "
                "softness is demanded"
            ),
        },
        "summary": {
            "mediator_sectors_inspected": len(blocks),
            "finite_vector_branches_needing_completion": len(requiring_completion),
            "unbroken_vector_branches_bookkeeping_ready": len(unbroken_ready),
            "scalar_or_auxiliary_interpretations_without_longitudinal_vector_risk": len(
                scalar_rows
            ),
            "contact_only_components_needing_uv_completion": p2["summary"][
                "allowed_symbolic_pole_components"
            ],
            "broken_transitions_needing_completion": p4["summary"][
                "broken_transitions_needing_completion"
            ],
            "tower_like_uv_forced_now": False,
            "tower_like_uv_status": (
                "not forced at P6; retained as a comparison class if finite-pole "
                "completion fails"
            ),
        },
        "transition_growth_summary": transition_growth_summary(p4),
        "diagnostic_seed_counts": diagnostic_seed_counts(p2),
        "vector_branch_audits": vector_rows,
        "scalar_auxiliary_branch_audits": scalar_rows,
        "contact_growth_summary": contact_growth_summary(p2),
        "candidate_implications": candidate_implications(p5),
    }
    return payload


def build_markdown(payload: dict) -> str:
    s = payload["summary"]
    lines = [
        "# Route-C P6 Controlled High-Energy-Growth Audit",
        "",
        "P6 consumes the P3 residue ledger and P4 Ward-identity bookkeeping.",
        "It does not compute complete amplitudes.  It classifies which symbolic",
        "branches have high-energy-growth risk unless broken-vector",
        "Goldstone/Higgs/source completion is supplied, and when a tower-like UV",
        "completion should become a serious comparison point.",
        "",
        "## High-Energy Targets",
        "",
        "For finite massive vector exchange, the dangerous terms are the",
        "longitudinal pieces.  The conservative target is",
        "",
        "```text",
        payload["high_energy_targets"]["finite_vector_target"],
        "```",
        "",
        "For a broken vector the required generalized Ward relation is",
        "",
        "```text",
        payload["high_energy_targets"]["broken_vector_identity"],
        "```",
        "",
        "A contact term can describe the low-energy EFT, but it is not a UV",
        "softness proof.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| mediator sectors inspected | {s['mediator_sectors_inspected']} |",
        f"| finite-vector branches needing completion | {s['finite_vector_branches_needing_completion']} |",
        f"| unbroken-vector branches bookkeeping-ready | {s['unbroken_vector_branches_bookkeeping_ready']} |",
        f"| scalar/auxiliary interpretations without longitudinal-vector risk | {s['scalar_or_auxiliary_interpretations_without_longitudinal_vector_risk']} |",
        f"| contact-only components needing UV completion | {s['contact_only_components_needing_uv_completion']} |",
        f"| broken transitions needing completion | {s['broken_transitions_needing_completion']} |",
        f"| tower-like UV forced now | {s['tower_like_uv_forced_now']} |",
        "",
        f"Tower status: {s['tower_like_uv_status']}.",
        "",
        "## Broken-Current Growth Summary",
        "",
        "| current class | count |",
        "| --- | ---: |",
    ]
    for cls, count in payload["transition_growth_summary"][
        "broken_transition_classes"
    ].items():
        lines.append(f"| {cls} | {count} |")
    lines.extend(
        [
            "",
            "## Vector-Branch Audit",
            "",
            "| mediator | channels | Ward target | P6 verdict | tower-like UV hint |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )
    for row in payload["vector_branch_audits"]:
        lines.append(
            "| `{mediator}` | {channels} | `{target}` | {verdict} | {tower} |".format(
                mediator=row["mediator"],
                channels=row["channel_count"],
                target=row["ward_target"],
                verdict=row["p6_verdict"],
                tower=row["tower_like_uv_hint"],
            )
        )

    lines.extend(
        [
            "",
            "## Scalar or Auxiliary Interpretation",
            "",
            "If a symbolic mediator is interpreted as a scalar or auxiliary",
            "exchange, the broken-vector longitudinal Ward problem is absent.",
            "This is not a phenomenological success claim: spin-statistics, action",
            "existence, matching, and perturbative-unitarity checks remain.",
            "",
            "| quantity | value |",
            "| --- | ---: |",
            f"| scalar/auxiliary rows | {len(payload['scalar_auxiliary_branch_audits'])} |",
            "",
            "## Contact-Only Boundary",
            "",
            "| allowed symbolic pole components | contact verdict | required resolution |",
            "| ---: | --- | --- |",
            "| {allowed} | {verdict} | {resolution} |".format(
                allowed=payload["contact_growth_summary"][
                    "allowed_symbolic_pole_components"
                ],
                verdict=payload["contact_growth_summary"]["verdict"],
                resolution=payload["contact_growth_summary"]["required_resolution"],
            ),
            "",
            "## Candidate Implications",
            "",
            "| candidate | P6 implication | tower-like UV role |",
            "| --- | --- | --- |",
        ]
    )
    for row in payload["candidate_implications"]:
        lines.append(
            "| {candidate} | {implication} | {tower} |".format(
                candidate=row["candidate"],
                implication=row["p6_implication"],
                tower=row["tower_like_uv_role"],
            )
        )

    lines.extend(
        [
            "",
            "## P6 Boundary",
            "",
            payload["boundary"],
            "",
            "P6 therefore does not say that superstring theory is required.  It says",
            "that a string-like tower becomes a meaningful comparison class only if",
            "the finite broken-vector plus Goldstone/Higgs/source completion cannot",
            "cancel the high-energy growth.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "high_energy_growth_audit.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (OUT / "high_energy_growth_audit.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )
    print(json.dumps(payload["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
