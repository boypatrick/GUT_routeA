#!/usr/bin/env python3
"""Route-C P8 comparative theory scorecard.

P8 is the first pass that compares candidate branches using the explicit
P1--P7 gates rather than anomaly cancellation or minimality alone.  It is a
hard-gate ledger, not a subjective weighted score.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P0_JSON = OUT / "candidate_ledger.json"
P1_JSON = OUT / "external_state_charges.json"
P2_JSON = OUT / "four_point_pole_ansatz.json"
P3_JSON = OUT / "residue_positivity.json"
P4_JSON = OUT / "ward_identity_ledger.json"
P5_JSON = OUT / "anomaly_chiral_consistency.json"
P6_JSON = OUT / "high_energy_growth_audit.json"
P7_JSON = OUT / "low_energy_matching_proton_report.json"


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def by_name(rows: list[dict], key: str) -> dict[str, dict]:
    return {row[key]: row for row in rows}


def gate(status: str, note: str) -> dict[str, str]:
    return {"status": status, "note": note}


def final_verdict(name: str) -> tuple[str, str, list[str]]:
    if name == "Spin(10)":
        return (
            "leading_conditional_route_c_branch",
            (
                "Only field-theory candidate passing the current minimal "
                "single-object filter, but still blocked from physical "
                "high-energy/proton claims by P6/P7 completion debt."
            ),
            [
                "explicit broken D5 generator matrices",
                "Spin(10)-breaking Higgs/Goldstone/source sector",
                "mediator masses and couplings",
                "flavor-basis Wilson tensors and proton bounds",
            ],
        )
    if name == "SU(5)":
        return (
            "conditional_branch_if_single_object_filter_is_relaxed",
            (
                "Anomaly-consistent and simple, but the family is split as "
                "10 + bar5 + 1 rather than one irreducible six-face object."
            ),
            [
                "explicit statement that Route-A single-object criterion is relaxed",
                "SU(5)-specific broken-vector completion",
                "classic X/Y proton matching and thresholds",
            ],
        )
    if name == "Pati-Salam":
        return (
            "conditional_intermediate_face_not_final_simple_object",
            (
                "Anomaly-consistent and physically natural, but product-group "
                "rather than simple single-object unification."
            ),
            [
                "explicit product-group bootstrap branch",
                "Pati-Salam leptoquark and SU(2)_R broken Ward completion",
                "proton and threshold audit if used as more than an intermediate face",
            ],
        )
    if name == "E6":
        return (
            "conditional_branch_with_exotics_lifting_debt",
            (
                "Simple and contains a Spin(10) 16 inside 27, but the extra "
                "10 + 1 states must be lifted or audited."
            ),
            [
                "mass/lifting mechanism for 10 + 1",
                "check no light chiral exotics remain",
                "E6 threshold/proton audit including lifted states",
            ],
        )
    return (
        "uv_completion_class_not_forced_by_p1_p7",
        (
            "Not a single finite four-dimensional GUT candidate.  Relevant as "
            "a comparison class if finite-pole GUT completions fail high-energy "
            "softness or require a tower."
        ),
        [
            "explicit compactification and massless spectrum",
            "Green-Schwarz/tadpole/modular consistency",
            "chiral index, exotics lifting, and matching to the Route-A skeleton",
        ],
    )


def candidate_gate_row(candidate: dict, p0_map: dict, summaries: dict) -> dict:
    name = candidate["candidate"]
    p0 = p0_map.get(name)
    final, final_note, next_audits = final_verdict(name)
    is_uv_class = name == "Superstring-derived completions"
    single_filter_pass = (
        bool(p0 and p0["survives_minimal_single_object_filter"])
        if p0 is not None
        else "not_applicable_uv_class"
    )
    p5_status = candidate["anomaly_status"]
    if p5_status is True:
        anomaly_gate = "pass"
    elif isinstance(p5_status, str) and "conditional" in p5_status:
        anomaly_gate = "conditional"
    else:
        anomaly_gate = "construction_dependent" if is_uv_class else "conditional"

    return {
        "candidate": name,
        "category": candidate["category"],
        "family_realization": candidate["family_realization"],
        "gates": {
            "C0_representation_filter": gate(
                "pass" if single_filter_pass is True else ("not_applicable" if is_uv_class else "fail_or_relax"),
                (
                    "passes minimal single-object filter"
                    if single_filter_pass is True
                    else (
                        "UV class, not directly comparable as one finite group"
                        if is_uv_class
                        else candidate["route_c_status"]
                    )
                ),
            ),
            "P1_sm_family_baseline": gate(
                "pass" if summaries["p1_all_checks_pass"] else "fail",
                "one-family SM + nu^c charge/anomaly baseline passes exactly",
            ),
            "P2_charge_filtered_poles": gate(
                "shared_prebranch_ledger",
                (
                    f"{summaries['p2_allowed_poles']} SM-face symbolic pole "
                    "components are available before candidate-specific generator data"
                ),
            ),
            "P3_residue_factorization": gate(
                "conditional_pass",
                (
                    f"{summaries['p3_factorization_passed']} symbolic poles "
                    "factorize with PSD residues if mediators have positive norm"
                ),
            ),
            "P4_ward_completion": gate(
                "completion_required",
                (
                    f"{summaries['p4_broken_transitions']} broken/off-face "
                    "transitions need Goldstone/Higgs/source completion"
                ),
            ),
            "P5_anomaly_chiral": gate(
                anomaly_gate,
                candidate["route_c_status"],
            ),
            "P6_high_energy_growth": gate(
                "completion_required",
                (
                    f"{summaries['p6_vector_completion']} finite-vector branches "
                    "need completion; tower-like UV is not forced now"
                ),
            ),
            "P7_proton_bound_readiness": gate(
                "not_evaluable_yet",
                (
                    f"{summaries['p7_bnv_poles']} BNV symbolic poles found; "
                    "0 physical proton bounds evaluable without completion/flavor/hadronic inputs"
                ),
            ),
        },
        "final_status": final,
        "final_note": final_note,
        "next_required_audits": next_audits,
    }


def build_payload() -> dict:
    p0 = load_json(P0_JSON)
    p1 = load_json(P1_JSON)
    p2 = load_json(P2_JSON)
    p3 = load_json(P3_JSON)
    p4 = load_json(P4_JSON)
    p5 = load_json(P5_JSON)
    p6 = load_json(P6_JSON)
    p7 = load_json(P7_JSON)

    p0_map = by_name(p0["candidates"], "name")
    summaries = {
        "p1_all_checks_pass": p1["all_checks_pass"],
        "p2_allowed_poles": p2["summary"]["allowed_symbolic_pole_components"],
        "p3_factorization_passed": p3["summary"]["factorization_checks_passed"],
        "p4_broken_transitions": p4["summary"]["broken_transitions_needing_completion"],
        "p6_vector_completion": p6["summary"][
            "finite_vector_branches_needing_completion"
        ],
        "p7_bnv_poles": p7["summary"]["baryon_violating_symbolic_poles"],
    }
    rows = [candidate_gate_row(c, p0_map, summaries) for c in p5["candidates"]]
    status_counts: dict[str, int] = {}
    for row in rows:
        status_counts[row["final_status"]] = status_counts.get(row["final_status"], 0) + 1
    return {
        "description": "Route-C P8 comparative theory scorecard",
        "boundary": (
            "P8 is a hard-gate scorecard over existing Route-C ledgers.  It does "
            "not upgrade conditional branches into completed GUTs."
        ),
        "summary": {
            "candidates_compared": len(rows),
            "p1_p7_gates_used": [
                "C0 representation pre-filter",
                "P1 charge/anomaly baseline",
                "P2 charge-filtered pole ledger",
                "P3 residue factorization/PSD",
                "P4 Ward-completion readiness",
                "P5 anomaly/chiral/exotics ledger",
                "P6 high-energy-growth completion gate",
                "P7 proton-bound readiness",
            ],
            "final_status_counts": status_counts,
            "leading_conditional_field_theory_branch": "Spin(10)",
            "tower_like_uv_forced_now": False,
        },
        "scorecard": rows,
        "route_c_decision": (
            "Spin(10) is the leading conditional finite field-theory branch under "
            "the current minimal single-object filter.  This is not yet a full "
            "action-level or phenomenological proof.  Superstring-derived "
            "completions remain a UV comparison class, not a P8 conclusion."
        ),
    }


def build_markdown(payload: dict) -> str:
    s = payload["summary"]
    lines = [
        "# Route-C P8 Comparative Theory Scorecard",
        "",
        "P8 compares candidate branches using the explicit P1--P7 gates.  It is a",
        "hard-gate ledger, not a weighted score and not a completed GUT proof.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| candidates compared | {s['candidates_compared']} |",
        f"| leading conditional field-theory branch | {s['leading_conditional_field_theory_branch']} |",
        f"| tower-like UV forced now | {s['tower_like_uv_forced_now']} |",
        "",
        "Final status counts:",
        "",
        "| final status | count |",
        "| --- | ---: |",
    ]
    for status, count in s["final_status_counts"].items():
        lines.append(f"| {status} | {count} |")

    lines.extend(
        [
            "",
            "## Candidate Gate Matrix",
            "",
            "| candidate | C0 | P5 anomaly/chiral | P6 high energy | P7 proton readiness | final status |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload["scorecard"]:
        gates = row["gates"]
        lines.append(
            "| {candidate} | {c0} | {p5} | {p6} | {p7} | {final} |".format(
                candidate=row["candidate"],
                c0=gates["C0_representation_filter"]["status"],
                p5=gates["P5_anomaly_chiral"]["status"],
                p6=gates["P6_high_energy_growth"]["status"],
                p7=gates["P7_proton_bound_readiness"]["status"],
                final=row["final_status"],
            )
        )

    lines.extend(
        [
            "",
            "## Candidate Notes",
            "",
        ]
    )
    for row in payload["scorecard"]:
        lines.extend(
            [
                f"### {row['candidate']}",
                "",
                f"- family realization: `{row['family_realization']}`",
                f"- final status: `{row['final_status']}`",
                f"- note: {row['final_note']}",
                "- next required audits:",
            ]
        )
        for item in row["next_required_audits"]:
            lines.append(f"  - {item}")
        lines.append("")

    lines.extend(
        [
            "## P8 Decision Boundary",
            "",
            payload["route_c_decision"],
            "",
            payload["boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "theory_comparison_scorecard.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (OUT / "theory_comparison_scorecard.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )
    print(json.dumps(payload["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
