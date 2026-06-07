#!/usr/bin/env python3
"""Route-C P4 Ward-identity bookkeeping ledger.

P4 does not yet compute complete amplitudes.  It separates the current sectors
that are ready for unbroken SM Ward checks from broken/off-face currents that
must be completed by massive-vector, Goldstone, and Higgs/source data before a
real Ward-identity test can be claimed.
"""

from __future__ import annotations

import json
import re
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P1_JSON = OUT / "external_state_charges.json"
P2_JSON = OUT / "four_point_pole_ansatz.json"
P3_JSON = OUT / "residue_positivity.json"


UNBROKEN_ADJOINTS = {
    ("8", "1", Fraction(0), Fraction(0)): "SU(3)_C adjoint",
    ("1", "3", Fraction(0), Fraction(0)): "SU(2)_L adjoint",
    ("1", "1", Fraction(0), Fraction(0)): "U(1)_Y neutral current",
}

PATI_SALAM_LEPTOQUARK_PAIRS = {
    frozenset(("Q", "L")),
    frozenset(("u^c", "nu^c")),
    frozenset(("d^c", "e^c")),
}

SU2R_PAIRS = {
    frozenset(("u^c", "d^c")),
    frozenset(("nu^c", "e^c")),
}


def parse_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def frac_json(value: Fraction) -> dict[str, int | str | float]:
    return {
        "fraction": frac(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "float": float(value),
    }


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_transition(current: str) -> tuple[str, str]:
    match = re.fullmatch(r"\((.+)\)\^dagger -> (.+)", current)
    if not match:
        raise ValueError(f"cannot parse current transition {current!r}")
    return match.group(1), match.group(2)


def classify_transition(source: str, target: str) -> tuple[str, str, bool]:
    if source == target:
        return (
            "unbroken diagonal SM-face current",
            "epsilon_mu -> p_mu Ward check is bookkeeping-ready in the massless/canonical current limit",
            False,
        )
    pair = frozenset((source, target))
    if pair in PATI_SALAM_LEPTOQUARK_PAIRS:
        return (
            "broken Pati-Salam leptoquark current",
            "requires massive leptoquark vector plus Goldstone/Higgs completion before Ward test",
            True,
        )
    if pair in SU2R_PAIRS:
        return (
            "broken SU(2)_R charged current",
            "requires SU(2)_R-breaking vector plus Goldstone/Higgs completion before Ward test",
            True,
        )
    return (
        "broken Spin(10) off-face current",
        "requires declared D5 generator and broken-sector completion before Ward test",
        True,
    )


def unbroken_current_checks(p1: dict) -> list[dict]:
    checks = []
    for sector in p1["current_sectors"]:
        if "SM-preserving" not in sector["classification"]:
            continue
        name = sector["name"]
        if name.startswith("SU(3)"):
            status = (
                "bookkeeping-ready: preserves multiplet label and acts inside color components; "
                "full check needs component-level SU(3) generators"
            )
        elif name.startswith("SU(2)"):
            status = (
                "bookkeeping-ready: preserves multiplet label and acts inside weak doublets; "
                "full check needs component-level SU(2) generators"
            )
        else:
            status = (
                "bookkeeping-ready: diagonal hypercharge current with exact P1 charge conservation"
            )
        checks.append(
            {
                "name": name,
                "representative_current": sector["representative_current"],
                "acts_on": sector["acts_on"],
                "ward_replacement": "epsilon_mu -> p_mu",
                "bookkeeping_pass": True,
                "full_component_ward_proved": False,
                "status": status,
            }
        )
    return checks


def transition_checks(p2: dict) -> list[dict]:
    checks = []
    for transition in p2["current_transition_seed"]:
        source, target = parse_transition(transition["current"])
        y = parse_fraction(transition["charge_carried_by_current"]["Y"])
        bl = parse_fraction(transition["charge_carried_by_current"]["B_minus_L"])
        classification, status, requires_completion = classify_transition(source, target)
        checks.append(
            {
                "current": transition["current"],
                "source": source,
                "target": target,
                "charge_carried_by_current": {
                    "Y": frac_json(y),
                    "B_minus_L": frac_json(bl),
                },
                "classification": classification,
                "requires_broken_sector_completion": requires_completion,
                "massless_ward_bookkeeping_pass": (
                    not requires_completion and y == 0 and bl == 0
                ),
                "status": status,
            }
        )
    return checks


def mediator_key(mediator: dict) -> tuple[str, str, Fraction, Fraction]:
    return (
        mediator["SU3_C"],
        mediator["SU2_L"],
        parse_fraction(mediator["Y"]),
        parse_fraction(mediator["B_minus_L"]),
    )


def mediator_ward_requirements(p3: dict) -> list[dict]:
    rows = []
    for block in p3["mediator_blocks"]:
        med = block["mediator"]
        key = mediator_key(med)
        if key in UNBROKEN_ADJOINTS:
            status = (
                "unbroken adjoint-compatible vector sector; ordinary massless Ward check can be attempted"
            )
            requires_completion = False
            interpretation = UNBROKEN_ADJOINTS[key]
        else:
            status = (
                "if interpreted as a vector mediator, requires massive-vector "
                "Goldstone/Higgs completion and a generalized Ward identity"
            )
            requires_completion = True
            interpretation = "broken or non-adjoint mediator sector"
        rows.append(
            {
                "mediator": med["label"],
                "channel_count": block["channel_count"],
                "vector_interpretation": interpretation,
                "requires_goldstone_higgs_completion_if_vector": requires_completion,
                "ward_identity_target": (
                    "p_mu M^mu(V_X)=0"
                    if not requires_completion
                    else "p_mu M^mu(V_X)=M_X M(phi_X)"
                ),
                "status": status,
            }
        )
    return rows


def build_markdown(payload: dict) -> str:
    summary = payload["summary"]
    lines = [
        "# Route-C P4 Ward-Identity Bookkeeping Ledger",
        "",
        "This is a bookkeeping layer for Ward identities.  It does not yet",
        "compute complete amplitudes.  It identifies which current sectors are",
        "ready for unbroken SM Ward checks and which broken/off-face sectors need",
        "Goldstone/Higgs completion before a real Ward test can be claimed.",
        "",
        "## Ward Targets",
        "",
        "For an unbroken gauge boson, the target check is",
        "",
        "```text",
        "epsilon_mu -> p_mu,    p_mu M^mu = 0.",
        "```",
        "",
        "For a broken massive vector, the target is instead the generalized",
        "Ward/Slavnov-Taylor relation",
        "",
        "```text",
        "p_mu M^mu(V_X) = M_X M(phi_X),",
        "```",
        "",
        "so the Goldstone/Higgs/source sector must be supplied before Route C can",
        "test the amplitude.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| unbroken SM current sectors | {summary['unbroken_sm_current_sectors']} |",
        f"| current transitions checked | {summary['current_transitions_checked']} |",
        f"| massless diagonal transitions ready | {summary['massless_diagonal_transitions_ready']} |",
        f"| broken/off-face transitions needing completion | {summary['broken_transitions_needing_completion']} |",
        f"| mediator sectors inspected | {summary['mediator_sectors_inspected']} |",
        f"| mediator sectors needing completion if vector | {summary['mediator_sectors_needing_completion_if_vector']} |",
        "",
        "## Unbroken SM Current Readiness",
        "",
        "| current sector | acts on | bookkeeping pass | full component Ward proved | status |",
        "| --- | --- | --- | --- | --- |",
    ]
    for check in payload["unbroken_current_checks"]:
        lines.append(
            "| {name} | `{acts}` | {book} | {full} | {status} |".format(
                name=check["name"],
                acts=", ".join(check["acts_on"]),
                book=check["bookkeeping_pass"],
                full=check["full_component_ward_proved"],
                status=check["status"],
            )
        )

    lines.extend(
        [
            "",
            "## Current Transition Classification",
            "",
            "| transition | Y carried | B-L carried | classification | completion needed |",
            "| --- | ---: | ---: | --- | --- |",
        ]
    )
    for check in payload["transition_checks"]:
        lines.append(
            "| `{current}` | {y} | {bl} | {classification} | {need} |".format(
                current=check["current"],
                y=check["charge_carried_by_current"]["Y"]["fraction"],
                bl=check["charge_carried_by_current"]["B_minus_L"]["fraction"],
                classification=check["classification"],
                need=check["requires_broken_sector_completion"],
            )
        )

    lines.extend(
        [
            "",
            "## P3 Mediator Ward Requirements",
            "",
            "| mediator | channels | Ward target | completion needed if vector |",
            "| --- | ---: | --- | --- |",
        ]
    )
    for row in payload["mediator_ward_requirements"]:
        lines.append(
            "| `{mediator}` | {n} | `{target}` | {need} |".format(
                mediator=row["mediator"],
                n=row["channel_count"],
                target=row["ward_identity_target"],
                need=row["requires_goldstone_higgs_completion_if_vector"],
            )
        )

    lines.extend(
        [
            "",
            "## P4 Boundary",
            "",
            "P4 is not a complete Ward-identity proof.  It verifies that the",
            "unbroken SM current sectors are bookkeeping-ready and reports all",
            "broken/off-face sectors that need explicit generator, mass,",
            "Goldstone, Higgs, or source data.  The actual cancellations belong",
            "to later P4 refinements and P6 high-energy-growth checks.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    p1 = load_json(P1_JSON)
    p2 = load_json(P2_JSON)
    p3 = load_json(P3_JSON)

    unbroken = unbroken_current_checks(p1)
    transitions = transition_checks(p2)
    mediators = mediator_ward_requirements(p3)

    diagonal_ready = sum(
        1 for row in transitions if row["massless_ward_bookkeeping_pass"]
    )
    broken_transitions = sum(
        1 for row in transitions if row["requires_broken_sector_completion"]
    )
    mediator_completion = sum(
        1
        for row in mediators
        if row["requires_goldstone_higgs_completion_if_vector"]
    )

    payload = {
        "description": "Route-C P4 Ward-identity bookkeeping ledger",
        "boundary": (
            "This is a Ward-identity readiness ledger, not a complete amplitude "
            "cancellation proof."
        ),
        "ward_targets": {
            "unbroken": "epsilon_mu -> p_mu gives p_mu M^mu = 0",
            "broken_massive_vector": "p_mu M^mu(V_X)=M_X M(phi_X)",
        },
        "summary": {
            "unbroken_sm_current_sectors": len(unbroken),
            "current_transitions_checked": len(transitions),
            "massless_diagonal_transitions_ready": diagonal_ready,
            "broken_transitions_needing_completion": broken_transitions,
            "mediator_sectors_inspected": len(mediators),
            "mediator_sectors_needing_completion_if_vector": mediator_completion,
            "all_broken_transitions_flagged": diagonal_ready + broken_transitions
            == len(transitions),
        },
        "unbroken_current_checks": unbroken,
        "transition_checks": transitions,
        "mediator_ward_requirements": mediators,
        "missing_action_level_data": [
            "component-level SU(3)_C and SU(2)_L generator matrices",
            "explicit broken Spin(10)/Pati-Salam generator matrices",
            "mass matrix for broken vector mediators",
            "Goldstone/Higgs/source sector defining p_mu M^mu(V)=M M(phi)",
            "spin-statistics and numerator structures for the actual amplitudes",
        ],
    }

    (OUT / "ward_identity_ledger.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "ward_identity_ledger.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )

    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
