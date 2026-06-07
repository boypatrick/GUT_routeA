#!/usr/bin/env python3
"""Route-C P7 low-energy matching and proton-bound readiness report.

P7 consumes the P2 symbolic pole ledger, P3 factorization ledger, and P6
high-energy-growth audit.  It builds a conditional dimension-six matching
ledger.  It does not compute physical proton lifetimes, because the current
Route-C branch still lacks completion data, mediator masses/couplings, flavor
rotations, RG factors, and hadronic matrix elements.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P2_JSON = OUT / "four_point_pole_ansatz.json"
P3_JSON = OUT / "residue_positivity.json"
P6_JSON = OUT / "high_energy_growth_audit.json"


BARYON = {
    "Q": Fraction(1, 3),
    "L": Fraction(0),
    "u^c": Fraction(-1, 3),
    "d^c": Fraction(-1, 3),
    "nu^c": Fraction(0),
    "e^c": Fraction(0),
}

LEPTON = {
    "Q": Fraction(0),
    "L": Fraction(1),
    "u^c": Fraction(0),
    "d^c": Fraction(0),
    "nu^c": Fraction(-1),
    "e^c": Fraction(-1),
}


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


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


def tokens(external_multiset: str) -> list[str]:
    return external_multiset.split()


def quantum_numbers(external_multiset: str) -> dict:
    ts = tokens(external_multiset)
    b = sum(BARYON[t] for t in ts)
    l = sum(LEPTON[t] for t in ts)
    return {
        "B": frac_json(b),
        "L": frac_json(l),
        "B_minus_L": frac_json(b - l),
        "violates_B": b != 0,
        "violates_L": l != 0,
        "preserves_B_minus_L": b - l == 0,
    }


def canonical_operator_name(external_multiset: str) -> str:
    count = Counter(tokens(external_multiset))
    if count == Counter({"Q": 3, "L": 1}):
        return "QQQL"
    if count == Counter({"u^c": 2, "d^c": 1, "e^c": 1}):
        return "u^c u^c d^c e^c"
    if count == Counter({"u^c": 1, "d^c": 2, "nu^c": 1}):
        return "u^c d^c d^c nu^c"
    if count == Counter({"Q": 2, "u^c": 1, "d^c": 1}):
        return "QQ u^c d^c"
    if count == Counter({"Q": 1, "L": 1, "u^c": 1, "e^c": 1}):
        return "QL u^c e^c"
    if count == Counter({"Q": 1, "L": 1, "d^c": 1, "nu^c": 1}):
        return "QL d^c nu^c"
    if count == Counter({"L": 2, "e^c": 1, "nu^c": 1}):
        return "LL e^c nu^c"
    return external_multiset


def operator_status(external_multiset: str) -> tuple[str, str, list[str]]:
    q = quantum_numbers(external_multiset)
    name = canonical_operator_name(external_multiset)
    if q["violates_B"] and q["violates_L"]:
        if "nu^c" in tokens(external_multiset):
            return (
                "BNV_with_sterile_neutrino",
                "B and L violating but B-L conserving; requires nu^c mass/seesaw or light-sterile treatment before proton-channel interpretation",
                ["sterile-neutrino treatment", "flavor basis", "hadronic/RG inputs"],
            )
        return (
            "standard_BNV_seed",
            "B and L violating but B-L conserving; can seed proton-decay Wilson operators after completion and flavor rotation",
            ["completion status", "mediator mass/couplings", "flavor basis", "hadronic/RG inputs"],
        )
    if not q["violates_B"] and not q["violates_L"]:
        return (
            "B_and_L_conserving",
            "not a proton-decay seed at the four-fermion level",
            ["ordinary EFT matching only if phenomenologically needed"],
        )
    return (
        "nonstandard_global_charge_pattern",
        f"{name} has nonstandard B/L bookkeeping and must be inspected manually",
        ["manual B/L audit"],
    )


def p6_by_mediator(p6: dict) -> dict[str, dict]:
    return {row["mediator"]: row for row in p6["vector_branch_audits"]}


def p3_by_mediator(p3: dict) -> dict[str, dict]:
    return {row["mediator"]["label"]: row for row in p3["mediator_blocks"]}


def coupling_symbol(pair: list[str]) -> str:
    return "g_{" + "".join(pair) + "X}"


def matching_row(pole: dict, p3_map: dict, p6_map: dict) -> dict:
    external = pole["external_multiset"]
    mediator = pole["mediator"]["label"]
    left_pair = pole["left_pair"]["pair"]
    right_pair = pole["right_pair"]["pair"]
    status, interpretation, missing = operator_status(external)
    q = quantum_numbers(external)
    p6 = p6_map[mediator]
    p3 = p3_map[mediator]
    coefficient = (
        f"C6[{canonical_operator_name(external)}; {mediator}] = "
        f"{coupling_symbol(left_pair)} {coupling_symbol(right_pair)}^*/M_{mediator}^2"
    )
    if q["violates_B"]:
        p6_block = (
            "blocked_as_vector_until_completion"
            if p6["requires_completion_if_vector"]
            else "vector_bookkeeping_ready"
        )
        bound_status = "not_evaluable_yet"
    else:
        p6_block = "not_a_proton_seed"
        bound_status = "not_applicable"

    return {
        "external_multiset": external,
        "canonical_operator": canonical_operator_name(external),
        "operator_status": status,
        "operator_interpretation": interpretation,
        "global_charges": q,
        "mediator": mediator,
        "left_pair": left_pair,
        "right_pair": right_pair,
        "symbolic_matching_coefficient": coefficient,
        "operator_dimension": pole["operator_dimension"],
        "p3_factorization_available": True,
        "p3_channel_count_for_mediator": p3["channel_count"],
        "p6_vector_completion_status": p6["p6_verdict"],
        "p6_matching_gate": p6_block,
        "physical_proton_bound_status": bound_status,
        "physical_bound_evaluable": False,
        "missing_inputs_before_physical_bound": missing
        + [
            "P6 finite-vector completion if using vector interpretation",
            "physical Wilson basis",
            "RG evolution",
            "lattice/chiral hadronic matrix elements",
            "experimental lifetime limit for the selected channel",
        ],
    }


def bound_interface() -> dict:
    return {
        "dimension_six_width_template": (
            "Gamma(p -> channel a) = sum_ij C_i^phys C_j^{phys *} H_ij^(a)"
        ),
        "bound_template": (
            "tau_a^exp > 1/Gamma_a implies C^phys dagger H^(a) C^phys < 1/tau_a^exp"
        ),
        "single_coefficient_schematic": "|C_i^phys| < sqrt(1/(tau_a^exp H_ii^(a)))",
        "required_inputs": [
            "completed high-energy branch passing P6",
            "mediator masses and couplings",
            "flavor rotations to physical fermion basis",
            "operator basis and chiral contractions",
            "RG factors from matching scale to hadronic scale",
            "lattice/chiral matrix elements",
            "current experimental lower limit for the selected proton channel",
        ],
        "boundary": (
            "P7 records the interface to proton bounds.  It does not insert "
            "numerical experimental limits in this bootstrap ledger."
        ),
    }


def summarize(rows: list[dict]) -> dict:
    by_status = Counter(row["operator_status"] for row in rows)
    bnv = [row for row in rows if row["global_charges"]["violates_B"]]
    blocked_vector = [
        row for row in bnv if row["p6_matching_gate"] == "blocked_as_vector_until_completion"
    ]
    return {
        "allowed_symbolic_poles_processed": len(rows),
        "operator_status_counts": dict(by_status),
        "baryon_violating_symbolic_poles": len(bnv),
        "baryon_violating_vector_poles_blocked_by_p6_completion": len(blocked_vector),
        "physical_proton_bounds_evaluable_now": 0,
        "all_physical_bounds_conditional": True,
    }


def group_rows(rows: list[dict]) -> list[dict]:
    grouped: dict[str, dict] = {}
    for row in rows:
        key = row["canonical_operator"]
        entry = grouped.setdefault(
            key,
            {
                "canonical_operator": key,
                "operator_status": row["operator_status"],
                "global_charges": row["global_charges"],
                "symbolic_pole_count": 0,
                "mediators": [],
                "physical_bound_evaluable": False,
            },
        )
        entry["symbolic_pole_count"] += 1
        entry["mediators"].append(row["mediator"])
    return list(grouped.values())


def build_payload() -> dict:
    p2 = load_json(P2_JSON)
    p3 = load_json(P3_JSON)
    p6 = load_json(P6_JSON)
    p3_map = p3_by_mediator(p3)
    p6_map = p6_by_mediator(p6)
    rows = [matching_row(pole, p3_map, p6_map) for pole in p2["allowed_poles"]]
    return {
        "description": "Route-C P7 low-energy matching and proton-bound readiness report",
        "boundary": (
            "P7 is a conditional matching ledger.  It does not compute physical "
            "proton lifetimes because P6 completion, flavor rotations, mediator "
            "masses/couplings, RG factors, and hadronic matrix elements are not "
            "yet supplied."
        ),
        "summary": summarize(rows),
        "operator_groups": group_rows(rows),
        "matching_rows": rows,
        "proton_bound_interface": bound_interface(),
    }


def build_markdown(payload: dict) -> str:
    s = payload["summary"]
    lines = [
        "# Route-C P7 Low-Energy Matching and Proton-Bound Readiness Report",
        "",
        "P7 consumes the P2 symbolic pole ledger, P3 factorization ledger, and P6",
        "high-energy-growth audit.  It builds conditional dimension-six matching",
        "entries and marks every physical proton bound as conditional until the",
        "missing completion/flavor/hadronic inputs are supplied.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| allowed symbolic poles processed | {s['allowed_symbolic_poles_processed']} |",
        f"| baryon-violating symbolic poles | {s['baryon_violating_symbolic_poles']} |",
        f"| BNV vector poles blocked by P6 completion | {s['baryon_violating_vector_poles_blocked_by_p6_completion']} |",
        f"| physical proton bounds evaluable now | {s['physical_proton_bounds_evaluable_now']} |",
        f"| all physical bounds conditional | {s['all_physical_bounds_conditional']} |",
        "",
        "Operator status counts:",
        "",
        "| status | count |",
        "| --- | ---: |",
    ]
    for status, count in s["operator_status_counts"].items():
        lines.append(f"| {status} | {count} |")

    lines.extend(
        [
            "",
            "## Operator Groups",
            "",
            "| operator | status | B | L | B-L | symbolic poles | bound evaluable |",
            "| --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in payload["operator_groups"]:
        q = row["global_charges"]
        lines.append(
            "| `{op}` | {status} | {b} | {l} | {bl} | {n} | {bound} |".format(
                op=row["canonical_operator"],
                status=row["operator_status"],
                b=q["B"]["fraction"],
                l=q["L"]["fraction"],
                bl=q["B_minus_L"]["fraction"],
                n=row["symbolic_pole_count"],
                bound=row["physical_bound_evaluable"],
            )
        )

    lines.extend(
        [
            "",
            "## Conditional Matching Rows",
            "",
            "| external multiset | mediator | operator status | P6 gate | symbolic coefficient |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload["matching_rows"]:
        lines.append(
            "| `{external}` | `{mediator}` | {status} | {gate} | `{coef}` |".format(
                external=row["external_multiset"],
                mediator=row["mediator"],
                status=row["operator_status"],
                gate=row["p6_matching_gate"],
                coef=row["symbolic_matching_coefficient"],
            )
        )

    interface = payload["proton_bound_interface"]
    lines.extend(
        [
            "",
            "## Proton-Bound Interface",
            "",
            "The physical width template is",
            "",
            "```text",
            interface["dimension_six_width_template"],
            "```",
            "",
            "The corresponding bound template is",
            "",
            "```text",
            interface["bound_template"],
            "```",
            "",
            "Required inputs before numerical proton bounds can be quoted:",
            "",
        ]
    )
    for item in interface["required_inputs"]:
        lines.append(f"- {item}")

    lines.extend(
        [
            "",
            "## P7 Boundary",
            "",
            payload["boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "low_energy_matching_proton_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (OUT / "low_energy_matching_proton_report.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )
    print(json.dumps(payload["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
