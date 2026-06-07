#!/usr/bin/env python3
"""Route-C P1 external-state and charge ledger.

The table uses left-handed Weyl conventions throughout:

    Q, L, u^c, d^c, nu^c, e^c.

It verifies the one-family hypercharge normalization and anomaly cancellations
that will be used as input for later amplitude-bootstrap stages.  This is still
pre-amplitude bookkeeping; it does not construct scattering amplitudes.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"


@dataclass(frozen=True)
class State:
    name: str
    su3: str
    su3_dim: int
    su3_cubic_index: int
    su2: str
    su2_dim: int
    y: Fraction
    b_minus_l: Fraction
    chirality: str
    spin10_origin: str
    components: tuple[tuple[str, Fraction], ...]

    @property
    def multiplicity(self) -> int:
        return self.su3_dim * self.su2_dim


STATES = [
    State(
        name="Q",
        su3="3",
        su3_dim=3,
        su3_cubic_index=1,
        su2="2",
        su2_dim=2,
        y=Fraction(1, 6),
        b_minus_l=Fraction(1, 3),
        chirality="left-handed Weyl",
        spin10_origin="(4,2,1) under Pati-Salam",
        components=(("u_L", Fraction(2, 3)), ("d_L", Fraction(-1, 3))),
    ),
    State(
        name="L",
        su3="1",
        su3_dim=1,
        su3_cubic_index=0,
        su2="2",
        su2_dim=2,
        y=Fraction(-1, 2),
        b_minus_l=Fraction(-1, 1),
        chirality="left-handed Weyl",
        spin10_origin="(4,2,1) under Pati-Salam",
        components=(("nu_L", Fraction(0)), ("e_L", Fraction(-1, 1))),
    ),
    State(
        name="u^c",
        su3="bar3",
        su3_dim=3,
        su3_cubic_index=-1,
        su2="1",
        su2_dim=1,
        y=Fraction(-2, 3),
        b_minus_l=Fraction(-1, 3),
        chirality="left-handed conjugate Weyl",
        spin10_origin="(bar4,1,2) under Pati-Salam",
        components=(("u^c", Fraction(-2, 3)),),
    ),
    State(
        name="d^c",
        su3="bar3",
        su3_dim=3,
        su3_cubic_index=-1,
        su2="1",
        su2_dim=1,
        y=Fraction(1, 3),
        b_minus_l=Fraction(-1, 3),
        chirality="left-handed conjugate Weyl",
        spin10_origin="(bar4,1,2) under Pati-Salam",
        components=(("d^c", Fraction(1, 3)),),
    ),
    State(
        name="nu^c",
        su3="1",
        su3_dim=1,
        su3_cubic_index=0,
        su2="1",
        su2_dim=1,
        y=Fraction(0),
        b_minus_l=Fraction(1, 1),
        chirality="left-handed conjugate Weyl",
        spin10_origin="(bar4,1,2) under Pati-Salam",
        components=(("nu^c", Fraction(0)),),
    ),
    State(
        name="e^c",
        su3="1",
        su3_dim=1,
        su3_cubic_index=0,
        su2="1",
        su2_dim=1,
        y=Fraction(1, 1),
        b_minus_l=Fraction(1, 1),
        chirality="left-handed conjugate Weyl",
        spin10_origin="(bar4,1,2) under Pati-Salam",
        components=(("e^c", Fraction(1, 1)),),
    ),
]


CURRENT_SECTORS = [
    {
        "name": "SU(3)_C current",
        "representative_current": "bar(psi) gamma_mu T_C psi",
        "acts_on": ["Q", "u^c", "d^c"],
        "classification": "SM-preserving",
        "bootstrap_use": "unbroken gauge Ward identities",
    },
    {
        "name": "SU(2)_L current",
        "representative_current": "bar(psi) gamma_mu T_L psi",
        "acts_on": ["Q", "L"],
        "classification": "SM-preserving",
        "bootstrap_use": "unbroken gauge Ward identities",
    },
    {
        "name": "U(1)_Y current",
        "representative_current": "sum_i Y_i bar(psi_i) gamma_mu psi_i",
        "acts_on": [s.name for s in STATES],
        "classification": "SM-preserving",
        "bootstrap_use": "hypercharge Ward identities and normalization",
    },
    {
        "name": "Pati-Salam leptoquark currents",
        "representative_current": "Q <-> L, u^c <-> nu^c, d^c <-> e^c",
        "acts_on": ["Q/L", "u^c/nu^c", "d^c/e^c"],
        "classification": "broken-GUT schematic",
        "bootstrap_use": "massive mediator pole and proton-decay matching audit",
    },
    {
        "name": "SU(2)_R charged currents",
        "representative_current": "u^c <-> d^c, nu^c <-> e^c",
        "acts_on": ["u^c/d^c", "nu^c/e^c"],
        "classification": "broken-GUT schematic after Y selection",
        "bootstrap_use": "broken-sector Ward and Goldstone completion audit",
    },
    {
        "name": "Spin(10) off-face currents",
        "representative_current": "16 weights connected by broken D5 roots",
        "acts_on": [s.name for s in STATES],
        "classification": "broken-GUT schematic",
        "bootstrap_use": "candidate mediator ledger before explicit D5 generator matrices",
    },
]


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


def state_json(state: State) -> dict:
    return {
        "name": state.name,
        "SU3_C": state.su3,
        "SU3_dimension": state.su3_dim,
        "SU3_cubic_index_convention": state.su3_cubic_index,
        "SU2_L": state.su2,
        "SU2_dimension": state.su2_dim,
        "multiplicity": state.multiplicity,
        "Y": frac_json(state.y),
        "B_minus_L": frac_json(state.b_minus_l),
        "chirality": state.chirality,
        "Spin10_origin": state.spin10_origin,
        "components": [
            {"name": name, "electric_charge_Qem": frac_json(q)}
            for name, q in state.components
        ],
    }


def sum_over_states(values: Iterable[Fraction]) -> Fraction:
    total = Fraction(0)
    for value in values:
        total += value
    return total


def mixed_su3_su3_u1(charge: str) -> Fraction:
    # T(fundamental) = T(antifundamental) = 1/2.
    return sum_over_states(
        Fraction(state.su2_dim, 2) * getattr(state, charge)
        for state in STATES
        if state.su3 in {"3", "bar3"}
    )


def mixed_su2_su2_u1(charge: str) -> Fraction:
    # T(SU(2) doublet) = 1/2.
    return sum_over_states(
        Fraction(state.su3_dim, 2) * getattr(state, charge)
        for state in STATES
        if state.su2 == "2"
    )


def gravitational_u1(charge: str) -> Fraction:
    return sum_over_states(state.multiplicity * getattr(state, charge) for state in STATES)


def cubic_u1(charge: str) -> Fraction:
    return sum_over_states(
        state.multiplicity * getattr(state, charge) ** 3 for state in STATES
    )


def mixed_u1_u1_u1(first: str, second: str, third: str) -> Fraction:
    return sum_over_states(
        state.multiplicity
        * getattr(state, first)
        * getattr(state, second)
        * getattr(state, third)
        for state in STATES
    )


def compute_checks() -> dict:
    tr_y2 = sum_over_states(state.multiplicity * state.y**2 for state in STATES)
    # Fundamental doublet has Tr(T3^2) = 1/2; color copies multiply the trace.
    tr_t3l2 = sum_over_states(
        Fraction(state.su3_dim, 2) for state in STATES if state.su2 == "2"
    )
    hypercharge_ratio = tr_y2 / tr_t3l2

    su2_doublets = sum(state.su3_dim for state in STATES if state.su2 == "2")

    anomaly_checks = {
        "SU3_cubed": sum_over_states(
            state.su2_dim * state.su3_cubic_index for state in STATES
        ),
        "SU3_squared_U1Y": mixed_su3_su3_u1("y"),
        "SU2_squared_U1Y": mixed_su2_su2_u1("y"),
        "grav_squared_U1Y": gravitational_u1("y"),
        "U1Y_cubed": cubic_u1("y"),
        "SU3_squared_BminusL": mixed_su3_su3_u1("b_minus_l"),
        "SU2_squared_BminusL": mixed_su2_su2_u1("b_minus_l"),
        "grav_squared_BminusL": gravitational_u1("b_minus_l"),
        "BminusL_cubed": cubic_u1("b_minus_l"),
        "U1Y_squared_BminusL": mixed_u1_u1_u1("y", "y", "b_minus_l"),
        "U1Y_BminusL_squared": mixed_u1_u1_u1("y", "b_minus_l", "b_minus_l"),
    }

    return {
        "hypercharge_normalization": {
            "Tr_Y_squared": tr_y2,
            "Tr_T3L_squared": tr_t3l2,
            "ratio": hypercharge_ratio,
            "passes": hypercharge_ratio == Fraction(5, 3),
        },
        "anomalies": {
            name: {"value": value, "passes": value == 0}
            for name, value in anomaly_checks.items()
        },
        "witten_SU2": {
            "number_of_left_handed_doublets": su2_doublets,
            "passes_even_doublet_count": su2_doublets % 2 == 0,
        },
    }


def serializable_checks(checks: dict) -> dict:
    out = {
        "hypercharge_normalization": {
            key: frac_json(value) if isinstance(value, Fraction) else value
            for key, value in checks["hypercharge_normalization"].items()
        },
        "anomalies": {},
        "witten_SU2": checks["witten_SU2"],
    }
    for name, payload in checks["anomalies"].items():
        out["anomalies"][name] = {
            "value": frac_json(payload["value"]),
            "passes": payload["passes"],
        }
    return out


def build_markdown(payload: dict, checks: dict) -> str:
    lines = [
        "# Route-C P1 External State and Charge Ledger",
        "",
        "This is a pre-amplitude ledger for one left-handed Standard Model",
        "family plus `nu^c`.  It fixes the charge conventions used by later",
        "pole-factorization, Ward-identity, anomaly, and proton-matching",
        "checks.",
        "",
        "## State Table",
        "",
        "| state | SU(3)_C | SU(2)_L | Y | B-L | multiplicity | convention |",
        "| --- | --- | --- | --- | --- | ---: | --- |",
    ]
    for state in STATES:
        lines.append(
            "| {name} | `{su3}` | `{su2}` | {y} | {bl} | {mult} | {conv} |".format(
                name=state.name,
                su3=state.su3,
                su2=state.su2,
                y=frac(state.y),
                bl=frac(state.b_minus_l),
                mult=state.multiplicity,
                conv=state.chirality,
            )
        )

    lines.extend(
        [
            "",
            "All entries use left-handed Weyl conventions.  The fields",
            "`u^c`, `d^c`, `nu^c`, and `e^c` are left-handed conjugate fields,",
            "so their listed charges are the charges of those conjugate Weyl",
            "fields.",
            "",
            "## Hypercharge Normalization",
            "",
            "For one full family plus `nu^c`,",
            "",
            "```text",
            f"Tr Y^2       = {frac(checks['hypercharge_normalization']['Tr_Y_squared'])}",
            f"Tr T3L^2     = {frac(checks['hypercharge_normalization']['Tr_T3L_squared'])}",
            f"Tr Y^2 / Tr T3L^2 = {frac(checks['hypercharge_normalization']['ratio'])}",
            "```",
            "",
            "The normalization check passes exactly:",
            "",
            "```text",
            f"{checks['hypercharge_normalization']['passes']}",
            "```",
            "",
            "## Anomaly Checks",
            "",
            "| check | exact value | passes |",
            "| --- | ---: | --- |",
        ]
    )
    for name, result in checks["anomalies"].items():
        lines.append(f"| `{name}` | {frac(result['value'])} | {result['passes']} |")

    lines.extend(
        [
            "",
            "The Witten `SU(2)` condition also passes:",
            "",
            "```text",
            "number of left-handed SU(2) doublets = "
            f"{checks['witten_SU2']['number_of_left_handed_doublets']}",
            "even doublet count = "
            f"{checks['witten_SU2']['passes_even_doublet_count']}",
            "```",
            "",
            "## Current Sectors",
            "",
            "The following current ledger is schematic.  It separates currents",
            "that preserve the SM face from currents that leave the SM face and",
            "therefore require broken-GUT mediator and Goldstone/Higgs-sector",
            "data in later Route-C stages.",
            "",
            "| current sector | representative current | classification | later use |",
            "| --- | --- | --- | --- |",
        ]
    )
    for current in CURRENT_SECTORS:
        lines.append(
            "| {name} | `{rep}` | {classification} | {use} |".format(
                name=current["name"],
                rep=current["representative_current"],
                classification=current["classification"],
                use=current["bootstrap_use"],
            )
        )

    lines.extend(
        [
            "",
            "## P1 Boundary",
            "",
            "This file verifies charge and anomaly bookkeeping.  It does not yet",
            "construct four-point amplitudes, residues, or Ward-identity",
            "cancellations.  Those begin in P2 and P4.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = compute_checks()

    all_pass = (
        checks["hypercharge_normalization"]["passes"]
        and all(item["passes"] for item in checks["anomalies"].values())
        and checks["witten_SU2"]["passes_even_doublet_count"]
    )

    payload = {
        "description": "Route-C P1 external-state and charge ledger",
        "convention": "one left-handed SM family plus nu^c",
        "states": [state_json(state) for state in STATES],
        "current_sectors": CURRENT_SECTORS,
        "checks": serializable_checks(checks),
        "all_checks_pass": all_pass,
    }

    (OUT / "external_state_charges.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "external_state_charges.md").write_text(
        build_markdown(payload, checks), encoding="utf-8"
    )

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
