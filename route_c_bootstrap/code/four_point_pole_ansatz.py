#!/usr/bin/env python3
"""Route-C P2 symbolic four-point pole ansatz ledger.

This script consumes the P1 external-state charge table and builds the first
charge-filtered pole ledger for all-incoming left-handed Weyl bilinears.

It does not assert that a candidate action contains any listed mediator.  It
only records which symbolic poles are not forbidden by the SM-face quantum
numbers already fixed in P1.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P1_JSON = OUT / "external_state_charges.json"


SU3_PRODUCTS = {
    ("1", "1"): ["1"],
    ("1", "3"): ["3"],
    ("1", "bar3"): ["bar3"],
    ("3", "1"): ["3"],
    ("bar3", "1"): ["bar3"],
    ("3", "3"): ["6", "bar3"],
    ("bar3", "bar3"): ["bar6", "3"],
    ("3", "bar3"): ["1", "8"],
    ("bar3", "3"): ["1", "8"],
}

SU2_PRODUCTS = {
    ("1", "1"): ["1"],
    ("1", "2"): ["2"],
    ("2", "1"): ["2"],
    ("2", "2"): ["1", "3"],
}

SU3_CONJUGATE_REPS = {
    "1": "1",
    "3": "bar3",
    "bar3": "3",
    "6": "bar6",
    "bar6": "6",
    "8": "8",
}

SU2_CONJUGATE_REPS = {
    "1": "1",
    "2": "2",
    "3": "3",
}


@dataclass(frozen=True)
class State:
    name: str
    su3: str
    su2: str
    y: Fraction
    b_minus_l: Fraction


@dataclass(frozen=True)
class PairComponent:
    pair: tuple[str, str]
    su3_product: str
    su2_product: str
    y_total: Fraction
    b_minus_l_total: Fraction

    @property
    def product_key(self) -> tuple[str, str, Fraction, Fraction]:
        return (
            self.su3_product,
            self.su2_product,
            self.y_total,
            self.b_minus_l_total,
        )

    @property
    def conjugate_product_key(self) -> tuple[str, str, Fraction, Fraction]:
        return (
            conjugate_su3_rep(self.su3_product),
            conjugate_su2_rep(self.su2_product),
            -self.y_total,
            -self.b_minus_l_total,
        )

    @property
    def mediator_key_for_vertex(self) -> tuple[str, str, Fraction, Fraction]:
        return self.conjugate_product_key


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


def conjugate_su3_rep(rep: str) -> str:
    if rep not in SU3_CONJUGATE_REPS:
        raise KeyError(f"missing SU(3) conjugate map for representation {rep!r}")
    return SU3_CONJUGATE_REPS[rep]


def conjugate_su2_rep(rep: str) -> str:
    if rep not in SU2_CONJUGATE_REPS:
        raise KeyError(f"missing SU(2) conjugate map for representation {rep!r}")
    return SU2_CONJUGATE_REPS[rep]


def rep_product(group: str, rep_a: str, rep_b: str) -> list[str]:
    table = SU3_PRODUCTS if group == "SU3" else SU2_PRODUCTS
    key = (rep_a, rep_b)
    if key not in table:
        raise KeyError(f"missing {group} product for {rep_a} x {rep_b}")
    return table[key]


def load_states() -> list[State]:
    payload = json.loads(P1_JSON.read_text(encoding="utf-8"))
    states = []
    for row in payload["states"]:
        states.append(
            State(
                name=row["name"],
                su3=row["SU3_C"],
                su2=row["SU2_L"],
                y=parse_fraction(row["Y"]),
                b_minus_l=parse_fraction(row["B_minus_L"]),
            )
        )
    return states


def build_pair_components(states: list[State]) -> list[PairComponent]:
    by_name = {state.name: state for state in states}
    components: list[PairComponent] = []
    for name_a, name_b in combinations_with_replacement(by_name, 2):
        a = by_name[name_a]
        b = by_name[name_b]
        for su3_prod in rep_product("SU3", a.su3, b.su3):
            for su2_prod in rep_product("SU2", a.su2, b.su2):
                components.append(
                    PairComponent(
                        pair=(name_a, name_b),
                        su3_product=su3_prod,
                        su2_product=su2_prod,
                        y_total=a.y + b.y,
                        b_minus_l_total=a.b_minus_l + b.b_minus_l,
                    )
                )
    return components


def pair_component_json(component: PairComponent) -> dict:
    med = component.mediator_key_for_vertex
    return {
        "pair": list(component.pair),
        "pair_product": {
            "SU3_C": component.su3_product,
            "SU2_L": component.su2_product,
            "Y": frac_json(component.y_total),
            "B_minus_L": frac_json(component.b_minus_l_total),
        },
        "mediator_for_chiral_vertex": {
            "SU3_C": med[0],
            "SU2_L": med[1],
            "Y": frac_json(med[2]),
            "B_minus_L": frac_json(med[3]),
        },
    }


def mediator_label(key: tuple[str, str, Fraction, Fraction]) -> str:
    su3, su2, y, bl = key
    return f"X[({su3},{su2});Y={frac(y)};B-L={frac(bl)}]"


def pair_label(component: PairComponent) -> str:
    a, b = component.pair
    return (
        f"({a},{b})_"
        f"({component.su3_product},{component.su2_product};"
        f"Y={frac(component.y_total)};B-L={frac(component.b_minus_l_total)})"
    )


def external_multiset(left: PairComponent, right: PairComponent) -> str:
    states = sorted([*left.pair, *right.pair])
    return " ".join(states)


def build_allowed_poles(components: list[PairComponent]) -> list[dict]:
    allowed = []
    for idx, left in enumerate(components):
        for right in components[idx:]:
            if left.conjugate_product_key != right.product_key:
                continue
            mediator = left.mediator_key_for_vertex
            channel = {
                "left_pair": pair_component_json(left),
                "right_pair": pair_component_json(right),
                "external_multiset": external_multiset(left, right),
                "mediator": {
                    "label": mediator_label(mediator),
                    "SU3_C": mediator[0],
                    "SU2_L": mediator[1],
                    "Y": frac_json(mediator[2]),
                    "B_minus_L": frac_json(mediator[3]),
                },
                "symbolic_pole": (
                    "R_s[{mediator}]/(s-M_{mediator}^2), "
                    "R_s = g_{left} g_{right_bar}"
                ).format(
                    mediator=mediator_label(mediator),
                    left="".join(left.pair),
                    right_bar="".join(right.pair) + "_barX",
                ),
                "contact_term": "C6({left};{right})".format(
                    left="+".join(left.pair), right="+".join(right.pair)
                ),
                "operator_dimension": 6,
                "filter_status": "allowed by SM-face charge and conjugacy only",
            }
            allowed.append(channel)
    return allowed


def build_forbidden_summary(components: list[PairComponent]) -> dict:
    total = 0
    abelian_forbidden = 0
    nonabelian_forbidden = 0
    for idx, left in enumerate(components):
        for right in components[idx:]:
            total += 1
            if left.y_total + right.y_total != 0 or left.b_minus_l_total + right.b_minus_l_total != 0:
                abelian_forbidden += 1
            elif left.conjugate_product_key != right.product_key:
                nonabelian_forbidden += 1
    return {
        "pair_component_pair_tests": total,
        "forbidden_by_abelian_charge": abelian_forbidden,
        "forbidden_by_nonabelian_conjugacy_after_abelian_pass": nonabelian_forbidden,
    }


def find_seed(channels: list[dict], required_multiset: str) -> list[dict]:
    target = " ".join(sorted(required_multiset.split()))
    return [channel for channel in channels if channel["external_multiset"] == target]


def current_transition_seed(states: list[State]) -> list[dict]:
    """Build a schematic current transition ledger.

    A current psi_a^\dagger psi_b carries additive charges q_b - q_a.  This is
    only a charge ledger for later vector-exchange Ward checks; non-Abelian
    generator matrices are deferred.
    """

    transitions = []
    for source in states:
        for target in states:
            if source.name == target.name:
                classification = "diagonal SM-face or Cartan current"
            else:
                classification = (
                    "off-diagonal candidate current; requires declared GUT generator"
                )
            transitions.append(
                {
                    "current": f"({source.name})^dagger -> {target.name}",
                    "charge_carried_by_current": {
                        "Y": frac_json(target.y - source.y),
                        "B_minus_L": frac_json(
                            target.b_minus_l - source.b_minus_l
                        ),
                    },
                    "classification": classification,
                }
            )
    return transitions


def build_markdown(payload: dict) -> str:
    summary = payload["summary"]
    lines = [
        "# Route-C P2 Symbolic Four-Point Pole Ansatz",
        "",
        "This is the first charge-filtered pole ledger.  It is not yet a",
        "unitarity, positivity, or Ward-identity proof.  A pole listed here is",
        "only a symbolic term not forbidden by the P1 SM-face charge table.",
        "",
        "## Ansatz",
        "",
        "For all-incoming left-handed Weyl bilinears, the seed ansatz is",
        "",
        "```text",
        "A_4(12|34) = sum_X R_s(X)/(s-M_X^2) + C6(12|34),",
        "R_s(X) = g_{12X} g_{34Xbar}.",
        "```",
        "",
        "The script declares the quantum numbers of `X` from the pair product",
        "and keeps every residue symbolic.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| pair types | {summary['pair_types']} |",
        f"| pair-product components | {summary['pair_product_components']} |",
        f"| tested pair-component pairings | {summary['pair_component_pair_tests']} |",
        f"| allowed symbolic pole components | {summary['allowed_symbolic_pole_components']} |",
        f"| forbidden by Abelian charge | {summary['forbidden_by_abelian_charge']} |",
        "| forbidden by non-Abelian conjugacy after Abelian pass | "
        f"{summary['forbidden_by_nonabelian_conjugacy_after_abelian_pass']} |",
        "",
        "## Representative Allowed Poles",
        "",
        "| left pair | right pair | mediator | external multiset |",
        "| --- | --- | --- | --- |",
    ]

    for channel in payload["allowed_poles"][:40]:
        left = channel["left_pair"]
        right = channel["right_pair"]
        lines.append(
            "| `{left_pair}` in `{left_rep}` | `{right_pair}` in `{right_rep}` | `{med}` | `{multiset}` |".format(
                left_pair="+".join(left["pair"]),
                left_rep="({},{})".format(
                    left["pair_product"]["SU3_C"], left["pair_product"]["SU2_L"]
                ),
                right_pair="+".join(right["pair"]),
                right_rep="({},{})".format(
                    right["pair_product"]["SU3_C"], right["pair_product"]["SU2_L"]
                ),
                med=channel["mediator"]["label"],
                multiset=channel["external_multiset"],
            )
        )

    if len(payload["allowed_poles"]) > 40:
        lines.append(
            f"| ... | ... | ... | {len(payload['allowed_poles']) - 40} more in JSON |"
        )

    lines.extend(
        [
            "",
            "## Proton-Operator Seeds",
            "",
            "These entries are useful later for low-energy matching.  Their",
            "presence here means only that a charge-allowed symbolic pole exists;",
            "P3--P7 must still check residues, Ward identities, high-energy",
            "behavior, and bounds.",
            "",
            "| external multiset | allowed pole components |",
            "| --- | ---: |",
        ]
    )

    for name, channels in payload["diagnostic_seeds"].items():
        lines.append(f"| `{name}` | {len(channels)} |")

    lines.extend(
        [
            "",
            "## Current-Transition Charge Ledger",
            "",
            "The current ledger records additive charges for schematic",
            "`psi_a^dagger psi_b` transitions.  It is included for later",
            "vector-exchange Ward checks, but it does not yet provide explicit",
            "GUT generator matrices.",
            "",
            "| transition | Y current charge | B-L current charge | classification |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    for transition in payload["current_transition_seed"]:
        y = transition["charge_carried_by_current"]["Y"]["fraction"]
        bl = transition["charge_carried_by_current"]["B_minus_L"]["fraction"]
        lines.append(
            "| `{current}` | {y} | {bl} | {classification} |".format(
                current=transition["current"],
                y=y,
                bl=bl,
                classification=transition["classification"],
            )
        )

    lines.extend(
        [
            "",
            "## P2 Boundary",
            "",
            "Allowed here means allowed by SM-face charge conservation and",
            "non-Abelian conjugacy of the pair products.  It does not mean the",
            "mediator exists in a candidate action, that the residue has positive",
            "norm, or that Ward identities cancel.  Those are P3--P6 tasks.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    if not P1_JSON.exists():
        raise SystemExit(
            f"missing {P1_JSON}; run external_state_charge_table.py first"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    states = load_states()
    components = build_pair_components(states)
    allowed = build_allowed_poles(components)
    forbidden = build_forbidden_summary(components)

    pair_types = {
        tuple(component.pair)
        for component in components
    }
    diagnostic_seeds = {
        "Q Q Q L": find_seed(allowed, "Q Q Q L"),
        "u^c u^c d^c e^c": find_seed(allowed, "u^c u^c d^c e^c"),
        "Q Q u^c e^c": find_seed(allowed, "Q Q u^c e^c"),
        "Q L u^c d^c": find_seed(allowed, "Q L u^c d^c"),
    }
    payload = {
        "description": "Route-C P2 symbolic four-point pole ansatz ledger",
        "boundary": (
            "Allowed poles are charge-filtered symbolic candidates only; "
            "mediator existence, residue positivity, and Ward identities are "
            "deferred."
        ),
        "ansatz": {
            "formula": "A4(12|34)=sum_X R_s(X)/(s-M_X^2)+C6(12|34)",
            "residue": "R_s(X)=g_12X g_34Xbar",
            "contact_operator_dimension": 6,
        },
        "summary": {
            "pair_types": len(pair_types),
            "pair_product_components": len(components),
            "allowed_symbolic_pole_components": len(allowed),
            **forbidden,
        },
        "pair_components": [pair_component_json(component) for component in components],
        "allowed_poles": allowed,
        "diagnostic_seeds": diagnostic_seeds,
        "current_transition_seed": current_transition_seed(states),
    }

    (OUT / "four_point_pole_ansatz.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "four_point_pole_ansatz.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )

    print(json.dumps(payload["summary"], indent=2))
    print("diagnostic seeds:")
    for name, channels in diagnostic_seeds.items():
        print(f"  {name}: {len(channels)}")


if __name__ == "__main__":
    main()
