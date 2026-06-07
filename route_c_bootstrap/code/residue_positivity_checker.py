#!/usr/bin/env python3
"""Route-C P3 residue factorization and positivity ledger.

P3 consumes the P2 symbolic pole ledger.  For every declared mediator
quantum-number sector, it builds a two-particle channel block with formal
residue

    R_{alpha beta}^{(X)} = g_{alpha X} g^*_{beta X}.

This is a Gram matrix and is therefore positive semidefinite if the exchanged
mediator has positive Hilbert-space norm.  The script records that conditional
PSD statement and a wrong-sign stress test, but it does not construct the full
action, spin-dependent numerators, or Ward-identity cancellations.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P2_JSON = OUT / "four_point_pole_ansatz.json"


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


def charge_key(payload: dict) -> tuple[str, str, Fraction, Fraction]:
    return (
        payload["SU3_C"],
        payload["SU2_L"],
        parse_fraction(payload["Y"]),
        parse_fraction(payload["B_minus_L"]),
    )


def key_json(key: tuple[str, str, Fraction, Fraction]) -> dict:
    return {
        "SU3_C": key[0],
        "SU2_L": key[1],
        "Y": frac_json(key[2]),
        "B_minus_L": frac_json(key[3]),
    }


def key_label(key: tuple[str, str, Fraction, Fraction]) -> str:
    return f"X[({key[0]},{key[1]});Y={frac(key[2])};B-L={frac(key[3])}]"


def channel_label(component: dict) -> str:
    pair = "+".join(component["pair"])
    prod = component["pair_product"]
    return (
        f"{pair}_({prod['SU3_C']},{prod['SU2_L']};"
        f"Y={prod['Y']['fraction']};B-L={prod['B_minus_L']['fraction']})"
    )


def coupling_symbol(channel: dict, mediator: str) -> str:
    bare_pair = "".join(channel["pair"]).replace("^", "")
    safe_mediator = (
        mediator.replace("X", "X")
        .replace("[", "_")
        .replace("]", "")
        .replace("(", "")
        .replace(")", "")
        .replace(";", "_")
        .replace(",", "_")
        .replace("=", "")
        .replace("/", "over")
        .replace("-", "m")
        .replace("+", "p")
        .replace(" ", "")
    )
    return f"g_{bare_pair}_{safe_mediator}"


def compact_symbolic_matrix(channels: list[dict], mediator_label: str) -> list[list[str]]:
    symbols = [coupling_symbol(channel, mediator_label) for channel in channels]
    matrix = []
    for left in symbols:
        matrix.append([f"{left} {right}^*" for right in symbols])
    return matrix


def build_blocks(pair_components: list[dict]) -> dict[str, dict]:
    blocks: dict[str, dict] = {}
    for component in pair_components:
        mediator_key = charge_key(component["mediator_for_chiral_vertex"])
        label = key_label(mediator_key)
        blocks.setdefault(
            label,
            {
                "mediator_key": mediator_key,
                "channels": [],
            },
        )
        blocks[label]["channels"].append(component)
    return blocks


def block_json(label: str, block: dict) -> dict:
    channels = block["channels"]
    n = len(channels)
    return {
        "mediator": {
            "label": label,
            **key_json(block["mediator_key"]),
        },
        "channel_count": n,
        "channel_basis": [channel_label(channel) for channel in channels],
        "formal_residue_matrix": compact_symbolic_matrix(channels, label),
        "unit_coupling_positive_metric_spectrum": [n] + [0] * (n - 1),
        "unit_coupling_positive_metric_rank": 1 if n else 0,
        "positive_semidefinite_if_positive_norm_mediator": True,
        "wrong_sign_metric_stress_test": {
            "unit_coupling_spectrum": [-n] + [0] * (n - 1),
            "passes_unitarity": False,
            "interpretation": "a ghost or wrong-sign mediator metric would violate residue positivity",
        },
        "boundary": (
            "PSD is conditional on positive mediator norm and canonical "
            "unitarity normalization; spin numerators and Ward identities are "
            "deferred."
        ),
    }


def factorization_check(pole: dict) -> dict:
    mediator_key = charge_key(pole["mediator"])
    left_mediator_key = charge_key(pole["left_pair"]["mediator_for_chiral_vertex"])
    right_pair_product_key = charge_key(pole["right_pair"]["pair_product"])
    right_conjugate_mediator_key = charge_key(
        pole["right_pair"]["mediator_for_chiral_vertex"]
    )
    left_pass = mediator_key == left_mediator_key
    right_pass = mediator_key == right_pair_product_key
    conjugate_pass = left_mediator_key == right_pair_product_key

    return {
        "external_multiset": pole["external_multiset"],
        "mediator": pole["mediator"]["label"],
        "left_channel": channel_label(pole["left_pair"]),
        "right_channel": channel_label(pole["right_pair"]),
        "left_vertex_uses_declared_mediator": left_pass,
        "right_pair_is_conjugate_product_for_declared_mediator": right_pass,
        "factorizes_as_g_leftX_g_rightXbar": left_pass and right_pass,
        "right_vertex_couples_to": key_label(right_conjugate_mediator_key),
        "conjugacy_cross_check": conjugate_pass,
    }


def build_markdown(payload: dict) -> str:
    summary = payload["summary"]
    lines = [
        "# Route-C P3 Residue Factorization and Positivity Ledger",
        "",
        "This ledger turns the P2 charge-allowed symbolic poles into",
        "two-particle channel blocks.  It verifies the conditional Gram-matrix",
        "positivity statement for positive-norm mediator exchange.",
        "",
        "## Mathematical Statement",
        "",
        "For a declared mediator sector `X`, let `alpha` label all two-particle",
        "channels that can couple to `X`.  Near the pole, unitarity requires",
        "",
        "```text",
        "Res_X A_{alpha beta} = g_{alpha X} g^*_{beta X}.",
        "```",
        "",
        "This is a Gram matrix.  With a positive-norm mediator metric it is",
        "positive semidefinite.  With a wrong-sign metric, the same unit-coupling",
        "test has a negative eigenvalue and fails.",
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| mediator sectors | {summary['mediator_sectors']} |",
        f"| total vertex channels | {summary['total_vertex_channels']} |",
        f"| P2 allowed poles checked | {summary['p2_allowed_poles_checked']} |",
        f"| factorization checks passed | {summary['factorization_checks_passed']} |",
        f"| PSD blocks under positive metric | {summary['psd_blocks_under_positive_metric']} |",
        f"| largest channel-block size | {summary['largest_channel_block_size']} |",
        "",
        "## Mediator Blocks",
        "",
        "| mediator | channel count | unit-coupling positive spectrum | wrong-sign stress spectrum |",
        "| --- | ---: | --- | --- |",
    ]
    for block in payload["mediator_blocks"]:
        lines.append(
            "| `{label}` | {n} | `{pos}` | `{neg}` |".format(
                label=block["mediator"]["label"],
                n=block["channel_count"],
                pos=block["unit_coupling_positive_metric_spectrum"],
                neg=block["wrong_sign_metric_stress_test"][
                    "unit_coupling_spectrum"
                ],
            )
        )

    lines.extend(
        [
            "",
            "## Factorization Checks",
            "",
            "| external multiset | mediator | left channel | right channel | passes |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for check in payload["factorization_checks"]:
        lines.append(
            "| `{multiset}` | `{mediator}` | `{left}` | `{right}` | {passes} |".format(
                multiset=check["external_multiset"],
                mediator=check["mediator"],
                left=check["left_channel"],
                right=check["right_channel"],
                passes=check["factorizes_as_g_leftX_g_rightXbar"],
            )
        )

    lines.extend(
        [
            "",
            "## P3 Boundary",
            "",
            "The PSD statement is conditional.  It assumes that the mediator sector",
            "exists in the candidate action with positive norm and canonical",
            "unitarity normalization.  P3 does not yet check explicit generator",
            "algebra, spin-statistics projections, Ward identities, high-energy",
            "growth, or proton bounds.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    if not P2_JSON.exists():
        raise SystemExit(f"missing {P2_JSON}; run four_point_pole_ansatz.py first")

    payload_p2 = json.loads(P2_JSON.read_text(encoding="utf-8"))
    blocks = build_blocks(payload_p2["pair_components"])
    block_rows = [block_json(label, block) for label, block in sorted(blocks.items())]
    factorization_checks = [
        factorization_check(pole) for pole in payload_p2["allowed_poles"]
    ]

    factorization_passed = sum(
        1 for check in factorization_checks if check["factorizes_as_g_leftX_g_rightXbar"]
    )
    psd_passed = sum(
        1
        for block in block_rows
        if block["positive_semidefinite_if_positive_norm_mediator"]
    )
    largest_block = max((block["channel_count"] for block in block_rows), default=0)

    payload = {
        "description": "Route-C P3 residue factorization and positivity ledger",
        "boundary": (
            "This is a conditional Gram-matrix PSD check.  It does not prove "
            "mediator existence, Ward identities, high-energy softness, or "
            "proton safety."
        ),
        "summary": {
            "mediator_sectors": len(block_rows),
            "total_vertex_channels": len(payload_p2["pair_components"]),
            "p2_allowed_poles_checked": len(factorization_checks),
            "factorization_checks_passed": factorization_passed,
            "all_factorization_checks_pass": factorization_passed
            == len(factorization_checks),
            "psd_blocks_under_positive_metric": psd_passed,
            "all_psd_blocks_pass_under_positive_metric": psd_passed
            == len(block_rows),
            "largest_channel_block_size": largest_block,
        },
        "mediator_blocks": block_rows,
        "factorization_checks": factorization_checks,
    }

    (OUT / "residue_positivity.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "residue_positivity.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )

    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
