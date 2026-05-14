#!/usr/bin/env python3
"""Verify the D5 half-spin weight-hypercube branching used in the lean note.

The script uses integer signs s_i = +/-1 for twice the usual half-spin weights
  w = (1/2)(s_1 e_1 + ... + s_5 e_5).
The positive half-spin representation is chosen as the even-minus sector.

The Pati-Salam face is the split
  D5 -> D3 + D2 ~= so(6) + so(4)
with (s_1,s_2,s_3) carrying SU(4)_C spinor data and (s_4,s_5) carrying
SU(2)_L x SU(2)_R spinor data.
"""

from __future__ import annotations

import csv
import json
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "d5_half_spin"


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def su4_label(color_minus: int) -> tuple[str, str, Fraction]:
    """Return SU(3) label, SU(4) component label, and B-L."""
    if color_minus % 2 == 0:
        if color_minus == 0:
            return "1", "lepton in 4", Fraction(-1, 1)
        return "3", "quark in 4", Fraction(1, 3)
    if color_minus == 3:
        return "1", "anti-lepton in bar4", Fraction(1, 1)
    return "bar3", "anti-quark in bar4", Fraction(-1, 3)


def weak_label(weak_minus: int, s4: int, s5: int) -> tuple[str, Fraction, Fraction]:
    """Return weak representation, T3L, T3R."""
    if weak_minus % 2 == 0:
        return "2_L", Fraction(s4 + s5, 4), Fraction(0, 1)
    return "2_R", Fraction(0, 1), Fraction(s4 - s5, 4)


def field_label(su3: str, weak: str, bminusl: Fraction, t3r: Fraction) -> str:
    y = t3r + bminusl / 2
    if weak == "2_L" and su3 == "3" and y == Fraction(1, 6):
        return "Q"
    if weak == "2_L" and su3 == "1" and y == Fraction(-1, 2):
        return "L"
    if weak == "2_R" and su3 == "bar3" and y == Fraction(-2, 3):
        return "u^c"
    if weak == "2_R" and su3 == "bar3" and y == Fraction(1, 3):
        return "d^c"
    if weak == "2_R" and su3 == "1" and y == Fraction(0, 1):
        return "nu^c"
    if weak == "2_R" and su3 == "1" and y == Fraction(1, 1):
        return "e^c"
    return "unmatched"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    rows = []
    for signs in product([1, -1], repeat=5):
        minus_total = sum(1 for s in signs if s < 0)
        if minus_total % 2 != 0:
            continue

        color_signs = signs[:3]
        weak_signs = signs[3:]
        color_minus = sum(1 for s in color_signs if s < 0)
        weak_minus = sum(1 for s in weak_signs if s < 0)

        su3, su4_component, bminusl = su4_label(color_minus)
        weak, t3l, t3r = weak_label(weak_minus, signs[3], signs[4])
        y = t3r + bminusl / 2
        field = field_label(su3, weak, bminusl, t3r)

        if color_minus % 2 == 0 and weak_minus % 2 == 0:
            ps = "(4,2,1)"
        elif color_minus % 2 == 1 and weak_minus % 2 == 1:
            ps = "(bar4,1,2)"
        else:
            ps = "parity-mismatch"

        rows.append(
            {
                "signs": "".join("+" if s > 0 else "-" for s in signs),
                "minus_total": minus_total,
                "color_minus": color_minus,
                "weak_minus": weak_minus,
                "pati_salam": ps,
                "su4_component": su4_component,
                "su3": su3,
                "weak": weak,
                "B_minus_L": fmt_frac(bminusl),
                "T3L": fmt_frac(t3l),
                "T3R": fmt_frac(t3r),
                "Y": fmt_frac(y),
                "field": field,
            }
        )

    field_counts: dict[str, int] = {}
    ps_counts: dict[str, int] = {}
    for row in rows:
        field_counts[row["field"]] = field_counts.get(row["field"], 0) + 1
        ps_counts[row["pati_salam"]] = ps_counts.get(row["pati_salam"], 0) + 1

    expected_field_counts = {
        "Q": 6,
        "L": 2,
        "u^c": 3,
        "d^c": 3,
        "nu^c": 1,
        "e^c": 1,
    }
    expected_ps_counts = {"(4,2,1)": 8, "(bar4,1,2)": 8}

    trace_y2 = sum(Fraction(row["Y"]) ** 2 for row in rows)
    trace_t3l2 = sum(Fraction(row["T3L"]) ** 2 for row in rows)

    summary = {
        "num_half_spin_weights": len(rows),
        "chirality": "even-minus D5 half-spin sector",
        "pati_salam_counts": ps_counts,
        "field_counts": field_counts,
        "expected_field_counts": expected_field_counts,
        "field_counts_match": field_counts == expected_field_counts,
        "pati_salam_counts_match": ps_counts == expected_ps_counts,
        "unmatched_fields": [row for row in rows if row["field"] == "unmatched"],
        "trace_Y2": fmt_frac(trace_y2),
        "trace_T3L2": fmt_frac(trace_t3l2),
        "kY": fmt_frac(trace_y2 / trace_t3l2),
    }

    csv_path = OUT / "d5_half_spin_weights.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    json_path = OUT / "d5_half_spin_summary.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report_path = OUT / "d5_half_spin_report.md"
    report_path.write_text(
        "\n".join(
            [
                "# D5 Half-Spin Hypercube Verification",
                "",
                f"Half-spin weights: `{summary['num_half_spin_weights']}`.",
                f"Pati-Salam counts: `{ps_counts}`.",
                f"Field counts: `{field_counts}`.",
                f"`Tr Y^2 = {summary['trace_Y2']}`.",
                f"`Tr T3L^2 = {summary['trace_T3L2']}`.",
                f"`k_Y = {summary['kY']}`.",
                f"All field counts match: `{summary['field_counts_match']}`.",
                f"All Pati-Salam counts match: `{summary['pati_salam_counts_match']}`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
