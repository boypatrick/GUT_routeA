#!/usr/bin/env python3
"""Combined charge-table check for Route-B and source stabilization.

No web lookup is used.  This script merges the Route-B mediator grading with
the source-stabilization sector into a single selection-rule table.

The nontrivial point is that U(1)_M alone cannot distinguish the neutral AA
link from the neutral BC link.  A minimal Z3 edge grading does:

    e(A)=0, e(B)=e(C)=1 mod 3.

Then AA has edge charge 0, BC has edge charge 2, and the corresponding link
fields carry 0 and 1.  This keeps L_AA from coupling to BC while preserving
the required link terms.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "combined_charge_table"


@dataclass(frozen=True)
class Charge:
    r: Fraction
    m: int
    e: int

    def as_json(self) -> dict[str, Any]:
        return {
            "R": str(self.r),
            "U1M": self.m,
            "Z3E": self.e % 3,
        }


SUPERPOTENTIAL_R = Fraction(2, 1)


def ch(r: int | Fraction, m: int, e: int) -> Charge:
    return Charge(Fraction(r), m, e % 3)


def field_charges() -> dict[str, Charge]:
    q: dict[str, Charge] = {}

    # Mediator nodes.
    q.update(
        {
            "A": ch(1, 0, 0),
            "B": ch(1, 1, 1),
            "C": ch(1, -1, 1),
        }
    )

    # Endpoint link fields.  Z3E link charge is minus the edge charge of the
    # corresponding bilinear.
    q.update(
        {
            "L_AA_1": ch(0, 0, 0),
            "L_AA_54": ch(0, 0, 0),
            "L_AA_210": ch(0, 0, 0),
            "L_AB_1": ch(0, -1, 2),
            "L_AB_54": ch(0, -1, 2),
            "L_AC_1": ch(0, 1, 2),
            "L_AC_54": ch(0, 1, 2),
            "L_BC_1": ch(0, 0, 1),
        }
    )

    # Endpoint source copies.  A single common S_54 is incompatible with U(1)_M
    # for AA/AB/AC simultaneously, so the combined completion uses endpoint
    # source copies with identical Spin(10) orbit constraints.
    q.update(
        {
            "S_AA_54": ch(0, 0, 0),
            "S_AB_54": ch(0, -1, 2),
            "S_AC_54": ch(0, 1, 2),
            "Omega_AA_210": ch(0, 0, 0),
            "Phi_C_210": ch(Fraction(2, 3), 0, 0),
        }
    )

    # Link/source drivers.  Their charge is the inverse of the linked field and
    # R=2, so D(L-S) and X(L-v) are allowed.
    q.update(
        {
            "D_AA_54": ch(2, 0, 0),
            "D_AB_54": ch(2, 1, 1),
            "D_AC_54": ch(2, -1, 1),
            "D_AA_210": ch(2, 0, 0),
            "X_AA_1": ch(2, 0, 0),
            "X_AA_54": ch(2, 0, 0),
            "X_AA_210": ch(2, 0, 0),
            "X_AB_n": ch(2, 1, 1),
            "X_AB_r": ch(2, 1, 1),
            "X_AC_n": ch(2, -1, 1),
            "X_AC_r": ch(2, -1, 1),
            "X_BC_1": ch(2, 0, 2),
            "Y_AB_54": ch(2, 1, 1),
            "Y_AC_54": ch(2, -1, 1),
        }
    )

    # Link-vev spurions used in linear driving constraints.  These are not
    # inserted into mediator bilinears; they are the charge-carrying constants
    # that make X(L-v) invariant.
    q.update(
        {
            "v_AA_1": q["L_AA_1"],
            "v_AA_54": q["L_AA_54"],
            "v_AA_210": q["L_AA_210"],
            "v_AB_54": q["L_AB_54"],
            "v_AC_54": q["L_AC_54"],
            "v_BC_1": q["L_BC_1"],
        }
    )

    # Dimensionful/color-cubic spurion.  M_C Phi_C^2 and Phi_C^3 both have R=2.
    q["M_C"] = ch(Fraction(2, 3), 0, 0)

    # Matter/Yukawa/Majorana sector.
    q.update(
        {
            "Psi16": ch(1, 0, 0),
            "H10": ch(0, 0, 0),
            "H120": ch(0, 0, 0),
            "H126bar": ch(0, 0, 0),
            "H16bar": ch(0, 0, 0),
            "Lambda_inv": ch(0, 0, 0),
        }
    )

    # Optional harmless cubic mediator spurion.
    q["rho_ABC"] = ch(-1, 0, 1)
    return q


def term_charge(term: list[str], charges: dict[str, Charge]) -> Charge:
    total_r = sum((charges[name].r for name in term), Fraction(0))
    total_m = sum(charges[name].m for name in term)
    total_e = sum(charges[name].e for name in term) % 3
    return Charge(total_r, total_m, total_e)


def is_allowed(term: list[str], charges: dict[str, Charge]) -> bool:
    c = term_charge(term, charges)
    return c.r == SUPERPOTENTIAL_R and c.m == 0 and c.e % 3 == 0


def term(name: str, factors: list[str], expected: bool, sector: str) -> dict[str, Any]:
    return {"name": name, "factors": factors, "expected": expected, "sector": sector}


def required_and_forbidden_terms() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    # W_link.
    rows += [
        term("L_AA^1 A A", ["L_AA_1", "A", "A"], True, "W_link"),
        term("L_AA^54 A A", ["L_AA_54", "A", "A"], True, "W_link"),
        term("L_AA^210 A A", ["L_AA_210", "A", "A"], True, "W_link"),
        term("L_AB^1 A B", ["L_AB_1", "A", "B"], True, "W_link"),
        term("L_AB^54 A B", ["L_AB_54", "A", "B"], True, "W_link"),
        term("L_AC^1 A C", ["L_AC_1", "A", "C"], True, "W_link"),
        term("L_AC^54 A C", ["L_AC_54", "A", "C"], True, "W_link"),
        term("L_BC^1 B C", ["L_BC_1", "B", "C"], True, "W_link"),
        term("rho A[B,C]", ["rho_ABC", "A", "B", "C"], True, "W_link_optional"),
    ]

    # W_drive.
    rows += [
        term("X_AA^1 L_AA^1", ["X_AA_1", "L_AA_1"], True, "W_drive"),
        term("X_AA^1 v_AA^1", ["X_AA_1", "v_AA_1"], True, "W_drive"),
        term("X_AA^54 L_AA^54", ["X_AA_54", "L_AA_54"], True, "W_drive"),
        term("X_AA^210 L_AA^210", ["X_AA_210", "L_AA_210"], True, "W_drive"),
        term("X_AB^n L_AB^54", ["X_AB_n", "L_AB_54"], True, "W_drive"),
        term("X_AB^n v_AB^54", ["X_AB_n", "v_AB_54"], True, "W_drive"),
        term("X_AB^r L_AB^1", ["X_AB_r", "L_AB_1"], True, "W_drive"),
        term("X_AB^r L_AB^54", ["X_AB_r", "L_AB_54"], True, "W_drive"),
        term("X_AC^n L_AC^54", ["X_AC_n", "L_AC_54"], True, "W_drive"),
        term("X_AC^n v_AC^54", ["X_AC_n", "v_AC_54"], True, "W_drive"),
        term("X_AC^r L_AC^1", ["X_AC_r", "L_AC_1"], True, "W_drive"),
        term("X_AC^r L_AC^54", ["X_AC_r", "L_AC_54"], True, "W_drive"),
        term("X_BC^1 L_BC^1", ["X_BC_1", "L_BC_1"], True, "W_drive"),
        term("X_BC^1 v_BC^1", ["X_BC_1", "v_BC_1"], True, "W_drive"),
    ]

    # W_src.
    rows += [
        term("D_AA^54 L_AA^54", ["D_AA_54", "L_AA_54"], True, "W_src"),
        term("D_AA^54 S_AA^54", ["D_AA_54", "S_AA_54"], True, "W_src"),
        term("D_AB^54 L_AB^54", ["D_AB_54", "L_AB_54"], True, "W_src"),
        term("D_AB^54 S_AB^54", ["D_AB_54", "S_AB_54"], True, "W_src"),
        term("D_AC^54 L_AC^54", ["D_AC_54", "L_AC_54"], True, "W_src"),
        term("D_AC^54 S_AC^54", ["D_AC_54", "S_AC_54"], True, "W_src"),
        term("D_AA^210 L_AA^210", ["D_AA_210", "L_AA_210"], True, "W_src"),
        term("D_AA^210 Omega_AA^210", ["D_AA_210", "Omega_AA_210"], True, "W_src"),
        term("Phi_C^3", ["Phi_C_210", "Phi_C_210", "Phi_C_210"], True, "W_src"),
        term("M_C Phi_C^2", ["M_C", "Phi_C_210", "Phi_C_210"], True, "W_src"),
    ]

    # Minimal relative-alignment driver for endpoint 54 source copies.
    rows += [
        term("Y_AB^54 S_AB^54", ["Y_AB_54", "S_AB_54"], True, "W_rel"),
        term("Y_AB^54 v_AB^54 S_AA^54", ["Y_AB_54", "v_AB_54", "S_AA_54"], True, "W_rel"),
        term("Y_AC^54 S_AC^54", ["Y_AC_54", "S_AC_54"], True, "W_rel"),
        term("Y_AC^54 v_AC^54 S_AA^54", ["Y_AC_54", "v_AC_54", "S_AA_54"], True, "W_rel"),
        term("Y_AB^54 S_AC^54", ["Y_AB_54", "S_AC_54"], False, "forbidden_rel_cross"),
        term("Y_AC^54 S_AB^54", ["Y_AC_54", "S_AB_54"], False, "forbidden_rel_cross"),
        term("Y_AB^54 S_AA^54", ["Y_AB_54", "S_AA_54"], False, "forbidden_rel_cross"),
        term("Y_AC^54 S_AA^54", ["Y_AC_54", "S_AA_54"], False, "forbidden_rel_cross"),
    ]

    # Yukawa/Majorana.
    rows += [
        term("16 16 10_H", ["Psi16", "Psi16", "H10"], True, "Yukawa"),
        term("16 16 120_H", ["Psi16", "Psi16", "H120"], True, "Yukawa"),
        term("16 16 overline126_H", ["Psi16", "Psi16", "H126bar"], True, "Yukawa/Majorana"),
        term(
            "16 16 overline16_H overline16_H / Lambda",
            ["Psi16", "Psi16", "H16bar", "H16bar", "Lambda_inv"],
            True,
            "Majorana_dim5",
        ),
    ]

    # Dangerous generic mediator/source terms.
    rows += [
        term("bare B B", ["B", "B"], False, "forbidden_generic"),
        term("bare C C", ["C", "C"], False, "forbidden_generic"),
        term("bare B C", ["B", "C"], False, "forbidden_generic"),
        term("L_AA^54 B C", ["L_AA_54", "B", "C"], False, "forbidden_wrong_portal"),
        term("L_BC^1 A A", ["L_BC_1", "A", "A"], False, "forbidden_wrong_portal"),
        term("L_AB^54 A C", ["L_AB_54", "A", "C"], False, "forbidden_wrong_portal"),
        term("L_AC^54 A B", ["L_AC_54", "A", "B"], False, "forbidden_wrong_portal"),
        term("S_AA^54 B C", ["S_AA_54", "B", "C"], False, "forbidden_source_leakage"),
        term("S_AB^54 B B", ["S_AB_54", "B", "B"], False, "forbidden_source_leakage"),
        term("S_AC^54 C C", ["S_AC_54", "C", "C"], False, "forbidden_source_leakage"),
        term("Omega_W^3", ["Omega_AA_210", "Omega_AA_210", "Omega_AA_210"], False, "forbidden_source_cubic"),
        term("S_AA Omega_W^2", ["S_AA_54", "Omega_AA_210", "Omega_AA_210"], False, "forbidden_source_cubic"),
        term("Omega_W Phi_C^2", ["Omega_AA_210", "Phi_C_210", "Phi_C_210"], False, "forbidden_source_cubic"),
    ]
    return rows


def evaluate_terms(charges: dict[str, Charge]) -> list[dict[str, Any]]:
    evaluated = []
    for row in required_and_forbidden_terms():
        c = term_charge(row["factors"], charges)
        actual = is_allowed(row["factors"], charges)
        evaluated.append(
            {
                **row,
                "actual": actual,
                "pass": actual == row["expected"],
                "total_charge": c.as_json(),
            }
        )
    return evaluated


def summarize(evaluated: list[dict[str, Any]]) -> dict[str, Any]:
    by_sector: dict[str, dict[str, Any]] = {}
    for row in evaluated:
        sec = row["sector"]
        by_sector.setdefault(sec, {"terms": 0, "passes": 0, "failures": []})
        by_sector[sec]["terms"] += 1
        by_sector[sec]["passes"] += int(row["pass"])
        if not row["pass"]:
            by_sector[sec]["failures"].append(row)
    return {
        "terms_checked": len(evaluated),
        "passes": sum(1 for row in evaluated if row["pass"]),
        "failures": [row for row in evaluated if not row["pass"]],
        "all_pass": all(row["pass"] for row in evaluated),
        "by_sector": by_sector,
    }


def compatibility_notes() -> list[str]:
    return [
        "The original source-parity assignment with odd link fields is not kept, because it would forbid W_link terms L_ij X_i X_j. R-sequestering alone forbids Omega_W^3, S Omega_W^2, and Omega_W Phi_C^2.",
        "A single common S_54 cannot couple to AA, AB, and AC links while preserving U(1)_M. The combined table uses endpoint source copies S_AA, S_AB, S_AC with identical Spin(10) orbit constraints.",
        "Endpoint source copies leave relative 54-orientation moduli unless one adds the allowed drivers Y_AB and Y_AC; these enforce S_AB proportional to S_AA and S_AC proportional to S_AA without opening cross-alignment terms.",
        "U(1)_M alone cannot distinguish neutral AA and BC channels. The Z3 edge grading with e(A)=0, e(B)=e(C)=1 forbids L_AA BC and L_BC AA.",
        "The Yukawa/Majorana sector is neutral under U(1)_M and Z3E, with R(Psi16)=1 and R(H)=0, so standard Spin(10) Yukawa and Majorana operators remain allowed.",
    ]


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Combined charge-table verification")
    lines.append("")
    lines.append("No web lookup was used.  This combines the Route-B mediator grading")
    lines.append("with the source-stabilization sector and the Yukawa/Majorana sector.")
    lines.append("")
    lines.append("## Symmetry")
    lines.append("")
    lines.append("```text")
    lines.append("G_sel = U(1)_R x U(1)_M x Z3_E,")
    lines.append("R(W)=2, e(A)=0, e(B)=e(C)=1 mod 3.")
    lines.append("```")
    lines.append("")
    lines.append("## Key compatibility fixes")
    lines.append("")
    for note in payload["compatibility_notes"]:
        lines.append(f"- {note}")
    lines.append("")
    lines.append("## Field charges")
    lines.append("")
    lines.append("| field | R | U1M | Z3E |")
    lines.append("|---|---:|---:|---:|")
    for name, c in payload["field_charges"].items():
        lines.append(f"| `{name}` | {c['R']} | {c['U1M']} | {c['Z3E']} |")
    lines.append("")
    lines.append("## Term check summary")
    lines.append("")
    summary = payload["summary"]
    lines.append("```text")
    lines.append(f"terms checked = {summary['terms_checked']}")
    lines.append(f"passes = {summary['passes']}")
    lines.append(f"failures = {len(summary['failures'])}")
    lines.append(f"all pass = {summary['all_pass']}")
    lines.append("```")
    lines.append("")
    lines.append("| sector | terms | passes |")
    lines.append("|---|---:|---:|")
    for sec, row in summary["by_sector"].items():
        lines.append(f"| `{sec}` | {row['terms']} | {row['passes']} |")
    lines.append("")
    lines.append("## Checked terms")
    lines.append("")
    lines.append("| term | expected | actual | total charge |")
    lines.append("|---|---:|---:|---|")
    for row in payload["evaluated_terms"]:
        charge = row["total_charge"]
        lines.append(
            f"| `{row['name']}` | {row['expected']} | {row['actual']} | "
            f"R={charge['R']}, M={charge['U1M']}, E={charge['Z3E']} |"
        )
    (OUT / "combined_charge_table_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    charges = field_charges()
    evaluated = evaluate_terms(charges)
    summary = summarize(evaluated)
    payload = {
        "note": "No web lookup used. Combined Route-B/source-stabilization/Yukawa charge-table check.",
        "symmetry": {
            "group": "U(1)_R x U(1)_M x Z3_E",
            "superpotential_R": "2",
            "edge_grading": {"A": 0, "B": 1, "C": 1, "modulus": 3},
        },
        "compatibility_notes": compatibility_notes(),
        "field_charges": {name: c.as_json() for name, c in charges.items()},
        "evaluated_terms": evaluated,
        "summary": summary,
        "passes": {
            "all_required_and_forbidden_terms_match": summary["all_pass"],
            "Yukawa_sector_allowed": all(
                row["pass"] for row in evaluated if row["sector"] in {"Yukawa", "Yukawa/Majorana", "Majorana_dim5"}
            ),
            "dangerous_terms_forbidden": all(
                row["pass"] for row in evaluated if row["sector"].startswith("forbidden")
            ),
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "combined_charge_table_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)

    print("Combined charge-table verification")
    print(f"  terms checked: {summary['terms_checked']}")
    print(f"  failures: {len(summary['failures'])}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
