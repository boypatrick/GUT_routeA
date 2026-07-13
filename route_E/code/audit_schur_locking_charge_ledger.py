#!/usr/bin/env python3
"""Audit a finite charge ledger for the polynomial Schur-locking sector.

The previous path-sum/determinant-locking audit showed that

    M = [[A,U],[V,G0 + V A^{-1} U]]

locks threshold determinants while leaving the AA inverse block sensitive to
open paths.  The polynomial completion introduces auxiliary fields B,C and
drivers X,Y,Z,

    W = Tr X(AB - I) + Tr Z(C - BU) + Tr Y(G - G0 - VC).

This script searches for a small finite grading, together with the usual
R(W)=2 assignment, that permits those terms while forbidding representative
dangerous lower-path and doublet-bridge terms.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


OUT = Path("output/schur_locking_charge_ledger")
FIELDS = (
    "A",
    "B",
    "C",
    "G",
    "G0",
    "U_D",
    "V_D",
    "U_L",
    "V_L",
    "I_A",
    "X",
    "Y",
    "Z",
    "H120",
    "H_u",
    "H_d",
    "M_mix",
)

R_CHARGE = {field: 0 for field in FIELDS}
R_CHARGE.update({"X": 2, "Y": 2, "Z": 2, "M_mix": 2})


@dataclass(frozen=True)
class Term:
    name: str
    fields: Tuple[str, ...]
    required: bool
    reason: str


TERMS: Tuple[Term, ...] = (
    Term("XAB", ("X", "A", "B"), True, "driver enforces AB=I_A"),
    Term("XI_A", ("X", "I_A"), True, "driver source for inverse identity"),
    Term("ZC", ("Z", "C"), True, "driver enforces C=BU_D"),
    Term("ZBU_D", ("Z", "B", "U_D"), True, "triplet bridge path"),
    Term("YG", ("Y", "G"), True, "driver source for sterile block"),
    Term("YG0", ("Y", "G0"), True, "unperturbed sterile block"),
    Term("YV_DC", ("Y", "V_D", "C"), True, "triplet return path"),
    Term("XA", ("X", "A"), False, "would force A=0 without inverse source"),
    Term("XB", ("X", "B"), False, "would force B=0 without inverse source"),
    Term("ZU_D", ("Z", "U_D"), False, "would collapse C=BU_D to C=U_D"),
    Term("ZBU_L", ("Z", "B", "U_L"), False, "doublet bridge analogue"),
    Term("YV_LC", ("Y", "V_L", "C"), False, "doublet return analogue"),
    Term("YV_D", ("Y", "V_D"), False, "open return tadpole"),
    Term("YC", ("Y", "C"), False, "sterile/source cross tadpole"),
    Term("M_H120_U_D", ("M_mix", "H120", "U_D"), False, "physical Higgs-triplet bridge mixing"),
    Term("M_H120_V_D", ("M_mix", "H120", "V_D"), False, "physical Higgs-triplet return mixing"),
    Term("M_H120_G", ("M_mix", "H120", "G"), False, "physical Higgs/source mass mixing"),
    Term("M_Hu_U_L", ("M_mix", "H_u", "U_L"), False, "doublet bridge mixes with MSSM H_u"),
    Term("M_Hd_V_L", ("M_mix", "H_d", "V_L"), False, "doublet bridge mixes with MSSM H_d"),
    Term("M_UD_UL", ("M_mix", "U_D", "U_L"), False, "triplet-doublet bridge leakage"),
    Term("M_VD_VL", ("M_mix", "V_D", "V_L"), False, "return-path doublet leakage"),
)


def mod_sum(charges: Dict[str, int], fields: Iterable[str], modulus: int) -> int:
    return sum(charges[f] for f in fields) % modulus


def r_sum(fields: Iterable[str]) -> int:
    return sum(R_CHARGE[f] for f in fields)


def is_allowed(
    term: Term,
    z_n: Dict[str, int],
    n: int,
    z2: Dict[str, int],
) -> bool:
    return (
        mod_sum(z_n, term.fields, n) == 0
        and mod_sum(z2, term.fields, 2) == 0
        and r_sum(term.fields) == 2
    )


def build_charges(
    n: int,
    a: int,
    b: int,
    u: int,
    g: int,
    h120: int,
    u_l: int,
    v_l: int,
) -> Dict[str, int]:
    """Solve the required Schur-locking charge equations over Z_n.

    The identity I_A is treated as a source spurion with the charge required
    by XAB.  In the minimal Z_2 solution it is neutral because q(A)+q(B)=0,
    but A and B themselves are nontrivially graded; that is what separates
    XAB from XA/XB and ZBU_D from ZU_D.
    """

    charges = {field: 0 for field in FIELDS}
    charges["A"] = a % n
    charges["B"] = b % n
    charges["U_D"] = u % n
    charges["G"] = g % n
    charges["G0"] = g % n
    charges["H120"] = h120 % n
    charges["U_L"] = u_l % n
    charges["V_L"] = v_l % n
    charges["H_u"] = 0
    charges["H_d"] = 0
    charges["M_mix"] = 0

    charges["X"] = (-(a + b)) % n
    charges["I_A"] = (a + b) % n
    charges["C"] = (b + u) % n
    charges["Z"] = (-(b + u)) % n
    charges["Y"] = (-g) % n
    charges["V_D"] = (g - charges["C"]) % n
    return charges


def term_rows(z_n: Dict[str, int], n: int, z2: Dict[str, int]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for term in TERMS:
        allowed = is_allowed(term, z_n, n, z2)
        rows.append(
            {
                "term": term.name,
                "fields": " ".join(term.fields),
                "required": term.required,
                "allowed": allowed,
                "zN_sum": mod_sum(z_n, term.fields, n),
                "z2_sum": mod_sum(z2, term.fields, 2),
                "R_sum": r_sum(term.fields),
                "status": "PASS" if allowed == term.required else "FAIL",
                "reason": term.reason,
            }
        )
    return rows


def valid_solution(z_n: Dict[str, int], n: int, z2: Dict[str, int]) -> bool:
    return all(row["status"] == "PASS" for row in term_rows(z_n, n, z2))


def score_solution(z_n: Dict[str, int], n: int, z2: Dict[str, int]) -> Tuple[int, int, int]:
    nonzero = sum(1 for f in FIELDS if z_n[f] % n != 0) + sum(1 for f in FIELDS if z2[f] % 2 != 0)
    bridge_weight = (z_n["U_D"] == 0) + (z_n["V_D"] == 0) + (z2["U_D"] == 0) + (z2["V_D"] == 0)
    return (n, bridge_weight, nonzero)


def find_solution() -> Tuple[int, Dict[str, int], Dict[str, int]]:
    candidates: List[Tuple[Tuple[int, int, int], int, Dict[str, int], Dict[str, int]]] = []
    for n in range(2, 9):
        for a in range(n):
            for b in range(n):
                for u in range(n):
                    for g in range(n):
                        for h120 in range(n):
                            for u_l in range(n):
                                for v_l in range(n):
                                    z_n = build_charges(n, a, b, u, g, h120, u_l, v_l)
                                    if z_n["B"] == 0:
                                        continue
                                    for a2 in range(2):
                                        for b2 in range(2):
                                            for u2 in range(2):
                                                for g2 in range(2):
                                                    for h120_2 in range(2):
                                                        for u_l2 in range(2):
                                                            for v_l2 in range(2):
                                                                z2 = build_charges(
                                                                    2,
                                                                    a2,
                                                                    b2,
                                                                    u2,
                                                                    g2,
                                                                    h120_2,
                                                                    u_l2,
                                                                    v_l2,
                                                                )
                                                                if valid_solution(z_n, n, z2):
                                                                    candidates.append(
                                                                        (score_solution(z_n, n, z2), n, z_n, z2)
                                                                    )
        if candidates:
            break
    if not candidates:
        raise RuntimeError("No finite Schur-locking charge ledger found up to Z_8 x Z_2.")
    candidates.sort(key=lambda item: item[0])
    _, n, z_n, z2 = candidates[0]
    return n, z_n, z2


def write_outputs(n: int, z_n: Dict[str, int], z2: Dict[str, int]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = term_rows(z_n, n, z2)

    with (OUT / "term_ledger.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    charge_rows = [
        {
            "field": field,
            f"Z{n}_lock": z_n[field],
            "Z2_path": z2[field],
            "R": R_CHARGE[field],
        }
        for field in FIELDS
    ]
    with (OUT / "charge_table.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(charge_rows[0].keys()))
        writer.writeheader()
        writer.writerows(charge_rows)

    summary = {
        "finite_group": f"Z_{n} x Z_2^path x U(1)_R",
        "all_required_allowed": all(row["allowed"] for row in rows if row["required"]),
        "all_forbidden_rejected": all(not row["allowed"] for row in rows if not row["required"]),
        "required_terms": [row for row in rows if row["required"]],
        "forbidden_terms": [row for row in rows if not row["required"]],
        "charge_table": charge_rows,
        "interpretation": {
            "inverse_block_grading": (
                "A and B carry a nontrivial inverse-block grading.  The source I_A "
                "has the charge of AB, which is neutral in the minimal Z_2 solution. "
                "This prevents the lower-path invariants XA, XB and ZU_D without "
                "changing the Schur F-term solution."
            ),
            "limitation": (
                "The ledger controls the Schur-locking auxiliary sector.  It does "
                "not by itself prove that every Spin(10) component of the 54/210 "
                "source sector is scalar-degenerate; that remains a spectrum audit."
            ),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = [
        "# Schur-locking charge ledger",
        "",
        f"Finite group: `Z_{n} x Z_2^path x U(1)_R`.",
        "",
        "## Charge table",
        "",
        f"| field | Z_{n} | Z_2^path | R |",
        "|---|---:|---:|---:|",
    ]
    for row in charge_rows:
        report.append(f"| {row['field']} | {row[f'Z{n}_lock']} | {row['Z2_path']} | {row['R']} |")
    report.extend(
        [
            "",
            "## Term ledger",
            "",
            "| term | required | allowed | Z_N | Z_2 | R | status |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in rows:
        report.append(
            f"| {row['term']} | {row['required']} | {row['allowed']} | "
            f"{row['zN_sum']} | {row['z2_sum']} | {row['R_sum']} | {row['status']} |"
        )
    report.extend(
        [
            "",
            "## Result",
            "",
            "- Required Schur-locking terms are all allowed.",
            "- Representative lower-path, physical-Higgs, and doublet-bridge leakages are all forbidden.",
            "- The nontrivial point is the inverse-block grading of `A,B`: if `B` is neutral, additive charges cannot distinguish `ZBU_D` from the lower-path term `ZU_D`.",
        ]
    )
    (OUT / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    n, z_n, z2 = find_solution()
    write_outputs(n, z_n, z2)
    rows = term_rows(z_n, n, z2)
    print(f"finite_group=Z_{n} x Z_2^path x U(1)_R")
    print(f"required_allowed={all(row['allowed'] for row in rows if row['required'])}")
    print(f"forbidden_rejected={all(not row['allowed'] for row in rows if not row['required'])}")
    print(f"outputs={OUT}")


if __name__ == "__main__":
    main()
