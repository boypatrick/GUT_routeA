#!/usr/bin/env python3
"""Low-Dynkin-index projector/filter enumerator.

No web lookup is used.  This script tests the UV-perturbativity escape route
suggested by the Landau-pole audit: can the X_(6,2,2) projector be reproduced
with a small propagating mediator sector if the large 54/210 alignment tensors
are treated as spurionic backgrounds rather than propagating multiplets?

The four relevant 45 fragments carry

    Sigma_L : F54= 2,    D210=+1,  P_X=0
    Sigma_R : F54= 2,    D210=-1,  P_X=0
    Sigma_8 : F54=-4/3, D210= 0,  P_X=0
    X_622   : F54= 1/3, D210= 0,  P_X=1

The original projector is

    P_X = -(9/25)(F54-2)(F54+4/3).

This is already a quadratic polynomial in the 54 spurion alone.  The important
question is therefore not algebraic existence but UV interpretation: whether
the 54/210 fields that generate the Clebsch filter must be propagating in the
SO(10) interval above M_G.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "low_index_projector_filters"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"

T_INDEX = {
    "10": 1.0,
    "16": 2.0,
    "16bar": 2.0,
    "45": 8.0,
    "54": 12.0,
    "210": 56.0,
}
GAUGE_CONTRIBUTION = 24.0

FRAGMENTS = [
    {"name": "Sigma_L", "F54": 2.0, "D210": 1.0, "target": 0.0},
    {"name": "Sigma_R", "F54": 2.0, "D210": -1.0, "target": 0.0},
    {"name": "Sigma8", "F54": -4.0 / 3.0, "D210": 0.0, "target": 0.0},
    {"name": "X_622", "F54": 1.0 / 3.0, "D210": 0.0, "target": 1.0},
]


def basis_value(label: str, f: float, d: float) -> float:
    return {
        "1": 1.0,
        "F": f,
        "D": d,
        "F2": f * f,
        "FD": f * d,
        "D2": d * d,
        "PX_exact": -(9.0 / 25.0) * (f - 2.0) * (f + 4.0 / 3.0),
    }[label]


CANDIDATES = [
    {
        "name": "constant_only",
        "basis": ["1"],
        "propagating_counts": {"45": 3},
        "spurions": [],
    },
    {
        "name": "linear_F_spurion",
        "basis": ["1", "F"],
        "propagating_counts": {"45": 3},
        "spurions": ["F54"],
    },
    {
        "name": "quadratic_F_spurion",
        "basis": ["1", "F", "F2"],
        "propagating_counts": {"45": 3},
        "spurions": ["F54"],
    },
    {
        "name": "quadratic_F_with_one_propagating_54",
        "basis": ["1", "F", "F2"],
        "propagating_counts": {"45": 3, "54": 1},
        "spurions": [],
    },
    {
        "name": "quadratic_FD_spurions",
        "basis": ["1", "F", "D", "F2", "FD", "D2"],
        "propagating_counts": {"45": 3},
        "spurions": ["F54", "D210"],
    },
    {
        "name": "quadratic_FD_with_one_54_one_210",
        "basis": ["1", "F", "D", "F2", "FD", "D2"],
        "propagating_counts": {"45": 3, "54": 1, "210": 1},
        "spurions": [],
    },
    {
        "name": "exact_PX_spurionized_operator",
        "basis": ["PX_exact"],
        "propagating_counts": {"45": 3},
        "spurions": ["F54"],
    },
    {
        "name": "documented_large_alignment_sector",
        "basis": ["PX_exact"],
        "propagating_counts": {"45": 3, "54": 6, "210": 2},
        "spurions": [],
    },
]


def load_vacuum_rows() -> list[dict[str, Any]]:
    payload = json.loads(VACUUM.read_text(encoding="utf-8"))
    return payload["replay_rows"]


def sum_dynkin(counts: dict[str, int]) -> float:
    return float(sum(T_INDEX[rep] * count for rep, count in counts.items()))


def landau_ratio(alpha_inv: float, b10: float) -> float:
    if b10 <= 0.0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def fit_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    basis = candidate["basis"]
    matrix = np.array(
        [[basis_value(label, frag["F54"], frag["D210"]) for label in basis] for frag in FRAGMENTS],
        dtype=float,
    )
    target = np.array([frag["target"] for frag in FRAGMENTS], dtype=float)
    coeffs, *_ = np.linalg.lstsq(matrix, target, rcond=None)
    predicted = matrix @ coeffs
    residual = predicted - target
    rank = int(np.linalg.matrix_rank(matrix))
    return {
        "basis": basis,
        "rank": rank,
        "coefficients": {label: float(value) for label, value in zip(basis, coeffs)},
        "predicted": {frag["name"]: float(value) for frag, value in zip(FRAGMENTS, predicted)},
        "max_abs_projector_error": float(np.max(np.abs(residual))),
        "residual_l2": float(np.linalg.norm(residual)),
        "exact_projector": bool(np.max(np.abs(residual)) < 1.0e-12),
    }


def audit_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    fit = fit_candidate(candidate)
    sum_t = sum_dynkin(candidate["propagating_counts"])
    b10 = sum_t - GAUGE_CONTRIBUTION
    landau_rows = []
    for vac in load_vacuum_rows():
        alpha_inv = float(vac["alphaG_inv"])
        r_target = float(vac["R"])
        ratio_lp = landau_ratio(alpha_inv, b10)
        alpha_at_r = alpha_inv - b10 * math.log(r_target) / (2.0 * math.pi)
        landau_rows.append(
            {
                "R": r_target,
                "alphaG_inv_MG": alpha_inv,
                "alpha_inv_at_RMG": alpha_at_r,
                "landau_ratio": ratio_lp,
                "passes": alpha_at_r > 0.0,
            }
        )
    exact = fit["exact_projector"]
    landau_r50 = next(row for row in landau_rows if row["R"] == 50.0)["passes"]
    landau_r200 = next(row for row in landau_rows if row["R"] == 200.0)["passes"]
    sigma_constraints = (
        abs(fit["predicted"]["Sigma_L"]) < 1.0e-12
        and abs(fit["predicted"]["Sigma_R"]) < 1.0e-12
        and abs(fit["predicted"]["Sigma8"]) < 1.0e-12
        and abs(fit["predicted"]["X_622"] - 1.0) < 1.0e-12
    )
    return {
        "name": candidate["name"],
        "basis": candidate["basis"],
        "spurions": candidate["spurions"],
        "propagating_counts": candidate["propagating_counts"],
        "sum_T": sum_t,
        "b10": b10,
        **fit,
        "landau_rows": landau_rows,
        "passes_R50": bool(exact and sigma_constraints and landau_r50),
        "passes_R200": bool(exact and sigma_constraints and landau_r200),
        "interpretation": interpretation(candidate, fit, b10, landau_r200),
    }


def interpretation(candidate: dict[str, Any], fit: dict[str, Any], b10: float, landau_r200: bool) -> str:
    if not fit["exact_projector"]:
        return "reject: does not reproduce the P_X eigenvalue target."
    if not landau_r200:
        return "reject as high-R propagating completion: Landau pole below R=200 M_G."
    if candidate["spurions"]:
        return "passes as a spurionized low-index filter, not as a fully dynamical source-sector derivation."
    return "passes as a propagating low-index candidate on this four-fragment projector test."


def build_payload() -> dict[str, Any]:
    rows = [audit_candidate(candidate) for candidate in CANDIDATES]
    best = [row for row in rows if row["passes_R200"]]
    return {
        "note": "No web lookup used. Low-index projector/filter enumeration.",
        "target_fragments": FRAGMENTS,
        "rows": rows,
        "verdict": {
            "has_R200_candidate": bool(best),
            "R200_candidates": [row["name"] for row in best],
            "minimal_escape": (
                "quadratic_F_spurion or exact_PX_spurionized_operator reproduces the projector "
                "with only three propagating 45 fields counted in b10. This is UV-safe only if "
                "the F54 source is treated as a fixed spurion/constrained background."
            ),
            "no_fully_dynamic_large_sector": (
                "The documented large 54/210 alignment sector remains rejected as a high-R "
                "propagating completion by the Landau audit."
            ),
        },
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "name",
        "basis",
        "spurions",
        "propagating_counts",
        "sum_T",
        "b10",
        "rank",
        "max_abs_projector_error",
        "passes_R50",
        "passes_R200",
        "interpretation",
    ]
    with (OUT / "low_index_projector_filters.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Low-index projector/filter enumerator")
    lines.append("")
    lines.append("No web lookup was used.  The target is the four-fragment projector")
    lines.append("`P_X=(0,0,0,1)` on `Sigma_L,Sigma_R,Sigma8,X_622`.")
    lines.append("")
    lines.append("## Projector algebra")
    lines.append("")
    lines.append("The exact filter is")
    lines.append("")
    lines.append("```text")
    lines.append("P_X = -(9/25)(F54-2)(F54+4/3).")
    lines.append("```")
    lines.append("")
    lines.append("It uses only the 54 Clebsch eigenvalue as a spurionic background.")
    lines.append("")
    lines.append("## Candidate table")
    lines.append("")
    lines.append("| candidate | basis | sum T | b10 | error | R=50 | R=200 | interpretation |")
    lines.append("|---|---|---:|---:|---:|---|---|---|")
    for row in payload["rows"]:
        lines.append(
            f"| `{row['name']}` | `{','.join(row['basis'])}` | {row['sum_T']:.0f} | {row['b10']:.0f} | "
            f"{row['max_abs_projector_error']:.3e} | {row['passes_R50']} | {row['passes_R200']} | "
            f"{row['interpretation']} |"
        )
    lines.append("")
    lines.append("## Key conclusion")
    lines.append("")
    lines.append("The projector itself is not the Landau problem: it is exactly reproduced by")
    lines.append("a quadratic `F54` filter.  The Landau problem comes from making the full")
    lines.append("54/210 alignment/source sector a propagating chiral sector above `M_G`.")
    lines.append("")
    lines.append("Therefore the high-R branch has a consistent low-index interpretation only")
    lines.append("if the large alignment tensors are fixed spurions or constrained background")
    lines.append("fields.  A fully dynamical large-representation derivation still fails the")
    lines.append("one-loop UV perturbativity test.")
    (OUT / "low_index_projector_filters_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    (OUT / "low_index_projector_filters_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(payload["rows"])
    write_report(payload)
    print("Low-index projector/filter enumerator")
    for row in payload["rows"]:
        if row["passes_R200"] or row["name"] == "documented_large_alignment_sector":
            print(
                f"  {row['name']}: error={row['max_abs_projector_error']:.3e}, "
                f"sumT={row['sum_T']:.0f}, b10={row['b10']:.0f}, R200={row['passes_R200']}"
            )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
