#!/usr/bin/env python3
"""Component Hessian audit for the constrained 54 orbit superpotential.

The goal is to compare multiplier choices for

    W_orbit = <Xi, S^2 - S - 6 I>

around S0=diag(-2^6,3^4).  The script is local and self-contained: no web
lookup is used.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "orbit_superpotential_hessian"


def inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a.T @ b))


def symmetric_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        m = np.zeros((n, n), dtype=float)
        m[i, i] = 1.0
        basis.append(m)
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(m)
    assert len(basis) == 55
    return basis


def traceless_symmetric_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(m)
    for k in range(n - 1):
        m = np.zeros((n, n), dtype=float)
        coeff = 1.0 / math.sqrt((k + 1) * (k + 2))
        for i in range(k + 1):
            m[i, i] = coeff
        m[k + 1, k + 1] = -(k + 1) * coeff
        basis.append(m)
    assert len(basis) == 54
    return basis


def so10_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = -1.0 / math.sqrt(2.0)
            basis.append(m)
    assert len(basis) == 45
    return basis


def linearized_constraint(s0: np.ndarray, x: np.ndarray) -> np.ndarray:
    return s0 @ x + x @ s0 - x


def map_matrix(codomain: list[np.ndarray], domain: list[np.ndarray], s0: np.ndarray) -> np.ndarray:
    j = np.zeros((len(codomain), len(domain)))
    for a, e in enumerate(domain):
        image = linearized_constraint(s0, e)
        for b, c in enumerate(codomain):
            j[b, a] = inner(c, image)
    return j


def hessian_stats(j: np.ndarray, xi_dim: int, s_dim: int) -> dict[str, Any]:
    singular = np.linalg.svd(j, compute_uv=False)
    rank = int(np.sum(singular > 1.0e-10))
    total_dim = xi_dim + s_dim
    h_rank = 2 * rank
    return {
        "constraint_rank": rank,
        "constraint_nullity_S": s_dim - rank,
        "multiplier_nullity": xi_dim - rank,
        "hessian_dimension": total_dim,
        "hessian_rank": h_rank,
        "hessian_zero_modes": total_dim - h_rank,
        "nonzero_singular_min": float(min(s for s in singular if s > 1.0e-10)),
        "nonzero_singular_max": float(max(singular)),
        "singular_values": singular.tolist(),
    }


def tangent_rank(s0: np.ndarray) -> int:
    tangents = [t @ s0 - s0 @ t for t in so10_basis()]
    gram = np.array([[inner(a, b) for b in tangents] for a in tangents])
    return int(np.sum(np.linalg.eigvalsh(gram) > 1.0e-10))


def audit() -> dict[str, Any]:
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    s_basis = traceless_symmetric_basis()
    sym_basis = symmetric_basis()
    ts_basis = traceless_symmetric_basis()

    j_full = map_matrix(sym_basis, s_basis, s0)
    j_projected54 = map_matrix(ts_basis, s_basis, s0)

    # The normal-bundle multiplier keeps only the nonzero image directions.
    # SVD gives an orthonormal basis for im(J) without choosing a PS component
    # coordinate system by hand.
    u, singular, vh = np.linalg.svd(j_full, full_matrices=True)
    rank = int(np.sum(singular > 1.0e-10))
    j_normal = np.diag(singular[:rank]) @ vh[:rank, :]

    b54_phys = np.array([34.0 / 5.0, 6.0, 8.0])
    cases = [
        (
            "unprojected_symmetric_Xi55",
            j_full,
            len(sym_basis),
            "dynamical unprojected Xi leaves 25 massless multiplier modes",
            "reject_as_dynamical; acceptable only if Xi is purely auxiliary",
        ),
        (
            "projected_traceless_Xi54",
            j_projected54,
            len(ts_basis),
            "projected Xi removes the identity multiplier but still leaves 24 multiplier zero modes",
            "reject_as_dynamical; acceptable only if Xi is auxiliary or further gauge-fixed",
        ),
        (
            "normal_bundle_Xi30",
            j_normal,
            rank,
            "normal-bundle Xi pairs all 30 normal S directions and leaves only 24 Goldstone modes",
            "best_action_level_candidate",
        ),
    ]
    rows = []
    for name, j, xi_dim, comment, verdict in cases:
        stats = hessian_stats(j, xi_dim=xi_dim, s_dim=len(s_basis))
        dynamical_threshold = (2.0 * b54_phys).tolist() if name == "normal_bundle_Xi30" else "model dependent"
        if verdict.startswith("reject"):
            dynamical_threshold = "bad unless multiplier zero modes are auxiliary/removed"
        rows.append(
            {
                "case": name,
                **{k: v for k, v in stats.items() if k != "singular_values"},
                "goldstone_modes_expected": 24,
                "dynamical_threshold_vector_if_propagating": dynamical_threshold,
                "comment": comment,
                "verdict": verdict,
            }
        )
    return {
        "S0": [-2.0] * 6 + [3.0] * 4,
        "S_dimension": len(s_basis),
        "orbit_tangent_rank": tangent_rank(s0),
        "cases": rows,
        "singular_spectra": {
            "unprojected_symmetric_Xi55": hessian_stats(j_full, 55, 54)["singular_values"],
            "projected_traceless_Xi54": hessian_stats(j_projected54, 54, 54)["singular_values"],
            "normal_bundle_Xi30": hessian_stats(j_normal, 30, 54)["singular_values"],
        },
        "threshold_vectors": {
            "single_54_physical": b54_phys.tolist(),
            "dynamical_normal_pair": (2.0 * b54_phys).tolist(),
            "auxiliary_normal_bundle": [0.0, 0.0, 0.0],
        },
        "verdict": {
            "recommended": "normal_bundle_Xi30_as_auxiliary_or_composite_constraint",
            "reason": "It leaves exactly the 24 Goldstone directions and contributes no physical single-54 threshold if the multiplier is auxiliary/composite.",
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_spectra(payload: dict[str, Any]) -> None:
    for name, spectrum in payload["singular_spectra"].items():
        rows = [{"index": i, "singular_value": value} for i, value in enumerate(spectrum)]
        write_csv(OUT / f"{name}_singular_values.csv", rows)


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Orbit superpotential Hessian audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "The normal-bundle multiplier is the only clean action-level candidate: "
        "it pairs the 30 normal 54 directions and leaves exactly the 24 Goldstone directions."
    )
    lines.append("")
    lines.append("## Cases")
    lines.append("")
    lines.append("| case | Xi dim | rank(J) | S nullity | Xi nullity | Hessian zeros | verdict |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for row in payload["cases"]:
        lines.append(
            f"| {row['case']} | {row['hessian_dimension'] - payload['S_dimension']} | "
            f"{row['constraint_rank']} | {row['constraint_nullity_S']} | "
            f"{row['multiplier_nullity']} | {row['hessian_zero_modes']} | {row['verdict']} |"
        )
    lines.append("")
    lines.append("## Threshold vectors")
    lines.append("")
    lines.append(f"single_54_physical = `{payload['threshold_vectors']['single_54_physical']}`")
    lines.append(f"dynamical_normal_pair = `{payload['threshold_vectors']['dynamical_normal_pair']}`")
    lines.append(f"auxiliary_normal_bundle = `{payload['threshold_vectors']['auxiliary_normal_bundle']}`")
    lines.append("")
    (OUT / "orbit_superpotential_hessian_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = audit()
    write_csv(OUT / "orbit_superpotential_hessian_cases.csv", payload["cases"])
    write_spectra(payload)
    (OUT / "orbit_superpotential_hessian_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Orbit superpotential Hessian audit complete")
    print(f"  orbit tangent rank={payload['orbit_tangent_rank']}")
    for row in payload["cases"]:
        print(
            f"  {row['case']}: rank={row['constraint_rank']} "
            f"S_null={row['constraint_nullity_S']} Xi_null={row['multiplier_nullity']} "
            f"H_zero={row['hessian_zero_modes']}"
        )
    print(f"  recommended={payload['verdict']['recommended']}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
