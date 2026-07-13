#!/usr/bin/env python3
"""Audit candidate UV mechanisms for the kappa54 locked branch.

No web input is used.  The script quantifies three possibilities:

1. quasi-extended Higgsing: exact mass locking is assumed by an extended
   gauge-Higgs relation, and we check the extra complete-multiplet threshold;
2. constrained/composite 54 order parameter: impose the spectral orbit
   S^2 - S - 6 I = 0 and verify that only the 24 Goldstone tangent directions
   remain;
3. conformal threshold sector: solve the anomalous-dimension coefficient
   needed to make lambda/g=sqrt(2) a fixed ratio despite the ordinary N=1
   obstruction.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "uv_locking_mechanisms"
KAPPA_CARDS = ROOT / "output" / "kappa54_global_scan" / "kappa54_benchmark_cards.json"


def symmetric_traceless_basis(n: int = 10) -> list[np.ndarray]:
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


def inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a.T @ b))


def constrained_orbit_audit() -> dict[str, Any]:
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    residual = s0 @ s0 - s0 - 6.0 * np.eye(10)
    domain = symmetric_traceless_basis()
    codomain = symmetric_basis()

    # Linearized spectral constraint: delta(S^2-S-6I)=S0 dS+dS S0-dS.
    jac = np.zeros((len(codomain), len(domain)))
    for a, e in enumerate(domain):
        image = s0 @ e + e @ s0 - e
        for b, c in enumerate(codomain):
            jac[b, a] = inner(c, image)
    singular_values = np.linalg.svd(jac, compute_uv=False)
    rank = int(np.sum(singular_values > 1.0e-10))
    nullity = len(domain) - rank

    # Orbit tangent from commutators [T,S0].
    tangents = []
    for t in so10_basis():
        tangent = t @ s0 - s0 @ t
        tangents.append(tangent)
    gram = np.array([[inner(a, b) for b in tangents] for a in tangents])
    tangent_eigs = np.linalg.eigvalsh(gram)
    tangent_rank = int(np.sum(tangent_eigs > 1.0e-10))

    # Classify the linearized constraint by the three block types.
    coeffs = {
        "cc_symmetric_block": -2.0 + -2.0 - 1.0,
        "ww_symmetric_block": 3.0 + 3.0 - 1.0,
        "cw_cross_block": -2.0 + 3.0 - 1.0,
    }
    multiplicities = {
        "cc_symmetric": 6 * 7 // 2,
        "ww_symmetric": 4 * 5 // 2,
        "trace_relation_removed": -1,
        "cw_cross_goldstone": 6 * 4,
    }

    return {
        "minimal_polynomial": "S^2 - S - 6 I = 0",
        "S0_trace": float(np.trace(s0)),
        "S0_trace_square": float(np.trace(s0 @ s0)),
        "polynomial_residual_norm": float(np.linalg.norm(residual)),
        "constraint_jacobian_rank": rank,
        "constraint_jacobian_nullity": nullity,
        "orbit_tangent_rank": tangent_rank,
        "normal_modes_removed": rank,
        "goldstone_modes": nullity,
        "block_linear_coefficients": coeffs,
        "block_multiplicities": multiplicities,
        "nonuniversal_threshold_vector": [0.0, 0.0, 0.0],
        "passes": residual.dot(residual).sum() < 1.0e-20 and rank == 30 and nullity == 24 and tangent_rank == 24,
        "singular_values": singular_values.tolist(),
    }


def ordinary_fixedline_rows(g_value: float) -> list[dict[str, Any]]:
    scenarios = [
        ("single_54_only", -12.0, "ordinary elementary 54 cubic"),
        ("54_plus_three_45", 12.0, "low-index projector mediator"),
        ("54_plus_minimal_flavor_higgs", 58.0, "minimal flavor/Higgs sector"),
        ("54_plus_projector_plus_flavor", 82.0, "low-index projector plus flavor/Higgs"),
        ("large_54_210_alignment_tower", 184.0, "large propagating source/driver tower"),
    ]
    rows: list[dict[str, Any]] = []
    target_beta_lambda_coeff = (84.0 / 5.0) * 2.0 - 60.0
    loop = 16.0 * math.pi * math.pi
    for name, b10, note in scenarios:
        extra_coeff = b10 - target_beta_lambda_coeff
        gamma_cft_coeff = extra_coeff / 3.0
        gamma_cft_value = gamma_cft_coeff * g_value * g_value / loop
        fixed_sq = (60.0 + b10) / (84.0 / 5.0)
        rows.append(
            {
                "scenario": name,
                "b10": b10,
                "ordinary_fixed_lambda_over_g": math.sqrt(fixed_sq) if fixed_sq > 0 else float("nan"),
                "required_extra_beta_lambda_coefficient": extra_coeff,
                "required_gamma_cft_coefficient": gamma_cft_coeff,
                "required_gamma_cft_at_R200_g": gamma_cft_value,
                "note": note,
            }
        )
    return rows


def mechanism_scorecard(constrained: dict[str, Any], gamma_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    qes_extra_index = 20.0  # diagonal 10 x 10 link = 1 + 45 + 54 under the diagonal SO(10)
    pure_qes_b10 = qes_extra_index - 24.0
    single_54_gamma = next(row for row in gamma_rows if row["scenario"] == "single_54_only")
    return [
        {
            "mechanism": "quasi_extended_higgsing_link",
            "core_equation": "N=2-like gauge-Higgs locking gives M_link=M_V, hence kappa54=1 if the 54 is a link fragment.",
            "lambda_over_g_prediction": math.sqrt(2.0),
            "evades_b10_obstruction": True,
            "nonuniversal_threshold": 0.0,
            "extra_dynkin_index": qes_extra_index,
            "pure_link_b10": pure_qes_b10,
            "status": "conditional",
            "main_caveat": "An explicit Spin(10) link construction must show that the 54 fragment, not only adjoint fields, inherits the extended mass locking.",
        },
        {
            "mechanism": "constrained_composite_54_orbit",
            "core_equation": "S^2-S-6I=0, Tr S=0; the only fluctuations are the 24 Spin(10)/PS orbit tangents.",
            "lambda_over_g_prediction": "not needed",
            "evades_b10_obstruction": True,
            "nonuniversal_threshold": 0.0,
            "constraint_rank": constrained["constraint_jacobian_rank"],
            "goldstone_modes": constrained["goldstone_modes"],
            "status": "best_current_route",
            "main_caveat": "It is a constrained/composite order parameter rather than a fully elementary renormalizable 54_H cubic.",
        },
        {
            "mechanism": "conformal_threshold_fixed_line",
            "core_equation": "Add gamma_CFT so beta_lambda/lambda=beta_g/g at lambda/g=sqrt(2).",
            "lambda_over_g_prediction": math.sqrt(2.0),
            "evades_b10_obstruction": True,
            "required_gamma_cft_at_R200_g_single54": single_54_gamma["required_gamma_cft_at_R200_g"],
            "required_extra_beta_coeff_single54": single_54_gamma["required_extra_beta_lambda_coefficient"],
            "nonuniversal_threshold": "model dependent",
            "status": "possible_but_unbuilt",
            "main_caveat": "The anomalous dimension must track g^2 and the charged CFT threshold must be added to the RGE/proton scan.",
        },
    ]


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


def write_singular_values(path: Path, values: list[float]) -> None:
    rows = [{"index": i, "singular_value": value} for i, value in enumerate(values)]
    write_csv(path, rows)


def write_report(payload: dict[str, Any]) -> None:
    constrained = payload["constrained_orbit"]
    lines: list[str] = []
    lines.append("# UV locking mechanism audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Recommended route")
    lines.append("")
    lines.append(
        "The constrained/composite 54 orbit is the cleanest current route: it evades the ordinary "
        "N=1 fixed-line obstruction by removing the propagating non-Goldstone 54 fragments."
    )
    lines.append("")
    lines.append("## Constrained 54 orbit")
    lines.append("")
    lines.append(f"Polynomial residual norm: `{constrained['polynomial_residual_norm']:.3e}`")
    lines.append(f"Constraint Jacobian rank/nullity: `{constrained['constraint_jacobian_rank']}/{constrained['constraint_jacobian_nullity']}`")
    lines.append(f"Orbit tangent rank: `{constrained['orbit_tangent_rank']}`")
    lines.append(f"Passes: `{constrained['passes']}`")
    lines.append("")
    lines.append("Linearized block coefficients for delta(S^2-S-6I):")
    for key, value in constrained["block_linear_coefficients"].items():
        lines.append(f"- `{key}`: `{value:.1f}`")
    lines.append("")
    lines.append("## Conformal anomalous-dimension requirement")
    lines.append("")
    lines.append("| scenario | b10 | ordinary fixed lambda/g | required extra coeff | required gamma_CFT at R200 |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in payload["conformal_required_gamma"]:
        lines.append(
            f"| {row['scenario']} | {row['b10']:.1f} | "
            f"{row['ordinary_fixed_lambda_over_g']:.6f} | "
            f"{row['required_extra_beta_lambda_coefficient']:.6f} | "
            f"{row['required_gamma_cft_at_R200_g']:.6e} |"
        )
    lines.append("")
    lines.append("## Scorecard")
    lines.append("")
    lines.append("| mechanism | status | nonuniversal threshold | main caveat |")
    lines.append("|---|---|---:|---|")
    for row in payload["scorecard"]:
        lines.append(
            f"| {row['mechanism']} | {row['status']} | {row['nonuniversal_threshold']} | {row['main_caveat']} |"
        )
    lines.append("")
    (OUT / "uv_locking_mechanism_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cards = json.loads(KAPPA_CARDS.read_text(encoding="utf-8"))
    r200 = cards["R200_kappa54_locked"]
    g_r200 = float(r200["g_G"])
    constrained = constrained_orbit_audit()
    gamma_rows = ordinary_fixedline_rows(g_r200)
    scorecard = mechanism_scorecard(constrained, gamma_rows)
    payload = {
        "note": "No web lookup used. Audit of UV mechanisms for kappa54 locking.",
        "locked_R200_inputs": {
            "g_G": g_r200,
            "lambda_locked": r200["lambda_value"],
            "kappa_54": r200["kappa_54"],
            "M_Sigma3_GeV": r200["best"]["M_Sigma3_GeV"],
            "tau_dim6_years": r200["best"]["tau_dim6_years"],
        },
        "constrained_orbit": constrained,
        "conformal_required_gamma": gamma_rows,
        "scorecard": scorecard,
        "verdict": {
            "recommended_route": "constrained_composite_54_orbit",
            "reason": (
                "It gives a first-principles algebraic orbit lock with zero nonuniversal "
                "single-54 threshold and does not rely on an impossible ordinary N=1 fixed line."
            ),
        },
    }
    write_singular_values(OUT / "constrained_54_jacobian_spectrum.csv", constrained["singular_values"])
    write_csv(OUT / "conformal_required_gamma.csv", gamma_rows)
    write_csv(OUT / "mechanism_scorecard.csv", scorecard)
    (OUT / "uv_locking_mechanism_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("UV locking mechanism audit complete")
    print(f"  constrained orbit passes={constrained['passes']}")
    print(f"  constrained rank/nullity={constrained['constraint_jacobian_rank']}/{constrained['constraint_jacobian_nullity']}")
    print(f"  recommended={payload['verdict']['recommended_route']}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
