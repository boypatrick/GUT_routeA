#!/usr/bin/env python3
"""Vector multiplet mass normalization for the single-54_H vacuum.

No web lookup is used.  This audit fixes the normalization question left by
the single-54_H Hessian scan:

    kappa_54 = 5 |m| / M_V.

The conventions are the same as in the Hessian script:

* S is a complex traceless-symmetric 10x10 chiral field.
* K = Tr(S^\dagger S).
* S0 = diag(-2^6, 3^4), and <S> = v S0 with v=-m/lambda.
* SO(10) vector generators are normalized as Tr_fund(T_a T_b)=delta_ab,
  i.e. T(vector)=1 and C2(adj)=8.

For a broken generator T, the kinetic term gives

    L_mass = g^2 A_mu^2 ||T . <S>||^2
           = (1/2) M_V^2 A_mu^2,

therefore M_V^2 = 2 g^2 ||T . <S>||^2.  For every mixed 6x4 broken generator,
||T . S0|| = 5, hence

    M_V = 5 sqrt(2) g |v| = 5 sqrt(2) g |m|/|lambda|,
    kappa_54 = 5 |m| / M_V = |lambda|/(sqrt(2) g).
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "single_54_vector_mass"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
HESSIAN = ROOT / "output" / "single_54_hessian" / "single_54_hessian_summary.json"

S0 = np.diag([-2.0] * 6 + [3.0] * 4)
N = 10


def normalized_so_generator(i: int, j: int) -> np.ndarray:
    """Real antisymmetric generator with Tr(T^T T)=1."""
    mat = np.zeros((N, N), dtype=float)
    mat[i, j] = 1.0 / math.sqrt(2.0)
    mat[j, i] = -1.0 / math.sqrt(2.0)
    return mat


def tangent(generator: np.ndarray, source: np.ndarray) -> np.ndarray:
    # For S -> R S R^T and antisymmetric generator A:
    # delta S = A S + S A^T = A S - S A = [A,S].
    return generator @ source - source @ generator


def generator_audit_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for i in range(N):
        for j in range(i + 1, N):
            gen = normalized_so_generator(i, j)
            tan = tangent(gen, S0)
            gen_norm = float(np.trace(gen.T @ gen))
            tan_norm2 = float(np.trace(tan.T @ tan))
            if i < 6 and j < 6:
                block = "so6_stabilizer"
            elif i >= 6 and j >= 6:
                block = "so4_stabilizer"
            else:
                block = "broken_mixed"
            rows.append(
                {
                    "i": i + 1,
                    "j": j + 1,
                    "block": block,
                    "generator_norm_TrTT": gen_norm,
                    "tangent_norm2_Tr": tan_norm2,
                    "tangent_norm": math.sqrt(max(tan_norm2, 0.0)),
                    "mass_coefficient_MV2_over_g2v2": 2.0 * tan_norm2,
                    "mass_coefficient_MV_over_gv": math.sqrt(2.0 * tan_norm2),
                }
            )
    return rows


def load_vacuum() -> dict[str, Any]:
    return json.loads(VACUUM.read_text(encoding="utf-8"))


def load_hessian() -> dict[str, Any]:
    return json.loads(HESSIAN.read_text(encoding="utf-8"))


def coupling_rows() -> list[dict[str, Any]]:
    payload = load_vacuum()
    rows: list[dict[str, Any]] = []
    for row in payload["replay_rows"]:
        alpha_inv = float(row["alphaG_inv"])
        g = math.sqrt(4.0 * math.pi / alpha_inv)
        lambda_for_kappa1 = math.sqrt(2.0) * g
        rows.append(
            {
                "R": float(row["R"]),
                "alphaG_inv": alpha_inv,
                "g_G": g,
                "lambda_for_kappa54_eq_1": lambda_for_kappa1,
                "lambda_over_g_for_kappa54_eq_1": math.sqrt(2.0),
                "lambda_over_4pi_sqrt": lambda_for_kappa1 / math.sqrt(4.0 * math.pi),
            }
        )
    return rows


def kappa_to_coupling_rows() -> list[dict[str, Any]]:
    hessian = load_hessian()
    r200 = next(row for row in coupling_rows() if abs(row["R"] - 200.0) < 1.0e-12)
    g = r200["g_G"]
    rows = []
    for row in hessian["kappa_replay_rows"]:
        kappa = float(row["kappa_54"])
        lam_over_g = math.sqrt(2.0) * kappa
        lam = lam_over_g * g
        rows.append(
            {
                "kappa_54": kappa,
                "lambda_over_g": lam_over_g,
                "lambda_at_R200_g": lam,
                "safe_points": int(row["safe_points"]),
                "alphaG_inv": float(row["best"]["alphaG_inv"]),
                "M_Sigma3_GeV": float(row["best"]["M_Sigma3_GeV"]),
                "projected_l2_delta54": float(row["projected_l2_delta54"]),
                "perturbative_lambda_lt_sqrt4pi": lam < math.sqrt(4.0 * math.pi),
            }
        )
    return rows


def tuning_budget() -> dict[str, Any]:
    hessian = load_hessian()
    unit = float(hessian["threshold_tolerances"]["unit_projected_l2_per_abs_log_kappa"])
    budgets = [5.0e-4, 1.0e-3, 1.0e-2, 5.0e-2, 1.0e-1]
    return {
        "unit_projected_l2_per_abs_log_kappa": unit,
        "budgets": {
            f"{budget:.1e}": {
                "abs_log_kappa_max": budget / unit,
                "kappa_range": [math.exp(-budget / unit), math.exp(budget / unit)],
                "lambda_over_g_range": [
                    math.sqrt(2.0) * math.exp(-budget / unit),
                    math.sqrt(2.0) * math.exp(budget / unit),
                ],
                "fractional_half_width_around_kappa1": math.exp(budget / unit) - 1.0,
            }
            for budget in budgets
        },
    }


def build_payload() -> dict[str, Any]:
    gen_rows = generator_audit_rows()
    broken = [row for row in gen_rows if row["block"] == "broken_mixed"]
    stab = [row for row in gen_rows if row["block"] != "broken_mixed"]
    tangent_norms = sorted({round(row["tangent_norm"], 12) for row in gen_rows})
    coupling = coupling_rows()
    kappa_rows = kappa_to_coupling_rows()
    budget = tuning_budget()
    return {
        "note": "No web lookup used. Single-54_H vector mass normalization audit.",
        "conventions": {
            "Kahler": "K = Tr(S^dagger S)",
            "generator_normalization": "Tr_fund(T_a T_b)=delta_ab; T(vector)=1",
            "mass_term_matching": "g^2 A^2 ||T<S>||^2 = (1/2) M_V^2 A^2",
        },
        "generator_summary": {
            "total_generators": len(gen_rows),
            "broken_generators": len(broken),
            "stabilizer_generators": len(stab),
            "unique_tangent_norms": tangent_norms,
            "broken_tangent_norm": broken[0]["tangent_norm"],
            "broken_MV_over_gv": broken[0]["mass_coefficient_MV_over_gv"],
            "passes": len(broken) == 24
            and len(stab) == 21
            and all(abs(row["tangent_norm"] - 5.0) < 1.0e-12 for row in broken)
            and all(abs(row["tangent_norm"]) < 1.0e-12 for row in stab),
        },
        "formula": {
            "MV": "5 sqrt(2) g |v|",
            "v": "|m|/|lambda|",
            "M54_phys": "5 |m|",
            "kappa54": "|lambda|/(sqrt(2) g)",
            "kappa54_eq_1": "|lambda|/g = sqrt(2)",
            "real_scalar_half_trace_convention_note": (
                "If instead one uses a purely real scalar kinetic term (1/2)Tr(DS DS), "
                "the sqrt(2) is absorbed and kappa54=|lambda|/g. The SUSY chiral "
                "normalization used here is K=Tr(S^dagger S)."
            ),
        },
        "coupling_rows": coupling,
        "kappa_to_coupling_rows_R200": kappa_rows,
        "tuning_budget": budget,
        "verdict": {
            "kappa54_eq_1_is_mild": True,
            "lambda_over_g_required": math.sqrt(2.0),
            "R200_g": next(row["g_G"] for row in coupling if abs(row["R"] - 200.0) < 1.0e-12),
            "R200_lambda_required": next(
                row["lambda_for_kappa54_eq_1"] for row in coupling if abs(row["R"] - 200.0) < 1.0e-12
            ),
            "interpretation": (
                "The clean kappa54=1 benchmark is not a small-number tuning: it requires "
                "|lambda|/g=sqrt(2), i.e. lambda≈0.780 at the R=200 point. It is an O(1) "
                "gauge-Yukawa locking relation. Without a symmetry or fixed-line argument it "
                "should be stated as a boundary condition; if not imposed, kappa54 must remain "
                "an explicit heavy-threshold scan parameter."
            ),
        },
    }


def write_csvs(payload: dict[str, Any]) -> None:
    gen_fields = [
        "i",
        "j",
        "block",
        "generator_norm_TrTT",
        "tangent_norm2_Tr",
        "tangent_norm",
        "mass_coefficient_MV2_over_g2v2",
        "mass_coefficient_MV_over_gv",
    ]
    with (OUT / "vector_generator_norms.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=gen_fields)
        writer.writeheader()
        writer.writerows(generator_audit_rows())

    coupling_fields = [
        "R",
        "alphaG_inv",
        "g_G",
        "lambda_for_kappa54_eq_1",
        "lambda_over_g_for_kappa54_eq_1",
        "lambda_over_4pi_sqrt",
    ]
    with (OUT / "vector_coupling_relation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=coupling_fields)
        writer.writeheader()
        writer.writerows(payload["coupling_rows"])

    kappa_fields = [
        "kappa_54",
        "lambda_over_g",
        "lambda_at_R200_g",
        "safe_points",
        "alphaG_inv",
        "M_Sigma3_GeV",
        "projected_l2_delta54",
        "perturbative_lambda_lt_sqrt4pi",
    ]
    with (OUT / "kappa54_lambda_over_g_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=kappa_fields)
        writer.writeheader()
        writer.writerows(payload["kappa_to_coupling_rows_R200"])


def write_report(payload: dict[str, Any]) -> None:
    summary = payload["generator_summary"]
    verdict = payload["verdict"]
    lines: list[str] = []
    lines.append("# Single-54_H vector mass normalization audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Convention")
    lines.append("")
    lines.append("```text")
    lines.append("K = Tr(S^dagger S)")
    lines.append("Tr_fund(T_a T_b)=delta_ab")
    lines.append("g^2 A^2 ||T<S>||^2 = (1/2) M_V^2 A^2")
    lines.append("```")
    lines.append("")
    lines.append("## Generator audit")
    lines.append("")
    lines.append(f"Broken generators: {summary['broken_generators']}; stabilizer generators: {summary['stabilizer_generators']}.")
    lines.append(f"For every broken mixed generator, `||T S0|| = {summary['broken_tangent_norm']:.6g}`.")
    lines.append(f"Thus `M_V/(g |v|) = {summary['broken_MV_over_gv']:.6g}`.")
    lines.append("")
    lines.append("## Coupling relation")
    lines.append("")
    lines.append("```text")
    lines.append("M_V = 5 sqrt(2) g |v|")
    lines.append("|v| = |m|/|lambda|")
    lines.append("M_54,phys = 5 |m|")
    lines.append("kappa_54 = M_54,phys/M_V = |lambda|/(sqrt(2) g)")
    lines.append("```")
    lines.append("")
    lines.append("For `kappa_54=1`, `|lambda|/g=sqrt(2)`.")
    lines.append("")
    lines.append("## R=200 benchmark")
    lines.append("")
    lines.append(f"`g_G = {verdict['R200_g']:.9f}` and `lambda_required = {verdict['R200_lambda_required']:.9f}`.")
    lines.append("")
    lines.append("## Threshold tolerance")
    lines.append("")
    budget = payload["tuning_budget"]["budgets"]["1.0e-02"]
    lines.append(
        "Keeping the 54 threshold below `||P Delta54||<1e-2` allows "
        f"`kappa_54` in [{budget['kappa_range'][0]:.6f}, {budget['kappa_range'][1]:.6f}], "
        f"or `lambda/g` in [{budget['lambda_over_g_range'][0]:.6f}, {budget['lambda_over_g_range'][1]:.6f}]."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(verdict["interpretation"])
    lines.append("")
    (OUT / "single_54_vector_mass_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    write_csvs(payload)
    (OUT / "single_54_vector_mass_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Single-54_H vector mass normalization audit complete")
    print(f"  generator_passes={payload['generator_summary']['passes']}")
    print(f"  broken_tangent_norm={payload['generator_summary']['broken_tangent_norm']:.6g}")
    print(f"  MV/(g|v|)={payload['generator_summary']['broken_MV_over_gv']:.6g}")
    print(f"  kappa54=lambda/(sqrt(2)g); lambda/g for kappa=1={math.sqrt(2.0):.9f}")
    print(f"  R200 g={payload['verdict']['R200_g']:.9f}, lambda_req={payload['verdict']['R200_lambda_required']:.9f}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
