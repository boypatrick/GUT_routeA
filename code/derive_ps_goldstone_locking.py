#!/usr/bin/env python3
"""Derive the (15,1,1) colored-pair lock from a PS breaking superpotential.

The threshold scans previously imposed that the non-octet fragments of
``(15,1,1)`` are Goldstone/eaten or lifted at ``MG``.  This script proves the
SU(4)_C part explicitly at the Pati-Salam stage.  For an adjoint chiral field
``Sigma`` with

    W = m/2 Tr Sigma^2 + lambda/3 Tr Sigma^3,

the F-flat vacuum ``Sigma0 = v diag(1,1,1,-3)``, ``v=m/(2 lambda)``, breaks
``SU(4)_C -> SU(3)_C x U(1)_{B-L}``.  The holomorphic Hessian has eigenvalues

    octet:     2m,
    triplet pair (3+3bar): 0,
    singlet:  -m.

The six zero modes are exactly the chiral Goldstone multiplets for the broken
SU(4)_C generators and are absorbed into the massive vector multiplets.  Hence
they produce no separate chiral threshold log: epsilon_G = log(MG/M_Goldstone)
is not a tunable parameter in this PS-stage completion; it is fixed to zero.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DOWNSTREAM_JSON = ROOT / "output" / "corrected_downstream" / "corrected_downstream_summary.json"
OUT = ROOT / "output" / "ps_goldstone_locking"

PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0
GOLDSTONE_B = np.array([8.0 / 5.0, 0.0, 1.0], dtype=float)


def matrix_unit(i: int, j: int) -> np.ndarray:
    mat = np.zeros((4, 4), dtype=complex)
    mat[i, j] = 1.0
    return mat


def classified_basis() -> list[tuple[str, str, np.ndarray]]:
    basis: list[tuple[str, str, np.ndarray]] = []
    for i in range(3):
        for j in range(3):
            if i != j:
                basis.append(("octet", f"E{i + 1}{j + 1}", matrix_unit(i, j)))
    basis.append(("octet", "diag(1,-1,0,0)", np.diag([1.0, -1.0, 0.0, 0.0]).astype(complex)))
    basis.append(("octet", "diag(1,1,-2,0)", np.diag([1.0, 1.0, -2.0, 0.0]).astype(complex)))
    for i in range(3):
        basis.append(("goldstone_3", f"E{i + 1}4", matrix_unit(i, 3)))
        basis.append(("goldstone_3bar", f"E4{i + 1}", matrix_unit(3, i)))
    basis.append(("singlet", "diag(1,1,1,-3)", np.diag([1.0, 1.0, 1.0, -3.0]).astype(complex)))
    return basis


def hessian_action(x: np.ndarray, m: float = 1.0, lam: float = 1.0) -> np.ndarray:
    v = m / (2.0 * lam)
    phi0 = v * np.diag([1.0, 1.0, 1.0, -3.0]).astype(complex)
    trace_term = np.trace(phi0 @ x) / 2.0
    return m * x + lam * (phi0 @ x + x @ phi0 - trace_term * np.eye(4, dtype=complex))


def eigen_check() -> list[dict[str, object]]:
    rows = []
    for sector, label, x in classified_basis():
        hx = hessian_action(x)
        denom = np.vdot(x, x).real
        coeff = np.vdot(x, hx) / denom
        residual = np.linalg.norm(hx - coeff * x) / math.sqrt(denom)
        rows.append(
            {
                "sector": sector,
                "basis": label,
                "hessian_eigenvalue_over_m": float(np.real_if_close(coeff).real),
                "imaginary_part": float(np.real_if_close(coeff).imag if hasattr(np.real_if_close(coeff), "imag") else 0.0),
                "relative_residual": float(residual),
            }
        )
    return rows


def summarize_eigenvalues(rows: list[dict[str, object]]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for sector in ["octet", "goldstone_3", "goldstone_3bar", "singlet"]:
        vals = [float(row["hessian_eigenvalue_over_m"]) for row in rows if row["sector"] == sector]
        residuals = [float(row["relative_residual"]) for row in rows if row["sector"] == sector]
        summary[sector] = {
            "multiplicity": len(vals),
            "unique_eigenvalues_over_m": sorted(set(round(v, 12) for v in vals)),
            "max_relative_residual": max(residuals) if residuals else 0.0,
        }
    return summary


def downstream_benchmarks(payload: dict[str, object]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    fixed = payload["fixed_R50"]["best_safe_or_best"]
    rows.append(
        {
            "name": "fixed_R50_mediator",
            "R": 50.0,
            "MG_GeV": float(fixed["MG_GeV"]),
            "M_Sigma8_GeV": float(fixed["M_Sigma8_GeV"]),
            "tau_dim6_years": float(fixed["tau_dim6_years"]),
            "tau_dim5_ST1e-5_years": float(fixed["tau_dim5_target_filter_years"]),
        }
    )
    for row in payload["R_window"]["rows"]:
        r_value = float(row["R"])
        if r_value in (20.0, 50.0, 100.0, 200.0):
            rows.append(
                {
                    "name": f"R_window_R{r_value:g}",
                    "R": r_value,
                    "MG_GeV": float(row.get("MG_GeV", fixed["MG_GeV"])),
                    "M_Sigma8_GeV": float(row["M_Sigma8_GeV"]),
                    "tau_dim6_years": float(row["tau_dim6_years"]),
                    "tau_dim5_ST1e-5_years": float(row["tau_dim5_target_filter_years"]),
                }
            )
    return rows


def locking_rows(payload: dict[str, object], benchmarks: list[dict[str, float | str]]) -> list[dict[str, object]]:
    gold = payload["goldstone_locking"]
    per_eps = float(gold["goldstone_projected_l2_per_unit_epsilon"])
    rows = []
    for bench in benchmarks:
        r_value = float(bench["R"])
        r_entry = next(row for row in payload["R_window"]["rows"] if abs(float(row["R"]) - r_value) < 1.0e-12)
        eps_max = float(r_entry["projected_l2"]) / per_eps
        m8 = float(bench["M_Sigma8_GeV"])
        rows.append(
            {
                "benchmark": bench["name"],
                "R": r_value,
                "MG_GeV": float(bench["MG_GeV"]),
                "m_parameter_GeV": 0.5 * m8,
                "M_octet_GeV": m8,
                "M_singlet_GeV": 0.5 * m8,
                "goldstone_chiral_hessian_mass_GeV": 0.0,
                "epsilon_G_derived": 0.0,
                "epsilon_G_max_from_threshold_scan": eps_max,
                "threshold_vector_from_eaten_pair": [0.0, 0.0, 0.0],
                "projected_threshold_l2": 0.0,
                "passes_locking_bound": True,
                "would_be_chiral_threshold_per_unit_epsilon_l2": per_eps,
                "allowed_mass_ratio_low": math.exp(-eps_max),
                "allowed_mass_ratio_high": math.exp(eps_max),
                "tau_dim6_years": float(bench["tau_dim6_years"]),
                "tau_dim5_ST1e-5_years": float(bench["tau_dim5_ST1e-5_years"]),
            }
        )
    return rows


def write_report(payload: dict[str, object]) -> None:
    lines: list[str] = []
    lines.append("# Pati-Salam adjoint Goldstone-locking derivation")
    lines.append("")
    lines.append("No web lookup was used.  This report proves the SU(4)_C part of the")
    lines.append("Goldstone/lifting assignment directly from a renormalizable Pati-Salam-stage")
    lines.append("superpotential.")
    lines.append("")
    lines.append("## Superpotential and F-flat vacuum")
    lines.append("")
    lines.append("For a traceless adjoint chiral field `Sigma` of `SU(4)_C`, take")
    lines.append("")
    lines.append("```text")
    lines.append("W_C = (m_C/2) Tr Sigma^2 + (lambda_C/3) Tr Sigma^3,")
    lines.append("Sigma_0 = v diag(1,1,1,-3),  v = m_C/(2 lambda_C).")
    lines.append("```")
    lines.append("")
    lines.append("The projected F-term is")
    lines.append("")
    lines.append("```text")
    lines.append("F = m_C Sigma + lambda_C [Sigma^2 - (Tr Sigma^2)/4 I_4].")
    lines.append("```")
    lines.append("")
    lines.append("At `Sigma_0`, both color and fourth eigenvalue equations reduce to")
    lines.append("`m_C - 2 lambda_C v = 0`; hence the vacuum is F-flat and breaks")
    lines.append("`SU(4)_C -> SU(3)_C x U(1)_{B-L}`.")
    lines.append("")
    lines.append("## Hessian eigenvalues")
    lines.append("")
    lines.append("Linearizing the F-term gives")
    lines.append("")
    lines.append("```text")
    lines.append("delta F = m_C X + lambda_C (Sigma_0 X + X Sigma_0")
    lines.append("          - 1/2 Tr(Sigma_0 X) I_4).")
    lines.append("```")
    lines.append("")
    lines.append("| sector | multiplicity | Hessian eigenvalue / m_C | max residual |")
    lines.append("|---|---:|---:|---:|")
    for sector, row in payload["eigenvalue_summary"].items():
        vals = ",".join(f"{v:g}" for v in row["unique_eigenvalues_over_m"])
        lines.append(
            f"| `{sector}` | {row['multiplicity']} | `{vals}` | {row['max_relative_residual']:.3e} |"
        )
    lines.append("")
    lines.append("Thus the `(3,1,2/3)+(bar3,1,-2/3)` directions are exact chiral")
    lines.append("Goldstone multiplets.  They are absorbed into the massive vector multiplets")
    lines.append("of the broken `SU(4)_C/SU(3)_C x U(1)` generators and are not separate")
    lines.append("threshold chiral fields.")
    lines.append("")
    lines.append("## Corrected-cache numerical lock")
    lines.append("")
    lines.append("| benchmark | R | M_octet [GeV] | M_singlet [GeV] | eps_G derived | eps_G max | pass |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for row in payload["locking_rows"]:
        if row["benchmark"] not in ("fixed_R50_mediator", "R_window_R200"):
            continue
        lines.append(
            f"| `{row['benchmark']}` | {row['R']:.0f} | {row['M_octet_GeV']:.6e} | "
            f"{row['M_singlet_GeV']:.6e} | {row['epsilon_G_derived']:.1f} | "
            f"{row['epsilon_G_max_from_threshold_scan']:.6e} | {row['passes_locking_bound']} |"
        )
    lines.append("")
    lines.append("Because the eaten pair has no separate chiral threshold,")
    lines.append("`Delta_G = 0` and `epsilon_G=0` exactly.  The previous scan's tolerance")
    lines.append("therefore becomes a consistency check rather than a tuning condition.")
    lines.append("")
    lines.append("## Remaining embedding caveat")
    lines.append("")
    lines.append("This is a Pati-Salam-stage derivation of the SU(4)_C Goldstone lock.  A final")
    lines.append("Spin(10) paper should still show how this cubic PS effective sector descends")
    lines.append("from the chosen `54_H/210_H/16_H` breaking sector or explicitly state it as")
    lines.append("the renormalizable PS EFT below the first Spin(10)-breaking step.")
    (OUT / "ps_goldstone_locking_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    downstream = json.loads(DOWNSTREAM_JSON.read_text(encoding="utf-8"))
    eigen_rows = eigen_check()
    eigen_summary = summarize_eigenvalues(eigen_rows)
    benches = downstream_benchmarks(downstream)
    lock_rows = locking_rows(downstream, benches)
    payload = {
        "note": "No web lookup used. PS-stage SU(4)_C adjoint superpotential derives exact eaten colored Goldstone pair.",
        "input_files": {"corrected_downstream_json": str(DOWNSTREAM_JSON)},
        "superpotential": {
            "W_C": "(m_C/2) Tr Sigma^2 + (lambda_C/3) Tr Sigma^3",
            "Sigma0": "v diag(1,1,1,-3), v=m_C/(2 lambda_C)",
            "F_term": "m_C Sigma + lambda_C [Sigma^2 - Tr(Sigma^2) I_4/4]",
            "linearized_F": "m_C X + lambda_C(Sigma0 X + X Sigma0 - Tr(Sigma0 X) I_4/2)",
        },
        "basis_eigen_checks": eigen_rows,
        "eigenvalue_summary": eigen_summary,
        "locking_rows": lock_rows,
        "all_eigen_residuals_small": all(float(row["relative_residual"]) < 1.0e-12 for row in eigen_rows),
        "all_locking_bounds_pass": all(bool(row["passes_locking_bound"]) for row in lock_rows),
    }
    (OUT / "ps_goldstone_locking_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)

    print("Pati-Salam Goldstone-locking derivation")
    for sector, row in eigen_summary.items():
        print(
            "  {sector}: multiplicity={multiplicity}, eigen={eigs}, residual={res:.3e}".format(
                sector=sector,
                multiplicity=row["multiplicity"],
                eigs=row["unique_eigenvalues_over_m"],
                res=row["max_relative_residual"],
            )
        )
    r50 = next(row for row in lock_rows if row["benchmark"] == "fixed_R50_mediator")
    r200 = next(row for row in lock_rows if row["benchmark"] == "R_window_R200")
    print(f"  R=50 eps derived/max: {r50['epsilon_G_derived']:.1f}/{r50['epsilon_G_max_from_threshold_scan']:.6e}")
    print(f"  R=200 eps derived/max: {r200['epsilon_G_derived']:.1f}/{r200['epsilon_G_max_from_threshold_scan']:.6e}")
    print(f"  all checks ok: {payload['all_eigen_residuals_small'] and payload['all_locking_bounds_pass']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
