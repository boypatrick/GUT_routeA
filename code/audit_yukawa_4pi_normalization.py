#!/usr/bin/env python3
"""Audit the missing 4pi factor in the two-loop Yukawa gauge terms.

No web lookup is used.  This script leaves the original benchmark files in
place and reruns the same cached-input scan with the convention

    alpha_y = Tr(Y^dagger Y)/(4 pi)

inside the two-loop gauge beta function.  It quantifies how much the
Yukawa-refined threshold branch moves when ordinary Yukawa squares are
converted to alpha_y in the alpha_i RGE.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "yukawa_4pi_audit"
ORIGINAL_SUMMARY = ROOT / "output" / "yukawa_superpotential_rge" / "yukawa_superpotential_rge_summary.json"


def beta_mssm_state_yukawa_alpha(state: np.ndarray, mu: float, ratios: dict[str, list[float]], mr: list[float]) -> np.ndarray:
    alpha = state[:3]
    yt, yb, ytau, ynu = state[3:]

    tr_u = yt * yt * sum(r * r for r in ratios["up"])
    tr_d = yb * yb * sum(r * r for r in ratios["down"])
    tr_e = ytau * ytau * sum(r * r for r in ratios["charged_lepton"])
    tr_nu = ynu * ynu * base.active_nu_trace_factor(mu, ratios["neutrino_dirac"], mr)

    ordinary_yukawa_drag = (
        base.MSSM_C_TOP * tr_u
        + base.MSSM_C_BOTTOM * tr_d
        + base.MSSM_C_TAU * tr_e
        + base.MSSM_C_NU * tr_nu
    )
    yukawa_drag_alpha = ordinary_yukawa_drag / (4.0 * math.pi)
    dalpha = (
        base.MSSM_B1 * alpha**2 / (2.0 * math.pi)
        + alpha**2 * ((base.MSSM_B2 @ alpha) - yukawa_drag_alpha) / (8.0 * math.pi**2)
    )

    g1sq, g2sq, g3sq = 4.0 * math.pi * alpha
    nu_on = 1.0 if mu >= max(mr) else 0.0
    ynu2 = ynu * ynu * nu_on

    dyt = yt * (6.0 * yt * yt + yb * yb + ynu2 - (13.0 / 15.0) * g1sq - 3.0 * g2sq - (16.0 / 3.0) * g3sq)
    dyb = yb * (yt * yt + 6.0 * yb * yb + ytau * ytau - (7.0 / 15.0) * g1sq - 3.0 * g2sq - (16.0 / 3.0) * g3sq)
    dytau = ytau * (3.0 * yb * yb + 4.0 * ytau * ytau + ynu2 - (9.0 / 5.0) * g1sq - 3.0 * g2sq)
    dynu = ynu * nu_on * (3.0 * yt * yt + 4.0 * ynu * ynu + ytau * ytau - (3.0 / 5.0) * g1sq - 3.0 * g2sq)
    dy = np.array([dyt, dyb, dytau, dynu], dtype=float) / (16.0 * math.pi**2)
    return np.concatenate([dalpha, dy])


def write_csv(rows: list[dict[str, float | bool]]) -> None:
    with (OUT / "yukawa_4pi_corrected_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def compare_best(original: dict[str, Any], corrected: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "MG_GeV",
        "alphaG_inv",
        "lambda_T",
        "lambda_S",
        "chi",
        "M_HC_GeV",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "log_HC",
        "log_Sigma3",
        "log_Sigma8",
    ]
    rows = {}
    for key in keys:
        old = float(original[key])
        new = float(corrected[key])
        rows[key] = {
            "original": old,
            "corrected": new,
            "delta": new - old,
            "ratio": new / old if old != 0 else math.nan,
        }
    return rows


def write_report(summary: dict[str, Any]) -> None:
    comparison = summary["best_comparison"]
    original = summary["original_summary"]
    corrected = summary["corrected_summary"]
    old_best = original["best"]
    new_best = corrected["best"]
    lines: list[str] = []
    lines.append("# Yukawa 4pi normalization audit")
    lines.append("")
    lines.append("No web lookup was used.  This audit reruns the Yukawa-refined two-loop")
    lines.append("scan with ordinary Yukawa squares converted to")
    lines.append("`alpha_y = Tr(Y^dagger Y)/(4 pi)` in the gauge beta function.")
    lines.append("")
    lines.append("## Convention")
    lines.append("")
    lines.append("The alpha-form two-loop term should be interpreted as")
    lines.append("")
    lines.append("```text")
    lines.append("d alpha_i/dt = b_i alpha_i^2/(2 pi)")
    lines.append("             + alpha_i^2/(8 pi^2)")
    lines.append("               [sum_j B_ij alpha_j - sum_a c_i^a Tr(Y_a^dag Y_a)/(4 pi)]")
    lines.append("```")
    lines.append("")
    lines.append("when `Y` denotes ordinary Yukawa matrices.")
    lines.append("")
    lines.append("## Best-point comparison")
    lines.append("")
    lines.append("| quantity | original | 4pi-corrected | ratio |")
    lines.append("|---|---:|---:|---:|")
    for key, row in comparison.items():
        lines.append(f"| `{key}` | {row['original']:.6e} | {row['corrected']:.6e} | {row['ratio']:.6e} |")
    lines.append("")
    lines.append("## Scan diagnostics")
    lines.append("")
    old_diag = original["scan_diagnostics"]
    new_diag = corrected["scan_diagnostics"]
    lines.append("| diagnostic | original | 4pi-corrected |")
    lines.append("|---|---:|---:|")
    for key in ["total_points", "proton_filter_perturbative_safe_points", "safe_single_scale_factor2_points"]:
        lines.append(f"| `{key}` | {old_diag[key]} | {new_diag[key]} |")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("The branch is not destroyed by the 4pi convention correction on this grid,")
    lines.append("but the best point and lifetimes move.  Therefore every downstream")
    lines.append("mediator-threshold and Goldstone-locking table derived from the old cached")
    lines.append("Yukawa scan must be regenerated before paper claims are final.")
    lines.append("")
    lines.append("Original best:")
    lines.append("")
    lines.append("```text")
    lines.append(f"MG = {old_best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {old_best['alphaG_inv']:.6f}")
    lines.append(f"M_Sigma3 = {old_best['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"tau_dim6 = {old_best['tau_dim6_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("Corrected best:")
    lines.append("")
    lines.append("```text")
    lines.append(f"MG = {new_best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {new_best['alphaG_inv']:.6f}")
    lines.append(f"M_Sigma3 = {new_best['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"tau_dim6 = {new_best['tau_dim6_years']:.6e} yr")
    lines.append("```")
    (OUT / "yukawa_4pi_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    original_beta = base.beta_mssm_state
    try:
        base.beta_mssm_state = beta_mssm_state_yukawa_alpha
        rows, corrected_summary = base.scan()
    finally:
        base.beta_mssm_state = original_beta

    original_summary = json.loads(ORIGINAL_SUMMARY.read_text(encoding="utf-8"))
    summary = {
        "note": "No web lookup used. Corrects the alpha-form Yukawa drag by using Tr(Y^dag Y)/(4pi).",
        "original_summary_path": str(ORIGINAL_SUMMARY),
        "corrected_summary": corrected_summary,
        "original_summary": original_summary,
        "best_comparison": compare_best(original_summary["best"], corrected_summary["best"]),
    }
    write_csv(rows)
    (OUT / "yukawa_4pi_audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)

    old_best = original_summary["best"]
    new_best = corrected_summary["best"]
    print("Yukawa 4pi normalization audit")
    print(f"  original MG: {old_best['MG_GeV']:.6e} GeV")
    print(f"  corrected MG: {new_best['MG_GeV']:.6e} GeV")
    print(f"  original alphaG^-1: {old_best['alphaG_inv']:.6f}")
    print(f"  corrected alphaG^-1: {new_best['alphaG_inv']:.6f}")
    print(f"  original safe points: {original_summary['scan_diagnostics']['proton_filter_perturbative_safe_points']}")
    print(f"  corrected safe points: {corrected_summary['scan_diagnostics']['proton_filter_perturbative_safe_points']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
