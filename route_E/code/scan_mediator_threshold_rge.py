#!/usr/bin/env python3
"""Two-loop Yukawa/RGE scan with finite mediator thresholds included.

This script takes the renormalizable 45-mediator UV completion seriously at the
matching scale.  The two heavy mediator eigenstates in each 45 block produce a
finite threshold vector Delta_med.  The matching equation is therefore

    alpha_i^{-1}(MG) = alpha_G^{-1} + Delta_med_i
                     + sum_r b_i^r log(MG/M_r)/(2 pi).

The remaining fit solves for the H_C, Sigma_3, and Sigma_8 logarithms exactly
in the traceless plane, then maps them back to the breaking-chain
superpotential.  The script also records the full (15,1,1) Goldstone/lifting
assignment: only the SM octet remains at M_Sigma8; the colored triplet pair is
eaten/lifted at MG, and the singlet carries no one-loop SM threshold.
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


UV_JSON = ROOT / "output" / "uv_projector_mediator" / "uv_projector_mediator_summary.json"
BASE_JSON = ROOT / "output" / "yukawa_superpotential_rge" / "yukawa_superpotential_rge_summary.json"
BASE_SCAN_CSV = ROOT / "output" / "yukawa_superpotential_rge" / "yukawa_superpotential_rge_scan.csv"
OUT = ROOT / "output" / "mediator_threshold_rge"

SIGMA8_OCTET_B = np.array([0.0, 0.0, 3.0], dtype=float)
SIGMA8_COLORED_PAIR_B = np.array([8.0 / 5.0, 0.0, 1.0], dtype=float)
SIGMA8_SINGLET_B = np.array([0.0, 0.0, 0.0], dtype=float)


def load_uv_threshold() -> dict[str, Any]:
    payload = json.loads(UV_JSON.read_text(encoding="utf-8"))
    heavy = payload["benchmark"]["heavy_threshold"]
    return {
        "source": str(UV_JSON),
        "R": float(payload["benchmark"]["r_med"]),
        "eta": float(payload["benchmark"]["parameters"]["eta"]),
        "delta_full": np.array(heavy["delta_full"], dtype=float),
        "projected_delta": np.array(heavy["projected_delta"], dtype=float),
        "projected_l2": float(heavy["projected_l2"]),
        "raw_payload": payload,
    }


def exact_threshold_logs_with_mediator(alpha_inv: np.ndarray, delta_med: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    mat = base.PROJECTOR @ base.HEAVY_BASIS / (2.0 * math.pi)
    target = base.PROJECTOR @ (alpha_inv - delta_med)
    logs, *_ = np.linalg.lstsq(mat, target, rcond=None)
    threshold = base.HEAVY_BASIS @ logs / (2.0 * math.pi)
    residual_vec = base.PROJECTOR @ (alpha_inv - delta_med - threshold)
    residual = float(np.linalg.norm(residual_vec))
    alpha_g_inv = float(np.mean(alpha_inv - delta_med - threshold))
    return logs, threshold, residual, alpha_g_inv


def analyze_point_with_mediator(
    mg: float,
    msusy: float,
    tan_beta: float,
    state: np.ndarray,
    prefactor: float,
    delta_med: np.ndarray,
) -> dict[str, float | bool]:
    alpha_inv = 1.0 / state[:3]
    logs, threshold, residual, alpha_g_inv = exact_threshold_logs_with_mediator(alpha_inv, delta_med)
    params = base.superpotential_from_logs(logs)
    mx_min = base.mx_required(alpha_g_inv, prefactor)
    yt_mg, yb_mg, ytau_mg, ynu_mg = state[3:]
    m_t = params["lambda_T"] * mg
    loop = (1.0 / 25.0) / (4.0 * math.pi)
    c6_limit = base.c6_max_for_target(prefactor)
    denominator = (yt_mg * yb_mg / m_t) * loop * base.MWINO_GEV / (msusy * msusy)
    filter_required = c6_limit / max(denominator, 1.0e-300)
    max_abs_log = float(np.max(np.abs(logs)))
    proton_penalty = max(0.0, math.log(mx_min / mg))
    filter_penalty = max(0.0, math.log(base.TRIPLET_FILTER_TARGET / filter_required))
    perturbative_penalty = 0.0 if params["perturbative"] else 10.0
    mediator_projected_l2 = float(np.linalg.norm(base.PROJECTOR @ delta_med))
    score = (
        max_abs_log
        + 0.10 * abs(params["chi"])
        + 0.10 * max(0.0, -math.log(max(params["cancellation_index"], 1.0e-12)))
        + 20.0 * proton_penalty * proton_penalty
        + 20.0 * filter_penalty * filter_penalty
        + perturbative_penalty
        + 1000.0 * residual
        + 0.25 * mediator_projected_l2
    )

    c6_gauge = (4.0 * math.pi / alpha_g_inv) / (mg * mg)
    c6_d5 = base.TRIPLET_FILTER_TARGET * yt_mg * yb_mg / m_t * loop * base.MWINO_GEV / (msusy * msusy)

    return {
        "score": float(score),
        "tan_beta": float(tan_beta),
        "MSUSY_GeV": float(msusy),
        "MG_GeV": float(mg),
        "alpha1_inv_MG_raw": float(alpha_inv[0]),
        "alpha2_inv_MG_raw": float(alpha_inv[1]),
        "alpha3_inv_MG_raw": float(alpha_inv[2]),
        "mediator_delta1": float(delta_med[0]),
        "mediator_delta2": float(delta_med[1]),
        "mediator_delta3": float(delta_med[2]),
        "mediator_projected_l2": mediator_projected_l2,
        "alphaG_inv": alpha_g_inv,
        "yt_MG": float(yt_mg),
        "yb_MG": float(yb_mg),
        "ytau_MG": float(ytau_mg),
        "ynu_MG": float(ynu_mg),
        "lambda_T": float(params["lambda_T"]),
        "lambda_S": float(params["lambda_S"]),
        "chi": float(params["chi"]),
        "kappa_3": float(params["kappa_3"]),
        "kappa_8": float(params["kappa_8"]),
        "cancellation_index": float(params["cancellation_index"]),
        "perturbative_superpotential": bool(params["perturbative"]),
        "single_scale_factor2": bool(params["single_scale_factor2"]),
        "log_HC": float(logs[0]),
        "log_Sigma3": float(logs[1]),
        "log_Sigma8": float(logs[2]),
        "M_HC_GeV": float(params["lambda_T"] * mg),
        "M_Sigma3_GeV": float(params["kappa_3"] * mg),
        "M_Sigma8_GeV": float(params["kappa_8"] * mg),
        "M_X_GeV": float(mg),
        "M_X_min_GeV": float(mx_min),
        "rho_X": float(max(1.0, mx_min / mg)),
        "residual_l2_after_mediator": float(residual),
        "max_abs_log": float(max_abs_log),
        "triplet_filter_required": float(filter_required),
        "triplet_filter_target": base.TRIPLET_FILTER_TARGET,
        "triplet_filter_safe": bool(base.TRIPLET_FILTER_TARGET <= filter_required),
        "tau_dim6_years": float(base.lifetime_from_c6(c6_gauge, prefactor)),
        "tau_dim5_target_filter_years": float(base.lifetime_from_c6(c6_d5, prefactor)),
    }


def analyze_cached_row_with_mediator(
    source: dict[str, str],
    prefactor: float,
    delta_med: np.ndarray,
) -> dict[str, float | bool]:
    mg = float(source["MG_GeV"])
    msusy = float(source["MSUSY_GeV"])
    tan_beta = float(source["tan_beta"])
    alpha_inv = np.array(
        [
            float(source["alpha1_inv_MG"]),
            float(source["alpha2_inv_MG"]),
            float(source["alpha3_inv_MG"]),
        ],
        dtype=float,
    )
    yt_mg = float(source["yt_MG"])
    yb_mg = float(source["yb_MG"])
    ytau_mg = float(source["ytau_MG"])
    ynu_mg = float(source["ynu_MG"])

    logs, threshold, residual, alpha_g_inv = exact_threshold_logs_with_mediator(alpha_inv, delta_med)
    params = base.superpotential_from_logs(logs)
    mx_min = base.mx_required(alpha_g_inv, prefactor)
    m_t = params["lambda_T"] * mg
    loop = (1.0 / 25.0) / (4.0 * math.pi)
    c6_limit = base.c6_max_for_target(prefactor)
    denominator = (yt_mg * yb_mg / m_t) * loop * base.MWINO_GEV / (msusy * msusy)
    filter_required = c6_limit / max(denominator, 1.0e-300)
    max_abs_log = float(np.max(np.abs(logs)))
    proton_penalty = max(0.0, math.log(mx_min / mg))
    filter_penalty = max(0.0, math.log(base.TRIPLET_FILTER_TARGET / filter_required))
    perturbative_penalty = 0.0 if params["perturbative"] else 10.0
    mediator_projected_l2 = float(np.linalg.norm(base.PROJECTOR @ delta_med))
    score = (
        max_abs_log
        + 0.10 * abs(params["chi"])
        + 0.10 * max(0.0, -math.log(max(params["cancellation_index"], 1.0e-12)))
        + 20.0 * proton_penalty * proton_penalty
        + 20.0 * filter_penalty * filter_penalty
        + perturbative_penalty
        + 1000.0 * residual
        + 0.25 * mediator_projected_l2
    )

    c6_gauge = (4.0 * math.pi / alpha_g_inv) / (mg * mg)
    c6_d5 = base.TRIPLET_FILTER_TARGET * yt_mg * yb_mg / m_t * loop * base.MWINO_GEV / (msusy * msusy)

    return {
        "score": float(score),
        "tan_beta": tan_beta,
        "MSUSY_GeV": msusy,
        "MG_GeV": mg,
        "alpha1_inv_MG_raw": float(alpha_inv[0]),
        "alpha2_inv_MG_raw": float(alpha_inv[1]),
        "alpha3_inv_MG_raw": float(alpha_inv[2]),
        "mediator_delta1": float(delta_med[0]),
        "mediator_delta2": float(delta_med[1]),
        "mediator_delta3": float(delta_med[2]),
        "mediator_projected_l2": mediator_projected_l2,
        "alphaG_inv": alpha_g_inv,
        "yt_MG": yt_mg,
        "yb_MG": yb_mg,
        "ytau_MG": ytau_mg,
        "ynu_MG": ynu_mg,
        "lambda_T": float(params["lambda_T"]),
        "lambda_S": float(params["lambda_S"]),
        "chi": float(params["chi"]),
        "kappa_3": float(params["kappa_3"]),
        "kappa_8": float(params["kappa_8"]),
        "cancellation_index": float(params["cancellation_index"]),
        "perturbative_superpotential": bool(params["perturbative"]),
        "single_scale_factor2": bool(params["single_scale_factor2"]),
        "log_HC": float(logs[0]),
        "log_Sigma3": float(logs[1]),
        "log_Sigma8": float(logs[2]),
        "M_HC_GeV": float(params["lambda_T"] * mg),
        "M_Sigma3_GeV": float(params["kappa_3"] * mg),
        "M_Sigma8_GeV": float(params["kappa_8"] * mg),
        "M_X_GeV": float(mg),
        "M_X_min_GeV": float(mx_min),
        "rho_X": float(max(1.0, mx_min / mg)),
        "residual_l2_after_mediator": float(residual),
        "max_abs_log": float(max_abs_log),
        "triplet_filter_required": float(filter_required),
        "triplet_filter_target": base.TRIPLET_FILTER_TARGET,
        "triplet_filter_safe": bool(base.TRIPLET_FILTER_TARGET <= filter_required),
        "tau_dim6_years": float(base.lifetime_from_c6(c6_gauge, prefactor)),
        "tau_dim5_target_filter_years": float(base.lifetime_from_c6(c6_d5, prefactor)),
    }


def goldstone_lifting_assignment(log_sigma8: float) -> dict[str, Any]:
    octet = SIGMA8_OCTET_B * log_sigma8 / (2.0 * math.pi)
    colored_if_light = SIGMA8_COLORED_PAIR_B * log_sigma8 / (2.0 * math.pi)
    singlet = SIGMA8_SINGLET_B * log_sigma8 / (2.0 * math.pi)
    extra_residual = base.PROJECTOR @ colored_if_light
    return {
        "decomposition": [
            {
                "component": "(8,1,0)",
                "role": "physical light octet Sigma8",
                "b": SIGMA8_OCTET_B.tolist(),
                "log_MG_over_M": log_sigma8,
                "threshold": octet.tolist(),
            },
            {
                "component": "(3,1,2/3)+(bar3,1,-2/3)",
                "role": "Goldstone pair for SU(4)_C -> SU(3)_C x U(1)_{B-L}, eaten or lifted at MG",
                "b": SIGMA8_COLORED_PAIR_B.tolist(),
                "log_MG_over_M": 0.0,
                "threshold": [0.0, 0.0, 0.0],
                "if_left_at_MSigma8_threshold": colored_if_light.tolist(),
            },
            {
                "component": "(1,1,0)",
                "role": "singlet/radial mode, lifted at MG; no one-loop SM gauge threshold",
                "b": SIGMA8_SINGLET_B.tolist(),
                "log_MG_over_M": 0.0,
                "threshold": singlet.tolist(),
            },
        ],
        "assigned_total_threshold": octet.tolist(),
        "colored_pair_unlifted_projected_vector": extra_residual.tolist(),
        "colored_pair_unlifted_projected_l2": float(np.linalg.norm(extra_residual)),
    }


def scan() -> tuple[list[dict[str, float | bool]], dict[str, Any]]:
    prefactor = base.proton_prefactor()
    uv = load_uv_threshold()
    delta_med = uv["delta_full"]

    rows: list[dict[str, float | bool]] = []
    best: dict[str, float | bool] | None = None
    with BASE_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for source_row in reader:
            row = analyze_cached_row_with_mediator(source_row, prefactor, delta_med)
            rows.append(row)
            if best is None or float(row["score"]) < float(best["score"]):
                best = row
    assert best is not None

    safe_rows = [
        row
        for row in rows
        if bool(row["triplet_filter_safe"])
        and float(row["rho_X"]) <= 1.000001
        and bool(row["perturbative_superpotential"])
    ]
    factor2_safe_rows = [row for row in safe_rows if bool(row["single_scale_factor2"])]
    best_factor2 = min(factor2_safe_rows, key=lambda row: float(row["score"]), default=None)
    base_summary = json.loads(BASE_JSON.read_text(encoding="utf-8"))
    best_logs = np.array([best["log_HC"], best["log_Sigma3"], best["log_Sigma8"]], dtype=float)
    heavy_threshold = base.HEAVY_BASIS @ best_logs / (2.0 * math.pi)
    summary = {
        "note": "No web lookup used.  Includes the finite 45-mediator threshold vector and a complete (15,1,1) Goldstone/lifting assignment.",
        "mediator_threshold": {
            "source": uv["source"],
            "R": uv["R"],
            "eta": uv["eta"],
            "delta_full": uv["delta_full"].tolist(),
            "projected_delta": uv["projected_delta"].tolist(),
            "projected_l2": uv["projected_l2"],
            "convention": "Delta_i=sum b_i log(MG/Mheavy)/(2pi), included before solving the H_C/Sigma3/Sigma8 logs",
        },
        "base_best_for_comparison": base_summary["best"],
        "scan_diagnostics": {
            "total_points": len(rows),
            "proton_filter_perturbative_safe_points": len(safe_rows),
            "safe_single_scale_factor2_points": len(factor2_safe_rows),
            "best_factor2_safe_point": best_factor2,
        },
        "best": best,
        "best_threshold_vector_without_mediator": heavy_threshold.tolist(),
        "best_total_threshold_vector": (heavy_threshold + uv["delta_full"]).tolist(),
        "goldstone_lifting_assignment": goldstone_lifting_assignment(float(best["log_Sigma8"])),
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, float | bool]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    best = summary["best"]
    base_best = summary["base_best_for_comparison"]
    med = summary["mediator_threshold"]
    gold = summary["goldstone_lifting_assignment"]
    lines: list[str] = []
    lines.append("# Two-loop RGE scan with finite mediator threshold")
    lines.append("")
    lines.append("No web lookup was used.  This scan adds the finite 45-mediator threshold")
    lines.append("from the renormalizable projector sector before solving the heavy spectrum.")
    lines.append("")
    lines.append("## Matching equation")
    lines.append("")
    lines.append("```text")
    lines.append("alpha_i^-1(MG) = alpha_G^-1 + Delta_med_i")
    lines.append("               + sum_r b_i^r log(MG/M_r)/(2 pi)")
    lines.append("r = H_C, Sigma_3, Sigma_8(octet)")
    lines.append("```")
    lines.append("")
    lines.append("Finite mediator threshold:")
    lines.append("")
    lines.append("```text")
    lines.append(f"R = {med['R']:.6g}, eta = {med['eta']:.9f}")
    lines.append(f"Delta_med = ({med['delta_full'][0]:+.9f}, {med['delta_full'][1]:+.9f}, {med['delta_full'][2]:+.9f})")
    lines.append(f"P Delta_med = ({med['projected_delta'][0]:+.9e}, {med['projected_delta'][1]:+.9e}, {med['projected_delta'][2]:+.9e})")
    lines.append(f"||P Delta_med||_2 = {med['projected_l2']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("## Best benchmark after mediator threshold")
    lines.append("")
    lines.append("```text")
    lines.append(f"tan beta = {best['tan_beta']:.2f}")
    lines.append(f"MSUSY = {best['MSUSY_GeV']:.6e} GeV")
    lines.append(f"MG = {best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {best['alphaG_inv']:.6f}")
    lines.append(f"lambda_T = {best['lambda_T']:.9f}")
    lines.append(f"lambda_S = {best['lambda_S']:.9f}")
    lines.append(f"chi = {best['chi']:.9f}")
    lines.append(f"logs = ({best['log_HC']:+.9f}, {best['log_Sigma3']:+.9f}, {best['log_Sigma8']:+.9f})")
    lines.append(f"residual_l2_after_mediator = {best['residual_l2_after_mediator']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("| quantity | before mediator | after mediator |")
    lines.append("|---|---:|---:|")
    for key in ["MG_GeV", "alphaG_inv", "M_HC_GeV", "M_Sigma3_GeV", "M_Sigma8_GeV", "tau_dim6_years", "tau_dim5_target_filter_years"]:
        before = base_best[key]
        after = best[key]
        lines.append(f"| `{key}` | {before:.6e} | {after:.6e} |")
    lines.append("")
    lines.append("The non-universal part of the mediator threshold shifts the required")
    lines.append("Clebsch-spectrum logs only mildly; the large common part is absorbed into")
    lines.append("the fitted unified coupling.")
    lines.append("")
    lines.append("## Complete (15,1,1) assignment")
    lines.append("")
    lines.append("Under `SU(4)_C -> SU(3)_C x U(1)_{B-L}`,")
    lines.append("")
    lines.append("```text")
    lines.append("(15,1,1) -> (8,1,0) + (3,1,2/3) + (bar3,1,-2/3) + (1,1,0).")
    lines.append("```")
    lines.append("")
    lines.append("| component | role | b-vector | log(MG/M) |")
    lines.append("|---|---|---:|---:|")
    for comp in gold["decomposition"]:
        b = comp["b"]
        lines.append(
            f"| `{comp['component']}` | {comp['role']} | ({b[0]:.6g}, {b[1]:.6g}, {b[2]:.6g}) | {comp['log_MG_over_M']:.9f} |"
        )
    lines.append("")
    lines.append("Only the octet contributes to the fitted `Sigma8` threshold.  If the colored")
    lines.append("triplet pair were accidentally left at `M_Sigma8`, it would add a projected")
    lines.append(f"threshold of norm `{gold['colored_pair_unlifted_projected_l2']:.6e}`, so the")
    lines.append("Goldstone/lifting assignment is a genuine consistency condition.")
    lines.append("")
    lines.append("## Scan verdict")
    lines.append("")
    diag = summary["scan_diagnostics"]
    lines.append(f"- Total scanned RGE points: `{diag['total_points']}`.")
    lines.append(f"- Proton-safe, dimension-5-safe, perturbative points: `{diag['proton_filter_perturbative_safe_points']}`.")
    lines.append(f"- Safe single-scale factor-2 points: `{diag['safe_single_scale_factor2_points']}`.")
    lines.append("")
    lines.append("The viable branch survives the finite mediator threshold.  It remains an")
    lines.append("intermediate-`Sigma3` branch, not a single-scale branch.")
    (OUT / "mediator_threshold_rge_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = scan()
    write_csv(OUT / "mediator_threshold_rge_scan.csv", rows)
    (OUT / "mediator_threshold_rge_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    best = summary["best"]
    print("Two-loop RGE scan with finite mediator threshold")
    print(f"  tan beta: {best['tan_beta']:.2f}")
    print(f"  MSUSY: {best['MSUSY_GeV']:.6e} GeV")
    print(f"  MG: {best['MG_GeV']:.6e} GeV")
    print(f"  alphaG^-1: {best['alphaG_inv']:.6f}")
    print(f"  lambda_T, lambda_S, chi: {best['lambda_T']:.6f}, {best['lambda_S']:.6f}, {best['chi']:.6f}")
    print(f"  residual after mediator: {best['residual_l2_after_mediator']:.3e}")
    print(f"  tau d=6: {best['tau_dim6_years']:.6e} yr")
    print(f"  tau d=5 at S_T=1e-5: {best['tau_dim5_target_filter_years']:.6e} yr")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
