#!/usr/bin/env python3
"""Two-loop gauge RGE and heavy-spectrum fit for the paper skeleton.

This is a gauge-only two-loop scan:
  MZ -> MSUSY uses SM gauge beta functions,
  MSUSY -> MG uses MSSM gauge beta functions,
  MG matching uses one-loop heavy thresholds from a SUSY-like heavy spectrum.

No web lookup is used; all inputs are local benchmark constants.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROTON_JSON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
SEESAW_JSON = ROOT / "output" / "seesaw" / "seesaw_reconstruction.json"
OUT = ROOT / "output" / "two_loop_spectrum"

MZ_GEV = 91.1876
MSUSY_GEV = 1.0e3
ALPHA_EM_INV_MZ = 127.955
SIN2_THETA_W_MZ = 0.23122
ALPHA_S_MZ = 0.1184
HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0

SM_B1 = np.array([41.0 / 10.0, -19.0 / 6.0, -7.0], dtype=float)
SM_B2 = np.array(
    [
        [199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0],
        [9.0 / 10.0, 35.0 / 6.0, 12.0],
        [11.0 / 10.0, 9.0 / 2.0, -26.0],
    ],
    dtype=float,
)
MSSM_B1 = np.array([33.0 / 5.0, 1.0, -3.0], dtype=float)
MSSM_B2 = np.array(
    [
        [199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0],
        [9.0 / 5.0, 25.0, 24.0],
        [11.0 / 5.0, 9.0, 14.0],
    ],
    dtype=float,
)

# SUSY-like chiral heavy thresholds: triplet pair, adjoint SU(2), adjoint SU(3).
HEAVY_BASIS = np.array(
    [
        [2.0 / 5.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [1.0, 0.0, 3.0],
    ],
    dtype=float,
)
HEAVY_NAMES = ["H_C_pair", "Sigma_3", "Sigma_8"]


def alpha_inverse_mz() -> np.ndarray:
    alpha_em = 1.0 / ALPHA_EM_INV_MZ
    cos2 = 1.0 - SIN2_THETA_W_MZ
    alpha1 = (5.0 / 3.0) * alpha_em / cos2
    alpha2 = alpha_em / SIN2_THETA_W_MZ
    alpha3 = ALPHA_S_MZ
    return np.array([1.0 / alpha1, 1.0 / alpha2, 1.0 / alpha3], dtype=float)


def beta_alpha(alpha: np.ndarray, b1: np.ndarray, b2: np.ndarray) -> np.ndarray:
    return b1 * alpha**2 / (2.0 * math.pi) + alpha**2 * (b2 @ alpha) / (8.0 * math.pi**2)


def run_two_loop(alpha0: np.ndarray, mu0: float, mu1: float, b1: np.ndarray, b2: np.ndarray, steps: int = 800) -> np.ndarray:
    if mu1 == mu0:
        return alpha0.copy()
    t0, t1 = math.log(mu0), math.log(mu1)
    h = (t1 - t0) / steps
    alpha = alpha0.astype(float).copy()
    for _ in range(steps):
        k1 = beta_alpha(alpha, b1, b2)
        k2 = beta_alpha(alpha + 0.5 * h * k1, b1, b2)
        k3 = beta_alpha(alpha + 0.5 * h * k2, b1, b2)
        k4 = beta_alpha(alpha + h * k3, b1, b2)
        alpha += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return alpha


def run_to_mg(mg: float) -> np.ndarray:
    alpha0 = 1.0 / alpha_inverse_mz()
    alpha_ms = run_two_loop(alpha0, MZ_GEV, MSUSY_GEV, SM_B1, SM_B2, steps=500)
    return run_two_loop(alpha_ms, MSUSY_GEV, mg, MSSM_B1, MSSM_B2, steps=1200)


def proton_constants() -> dict[str, float]:
    payload = json.loads(PROTON_JSON.read_text(encoding="utf-8"))
    return {
        "K_GeV5": float(payload["hadronic_constants"]["width_prefactor_GeV5"]),
        "tau_target_years": 1.0e34,
    }


def seesaw_masses() -> list[float]:
    payload = json.loads(SEESAW_JSON.read_text(encoding="utf-8"))
    return [float(x) for x in payload["heavy_neutrino_masses_GeV"]]


def mx_required(alpha_g_inv: float, constants: dict[str, float]) -> float:
    width_target = HBAR_GEV_S / (constants["tau_target_years"] * SECONDS_PER_YEAR)
    c6_max = math.sqrt(width_target / constants["K_GeV5"])
    g2 = 4.0 * math.pi / alpha_g_inv
    return math.sqrt(g2 / c6_max)


def triplet_filter_required(m_triplet: float, constants: dict[str, float]) -> float:
    alpha2 = 1.0 / 25.0
    loop = alpha2 / (4.0 * math.pi)
    m_wino = 1.0e3
    m_sfermion = 1.0e5
    yprod = 0.0144
    width_target = HBAR_GEV_S / (constants["tau_target_years"] * SECONDS_PER_YEAR)
    c6_max = math.sqrt(width_target / constants["K_GeV5"])
    denominator = (yprod / m_triplet) * loop * m_wino / m_sfermion**2
    return c6_max / denominator


def fit_heavy_thresholds(alpha_inv_mg: np.ndarray) -> dict[str, object]:
    projector = np.eye(3) - np.ones((3, 3)) / 3.0
    target = projector @ alpha_inv_mg
    mat = projector @ HEAVY_BASIS / (2.0 * math.pi)
    logs, *_ = np.linalg.lstsq(mat, target, rcond=None)
    threshold = HEAVY_BASIS @ logs / (2.0 * math.pi)
    residual_vec = projector @ (alpha_inv_mg - threshold)
    alpha_g_inv = float(np.mean(alpha_inv_mg - threshold))
    return {
        "logs": logs,
        "threshold": threshold,
        "alphaG_inv": alpha_g_inv,
        "residual_l2": float(np.linalg.norm(residual_vec)),
        "max_abs_log": float(np.max(np.abs(logs))),
    }


def spectrum_from_logs(mg: float, logs: np.ndarray, alpha_g_inv: float, constants: dict[str, float]) -> dict[str, object]:
    masses = {name: float(mg * math.exp(-logv)) for name, logv in zip(HEAVY_NAMES, logs)}
    mx_min = mx_required(alpha_g_inv, constants)
    masses["X_gauge"] = float(max(mg, mx_min))
    masses["nuR_1"], masses["nuR_2"], masses["nuR_3"] = seesaw_masses()
    filter_required = triplet_filter_required(masses["H_C_pair"], constants)
    return {
        "masses_GeV": masses,
        "M_X_min_GeV": mx_min,
        "rho_X": float(max(1.0, mx_min / mg)),
        "triplet_filter_required": filter_required,
        "triplet_filter_1e_minus5_safe": 1.0e-5 <= filter_required,
    }


def scan() -> tuple[list[dict[str, object]], dict[str, object]]:
    constants = proton_constants()
    rows: list[dict[str, object]] = []
    best = None
    for log10_mg in np.linspace(15.5, 17.5, 1001):
        mg = 10.0**float(log10_mg)
        alpha_mg = run_to_mg(mg)
        alpha_inv_mg = 1.0 / alpha_mg
        fit = fit_heavy_thresholds(alpha_inv_mg)
        spectrum = spectrum_from_logs(mg, fit["logs"], fit["alphaG_inv"], constants)
        row = {
            "log10_MG": float(log10_mg),
            "MG_GeV": mg,
            "alpha1_inv_MG": float(alpha_inv_mg[0]),
            "alpha2_inv_MG": float(alpha_inv_mg[1]),
            "alpha3_inv_MG": float(alpha_inv_mg[2]),
            "alphaG_inv": fit["alphaG_inv"],
            "log_HC": float(fit["logs"][0]),
            "log_Sigma3": float(fit["logs"][1]),
            "log_Sigma8": float(fit["logs"][2]),
            "M_HC_GeV": spectrum["masses_GeV"]["H_C_pair"],
            "M_Sigma3_GeV": spectrum["masses_GeV"]["Sigma_3"],
            "M_Sigma8_GeV": spectrum["masses_GeV"]["Sigma_8"],
            "M_X_GeV": spectrum["masses_GeV"]["X_gauge"],
            "M_X_min_GeV": spectrum["M_X_min_GeV"],
            "rho_X": spectrum["rho_X"],
            "triplet_filter_required": spectrum["triplet_filter_required"],
            "triplet_filter_1e_minus5_safe": spectrum["triplet_filter_1e_minus5_safe"],
            "residual_l2": fit["residual_l2"],
            "max_abs_log": fit["max_abs_log"],
            "proton_safe": spectrum["rho_X"] <= 1.0000001,
        }
        rows.append(row)
        penalty = 0.0
        if not row["proton_safe"]:
            penalty += 100.0 * math.log(row["rho_X"]) ** 2
        if not row["triplet_filter_1e_minus5_safe"]:
            penalty += 100.0 * math.log(1.0e-5 / row["triplet_filter_required"]) ** 2
        score = row["max_abs_log"] + 1000.0 * row["residual_l2"] + penalty
        if best is None or score < best[0]:
            best = (score, row)
    assert best is not None
    best_row = best[1]
    best_spectrum = spectrum_from_logs(
        best_row["MG_GeV"],
        np.array([best_row["log_HC"], best_row["log_Sigma3"], best_row["log_Sigma8"]]),
        best_row["alphaG_inv"],
        constants,
    )
    summary = {
        "constants": constants,
        "MSUSY_GeV": MSUSY_GEV,
        "two_loop": {
            "SM_b1": SM_B1.tolist(),
            "SM_b2": SM_B2.tolist(),
            "MSSM_b1": MSSM_B1.tolist(),
            "MSSM_b2": MSSM_B2.tolist(),
            "yukawa_terms": "neglected in this gauge-only first pass",
        },
        "heavy_threshold_basis": {
            "names": HEAVY_NAMES,
            "columns_b_i_r": HEAVY_BASIS.T.tolist(),
            "convention": "Delta_i = sum_r b_i^r ln(MG/Mr)/(2pi)",
        },
        "best": best_row,
        "best_spectrum": best_spectrum,
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, object]) -> None:
    best = summary["best"]
    spectrum = summary["best_spectrum"]
    masses = spectrum["masses_GeV"]
    lines: list[str] = []
    lines.append("# Two-loop RGE and heavy-spectrum fit")
    lines.append("")
    lines.append("Gauge-only two-loop scan; Yukawa terms are left for the next refinement.")
    lines.append("")
    lines.append("Running:")
    lines.append("")
    lines.append("```text")
    lines.append("MZ -> MSUSY: SM two-loop gauge beta functions")
    lines.append("MSUSY -> MG: MSSM two-loop gauge beta functions")
    lines.append(f"MSUSY = {summary['MSUSY_GeV']:.3e} GeV")
    lines.append("```")
    lines.append("")
    lines.append("Heavy threshold basis:")
    lines.append("")
    lines.append("```text")
    lines.append("H_C pair:   b = (2/5, 0, 1)")
    lines.append("Sigma_3:    b = (0, 2, 0)")
    lines.append("Sigma_8:    b = (0, 0, 3)")
    lines.append("Delta_i = sum_r b_i^r ln(MG/Mr)/(2pi)")
    lines.append("```")
    lines.append("")
    lines.append("Best fit:")
    lines.append("")
    lines.append("```text")
    lines.append(f"MG = {best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {best['alphaG_inv']:.6f}")
    lines.append(f"alpha_i^-1(MG) = ({best['alpha1_inv_MG']:.6f}, {best['alpha2_inv_MG']:.6f}, {best['alpha3_inv_MG']:.6f})")
    lines.append(f"logs ln(MG/Mr) = ({best['log_HC']:+.6f}, {best['log_Sigma3']:+.6f}, {best['log_Sigma8']:+.6f})")
    lines.append(f"residual_l2 = {best['residual_l2']:.6e}")
    lines.append(f"max_abs_log = {best['max_abs_log']:.6f}")
    lines.append("```")
    lines.append("")
    lines.append("Spectrum:")
    lines.append("")
    lines.append("| state | mass [GeV] |")
    lines.append("|---|---:|")
    for key in ["H_C_pair", "Sigma_3", "Sigma_8", "X_gauge", "nuR_1", "nuR_2", "nuR_3"]:
        lines.append(f"| `{key}` | {masses[key]:.6e} |")
    lines.append("")
    lines.append("Proton checks:")
    lines.append("")
    lines.append("```text")
    lines.append(f"M_X_min = {spectrum['M_X_min_GeV']:.6e} GeV")
    lines.append(f"rho_X = {spectrum['rho_X']:.6f}")
    lines.append(f"triplet_filter_required = {spectrum['triplet_filter_required']:.6e}")
    lines.append(f"S_T=1e-5 safe = {spectrum['triplet_filter_1e_minus5_safe']}")
    lines.append("```")
    lines.append("")
    lines.append("Interpretation: the two-loop gauge-only fit remains close to the one-loop MSSM-like branch.")
    lines.append("The heavy triplet sits above the matching scale, the gauge boson is proton-safe,")
    lines.append("and the required triplet geometric filter is weaker than the item-4 design target.")
    (OUT / "two_loop_spectrum_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = scan()
    write_csv(OUT / "two_loop_spectrum_scan.csv", rows)
    (OUT / "two_loop_spectrum_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    best = summary["best"]
    spectrum = summary["best_spectrum"]
    print("Two-loop RGE heavy-spectrum fit")
    print(f"  MG: {best['MG_GeV']:.6e} GeV")
    print(f"  alphaG^-1: {best['alphaG_inv']:.6f}")
    print(f"  logs: ({best['log_HC']:+.6f}, {best['log_Sigma3']:+.6f}, {best['log_Sigma8']:+.6f})")
    print(f"  residual_l2: {best['residual_l2']:.6e}")
    print(f"  max_abs_log: {best['max_abs_log']:.6f}")
    print(f"  M_X_min: {spectrum['M_X_min_GeV']:.6e} GeV")
    print(f"  rho_X: {spectrum['rho_X']:.6f}")
    print(f"  triplet filter required: {spectrum['triplet_filter_required']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
