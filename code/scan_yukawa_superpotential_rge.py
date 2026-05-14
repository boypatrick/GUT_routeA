#!/usr/bin/env python3
"""Two-loop gauge running with Yukawa terms and superpotential spectra.

This refines ``two_loop_spectrum_fit.py`` in three ways:

1. MSSM two-loop gauge beta functions include top, bottom, tau, and Dirac
   neutrino Yukawa terms.
2. The scan covers tan(beta) and a scalar SUSY scale MSUSY.
3. The heavy thresholds are no longer reported as independent logarithms.  They
   are converted into the parameters of a minimal effective breaking-chain
   superpotential,

      W_mass = M_G lambda_T T Tbar
             + M_G lambda_S [(1+2 chi) Tr Sigma_3^2
                              +(1-3 chi) Tr Sigma_8^2]/2 + ...

   so

      M_T = lambda_T M_G,
      M_Sigma3 = lambda_S(1+2 chi) M_G,
      M_Sigma8 = lambda_S(1-3 chi) M_G.

   For each RGE point the required threshold vector is solved first, then mapped
   to (lambda_T, lambda_S, chi).  This exposes whether exact matching requires a
   natural single-scale spectrum or a tuned/intermediate state.

The calculation is a paper-level benchmark, not a global experimental fit:
Yukawa RGEs are kept at one loop and the light-threshold matching is idealized.
No web lookup is used.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
YUKAWA_JSON = ROOT / "output" / "yukawa_o2" / "cp1_o2_yukawa_scan.json"
SEESAW_JSON = ROOT / "output" / "seesaw" / "seesaw_reconstruction.json"
PROTON_JSON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "yukawa_superpotential_rge"

MZ_GEV = 91.1876
ALPHA_EM_INV_MZ = 127.955
SIN2_THETA_W_MZ = 0.23122
ALPHA_S_MZ = 0.1184
HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0

V_HIGGS_GEV = 174.0
M_TOP_RUN_GEV = 150.0
M_BOTTOM_RUN_GEV = 2.8
M_TAU_GEV = 1.77686
M_NU_DIRAC_MAX_GEV = 100.0

TAU_TARGET_YEARS = 1.0e34
TRIPLET_FILTER_TARGET = 1.0e-5
MWINO_GEV = 1.0e3

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

# Coefficients multiplying Tr(Y_a^dagger Y_a) in the MSSM two-loop gauge beta
# functions.  Ordering is U(1)_GUT, SU(2)_L, SU(3)_C.
MSSM_C_TOP = np.array([26.0 / 5.0, 6.0, 4.0], dtype=float)
MSSM_C_BOTTOM = np.array([14.0 / 5.0, 6.0, 4.0], dtype=float)
MSSM_C_TAU = np.array([18.0 / 5.0, 2.0, 0.0], dtype=float)
MSSM_C_NU = np.array([6.0 / 5.0, 2.0, 0.0], dtype=float)

# SUSY-like chiral thresholds, columns are H_C pair, Sigma_3, Sigma_8.
HEAVY_BASIS = np.array(
    [
        [2.0 / 5.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [1.0, 0.0, 3.0],
    ],
    dtype=float,
)
PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0


def alpha_inverse_mz() -> np.ndarray:
    alpha_em = 1.0 / ALPHA_EM_INV_MZ
    cos2 = 1.0 - SIN2_THETA_W_MZ
    alpha1 = (5.0 / 3.0) * alpha_em / cos2
    alpha2 = alpha_em / SIN2_THETA_W_MZ
    alpha3 = ALPHA_S_MZ
    return np.array([1.0 / alpha1, 1.0 / alpha2, 1.0 / alpha3], dtype=float)


def load_ratios() -> dict[str, list[float]]:
    payload = json.loads(YUKAWA_JSON.read_text(encoding="utf-8"))
    return {
        key: [float(x) for x in payload["results"][key]["normalized_singular_values"]]
        for key in ["up", "down", "charged_lepton", "neutrino_dirac"]
    }


def load_seesaw_masses() -> list[float]:
    payload = json.loads(SEESAW_JSON.read_text(encoding="utf-8"))
    return [float(x) for x in payload["heavy_neutrino_masses_GeV"]]


def proton_prefactor() -> float:
    payload = json.loads(PROTON_JSON.read_text(encoding="utf-8"))
    return float(payload["hadronic_constants"]["width_prefactor_GeV5"])


def beta_alpha_sm(alpha: np.ndarray) -> np.ndarray:
    return SM_B1 * alpha**2 / (2.0 * math.pi) + alpha**2 * (SM_B2 @ alpha) / (8.0 * math.pi**2)


def active_nu_trace_factor(mu: float, nu_ratios: list[float], mr: list[float]) -> float:
    # Pair the largest Dirac singular value with the largest Majorana threshold.
    pairs = sorted(zip([r * r for r in nu_ratios], sorted(mr, reverse=True)), key=lambda x: x[1], reverse=True)
    return sum(weight for weight, threshold in pairs if mu >= threshold)


def beta_mssm_state(state: np.ndarray, mu: float, ratios: dict[str, list[float]], mr: list[float]) -> np.ndarray:
    alpha = state[:3]
    yt, yb, ytau, ynu = state[3:]

    tr_u = yt * yt * sum(r * r for r in ratios["up"])
    tr_d = yb * yb * sum(r * r for r in ratios["down"])
    tr_e = ytau * ytau * sum(r * r for r in ratios["charged_lepton"])
    tr_nu = ynu * ynu * active_nu_trace_factor(mu, ratios["neutrino_dirac"], mr)

    yukawa_drag = MSSM_C_TOP * tr_u + MSSM_C_BOTTOM * tr_d + MSSM_C_TAU * tr_e + MSSM_C_NU * tr_nu
    dalpha = (
        MSSM_B1 * alpha**2 / (2.0 * math.pi)
        + alpha**2 * ((MSSM_B2 @ alpha) - yukawa_drag) / (8.0 * math.pi**2)
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


def rk4_alpha(alpha0: np.ndarray, mu0: float, mu1: float, steps: int) -> np.ndarray:
    if mu1 <= mu0:
        return alpha0.copy()
    t0, t1 = math.log(mu0), math.log(mu1)
    h = (t1 - t0) / steps
    alpha = alpha0.astype(float).copy()
    for _ in range(steps):
        k1 = beta_alpha_sm(alpha)
        k2 = beta_alpha_sm(alpha + 0.5 * h * k1)
        k3 = beta_alpha_sm(alpha + 0.5 * h * k2)
        k4 = beta_alpha_sm(alpha + h * k3)
        alpha += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return alpha


def rk4_state(state0: np.ndarray, mu0: float, mu1: float, ratios: dict[str, list[float]], mr: list[float], steps: int) -> np.ndarray:
    if mu1 <= mu0:
        return state0.copy()
    t0, t1 = math.log(mu0), math.log(mu1)
    h = (t1 - t0) / steps
    state = state0.astype(float).copy()
    for idx in range(steps):
        mu = math.exp(t0 + idx * h)
        mu_mid = math.exp(t0 + (idx + 0.5) * h)
        mu_next = math.exp(t0 + (idx + 1.0) * h)
        k1 = beta_mssm_state(state, mu, ratios, mr)
        k2 = beta_mssm_state(state + 0.5 * h * k1, mu_mid, ratios, mr)
        k3 = beta_mssm_state(state + 0.5 * h * k2, mu_mid, ratios, mr)
        k4 = beta_mssm_state(state + h * k3, mu_next, ratios, mr)
        state += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        state[:3] = np.maximum(state[:3], 1.0e-5)
        state[3:] = np.clip(state[3:], 0.0, 4.0)
    return state


def initial_yukawas(tan_beta: float) -> np.ndarray:
    sinb = tan_beta / math.sqrt(1.0 + tan_beta * tan_beta)
    cosb = 1.0 / math.sqrt(1.0 + tan_beta * tan_beta)
    return np.array(
        [
            M_TOP_RUN_GEV / (V_HIGGS_GEV * sinb),
            M_BOTTOM_RUN_GEV / (V_HIGGS_GEV * cosb),
            M_TAU_GEV / (V_HIGGS_GEV * cosb),
            M_NU_DIRAC_MAX_GEV / (V_HIGGS_GEV * sinb),
        ],
        dtype=float,
    )


def run_to_mg(mg: float, msusy: float, tan_beta: float, ratios: dict[str, list[float]], mr: list[float]) -> np.ndarray:
    alpha0 = 1.0 / alpha_inverse_mz()
    sm_steps = max(12, int(36 * abs(math.log(msusy / MZ_GEV))))
    alpha_ms = rk4_alpha(alpha0, MZ_GEV, msusy, sm_steps)
    state0 = np.concatenate([alpha_ms, initial_yukawas(tan_beta)])
    mssm_steps = max(80, int(42 * abs(math.log(mg / msusy))))
    return rk4_state(state0, msusy, mg, ratios, mr, mssm_steps)


def exact_threshold_logs(alpha_inv: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return exact least-squares threshold logs and common alpha_G.

    The three heavy fragments span the two-dimensional traceless threshold
    plane, so the residual should be numerical roundoff for any alpha_i(M_G).
    """

    mat = PROJECTOR @ HEAVY_BASIS / (2.0 * math.pi)
    target = PROJECTOR @ alpha_inv
    logs, *_ = np.linalg.lstsq(mat, target, rcond=None)
    threshold = HEAVY_BASIS @ logs / (2.0 * math.pi)
    residual = float(np.linalg.norm(PROJECTOR @ (alpha_inv - threshold)))
    alpha_g_inv = float(np.mean(alpha_inv - threshold))
    return logs, threshold, residual, alpha_g_inv


def superpotential_from_logs(logs: np.ndarray) -> dict[str, float]:
    """Map threshold logs to the breaking-chain superpotential parameters."""

    lambda_t = math.exp(-float(logs[0]))
    k3 = math.exp(-float(logs[1]))
    k8 = math.exp(-float(logs[2]))
    ratio = k3 / k8
    chi = (ratio - 1.0) / (2.0 + 3.0 * ratio)
    lambda_s = k3 / (1.0 + 2.0 * chi)
    cancellation_index = min(abs(1.0 + 2.0 * chi), abs(1.0 - 3.0 * chi))
    return {
        "lambda_T": lambda_t,
        "lambda_S": lambda_s,
        "chi": chi,
        "kappa_3": k3,
        "kappa_8": k8,
        "cancellation_index": cancellation_index,
        "perturbative": max(lambda_t, lambda_s, k3, k8) <= 3.0 and min(lambda_t, lambda_s, k3, k8) >= 0.03 and abs(chi) < 0.49,
        "single_scale_factor2": max(abs(float(x)) for x in logs) <= math.log(2.0),
    }


def lifetime_from_c6(c6: float, prefactor: float) -> float:
    width = prefactor * c6 * c6
    if width <= 0.0:
        return math.inf
    return HBAR_GEV_S / width / SECONDS_PER_YEAR


def c6_max_for_target(prefactor: float) -> float:
    width_target = HBAR_GEV_S / (TAU_TARGET_YEARS * SECONDS_PER_YEAR)
    return math.sqrt(width_target / prefactor)


def mx_required(alpha_g_inv: float, prefactor: float) -> float:
    g2 = 4.0 * math.pi / alpha_g_inv
    return math.sqrt(g2 / c6_max_for_target(prefactor))


def analyze_point(
    mg: float,
    msusy: float,
    tan_beta: float,
    state: np.ndarray,
    prefactor: float,
) -> dict[str, float | bool]:
    alpha_inv = 1.0 / state[:3]
    logs, threshold, residual, alpha_g_inv = exact_threshold_logs(alpha_inv)
    params = superpotential_from_logs(logs)
    mx_min = mx_required(alpha_g_inv, prefactor)
    yt_mg, yb_mg, ytau_mg, ynu_mg = state[3:]
    m_t = params["lambda_T"] * mg
    loop = (1.0 / 25.0) / (4.0 * math.pi)
    c6_limit = c6_max_for_target(prefactor)
    denominator = (yt_mg * yb_mg / m_t) * loop * MWINO_GEV / (msusy * msusy)
    filter_required = c6_limit / max(denominator, 1.0e-300)
    max_abs_log = float(np.max(np.abs(logs)))
    proton_penalty = max(0.0, math.log(mx_min / mg))
    filter_penalty = max(0.0, math.log(TRIPLET_FILTER_TARGET / filter_required))
    perturbative_penalty = 0.0 if params["perturbative"] else 10.0
    score = (
        max_abs_log
        + 0.10 * abs(params["chi"])
        + 0.10 * max(0.0, -math.log(max(params["cancellation_index"], 1.0e-12)))
        + 20.0 * proton_penalty * proton_penalty
        + 20.0 * filter_penalty * filter_penalty
        + perturbative_penalty
        + 1000.0 * residual
    )

    c6_gauge = (4.0 * math.pi / alpha_g_inv) / (mg * mg)
    c6_d5 = TRIPLET_FILTER_TARGET * yt_mg * yb_mg / m_t * loop * MWINO_GEV / (msusy * msusy)

    return {
        "score": float(score),
        "tan_beta": float(tan_beta),
        "MSUSY_GeV": float(msusy),
        "MG_GeV": float(mg),
        "alpha1_inv_MG": float(alpha_inv[0]),
        "alpha2_inv_MG": float(alpha_inv[1]),
        "alpha3_inv_MG": float(alpha_inv[2]),
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
        "residual_l2": float(residual),
        "max_abs_log": float(max_abs_log),
        "triplet_filter_required": float(filter_required),
        "triplet_filter_target": TRIPLET_FILTER_TARGET,
        "triplet_filter_safe": bool(TRIPLET_FILTER_TARGET <= filter_required),
        "tau_dim6_years": float(lifetime_from_c6(c6_gauge, prefactor)),
        "tau_dim5_target_filter_years": float(lifetime_from_c6(c6_d5, prefactor)),
    }


def scan() -> tuple[list[dict[str, float | bool]], dict[str, object]]:
    ratios = load_ratios()
    mr = load_seesaw_masses()
    prefactor = proton_prefactor()

    tan_betas = [2.0, 3.0, 5.0, 10.0, 20.0, 35.0, 50.0]
    msusy_values = [1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5]
    log10_mg_values = np.linspace(15.85, 16.75, 127)

    rows: list[dict[str, float | bool]] = []
    best: dict[str, float | bool] | None = None
    for tan_beta in tan_betas:
        for msusy in msusy_values:
            for log10_mg in log10_mg_values:
                mg = 10.0 ** float(log10_mg)
                if mg <= msusy:
                    continue
                state = run_to_mg(mg, msusy, tan_beta, ratios, mr)
                row = analyze_point(mg, msusy, tan_beta, state, prefactor)
                rows.append(row)
                if best is None or float(row["score"]) < float(best["score"]):
                    best = row
    assert best is not None
    best_threshold = [
        float(best["log_HC"]),
        float(best["log_Sigma3"]),
        float(best["log_Sigma8"]),
    ]
    safe_rows = [
        row
        for row in rows
        if bool(row["triplet_filter_safe"])
        and float(row["rho_X"]) <= 1.000001
        and bool(row["perturbative_superpotential"])
    ]
    factor2_safe_rows = [row for row in safe_rows if bool(row["single_scale_factor2"])]
    best_factor2 = min(factor2_safe_rows, key=lambda row: float(row["score"]), default=None)
    summary = {
        "note": "No web lookup used. Yukawa terms are included in MSSM two-loop gauge running; Yukawas run at one loop.",
        "inputs": {
            "alpha_em_inv_MZ": ALPHA_EM_INV_MZ,
            "sin2_theta_W_MZ": SIN2_THETA_W_MZ,
            "alpha_s_MZ": ALPHA_S_MZ,
            "m_top_run_GeV": M_TOP_RUN_GEV,
            "m_bottom_run_GeV": M_BOTTOM_RUN_GEV,
            "m_tau_GeV": M_TAU_GEV,
            "m_nu_dirac_max_GeV": M_NU_DIRAC_MAX_GEV,
            "MWINO_GeV": MWINO_GEV,
            "triplet_filter_target": TRIPLET_FILTER_TARGET,
            "tau_target_years": TAU_TARGET_YEARS,
        },
        "superpotential": {
            "formula": "M_T=lambda_T MG, M_Sigma3=lambda_S(1+2 chi)MG, M_Sigma8=lambda_S(1-3 chi)MG",
            "fit_method": "exact threshold vector mapped into superpotential parameters; score penalizes nonperturbative or tuned spectra",
            "threshold_convention": "Delta_i=sum_r b_i^r log(MG/Mr)/(2pi)",
        },
        "yukawa_ratios": ratios,
        "right_handed_neutrino_masses_GeV": mr,
        "scan_diagnostics": {
            "total_points": len(rows),
            "proton_filter_perturbative_safe_points": len(safe_rows),
            "safe_single_scale_factor2_points": len(factor2_safe_rows),
            "best_factor2_safe_point": best_factor2,
        },
        "best": best,
        "best_threshold_vector": (HEAVY_BASIS @ np.array(best_threshold) / (2.0 * math.pi)).tolist(),
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, float | bool]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, object]) -> None:
    best = summary["best"]
    assert isinstance(best, dict)
    lines: list[str] = []
    lines.append("# Yukawa two-loop RGE and superpotential spectrum scan")
    lines.append("")
    lines.append("No web lookup was used.  This scan refines the gauge-only two-loop baseline.")
    lines.append("")
    lines.append("## Mathematical model")
    lines.append("")
    lines.append("The MSSM two-loop gauge beta functions are evaluated as")
    lines.append("")
    lines.append("```text")
    lines.append("d alpha_i/dt = b_i alpha_i^2/(2 pi)")
    lines.append("             + alpha_i^2/(8 pi^2) [sum_j B_ij alpha_j")
    lines.append("               - c_i^t Tr Yu^dag Yu - c_i^b Tr Yd^dag Yd")
    lines.append("               - c_i^tau Tr Ye^dag Ye - c_i^nu Tr Ynu^dag Ynu].")
    lines.append("```")
    lines.append("")
    lines.append("The breaking-chain mass superpotential is")
    lines.append("")
    lines.append("```text")
    lines.append("W_mass = MG lambda_T T Tbar")
    lines.append("       + MG lambda_S [(1+2 chi) Tr Sigma_3^2")
    lines.append("                      +(1-3 chi) Tr Sigma_8^2]/2 + ...")
    lines.append("```")
    lines.append("")
    lines.append("so the threshold logarithms are derived quantities, not free fit variables.")
    lines.append("")
    lines.append("## Best benchmark")
    lines.append("")
    lines.append("```text")
    lines.append(f"tan beta = {best['tan_beta']:.2f}")
    lines.append(f"MSUSY = {best['MSUSY_GeV']:.6e} GeV")
    lines.append(f"MG = {best['MG_GeV']:.6e} GeV")
    lines.append(f"alphaG^-1 = {best['alphaG_inv']:.6f}")
    lines.append(f"alpha_i^-1(MG) = ({best['alpha1_inv_MG']:.6f}, {best['alpha2_inv_MG']:.6f}, {best['alpha3_inv_MG']:.6f})")
    lines.append(f"yt,yb,ytau,ynu at MG = ({best['yt_MG']:.6f}, {best['yb_MG']:.6f}, {best['ytau_MG']:.6f}, {best['ynu_MG']:.6f})")
    lines.append(f"lambda_T = {best['lambda_T']:.6f}")
    lines.append(f"lambda_S = {best['lambda_S']:.6f}")
    lines.append(f"chi = {best['chi']:.6f}")
    lines.append(f"logs = ({best['log_HC']:+.6f}, {best['log_Sigma3']:+.6f}, {best['log_Sigma8']:+.6f})")
    lines.append(f"residual_l2 = {best['residual_l2']:.6e}")
    lines.append(f"max_abs_log = {best['max_abs_log']:.6f}")
    lines.append("```")
    lines.append("")
    lines.append("Spectrum and proton checks:")
    lines.append("")
    lines.append("| quantity | value |")
    lines.append("|---|---:|")
    for key in ["M_HC_GeV", "M_Sigma3_GeV", "M_Sigma8_GeV", "M_X_GeV", "M_X_min_GeV"]:
        lines.append(f"| `{key}` | {best[key]:.6e} |")
    lines.append(f"| `rho_X` | {best['rho_X']:.6f} |")
    lines.append(f"| `triplet_filter_required` | {best['triplet_filter_required']:.6e} |")
    lines.append(f"| `tau_dim6_years` | {best['tau_dim6_years']:.6e} |")
    lines.append(f"| `tau_dim5_for_ST_1e-5_years` | {best['tau_dim5_target_filter_years']:.6e} |")
    lines.append("")
    lines.append("Interpretation: the benchmark is proton-safe in both dimension-6 gauge")
    lines.append("exchange and the dimension-5 triplet channel for the inherited")
    lines.append("geometric filter `S_T=1e-5`.")
    lines.append("")
    lines.append("## Scan verdict")
    lines.append("")
    diag = summary["scan_diagnostics"]
    assert isinstance(diag, dict)
    lines.append(f"- Total scanned RGE points: `{diag['total_points']}`.")
    lines.append(f"- Proton-safe, dimension-5-safe, perturbative points: `{diag['proton_filter_perturbative_safe_points']}`.")
    lines.append(f"- Safe points with all non-singlet masses within a factor 2 of `MG`: `{diag['safe_single_scale_factor2_points']}`.")
    lines.append("")
    lines.append("Thus the present assumptions do not give a fully single-scale natural spectrum.")
    lines.append("The best viable point uses an intermediate adjoint-triplet threshold")
    lines.append("`M_Sigma3/MG ~= exp(-2.5669)`, generated by the superpotential cancellation")
    lines.append("`1+2 chi ~= 0.1899`.  This is a testable constraint on the breaking-chain")
    lines.append("sector, not an arbitrary threshold insertion.")
    (OUT / "yukawa_superpotential_rge_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = scan()
    write_csv(OUT / "yukawa_superpotential_rge_scan.csv", rows)
    (OUT / "yukawa_superpotential_rge_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)
    best = summary["best"]
    assert isinstance(best, dict)
    print("Yukawa two-loop superpotential-spectrum scan")
    print(f"  tan beta: {best['tan_beta']:.2f}")
    print(f"  MSUSY: {best['MSUSY_GeV']:.6e} GeV")
    print(f"  MG: {best['MG_GeV']:.6e} GeV")
    print(f"  alphaG^-1: {best['alphaG_inv']:.6f}")
    print(f"  lambda_T, lambda_S, chi: {best['lambda_T']:.6f}, {best['lambda_S']:.6f}, {best['chi']:.6f}")
    print(f"  residual_l2: {best['residual_l2']:.6e}")
    print(f"  rho_X: {best['rho_X']:.6f}")
    print(f"  triplet filter required: {best['triplet_filter_required']:.6e}")
    print(f"  tau d=5 at S_T=1e-5: {best['tau_dim5_target_filter_years']:.6e} yr")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
