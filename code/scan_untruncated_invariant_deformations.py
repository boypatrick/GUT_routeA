#!/usr/bin/env python3
"""Route-A scan of allowed untruncated Spin(10) invariant deformations.

No web lookup is used.  This scan starts from the R=200 component Hessian and
adds the allowed renormalizable perturbations found by the invariant audit:

    delta W = 1/2 X_i (dm_ij + ds_ij F54 + dt_ij H210) X_j
            + rho Tr A [B,C].

Here X_i=(A,B,C), the coefficient matrices dm,ds,dt are real symmetric
3-by-3 matrices, and H210 is the aligned weak-volume plus color-fourform
operator.  The scan quantifies how quickly these allowed terms spoil

* the intermediate Sigma_3 threshold;
* X_(6,2,2) lifting at M_G;
* colored Goldstone locking;
* the proton-safe two-loop/RGE branch.

The output is a robustness audit, not a final global vacuum fit.
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

import audit_untruncated_spin10_invariants as inv  # noqa: E402
import compute_210_cubic_matching as cubic210  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "untruncated_invariant_deformation_scan"
BASE_SCAN = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"

BETAS = {
    "Sigma_L": np.array([0.0, 2.0, 0.0], dtype=float),
    "Sigma_R": np.array([6.0 / 5.0, 0.0, 0.0], dtype=float),
    "X_622": np.array([26.0 / 5.0, 6.0, 4.0], dtype=float),
    "Sigma_8_octet": np.array([0.0, 0.0, 3.0], dtype=float),
    "colored_Goldstone_pair": np.array([8.0 / 5.0, 0.0, 1.0], dtype=float),
    "radial_singlet": np.array([0.0, 0.0, 0.0], dtype=float),
    "H_C+H_Cbar": np.array([2.0 / 5.0, 0.0, 1.0], dtype=float),
}

REPS = [
    ("Sigma_L", 3, 0.13012246874641864),
    ("Sigma_R", 3, 1.0),
    ("X_622", 24, 1.0),
    ("Sigma_8_octet", 8, 0.5491736620100269),
    ("colored_Goldstone_pair", 6, 0.0),
    ("radial_singlet", 1, -0.27458683100501345),
]


def color_hessian_full() -> np.ndarray:
    """Embed the verified 210_H^3 PS Hessian into the 45 pair basis."""

    x0 = np.zeros(len(cubic210.COLOR_PAIRS), dtype=float)
    for idx, pair in enumerate(cubic210.COLOR_PAIRS):
        if pair in [(0, 1), (2, 3), (4, 5)]:
            x0[idx] = 1.0

    def w_of_x(x: np.ndarray) -> float:
        a = cubic210.vector_to_twoform(x)
        return 0.5 * float(np.dot(x, x)) - cubic210.invariant_i3(a) / 6.0

    h = 1.0e-4
    n = len(x0)
    hess = np.zeros((n, n), dtype=float)
    w0 = w_of_x(x0)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h
        hess[i, i] = (w_of_x(x0 + ei) - 2.0 * w0 + w_of_x(x0 - ei)) / (h * h)
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = h
            val = (
                w_of_x(x0 + ei + ej)
                - w_of_x(x0 + ei - ej)
                - w_of_x(x0 - ei + ej)
                + w_of_x(x0 - ei - ej)
            ) / (4.0 * h * h)
            hess[i, j] = val
            hess[j, i] = val

    full = np.zeros((45, 45), dtype=float)
    for a, pair_a in enumerate(cubic210.COLOR_PAIRS):
        ia = inv.PAIR_INDEX[pair_a]
        for b, pair_b in enumerate(cubic210.COLOR_PAIRS):
            ib = inv.PAIR_INDEX[pair_b]
            full[ia, ib] = hess[a, b]
    return 0.5 * (full + full.T)


def rep_projectors(hcolor: np.ndarray) -> dict[str, np.ndarray]:
    broad = inv.block_projectors()
    p_color = broad["P_color_15_1_1"]
    color_basis = basis_from_projector(p_color)
    h_sub = color_basis.T @ hcolor @ color_basis
    evals, evecs_sub = np.linalg.eigh(0.5 * (h_sub + h_sub.T))

    def spectral(target: float, atol: float = 2.0e-6) -> np.ndarray:
        cols_sub = evecs_sub[:, np.isclose(evals, target, atol=atol)]
        cols = color_basis @ cols_sub
        return cols @ cols.T

    return {
        "Sigma_L": broad["P_L_1_3_1"],
        "Sigma_R": broad["P_R_1_1_3"],
        "X_622": broad["P_X_6_2_2"],
        "Sigma_8_octet": spectral(2.0),
        "colored_Goldstone_pair": spectral(0.0),
        "radial_singlet": spectral(-1.0),
    }


def block_diag3(op: np.ndarray) -> np.ndarray:
    out = np.zeros((135, 135), dtype=float)
    for idx in range(3):
        start = 45 * idx
        out[start : start + 45, start : start + 45] = op
    return out


def insert_block(big: np.ndarray, i: int, j: int, op: np.ndarray) -> None:
    si = 45 * i
    sj = 45 * j
    big[si : si + 45, sj : sj + 45] += op
    if i != j:
        big[sj : sj + 45, si : si + 45] += op.T


def baseline_operators(card: dict[str, Any]) -> dict[str, Any]:
    pars = card["vev_parameters_and_couplings"]
    states = card["light_and_mediator_eigenvalues"]
    sigma8_matrix = np.array(states["Sigma8_block"]["matrix_over_MG"], dtype=float)
    q8 = float(sigma8_matrix[0, 0])
    eta = float(pars["eta"])
    r_med = float(card["R"])

    i45 = np.eye(45)
    f54 = inv.f54_operator()
    hweak = inv.d_matrix_from_fourform(inv.weak_volume_fourform())
    hcolor = color_hessian_full()
    h210 = hweak + hcolor
    p_color = inv.block_projectors()["P_color_15_1_1"]
    p_noncolor = i45 - p_color

    aa = p_noncolor @ (
        float(pars["mu"]) * i45 + float(pars["a54"]) * f54 + float(pars["b210"]) * hweak
    ) @ p_noncolor
    aa += p_color @ (0.5 * q8 * hcolor) @ p_color
    ab = eta * (f54 - 2.0 * i45)
    ac = eta * (f54 + (4.0 / 3.0) * i45)
    bc = r_med * i45

    h0 = np.zeros((135, 135), dtype=float)
    insert_block(h0, 0, 0, aa)
    insert_block(h0, 0, 1, ab)
    insert_block(h0, 0, 2, ac)
    insert_block(h0, 1, 2, bc)

    c0_over_a0 = -eta * (-4.0 / 3.0 - 2.0) / r_med
    ad_a0 = inv.adjoint_action_operator(inv.color_adj_vev_10d())

    return {
        "I45": i45,
        "F54": f54,
        "H210": h210,
        "Hcolor": hcolor,
        "H0": 0.5 * (h0 + h0.T),
        "rep_projectors": rep_projectors(hcolor),
        "ad_A0": ad_a0,
        "C0_over_A0": c0_over_a0,
        "q8": q8,
    }


def symmetric_unit(rng: np.random.Generator) -> np.ndarray:
    raw = rng.normal(size=(3, 3))
    mat = 0.5 * (raw + raw.T)
    norm = np.linalg.norm(mat)
    return mat / max(norm, 1.0e-300)


def deformation_matrix(
    ops: dict[str, Any],
    dm: np.ndarray,
    ds: np.ndarray,
    dt: np.ndarray,
    rho: float,
) -> np.ndarray:
    h = np.array(ops["H0"], copy=True)
    for i in range(3):
        for j in range(i, 3):
            op = dm[i, j] * ops["I45"] + ds[i, j] * ops["F54"] + dt[i, j] * ops["H210"]
            insert_block(h, i, j, op)

    # rho Tr A[B,C].  The A0 vev deforms the B-C block.  The Schur-flat
    # C0 parallel vev gives a smaller A-B block of the same adjoint action.
    ad = ops["ad_A0"]
    insert_block(h, 1, 2, rho * ad)
    insert_block(h, 0, 1, rho * float(ops["C0_over_A0"]) * ad)
    return 0.5 * (h + h.T)


def rep_weights(vec: np.ndarray, projectors: dict[str, np.ndarray]) -> dict[str, float]:
    weights = {}
    for name, p in projectors.items():
        total = 0.0
        for block in range(3):
            part = vec[45 * block : 45 * (block + 1)]
            total += float(part @ p @ part)
        weights[name] = total
    return weights


def basis_from_projector(projector: np.ndarray, atol: float = 1.0e-8) -> np.ndarray:
    evals, evecs = np.linalg.eigh(projector)
    return evecs[:, evals > 1.0 - atol]


def lift_basis_to_abc(basis45: np.ndarray) -> np.ndarray:
    dim = basis45.shape[1]
    out = np.zeros((135, 3 * dim), dtype=float)
    for block in range(3):
        out[45 * block : 45 * (block + 1), block * dim : (block + 1) * dim] = basis45
    return out


def restricted_eigenvalues(h: np.ndarray, projector45: np.ndarray) -> np.ndarray:
    basis = lift_basis_to_abc(basis_from_projector(projector45))
    restricted = basis.T @ h @ basis
    return np.linalg.eigvalsh(0.5 * (restricted + restricted.T))


def fine_offblock_norm(h: np.ndarray, projectors45: dict[str, np.ndarray]) -> dict[str, float]:
    projectors135 = {name: block_diag3(p) for name, p in projectors45.items()}
    off = np.zeros_like(h)
    names = list(projectors135)
    for a, name_a in enumerate(names):
        for b, name_b in enumerate(names):
            if a == b:
                continue
            off += projectors135[name_a] @ h @ projectors135[name_b]
    return {
        "max_abs": float(np.max(np.abs(off))),
        "frobenius": float(np.linalg.norm(off)),
    }


def classify_eigenstates(evals: np.ndarray, evecs: np.ndarray, projectors: dict[str, np.ndarray]) -> dict[str, list[tuple[float, int]]]:
    classes = {name: [] for name in projectors}
    for idx, val in enumerate(evals):
        weights = rep_weights(evecs[:, idx], projectors)
        name = max(weights, key=weights.get)
        classes[name].append((float(val), idx))
    return classes


def select_light(values: list[tuple[float, int]], count: int, target: float) -> list[tuple[float, int]]:
    return sorted(values, key=lambda row: abs(row[0] - target))[:count]


def geom_log_from_kappas(kappas: list[float], floor: float = 1.0e-16) -> float:
    return float(np.mean([-math.log(max(abs(k), floor)) for k in kappas]))


def spectrum_metrics(h: np.ndarray, ops: dict[str, Any]) -> dict[str, Any]:
    reps_out: dict[str, Any] = {}
    heavy_threshold = np.zeros(3, dtype=float)
    fine_off = fine_offblock_norm(h, ops["rep_projectors"])

    for name, dim, target in REPS:
        values = [(float(v), idx) for idx, v in enumerate(restricted_eigenvalues(h, ops["rep_projectors"][name]))]
        light = select_light(values, dim, target)
        light_indices = {idx for _val, idx in light}
        light_vals = [val for val, _idx in light]
        heavy_vals = [val for val, idx in values if idx not in light_indices]
        light_log = geom_log_from_kappas(light_vals)
        heavy_log_sum = sum(-math.log(max(abs(v), 1.0e-16)) for v in heavy_vals)
        heavy_threshold += BETAS[name] * (heavy_log_sum / max(dim, 1)) / (2.0 * math.pi)
        reps_out[name] = {
            "count_total": len(values),
            "light_values": light_vals,
            "light_abs_mean": float(np.mean(np.abs(light_vals))) if light_vals else float("nan"),
            "light_abs_max": float(np.max(np.abs(light_vals))) if light_vals else float("nan"),
            "light_signed_mean": float(np.mean(light_vals)) if light_vals else float("nan"),
            "light_log_mean": light_log,
            "heavy_log_sum_per_component": heavy_log_sum / max(dim, 1),
        }

    k3 = reps_out["Sigma_L"]["light_abs_mean"]
    k8 = reps_out["Sigma_8_octet"]["light_abs_mean"]
    extra_unwanted = (
        BETAS["Sigma_R"] * reps_out["Sigma_R"]["light_log_mean"]
        + BETAS["X_622"] * reps_out["X_622"]["light_log_mean"]
    ) / (2.0 * math.pi)
    gold_k = reps_out["colored_Goldstone_pair"]["light_abs_max"]
    gold_damage = BETAS["colored_Goldstone_pair"] * (-math.log(max(gold_k, 1.0e-16))) / (
        2.0 * math.pi
    )
    return {
        "reps": reps_out,
        "heavy_threshold_delta": heavy_threshold.tolist(),
        "heavy_threshold_projected_l2": float(np.linalg.norm(base.PROJECTOR @ heavy_threshold)),
        "kappa3_eff": k3,
        "kappa8_eff": k8,
        "log_sigma3_eff": reps_out["Sigma_L"]["light_log_mean"],
        "log_sigma8_eff": reps_out["Sigma_8_octet"]["light_log_mean"],
        "x622_unwanted_projected_l2": float(
            np.linalg.norm(base.PROJECTOR @ (BETAS["X_622"] * reps_out["X_622"]["light_log_mean"] / (2.0 * math.pi)))
        ),
        "sigmaR_unwanted_projected_l2": float(
            np.linalg.norm(base.PROJECTOR @ (BETAS["Sigma_R"] * reps_out["Sigma_R"]["light_log_mean"] / (2.0 * math.pi)))
        ),
        "unwanted_threshold_delta": extra_unwanted.tolist(),
        "unwanted_threshold_projected_l2": float(np.linalg.norm(base.PROJECTOR @ extra_unwanted)),
        "goldstone_max_abs_kappa": gold_k,
        "goldstone_if_physical_projected_l2": float(np.linalg.norm(base.PROJECTOR @ gold_damage)),
        "fine_rep_offblock_max_abs": fine_off["max_abs"],
        "fine_rep_offblock_frobenius": fine_off["frobenius"],
    }


def load_cached_grid() -> dict[str, np.ndarray]:
    rows: list[dict[str, str]] = []
    with BASE_SCAN.open(newline="", encoding="utf-8") as handle:
        rows.extend(csv.DictReader(handle))
    return {
        "MG_GeV": np.array([float(r["MG_GeV"]) for r in rows]),
        "MSUSY_GeV": np.array([float(r["MSUSY_GeV"]) for r in rows]),
        "tan_beta": np.array([float(r["tan_beta"]) for r in rows]),
        "alpha_inv": np.array(
            [
                [float(r["alpha1_inv_MG"]), float(r["alpha2_inv_MG"]), float(r["alpha3_inv_MG"])]
                for r in rows
            ],
            dtype=float,
        ),
        "yt": np.array([float(r["yt_MG"]) for r in rows]),
        "yb": np.array([float(r["yb_MG"]) for r in rows]),
    }


def fixed_spectrum_rge_scan(metrics: dict[str, Any], grid: dict[str, np.ndarray], prefactor: float) -> dict[str, Any]:
    alpha_inv = grid["alpha_inv"]
    delta_med = np.array(metrics["heavy_threshold_delta"], dtype=float)
    unwanted = np.array(metrics["unwanted_threshold_delta"], dtype=float)
    fixed_chiral = (
        BETAS["Sigma_L"] * float(metrics["log_sigma3_eff"])
        + BETAS["Sigma_8_octet"] * float(metrics["log_sigma8_eff"])
    ) / (2.0 * math.pi)
    fixed = delta_med + unwanted + fixed_chiral
    target = (alpha_inv - fixed) @ base.PROJECTOR.T
    mat = base.PROJECTOR @ BETAS["H_C+H_Cbar"] / (2.0 * math.pi)
    denom = float(mat @ mat)
    log_hc = target @ mat / denom
    hc_threshold = np.outer(log_hc, BETAS["H_C+H_Cbar"] / (2.0 * math.pi))
    residual_vec = (alpha_inv - fixed - hc_threshold) @ base.PROJECTOR.T
    residual = np.linalg.norm(residual_vec, axis=1)
    alpha_g_inv = np.mean(alpha_inv - fixed - hc_threshold, axis=1)

    lambda_t = np.exp(-log_hc)
    k3 = max(float(metrics["kappa3_eff"]), 1.0e-300)
    k8 = max(float(metrics["kappa8_eff"]), 1.0e-300)
    ratio = k3 / k8
    chi = (ratio - 1.0) / (2.0 + 3.0 * ratio)
    lambda_s = k3 / (1.0 + 2.0 * chi)
    cancellation = min(abs(1.0 + 2.0 * chi), abs(1.0 - 3.0 * chi))
    perturbative = (
        max(k3, k8, lambda_s, float(np.min(lambda_t))) <= 3.0
        and min(k3, k8, lambda_s, float(np.max(lambda_t))) >= 0.03
        and abs(chi) < 0.49
    )

    mx_min = np.array([base.mx_required(float(a), prefactor) for a in alpha_g_inv])
    c6_gauge = (4.0 * math.pi / alpha_g_inv) / (grid["MG_GeV"] * grid["MG_GeV"])
    tau_d6 = np.array([base.lifetime_from_c6(float(c), prefactor) for c in c6_gauge])
    m_t = lambda_t * grid["MG_GeV"]
    loop = (1.0 / 25.0) / (4.0 * math.pi)
    c6_limit = base.c6_max_for_target(prefactor)
    denom_d5 = (grid["yt"] * grid["yb"] / m_t) * loop * base.MWINO_GEV / (
        grid["MSUSY_GeV"] * grid["MSUSY_GeV"]
    )
    filter_required = c6_limit / np.maximum(denom_d5, 1.0e-300)
    c6_d5 = base.TRIPLET_FILTER_TARGET * grid["yt"] * grid["yb"] / m_t * loop * base.MWINO_GEV / (
        grid["MSUSY_GeV"] * grid["MSUSY_GeV"]
    )
    tau_d5 = np.array([base.lifetime_from_c6(float(c), prefactor) for c in c6_d5])

    proton_safe = grid["MG_GeV"] >= mx_min
    d5_safe = base.TRIPLET_FILTER_TARGET <= filter_required
    residual_ok = residual < 1.0e-3
    all_safe = proton_safe & d5_safe & residual_ok & perturbative
    score = (
        residual * 1000.0
        + np.maximum(0.0, np.log(mx_min / grid["MG_GeV"])) ** 2 * 20.0
        + np.maximum(0.0, np.log(base.TRIPLET_FILTER_TARGET / filter_required)) ** 2 * 20.0
        + np.abs(log_hc)
        + 0.1 * abs(chi)
    )
    idx = int(np.argmin(score))
    safe_count = int(np.sum(all_safe))
    best_safe_idx = int(np.argmin(np.where(all_safe, score, np.inf))) if safe_count else -1
    pick = best_safe_idx if safe_count else idx
    return {
        "best_index": idx,
        "best_safe_index": best_safe_idx,
        "safe_count_residual_lt_1e_3": safe_count,
        "best_score": float(score[idx]),
        "best_residual_l2": float(residual[idx]),
        "best_safe_residual_l2": float(residual[pick]),
        "best_MG_GeV": float(grid["MG_GeV"][pick]),
        "best_alphaG_inv": float(alpha_g_inv[pick]),
        "best_log_HC": float(log_hc[pick]),
        "best_M_HC_GeV": float(lambda_t[pick] * grid["MG_GeV"][pick]),
        "best_tau_dim6_years": float(tau_d6[pick]),
        "best_tau_dim5_ST_1e_5_years": float(tau_d5[pick]),
        "best_filter_required": float(filter_required[pick]),
        "proton_safe_at_pick": bool(proton_safe[pick]),
        "d5_safe_at_pick": bool(d5_safe[pick]),
        "residual_ok_at_pick": bool(residual_ok[pick]),
        "perturbative_fixed_spectrum": bool(perturbative),
        "chi_from_fixed_kappas": float(chi),
        "lambda_S_from_fixed_kappas": float(lambda_s),
        "cancellation_index": float(cancellation),
    }


def sample_deformations() -> list[dict[str, Any]]:
    rng = np.random.default_rng(20260506)
    amplitudes = [
        0.0,
        1.0e-6,
        3.0e-6,
        1.0e-5,
        3.0e-5,
        1.0e-4,
        3.0e-4,
        1.0e-3,
        3.0e-3,
        1.0e-2,
        3.0e-2,
        1.0e-1,
        3.0e-1,
        1.0,
    ]
    rows: list[dict[str, Any]] = []
    rows.append({"class": "baseline", "amplitude": 0.0, "dm": np.zeros((3, 3)), "ds": np.zeros((3, 3)), "dt": np.zeros((3, 3)), "rho": 0.0})
    for amp in amplitudes[1:]:
        for cls in ["m", "sigma", "tau", "rho", "all"]:
            samples = 24 if cls != "rho" else 8
            for sample in range(samples):
                dm = np.zeros((3, 3))
                ds = np.zeros((3, 3))
                dt = np.zeros((3, 3))
                rho = 0.0
                if cls in {"m", "all"}:
                    dm = amp * symmetric_unit(rng)
                if cls in {"sigma", "all"}:
                    ds = amp * symmetric_unit(rng)
                if cls in {"tau", "all"}:
                    dt = amp * symmetric_unit(rng)
                if cls in {"rho", "all"}:
                    rho = amp * float(rng.choice([-1.0, 1.0]))
                rows.append({"class": cls, "amplitude": amp, "sample": sample, "dm": dm, "ds": ds, "dt": dt, "rho": rho})
    return rows


def evaluate_sample(sample: dict[str, Any], ops: dict[str, Any], grid: dict[str, np.ndarray], prefactor: float) -> dict[str, Any]:
    h = deformation_matrix(ops, sample["dm"], sample["ds"], sample["dt"], float(sample["rho"]))
    metrics = spectrum_metrics(h, ops)
    if "baseline_raw_heavy_delta" in ops:
        raw = np.array(metrics["heavy_threshold_delta"], dtype=float)
        calibrated = (
            np.array(ops["card_heavy_delta"], dtype=float)
            + raw
            - np.array(ops["baseline_raw_heavy_delta"], dtype=float)
        )
        metrics["heavy_threshold_delta"] = calibrated.tolist()
        metrics["heavy_threshold_projected_l2"] = float(
            np.linalg.norm(base.PROJECTOR @ calibrated)
        )
    rge = fixed_spectrum_rge_scan(metrics, grid, prefactor)
    base_k3 = REPS[0][2]
    base_k8 = REPS[3][2]
    safe = (
        metrics["goldstone_max_abs_kappa"] < 1.0e-6
        and metrics["x622_unwanted_projected_l2"] < 5.022738709841473e-4
        and metrics["fine_rep_offblock_max_abs"] < 1.0e-10
        and rge["safe_count_residual_lt_1e_3"] > 0
    )
    return {
        "class": sample["class"],
        "amplitude": float(sample["amplitude"]),
        "sample": int(sample.get("sample", 0)),
        "rho": float(sample["rho"]),
        "kappa3_eff": float(metrics["kappa3_eff"]),
        "kappa8_eff": float(metrics["kappa8_eff"]),
        "sigma3_rel_shift": abs(float(metrics["kappa3_eff"]) - base_k3) / base_k3,
        "sigma8_rel_shift": abs(float(metrics["kappa8_eff"]) - base_k8) / base_k8,
        "x622_unwanted_projected_l2": float(metrics["x622_unwanted_projected_l2"]),
        "sigmaR_unwanted_projected_l2": float(metrics["sigmaR_unwanted_projected_l2"]),
        "goldstone_max_abs_kappa": float(metrics["goldstone_max_abs_kappa"]),
        "goldstone_if_physical_projected_l2": float(metrics["goldstone_if_physical_projected_l2"]),
        "fine_rep_offblock_max_abs": float(metrics["fine_rep_offblock_max_abs"]),
        "heavy_threshold_projected_l2": float(metrics["heavy_threshold_projected_l2"]),
        "rge_best_residual_l2": float(rge["best_residual_l2"]),
        "rge_safe_count": int(rge["safe_count_residual_lt_1e_3"]),
        "best_MG_GeV": float(rge["best_MG_GeV"]),
        "best_tau_dim6_years": float(rge["best_tau_dim6_years"]),
        "best_tau_dim5_ST_1e_5_years": float(rge["best_tau_dim5_ST_1e_5_years"]),
        "safe_routeA": bool(safe),
    }


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for cls in sorted({row["class"] for row in rows}):
        cls_rows = [row for row in rows if row["class"] == cls]
        amp_rows = []
        for amp in sorted({row["amplitude"] for row in cls_rows}):
            selected = [row for row in cls_rows if row["amplitude"] == amp]
            amp_rows.append(
                {
                    "amplitude": amp,
                    "samples": len(selected),
                    "safe_fraction": float(np.mean([row["safe_routeA"] for row in selected])),
                    "median_goldstone_kappa": float(np.median([row["goldstone_max_abs_kappa"] for row in selected])),
                    "median_x622_l2": float(np.median([row["x622_unwanted_projected_l2"] for row in selected])),
                    "median_best_residual": float(np.median([row["rge_best_residual_l2"] for row in selected])),
                    "median_sigma3_rel_shift": float(np.median([row["sigma3_rel_shift"] for row in selected])),
                    "median_heavy_threshold_l2": float(np.median([row["heavy_threshold_projected_l2"] for row in selected])),
                    "median_fine_offblock": float(np.median([row["fine_rep_offblock_max_abs"] for row in selected])),
                }
            )
        safe_amps = [row["amplitude"] for row in amp_rows if row["safe_fraction"] >= 0.5]
        summary[cls] = {
            "by_amplitude": amp_rows,
            "largest_amplitude_with_safe_fraction_ge_0p5": max(safe_amps) if safe_amps else None,
        }
    return summary


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "class",
        "amplitude",
        "sample",
        "rho",
        "kappa3_eff",
        "kappa8_eff",
        "sigma3_rel_shift",
        "sigma8_rel_shift",
        "x622_unwanted_projected_l2",
        "sigmaR_unwanted_projected_l2",
        "goldstone_max_abs_kappa",
        "goldstone_if_physical_projected_l2",
        "fine_rep_offblock_max_abs",
        "heavy_threshold_projected_l2",
        "rge_best_residual_l2",
        "rge_safe_count",
        "best_MG_GeV",
        "best_tau_dim6_years",
        "best_tau_dim5_ST_1e_5_years",
        "safe_routeA",
    ]
    with (OUT / "untruncated_invariant_deformation_scan.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{key: row[key] for key in fields} for row in rows])


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Route-A scan of untruncated invariant deformations")
    lines.append("")
    lines.append("No web lookup was used.  The scan adds the allowed coefficients")
    lines.append("`dm_ij`, `ds_ij`, `dt_ij`, and `rho_ABC` to the R=200 component Hessian,")
    lines.append("then replays a fixed-spectrum threshold fit on the cached two-loop grid.")
    lines.append("")
    lines.append("## Deformation model")
    lines.append("")
    lines.append("```text")
    lines.append("delta W = 1/2 X_i (dm_ij + ds_ij F54 + dt_ij H210) X_j")
    lines.append("        + rho Tr A[B,C],     X_i=(A,B,C).")
    lines.append("```")
    lines.append("")
    lines.append("Here `H210` includes the weak-volume 210 operator and the color-fourform")
    lines.append("PS Hessian with eigenvalues `-1,0,2` on radial/Goldstone/octet.")
    lines.append("")
    lines.append("## Robustness summary")
    lines.append("")
    lines.append("A sample is counted as Route-A safe when")
    lines.append("")
    lines.append("```text")
    lines.append("max |kappa_Goldstone| < 1e-6,")
    lines.append("||P Delta_X622|| < 5.0227387e-4,")
    lines.append("fine-representation off-block max < 1e-10,")
    lines.append("and at least one cached two-loop point has residual < 1e-3")
    lines.append("with d=6 and d=5 proton bounds satisfied.")
    lines.append("```")
    lines.append("")
    lines.append("| class | largest amp with safe fraction >= 1/2 | first unsafe median signal |")
    lines.append("|---|---:|---|")
    for cls, data in payload["summary"].items():
        largest = data["largest_amplitude_with_safe_fraction_ge_0p5"]
        largest_text = "none" if largest is None else f"{largest:.1e}"
        first_bad = next((row for row in data["by_amplitude"] if row["safe_fraction"] < 0.5), None)
        if first_bad is None:
            signal = "not reached"
        else:
            signal = (
                f"amp={first_bad['amplitude']:.1e}, "
                f"gold={first_bad['median_goldstone_kappa']:.2e}, "
                f"X={first_bad['median_x622_l2']:.2e}, "
                f"off={first_bad['median_fine_offblock']:.2e}, "
                f"res={first_bad['median_best_residual']:.2e}"
            )
        lines.append(f"| `{cls}` | {largest_text} | {signal} |")
    lines.append("")
    lines.append("## Amplitude table")
    lines.append("")
    for cls, data in payload["summary"].items():
        lines.append(f"### {cls}")
        lines.append("")
        lines.append("| amp | safe frac | median Goldstone kappa | median X622 l2 | median RGE residual | median Sigma3 shift |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        for row in data["by_amplitude"]:
            lines.append(
                f"| {row['amplitude']:.1e} | {row['safe_fraction']:.2f} | "
                f"{row['median_goldstone_kappa']:.3e} | {row['median_x622_l2']:.3e} | "
                f"{row['median_best_residual']:.3e} | {row['median_sigma3_rel_shift']:.3e} |"
            )
        lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("The scan should be read as a local robustness bound around the R=200 card.")
    lines.append("If a coefficient class fails at amplitude epsilon, it does not prove the")
    lines.append("model impossible; it means the full vacuum fit must either tune that")
    lines.append("coefficient below epsilon, correlate it with other entries, or forbid it by")
    lines.append("a mediator grading/R-symmetry.")
    (OUT / "untruncated_invariant_deformation_scan_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = comp.load_card()
    ops = baseline_operators(card)
    baseline_raw = spectrum_metrics(ops["H0"], ops)
    ops["baseline_raw_heavy_delta"] = baseline_raw["heavy_threshold_delta"]
    ops["card_heavy_delta"] = card["mediator_heavy_threshold"]["delta_full"]
    grid = load_cached_grid()
    prefactor = base.proton_prefactor()
    rows = [evaluate_sample(sample, ops, grid, prefactor) for sample in sample_deformations()]
    payload = {
        "note": "No web lookup used. Route-A scan of allowed untruncated Spin(10) invariant deformations.",
        "source_card": str(comp.CARD),
        "samples": len(rows),
        "safety_definition": {
            "goldstone_max_abs_kappa": "< 1e-6",
            "X622_projected_l2": "< baseline finite mediator projected l2 = 5.0227387e-4",
            "fine_rep_offblock_max_abs": "< 1e-10",
            "RGE": "at least one cached two-loop point with residual < 1e-3 and proton d=6/d=5 safe",
        },
        "baseline_row": rows[0],
        "summary": summarize(rows),
    }
    (OUT / "untruncated_invariant_deformation_scan_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(payload)

    print("Route-A untruncated invariant deformation scan")
    print(f"  samples: {len(rows)}")
    for cls, data in payload["summary"].items():
        largest = data["largest_amplitude_with_safe_fraction_ge_0p5"]
        print(f"  {cls}: largest amp with safe fraction >= 1/2 = {largest}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
