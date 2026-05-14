#!/usr/bin/env python3
"""Scan a minimal 120_H-like antisymmetric correction for CKM misalignment.

No web lookup is used.  This is an effective flavor audit, not yet a full
Spin(10) Yukawa fit: a complex antisymmetric family matrix A is added only to
the down-sector benchmark

    Y_d -> Y_d + A,        A^T = -A,

while Y_u is held fixed.  The question is whether the exact O(-4) symmetric
core can be repaired by the smallest 120_H-type left-misalignment ingredient
without destroying the down-sector hierarchy.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
OUT = ROOT / "output" / "flavor_120h_scan"


TARGET_MASS = {"down": (1.0e-3, 2.0e-2)}
TARGET_CKM = {"Vus": 0.22499845986972888, "Vcb": 0.04099971935403949, "Vub": 0.0037, "J": 3.097061593767482e-05}


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in m]


def matrix_abs_json(m: np.ndarray) -> list[list[float]]:
    return [[float(abs(z)) for z in row] for row in m]


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def antisymmetric_from_x(x: np.ndarray) -> np.ndarray:
    a12 = x[0] + 1j * x[1]
    a13 = x[2] + 1j * x[3]
    a23 = x[4] + 1j * x[5]
    return np.array(
        [
            [0.0, a12, a13],
            [-a12, 0.0, a23],
            [-a13, -a23, 0.0],
        ],
        dtype=complex,
    )


def ckm_observables(u_left: np.ndarray, d_left: np.ndarray) -> dict[str, float]:
    v = u_left.conjugate().T @ d_left
    j = float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))
    return {
        "Vus": float(abs(v[0, 1])),
        "Vcb": float(abs(v[1, 2])),
        "Vub": float(abs(v[0, 2])),
        "J": abs(j),
    }


def signed_jarlskog(u_left: np.ndarray, d_left: np.ndarray) -> float:
    v = u_left.conjugate().T @ d_left
    return float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))


def ratios(s: np.ndarray) -> tuple[float, float]:
    return float(s[0] / s[2]), float(s[1] / s[2])


def log10_ratio(found: float, target: float) -> float:
    return math.log10(max(found, 1.0e-30) / target)


def objective(x: np.ndarray, y_u: np.ndarray, y_d: np.ndarray, weights: dict[str, float]) -> np.ndarray:
    u_left, _ = left_rotation(y_u)
    d_left, d_s = left_rotation(y_d + antisymmetric_from_x(x))
    obs = ckm_observables(u_left, d_left)
    small, mid = ratios(d_s)
    return np.array(
        [
            weights["ckm"] * log10_ratio(obs["Vus"], TARGET_CKM["Vus"]),
            weights["ckm"] * log10_ratio(obs["Vcb"], TARGET_CKM["Vcb"]),
            weights["ckm"] * log10_ratio(obs["Vub"], TARGET_CKM["Vub"]),
            weights["jarlskog"] * log10_ratio(obs["J"], TARGET_CKM["J"]),
            weights["mass"] * log10_ratio(small, TARGET_MASS["down"][0]),
            weights["mass"] * log10_ratio(mid, TARGET_MASS["down"][1]),
            weights["reg"] * float(np.linalg.norm(x)),
        ],
        dtype=float,
    )


def evaluate(x: np.ndarray, y_u: np.ndarray, y_d: np.ndarray) -> dict[str, Any]:
    a = antisymmetric_from_x(x)
    y_new = y_d + a
    u_left, u_s = left_rotation(y_u)
    d_left, d_s = left_rotation(y_new)
    v = u_left.conjugate().T @ d_left
    obs = ckm_observables(u_left, d_left)
    small, mid = ratios(d_s)
    mass_score = log10_ratio(small, TARGET_MASS["down"][0]) ** 2 + log10_ratio(mid, TARGET_MASS["down"][1]) ** 2
    ckm_score = sum(log10_ratio(obs[key], TARGET_CKM[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = log10_ratio(obs["J"], TARGET_CKM["J"]) ** 2
    return {
        "x": [float(vv) for vv in x],
        "A_120": matrix_json(a),
        "Yd_corrected": matrix_json(y_new),
        "CKM_abs": matrix_abs_json(v),
        "CKM_observables": obs,
        "jarlskog_signed": signed_jarlskog(u_left, d_left),
        "down_singular_values": [float(vv) for vv in d_s],
        "down_mass_ratios": {"small": small, "mid": mid},
        "scores": {
            "mass_log_score": float(mass_score),
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "total_unweighted": float(mass_score + ckm_score + j_score),
        },
        "norms": {
            "A_fro_over_Yd_fro": float(np.linalg.norm(a) / np.linalg.norm(y_d)),
            "A_spectral_over_Yd_spectral": float(np.linalg.norm(a, ord=2) / np.linalg.norm(y_d, ord=2)),
            "A_max_abs": float(np.max(np.abs(a))),
            "Yd_shift_fro": float(np.linalg.norm(a)),
        },
    }


def scan_bound(bound: float, y_u: np.ndarray, y_d: np.ndarray, rng: np.random.Generator, starts: int = 48) -> dict[str, Any]:
    weights = {"ckm": 1.0, "jarlskog": 0.35, "mass": 2.0, "reg": 0.02 / max(bound, 1.0e-9)}
    seeds = [np.zeros(6)]
    seeds.extend(rng.uniform(-bound, bound, size=6) for _ in range(starts - 1))
    best = None
    for seed in seeds:
        result = least_squares(
            objective,
            seed,
            bounds=(-bound * np.ones(6), bound * np.ones(6)),
            args=(y_u, y_d, weights),
            xtol=1.0e-11,
            ftol=1.0e-11,
            gtol=1.0e-11,
            max_nfev=3500,
        )
        evald = evaluate(result.x, y_u, y_d)
        weighted_cost = float(np.sum(objective(result.x, y_u, y_d, weights) ** 2))
        evald["optimizer"] = {
            "success": bool(result.success),
            "cost_weighted_sum_squares": weighted_cost,
            "nfev": int(result.nfev),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(result.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or weighted_cost < best["optimizer"]["cost_weighted_sum_squares"]:
            best = evald
    assert best is not None
    best["component_bound"] = bound
    return best


def write_csv(rows: list[dict[str, Any]]) -> None:
    flat = []
    for row in rows:
        obs = row["CKM_observables"]
        mass = row["down_mass_ratios"]
        score = row["scores"]
        norms = row["norms"]
        opt = row["optimizer"]
        flat.append(
            {
                "component_bound": row["component_bound"],
                "Vus": obs["Vus"],
                "Vcb": obs["Vcb"],
                "Vub": obs["Vub"],
                "J_abs": obs["J"],
                "down_small": mass["small"],
                "down_mid": mass["mid"],
                "ckm_score": score["ckm_magnitude_log_score"],
                "mass_score": score["mass_log_score"],
                "j_score": score["jarlskog_log_score"],
                "A_fro_over_Yd_fro": norms["A_fro_over_Yd_fro"],
                "A_spectral_over_Yd_spectral": norms["A_spectral_over_Yd_spectral"],
                "A_max_abs": norms["A_max_abs"],
                "weighted_cost": opt["cost_weighted_sum_squares"],
                "active_bound_fraction": opt["active_bound_fraction"],
            }
        )
    with (OUT / "scan_120h_ckm_bounds.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0].keys()))
        writer.writeheader()
        writer.writerows(flat)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# 120H CKM misalignment scan")
    lines.append("")
    lines.append("No web lookup was used.  This is an effective down-sector 120H audit.")
    lines.append("")
    lines.append("The scan adds a complex antisymmetric correction `A` to the exact card:")
    lines.append("")
    lines.append("```text")
    lines.append("Y_d -> Y_d + A,  A^T = -A.")
    lines.append("```")
    lines.append("")
    lines.append("| bound | Vus | Vcb | Vub | J | down small | down mid | CKM score | A_F/Yd_F |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["bound_scan"]:
        obs = row["CKM_observables"]
        mass = row["down_mass_ratios"]
        lines.append(
            f"| {row['component_bound']:.2f} | {obs['Vus']:.4e} | {obs['Vcb']:.4e} | "
            f"{obs['Vub']:.4e} | {obs['J']:.4e} | {mass['small']:.4e} | "
            f"{mass['mid']:.4e} | {row['scores']['ckm_magnitude_log_score']:.3e} | "
            f"{row['norms']['A_fro_over_Yd_fro']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
    mass = best["down_mass_ratios"]
    lines.append("")
    lines.append("## Best point")
    lines.append("")
    lines.append(f"Component bound: `{best['component_bound']:.3g}`.")
    lines.append(
        "CKM: "
        f"`|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, |J|={obs['J']:.6e}`."
    )
    lines.append(
        "Down ratios: "
        f"`small={mass['small']:.6e}, mid={mass['mid']:.6e}`."
    )
    lines.append(
        "Shift size: "
        f"`||A||_F/||Y_d||_F={best['norms']['A_fro_over_Yd_fro']:.6e}`, "
        f"`||A||_2/||Y_d||_2={best['norms']['A_spectral_over_Yd_spectral']:.6e}`."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "scan_120h_ckm_misalignment_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = json.loads(CARD.read_text(encoding="utf-8"))
    y_u = cmat(card["yukawa_sectors"]["up"]["Y_matrix_normalized"])
    y_d = cmat(card["yukawa_sectors"]["down"]["Y_matrix_normalized"])
    rng = np.random.default_rng(20260507)
    bounds = [0.02, 0.05, 0.10, 0.20, 0.50, 1.00]
    rows = [scan_bound(bound, y_u, y_d, rng) for bound in bounds]
    best = min(rows, key=lambda row: row["optimizer"]["cost_weighted_sum_squares"])
    verdict = (
        "A down-sector antisymmetric 120H-like correction alone does not repair "
        "the CKM failure.  It preserves the down hierarchy and can tune the "
        "Jarlskog size, but it gets stuck with large |Vcb| and |Vub| even after "
        "a moderate ||A||_F/||Y_d||_F ~= 0.26 deformation.  Therefore the next "
        "test should not be a larger one-sector 120H scan.  It must either impose "
        "SO(10) Clebsch-correlated 120H contributions across u,d,e,nu, or try a "
        "Toeplitz/Kahler canonical-normalization metric before replaying "
        "dimension-five proton decay with fitted rotations."
    )
    payload = {
        "note": "No web lookup used. Effective 120H antisymmetric down-sector CKM scan.",
        "input_card": str(CARD),
        "target_CKM": TARGET_CKM,
        "target_down_mass_ratios": TARGET_MASS["down"],
        "bound_scan": rows,
        "best": best,
        "checks": {
            "scan_completed_with_finite_scores": all(
                np.isfinite(row["scores"]["ckm_magnitude_log_score"]) and np.isfinite(row["scores"]["mass_log_score"])
                for row in rows
            ),
            "best_preserves_down_mass_score_below_0p05": best["scores"]["mass_log_score"] < 0.05,
            "best_still_fails_ckm_score_above_1": best["scores"]["ckm_magnitude_log_score"] > 1.0,
            "best_shift_is_moderate": 0.1 < best["norms"]["A_fro_over_Yd_fro"] < 0.5,
        },
        "verdict": verdict,
    }
    payload["all_checks_ok"] = all(payload["checks"].values())
    (OUT / "scan_120h_ckm_misalignment.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(rows)
    write_report(payload)
    print("120H CKM misalignment scan")
    for row in rows:
        obs = row["CKM_observables"]
        print(
            f"  bound={row['component_bound']:.2f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"ckm_score={row['scores']['ckm_magnitude_log_score']:.3e} "
            f"A/Yd={row['norms']['A_fro_over_Yd_fro']:.3e}"
        )
    print(f"  best bound: {best['component_bound']:.2f}")
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
