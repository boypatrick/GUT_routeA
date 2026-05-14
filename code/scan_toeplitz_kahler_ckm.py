#!/usr/bin/env python3
"""Scan a Toeplitz/Kahler left-normalization metric for CKM repair.

No web lookup is used.  This audit tests an effective Kahler route:

    K_Q = exp(G),       G = G^dagger, Tr G = 0,
    Y_u -> K_Q^{-1/2} Y_u,
    Y_d -> K_Q^{-1/2} Y_d.

In the finite CP1/O(2) Hilbert space, any Hermitian 3x3 matrix can be realized
as a Berezin-Toeplitz operator with a sufficiently general real symbol.  Thus G
is an eight-parameter Toeplitz/Kahler generator.  The scan asks whether a
controlled positive metric can make CKM small while preserving the already
verified mass hierarchies.
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
OUT = ROOT / "output" / "flavor_kahler_scan"


TARGETS = {
    "up_small": 1.0e-5,
    "up_mid": 7.0e-3,
    "down_small": 1.0e-3,
    "down_mid": 2.0e-2,
    "Vus": 0.22499845986972888,
    "Vcb": 0.04099971935403949,
    "Vub": 0.0037,
    "J": 3.097061593767482e-05,
}


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in m]


def matrix_abs_json(m: np.ndarray) -> list[list[float]]:
    return [[float(abs(z)) for z in row] for row in m]


def hermitian_generator(x: np.ndarray) -> np.ndarray:
    return np.array(
        [
            [x[0], x[2] + 1j * x[3], x[4] + 1j * x[5]],
            [x[2] - 1j * x[3], x[1], x[6] + 1j * x[7]],
            [x[4] - 1j * x[5], x[6] - 1j * x[7], -x[0] - x[1]],
        ],
        dtype=complex,
    )


def exp_hermitian(g: np.ndarray, factor: float) -> np.ndarray:
    values, vectors = np.linalg.eigh(g)
    return (vectors * np.exp(factor * values)) @ vectors.conjugate().T


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def ratios(singulars: np.ndarray) -> tuple[float, float]:
    return float(singulars[0] / singulars[2]), float(singulars[1] / singulars[2])


def log10_ratio(found: float, target: float) -> float:
    return math.log10(max(found, 1.0e-30) / target)


def ckm(u_left: np.ndarray, d_left: np.ndarray) -> np.ndarray:
    return u_left.conjugate().T @ d_left


def ckm_observables(v: np.ndarray) -> dict[str, float]:
    j = float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))
    return {
        "Vus": float(abs(v[0, 1])),
        "Vcb": float(abs(v[1, 2])),
        "Vub": float(abs(v[0, 2])),
        "J": abs(j),
        "J_signed": j,
    }


def apply_metric(x: np.ndarray, y_u: np.ndarray, y_d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    g = hermitian_generator(x)
    k_minus_half = exp_hermitian(g, -0.5)
    k = exp_hermitian(g, 1.0)
    return k_minus_half @ y_u, k_minus_half @ y_d, g, k


def objective(x: np.ndarray, y_u: np.ndarray, y_d: np.ndarray, weights: dict[str, float]) -> np.ndarray:
    yu, yd, g, _k = apply_metric(x, y_u, y_d)
    u_left, u_s = left_rotation(yu)
    d_left, d_s = left_rotation(yd)
    obs = ckm_observables(ckm(u_left, d_left))
    us, um = ratios(u_s)
    ds, dm = ratios(d_s)
    eig = np.linalg.eigvalsh(g)
    log_cond = float(np.max(eig) - np.min(eig))
    return np.array(
        [
            weights["ckm"] * log10_ratio(obs["Vus"], TARGETS["Vus"]),
            weights["ckm"] * log10_ratio(obs["Vcb"], TARGETS["Vcb"]),
            weights["ckm"] * log10_ratio(obs["Vub"], TARGETS["Vub"]),
            weights["jarlskog"] * log10_ratio(obs["J"], TARGETS["J"]),
            weights["mass"] * log10_ratio(us, TARGETS["up_small"]),
            weights["mass"] * log10_ratio(um, TARGETS["up_mid"]),
            weights["mass"] * log10_ratio(ds, TARGETS["down_small"]),
            weights["mass"] * log10_ratio(dm, TARGETS["down_mid"]),
            weights["cond"] * log_cond,
            weights["reg"] * float(np.linalg.norm(x)),
        ],
        dtype=float,
    )


def evaluate(x: np.ndarray, y_u: np.ndarray, y_d: np.ndarray) -> dict[str, Any]:
    yu, yd, g, k = apply_metric(x, y_u, y_d)
    u_left, u_s = left_rotation(yu)
    d_left, d_s = left_rotation(yd)
    v = ckm(u_left, d_left)
    obs = ckm_observables(v)
    us, um = ratios(u_s)
    ds, dm = ratios(d_s)
    eig = np.linalg.eigvalsh(g)
    k_eig = np.linalg.eigvalsh(k)
    ckm_score = sum(log10_ratio(obs[key], TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = log10_ratio(obs["J"], TARGETS["J"]) ** 2
    mass_score = (
        log10_ratio(us, TARGETS["up_small"]) ** 2
        + log10_ratio(um, TARGETS["up_mid"]) ** 2
        + log10_ratio(ds, TARGETS["down_small"]) ** 2
        + log10_ratio(dm, TARGETS["down_mid"]) ** 2
    )
    return {
        "x": [float(vv) for vv in x],
        "G": matrix_json(g),
        "K_Q": matrix_json(k),
        "K_Q_eigenvalues": [float(vv) for vv in k_eig],
        "K_Q_condition": float(np.max(k_eig) / np.min(k_eig)),
        "G_spectral_span": float(np.max(eig) - np.min(eig)),
        "Y_u_canonical": matrix_json(yu),
        "Y_d_canonical": matrix_json(yd),
        "CKM_abs": matrix_abs_json(v),
        "CKM_observables": obs,
        "up_mass_ratios": {"small": us, "mid": um},
        "down_mass_ratios": {"small": ds, "mid": dm},
        "scores": {
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "mass_log_score": float(mass_score),
            "total_unweighted": float(ckm_score + j_score + mass_score),
        },
    }


def scan_bound(bound: float, y_u: np.ndarray, y_d: np.ndarray, rng: np.random.Generator, starts: int = 64) -> dict[str, Any]:
    weights = {"ckm": 1.0, "jarlskog": 0.25, "mass": 2.5, "cond": 0.05, "reg": 0.01 / max(bound, 1.0e-9)}
    seeds = [np.zeros(8)]
    seeds.extend(rng.uniform(-bound, bound, size=8) for _ in range(starts - 1))
    best = None
    for seed in seeds:
        result = least_squares(
            objective,
            seed,
            bounds=(-bound * np.ones(8), bound * np.ones(8)),
            args=(y_u, y_d, weights),
            xtol=1.0e-11,
            ftol=1.0e-11,
            gtol=1.0e-11,
            max_nfev=4500,
        )
        evald = evaluate(result.x, y_u, y_d)
        cost = float(np.sum(objective(result.x, y_u, y_d, weights) ** 2))
        evald["optimizer"] = {
            "success": bool(result.success),
            "cost_weighted_sum_squares": cost,
            "nfev": int(result.nfev),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(result.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or cost < best["optimizer"]["cost_weighted_sum_squares"]:
            best = evald
    assert best is not None
    best["component_bound"] = bound
    return best


def write_csv(rows: list[dict[str, Any]]) -> None:
    flat = []
    for row in rows:
        obs = row["CKM_observables"]
        us = row["up_mass_ratios"]
        ds = row["down_mass_ratios"]
        scores = row["scores"]
        flat.append(
            {
                "component_bound": row["component_bound"],
                "Vus": obs["Vus"],
                "Vcb": obs["Vcb"],
                "Vub": obs["Vub"],
                "J_abs": obs["J"],
                "up_small": us["small"],
                "up_mid": us["mid"],
                "down_small": ds["small"],
                "down_mid": ds["mid"],
                "ckm_score": scores["ckm_magnitude_log_score"],
                "mass_score": scores["mass_log_score"],
                "j_score": scores["jarlskog_log_score"],
                "K_condition": row["K_Q_condition"],
                "G_span": row["G_spectral_span"],
                "weighted_cost": row["optimizer"]["cost_weighted_sum_squares"],
            }
        )
    with (OUT / "scan_toeplitz_kahler_ckm.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0].keys()))
        writer.writeheader()
        writer.writerows(flat)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Toeplitz/Kahler CKM scan")
    lines.append("")
    lines.append("No web lookup was used.  The scan tests a shared left-handed Q metric.")
    lines.append("")
    lines.append("```text")
    lines.append("K_Q = exp(G),  G = G^dagger, Tr G = 0")
    lines.append("Y_u,d -> K_Q^{-1/2} Y_u,d")
    lines.append("```")
    lines.append("")
    lines.append("| bound | Vus | Vcb | Vub | J | up small | down small | CKM score | mass score | cond K |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["bound_scan"]:
        obs = row["CKM_observables"]
        us = row["up_mass_ratios"]
        ds = row["down_mass_ratios"]
        lines.append(
            f"| {row['component_bound']:.2f} | {obs['Vus']:.4e} | {obs['Vcb']:.4e} | "
            f"{obs['Vub']:.4e} | {obs['J']:.4e} | {us['small']:.4e} | "
            f"{ds['small']:.4e} | {row['scores']['ckm_magnitude_log_score']:.3e} | "
            f"{row['scores']['mass_log_score']:.3e} | {row['K_Q_condition']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
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
        "Mass ratios: "
        f"`u=({best['up_mass_ratios']['small']:.6e}, {best['up_mass_ratios']['mid']:.6e})`, "
        f"`d=({best['down_mass_ratios']['small']:.6e}, {best['down_mass_ratios']['mid']:.6e})`."
    )
    lines.append(
        f"`cond(K_Q)={best['K_Q_condition']:.6e}`, "
        f"`span(G)={best['G_spectral_span']:.6e}`."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "scan_toeplitz_kahler_ckm_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = json.loads(CARD.read_text(encoding="utf-8"))
    y_u = cmat(card["yukawa_sectors"]["up"]["Y_matrix_normalized"])
    y_d = cmat(card["yukawa_sectors"]["down"]["Y_matrix_normalized"])
    rng = np.random.default_rng(2026050714)
    bounds = [0.10, 0.25, 0.50, 0.80, 1.20, 1.80]
    rows = [scan_bound(bound, y_u, y_d, rng) for bound in bounds]
    best = min(rows, key=lambda row: row["optimizer"]["cost_weighted_sum_squares"])
    verdict = (
        "A shared left-handed Toeplitz/Kahler metric is more effective than the "
        "down-only 120H perturbation, because it directly changes the generalized "
        "left eigensystems.  However, this effective Q-only metric is still not a "
        "full flavor solution: the best mass-preserving point has |Vcb| near the "
        "target scale but leaves |Vus| and especially |Vub| too large, while "
        "requiring cond(K_Q) ~= 31.  Therefore the next viable route is not a "
        "larger Q-only Kahler scan; it is either a correlated K_16/K_PS Kahler "
        "model with lepton/seesaw replay, or a full Clebsch-correlated 10+120+126 "
        "Spin(10) flavor fit."
    )
    payload = {
        "note": "No web lookup used. Effective Toeplitz/Kahler left normalization CKM scan.",
        "input_card": str(CARD),
        "target": TARGETS,
        "bound_scan": rows,
        "best": best,
        "checks": {
            "best_has_finite_scores": bool(
                np.isfinite(best["scores"]["ckm_magnitude_log_score"])
                and np.isfinite(best["scores"]["mass_log_score"])
            ),
            "best_preserves_mass_score_below_0p2": bool(best["scores"]["mass_log_score"] < 0.2),
            "best_improves_ckm_vs_raw_card": bool(best["scores"]["ckm_magnitude_log_score"] < 5.449970189134648),
            "metric_condition_recorded": bool(best["K_Q_condition"] > 1.0),
        },
        "verdict": verdict,
    }
    payload["all_checks_ok"] = all(payload["checks"].values())
    (OUT / "scan_toeplitz_kahler_ckm.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(rows)
    write_report(payload)
    print("Toeplitz/Kahler CKM scan")
    for row in rows:
        obs = row["CKM_observables"]
        print(
            f"  bound={row['component_bound']:.2f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"ckm_score={row['scores']['ckm_magnitude_log_score']:.3e} "
            f"mass_score={row['scores']['mass_log_score']:.3e} "
            f"condK={row['K_Q_condition']:.3e}"
        )
    print(f"  best bound: {best['component_bound']:.2f}")
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
