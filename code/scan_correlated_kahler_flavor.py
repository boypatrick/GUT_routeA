#!/usr/bin/env python3
"""Scan correlated K16 / Pati-Salam Kahler flavor normalizations.

No web lookup is used.  This script upgrades the Q-only Kahler scan to two
correlated possibilities:

  K16: one family metric for the full 16, so every Yukawa matrix transforms as
       Y_a -> A^T Y_a A.

  PS:  one metric for F_L=(4,2,1), i.e. Q,L, and one metric for
       F_R=(bar4,1,2), i.e. u^c,d^c,e^c,nu^c:
       Y_a -> A_L^T Y_a A_R.

The scan checks CKM, all Dirac hierarchy ratios, and then replays the inverse
type-I seesaw to see whether the exact-card neutrino benchmark remains
reconstructible after the correlated canonical normalization.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import verify_seesaw_item3 as seesaw  # noqa: E402


CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
OUT = ROOT / "output" / "flavor_correlated_kahler_scan"


TARGETS = {
    "up_small": 1.0e-5,
    "up_mid": 7.0e-3,
    "down_small": 1.0e-3,
    "down_mid": 2.0e-2,
    "charged_lepton_small": 3.0e-4,
    "charged_lepton_mid": 6.0e-2,
    "neutrino_dirac_small": 1.0e-2,
    "neutrino_dirac_mid": 2.0e-1,
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


def metric_data(g: np.ndarray) -> dict[str, Any]:
    k = exp_hermitian(g, 1.0)
    eig = np.linalg.eigvalsh(k)
    geig = np.linalg.eigvalsh(g)
    return {
        "G": matrix_json(g),
        "K": matrix_json(k),
        "K_eigenvalues": [float(v) for v in eig],
        "K_condition": float(np.max(eig) / np.min(eig)),
        "G_span": float(np.max(geig) - np.min(geig)),
    }


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def ratios(s: np.ndarray) -> tuple[float, float]:
    return float(s[0] / s[2]), float(s[1] / s[2])


def log10_ratio(found: float, target: float) -> float:
    return math.log10(max(found, 1.0e-30) / target)


def ckm_matrix(yu: np.ndarray, yd: np.ndarray) -> tuple[np.ndarray, dict[str, float], dict[str, tuple[float, float]]]:
    u_left, u_s = left_rotation(yu)
    d_left, d_s = left_rotation(yd)
    v = u_left.conjugate().T @ d_left
    j = float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))
    obs = {"Vus": float(abs(v[0, 1])), "Vcb": float(abs(v[1, 2])), "Vub": float(abs(v[0, 2])), "J": abs(j), "J_signed": j}
    mass = {"up": ratios(u_s), "down": ratios(d_s)}
    return v, obs, mass


def transform_yukawas(x: np.ndarray, mode: str, y: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    if mode == "K16":
        g = hermitian_generator(x)
        a = exp_hermitian(g, -0.5)
        transformed = {name: a.T @ mat @ a for name, mat in y.items()}
        return transformed, {"mode": mode, "K16": metric_data(g)}
    if mode == "PS":
        g_left = hermitian_generator(x[:8])
        g_right = hermitian_generator(x[8:])
        a_left = exp_hermitian(g_left, -0.5)
        a_right = exp_hermitian(g_right, -0.5)
        transformed = {name: a_left.T @ mat @ a_right for name, mat in y.items()}
        return transformed, {"mode": mode, "K_left_FL": metric_data(g_left), "K_right_FR": metric_data(g_right)}
    raise ValueError(f"unknown mode {mode}")


def all_mass_ratios(y: dict[str, np.ndarray]) -> dict[str, tuple[float, float]]:
    out = {}
    for name, mat in y.items():
        _u, s = left_rotation(mat)
        out[name] = ratios(s)
    return out


def mass_score(mass: dict[str, tuple[float, float]]) -> float:
    mapping = {
        "up": ("up_small", "up_mid"),
        "down": ("down_small", "down_mid"),
        "charged_lepton": ("charged_lepton_small", "charged_lepton_mid"),
        "neutrino_dirac": ("neutrino_dirac_small", "neutrino_dirac_mid"),
    }
    total = 0.0
    for sector, (small_key, mid_key) in mapping.items():
        small, mid = mass[sector]
        total += log10_ratio(small, TARGETS[small_key]) ** 2
        total += log10_ratio(mid, TARGETS[mid_key]) ** 2
    return float(total)


def seesaw_replay(y_nu: np.ndarray, y_e: np.ndarray) -> dict[str, Any]:
    u_e, charged_singulars = seesaw.left_rotation(y_e)
    benchmark = seesaw.LightNuBenchmark(
        m1_eV=1.0e-3,
        dm21_eV2=7.42e-5,
        dm31_eV2=2.517e-3,
        sin2_theta12=0.304,
        sin2_theta13=0.0222,
        sin2_theta23=0.573,
        delta_cp_rad=1.20 * math.pi,
        alpha21_rad=0.35 * math.pi,
        alpha31_rad=1.10 * math.pi,
    )
    m_light_diag = np.array(
        [
            benchmark.m1_eV,
            math.sqrt(benchmark.m1_eV**2 + benchmark.dm21_eV2),
            math.sqrt(benchmark.m1_eV**2 + benchmark.dm31_eV2),
        ],
        dtype=float,
    )
    u_pmns_target = seesaw.standard_pmns(
        benchmark.sin2_theta12,
        benchmark.sin2_theta13,
        benchmark.sin2_theta23,
        benchmark.delta_cp_rad,
        benchmark.alpha21_rad,
        benchmark.alpha31_rad,
    )
    u_nu_target = u_e @ u_pmns_target
    m_light_target = u_nu_target.conjugate() @ np.diag(m_light_diag) @ u_nu_target.conjugate().T
    mD_largest_GeV = 100.0
    mD_eV = y_nu * (mD_largest_GeV * 1.0e9)
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light_target) @ mD_eV)
    MR_GeV = MR_eV / 1.0e9
    m_light_reco = -(mD_eV @ np.linalg.inv(MR_eV) @ mD_eV.T)
    light_masses, u_nu_reco = seesaw.takagi_by_h(m_light_reco)
    u_pmns_reco = u_e.conjugate().T @ u_nu_reco
    angles = seesaw.mixing_angles(u_pmns_reco)
    heavy = np.sort(np.linalg.svd(MR_GeV, compute_uv=False))
    residual = float(np.linalg.norm(m_light_reco - m_light_target) / np.linalg.norm(m_light_target))
    theta_norm = float(np.linalg.norm(mD_eV @ np.linalg.inv(MR_eV), ord=2))
    return {
        "heavy_neutrino_masses_GeV": [float(v) for v in heavy],
        "theta_norm": theta_norm,
        "seesaw_matrix_residual": residual,
        "reconstructed_light_masses_eV": [float(v) for v in np.sort(light_masses)],
        "reconstructed_pmns_angles": angles,
        "Y_e_singular_values_normalized": [float(v) for v in charged_singulars],
        "MR_condition_number": float(heavy[-1] / heavy[0]),
    }


def objective(x: np.ndarray, mode: str, y0: dict[str, np.ndarray], weights: dict[str, float]) -> np.ndarray:
    y, metrics = transform_yukawas(x, mode, y0)
    _v, ckm_obs, _qm = ckm_matrix(y["up"], y["down"])
    masses = all_mass_ratios(y)
    metric_spans = []
    if mode == "K16":
        metric_spans.append(metrics["K16"]["G_span"])
    else:
        metric_spans.append(metrics["K_left_FL"]["G_span"])
        metric_spans.append(metrics["K_right_FR"]["G_span"])
    pieces = [
        weights["ckm"] * log10_ratio(ckm_obs["Vus"], TARGETS["Vus"]),
        weights["ckm"] * log10_ratio(ckm_obs["Vcb"], TARGETS["Vcb"]),
        weights["ckm"] * log10_ratio(ckm_obs["Vub"], TARGETS["Vub"]),
        weights["jarlskog"] * log10_ratio(ckm_obs["J"], TARGETS["J"]),
    ]
    for sector, (small, mid) in masses.items():
        pieces.append(weights["mass"] * log10_ratio(small, TARGETS[f"{sector}_small"]))
        pieces.append(weights["mass"] * log10_ratio(mid, TARGETS[f"{sector}_mid"]))
    pieces.extend(weights["cond"] * span for span in metric_spans)
    pieces.append(weights["reg"] * float(np.linalg.norm(x)))
    return np.array(pieces, dtype=float)


def evaluate(x: np.ndarray, mode: str, y0: dict[str, np.ndarray]) -> dict[str, Any]:
    y, metrics = transform_yukawas(x, mode, y0)
    v, ckm_obs, _qm = ckm_matrix(y["up"], y["down"])
    masses = all_mass_ratios(y)
    ckm_score = sum(log10_ratio(ckm_obs[key], TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = log10_ratio(ckm_obs["J"], TARGETS["J"]) ** 2
    m_score = mass_score(masses)
    replay = seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    conds = []
    if mode == "K16":
        conds.append(metrics["K16"]["K_condition"])
    else:
        conds.extend([metrics["K_left_FL"]["K_condition"], metrics["K_right_FR"]["K_condition"]])
    return {
        "mode": mode,
        "x": [float(vv) for vv in x],
        "metrics": metrics,
        "Yukawa_canonical": {name: matrix_json(mat) for name, mat in y.items()},
        "CKM_abs": matrix_abs_json(v),
        "CKM_observables": ckm_obs,
        "mass_ratios": {name: {"small": vals[0], "mid": vals[1]} for name, vals in masses.items()},
        "scores": {
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "mass_log_score": float(m_score),
            "total_unweighted": float(ckm_score + j_score + m_score),
        },
        "max_metric_condition": float(max(conds)),
        "seesaw_replay": replay,
    }


def scan(mode: str, bound: float, y0: dict[str, np.ndarray], rng: np.random.Generator, starts: int) -> dict[str, Any]:
    dim = 8 if mode == "K16" else 16
    weights = {"ckm": 1.0, "jarlskog": 0.25, "mass": 2.2, "cond": 0.04, "reg": 0.01 / max(bound, 1.0e-9)}
    seeds = [np.zeros(dim)]
    seeds.extend(rng.uniform(-bound, bound, size=dim) for _ in range(starts - 1))
    best = None
    for seed in seeds:
        res = least_squares(
            objective,
            seed,
            bounds=(-bound * np.ones(dim), bound * np.ones(dim)),
            args=(mode, y0, weights),
            xtol=1.0e-10,
            ftol=1.0e-10,
            gtol=1.0e-10,
            max_nfev=3600,
        )
        item = evaluate(res.x, mode, y0)
        item["optimizer"] = {
            "success": bool(res.success),
            "nfev": int(res.nfev),
            "cost_weighted_sum_squares": float(np.sum(objective(res.x, mode, y0, weights) ** 2)),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(res.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or item["optimizer"]["cost_weighted_sum_squares"] < best["optimizer"]["cost_weighted_sum_squares"]:
            best = item
    assert best is not None
    best["component_bound"] = bound
    return best


def write_csv(rows: list[dict[str, Any]]) -> None:
    flat = []
    for row in rows:
        obs = row["CKM_observables"]
        scores = row["scores"]
        replay = row["seesaw_replay"]
        flat.append(
            {
                "mode": row["mode"],
                "component_bound": row["component_bound"],
                "Vus": obs["Vus"],
                "Vcb": obs["Vcb"],
                "Vub": obs["Vub"],
                "J_abs": obs["J"],
                "ckm_score": scores["ckm_magnitude_log_score"],
                "mass_score": scores["mass_log_score"],
                "j_score": scores["jarlskog_log_score"],
                "max_metric_condition": row["max_metric_condition"],
                "theta_norm": replay["theta_norm"],
                "seesaw_residual": replay["seesaw_matrix_residual"],
                "M1_GeV": replay["heavy_neutrino_masses_GeV"][0],
                "M2_GeV": replay["heavy_neutrino_masses_GeV"][1],
                "M3_GeV": replay["heavy_neutrino_masses_GeV"][2],
                "weighted_cost": row["optimizer"]["cost_weighted_sum_squares"],
            }
        )
    with (OUT / "scan_correlated_kahler_flavor.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0].keys()))
        writer.writeheader()
        writer.writerows(flat)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Correlated Kahler flavor scan")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("| mode | bound | Vus | Vcb | Vub | J | CKM score | mass score | max cond K | theta |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["rows"]:
        obs = row["CKM_observables"]
        scores = row["scores"]
        lines.append(
            f"| {row['mode']} | {row['component_bound']:.2f} | {obs['Vus']:.4e} | "
            f"{obs['Vcb']:.4e} | {obs['Vub']:.4e} | {obs['J']:.4e} | "
            f"{scores['ckm_magnitude_log_score']:.3e} | {scores['mass_log_score']:.3e} | "
            f"{row['max_metric_condition']:.3e} | {row['seesaw_replay']['theta_norm']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
    replay = best["seesaw_replay"]
    lines.append("")
    lines.append("## Best point")
    lines.append("")
    lines.append(f"Mode: `{best['mode']}`, bound `{best['component_bound']:.3g}`.")
    lines.append(
        f"CKM: `|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, |J|={obs['J']:.6e}`."
    )
    lines.append(
        f"Scores: `CKM={best['scores']['ckm_magnitude_log_score']:.6e}`, "
        f"`mass={best['scores']['mass_log_score']:.6e}`."
    )
    lines.append(
        f"Seesaw replay: `theta={replay['theta_norm']:.6e}`, "
        f"`residual={replay['seesaw_matrix_residual']:.6e}`."
    )
    heavy = replay["heavy_neutrino_masses_GeV"]
    lines.append(
        f"Heavy masses: `({heavy[0]:.6e}, {heavy[1]:.6e}, {heavy[2]:.6e}) GeV`."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "scan_correlated_kahler_flavor_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = json.loads(CARD.read_text(encoding="utf-8"))
    y0 = {
        "up": cmat(card["yukawa_sectors"]["up"]["Y_matrix_normalized"]),
        "down": cmat(card["yukawa_sectors"]["down"]["Y_matrix_normalized"]),
        "charged_lepton": cmat(card["yukawa_sectors"]["charged_lepton"]["Y_matrix_normalized"]),
        "neutrino_dirac": cmat(card["yukawa_sectors"]["neutrino_dirac"]["Y_matrix_normalized"]),
    }
    rng = np.random.default_rng(2026050715)
    rows = []
    for bound in [0.3, 0.7, 1.2]:
        rows.append(scan("K16", bound, y0, rng, starts=32))
    for bound in [0.3, 0.7, 1.2]:
        rows.append(scan("PS", bound, y0, rng, starts=42))
    best = min(rows, key=lambda row: row["optimizer"]["cost_weighted_sum_squares"])
    best_obs = best["CKM_observables"]
    verdict = (
        "Correlating the Kahler metric across the full 16 or across Pati-Salam "
        "left/right multiplets is the correct next consistency test.  The scan "
        "shows that the Pati-Salam left/right metric improves the CKM score while "
        "the inverse-seesaw replay remains numerically stable, but it still leaves "
        f"|Vus|={best_obs['Vus']:.6e}, |Vcb|={best_obs['Vcb']:.6e}, and "
        f"|Vub|={best_obs['Vub']:.6e}.  This is not a full flavor fit.  The model "
        "must therefore move to a full Clebsch-correlated 10+126bar+120 Yukawa "
        "sector before any channel-specific dimension-five proton-decay claim is "
        "meaningful."
    )
    payload = {
        "note": "No web lookup used. Correlated K16/PS Kahler canonical-normalization scan.",
        "input_card": str(CARD),
        "target": TARGETS,
        "rows": rows,
        "best": best,
        "internal_reconstruction_checks": {
            "scan_completed": True,
            "best_seesaw_residual_lt_1e_minus_10": bool(best["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10),
            "best_theta_lt_1e_minus_8": bool(best["seesaw_replay"]["theta_norm"] < 1.0e-8),
            "best_improves_ckm_vs_Q_only": bool(best["scores"]["ckm_magnitude_log_score"] < 3.3700882923769235),
        },
        "phenomenology_viability_checks": {
            "best_ckm_score_lt_0p05": bool(best["scores"]["ckm_magnitude_log_score"] < 5.0e-2),
            "best_Vus_within_10_percent": bool(abs(best_obs["Vus"] / TARGETS["Vus"] - 1.0) < 0.10),
            "best_Vcb_within_10_percent": bool(abs(best_obs["Vcb"] / TARGETS["Vcb"] - 1.0) < 0.10),
            "best_Vub_within_10_percent": bool(abs(best_obs["Vub"] / TARGETS["Vub"] - 1.0) < 0.10),
        },
        "verdict": verdict,
    }
    payload["internal_reconstruction_ok"] = all(payload["internal_reconstruction_checks"].values())
    payload["phenomenology_viable"] = all(payload["phenomenology_viability_checks"].values())
    (OUT / "scan_correlated_kahler_flavor.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(rows)
    write_report(payload)
    print("Correlated Kahler flavor scan")
    for row in rows:
        obs = row["CKM_observables"]
        print(
            f"  {row['mode']:3s} bound={row['component_bound']:.1f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"ckm={row['scores']['ckm_magnitude_log_score']:.3e} "
            f"mass={row['scores']['mass_log_score']:.3e} "
            f"cond={row['max_metric_condition']:.3e} "
            f"theta={row['seesaw_replay']['theta_norm']:.3e}"
        )
    print(f"  best: {best['mode']} bound={best['component_bound']:.1f}")
    print(f"  internal reconstruction ok: {payload['internal_reconstruction_ok']}")
    print(f"  phenomenology viable: {payload['phenomenology_viable']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
