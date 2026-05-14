#!/usr/bin/env python3
"""Two-120-direction Clebsch flavor scan.

No web lookup is used.  This is the representation-controlled enlargement
after the single-120_H homotopy scan: keep H,F in the CP1/O(-4) Veronese
symmetric subspace, but allow two independent antisymmetric 120_H family
directions G_A,G_B.

The ansatz is

  Y_u  = H + r_u F + a_u G_A + b_u G_B,
  Y_nu = H - 3 r_u F + a_nu G_A + b_nu G_B,
  Y_d  = r_d (H + F + a_d G_A + b_d G_B),
  Y_e  = r_d (H - 3 F + a_e G_A + b_e G_B).

This is not a fully arbitrary flavor fit: the symmetric tensors remain
geometric, and the new freedom is precisely one additional antisymmetric
Spin(10) 120-like direction.
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

import scan_clebsch_flavor_fit as fit  # noqa: E402


HOMOTOPY = ROOT / "output" / "flavor_clebsch_homotopy" / "scan_clebsch_flavor_homotopy.json"
OUT = ROOT / "output" / "flavor_clebsch_two120"


STAGES = [
    {"name": "ckm_lock", "ckm": 4.5, "jarlskog": 0.25, "mass": 1.0, "reg": 0.003},
    {"name": "mass_2", "ckm": 4.5, "jarlskog": 0.25, "mass": 2.0, "reg": 0.003},
    {"name": "mass_4", "ckm": 4.8, "jarlskog": 0.25, "mass": 4.0, "reg": 0.004},
    {"name": "mass_8", "ckm": 5.0, "jarlskog": 0.25, "mass": 8.0, "reg": 0.005},
    {"name": "mass_12", "ckm": 5.5, "jarlskog": 0.25, "mass": 12.0, "reg": 0.006},
]


def load_basis() -> np.ndarray:
    card = json.loads(fit.CARD.read_text(encoding="utf-8"))
    return np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )


def unpack_complex(x: np.ndarray, start: int, count: int) -> tuple[np.ndarray, int]:
    vals = x[start : start + 2 * count]
    return vals[:count] + 1j * vals[count:], start + 2 * count


def pack_complex(v: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(v), np.imag(v)]).astype(float)


def decode_single_to_two120(x_single: np.ndarray) -> np.ndarray:
    pos = 0
    h_coeff, pos = unpack_complex(x_single, pos, 5)
    f_coeff, pos = unpack_complex(x_single, pos, 5)
    g_coeff, pos = unpack_complex(x_single, pos, 3)
    mix, pos = unpack_complex(x_single, pos, 6)
    assert pos == len(x_single)
    r_u, r_d, g_u, g_d, g_e, g_nu = mix
    g2 = np.array([1.0e-3, -1.0e-3j, 1.0e-3 + 1.0e-3j], dtype=complex)
    mix2 = np.array(
        [
            r_u,
            r_d,
            g_u,
            0.0 + 0.0j,
            g_d,
            0.0 + 0.0j,
            g_e,
            0.0 + 0.0j,
            g_nu,
            0.0 + 0.0j,
        ],
        dtype=complex,
    )
    return np.concatenate(
        [
            pack_complex(h_coeff),
            pack_complex(f_coeff),
            pack_complex(g_coeff),
            pack_complex(g2),
            pack_complex(mix2),
        ]
    )


def load_seed() -> np.ndarray:
    payload = json.loads(HOMOTOPY.read_text(encoding="utf-8"))
    return decode_single_to_two120(np.array(payload["best_pareto"]["x"], dtype=float))


def decode(x: np.ndarray, basis: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    pos = 0
    h_coeff, pos = unpack_complex(x, pos, 5)
    f_coeff, pos = unpack_complex(x, pos, 5)
    ga_coeff, pos = unpack_complex(x, pos, 3)
    gb_coeff, pos = unpack_complex(x, pos, 3)
    mix, pos = unpack_complex(x, pos, 10)
    assert pos == len(x)

    h = fit.symmetric_from_coeff(h_coeff, basis)
    f = fit.symmetric_from_coeff(f_coeff, basis)
    ga = fit.antisym_from_coeff(ga_coeff)
    gb = fit.antisym_from_coeff(gb_coeff)
    r_u, r_d, a_u, b_u, a_d, b_d, a_e, b_e, a_nu, b_nu = mix
    y = {
        "up": h + r_u * f + a_u * ga + b_u * gb,
        "down": r_d * (h + f + a_d * ga + b_d * gb),
        "charged_lepton": r_d * (h - 3.0 * f + a_e * ga + b_e * gb),
        "neutrino_dirac": h - 3.0 * r_u * f + a_nu * ga + b_nu * gb,
    }
    data = {
        "H_coefficients": fit.complex_list_json(h_coeff),
        "F_coefficients": fit.complex_list_json(f_coeff),
        "G_A_coefficients": fit.complex_list_json(ga_coeff),
        "G_B_coefficients": fit.complex_list_json(gb_coeff),
        "mixing_coefficients": {
            "r_u": fit.cjson(r_u),
            "r_d": fit.cjson(r_d),
            "a_u": fit.cjson(a_u),
            "b_u": fit.cjson(b_u),
            "a_d": fit.cjson(a_d),
            "b_d": fit.cjson(b_d),
            "a_e": fit.cjson(a_e),
            "b_e": fit.cjson(b_e),
            "a_nu": fit.cjson(a_nu),
            "b_nu": fit.cjson(b_nu),
        },
        "norms": {
            "H_fro": float(np.linalg.norm(h)),
            "F_fro": float(np.linalg.norm(f)),
            "G_A_fro": float(np.linalg.norm(ga)),
            "G_B_fro": float(np.linalg.norm(gb)),
            "parameter_l2": float(np.linalg.norm(x)),
        },
    }
    return y, data


def objective(x: np.ndarray, basis: np.ndarray, weights: dict[str, float]) -> np.ndarray:
    y, _data = decode(x, basis)
    _v, ckm_obs, _qm = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    pieces = [
        weights["ckm"] * fit.log10_ratio(ckm_obs["Vus"], fit.TARGETS["Vus"]),
        weights["ckm"] * fit.log10_ratio(ckm_obs["Vcb"], fit.TARGETS["Vcb"]),
        weights["ckm"] * fit.log10_ratio(ckm_obs["Vub"], fit.TARGETS["Vub"]),
        weights["jarlskog"] * fit.log10_ratio(ckm_obs["J"], fit.TARGETS["J"]),
    ]
    for sector, (small, mid) in masses.items():
        pieces.append(weights["mass"] * fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]))
        pieces.append(weights["mass"] * fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]))
    pieces.append(weights["reg"] * float(np.linalg.norm(x)))
    return np.array(pieces, dtype=float)


def evaluate(x: np.ndarray, basis: np.ndarray, stage: str, bound: float) -> dict[str, Any]:
    y, data = decode(x, basis)
    v, ckm_obs, _qm = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm_obs[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = fit.log10_ratio(ckm_obs["J"], fit.TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    replay = fit.seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    return {
        "stage": stage,
        "component_bound": bound,
        "x": [float(vv) for vv in x],
        "model_data": data,
        "Yukawa_fit": {name: fit.matrix_json(mat) for name, mat in y.items()},
        "CKM_abs": fit.matrix_abs_json(v),
        "CKM_observables": ckm_obs,
        "mass_ratios": {name: {"small": vals[0], "mid": vals[1]} for name, vals in masses.items()},
        "scores": {
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "mass_log_score": float(mass_score),
            "total_unweighted": float(ckm_score + j_score + mass_score),
        },
        "seesaw_replay": replay,
    }


def run_stage(
    stage: dict[str, float | str],
    seed: np.ndarray,
    basis: np.ndarray,
    rng: np.random.Generator,
    bound: float,
) -> dict[str, Any]:
    weights = {
        "ckm": float(stage["ckm"]),
        "jarlskog": float(stage["jarlskog"]),
        "mass": float(stage["mass"]),
        "reg": float(stage["reg"]) / bound,
    }
    candidates = [np.clip(seed, -0.95 * bound, 0.95 * bound)]
    for scale in [0.025, 0.060]:
        candidates.append(np.clip(seed + rng.normal(0.0, scale, size=seed.shape), -0.95 * bound, 0.95 * bound))

    best = None
    for candidate in candidates:
        res = least_squares(
            objective,
            candidate,
            bounds=(-bound * np.ones_like(seed), bound * np.ones_like(seed)),
            args=(basis, weights),
            xtol=3.0e-10,
            ftol=3.0e-10,
            gtol=3.0e-10,
            max_nfev=900,
        )
        item = evaluate(res.x, basis, str(stage["name"]), bound)
        item["stage_weights"] = weights
        item["optimizer"] = {
            "success": bool(res.success),
            "nfev": int(res.nfev),
            "weighted_cost": float(np.sum(objective(res.x, basis, weights) ** 2)),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(res.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or item["optimizer"]["weighted_cost"] < best["optimizer"]["weighted_cost"]:
            best = item
    assert best is not None
    return best


def flat_row(item: dict[str, Any]) -> dict[str, Any]:
    obs = item["CKM_observables"]
    scores = item["scores"]
    replay = item["seesaw_replay"]
    return {
        "stage": item["stage"],
        "mass_weight": item["stage_weights"]["mass"],
        "ckm_weight": item["stage_weights"]["ckm"],
        "Vus": obs["Vus"],
        "Vcb": obs["Vcb"],
        "Vub": obs["Vub"],
        "J_abs": obs["J"],
        "ckm_score": scores["ckm_magnitude_log_score"],
        "mass_score": scores["mass_log_score"],
        "j_score": scores["jarlskog_log_score"],
        "theta_norm": replay["theta_norm"],
        "seesaw_residual": replay["seesaw_matrix_residual"],
        "M1_GeV": replay["heavy_neutrino_masses_GeV"][0],
        "M2_GeV": replay["heavy_neutrino_masses_GeV"][1],
        "M3_GeV": replay["heavy_neutrino_masses_GeV"][2],
        "weighted_cost": item["optimizer"]["weighted_cost"],
    }


def write_csv(items: list[dict[str, Any]]) -> None:
    rows = [flat_row(item) for item in items]
    with (OUT / "scan_clebsch_flavor_two120.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Two-120-direction Clebsch flavor scan")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Ansatz")
    lines.append("")
    lines.append("```text")
    lines.append("Y_u  = H + r_u F + a_u G_A + b_u G_B")
    lines.append("Y_nu = H - 3 r_u F + a_nu G_A + b_nu G_B")
    lines.append("Y_d  = r_d (H + F + a_d G_A + b_d G_B)")
    lines.append("Y_e  = r_d (H - 3 F + a_e G_A + b_e G_B)")
    lines.append("H,F in the CP1 O(-4) Veronese symmetric subspace")
    lines.append("G_A,G_B in Lambda^2 C^3")
    lines.append("```")
    lines.append("")
    lines.append("| stage | mass weight | Vus | Vcb | Vub | J | CKM score | mass score | theta |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for item in payload["stages"]:
        row = flat_row(item)
        lines.append(
            f"| {row['stage']} | {row['mass_weight']:.2f} | {row['Vus']:.4e} | "
            f"{row['Vcb']:.4e} | {row['Vub']:.4e} | {row['J_abs']:.4e} | "
            f"{row['ckm_score']:.3e} | {row['mass_score']:.3e} | {row['theta_norm']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
    scores = best["scores"]
    replay = best["seesaw_replay"]
    lines.append("")
    lines.append("## Best point")
    lines.append("")
    lines.append(f"Stage: `{best['stage']}`.")
    lines.append(
        f"CKM: `|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, |J|={obs['J']:.6e}`."
    )
    lines.append(
        f"Scores: `CKM={scores['ckm_magnitude_log_score']:.6e}`, "
        f"`mass={scores['mass_log_score']:.6e}`."
    )
    heavy = replay["heavy_neutrino_masses_GeV"]
    lines.append(
        f"Seesaw: `theta={replay['theta_norm']:.6e}`, "
        f"`residual={replay['seesaw_matrix_residual']:.6e}`, "
        f"`M=({heavy[0]:.6e}, {heavy[1]:.6e}, {heavy[2]:.6e}) GeV`."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "scan_clebsch_flavor_two120_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    basis = load_basis()
    seed = load_seed()
    rng = np.random.default_rng(2026050721)
    bound = 3.2
    stages = []
    current = seed
    print("Two-120-direction Clebsch flavor scan", flush=True)
    for stage in STAGES:
        item = run_stage(stage, current, basis, rng, bound)
        stages.append(item)
        current = np.array(item["x"], dtype=float)
        obs = item["CKM_observables"]
        scores = item["scores"]
        print(
            f"  {stage['name']:8s} mass_w={stage['mass']:.1f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"ckm={scores['ckm_magnitude_log_score']:.3e} "
            f"mass={scores['mass_log_score']:.3e}",
            flush=True,
        )

    viable = [
        item
        for item in stages
        if item["scores"]["ckm_magnitude_log_score"] < 5.0e-2
        and item["scores"]["mass_log_score"] < 2.0e-1
    ]
    if viable:
        best = min(viable, key=lambda item: item["scores"]["mass_log_score"])
    else:
        best = min(stages, key=lambda item: item["scores"]["mass_log_score"] + 4.0 * item["scores"]["ckm_magnitude_log_score"])

    obs = best["CKM_observables"]
    scores = best["scores"]
    verdict = (
        "The two-120-direction scan tests the minimal representation-controlled "
        "antisymmetric enlargement.  The best point has "
        f"|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, CKM score={scores['ckm_magnitude_log_score']:.6e}, "
        f"and mass score={scores['mass_log_score']:.6e}.  A joint viable point "
        "exists only if both CKM and mass scores pass the stored thresholds."
    )
    payload = {
        "note": "No web lookup used. Two independent antisymmetric 120-like Clebsch directions.",
        "input_homotopy": str(HOMOTOPY),
        "stages": stages,
        "best": best,
        "checks": {
            "scan_completed": True,
            "any_ckm_score_lt_0p05": any(item["scores"]["ckm_magnitude_log_score"] < 5.0e-2 for item in stages),
            "any_mass_score_lt_0p20": any(item["scores"]["mass_log_score"] < 2.0e-1 for item in stages),
            "any_joint_viable": bool(viable),
            "best_seesaw_residual_lt_1e_minus_10": bool(best["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10),
            "best_theta_lt_1e_minus_8": bool(best["seesaw_replay"]["theta_norm"] < 1.0e-8),
            "improves_mass_vs_single120_homotopy": bool(best["scores"]["mass_log_score"] < 0.7952350183081083),
        },
        "verdict": verdict,
    }
    payload["phenomenology_viable"] = payload["checks"]["any_joint_viable"]
    (OUT / "scan_clebsch_flavor_two120.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(stages)
    write_report(payload)
    print(f"  best: {best['stage']}", flush=True)
    print(f"  phenomenology viable: {payload['phenomenology_viable']}", flush=True)
    print(f"  wrote: {OUT}", flush=True)


if __name__ == "__main__":
    main()
