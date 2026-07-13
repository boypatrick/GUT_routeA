#!/usr/bin/env python3
"""Orthogonal-symmetric relaxation diagnostic for the Clebsch flavor scan.

No web lookup is used.  This is a deliberately controlled test of the flavor
blocker found by the two-120-direction audit.  The CP1/O(-4) symmetric product
map spans a five-complex-dimensional subspace of Sym^2(C^3).  Since
dim_C Sym^2(C^3)=6, there is exactly one normalized symmetric direction
S_perp orthogonal to the Veronese image.  We add only this direction,
separately to H or F,

  H -> H + eps_H S_perp,     or     F -> F + eps_F S_perp,

while keeping the two antisymmetric 120-like directions G_A,G_B.  The scan
measures whether a small non-Veronese symmetric leakage is enough to recover a
joint CKM/mass fit, and replays the inverse type-I seesaw.
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
import scan_clebsch_flavor_two120 as two120  # noqa: E402


TWO120 = ROOT / "output" / "flavor_clebsch_two120" / "scan_clebsch_flavor_two120.json"
OUT = ROOT / "output" / "flavor_symmetric_relaxation"

BASE_BOUND = 3.2
EPS_BUDGETS = [0.03, 0.10, 0.30, 1.00]
STAGES = [
    {"name": "mass_12", "ckm": 5.5, "jarlskog": 0.25, "mass": 12.0, "reg": 0.006},
    {"name": "mass_24", "ckm": 6.0, "jarlskog": 0.25, "mass": 24.0, "reg": 0.007},
    {"name": "mass_40", "ckm": 7.0, "jarlskog": 0.25, "mass": 40.0, "reg": 0.008},
]


def load_basis() -> np.ndarray:
    card = json.loads(fit.CARD.read_text(encoding="utf-8"))
    return np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )


def sym_basis() -> list[np.ndarray]:
    out = []
    for i in range(3):
        m = np.zeros((3, 3), dtype=complex)
        m[i, i] = 1.0
        out.append(m)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        m = np.zeros((3, 3), dtype=complex)
        m[i, j] = 1.0 / math.sqrt(2.0)
        m[j, i] = 1.0 / math.sqrt(2.0)
        out.append(m)
    return out


def sym_coordinates(m: np.ndarray, basis6: list[np.ndarray]) -> np.ndarray:
    return np.array([np.vdot(b, m) for b in basis6], dtype=complex)


def matrix_from_sym_coordinates(c: np.ndarray, basis6: list[np.ndarray]) -> np.ndarray:
    out = np.zeros((3, 3), dtype=complex)
    for coeff, b in zip(c, basis6):
        out += coeff * b
    return out


def orthogonal_symmetric_direction(veronese_basis: np.ndarray) -> dict[str, Any]:
    basis6 = sym_basis()
    coords = np.column_stack([sym_coordinates(veronese_basis[:, :, k], basis6) for k in range(5)])
    _u, singulars, vh = np.linalg.svd(coords.conjugate().T, full_matrices=True)
    null_vec = vh[-1].conjugate()
    null_vec = null_vec / np.linalg.norm(null_vec)
    s_perp = matrix_from_sym_coordinates(null_vec, basis6)
    overlaps = [np.vdot(veronese_basis[:, :, k], s_perp) for k in range(5)]
    symmetry_error = float(np.linalg.norm(s_perp - s_perp.T))
    fro_norm = float(np.linalg.norm(s_perp))
    return {
        "S_perp": s_perp,
        "coordinates": null_vec,
        "veronese_rank": int(np.linalg.matrix_rank(coords, tol=1.0e-12)),
        "singular_values": [float(v) for v in singulars],
        "orthogonality_max_abs": float(max(abs(z) for z in overlaps)),
        "symmetry_error": symmetry_error,
        "fro_norm": fro_norm,
    }


def unpack_complex(x: np.ndarray, start: int, count: int) -> tuple[np.ndarray, int]:
    vals = x[start : start + 2 * count]
    return vals[:count] + 1j * vals[count:], start + 2 * count


def pack_complex(v: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(v), np.imag(v)]).astype(float)


def load_seed() -> np.ndarray:
    payload = json.loads(TWO120.read_text(encoding="utf-8"))
    return np.array(payload["best"]["x"], dtype=float)


def decode(
    x: np.ndarray,
    veronese_basis: np.ndarray,
    s_perp: np.ndarray,
    branch: str,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    pos = 0
    h_coeff, pos = unpack_complex(x, pos, 5)
    f_coeff, pos = unpack_complex(x, pos, 5)
    ga_coeff, pos = unpack_complex(x, pos, 3)
    gb_coeff, pos = unpack_complex(x, pos, 3)
    mix, pos = unpack_complex(x, pos, 10)
    eps_arr, pos = unpack_complex(x, pos, 1)
    assert pos == len(x)
    eps = eps_arr[0]

    h_geo = fit.symmetric_from_coeff(h_coeff, veronese_basis)
    f_geo = fit.symmetric_from_coeff(f_coeff, veronese_basis)
    h = h_geo + (eps * s_perp if branch == "H" else 0.0)
    fmat = f_geo + (eps * s_perp if branch == "F" else 0.0)
    ga = fit.antisym_from_coeff(ga_coeff)
    gb = fit.antisym_from_coeff(gb_coeff)
    r_u, r_d, a_u, b_u, a_d, b_d, a_e, b_e, a_nu, b_nu = mix
    y = {
        "up": h + r_u * fmat + a_u * ga + b_u * gb,
        "down": r_d * (h + fmat + a_d * ga + b_d * gb),
        "charged_lepton": r_d * (h - 3.0 * fmat + a_e * ga + b_e * gb),
        "neutrino_dirac": h - 3.0 * r_u * fmat + a_nu * ga + b_nu * gb,
    }
    reference_norm = np.linalg.norm(h_geo if branch == "H" else f_geo)
    leak_norm = abs(eps) * np.linalg.norm(s_perp)
    data = {
        "branch": branch,
        "epsilon": fit.cjson(eps),
        "epsilon_abs": float(abs(eps)),
        "relative_leakage_norm": float(leak_norm / max(reference_norm, 1.0e-30)),
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
            "H_geo_fro": float(np.linalg.norm(h_geo)),
            "F_geo_fro": float(np.linalg.norm(f_geo)),
            "H_total_fro": float(np.linalg.norm(h)),
            "F_total_fro": float(np.linalg.norm(fmat)),
            "S_perp_fro": float(np.linalg.norm(s_perp)),
            "G_A_fro": float(np.linalg.norm(ga)),
            "G_B_fro": float(np.linalg.norm(gb)),
            "parameter_l2": float(np.linalg.norm(x)),
        },
    }
    return y, data


def objective(
    x: np.ndarray,
    veronese_basis: np.ndarray,
    s_perp: np.ndarray,
    branch: str,
    weights: dict[str, float],
) -> np.ndarray:
    y, data = decode(x, veronese_basis, s_perp, branch)
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
    pieces.append(weights["reg"] * float(np.linalg.norm(x[:-2])))
    pieces.append(weights["epsilon"] * data["epsilon_abs"])
    return np.array(pieces, dtype=float)


def evaluate(
    x: np.ndarray,
    veronese_basis: np.ndarray,
    s_perp: np.ndarray,
    branch: str,
    stage: str,
    eps_budget: float,
) -> dict[str, Any]:
    y, data = decode(x, veronese_basis, s_perp, branch)
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
        "branch": branch,
        "stage": stage,
        "epsilon_component_budget": eps_budget,
        "base_component_bound": BASE_BOUND,
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


def bounds_for(seed: np.ndarray, eps_budget: float) -> tuple[np.ndarray, np.ndarray]:
    lower = -BASE_BOUND * np.ones_like(seed)
    upper = BASE_BOUND * np.ones_like(seed)
    lower[-2:] = -eps_budget
    upper[-2:] = eps_budget
    return lower, upper


def run_stage(
    stage: dict[str, float | str],
    seed: np.ndarray,
    veronese_basis: np.ndarray,
    s_perp: np.ndarray,
    branch: str,
    eps_budget: float,
    rng: np.random.Generator,
) -> dict[str, Any]:
    lower, upper = bounds_for(seed, eps_budget)
    weights = {
        "ckm": float(stage["ckm"]),
        "jarlskog": float(stage["jarlskog"]),
        "mass": float(stage["mass"]),
        "reg": float(stage["reg"]) / BASE_BOUND,
        "epsilon": 0.012 / max(eps_budget, 1.0e-12),
    }
    clipped = np.clip(seed, lower * 0.95, upper * 0.95)
    candidates = [clipped]
    jitter = np.zeros_like(seed)
    jitter[:-2] = rng.normal(0.0, 0.025, size=len(seed) - 2)
    jitter[-2:] = rng.normal(0.0, 0.12 * eps_budget, size=2)
    candidates.append(np.clip(clipped + jitter, lower * 0.95, upper * 0.95))
    best = None
    for candidate in candidates:
        try:
            res = least_squares(
                objective,
                candidate,
                bounds=(lower, upper),
                args=(veronese_basis, s_perp, branch, weights),
                xtol=5.0e-10,
                ftol=5.0e-10,
                gtol=5.0e-10,
                max_nfev=900,
            )
            x_opt = res.x
            optimizer = {
                "success": bool(res.success),
                "nfev": int(res.nfev),
                "weighted_cost": float(np.sum(objective(res.x, veronese_basis, s_perp, branch, weights) ** 2)),
                "active_base_bound_fraction": float(np.mean(np.isclose(np.abs(res.x[:-2]), BASE_BOUND, rtol=0.0, atol=1.0e-5))),
                "active_epsilon_bound_fraction": float(np.mean(np.isclose(np.abs(res.x[-2:]), eps_budget, rtol=0.0, atol=1.0e-5))),
                "fallback_used": False,
            }
        except Exception as exc:  # scipy may occasionally fail an internal SVD.
            x_opt = candidate
            optimizer = {
                "success": False,
                "nfev": 0,
                "weighted_cost": float(np.sum(objective(candidate, veronese_basis, s_perp, branch, weights) ** 2)),
                "active_base_bound_fraction": float(np.mean(np.isclose(np.abs(candidate[:-2]), BASE_BOUND, rtol=0.0, atol=1.0e-5))),
                "active_epsilon_bound_fraction": float(np.mean(np.isclose(np.abs(candidate[-2:]), eps_budget, rtol=0.0, atol=1.0e-5))),
                "fallback_used": True,
                "fallback_reason": repr(exc),
            }
        item = evaluate(x_opt, veronese_basis, s_perp, branch, str(stage["name"]), eps_budget)
        item["stage_weights"] = weights
        item["optimizer"] = optimizer
        if best is None or item["optimizer"]["weighted_cost"] < best["optimizer"]["weighted_cost"]:
            best = item
    assert best is not None
    return best


def flat_row(item: dict[str, Any]) -> dict[str, Any]:
    obs = item["CKM_observables"]
    scores = item["scores"]
    replay = item["seesaw_replay"]
    data = item["model_data"]
    return {
        "branch": item["branch"],
        "eps_budget": item["epsilon_component_budget"],
        "stage": item["stage"],
        "mass_weight": item["stage_weights"]["mass"],
        "Vus": obs["Vus"],
        "Vcb": obs["Vcb"],
        "Vub": obs["Vub"],
        "J_abs": obs["J"],
        "ckm_score": scores["ckm_magnitude_log_score"],
        "mass_score": scores["mass_log_score"],
        "j_score": scores["jarlskog_log_score"],
        "epsilon_abs": data["epsilon_abs"],
        "relative_leakage_norm": data["relative_leakage_norm"],
        "theta_norm": replay["theta_norm"],
        "seesaw_residual": replay["seesaw_matrix_residual"],
        "M1_GeV": replay["heavy_neutrino_masses_GeV"][0],
        "M2_GeV": replay["heavy_neutrino_masses_GeV"][1],
        "M3_GeV": replay["heavy_neutrino_masses_GeV"][2],
        "weighted_cost": item["optimizer"]["weighted_cost"],
    }


def write_csv(items: list[dict[str, Any]]) -> None:
    rows = [flat_row(item) for item in items]
    with (OUT / "scan_clebsch_flavor_symmetric_relaxation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Orthogonal-symmetric relaxation diagnostic")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Geometric split")
    lines.append("")
    lines.append("```text")
    lines.append("Sym^2(C^3) = Veronese_5 \\oplus C S_perp")
    lines.append("H -> H + eps_H S_perp, or F -> F + eps_F S_perp")
    lines.append("G_A,G_B remain in Lambda^2 C^3")
    lines.append("```")
    lines.append("")
    diag = payload["S_perp_diagnostics"]
    lines.append(
        f"`rank(Veronese)={diag['veronese_rank']}`, "
        f"`max |<C_m,S_perp>|={diag['orthogonality_max_abs']:.3e}`, "
        f"`||S_perp-S_perp^T||={diag['symmetry_error']:.3e}`."
    )
    lines.append("")
    lines.append("| branch | eps budget | stage | Vus | Vcb | Vub | CKM score | mass score | |eps| | rel leak | theta |")
    lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for item in payload["rows"]:
        row = flat_row(item)
        lines.append(
            f"| {row['branch']} | {row['eps_budget']:.2f} | {row['stage']} | "
            f"{row['Vus']:.4e} | {row['Vcb']:.4e} | {row['Vub']:.4e} | "
            f"{row['ckm_score']:.3e} | {row['mass_score']:.3e} | "
            f"{row['epsilon_abs']:.3e} | {row['relative_leakage_norm']:.3e} | "
            f"{row['theta_norm']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
    scores = best["scores"]
    replay = best["seesaw_replay"]
    data = best["model_data"]
    heavy = replay["heavy_neutrino_masses_GeV"]
    lines.append("")
    lines.append("## Best point")
    lines.append("")
    lines.append(f"Branch: `{best['branch']}`, stage `{best['stage']}`, eps budget `{best['epsilon_component_budget']:.3g}`.")
    lines.append(
        f"CKM: `|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, |J|={obs['J']:.6e}`."
    )
    lines.append(
        f"Scores: `CKM={scores['ckm_magnitude_log_score']:.6e}`, "
        f"`mass={scores['mass_log_score']:.6e}`."
    )
    lines.append(
        f"Leakage: `|eps|={data['epsilon_abs']:.6e}`, "
        f"`relative={data['relative_leakage_norm']:.6e}`."
    )
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
    (OUT / "scan_clebsch_flavor_symmetric_relaxation_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(__doc__.strip())
        return
    OUT.mkdir(parents=True, exist_ok=True)
    veronese_basis = load_basis()
    ortho = orthogonal_symmetric_direction(veronese_basis)
    s_perp = ortho["S_perp"]
    seed_base = load_seed()
    rng = np.random.default_rng(2026050722)
    rows = []
    print("Orthogonal-symmetric relaxation diagnostic", flush=True)
    print(
        f"  Veronese rank={ortho['veronese_rank']} "
        f"max_overlap={ortho['orthogonality_max_abs']:.3e}",
        flush=True,
    )
    for branch in ["F", "H"]:
        current = np.concatenate([seed_base, np.zeros(2, dtype=float)])
        for eps_budget in EPS_BUDGETS:
            for stage in STAGES:
                item = run_stage(stage, current, veronese_basis, s_perp, branch, eps_budget, rng)
                rows.append(item)
                current = np.array(item["x"], dtype=float)
                obs = item["CKM_observables"]
                scores = item["scores"]
                data = item["model_data"]
                print(
                    f"  {branch} eps={eps_budget:.2f} {stage['name']:7s} "
                    f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
                    f"ckm={scores['ckm_magnitude_log_score']:.3e} "
                    f"mass={scores['mass_log_score']:.3e} "
                    f"|eps|={data['epsilon_abs']:.3e} rel={data['relative_leakage_norm']:.3e}",
                    flush=True,
                )

    viable = [
        item
        for item in rows
        if item["scores"]["ckm_magnitude_log_score"] < 5.0e-2
        and item["scores"]["mass_log_score"] < 2.0e-1
        and item["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10
        and item["seesaw_replay"]["theta_norm"] < 1.0e-8
    ]
    if viable:
        best = min(
            viable,
            key=lambda item: (
                item["model_data"]["relative_leakage_norm"],
                item["scores"]["mass_log_score"],
                item["scores"]["ckm_magnitude_log_score"],
            ),
        )
    else:
        best = min(rows, key=lambda item: item["scores"]["mass_log_score"] + 4.0 * item["scores"]["ckm_magnitude_log_score"])

    best_obs = best["CKM_observables"]
    best_scores = best["scores"]
    best_data = best["model_data"]
    if viable:
        verdict = (
            "A joint CKM/mass/seesaw point appears after adding the single "
            "orthogonal symmetric direction.  The minimum viable relative "
            f"leakage found is {best_data['relative_leakage_norm']:.6e} in the "
            f"{best['branch']} branch, with CKM score "
            f"{best_scores['ckm_magnitude_log_score']:.6e} and mass score "
            f"{best_scores['mass_log_score']:.6e}.  This turns the previous "
            "strict-Veronese no-go into a controlled deformation target that "
            "must be derived from a covariant geometric correction."
        )
    else:
        verdict = (
            "No joint viable point was found even after adding the single "
            "orthogonal symmetric direction within the scanned budgets.  The "
            f"best scalarized point has |Vus|={best_obs['Vus']:.6e}, "
            f"|Vcb|={best_obs['Vcb']:.6e}, |Vub|={best_obs['Vub']:.6e}, CKM "
            f"score={best_scores['ckm_magnitude_log_score']:.6e}, mass "
            f"score={best_scores['mass_log_score']:.6e}, and relative leakage "
            f"{best_data['relative_leakage_norm']:.6e}.  The obstruction is "
            "therefore stronger than a one-direction symmetric relaxation."
        )

    payload = {
        "note": "No web lookup used. One orthogonal symmetric non-Veronese direction added to H or F.",
        "input_two120": str(TWO120),
        "S_perp": fit.matrix_json(s_perp),
        "S_perp_coordinates_in_orthonormal_sym_basis": fit.complex_list_json(ortho["coordinates"]),
        "S_perp_diagnostics": {
            "veronese_rank": ortho["veronese_rank"],
            "singular_values": ortho["singular_values"],
            "orthogonality_max_abs": ortho["orthogonality_max_abs"],
            "symmetry_error": ortho["symmetry_error"],
            "fro_norm": ortho["fro_norm"],
        },
        "rows": rows,
        "best": best,
        "checks": {
            "scan_completed": True,
            "S_perp_rank_check_pass": bool(ortho["veronese_rank"] == 5),
            "S_perp_orthogonal_lt_1e_minus_10": bool(ortho["orthogonality_max_abs"] < 1.0e-10),
            "S_perp_symmetric_lt_1e_minus_12": bool(ortho["symmetry_error"] < 1.0e-12),
            "any_ckm_score_lt_0p05": any(item["scores"]["ckm_magnitude_log_score"] < 5.0e-2 for item in rows),
            "any_mass_score_lt_0p20": any(item["scores"]["mass_log_score"] < 2.0e-1 for item in rows),
            "any_joint_viable": bool(viable),
            "best_seesaw_residual_lt_1e_minus_10": bool(best["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10),
            "best_theta_lt_1e_minus_8": bool(best["seesaw_replay"]["theta_norm"] < 1.0e-8),
            "improves_mass_vs_two120": bool(best["scores"]["mass_log_score"] < 0.6440132),
        },
        "verdict": verdict,
    }
    payload["phenomenology_viable"] = payload["checks"]["any_joint_viable"]
    (OUT / "scan_clebsch_flavor_symmetric_relaxation.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(rows)
    write_report(payload)
    print(f"  best: {best['branch']} {best['stage']} eps={best['epsilon_component_budget']:.2f}", flush=True)
    print(f"  phenomenology viable: {payload['phenomenology_viable']}", flush=True)
    print(f"  wrote: {OUT}", flush=True)


if __name__ == "__main__":
    main()
