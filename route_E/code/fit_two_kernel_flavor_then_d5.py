#!/usr/bin/env python3
"""Fit a CP1 two-kernel flavor ansatz, then derive d=5 triplet tensors.

No web lookup is used.  This audit is deliberately stricter than the previous
matrix-deformation scans: the doublet Yukawa matrices and the colored-triplet
Wilson tensors are derived from the same operator-level ingredients.

The symmetric family tensors are

    H = V[h_H] + eps_H K,
    F = V[h_F] + eps_F K,

where V is the CP1 O(-4) Veronese dual-density product map and

    K = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]]

is the spin-zero second transvectant/contact functional.  The doublet sector is
the same Spin(10)-Clebsch ansatz used in the previous successful flavor card,
with two antisymmetric 120-like tensors G_A,G_B.  For proton decay, the triplet
QQ and QL couplings are not arbitrary matrix perturbations: a finite list of
Clebsch-phase profiles reuses the same H,F,G_A,G_B and doublet-mixing
coefficients.

The goal is not to publish a final flavor fit; it is to answer a sharper local
question: can an operator-level two-kernel ansatz simultaneously satisfy the
flavor scores and the calibrated Knu future-stress proxy?
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_clebsch_flavor_symmetric_relaxation as relax  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


RELAX = ROOT / "output" / "flavor_symmetric_relaxation" / "scan_clebsch_flavor_symmetric_relaxation.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "two_kernel_flavor_then_d5"

RNG_SEED = 202605082139
BASE_BOUND = 3.2
EPS_BOUND = 0.30
STRICT_CKM = 1.0e-3
STRICT_MASS_FACTOR = 1.25
LOOSE_CKM = 5.0e-2
LOOSE_MASS = 2.0e-1


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def analytic_k() -> np.ndarray:
    return np.array(
        [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        dtype=complex,
    ) / math.sqrt(3.0)


def unpack_complex(x: np.ndarray, start: int, count: int) -> tuple[np.ndarray, int]:
    vals = x[start : start + 2 * count]
    return vals[:count] + 1j * vals[count:], start + 2 * count


def pack_complex(v: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(v), np.imag(v)]).astype(float)


def load_veronese_basis() -> np.ndarray:
    card = read_json(fit.CARD)
    return np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )


def seed_from_relaxation() -> np.ndarray:
    payload = read_json(RELAX)
    best = payload["best"]
    s_raw = np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in payload["S_perp"]], dtype=complex)
    phase_to_analytic_k = np.vdot(analytic_k(), s_raw)
    x_old = np.array(best["x"], dtype=float)
    pos = 0
    h_coeff, pos = unpack_complex(x_old, pos, 5)
    f_coeff, pos = unpack_complex(x_old, pos, 5)
    ga_coeff, pos = unpack_complex(x_old, pos, 3)
    gb_coeff, pos = unpack_complex(x_old, pos, 3)
    mix, pos = unpack_complex(x_old, pos, 10)
    eps_arr, pos = unpack_complex(x_old, pos, 1)
    assert pos == len(x_old)
    eps_in_analytic_k_basis = eps_arr[0] * phase_to_analytic_k
    eps_h = eps_in_analytic_k_basis if best["branch"] == "H" else 0.0 + 0.0j
    eps_f = eps_in_analytic_k_basis if best["branch"] == "F" else 0.0 + 0.0j
    return np.concatenate(
        [
            pack_complex(h_coeff),
            pack_complex(f_coeff),
            pack_complex(ga_coeff),
            pack_complex(gb_coeff),
            pack_complex(mix),
            pack_complex(np.array([eps_h, eps_f], dtype=complex)),
        ]
    )


@dataclass(frozen=True)
class TripletProfile:
    name: str
    qf: complex
    lf: complex
    qa: complex
    qb: complex
    la: complex
    lb: complex


PROFILES = [
    TripletProfile("doublet_like", 1, 1, 1, 1, 1, 1),
    TripletProfile("ql_F_minus", 1, -1, 1, 1, 1, 1),
    TripletProfile("both_F_minus", -1, -1, 1, 1, 1, 1),
    TripletProfile("anti_flip_ql", 1, 1, 1, 1, -1, -1),
    TripletProfile("F_quarter_phase", 1j, -1j, 1, 1, 1, 1),
]


def profile_json(profile: TripletProfile) -> dict[str, Any]:
    return {
        "name": profile.name,
        "qf": cjson(profile.qf),
        "lf": cjson(profile.lf),
        "qa": cjson(profile.qa),
        "qb": cjson(profile.qb),
        "la": cjson(profile.la),
        "lb": cjson(profile.lb),
    }


def decode(x: np.ndarray, basis: np.ndarray, k: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    pos = 0
    h_coeff, pos = unpack_complex(x, pos, 5)
    f_coeff, pos = unpack_complex(x, pos, 5)
    ga_coeff, pos = unpack_complex(x, pos, 3)
    gb_coeff, pos = unpack_complex(x, pos, 3)
    mix, pos = unpack_complex(x, pos, 10)
    eps, pos = unpack_complex(x, pos, 2)
    assert pos == len(x)

    h_geo = fit.symmetric_from_coeff(h_coeff, basis)
    f_geo = fit.symmetric_from_coeff(f_coeff, basis)
    h = h_geo + eps[0] * k
    fmat = f_geo + eps[1] * k
    ga = fit.antisym_from_coeff(ga_coeff)
    gb = fit.antisym_from_coeff(gb_coeff)
    r_u, r_d, a_u, b_u, a_d, b_d, a_e, b_e, a_nu, b_nu = mix
    y = {
        "up": h + r_u * fmat + a_u * ga + b_u * gb,
        "down": r_d * (h + fmat + a_d * ga + b_d * gb),
        "charged_lepton": r_d * (h - 3.0 * fmat + a_e * ga + b_e * gb),
        "neutrino_dirac": h - 3.0 * r_u * fmat + a_nu * ga + b_nu * gb,
    }
    data = {
        "H_geo": h_geo,
        "F_geo": f_geo,
        "H": h,
        "F": fmat,
        "G_A": ga,
        "G_B": gb,
        "mix": {
            "r_u": r_u,
            "r_d": r_d,
            "a_u": a_u,
            "b_u": b_u,
            "a_d": a_d,
            "b_d": b_d,
            "a_e": a_e,
            "b_e": b_e,
            "a_nu": a_nu,
            "b_nu": b_nu,
        },
        "json": {
            "H_coefficients": fit.complex_list_json(h_coeff),
            "F_coefficients": fit.complex_list_json(f_coeff),
            "G_A_coefficients": fit.complex_list_json(ga_coeff),
            "G_B_coefficients": fit.complex_list_json(gb_coeff),
            "epsilon_H": fit.cjson(eps[0]),
            "epsilon_F": fit.cjson(eps[1]),
            "epsilon_H_abs": float(abs(eps[0])),
            "epsilon_F_abs": float(abs(eps[1])),
            "mixing_coefficients": {key: fit.cjson(val) for key, val in {
                "r_u": r_u, "r_d": r_d, "a_u": a_u, "b_u": b_u, "a_d": a_d,
                "b_d": b_d, "a_e": a_e, "b_e": b_e, "a_nu": a_nu, "b_nu": b_nu,
            }.items()},
            "norms": {
                "H_geo_fro": float(np.linalg.norm(h_geo)),
                "F_geo_fro": float(np.linalg.norm(f_geo)),
                "H_contact_relative": float(abs(eps[0]) * np.linalg.norm(k) / max(np.linalg.norm(h_geo), 1.0e-30)),
                "F_contact_relative": float(abs(eps[1]) * np.linalg.norm(k) / max(np.linalg.norm(f_geo), 1.0e-30)),
                "parameter_l2": float(np.linalg.norm(x)),
            },
        },
    }
    return y, data


def left_right_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, s, vh = np.linalg.svd(y, full_matrices=True)
    order = np.argsort(s)
    s = s[order]
    left = u[:, order]
    right = vh.conjugate().T[:, order]
    diag = left.conjugate().T @ y @ right
    phases = np.array([1.0 if abs(diag[i, i]) < 1.0e-30 else diag[i, i] / abs(diag[i, i]) for i in range(3)], dtype=complex)
    right = right @ np.diag(phases.conjugate())
    return left, right, s


def pmns_target(u_e: np.ndarray) -> np.ndarray:
    import verify_seesaw_item3 as seesaw

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
    u_pmns = seesaw.standard_pmns(
        benchmark.sin2_theta12,
        benchmark.sin2_theta13,
        benchmark.sin2_theta23,
        benchmark.delta_cp_rad,
        benchmark.alpha21_rad,
        benchmark.alpha31_rad,
    )
    return u_e @ u_pmns


def score_flavor(y: dict[str, np.ndarray]) -> dict[str, Any]:
    v, ckm, _mass_quark = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = fit.log10_ratio(ckm["J"], fit.TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    return {
        "CKM_abs": fit.matrix_abs_json(v),
        "CKM_observables": ckm,
        "mass_ratios": {name: {"small": vals[0], "mid": vals[1]} for name, vals in masses.items()},
        "scores": {
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "mass_log_score": float(mass_score),
            "total_unweighted": float(ckm_score + j_score + mass_score),
        },
    }


def rotations_for(y: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    ul: dict[str, np.ndarray] = {}
    ur: dict[str, np.ndarray] = {}
    for name, mat in y.items():
        left, right, _s = left_right_rotation(mat)
        ul[name] = left
        ur[name] = right
    return ul, ur, pmns_target(ul["charged_lepton"])


def triplet_matrices(data: dict[str, Any], profile: TripletProfile) -> tuple[np.ndarray, np.ndarray]:
    h = data["H"]
    fmat = data["F"]
    ga = data["G_A"]
    gb = data["G_B"]
    mix = data["mix"]
    qq_raw = h + profile.qf * mix["r_u"] * fmat + profile.qa * mix["a_u"] * ga + profile.qb * mix["b_u"] * gb
    ql_raw = mix["r_d"] * (h + profile.lf * fmat + profile.la * mix["a_d"] * ga + profile.lb * mix["b_d"] * gb)
    return d5.scale_to_largest_singular(qq_raw, 0.60), d5.scale_to_largest_singular(ql_raw, 0.024)


def knu_for(y: dict[str, np.ndarray], data: dict[str, Any], profile: TripletProfile) -> dict[str, Any]:
    ul, _ur, u_nu = rotations_for(y)
    y_qq, y_ql = triplet_matrices(data, profile)
    tensor = d5.tensor_llll(d5.sym(y_qq), y_ql, ul["up"], ul["up"], ul["down"], u_nu)
    amp, idx = d5.max_perm_entry(tensor, (0, 0, 1), None)
    return {"amplitude": float(amp), "index": "".join(str(i) for i in idx)}


def best_knu_profile(y: dict[str, np.ndarray], data: dict[str, Any], baseline_amp: float, future_margin0: float) -> dict[str, Any]:
    rows = []
    for profile in PROFILES:
        item = knu_for(y, data, profile)
        ratio = item["amplitude"] / max(baseline_amp, 1.0e-30)
        rows.append(
            {
                "profile": profile.name,
                "amplitude": item["amplitude"],
                "selected_index": item["index"],
                "amplitude_ratio_to_baseline": ratio,
                "future_margin_1e35": future_margin0 / max(ratio * ratio, 1.0e-30),
            }
        )
    return min(rows, key=lambda row: row["amplitude"])


def objective(
    x: np.ndarray,
    basis: np.ndarray,
    k: np.ndarray,
    weights: dict[str, float],
    profile: TripletProfile | None,
    baseline_amp: float,
    future_margin0: float,
) -> np.ndarray:
    y, data = decode(x, basis, k)
    flv = score_flavor(y)
    ckm = flv["CKM_observables"]
    masses = flv["mass_ratios"]
    pieces = [
        weights["ckm"] * fit.log10_ratio(ckm["Vus"], fit.TARGETS["Vus"]),
        weights["ckm"] * fit.log10_ratio(ckm["Vcb"], fit.TARGETS["Vcb"]),
        weights["ckm"] * fit.log10_ratio(ckm["Vub"], fit.TARGETS["Vub"]),
        weights["jarlskog"] * fit.log10_ratio(ckm["J"], fit.TARGETS["J"]),
    ]
    for sector, vals in masses.items():
        pieces.append(weights["mass"] * fit.log10_ratio(vals["small"], fit.TARGETS[f"{sector}_small"]))
        pieces.append(weights["mass"] * fit.log10_ratio(vals["mid"], fit.TARGETS[f"{sector}_mid"]))
    pieces.append(weights["reg"] * float(np.linalg.norm(x[:-4])))
    pieces.append(weights["epsilon"] * float(np.linalg.norm(x[-4:])))
    if profile is not None and weights["d5"] > 0.0:
        amp = knu_for(y, data, profile)["amplitude"]
        amp_ratio = amp / max(baseline_amp, 1.0e-30)
        target_ratio = math.sqrt(future_margin0)
        pieces.append(weights["d5"] * max(0.0, math.log10(max(amp_ratio, 1.0e-30) / target_ratio)))
    return np.array(pieces, dtype=float)


def bounds_for(seed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lower = -BASE_BOUND * np.ones_like(seed)
    upper = BASE_BOUND * np.ones_like(seed)
    lower[-4:] = -EPS_BOUND
    upper[-4:] = EPS_BOUND
    return lower, upper


def evaluate(
    label: str,
    x: np.ndarray,
    basis: np.ndarray,
    k: np.ndarray,
    profile_hint: str,
    baseline_amp: float,
    future_margin0: float,
) -> dict[str, Any]:
    y, data = decode(x, basis, k)
    flv = score_flavor(y)
    replay = fit.seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    best_profile = best_knu_profile(y, data, baseline_amp, future_margin0)
    return {
        "label": label,
        "profile_hint": profile_hint,
        "x": [float(v) for v in x],
        "model_data": data["json"],
        "Yukawa_fit": {name: fit.matrix_json(mat) for name, mat in y.items()},
        **flv,
        "seesaw_replay": replay,
        "best_triplet_profile": best_profile,
    }


def run_fit(
    label: str,
    seed: np.ndarray,
    basis: np.ndarray,
    k: np.ndarray,
    weights: dict[str, float],
    profile: TripletProfile | None,
    baseline_amp: float,
    future_margin0: float,
    rng: np.random.Generator,
) -> dict[str, Any]:
    lower, upper = bounds_for(seed)
    # Deterministic local audit: the random generator is retained for metadata
    # continuity, but the heartbeat run uses one seeded refinement per profile.
    del rng
    candidates = [np.clip(seed, 0.95 * lower, 0.95 * upper)]
    best = None
    for idx, candidate in enumerate(candidates):
        res = least_squares(
            objective,
            candidate,
            bounds=(lower, upper),
            args=(basis, k, weights, profile, baseline_amp, future_margin0),
            xtol=4.0e-10,
            ftol=4.0e-10,
            gtol=4.0e-10,
            max_nfev=260,
        )
        item = evaluate(f"{label}_{idx}", res.x, basis, k, profile.name if profile else "flavor_only", baseline_amp, future_margin0)
        item["optimizer"] = {
            "success": bool(res.success),
            "nfev": int(res.nfev),
            "weighted_cost": float(np.sum(objective(res.x, basis, k, weights, profile, baseline_amp, future_margin0) ** 2)),
            "active_base_bound_fraction": float(np.mean(np.isclose(np.abs(res.x[:-4]), BASE_BOUND, rtol=0.0, atol=1.0e-5))),
            "active_epsilon_bound_fraction": float(np.mean(np.isclose(np.abs(res.x[-4:]), EPS_BOUND, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or item["optimizer"]["weighted_cost"] < best["optimizer"]["weighted_cost"]:
            best = item
    assert best is not None
    return best


def flat_row(item: dict[str, Any], baseline_mass: float) -> dict[str, Any]:
    scores = item["scores"]
    obs = item["CKM_observables"]
    prof = item["best_triplet_profile"]
    replay = item["seesaw_replay"]
    return {
        "label": item["label"],
        "profile_hint": item["profile_hint"],
        "best_triplet_profile": prof["profile"],
        "Vus": obs["Vus"],
        "Vcb": obs["Vcb"],
        "Vub": obs["Vub"],
        "J_abs": obs["J"],
        "ckm_score": scores["ckm_magnitude_log_score"],
        "j_score": scores["jarlskog_log_score"],
        "mass_score": scores["mass_log_score"],
        "mass_score_ratio": scores["mass_log_score"] / max(baseline_mass, 1.0e-30),
        "knu_amp_ratio": prof["amplitude_ratio_to_baseline"],
        "future_margin": prof["future_margin_1e35"],
        "epsilon_H_abs": item["model_data"]["epsilon_H_abs"],
        "epsilon_F_abs": item["model_data"]["epsilon_F_abs"],
        "H_contact_relative": item["model_data"]["norms"]["H_contact_relative"],
        "F_contact_relative": item["model_data"]["norms"]["F_contact_relative"],
        "theta_norm": replay["theta_norm"],
        "seesaw_residual": replay["seesaw_matrix_residual"],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    basis = load_veronese_basis()
    k = analytic_k()
    seed = seed_from_relaxation()
    future_margin0 = float(read_json(KNU_TARGET)["verdict"]["future_margin"])

    y0, data0 = decode(seed, basis, k)
    flv0 = score_flavor(y0)
    baseline_amp = knu_for(y0, data0, PROFILES[0])["amplitude"]
    baseline_mass = flv0["scores"]["mass_log_score"]

    rows = [evaluate("seed_relaxation", seed, basis, k, "seed", baseline_amp, future_margin0)]

    flavor_weights = {"ckm": 7.0, "jarlskog": 0.25, "mass": 38.0, "reg": 0.0025 / BASE_BOUND, "epsilon": 0.012 / EPS_BOUND, "d5": 0.0}
    flavor_best = run_fit("flavor_refit", seed, basis, k, flavor_weights, None, baseline_amp, future_margin0, rng)
    rows.append(flavor_best)

    # The d5 refinement uses finite Clebsch profiles only; no free triplet
    # matrices are introduced.
    current_seed = np.array(flavor_best["x"], dtype=float)
    d5_weights = {"ckm": 7.0, "jarlskog": 0.25, "mass": 38.0, "reg": 0.0025 / BASE_BOUND, "epsilon": 0.010 / EPS_BOUND, "d5": 4.5}
    for profile in PROFILES:
        rows.append(run_fit(f"d5_{profile.name}", current_seed, basis, k, d5_weights, profile, baseline_amp, future_margin0, rng))

    flat = [flat_row(item, baseline_mass) for item in rows]
    with (OUT / "two_kernel_summary_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0].keys()))
        writer.writeheader()
        writer.writerows(flat)

    strict = [
        row for row in flat
        if row["ckm_score"] < STRICT_CKM
        and row["mass_score"] <= STRICT_MASS_FACTOR * baseline_mass
        and row["future_margin"] >= 1.0
    ]
    loose = [
        row for row in flat
        if row["ckm_score"] < LOOSE_CKM
        and row["mass_score"] <= LOOSE_MASS
        and row["future_margin"] >= 1.0
    ]
    best_ckm = min(flat, key=lambda row: row["ckm_score"])
    best_mass = min(flat, key=lambda row: row["mass_score"])
    best_knu = max(flat, key=lambda row: row["future_margin"])
    best_balanced = min(flat, key=lambda row: row["ckm_score"] / STRICT_CKM + row["mass_score"] / max(baseline_mass, 1.0e-30) + max(0.0, 1.0 - row["future_margin"]))

    keep_labels = {best_ckm["label"], best_mass["label"], best_knu["label"], best_balanced["label"], "seed_relaxation"}
    item_by_label = {item["label"]: item for item in rows}
    matrix_cards = {
        label: {
            "Yukawa_fit": item_by_label[label]["Yukawa_fit"],
            "model_data": item_by_label[label]["model_data"],
            "best_triplet_profile": item_by_label[label]["best_triplet_profile"],
        }
        for label in sorted(keep_labels)
        if label in item_by_label
    }

    verdict = {
        "rows": len(flat),
        "strict_closure_count": len(strict),
        "loose_closure_count": len(loose),
        "baseline_ckm_score": flv0["scores"]["ckm_magnitude_log_score"],
        "baseline_mass_score": baseline_mass,
        "baseline_future_margin": future_margin0,
        "best_ckm": best_ckm,
        "best_mass": best_mass,
        "best_knu": best_knu,
        "best_balanced": best_balanced,
        "interpretation": (
            "strict_closure_count tests the publication-grade local target "
            "CKM<1e-3, mass within 25 percent of the seed, and future Knu "
            "margin>1.  loose_closure_count tests whether the earlier "
            "phenomenology-viable flavor threshold can be kept while making "
            "the derived triplet tensor future-safe.  Because triplet tensors "
            "are finite Clebsch reuses of H,F,G_A,G_B, any passing row would be "
            "an operator-level candidate rather than an a posteriori matrix "
            "filter."
        ),
    }
    payload = {
        "note": "No web lookup used. CP1 two-kernel operator-level flavor fit followed by derived d=5 triplet tensor audit.",
        "mathematical_setup": {
            "H": "V[h_H] + eps_H K",
            "F": "V[h_F] + eps_F K",
            "K": fit.matrix_json(k),
            "triplet_profiles": [profile_json(profile) for profile in PROFILES],
        },
        "config": {
            "rng_seed": RNG_SEED,
            "base_bound": BASE_BOUND,
            "eps_bound": EPS_BOUND,
            "strict_ckm": STRICT_CKM,
            "strict_mass_factor": STRICT_MASS_FACTOR,
            "loose_ckm": LOOSE_CKM,
            "loose_mass": LOOSE_MASS,
        },
        "rows": flat,
        "matrix_cards": matrix_cards,
        "verdict": verdict,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Two-kernel flavor then d=5 audit",
        "",
        "No web lookup was used.",
        "",
        "Operator-level ansatz:",
        "",
        "```text",
        "H = V[h_H] + eps_H K",
        "F = V[h_F] + eps_F K",
        "Y_QQ^T, Y_QL^T derived from the same H,F,G_A,G_B with finite Clebsch phases.",
        "```",
        "",
        "| role | profile hint | best triplet | CKM score | mass score | Knu margin | epsH | epsF |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    selected = [
        ("seed", flat[0]),
        ("best_ckm", best_ckm),
        ("best_mass", best_mass),
        ("best_knu", best_knu),
        ("best_balanced", best_balanced),
    ]
    seen: set[str] = set()
    for role, row in selected:
        key = role + row["label"]
        if key in seen:
            continue
        seen.add(key)
        report.append(
            f"| `{role}` | `{row['profile_hint']}` | `{row['best_triplet_profile']}` | "
            f"{row['ckm_score']:.6e} | {row['mass_score']:.6e} | "
            f"{row['future_margin']:.6e} | {row['epsilon_H_abs']:.3e} | {row['epsilon_F_abs']:.3e} |"
        )
    report += [
        "",
        "## Verdict",
        "",
        f"Strict closures: `{len(strict)}`.",
        f"Loose closures: `{len(loose)}`.",
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Two-kernel flavor then d=5 audit")
    print(f"  rows: {len(flat)}")
    print(f"  strict closures: {len(strict)}")
    print(f"  loose closures: {len(loose)}")
    print(f"  best CKM score: {best_ckm['ckm_score']:.6e}")
    print(f"  best mass score: {best_mass['mass_score']:.6e}")
    print(f"  best Knu margin: {best_knu['future_margin']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
