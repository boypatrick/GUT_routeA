#!/usr/bin/env python3
"""Scan a CP1/O(-4)-restricted Spin(10) Clebsch flavor fit.

No web lookup is used.  This is the first full-flavor baseline after the
one-sector 120_H, Q-only Kahler, and correlated Kahler scans failed to give a
phenomenological CKM fit.

The family tensors are restricted to the geometric spaces already used in the
paper:

  H,F in Sym^2 H^0(CP1,O(2)) projected through the O(-4) dual-density product
      map, hence a five-complex-dimensional Veronese subspace.

  G in Lambda^2 C^3, the antisymmetric 120_H family tensor.

The renormalizable Spin(10) Clebsch ansatz is

  Y_u  = H + r_u F + g_u G,
  Y_nu = H - 3 r_u F + g_nu G,
  Y_d  = r_d (H + F + g_d G),
  Y_e  = r_d (H - 3 F + g_e G).

The coefficients are complex doublet-mixing parameters.  The scan optimizes
mass ratios, CKM magnitudes, and the Jarlskog invariant, then replays the
inverse type-I seesaw for the resulting Y_nu,Y_e.
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
OUT = ROOT / "output" / "flavor_clebsch_scan"

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


def complex_list_json(v: np.ndarray) -> list[dict[str, float]]:
    return [cjson(z) for z in v]


def unpack_complex(x: np.ndarray, start: int, count: int) -> tuple[np.ndarray, int]:
    vals = x[start : start + 2 * count]
    return vals[:count] + 1j * vals[count:], start + 2 * count


def pack_complex(v: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(v), np.imag(v)]).astype(float)


def antisym_from_coeff(c: np.ndarray) -> np.ndarray:
    g = np.zeros((3, 3), dtype=complex)
    g[0, 1] = c[0]
    g[1, 0] = -c[0]
    g[0, 2] = c[1]
    g[2, 0] = -c[1]
    g[1, 2] = c[2]
    g[2, 1] = -c[2]
    return g


def coeff_from_symmetric_basis(m: np.ndarray, basis: np.ndarray) -> np.ndarray:
    design = basis.reshape(9, 5)
    coeff, *_ = np.linalg.lstsq(design, m.reshape(9), rcond=None)
    return coeff


def symmetric_from_coeff(c: np.ndarray, basis: np.ndarray) -> np.ndarray:
    return np.einsum("ijm,m->ij", basis, c)


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def ratios(s: np.ndarray) -> tuple[float, float]:
    largest = max(float(s[2]), 1.0e-30)
    return float(s[0] / largest), float(s[1] / largest)


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


def all_mass_ratios(y: dict[str, np.ndarray]) -> dict[str, tuple[float, float]]:
    out = {}
    for name, mat in y.items():
        _u, s = left_rotation(mat)
        out[name] = ratios(s)
    return out


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
    heavy = np.sort(np.linalg.svd(MR_GeV, compute_uv=False))
    residual = float(np.linalg.norm(m_light_reco - m_light_target) / np.linalg.norm(m_light_target))
    theta_norm = float(np.linalg.norm(mD_eV @ np.linalg.inv(MR_eV), ord=2))
    return {
        "heavy_neutrino_masses_GeV": [float(v) for v in heavy],
        "theta_norm": theta_norm,
        "seesaw_matrix_residual": residual,
        "Y_e_singular_values_normalized": [float(v) for v in charged_singulars],
        "MR_condition_number": float(heavy[-1] / max(heavy[0], 1.0e-30)),
    }


def decode(x: np.ndarray, basis: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    pos = 0
    h_coeff, pos = unpack_complex(x, pos, 5)
    f_coeff, pos = unpack_complex(x, pos, 5)
    g_coeff, pos = unpack_complex(x, pos, 3)
    mix, pos = unpack_complex(x, pos, 6)
    assert pos == len(x)

    h = symmetric_from_coeff(h_coeff, basis)
    f = symmetric_from_coeff(f_coeff, basis)
    g = antisym_from_coeff(g_coeff)
    r_u, r_d, g_u, g_d, g_e, g_nu = mix
    y = {
        "up": h + r_u * f + g_u * g,
        "down": r_d * (h + f + g_d * g),
        "charged_lepton": r_d * (h - 3.0 * f + g_e * g),
        "neutrino_dirac": h - 3.0 * r_u * f + g_nu * g,
    }
    data = {
        "H_coefficients": complex_list_json(h_coeff),
        "F_coefficients": complex_list_json(f_coeff),
        "G_antisymmetric_coefficients": complex_list_json(g_coeff),
        "mixing_coefficients": {
            "r_u": cjson(r_u),
            "r_d": cjson(r_d),
            "g_u": cjson(g_u),
            "g_d": cjson(g_d),
            "g_e": cjson(g_e),
            "g_nu": cjson(g_nu),
        },
        "H_matrix": matrix_json(h),
        "F_matrix": matrix_json(f),
        "G_matrix": matrix_json(g),
        "norms": {
            "H_fro": float(np.linalg.norm(h)),
            "F_fro": float(np.linalg.norm(f)),
            "G_fro": float(np.linalg.norm(g)),
            "parameter_l2": float(np.linalg.norm(x)),
        },
    }
    return y, data


def objective(x: np.ndarray, basis: np.ndarray, weights: dict[str, float]) -> np.ndarray:
    y, _data = decode(x, basis)
    _v, ckm_obs, _qm = ckm_matrix(y["up"], y["down"])
    masses = all_mass_ratios(y)
    pieces = [
        weights["ckm"] * log10_ratio(ckm_obs["Vus"], TARGETS["Vus"]),
        weights["ckm"] * log10_ratio(ckm_obs["Vcb"], TARGETS["Vcb"]),
        weights["ckm"] * log10_ratio(ckm_obs["Vub"], TARGETS["Vub"]),
        weights["jarlskog"] * log10_ratio(ckm_obs["J"], TARGETS["J"]),
    ]
    for sector, (small, mid) in masses.items():
        pieces.append(weights["mass"] * log10_ratio(small, TARGETS[f"{sector}_small"]))
        pieces.append(weights["mass"] * log10_ratio(mid, TARGETS[f"{sector}_mid"]))
    pieces.append(weights["reg"] * float(np.linalg.norm(x)))
    return np.array(pieces, dtype=float)


def evaluate(x: np.ndarray, basis: np.ndarray, profile: str, bound: float) -> dict[str, Any]:
    y, data = decode(x, basis)
    v, ckm_obs, _qm = ckm_matrix(y["up"], y["down"])
    masses = all_mass_ratios(y)
    ckm_score = sum(log10_ratio(ckm_obs[key], TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = log10_ratio(ckm_obs["J"], TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += log10_ratio(small, TARGETS[f"{sector}_small"]) ** 2
        mass_score += log10_ratio(mid, TARGETS[f"{sector}_mid"]) ** 2
    replay = seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    return {
        "profile": profile,
        "component_bound": bound,
        "x": [float(vv) for vv in x],
        "model_data": data,
        "Yukawa_fit": {name: matrix_json(mat) for name, mat in y.items()},
        "CKM_abs": matrix_abs_json(v),
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


def seed_from_card(card: dict[str, Any], basis: np.ndarray) -> np.ndarray:
    sectors = card["yukawa_sectors"]
    y_u = cmat(sectors["up"]["Y_matrix_normalized"])
    y_d = cmat(sectors["down"]["Y_matrix_normalized"])
    y_e = cmat(sectors["charged_lepton"]["Y_matrix_normalized"])

    h0 = (3.0 * y_d + y_e) / 4.0
    f0 = (y_d - y_e) / 4.0
    h_coeff = coeff_from_symmetric_basis(h0, basis)
    f_coeff = coeff_from_symmetric_basis(f0, basis)
    f_norm = np.vdot(f0, f0)
    r_u = np.vdot(f0, y_u - h0) / f_norm if abs(f_norm) > 1.0e-30 else 0.0
    g_coeff = np.array([1.0e-3j, -1.0e-3, 1.0e-3 + 0.0j], dtype=complex)
    mix = np.array([r_u, 1.0 + 0.0j, 0.1 + 0.0j, 0.1 + 0.0j, 0.1 + 0.0j, 0.1 + 0.0j], dtype=complex)
    return np.concatenate([pack_complex(h_coeff), pack_complex(f_coeff), pack_complex(g_coeff), pack_complex(mix)])


def scan(
    profile: str,
    bound: float,
    basis: np.ndarray,
    base_seed: np.ndarray,
    rng: np.random.Generator,
    starts: int,
) -> dict[str, Any]:
    profiles = {
        "balanced": {"ckm": 1.3, "jarlskog": 0.25, "mass": 2.0, "reg": 0.008 / max(bound, 1.0e-9)},
        "ckm_heavy": {"ckm": 2.8, "jarlskog": 0.20, "mass": 1.0, "reg": 0.006 / max(bound, 1.0e-9)},
    }
    weights = profiles[profile]
    dim = len(base_seed)
    seeds = [np.clip(base_seed, -0.85 * bound, 0.85 * bound)]
    seeds.extend(rng.uniform(-bound, bound, size=dim) for _ in range(starts - 1))
    best = None
    for seed in seeds:
        res = least_squares(
            objective,
            seed,
            bounds=(-bound * np.ones(dim), bound * np.ones(dim)),
            args=(basis, weights),
            xtol=2.0e-10,
            ftol=2.0e-10,
            gtol=2.0e-10,
            max_nfev=900,
        )
        item = evaluate(res.x, basis, profile, bound)
        item["optimizer"] = {
            "success": bool(res.success),
            "nfev": int(res.nfev),
            "cost_weighted_sum_squares": float(np.sum(objective(res.x, basis, weights) ** 2)),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(res.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or item["optimizer"]["cost_weighted_sum_squares"] < best["optimizer"]["cost_weighted_sum_squares"]:
            best = item
    assert best is not None
    return best


def write_csv(rows: list[dict[str, Any]]) -> None:
    flat = []
    for row in rows:
        obs = row["CKM_observables"]
        scores = row["scores"]
        replay = row["seesaw_replay"]
        flat.append(
            {
                "profile": row["profile"],
                "component_bound": row["component_bound"],
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
                "weighted_cost": row["optimizer"]["cost_weighted_sum_squares"],
            }
        )
    with (OUT / "scan_clebsch_flavor_fit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0].keys()))
        writer.writeheader()
        writer.writerows(flat)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Clebsch-correlated Spin(10) flavor scan")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Ansatz")
    lines.append("")
    lines.append("```text")
    lines.append("Y_u  = H + r_u F + g_u G")
    lines.append("Y_nu = H - 3 r_u F + g_nu G")
    lines.append("Y_d  = r_d (H + F + g_d G)")
    lines.append("Y_e  = r_d (H - 3 F + g_e G)")
    lines.append("H,F in the five-dimensional CP1 O(-4) symmetric subspace")
    lines.append("G in Lambda^2 C^3")
    lines.append("```")
    lines.append("")
    lines.append("| profile | bound | Vus | Vcb | Vub | J | CKM score | mass score | theta |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["rows"]:
        obs = row["CKM_observables"]
        scores = row["scores"]
        lines.append(
            f"| {row['profile']} | {row['component_bound']:.2f} | {obs['Vus']:.4e} | "
            f"{obs['Vcb']:.4e} | {obs['Vub']:.4e} | {obs['J']:.4e} | "
            f"{scores['ckm_magnitude_log_score']:.3e} | {scores['mass_log_score']:.3e} | "
            f"{row['seesaw_replay']['theta_norm']:.3e} |"
        )
    best = payload["best"]
    obs = best["CKM_observables"]
    replay = best["seesaw_replay"]
    lines.append("")
    lines.append("## Best point")
    lines.append("")
    lines.append(f"Profile: `{best['profile']}`, bound `{best['component_bound']:.3g}`.")
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
    (OUT / "scan_clebsch_flavor_fit_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = json.loads(CARD.read_text(encoding="utf-8"))
    basis = np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )
    base_seed = seed_from_card(card, basis)
    rng = np.random.default_rng(2026050719)
    rows = []
    # Heartbeat-safe baseline: this is deliberately a finite reproducible scan,
    # not a final exhaustive global optimizer.  The output tells us whether the
    # Clebsch branch is promising enough to justify a longer run.
    for bound in [1.2, 2.4]:
        rows.append(scan("balanced", bound, basis, base_seed, rng, starts=5))
    for bound in [1.2, 2.4]:
        rows.append(scan("ckm_heavy", bound, basis, base_seed, rng, starts=5))
    best = min(rows, key=lambda row: row["optimizer"]["cost_weighted_sum_squares"])
    best_obs = best["CKM_observables"]
    verdict = (
        "The full Clebsch-correlated ansatz is the first branch with enough "
        "structure to test flavor and proton decay together.  This heartbeat-safe "
        f"baseline scan finds |Vus|={best_obs['Vus']:.6e}, "
        f"|Vcb|={best_obs['Vcb']:.6e}, |Vub|={best_obs['Vub']:.6e}, and "
        f"|J|={best_obs['J']:.6e}, so CKM can be made close.  However the mass "
        f"log-score is {best['scores']['mass_log_score']:.6e}, so the point is "
        "not yet a phenomenological flavor fit.  The next step is a Pareto or "
        "homotopy scan that starts from the CKM point and adiabatically raises "
        "the mass-hierarchy weights before any p -> K+ nu_bar claim is made."
    )
    payload = {
        "note": "No web lookup used. CP1/O(-4)-restricted 10H+126barH+120H Clebsch flavor scan.",
        "input_card": str(CARD),
        "targets": TARGETS,
        "ansatz": {
            "Y_u": "H + r_u F + g_u G",
            "Y_nu": "H - 3 r_u F + g_nu G",
            "Y_d": "r_d (H + F + g_d G)",
            "Y_e": "r_d (H - 3 F + g_e G)",
            "H_F_space": "five-dimensional CP1 O(-4) symmetric product subspace",
            "G_space": "Lambda^2 C^3 antisymmetric 120_H family tensor",
        },
        "rows": rows,
        "best": best,
        "internal_reconstruction_checks": {
            "scan_completed": True,
            "best_seesaw_residual_lt_1e_minus_10": bool(best["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10),
            "best_theta_lt_1e_minus_8": bool(best["seesaw_replay"]["theta_norm"] < 1.0e-8),
            "best_improves_ckm_vs_correlated_kahler": bool(best["scores"]["ckm_magnitude_log_score"] < 2.682300),
        },
        "phenomenology_viability_checks": {
            "best_ckm_score_lt_0p05": bool(best["scores"]["ckm_magnitude_log_score"] < 5.0e-2),
            "best_mass_score_lt_0p20": bool(best["scores"]["mass_log_score"] < 2.0e-1),
            "best_Vus_within_10_percent": bool(abs(best_obs["Vus"] / TARGETS["Vus"] - 1.0) < 0.10),
            "best_Vcb_within_10_percent": bool(abs(best_obs["Vcb"] / TARGETS["Vcb"] - 1.0) < 0.10),
            "best_Vub_within_10_percent": bool(abs(best_obs["Vub"] / TARGETS["Vub"] - 1.0) < 0.10),
        },
        "verdict": verdict,
    }
    payload["internal_reconstruction_ok"] = all(payload["internal_reconstruction_checks"].values())
    payload["phenomenology_viable"] = all(payload["phenomenology_viability_checks"].values())
    (OUT / "scan_clebsch_flavor_fit.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(rows)
    write_report(payload)
    print("Clebsch-correlated flavor scan")
    for row in rows:
        obs = row["CKM_observables"]
        print(
            f"  {row['profile']:9s} bound={row['component_bound']:.1f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"J={obs['J']:.4e} ckm={row['scores']['ckm_magnitude_log_score']:.3e} "
            f"mass={row['scores']['mass_log_score']:.3e} "
            f"theta={row['seesaw_replay']['theta_norm']:.3e}"
        )
    print(f"  best: {best['profile']} bound={best['component_bound']:.1f}")
    print(f"  internal reconstruction ok: {payload['internal_reconstruction_ok']}")
    print(f"  phenomenology viable: {payload['phenomenology_viable']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
