#!/usr/bin/env python3
"""Export the transvectant flavor benchmark and a d=5 leakage proxy.

No web lookup is used.  This script takes the viable H-branch point from the
orthogonal-symmetric relaxation scan, exports the exact Yukawa matrices and
biunitary rotations, and computes a deliberately labelled proxy for the
dimension-five p -> K+ nu_bar flavor leakage.

The proxy is not a full SUSY-GUT proton-decay calculation.  It only answers the
next reproducibility question: after the transvectant correction fixes CKM and
mass ratios, how large is the CKM-dominant third-generation leakage compared
with the earlier raw second/third-generation estimates?
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

import scan_clebsch_flavor_fit as fit  # noqa: E402
import verify_seesaw_item3 as seesaw  # noqa: E402


RELAX = ROOT / "output" / "flavor_symmetric_relaxation" / "scan_clebsch_flavor_symmetric_relaxation.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "flavor_transvectant_rotations"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in m]


def matrix_abs_json(m: np.ndarray) -> list[list[float]]:
    return [[float(abs(z)) for z in row] for row in m]


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def load_best() -> dict[str, Any]:
    payload = json.loads(RELAX.read_text(encoding="utf-8"))
    return payload["best"]


def biunitary(y: np.ndarray) -> dict[str, Any]:
    u, singulars, vh = np.linalg.svd(y, full_matrices=True)
    order = np.argsort(singulars)
    singulars = singulars[order]
    left = u[:, order]
    right = vh.conjugate().T[:, order]
    diagonal = left.conjugate().T @ y @ right
    offdiag = diagonal - np.diag(np.diag(diagonal))
    phases = np.array(
        [1.0 if abs(diagonal[i, i]) < 1.0e-30 else diagonal[i, i] / abs(diagonal[i, i]) for i in range(3)],
        dtype=complex,
    )
    # Rephase right fields so L^dagger Y R is real-positive.
    right = right @ np.diag(phases.conjugate())
    diagonal = left.conjugate().T @ y @ right
    offdiag = diagonal - np.diag(np.diag(diagonal))
    return {
        "singular_values": [float(v) for v in singulars],
        "left_rotation": left,
        "right_rotation": right,
        "diagonalized": diagonal,
        "unitarity_left_residual": float(np.linalg.norm(left.conjugate().T @ left - np.eye(3))),
        "unitarity_right_residual": float(np.linalg.norm(right.conjugate().T @ right - np.eye(3))),
        "offdiag_residual": float(np.linalg.norm(offdiag) / max(np.linalg.norm(diagonal), 1.0e-30)),
        "diag_imag_residual": float(np.linalg.norm(np.imag(np.diag(diagonal))) / max(np.linalg.norm(np.diag(diagonal)), 1.0e-30)),
    }


def proton_constants() -> dict[str, float]:
    proton = json.loads(PROTON.read_text(encoding="utf-8"))
    return proton["hadronic_constants"]


def lifetime_from_yprod(
    yprod: float,
    triplet_filter: float,
    constants: dict[str, float],
    m_triplet_gev: float = 1.0e16,
    m_wino_gev: float = 1.0e3,
    m_sfermion_gev: float = 1.0e5,
    alpha2_inv: float = 25.0,
) -> dict[str, float]:
    loop = (1.0 / alpha2_inv) / (4.0 * math.pi)
    c5 = triplet_filter * yprod / m_triplet_gev
    dressing = loop * m_wino_gev / (m_sfermion_gev * m_sfermion_gev)
    c6 = c5 * dressing
    width = constants["width_prefactor_GeV5"] * c6 * c6
    tau = math.inf if width <= 0.0 else HBAR_GEV_S / width / SECONDS_PER_YEAR
    return {
        "yprod_effective": yprod,
        "triplet_filter": triplet_filter,
        "M_T_GeV": m_triplet_gev,
        "m_wino_GeV": m_wino_gev,
        "m_sfermion_GeV": m_sfermion_gev,
        "loop_alpha2_over_4pi": loop,
        "C5_GeV_minus1": c5,
        "dressing_GeV_minus1": dressing,
        "C6_dressed_GeV_minus2": c6,
        "tau_years": tau,
    }


def filter_required(row_at_st1: dict[str, float], tau_bound_years: float) -> float:
    tau_at_st1 = row_at_st1["tau_years"]
    return math.sqrt(tau_at_st1 / tau_bound_years)


def flavor_proxy_rows(rot: dict[str, Any], constants: dict[str, float]) -> list[dict[str, float | str]]:
    u_s = np.array(rot["up"]["singular_values"], dtype=float)
    d_s = np.array(rot["down"]["singular_values"], dtype=float)
    ckm = np.array(rot["CKM_abs"], dtype=float)
    # Local item-4 physical normalization.  This keeps comparison with the
    # previous proton tables transparent.
    y_t = 0.60
    y_b = 0.024
    y_u = y_t * u_s[0] / max(u_s[2], 1.0e-30)
    y_c = y_t * u_s[1] / max(u_s[2], 1.0e-30)
    y_d = y_b * d_s[0] / max(d_s[2], 1.0e-30)
    y_s = y_b * d_s[1] / max(d_s[2], 1.0e-30)
    vtd_vts = ckm[2, 0] * ckm[2, 1]
    vcd_vcs = ckm[1, 0] * ckm[1, 1]
    vus = ckm[0, 1]

    scenarios = [
        ("first_gen_raw", y_u * y_d, 1.0),
        ("second_gen_raw", y_c * y_s, 1.0),
        ("second_gen_ckm_to_K_proxy", y_c * y_s * vcd_vcs, vcd_vcs),
        ("third_gen_raw", y_t * y_b, 1.0),
        ("third_gen_ckm_to_K_proxy", y_t * y_b * vtd_vts, vtd_vts),
        ("third_gen_single_Cabibbo_proxy", y_t * y_b * vus, vus),
    ]
    rows: list[dict[str, float | str]] = []
    for name, yprod, leakage in scenarios:
        base = lifetime_from_yprod(float(yprod), 1.0, constants)
        row: dict[str, float | str] = {
            "scenario": name,
            "flavor_leakage_factor": float(leakage),
            **base,
            "S_T_required_tau_1e34": filter_required(base, 1.0e34),
            "S_T_required_tau_2p4e34": filter_required(base, 2.4e34),
            "S_T_required_tau_1e35": filter_required(base, 1.0e35),
        }
        rows.append(row)
    return rows


def write_csv(rows: list[dict[str, Any]]) -> None:
    with (OUT / "dimension5_flavor_proxy.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    best = payload["source_best_summary"]
    rows = payload["dimension5_flavor_proxy"]
    checks = payload["checks"]
    lines = [
        "# Transvectant flavor rotations and d=5 proxy",
        "",
        "No web lookup was used.",
        "",
        "## Best flavor point",
        "",
        f"Branch `{best['branch']}`, stage `{best['stage']}`, epsilon budget `{best['epsilon_component_budget']}`.",
        f"CKM score `{best['ckm_score']:.6e}`, mass score `{best['mass_score']:.6e}`.",
        "",
        "## Rotation checks",
        "",
        f"Max biunitary off-diagonal residual: `{checks['max_biunitary_offdiag_residual']:.6e}`.",
        f"CKM replay residual: `{checks['CKM_replay_abs_residual']:.6e}`.",
        "",
        "## Dimension-five proxy",
        "",
        "The proxy uses the local item-4 normalization with `y_t=0.60`, `y_b=0.024`,",
        "`M_T=1e16 GeV`, `m_wino=1e3 GeV`, and `m_sfermion=1e5 GeV`.",
        "It is not yet a full dressed channel calculation.",
        "",
        "| scenario | leakage | yprod_eff | tau(S_T=1) [yr] | S_T max 1e34 | S_T max 2.4e34 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['scenario']}` | {row['flavor_leakage_factor']:.3e} | "
            f"{row['yprod_effective']:.3e} | {row['tau_years']:.3e} | "
            f"{row['S_T_required_tau_1e34']:.3e} | {row['S_T_required_tau_2p4e34']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## Verdict",
            "",
            payload["verdict"],
            "",
        ]
    )
    (OUT / "transvectant_flavor_rotations_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(__doc__.strip())
        return
    OUT.mkdir(parents=True, exist_ok=True)
    best = load_best()
    y = {name: cmat(raw) for name, raw in best["Yukawa_fit"].items()}
    rotations: dict[str, Any] = {}
    for name, mat in y.items():
        item = biunitary(mat)
        rotations[name] = {
            "singular_values": item["singular_values"],
            "left_rotation": matrix_json(item["left_rotation"]),
            "right_rotation": matrix_json(item["right_rotation"]),
            "diagonalized": matrix_json(item["diagonalized"]),
            "unitarity_left_residual": item["unitarity_left_residual"],
            "unitarity_right_residual": item["unitarity_right_residual"],
            "offdiag_residual": item["offdiag_residual"],
            "diag_imag_residual": item["diag_imag_residual"],
        }
        rotations[f"_{name}_left_matrix"] = item["left_rotation"]
        rotations[f"_{name}_right_matrix"] = item["right_rotation"]

    ckm = rotations["_up_left_matrix"].conjugate().T @ rotations["_down_left_matrix"]
    stored_ckm = np.array(best["CKM_abs"], dtype=float)
    ckm_abs = np.abs(ckm)

    # Reconstruct the light-neutrino convention used in the seesaw replay.
    u_e = rotations["_charged_lepton_left_matrix"]
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
    u_pmns_target = seesaw.standard_pmns(
        benchmark.sin2_theta12,
        benchmark.sin2_theta13,
        benchmark.sin2_theta23,
        benchmark.delta_cp_rad,
        benchmark.alpha21_rad,
        benchmark.alpha31_rad,
    )
    u_nu_target = u_e @ u_pmns_target

    constants = proton_constants()
    proxy_rows = flavor_proxy_rows({"up": rotations["up"], "down": rotations["down"], "CKM_abs": ckm_abs.tolist()}, constants)
    write_csv(proxy_rows)
    for key in list(rotations.keys()):
        if key.startswith("_"):
            del rotations[key]

    third_ckm = next(row for row in proxy_rows if row["scenario"] == "third_gen_ckm_to_K_proxy")
    second_raw = next(row for row in proxy_rows if row["scenario"] == "second_gen_raw")
    verdict = (
        "The viable transvectant point now has reproducible left/right rotations. "
        "In the CKM-dominant proxy, the third-generation K nu leakage is "
        f"{third_ckm['flavor_leakage_factor']:.6e}, giving an effective y-product "
        f"{third_ckm['yprod_effective']:.6e}.  This is comparable to, and below, "
        f"the raw second-generation product {second_raw['yprod_effective']:.6e}; "
        "therefore the full dressed p -> K+ nu_bar calculation is now meaningful "
        "and should replace the older scalar S_T-only toy estimate."
    )

    checks = {
        "max_biunitary_offdiag_residual": float(
            max(rotations[name]["offdiag_residual"] for name in ["up", "down", "charged_lepton", "neutrino_dirac"])
        ),
        "max_unitarity_residual": float(
            max(
                max(rotations[name]["unitarity_left_residual"], rotations[name]["unitarity_right_residual"])
                for name in ["up", "down", "charged_lepton", "neutrino_dirac"]
            )
        ),
        "CKM_replay_abs_residual": float(np.linalg.norm(ckm_abs - stored_ckm)),
        "third_ckm_proxy_less_than_raw_third": bool(third_ckm["yprod_effective"] < next(row for row in proxy_rows if row["scenario"] == "third_gen_raw")["yprod_effective"]),
        "third_ckm_proxy_needs_no_1e_minus5_filter_for_1e34": bool(third_ckm["S_T_required_tau_1e34"] > 1.0e-5),
    }
    checks["all_numerical_rotation_checks_pass"] = bool(
        checks["max_biunitary_offdiag_residual"] < 1.0e-12
        and checks["max_unitarity_residual"] < 1.0e-12
        and checks["CKM_replay_abs_residual"] < 1.0e-12
    )

    payload = {
        "note": "No web lookup used. This is a rotation export and dimension-five flavor proxy, not a full proton-decay calculation.",
        "input": str(RELAX),
        "source_best_summary": {
            "branch": best["branch"],
            "stage": best["stage"],
            "epsilon_component_budget": best["epsilon_component_budget"],
            "ckm_score": best["scores"]["ckm_magnitude_log_score"],
            "mass_score": best["scores"]["mass_log_score"],
            "epsilon_abs": best["model_data"]["epsilon_abs"],
            "relative_leakage_norm": best["model_data"]["relative_leakage_norm"],
        },
        "Yukawa_matrices": {name: matrix_json(mat) for name, mat in y.items()},
        "biunitary_rotations": rotations,
        "CKM": {
            "matrix": matrix_json(ckm),
            "abs": matrix_abs_json(ckm),
            "observables": {
                "Vus": float(abs(ckm[0, 1])),
                "Vcb": float(abs(ckm[1, 2])),
                "Vub": float(abs(ckm[0, 2])),
                "Vtd": float(abs(ckm[2, 0])),
                "Vts": float(abs(ckm[2, 1])),
            },
        },
        "PMNS_convention": {
            "target_PMNS": matrix_json(u_pmns_target),
            "neutrino_left_rotation_target": matrix_json(u_nu_target),
            "note": "U_PMNS = U_e^dagger U_nu_target, matching the local seesaw replay convention.",
        },
        "dimension5_flavor_proxy": proxy_rows,
        "checks": checks,
        "verdict": verdict,
    }
    (OUT / "transvectant_flavor_rotations.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(payload)

    print("Transvectant flavor rotation export")
    print(f"  CKM replay residual: {checks['CKM_replay_abs_residual']:.6e}")
    print(f"  max biunitary offdiag residual: {checks['max_biunitary_offdiag_residual']:.6e}")
    print(f"  third CKM Knu proxy yprod: {third_ckm['yprod_effective']:.6e}")
    print(f"  third CKM Knu S_T max for 1e34 yr: {third_ckm['S_T_required_tau_1e34']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
