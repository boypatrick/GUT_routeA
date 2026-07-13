#!/usr/bin/env python3
"""Export an action-level CP1/O(2) flavor benchmark card.

No web lookup is used.  This script turns the covariant O(-4) dual-density
interpretation of the legacy Yukawa scan into machine-readable benchmark data:

  * the O(4) product Clebsch coefficients C_ij^m,
  * the dual-density coefficients a_m for each Yukawa sector,
  * normalized Yukawa matrices Y_u,Y_d,Y_e,Y_nu,
  * the type-I seesaw reconstruction from the exact Y_nu and Y_e matrices.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_cp1_yukawa_covariant as cov  # noqa: E402
import scan_cp1_o2_yukawa as legacy  # noqa: E402
import verify_seesaw_item3 as seesaw  # noqa: E402


OUT = ROOT / "output" / "flavor_benchmark"
LEGACY_JSON = ROOT / "output" / "yukawa_o2" / "cp1_o2_yukawa_scan.json"


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in m]


def tensor3_json(t: np.ndarray) -> list[list[list[dict[str, float]]]]:
    return [[[cjson(z) for z in row] for row in plane] for plane in t]


def complex_list_json(v: np.ndarray) -> list[dict[str, float]]:
    return [cjson(z) for z in v]


def matrix_markdown(m: np.ndarray) -> str:
    lines = []
    for row in m:
        lines.append("[ " + ", ".join(f"{z.real:+.6e}{z.imag:+.6e}i" for z in row) + " ]")
    return "\n".join(lines)


def load_legacy_payload() -> dict[str, Any]:
    return json.loads(LEGACY_JSON.read_text(encoding="utf-8"))


def kernel_params(item: dict[str, Any]) -> legacy.KernelParams:
    data = {key: float(value) for key, value in item["params"].items() if key != "center_separation"}
    return legacy.KernelParams(**data)


def export_sector(
    grid: legacy.CP1O2Grid,
    chi: np.ndarray,
    coeff: np.ndarray,
    name: str,
    item: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray]:
    params = kernel_params(item)
    h = grid.higgs_profile(params)
    direct = grid.y_matrix(params)
    singular_values, norm = legacy.singular_data(direct)
    scale = float(singular_values[0])

    dual_raw = np.array([np.sum(grid.weight * chi[m] * h) for m in range(5)], dtype=complex)
    dual_direct = np.einsum("ijm,m->ij", coeff, dual_raw)
    y_normalized = direct / scale
    dual_normalized = dual_raw / scale
    y_from_normalized_dual = np.einsum("ijm,m->ij", coeff, dual_normalized)

    sector_card = {
        "target": item["target"],
        "params": item["params"],
        "normalization": {
            "largest_singular_value_raw": scale,
            "Y_matrix_is_normalized_to_largest_singular_value_one": True,
        },
        "dual_coefficients_raw": complex_list_json(dual_raw),
        "dual_coefficients_normalized": complex_list_json(dual_normalized),
        "singular_values_raw": [float(x) for x in singular_values],
        "normalized_singular_values": [float(x) for x in norm],
        "Y_matrix_normalized": matrix_json(y_normalized),
        "Y_matrix_markdown": matrix_markdown(y_normalized),
        "checks": {
            "direct_vs_dual_raw_max_abs": float(np.max(np.abs(direct - dual_direct))),
            "direct_vs_dual_normalized_max_abs": float(np.max(np.abs(y_normalized - y_from_normalized_dual))),
            "veronese_relation_abs_normalized": float(abs(y_normalized[1, 1] - 2.0 * y_normalized[0, 2])),
        },
    }
    return sector_card, y_normalized


def seesaw_card_from_yukawas(y_nu: np.ndarray, y_e: np.ndarray) -> dict[str, Any]:
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
    mD_GeV = y_nu * mD_largest_GeV
    mD_eV = mD_GeV * 1.0e9
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light_target) @ mD_eV)
    MR_GeV = MR_eV / 1.0e9

    m_light_reco = -(mD_eV @ np.linalg.inv(MR_eV) @ mD_eV.T)
    light_masses, u_nu_reco = seesaw.takagi_by_h(m_light_reco)
    u_pmns_reco = u_e.conjugate().T @ u_nu_reco
    angles_reco = seesaw.mixing_angles(u_pmns_reco)
    angle_residual = {
        key: float(angles_reco[key] - getattr(benchmark, key))
        for key in ["sin2_theta12", "sin2_theta13", "sin2_theta23"]
    }

    heavy_masses_GeV = np.sort(np.linalg.svd(MR_GeV, compute_uv=False))
    theta_norm = float(np.linalg.norm(mD_eV @ np.linalg.inv(MR_eV), ord=2))
    seesaw_residual = float(np.linalg.norm(m_light_reco - m_light_target) / np.linalg.norm(m_light_target))
    mass_residual = float(np.linalg.norm(np.sort(light_masses) - m_light_diag) / np.linalg.norm(m_light_diag))
    MR_scale = float(np.linalg.svd(MR_GeV, compute_uv=False)[-1])
    MR_normalized = MR_GeV / MR_scale
    scalar_fit = seesaw.veronese_kernel_fit(MR_normalized)
    lifted_fit = seesaw.trace_lift_kernel_fit(MR_normalized)

    checks = {
        "seesaw_matrix_residual_lt_1e-10": seesaw_residual < 1.0e-10,
        "light_mass_residual_lt_1e-10": mass_residual < 1.0e-10,
        "theta_norm_lt_1e-10": theta_norm < 1.0e-10,
        "pmns_angle_max_residual_lt_1e-10": max(abs(x) for x in angle_residual.values()) < 1.0e-10,
        "MR_symmetric": bool(np.linalg.norm(MR_GeV - MR_GeV.T) / np.linalg.norm(MR_GeV) < 1.0e-12),
        "trace_lift_residual_lt_1e-10": lifted_fit["relative_residual"] < 1.0e-10,
    }

    return {
        "basis_convention": {
            "family_basis": "normalized CP1 O(2) holomorphic basis psi_0,psi_1,psi_2",
            "charged_lepton_left_rotation": "U_e diagonalizes Y_e Y_e^dagger with ascending singular values",
            "pmns_definition": "U_PMNS = U_e^dagger U_nu",
            "majorana_convention": "m_light = U_nu^* diag(m_i) U_nu^dagger",
        },
        "benchmark": asdict(benchmark),
        "mD_largest_GeV": mD_largest_GeV,
        "Y_nu_singular_values_normalized": [float(x) for x in np.linalg.svd(y_nu, compute_uv=False)],
        "Y_e_singular_values_normalized": [float(x) for x in charged_singulars],
        "U_e": matrix_json(u_e),
        "U_PMNS_target": matrix_json(u_pmns_target),
        "U_PMNS_reconstructed": matrix_json(u_pmns_reco),
        "mD_GeV": matrix_json(mD_GeV),
        "m_light_target_eV": matrix_json(m_light_target),
        "m_light_reconstructed_eV": matrix_json(m_light_reco),
        "target_light_masses_eV": [float(x) for x in m_light_diag],
        "reconstructed_light_masses_eV": [float(x) for x in np.sort(light_masses)],
        "reconstructed_pmns_angles": angles_reco,
        "pmns_angle_residuals": angle_residual,
        "MR_GeV": matrix_json(MR_GeV),
        "MR_normalized_by_smallest_singular": matrix_json(MR_normalized),
        "MR_normalized_markdown": matrix_markdown(MR_normalized),
        "heavy_neutrino_masses_GeV": [float(x) for x in heavy_masses_GeV],
        "MR_condition_number": float(heavy_masses_GeV[-1] / heavy_masses_GeV[0]),
        "theta_norm": theta_norm,
        "seesaw_matrix_residual": seesaw_residual,
        "light_mass_residual": mass_residual,
        "scalar_veronese_kernel_fit": scalar_fit,
        "trace_lift_kernel_fit": lifted_fit,
        "checks": checks,
        "all_checks_ok": all(checks.values()),
    }


def write_dual_csv(sectors: dict[str, Any]) -> None:
    rows = []
    for sector, item in sectors.items():
        for idx, value in enumerate(item["dual_coefficients_normalized"]):
            rows.append(
                {
                    "sector": sector,
                    "m": idx,
                    "a_m_normalized_re": value["re"],
                    "a_m_normalized_im": value["im"],
                    "a_m_normalized_abs": math.hypot(value["re"], value["im"]),
                    "largest_singular_value_raw": item["normalization"]["largest_singular_value_raw"],
                }
            )
    with (OUT / "dual_density_coefficients.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Flavor benchmark card")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Covariant Yukawa data")
    lines.append("")
    lines.append("The superpotential branch is")
    lines.append("")
    lines.append("```text")
    lines.append("Y_ij^(a) = <H_a, psi_i psi_j> = sum_m a_m^(a) C_ij^m")
    lines.append("H_a in Gamma_sm(O(-4)) tensor Omega^(1,1)")
    lines.append("```")
    lines.append("")
    lines.append("| sector | small | mid | direct-dual error | Veronese error |")
    lines.append("|---|---:|---:|---:|---:|")
    for sector, item in payload["yukawa_sectors"].items():
        norm = item["normalized_singular_values"]
        checks = item["checks"]
        lines.append(
            f"| {sector} | {norm[2]:.6e} | {norm[1]:.6e} | "
            f"{checks['direct_vs_dual_normalized_max_abs']:.3e} | "
            f"{checks['veronese_relation_abs_normalized']:.3e} |"
        )
    lines.append("")
    lines.append("## Seesaw reproducibility")
    lines.append("")
    s = payload["seesaw"]
    heavy = s["heavy_neutrino_masses_GeV"]
    light = s["reconstructed_light_masses_eV"]
    angles = s["reconstructed_pmns_angles"]
    lines.append(
        "Heavy Majorana masses: "
        f"`({heavy[0]:.6e}, {heavy[1]:.6e}, {heavy[2]:.6e}) GeV`."
    )
    lines.append(
        "Light masses: "
        f"`({light[0]:.6e}, {light[1]:.6e}, {light[2]:.6e}) eV`."
    )
    lines.append(
        "Reconstructed PMNS sin^2 angles: "
        f"`({angles['sin2_theta12']:.6f}, {angles['sin2_theta13']:.6f}, "
        f"{angles['sin2_theta23']:.6f})`."
    )
    lines.append("")
    lines.append("Checks:")
    lines.append("")
    for key, value in s["checks"].items():
        lines.append(f"- `{key}` = `{value}`")
    lines.append("")
    lines.append(f"`theta_norm = {s['theta_norm']:.6e}`.")
    lines.append(f"`seesaw_matrix_residual = {s['seesaw_matrix_residual']:.6e}`.")
    lines.append(f"`trace_lift_relative_residual = {s['trace_lift_kernel_fit']['relative_residual']:.6e}`.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `output/flavor_benchmark/flavor_benchmark_card.json`")
    lines.append("- `output/flavor_benchmark/dual_density_coefficients.csv`")
    lines.append("")
    (OUT / "flavor_benchmark_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    legacy.require_verified_fields()
    OUT.mkdir(parents=True, exist_ok=True)
    grid = legacy.CP1O2Grid(n_theta=54, n_phi=108)
    chi = cov.o4_basis(grid)
    pmap = cov.product_map(grid, chi)
    coeff = pmap["coefficients"]
    payload_in = load_legacy_payload()

    sectors = {}
    matrices = {}
    for name, item in payload_in["results"].items():
        sectors[name], matrices[name] = export_sector(grid, chi, coeff, name, item)

    seesaw_payload = seesaw_card_from_yukawas(matrices["neutrino_dirac"], matrices["charged_lepton"])
    bundle_checks = {
        "O2_gram_max_error": float(np.max(np.abs(grid.gram() - np.eye(3)))),
        "O4_gram_max_error": float(np.max(np.abs((chi.conjugate() * grid.weight) @ chi.T - np.eye(5)))),
        "product_reconstruction_max_error": pmap["max_product_reconstruction_error"],
        "veronese_relation_max_abs": pmap["veronese_relation_max_abs"],
    }

    payload = {
        "note": "No web lookup used. Action-level CP1/O(2) flavor benchmark card.",
        "input": {
            "legacy_yukawa_json": str(LEGACY_JSON),
            "grid": payload_in["grid"],
        },
        "basis": {
            "psi": payload_in["basis"],
            "chi_m": "sqrt(5*binom(4,m))*exp(i*m*phi)*sin(theta/2)^m*cos(theta/2)^(4-m)",
        },
        "bundle_checks": bundle_checks,
        "product_clebsch_C_ij_m": tensor3_json(coeff),
        "yukawa_sectors": sectors,
        "seesaw": seesaw_payload,
        "all_checks_ok": bool(
            seesaw_payload["all_checks_ok"]
            and all(
                item["checks"]["direct_vs_dual_normalized_max_abs"] < 1.0e-12
                for item in sectors.values()
            )
        ),
    }

    (OUT / "flavor_benchmark_card.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_dual_csv(sectors)
    write_report(payload)

    print("Flavor benchmark card export")
    print(f"  O2 gram error: {bundle_checks['O2_gram_max_error']:.3e}")
    print(f"  O4 gram error: {bundle_checks['O4_gram_max_error']:.3e}")
    for name, item in sectors.items():
        norm = item["normalized_singular_values"]
        err = item["checks"]["direct_vs_dual_normalized_max_abs"]
        print(f"  {name:15s}: small={norm[2]:.6e} mid={norm[1]:.6e} dual_err={err:.3e}")
    heavy = seesaw_payload["heavy_neutrino_masses_GeV"]
    print("  heavy masses GeV: " + ", ".join(f"{x:.6e}" for x in heavy))
    print(f"  seesaw residual: {seesaw_payload['seesaw_matrix_residual']:.3e}")
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")
    if not payload["all_checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
