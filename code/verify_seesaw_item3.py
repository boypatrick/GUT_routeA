#!/usr/bin/env python3
"""Item 3: type-I seesaw verification using the item-2 nu_D texture.

The script reads the CP1/O(2) Yukawa textures from item 2, reconstructs the
Majorana matrix required for a normal-ordering light-neutrino benchmark, verifies
the seesaw relation, and checks that the resulting complex symmetric M_R can be
represented by a finite CP1/O(2) Veronese moment kernel.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
YUKAWA_JSON = ROOT / "output" / "yukawa_o2" / "cp1_o2_yukawa_scan.json"
OUT = ROOT / "output" / "seesaw"


@dataclass(frozen=True)
class LightNuBenchmark:
    m1_eV: float
    dm21_eV2: float
    dm31_eV2: float
    sin2_theta12: float
    sin2_theta13: float
    sin2_theta23: float
    delta_cp_rad: float
    alpha21_rad: float
    alpha31_rad: float


def load_matrix(sector: str) -> np.ndarray:
    payload = json.loads(YUKAWA_JSON.read_text(encoding="utf-8"))
    raw = payload["results"][sector]["matrix"]
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def standard_pmns(
    sin2_theta12: float,
    sin2_theta13: float,
    sin2_theta23: float,
    delta: float,
    alpha21: float,
    alpha31: float,
) -> np.ndarray:
    s12, s13, s23 = math.sqrt(sin2_theta12), math.sqrt(sin2_theta13), math.sqrt(sin2_theta23)
    c12, c13, c23 = math.sqrt(1.0 - sin2_theta12), math.sqrt(1.0 - sin2_theta13), math.sqrt(1.0 - sin2_theta23)
    eid = np.exp(1j * delta)
    emid = np.exp(-1j * delta)
    u = np.array(
        [
            [c12 * c13, s12 * c13, s13 * emid],
            [-s12 * c23 - c12 * s23 * s13 * eid, c12 * c23 - s12 * s23 * s13 * eid, s23 * c13],
            [s12 * s23 - c12 * c23 * s13 * eid, -c12 * s23 - s12 * c23 * s13 * eid, c23 * c13],
        ],
        dtype=complex,
    )
    majorana = np.diag([1.0, np.exp(0.5j * alpha21), np.exp(0.5j * alpha31)])
    return u @ majorana


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def takagi_by_h(m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = m.conjugate().T @ m
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    masses = np.sqrt(np.maximum(values, 0.0))
    return masses, vectors


def mixing_angles(u: np.ndarray) -> dict[str, float]:
    uabs2 = np.abs(u) ** 2
    s13 = float(uabs2[0, 2])
    s12 = float(uabs2[0, 1] / max(1.0 - s13, 1e-30))
    s23 = float(uabs2[1, 2] / max(1.0 - s13, 1e-30))
    return {"sin2_theta12": s12, "sin2_theta13": s13, "sin2_theta23": s23}


def cp1_o2_eval(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [
            math.sqrt(3.0) * math.cos(theta / 2.0) ** 2,
            math.sqrt(6.0) * np.exp(1j * phi) * math.sin(theta / 2.0) * math.cos(theta / 2.0),
            math.sqrt(3.0) * np.exp(2j * phi) * math.sin(theta / 2.0) ** 2,
        ],
        dtype=complex,
    )


def sym_entries(m: np.ndarray) -> np.ndarray:
    return np.array([m[0, 0], m[0, 1], m[0, 2], m[1, 1], m[1, 2], m[2, 2]], dtype=complex)


def unsym_entries(entries: np.ndarray) -> np.ndarray:
    out = np.zeros((3, 3), dtype=complex)
    positions = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    for value, (i, j) in zip(entries, positions):
        out[i, j] = value
        out[j, i] = value
    return out


def veronese_kernel_fit(m: np.ndarray) -> dict[str, object]:
    # Six generic CP1 points.  Their Veronese images span Sym^2(C^3).
    points = [
        (0.47, 0.10),
        (0.91, 1.17),
        (1.23, 2.41),
        (1.72, 3.04),
        (2.21, 4.36),
        (2.63, 5.51),
    ]
    columns = []
    for theta, phi in points:
        v = cp1_o2_eval(theta, phi)
        columns.append(sym_entries(np.outer(v, v)))
    basis = np.stack(columns, axis=1)
    target = sym_entries(m)
    weights, *_ = np.linalg.lstsq(basis, target, rcond=None)
    reconstructed = unsym_entries(basis @ weights)
    residual = np.linalg.norm(reconstructed - m) / max(np.linalg.norm(m), 1e-30)
    condition = np.linalg.cond(basis)
    rows = []
    for (theta, phi), weight in zip(points, weights):
        rows.append(
            {
                "theta": float(theta),
                "phi": float(phi),
                "weight_re": float(weight.real),
                "weight_im": float(weight.imag),
                "abs_weight": float(abs(weight)),
            }
        )
    return {
        "points": rows,
        "basis_condition": float(condition),
        "relative_residual": float(residual),
        "constraint": "scalar O(2) moment matrices obey M_11 = 2 M_02 in the psi_0,psi_1,psi_2 basis",
        "constraint_violation_abs": float(abs(m[1, 1] - 2.0 * m[0, 2])),
    }


def trace_lift_kernel_fit(m: np.ndarray) -> dict[str, object]:
    # The scalar multiplication map Sym^2 H^0(O(2)) -> H^0(O(4)) is five-dimensional:
    # psi_1^2 = 2 psi_0 psi_2.  The missing spin-0 bilinear can be represented by
    # the SU(2)-invariant charge-conjugation/contact matrix below.
    points = [
        (0.47, 0.10),
        (0.91, 1.17),
        (1.23, 2.41),
        (1.72, 3.04),
        (2.21, 4.36),
    ]
    scalar_columns = []
    for theta, phi in points:
        v = cp1_o2_eval(theta, phi)
        scalar_columns.append(sym_entries(np.outer(v, v)))
    trace_lift = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=complex,
    )
    columns = scalar_columns + [sym_entries(trace_lift)]
    basis = np.stack(columns, axis=1)
    target = sym_entries(m)
    weights = np.linalg.solve(basis, target)
    reconstructed = unsym_entries(basis @ weights)
    residual = np.linalg.norm(reconstructed - m) / max(np.linalg.norm(m), 1e-30)
    rows = []
    for (theta, phi), weight in zip(points, weights[:5]):
        rows.append(
            {
                "theta": float(theta),
                "phi": float(phi),
                "weight_re": float(weight.real),
                "weight_im": float(weight.imag),
                "abs_weight": float(abs(weight)),
            }
        )
    zeta = weights[5]
    return {
        "scalar_points": rows,
        "trace_lift_weight": {"re": float(zeta.real), "im": float(zeta.imag), "abs": float(abs(zeta))},
        "basis_condition": float(np.linalg.cond(basis)),
        "relative_residual": float(residual),
        "trace_lift_matrix": complex_matrix_json(trace_lift),
    }


def complex_matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[{"re": float(z.real), "im": float(z.imag)} for z in row] for row in m]


def matrix_markdown(m: np.ndarray) -> str:
    lines = []
    for row in m:
        lines.append("[ " + ", ".join(f"{z.real:+.4e}{z.imag:+.4e}i" for z in row) + " ]")
    return "\n".join(lines)


def write_veronese_csv(fit: dict[str, object]) -> None:
    with (OUT / "majorana_veronese_kernel.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["theta", "phi", "weight_re", "weight_im", "abs_weight"])
        writer.writeheader()
        writer.writerows(fit["points"])


def write_trace_lift_csv(fit: dict[str, object]) -> None:
    with (OUT / "majorana_trace_lift_kernel.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["theta", "phi", "weight_re", "weight_im", "abs_weight"])
        writer.writeheader()
        writer.writerows(fit["scalar_points"])


def write_report(payload: dict[str, object]) -> None:
    b = payload["benchmark"]
    checks = payload["checks"]
    heavy = payload["heavy_neutrino_masses_GeV"]
    light = payload["reconstructed_light_masses_eV"]
    angles = payload["reconstructed_pmns_angles"]
    target = payload["target_pmns_angles"]
    fit = payload["scalar_veronese_kernel_fit"]
    lifted = payload["trace_lift_kernel_fit"]

    lines: list[str] = []
    lines.append("# Item 3 seesaw verification")
    lines.append("")
    lines.append("Input textures are read from item 2:")
    lines.append("")
    lines.append("- `Y_nuD`: `L nu^c H_u/nu`")
    lines.append("- `Y_e`: `L e^c H_d/e`, used for the charged-lepton left rotation")
    lines.append("")
    lines.append("## Derivation")
    lines.append("")
    lines.append("In Pati-Salam language, `nu^c` has `B-L=+1`, so the bilinear `nu^c nu^c` has")
    lines.append("`B-L=+2`. A Majorana mass therefore requires a `B-L=-2` order parameter,")
    lines.append("realized in a Spin(10) `bar{126}_H` channel or by the effective operator")
    lines.append("`(16_i 16_j bar{16}_H bar{16}_H)/Lambda`.")
    lines.append("")
    lines.append("The neutral-fermion mass matrix is")
    lines.append("")
    lines.append("```text")
    lines.append("M_N = [ 0      m_D ]")
    lines.append("      [ m_D^T  M_R ]")
    lines.append("```")
    lines.append("")
    lines.append("For `||m_D M_R^{-1}|| << 1`, Schur complement block diagonalization gives")
    lines.append("")
    lines.append("```text")
    lines.append("m_light = -m_D M_R^{-1} m_D^T + O(m_D^4/M_R^3).")
    lines.append("```")
    lines.append("")
    lines.append("Conversely, if `m_D` and `m_light` are nonsingular, the required Majorana")
    lines.append("matrix is uniquely")
    lines.append("")
    lines.append("```text")
    lines.append("M_R = -m_D^T m_light^{-1} m_D.")
    lines.append("```")
    lines.append("")
    lines.append("## Benchmark")
    lines.append("")
    lines.append("No web lookup was used. The benchmark is an approximate normal-ordering target:")
    lines.append("")
    lines.append("```text")
    lines.append(f"m1 = {b['m1_eV']:.4e} eV")
    lines.append(f"Delta m21^2 = {b['dm21_eV2']:.4e} eV^2")
    lines.append(f"Delta m31^2 = {b['dm31_eV2']:.4e} eV^2")
    lines.append(
        "sin^2(theta12,theta13,theta23) = "
        f"({target['sin2_theta12']:.4f}, {target['sin2_theta13']:.4f}, {target['sin2_theta23']:.4f})"
    )
    lines.append(f"m_D largest singular value = {payload['mD_largest_GeV']:.3g} GeV")
    lines.append("```")
    lines.append("")
    lines.append("## Numerical result")
    lines.append("")
    lines.append("Reconstructed light masses:")
    lines.append("")
    lines.append("```text")
    lines.append(f"m1,m2,m3 = ({light[0]:.6e}, {light[1]:.6e}, {light[2]:.6e}) eV")
    lines.append("```")
    lines.append("")
    lines.append("Reconstructed PMNS angles:")
    lines.append("")
    lines.append("```text")
    lines.append(
        "sin^2(theta12,theta13,theta23) = "
        f"({angles['sin2_theta12']:.6f}, {angles['sin2_theta13']:.6f}, {angles['sin2_theta23']:.6f})"
    )
    lines.append("```")
    lines.append("")
    lines.append("Heavy Majorana Takagi masses:")
    lines.append("")
    lines.append("```text")
    lines.append(f"M1,M2,M3 = ({heavy[0]:.6e}, {heavy[1]:.6e}, {heavy[2]:.6e}) GeV")
    lines.append("```")
    lines.append("")
    lines.append("Dimensionless right-handed mixing parameter:")
    lines.append("")
    lines.append("```text")
    lines.append(f"||Theta||_2 = ||m_D M_R^-1||_2 = {payload['theta_norm']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("Checks:")
    lines.append("")
    for key, value in checks.items():
        lines.append(f"- `{key}` = `{value}`")
    lines.append("")
    lines.append("## Geometric realization of M_R")
    lines.append("")
    lines.append("A pure scalar `CP1/O(2)` kernel is not enough for a general `M_R`. Because")
    lines.append("`psi_1^2 = 2 psi_0 psi_2`, scalar moment matrices obey")
    lines.append("")
    lines.append("```text")
    lines.append("(M_scalar)_11 = 2 (M_scalar)_02.")
    lines.append("```")
    lines.append("")
    lines.append("The reconstructed `M_R` violates this scalar Veronese constraint, which is")
    lines.append("a useful no-go for the minimal Majorana kernel:")
    lines.append("")
    lines.append("```text")
    lines.append(f"scalar basis condition = {fit['basis_condition']:.6e}")
    lines.append(f"scalar relative residual = {fit['relative_residual']:.6e}")
    lines.append(f"|M11 - 2 M02| = {fit['constraint_violation_abs']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("The fix is a Majorana-only trace-lift/contact channel, representing the missing")
    lines.append("spin-0 part in `Sym^2 H^0(O(2)) = spin-2 plus spin-0`:")
    lines.append("")
    lines.append("```text")
    lines.append("C_0 = [[0,0,1],[0,-1,0],[1,0,0]]")
    lines.append("M_R/M_scale = sum_{a=1}^5 w_a v(x_a)v(x_a)^T + zeta C_0.")
    lines.append("```")
    lines.append("")
    lines.append("With this trace-lift channel:")
    lines.append("")
    lines.append("```text")
    lines.append(f"lifted basis condition = {lifted['basis_condition']:.6e}")
    lines.append(f"lifted relative residual = {lifted['relative_residual']:.6e}")
    lines.append(
        "zeta = "
        f"{lifted['trace_lift_weight']['re']:+.6e}{lifted['trace_lift_weight']['im']:+.6e}i"
    )
    lines.append("```")
    lines.append("")
    lines.append("Normalized reconstructed `M_R` matrix:")
    lines.append("")
    lines.append("```text")
    lines.append(payload["MR_normalized_markdown"])
    lines.append("```")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- The seesaw scale is automatically in the intermediate/GUT neighborhood for a top-like `m_D`.")
    lines.append("- The construction is not yet predictive because inverse reconstruction uses light-neutrino targets.")
    lines.append("- The scalar-kernel no-go is a genuine constraint: the Majorana sector needs the trace-lift channel or an equivalent derivative/curvature insertion.")
    lines.append("- The next stricter test is to restrict the lifted `M_R` kernel to a small number of PSLT/Ed centers and scan whether the same benchmark remains reachable.")
    (OUT / "seesaw_item3_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    y_nu = load_matrix("neutrino_dirac")
    y_e = load_matrix("charged_lepton")
    u_e, charged_singulars = left_rotation(y_e)

    benchmark = LightNuBenchmark(
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

    u_pmns_target = standard_pmns(
        benchmark.sin2_theta12,
        benchmark.sin2_theta13,
        benchmark.sin2_theta23,
        benchmark.delta_cp_rad,
        benchmark.alpha21_rad,
        benchmark.alpha31_rad,
    )
    # Gauge/family-basis neutrino Takagi matrix must satisfy U_PMNS = U_e^dagger U_nu.
    u_nu_target = u_e @ u_pmns_target
    m_light_target = u_nu_target.conjugate() @ np.diag(m_light_diag) @ u_nu_target.conjugate().T

    mD_largest_GeV = 100.0
    mD_eV = y_nu * (mD_largest_GeV * 1.0e9)
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light_target) @ mD_eV)
    MR_GeV = MR_eV / 1.0e9

    m_light_reco = -(mD_eV @ np.linalg.inv(MR_eV) @ mD_eV.T)
    light_masses, u_nu_reco = takagi_by_h(m_light_reco)
    u_pmns_reco = u_e.conjugate().T @ u_nu_reco
    angles_reco = mixing_angles(u_pmns_reco)

    heavy_masses_GeV = np.linalg.svd(MR_GeV, compute_uv=False)
    heavy_masses_GeV = np.sort(heavy_masses_GeV)
    theta_norm = float(np.linalg.norm(mD_eV @ np.linalg.inv(MR_eV), ord=2))
    seesaw_residual = float(np.linalg.norm(m_light_reco - m_light_target) / np.linalg.norm(m_light_target))
    mass_residual = float(np.linalg.norm(np.sort(light_masses) - m_light_diag) / np.linalg.norm(m_light_diag))
    angle_residual = {
        key: float(angles_reco[key] - getattr(benchmark, key))
        for key in ["sin2_theta12", "sin2_theta13", "sin2_theta23"]
    }

    MR_scale = float(np.linalg.svd(MR_GeV, compute_uv=False)[-1])
    MR_normalized = MR_GeV / MR_scale
    fit = veronese_kernel_fit(MR_normalized)
    lifted_fit = trace_lift_kernel_fit(MR_normalized)
    write_veronese_csv(fit)
    write_trace_lift_csv(lifted_fit)

    checks = {
        "seesaw_matrix_residual_lt_1e-10": seesaw_residual < 1e-10,
        "light_mass_residual_lt_1e-10": mass_residual < 1e-10,
        "theta_norm_lt_1e-10": theta_norm < 1e-10,
        "scalar_veronese_residual_gt_1e-3": fit["relative_residual"] > 1e-3,
        "trace_lift_residual_lt_1e-10": lifted_fit["relative_residual"] < 1e-10,
        "MR_symmetric": bool(np.linalg.norm(MR_GeV - MR_GeV.T) / np.linalg.norm(MR_GeV) < 1e-12),
    }

    payload: dict[str, object] = {
        "input": {
            "yukawa_json": str(YUKAWA_JSON),
            "Y_nuD_source": "neutrino_dirac",
            "Y_e_source": "charged_lepton",
        },
        "benchmark": {
            "m1_eV": benchmark.m1_eV,
            "dm21_eV2": benchmark.dm21_eV2,
            "dm31_eV2": benchmark.dm31_eV2,
            "delta_cp_rad": benchmark.delta_cp_rad,
            "alpha21_rad": benchmark.alpha21_rad,
            "alpha31_rad": benchmark.alpha31_rad,
        },
        "target_pmns_angles": {
            "sin2_theta12": benchmark.sin2_theta12,
            "sin2_theta13": benchmark.sin2_theta13,
            "sin2_theta23": benchmark.sin2_theta23,
        },
        "mD_largest_GeV": mD_largest_GeV,
        "Y_nu_singular_values_normalized": [float(x) for x in np.linalg.svd(y_nu, compute_uv=False)],
        "Y_e_singular_values_normalized": [float(x) for x in charged_singulars],
        "target_light_masses_eV": [float(x) for x in m_light_diag],
        "reconstructed_light_masses_eV": [float(x) for x in np.sort(light_masses)],
        "reconstructed_pmns_angles": angles_reco,
        "pmns_angle_residuals": angle_residual,
        "heavy_neutrino_masses_GeV": [float(x) for x in heavy_masses_GeV],
        "MR_condition_number": float(heavy_masses_GeV[-1] / heavy_masses_GeV[0]),
        "theta_norm": theta_norm,
        "seesaw_matrix_residual": seesaw_residual,
        "light_mass_residual": mass_residual,
        "MR_GeV": complex_matrix_json(MR_GeV),
        "MR_normalized_by_smallest_singular": complex_matrix_json(MR_normalized),
        "MR_normalized_markdown": matrix_markdown(MR_normalized),
        "scalar_veronese_kernel_fit": fit,
        "trace_lift_kernel_fit": lifted_fit,
        "checks": checks,
        "all_checks_ok": all(checks.values()),
    }
    (OUT / "seesaw_reconstruction.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)

    print("Item 3 seesaw verification")
    print(f"  input: {YUKAWA_JSON}")
    print(f"  mD largest singular value: {mD_largest_GeV:.3g} GeV")
    print(
        "  light masses eV: "
        + ", ".join(f"{x:.6e}" for x in payload["reconstructed_light_masses_eV"])
    )
    print(
        "  PMNS sin^2: "
        + ", ".join(f"{angles_reco[k]:.6f}" for k in ["sin2_theta12", "sin2_theta13", "sin2_theta23"])
    )
    print("  heavy masses GeV: " + ", ".join(f"{x:.6e}" for x in payload["heavy_neutrino_masses_GeV"]))
    print(f"  theta norm: {theta_norm:.6e}")
    print(f"  seesaw residual: {seesaw_residual:.6e}")
    print(f"  scalar Veronese residual: {fit['relative_residual']:.6e}")
    print(f"  trace-lift residual: {lifted_fit['relative_residual']:.6e}")
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")

    if not payload["all_checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
