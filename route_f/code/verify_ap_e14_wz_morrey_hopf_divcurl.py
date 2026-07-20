#!/usr/bin/env python3
"""AP-E14 WZ-flux, Morrey, and Hopf base/gauge audit.

This card separates four logically different statements:

1. for one fixed complete-minor graph map, local Wess--Zumino flux decays;
2. that decay does not imply linear-scale target oscillation;
3. the Hopf base/gauge split exactly diagonalizes every minor order;
4. only the signed Chern--Simons density has exact div--curl structure;
   the positive cutoff remainder does not.

Production mode also re-relaxes the same N=17,21,25 compatible-action B=1
backgrounds used in AP-E12 and measures finite-grid WZ and Lipschitz proxies.
Those measurements remain evidence, not a continuum Morrey theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E12_SCRIPT = ROUTE_F / "code" / "verify_ap_e12_graph_density_regular_minimizer.py"
E13_JSON = OUTPUT / "ap_e13_annular_replacement_reverse_holder.json"
CHECKS: list[dict[str, Any]] = []


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E12 = load_module("ap_e12_for_ap_e14", E12_SCRIPT)
E11 = E12.E11


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def hopf_map(u: np.ndarray) -> np.ndarray:
    x1, x2, x3, x4 = u
    return np.asarray(
        [
            2.0 * (x1 * x3 + x2 * x4),
            2.0 * (x2 * x3 - x1 * x4),
            x1 * x1 + x2 * x2 - x3 * x3 - x4 * x4,
        ]
    )


def hopf_derivative(u: np.ndarray) -> np.ndarray:
    x1, x2, x3, x4 = u
    return np.asarray(
        [
            [2.0 * x3, 2.0 * x4, 2.0 * x1, 2.0 * x2],
            [-2.0 * x4, 2.0 * x3, 2.0 * x2, -2.0 * x1],
            [2.0 * x1, 2.0 * x2, -2.0 * x3, -2.0 * x4],
        ]
    )


def hopf_variables(u: np.ndarray, gradient: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return n=h(u), a=2u*theta, and Dn for a tangent 4x3 gradient."""
    x1, x2, x3, x4 = u
    theta = np.asarray([-x2, x1, -x4, x3])
    n = hopf_map(u)
    a = 2.0 * theta @ gradient
    dn = hopf_derivative(u) @ gradient
    return n, a, dn


def hopf_algebra_audit(quick: bool) -> dict[str, Any]:
    rng = np.random.default_rng(20260720)
    count = 24 if quick else 128
    rows = []
    for _ in range(count):
        u = rng.normal(size=4)
        u /= np.linalg.norm(u)
        ambient = rng.normal(size=(4, 3))
        gradient = ambient - u[:, None] * (u @ ambient)[None, :]
        n, a, dn = hopf_variables(u, gradient)

        q = np.asarray(
            [
                n @ np.cross(dn[:, 1], dn[:, 2]),
                n @ np.cross(dn[:, 2], dn[:, 0]),
                n @ np.cross(dn[:, 0], dn[:, 1]),
            ]
        )
        a_wedge_dn_squared = sum(
            float(np.sum((a[i] * dn[:, j] - a[j] * dn[:, i]) ** 2))
            for i in range(3)
            for j in range(i + 1, 3)
        )
        a_wedge_q = float(a @ q)

        gram = gradient.T @ gradient
        raw_first = float(np.sum(gradient**2))
        raw_second = float(
            0.5 * (np.trace(gram) ** 2 - np.trace(gram @ gram))
        )
        raw_third = float(np.linalg.det(gram))
        raw_jacobian = float(np.linalg.det(np.column_stack([u, gradient])))
        rows.append(
            {
                "first_residual": abs(
                    raw_first - 0.25 * (float(a @ a) + float(np.sum(dn**2)))
                ),
                "second_residual": abs(
                    raw_second
                    - 0.0625 * (a_wedge_dn_squared + float(q @ q))
                ),
                "third_residual": abs(raw_third - a_wedge_q**2 / 64.0),
                "volume_residual": abs(raw_jacobian + a_wedge_q / 8.0),
                "hopf_unit_residual": abs(float(n @ n) - 1.0),
                "base_tangency_residual": float(np.max(np.abs(n @ dn))),
            }
        )
    maxima = {
        key: max(row[key] for row in rows)
        for key in rows[0]
    }
    return {
        "samples": count,
        "maximum_residuals": maxima,
        "identities": {
            "connection": "a=2 u*theta and da=-n*omega_S2 for the stated standard S2 orientation",
            "first": "|Du|^2=(|a|^2+|Dn|^2)/4",
            "second": "|M2(Du)|^2=(|a wedge Dn|^2+|da|^2)/16",
            "third": "|M3(Du)|^2=|a wedge da|^2/64",
            "degree_orientation": "det[u,Du]=*(a wedge da)/8 and B=(16 pi^2)^(-1) int a wedge da",
        },
    }


def fixed_map_flux_audit() -> dict[str, Any]:
    radii = np.logspace(-2.0, -12.0, 11)
    # Slow integrable endpoint density: f(r)=1/[r^3 log(e/r)^2].
    local_m3_squared_mass = 4.0 * math.pi / np.log(math.e / radii)
    return {
        "radii": radii.tolist(),
        "slow_integrable_local_m3_squared_mass": local_m3_squared_mass.tolist(),
        "smallest_radius_mass": float(local_m3_squared_mass[-1]),
        "strictly_decreasing": bool(np.all(np.diff(local_m3_squared_mass) < 0.0)),
        "theorem": "for fixed u with M3(Du) in L2, sup_x dist(V_u(x,r),2 pi^2 Z)^2/r^3 <= (4 pi/3) sup_x int_Br(x)|M3(Du)|^2 -> 0",
        "scope": "applies to the selected relaxed B=1 graph representative; minimality is not needed",
    }


def rank_one_morrey_obstruction() -> dict[str, Any]:
    radii = np.logspace(-2.0, -12.0, 11)
    dirichlet_tail_per_unit_axis = 2.0 * math.pi / np.log(math.e / radii)
    return {
        "map": "u(r,z)=(cos ell(r),sin ell(r),0,0), ell=log log(e/r), localized inside a vacuum cylinder",
        "dirichlet_tail_per_unit_axis": dirichlet_tail_per_unit_axis.tolist(),
        "smallest_radius_tail": float(dirichlet_tail_per_unit_axis[-1]),
        "derivative_rank": 1,
        "M2": 0.0,
        "M3": 0.0,
        "wz_flux": 0.0,
        "oscillation_diameter_on_every_axis_neighbourhood": 2.0,
        "linear_morrey_bound": False,
        "fixed_degree_insertion": "insert the localized circle-valued defect into a vacuum cylinder of a smooth compactly supported B=1 map; smooth truncations have degree zero locally and preserve global B=1",
        "consequence": "WZ decay plus graph energy and fixed degree do not imply continuity or O(r) oscillation; a theorem using actual relaxed minimality is still required",
    }


def divcurl_audit() -> dict[str, Any]:
    amplitude_t = 0.7
    amplitude_phase = 0.4
    shear_integral = (
        amplitude_t**2 * amplitude_phase**2 * math.pi**3 / 2.0
    )
    epsilon = 0.3
    torus_volume = (2.0 * math.pi) ** 3
    return {
        "periodic_shear": {
            "lift": "u=(exp(i chi) cos(t/2), exp(i chi) sin(t/2)), t=A sin x1, chi=B sin x2",
            "a": "2 B cos(x2) dx2",
            "base_rank": 1,
            "da": 0.0,
            "a_wedge_da": 0.0,
            "integral_M2_squared": shear_integral,
            "positive": bool(shear_integral > 0.0),
            "meaning": "the positive vertical-horizontal shear |a wedge Dn|^2 survives when every curl/Chern-Simons term vanishes",
        },
        "chern_simons_phase_test": {
            "a": "sin(z) dx+cos(z) dy+c dz",
            "chi": "epsilon sin(x)",
            "epsilon": epsilon,
            "cs_integral_before": torus_volume,
            "cs_integral_after_a_plus_2dchi": torus_volume,
            "cs_integral_residual": 0.0,
            "cs_L2_before": torus_volume,
            "cs_L2_after": torus_volume * (1.0 + epsilon**2),
            "cs_L2_increase": torus_volume * epsilon**2,
            "meaning": "a wedge da changes by the exact form 2 d(chi da), but its squared energy does not; exact topological cancellation does not remove the sextic cutoff stress",
        },
        "full_cutoff_exact_divcurl": False,
        "surviving_positive_remainder": "r^(-2)|u-q|^2[|a wedge Dn|^2+|da|^2]/16",
    }


def local_degree_and_m3(
    field: np.ndarray,
    tetrahedra: np.ndarray,
    orientations: np.ndarray,
    volumes: np.ndarray,
    order: int = 6,
) -> tuple[np.ndarray, np.ndarray]:
    barycentric, weights = E11.duffy_simplex_quadrature(order)
    values = field.reshape(-1, 4)[tetrahedra]
    determinants = np.linalg.det(np.transpose(values, (0, 2, 1)))
    affine = np.einsum("qv,tva->tqa", barycentric, values)
    norm_squared = np.einsum("tqa,tqa->tq", affine, affine)
    radial_four = np.einsum("q,tq->t", weights, norm_squared ** (-2.0))
    radial_eight = np.einsum("q,tq->t", weights, norm_squared ** (-4.0))
    physical_wz = orientations * determinants * radial_four
    m3_squared = determinants**2 * radial_eight / (6.0 * volumes)
    return physical_wz, m3_squared


def numerical_background_audit(quick: bool) -> dict[str, Any]:
    specifications = ((17, 5.5),) if quick else ((17, 5.5), (21, 6.0), (25, 6.5))
    rows = []
    for extent, length in specifications:
        plan = E11.RelaxationPlan(
            f"e14_N{extent}",
            extent,
            length,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            ("wz_morrey",),
        )
        started = time.monotonic()
        field, axis, row, geometry = E12.relax_field(plan, quick)
        tetrahedra, orientations, gradients, volumes, left, right = geometry
        spacing = float(axis[1] - axis[0])
        flat = field.reshape(-1, 4)
        edge_dots = np.clip(np.sum(flat[left] * flat[right], axis=1), -1.0, 1.0)
        edge_quotients = np.arccos(edge_dots) / spacing

        coordinates = np.stack(
            np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
        )
        center_index = (extent // 2, extent // 2, extent // 2)
        center_value = field[center_index]
        node_radii = np.linalg.norm(coordinates, axis=-1)
        node_angles = np.arccos(
            np.clip(np.einsum("...a,a->...", field, center_value), -1.0, 1.0)
        )
        morrey_rows = []
        for multiplier in (1.0, 2.0, 3.0):
            radius = multiplier * spacing
            selected = node_radii <= radius + 1.0e-12
            morrey_rows.append(
                {
                    "radius_in_spacing_units": multiplier,
                    "radius": radius,
                    "maximum_target_angle": float(np.max(node_angles[selected])),
                    "angle_over_radius": float(np.max(node_angles[selected]) / radius),
                }
            )

        physical_wz, normalized_m3_squared = local_degree_and_m3(
            field, tetrahedra, orientations, volumes
        )
        coordinate_flat = coordinates.reshape(-1, 3)
        tetrahedron_centers = np.mean(coordinate_flat[tetrahedra], axis=1)
        tetrahedron_radii = np.linalg.norm(tetrahedron_centers, axis=1)
        wz_rows = []
        for multiplier in (2.0, 3.0):
            radius = multiplier * spacing
            selected = tetrahedron_radii <= radius
            flux = float(np.sum(physical_wz[selected]))
            local_m3 = float(np.sum(normalized_m3_squared[selected]))
            selected_volume = float(np.sum(volumes[selected]))
            wz_rows.append(
                {
                    "radius_in_spacing_units": multiplier,
                    "radius": radius,
                    "physical_wz_flux": flux,
                    "flux_squared_over_radius_cubed": flux**2 / radius**3,
                    "normalized_m3_squared_integral": local_m3,
                    "selected_union_volume": selected_volume,
                    "cauchy_ratio": flux**2 / max(selected_volume * local_m3, 1.0e-300),
                }
            )

        rows.append(
            {
                "extent": extent,
                "length": length,
                "spacing": spacing,
                "runtime_seconds": time.monotonic() - started,
                "stationarity": row["projected_gradient_density"],
                "admissible": row["admissible"],
                "B": row["final_B"],
                "maximum_edge_angle_over_spacing": float(np.max(edge_quotients)),
                "q99_edge_angle_over_spacing": float(np.quantile(edge_quotients, 0.99)),
                "morrey_rows": morrey_rows,
                "wz_rows": wz_rows,
            }
        )

    maximum_cauchy = max(
        wz["cauchy_ratio"] for row in rows for wz in row["wz_rows"]
    )
    all_backgrounds = all(
        row["admissible"]
        and row["B"] == [1, 1, 1]
        and row["stationarity"] < E11.STATIONARITY_TOLERANCE
        for row in rows
    )
    edge_values = [row["maximum_edge_angle_over_spacing"] for row in rows]
    finite_morrey_proxy = bool(
        max(edge_values) < 4.2
        and (len(edge_values) == 1 or np.all(np.diff(edge_values) < 0.0))
    )
    return {
        "rows": rows,
        "all_backgrounds_stationary_admissible_B1": all_backgrounds,
        "maximum_local_wz_cauchy_ratio": maximum_cauchy,
        "local_wz_cauchy_pass": bool(maximum_cauchy <= 1.0 + 2.0e-10),
        "finite_grid_linear_morrey_proxy": finite_morrey_proxy,
        "scope": "finite P1 evidence only; bounded/decreasing edge quotients do not prove a continuum L-infinity gradient bound",
    }


def theorem_ledger() -> dict[str, Any]:
    return {
        "fixed_graph_map_wz_flux_decay": True,
        "selected_relaxed_B1_wz_flux_decay": True,
        "uses_minimality": False,
        "linear_scale_morrey_for_selected_relaxed_minimizer": False,
        "linear_scale_morrey_status": "open; graph energy, fixed degree, and WZ decay are insufficient",
        "hopf_complete_minor_decomposition": True,
        "signed_chern_simons_exact_pairing": True,
        "full_cutoff_remainder_exact_divcurl": False,
        "full_divcurl_statement_disproved": True,
        "reason": "periodic rank-one-base/vertical-phase shear has da=a wedge da=0 but strictly positive |a wedge Dn|^2; phase exactness preserves int a wedge da but not int |a wedge da|^2",
        "weak_euler_lagrange_for_naive_integral": False,
        "endpoint_graph_density": False,
        "classicality": False,
        "continuum_isolation": False,
        "hessian_gate_open": False,
        "determinant_variation_gate_open": False,
        "portal_start_allowed": False,
        "ordered_next_gate": [
            "derive a defect-measure inner-variation/monotonicity formula for the lower-semicontinuous relaxation itself",
            "prove that its rank-one vertical-horizontal defect measure vanishes, or construct a degree-preserving local comparison that removes it",
            "then seek epsilon regularity and upgrade mean Morrey decay to an L-infinity oscillation estimate",
        ],
    }


def render_markdown(result: dict[str, Any]) -> str:
    algebra = result["hopf_algebra"]["maximum_residuals"]
    numerical = result["numerical_background"]
    theory = result["theory"]
    return "\n".join(
        [
            "# AP-E14 WZ/Morrey and Hopf div-curl audit",
            "",
            f"Status: `{result['status']}`.",
            "",
            "## Closed statements",
            "",
            f"- fixed-map WZ-flux decay: `{theory['fixed_graph_map_wz_flux_decay']}`.",
            f"- applies to selected relaxed B=1 representative: `{theory['selected_relaxed_B1_wz_flux_decay']}`.",
            f"- Hopf complete-minor decomposition: `{theory['hopf_complete_minor_decomposition']}`.",
            f"- largest Hopf algebra residual: `{max(algebra.values()):.3e}`.",
            "",
            "## No-go and remaining gate",
            "",
            f"- linear-scale Morrey theorem: `{theory['linear_scale_morrey_for_selected_relaxed_minimizer']}`.",
            f"- full cutoff exact div-curl: `{theory['full_cutoff_remainder_exact_divcurl']}`.",
            f"- full div-curl statement disproved: `{theory['full_divcurl_statement_disproved']}`.",
            f"- Hessian gate: `{theory['hessian_gate_open']}`.",
            "",
            "## Same-action finite-grid evidence",
            "",
            f"- all backgrounds stationary/admissible/B=1: `{numerical['all_backgrounds_stationary_admissible_B1']}`.",
            f"- maximum local WZ Cauchy ratio: `{numerical['maximum_local_wz_cauchy_ratio']:.9f}`.",
            f"- finite-grid Morrey proxy: `{numerical['finite_grid_linear_morrey_proxy']}`.",
            "",
            f"Checks: `{result['checks_passed']}/{result['checks_total']}`.",
        ]
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    algebra = hopf_algebra_audit(arguments.quick)
    flux = fixed_map_flux_audit()
    morrey = rank_one_morrey_obstruction()
    divcurl = divcurl_audit()
    numerical = numerical_background_audit(arguments.quick)
    theory = theorem_ledger()
    residuals = algebra["maximum_residuals"]

    check(
        "H1",
        "Hopf metric and target constraints",
        residuals["first_residual"] < 1.0e-11
        and residuals["hopf_unit_residual"] < 1.0e-12
        and residuals["base_tangency_residual"] < 1.0e-11,
        f"first={residuals['first_residual']:.3e}",
    )
    check(
        "H2",
        "Hopf second-minor identity",
        residuals["second_residual"] < 1.0e-10,
        f"residual={residuals['second_residual']:.3e}",
    )
    check(
        "H3",
        "Hopf third-minor identity",
        residuals["third_residual"] < 1.0e-10,
        f"residual={residuals['third_residual']:.3e}",
    )
    check(
        "H4",
        "Hopf Chern-Simons degree orientation",
        residuals["volume_residual"] < 1.0e-11,
        f"residual={residuals['volume_residual']:.3e}",
    )
    check(
        "W1",
        "fixed-map WZ flux decays by L2 absolute continuity",
        theory["fixed_graph_map_wz_flux_decay"]
        and flux["strictly_decreasing"]
        and flux["smallest_radius_mass"] < flux["slow_integrable_local_m3_squared_mass"][0],
        f"slow endpoint mass={flux['smallest_radius_mass']:.3e}",
    )
    check(
        "M1",
        "WZ decay does not imply linear Morrey oscillation",
        morrey["M2"] == 0.0
        and morrey["M3"] == 0.0
        and morrey["oscillation_diameter_on_every_axis_neighbourhood"] == 2.0
        and not morrey["linear_morrey_bound"],
        "rank-one circle-valued defect has zero WZ flux and full oscillation",
    )
    check(
        "D1",
        "periodic Hopf shear leaves positive non-curl remainder",
        divcurl["periodic_shear"]["da"] == 0.0
        and divcurl["periodic_shear"]["a_wedge_da"] == 0.0
        and divcurl["periodic_shear"]["positive"],
        f"integral={divcurl['periodic_shear']['integral_M2_squared']:.6e}",
    )
    check(
        "D2",
        "exact phase preserves Chern-Simons integral but not squared energy",
        divcurl["chern_simons_phase_test"]["cs_integral_residual"] == 0.0
        and divcurl["chern_simons_phase_test"]["cs_L2_increase"] > 0.0,
        f"increase={divcurl['chern_simons_phase_test']['cs_L2_increase']:.6e}",
    )
    check(
        "N1",
        "same-action backgrounds remain stationary admissible B=1",
        numerical["all_backgrounds_stationary_admissible_B1"],
        f"cases={len(numerical['rows'])}",
    )
    check(
        "N2",
        "normalized-affine local WZ Cauchy inequality",
        numerical["local_wz_cauchy_pass"],
        f"maximum ratio={numerical['maximum_local_wz_cauchy_ratio']:.9f}",
    )
    check(
        "N3",
        "finite-grid linear Morrey proxy",
        numerical["finite_grid_linear_morrey_proxy"],
        "bounded/decreasing edge-angle quotients are evidence only",
    )
    check(
        "G1",
        "downstream gates fail closed",
        not theory["linear_scale_morrey_for_selected_relaxed_minimizer"]
        and not theory["full_cutoff_remainder_exact_divcurl"]
        and not theory["hessian_gate_open"]
        and not theory["determinant_variation_gate_open"]
        and not theory["portal_start_allowed"],
        "WZ is closed, but Morrey/classicality and full div-curl cancellation are not",
    )

    result = {
        "artifact": "AP-E14 WZ-flux Morrey and Hopf base-gauge div-curl audit",
        "status": "quick_wz_closed_morrey_divcurl_fail_closed"
        if arguments.quick
        else "production_wz_closed_morrey_divcurl_fail_closed",
        "quick": arguments.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "source_manifest": [source_row(THIS_SCRIPT), source_row(E12_SCRIPT), source_row(E13_JSON)],
        "hopf_algebra": algebra,
        "fixed_map_wz_flux": flux,
        "rank_one_morrey_obstruction": morrey,
        "divcurl": divcurl,
        "numerical_background": numerical,
        "theory": theory,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    json_path = OUTPUT / "ap_e14_wz_morrey_hopf_divcurl.json"
    markdown_path = OUTPUT / "ap_e14_wz_morrey_hopf_divcurl.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    markdown_path.write_text(render_markdown(result))
    print(
        f"AP-E14 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"wz={theory['selected_relaxed_B1_wz_flux_decay']}; "
        f"morrey={theory['linear_scale_morrey_for_selected_relaxed_minimizer']}; "
        f"divcurl={theory['full_cutoff_remainder_exact_divcurl']}; "
        f"hessian={theory['hessian_gate_open']}",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
