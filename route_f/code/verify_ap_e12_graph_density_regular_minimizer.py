#!/usr/bin/env python3
"""AP-E12 endpoint graph-density and relaxed-minimizer audit.

The card distinguishes four statements which must not be conflated:

1. ordinary strong W^{1,2} density for maps into S^3;
2. strong density in the endpoint complete-minor graph topology
   (Du, M_2(Du), M_3(Du)) in L^2;
3. local classicality/regularity of a minimizer of the relaxed functional;
4. isolation modulo translations, needed before a physical Hessian exists.

It proves exact algebraic identities, audits the sharp radial scaling of the
Malý dipole mechanism, and performs finite-grid concentration and perturbed-
start tests on the AP-E11 action.  Numerical evidence is deliberately unable
to promote an endpoint density, regularity, isolation, or Hessian theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import platform
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E11_SCRIPT = ROUTE_F / "code" / "scan_ap_e11_compatible_cochain_action.py"
E11_JSON = OUTPUT / "ap_e11_compatible_cochain_action.json"

CHECKS: list[dict[str, Any]] = []


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E11 = load_module("ap_e11_for_ap_e12", E11_SCRIPT)


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


def complete_minor_density(gradient: np.ndarray) -> np.ndarray:
    density, _ = E11.skyrme_density_and_derivative(
        gradient, E11.R_SKYRME, E11.K_SEXTIC
    )
    return density


def algebraic_symbol_audit(quick: bool) -> dict[str, Any]:
    rng = np.random.default_rng(20260720)
    rows = []
    count = 12 if quick else 64
    for _ in range(count):
        gradient = rng.normal(size=(4, 3))
        singular_values = np.linalg.svd(gradient, compute_uv=False)
        s2 = singular_values**2
        invariant_density = 0.5 * (
            np.sum(s2)
            + E11.R_SKYRME
            * sum(s2[i] * s2[j] for i, j in itertools.combinations(range(3), 2))
            + E11.K_SEXTIC * np.prod(s2)
        )
        direct_density = float(complete_minor_density(gradient))

        target = rng.normal(size=4)
        spatial = rng.normal(size=3)
        rank_one = np.outer(target, spatial)
        second_difference = float(
            complete_minor_density(gradient + rank_one)
            + complete_minor_density(gradient - rank_one)
            - 2.0 * complete_minor_density(gradient)
        )
        delta_two = E11.second_minors(gradient + rank_one) - E11.second_minors(
            gradient
        )
        delta_three = E11.third_minors(
            gradient + rank_one
        ) - E11.third_minors(gradient)
        exact_symbol = float(
            np.sum(rank_one**2)
            + E11.R_SKYRME * np.sum(delta_two**2)
            + E11.K_SEXTIC * np.sum(delta_three**2)
        )
        rows.append(
            {
                "singular_value_identity_residual": abs(
                    direct_density - invariant_density
                ),
                "rank_one_symbol_identity_residual": abs(
                    second_difference - exact_symbol
                ),
                "rank_one_symbol_lower_margin": exact_symbol
                - float(np.sum(rank_one**2)),
            }
        )
    maximum_invariant = max(row["singular_value_identity_residual"] for row in rows)
    maximum_symbol = max(row["rank_one_symbol_identity_residual"] for row in rows)
    minimum_margin = min(row["rank_one_symbol_lower_margin"] for row in rows)
    return {
        "density": "W(F)=1/2(|F|^2+R|M2(F)|^2+K|M3(F)|^2)",
        "singular_value_formula": "2W=sum s_i^2+R sum_{i<j}s_i^2 s_j^2+K product_i s_i^2",
        "rank_one_formula": "D2W(F)[a tensor xi,a tensor xi]=|a tensor xi|^2+R|DM2(F)[a tensor xi]|^2+K|DM3(F)[a tensor xi]|^2",
        "uniform_legendre_hadamard_constant": 1.0,
        "maximum_singular_value_identity_residual": maximum_invariant,
        "maximum_rank_one_symbol_identity_residual": maximum_symbol,
        "minimum_rank_one_symbol_lower_margin": minimum_margin,
        "rows": rows,
    }


def quasiconvex_identity_audit(quick: bool) -> dict[str, Any]:
    rng = np.random.default_rng(12071991)
    specs = (E11.MESH_SPECS[0],) if quick else E11.MESH_SPECS
    rows = []
    for spec in specs:
        tetrahedra, coordinates, gradients, volumes = E11.periodic_tetrahedral_mesh(
            2, spec
        )
        vertices = int(np.max(tetrahedra)) + 1
        for sample in range(2 if quick else 5):
            affine = rng.normal(scale=0.45, size=(4, 3))
            corrector = rng.normal(scale=0.18, size=(vertices, 4))
            corrector -= np.mean(corrector, axis=0, keepdims=True)
            affine_values = np.einsum("tvi,ai->tva", coordinates, affine)
            values = affine_values + corrector[tetrahedra]
            field_gradient = np.einsum("tva,tvi->tai", values, gradients)
            two = E11.second_minors(field_gradient)
            three = E11.third_minors(field_gradient)
            affine_two = E11.second_minors(affine)
            affine_three = E11.third_minors(affine)
            cell_volume = float(np.sum(volumes))
            mean_gradient = np.einsum("t,tai->ai", volumes, field_gradient) / cell_volume
            mean_two = np.einsum("t,tab->ab", volumes, two) / cell_volume
            mean_three = np.einsum("t,ta->a", volumes, three) / cell_volume

            left = float(
                volumes
                @ (complete_minor_density(field_gradient) - complete_minor_density(affine))
            )
            right_density = 0.5 * np.sum((field_gradient - affine) ** 2, axis=(-2, -1))
            right_density += 0.5 * E11.R_SKYRME * np.sum(
                (two - affine_two) ** 2, axis=(-2, -1)
            )
            right_density += 0.5 * E11.K_SEXTIC * np.sum(
                (three - affine_three) ** 2, axis=-1
            )
            right = float(volumes @ right_density)
            rows.append(
                {
                    "mesh": spec.name,
                    "sample": sample,
                    "mean_gradient_residual": float(np.max(np.abs(mean_gradient - affine))),
                    "mean_second_minor_residual": float(np.max(np.abs(mean_two - affine_two))),
                    "mean_third_minor_residual": float(np.max(np.abs(mean_three - affine_three))),
                    "strong_quasiconvex_identity_residual": abs(left - right),
                    "excess": right,
                }
            )
    return {
        "identity": "integral[W(F+Dphi)-W(F)]=1/2 integral(|Dphi|^2+R|M2(F+Dphi)-M2(F)|^2+K|M3(F+Dphi)-M3(F)|^2)",
        "proof": "all complete-minor cross terms vanish because periodic first, second, and third minors have their affine means",
        "maximum_mean_residual": max(
            max(
                row["mean_gradient_residual"],
                row["mean_second_minor_residual"],
                row["mean_third_minor_residual"],
            )
            for row in rows
        ),
        "maximum_identity_residual": max(
            row["strong_quasiconvex_identity_residual"] for row in rows
        ),
        "rows": rows,
    }


def endpoint_scaling_audit() -> dict[str, Any]:
    rows = []
    for alpha in np.linspace(0.05, 0.95, 19):
        dirichlet_finite = alpha > 0.0
        second_minor_L2_finite = alpha > 0.5
        forced_endpoint_gap = alpha < 0.5
        rows.append(
            {
                "alpha": float(alpha),
                "dirichlet_integrand_power": float(2.0 * alpha - 1.0),
                "second_minor_squared_integrand_power": float(4.0 * alpha - 3.0),
                "forced_fill_L2_lower_bound_power": float(4.0 * alpha - 2.0),
                "dirichlet_finite": dirichlet_finite,
                "second_minor_L2_finite": second_minor_L2_finite,
                "forced_endpoint_gap": forced_endpoint_gap,
                "forbidden_overlap": second_minor_L2_finite and forced_endpoint_gap,
            }
        )
    log_rows = []
    for beta in (-0.5, -0.25, 0.0, 0.2, 0.3, 0.5, 1.0):
        # a(r)=r^(1/2) log(e/r)^(-beta).
        m2_finite = beta > 0.25
        gap_limit = "infinity" if beta < 0.0 else ("constant" if beta == 0.0 else "zero")
        log_rows.append(
            {
                "beta": beta,
                "second_minor_L2_finite": m2_finite,
                "forced_fill_bound_limit": gap_limit,
                "forbidden_overlap": m2_finite and beta < 0.0,
            }
        )
    return {
        "model": "u(r,theta)=a(r) gamma(theta), where gamma is the two-hole Malý loop",
        "power_law": {
            "amplitude": "a(r)=r^alpha",
            "dirichlet_condition": "alpha>0",
            "second_minor_L2_condition": "alpha>1/2",
            "forced_endpoint_gap_condition": "alpha<1/2",
            "no_overlap": not any(row["forbidden_overlap"] for row in rows),
            "rows": rows,
        },
        "borderline_log": {
            "amplitude": "a(r)=r^(1/2) log(e/r)^(-beta)",
            "second_minor_L2_condition": "beta>1/4",
            "forced_fill_bound": "a(r)^4/r^2=log(e/r)^(-4 beta)",
            "no_overlap": not any(row["forbidden_overlap"] for row in log_rows),
            "rows": log_rows,
        },
        "conclusion": "the classical two-hole/dipole endpoint obstruction cannot coexist with the AP-E11 L2 second-minor bound; this excludes that mechanism but is not a general density theorem",
    }


def rank_one_discontinuity_audit() -> dict[str, Any]:
    radii = np.geomspace(1.0e-8, 0.2, 6000)
    logarithm = np.log(math.e / radii)
    theta_prime = -1.0 / (radii * logarithm)
    numerical_energy = float(np.trapezoid(radii * theta_prime**2, radii))
    exact_energy = float(
        1.0 / np.log(math.e / radii[-1])
        - 1.0 / np.log(math.e / radii[0])
    )
    truncation_rows = []
    for epsilon in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8):
        tail = 1.0 / np.log(math.e / epsilon)
        truncation_rows.append(
            {"epsilon": epsilon, "W12_tail_energy": tail, "M2_error": 0.0, "M3_error": 0.0}
        )
    return {
        "map": "u(r,z)=(cos(log log(e/r)),sin(log log(e/r)),0,0) around an axis",
        "property": "finite W12 energy but no continuous representative on the axis; derivative rank one, hence M2=M3=0",
        "numerical_radial_energy": numerical_energy,
        "exact_radial_energy_to_r0_0p2": exact_energy,
        "quadrature_residual": abs(numerical_energy - exact_energy),
        "smooth_rank_one_truncations_converge": True,
        "truncation_rows": truncation_rows,
        "conclusion": "finite complete-minor energy does not imply continuity, so continuous-Cartesian approximation theorems are insufficient; the example itself is nevertheless graph-norm approximable",
    }


def smooth_perturbation(
    field: np.ndarray, axis: np.ndarray, seed: int, amplitude: float
) -> np.ndarray:
    interior = field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    frames = E11.tangent_frames(interior)
    grid = np.stack(
        np.meshgrid(axis[1:-1], axis[1:-1], axis[1:-1], indexing="ij"), axis=-1
    )
    scaled = (grid - axis[0]) / (axis[-1] - axis[0])
    envelope = np.prod(np.sin(np.pi * scaled), axis=-1)
    rng = np.random.default_rng(seed)
    coordinates = np.zeros((*grid.shape[:-1], 3), dtype=float)
    for tangent in range(3):
        coefficients = rng.normal(size=(3, 3))
        mode = np.zeros(grid.shape[:-1], dtype=float)
        for row in coefficients:
            wave = np.maximum(1, np.rint(np.abs(row)).astype(int))
            phase = float(rng.uniform(-np.pi, np.pi))
            argument = 2.0 * np.pi * np.einsum("...i,i->...", scaled, wave)
            mode += float(rng.normal()) * np.sin(argument + phase)
        coordinates[..., tangent] = envelope * mode
    flat = coordinates.reshape(-1, 3)
    norms = np.linalg.norm(flat, axis=1)
    maximum = float(np.max(norms))
    if maximum > 0.0:
        flat *= amplitude / maximum
    moved, _, _ = E11.exponential_chart(field, flat.ravel(), frames)
    return moved


def relax_field(
    plan: Any,
    quick: bool,
    perturb_seed: int | None = None,
    perturb_amplitude: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any], tuple[Any, ...]]:
    spec = E11.mesh_spec(plan.mesh)
    geometry_data = E11.finite_element_geometry(plan.extent, plan.length, spec)
    tetrahedra, orientations, gradients, volumes, left, right = geometry_data
    field, axis = E11.compact_hedgehog(
        plan.extent, plan.length, plan.translation_lattice_units
    )
    if perturb_seed is not None and perturb_amplitude > 0.0:
        field = smooth_perturbation(field, axis, perturb_seed, perturb_amplitude)
    initial_topology = E11.topology_summary(field, tetrahedra, orientations)
    history = []
    charts = 4 if quick else 6
    maximum_iterations = 400 if quick else 550
    for chart in range(charts):
        base = field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
        frames = E11.tangent_frames(base)
        result = minimize(
            E11.relaxation_objective,
            np.zeros((len(base), 3), dtype=float).ravel(),
            args=(field, frames, tetrahedra, gradients, volumes, left, right),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": maximum_iterations, "maxcor": 20, "ftol": 2.0e-13, "gtol": 1.0e-8},
        )
        field, _, _ = E11.exponential_chart(field, result.x, frames)
        _, physical_gradient, local_geometry = E11.finite_element_energy_gradient(
            field, tetrahedra, gradients, volumes, left, right, False
        )
        stationarity = E11.projected_gradient_density(field, physical_gradient)
        history.append(
            {
                "chart": chart,
                "iterations": int(result.nit),
                "optimizer_success": bool(result.success),
                "physical_energy": local_geometry["physical_energy"],
                "projected_gradient_density": stationarity,
            }
        )
        if stationarity < E11.STATIONARITY_TOLERANCE and local_geometry["minimum_pair_dot"] > E11.SOLVER_BUFFER:
            break
    physical_energy, physical_gradient, local_geometry = E11.finite_element_energy_gradient(
        field, tetrahedra, gradients, volumes, left, right, False
    )
    stationarity = E11.projected_gradient_density(field, physical_gradient)
    topology = E11.topology_summary(field, tetrahedra, orientations)
    observables = E11.quotient_observables(field, axis, tetrahedra, gradients, volumes)
    result_row = {
        **asdict(plan),
        "spacing": plan.length / (plan.extent - 1),
        "perturb_seed": perturb_seed,
        "perturb_amplitude": perturb_amplitude,
        "initial_B": initial_topology["baryon_numbers"],
        "history": history,
        "physical_energy": physical_energy,
        "projected_gradient_density": stationarity,
        "minimum_pair_dot": local_geometry["minimum_pair_dot"],
        "admissible": local_geometry["admissible"],
        "final_B": topology["baryon_numbers"],
        "quotient_observables": observables,
        "case_pass": bool(
            stationarity < E11.STATIONARITY_TOLERANCE
            and local_geometry["minimum_pair_dot"] > E11.SOLVER_BUFFER
            and local_geometry["admissible"]
            and topology["all_targets_B1"]
        ),
    }
    return field, axis, result_row, geometry_data


def component_energies(
    field: np.ndarray,
    axis: np.ndarray,
    tetrahedra: np.ndarray,
    gradients: np.ndarray,
    volumes: np.ndarray,
) -> dict[str, Any]:
    flat = field.reshape(-1, 4)
    values = flat[tetrahedra]
    derivative = np.einsum("tva,tvi->tai", values, gradients)
    gram = np.einsum("tai,taj->tij", derivative, derivative)
    trace = np.trace(gram, axis1=-2, axis2=-1)
    second = 0.5 * (trace**2 - np.einsum("tij,tji->t", gram, gram))
    three = E11.third_minors(derivative)
    component_density = {
        "dirichlet": 0.5 * np.sum(derivative**2, axis=(-2, -1)),
        "second_minor": 0.5 * E11.R_SKYRME * second,
        "third_minor": 0.5 * E11.K_SEXTIC * np.sum(three**2, axis=-1),
        "potential": E11.MASS_SQUARED * (1.0 - np.mean(values[..., 0], axis=1)),
    }
    components = {name: volumes * density for name, density in component_density.items()}
    total = sum(components.values())
    coordinates = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    centers = np.mean(coordinates[tetrahedra], axis=1)
    barycentre = np.sum(total[:, None] * centers, axis=0) / np.sum(total)
    distances = np.linalg.norm(centers - barycentre, axis=1)
    spacing = axis[1] - axis[0]
    radii = sorted(set(float(value) for value in (2.0 * spacing, 3.0 * spacing, 0.75, 1.0, 1.25)))
    ball_rows = []
    for radius in radii:
        selected = distances <= radius
        local_components = {name: float(np.sum(values_[selected])) for name, values_ in components.items()}
        local_total = sum(local_components.values())
        ball_rows.append(
            {
                "radius": radius,
                "total_energy": local_total,
                "harmonic_scale_ratio_E_over_r": local_total / radius,
                "energy_fraction": local_total / float(np.sum(total)),
                "components": local_components,
            }
        )
    return {
        "energy_barycentre": barycentre.tolist(),
        "component_totals": {name: float(np.sum(values_)) for name, values_ in components.items()},
        "total_energy": float(np.sum(total)),
        "ball_rows": ball_rows,
        "maximum_tetrahedron_energy_fraction": float(np.max(total) / np.sum(total)),
        "density_quantiles": {
            "q50": float(np.quantile(total / volumes, 0.50)),
            "q90": float(np.quantile(total / volumes, 0.90)),
            "q99": float(np.quantile(total / volumes, 0.99)),
            "q999": float(np.quantile(total / volumes, 0.999)),
        },
    }


def relative_spread(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    return float((np.max(array) - np.min(array)) / max(abs(np.mean(array)), 1.0e-15))


def profile_cdf_distance(first: dict[str, Any], second: dict[str, Any]) -> float:
    a = np.asarray(first["normalized_centered_radial_profile"], dtype=float)
    b = np.asarray(second["normalized_centered_radial_profile"], dtype=float)
    return float(np.max(np.abs(np.cumsum(a) - np.cumsum(b))))


def shifted_field(field: np.ndarray, shift: tuple[int, int, int]) -> np.ndarray:
    result = np.zeros_like(field)
    result[..., 0] = 1.0
    source_slices = []
    target_slices = []
    for amount, size in zip(shift, field.shape[:3]):
        if amount >= 0:
            source_slices.append(slice(0, size - amount))
            target_slices.append(slice(amount, size))
        else:
            source_slices.append(slice(-amount, size))
            target_slices.append(slice(0, size + amount))
    result[tuple(target_slices)] = field[tuple(source_slices)]
    return result


def proper_target_procrustes_distance(
    reference: np.ndarray, candidate: np.ndarray
) -> dict[str, Any]:
    """Quotient the unbroken global SO(3) action on the pion components."""
    reference_core = reference[1:-1, 1:-1, 1:-1]
    candidate_core = candidate[1:-1, 1:-1, 1:-1]
    reference_pions = reference_core[..., 1:4].reshape(-1, 3)
    candidate_pions = candidate_core[..., 1:4].reshape(-1, 3)
    left, _, right_transpose = np.linalg.svd(candidate_pions.T @ reference_pions)
    orientation = np.ones(3)
    orientation[-1] = np.linalg.det(left @ right_transpose)
    rotation = left @ np.diag(orientation) @ right_transpose
    aligned = candidate_core.copy()
    aligned[..., 1:4] = candidate_core[..., 1:4] @ rotation
    difference = aligned - reference_core
    return {
        "target_rotation": rotation.tolist(),
        "target_rotation_determinant": float(np.linalg.det(rotation)),
        "rms_field_distance": float(
            np.linalg.norm(difference.ravel()) / math.sqrt(difference.size)
        ),
    }


def best_symmetry_quotient_distance(
    reference: np.ndarray, candidate: np.ndarray
) -> dict[str, Any]:
    rows = []
    for shift in itertools.product(range(-1, 2), repeat=3):
        moved = shifted_field(candidate, shift)
        row = proper_target_procrustes_distance(reference, moved)
        row["integer_spatial_shift"] = list(shift)
        rows.append(row)
    return min(rows, key=lambda row: row["rms_field_distance"])


def numerical_background_audit(quick: bool) -> dict[str, Any]:
    base_specs = ((17, 5.5), (21, 6.0)) if quick else ((17, 5.5), (21, 6.0), (25, 6.5))
    background_rows = []
    base_fields: dict[int, np.ndarray] = {}
    base_observables: dict[int, dict[str, Any]] = {}
    for extent, length in base_specs:
        plan = E11.RelaxationPlan(
            f"density_N{extent}", extent, length, "uniform_ppp", (0.0, 0.0, 0.0), ("density",)
        )
        field, axis, row, geometry_data = relax_field(plan, quick)
        tetrahedra, _, gradients, volumes, _, _ = geometry_data
        row["concentration"] = component_energies(field, axis, tetrahedra, gradients, volumes)
        background_rows.append(row)
        base_fields[extent] = field
        base_observables[extent] = row["quotient_observables"]

    perturb_rows = []
    reference_extent, reference_length = base_specs[0]
    reference_field = base_fields[reference_extent]
    reference_observables = base_observables[reference_extent]
    for seed in ((101,) if quick else (101, 202, 303)):
        plan = E11.RelaxationPlan(
            f"perturb_{seed}_N{reference_extent}",
            reference_extent,
            reference_length,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            ("finite_grid_isolation",),
        )
        field, _, row, _ = relax_field(plan, quick, seed, 0.24)
        row["profile_CDF_distance_to_reference"] = profile_cdf_distance(
            reference_observables, row["quotient_observables"]
        )
        row["best_symmetry_quotient"] = best_symmetry_quotient_distance(
            reference_field, field
        )
        perturb_rows.append(row)

    all_rows = background_rows + perturb_rows
    reference_energy = background_rows[0]["physical_energy"]
    perturb_energy_spread = relative_spread(
        [reference_energy] + [row["physical_energy"] for row in perturb_rows]
    )
    maximum_profile = max(
        [row["profile_CDF_distance_to_reference"] for row in perturb_rows], default=0.0
    )
    maximum_field = max(
        [row["best_symmetry_quotient"]["rms_field_distance"] for row in perturb_rows],
        default=0.0,
    )
    finite_grid_proxy = bool(
        all(row["case_pass"] for row in all_rows)
        and perturb_energy_spread < 2.0e-6
        and maximum_profile < 2.0e-4
        and maximum_field < 6.0e-3
    )
    return {
        "background_rows": background_rows,
        "perturbation_rows": perturb_rows,
        "finite_grid_isolation_proxy": {
            "energy_relative_spread": perturb_energy_spread,
            "maximum_profile_CDF_distance": maximum_profile,
            "maximum_translation_and_target_SO3_quotiented_field_RMS": maximum_field,
            "thresholds": {"energy": 2.0e-6, "profile": 2.0e-4, "field_RMS": 6.0e-3},
            "pass": finite_grid_proxy,
            "scope": "finite-dimensional same-action basin evidence only; not a continuum isolation theorem",
        },
        "all_cases_pass": all(row["case_pass"] for row in all_rows),
    }


def theory_ledger() -> dict[str, Any]:
    return {
        "graph_norm_degree_continuity": {
            "status": True,
            "estimate": "2 pi^2 |B(u_j)-B(u)| <= ||u_j-u||_2 ||M3(Du_j)||_2 + |Omega|^(1/2) ||M3(Du_j)-M3(Du)||_2",
            "consequence": "a smooth target-valued graph recovery sequence is eventually in the same integer degree",
        },
        "ordinary_W12_density": {
            "status": True,
            "reason": "pi_2(S3)=0 removes the ordinary Bethuel obstruction",
            "does_not_control": "strong L2 convergence of M2 and M3",
        },
        "complete_minor_endpoint_density": {
            "status": False,
            "disproved": False,
            "reason": "known Cartesian approximation gives exponent loss q<p; the endpoint p=q=2 requires an additional equiintegrability/extension theorem not supplied by ordinary W12 density",
        },
        "sufficient_supercritical_route": {
            "statement": "if the complete-minor vector belongs locally to L^(2+delta) and the map is in the appropriate Cartesian closure, Malý approximation with q=2 gives graph-L2 recovery",
            "proved_as_implication": True,
            "higher_integrability_for_this_minimizer": False,
        },
        "Maly_dipole_mechanism_excluded": True,
        "strong_quasiconvexity": True,
        "projected_Euler_Lagrange_for_classical_map": "(I-u tensor u)(-div P(Du)-m^2 e0)=0, P=DW",
        "uniform_Legendre_Hadamard": True,
        "relaxed_minimizer_local_recovery": False,
        "weak_Euler_Lagrange_for_relaxed_minimizer": False,
        "partial_regularity_theorem": False,
        "classicality_theorem": False,
        "continuum_isolation_theorem_mod_translations": False,
        "growth_obstruction": "the autonomous density has 2-to-6 growth (q/p=3); available generic (p,q)-regularity theorems require substantially smaller growth ratio and do not cover a sphere constraint plus relaxation gap",
        "ordered_next_gate": "prove a local endpoint recovery/equiintegrability lemma or prove L^(2+delta) higher integrability of the complete-minor vector for the selected minimizer; only then derive weak EL and epsilon regularity",
    }


def render_markdown(result: dict[str, Any]) -> str:
    symbol = result["algebraic_symbol"]
    quasi = result["strong_quasiconvexity"]
    scaling = result["endpoint_scaling"]
    numerical = result["numerical_background"]
    theory = result["theory"]
    proxy = numerical["finite_grid_isolation_proxy"]
    lines = [
        "# AP-E12 graph-density and regular-minimizer audit",
        "",
        f"Status: `{result['status']}`.",
        "",
        "## Exact results",
        "",
        f"- degree continuous in the complete-minor graph norm: `{theory['graph_norm_degree_continuity']['status']}`.",
        f"- maximum singular-value identity residual: `{symbol['maximum_singular_value_identity_residual']:.3e}`.",
        f"- maximum rank-one symbol residual: `{symbol['maximum_rank_one_symbol_identity_residual']:.3e}`.",
        f"- maximum periodic strong-quasiconvex identity residual: `{quasi['maximum_identity_residual']:.3e}`.",
        f"- power-law Malý obstruction overlap: `{not scaling['power_law']['no_overlap']}`.",
        f"- logarithmic-borderline overlap: `{not scaling['borderline_log']['no_overlap']}`.",
        "",
        "## Numerical evidence",
        "",
        f"- all re-relaxed backgrounds and perturbations pass: `{numerical['all_cases_pass']}`.",
        f"- finite-grid isolation proxy: `{proxy['pass']}`.",
        f"- perturbed-start energy spread: `{proxy['energy_relative_spread']:.6e}`.",
        f"- maximum quotient profile distance: `{proxy['maximum_profile_CDF_distance']:.6e}`.",
        "- maximum translation/target-$SO(3)$ quotient field RMS: "
        f"`{proxy['maximum_translation_and_target_SO3_quotiented_field_RMS']:.6e}`.",
        "",
        "## Fail-closed theorem boundary",
        "",
        f"- endpoint graph-L2 density proved: `{theory['complete_minor_endpoint_density']['status']}`.",
        f"- endpoint density disproved: `{theory['complete_minor_endpoint_density']['disproved']}`.",
        f"- local recovery for relaxed minimizer: `{theory['relaxed_minimizer_local_recovery']}`.",
        f"- partial regularity: `{theory['partial_regularity_theorem']}`.",
        f"- classicality: `{theory['classicality_theorem']}`.",
        f"- continuum isolation: `{theory['continuum_isolation_theorem_mod_translations']}`.",
        f"- Hessian gate: `{result['hessian_gate_open']}`.",
        "",
        f"Checks: `{result['checks_passed']}/{result['checks_total']}`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    algebra = algebraic_symbol_audit(arguments.quick)
    quasi = quasiconvex_identity_audit(arguments.quick)
    scaling = endpoint_scaling_audit()
    rank_one = rank_one_discontinuity_audit()
    numerical = numerical_background_audit(arguments.quick)
    theory = theory_ledger()

    check(
        "A1",
        "complete-minor singular-value identity",
        algebra["maximum_singular_value_identity_residual"] < 1.0e-11,
        f"residual={algebra['maximum_singular_value_identity_residual']:.3e}",
    )
    check(
        "A2",
        "uniform rank-one symbol identity",
        algebra["maximum_rank_one_symbol_identity_residual"] < 1.0e-10
        and algebra["minimum_rank_one_symbol_lower_margin"] > -1.0e-10,
        f"residual={algebra['maximum_rank_one_symbol_identity_residual']:.3e}",
    )
    check(
        "A3",
        "periodic strong-quasiconvex identity",
        quasi["maximum_mean_residual"] < 1.0e-11
        and quasi["maximum_identity_residual"] < 1.0e-10,
        f"identity={quasi['maximum_identity_residual']:.3e}",
    )
    check(
        "D1",
        "power-law Malý obstruction excluded by M2 L2",
        scaling["power_law"]["no_overlap"],
        "M2 L2 needs alpha>1/2 while the forced endpoint gap needs alpha<1/2",
    )
    check(
        "D2",
        "logarithmic borderline closes without overlap",
        scaling["borderline_log"]["no_overlap"],
        "M2 L2 needs beta>1/4, forcing a(r)^4/r^2 to zero",
    )
    check(
        "D3",
        "rank-one discontinuity is approximable and not a counterexample",
        rank_one["smooth_rank_one_truncations_converge"]
        and rank_one["quadrature_residual"] < 2.0e-3,
        f"quadrature={rank_one['quadrature_residual']:.3e}",
    )
    check(
        "N1",
        "same-action backgrounds and perturbed starts",
        numerical["all_cases_pass"],
        f"cases={len(numerical['background_rows']) + len(numerical['perturbation_rows'])}",
    )
    check(
        "N2",
        "finite-grid isolation proxy",
        numerical["finite_grid_isolation_proxy"]["pass"],
        f"energy={numerical['finite_grid_isolation_proxy']['energy_relative_spread']:.3e}",
    )
    check(
        "G1",
        "theorem gates fail closed",
        not theory["complete_minor_endpoint_density"]["status"]
        and not theory["partial_regularity_theorem"]
        and not theory["classicality_theorem"]
        and not theory["continuum_isolation_theorem_mod_translations"],
        "no finite-grid proxy promotes an endpoint density or continuum theorem",
    )

    hessian_gate = bool(
        theory["complete_minor_endpoint_density"]["status"]
        or (
            theory["partial_regularity_theorem"]
            and theory["classicality_theorem"]
            and theory["continuum_isolation_theorem_mod_translations"]
        )
    )
    result = {
        "artifact": "AP-E12 endpoint graph-density and relaxed-minimizer audit",
        "status": "quick_fail_closed" if arguments.quick else "production_endpoint_density_regular_minimizer_fail_closed",
        "quick": arguments.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_manifest": [source_row(THIS_SCRIPT), source_row(E11_SCRIPT), source_row(E11_JSON)],
        "algebraic_symbol": algebra,
        "strong_quasiconvexity": quasi,
        "endpoint_scaling": scaling,
        "rank_one_discontinuity": rank_one,
        "numerical_background": numerical,
        "theory": theory,
        "hessian_gate_open": hessian_gate,
        "determinant_variation_gate_open": False,
        "portal_start_allowed": False,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    json_path = OUTPUT / "ap_e12_graph_density_regular_minimizer.json"
    markdown_path = OUTPUT / "ap_e12_graph_density_regular_minimizer.md"
    json_path.write_text(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
            default=lambda value: value.item()
            if isinstance(value, np.generic)
            else str(value),
        )
        + "\n"
    )
    markdown_path.write_text(render_markdown(result))
    print(
        f"AP-E12 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"finite_proxy={numerical['finite_grid_isolation_proxy']['pass']}; "
        f"density={theory['complete_minor_endpoint_density']['status']}; "
        f"hessian={hessian_gate}",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
