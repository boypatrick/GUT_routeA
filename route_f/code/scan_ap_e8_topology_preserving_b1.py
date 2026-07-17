#!/usr/bin/env python3
"""AP-E8 topology-preserving finite-grid B=1 construction.

This card changes the AP-E7 finite-site obstruction by adding an infinite
admissibility barrier on the unique vertex pairs of a fixed Freudenthal
triangulation.  It proves finite-grid sector preservation and existence of an
unanchored minimizer in every nonempty component, and it computes independent
multi-spacing and multi-volume stationary representatives.

The result is deliberately finite-grid.  It is not a Gamma-convergence proof,
not lattice QC2D, and not a physical gauge--meson--ghost--fermion Hessian.
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
from typing import Any

import numpy as np
import scipy
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E7_SCRIPT = ROUTE_F / "code" / "scan_ap_e7_discrete_rerelax_superhessian.py"
E6_SCRIPT = ROUTE_F / "code" / "scan_ap_e6_relaxed_b1_multigrid_hessian.py"

EPSILON = 0.01
GAMMA = 1.0
STATIONARITY_TOLERANCE = 2.0e-6
CHECKS: list[dict[str, Any]] = []


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E7 = load_module("ap_e7_for_ap_e8", E7_SCRIPT)
E6 = E7.E6


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


def freudenthal_unique_edges(extent: int) -> tuple[np.ndarray, np.ndarray]:
    """Unique unordered pairs co-occurring in a Freudenthal tetrahedron."""

    edges: set[tuple[int, int]] = set()
    shape = (extent,) * 3
    for origin_tuple in np.ndindex((extent - 1,) * 3):
        origin = np.asarray(origin_tuple, dtype=np.int64)
        for permutation in itertools.permutations((0, 1, 2)):
            vertices = [origin.copy()]
            running = origin.copy()
            for direction in permutation:
                running = running.copy()
                running[direction] += 1
                vertices.append(running)
            linear = [
                int(np.ravel_multi_index(tuple(vertex), shape))
                for vertex in vertices
            ]
            for first, second in itertools.combinations(linear, 2):
                edges.add((min(first, second), max(first, second)))
    array = np.asarray(sorted(edges), dtype=np.int64)
    return array[:, 0], array[:, 1]


def barrier_value_prime_second(
    dots: np.ndarray, epsilon: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact barrier calculus on the admissible open set dots > epsilon."""

    if not np.all(dots > epsilon):
        raise ValueError("barrier evaluated outside its finite-energy domain")
    denominator = 1.0 - epsilon
    shifted = dots - epsilon
    values = -np.log(shifted / denominator) + (dots - 1.0) / denominator
    prime = -1.0 / shifted + 1.0 / denominator
    second = 1.0 / shifted**2
    return values, prime, second


def barrier_energy_gradient(
    field: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    spacing: float,
    gamma: float,
    epsilon: float,
    solver_extension: bool = False,
) -> tuple[float, np.ndarray, dict[str, Any]]:
    """Barrier energy and ambient gradient.

    The declared action is infinite for a dot product <= epsilon.  Optimizer
    trial points can optionally receive a large differentiable rejection
    extension; no result at such a point is ever reported as physical.
    """

    flat = field.reshape(-1, 4)
    dots = np.einsum("ij,ij->i", flat[left], flat[right])
    weight = gamma * spacing
    invalid = dots <= epsilon
    if solver_extension:
        # Solver-only C1 rejection continuation.  Match the exact function,
        # first derivative, and second derivative at scaled=cutoff, then add a
        # C1 quartic penalty below zero.  Energy and returned gradient are an
        # exact pair, unlike clipping the logarithm while retaining its pole.
        # Reported fields are always re-evaluated with solver_extension=False.
        denominator = 1.0 - epsilon
        scaled = (dots - epsilon) / denominator
        cutoff = 1.0e-4
        negative = np.minimum(scaled, 0.0)
        exact_mask = scaled >= cutoff
        values = np.empty_like(scaled)
        derivative_scaled = np.empty_like(scaled)
        values[exact_mask] = (
            -np.log(scaled[exact_mask]) + scaled[exact_mask] - 1.0
        )
        derivative_scaled[exact_mask] = -1.0 / scaled[exact_mask] + 1.0
        displacement = scaled[~exact_mask] - cutoff
        value_cutoff = -math.log(cutoff) + cutoff - 1.0
        prime_cutoff = -1.0 / cutoff + 1.0
        second_cutoff = 1.0 / cutoff**2
        values[~exact_mask] = (
            value_cutoff
            + prime_cutoff * displacement
            + 0.5 * second_cutoff * displacement**2
        )
        derivative_scaled[~exact_mask] = (
            prime_cutoff + second_cutoff * displacement
        )
        quartic_strength = 1.0e6
        values += quartic_strength * negative**4
        derivative_scaled += 4.0 * quartic_strength * negative**3
        energy = weight * float(np.sum(values))
        coefficients = weight * derivative_scaled / denominator
        gradient_flat = np.column_stack(
            [
                np.bincount(
                    left,
                    weights=coefficients * flat[right, component],
                    minlength=len(flat),
                )
                + np.bincount(
                    right,
                    weights=coefficients * flat[left, component],
                    minlength=len(flat),
                )
                for component in range(4)
            ]
        )
        return energy, gradient_flat.reshape(field.shape), {
            "minimum_pair_dot": float(np.min(dots)),
            "minimum_pair_margin": float(np.min(dots) - epsilon),
            "invalid_edges": int(np.count_nonzero(invalid)),
            "finite_declared_action": not bool(np.any(invalid)),
        }
    if np.any(invalid):
        raise ValueError("infinite topology barrier")
    values, prime, _ = barrier_value_prime_second(dots, epsilon)
    coefficients = weight * prime
    gradient_flat = np.column_stack(
        [
            np.bincount(
                left,
                weights=coefficients * flat[right, component],
                minlength=len(flat),
            )
            + np.bincount(
                right,
                weights=coefficients * flat[left, component],
                minlength=len(flat),
            )
            for component in range(4)
        ]
    )
    return weight * float(np.sum(values)), gradient_flat.reshape(field.shape), {
        "minimum_pair_dot": float(np.min(dots)),
        "minimum_pair_margin": float(np.min(dots) - epsilon),
        "invalid_edges": 0,
        "finite_declared_action": True,
    }


def total_energy_gradient(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    epsilon: float,
    solver_extension: bool = False,
) -> tuple[float, np.ndarray, np.ndarray, dict[str, Any]]:
    base_energy, gradient_n, gradient_s = E7.energy_embedding_gradient(
        field, breathing, spacing, model
    )
    barrier_energy, barrier_gradient, barrier_geometry = barrier_energy_gradient(
        field,
        left,
        right,
        spacing,
        gamma,
        epsilon,
        solver_extension=solver_extension,
    )
    geometry = dict(barrier_geometry)
    geometry["base_energy"] = float(base_energy)
    geometry["barrier_energy"] = float(barrier_energy)
    return (
        float(base_energy + barrier_energy),
        gradient_n + barrier_gradient,
        gradient_s,
        geometry,
    )


def chart_objective(
    coordinates: np.ndarray,
    base_field: np.ndarray,
    base_breathing: np.ndarray,
    base: np.ndarray,
    scalar: np.ndarray,
    frames: np.ndarray,
    spacing: float,
    model: Any,
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    epsilon: float,
) -> tuple[float, np.ndarray]:
    field, breathing, z, radius = E7.chart_fields(
        coordinates, base_field, base_breathing, base, scalar, frames
    )
    energy, gradient_n, gradient_s, _ = total_energy_gradient(
        field,
        breathing,
        spacing,
        model,
        left,
        right,
        gamma,
        epsilon,
        solver_extension=True,
    )
    embedding = gradient_n[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    variables = coordinates.reshape(-1, 4)
    derivative = np.zeros_like(variables)
    small = radius <= 1.0e-8
    if np.any(small):
        derivative[small, :3] = np.einsum(
            "nij,ni->nj", frames[small], embedding[small]
        )
    if np.any(~small):
        selected = ~small
        r = radius[selected]
        unit = z[selected] / r[:, None]
        image_unit = np.einsum("nij,nj->ni", frames[selected], unit)
        radial_image = (
            -np.sin(r)[:, None] * base[selected]
            + np.cos(r)[:, None] * image_unit
        )
        radial = np.sum(radial_image * embedding[selected], axis=1)
        pulled = np.einsum("nij,ni->nj", frames[selected], embedding[selected])
        transverse = pulled - unit * np.sum(unit * pulled, axis=1)[:, None]
        derivative[selected, :3] = (
            unit * radial[:, None]
            + (np.sin(r) / r)[:, None] * transverse
        )
    derivative[:, 3] = gradient_s[1:-1, 1:-1, 1:-1].ravel()
    return energy, derivative.ravel()


def direct_stationarity(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    epsilon: float,
) -> dict[str, Any]:
    energy, gradient_n, gradient_s, geometry = total_energy_gradient(
        field,
        breathing,
        spacing,
        model,
        left,
        right,
        gamma,
        epsilon,
    )
    interior_field = field[1:-1, 1:-1, 1:-1]
    interior_gradient = gradient_n[1:-1, 1:-1, 1:-1]
    radial = np.sum(interior_gradient * interior_field, axis=-1, keepdims=True)
    tangent = interior_gradient - radial * interior_field
    scalar_gradient = gradient_s[1:-1, 1:-1, 1:-1]
    dofs = 4 * tangent.shape[0] * tangent.shape[1] * tangent.shape[2]
    norm = math.sqrt(float(np.sum(tangent**2) + np.sum(scalar_gradient**2)))
    return {
        "total_energy": energy,
        "base_energy": geometry["base_energy"],
        "barrier_energy": geometry["barrier_energy"],
        "direct_tangent_plus_scalar_gradient_l2": norm,
        "direct_projected_gradient_l2_density": norm / math.sqrt(dofs) / spacing**3,
        "minimum_pair_dot": geometry["minimum_pair_dot"],
        "minimum_pair_margin": geometry["minimum_pair_margin"],
        "invalid_edges": geometry["invalid_edges"],
    }


def soliton_radius(field: np.ndarray, axis: np.ndarray) -> dict[str, float]:
    coordinates = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    )
    weight = np.maximum(1.0 - field[..., 0], 0.0)
    total = float(np.sum(weight))
    radius_sq = np.sum(coordinates**2, axis=-1)
    rms = math.sqrt(float(np.sum(weight * radius_sq)) / max(total, 1.0e-300))
    spacing = float(axis[1] - axis[0])
    return {
        "one_minus_n0_weight": total,
        "rms_radius": rms,
        "rms_radius_over_a": rms / spacing,
    }


def topology_summary(field: np.ndarray) -> dict[str, Any]:
    topology = E7.geometric_preimage_degree(field)
    baryon_rows = [row["baryon_number"] for row in topology["target_rows"]]
    return {
        "method": topology["method"],
        "targets_agree": topology["targets_agree"],
        "baryon_numbers": baryon_rows,
        "all_targets_B1": topology["targets_agree"] and baryon_rows == [1, 1, 1],
        "target_rows": topology["target_rows"],
        "minimum_nonzero_absolute_tetrahedron_determinant": topology[
            "minimum_nonzero_absolute_tetrahedron_determinant"
        ],
    }


def relax_case(
    solution: Any,
    model: Any,
    extent: int,
    length: float,
    gamma: float = GAMMA,
    epsilon: float = EPSILON,
    seed: int | None = None,
    maximum_recentres: int = 5,
) -> dict[str, Any]:
    field, breathing, axis, _ = E6.sampled_background(solution, extent, length)
    field, breathing = E7.exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    left, right = freudenthal_unique_edges(extent)
    if seed is not None:
        rng = np.random.default_rng(seed)
        base, scalar, frames = E7.fixed_chart(field, breathing)
        coordinates = np.zeros(4 * len(base))
        perturbation = coordinates.reshape(-1, 4)
        perturbation[:, :3] = 2.0e-3 * rng.normal(size=(len(base), 3))
        perturbation[:, 3] = 2.0e-3 * rng.normal(size=len(base))
        field, breathing, _, _ = E7.chart_fields(
            coordinates, field, breathing, base, scalar, frames
        )
    initial = direct_stationarity(
        field, breathing, spacing, model, left, right, gamma, epsilon
    )
    initial_topology = topology_summary(field)
    history: list[dict[str, Any]] = []
    success = False
    for recenter in range(maximum_recentres):
        base, scalar, frames = E7.fixed_chart(field, breathing)

        def objective(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
            return chart_objective(
                coordinates,
                field,
                breathing,
                base,
                scalar,
                frames,
                spacing,
                model,
                left,
                right,
                gamma,
                epsilon,
            )

        result = minimize(
            objective,
            np.zeros(4 * len(base)),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": 700,
                "maxcor": 30,
                "maxls": 80,
                "ftol": 1.0e-15,
                "gtol": 1.0e-10,
            },
        )
        field, breathing, _, _ = E7.chart_fields(
            result.x, field, breathing, base, scalar, frames
        )
        stationarity = direct_stationarity(
            field, breathing, spacing, model, left, right, gamma, epsilon
        )
        history.append(
            {
                "recenter": recenter,
                "optimizer_success": bool(result.success),
                "optimizer_message": str(result.message),
                "iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "direct_projected_gradient_l2_density": stationarity[
                    "direct_projected_gradient_l2_density"
                ],
                "energy": stationarity["total_energy"],
                "minimum_pair_dot": stationarity["minimum_pair_dot"],
            }
        )
        if (
            stationarity["direct_projected_gradient_l2_density"]
            < STATIONARITY_TOLERANCE
        ):
            success = True
            break
    final = direct_stationarity(
        field, breathing, spacing, model, left, right, gamma, epsilon
    )
    final_topology = topology_summary(field)
    lower_bound = math.sqrt((1.0 + 3.0 * epsilon) / 4.0)
    return {
        "N": extent,
        "L": length,
        "a": spacing,
        "gamma": gamma,
        "epsilon": epsilon,
        "unique_freudenthal_pair_edges": int(len(left)),
        "core_or_centre_anchor_used": False,
        "deterministic_perturbation_seed": seed,
        "initial": initial,
        "initial_topology": initial_topology,
        "recentering_history": history,
        "stationary_by_direct_tangent_test": success,
        "final": final,
        "final_topology": final_topology,
        "analytic_normalized_affine_norm_lower_bound": lower_bound,
        "radius": soliton_radius(field, axis),
    }


def gradient_audit(solution: Any, model: Any) -> dict[str, Any]:
    extent, length = 15, 6.0
    field, breathing, axis, _ = E6.sampled_background(solution, extent, length)
    field, breathing = E7.exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    left, right = freudenthal_unique_edges(extent)
    base, scalar, frames = E7.fixed_chart(field, breathing)
    rng = np.random.default_rng(20260717)
    coordinates = 1.0e-3 * rng.normal(size=4 * len(base))
    direction = rng.normal(size=len(coordinates))
    direction /= np.linalg.norm(direction)
    energy, gradient = chart_objective(
        coordinates,
        field,
        breathing,
        base,
        scalar,
        frames,
        spacing,
        model,
        left,
        right,
        GAMMA,
        EPSILON,
    )
    predicted = float(gradient @ direction)
    rows = []
    for step in (2.0e-5, 1.0e-5, 5.0e-6):
        plus = chart_objective(
            coordinates + step * direction,
            field,
            breathing,
            base,
            scalar,
            frames,
            spacing,
            model,
            left,
            right,
            GAMMA,
            EPSILON,
        )[0]
        minus = chart_objective(
            coordinates - step * direction,
            field,
            breathing,
            base,
            scalar,
            frames,
            spacing,
            model,
            left,
            right,
            GAMMA,
            EPSILON,
        )[0]
        finite = (plus - minus) / (2.0 * step)
        rows.append(
            {
                "step": step,
                "finite_difference": finite,
                "analytic": predicted,
                "relative_residual": abs(finite - predicted)
                / max(abs(predicted), 1.0),
            }
        )
    return {
        "energy": energy,
        "analytic_directional_derivative": predicted,
        "rows": rows,
        "best_relative_residual": min(row["relative_residual"] for row in rows),
    }


def solver_continuation_gradient_audit() -> dict[str, Any]:
    """Check the optimizer-only continuation at an inadmissible one-edge point."""

    q = EPSILON - 0.02
    angle = math.acos(q)
    field = np.zeros((2, 1, 1, 4))
    field[0, 0, 0] = np.array([1.0, 0.0, 0.0, 0.0])
    field[1, 0, 0] = np.array([math.cos(angle), math.sin(angle), 0.0, 0.0])
    tangent = np.array([-math.sin(angle), math.cos(angle), 0.0, 0.0])
    left = np.array([0], dtype=np.int64)
    right = np.array([1], dtype=np.int64)
    energy, gradient, geometry = barrier_energy_gradient(
        field,
        left,
        right,
        spacing=1.0,
        gamma=1.0,
        epsilon=EPSILON,
        solver_extension=True,
    )
    predicted = float(gradient[1, 0, 0] @ tangent)

    def value(displacement: float) -> float:
        moved = field.copy()
        moved[1, 0, 0] = np.array(
            [
                math.cos(angle + displacement),
                math.sin(angle + displacement),
                0.0,
                0.0,
            ]
        )
        return barrier_energy_gradient(
            moved,
            left,
            right,
            spacing=1.0,
            gamma=1.0,
            epsilon=EPSILON,
            solver_extension=True,
        )[0]

    rows = []
    for step in (2.0e-6, 1.0e-6, 5.0e-7):
        finite = (value(step) - value(-step)) / (2.0 * step)
        rows.append(
            {
                "step": step,
                "finite_difference": finite,
                "analytic": predicted,
                "relative_residual": abs(finite - predicted)
                / max(abs(predicted), 1.0),
            }
        )
    return {
        "test_pair_dot": q,
        "invalid_edges": geometry["invalid_edges"],
        "energy": energy,
        "analytic_directional_derivative": predicted,
        "rows": rows,
        "best_relative_residual": min(row["relative_residual"] for row in rows),
    }


def barrier_hessian_audit(solution: Any) -> dict[str, Any]:
    extent, length = 15, 6.0
    field, breathing, axis, _ = E6.sampled_background(solution, extent, length)
    field, breathing = E7.exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    left, right = freudenthal_unique_edges(extent)
    flat = field.reshape(-1, 4)
    dots = np.einsum("ij,ij->i", flat[left], flat[right])
    _, prime, second = barrier_value_prime_second(dots, EPSILON)
    rng = np.random.default_rng(20260718)
    tangent = rng.normal(size=field.shape)
    tangent -= np.sum(tangent * field, axis=-1, keepdims=True) * field
    boundary = np.zeros(field.shape[:3], dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    tangent[boundary] = 0.0
    tangent /= np.linalg.norm(tangent)
    tangent_flat = tangent.reshape(-1, 4)
    delta_dot = (
        np.einsum("ij,ij->i", tangent_flat[left], flat[right])
        + np.einsum("ij,ij->i", flat[left], tangent_flat[right])
    )
    edge_second = (
        second * delta_dot**2
        + prime
        * (
            2.0
            * np.einsum("ij,ij->i", tangent_flat[left], tangent_flat[right])
            - dots
            * (
                np.sum(tangent_flat[left] ** 2, axis=1)
                + np.sum(tangent_flat[right] ** 2, axis=1)
            )
        )
    )
    predicted = GAMMA * spacing * float(np.sum(edge_second))
    centre = barrier_energy_gradient(
        field, left, right, spacing, GAMMA, EPSILON
    )[0]

    def move(step: float) -> np.ndarray:
        norms = np.linalg.norm(tangent, axis=-1)
        moved = field.copy()
        nonzero = norms > 0.0
        moved[nonzero] = (
            np.cos(step * norms[nonzero])[:, None] * field[nonzero]
            + (
                np.sin(step * norms[nonzero]) / norms[nonzero]
            )[:, None]
            * tangent[nonzero]
        )
        return moved

    rows = []
    for step in (2.0e-3, 1.0e-3, 5.0e-4):
        plus = barrier_energy_gradient(
            move(step), left, right, spacing, GAMMA, EPSILON
        )[0]
        minus = barrier_energy_gradient(
            move(-step), left, right, spacing, GAMMA, EPSILON
        )[0]
        finite = (plus - 2.0 * centre + minus) / step**2
        rows.append(
            {
                "step": step,
                "finite_difference": finite,
                "analytic": predicted,
                "relative_residual": abs(finite - predicted)
                / max(abs(predicted), 1.0),
            }
        )
    return {
        "analytic_intrinsic_quadratic_form": predicted,
        "rows": rows,
        "best_relative_residual": min(row["relative_residual"] for row in rows),
    }


def compact_supported_degree_one_field(
    extent: int, length: float = 6.0, support_radius: float = 2.0
) -> tuple[np.ndarray, np.ndarray]:
    """C2 compact-support hedgehog compatible with the exact vacuum boundary."""

    axis = np.linspace(-0.5 * length, 0.5 * length, extent)
    coordinates = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    radius = np.linalg.norm(coordinates, axis=-1)
    z = radius / support_radius
    profile = np.zeros_like(radius)
    inside = z < 1.0
    zi = z[inside]
    # h(0)=1, h'(0)=-1, h''(0)=0 and h(1)=h'(1)=h''(1)=0.
    profile[inside] = 1.0 - zi - 4.0 * zi**3 + 7.0 * zi**4 - 3.0 * zi**5
    angle = math.pi * profile
    field = np.zeros(coordinates.shape[:-1] + (4,))
    field[..., 0] = np.cos(angle)
    nonzero = radius > 0.0
    field[nonzero, 1:] = (
        np.sin(angle[nonzero]) / radius[nonzero]
    )[:, None] * coordinates[nonzero]
    return field, axis


def smooth_field_scaling() -> dict[str, Any]:
    rows = []
    for extent in (17, 21, 25, 29, 33):
        field, axis = compact_supported_degree_one_field(extent)
        spacing = float(axis[1] - axis[0])
        left, right = freudenthal_unique_edges(extent)
        energy, _, geometry = barrier_energy_gradient(
            field, left, right, spacing, GAMMA, EPSILON
        )
        rows.append(
            {
                "N": extent,
                "a": spacing,
                "barrier_energy": energy,
                "minimum_pair_dot": geometry["minimum_pair_dot"],
            }
        )
    log_a = np.log([row["a"] for row in rows])
    log_energy = np.log([row["barrier_energy"] for row in rows])
    slope, intercept = np.polyfit(log_a, log_energy, 1)
    predicted = slope * log_a + intercept
    residual = log_energy - predicted
    total = log_energy - np.mean(log_energy)
    r_squared = 1.0 - float(residual @ residual) / float(total @ total)
    local_slopes = [
        float(
            (log_energy[index + 1] - log_energy[index])
            / (log_a[index + 1] - log_a[index])
        )
        for index in range(len(rows) - 1)
    ]
    return {
        "rows": rows,
        "fitted_log_log_slope": float(slope),
        "fitted_intercept": float(intercept),
        "fit_r_squared": r_squared,
        "local_log_log_slopes": local_slopes,
        "profile": "C2 degree-one hedgehog with support radius 2 inside the L=6 box and exact vacuum boundary",
        "scope": "fixed compact-supported boundary-compatible test field only; not a uniform bound on minimizers",
    }


def clean_case(row: dict[str, Any]) -> dict[str, Any]:
    return row


def markdown(result: dict[str, Any]) -> str:
    lines = [
        "# AP-E8 topology-preserving finite-grid B=1 card",
        "",
        f"Status: `{result['status']}`",
        "",
        "## Main unanchored stationary cases",
        "",
        "| N | L | a | E | gradient density | min dot | B rows | R_rms/a |",
        "|---:|---:|---:|---:|---:|---:|:---:|---:|",
    ]
    for row in result["main_stationary_cases"]:
        lines.append(
            "| {N} | {L:.1f} | {a:.6f} | {energy:.9f} | {gradient:.3e} | "
            "{dot:.6f} | {baryon} | {radius:.4f} |".format(
                N=row["N"],
                L=row["L"],
                a=row["a"],
                energy=row["final"]["total_energy"],
                gradient=row["final"]["direct_projected_gradient_l2_density"],
                dot=row["final"]["minimum_pair_dot"],
                baryon=row["final_topology"]["baryon_numbers"],
                radius=row["radius"]["rms_radius_over_a"],
            )
        )
    lines.extend(
        [
            "",
            "## Theorem-level result",
            "",
            "For each fixed grid and each nonempty admissible B=1 component, "
            "the infinite barrier makes bounded sublevels compactly contained "
            "in that component. With the coercive breathing-field mass term, "
            "the action attains an interior unanchored minimum.",
            "",
            "This is a grid-indexed finite-a stationary family, not a proof of "
            "Gamma convergence, regulator independence, or a continuum soliton.",
            "",
            f"Checks: **{result['checks_passed']}/{result['checks_total']}**.",
            f"Physics promotion: `{result['physics_promotion_allowed']}`.",
            f"Portal start: `{result['portal_start_allowed']}`.",
            "",
            "## Checks",
            "",
        ]
    )
    for row in result["checks"]:
        lines.append(
            f"- `[{'PASS' if row['pass'] else 'FAIL'}]` {row['group']}: "
            f"{row['name']} — {row['detail']}"
        )
    return "\n".join(lines) + "\n"


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e8_topology_preserving_b1.json"
    md_path = OUTPUT / "ap_e8_topology_preserving_b1.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    md_path.write_text(markdown(result), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run a smaller smoke grid while retaining fail-closed status",
    )
    args = parser.parse_args()
    CHECKS.clear()
    model = E6.Model()
    profile = E6.solve_relaxed_profile(
        box=16.0, model=model, sample_points=4097
    )
    checksum = profile["profile_sha256_r_F_s_little_endian_float64"]
    solution = profile["solution"]
    source_manifest = [source_row(path) for path in (THIS_SCRIPT, E7_SCRIPT, E6_SCRIPT)]

    check(
        "P0_provenance",
        "the canonical AP-E6 continuum seed has the frozen checksum",
        checksum == "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59",
        f"checksum={checksum}",
    )
    values, prime, second = barrier_value_prime_second(
        np.array([1.0, 0.8, 0.2]), EPSILON
    )
    check(
        "P1_action",
        "the subtracted logarithmic barrier is nonnegative, vacuum-flat, convex, and divergent at the floor",
        abs(values[0]) < 1.0e-15
        and abs(prime[0]) < 1.0e-15
        and np.all(values >= -1.0e-15)
        and np.all(second > 0.0)
        and barrier_value_prime_second(np.array([EPSILON + 1.0e-12]), EPSILON)[0][0] > 20.0,
        f"b(1)={values[0]:.3e}; b'(1)={prime[0]:.3e}; min b''={np.min(second):.3e}",
    )
    lower_bound = math.sqrt((1.0 + 3.0 * EPSILON) / 4.0)
    check(
        "P1_action",
        "pairwise admissibility gives a uniform no-zero bound for normalized affine interpolation",
        lower_bound > 0.5,
        f"sqrt((1+3 epsilon)/4)={lower_bound:.12f}",
    )
    derivative_audit = gradient_audit(solution, model)
    check(
        "P1_calculus",
        "the full AP-E8 chart gradient agrees with a centred directional derivative",
        derivative_audit["best_relative_residual"] < 2.0e-8,
        f"best relative residual={derivative_audit['best_relative_residual']:.3e}",
    )
    continuation_audit = solver_continuation_gradient_audit()
    check(
        "P1_calculus",
        "the optimizer-only invalid-domain continuation returns a derivative consistent with its energy",
        continuation_audit["invalid_edges"] == 1
        and continuation_audit["best_relative_residual"] < 2.0e-7,
        f"best relative residual={continuation_audit['best_relative_residual']:.3e}",
    )
    hessian_audit = barrier_hessian_audit(solution)
    check(
        "P1_calculus",
        "the intrinsic barrier second variation agrees with a centred geodesic second difference",
        hessian_audit["best_relative_residual"] < 2.0e-7,
        f"best relative residual={hessian_audit['best_relative_residual']:.3e}",
    )

    if args.quick:
        main_specs = [(15, 6.0), (17, 6.0), (16, 6.0)]
    else:
        main_specs = [
            (15, 6.0),
            (17, 6.0),
            (19, 6.0),
            (21, 6.0),
            (16, 6.0),
            (21, 8.0),
            (26, 10.0),
        ]
    main_cases = [
        relax_case(solution, model, extent, length)
        for extent, length in main_specs
    ]
    main_stationary = all(
        row["stationary_by_direct_tangent_test"]
        and row["final"]["minimum_pair_margin"] > 0.0
        and row["final_topology"]["all_targets_B1"]
        and not row["core_or_centre_anchor_used"]
        for row in main_cases
    )
    check(
        "P2_stationarity",
        "every main final field is unanchored, admissible, three-target B=1, and stationary by the direct tangent test",
        main_stationary,
        "max gradient={:.3e}; min pair margin={:.3e}; B rows={}".format(
            max(
                row["final"]["direct_projected_gradient_l2_density"]
                for row in main_cases
            ),
            min(row["final"]["minimum_pair_margin"] for row in main_cases),
            [row["final_topology"]["baryon_numbers"] for row in main_cases],
        ),
    )
    fixed_box = [
        row
        for row in main_cases
        if row["L"] == 6.0 and row["N"] in (15, 17, 19, 21)
    ]
    fixed_spacing = [
        row
        for row in main_cases
        if abs(row["a"] - 0.4) < 1.0e-12 and row["N"] in (16, 21, 26)
    ]
    check(
        "P2_sequence",
        "the fixed-box sequence contains four independently re-relaxed spacings",
        args.quick or (len(fixed_box) == 4 and all(row["stationary_by_direct_tangent_test"] for row in fixed_box)),
        f"(N,a)={[(row['N'], row['a']) for row in fixed_box]}",
    )
    volume_relative_spread = (
        (max(row["final"]["total_energy"] for row in fixed_spacing)
         - min(row["final"]["total_energy"] for row in fixed_spacing))
        / max(row["final"]["total_energy"] for row in fixed_spacing)
        if fixed_spacing
        else math.inf
    )
    check(
        "P2_sequence",
        "the fixed-spacing L=6,8,10 sequence is stationary with controlled finite-volume energy spread",
        args.quick
        or (
            len(fixed_spacing) == 3
            and all(row["stationary_by_direct_tangent_test"] for row in fixed_spacing)
            and volume_relative_spread < 0.01
        ),
        f"relative energy spread={volume_relative_spread:.6e}",
    )

    if args.quick:
        parameter_controls: list[dict[str, Any]] = []
    else:
        parameter_controls = [
            relax_case(solution, model, 15, 6.0, gamma=0.5, epsilon=0.01),
            relax_case(solution, model, 15, 6.0, gamma=2.0, epsilon=0.01),
            relax_case(solution, model, 15, 6.0, gamma=1.0, epsilon=0.005),
            relax_case(solution, model, 15, 6.0, gamma=1.0, epsilon=0.02),
        ]
    parameter_pass = args.quick or all(
        row["stationary_by_direct_tangent_test"]
        and row["final_topology"]["all_targets_B1"]
        and row["final"]["minimum_pair_margin"] > 0.0
        for row in parameter_controls
    )
    check(
        "P3_regulator_controls",
        "nearby gamma and epsilon choices retain admissible unanchored B=1 stationary representatives",
        parameter_pass,
        "controls="
        + str(
            [
                (
                    row["gamma"],
                    row["epsilon"],
                    row["final"]["direct_projected_gradient_l2_density"],
                    row["final"]["minimum_pair_margin"],
                )
                for row in parameter_controls
            ]
        ),
    )
    if args.quick:
        perturbed_start_control: dict[str, Any] | None = None
    else:
        perturbed_start_control = relax_case(
            solution, model, 15, 6.0, seed=20260719
        )
    baseline_energy = main_cases[0]["final"]["total_energy"]
    basin_energy_relative_difference = (
        abs(perturbed_start_control["final"]["total_energy"] - baseline_energy)
        / baseline_energy
        if perturbed_start_control is not None
        else None
    )
    check(
        "P3_basin_control",
        "a deterministic tangent-plus-scalar perturbation returns to a numerically indistinguishable same-energy admissible stationary B=1 representative",
        args.quick
        or (
            perturbed_start_control is not None
            and perturbed_start_control["stationary_by_direct_tangent_test"]
            and perturbed_start_control["final_topology"]["all_targets_B1"]
            and perturbed_start_control["final"]["minimum_pair_margin"] > 0.0
            and basin_energy_relative_difference is not None
            and basin_energy_relative_difference < 1.0e-10
        ),
        "seed=20260719; relative energy difference={}".format(
            basin_energy_relative_difference
        ),
    )
    scaling = smooth_field_scaling()
    check(
        "P3_smooth_field",
        "the barrier has a stable approximately quadratic power on a compact-supported boundary-compatible profile",
        1.7 < scaling["fitted_log_log_slope"] < 2.4
        and scaling["fit_r_squared"] > 0.995
        and all(1.7 < value < 2.5 for value in scaling["local_log_log_slopes"][-2:]),
        "slope={:.6f}; R2={:.6f}; local slopes={}".format(
            scaling["fitted_log_log_slope"],
            scaling["fit_r_squared"],
            scaling["local_log_log_slopes"],
        ),
    )

    vacuum = np.zeros((4, 4, 4, 4))
    vacuum[..., 0] = 1.0
    vacuum_s = np.ones((4, 4, 4))
    vacuum_left, vacuum_right = freudenthal_unique_edges(4)
    vacuum_barrier = barrier_energy_gradient(
        vacuum, vacuum_left, vacuum_right, 1.0, GAMMA, EPSILON
    )[0]
    vacuum_topology = topology_summary(vacuum)
    check(
        "P4_controls",
        "the exact vacuum has zero barrier and B=0",
        abs(vacuum_barrier) < 1.0e-14
        and vacuum_topology["targets_agree"]
        and vacuum_topology["baryon_numbers"] == [0, 0, 0],
        f"barrier={vacuum_barrier:.3e}; B={vacuum_topology['baryon_numbers']}",
    )

    theorem_flags = {
        "topology_preserving_finite_grid_action_defined": True,
        "finite_energy_sector_preservation_proven": True,
        "finite_grid_B1_minimizer_exists_in_each_nonempty_component": True,
        "grid_indexed_unanchored_B1_stationary_family_exists": True,
        "numerical_finite_multigrid_unanchored_B1_stationary_family_found": main_stationary and not args.quick,
        "multi_spacing_and_multi_volume_representatives_computed": main_stationary and not args.quick,
    }
    fail_closed = {
        "complete_line_search_segments_certified_admissible": False,
        "global_minimizer_numerically_certified": False,
        "minimizer_sequence_equicoercive": False,
        "gamma_convergence_proven": False,
        "continuum_stationary_limit_constructed": False,
        "barrier_regulator_independence_established": False,
        "triangulation_independence_established": False,
        "physical_qc2d_action_derived": False,
        "physical_gauge_meson_ghost_fermion_superhessian_complete": False,
        "four_dimensional_dynamics_performed": False,
        "quantum_continuum_limit_proven": False,
    }
    result = {
        "artifact": "AP-E8 topology-preserving finite-grid B=1 construction",
        "status": "mechanical_pass_finite_grid_blocker_changed_continuum_fail_closed",
        "generated_utc": "deterministic-no-wall-clock",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "model": asdict(model),
        "canonical_profile_checksum": checksum,
        "source_manifest": source_manifest,
        "action_definition": {
            "base_action": "AP-E7 declared finite-difference n+s action",
            "triangulation": "fixed six-tetrahedron Freudenthal split",
            "edge_measure": "unique unordered vertex pairs co-occurring in any Freudenthal tetrahedron",
            "barrier": "b_epsilon(t)=-log((t-epsilon)/(1-epsilon))+(t-1)/(1-epsilon), t>epsilon; +infinity otherwise",
            "coefficient": "gamma*a",
            "epsilon": EPSILON,
            "gamma": GAMMA,
            "normalized_affine_lower_bound": lower_bound,
        },
        "proof_contract": {
            "sector_preservation": "all six pair dots in every tetrahedron exceed epsilon, so the normalized affine map has norm at least sqrt((1+3 epsilon)/4); its degree is locally constant and the barrier forbids a finite-energy crossing",
            "finite_grid_existence": "compact n-space plus coercive m_s>0 breathing potential and barrier-bounded sublevels yield an interior minimizer in every nonempty admissible component",
            "grid_sequence_nonemptiness": "sampling a smooth compactly supported degree-one continuum field gives pair dots>epsilon on every sufficiently fine grid by uniform continuity",
            "continuum_boundary": "the theorem supplies finite-grid minimizers indexed by a; it gives neither equicoercivity nor convergence of those minimizers",
        },
        "gradient_audit": derivative_audit,
        "solver_continuation_gradient_audit": continuation_audit,
        "barrier_intrinsic_hessian_audit": hessian_audit,
        "main_case_specifications": main_specs,
        "main_stationary_cases": [clean_case(row) for row in main_cases],
        "parameter_controls": [clean_case(row) for row in parameter_controls],
        "perturbed_start_basin_control": clean_case(perturbed_start_control)
        if perturbed_start_control is not None
        else None,
        "smooth_compact_supported_profile_barrier_scaling": scaling,
        "vacuum_control": {
            "barrier_energy": vacuum_barrier,
            "topology": vacuum_topology,
        },
        "stationarity_tolerance": STATIONARITY_TOLERANCE,
        "theorem_flags": theorem_flags,
        "fail_closed": fail_closed,
        "physics_promotion_allowed": False,
        "portal_start_allowed": False,
        "lane_closed": False,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    if result["checks_passed"] != result["checks_total"] or args.quick:
        result["status"] = (
            "quick_smoke_finite_grid_only_physics_fail_closed"
            if args.quick and result["checks_passed"] == result["checks_total"]
            else "mechanical_failure_physics_fail_closed"
        )
        if args.quick:
            result["theorem_flags"]["numerical_finite_multigrid_unanchored_B1_stationary_family_found"] = False
            result["theorem_flags"]["multi_spacing_and_multi_volume_representatives_computed"] = False
    write_outputs(result)
    print(
        f"AP-E8 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"status={result['status']}; lane_closed={result['lane_closed']}"
    )
    for row in main_cases:
        print(
            "N={N} L={L:g} a={a:.6g} E={E:.9f} grad={g:.3e} "
            "min_dot={d:.6f} B={B} R/a={r:.4f}".format(
                N=row["N"],
                L=row["L"],
                a=row["a"],
                E=row["final"]["total_energy"],
                g=row["final"]["direct_projected_gradient_l2_density"],
                d=row["final"]["minimum_pair_dot"],
                B=row["final_topology"]["baryon_numbers"],
                r=row["radius"]["rms_radius_over_a"],
            )
        )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
