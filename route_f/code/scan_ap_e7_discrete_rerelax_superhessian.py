#!/usr/bin/env python3
"""AP-E7 discrete re-relaxation and same-grid operator audit.

This calculation imports the AP-E6 continuum B=1 solution only as an initial
condition.  Every reported cubic background is then minimized against the
same declared finite-difference n+s energy with its sampled Dirichlet boundary
held fixed.  A geometric preimage-counting degree, an exact tangent n+s
Hessian, charge-two diquark, gauge, ghost, and a declared Wilson/Yukawa normal
operator are built on each identical interior grid.

The Wilson/Yukawa normal operator is not the second derivative of a fermion
determinant with respect to the bosonic fields.  Consequently the physical
gauge-meson-ghost-fermion aggregate gate remains false.  Nothing here is HMC,
four-dimensional dynamics, or a quantum continuum limit.
"""

from __future__ import annotations

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
from scipy import sparse
from scipy.optimize import minimize
from scipy.sparse.linalg import eigsh


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E6_SCRIPT = ROUTE_F / "code" / "scan_ap_e6_relaxed_b1_multigrid_hessian.py"
TEX = ROUTE_F / "tex" / "ap_e7_discrete_rerelax_superhessian.tex"
BIB = ROUTE_F / "tex" / "ap_e7_discrete_rerelax_superhessian.bib"

CHECKS: list[dict[str, Any]] = []
TARGETS = (
    np.array([0.0, 1.0, 2.0, 3.0]) / math.sqrt(14.0),
    np.array([0.25, -0.45, 0.70, 0.4993746089]),
    np.array([-0.33, 0.61, 0.22, -0.684]),
)
TARGETS = tuple(target / np.linalg.norm(target) for target in TARGETS)


def load_e6() -> Any:
    spec = importlib.util.spec_from_file_location("ap_e6_relaxed", E6_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import AP-E6 source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E6 = load_e6()


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


def energy_embedding_gradient(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Vectorized copy of the exact AP-E6 declared energy and first derivative."""

    extent = field.shape[0]
    gradient_n = np.zeros_like(field)
    gradient_s = np.zeros_like(breathing)
    energy = 0.0

    for direction in range(3):
        difference_n = np.diff(field, axis=direction)
        difference_s = np.diff(breathing, axis=direction)
        energy += 0.5 * spacing * float(
            np.sum(difference_n**2) + np.sum(difference_s**2)
        )
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        lower[direction] = slice(0, extent - 1)
        upper[direction] = slice(1, extent)
        gradient_n[tuple(lower)] -= spacing * difference_n
        gradient_n[tuple(upper)] += spacing * difference_n
        gradient_s[tuple(lower)] -= spacing * difference_s
        gradient_s[tuple(upper)] += spacing * difference_s

    base = field[:-1, :-1, :-1]
    differences = [
        (field[1:, :-1, :-1] - base) / spacing,
        (field[:-1, 1:, :-1] - base) / spacing,
        (field[:-1, :-1, 1:] - base) / spacing,
    ]
    centre_slice = (slice(0, extent - 1),) * 3
    neighbour_slices = [
        (slice(1, extent), slice(0, extent - 1), slice(0, extent - 1)),
        (slice(0, extent - 1), slice(1, extent), slice(0, extent - 1)),
        (slice(0, extent - 1), slice(0, extent - 1), slice(1, extent)),
    ]
    for first, second in ((0, 1), (0, 2), (1, 2)):
        u = differences[first]
        v = differences[second]
        u_sq = np.sum(u * u, axis=-1)
        v_sq = np.sum(v * v, axis=-1)
        uv = np.sum(u * v, axis=-1)
        energy += 0.5 * model.skyrme_coefficient * spacing**3 * float(
            np.sum(u_sq * v_sq - uv**2)
        )
        grad_u = model.skyrme_coefficient * (
            v_sq[..., None] * u - uv[..., None] * v
        )
        grad_v = model.skyrme_coefficient * (
            u_sq[..., None] * v - uv[..., None] * u
        )
        force_u = spacing**2 * grad_u
        force_v = spacing**2 * grad_v
        gradient_n[centre_slice] -= force_u + force_v
        gradient_n[neighbour_slices[first]] += force_u
        gradient_n[neighbour_slices[second]] += force_v

    n_zero = field[1:-1, 1:-1, 1:-1, 0]
    scalar = breathing[1:-1, 1:-1, 1:-1]
    energy += spacing**3 * float(
        np.sum(
            0.5 * model.breathing_mass**2 * (scalar - 1.0) ** 2
            + model.beta**2 * scalar**2 * (1.0 - n_zero)
        )
    )
    gradient_n[1:-1, 1:-1, 1:-1, 0] -= (
        spacing**3 * model.beta**2 * scalar**2
    )
    gradient_s[1:-1, 1:-1, 1:-1] += spacing**3 * (
        model.breathing_mass**2 * (scalar - 1.0)
        + 2.0 * model.beta**2 * scalar * (1.0 - n_zero)
    )
    return energy, gradient_n, gradient_s


def fixed_chart(
    base_field: np.ndarray, base_breathing: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    base = base_field[1:-1, 1:-1, 1:-1].reshape(-1, 4).copy()
    scalar = base_breathing[1:-1, 1:-1, 1:-1].ravel().copy()
    frames = np.stack([E6.tangent_frame(vector) for vector in base])
    return base, scalar, frames


def chart_fields(
    coordinates: np.ndarray,
    base_field: np.ndarray,
    base_breathing: np.ndarray,
    base: np.ndarray,
    scalar: np.ndarray,
    frames: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    variables = coordinates.reshape(-1, 4)
    z = variables[:, :3]
    radius = np.linalg.norm(z, axis=1)
    tangent = np.einsum("nij,nj->ni", frames, z)
    sinc = np.ones_like(radius)
    nonzero = radius > 1.0e-10
    sinc[nonzero] = np.sin(radius[nonzero]) / radius[nonzero]
    sinc[~nonzero] = 1.0 - radius[~nonzero] ** 2 / 6.0
    moved = np.cos(radius)[:, None] * base + sinc[:, None] * tangent
    field = base_field.copy()
    breathing = base_breathing.copy()
    field[1:-1, 1:-1, 1:-1] = moved.reshape(field.shape[0] - 2, field.shape[0] - 2, field.shape[0] - 2, 4)
    breathing[1:-1, 1:-1, 1:-1] = (
        scalar + variables[:, 3]
    ).reshape(breathing.shape[0] - 2, breathing.shape[0] - 2, breathing.shape[0] - 2)
    return field, breathing, z, radius


def chart_objective(
    coordinates: np.ndarray,
    base_field: np.ndarray,
    base_breathing: np.ndarray,
    base: np.ndarray,
    scalar: np.ndarray,
    frames: np.ndarray,
    spacing: float,
    model: Any,
) -> tuple[float, np.ndarray]:
    field, breathing, z, radius = chart_fields(
        coordinates, base_field, base_breathing, base, scalar, frames
    )
    energy, gradient_n, gradient_s = energy_embedding_gradient(
        field, breathing, spacing, model
    )
    g = gradient_n[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    variables = coordinates.reshape(-1, 4)
    derivative = np.zeros_like(variables)
    small = radius <= 1.0e-8
    if np.any(small):
        derivative[small, :3] = np.einsum(
            "nij,ni->nj", frames[small], g[small]
        )
    if np.any(~small):
        selected = ~small
        r = radius[selected]
        unit = z[selected] / r[:, None]
        eunit = np.einsum("nij,nj->ni", frames[selected], unit)
        radial_image = (
            -np.sin(r)[:, None] * base[selected]
            + np.cos(r)[:, None] * eunit
        )
        radial = np.sum(radial_image * g[selected], axis=1)
        projected_g = np.einsum("nij,ni->nj", frames[selected], g[selected])
        transverse = projected_g - unit * np.sum(unit * projected_g, axis=1)[:, None]
        derivative[selected, :3] = (
            unit * radial[:, None]
            + (np.sin(r) / r)[:, None] * transverse
        )
    derivative[:, 3] = gradient_s[1:-1, 1:-1, 1:-1].ravel()
    return energy, derivative.ravel()


def gradient_audit(
    base_field: np.ndarray,
    base_breathing: np.ndarray,
    spacing: float,
    model: Any,
) -> dict[str, float]:
    base, scalar, frames = fixed_chart(base_field, base_breathing)
    rng = np.random.default_rng(20260717)
    coordinates = 0.025 * rng.normal(size=4 * len(base))
    direction = rng.normal(size=len(coordinates))
    direction /= np.linalg.norm(direction)
    energy, gradient = chart_objective(
        coordinates, base_field, base_breathing, base, scalar, frames, spacing, model
    )
    predicted = float(gradient @ direction)
    rows = []
    for step in (2.0e-5, 1.0e-5, 5.0e-6):
        plus = chart_objective(
            coordinates + step * direction,
            base_field,
            base_breathing,
            base,
            scalar,
            frames,
            spacing,
            model,
        )[0]
        minus = chart_objective(
            coordinates - step * direction,
            base_field,
            base_breathing,
            base,
            scalar,
            frames,
            spacing,
            model,
        )[0]
        finite = (plus - minus) / (2.0 * step)
        rows.append(abs(finite - predicted) / max(abs(predicted), 1.0))
    return {
        "energy": energy,
        "analytic_directional_derivative": predicted,
        "best_relative_residual": min(rows),
    }


def permutation_sign(permutation: tuple[int, int, int]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(3)
        for j in range(i + 1, 3)
    )
    return -1 if inversions % 2 else 1


def geometric_preimage_degree(
    field: np.ndarray, targets: tuple[np.ndarray, ...] = TARGETS
) -> dict[str, Any]:
    """Degree of normalized affine tetrahedral interpolation by regular values.

    A target q is covered by a tetrahedron with vertex matrix N precisely when
    every component of N^{-1}q is positive.  The local map degree is
    sign(det N)/sign(det D).  Our baryon convention is minus the map degree.
    """

    rows = []
    minimum_nonzero_abs_det = math.inf
    boundary = np.concatenate(
        [
            field[0].reshape(-1, 4),
            field[-1].reshape(-1, 4),
            field[:, 0].reshape(-1, 4),
            field[:, -1].reshape(-1, 4),
            field[:, :, 0].reshape(-1, 4),
            field[:, :, -1].reshape(-1, 4),
        ]
    )
    for target in targets:
        degree = 0
        hits = 0
        hit_abs_determinants: list[float] = []
        hit_cone_margins: list[float] = []
        closest_boundary = float(np.min(np.linalg.norm(boundary - target, axis=1)))
        for origin in np.ndindex((field.shape[0] - 1,) * 3):
            origin_array = np.array(origin)
            for permutation in itertools.permutations((0, 1, 2)):
                vertices = [origin_array.copy()]
                running = origin_array.copy()
                for direction in permutation:
                    running = running.copy()
                    running[direction] += 1
                    vertices.append(running)
                matrix = np.column_stack(
                    [field[tuple(vertex)] for vertex in vertices]
                )
                determinant = float(np.linalg.det(matrix))
                if abs(determinant) < 1.0e-13:
                    continue
                minimum_nonzero_abs_det = min(
                    minimum_nonzero_abs_det, abs(determinant)
                )
                coefficients = np.linalg.solve(matrix, target)
                margin = float(np.min(coefficients))
                if margin > 1.0e-11:
                    hits += 1
                    hit_abs_determinants.append(abs(determinant))
                    hit_cone_margins.append(margin)
                    degree += int(math.copysign(1, determinant)) * permutation_sign(permutation)
        rows.append(
            {
                "target": target.tolist(),
                "preimage_tetrahedra": hits,
                "map_degree": degree,
                "baryon_number": -degree,
                "closest_boundary_vertex_distance": closest_boundary,
                "minimum_hit_absolute_determinant": min(
                    hit_abs_determinants, default=None
                ),
                "minimum_hit_positive_cone_margin": min(
                    hit_cone_margins, default=None
                ),
            }
        )
    return {
        "method": "Freudenthal tetrahedra plus normalized-affine regular-value preimage count",
        "target_rows": rows,
        "targets_agree": len({row["map_degree"] for row in rows}) == 1,
        "minimum_nonzero_absolute_tetrahedron_determinant": (
            minimum_nonzero_abs_det
            if math.isfinite(minimum_nonzero_abs_det)
            else None
        ),
    }


def tetrahedral_admissibility(field: np.ndarray) -> dict[str, Any]:
    """Check that normalized affine interpolation never meets the origin."""

    vertex_rows: list[np.ndarray] = []
    for origin in np.ndindex((field.shape[0] - 1,) * 3):
        origin_array = np.array(origin)
        for permutation in itertools.permutations((0, 1, 2)):
            vertices = [origin_array.copy()]
            running = origin_array.copy()
            for direction in permutation:
                running = running.copy()
                running[direction] += 1
                vertices.append(running)
            vertex_rows.append(np.stack([field[tuple(vertex)] for vertex in vertices]))
    values = np.stack(vertex_rows)  # tetrahedron, vertex, R4
    tetrahedra = len(values)
    distances = np.ones(tetrahedra)

    # For two unit vectors the closest point of their segment is the midpoint.
    for first, second in itertools.combinations(range(4), 2):
        dot = np.sum(values[:, first] * values[:, second], axis=1)
        distances = np.minimum(
            distances, np.sqrt(np.maximum(0.5 * (1.0 + dot), 0.0))
        )

    # Interior stationary point on every triangular face and on the tetrahedron.
    # The augmented KKT system is essential: G can be singular precisely when
    # zero lies in the relative interior, in which case G^+ 1 / (1^T G^+ 1)
    # is invalid (its denominator can vanish).  Boundary optima were already
    # included at lower active-set cardinality.
    for size in (3, 4):
        for subset in itertools.combinations(range(4), size):
            selected = values[:, subset, :]
            gram = np.einsum("tij,tkj->tik", selected, selected)
            kkt = np.zeros((tetrahedra, size + 1, size + 1))
            kkt[:, :size, :size] = gram
            kkt[:, :size, size] = 1.0
            kkt[:, size, :size] = 1.0
            rhs = np.zeros((tetrahedra, size + 1))
            rhs[:, size] = 1.0
            solution = np.einsum(
                "tij,tj->ti", np.linalg.pinv(kkt, rcond=1.0e-13), rhs
            )
            weights = solution[:, :size]
            residual = np.linalg.norm(
                np.einsum("tij,tj->ti", kkt, solution) - rhs, axis=1
            )
            valid = (residual < 1.0e-9) & (
                np.min(weights, axis=1) >= -1.0e-10
            )
            candidate = np.full(tetrahedra, np.inf)
            image = np.einsum("ti,tij->tj", weights, selected)
            candidate[valid] = np.linalg.norm(image[valid], axis=1)
            distances = np.minimum(distances, candidate)

    minimum = float(np.min(distances))
    below_threshold = int(np.count_nonzero(distances <= 1.0e-7))
    return {
        "definition": "minimum Euclidean distance from zero to every vertex convex hull",
        "tetrahedra": tetrahedra,
        "minimum_convex_hull_origin_distance": minimum,
        "tetrahedra_below_1e_minus_7": below_threshold,
        "normalized_affine_interpolation_admissible": below_threshold == 0,
    }


def regular_tetrahedron_dislocation_regression() -> dict[str, Any]:
    """A four-vector convex hull containing zero must have zero margin."""

    field = np.zeros((2, 2, 2, 4))
    field[..., 0] = 1.0
    vertices = np.array(
        [
            [1.0, 1.0, 1.0, 0.0],
            [1.0, -1.0, -1.0, 0.0],
            [-1.0, 1.0, -1.0, 0.0],
            [-1.0, -1.0, 1.0, 0.0],
        ]
    ) / math.sqrt(3.0)
    for coordinate, vector in zip(
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)), vertices
    ):
        field[coordinate] = vector
    result = tetrahedral_admissibility(field)
    result["expected_equal_weight_convex_combination_norm"] = float(
        np.linalg.norm(np.mean(vertices, axis=0))
    )
    return result


def exact_vacuum_boundary(
    field: np.ndarray, breathing: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    field = field.copy()
    breathing = breathing.copy()
    boundary = np.zeros(field.shape[:3], dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    field[boundary] = np.array([1.0, 0.0, 0.0, 0.0])
    breathing[boundary] = 1.0
    return field, breathing


def core_anchor_mask(extent: int) -> np.ndarray:
    """Central 3^3 n-core defining the declared discrete B=1 sector."""

    interior = extent - 2
    centre = interior // 2
    anchor = np.zeros((interior, interior, interior), dtype=bool)
    anchor[
        centre - 1 : centre + 2,
        centre - 1 : centre + 2,
        centre - 1 : centre + 2,
    ] = True
    return anchor.ravel()


def unconstrained_unwind_case(
    solution: Any, model: Any, extent: int, length: float, centre_pin: bool = False
) -> dict[str, Any]:
    """Unconstrained negative control: minimize the declared finite-site action."""

    field0, breathing0, axis, _ = E6.sampled_background(solution, extent, length)
    field0, breathing0 = exact_vacuum_boundary(field0, breathing0)
    spacing = float(axis[1] - axis[0])
    base, scalar, frames = fixed_chart(field0, breathing0)
    nodes = len(base)

    def objective(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        return chart_objective(
            coordinates,
            field0,
            breathing0,
            base,
            scalar,
            frames,
            spacing,
            model,
        )

    bounds: list[tuple[float | None, float | None]] | None = None
    free_mask = np.ones((nodes, 4), dtype=bool)
    if centre_pin:
        centre = (extent - 2) // 2
        centre_node = int(
            np.ravel_multi_index((centre, centre, centre), (extent - 2,) * 3)
        )
        bounds = []
        for node in range(nodes):
            bounds.extend(
                [(0.0, 0.0)] * 3 if node == centre_node else [(None, None)] * 3
            )
            bounds.append((None, None))
        free_mask[centre_node, :3] = False
    result = minimize(
        objective,
        np.zeros(4 * nodes),
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={
            "maxiter": 1800,
            "maxcor": 30,
            "maxls": 60,
            "ftol": 2.0e-15,
            "gtol": 2.0e-10,
        },
    )
    field, breathing, _, _ = chart_fields(
        result.x, field0, breathing0, base, scalar, frames
    )
    assembled = E6.assemble_declared_hessian(field, breathing, spacing, model)
    gradient = assembled["gradient"].reshape(-1, 4)
    free_density = float(
        np.linalg.norm(gradient[free_mask])
        / math.sqrt(int(np.count_nonzero(free_mask)))
        / spacing**3
    )
    initial_energy = objective(np.zeros(4 * nodes))[0]
    return {
        "N": extent,
        "L": length,
        "a": spacing,
        "centre_n_minus_one_pinned": centre_pin,
        "success": bool(result.success),
        "iterations": int(result.nit),
        "initial_energy": float(initial_energy),
        "final_energy": float(result.fun),
        "relative_energy_decrease": float((initial_energy - result.fun) / initial_energy),
        "final_free_projected_gradient_l2_density": free_density,
        "topology_initial": geometric_preimage_degree(field0),
        "topology_final": geometric_preimage_degree(field),
        "admissibility_final": tetrahedral_admissibility(field),
    }


def guarded_admissible_search(
    solution: Any,
    model: Any,
    extent: int = 7,
    length: float = 6.0,
    admissibility_floor: float = 0.05,
    maximum_steps: int = 180,
) -> dict[str, Any]:
    """Topology-preserving descent with an explicit normalized-affine guard."""

    field, breathing, axis, _ = E6.sampled_background(solution, extent, length)
    field, breathing = exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    initial_topology = geometric_preimage_degree(field)
    initial_admissibility = tetrahedral_admissibility(field)
    energy, gradient_n, gradient_s = energy_embedding_gradient(
        field, breathing, spacing, model
    )
    trace = [
        {
            "accepted_step": 0,
            "energy": float(energy),
            "admissibility_margin": initial_admissibility[
                "minimum_convex_hull_origin_distance"
            ],
            "baryon_numbers_at_iterate": [
                row["baryon_number"] for row in initial_topology["target_rows"]
            ],
        }
    ]
    rejected_guard = 0
    rejected_armijo = 0
    termination = "maximum_steps"
    last_trial_margin = None
    for accepted in range(1, maximum_steps + 1):
        base, scalar, frames = fixed_chart(field, breathing)
        tangent_gradient = np.einsum(
            "nij,ni->nj",
            frames,
            gradient_n[1:-1, 1:-1, 1:-1].reshape(-1, 4),
        )
        scalar_gradient = gradient_s[1:-1, 1:-1, 1:-1].ravel()
        direction = np.column_stack(
            (-tangent_gradient, -scalar_gradient)
        ).ravel()
        gradient_norm_sq = float(direction @ direction)
        gradient_density = float(
            math.sqrt(gradient_norm_sq / len(direction)) / spacing**3
        )
        if gradient_density < 2.0e-6:
            termination = "interior_stationary"
            break
        maximum_site_norm = float(
            np.max(np.linalg.norm(direction.reshape(-1, 4)[:, :3], axis=1))
        )
        step = min(0.12, 0.04 / max(maximum_site_norm, 1.0e-12))
        accepted_candidate = False
        for _ in range(40):
            if step < 1.0e-9:
                break
            coordinates = step * direction
            candidate_field, candidate_s, _, _ = chart_fields(
                coordinates, field, breathing, base, scalar, frames
            )
            candidate_energy = energy_embedding_gradient(
                candidate_field, candidate_s, spacing, model
            )[0]
            if candidate_energy > energy - 1.0e-4 * step * gradient_norm_sq:
                rejected_armijo += 1
                step *= 0.5
                continue
            topology = geometric_preimage_degree(candidate_field)
            admissibility = tetrahedral_admissibility(candidate_field)
            last_trial_margin = admissibility[
                "minimum_convex_hull_origin_distance"
            ]
            topology_ok = topology["targets_agree"] and all(
                row["baryon_number"] == 1 for row in topology["target_rows"]
            )
            admissible = (
                admissibility["normalized_affine_interpolation_admissible"]
                and last_trial_margin >= admissibility_floor
            )
            if not topology_ok or not admissible:
                rejected_guard += 1
                step *= 0.5
                continue
            field = candidate_field
            breathing = candidate_s
            energy, gradient_n, gradient_s = energy_embedding_gradient(
                field, breathing, spacing, model
            )
            trace.append(
                {
                    "accepted_step": accepted,
                    "energy": float(energy),
                    "step_size": step,
                    "projected_gradient_l2_density_before_step": gradient_density,
                    "admissibility_margin": last_trial_margin,
                    "baryon_numbers_at_iterate": [
                        row["baryon_number"] for row in topology["target_rows"]
                    ],
                }
            )
            accepted_candidate = True
            break
        if not accepted_candidate:
            termination = "admissibility_boundary_blocked_descent"
            break
    assembled = E6.assemble_declared_hessian(field, breathing, spacing, model)
    final_topology = geometric_preimage_degree(field)
    final_admissibility = tetrahedral_admissibility(field)
    return {
        "N": extent,
        "L": length,
        "a": spacing,
        "admissibility_floor": admissibility_floor,
        "termination": termination,
        "accepted_steps": len(trace) - 1,
        "rejected_by_guard": rejected_guard,
        "rejected_by_armijo": rejected_armijo,
        "initial_energy": trace[0]["energy"],
        "final_energy": float(assembled["energy"]),
        "relative_energy_decrease": float(
            (trace[0]["energy"] - assembled["energy"]) / trace[0]["energy"]
        ),
        "final_projected_gradient_l2_density": assembled[
            "gradient_l2_density"
        ],
        "initial_topology": initial_topology,
        "final_topology": final_topology,
        "initial_admissibility": initial_admissibility,
        "final_admissibility": final_admissibility,
        "last_rejected_trial_margin": last_trial_margin,
        "energy_monotone": all(
            trace[index + 1]["energy"] <= trace[index]["energy"] + 1.0e-12
            for index in range(len(trace) - 1)
        ),
        "trace": trace,
    }


def restricted_hessian_directional_audit(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
    assembled: dict[str, Any],
    free_mask: np.ndarray,
) -> dict[str, Any]:
    """Centred geodesic check with every fixed-core tangent entry set to zero."""

    rng = np.random.default_rng(20260718)
    direction = np.zeros(assembled["hessian"].shape[0])
    direction[free_mask.ravel()] = rng.normal(size=int(np.count_nonzero(free_mask)))
    direction /= np.linalg.norm(direction)
    predicted = float(direction @ (assembled["hessian"] @ direction))
    centre = E6.declared_lattice_energy(field, breathing, spacing, model)
    rows = []
    for step in (2.0e-3, 1.0e-3, 5.0e-4):
        plus_n, plus_s = E6.exponential_map_direction(
            field, breathing, assembled, direction, step
        )
        minus_n, minus_s = E6.exponential_map_direction(
            field, breathing, assembled, direction, -step
        )
        finite = (
            E6.declared_lattice_energy(plus_n, plus_s, spacing, model)
            - 2.0 * centre
            + E6.declared_lattice_energy(minus_n, minus_s, spacing, model)
        ) / step**2
        rows.append(
            {
                "step": step,
                "finite_difference": finite,
                "analytic_quadratic_form": predicted,
                "relative_residual": abs(finite - predicted)
                / max(abs(predicted), 1.0),
            }
        )
    return {
        "method": "sitewise S3 exponential map with fixed-core tangent coordinates identically zero",
        "rows": rows,
        "best_relative_residual": min(row["relative_residual"] for row in rows),
    }


def rerelax_case(
    solution: Any,
    model: Any,
    extent: int,
    length: float,
    maxiter: int = 1200,
) -> dict[str, Any]:
    field0, breathing0, axis, _ = E6.sampled_background(solution, extent, length)
    field0, breathing0 = exact_vacuum_boundary(field0, breathing0)
    spacing = float(axis[1] - axis[0])
    base, scalar, frames = fixed_chart(field0, breathing0)
    anchor = core_anchor_mask(extent)
    free_mask = np.ones((len(base), 4), dtype=bool)
    free_mask[anchor, :3] = False
    initial_energy, initial_gradient = chart_objective(
        np.zeros(4 * len(base)),
        field0,
        breathing0,
        base,
        scalar,
        frames,
        spacing,
        model,
    )
    trace: list[float] = [float(initial_energy)]

    def objective(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        return chart_objective(
            coordinates,
            field0,
            breathing0,
            base,
            scalar,
            frames,
            spacing,
            model,
        )

    def callback(coordinates: np.ndarray) -> None:
        trace.append(float(objective(coordinates)[0]))

    bounds: list[tuple[float | None, float | None]] = []
    for is_anchor in anchor:
        bounds.extend([(0.0, 0.0)] * 3 if is_anchor else [(None, None)] * 3)
        bounds.append((None, None))

    result = minimize(
        objective,
        np.zeros(4 * len(base)),
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        callback=callback,
        options={
            "maxiter": maxiter,
            "maxcor": 30,
            "maxls": 60,
            "ftol": 2.0e-15,
            "gtol": 2.0e-10,
        },
    )
    field, breathing, _, radii = chart_fields(
        result.x, field0, breathing0, base, scalar, frames
    )
    assembled = E6.assemble_declared_hessian(field, breathing, spacing, model)
    hessian_directional_audit = restricted_hessian_directional_audit(
        field, breathing, spacing, model, assembled, free_mask
    )
    gradient_blocks = assembled["gradient"].reshape(-1, 4)
    free_gradient_density = float(
        np.linalg.norm(gradient_blocks[free_mask])
        / math.sqrt(int(np.count_nonzero(free_mask)))
        / spacing**3
    )
    hessian_density = (assembled["hessian"] / spacing**3).tocsr()
    constrained_hessian = hessian_density[free_mask.ravel()][
        :, free_mask.ravel()
    ].tocsr()
    constrained_values, constrained_vectors = eigsh(
        constrained_hessian,
        k=min(6, constrained_hessian.shape[0] - 2),
        which="SA",
        tol=2.0e-8,
        maxiter=150000,
        v0=np.linspace(1.0, 2.0, constrained_hessian.shape[0]),
    )
    order = np.argsort(constrained_values)
    constrained_values = constrained_values[order]
    constrained_vectors = constrained_vectors[:, order]
    reconstructed = np.zeros(4 * len(base))
    reconstructed[free_mask.ravel()] = constrained_vectors[:, 0]
    node_amplitudes = np.linalg.norm(reconstructed.reshape(-1, 4), axis=1)
    inverse_participation = float(
        np.sum(node_amplitudes**4) / np.sum(node_amplitudes**2) ** 2
    )
    diquark_one = E6.scalar_dirichlet_laplacian(extent, spacing) + sparse.diags(
        model.diquark_mass_sq
        + model.diquark_core_coupling
        * (1.0 - field[1:-1, 1:-1, 1:-1, 0]).ravel(),
        format="csr",
    )
    diquark = sparse.kron(sparse.eye(2, format="csr"), diquark_one).tocsr()
    diquark_values = np.sort(
        eigsh(
            diquark_one,
            k=min(4, diquark_one.shape[0] - 2),
            which="SA",
            v0=np.linspace(1.0, 2.0, diquark_one.shape[0]),
            return_eigenvectors=False,
        )
    )
    laplacian = E6.scalar_dirichlet_laplacian(extent, spacing)
    gauge = sparse.kron(sparse.eye(3, format="csr"), laplacian).tocsr()
    ghost = laplacian.tocsr()
    gauge_min = float(
        eigsh(
            laplacian,
            k=1,
            which="SA",
            v0=np.linspace(1.0, 2.0, laplacian.shape[0]),
            return_eigenvectors=False,
        )[0]
    )
    fermion = declared_wilson_yukawa(field, spacing)
    normal = (fermion.getH() @ fermion).tocsr()
    fermion_min = float(
        eigsh(
            normal,
            k=1,
            which="SA",
            tol=2.0e-8,
            maxiter=100000,
            v0=np.linspace(1.0, 2.0, normal.shape[0]),
            return_eigenvectors=False,
        )[0]
    )
    topology_initial = geometric_preimage_degree(field0)
    topology_final = geometric_preimage_degree(field)
    admissibility_initial = tetrahedral_admissibility(field0)
    admissibility_final = tetrahedral_admissibility(field)
    final_energy = float(assembled["energy"])
    monotonic_tolerance = 2.0e-11 * max(abs(final_energy), 1.0)
    return {
        "N": extent,
        "L": length,
        "a": spacing,
        "boundary_condition": "exact vacuum n=(1,0,0,0), s=1 Dirichlet data on all six faces",
        "sector_constraint": {
            "definition": "central 3x3x3 AP-E6 hedgehog n values fixed; s remains dynamical there",
            "anchored_n_sites": int(np.count_nonzero(anchor)),
            "anchored_n_tangent_dofs": int(3 * np.count_nonzero(anchor)),
            "free_constrained_dofs": int(np.count_nonzero(free_mask)),
            "full_unanchored_stationarity_tested": False,
        },
        "optimizer": {
            "method": "fixed-base sitewise exponential chart plus analytic-gradient L-BFGS-B",
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
            "initial_energy": float(initial_energy),
            "final_energy": final_energy,
            "relative_energy_decrease": float((initial_energy - final_energy) / initial_energy),
            "initial_chart_gradient_l2": float(np.linalg.norm(initial_gradient)),
            "final_free_projected_gradient_l2_density": free_gradient_density,
            "full_unanchored_projected_gradient_l2_density": assembled[
                "gradient_l2_density"
            ],
            "maximum_chart_radius": float(np.max(radii)),
            "trace_samples": len(trace),
            "maximum_trace_energy_increase": float(
                max([trace[index + 1] - trace[index] for index in range(len(trace) - 1)] + [0.0])
            ),
            "monotonic_with_tolerance": bool(
                all(trace[index + 1] <= trace[index] + monotonic_tolerance for index in range(len(trace) - 1))
            ),
        },
        "topology_initial": topology_initial,
        "topology_final": topology_final,
        "admissibility_initial": admissibility_initial,
        "admissibility_final": admissibility_final,
        "continuum_derivative_B_after_relax": E6.lattice_baryon_number(field, axis),
        "n_s": {
            "unanchored_reduced_dofs": assembled["reduced_dofs"],
            "constraint_projected_dofs": int(constrained_hessian.shape[0]),
            "embedding_dofs": assembled["embedding_dofs"],
            "unanchored_nnz": assembled["nnz"],
            "constraint_projected_nnz": int(constrained_hessian.nnz),
            "symmetry_residual": assembled["symmetry_residual"],
            "constraint_projected_lowest": constrained_values.tolist(),
            "lowest_mode_inverse_participation_ratio": inverse_participation,
            "collective_projection_applied": False,
            "reason": "exact-vacuum boundary and fixed core explicitly break translations and isorotations",
            "directional_second_difference_audit": hessian_directional_audit,
        },
        "diquark": {
            "real_dofs": int(diquark.shape[0]),
            "nnz": int(diquark.nnz),
            "lowest": diquark_values.tolist(),
        },
        "gauge": {
            "gauge_condition": "xi=1 Feynman gauge",
            "real_dofs": int(gauge.shape[0]),
            "nnz": int(gauge.nnz),
            "lowest": gauge_min,
        },
        "ghost": {
            "complex_dofs": int(ghost.shape[0]),
            "nnz": int(ghost.nnz),
            "lowest": gauge_min,
        },
        "fermion_normal": {
            "complex_dofs": int(normal.shape[0]),
            "nnz": int(normal.nnz),
            "minimum_eigenvalue": fermion_min,
            "minimum_singular_value": math.sqrt(max(fermion_min, 0.0)),
            "operator": "K_F=D_W^dagger D_W; not the bosonic Hessian of -log det D_W",
        },
        "same_grid_declared_direct_sum_dofs": int(
            constrained_hessian.shape[0]
            + diquark.shape[0]
            + gauge.shape[0]
            + ghost.shape[0]
            + normal.shape[0]
        ),
        "same_grid_declared_direct_sum_nnz": int(
            constrained_hessian.nnz
            + diquark.nnz
            + gauge.nnz
            + ghost.nnz
            + normal.nnz
        ),
        "_field": field,
        "_breathing": breathing,
    }


def pauli() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    )


def declared_wilson_yukawa(field: np.ndarray, spacing: float) -> sparse.csr_matrix:
    """Declared 3D Dirichlet Wilson/Yukawa operator on C4 spin x C2 iso."""

    sigma = pauli()
    identity_two = np.eye(2, dtype=complex)
    gamma = [
        np.block([[np.zeros((2, 2)), -1j * item], [1j * item, np.zeros((2, 2))]])
        for item in sigma
    ]
    gamma5 = np.block(
        [[np.eye(2), np.zeros((2, 2))], [np.zeros((2, 2)), -np.eye(2)]]
    ).astype(complex)
    identity_four = np.eye(4, dtype=complex)
    pr = 0.5 * (identity_four + gamma5)
    pl = 0.5 * (identity_four - gamma5)
    extent = field.shape[0]
    interior = extent - 2
    nodes = interior**3
    one_shift_plus = sparse.diags(np.ones(interior - 1), 1, shape=(interior, interior), format="csr")
    one_shift_minus = one_shift_plus.T.tocsr()
    identity = sparse.eye(interior, format="csr")
    shifts_plus = [
        sparse.kron(sparse.kron(one_shift_plus, identity), identity, format="csr"),
        sparse.kron(sparse.kron(identity, one_shift_plus), identity, format="csr"),
        sparse.kron(sparse.kron(identity, identity), one_shift_plus, format="csr"),
    ]
    shifts_minus = [operator.T.tocsr() for operator in shifts_plus]
    node_identity = sparse.eye(nodes, format="csr")
    internal_identity = sparse.eye(8, format="csr", dtype=complex)
    mass = 0.60
    yukawa = 0.80
    wilson_r = 1.0
    operator = sparse.kron(node_identity, mass * internal_identity, format="csr")
    for direction in range(3):
        derivative = (shifts_plus[direction] - shifts_minus[direction]) / (2.0 * spacing)
        laplace_piece = (2.0 * node_identity - shifts_plus[direction] - shifts_minus[direction]) / spacing
        operator += sparse.kron(
            derivative,
            sparse.csr_matrix(np.kron(gamma[direction], identity_two)),
            format="csr",
        )
        operator += sparse.kron(
            laplace_piece,
            0.5 * wilson_r * internal_identity,
            format="csr",
        )
    interior_field = field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    for node, vector in enumerate(interior_field):
        unitary = vector[0] * identity_two
        for component in range(3):
            unitary += 1j * vector[component + 1] * sigma[component]
        local = yukawa * (np.kron(pr, unitary) + np.kron(pl, unitary.conj().T))
        for row, column in zip(*np.nonzero(np.abs(local) > 1.0e-15)):
            rows.append(8 * node + int(row))
            cols.append(8 * node + int(column))
            data.append(complex(local[row, column]))
    operator += sparse.coo_matrix(
        (data, (rows, cols)), shape=(8 * nodes, 8 * nodes), dtype=complex
    ).tocsr()
    return operator.tocsr()


def clean_case(case: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in case.items() if not key.startswith("_")}


def deterministic_digest(cases: list[dict[str, Any]]) -> str:
    payload = json.dumps(
        [clean_case(case) for case in cases],
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e7_discrete_rerelax_superhessian.json"
    md_path = OUTPUT / "ap_e7_discrete_rerelax_superhessian.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    unwind_rows = "\n".join(
        "| {N} | {L:.1f} | {a:.3f} | {energy:.8f} | {gradient:.3e} | {B:d} | {margin:.3e} |".format(
            N=row["N"],
            L=row["L"],
            a=row["a"],
            energy=row["final_energy"],
            gradient=row["final_free_projected_gradient_l2_density"],
            B=row["topology_final"]["target_rows"][0]["baryon_number"],
            margin=row["admissibility_final"]["minimum_convex_hull_origin_distance"],
        )
        for row in result["unconstrained_unwind_cases"]
    )
    guard_rows = "\n".join(
        "| {floor:.3f} | {steps} | {energy:.8f} | {gradient:.3e} | {margin:.9f} | {B:d} | {termination} |".format(
            floor=row["admissibility_floor"],
            steps=row["accepted_steps"],
            energy=row["final_energy"],
            gradient=row["final_projected_gradient_l2_density"],
            margin=row["final_admissibility"]["minimum_convex_hull_origin_distance"],
            B=row["final_topology"]["target_rows"][0]["baryon_number"],
            termination=row["termination"],
        )
        for row in result["admissibility_guard_sequence"]
    )
    anchor_rows = "\n".join(
        "| {N} | {L:.1f} | {a:.3f} | {iterations} | {energy:.8f} | {gradient:.3e} | {B:d} | {derivative_B:.6f} | {lambda0:.7f} | {delta0:.7f} | {sigma:.7f} |".format(
            N=row["N"],
            L=row["L"],
            a=row["a"],
            iterations=row["optimizer"]["iterations"],
            energy=row["optimizer"]["final_energy"],
            gradient=row["optimizer"]["final_free_projected_gradient_l2_density"],
            B=row["topology_final"]["target_rows"][0]["baryon_number"],
            derivative_B=row["continuum_derivative_B_after_relax"],
            lambda0=row["n_s"]["constraint_projected_lowest"][0],
            delta0=row["diquark"]["lowest"][0],
            sigma=row["fermion_normal"]["minimum_singular_value"],
        )
        for row in result["artificial_core_control_cases"]
    )
    checks = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - {row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    md_path.write_text(
        f"""# AP-E7 discrete re-relaxation and same-grid operator audit

- Status: `{result['status']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Finite-site configuration space path-connected: `{str(result['finite_site_configuration_space_path_connected']).lower()}`
- Exact unconstrained lattice B sector exists: `{str(result['exact_unconstrained_lattice_B_sector_exists']).lower()}`
- Unconstrained minimizers unwind to B=0: `{str(result['unconstrained_minimizers_unwind_to_B0']).lower()}`
- Admissible interior B=1 stationary point found: `{str(result['admissible_interior_B1_stationary_point_found']).lower()}`
- Artificial 3x3x3-core constrained solver converged: `{str(result['artificial_core_constrained_solver_converged']).lower()}`
- Exact n+s sparse Hessian in the artificial constrained sector: `{str(result['exact_n_s_sparse_hessian_artificial_sector']).lower()}`
- Same-grid declared proxy blocks assembled: `{str(result['same_grid_declared_proxy_blocks_assembled']).lower()}`
- Physical SU(2) gauge-ghost blocks assembled: `{str(result['physical_su2_gauge_ghost_blocks_assembled']).lower()}`
- Physical aggregate super-Hessian complete: `{str(result['physical_aggregate_superhessian_complete']).lower()}`
- 4D dynamics performed: `{str(result['four_dimensional_dynamics_performed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`
- Portal start allowed: `{str(result['portal_start_allowed']).lower()}`
- Lane closed: `{str(result['lane_closed']).lower()}`

## Unconstrained same-action minimization

| N | L | a | E_final | free gradient density | B_geom | admissibility margin |
|---:|---:|---:|---:|---:|---:|---:|
{unwind_rows}

## Admissible-component guard at N=7, L=6

| floor | accepted steps | E_final | gradient density | final margin | B_geom | termination |
|---:|---:|---:|---:|---:|---:|:---|
{guard_rows}

## Artificial fixed-core control only

| N | L | a | iters | E | free gradient density | B_geom | B_derivative | restricted n+s lambda0 | Delta lambda0 | sigma_min(D_W) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{anchor_rows}

For fixed boundary, the finite-site space is a finite product of connected
spaces.  It therefore has no exact pi_3 sector.  The regular-value degree is
locally constant only after excluding tetrahedra whose convex hull meets the
origin.  Unconstrained minimization crosses precisely that dislocation set.
Every accepted guard endpoint has B=1, but each search reaches its imposed
admissibility floor with a nonzero gradient.  The fixed 3x3x3 core is only an artificial
solver and sparse-operator control, never a genuine unanchored B=1 solution.
The Wilson/Yukawa normal operator is not the bosonic second variation of
`-log det D_W`; the physical aggregate and four-dimensional gates remain false.

## Mechanical checks

{checks}

## Fail-closed boundary

"""
        + "\n".join(
            f"- `{name}`: `{str(value).lower()}`" for name, value in result["fail_closed"].items()
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    model = E6.Model()
    source_manifest = [source_row(path) for path in (THIS_SCRIPT, E6_SCRIPT, TEX, BIB)]
    check(
        "P0_provenance",
        "AP-E7 source, imported AP-E6 action, TeX, and bibliography are present and hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )
    profile = E6.solve_relaxed_profile(box=16.0, model=model, sample_points=4097)
    checksum = profile["profile_sha256_r_F_s_little_endian_float64"]
    check(
        "P0_provenance",
        "the AP-E6 canonical continuum profile checksum is reproduced",
        checksum == "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59",
        f"checksum={checksum}",
    )
    solution = profile["solution"]
    audit_field, audit_s, audit_axis, _ = E6.sampled_background(solution, 7, 6.0)
    audit_field, audit_s = exact_vacuum_boundary(audit_field, audit_s)
    gradient_check = gradient_audit(audit_field, audit_s, float(audit_axis[1] - audit_axis[0]), model)
    check(
        "P1_relaxation",
        "the independent chart gradient agrees with a centred directional derivative",
        gradient_check["best_relative_residual"] < 2.0e-8,
        f"relative residual={gradient_check['best_relative_residual']:.3e}",
    )
    dislocation_regression = regular_tetrahedron_dislocation_regression()
    check(
        "P1_geometry",
        "a regular tetrahedron whose equal-weight convex combination is zero is rejected as non-admissible",
        dislocation_regression["minimum_convex_hull_origin_distance"] < 1.0e-12
        and not dislocation_regression[
            "normalized_affine_interpolation_admissible"
        ],
        "margin={:.3e}; admissible={}".format(
            dislocation_regression["minimum_convex_hull_origin_distance"],
            dislocation_regression["normalized_affine_interpolation_admissible"],
        ),
    )

    finite_site_path_connected = True
    check(
        "P2_configuration_space",
        "the fixed-boundary finite-site configuration space is path-connected",
        finite_site_path_connected,
        "C_N=(S^3)^(N_int) x R^(N_int); finite products of path-connected spaces are path-connected",
    )

    case_specs = [(7, 6.0), (9, 6.0), (11, 6.0), (9, 8.0), (11, 10.0)]
    unwind_cases = [
        unconstrained_unwind_case(solution, model, extent, length)
        for extent, length in case_specs
    ]
    check(
        "P2_unwind",
        "unconstrained same-action minimization converges but unwinds every tested B=1 sample to B=0",
        all(
            row["success"]
            and row["final_free_projected_gradient_l2_density"] < 2.0e-6
            and row["topology_final"]["targets_agree"]
            and all(target["baryon_number"] == 0 for target in row["topology_final"]["target_rows"])
            for row in unwind_cases
        ),
        "max gradient={:.3e}; final B rows={}".format(
            max(row["final_free_projected_gradient_l2_density"] for row in unwind_cases),
            [row["topology_final"]["target_rows"][0]["baryon_number"] for row in unwind_cases],
        ),
    )
    centre_pin_control = unconstrained_unwind_case(
        solution, model, 7, 6.0, centre_pin=True
    )
    check(
        "P2_unwind",
        "pinning only the centre to n=-1 does not define a B=1 lattice sector",
        centre_pin_control["success"]
        and all(
            target["baryon_number"] == 0
            for target in centre_pin_control["topology_final"]["target_rows"]
        ),
        "final B=" + str(
            [target["baryon_number"] for target in centre_pin_control["topology_final"]["target_rows"]]
        ),
    )

    guard_sequence = [
        guarded_admissible_search(
            solution, model, admissibility_floor=floor, maximum_steps=120
        )
        for floor in (0.10, 0.05, 0.02)
    ]
    check(
        "P3_admissible_guard",
        "three regular targets give B=1 at every accepted endpoint-guard iterate",
        all(
            row["final_topology"]["targets_agree"]
            and all(target["baryon_number"] == 1 for target in row["final_topology"]["target_rows"])
            and row["energy_monotone"]
            and all(
                trace_row["baryon_numbers_at_iterate"] == [1, 1, 1]
                for trace_row in row["trace"]
            )
            for row in guard_sequence
        ),
        "B rows=" + str([[target["baryon_number"] for target in row["final_topology"]["target_rows"]] for row in guard_sequence]),
    )
    check(
        "P3_admissible_guard",
        "each guarded descent reaches its admissibility floor before an interior stationary point",
        all(
            row["termination"] == "admissibility_boundary_blocked_descent"
            and row["final_projected_gradient_l2_density"] > 1.0e-2
            and abs(
                row["final_admissibility"]["minimum_convex_hull_origin_distance"]
                - row["admissibility_floor"]
            )
            < 2.0e-7
            for row in guard_sequence
        ),
        "(floor,margin,gradient)="
        + str(
            [
                (
                    row["admissibility_floor"],
                    row["final_admissibility"]["minimum_convex_hull_origin_distance"],
                    row["final_projected_gradient_l2_density"],
                )
                for row in guard_sequence
            ]
        ),
    )

    anchored_cases = [
        rerelax_case(solution, model, extent, length)
        for extent, length in case_specs
    ]
    clean_cases = [clean_case(case) for case in anchored_cases]
    check(
        "P4_artificial_control",
        "each fixed-3x3x3-core control converges in its explicitly restricted sector",
        all(
            row["optimizer"]["success"]
            and row["optimizer"]["final_free_projected_gradient_l2_density"] < 2.0e-6
            and row["optimizer"]["monotonic_with_tolerance"]
            for row in clean_cases
        ),
        "max free gradient={:.3e}".format(
            max(row["optimizer"]["final_free_projected_gradient_l2_density"] for row in clean_cases)
        ),
    )
    check(
        "P4_artificial_control",
        "the artificial controls retain the three-target B=1 preimage count",
        all(
            row["topology_final"]["targets_agree"]
            and all(target["baryon_number"] == 1 for target in row["topology_final"]["target_rows"])
            for row in clean_cases
        ),
        "B rows=" + str([[target["baryon_number"] for target in row["topology_final"]["target_rows"]] for row in clean_cases]),
    )
    check(
        "P5_sparse_hessian",
        "each artificial-sector exact tangent n+s Hessian is symmetric and positive",
        all(
            row["n_s"]["symmetry_residual"] < 1.0e-12
            and row["n_s"]["constraint_projected_lowest"][0] > 0.0
            for row in clean_cases
        ),
        "lowest restricted values=" + str([row["n_s"]["constraint_projected_lowest"][0] for row in clean_cases]),
    )
    check(
        "P5_sparse_hessian",
        "a core-compatible geodesic second difference reproduces every analytic n+s Hessian",
        all(
            row["n_s"]["directional_second_difference_audit"]["best_relative_residual"]
            < 2.0e-7
            for row in clean_cases
        ),
        "best residuals="
        + str(
            [
                row["n_s"]["directional_second_difference_audit"]["best_relative_residual"]
                for row in clean_cases
            ]
        ),
    )
    check(
        "P5_same_grid_blocks",
        "diquark, gauge, ghost, and declared Wilson/Yukawa normal operators are positive on every identical grid",
        all(
            row["diquark"]["lowest"][0] > 0.0
            and row["gauge"]["lowest"] > 0.0
            and row["ghost"]["lowest"] > 0.0
            and row["fermion_normal"]["minimum_eigenvalue"] > -1.0e-10
            for row in clean_cases
        ),
        "minima Delta={:.6g}, gauge={:.6g}, sigma(D)={:.6g}".format(
            min(row["diquark"]["lowest"][0] for row in clean_cases),
            min(row["gauge"]["lowest"] for row in clean_cases),
            min(row["fermion_normal"]["minimum_singular_value"] for row in clean_cases),
        ),
    )
    check(
        "P5_same_grid_blocks",
        "every declared proxy has the advertised same-grid multiplicity",
        all(
            row["diquark"]["real_dofs"] == 2 * (row["N"] - 2) ** 3
            and row["gauge"]["real_dofs"] == 3 * (row["N"] - 2) ** 3
            and row["ghost"]["complex_dofs"] == (row["N"] - 2) ** 3
            and row["fermion_normal"]["complex_dofs"]
            == 8 * (row["N"] - 2) ** 3
            for row in clean_cases
        ),
        "multiplicities=(Delta 2M, abelian spatial gauge 3M, abelian ghost M, fermion 8M)",
    )

    digest_first = deterministic_digest(anchored_cases)
    repeat = rerelax_case(solution, model, 7, 6.0)
    digest_case_first = deterministic_digest([anchored_cases[0]])
    digest_case_repeat = deterministic_digest([repeat])
    check(
        "P5_determinism",
        "an independent repeat of the artificial N=7,L=6 control is byte-stable after canonical JSON serialization",
        digest_case_first == digest_case_repeat,
        f"first={digest_case_first}; repeat={digest_case_repeat}",
    )

    artificial_converged = all(
        row["optimizer"]["success"]
        and row["optimizer"]["final_free_projected_gradient_l2_density"] < 2.0e-6
        for row in clean_cases
    )
    unwind_to_b0 = all(
        row["success"]
        and row["final_free_projected_gradient_l2_density"] < 2.0e-6
        and row["topology_final"]["targets_agree"]
        and all(
            target["baryon_number"] == 0
            for target in row["topology_final"]["target_rows"]
        )
        for row in unwind_cases
    )
    guard_iterates_have_b1 = all(
        row["energy_monotone"]
        and all(
            trace_row["baryon_numbers_at_iterate"] == [1, 1, 1]
            for trace_row in row["trace"]
        )
        for row in guard_sequence
    )
    exact_ns = all(
        row["n_s"]["symmetry_residual"] < 1.0e-12
        and row["n_s"]["directional_second_difference_audit"][
            "best_relative_residual"
        ]
        < 2.0e-7
        for row in clean_cases
    )
    same_grid_proxies = all(
        row["diquark"]["real_dofs"] == 2 * (row["N"] - 2) ** 3
        and row["gauge"]["real_dofs"] == 3 * (row["N"] - 2) ** 3
        and row["ghost"]["complex_dofs"] == (row["N"] - 2) ** 3
        and row["fermion_normal"]["complex_dofs"] == 8 * (row["N"] - 2) ** 3
        and row["same_grid_declared_direct_sum_dofs"]
        == row["n_s"]["constraint_projected_dofs"]
        + row["diquark"]["real_dofs"]
        + row["gauge"]["real_dofs"]
        + row["ghost"]["complex_dofs"]
        + row["fermion_normal"]["complex_dofs"]
        and row["same_grid_declared_direct_sum_nnz"]
        == row["n_s"]["constraint_projected_nnz"]
        + row["diquark"]["nnz"]
        + row["gauge"]["nnz"]
        + row["ghost"]["nnz"]
        + row["fermion_normal"]["nnz"]
        and row["diquark"]["lowest"][0] > 0.0
        and row["gauge"]["lowest"] > 0.0
        and row["ghost"]["lowest"] > 0.0
        and row["fermion_normal"]["minimum_eigenvalue"] > -1.0e-10
        for row in clean_cases
    )
    fail_closed = {
        "genuine_unanchored_B1_lattice_stationary_point_found": False,
        "full_unanchored_discrete_stationarity_achieved": False,
        "admissibility_independent_continuum_B1_solution_constructed": False,
        "continuous_line_search_segments_certified_admissible": False,
        "physical_su2_gauge_ghost_blocks_assembled": False,
        "fermion_determinant_bosonic_second_variation_computed": False,
        "interacting_gauge_meson_cross_blocks_computed": False,
        "brst_superdeterminant_on_the_relaxed_soliton_computed": False,
        "physical_aggregate_superhessian_complete": False,
        "projected_stability_continuum_extrapolated": False,
        "four_dimensional_dynamics_performed": False,
        "importance_sampling_or_hmc_performed": False,
        "quantum_continuum_limit_proven": False,
    }
    lane_closed = False
    result = {
        "artifact": "AP-E7 discrete re-relaxation and same-grid operator audit",
        "status": "mechanical_pass_physics_fail_closed",
        "generated_utc": "deterministic-no-wall-clock",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "model": asdict(model),
        "source_manifest": source_manifest,
        "canonical_profile_checksum": checksum,
        "gradient_audit": gradient_check,
        "regular_tetrahedron_dislocation_regression": dislocation_regression,
        "case_specifications": case_specs,
        "unconstrained_unwind_cases": unwind_cases,
        "single_centre_pin_negative_control": centre_pin_control,
        "admissibility_guard_sequence": guard_sequence,
        "artificial_core_control_cases": clean_cases,
        "deterministic_all_cases_digest": digest_first,
        "deterministic_repeat_digest": digest_case_repeat,
        "projected_gradient_density_tolerance": 2.0e-6,
        "each_grid_independently_relaxed_from_canonical_continuum_initial_data": True,
        "common_exact_vacuum_dirichlet_boundary": True,
        "finite_site_configuration_space_path_connected": finite_site_path_connected,
        "exact_unconstrained_lattice_B_sector_exists": False,
        "unconstrained_minimizers_unwind_to_B0": unwind_to_b0,
        "single_centre_pin_preserves_B1": False,
        "admissibility_guard_accepted_iterates_have_B1": guard_iterates_have_b1,
        "admissible_interior_B1_stationary_point_found": False,
        "artificial_core_constrained_solver_converged": artificial_converged,
        "exact_n_s_sparse_hessian_artificial_sector": exact_ns,
        "same_grid_declared_proxy_blocks_assembled": same_grid_proxies,
        "same_grid_abelian_gauge_ghost_proxy_blocks_assembled": same_grid_proxies,
        "physical_su2_gauge_ghost_blocks_assembled": False,
        "declared_same_grid_wilson_yukawa_normal_operator_assembled": all(
            row["fermion_normal"]["complex_dofs"] == 8 * (row["N"] - 2) ** 3
            and row["fermion_normal"]["nnz"] > 0
            for row in clean_cases
        ),
        "operator_classification": {
            "n_plus_s": "exact constrained Hessian of the declared lattice energy, only after deleting fixed-core tangent rows and columns",
            "diquark": "exact Hessian of the separately declared charge-two quadratic Delta sector on the fixed background",
            "gauge": "free abelian xi=1 spatial-triplet proxy with 3M real dofs; not the 9M spatial-adjoint SU(2) block and no interacting cross blocks",
            "ghost": "free abelian Faddeev-Popov Grassmann proxy with M complex dofs; not the 3M-complex SU(2)-adjoint ghost",
            "fermion": "positive normal proxy K_F=D_W^dagger D_W; not the bosonic second variation of -log det D_W",
            "aggregate": "declared direct sum only; not a BRST super-Hessian and not a physical QC2D Hessian",
        },
        "physical_aggregate_superhessian_complete": False,
        "four_dimensional_dynamics_performed": False,
        "physics_promotion_allowed": False,
        "portal_start_allowed": False,
        "lane_closed": lane_closed,
        "fail_closed": fail_closed,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    if result["checks_passed"] != result["checks_total"]:
        result["status"] = "mechanical_failure_physics_fail_closed"
    write_outputs(result)
    print(
        f"AP-E7 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"status={result['status']}; lane_closed={lane_closed}"
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
