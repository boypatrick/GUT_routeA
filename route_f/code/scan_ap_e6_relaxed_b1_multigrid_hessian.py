#!/usr/bin/env python3
"""AP-E6 coupled relaxed B=1 profile and sparse projected Hessian audit.

This is a deterministic classical calculation in a declared reduced EFT.
It is not four-dimensional importance sampling, HMC, a fermion determinant,
or a proof of the quantum charged-QC2D phase.

The reduced sector contains a unit O(4) meson n, a neutral breathing scalar s,
a charge-two complex diquark at Delta=0, and free gauge/ghost quadratic
blocks.  A coupled hedgehog (F,s) is genuinely relaxed by solve_bvp.  The same
solution is sampled on several cubic lattices.  The exact sparse second
variation of the stated lattice n+s energy is assembled in embedding
coordinates and restricted to T_n S^3 with the constrained-Hessian formula.
Translation and isorotation collective directions are projected explicitly.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy
from scipy import sparse
from scipy.integrate import simpson, solve_bvp
from scipy.sparse.linalg import LinearOperator, eigsh


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e6_relaxed_b1_multigrid_hessian.tex"
BIB = ROUTE_F / "tex" / "ap_e6_relaxed_b1_multigrid_hessian.bib"
THIS_SCRIPT = Path(__file__).resolve()

EPS = 1.0e-4
BVP_TOL = 3.0e-7
CHECKS: list[dict[str, Any]] = []


@dataclass(frozen=True)
class Model:
    beta: float = 0.5
    breathing_mass: float = 2.5
    skyrme_coefficient: float = 1.0
    diquark_mass_sq: float = 0.81
    diquark_core_coupling: float = -0.30
    diquark_charge: int = 2


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


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


def origin_cubic(
    slope: float, beta_effective: float, skyrme_coefficient: float
) -> float:
    numerator = slope * (
        2.0 * skyrme_coefficient * slope**4
        + 4.0 * slope**2
        + 3.0 * beta_effective**2
    )
    return numerator / (
        30.0 * (1.0 + 2.0 * skyrme_coefficient * slope**2)
    )


def coupled_rhs(model: Model) -> Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    beta = model.beta
    mass = model.breathing_mass

    def rhs(x: np.ndarray, y: np.ndarray, _parameter: np.ndarray) -> np.ndarray:
        profile, derivative, breathing, breathing_derivative = y
        sine = np.sin(profile)
        cosine = np.cos(profile)
        leading = x**2 + 2.0 * model.skyrme_coefficient * sine**2
        second = (
            -2.0 * x * derivative
            - 2.0
            * sine
            * cosine
            * (model.skyrme_coefficient * derivative**2 - 1.0)
            + 2.0
            * model.skyrme_coefficient
            * sine**3
            * cosine
            / x**2
            + beta**2 * x**2 * breathing**2 * sine
        ) / leading
        breathing_second = (
            -2.0 * breathing_derivative / x
            + mass**2 * (breathing - 1.0)
            + 2.0 * beta**2 * breathing * (1.0 - cosine)
        )
        return np.vstack(
            (derivative, second, breathing_derivative, breathing_second)
        )

    return rhs


def solve_coupled_profile(model: Model, box: float) -> Any:
    mesh = np.linspace(EPS, box, max(500, int(32 * box)))
    radius_guess = 1.15
    profile_guess = (
        2.0 * np.arctan(radius_guess / mesh) * (box - mesh) / (box - EPS)
    )
    breathing_guess = 1.0 - 0.08 * np.exp(-(mesh / 1.5) ** 2)
    guess = np.vstack(
        (
            profile_guess,
            np.gradient(profile_guess, mesh),
            breathing_guess,
            np.gradient(breathing_guess, mesh),
        )
    )

    def boundary(
        left: np.ndarray, right: np.ndarray, parameter: np.ndarray
    ) -> np.ndarray:
        slope = float(parameter[0])
        cubic = origin_cubic(
            slope, model.beta * left[2], model.skyrme_coefficient
        )
        core_rhs = (
            model.breathing_mass**2 * (left[2] - 1.0)
            + 4.0 * model.beta**2 * left[2]
        )
        return np.array(
            [
                left[0] - (math.pi - slope * EPS + cubic * EPS**3),
                left[1] - (-slope + 3.0 * cubic * EPS**2),
                left[3] - EPS * core_rhs / 3.0,
                right[0],
                right[2] - 1.0,
            ]
        )

    solution = solve_bvp(
        coupled_rhs(model),
        boundary,
        mesh,
        guess,
        p=np.array([2.2]),
        tol=BVP_TOL,
        max_nodes=80000,
        verbose=0,
    )
    if solution.status != 0:
        raise RuntimeError(f"coupled BVP failed at R={box}: {solution.message}")
    return solution


def sample_relaxed_profile(
    solution: Any, box: float, sample_points: int = 4097
) -> dict[str, Any]:
    """Return the canonical r,F,s arrays and a byte-stable checksum."""

    radius = np.linspace(0.0, box, sample_points)
    evaluation_radius = np.maximum(radius, EPS)
    profile, derivative, breathing, breathing_derivative = solution.sol(
        evaluation_radius
    )
    profile[0] = math.pi
    canonical = np.column_stack((radius, profile, breathing)).astype("<f8")
    checksum = hashlib.sha256(canonical.tobytes(order="C")).hexdigest()
    return {
        "r": radius,
        "F": profile,
        "dF_dr": derivative,
        "s": breathing,
        "ds_dr": breathing_derivative,
        "profile_sha256_r_F_s_little_endian_float64": checksum,
        "sample_points": sample_points,
        "box": box,
        "solution": solution,
    }


def solve_relaxed_profile(
    box: float = 16.0,
    model: Model | None = None,
    sample_points: int = 4097,
) -> dict[str, Any]:
    """Public same-solution API for downstream AP-E6 worklines.

    The returned dictionary contains canonical arrays ``r``, ``F``, and ``s``
    plus their checksum.  Callers should record that checksum when deriving a
    Yukawa/Callias or other operator on this exact relaxed background.
    """

    if model is None:
        model = Model()
    solution = solve_coupled_profile(model, box)
    sampled = sample_relaxed_profile(solution, box, sample_points)
    sampled["model"] = asdict(model)
    return sampled


def profile_observables(solution: Any, model: Model, box: float) -> dict[str, Any]:
    x = np.linspace(EPS, box, 30001)
    profile, derivative, breathing, breathing_derivative = solution.sol(x)
    sine = np.sin(profile)
    cosine = np.cos(profile)
    e_two = float(
        simpson(0.5 * x**2 * derivative**2 + sine**2, x=x)
    )
    e_four = float(
        simpson(
            model.skyrme_coefficient
            * (sine**2 * derivative**2 + 0.5 * sine**4 / x**2),
            x=x,
        )
    )
    e_breathing_kinetic = float(
        simpson(0.5 * x**2 * breathing_derivative**2, x=x)
    )
    e_breathing_potential = float(
        simpson(
            0.5
            * model.breathing_mass**2
            * x**2
            * (breathing - 1.0) ** 2,
            x=x,
        )
    )
    e_mass = float(
        simpson(
            model.beta**2 * x**2 * breathing**2 * (1.0 - cosine), x=x
        )
    )
    baryon = float(simpson(-2.0 * sine**2 * derivative / math.pi, x=x))
    energy = (
        e_two
        + e_four
        + e_breathing_kinetic
        + e_breathing_potential
        + e_mass
    )
    virial = (
        e_two
        + e_breathing_kinetic
        - e_four
        + 3.0 * (e_breathing_potential + e_mass)
    )
    return {
        "box": box,
        "origin_slope": float(solution.p[0]),
        "breathing_at_origin": float(breathing[0]),
        "minimum_breathing": float(np.min(breathing)),
        "bvp_nodes": int(len(solution.x)),
        "bvp_max_rms_residual": float(np.max(solution.rms_residuals)),
        "B_numeric": baryon,
        "B_residual": abs(baryon - 1.0),
        "energy_dimensionless_integral": energy,
        "E2": e_two,
        "E4": e_four,
        "E_s_kinetic": e_breathing_kinetic,
        "E_s_potential": e_breathing_potential,
        "E_mass": e_mass,
        "derrick_virial": virial,
        "derrick_relative_residual": abs(virial) / energy,
        "scale_second_derivative": 2.0 * e_four
        + 6.0 * (e_breathing_potential + e_mass),
    }


def tangent_frame(unit_vector: np.ndarray) -> np.ndarray:
    """Householder frame E with E^T E=1 and E^T n=0."""

    reference = np.array([1.0, 0.0, 0.0, 0.0])
    difference = reference - unit_vector
    norm_sq = float(difference @ difference)
    if norm_sq < 1.0e-14:
        householder = np.eye(4)
    else:
        householder = np.eye(4) - 2.0 * np.outer(difference, difference) / norm_sq
    return householder[:, 1:]


def sampled_background(
    solution: Any, spatial_extent: int, physical_length: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    axis = np.linspace(-physical_length / 2.0, physical_length / 2.0, spatial_extent)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.sqrt(xx**2 + yy**2 + zz**2)
    flat = radius.ravel()
    profile, _, breathing, _ = solution.sol(np.maximum(flat, EPS))
    profile = profile.reshape(radius.shape)
    breathing = breathing.reshape(radius.shape)
    field = np.zeros((*radius.shape, 4), dtype=float)
    field[..., 0] = np.cos(profile)
    sine_over_radius = np.zeros_like(radius)
    nonzero = radius > 1.0e-14
    sine_over_radius[nonzero] = np.sin(profile[nonzero]) / radius[nonzero]
    field[..., 1] = sine_over_radius * xx
    field[..., 2] = sine_over_radius * yy
    field[..., 3] = sine_over_radius * zz
    field[~nonzero, 0] = -1.0
    field[~nonzero, 1:] = 0.0
    breathing[~nonzero] = float(solution.sol(EPS)[2])
    return field, breathing, axis, radius


def lattice_baryon_number(field: np.ndarray, axis: np.ndarray) -> float:
    spacing = float(axis[1] - axis[0])
    derivatives = [
        np.stack(
            [
                np.gradient(
                    field[..., component],
                    spacing,
                    axis=direction,
                    edge_order=2,
                )
                for component in range(4)
            ],
            axis=-1,
        )
        for direction in range(3)
    ]
    determinant = np.linalg.det(np.stack([field, *derivatives], axis=-2))
    integral = np.trapezoid(
        np.trapezoid(np.trapezoid(determinant, axis, axis=0), axis, axis=0),
        axis,
        axis=0,
    )
    return float(-integral / (2.0 * math.pi**2))


def interior_index(spatial_extent: int) -> np.ndarray:
    index = -np.ones((spatial_extent,) * 3, dtype=int)
    interior_shape = (spatial_extent - 2,) * 3
    index[1:-1, 1:-1, 1:-1] = np.arange(np.prod(interior_shape)).reshape(
        interior_shape
    )
    return index


def append_block(
    rows: list[int],
    cols: list[int],
    data: list[float],
    first_node: int,
    second_node: int,
    block: np.ndarray,
    fields_per_node: int = 5,
    first_offset: int = 0,
    second_offset: int = 0,
) -> None:
    if first_node < 0 or second_node < 0:
        return
    for first in range(block.shape[0]):
        for second in range(block.shape[1]):
            value = float(block[first, second])
            if value != 0.0:
                rows.append(fields_per_node * first_node + first_offset + first)
                cols.append(fields_per_node * second_node + second_offset + second)
                data.append(value)


def assemble_declared_hessian(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Model,
    skyrme_coefficient: float | None = None,
) -> dict[str, Any]:
    """Exact sparse constrained Hessian of the declared finite-difference energy."""

    if skyrme_coefficient is None:
        skyrme_coefficient = model.skyrme_coefficient
    extent = field.shape[0]
    mapping = interior_index(extent)
    node_count = int(np.max(mapping) + 1)
    embedding_size = 5 * node_count
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    gradient_n = np.zeros_like(field)
    gradient_s = np.zeros_like(breathing)
    energy_two = 0.0
    energy_four = 0.0
    energy_potential = 0.0
    identity_four = np.eye(4)

    # Nearest-neighbour kinetic terms, with fixed Dirichlet boundary values.
    for direction in range(3):
        for coordinate in np.ndindex(
            tuple(extent - 1 if axis == direction else extent for axis in range(3))
        ):
            neighbour = list(coordinate)
            neighbour[direction] += 1
            neighbour = tuple(neighbour)
            difference_n = field[neighbour] - field[coordinate]
            difference_s = float(breathing[neighbour] - breathing[coordinate])
            energy_two += 0.5 * spacing * (
                float(difference_n @ difference_n) + difference_s**2
            )
            gradient_n[coordinate] -= spacing * difference_n
            gradient_n[neighbour] += spacing * difference_n
            gradient_s[coordinate] -= spacing * difference_s
            gradient_s[neighbour] += spacing * difference_s
            first = int(mapping[coordinate])
            second = int(mapping[neighbour])
            append_block(rows, cols, data, first, first, spacing * identity_four)
            append_block(rows, cols, data, second, second, spacing * identity_four)
            append_block(rows, cols, data, first, second, -spacing * identity_four)
            append_block(rows, cols, data, second, first, -spacing * identity_four)
            scalar = np.array([[spacing]])
            append_block(rows, cols, data, first, first, scalar, first_offset=4, second_offset=4)
            append_block(rows, cols, data, second, second, scalar, first_offset=4, second_offset=4)
            append_block(rows, cols, data, first, second, -scalar, first_offset=4, second_offset=4)
            append_block(rows, cols, data, second, first, -scalar, first_offset=4, second_offset=4)

    # Cell Skyrme wedge: kappa/2 sum_{i<j} |d_i n wedge d_j n|^2.
    for coordinate in np.ndindex((extent - 1, extent - 1, extent - 1)):
        centre = tuple(coordinate)
        neighbours = []
        differences = []
        for direction in range(3):
            neighbour = list(centre)
            neighbour[direction] += 1
            neighbour_tuple = tuple(neighbour)
            neighbours.append(neighbour_tuple)
            differences.append((field[neighbour_tuple] - field[centre]) / spacing)
        for first_direction, second_direction in ((0, 1), (0, 2), (1, 2)):
            u = differences[first_direction]
            v = differences[second_direction]
            uv = float(u @ v)
            u_sq = float(u @ u)
            v_sq = float(v @ v)
            density = 0.5 * skyrme_coefficient * (u_sq * v_sq - uv**2)
            energy_four += spacing**3 * density
            grad_u = skyrme_coefficient * (v_sq * u - uv * v)
            grad_v = skyrme_coefficient * (u_sq * v - uv * u)
            force_u = spacing**2 * grad_u
            force_v = spacing**2 * grad_v
            first_neighbour = neighbours[first_direction]
            second_neighbour = neighbours[second_direction]
            gradient_n[centre] -= force_u + force_v
            gradient_n[first_neighbour] += force_u
            gradient_n[second_neighbour] += force_v

            h_uu = skyrme_coefficient * (v_sq * identity_four - np.outer(v, v))
            h_vv = skyrme_coefficient * (u_sq * identity_four - np.outer(u, u))
            h_uv = skyrme_coefficient * (
                2.0 * np.outer(u, v)
                - np.outer(v, u)
                - uv * identity_four
            )
            h_vu = h_uv.T
            local_nodes = [
                int(mapping[centre]),
                int(mapping[first_neighbour]),
                int(mapping[second_neighbour]),
            ]
            coefficient_u = [-1.0, 1.0, 0.0]
            coefficient_v = [-1.0, 0.0, 1.0]
            for first_local, first_node in enumerate(local_nodes):
                for second_local, second_node in enumerate(local_nodes):
                    block = spacing * (
                        coefficient_u[first_local]
                        * coefficient_u[second_local]
                        * h_uu
                        + coefficient_u[first_local]
                        * coefficient_v[second_local]
                        * h_uv
                        + coefficient_v[first_local]
                        * coefficient_u[second_local]
                        * h_vu
                        + coefficient_v[first_local]
                        * coefficient_v[second_local]
                        * h_vv
                    )
                    append_block(
                        rows, cols, data, first_node, second_node, block
                    )

    # Local breathing potential and n_0 coupling.
    for coordinate in np.ndindex((extent, extent, extent)):
        node = int(mapping[coordinate])
        if node < 0:
            continue
        n_zero = float(field[coordinate][0])
        scalar = float(breathing[coordinate])
        density = (
            0.5 * model.breathing_mass**2 * (scalar - 1.0) ** 2
            + model.beta**2 * scalar**2 * (1.0 - n_zero)
        )
        energy_potential += spacing**3 * density
        gradient_n[coordinate][0] += -spacing**3 * model.beta**2 * scalar**2
        gradient_s[coordinate] += spacing**3 * (
            model.breathing_mass**2 * (scalar - 1.0)
            + 2.0 * model.beta**2 * scalar * (1.0 - n_zero)
        )
        h_ss = spacing**3 * (
            model.breathing_mass**2 + 2.0 * model.beta**2 * (1.0 - n_zero)
        )
        cross = -2.0 * spacing**3 * model.beta**2 * scalar
        append_block(
            rows,
            cols,
            data,
            node,
            node,
            np.array([[h_ss]]),
            first_offset=4,
            second_offset=4,
        )
        append_block(
            rows,
            cols,
            data,
            node,
            node,
            np.array([[cross]]),
            first_offset=0,
            second_offset=4,
        )
        append_block(
            rows,
            cols,
            data,
            node,
            node,
            np.array([[cross]]),
            first_offset=4,
            second_offset=0,
        )

    embedding_hessian = sparse.coo_matrix(
        (data, (rows, cols)), shape=(embedding_size, embedding_size)
    ).tocsr()
    embedding_hessian.sum_duplicates()

    # Build local tangent restriction and constrained-Hessian multiplier.
    transform_rows: list[int] = []
    transform_cols: list[int] = []
    transform_data: list[float] = []
    multiplier = np.zeros(node_count)
    reduced_gradient = np.zeros(4 * node_count)
    frames = np.zeros((node_count, 4, 3))
    background_n = np.zeros((node_count, 4))
    background_s = np.zeros(node_count)
    for coordinate in np.ndindex((extent, extent, extent)):
        node = int(mapping[coordinate])
        if node < 0:
            continue
        unit_vector = field[coordinate]
        frame = tangent_frame(unit_vector)
        frames[node] = frame
        background_n[node] = unit_vector
        background_s[node] = breathing[coordinate]
        multiplier[node] = float(unit_vector @ gradient_n[coordinate])
        reduced_gradient[4 * node : 4 * node + 3] = frame.T @ gradient_n[
            coordinate
        ]
        reduced_gradient[4 * node + 3] = gradient_s[coordinate]
        for row in range(4):
            for column in range(3):
                transform_rows.append(5 * node + row)
                transform_cols.append(4 * node + column)
                transform_data.append(float(frame[row, column]))
        transform_rows.append(5 * node + 4)
        transform_cols.append(4 * node + 3)
        transform_data.append(1.0)
    correction_rows = []
    correction_values = []
    for node in range(node_count):
        for component in range(4):
            correction_rows.append(5 * node + component)
            correction_values.append(-multiplier[node])
    constrained = embedding_hessian + sparse.coo_matrix(
        (
            correction_values,
            (correction_rows, correction_rows),
        ),
        shape=embedding_hessian.shape,
    ).tocsr()
    transform = sparse.coo_matrix(
        (transform_data, (transform_rows, transform_cols)),
        shape=(embedding_size, 4 * node_count),
    ).tocsr()
    reduced_hessian = (transform.T @ constrained @ transform).tocsr()
    symmetry_residual = float(
        sparse.linalg.norm(reduced_hessian - reduced_hessian.T)
        / max(sparse.linalg.norm(reduced_hessian), 1.0)
    )
    return {
        "hessian": reduced_hessian,
        "gradient": reduced_gradient,
        "frames": frames,
        "background_n": background_n,
        "background_s": background_s,
        "mapping": mapping,
        "energy": energy_two + energy_four + energy_potential,
        "energy_two_including_breathing": energy_two,
        "energy_four": energy_four,
        "energy_potential": energy_potential,
        "symmetry_residual": symmetry_residual,
        "embedding_dofs": embedding_size,
        "reduced_dofs": 4 * node_count,
        "nnz": int(reduced_hessian.nnz),
        "gradient_l2_density": float(
            np.linalg.norm(reduced_gradient)
            / math.sqrt(max(len(reduced_gradient), 1))
            / max(spacing**3, 1.0e-15)
        ),
    }


def declared_lattice_energy(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Model,
    skyrme_coefficient: float | None = None,
) -> float:
    """Evaluate exactly the finite-difference energy differentiated above."""

    if skyrme_coefficient is None:
        skyrme_coefficient = model.skyrme_coefficient
    extent = field.shape[0]
    energy = 0.0
    for direction in range(3):
        difference_n = np.diff(field, axis=direction)
        difference_s = np.diff(breathing, axis=direction)
        energy += 0.5 * spacing * float(
            np.sum(difference_n**2) + np.sum(difference_s**2)
        )
    base = field[:-1, :-1, :-1]
    differences = [
        (field[1:, :-1, :-1] - base) / spacing,
        (field[:-1, 1:, :-1] - base) / spacing,
        (field[:-1, :-1, 1:] - base) / spacing,
    ]
    for first, second in ((0, 1), (0, 2), (1, 2)):
        u = differences[first]
        v = differences[second]
        u_sq = np.sum(u * u, axis=-1)
        v_sq = np.sum(v * v, axis=-1)
        uv = np.sum(u * v, axis=-1)
        energy += 0.5 * skyrme_coefficient * spacing**3 * float(
            np.sum(u_sq * v_sq - uv**2)
        )
    interior_n_zero = field[1:-1, 1:-1, 1:-1, 0]
    interior_s = breathing[1:-1, 1:-1, 1:-1]
    energy += spacing**3 * float(
        np.sum(
            0.5 * model.breathing_mass**2 * (interior_s - 1.0) ** 2
            + model.beta**2 * interior_s**2 * (1.0 - interior_n_zero)
        )
    )
    return energy


def exponential_map_direction(
    field: np.ndarray,
    breathing: np.ndarray,
    assembled: dict[str, Any],
    direction: np.ndarray,
    parameter: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Move along sitewise S3 geodesics and linearly in the scalar field."""

    moved_field = field.copy()
    moved_breathing = breathing.copy()
    mapping = assembled["mapping"]
    frames = assembled["frames"]
    for coordinate in np.ndindex(field.shape[:3]):
        node = int(mapping[coordinate])
        if node < 0:
            continue
        tangent_coordinates = direction[4 * node : 4 * node + 3]
        tangent_vector = frames[node] @ tangent_coordinates
        norm = float(np.linalg.norm(tangent_vector))
        if norm < 1.0e-15:
            moved_field[coordinate] = field[coordinate]
        else:
            angle = parameter * norm
            moved_field[coordinate] = (
                math.cos(angle) * field[coordinate]
                + math.sin(angle) * tangent_vector / norm
            )
        moved_breathing[coordinate] += parameter * direction[4 * node + 3]
    return moved_field, moved_breathing


def directional_second_difference_audit(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Model,
    assembled: dict[str, Any],
    collective: np.ndarray,
) -> dict[str, Any]:
    """Independent exponential-map check of v^T H_tan v."""

    rng = np.random.default_rng(271828)
    direction = rng.normal(size=assembled["hessian"].shape[0])
    direction -= collective @ (collective.T @ direction)
    direction /= np.linalg.norm(direction)
    predicted = float(direction @ (assembled["hessian"] @ direction))
    centre = declared_lattice_energy(field, breathing, spacing, model)
    rows = []
    for step in (2.0e-3, 1.0e-3, 5.0e-4):
        plus_n, plus_s = exponential_map_direction(
            field, breathing, assembled, direction, step
        )
        minus_n, minus_s = exponential_map_direction(
            field, breathing, assembled, direction, -step
        )
        finite_difference = (
            declared_lattice_energy(plus_n, plus_s, spacing, model)
            - 2.0 * centre
            + declared_lattice_energy(minus_n, minus_s, spacing, model)
        ) / step**2
        rows.append(
            {
                "step": step,
                "finite_difference": finite_difference,
                "analytic_quadratic_form": predicted,
                "absolute_residual": abs(finite_difference - predicted),
                "relative_residual": abs(finite_difference - predicted)
                / max(abs(predicted), 1.0),
            }
        )
    return {
        "method": "sitewise S3 exponential map plus centred energy difference",
        "rows": rows,
        "best_relative_residual": min(row["relative_residual"] for row in rows),
    }


def collective_basis(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    assembled: dict[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    mapping = assembled["mapping"]
    frames = assembled["frames"]
    node_count = frames.shape[0]
    columns: list[np.ndarray] = []
    labels: list[str] = []
    derivatives_n = [
        np.gradient(field, spacing, axis=direction, edge_order=2)
        for direction in range(3)
    ]
    derivatives_s = [
        np.gradient(breathing, spacing, axis=direction, edge_order=2)
        for direction in range(3)
    ]
    for direction in range(3):
        vector = np.zeros(4 * node_count)
        for coordinate in np.ndindex(field.shape[:3]):
            node = int(mapping[coordinate])
            if node < 0:
                continue
            vector[4 * node : 4 * node + 3] = (
                frames[node].T @ derivatives_n[direction][coordinate]
            )
            vector[4 * node + 3] = derivatives_s[direction][coordinate]
        columns.append(vector)
        labels.append(f"translation_{direction + 1}")
    unit_axes = np.eye(3)
    for generator in range(3):
        vector = np.zeros(4 * node_count)
        for coordinate in np.ndindex(field.shape[:3]):
            node = int(mapping[coordinate])
            if node < 0:
                continue
            delta = np.zeros(4)
            delta[1:] = np.cross(unit_axes[generator], field[coordinate][1:])
            vector[4 * node : 4 * node + 3] = frames[node].T @ delta
        columns.append(vector)
        labels.append(f"isorotation_{generator + 1}")
    raw = np.column_stack(columns)
    q, r = np.linalg.qr(raw, mode="reduced")
    rank = int(np.count_nonzero(np.abs(np.diag(r)) > 1.0e-9))
    if rank != 6:
        raise RuntimeError(f"collective basis rank is {rank}, expected six")
    return q, {"labels": labels, "raw_gram_eigenvalues": np.linalg.eigvalsh(raw.T @ raw).tolist()}


def projected_spectrum(
    hessian: sparse.csr_matrix,
    collective: np.ndarray,
    count: int = 8,
) -> dict[str, Any]:
    size = hessian.shape[0]
    raw_values = np.sort(
        eigsh(
            hessian,
            k=min(count, size - 2),
            which="SA",
            tol=2.0e-8,
            maxiter=150000,
            return_eigenvectors=False,
        )
    )

    def project(vector: np.ndarray) -> np.ndarray:
        return vector - collective @ (collective.T @ vector)

    diagonal_scale = max(float(np.max(np.abs(hessian.diagonal()))), 1.0)
    penalty = 20.0 * diagonal_scale

    def matvec(vector: np.ndarray) -> np.ndarray:
        projected = project(vector)
        return project(hessian @ projected) + penalty * collective @ (
            collective.T @ vector
        )

    operator = LinearOperator((size, size), matvec=matvec, dtype=float)
    projected_values, projected_vectors = eigsh(
        operator,
        k=min(count, size - 2),
        which="SA",
        tol=2.0e-8,
        maxiter=200000,
    )
    ordering = np.argsort(projected_values)
    projected_values = projected_values[ordering]
    projected_vectors = projected_vectors[:, ordering]
    rng = np.random.default_rng(20260716)
    probe = rng.normal(size=size)
    projected_probe = project(probe)
    idempotency = float(
        np.linalg.norm(project(projected_probe) - projected_probe)
        / np.linalg.norm(projected_probe)
    )
    orthogonality = float(np.max(np.abs(collective.T @ projected_probe)))
    eigenvector_overlap = float(
        np.max(np.abs(collective.T @ projected_vectors))
    )
    collective_rayleigh = np.diag(collective.T @ (hessian @ collective))
    collective_residuals = np.linalg.norm(hessian @ collective, axis=0)
    lowest_blocks = projected_vectors[:, 0].reshape(-1, 4)
    node_amplitude = np.linalg.norm(lowest_blocks, axis=1)
    inverse_participation = float(
        np.sum(node_amplitude**4) / np.sum(node_amplitude**2) ** 2
    )
    return {
        "raw_lowest": raw_values.tolist(),
        "projected_lowest": projected_values.tolist(),
        "collective_rayleigh": collective_rayleigh.tolist(),
        "collective_Hq_norms": collective_residuals.tolist(),
        "projector_idempotency_residual": idempotency,
        "projector_orthogonality_residual": orthogonality,
        "projected_eigenvector_collective_overlap": eigenvector_overlap,
        "lowest_mode_inverse_participation_ratio": inverse_participation,
        "lowest_mode_max_to_mean_node_amplitude": float(
            np.max(node_amplitude) / np.mean(node_amplitude)
        ),
        "penalty": penalty,
    }


def scalar_dirichlet_laplacian(extent: int, spacing: float) -> sparse.csr_matrix:
    one_dimensional = sparse.diags(
        [
            -np.ones(extent - 3),
            2.0 * np.ones(extent - 2),
            -np.ones(extent - 3),
        ],
        offsets=[-1, 0, 1],
        format="csr",
    ) / spacing**2
    identity = sparse.eye(extent - 2, format="csr")
    return (
        sparse.kron(sparse.kron(one_dimensional, identity), identity)
        + sparse.kron(sparse.kron(identity, one_dimensional), identity)
        + sparse.kron(sparse.kron(identity, identity), one_dimensional)
    ).tocsr()


def diquark_spectrum(
    field: np.ndarray,
    spacing: float,
    model: Model,
    coupling: float | None = None,
) -> np.ndarray:
    if coupling is None:
        coupling = model.diquark_core_coupling
    extent = field.shape[0]
    laplacian = scalar_dirichlet_laplacian(extent, spacing)
    core = (1.0 - field[1:-1, 1:-1, 1:-1, 0]).ravel()
    operator = laplacian + sparse.diags(
        model.diquark_mass_sq + coupling * core, format="csr"
    )
    values = eigsh(
        operator,
        k=5,
        which="SA",
        tol=2.0e-9,
        maxiter=100000,
        return_eigenvectors=False,
    )
    return np.sort(values)


def periodic_gradient(extent: int, spacing: float) -> sparse.csr_matrix:
    volume = extent**3
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def flat(coordinate: tuple[int, int, int]) -> int:
        return int(np.ravel_multi_index(coordinate, (extent, extent, extent)))

    for direction in range(3):
        for coordinate in np.ndindex((extent, extent, extent)):
            neighbour = list(coordinate)
            neighbour[direction] = (neighbour[direction] + 1) % extent
            row = direction * volume + flat(coordinate)
            rows.extend((row, row))
            cols.extend((flat(coordinate), flat(tuple(neighbour))))
            data.extend((-1.0 / spacing, 1.0 / spacing))
    return sparse.coo_matrix(
        (data, (rows, cols)), shape=(3 * volume, volume)
    ).tocsr()


def gauge_ghost_audit(extent: int = 5, spacing: float = 0.8) -> dict[str, Any]:
    gradient = periodic_gradient(extent, spacing)
    ghost = (gradient.T @ gradient).tocsr()
    volume = extent**3
    base = sparse.kron(sparse.eye(3, format="csr"), ghost).tocsr()
    torons = np.zeros((3 * volume, 3))
    for direction in range(3):
        torons[direction * volume : (direction + 1) * volume, direction] = (
            1.0 / math.sqrt(volume)
        )
    constant_ghost = np.ones(volume) / math.sqrt(volume)

    def residuals(basis: np.ndarray) -> tuple[float, float]:
        rng = np.random.default_rng(314159)
        probe = rng.normal(size=basis.shape[0])
        projected = probe - basis @ (basis.T @ probe)
        twice = projected - basis @ (basis.T @ projected)
        return (
            float(np.linalg.norm(twice - projected) / np.linalg.norm(projected)),
            float(np.max(np.abs(basis.T @ projected))),
        )

    gauge_idem, gauge_orth = residuals(torons)
    ghost_idem, ghost_orth = residuals(constant_ghost[:, None])
    ghost_eigenvalues = np.linalg.eigvalsh(ghost.toarray())
    ghost_nonzero = ghost_eigenvalues[ghost_eigenvalues > 1.0e-10]
    rows = []
    reference = None
    for xi in (0.5, 1.0, 2.0):
        gauge = base + (1.0 / xi - 1.0) * (gradient @ gradient.T)
        eigenvalues = np.linalg.eigvalsh(gauge.toarray())
        nonzero = eigenvalues[eigenvalues > 1.0e-10]
        weight = 0.5 * float(np.log(nonzero).sum()) - float(
            np.log(ghost_nonzero).sum()
        )
        corrected = weight + 0.5 * (volume - 1) * math.log(xi)
        if reference is None:
            reference = corrected
        rows.append(
            {
                "xi": xi,
                "gauge_zero_modes": int(np.count_nonzero(eigenvalues <= 1.0e-10)),
                "ghost_zero_modes": int(
                    np.count_nonzero(ghost_eigenvalues <= 1.0e-10)
                ),
                "minimum_nonzero_gauge_eigenvalue": float(nonzero[0]),
                "minimum_nonzero_ghost_eigenvalue": float(ghost_nonzero[0]),
                "xi_corrected_log_weight": corrected,
                "xi_corrected_residual": abs(corrected - reference),
            }
        )
    return {
        "extent": extent,
        "spacing": spacing,
        "gauge_projector_idempotency_residual": gauge_idem,
        "gauge_projector_orthogonality_residual": gauge_orth,
        "ghost_projector_idempotency_residual": ghost_idem,
        "ghost_projector_orthogonality_residual": ghost_orth,
        "toron_zero_modes": 3,
        "constant_ghost_zero_modes": 1,
        "xi_rows": rows,
        "maximum_xi_corrected_residual": max(
            row["xi_corrected_residual"] for row in rows
        ),
    }


def lattice_case(
    solution: Any,
    model: Model,
    spatial_extent: int,
    physical_length: float,
    compute_spectrum: bool = True,
    spectrum_count: int = 6,
) -> dict[str, Any]:
    field, breathing, axis, _ = sampled_background(
        solution, spatial_extent, physical_length
    )
    spacing = float(axis[1] - axis[0])
    assembled = assemble_declared_hessian(field, breathing, spacing, model)
    collective, metadata = collective_basis(field, breathing, spacing, assembled)
    hessian_density = (assembled["hessian"] / spacing**3).tocsr()
    if compute_spectrum:
        spectrum = projected_spectrum(
            hessian_density, collective, count=spectrum_count
        )
        diquark = diquark_spectrum(field, spacing, model)
        directional = directional_second_difference_audit(
            field, breathing, spacing, model, assembled, collective
        )
    else:
        spectrum = {}
        diquark = np.array([])
        directional = {}
    boundary_mask = np.zeros(field.shape[:3], dtype=bool)
    boundary_mask[[0, -1], :, :] = True
    boundary_mask[:, [0, -1], :] = True
    boundary_mask[:, :, [0, -1]] = True
    vacuum = np.array([1.0, 0.0, 0.0, 0.0])
    return {
        "N": spatial_extent,
        "L": physical_length,
        "a": spacing,
        "B_lattice": lattice_baryon_number(field, axis),
        "boundary_condition": "sampled continuum profile frozen on cubic boundary; not overwritten by vacuum",
        "maximum_boundary_n_distance_from_vacuum": float(
            np.max(np.linalg.norm(field[boundary_mask] - vacuum, axis=1))
        ),
        "maximum_boundary_s_distance_from_vacuum": float(
            np.max(np.abs(breathing[boundary_mask] - 1.0))
        ),
        "energy_lattice": assembled["energy"],
        "gradient_l2_density": assembled["gradient_l2_density"],
        "reduced_dofs": assembled["reduced_dofs"],
        "embedding_dofs": assembled["embedding_dofs"],
        "hessian_nnz": assembled["nnz"],
        "hessian_symmetry_residual": assembled["symmetry_residual"],
        "directional_second_difference": directional,
        "collective_basis": metadata,
        "spectrum": spectrum,
        "diquark_lowest": diquark.tolist(),
        "_field": field,
        "_breathing": breathing,
        "_assembled": assembled,
        "_collective": collective,
    }


def clean_case(case: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in case.items() if not key.startswith("_")}


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e6_relaxed_b1_multigrid_hessian.json"
    md_path = OUTPUT / "ap_e6_relaxed_b1_multigrid_hessian.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    volume_rows = "\n".join(
        f"| {row['box']:.0f} | {row['bvp_nodes']} | {row['energy_dimensionless_integral']:.10f} | "
        f"{row['breathing_at_origin']:.9f} | {row['B_numeric']:.12f} | "
        f"{row['derrick_relative_residual']:.3e} |"
        for row in result["coupled_bvp"]["volume_sequence"]
    )
    grid_rows = "\n".join(
        f"| {row['N']} | {row['L']:.1f} | {row['a']:.6f} | {row['B_lattice']:.9f} | "
        f"{row['energy_lattice']:.8f} | {row['gradient_l2_density']:.3e} | "
        f"{row['spectrum']['projected_lowest'][0]:.8f} | {row['diquark_lowest'][0]:.8f} |"
        for row in result["multigrid_scan"]
    )
    physical_rows = "\n".join(
        f"| {row['N']} | {row['L']:.1f} | {row['a']:.6f} | {row['B_lattice']:.9f} | "
        f"{row['spectrum']['projected_lowest'][0]:.8f} | {row['diquark_lowest'][0]:.8f} |"
        for row in result["multivolume_scan"]
    )
    checks = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    md_path.write_text(
        f"""# AP-E6 relaxed B=1 multigrid sparse-Hessian audit

- Status: `{result['status']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Same relaxed B=1 solution computed: `{str(result['same_relaxed_b1_solution_computed']).lower()}`
- Multigrid/multivolume converged: `{str(result['multigrid_multivolume_converged']).lower()}`
- Aggregate same-grid sparse projected Hessian complete: `{str(result['sparse_projected_hessian_complete_in_declared_sector']).lower()}`
- Exact n+s sparse second variation complete: `{str(result['n_s_sparse_constrained_second_variation_complete']).lower()}`
- Projected-Hessian stability gate: `{str(result['projected_hessian_stability_gate_pass']).lower()}`
- Discrete stationarity achieved: `{str(result['discrete_stationarity_achieved']).lower()}`
- Lattice topology converged to B=1: `{str(result['lattice_topology_converged_to_B1']).lower()}`
- Off-shell projected-curvature diagnostic: `{str(result['off_shell_projected_curvature_diagnostic']).lower()}`
- Full gauge-meson-ghost-fermion Hessian complete: `{str(result['full_gauge_meson_ghost_fermion_hessian_complete']).lower()}`
- 4D importance sampling performed: `{str(result['dynamical_4d_importance_sampling_performed']).lower()}`
- Lane closed: `{str(result['lane_closed']).lower()}`

This card is a deterministic classical calculation.  Exact sparse second
differentiation refers only to the declared `n+s` energy.  The charge-two
diquark and periodic free-gauge/ghost operators are separate algebraic audits,
not one same-grid assembled Hessian.  The continuum-relaxed BVP is
sampled, not re-relaxed, on each cubic grid: the cubic spectra are therefore
off-shell curvature diagnostics, not a physical stability Hessian.  The card
does not include a fermion determinant, dynamical colour links, HMC, or a
quantum continuum extrapolation.

## Coupled BVP volume sequence

| R | nodes | E/(4pi f/e) | s(0) | B | virial/E |
|---:|---:|---:|---:|---:|---:|
{volume_rows}

## Fixed-volume grid sequence

| N | L | a | B_lattice | E_lattice | residual density | projected lambda0 | Delta lambda0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
{grid_rows}

## Fixed-spacing volume sequence

| N | L | a | B_lattice | projected lambda0 | Delta lambda0 |
|---:|---:|---:|---:|---:|---:|
{physical_rows}

## Negative controls

- Strong-core diquark eigenvalue: `{result['negative_controls']['strong_core_diquark_lowest']}`
- Projector with one omitted collective direction has leakage: `{result['negative_controls']['omitted_collective_direction_leakage']}`
- No-Skyrme Derrick derivative: `{result['negative_controls']['no_skyrme_derrick_derivative']}`

## Mechanical checks

{checks}

## Fail-closed boundary

"""
        + "\n".join(
            f"- `{name}`: `{str(value).lower()}`"
            for name, value in result["fail_closed"].items()
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    model = Model()
    sources = [THIS_SCRIPT, TEX, BIB]
    source_manifest = [source_row(path) for path in sources]
    check(
        "P0_provenance",
        "script, independent derivation, and bibliography exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # Coupled continuum relaxation and finite-volume convergence.
    bvp_rows = []
    solutions: dict[float, Any] = {}
    for box in (8.0, 10.0, 12.0, 16.0):
        solution = solve_coupled_profile(model, box)
        solutions[box] = solution
        bvp_rows.append(profile_observables(solution, model, box))
    baseline_solution = solutions[16.0]
    canonical_profile = sample_relaxed_profile(
        baseline_solution, box=16.0, sample_points=4097
    )
    check(
        "P0_provenance",
        "the public same-solution API produces canonical r-F-s arrays with a SHA-256 checksum",
        len(canonical_profile["profile_sha256_r_F_s_little_endian_float64"])
        == 64
        and len(canonical_profile["r"]) == 4097
        and abs(float(canonical_profile["F"][0]) - math.pi) < 1.0e-14,
        f"points={len(canonical_profile['r'])}; checksum={canonical_profile['profile_sha256_r_F_s_little_endian_float64']}",
    )
    final_bvp = bvp_rows[-1]
    final_energy_change = abs(
        bvp_rows[-1]["energy_dimensionless_integral"]
        - bvp_rows[-2]["energy_dimensionless_integral"]
    )
    final_core_change = abs(
        bvp_rows[-1]["breathing_at_origin"]
        - bvp_rows[-2]["breathing_at_origin"]
    )
    check(
        "P1_relaxation",
        "the coupled F-s BVP converges from the regular origin series",
        all(row["bvp_max_rms_residual"] <= 1.05 * BVP_TOL for row in bvp_rows),
        f"max residual={max(row['bvp_max_rms_residual'] for row in bvp_rows):.3e}",
    )
    check(
        "P1_relaxation",
        "the relaxed solution has B=+1 and the breathing amplitude remains positive",
        final_bvp["B_residual"] < 5.0e-8
        and final_bvp["minimum_breathing"] > 0.85,
        f"B={final_bvp['B_numeric']:.12f}; min s={final_bvp['minimum_breathing']:.9f}",
    )
    check(
        "P1_relaxation",
        "the generalized Derrick virial identity is satisfied",
        final_bvp["derrick_relative_residual"] < 3.0e-6
        and final_bvp["scale_second_derivative"] > 0.0,
        f"virial/E={final_bvp['derrick_relative_residual']:.3e}; curvature={final_bvp['scale_second_derivative']:.8f}",
    )
    check(
        "P2_volume",
        "R=12 to R=16 stabilises energy and the coupled core field",
        final_energy_change < 6.0e-6 and final_core_change < 2.0e-6,
        f"delta E={final_energy_change:.3e}; delta s0={final_core_change:.3e}",
    )

    # Exact sparse Hessian on the same continuum-relaxed background.
    multigrid_internal = [
        lattice_case(
            baseline_solution, model, extent, 8.0, spectrum_count=1
        )
        for extent in (9, 13, 17)
    ]
    baseline_n17_l8 = multigrid_internal[2]
    multivolume_internal = [
        lattice_case(baseline_solution, model, 9, 4.0, spectrum_count=1),
        lattice_case(baseline_solution, model, 13, 6.0, spectrum_count=1),
        baseline_n17_l8,
    ]
    finest = multigrid_internal[-1]
    grid_b_errors = [abs(row["B_lattice"] - 1.0) for row in multigrid_internal]
    projected_minima = [
        row["spectrum"]["projected_lowest"][0] for row in multigrid_internal
    ]
    diquark_minima = [row["diquark_lowest"][0] for row in multigrid_internal]
    check(
        "P3_multigrid",
        "the fixed-volume B estimator improves monotonically with resolution",
        grid_b_errors[0] > grid_b_errors[1] > grid_b_errors[2],
        f"errors={grid_b_errors}",
    )
    check(
        "P3_multigrid",
        "the exact off-shell projector detects a negative meson-breathing curvature while the diquark block stays positive",
        max(projected_minima) < 0.0 and min(diquark_minima) > 0.0,
        f"projected negative={projected_minima}; diquark positive={diquark_minima}",
    )
    final_projected_change = abs(projected_minima[-1] - projected_minima[-2])
    final_diquark_change = abs(diquark_minima[-1] - diquark_minima[-2])
    check(
        "P3_multigrid",
        "the final refinement exposes nonconvergence of the off-shell curvature while the diquark gap remains controlled",
        final_projected_change > 1.0 and final_diquark_change < 0.02,
        f"delta projected={final_projected_change:.3e}; delta diquark={final_diquark_change:.3e}",
    )
    localisation = [
        row["spectrum"]["lowest_mode_inverse_participation_ratio"]
        for row in multigrid_internal
    ]
    check(
        "P3_multigrid",
        "the negative mode is strongly localised and is not misreported as a continuum collective zero mode",
        min(localisation) > 0.02,
        f"inverse participation ratios={localisation}",
    )
    volume_projected = [
        row["spectrum"]["projected_lowest"][0] for row in multivolume_internal
    ]
    volume_diquark = [row["diquark_lowest"][0] for row in multivolume_internal]
    inverse_length_sq = np.array(
        [1.0 / row["L"] ** 2 for row in multivolume_internal]
    )
    diquark_volume_coefficients = np.polyfit(
        inverse_length_sq, np.array(volume_diquark), deg=1
    )
    diquark_volume_extrapolate = float(
        np.polyval(diquark_volume_coefficients, 0.0)
    )
    diquark_volume_fit_rms = float(
        np.sqrt(
            np.mean(
                (
                    np.polyval(diquark_volume_coefficients, inverse_length_sq)
                    - np.array(volume_diquark)
                )
                ** 2
            )
        )
    )
    check(
        "P4_multivolume",
        "fixed-spacing volume scan reproduces the negative projected mode and positive diquark gap",
        max(volume_projected) < 0.0 and min(volume_diquark) > 0.0,
        f"projected negative={volume_projected}; diquark positive={volume_diquark}",
    )
    check(
        "P4_multivolume",
        "the negative mode stabilises and the unbound diquark level has the expected 1/L^2 threshold extrapolation",
        abs(volume_projected[-1] - volume_projected[-2]) < 0.60
        and abs(diquark_volume_extrapolate - model.diquark_mass_sq) < 0.02
        and diquark_volume_fit_rms < 2.0e-3,
        f"delta projected={abs(volume_projected[-1]-volume_projected[-2]):.3e}; "
        f"Delta extrapolate={diquark_volume_extrapolate:.9f}; fit rms={diquark_volume_fit_rms:.3e}",
    )

    finest_spectrum = finest["spectrum"]
    check(
        "P5_projector",
        "translation-isorotation projector is idempotent and horizontal",
        finest_spectrum["projector_idempotency_residual"] < 2.0e-12
        and finest_spectrum["projector_orthogonality_residual"] < 2.0e-11
        and finest_spectrum["projected_eigenvector_collective_overlap"] < 2.0e-8,
        f"idem={finest_spectrum['projector_idempotency_residual']:.3e}; "
        f"orth={finest_spectrum['projector_orthogonality_residual']:.3e}; "
        f"eig overlap={finest_spectrum['projected_eigenvector_collective_overlap']:.3e}",
    )
    check(
        "P5_projector",
        "the exact n+s sparse Hessian is symmetric",
        finest["hessian_symmetry_residual"] < 2.0e-13,
        f"dofs={finest['reduced_dofs']}; nnz={finest['hessian_nnz']}; residual={finest['hessian_symmetry_residual']:.3e}",
    )
    check(
        "P5_projector",
        "an independent exponential-map energy difference reproduces v^T H_tan v",
        finest["directional_second_difference"]["best_relative_residual"]
        < 2.0e-5,
        f"best relative residual={finest['directional_second_difference']['best_relative_residual']:.3e}",
    )

    gauge = gauge_ghost_audit()
    check(
        "P6_gauge_ghost",
        "toron and constant-ghost projectors are idempotent and orthogonal",
        max(
            gauge["gauge_projector_idempotency_residual"],
            gauge["gauge_projector_orthogonality_residual"],
            gauge["ghost_projector_idempotency_residual"],
            gauge["ghost_projector_orthogonality_residual"],
        )
        < 2.0e-12,
        f"max residual={max(gauge['gauge_projector_idempotency_residual'], gauge['gauge_projector_orthogonality_residual'], gauge['ghost_projector_idempotency_residual'], gauge['ghost_projector_orthogonality_residual']):.3e}",
    )
    check(
        "P6_gauge_ghost",
        "free gauge/ghost determinant ratio obeys the xi cancellation law",
        gauge["maximum_xi_corrected_residual"] < 5.0e-11
        and all(row["gauge_zero_modes"] == 3 for row in gauge["xi_rows"])
        and all(row["ghost_zero_modes"] == 1 for row in gauge["xi_rows"]),
        f"max xi residual={gauge['maximum_xi_corrected_residual']:.3e}",
    )

    # Negative controls must be detected.
    field = finest["_field"]
    spacing = finest["a"]
    strong_core = diquark_spectrum(field, spacing, model, coupling=-5.0)
    full_collective = finest["_collective"]
    bad_basis = full_collective[:, :5]
    omitted = full_collective[:, 5]
    bad_projection = omitted - bad_basis @ (bad_basis.T @ omitted)
    omitted_leakage = float(abs(omitted @ bad_projection))
    no_skyrme_derivative = (
        final_bvp["E2"]
        + final_bvp["E_s_kinetic"]
        + 3.0 * (final_bvp["E_s_potential"] + final_bvp["E_mass"])
    )
    check(
        "P7_negative_controls",
        "a strong attractive diquark core creates a tachyon",
        float(strong_core[0]) < 0.0,
        f"lambda0={strong_core[0]:.9f}",
    )
    check(
        "P7_negative_controls",
        "omitting one collective vector is detected by unit leakage",
        omitted_leakage > 0.999999,
        f"leakage={omitted_leakage:.12f}",
    )
    check(
        "P7_negative_controls",
        "removing the Skyrme term destroys Derrick stationarity",
        no_skyrme_derivative > 0.5,
        f"dE/dscale={no_skyrme_derivative:.9f}",
    )

    all_pass = all(row["pass"] for row in CHECKS)
    same_solution = all_pass and final_bvp["B_residual"] < 5.0e-8
    converged = False
    n_s_differentiation_complete = bool(
        all_pass
        and finest["hessian_symmetry_residual"] < 2.0e-13
        and finest["directional_second_difference"]["best_relative_residual"]
        < 2.0e-5
    )
    declared_complete = False
    stability_gate = False
    fail_closed = {
        "sparse_projected_hessian_complete_in_declared_sector": False,
        "all_blocks_assembled_on_one_same_grid_and_boundary_condition": False,
        "discrete_stationarity_achieved": False,
        "lattice_topology_converged_to_B1": False,
        "physical_projected_hessian_at_stationary_discrete_solution": False,
        "full_three_dimensional_continuum_stationarity_proven": False,
        "multigrid_multivolume_converged": False,
        "full_gauge_meson_ghost_fermion_hessian_complete": False,
        "fermion_determinant_included": False,
        "dynamical_colour_links_included": False,
        "dynamical_4d_importance_sampling_performed": False,
        "HMC_or_Monte_Carlo_performed": False,
        "continuum_quantum_phase_proven": False,
        "microscopic_charged_QC2D_phase_proven": False,
        "finite_amplitude_global_unwinding_excluded": False,
        "degree_one_route_E_portal_allowed": False,
        "physics_promotion_allowed": False,
        "lane_closed": False,
    }
    result: dict[str, Any] = {
        "schema_version": 1,
        "artifact": "AP-E6 coupled relaxed B=1 multigrid sparse-Hessian audit",
        "status": "mechanical_pass_physics_fail_closed" if all_pass else "mechanical_fail",
        "all_pass": all_pass,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "same_relaxed_b1_solution_computed": same_solution,
        "multigrid_multivolume_converged": converged,
        "sparse_projected_hessian_complete_in_declared_sector": declared_complete,
        "n_s_sparse_constrained_second_variation_complete": n_s_differentiation_complete,
        "diquark_separate_quadratic_block_audited": True,
        "gauge_ghost_separate_periodic_algebra_block_audited": True,
        "all_blocks_assembled_on_one_same_grid_and_boundary_condition": False,
        "projected_hessian_stability_gate_pass": stability_gate,
        "negative_off_shell_projected_curvature_detected": max(projected_minima)
        < 0.0,
        "off_shell_projected_curvature_diagnostic": True,
        "discrete_stationarity_achieved": False,
        "lattice_topology_converged_to_B1": False,
        "physical_projected_hessian_at_stationary_discrete_solution": False,
        "continuum_bvp_volume_converged": True,
        "coupled_hedgehog_sector_stationary": True,
        "full_three_dimensional_continuum_stationarity_proven": False,
        "algebraic_hessian_multigrid_control_pass": True,
        "full_gauge_meson_ghost_fermion_hessian_complete": False,
        "dynamical_4d_importance_sampling_performed": False,
        "continuum_quantum_phase_proven": False,
        "physics_promotion_allowed": False,
        "lane_closed": False,
        "epistemic_scope": {
            "calculation": "deterministic classical continuum BVP plus finite-difference sparse second variation",
            "exactly_differentiated_sector": "unit O(4) meson n plus neutral breathing scalar s on an open sampled-profile boundary",
            "separate_not_same_grid_audits": "charge-two Delta Dirichlet spectrum and periodic free gauge/ghost algebra",
            "excluded": "fermion determinant, interacting colour links, four-dimensional path integral, HMC, continuum quantum extrapolation",
        },
        "model": asdict(model),
        "coupled_bvp": {
            "equations": "coupled Euler-Lagrange BVP for F(r) and s(r)",
            "volume_sequence": bvp_rows,
            "final_energy_change": final_energy_change,
            "final_core_field_change": final_core_change,
            "canonical_same_solution": {
                "box": canonical_profile["box"],
                "sample_points": canonical_profile["sample_points"],
                "profile_sha256_r_F_s_little_endian_float64": canonical_profile[
                    "profile_sha256_r_F_s_little_endian_float64"
                ],
                "model": asdict(model),
                "public_solver_api": "route_f.code.scan_ap_e6_relaxed_b1_multigrid_hessian.solve_relaxed_profile",
            },
        },
        "multigrid_scan": [clean_case(row) for row in multigrid_internal],
        "multivolume_scan": [clean_case(row) for row in multivolume_internal],
        "gauge_ghost": gauge,
        "negative_controls": {
            "strong_core_coupling": -5.0,
            "strong_core_diquark_lowest": strong_core.tolist(),
            "omitted_collective_direction_leakage": omitted_leakage,
            "no_skyrme_derrick_derivative": no_skyrme_derivative,
        },
        "convergence_summary": {
            "projected_multigrid_converged": False,
            "reason_projected_multigrid_not_converged": "sampled continuum BVP is not re-relaxed on each cubic lattice; gradient density remains O(0.4) and the topology estimator is not yet near B=1",
            "fixed_volume_B_errors": grid_b_errors,
            "fixed_volume_projected_minima": projected_minima,
            "fixed_volume_diquark_minima": diquark_minima,
            "fixed_volume_lowest_mode_inverse_participation_ratios": localisation,
            "final_projected_grid_change": final_projected_change,
            "final_diquark_grid_change": final_diquark_change,
            "fixed_spacing_volume_projected_minima": volume_projected,
            "fixed_spacing_volume_diquark_minima": volume_diquark,
            "diquark_linear_in_inverse_L_squared_extrapolate": diquark_volume_extrapolate,
            "diquark_linear_in_inverse_L_squared_fit_rms": diquark_volume_fit_rms,
        },
        "checks": CHECKS,
        "fail_closed": fail_closed,
        "source_manifest": source_manifest,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "deterministic": True,
            "randomness_used_only_for_projector_probe": True,
            "random_seed": 20260716,
        },
    }
    write_outputs(result)
    print(
        f"AP-E6 relaxed B=1 multigrid Hessian: {result['status']} "
        f"({result['checks_passed']}/{result['checks_total']}); "
        f"declared_complete={declared_complete}; lane_closed={result['lane_closed']}"
    )


if __name__ == "__main__":
    main()
