#!/usr/bin/env python3
"""AP-E10 compactness counterexample, cell formula, and quotient scan.

This card follows the AP-E9 fail-closed order.

1. It disproves compensated compactness for the one-corner forward Skyrme
   stencil by an exact two-periodic sphere-valued sequence.  The sequence
   converges uniformly to the vacuum but its forward topological current has
   a nonzero diffuse limit.  The same Skyrme bound does give an L^(4/3) bound
   on that current, so singular concentration is excluded; the surviving
   obstruction is an oscillatory stencil defect.
2. It constructs the finite-R periodic bond cell problem for the uniform
   six-tetrahedron and checkerboard five-tetrahedron graphs.  A direction-wise
   cycle-balance certificate and Jensen's inequality prove that the zero
   corrector is exact.  The raw AP-E9 quartics are therefore the homogenized
   densities and their nonproportionality survives homogenization.
3. It re-relaxes a smaller-a/larger-L sequence with a physical translation
   quotient and an optional dynamic-centering projection.  Numerical evidence
   cannot override the exact stencil counterexample, so the Hessian and
   determinant gates remain closed in every mode.
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
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy
from scipy import ndimage
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E9_SCRIPT = ROUTE_F / "code" / "scan_ap_e9_gamma_triangulation.py"

EPSILON = 0.01
STATIONARITY_TOLERANCE = 2.0e-6
PAIR_MARGIN_TOLERANCE = 1.0e-8
COUNTEREXAMPLE_ETA = 0.2
PROFILE_RADIUS_MAXIMUM = 3.5
PROFILE_BINS = 35
CHECKS: list[dict[str, Any]] = []


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E9 = load_module("ap_e9_for_ap_e10", E9_SCRIPT)
E8 = E9.E8
E7 = E9.E7
E6 = E9.E6
TRI = E9.TRI

TRIANGULATIONS = {
    spec.name: spec
    for spec in (TRI.UNIFORM_SPECS[0], *TRI.FIVE_TET_SPECS)
}


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


def wedge_square(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Squared exterior-product norm for arrays ending in R^4."""

    return np.sum(first * first, axis=-1) * np.sum(second * second, axis=-1) - np.sum(
        first * second, axis=-1
    ) ** 2


def forward_current(field: np.ndarray, spacing: float) -> dict[str, Any]:
    """One-corner AP-E7 derivatives, minors, and topological current."""

    base = field[:-1, :-1, :-1]
    derivatives = (
        (field[1:, :-1, :-1] - base) / spacing,
        (field[:-1, 1:, :-1] - base) / spacing,
        (field[:-1, :-1, 1:] - base) / spacing,
    )
    wedge_squares = [
        wedge_square(derivatives[first], derivatives[second])
        for first, second in ((0, 1), (0, 2), (1, 2))
    ]
    matrices = np.stack((base, *derivatives), axis=-1)
    current = np.linalg.det(matrices)
    cell_volume = spacing**3
    dirichlet = 0.5 * cell_volume * sum(
        float(np.sum(derivative**2)) for derivative in derivatives
    )
    skyrme = 0.5 * cell_volume * sum(
        float(np.sum(value)) for value in wedge_squares
    )
    current_l43_power = cell_volume * float(np.sum(np.abs(current) ** (4.0 / 3.0)))
    wedge_l2_power = cell_volume * sum(float(np.sum(value)) for value in wedge_squares)
    return {
        "dirichlet_energy_unit_coefficient": dirichlet,
        "skyrme_energy_unit_coefficient": skyrme,
        "forward_current_integral": cell_volume * float(np.sum(current)),
        "forward_current_minimum": float(np.min(current)),
        "forward_current_maximum": float(np.max(current)),
        "forward_current_L43_power": current_l43_power,
        "wedge2_L2_power": wedge_l2_power,
        "L43_bound_ratio": (
            3.0 * current_l43_power / wedge_l2_power if wedge_l2_power > 0.0 else 0.0
        ),
    }


def torus_counterexample(cells: int, eta: float = COUNTEREXAMPLE_ETA) -> dict[str, Any]:
    """Exact two-periodic counterexample on the unit three-torus."""

    if cells % 2:
        raise ValueError("the exact torus counterexample needs an even cell count")
    spacing = 1.0 / cells
    indices = np.indices((cells, cells, cells), dtype=np.int64)
    x, y, z = indices
    oscillation = np.stack(
        (
            (-1.0) ** (x + y),
            (-1.0) ** (y + z),
            (-1.0) ** (z + x),
        ),
        axis=-1,
    )
    field = np.empty((cells, cells, cells, 4), dtype=float)
    field[..., 0] = math.sqrt(1.0 - 3.0 * eta**2 * spacing**2)
    field[..., 1:] = eta * spacing * oscillation
    derivatives = np.stack(
        [
            (np.roll(field, -1, axis=direction) - field) / spacing
            for direction in range(3)
        ],
        axis=-1,
    )
    matrices = np.concatenate((field[..., :, None], derivatives), axis=-1)
    current = np.linalg.det(matrices)
    expected = -16.0 * eta**3 * field[0, 0, 0, 0]
    return {
        "cells": cells,
        "a": spacing,
        "maximum_distance_from_vacuum": float(
            np.max(np.linalg.norm(field - np.array([1.0, 0.0, 0.0, 0.0]), axis=-1))
        ),
        "minimum_n0": float(np.min(field[..., 0])),
        "hemisphere_degree_certificate": 0,
        "forward_current_integral": float(np.mean(current)),
        "exact_expected_forward_current": expected,
        "maximum_current_residual": float(np.max(np.abs(current - expected))),
    }


def torus_minor_counterexample(
    cells: int, eta: float = COUNTEREXAMPLE_ETA
) -> dict[str, Any]:
    """Exact three-periodic failure of weak continuity for a Skyrme 2-minor."""

    if cells % 3:
        raise ValueError("the exact minor counterexample needs a cell count divisible by three")
    spacing = 1.0 / cells
    x, y, _ = np.indices((cells, cells, cells), dtype=np.int64)
    phase = 2.0 * math.pi * (x + 2 * y) / 3.0
    field = np.zeros((cells, cells, cells, 4), dtype=float)
    field[..., 0] = math.sqrt(1.0 - eta**2 * spacing**2)
    field[..., 1] = eta * spacing * np.cos(phase)
    field[..., 2] = eta * spacing * np.sin(phase)
    derivative_x = (np.roll(field, -1, axis=0) - field) / spacing
    derivative_y = (np.roll(field, -1, axis=1) - field) / spacing
    minor_12 = (
        derivative_x[..., 1] * derivative_y[..., 2]
        - derivative_x[..., 2] * derivative_y[..., 1]
    )
    expected = 1.5 * math.sqrt(3.0) * eta**2
    return {
        "cells": cells,
        "a": spacing,
        "maximum_distance_from_vacuum": float(
            np.max(
                np.linalg.norm(
                    field - np.array([1.0, 0.0, 0.0, 0.0]), axis=-1
                )
            )
        ),
        "minor_12_mean": float(np.mean(minor_12)),
        "exact_expected_minor_12": expected,
        "maximum_minor_residual": float(np.max(np.abs(minor_12 - expected))),
        "continuum_limit_minor": 0.0,
        "hemisphere_degree_certificate": 0,
    }


def cutoff_counterexample(cells: int, eta: float = COUNTEREXAMPLE_ETA) -> dict[str, Any]:
    """Boundary-compatible localization of the exact periodic defect."""

    if cells % 2:
        raise ValueError("the localized counterexample needs an even cell count")
    spacing = 1.0 / cells
    axis = np.linspace(0.0, 1.0, cells + 1)
    x, y, z = np.indices((cells + 1, cells + 1, cells + 1), dtype=np.int64)
    mesh = np.meshgrid(axis, axis, axis, indexing="ij")
    cutoff = np.prod([np.sin(math.pi * coordinate) ** 2 for coordinate in mesh], axis=0)
    oscillation = np.stack(
        (
            (-1.0) ** (x + y),
            (-1.0) ** (y + z),
            (-1.0) ** (z + x),
        ),
        axis=-1,
    )
    tangent = eta * spacing * cutoff[..., None] * oscillation
    field = np.concatenate((np.ones((*cutoff.shape, 1)), tangent), axis=-1)
    field /= np.linalg.norm(field, axis=-1, keepdims=True)
    audit = forward_current(field, spacing)
    continuum_limit = -16.0 * eta**3 * (5.0 / 16.0) ** 3
    return {
        "cells": cells,
        "N": cells + 1,
        "a": spacing,
        "maximum_distance_from_vacuum": float(
            np.max(np.linalg.norm(field - np.array([1.0, 0.0, 0.0, 0.0]), axis=-1))
        ),
        "minimum_n0": float(np.min(field[..., 0])),
        "exact_vacuum_boundary": bool(
            np.all(field[0, ..., 0] == 1.0)
            and np.all(field[-1, ..., 0] == 1.0)
            and np.all(field[:, 0, ..., 0] == 1.0)
            and np.all(field[:, -1, ..., 0] == 1.0)
            and np.all(field[:, :, 0, 0] == 1.0)
            and np.all(field[:, :, -1, 0] == 1.0)
        ),
        "hemisphere_degree_certificate": 0,
        "predicted_diffuse_current_limit": continuum_limit,
        "current_error_from_predicted_limit": abs(
            audit["forward_current_integral"] - continuum_limit
        ),
        **audit,
    }


def compactness_counterexample(quick: bool) -> dict[str, Any]:
    torus_cells = (8, 16) if quick else (8, 16, 32, 64)
    minor_cells = (6, 12) if quick else (6, 12, 24, 48)
    cutoff_cells = (8, 16) if quick else (8, 12, 16, 24, 32, 48)
    torus_rows = [torus_counterexample(cells) for cells in torus_cells]
    minor_rows = [torus_minor_counterexample(cells) for cells in minor_cells]
    cutoff_rows = [cutoff_counterexample(cells) for cells in cutoff_cells]
    return {
        "sequence": (
            "u=( (-1)^(x+y), (-1)^(y+z), (-1)^(z+x) ); "
            "n_a=(sqrt(1-3 eta^2 a^2), eta a u) on the torus, localized by "
            "chi=prod_i sin^2(pi x_i) for the fixed-vacuum cube"
        ),
        "eta": COUNTEREXAMPLE_ETA,
        "exact_torus_identity": "J_forward=-16 eta^3 sqrt(1-3 eta^2 a^2)",
        "exact_minor_identity": (
            "for theta=2 pi (x+2y)/3 and (n1,n2)=eta a(cos theta,sin theta), "
            "M_12^(1,2)=(3 sqrt(3)/2) eta^2 although the continuum-limit minor is zero"
        ),
        "localized_limit": "J_forward converges weakly to -16 eta^3 chi^3 dx",
        "current_integrability_theorem": (
            "|J_a|^(4/3) <= |wedge^3 D^a n|^(4/3) "
            "<= (1/3)|wedge^2 D^a n|^2 pointwise"
        ),
        "verdict": (
            "singular concentration is excluded by the L^(4/3) bound, but "
            "compensated compactness fails through a nonzero diffuse oscillatory defect"
        ),
        "torus_rows": torus_rows,
        "minor_rows": minor_rows,
        "fixed_boundary_rows": cutoff_rows,
    }


def canonical_bond(
    left: Sequence[int], right: Sequence[int], period: int
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    left_tuple = tuple(int(value) for value in left)
    right_tuple = tuple(int(value) for value in right)
    displacement = tuple(right_tuple[i] - left_tuple[i] for i in range(3))
    first_nonzero = next(value for value in displacement if value != 0)
    if first_nonzero < 0:
        left_tuple, right_tuple = right_tuple, left_tuple
        displacement = tuple(-value for value in displacement)
    return tuple(value % period for value in left_tuple), displacement


def periodic_bonds(spec: Any, period: int = 2) -> list[tuple[Any, Any]]:
    """Unique infinite-graph bonds per periodic cell."""

    bonds: set[tuple[tuple[int, int, int], tuple[int, int, int]]] = set()
    for origin in itertools.product(range(period), repeat=3):
        if isinstance(spec, TRI.TriangulationSpec):
            tetrahedra, _ = TRI.cube_tetrahedra(origin, TRI.cube_signs(origin, spec))
        elif isinstance(spec, TRI.FiveTetrahedronSpec):
            central_parity = (spec.phase + sum(origin)) % 2
            tetrahedra, _ = TRI.cube_five_tetrahedra(origin, central_parity)
        else:
            raise TypeError("unsupported triangulation")
        for tetrahedron in tetrahedra:
            for left, right in itertools.combinations(tetrahedron, 2):
                bonds.add(canonical_bond(left, right, period))
    return sorted(bonds)


def direction_cycle_balance(spec: Any, period: int = 2) -> dict[str, Any]:
    nodes = list(itertools.product(range(period), repeat=3))
    rows = []
    bonds = periodic_bonds(spec, period)
    for direction in sorted(set(displacement for _, displacement in bonds)):
        incoming = {node: 0 for node in nodes}
        outgoing = {node: 0 for node in nodes}
        selected = [row for row in bonds if row[1] == direction]
        for left, displacement in selected:
            right = tuple((left[i] + displacement[i]) % period for i in range(3))
            outgoing[left] += 1
            incoming[right] += 1
        rows.append(
            {
                "direction": list(direction),
                "bonds_per_period_cell": len(selected),
                "nodewise_incoming_equals_outgoing": all(
                    incoming[node] == outgoing[node] for node in nodes
                ),
                "maximum_node_imbalance": max(
                    abs(incoming[node] - outgoing[node]) for node in nodes
                ),
            }
        )
    return {
        "mesh": spec.name,
        "period": period,
        "bonds": len(bonds),
        "directions": rows,
        "all_directions_cycle_balanced": all(
            row["nodewise_incoming_equals_outgoing"] for row in rows
        ),
    }


def raw_six_quartic(matrix: np.ndarray) -> float:
    columns = [matrix[:, index] for index in range(3)]
    values = [*columns]
    values.extend(columns[i] + columns[j] for i, j in ((0, 1), (0, 2), (1, 2)))
    values.append(columns[0] + columns[1] + columns[2])
    return float(sum((value @ value) ** 2 for value in values))


def raw_five_quartic(matrix: np.ndarray) -> float:
    columns = [matrix[:, index] for index in range(3)]
    value = sum((column @ column) ** 2 for column in columns)
    for first, second in ((0, 1), (0, 2), (1, 2)):
        plus = columns[first] + columns[second]
        minus = columns[first] - columns[second]
        value += 0.5 * ((plus @ plus) ** 2 + (minus @ minus) ** 2)
    return float(value)


def cell_energy_gradient(
    variables: np.ndarray,
    matrix: np.ndarray,
    bonds: Sequence[tuple[Any, Any]],
    period: int,
) -> tuple[float, np.ndarray]:
    nodes = list(itertools.product(range(period), repeat=3))
    node_index = {node: index for index, node in enumerate(nodes)}
    target_dimension = matrix.shape[0]
    corrector = np.zeros((len(nodes), target_dimension), dtype=float)
    corrector[1:] = variables.reshape(-1, target_dimension)
    gradient = np.zeros_like(corrector)
    energy = 0.0
    for left, displacement in bonds:
        right = tuple((left[i] + displacement[i]) % period for i in range(3))
        bond = (
            matrix @ np.asarray(displacement, dtype=float)
            + corrector[node_index[right]]
            - corrector[node_index[left]]
        )
        norm_square = float(bond @ bond)
        energy += norm_square**2
        derivative = 4.0 * norm_square * bond
        gradient[node_index[right]] += derivative
        gradient[node_index[left]] -= derivative
    normalization = float(period**3)
    return energy / normalization, gradient[1:].ravel() / normalization


def solve_cell_problem(
    spec: Any, matrix: np.ndarray, period: int, seed: int
) -> dict[str, Any]:
    bonds = periodic_bonds(spec, period)
    variables = 0.15 * np.random.default_rng(seed).normal(
        size=(period**3 - 1) * matrix.shape[0]
    )
    result = minimize(
        lambda value: cell_energy_gradient(value, matrix, bonds, period),
        variables,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 4000, "gtol": 1.0e-12, "ftol": 1.0e-15, "maxls": 80},
    )
    zero_variables = np.zeros_like(variables)
    cauchy_born = cell_energy_gradient(zero_variables, matrix, bonds, period)[0]
    formula = raw_six_quartic(matrix) if isinstance(spec, TRI.TriangulationSpec) else raw_five_quartic(matrix)
    return {
        "mesh": spec.name,
        "period": period,
        "matrix": matrix.tolist(),
        "bonds": len(bonds),
        "optimizer_success": bool(result.success),
        "optimizer_iterations": int(result.nit),
        "cell_minimum": float(result.fun),
        "zero_corrector_energy": cauchy_born,
        "closed_formula_energy": formula,
        "relative_cell_vs_formula_residual": abs(float(result.fun) - formula)
        / max(abs(formula), 1.0e-15),
        "relative_zero_vs_formula_residual": abs(cauchy_born - formula)
        / max(abs(formula), 1.0e-15),
        "corrector_maximum_absolute_value": float(np.max(np.abs(result.x))),
        "gradient_infinity_norm": float(np.max(np.abs(result.jac))),
    }


def homogenized_cell_audit(quick: bool) -> dict[str, Any]:
    matrices = {
        "rank_one_x": np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "equal_x_y": np.array([[1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "isometry": np.eye(3),
        "generic": np.array([[1.0, 0.3, -0.2], [0.4, -0.7, 0.1], [0.0, 0.2, 0.8]]),
    }
    if not quick:
        rng = np.random.default_rng(20260718)
        for index in range(4):
            matrices[f"random_{index}"] = rng.normal(size=(3, 3))
    specs = (TRI.UNIFORM_SPECS[0], *TRI.FIVE_TET_SPECS)
    balances = [direction_cycle_balance(spec, 2) for spec in specs]
    rows = []
    seed = 7300
    for spec in specs:
        for matrix in matrices.values():
            seed += 1
            rows.append(solve_cell_problem(spec, matrix, 2, seed))
    rank_one = matrices["rank_one_x"]
    equal = matrices["equal_x_y"]
    return {
        "cell_formula": (
            "Q_hom(A)=|Y|^(-1) inf_phi sum_(p,r in B_Y) "
            "|A r + phi(p+r)-phi(p)|^4"
        ),
        "proof_certificate": (
            "for every bond direction the periodic directed graph is cycle-balanced; "
            "the mean corrector increment is zero and convex Jensen gives phi=0 exactly"
        ),
        "six_tet_closed_formula": (
            "sum_i |A e_i|^4 + sum_(i<j)|A(e_i+e_j)|^4 "
            "+ |A(e1+e2+e3)|^4"
        ),
        "five_tet_closed_formula": (
            "sum_i |A e_i|^4 + 1/2 sum_(i<j)(|A(e_i+e_j)|^4+|A(e_i-e_j)|^4)"
        ),
        "direction_balances": balances,
        "optimization_rows": rows,
        "nonproportionality": {
            "six_over_five_rank_one_x": raw_six_quartic(rank_one) / raw_five_quartic(rank_one),
            "six_over_five_equal_x_y": raw_six_quartic(equal) / raw_five_quartic(equal),
            "homogenized_densities_scalar_multiples": False,
        },
        "verdict": (
            "the finite-R barrier density is computed, but triangulation dependence "
            "survives homogenization and therefore blocks regulator universality"
        ),
    }


@dataclass(frozen=True)
class QuotientPlan:
    case_id: str
    extent: int
    length: float
    triangulation: str = "uniform_ppp"
    translation_lattice_units: tuple[float, float, float] = (0.0, 0.0, 0.0)
    roles: tuple[str, ...] = ()


def add_plan(
    plans: dict[tuple[Any, ...], QuotientPlan],
    case_id: str,
    extent: int,
    length: float,
    triangulation: str,
    translation: Sequence[float],
    role: str,
) -> None:
    translation_tuple = tuple(float(value) for value in translation)
    key = (extent, float(length), triangulation, translation_tuple)
    if key in plans:
        previous = plans[key]
        plans[key] = QuotientPlan(
            previous.case_id,
            previous.extent,
            previous.length,
            previous.triangulation,
            previous.translation_lattice_units,
            (*previous.roles, role),
        )
    else:
        plans[key] = QuotientPlan(
            case_id, extent, length, triangulation, translation_tuple, (role,)
        )


def quotient_case_plan(quick: bool) -> list[QuotientPlan]:
    plans: dict[tuple[Any, ...], QuotientPlan] = {}
    joint = ((17, 6.0), (21, 6.5)) if quick else (
        (25, 6.0),
        (29, 6.5),
        (33, 7.0),
        (37, 7.5),
    )
    for extent, length in joint:
        add_plan(
            plans,
            f"joint_N{extent}_L{length:g}",
            extent,
            length,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            "joint_a_down_L_up",
        )
    comparison_extent, comparison_length = joint[1]
    offsets = (
        (0.37, -0.23, 0.41),
        (0.5, 0.5, 0.5),
    ) if quick else (
        (0.37, -0.23, 0.41),
        (0.5, 0.5, 0.5),
        (-0.41, 0.19, 0.33),
    )
    for index, offset in enumerate(offsets):
        add_plan(
            plans,
            f"translation_{index}_N{comparison_extent}",
            comparison_extent,
            comparison_length,
            "uniform_ppp",
            offset,
            "translation_quotient",
        )
    add_plan(
        plans,
        f"translation_center_N{comparison_extent}",
        comparison_extent,
        comparison_length,
        "uniform_ppp",
        (0.0, 0.0, 0.0),
        "translation_quotient",
    )
    meshes = ("five_tet_phase0",) if quick else ("five_tet_phase0", "five_tet_phase1")
    for mesh in meshes:
        for offset in ((0.0, 0.0, 0.0), (0.37, -0.23, 0.41)):
            add_plan(
                plans,
                f"mesh_{mesh}_{'center' if offset == (0.0, 0.0, 0.0) else 'translated'}",
                comparison_extent,
                comparison_length,
                mesh,
                offset,
                "triangulation_quotient",
            )
    add_plan(
        plans,
        f"mesh_uniform_center_N{comparison_extent}",
        comparison_extent,
        comparison_length,
        "uniform_ppp",
        (0.0, 0.0, 0.0),
        "triangulation_quotient",
    )
    add_plan(
        plans,
        f"mesh_uniform_translated_N{comparison_extent}",
        comparison_extent,
        comparison_length,
        "uniform_ppp",
        (0.37, -0.23, 0.41),
        "triangulation_quotient",
    )
    return list(plans.values())


def quotient_observables(field: np.ndarray, axis: np.ndarray) -> dict[str, Any]:
    mesh = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    density = np.maximum(1.0 - field[..., 0], 0.0)
    total = float(np.sum(density))
    if total <= 0.0:
        barycentre = np.zeros(3)
        covariance = np.zeros((3, 3))
    else:
        barycentre = np.sum(density[..., None] * mesh, axis=(0, 1, 2)) / total
        centered = mesh - barycentre
        covariance = np.sum(
            density[..., None, None]
            * centered[..., :, None]
            * centered[..., None, :],
            axis=(0, 1, 2),
        ) / total
    eigenvalues = np.linalg.eigvalsh(covariance)
    radius = np.linalg.norm(mesh - barycentre, axis=-1)
    edges = np.linspace(0.0, PROFILE_RADIUS_MAXIMUM, PROFILE_BINS + 1)
    profile, _ = np.histogram(radius, bins=edges, weights=density)
    profile = profile.astype(float)
    if np.sum(profile) > 0.0:
        profile /= np.sum(profile)
    return {
        "density": "1-n0",
        "barycentre": barycentre.tolist(),
        "barycentre_norm": float(np.linalg.norm(barycentre)),
        "centered_covariance_eigenvalues": eigenvalues.tolist(),
        "centered_rms": float(math.sqrt(max(float(np.trace(covariance)), 0.0))),
        "anisotropy_ratio": float(
            (eigenvalues[-1] - eigenvalues[0]) / max(float(np.sum(eigenvalues)), 1.0e-15)
        ),
        "radial_profile_edges": edges.tolist(),
        "normalized_centered_radial_profile": profile.tolist(),
    }


def shifted_to_barycentre(
    field: np.ndarray,
    breathing: np.ndarray,
    barycentre: np.ndarray,
    spacing: float,
    fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    voxel_shift = -fraction * barycentre / spacing
    shifted_field = np.empty_like(field)
    vacuum = (1.0, 0.0, 0.0, 0.0)
    for component in range(4):
        shifted_field[..., component] = ndimage.shift(
            field[..., component],
            shift=voxel_shift,
            order=1,
            mode="constant",
            cval=vacuum[component],
            prefilter=False,
        )
    norm = np.linalg.norm(shifted_field, axis=-1, keepdims=True)
    shifted_field /= np.maximum(norm, 1.0e-14)
    shifted_breathing = ndimage.shift(
        breathing,
        shift=voxel_shift,
        order=1,
        mode="constant",
        cval=1.0,
        prefilter=False,
    )
    return E7.exact_vacuum_boundary(shifted_field, shifted_breathing)


def local_relax(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    quick: bool,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    history = []
    maximum_charts = 2 if quick else 4
    maximum_iterations = 350 if quick else 650
    for chart in range(maximum_charts):
        base, scalar, frames = E7.fixed_chart(field, breathing)

        def objective(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
            return E8.chart_objective(
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
                EPSILON,
            )

        result = minimize(
            objective,
            np.zeros(4 * len(base)),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": maximum_iterations,
                "maxcor": 30,
                "maxls": 80,
                "ftol": 1.0e-15,
                "gtol": 1.0e-10,
            },
        )
        field, breathing, _, _ = E7.chart_fields(
            result.x, field, breathing, base, scalar, frames
        )
        stationarity = E9.stationarity_and_concentration(
            field, breathing, spacing, model, left, right, gamma, EPSILON
        )
        history.append(
            {
                "chart": chart,
                "optimizer_success": bool(result.success),
                "optimizer_message": str(result.message),
                "iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "energy": stationarity["total_energy"],
                "gradient_density": stationarity[
                    "direct_projected_gradient_l2_density"
                ],
                "minimum_pair_margin": stationarity["minimum_pair_margin"],
            }
        )
        if stationarity["direct_projected_gradient_l2_density"] < STATIONARITY_TOLERANCE:
            break
    return field, breathing, history


def relax_quotient_case(
    plan: QuotientPlan, solution: Any, model: Any, quick: bool
) -> dict[str, Any]:
    spec = TRIANGULATIONS[plan.triangulation]
    field, breathing, axis, _ = TRI.translated_radial_background(
        solution,
        plan.extent,
        plan.length,
        plan.translation_lattice_units,
        "lattice",
        radial_epsilon=E6.EPS,
    )
    field, breathing = E7.exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    gamma = E9.gamma_trajectory(spacing, 1.0)
    tetrahedra, orientations = TRI.tetrahedral_complex(plan.extent, spec)
    left, right = TRI.unique_edges(tetrahedra)
    initial_observables = quotient_observables(field, axis)
    all_relaxation_history: list[dict[str, Any]] = []
    centering_history: list[dict[str, Any]] = []
    maximum_centerings = 1 if quick else 2
    for centering_step in range(maximum_centerings + 1):
        field, breathing, local_history = local_relax(
            field, breathing, spacing, model, left, right, gamma, quick
        )
        for row in local_history:
            all_relaxation_history.append({"centering_step": centering_step, **row})
        observables = quotient_observables(field, axis)
        barycentre = np.asarray(observables["barycentre"], dtype=float)
        if (
            centering_step == maximum_centerings
            or np.linalg.norm(barycentre) <= 0.08 * spacing
        ):
            break
        accepted = False
        for fraction in (1.0, 0.5, 0.25):
            candidate_field, candidate_breathing = shifted_to_barycentre(
                field, breathing, barycentre, spacing, fraction
            )
            flat = candidate_field.reshape(-1, 4)
            dots = np.einsum("ij,ij->i", flat[left], flat[right])
            candidate_topology = E9.topology_summary(
                candidate_field, tetrahedra, orientations
            )
            if float(np.min(dots)) > EPSILON and candidate_topology["all_targets_B1"]:
                field, breathing = candidate_field, candidate_breathing
                centering_history.append(
                    {
                        "step": centering_step,
                        "fraction": fraction,
                        "barycentre_before": barycentre.tolist(),
                        "physical_shift": (-fraction * barycentre).tolist(),
                        "minimum_pair_margin_after_shift": float(np.min(dots) - EPSILON),
                        "degree_after_shift": candidate_topology["baryon_numbers"],
                    }
                )
                accepted = True
                break
        if not accepted:
            centering_history.append(
                {
                    "step": centering_step,
                    "accepted": False,
                    "reason": "no admissible degree-preserving interpolation fraction",
                }
            )
            break
    final = E9.stationarity_and_concentration(
        field, breathing, spacing, model, left, right, gamma, EPSILON
    )
    topology = E9.topology_summary(field, tetrahedra, orientations)
    observables = quotient_observables(field, axis)
    triangulation_audit = TRI.triangulation_audit(plan.extent, spec)
    passed = bool(
        final["direct_projected_gradient_l2_density"] < STATIONARITY_TOLERANCE
        and final["minimum_pair_margin"] > PAIR_MARGIN_TOLERANCE
        and topology["all_targets_B1"]
        and triangulation_audit["valid_conforming_triangulation"]
    )
    return {
        "case_id": plan.case_id,
        "roles": list(plan.roles),
        "N": plan.extent,
        "L": plan.length,
        **E9.scaling_numbers(spacing, 1.0),
        "triangulation": asdict(spec),
        "translation_lattice_units": list(plan.translation_lattice_units),
        "translation_physical_units": (
            spacing * np.asarray(plan.translation_lattice_units)
        ).tolist(),
        "translation_quotient": "barycentre-aligned 1-n0 density",
        "core_or_centre_anchor_used": False,
        "dynamic_centering_projection_used": bool(centering_history),
        "initial_quotient_observables": initial_observables,
        "centering_history": centering_history,
        "relaxation_history": all_relaxation_history,
        "final": final,
        "final_topology": topology,
        "quotient_observables": observables,
        "case_pass": passed,
    }


def relative_spread(values: Sequence[float]) -> float:
    values = [float(value) for value in values]
    scale = max(max(abs(value) for value in values), 1.0e-15)
    return (max(values) - min(values)) / scale


def profile_l1(first: dict[str, Any], second: dict[str, Any]) -> float:
    a = np.asarray(first["normalized_centered_radial_profile"], dtype=float)
    b = np.asarray(second["normalized_centered_radial_profile"], dtype=float)
    return 0.5 * float(np.sum(np.abs(a - b)))


def quotient_gate(cases: list[dict[str, Any]], quick: bool) -> dict[str, Any]:
    joint = sorted(
        [row for row in cases if "joint_a_down_L_up" in row["roles"]],
        key=lambda row: row["a"],
        reverse=True,
    )
    translations = [row for row in cases if "translation_quotient" in row["roles"]]
    triangulations = [
        row for row in cases if "triangulation_quotient" in row["roles"]
    ]
    joint_tail = joint[-3:] if len(joint) >= 3 else joint
    joint_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in joint_tail]
    )
    joint_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in joint_tail]
    )
    joint_profile_changes = [
        profile_l1(joint[index]["quotient_observables"], joint[index + 1]["quotient_observables"])
        for index in range(len(joint) - 1)
    ]
    translation_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in translations]
    )
    translation_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in translations]
    )
    translation_profiles = [
        profile_l1(
            translations[0]["quotient_observables"], row["quotient_observables"]
        )
        for row in translations[1:]
    ]
    triangulation_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in triangulations]
    )
    triangulation_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in triangulations]
    )
    maximum_center_residual_over_a = max(
        row["quotient_observables"]["barycentre_norm"] / row["a"] for row in cases
    )
    production_coverage = bool(
        not quick
        and len(joint) >= 4
        and min(row["a"] for row in joint) <= 0.21
        and max(row["L"] for row in joint) >= 7.5
        and len(translations) >= 4
        and len(triangulations) >= 6
    )
    translation_quotient_pass = bool(
        translation_energy_spread < 0.01
        and translation_radius_spread < 0.02
        and max(translation_profiles, default=0.0) < 0.03
    )
    dynamic_centering_pass = maximum_center_residual_over_a < 0.15
    numerical_pass = bool(
        all(row["case_pass"] for row in cases)
        and production_coverage
        and joint_energy_spread < 0.03
        and joint_radius_spread < 0.05
        and max(joint_profile_changes, default=0.0) < 0.05
        and translation_quotient_pass
        and triangulation_energy_spread < 0.03
        and triangulation_radius_spread < 0.05
    )
    return {
        "all_cases_pass": all(row["case_pass"] for row in cases),
        "production_coverage": production_coverage,
        "joint_tail_energy_relative_spread": joint_energy_spread,
        "joint_tail_energy_threshold": 0.03,
        "joint_tail_centered_rms_relative_spread": joint_radius_spread,
        "joint_tail_centered_rms_threshold": 0.05,
        "joint_maximum_successive_profile_L1": max(joint_profile_changes, default=0.0),
        "joint_profile_L1_threshold": 0.05,
        "translation_energy_relative_spread": translation_energy_spread,
        "translation_energy_threshold": 0.01,
        "translation_centered_rms_relative_spread": translation_radius_spread,
        "translation_centered_rms_threshold": 0.02,
        "translation_maximum_profile_L1": max(translation_profiles, default=0.0),
        "translation_profile_L1_threshold": 0.03,
        "translation_quotient_thresholds_pass": translation_quotient_pass,
        "maximum_barycentre_residual_over_a": maximum_center_residual_over_a,
        "barycentre_residual_over_a_threshold": 0.15,
        "dynamic_centering_threshold_pass": dynamic_centering_pass,
        "centering_interpretation": (
            "translation quotient is the primary gate; barycentre residual is a "
            "Peierls-locking diagnostic because re-relaxation may return a shifted "
            "field to a mesh-preferred representative"
        ),
        "triangulation_energy_relative_spread": triangulation_energy_spread,
        "triangulation_energy_threshold": 0.03,
        "triangulation_centered_rms_relative_spread": triangulation_radius_spread,
        "triangulation_centered_rms_threshold": 0.05,
        "numerical_quotient_gate_pass": numerical_pass,
    }


def markdown(result: dict[str, Any]) -> str:
    compactness = result["compactness_counterexample"]
    cell = result["finite_R_cell_problem"]
    gate = result["quotient_scan"]["gate"]
    lines = [
        "# AP-E10 compactness, cell formula, and translation quotient card",
        "",
        f"Status: `{result['status']}`",
        "",
        "## Exact compensated-compactness verdict",
        "",
        compactness["verdict"] + ".",
        "",
        "| cells | a | sup distance to vacuum | forward current | predicted limit | L4/3 ratio |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in compactness["fixed_boundary_rows"]:
        lines.append(
            "| {cells} | {a:.6f} | {distance:.3e} | {current:.8f} | {limit:.8f} | {ratio:.6f} |".format(
                cells=row["cells"],
                a=row["a"],
                distance=row["maximum_distance_from_vacuum"],
                current=row["forward_current_integral"],
                limit=row["predicted_diffuse_current_limit"],
                ratio=row["L43_bound_ratio"],
            )
        )
    lines.extend(
        [
            "",
            "## Finite-R periodic cell problem",
            "",
            cell["verdict"] + ".",
            "",
            "The exact six/five ratios are `{:.12g}` on rank-one x and `{:.12g}` on equal x-y.".format(
                cell["nonproportionality"]["six_over_five_rank_one_x"],
                cell["nonproportionality"]["six_over_five_equal_x_y"],
            ),
            "",
            "## Dynamic-centering and quotient cases",
            "",
            "| id | N | L | a | mesh | offset/a | E | grad | B | center/a | rms | pass |",
            "|:--|---:|---:|---:|:--|:--|---:|---:|:--|---:|---:|:--:|",
        ]
    )
    for row in result["quotient_scan"]["cases"]:
        lines.append(
            "| {id} | {N} | {L:.3f} | {a:.6f} | {mesh} | {offset} | {E:.8f} | {grad:.2e} | {B} | {center:.3e} | {rms:.5f} | {passed} |".format(
                id=row["case_id"],
                N=row["N"],
                L=row["L"],
                a=row["a"],
                mesh=row["triangulation"]["name"],
                offset=row["translation_lattice_units"],
                E=row["final"]["total_energy"],
                grad=row["final"]["direct_projected_gradient_l2_density"],
                B=row["final_topology"]["baryon_numbers"],
                center=row["quotient_observables"]["barycentre_norm"] / row["a"],
                rms=row["quotient_observables"]["centered_rms"],
                passed=row["case_pass"],
            )
        )
    lines.extend(
        [
            "",
            "## Fail-closed gate",
            "",
            f"Numerical quotient gate: `{gate['numerical_quotient_gate_pass']}`.",
            f"Discrete compensated compactness: `{result['fail_closed']['forward_skyrme_compensated_compactness']}`.",
            f"Finite-R cell density computed: `{result['fail_closed']['finite_R_periodic_cell_density_computed']}`.",
            f"Triangulation-independent finite-R density: `{result['fail_closed']['finite_R_density_triangulation_independent']}`.",
            f"Regulator-stable continuum background: `{result['fail_closed']['regulator_stable_continuum_background']}`.",
            f"Same-action Hessian allowed: `{result['hessian_gate_open']}`.",
            f"Regulated determinant variation allowed: `{result['determinant_variation_gate_open']}`.",
            "",
            f"Checks: **{result['checks_passed']}/{result['checks_total']}**.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT / "ap_e10_compactness_homogenization_centering"
    stem.with_suffix(".json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    stem.with_suffix(".md").write_text(markdown(result), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    arguments = parser.parse_args()
    CHECKS.clear()

    counterexample = compactness_counterexample(arguments.quick)
    torus_residual = max(
        row["maximum_current_residual"] for row in counterexample["torus_rows"]
    )
    minor_residual = max(
        row["maximum_minor_residual"] for row in counterexample["minor_rows"]
    )
    cutoff_rows = counterexample["fixed_boundary_rows"]
    check(
        "P1_compactness",
        "the exact three-periodic sphere-valued sequence disproves weak continuity of the one-corner forward Skyrme minors",
        minor_residual < 1.0e-12
        and max(
            row["maximum_distance_from_vacuum"]
            for row in counterexample["minor_rows"]
        )
        < 0.1
        and all(
            row["hemisphere_degree_certificate"] == 0
            for row in counterexample["minor_rows"]
        ),
        f"maximum exact-minor residual={minor_residual:.3e}",
    )
    check(
        "P1_compactness",
        "the exact two-periodic sequence gives a nonzero diffuse forward-current defect at degree zero",
        torus_residual < 1.0e-12
        and all(
            row["hemisphere_degree_certificate"] == 0
            for row in counterexample["torus_rows"]
        ),
        f"maximum exact-current residual={torus_residual:.3e}",
    )
    check(
        "P1_current",
        "the fixed-boundary localization remains energy-bounded and approaches a nonzero diffuse current defect",
        all(row["exact_vacuum_boundary"] for row in cutoff_rows)
        and max(row["dirichlet_energy_unit_coefficient"] for row in cutoff_rows) < 2.0
        and max(row["skyrme_energy_unit_coefficient"] for row in cutoff_rows) < 2.0
        and cutoff_rows[-1]["current_error_from_predicted_limit"]
        < cutoff_rows[0]["current_error_from_predicted_limit"],
        "current errors="
        + str([row["current_error_from_predicted_limit"] for row in cutoff_rows]),
    )
    check(
        "P1_current",
        "the Skyrme-minor bound enforces the pointwise L4/3 current estimate and excludes singular concentration",
        max(row["L43_bound_ratio"] for row in cutoff_rows) <= 1.0 + 1.0e-12,
        "maximum normalized inequality ratio="
        + str(max(row["L43_bound_ratio"] for row in cutoff_rows)),
    )

    cell = homogenized_cell_audit(arguments.quick)
    check(
        "P2_cell",
        "all six- and five-tet bond directions are cycle-balanced, proving zero corrector by Jensen",
        all(
            row["all_directions_cycle_balanced"]
            for row in cell["direction_balances"]
        ),
        "balances="
        + str(
            [
                (row["mesh"], row["all_directions_cycle_balanced"])
                for row in cell["direction_balances"]
            ]
        ),
    )
    check(
        "P2_cell",
        "random-start corrector solves reproduce the exact closed cell formulas",
        max(
            row["relative_cell_vs_formula_residual"]
            for row in cell["optimization_rows"]
        )
        < 1.0e-10,
        "maximum relative residual="
        + str(
            max(
                row["relative_cell_vs_formula_residual"]
                for row in cell["optimization_rows"]
            )
        ),
    )
    check(
        "P2_cell",
        "six- and five-tet homogenized quartics are nonproportional",
        abs(cell["nonproportionality"]["six_over_five_rank_one_x"] - 4.0 / 3.0)
        < 1.0e-14
        and abs(cell["nonproportionality"]["six_over_five_equal_x_y"] - 3.0)
        < 1.0e-14,
        str(cell["nonproportionality"]),
    )

    model = E6.Model()
    profile = E6.solve_relaxed_profile(box=16.0, model=model, sample_points=4097)
    checksum = profile["profile_sha256_r_F_s_little_endian_float64"]
    solution = profile["solution"]
    check(
        "P0_provenance",
        "the canonical AP-E6 seed checksum remains frozen",
        checksum == "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59",
        f"checksum={checksum}",
    )
    plans = quotient_case_plan(arguments.quick)
    cases = []
    for index, plan in enumerate(plans, start=1):
        print(
            f"AP-E10 case {index}/{len(plans)}: {plan.case_id} "
            f"N={plan.extent} L={plan.length:g} mesh={plan.triangulation} "
            f"offset/a={plan.translation_lattice_units}",
            flush=True,
        )
        cases.append(relax_quotient_case(plan, solution, model, arguments.quick))
    quotient = quotient_gate(cases, arguments.quick)
    check(
        "P3_quotient",
        "all quotient representatives are unanchored stationary admissible B=1 fields",
        all(row["case_pass"] for row in cases),
        "failed=" + str([row["case_id"] for row in cases if not row["case_pass"]]),
    )
    check(
        "P3_quotient",
        "the smaller-a/larger-L and translation/triangulation thresholds are evaluated without promoting failed mathematics",
        isinstance(quotient["numerical_quotient_gate_pass"], bool),
        f"numerical quotient gate={quotient['numerical_quotient_gate_pass']}; production={quotient['production_coverage']}",
    )

    fail_closed = {
        "forward_skyrme_compensated_compactness": False,
        "singular_topological_current_concentration_excluded_for_forward_current": True,
        "diffuse_oscillatory_topological_current_defect_present": True,
        "geometric_degree_current_controlled_by_forward_skyrme_stencil": False,
        "finite_R_periodic_cell_density_computed": True,
        "finite_R_zero_corrector_proven": True,
        "finite_R_density_triangulation_independent": False,
        "translation_quotient_implemented": True,
        "dynamic_centering_implemented": True,
        "full_action_gamma_liminf_proven": False,
        "regulator_stable_continuum_background": False,
        "same_action_riemann_hessian_allowed": False,
        "regulated_determinant_variation_allowed": False,
        "gw_domain_wall_parallel_lane_retained": True,
        "so3_mod_two_index_parallel_lane_retained": True,
    }
    result = {
        "artifact": "AP-E10 compactness, finite-R homogenization, and translation quotient",
        "status": "quick_fail_closed" if arguments.quick else "production_exact_stencil_no_go_fail_closed",
        "generated_utc": "deterministic-no-wall-clock",
        "quick": arguments.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "canonical_profile_checksum": checksum,
        "source_manifest": [
            source_row(path)
            for path in (
                THIS_SCRIPT,
                E9_SCRIPT,
                E9.TRI_SCRIPT,
                E9.E8_SCRIPT,
                E8.E7_SCRIPT,
                E8.E6_SCRIPT,
            )
        ],
        "compactness_counterexample": counterexample,
        "finite_R_cell_problem": cell,
        "quotient_scan": {
            "plan": [asdict(plan) for plan in plans],
            "cases": cases,
            "gate": quotient,
        },
        "fail_closed": fail_closed,
        "regulator_stable": False,
        "hessian_gate_open": False,
        "determinant_variation_gate_open": False,
        "portal_start_allowed": False,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    if result["checks_passed"] != result["checks_total"]:
        result["status"] = "mechanical_failure_all_physics_gates_fail_closed"
    write_outputs(result)
    print(
        f"AP-E10 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"quotient_gate={quotient['numerical_quotient_gate_pass']}; "
        "regulator_stable=False; hessian_gate=False",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
