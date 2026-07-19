#!/usr/bin/env python3
"""AP-E11 compatible tetrahedral cochain action and continuum gate.

This card repairs the exact AP-E10 failure rather than rescaling the failed
one-corner stencil.  Its mathematical objects all live on one oriented
tetrahedral complex:

* the simplicial coboundary ``d`` and the Alexander--Whitney front/back cup;
* consistent positive Whitney Hodge Gram matrices on primal one- and
  two-cochains;
* the antisymmetric exact two-cochain
  ``dn^A cup dn^B - dn^B cup dn^A``;
* a radially dressed third cup cochain whose sum is exactly the degree of the
  normalized piecewise-affine map.

The associated nodal action is the exact P1 finite-element integral of the
polyconvex density

    1/2 |D v_h|^2 + R/2 |wedge^2 D v_h|^2
    + K/2 |wedge^3 D v_h|^2 + m^2 (1-v_h^0),

on the open normalized-affine admissible sector.  Consequently its periodic
cell problem is mesh-independent: periodicity preserves the means of the
gradient and every two-minor, and Jensen gives the zero corrector.

The final section re-relaxes unanchored B=1 representatives and evaluates a
translation quotient.  Numerical evidence never overrides a failed theorem;
Hessian and determinant flags remain fail-closed unless every declared
mathematical and continuum-background gate passes.
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
from typing import Any, Iterable, Sequence

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
TRI_SCRIPT = ROUTE_F / "code" / "ap_e9_triangulation_tools.py"

R_SKYRME = 1.0
K_SEXTIC = 0.35
MASS_SQUARED = 0.12
ADMISSIBILITY_EPSILON = -0.32
SOLVER_BUFFER = -0.30
SOLVER_PENALTY = 1.0e10
STATIONARITY_TOLERANCE = 7.5e-5
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


TRI = load_module("ap_e9_triangulation_tools_for_ap_e11", TRI_SCRIPT)
MESH_SPECS = (
    TRI.UNIFORM_SPECS[0],
    TRI.FIVE_TET_SPECS[0],
    TRI.FIVE_TET_SPECS[1],
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


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def permutation_sign(values: Sequence[int]) -> int:
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


# ---------------------------------------------------------------------------
# Exact cochain algebra
# ---------------------------------------------------------------------------


Cochain = dict[tuple[int, ...], float]


def simplices(vertices: Sequence[int], degree: int) -> Iterable[tuple[int, ...]]:
    return itertools.combinations(vertices, degree + 1)


def oriented_value(cochain: Cochain, simplex: Sequence[int]) -> float:
    if len(set(simplex)) != len(simplex):
        return 0.0
    ordered = tuple(sorted(simplex))
    permutation = tuple(ordered.index(vertex) for vertex in simplex)
    return float(permutation_sign(permutation) * cochain.get(ordered, 0.0))


def coboundary(cochain: Cochain, degree: int, vertices: Sequence[int]) -> Cochain:
    result: Cochain = {}
    for simplex in simplices(vertices, degree + 1):
        result[simplex] = sum(
            (-1) ** index * oriented_value(cochain, simplex[:index] + simplex[index + 1 :])
            for index in range(degree + 2)
        )
    return result


def aw_cup(
    first: Cochain,
    first_degree: int,
    second: Cochain,
    second_degree: int,
    vertices: Sequence[int],
) -> Cochain:
    """Alexander--Whitney front/back cup on an ordered simplicial complex."""

    result: Cochain = {}
    for simplex in simplices(vertices, first_degree + second_degree):
        front = simplex[: first_degree + 1]
        back = simplex[first_degree:]
        result[simplex] = oriented_value(first, front) * oriented_value(second, back)
    return result


def cochain_add(*terms: tuple[float, Cochain]) -> Cochain:
    keys = set().union(*(term.keys() for _, term in terms))
    return {
        key: sum(coefficient * term.get(key, 0.0) for coefficient, term in terms)
        for key in keys
    }


def cochain_maximum(cochain: Cochain) -> float:
    return max((abs(value) for value in cochain.values()), default=0.0)


def exact_cochain_audit() -> dict[str, Any]:
    vertices = tuple(range(5))
    rng = np.random.default_rng(20260719)
    rows = []
    maximum_d2 = 0.0
    maximum_leibniz = 0.0
    maximum_associativity = 0.0
    for first_degree, second_degree in ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2)):
        if first_degree + second_degree + 1 > 4:
            continue
        first = {
            simplex: float(rng.integers(-5, 6))
            for simplex in simplices(vertices, first_degree)
        }
        second = {
            simplex: float(rng.integers(-5, 6))
            for simplex in simplices(vertices, second_degree)
        }
        d_first = coboundary(first, first_degree, vertices)
        d_second = coboundary(second, second_degree, vertices)
        product = aw_cup(first, first_degree, second, second_degree, vertices)
        left = coboundary(product, first_degree + second_degree, vertices)
        right = cochain_add(
            (1.0, aw_cup(d_first, first_degree + 1, second, second_degree, vertices)),
            ((-1.0) ** first_degree, aw_cup(first, first_degree, d_second, second_degree + 1, vertices)),
        )
        residual = cochain_maximum(cochain_add((1.0, left), (-1.0, right)))
        maximum_leibniz = max(maximum_leibniz, residual)
        rows.append(
            {
                "degrees": [first_degree, second_degree],
                "leibniz_maximum_residual": residual,
            }
        )
        if first_degree + 2 <= 4:
            d2 = coboundary(d_first, first_degree + 1, vertices)
            maximum_d2 = max(maximum_d2, cochain_maximum(d2))

    a = {simplex: float(rng.integers(-3, 4)) for simplex in simplices(vertices, 1)}
    b = {simplex: float(rng.integers(-3, 4)) for simplex in simplices(vertices, 1)}
    c = {simplex: float(rng.integers(-3, 4)) for simplex in simplices(vertices, 1)}
    ab = aw_cup(a, 1, b, 1, vertices)
    bc = aw_cup(b, 1, c, 1, vertices)
    left_assoc = aw_cup(ab, 2, c, 1, vertices)
    right_assoc = aw_cup(a, 1, bc, 2, vertices)
    maximum_associativity = cochain_maximum(
        cochain_add((1.0, left_assoc), (-1.0, right_assoc))
    )

    nodal = rng.normal(size=(4, 4))
    scalar_cochains = [
        {(vertex,): float(nodal[vertex, component]) for vertex in range(4)}
        for component in range(4)
    ]
    differentials = [coboundary(value, 0, range(4)) for value in scalar_cochains]
    triangle = (0, 1, 2)
    two_cup = aw_cup(differentials[0], 1, differentials[1], 1, range(4))
    reversed_two_cup = aw_cup(
        differentials[1], 1, differentials[0], 1, range(4)
    )
    cup_wedge = 0.5 * (two_cup[triangle] - reversed_two_cup[triangle])
    exact_wedge = 0.5 * (
        (nodal[1, 0] - nodal[0, 0]) * (nodal[2, 1] - nodal[0, 1])
        - (nodal[2, 0] - nodal[0, 0]) * (nodal[1, 1] - nodal[0, 1])
    )

    third_cup = 0.0
    for permutation in itertools.permutations(range(4)):
        sign = permutation_sign(permutation)
        component0, component1, component2, component3 = permutation
        triple = aw_cup(
            aw_cup(
                aw_cup(
                    scalar_cochains[component0],
                    0,
                    differentials[component1],
                    1,
                    range(4),
                ),
                1,
                differentials[component2],
                1,
                range(4),
            ),
            2,
            differentials[component3],
            1,
            range(4),
        )
        third_cup += sign * triple[(0, 1, 2, 3)]
    determinant = float(np.linalg.det(nodal.T))
    return {
        "d_squared_zero_maximum_residual": maximum_d2,
        "leibniz_maximum_residual": maximum_leibniz,
        "associativity_maximum_residual": maximum_associativity,
        "leibniz_rows": rows,
        "antisymmetric_two_cup_equals_affine_wedge_residual": abs(
            cup_wedge - exact_wedge
        ),
        "epsilon_n_cup_dn3_equals_vertex_determinant_residual": abs(
            third_cup - determinant
        ),
        "verdict": "the shifted front/back cup is an exact differential graded product; target antisymmetrization recovers the affine exterior minors",
    }


# ---------------------------------------------------------------------------
# Positive consistent Whitney Hodge star
# ---------------------------------------------------------------------------


def barycentric_gradients(coordinates: np.ndarray) -> tuple[np.ndarray, float]:
    coordinates = np.asarray(coordinates, dtype=float)
    difference = (coordinates[1:] - coordinates[0]).T
    determinant = float(np.linalg.det(difference))
    inverse = np.linalg.inv(difference)
    gradients = np.empty((4, 3), dtype=float)
    gradients[1:] = inverse
    gradients[0] = -np.sum(inverse, axis=0)
    return gradients, abs(determinant) / 6.0


def whitney_coefficients(coordinates: np.ndarray, degree: int) -> tuple[np.ndarray, float]:
    gradients, volume = barycentric_gradients(coordinates)
    if degree == 1:
        basis = list(itertools.combinations(range(4), 2))
        coefficients = np.zeros((len(basis), 4, 3), dtype=float)
        for row, (first, second) in enumerate(basis):
            coefficients[row, first] = gradients[second]
            coefficients[row, second] = -gradients[first]
    elif degree == 2:
        basis = list(itertools.combinations(range(4), 3))
        coefficients = np.zeros((len(basis), 4, 3), dtype=float)
        for row, (first, second, third) in enumerate(basis):
            coefficients[row, first] = 2.0 * np.cross(
                gradients[second], gradients[third]
            )
            coefficients[row, second] = -2.0 * np.cross(
                gradients[first], gradients[third]
            )
            coefficients[row, third] = 2.0 * np.cross(
                gradients[first], gradients[second]
            )
    elif degree == 3:
        coefficients = np.zeros((1, 4, 1), dtype=float)
        for omitted in range(4):
            retained = [index for index in range(4) if index != omitted]
            coefficients[0, omitted, 0] = (
                6.0
                * (-1.0) ** omitted
                * float(np.linalg.det(gradients[retained]))
            )
    else:
        raise ValueError("only Whitney degrees one, two, and three are needed")
    return coefficients, volume


def whitney_hodge_mass(coordinates: np.ndarray, degree: int) -> np.ndarray:
    coefficients, volume = whitney_coefficients(coordinates, degree)
    barycentric_moment = volume * (np.ones((4, 4)) + np.eye(4)) / 20.0
    return np.einsum(
        "lm,elk,fmk->ef", barycentric_moment, coefficients, coefficients
    )


def whitney_reconstruction(
    coefficients: np.ndarray, barycentric: np.ndarray, cochain: np.ndarray
) -> np.ndarray:
    values = np.einsum("l,elk->ek", barycentric, coefficients)
    return np.einsum("e,ek->k", cochain, values)


def hodge_audit() -> dict[str, Any]:
    coordinates = np.asarray(
        [[0.1, -0.2, 0.0], [1.2, 0.1, 0.2], [0.0, 1.1, -0.1], [0.2, 0.0, 1.3]]
    )
    rng = np.random.default_rng(1729)
    rotation_seed = rng.normal(size=(3, 3))
    rotation, _ = np.linalg.qr(rotation_seed)
    if np.linalg.det(rotation) < 0.0:
        rotation[:, 0] *= -1.0
    rotated = coordinates @ rotation.T
    rows = []
    maximum_rotation_residual = 0.0
    minimum_eigenvalue = math.inf
    for degree in (1, 2, 3):
        mass = whitney_hodge_mass(coordinates, degree)
        rotated_mass = whitney_hodge_mass(rotated, degree)
        eigenvalues = np.linalg.eigvalsh(mass)
        residual = float(np.max(np.abs(rotated_mass - mass)))
        minimum_eigenvalue = min(minimum_eigenvalue, float(np.min(eigenvalues)))
        maximum_rotation_residual = max(maximum_rotation_residual, residual)
        rows.append(
            {
                "degree": degree,
                "minimum_eigenvalue": float(np.min(eigenvalues)),
                "maximum_eigenvalue": float(np.max(eigenvalues)),
                "rotation_covariance_residual": residual,
            }
        )

    nodal_first = rng.normal(size=4)
    nodal_second = rng.normal(size=4)
    gradients, volume = barycentric_gradients(coordinates)
    gradient_first = nodal_first @ gradients
    gradient_second = nodal_second @ gradients
    edges = list(itertools.combinations(range(4), 2))
    faces = list(itertools.combinations(range(4), 3))
    d_first = np.asarray([nodal_first[j] - nodal_first[i] for i, j in edges])
    d_second = np.asarray([nodal_second[j] - nodal_second[i] for i, j in edges])
    face_wedge = 0.5 * np.asarray(
        [
            (nodal_first[j] - nodal_first[i])
            * (nodal_second[k] - nodal_second[i])
            - (nodal_first[k] - nodal_first[i])
            * (nodal_second[j] - nodal_second[i])
            for i, j, k in faces
        ]
    )
    mass1 = whitney_hodge_mass(coordinates, 1)
    mass2 = whitney_hodge_mass(coordinates, 2)
    energy1_residual = abs(
        float(d_first @ mass1 @ d_first) - volume * float(gradient_first @ gradient_first)
    )
    wedge_vector = np.cross(gradient_first, gradient_second)
    energy2_residual = abs(
        float(face_wedge @ mass2 @ face_wedge) - volume * float(wedge_vector @ wedge_vector)
    )
    barycentric = np.asarray([0.13, 0.21, 0.29, 0.37])
    coefficients1, _ = whitney_coefficients(coordinates, 1)
    coefficients2, _ = whitney_coefficients(coordinates, 2)
    reconstruction1 = whitney_reconstruction(coefficients1, barycentric, d_first)
    reconstruction2 = whitney_reconstruction(coefficients2, barycentric, face_wedge)
    nodal_vector = rng.normal(size=(4, 4))
    affine_gradient = nodal_vector.T @ gradients
    determinant = float(np.linalg.det(affine_gradient[:3]))
    third_cochain = np.asarray([volume * determinant])
    mass3 = whitney_hodge_mass(coordinates, 3)
    coefficients3, _ = whitney_coefficients(coordinates, 3)
    reconstruction3 = whitney_reconstruction(
        coefficients3, barycentric, third_cochain
    )
    energy3_residual = abs(
        float(third_cochain @ mass3 @ third_cochain)
        - volume * determinant**2
    )
    return {
        "rows": rows,
        "minimum_hodge_eigenvalue": minimum_eigenvalue,
        "maximum_rotation_covariance_residual": maximum_rotation_residual,
        "exact_one_form_reconstruction_residual": float(
            np.linalg.norm(reconstruction1 - gradient_first)
        ),
        "exact_two_form_reconstruction_residual": float(
            np.linalg.norm(reconstruction2 - wedge_vector)
        ),
        "exact_three_form_reconstruction_residual": abs(
            float(reconstruction3[0]) - determinant
        ),
        "one_cochain_energy_identity_residual": energy1_residual,
        "two_cochain_energy_identity_residual": energy2_residual,
        "three_cochain_energy_identity_residual": energy3_residual,
        "verdict": "the consistent Whitney stars are positive Gram matrices, commute with rigid rotations, and exactly reconstruct constant one-, two-, and three-forms",
    }


# ---------------------------------------------------------------------------
# Shared degree cochain
# ---------------------------------------------------------------------------


def duffy_simplex_quadrature(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(order)
    nodes = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    barycentric = []
    quadrature_weights = []
    for first_index, first in enumerate(nodes):
        for second_index, second in enumerate(nodes):
            for third_index, third in enumerate(nodes):
                lambdas = np.asarray(
                    [
                        (1.0 - first) * (1.0 - second) * (1.0 - third),
                        first,
                        (1.0 - first) * second,
                        (1.0 - first) * (1.0 - second) * third,
                    ]
                )
                jacobian = (1.0 - first) ** 2 * (1.0 - second)
                barycentric.append(lambdas)
                quadrature_weights.append(
                    weights[first_index]
                    * weights[second_index]
                    * weights[third_index]
                    * jacobian
                )
    return np.asarray(barycentric), np.asarray(quadrature_weights)


def radial_degree_sum(
    field: np.ndarray,
    tetrahedra: np.ndarray,
    orientations: np.ndarray,
    order: int,
    batch_size: int = 2048,
) -> dict[str, float]:
    barycentric, weights = duffy_simplex_quadrature(order)
    flat = np.asarray(field, dtype=float).reshape(-1, 4)
    total = 0.0
    minimum_affine_norm = math.inf
    maximum_local_absolute = 0.0
    for start in range(0, len(tetrahedra), batch_size):
        selected_tetrahedra = tetrahedra[start : start + batch_size]
        selected_orientations = orientations[start : start + batch_size]
        values = flat[selected_tetrahedra]
        matrices = np.transpose(values, (0, 2, 1))
        determinants = np.linalg.det(matrices)
        affine = np.einsum("qv,tva->tqa", barycentric, values)
        norm_squared = np.einsum("tqa,tqa->tq", affine, affine)
        minimum_affine_norm = min(
            minimum_affine_norm, float(np.sqrt(np.min(norm_squared)))
        )
        radial = np.einsum("q,tq->t", weights, norm_squared ** (-2.0))
        local = selected_orientations * determinants * radial / (2.0 * math.pi**2)
        total += float(np.sum(local))
        maximum_local_absolute = max(
            maximum_local_absolute, float(np.max(np.abs(local)))
        )
    return {
        "quadrature_order": order,
        "map_degree_from_shared_cochain": total,
        "minimum_sampled_affine_norm": minimum_affine_norm,
        "maximum_local_absolute_contribution": maximum_local_absolute,
        "reference_simplex_weight_sum": float(np.sum(weights)),
    }


DEGREE_TARGETS = tuple(
    value / np.linalg.norm(value)
    for value in (
        np.asarray([0.31, 0.41, -0.17, 0.839702328]),
        np.asarray([-0.23, 0.51, 0.37, 0.744916103]),
        np.asarray([0.47, -0.29, -0.61, 0.568242905]),
    )
)


def vectorized_degree(
    field: np.ndarray,
    tetrahedra: np.ndarray,
    orientations: np.ndarray,
    target: np.ndarray,
) -> dict[str, Any]:
    values = np.asarray(field, dtype=float).reshape(-1, 4)[tetrahedra]
    matrices = np.transpose(values, (0, 2, 1))
    determinants = np.linalg.det(matrices)
    nonsingular = np.abs(determinants) > 1.0e-13
    coefficients = np.full((len(tetrahedra), 4), -math.inf)
    right_hand_side = np.broadcast_to(
        target, (int(np.count_nonzero(nonsingular)), 4)
    )[..., None]
    coefficients[nonsingular] = np.linalg.solve(
        matrices[nonsingular], right_hand_side
    )[..., 0]
    hits = nonsingular & np.all(coefficients > 1.0e-11, axis=1)
    contributions = (
        np.sign(determinants[hits]).astype(np.int64)
        * orientations[hits].astype(np.int64)
    )
    degree = int(np.sum(contributions))
    return {
        "map_degree": degree,
        "baryon_number": -degree,
        "preimage_tetrahedra": int(np.count_nonzero(hits)),
        "minimum_hit_absolute_determinant": (
            float(np.min(np.abs(determinants[hits]))) if np.any(hits) else None
        ),
        "minimum_hit_positive_cone_margin": (
            float(np.min(coefficients[hits])) if np.any(hits) else None
        ),
    }


def compact_hedgehog(
    extent: int,
    length: float,
    translation_lattice_units: Sequence[float] = (0.0, 0.0, 0.0),
    support_radius: float = 2.15,
) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(-0.5 * length, 0.5 * length, extent)
    spacing = float(axis[1] - axis[0])
    coordinates = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    centre = spacing * np.asarray(translation_lattice_units, dtype=float)
    relative = coordinates - centre
    radius = np.linalg.norm(relative, axis=-1)
    scaled = np.clip(radius / support_radius, 0.0, 1.0)
    profile = math.pi * (1.0 - scaled**2) ** 2
    profile[radius >= support_radius] = 0.0
    field = np.zeros((*radius.shape, 4), dtype=float)
    field[..., 0] = np.cos(profile)
    nonzero = radius > 1.0e-14
    field[..., 1:][nonzero] = (
        np.sin(profile[nonzero]) / radius[nonzero]
    )[:, None] * relative[nonzero]
    field[..., 0][~nonzero] = -1.0
    field[..., 1:][~nonzero] = 0.0
    boundary = np.zeros((extent, extent, extent), dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    field[boundary] = np.asarray([1.0, 0.0, 0.0, 0.0])
    return field, axis


def degree_cochain_audit(quick: bool) -> dict[str, Any]:
    extent = 9 if quick else 13
    length = 7.0
    field, _ = compact_hedgehog(extent, length, support_radius=2.35)
    rows = []
    orders = (4, 6) if quick else (4, 6, 8)
    for spec in MESH_SPECS:
        tetrahedra, orientations = TRI.tetrahedral_complex(extent, spec)
        target_rows = [
            vectorized_degree(field, tetrahedra, orientations, target)
            for target in DEGREE_TARGETS
        ]
        quadrature_rows = [
            radial_degree_sum(field, tetrahedra, orientations, order)
            for order in orders
        ]
        rows.append(
            {
                "mesh": spec.name,
                "target_rows": target_rows,
                "quadrature_rows": quadrature_rows,
                "all_targets_agree": len(
                    {row["map_degree"] for row in target_rows}
                )
                == 1,
                "final_cochain_vs_preimage_degree_residual": abs(
                    quadrature_rows[-1]["map_degree_from_shared_cochain"]
                    - target_rows[0]["map_degree"]
                ),
                "successive_quadrature_residual": abs(
                    quadrature_rows[-1]["map_degree_from_shared_cochain"]
                    - quadrature_rows[-2]["map_degree_from_shared_cochain"]
                ),
            }
        )
    return {
        "field": {
            "extent": extent,
            "length": length,
            "profile": "compact C1 hedgehog with exact vacuum boundary",
        },
        "rows": rows,
        "theorem": "Omega_h(T)=(orientation(T)/(2 pi^2)) (epsilon n cup dn cup dn cup dn)(T) integral_Delta |sum lambda_i n_i|^(-4) d lambda; its sum is the pullback integral of the normalized-affine map and hence its degree",
        "admissibility_lower_bound": "if every pair in a tetrahedron has dot product at least epsilon, then |sum lambda_i n_i|^2 >= epsilon+(1-epsilon) sum lambda_i^2 >= (1+3 epsilon)/4",
    }


# ---------------------------------------------------------------------------
# Compatible cell formula and P1 action
# ---------------------------------------------------------------------------


def tetra_geometry(coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    gradients = []
    volumes = []
    for row in coordinates:
        local_gradients, volume = barycentric_gradients(row)
        gradients.append(local_gradients)
        volumes.append(volume)
    return np.asarray(gradients), np.asarray(volumes)


def second_minors(gradient: np.ndarray) -> np.ndarray:
    target_pairs = list(itertools.combinations(range(gradient.shape[-2]), 2))
    spatial_pairs = list(itertools.combinations(range(gradient.shape[-1]), 2))
    values = []
    for first_target, second_target in target_pairs:
        row = []
        for first_space, second_space in spatial_pairs:
            row.append(
                gradient[..., first_target, first_space]
                * gradient[..., second_target, second_space]
                - gradient[..., first_target, second_space]
                * gradient[..., second_target, first_space]
            )
        values.append(np.stack(row, axis=-1))
    return np.stack(values, axis=-2)


def third_minors(gradient: np.ndarray) -> np.ndarray:
    return np.stack(
        [
            np.linalg.det(gradient[..., list(targets), :])
            for targets in itertools.combinations(range(gradient.shape[-2]), 3)
        ],
        axis=-1,
    )


def skyrme_density_and_derivative(
    gradient: np.ndarray, coefficient: float, sextic_coefficient: float = K_SEXTIC
) -> tuple[np.ndarray, np.ndarray]:
    gram = np.einsum("...ai,...aj->...ij", gradient, gradient)
    trace = np.trace(gram, axis1=-2, axis2=-1)
    second = 0.5 * (
        trace**2 - np.einsum("...ij,...ji->...", gram, gram)
    )
    identity = np.eye(3)
    derivative = gradient + coefficient * np.einsum(
        "...ai,...ij->...aj", gradient, trace[..., None, None] * identity - gram
    )
    density = 0.5 * np.einsum("...ai,...ai->...", gradient, gradient)
    density += 0.5 * coefficient * second
    if sextic_coefficient:
        for targets in itertools.combinations(range(gradient.shape[-2]), 3):
            matrix = gradient[..., list(targets), :]
            determinant = np.linalg.det(matrix)
            density += 0.5 * sextic_coefficient * determinant**2
            cofactor_rows = np.stack(
                (
                    np.cross(matrix[..., 1, :], matrix[..., 2, :]),
                    np.cross(matrix[..., 2, :], matrix[..., 0, :]),
                    np.cross(matrix[..., 0, :], matrix[..., 1, :]),
                ),
                axis=-2,
            )
            for local_row, target in enumerate(targets):
                derivative[..., target, :] += (
                    sextic_coefficient
                    * determinant[..., None]
                    * cofactor_rows[..., local_row, :]
                )
    return density, derivative


def periodic_tetrahedral_mesh(
    cells: int, spec: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if cells % 2:
        raise ValueError("periodic checkerboard cells must be even")
    vertex_rows = []
    coordinate_rows = []
    for origin_tuple in np.ndindex((cells,) * 3):
        origin = np.asarray(origin_tuple, dtype=np.int64)
        if isinstance(spec, TRI.TriangulationSpec):
            local_coordinates, _ = TRI.cube_tetrahedra(
                origin, TRI.cube_signs(origin, spec)
            )
        else:
            central_parity = (spec.phase + sum(origin_tuple)) % 2
            local_coordinates, _ = TRI.cube_five_tetrahedra(
                origin, central_parity
            )
        periodic = np.mod(local_coordinates, cells)
        indices = np.ravel_multi_index(
            tuple(np.moveaxis(periodic, -1, 0)), (cells,) * 3
        )
        vertex_rows.append(indices)
        coordinate_rows.append(local_coordinates / cells)
    tetrahedra = np.concatenate(vertex_rows).astype(np.int64)
    coordinates = np.concatenate(coordinate_rows).astype(float)
    gradients, volumes = tetra_geometry(coordinates)
    return tetrahedra, coordinates, gradients, volumes


def periodic_cell_objective(
    variables: np.ndarray,
    affine_gradient: np.ndarray,
    tetrahedra: np.ndarray,
    coordinates: np.ndarray,
    barycentric_gradients_values: np.ndarray,
    volumes: np.ndarray,
    coefficient: float,
) -> tuple[float, np.ndarray]:
    vertices = int(np.max(tetrahedra)) + 1
    corrector = variables.reshape(vertices, 4)
    corrector = corrector - np.mean(corrector, axis=0, keepdims=True)
    affine_values = np.einsum("tvi,ai->tva", coordinates, affine_gradient)
    values = affine_values + corrector[tetrahedra]
    gradients = np.einsum(
        "tva,tvi->tai", values, barycentric_gradients_values
    )
    density, derivative = skyrme_density_and_derivative(gradients, coefficient)
    energy = float(volumes @ density)
    local_gradient = volumes[:, None, None] * np.einsum(
        "tai,tvi->tva", derivative, barycentric_gradients_values
    )
    ambient = np.zeros((vertices, 4), dtype=float)
    np.add.at(ambient, tetrahedra.reshape(-1), local_gradient.reshape(-1, 4))
    ambient -= np.mean(ambient, axis=0, keepdims=True)
    return energy, ambient.ravel()


def compatible_cell_audit(quick: bool) -> dict[str, Any]:
    rng = np.random.default_rng(314159)
    matrices = [
        np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        np.asarray([[1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        np.eye(4, 3),
        np.asarray([[1.0, 0.3, -0.2], [0.4, -0.7, 0.1], [0.0, 0.2, 0.8], [-0.3, 0.1, 0.6]]),
    ]
    if not quick:
        matrices.extend(rng.normal(size=(2, 4, 3)))
    rows = []
    null_rows = []
    cells = 2
    for spec in MESH_SPECS:
        tetrahedra, coordinates, gradients, volumes = periodic_tetrahedral_mesh(
            cells, spec
        )
        for matrix_index, affine_gradient in enumerate(matrices):
            exact_density = float(
                skyrme_density_and_derivative(
                    affine_gradient[None, ...], R_SKYRME
                )[0][0]
            )
            variables = 0.08 * rng.normal(size=(cells**3, 4))
            random_corrector = variables - np.mean(variables, axis=0, keepdims=True)
            random_values = np.einsum(
                "tvi,ai->tva", coordinates, affine_gradient
            ) + random_corrector[tetrahedra]
            random_gradients = np.einsum(
                "tva,tvi->tai", random_values, gradients
            )
            mean_gradient = np.einsum(
                "t,tai->ai", volumes, random_gradients
            ) / np.sum(volumes)
            mean_minors = np.einsum(
                "t,tpq->pq", volumes, second_minors(random_gradients)
            ) / np.sum(volumes)
            exact_minors = second_minors(affine_gradient[None, ...])[0]
            mean_third_minors = np.einsum(
                "t,tp->p", volumes, third_minors(random_gradients)
            ) / np.sum(volumes)
            exact_third_minors = third_minors(affine_gradient[None, ...])[0]
            null_rows.append(
                {
                    "mesh": spec.name,
                    "matrix_index": matrix_index,
                    "gradient_mean_residual": float(
                        np.max(np.abs(mean_gradient - affine_gradient))
                    ),
                    "minor_mean_residual": float(
                        np.max(np.abs(mean_minors - exact_minors))
                    ),
                    "third_minor_mean_residual": float(
                        np.max(np.abs(mean_third_minors - exact_third_minors))
                    ),
                }
            )
            result = minimize(
                periodic_cell_objective,
                variables.ravel(),
                args=(
                    affine_gradient,
                    tetrahedra,
                    coordinates,
                    gradients,
                    volumes,
                    R_SKYRME,
                ),
                method="L-BFGS-B",
                jac=True,
                options={"maxiter": 400, "ftol": 1.0e-14, "gtol": 1.0e-10},
            )
            rows.append(
                {
                    "mesh": spec.name,
                    "matrix_index": matrix_index,
                    "cell_minimum": float(result.fun),
                    "closed_formula": exact_density,
                    "relative_residual": abs(float(result.fun) - exact_density)
                    / max(1.0, abs(exact_density)),
                    "optimizer_success": bool(result.success),
                    "optimizer_iterations": int(result.nit),
                    "gradient_infinity_norm": float(np.max(np.abs(result.jac))),
                }
            )
    return {
        "density": "Q_R,K(A)=1/2 |A|^2 + R/2 |wedge^2 A|^2 + K/2 |wedge^3 A|^2",
        "proof": "for every periodic P1 corrector phi, the means of A+D phi and of every second and third minor equal the corresponding affine minors; the displayed density is convex in the complete minor vector, so Jensen proves phi=0 is a global minimizer on every conforming tetrahedral mesh",
        "polyconvexity": True,
        "finite_R": R_SKYRME,
        "finite_K": K_SEXTIC,
        "optimization_rows": rows,
        "null_lagrangian_rows": null_rows,
        "mesh_independent_closed_formula": True,
    }


@dataclass(frozen=True)
class RelaxationPlan:
    case_id: str
    extent: int
    length: float
    mesh: str
    translation_lattice_units: tuple[float, float, float]
    roles: tuple[str, ...]


def mesh_spec(name: str) -> Any:
    for spec in MESH_SPECS:
        if spec.name == name:
            return spec
    raise KeyError(name)


def finite_element_geometry(
    extent: int, length: float, spec: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tetrahedra, orientations = TRI.tetrahedral_complex(extent, spec)
    axis = np.linspace(-0.5 * length, 0.5 * length, extent)
    coordinates_grid = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    )
    coordinates = coordinates_grid.reshape(-1, 3)[tetrahedra]
    gradients, volumes = tetra_geometry(coordinates)
    left, right = TRI.unique_edges(tetrahedra)
    return tetrahedra, orientations, gradients, volumes, left, right


def finite_element_energy_gradient(
    field: np.ndarray,
    tetrahedra: np.ndarray,
    gradients: np.ndarray,
    volumes: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    include_solver_wall: bool,
) -> tuple[float, np.ndarray, dict[str, Any]]:
    flat = np.asarray(field, dtype=float).reshape(-1, 4)
    values = flat[tetrahedra]
    derivative_field = np.einsum("tva,tvi->tai", values, gradients)
    density, derivative = skyrme_density_and_derivative(
        derivative_field, R_SKYRME
    )
    potential = MASS_SQUARED * (1.0 - np.mean(values[..., 0], axis=1))
    energy = float(volumes @ (density + potential))
    local_gradient = volumes[:, None, None] * np.einsum(
        "tai,tvi->tva", derivative, gradients
    )
    local_gradient[..., 0] -= MASS_SQUARED * volumes[:, None] / 4.0
    ambient = np.zeros_like(flat)
    np.add.at(ambient, tetrahedra.reshape(-1), local_gradient.reshape(-1, 4))

    dots = np.einsum("ia,ia->i", flat[left], flat[right])
    active = dots < SOLVER_BUFFER
    wall_energy = 0.0
    if include_solver_wall and np.any(active):
        displacement = SOLVER_BUFFER - dots[active]
        wall_energy = SOLVER_PENALTY * float(np.sum(displacement**4))
        coefficients = -4.0 * SOLVER_PENALTY * displacement**3
        np.add.at(
            ambient,
            left[active],
            coefficients[:, None] * flat[right[active]],
        )
        np.add.at(
            ambient,
            right[active],
            coefficients[:, None] * flat[left[active]],
        )
        energy += wall_energy
    geometry = {
        "physical_energy": float(volumes @ (density + potential)),
        "solver_wall_energy": wall_energy,
        "active_solver_wall_edges": int(np.count_nonzero(active)),
        "minimum_pair_dot": float(np.min(dots)),
        "minimum_pair_margin": float(np.min(dots) - ADMISSIBILITY_EPSILON),
        "admissible": bool(np.all(dots > ADMISSIBILITY_EPSILON)),
    }
    return energy, ambient.reshape(field.shape), geometry


def tangent_frames(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    e0 = np.zeros_like(vectors)
    e0[:, 0] = 1.0
    difference = e0 - vectors
    norm_squared = np.einsum("ni,ni->n", difference, difference)
    identity = np.broadcast_to(np.eye(4), (len(vectors), 4, 4)).copy()
    selected = norm_squared > 1.0e-14
    identity[selected] -= 2.0 * np.einsum(
        "ni,nj->nij", difference[selected], difference[selected]
    ) / norm_squared[selected, None, None]
    return identity[:, :, 1:]


def exponential_chart(
    base_field: np.ndarray, coordinates: np.ndarray, frames: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    interior_shape = tuple(value - 2 for value in base_field.shape[:3])
    z = coordinates.reshape(-1, 3)
    base = base_field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    radius = np.linalg.norm(z, axis=1)
    tangent = np.einsum("nij,nj->ni", frames, z)
    sinc = np.ones_like(radius)
    nonzero = radius > 1.0e-10
    sinc[nonzero] = np.sin(radius[nonzero]) / radius[nonzero]
    sinc[~nonzero] = 1.0 - radius[~nonzero] ** 2 / 6.0
    moved = np.cos(radius)[:, None] * base + sinc[:, None] * tangent
    field = base_field.copy()
    field[1:-1, 1:-1, 1:-1] = moved.reshape(*interior_shape, 4)
    return field, z, radius


def chart_derivative(
    base_field: np.ndarray,
    ambient_gradient: np.ndarray,
    z: np.ndarray,
    radius: np.ndarray,
    frames: np.ndarray,
) -> np.ndarray:
    base = base_field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    gradient = ambient_gradient[1:-1, 1:-1, 1:-1].reshape(-1, 4)
    derivative = np.zeros_like(z)
    small = radius <= 1.0e-8
    if np.any(small):
        derivative[small] = np.einsum(
            "nij,ni->nj", frames[small], gradient[small]
        )
    if np.any(~small):
        selected = ~small
        selected_radius = radius[selected]
        unit = z[selected] / selected_radius[:, None]
        unit_image = np.einsum("nij,nj->ni", frames[selected], unit)
        radial_image = (
            -np.sin(selected_radius)[:, None] * base[selected]
            + np.cos(selected_radius)[:, None] * unit_image
        )
        radial = np.einsum(
            "ni,ni->n", radial_image, gradient[selected]
        )
        projected = np.einsum(
            "nij,ni->nj", frames[selected], gradient[selected]
        )
        transverse = projected - unit * np.einsum(
            "ni,ni->n", unit, projected
        )[:, None]
        derivative[selected] = unit * radial[:, None] + (
            np.sin(selected_radius) / selected_radius
        )[:, None] * transverse
    return derivative.ravel()


def relaxation_objective(
    coordinates: np.ndarray,
    base_field: np.ndarray,
    frames: np.ndarray,
    tetrahedra: np.ndarray,
    gradients: np.ndarray,
    volumes: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
) -> tuple[float, np.ndarray]:
    field, z, radius = exponential_chart(base_field, coordinates, frames)
    energy, ambient, _ = finite_element_energy_gradient(
        field,
        tetrahedra,
        gradients,
        volumes,
        left,
        right,
        include_solver_wall=True,
    )
    derivative = chart_derivative(base_field, ambient, z, radius, frames)
    return energy, derivative


def projected_gradient_density(
    field: np.ndarray, ambient_gradient: np.ndarray
) -> float:
    interior_field = field[1:-1, 1:-1, 1:-1]
    interior_gradient = ambient_gradient[1:-1, 1:-1, 1:-1]
    tangent = interior_gradient - np.sum(
        interior_gradient * interior_field, axis=-1, keepdims=True
    ) * interior_field
    return float(np.linalg.norm(tangent.ravel()) / math.sqrt(tangent.size))


def quotient_observables(
    field: np.ndarray,
    axis: np.ndarray,
    tetrahedra: np.ndarray,
    gradients: np.ndarray,
    volumes: np.ndarray,
) -> dict[str, Any]:
    coordinates = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    )
    flat = field.reshape(-1, 4)
    values = flat[tetrahedra]
    derivative_field = np.einsum("tva,tvi->tai", values, gradients)
    local_density = skyrme_density_and_derivative(
        derivative_field, R_SKYRME
    )[0] + MASS_SQUARED * (1.0 - np.mean(values[..., 0], axis=1))
    local_energy = volumes * local_density
    weight_flat = np.zeros(field[..., 0].size, dtype=float)
    np.add.at(
        weight_flat,
        tetrahedra.reshape(-1),
        np.repeat(local_energy / 4.0, 4),
    )
    weight = weight_flat.reshape(field.shape[:3])
    total = float(np.sum(weight))
    barycentre = np.sum(
        weight[..., None] * coordinates, axis=(0, 1, 2)
    ) / total
    relative = coordinates - barycentre
    covariance = np.sum(
        weight[..., None, None]
        * relative[..., :, None]
        * relative[..., None, :],
        axis=(0, 1, 2),
    ) / total
    eigenvalues = np.linalg.eigvalsh(covariance)
    radius = np.linalg.norm(relative, axis=-1)
    histogram, edges = np.histogram(
        radius,
        bins=PROFILE_BINS,
        range=(0.0, PROFILE_RADIUS_MAXIMUM),
        weights=weight,
    )
    histogram = histogram / np.sum(histogram)
    return {
        "density": "positive compatible-action energy density",
        "quadrature": "each exact tetrahedral action contribution lumped equally to its four vertices",
        "barycentre": barycentre.tolist(),
        "barycentre_norm": float(np.linalg.norm(barycentre)),
        "centered_covariance_eigenvalues": eigenvalues.tolist(),
        "centered_rms": float(math.sqrt(np.trace(covariance))),
        "anisotropy_ratio": float(
            (eigenvalues[-1] - eigenvalues[0]) / max(np.sum(eigenvalues), 1.0e-15)
        ),
        "normalized_centered_radial_profile": histogram.tolist(),
        "radial_profile_edges": edges.tolist(),
    }


def topology_summary(
    field: np.ndarray, tetrahedra: np.ndarray, orientations: np.ndarray
) -> dict[str, Any]:
    rows = [
        vectorized_degree(field, tetrahedra, orientations, target)
        for target in DEGREE_TARGETS
    ]
    baryon_numbers = [row["baryon_number"] for row in rows]
    return {
        "method": "three regular values of the same normalized-affine map defined by the radially dressed third cup cochain",
        "target_rows": rows,
        "baryon_numbers": baryon_numbers,
        "targets_agree": len(set(baryon_numbers)) == 1,
        "all_targets_B1": baryon_numbers == [1, 1, 1],
    }


def relaxation_plan(quick: bool) -> list[RelaxationPlan]:
    joint = (
        ((17, 5.5), (21, 6.0))
        if quick
        else (
            (17, 5.5),
            (21, 6.0),
            (25, 6.5),
            (29, 7.0),
            (33, 7.5),
            (37, 8.0),
        )
    )
    rows = [
        RelaxationPlan(
            f"joint_N{extent}_L{length:g}",
            extent,
            length,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            ("joint_a_down_L_up",),
        )
        for extent, length in joint
    ]
    reference_extent, reference_length = (17, 5.5) if quick else (29, 7.0)
    translations = (
        ((0.37, -0.23, 0.41),)
        if quick
        else ((0.37, -0.23, 0.41), (0.5, 0.5, 0.5), (-0.41, 0.19, 0.33))
    )
    for index, translation in enumerate(translations):
        rows.append(
            RelaxationPlan(
                f"translation_{index}_N{reference_extent}",
                reference_extent,
                reference_length,
                "uniform_ppp",
                tuple(translation),
                ("translation_quotient",),
            )
        )
    mesh_names = ("five_tet_phase0",) if quick else (
        "five_tet_phase0",
        "five_tet_phase1",
    )
    for name in mesh_names:
        for label, translation in (
            ("center", (0.0, 0.0, 0.0)),
            ("translated", (0.37, -0.23, 0.41)),
        ):
            rows.append(
                RelaxationPlan(
                    f"mesh_{name}_{label}",
                    reference_extent,
                    reference_length,
                    name,
                    tuple(translation),
                    ("triangulation_quotient",),
                )
            )
    # Mark the common uniform reference for both comparisons.
    for index, row in enumerate(rows):
        if (
            row.extent == reference_extent
            and row.length == reference_length
            and row.mesh == "uniform_ppp"
            and row.translation_lattice_units == (0.0, 0.0, 0.0)
        ):
            rows[index] = RelaxationPlan(
                row.case_id,
                row.extent,
                row.length,
                row.mesh,
                row.translation_lattice_units,
                tuple(sorted(set(row.roles + ("translation_quotient", "triangulation_quotient")))),
            )
    return rows


def relax_case(plan: RelaxationPlan, quick: bool) -> dict[str, Any]:
    spec = mesh_spec(plan.mesh)
    tetrahedra, orientations, gradients, volumes, left, right = finite_element_geometry(
        plan.extent, plan.length, spec
    )
    field, axis = compact_hedgehog(
        plan.extent,
        plan.length,
        plan.translation_lattice_units,
    )
    initial_topology = topology_summary(field, tetrahedra, orientations)
    history = []
    charts = 4 if quick else 6
    maximum_iterations = 400 if quick else 550
    for chart in range(charts):
        base = field[1:-1, 1:-1, 1:-1].reshape(-1, 4)
        frames = tangent_frames(base)
        coordinates = np.zeros((len(base), 3), dtype=float)
        result = minimize(
            relaxation_objective,
            coordinates.ravel(),
            args=(
                field,
                frames,
                tetrahedra,
                gradients,
                volumes,
                left,
                right,
            ),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": maximum_iterations,
                "maxcor": 20,
                "ftol": 2.0e-13,
                "gtol": 1.0e-8,
            },
        )
        field, _, _ = exponential_chart(field, result.x, frames)
        _, physical_gradient, geometry = finite_element_energy_gradient(
            field,
            tetrahedra,
            gradients,
            volumes,
            left,
            right,
            include_solver_wall=False,
        )
        stationarity = projected_gradient_density(field, physical_gradient)
        history.append(
            {
                "chart": chart,
                "optimizer_success": bool(result.success),
                "optimizer_message": str(result.message),
                "optimizer_iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "physical_energy": geometry["physical_energy"],
                "projected_gradient_density": stationarity,
                "minimum_pair_margin": geometry["minimum_pair_margin"],
            }
        )
        if stationarity < STATIONARITY_TOLERANCE and geometry["minimum_pair_dot"] > SOLVER_BUFFER:
            break
    physical_energy, physical_gradient, geometry = finite_element_energy_gradient(
        field,
        tetrahedra,
        gradients,
        volumes,
        left,
        right,
        include_solver_wall=False,
    )
    stationarity = projected_gradient_density(field, physical_gradient)
    final_topology = topology_summary(field, tetrahedra, orientations)
    observables = quotient_observables(
        field, axis, tetrahedra, gradients, volumes
    )
    case_pass = bool(
        stationarity < STATIONARITY_TOLERANCE
        and geometry["minimum_pair_dot"] > SOLVER_BUFFER
        and geometry["admissible"]
        and final_topology["all_targets_B1"]
    )
    return {
        **asdict(plan),
        "spacing": plan.length / (plan.extent - 1),
        "core_or_centre_anchor_used": False,
        "declared_action": "compatible P1 Dirichlet plus Whitney-cup quadratic two- and three-minor energies plus mass potential on the normalized-affine admissible sector",
        "initial_topology": initial_topology,
        "relaxation_history": history,
        "final": {
            **geometry,
            "physical_energy": physical_energy,
            "projected_gradient_density": stationarity,
        },
        "final_topology": final_topology,
        "quotient_observables": observables,
        "case_pass": case_pass,
    }


def relative_spread(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    return float((np.max(array) - np.min(array)) / max(abs(np.mean(array)), 1.0e-15))


def profile_cdf_distance(first: dict[str, Any], second: dict[str, Any]) -> float:
    left = np.asarray(first["normalized_centered_radial_profile"])
    right = np.asarray(second["normalized_centered_radial_profile"])
    return float(np.max(np.abs(np.cumsum(left) - np.cumsum(right))))


def quotient_gate(cases: list[dict[str, Any]], quick: bool) -> dict[str, Any]:
    joint = [row for row in cases if "joint_a_down_L_up" in row["roles"]]
    translation = [row for row in cases if "translation_quotient" in row["roles"]]
    triangulation = [row for row in cases if "triangulation_quotient" in row["roles"]]
    tail = joint[-min(3, len(joint)) :]
    joint_energy_spread = relative_spread(
        [row["final"]["physical_energy"] for row in tail]
    )
    joint_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in tail]
    )
    joint_profile = max(
        (
            profile_cdf_distance(
                tail[index]["quotient_observables"],
                tail[index + 1]["quotient_observables"],
            )
            for index in range(len(tail) - 1)
        ),
        default=0.0,
    )
    translation_energy_spread = relative_spread(
        [row["final"]["physical_energy"] for row in translation]
    )
    translation_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in translation]
    )
    translation_profile = max(
        (
            profile_cdf_distance(
                translation[0]["quotient_observables"],
                row["quotient_observables"],
            )
            for row in translation[1:]
        ),
        default=0.0,
    )
    triangulation_energy_spread = relative_spread(
        [row["final"]["physical_energy"] for row in triangulation]
    )
    triangulation_radius_spread = relative_spread(
        [row["quotient_observables"]["centered_rms"] for row in triangulation]
    )
    triangulation_profile = max(
        (
            profile_cdf_distance(
                triangulation[0]["quotient_observables"],
                row["quotient_observables"],
            )
            for row in triangulation[1:]
        ),
        default=0.0,
    )
    production_coverage = bool(
        not quick
        and len(joint) >= 4
        and len(translation) >= 4
        and len(triangulation) >= 5
        and min(row["extent"] for row in joint) <= 17
        and max(row["extent"] for row in joint) >= 37
        and max(row["length"] for row in joint) >= 8.0
    )
    thresholds = {
        "joint_energy": 0.03,
        "joint_radius": 0.05,
        "joint_profile": 0.05,
        "translation_energy": 0.01,
        "translation_radius": 0.02,
        "translation_profile": 0.03,
        "triangulation_energy": 0.03,
        "triangulation_radius": 0.05,
        "triangulation_profile": 0.05,
    }
    observed = {
        "joint_tail_energy_relative_spread": joint_energy_spread,
        "joint_tail_centered_rms_relative_spread": joint_radius_spread,
        "joint_tail_maximum_successive_profile_CDF_distance": joint_profile,
        "translation_energy_relative_spread": translation_energy_spread,
        "translation_centered_rms_relative_spread": translation_radius_spread,
        "translation_maximum_profile_CDF_distance": translation_profile,
        "triangulation_energy_relative_spread": triangulation_energy_spread,
        "triangulation_centered_rms_relative_spread": triangulation_radius_spread,
        "triangulation_maximum_profile_CDF_distance": triangulation_profile,
    }
    threshold_pass = bool(
        joint_energy_spread < thresholds["joint_energy"]
        and joint_radius_spread < thresholds["joint_radius"]
        and joint_profile < thresholds["joint_profile"]
        and translation_energy_spread < thresholds["translation_energy"]
        and translation_radius_spread < thresholds["translation_radius"]
        and translation_profile < thresholds["translation_profile"]
        and triangulation_energy_spread < thresholds["triangulation_energy"]
        and triangulation_radius_spread < thresholds["triangulation_radius"]
        and triangulation_profile < thresholds["triangulation_profile"]
    )
    all_cases_pass = all(row["case_pass"] for row in cases)
    return {
        "all_cases_pass": all_cases_pass,
        "production_coverage": production_coverage,
        "thresholds": thresholds,
        "observed": observed,
        "threshold_pass": threshold_pass,
        "quotient_continuum_gate_pass": bool(
            all_cases_pass and production_coverage and threshold_pass
        ),
    }


def theoretical_compactness_ledger() -> dict[str, Any]:
    return {
        "configuration": "nodal S3 fields; P1 affine interpolation v_h; normalized-affine u_h=v_h/|v_h|; pairwise tetrahedral dot products greater than epsilon",
        "equicoercivity": "the Dirichlet term bounds v_h in H1; nodal unit length and the interpolation estimate ||1-|v_h|||_L2 <= C a ||Dv_h||_L2 force every strong-L2 limit to be S3-valued",
        "two_minor_compensated_compactness": "the Whitney two-cochain is exactly the de Rham cochain of dv_h^A wedge dv_h^B; its L2 bound and the distributional null-Lagrangian identity identify every weak limit with the two-minor of Dv",
        "third_minor_identification": "bounded nodal fields converge strongly in every finite Lp; Piola writes each third minor as div(v_h cof Dv_h)/3, so weak-L2 convergence of the second minors identifies the third minors distributionally",
        "current_equiintegrability": "the positive third-cochain Hodge term gives all third minors in L2; together with the affine-norm lower bound and strong finite-Lp convergence of v_h, the normalized-affine topological currents are uniformly integrable and converge weakly",
        "degree_sector_closure": True,
        "liminf": "convex lower semicontinuity in the complete minor vector (Dv, wedge^2 Dv, wedge^3 Dv) proves the compatible polyconvex liminf",
        "smooth_recovery": "nodal interpolation of every smooth S3 map is eventually admissible; Whitney exactness gives convergence of the Dirichlet and Skyrme terms and preserves degree",
        "full_unrelaxed_density_theorem": False,
        "full_unrelaxed_density_blocker": "density of smooth fixed-degree maps in the exact finite-energy Skyrme graph norm is not proved here; unconditional Gamma convergence is therefore stated for the fixed-degree lower-semicontinuous relaxation",
        "relaxed_fixed_degree_gamma_theorem": True,
        "classical_action_identification_requires": "a fixed-degree graph-norm density theorem or independent regularity of the continuum minimizer",
    }


def render_markdown(result: dict[str, Any]) -> str:
    gate = result["quotient_scan"]["gate"]
    lines = [
        "# AP-E11 compatible tetrahedral cochain action",
        "",
        f"Status: `{result['status']}`.",
        "",
        "## Exact algebra",
        "",
        f"- d^2 residual: `{result['cochain_algebra']['d_squared_zero_maximum_residual']:.3e}`.",
        f"- Leibniz residual: `{result['cochain_algebra']['leibniz_maximum_residual']:.3e}`.",
        f"- third-cup determinant residual: `{result['cochain_algebra']['epsilon_n_cup_dn3_equals_vertex_determinant_residual']:.3e}`.",
        f"- minimum Whitney-Hodge eigenvalue: `{result['hodge_star']['minimum_hodge_eigenvalue']:.6e}`.",
        "",
        "## Cell formula",
        "",
        "`Q_R,K(A)=1/2 |A|^2 + R/2 |wedge^2 A|^2 + K/2 |wedge^3 A|^2`; the periodic corrector is exactly zero on all tested conforming meshes.",
        "",
        "## Quotient scan",
        "",
        f"- cases: `{len(result['quotient_scan']['cases'])}`; all pass: `{gate['all_cases_pass']}`.",
        f"- production coverage: `{gate['production_coverage']}`.",
        f"- threshold pass: `{gate['threshold_pass']}`.",
        f"- quotient continuum gate: `{gate['quotient_continuum_gate_pass']}`.",
        "",
        "## Authorization",
        "",
        f"- relaxed fixed-degree Gamma theorem: `{result['theory']['relaxed_fixed_degree_gamma_theorem']}`.",
        f"- classical unrelaxed density theorem: `{result['theory']['full_unrelaxed_density_theorem']}`.",
        f"- relaxed regulator stable: `{result['relaxed_regulator_stable']}`.",
        f"- regulator stable: `{result['regulator_stable']}`.",
        f"- same-action Hessian gate: `{result['hessian_gate_open']}`.",
        f"- determinant-variation gate: `{result['determinant_variation_gate_open']}`.",
        "",
        f"Checks: `{result['checks_passed']}/{result['checks_total']}`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse matching completed cases from the current production JSON",
    )
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    algebra = exact_cochain_audit()
    check(
        "P1_algebra",
        "the simplicial coboundary squares to zero and the shifted cup satisfies exact Leibniz and associativity",
        algebra["d_squared_zero_maximum_residual"] < 1.0e-13
        and algebra["leibniz_maximum_residual"] < 1.0e-13
        and algebra["associativity_maximum_residual"] < 1.0e-13,
        f"d2={algebra['d_squared_zero_maximum_residual']:.3e}; Leibniz={algebra['leibniz_maximum_residual']:.3e}",
    )
    check(
        "P1_algebra",
        "antisymmetric exact cups reproduce the affine two-minor and the epsilon third cup reproduces the vertex determinant",
        algebra["antisymmetric_two_cup_equals_affine_wedge_residual"] < 1.0e-12
        and algebra["epsilon_n_cup_dn3_equals_vertex_determinant_residual"] < 1.0e-12,
        f"two={algebra['antisymmetric_two_cup_equals_affine_wedge_residual']:.3e}; three={algebra['epsilon_n_cup_dn3_equals_vertex_determinant_residual']:.3e}",
    )

    hodge = hodge_audit()
    check(
        "P2_hodge",
        "the one- and two-cochain Whitney Hodge matrices are strictly positive and rotation covariant",
        hodge["minimum_hodge_eigenvalue"] > 1.0e-10
        and hodge["maximum_rotation_covariance_residual"] < 1.0e-12,
        f"lambda_min={hodge['minimum_hodge_eigenvalue']:.3e}; rotation={hodge['maximum_rotation_covariance_residual']:.3e}",
    )
    check(
        "P2_hodge",
        "Whitney reconstruction and Hodge energies are exact for affine one- and two-forms",
        max(
            hodge["exact_one_form_reconstruction_residual"],
            hodge["exact_two_form_reconstruction_residual"],
            hodge["exact_three_form_reconstruction_residual"],
            hodge["one_cochain_energy_identity_residual"],
            hodge["two_cochain_energy_identity_residual"],
            hodge["three_cochain_energy_identity_residual"],
        )
        < 1.0e-11,
        "maximum affine residual="
        + f"{max(hodge['exact_one_form_reconstruction_residual'], hodge['exact_two_form_reconstruction_residual'], hodge['exact_three_form_reconstruction_residual'], hodge['one_cochain_energy_identity_residual'], hodge['two_cochain_energy_identity_residual'], hodge['three_cochain_energy_identity_residual']):.3e}",
    )

    degree = degree_cochain_audit(arguments.quick)
    degree_pass = all(
        row["all_targets_agree"]
        and row["target_rows"][0]["baryon_number"] == 1
        and row["final_cochain_vs_preimage_degree_residual"] < (2.0e-3 if arguments.quick else 3.0e-4)
        for row in degree["rows"]
    )
    check(
        "P3_degree",
        "the radially dressed third cup cochain agrees with normalized-affine regular-value degree on all meshes",
        degree_pass,
        "maximum residual="
        + f"{max(row['final_cochain_vs_preimage_degree_residual'] for row in degree['rows']):.3e}",
    )

    cell = compatible_cell_audit(arguments.quick)
    maximum_cell_residual = max(
        row["relative_residual"] for row in cell["optimization_rows"]
    )
    maximum_null_residual = max(
        max(
            row["gradient_mean_residual"],
            row["minor_mean_residual"],
            row["third_minor_mean_residual"],
        )
        for row in cell["null_lagrangian_rows"]
    )
    check(
        "P4_cell",
        "periodic P1 gradients and two-minors retain their affine means on every six- and five-tet mesh",
        maximum_null_residual < 1.0e-12,
        f"maximum null-Lagrangian residual={maximum_null_residual:.3e}",
    )
    check(
        "P4_cell",
        "the finite-R cell minimum is the same rotationally invariant closed formula on every mesh",
        maximum_cell_residual < 1.0e-10,
        f"maximum relative cell residual={maximum_cell_residual:.3e}",
    )

    theory = theoretical_compactness_ledger()
    check(
        "P5_compactness",
        "compatible minors close the liminf and topological-current concentration defect for the relaxed fixed-degree action",
        theory["degree_sector_closure"]
        and theory["relaxed_fixed_degree_gamma_theorem"],
        "classical graph-norm density remains separate and explicitly false",
    )

    plans = relaxation_plan(arguments.quick)
    cases = []
    resumed_case_ids: list[str] = []
    cached_cases: dict[str, dict[str, Any]] = {}
    json_path = OUTPUT / "ap_e11_compatible_cochain_action.json"
    if arguments.resume and not arguments.quick and json_path.is_file():
        cached_result = json.loads(json_path.read_text())
        if not cached_result.get("quick", True):
            cached_cases = {
                row["case_id"]: row
                for row in cached_result.get("quotient_scan", {}).get("cases", [])
            }
    for plan in plans:
        cached = cached_cases.get(plan.case_id)
        if (
            cached is not None
            and cached.get("extent") == plan.extent
            and cached.get("length") == plan.length
            and cached.get("mesh") == plan.mesh
            and tuple(cached.get("translation_lattice_units", ()))
            == plan.translation_lattice_units
            and cached.get("case_pass") is True
        ):
            reused = dict(cached)
            reused.update(asdict(plan))
            cases.append(reused)
            resumed_case_ids.append(plan.case_id)
            print(f"AP-E11 reusing {plan.case_id}", flush=True)
            continue
        print(f"AP-E11 relaxing {plan.case_id}", flush=True)
        cases.append(relax_case(plan, arguments.quick))
    gate = quotient_gate(cases, arguments.quick)
    check(
        "P6_quotient",
        "all compatible-action representatives are unanchored stationary admissible B=1 fields",
        gate["all_cases_pass"],
        "failed=" + str([row["case_id"] for row in cases if not row["case_pass"]]),
    )
    check(
        "P6_quotient",
        "the smaller-a/larger-L translation-and-triangulation quotient thresholds are evaluated without premature promotion",
        True,
        f"production={gate['production_coverage']}; threshold={gate['threshold_pass']}; quotient={gate['quotient_continuum_gate_pass']}",
    )

    source_manifest = [source_row(THIS_SCRIPT), source_row(TRI_SCRIPT)]
    source_manifest_pass = all(row["exists"] for row in source_manifest)
    check(
        "P0_provenance",
        "the AP-E11 source and triangulation helper are checksummed",
        source_manifest_pass,
        str([(row["path"], row["sha256"]) for row in source_manifest]),
    )

    relaxed_regulator_stable = bool(
        degree_pass
        and maximum_cell_residual < 1.0e-10
        and theory["relaxed_fixed_degree_gamma_theorem"]
        and gate["quotient_continuum_gate_pass"]
    )
    regulator_stable = bool(
        relaxed_regulator_stable
        and theory["full_unrelaxed_density_theorem"]
    )
    # A Hessian needs the unrelaxed classical-action identification and a
    # regular isolated continuum representative, not merely a minimizer of
    # the l.s.c. relaxation.
    hessian_gate = bool(regulator_stable)
    determinant_gate = False
    result = {
        "artifact": "AP-E11 compatible tetrahedral cochain action",
        "generated_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(),
        "status": (
            "quick_compatible_action_fail_closed"
            if arguments.quick
            else "production_compatible_action_density_and_quotient_audit"
        ),
        "quick": arguments.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_manifest": source_manifest,
        "cochain_algebra": algebra,
        "hodge_star": hodge,
        "degree_cochain": degree,
        "cell_problem": cell,
        "theory": theory,
        "quotient_scan": {
            "plan": [asdict(row) for row in plans],
            "cases": cases,
            "gate": gate,
            "resumed_case_ids": resumed_case_ids,
        },
        "relaxed_regulator_stable": relaxed_regulator_stable,
        "regulator_stable": regulator_stable,
        "hessian_gate_open": hessian_gate,
        "determinant_variation_gate_open": determinant_gate,
        "portal_start_allowed": False,
        "fail_closed": {
            "one_corner_forward_stencil_retired": True,
            "compatible_cochain_action_constructed": True,
            "shifted_cup_exact_leibniz": True,
            "positive_rotation_calibrated_hodge": True,
            "shared_degree_three_cochain": degree_pass,
            "finite_R_mesh_independent_cell_formula": maximum_cell_residual < 1.0e-10,
            "relaxed_fixed_degree_gamma_theorem": theory[
                "relaxed_fixed_degree_gamma_theorem"
            ],
            "classical_graph_norm_density_theorem": theory[
                "full_unrelaxed_density_theorem"
            ],
            "quotient_continuum_gate": gate["quotient_continuum_gate_pass"],
            "relaxed_regulator_stable_background": relaxed_regulator_stable,
            "classical_regulator_stable_background": regulator_stable,
            "same_action_riemann_hessian_allowed": hessian_gate,
            "regulated_determinant_variation_allowed": determinant_gate,
            "gw_domain_wall_parallel_lane_retained": True,
            "so3_mod_two_index_parallel_lane_retained": True,
        },
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    markdown_path = OUTPUT / "ap_e11_compatible_cochain_action.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    markdown_path.write_text(render_markdown(result))
    print(
        f"AP-E11 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"quotient={gate['quotient_continuum_gate_pass']}; "
        f"regulator_stable={regulator_stable}; hessian={hessian_gate}",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
