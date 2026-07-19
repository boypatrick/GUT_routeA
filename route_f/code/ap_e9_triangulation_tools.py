#!/usr/bin/env python3
"""Conforming triangulations and translated radial starts for AP-E9.

The helpers in this module are deliberately independent of the AP-E8 driver.
They provide two families of six-tetrahedra-per-cube Kuhn triangulations and
a genuinely different five-tetrahedra-per-cube checkerboard triangulation:

* ``alternation=(0, 0, 0)`` selects one of the four uniform body diagonals;
* nonzero ``alternation`` entries flip the corresponding directed axis on
  alternate cube layers, giving genuinely periodic, non-uniform meshes.
* ``FiveTetrahedronSpec`` chooses the four vertices of one local parity for a
  central volume-1/3 tetrahedron and the other four vertices as volume-1/6
  corners.  The central parity is ``phase + sum(cube_origin) (mod 2)``.

The layer rule is conforming.  Across a face normal to direction ``k``, only
the normal sign is allowed to change.  The two tangential signs -- and hence
the diagonal induced on that square face -- agree on the adjacent cubes.

Run this file directly for exact combinatorial/orientation self-tests and a
translation-covariance test of the radial sampler.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence

import numpy as np


@dataclass(frozen=True)
class TriangulationSpec:
    """A globally conforming directed-axis Kuhn triangulation.

    ``base_signs[i]`` is ``+1`` or ``-1``.  If ``alternation[i] == 1``, its
    sign is multiplied by ``(-1)**origin[i]`` in the cube with lower corner
    ``origin``.  Dependence only on the matching cube coordinate is precisely
    what preserves the tangential signs on every shared face.
    """

    name: str = "uniform_ppp"
    base_signs: tuple[int, int, int] = (1, 1, 1)
    alternation: tuple[int, int, int] = (0, 0, 0)

    def __post_init__(self) -> None:
        if len(self.base_signs) != 3 or any(
            sign not in (-1, 1) for sign in self.base_signs
        ):
            raise ValueError("base_signs must contain exactly three +/-1 entries")
        if len(self.alternation) != 3 or any(
            value not in (0, 1) for value in self.alternation
        ):
            raise ValueError("alternation must contain exactly three 0/1 entries")


@dataclass(frozen=True)
class FiveTetrahedronSpec:
    """Conforming five-tetrahedra checkerboard mesh.

    On cube ``c``, the central tetrahedron contains the four local vertices of
    parity ``phase + c[0] + c[1] + c[2]`` modulo two.  Across a face normal to
    any axis, local parity changes because the normal local bit changes; the
    cube-origin parity changes by the same amount.  Thus both adjacent cubes
    induce the identical diagonal on their shared square.
    """

    name: str = "five_tet_phase0"
    phase: int = 0

    def __post_init__(self) -> None:
        if self.phase not in (0, 1):
            raise ValueError("phase must be zero or one")


UNIFORM_SPECS = (
    TriangulationSpec("uniform_ppp", (1, 1, 1), (0, 0, 0)),
    TriangulationSpec("uniform_ppm", (1, 1, -1), (0, 0, 0)),
    TriangulationSpec("uniform_pmp", (1, -1, 1), (0, 0, 0)),
    TriangulationSpec("uniform_pmm", (1, -1, -1), (0, 0, 0)),
)

PERIODIC_SPECS = (
    TriangulationSpec("layer_x", (1, 1, 1), (1, 0, 0)),
    TriangulationSpec("layer_xyz", (1, 1, 1), (1, 1, 1)),
)

FIVE_TET_SPECS = (
    FiveTetrahedronSpec("five_tet_phase0", 0),
    FiveTetrahedronSpec("five_tet_phase1", 1),
)


def permutation_sign(permutation: Sequence[int]) -> int:
    """Sign of a permutation of ``(0, 1, 2)``."""

    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(3)
        for j in range(i + 1, 3)
    )
    return -1 if inversions % 2 else 1


def cube_signs(
    origin: Sequence[int], spec: TriangulationSpec
) -> tuple[int, int, int]:
    """Directed coordinate signs in one cube."""

    if len(origin) != 3:
        raise ValueError("origin must have three coordinates")
    return tuple(
        int(spec.base_signs[i] * (-1) ** (spec.alternation[i] * int(origin[i])))
        for i in range(3)
    )


def cube_tetrahedra(
    origin: Sequence[int], signs: Sequence[int]
) -> tuple[np.ndarray, np.ndarray]:
    """Six Kuhn tetrahedra in one cube and their spatial orientations.

    The returned vertex array has shape ``(6, 4, 3)``.  If ``D`` has columns
    ``v_1-v_0, v_2-v_0, v_3-v_0``, the matching orientation entry is exactly
    ``sign(det D) = prod(signs) * sign(permutation)``.
    """

    origin_array = np.asarray(origin, dtype=np.int64)
    signs_array = np.asarray(signs, dtype=np.int64)
    if origin_array.shape != (3,):
        raise ValueError("origin must have shape (3,)")
    if signs_array.shape != (3,) or not np.all(np.isin(signs_array, (-1, 1))):
        raise ValueError("signs must have shape (3,) and entries +/-1")
    start = origin_array + (signs_array < 0).astype(np.int64)
    rows: list[np.ndarray] = []
    orientations: list[int] = []
    sign_product = int(np.prod(signs_array))
    for permutation in itertools.permutations((0, 1, 2)):
        running = start.copy()
        vertices = [running.copy()]
        for direction in permutation:
            running = running.copy()
            running[direction] += signs_array[direction]
            vertices.append(running)
        rows.append(np.asarray(vertices, dtype=np.int64))
        orientations.append(sign_product * permutation_sign(permutation))
    return np.stack(rows), np.asarray(orientations, dtype=np.int8)


def cube_five_tetrahedra(
    origin: Sequence[int], central_parity: int
) -> tuple[np.ndarray, np.ndarray]:
    """Central-plus-four-corners five-tetrahedra decomposition of one cube."""

    origin_array = np.asarray(origin, dtype=np.int64)
    if origin_array.shape != (3,):
        raise ValueError("origin must have shape (3,)")
    if central_parity not in (0, 1):
        raise ValueError("central_parity must be zero or one")
    bits = [
        np.asarray(vertex, dtype=np.int64)
        for vertex in itertools.product((0, 1), repeat=3)
    ]
    central = [vertex for vertex in bits if int(np.sum(vertex)) % 2 == central_parity]
    corners = [vertex for vertex in bits if int(np.sum(vertex)) % 2 != central_parity]
    rows = [np.stack(central) + origin_array]
    for corner in corners:
        adjacent = []
        for direction in range(3):
            neighbor = corner.copy()
            neighbor[direction] = 1 - neighbor[direction]
            adjacent.append(neighbor)
        rows.append(np.stack([corner, *adjacent]) + origin_array)
    coordinates = np.stack(rows)
    determinants = np.rint(
        np.linalg.det(
            np.moveaxis(coordinates[:, 1:] - coordinates[:, :1], 1, 2)
        )
    ).astype(np.int64)
    if sorted(np.abs(determinants).tolist()) != [1, 1, 1, 1, 2]:
        raise AssertionError("invalid five-tetrahedra cube decomposition")
    return coordinates, np.sign(determinants).astype(np.int8)


def freudenthal_tetrahedra(
    extent: int, spec: TriangulationSpec = UNIFORM_SPECS[0]
) -> tuple[np.ndarray, np.ndarray]:
    """Linear-index tetrahedra and exact spatial orientation signs."""

    if extent < 2:
        raise ValueError("extent must be at least two")
    shape = (extent,) * 3
    tetrahedra: list[np.ndarray] = []
    orientations: list[np.ndarray] = []
    for origin in np.ndindex((extent - 1,) * 3):
        coordinates, signs = cube_tetrahedra(origin, cube_signs(origin, spec))
        linear = np.ravel_multi_index(
            tuple(np.moveaxis(coordinates, -1, 0)), shape
        )
        tetrahedra.append(linear.astype(np.int64))
        orientations.append(signs)
    return np.concatenate(tetrahedra), np.concatenate(orientations)


def checkerboard_five_tetrahedra(
    extent: int, spec: FiveTetrahedronSpec = FIVE_TET_SPECS[0]
) -> tuple[np.ndarray, np.ndarray]:
    """Linear-index five-tet checkerboard complex and orientation signs."""

    if extent < 2:
        raise ValueError("extent must be at least two")
    shape = (extent,) * 3
    tetrahedra: list[np.ndarray] = []
    orientations: list[np.ndarray] = []
    for origin in np.ndindex((extent - 1,) * 3):
        central_parity = (spec.phase + sum(origin)) % 2
        coordinates, signs = cube_five_tetrahedra(origin, central_parity)
        linear = np.ravel_multi_index(
            tuple(np.moveaxis(coordinates, -1, 0)), shape
        )
        tetrahedra.append(linear.astype(np.int64))
        orientations.append(signs)
    return np.concatenate(tetrahedra), np.concatenate(orientations)


def tetrahedral_complex(
    extent: int, spec: Any
) -> tuple[np.ndarray, np.ndarray]:
    """Dispatch to a six-tet Kuhn or five-tet checkerboard complex."""

    if isinstance(spec, TriangulationSpec):
        return freudenthal_tetrahedra(extent, spec)
    if isinstance(spec, FiveTetrahedronSpec):
        return checkerboard_five_tetrahedra(extent, spec)
    raise TypeError("unsupported triangulation specification")


def unique_edges(tetrahedra: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unique unordered vertex pairs in a tetrahedral complex."""

    tetrahedra = np.asarray(tetrahedra, dtype=np.int64)
    if tetrahedra.ndim != 2 or tetrahedra.shape[1] != 4:
        raise ValueError("tetrahedra must have shape (number, 4)")
    edges: set[tuple[int, int]] = set()
    for row in tetrahedra:
        for left, right in itertools.combinations(row.tolist(), 2):
            edges.add((min(left, right), max(left, right)))
    values = np.asarray(sorted(edges), dtype=np.int64)
    return values[:, 0], values[:, 1]


def regular_value_degree(
    field: np.ndarray,
    tetrahedra: np.ndarray,
    orientations: np.ndarray,
    target: np.ndarray,
    determinant_tolerance: float = 1.0e-13,
    cone_tolerance: float = 1.0e-11,
) -> dict[str, Any]:
    """Degree of normalized-affine interpolation at one regular value.

    For a tetrahedron with field-vertex matrix ``N`` and spatial orientation
    ``sigma = sign(det D)``, a positive solution of ``N lambda = target``
    contributes ``sign(det N) * sigma`` to the map degree.  Route E's baryon
    convention is the negative of this degree.
    """

    flat = np.asarray(field, dtype=float).reshape(-1, 4)
    tetrahedra = np.asarray(tetrahedra, dtype=np.int64)
    orientations = np.asarray(orientations, dtype=np.int8)
    target = np.asarray(target, dtype=float)
    if tetrahedra.shape != (len(orientations), 4):
        raise ValueError("tetrahedra/orientations have inconsistent shapes")
    if target.shape != (4,):
        raise ValueError("target must have shape (4,)")
    degree = 0
    hits = 0
    minimum_abs_determinant = math.inf
    minimum_cone_margin = math.inf
    for vertices, orientation in zip(tetrahedra, orientations):
        matrix = flat[vertices].T
        determinant = float(np.linalg.det(matrix))
        if abs(determinant) < determinant_tolerance:
            continue
        coefficients = np.linalg.solve(matrix, target)
        margin = float(np.min(coefficients))
        if margin > cone_tolerance:
            hits += 1
            degree += int(math.copysign(1, determinant)) * int(orientation)
            minimum_abs_determinant = min(
                minimum_abs_determinant, abs(determinant)
            )
            minimum_cone_margin = min(minimum_cone_margin, margin)
    return {
        "map_degree": degree,
        "baryon_number": -degree,
        "preimage_tetrahedra": hits,
        "minimum_hit_absolute_determinant": (
            minimum_abs_determinant
            if math.isfinite(minimum_abs_determinant)
            else None
        ),
        "minimum_hit_positive_cone_margin": (
            minimum_cone_margin if math.isfinite(minimum_cone_margin) else None
        ),
    }


def triangulation_audit(extent: int, spec: Any) -> dict[str, Any]:
    """Exact volume, orientation, face-conformity, and Euler audit."""

    tetrahedra, orientations = tetrahedral_complex(extent, spec)
    coordinates = np.stack(np.unravel_index(tetrahedra, (extent,) * 3), axis=-1)
    determinants = np.rint(
        np.linalg.det(
            np.moveaxis(coordinates[:, 1:] - coordinates[:, :1], 1, 2)
        )
    ).astype(np.int64)
    orientation_ok = bool(
        np.array_equal(np.sign(determinants), orientations.astype(np.int64))
    )
    absolute_determinants = np.abs(determinants)
    cubes = (extent - 1) ** 3
    volume_sum_ok = bool(np.sum(absolute_determinants) == 6 * cubes)
    if isinstance(spec, TriangulationSpec):
        volume_distribution_ok = bool(np.all(absolute_determinants == 1))
        tetrahedra_per_cube = 6
    else:
        volume_distribution_ok = bool(
            np.count_nonzero(absolute_determinants == 2) == cubes
            and np.count_nonzero(absolute_determinants == 1) == 4 * cubes
            and np.all(np.isin(absolute_determinants, (1, 2)))
        )
        tetrahedra_per_cube = 5

    sorted_tetrahedra = [tuple(sorted(row.tolist())) for row in tetrahedra]
    no_duplicate_tetrahedra = len(set(sorted_tetrahedra)) == len(tetrahedra)

    face_counts: dict[tuple[int, int, int], int] = {}
    for tetrahedron in tetrahedra:
        for face in itertools.combinations(tetrahedron.tolist(), 3):
            key = tuple(sorted(face))
            face_counts[key] = face_counts.get(key, 0) + 1
    bad_boundary_faces = 0
    bad_interior_faces = 0
    boundary_faces = 0
    interior_faces = 0
    for face, count in face_counts.items():
        face_coordinates = np.stack(
            np.unravel_index(np.asarray(face), (extent,) * 3), axis=-1
        )
        on_boundary = any(
            np.all(face_coordinates[:, direction] == endpoint)
            for direction in range(3)
            for endpoint in (0, extent - 1)
        )
        if on_boundary:
            boundary_faces += 1
            bad_boundary_faces += int(count != 1)
        else:
            interior_faces += 1
            bad_interior_faces += int(count != 2)

    left, right = unique_edges(tetrahedra)
    vertices = extent**3
    faces = len(face_counts)
    euler_characteristic = vertices - len(left) + faces - len(tetrahedra)
    expected_tetrahedra = tetrahedra_per_cube * cubes
    expected_boundary_faces = 12 * (extent - 1) ** 2
    valid = all(
        (
            len(tetrahedra) == expected_tetrahedra,
            orientation_ok,
            volume_sum_ok,
            volume_distribution_ok,
            no_duplicate_tetrahedra,
            bad_boundary_faces == 0,
            bad_interior_faces == 0,
            boundary_faces == expected_boundary_faces,
            euler_characteristic == 1,
        )
    )
    return {
        "spec": asdict(spec),
        "extent": extent,
        "tetrahedra": len(tetrahedra),
        "expected_tetrahedra": expected_tetrahedra,
        "tetrahedra_per_cube": tetrahedra_per_cube,
        "unique_edges": len(left),
        "unique_faces": faces,
        "boundary_faces": boundary_faces,
        "expected_boundary_faces": expected_boundary_faces,
        "interior_faces": interior_faces,
        "bad_boundary_face_multiplicities": bad_boundary_faces,
        "bad_interior_face_multiplicities": bad_interior_faces,
        "orientation_formula_exact": orientation_ok,
        "six_times_total_tetrahedron_volume_per_cube_equals_six": (
            volume_sum_ok
        ),
        "tetrahedron_volume_distribution_exact": volume_distribution_ok,
        "no_duplicate_tetrahedra": no_duplicate_tetrahedra,
        "euler_characteristic": euler_characteristic,
        "valid_conforming_triangulation": valid,
    }


def translated_radial_background(
    solution: Any,
    spatial_extent: int,
    physical_length: float,
    translation: Sequence[float] = (0.0, 0.0, 0.0),
    translation_units: str = "physical",
    radial_epsilon: float = 1.0e-8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample a hedgehog centered at ``translation``.

    ``translation_units='lattice'`` interprets the three entries in units of
    the grid spacing.  The sampler does not alter the boundary; callers should
    apply the same exact-vacuum boundary operation as their action driver.
    """

    if spatial_extent < 2:
        raise ValueError("spatial_extent must be at least two")
    translation_array = np.asarray(translation, dtype=float)
    if translation_array.shape != (3,):
        raise ValueError("translation must have shape (3,)")
    axis = np.linspace(
        -physical_length / 2.0, physical_length / 2.0, spatial_extent
    )
    spacing = float(axis[1] - axis[0])
    if translation_units == "lattice":
        translation_array = spacing * translation_array
    elif translation_units != "physical":
        raise ValueError("translation_units must be 'physical' or 'lattice'")
    mesh = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    relative = mesh - translation_array
    radius = np.linalg.norm(relative, axis=-1)
    flat_radius = radius.ravel()
    profile, _, breathing, _ = solution.sol(
        np.maximum(flat_radius, radial_epsilon)
    )
    profile = np.asarray(profile).reshape(radius.shape)
    breathing = np.asarray(breathing).reshape(radius.shape)
    field = np.zeros((*radius.shape, 4), dtype=float)
    field[..., 0] = np.cos(profile)
    nonzero = radius > 1.0e-14
    sine_over_radius = np.zeros_like(radius)
    sine_over_radius[nonzero] = np.sin(profile[nonzero]) / radius[nonzero]
    field[..., 1:] = sine_over_radius[..., None] * relative
    field[~nonzero, 0] = -1.0
    field[~nonzero, 1:] = 0.0
    breathing[~nonzero] = float(
        np.asarray(solution.sol(radial_epsilon))[2]
    )
    return field, breathing, axis, radius


def run_self_tests() -> dict[str, Any]:
    """Run exact structural tests used before AP-E9 integration."""

    specs = (*UNIFORM_SPECS, *PERIODIC_SPECS, *FIVE_TET_SPECS)
    audits = [triangulation_audit(5, spec) for spec in specs]
    all_audits_pass = all(row["valid_conforming_triangulation"] for row in audits)

    uniform_complexes = []
    for spec in UNIFORM_SPECS:
        tetrahedra, _ = freudenthal_tetrahedra(3, spec)
        uniform_complexes.append(
            frozenset(tuple(sorted(row.tolist())) for row in tetrahedra)
        )
    four_uniform_body_diagonals_distinct = len(set(uniform_complexes)) == 4

    five_tet_complexes = []
    for spec in FIVE_TET_SPECS:
        tetrahedra, _ = checkerboard_five_tetrahedra(4, spec)
        five_tet_complexes.append(
            frozenset(tuple(sorted(row.tolist())) for row in tetrahedra)
        )
    five_tet_phases_distinct = len(set(five_tet_complexes)) == 2

    sign_reversal_checks = []
    for spec in (*UNIFORM_SPECS, *PERIODIC_SPECS):
        reversed_spec = TriangulationSpec(
            f"{spec.name}_reversed",
            tuple(-entry for entry in spec.base_signs),
            spec.alternation,
        )
        original, _ = freudenthal_tetrahedra(4, spec)
        reversed_rows, _ = freudenthal_tetrahedra(4, reversed_spec)
        original_set = frozenset(tuple(sorted(row.tolist())) for row in original)
        reversed_set = frozenset(
            tuple(sorted(row.tolist())) for row in reversed_rows
        )
        sign_reversal_checks.append(original_set == reversed_set)

    class MockSolution:
        @staticmethod
        def sol(radius: Any) -> np.ndarray:
            radius_array = np.asarray(radius, dtype=float)
            profile = math.pi * np.exp(-(radius_array**2))
            derivative = -2.0 * radius_array * profile
            breathing = 1.0 + 0.2 * np.exp(-(radius_array**2))
            breathing_derivative = -0.4 * radius_array * np.exp(
                -(radius_array**2)
            )
            return np.asarray(
                [profile, derivative, breathing, breathing_derivative]
            )

    centered = translated_radial_background(MockSolution(), 7, 6.0)
    spacing = float(centered[2][1] - centered[2][0])
    shifted_physical = translated_radial_background(
        MockSolution(), 7, 6.0, (spacing, 0.0, 0.0), "physical"
    )
    shifted_lattice = translated_radial_background(
        MockSolution(), 7, 6.0, (1.0, 0.0, 0.0), "lattice"
    )
    translation_covariant = bool(
        np.allclose(shifted_physical[0][1:], centered[0][:-1], atol=2.0e-14)
        and np.allclose(
            shifted_physical[1][1:], centered[1][:-1], atol=2.0e-14
        )
    )
    unit_conventions_agree = bool(
        np.allclose(shifted_physical[0], shifted_lattice[0], atol=2.0e-14)
        and np.allclose(
            shifted_physical[1], shifted_lattice[1], atol=2.0e-14
        )
    )
    maximum_norm_residual = float(
        np.max(np.abs(np.linalg.norm(shifted_physical[0], axis=-1) - 1.0))
    )

    extent = 7
    axis = np.linspace(-4.0, 4.0, extent)
    coordinates = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    )
    radius = np.linalg.norm(coordinates, axis=-1)
    profile = math.pi * np.clip(1.0 - radius / 2.5, 0.0, 1.0)
    hedgehog = np.zeros((*radius.shape, 4))
    hedgehog[..., 0] = np.cos(profile)
    nonzero = radius > 1.0e-14
    hedgehog[..., 1:][nonzero] = (
        np.sin(profile[nonzero]) / radius[nonzero]
    )[:, None] * coordinates[nonzero]
    hedgehog[..., 0][~nonzero] = -1.0
    targets = (
        np.asarray([0.31, 0.41, -0.17, 0.839702328]),
        np.asarray([-0.23, 0.51, 0.37, 0.744916103]),
    )
    targets = tuple(target / np.linalg.norm(target) for target in targets)
    degree_rows = []
    for spec in specs:
        tetrahedra, orientations = tetrahedral_complex(extent, spec)
        baryon_numbers = [
            regular_value_degree(
                hedgehog, tetrahedra, orientations, target
            )["baryon_number"]
            for target in targets
        ]
        degree_rows.append(
            {"spec": spec.name, "baryon_numbers": baryon_numbers}
        )
    orientation_aware_degree_consistent = all(
        row["baryon_numbers"] == [1, 1] for row in degree_rows
    )
    passed = all(
        (
            all_audits_pass,
            four_uniform_body_diagonals_distinct,
            five_tet_phases_distinct,
            all(sign_reversal_checks),
            translation_covariant,
            unit_conventions_agree,
            maximum_norm_residual < 1.0e-13,
            orientation_aware_degree_consistent,
        )
    )
    return {
        "pass": passed,
        "triangulation_audits": audits,
        "four_uniform_body_diagonals_distinct": (
            four_uniform_body_diagonals_distinct
        ),
        "five_tet_checkerboard_phases_distinct": five_tet_phases_distinct,
        "global_sign_reversal_same_unoriented_complex": all(
            sign_reversal_checks
        ),
        "translated_sampler_one_site_covariance": translation_covariant,
        "translated_sampler_unit_conventions_agree": unit_conventions_agree,
        "translated_sampler_maximum_S3_norm_residual": maximum_norm_residual,
        "orientation_aware_degree_B1_consistent": (
            orientation_aware_degree_consistent
        ),
        "orientation_aware_degree_rows": degree_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if not arguments.self_test:
        parser.error("pass --self-test")
    result = run_self_tests()
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
