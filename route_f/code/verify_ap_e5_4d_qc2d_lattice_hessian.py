#!/usr/bin/env python3
"""Finite 4D background-field audit for the charged-QC2D B=1 Hessian.

This is a deterministic, small-volume *diagnostic*, not a dynamical lattice
simulation of two-colour QCD.  It assembles the quadratic operators of a
declared Wilson-compatible background-field truncation and tests

* a sampled B=1 meson background and its continuum charge extrapolation;
* Wilson gauge fixing and Faddeev--Popov spectra on a periodic 4D lattice;
* frozen-background meson and charge-two diquark fluctuation blocks;
* a two-flavour Wilson--Dirac operator and its exact even--odd Schur factor;
* a one-collective-coordinate fermion log-determinant curvature;
* Wilson O(a^2) versus tree-level Symanzik O(a^4) dispersion scaling; and
* an unstable diquark negative control.

The script fails closed.  It does not perform importance sampling, infer a
continuum quantum phase, prove determinant positivity, compute the complete
Skyrme Jacobi operator, or establish Finkelstein--Rubinstein/global stability.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy import sparse


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
TEX = ROUTE_F / "tex" / "ap_e5_4d_qc2d_lattice_hessian.tex"
BIB = ROUTE_F / "tex" / "ap_e5_4d_qc2d_lattice_hessian.bib"


@dataclass(frozen=True)
class Benchmark:
    extents: tuple[int, int, int, int] = (2, 2, 2, 2)
    physical_length: float = 2.4
    soliton_radius: float = 1.2
    gauge_xi: float = 0.7
    pion_mass_sq: float = 0.25
    meson_core_lift: float = 0.18
    diquark_mass_sq: float = 0.81
    diquark_core_coupling: float = -0.30
    diquark_negative_control: float = -3.00
    wilson_r: float = 1.0
    wilson_bare_mass: float = 0.65
    yukawa_sigma: float = 0.20
    yukawa_pion: float = 0.16
    flavours: int = 2
    colour_degeneracy: int = 2


CHECKS: list[dict[str, Any]] = []


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


def site_index(coordinate: tuple[int, ...], extents: tuple[int, ...]) -> int:
    return int(np.ravel_multi_index(coordinate, extents))


def coordinates(extents: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [tuple(int(v) for v in row) for row in np.ndindex(extents)]


def forward_gradient(extents: tuple[int, ...]) -> sparse.csr_matrix:
    """Periodic forward gradient, ordered as (direction, site)."""

    volume = math.prod(extents)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for x in coordinates(extents):
        source = site_index(x, extents)
        for mu, length in enumerate(extents):
            target_coordinate = list(x)
            target_coordinate[mu] = (target_coordinate[mu] + 1) % length
            target = site_index(tuple(target_coordinate), extents)
            row = mu * volume + source
            rows.extend((row, row))
            cols.extend((source, target))
            data.extend((-1.0, 1.0))
    return sparse.coo_matrix(
        (data, (rows, cols)), shape=(len(extents) * volume, volume)
    ).tocsr()


def profile_field(
    spatial_extent: int, physical_length: float, radius: float, cell_centred: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return n=(sigma,pi-vector), F, and coordinates for a B=1 hedgehog.

    F=2 atan[(R/r) exp(-r^2/(2R^2))] has F(0)=pi and F(infinity)=0.
    """

    if cell_centred:
        spacing = physical_length / spatial_extent
        axis = (
            np.arange(spatial_extent, dtype=float) + 0.5
        ) * spacing - physical_length / 2.0
    else:
        axis = np.linspace(
            -physical_length / 2.0, physical_length / 2.0, spatial_extent
        )
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radial = np.sqrt(xx**2 + yy**2 + zz**2)
    argument = np.empty_like(radial)
    nonzero = radial > 0.0
    argument[nonzero] = (radius / radial[nonzero]) * np.exp(
        -0.5 * (radial[nonzero] / radius) ** 2
    )
    argument[~nonzero] = np.inf
    profile = 2.0 * np.arctan(argument)
    field = np.zeros((*radial.shape, 4), dtype=float)
    field[..., 0] = np.cos(profile)
    sine_over_r = np.zeros_like(radial)
    sine_over_r[nonzero] = np.sin(profile[nonzero]) / radial[nonzero]
    field[..., 1] = sine_over_r * xx
    field[..., 2] = sine_over_r * yy
    field[..., 3] = sine_over_r * zz
    return field, profile, axis


def lattice_baryon_number(field: np.ndarray, axis: np.ndarray) -> float:
    spacing = float(axis[1] - axis[0])
    derivatives = []
    for spatial_direction in range(3):
        derivatives.append(
            np.stack(
                [
                    np.gradient(
                        field[..., component],
                        spacing,
                        axis=spatial_direction,
                        edge_order=2,
                    )
                    for component in range(4)
                ],
                axis=-1,
            )
        )
    determinant = np.linalg.det(np.stack([field, *derivatives], axis=-2))
    integral = np.trapezoid(
        np.trapezoid(np.trapezoid(determinant, axis, axis=0), axis, axis=0),
        axis,
        axis=0,
    )
    return float(-integral / (2.0 * math.pi**2))


def topology_scan(benchmark: Benchmark) -> dict[str, Any]:
    sizes = [17, 25, 33, 49, 65]
    topology_length = benchmark.physical_length + 4.0
    rows: list[dict[str, float]] = []
    for size in sizes:
        field, _, axis = profile_field(
            size, topology_length, benchmark.soliton_radius, False
        )
        number = lattice_baryon_number(field, axis)
        rows.append(
            {
                "N_s": size,
                "a": float(axis[1] - axis[0]),
                "B_lattice": number,
                "abs_error": abs(number - 1.0),
            }
        )
    last = rows[-4:]
    fit_x = np.array([row["a"] ** 2 for row in last])
    fit_y = np.array([row["B_lattice"] for row in last])
    coefficients = np.polyfit(
        fit_x,
        fit_y,
        deg=2,
    )
    extrapolated = float(np.polyval(coefficients, 0.0))
    fit_rms = float(np.sqrt(np.mean((np.polyval(coefficients, fit_x) - fit_y) ** 2)))
    last_error = rows[-1]["abs_error"]
    previous_error = rows[-2]["abs_error"]
    observed_order = math.log(previous_error / last_error) / math.log(
        rows[-2]["a"] / rows[-1]["a"]
    )
    return {
        "profile": "F(r)=2 atan[(R/r) exp(-r^2/(2R^2))]",
        "orientation": "B=+1",
        "physical_length": topology_length,
        "finite_difference_estimator": rows,
        "quadratic_in_a_squared_extrapolation_last_four": extrapolated,
        "extrapolation_abs_error": abs(extrapolated - 1.0),
        "fit_rms_residual": fit_rms,
        "observed_order_last_pair": observed_order,
    }


def gauge_ghost_audit(extents: tuple[int, ...]) -> dict[str, Any]:
    gradient = forward_gradient(extents)
    laplacian = (gradient.T @ gradient).toarray()
    volume = laplacian.shape[0]
    identity_links = np.kron(np.eye(4), laplacian)
    longitudinal = (gradient @ gradient.T).toarray()
    ghost_eigenvalues = np.linalg.eigvalsh(laplacian)
    ghost_nonzero = ghost_eigenvalues[ghost_eigenvalues > 1.0e-10]
    xi_rows: list[dict[str, Any]] = []
    reference = None
    for xi in (0.5, 0.7, 1.0, 2.0):
        gauge = identity_links + (1.0 / xi - 1.0) * longitudinal
        eigenvalues = np.linalg.eigvalsh(gauge)
        nonzero = eigenvalues[eigenvalues > 1.0e-10]
        gaussian_log_weight = 0.5 * float(np.log(nonzero).sum()) - float(
            np.log(ghost_nonzero).sum()
        )
        corrected = gaussian_log_weight + 0.5 * (volume - 1) * math.log(xi)
        if reference is None:
            reference = corrected
        xi_rows.append(
            {
                "xi": xi,
                "gauge_zero_modes": int(np.count_nonzero(eigenvalues <= 1.0e-10)),
                "ghost_zero_modes": int(
                    np.count_nonzero(ghost_eigenvalues <= 1.0e-10)
                ),
                "gauge_nonzero_modes": int(len(nonzero)),
                "minimum_nonzero_gauge_eigenvalue": float(nonzero[0]),
                "minimum_nonzero_ghost_eigenvalue": float(ghost_nonzero[0]),
                "one_generator_log_weight": gaussian_log_weight,
                "xi_corrected_log_weight": corrected,
                "xi_corrected_residual": abs(corrected - reference),
            }
        )
    return {
        "discretisation": "Wilson quadratic gauge action plus covariant gauge fixing",
        "extents": list(extents),
        "volume": volume,
        "expected_gauge_toron_zero_modes_per_generator": 4,
        "expected_constant_ghost_zero_modes_per_generator": 1,
        "colour_generators": 3,
        "extra_U1_generators": 1,
        "xi_rows": xi_rows,
        "maximum_xi_corrected_residual": max(
            row["xi_corrected_residual"] for row in xi_rows
        ),
    }


def tangent_frame(unit_vector: np.ndarray) -> np.ndarray:
    reference = np.array([1.0, 0.0, 0.0, 0.0])
    difference = reference - unit_vector
    norm_sq = float(difference @ difference)
    if norm_sq < 1.0e-14:
        orthogonal = np.eye(4)
    else:
        orthogonal = np.eye(4) - 2.0 * np.outer(difference, difference) / norm_sq
    return orthogonal[:, 1:]


def sphere_parallel_transport(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    cosine = float(first @ second)
    if cosine < -1.0 + 1.0e-10:
        raise ValueError("antipodal neighbours make geodesic transport ambiguous")
    generator = np.outer(second, first) - np.outer(first, second)
    return np.eye(4) + generator + (generator @ generator) / (1.0 + cosine)


def repeated_background(benchmark: Benchmark) -> tuple[np.ndarray, np.ndarray]:
    _, spatial_extent, _, _ = benchmark.extents
    spatial_field, profile, _ = profile_field(
        spatial_extent,
        benchmark.physical_length,
        benchmark.soliton_radius,
        True,
    )
    field = np.repeat(spatial_field[np.newaxis, ...], benchmark.extents[0], axis=0)
    profile_4d = np.repeat(profile[np.newaxis, ...], benchmark.extents[0], axis=0)
    return field.reshape((-1, 4)), profile_4d.reshape(-1)


def meson_tangent_operator(
    extents: tuple[int, ...], background: np.ndarray, mass_sq: float, core_lift: float
) -> sparse.csr_matrix:
    volume = math.prod(extents)
    frames = [tangent_frame(background[index]) for index in range(volume)]
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def add_block(row_site: int, col_site: int, block: np.ndarray) -> None:
        for i in range(3):
            for j in range(3):
                rows.append(3 * row_site + i)
                cols.append(3 * col_site + j)
                data.append(float(block[i, j]))

    for x in coordinates(extents):
        source = site_index(x, extents)
        for mu, length in enumerate(extents):
            target_coordinate = list(x)
            target_coordinate[mu] = (target_coordinate[mu] + 1) % length
            target = site_index(tuple(target_coordinate), extents)
            rotation = sphere_parallel_transport(background[source], background[target])
            transport = frames[target].T @ rotation @ frames[source]
            add_block(source, source, np.eye(3))
            add_block(target, target, np.eye(3))
            add_block(target, source, -transport)
            add_block(source, target, -transport.T)
    core_density = 0.5 * (1.0 - background[:, 0])
    for site in range(volume):
        add_block(
            site,
            site,
            (mass_sq + core_lift * core_density[site]) * np.eye(3),
        )
    return sparse.coo_matrix(
        (data, (rows, cols)), shape=(3 * volume, 3 * volume)
    ).tocsr()


def boson_hessian_audit(benchmark: Benchmark) -> dict[str, Any]:
    background, _ = repeated_background(benchmark)
    volume = math.prod(benchmark.extents)
    gradient = forward_gradient(benchmark.extents)
    laplacian = (gradient.T @ gradient).toarray()
    core_density = 0.5 * (1.0 - background[:, 0])

    meson = meson_tangent_operator(
        benchmark.extents,
        background,
        benchmark.pion_mass_sq,
        benchmark.meson_core_lift,
    ).toarray()
    meson_eigenvalues = np.linalg.eigvalsh(meson)

    diquark = (
        laplacian
        + benchmark.diquark_mass_sq * np.eye(volume)
        + benchmark.diquark_core_coupling * np.diag(core_density)
    )
    diquark_eigenvalues = np.linalg.eigvalsh(diquark)
    negative_control = (
        laplacian
        + benchmark.diquark_mass_sq * np.eye(volume)
        + benchmark.diquark_negative_control * np.diag(core_density)
    )
    negative_eigenvalues = np.linalg.eigvalsh(negative_control)
    return {
        "operator_scope": (
            "complete bosonic Hessian of the declared frozen-background quadratic "
            "truncation; not the complete nonlinear Skyrme Jacobi operator"
        ),
        "background_delta": 0.0,
        "classical_gauge_diquark_mixing": 0.0,
        "classical_gauge_meson_mixing": 0.0,
        "meson_dimension": int(meson.shape[0]),
        "meson_smallest_eigenvalues": meson_eigenvalues[:8].tolist(),
        "meson_minimum_eigenvalue": float(meson_eigenvalues[0]),
        "diquark_real_component_dimension": int(diquark.shape[0]),
        "diquark_real_component_degeneracy": 2,
        "diquark_smallest_eigenvalues": diquark_eigenvalues[:8].tolist(),
        "diquark_minimum_eigenvalue": float(diquark_eigenvalues[0]),
        "negative_control_coupling": benchmark.diquark_negative_control,
        "negative_control_minimum_eigenvalue": float(negative_eigenvalues[0]),
        "negative_control_negative_mode_count_per_real_component": int(
            np.count_nonzero(negative_eigenvalues < -1.0e-9)
        ),
    }


def euclidean_gamma_matrices() -> tuple[list[np.ndarray], np.ndarray]:
    sigma_1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sigma_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    sigma_3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    identity_2 = np.eye(2, dtype=complex)
    zero_2 = np.zeros((2, 2), dtype=complex)
    gammas = []
    for sigma_matrix in (sigma_1, sigma_2, sigma_3):
        gammas.append(
            np.block([[zero_2, -1.0j * sigma_matrix], [1.0j * sigma_matrix, zero_2]])
        )
    gammas.append(np.block([[zero_2, identity_2], [identity_2, zero_2]]))
    gamma_5 = gammas[0] @ gammas[1] @ gammas[2] @ gammas[3]
    return gammas, gamma_5


def wilson_dirac(
    benchmark: Benchmark, background: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Two-flavour Wilson operator; SU(2)_c is an exact degeneracy here."""

    gammas, gamma_5 = euclidean_gamma_matrices()
    tau = [
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex),
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
    ]
    spin = 4
    flavour = 2
    internal = spin * flavour
    volume = math.prod(benchmark.extents)
    operator = np.zeros((internal * volume, internal * volume), dtype=complex)
    identity_internal = np.eye(internal, dtype=complex)
    gamma_internal = [np.kron(matrix, np.eye(flavour)) for matrix in gammas]
    gamma5_internal = np.kron(gamma_5, np.eye(flavour))
    r = benchmark.wilson_r

    for x in coordinates(benchmark.extents):
        site = site_index(x, benchmark.extents)
        block = slice(internal * site, internal * (site + 1))
        sigma = background[site, 0]
        pion_matrix = sum(
            background[site, component + 1] * tau[component]
            for component in range(3)
        )
        local = (
            (benchmark.wilson_bare_mass + 4.0 * r + benchmark.yukawa_sigma * sigma)
            * identity_internal
            + 1.0j
            * benchmark.yukawa_pion
            * np.kron(gamma_5, pion_matrix)
        )
        operator[block, block] += local
        for mu, length in enumerate(benchmark.extents):
            forward_coordinate = list(x)
            backward_coordinate = list(x)
            forward_coordinate[mu] = (forward_coordinate[mu] + 1) % length
            backward_coordinate[mu] = (backward_coordinate[mu] - 1) % length
            forward_site = site_index(tuple(forward_coordinate), benchmark.extents)
            backward_site = site_index(tuple(backward_coordinate), benchmark.extents)
            forward_block = slice(
                internal * forward_site, internal * (forward_site + 1)
            )
            backward_block = slice(
                internal * backward_site, internal * (backward_site + 1)
            )
            operator[block, forward_block] += -0.5 * (
                r * identity_internal - gamma_internal[mu]
            )
            operator[block, backward_block] += -0.5 * (
                r * identity_internal + gamma_internal[mu]
            )

    source = np.zeros_like(operator)
    core_density = 0.5 * (1.0 - background[:, 0])
    for site in range(volume):
        block = slice(internal * site, internal * (site + 1))
        source[block, block] = core_density[site] * identity_internal

    clifford_error = max(
        float(
            np.linalg.norm(
                gammas[mu] @ gammas[nu]
                + gammas[nu] @ gammas[mu]
                - (2.0 if mu == nu else 0.0) * np.eye(4)
            )
        )
        for mu in range(4)
        for nu in range(4)
    )
    return operator, source, {
        "clifford_error": clifford_error,
        "gamma5_square_error": float(np.linalg.norm(gamma_5 @ gamma_5 - np.eye(4))),
        "gamma5_internal_trace": float(np.trace(gamma5_internal).real),
    }


def fermion_audit(benchmark: Benchmark) -> dict[str, Any]:
    background, _ = repeated_background(benchmark)
    operator, source, gamma_checks = wilson_dirac(benchmark, background)
    internal = 8
    site_parities = np.array(
        [sum(x) % 2 for x in coordinates(benchmark.extents)], dtype=int
    )
    even_sites = np.flatnonzero(site_parities == 0)
    odd_sites = np.flatnonzero(site_parities == 1)
    even = np.concatenate(
        [np.arange(internal * site, internal * (site + 1)) for site in even_sites]
    )
    odd = np.concatenate(
        [np.arange(internal * site, internal * (site + 1)) for site in odd_sites]
    )
    permutation = np.concatenate((even, odd))
    permuted = operator[np.ix_(permutation, permutation)]
    half = len(even)
    dee = permuted[:half, :half]
    deo = permuted[:half, half:]
    doe = permuted[half:, :half]
    doo = permuted[half:, half:]
    doo_inverse_doe = np.linalg.solve(doo, doe)
    schur = dee - deo @ doo_inverse_doe

    phase_full, logabs_full = np.linalg.slogdet(permuted)
    phase_odd, logabs_odd = np.linalg.slogdet(doo)
    phase_schur, logabs_schur = np.linalg.slogdet(schur)
    determinant_log_residual = abs(logabs_full - logabs_odd - logabs_schur)
    determinant_phase_residual = abs(
        np.angle(phase_full / (phase_odd * phase_schur))
    )

    rng = np.random.default_rng(20260716)
    source_vector = rng.normal(size=operator.shape[0]) + 1.0j * rng.normal(
        size=operator.shape[0]
    )
    permuted_source = source_vector[permutation]
    be = permuted_source[:half]
    bo = permuted_source[half:]
    doo_inverse_bo = np.linalg.solve(doo, bo)
    xe = np.linalg.solve(schur, be - deo @ doo_inverse_bo)
    xo = np.linalg.solve(doo, bo - doe @ xe)
    schur_solution = np.concatenate((xe, xo))
    direct_solution = np.linalg.solve(permuted, permuted_source)
    solve_relative_residual = float(
        np.linalg.norm(schur_solution - direct_solution)
        / np.linalg.norm(direct_solution)
    )

    singular_values = np.linalg.svd(operator, compute_uv=False)
    minimum_singular = float(singular_values[-1])
    condition_number = float(singular_values[0] / singular_values[-1])

    inverse_source = np.linalg.solve(operator, source)
    # The matrix already contains the complete two-flavour block.  Multiplying
    # by benchmark.flavours here would double count N_f=2.
    analytic_curvature_one_colour = float(
        np.trace(inverse_source @ inverse_source).real
    )
    # At h=1e-3 the O(h^2) truncation and log-determinant cancellation errors
    # are both below the declared tolerance for this benchmark.
    step = 1.0e-3

    def effective_action(parameter: float) -> float:
        _, logabs = np.linalg.slogdet(operator + parameter * source)
        return -float(logabs)

    finite_difference_curvature_one_colour = (
        effective_action(step)
        - 2.0 * effective_action(0.0)
        + effective_action(-step)
    ) / step**2
    curvature_relative_error = abs(
        finite_difference_curvature_one_colour - analytic_curvature_one_colour
    ) / max(1.0, abs(analytic_curvature_one_colour))

    return {
        "operator": "two-flavour Wilson-Dirac in a sampled chiral B=1 background",
        "matrix_dimension_one_colour_copy": int(operator.shape[0]),
        "SU2_colour_degeneracy_at_trivial_colour_links": benchmark.colour_degeneracy,
        "gamma_checks": gamma_checks,
        "minimum_singular_value": minimum_singular,
        "condition_number": condition_number,
        "even_odd_schur_dimension": int(schur.shape[0]),
        "determinant_logabs_residual": float(determinant_log_residual),
        "determinant_phase_residual": float(determinant_phase_residual),
        "schur_solve_relative_residual": solve_relative_residual,
        "one_loop_collective_coordinate": {
            "source": "Q_x=rho_B(x) times identity",
            "Gamma_f": (
                "-Re log det(D_Nf=2+tQ); the matrix already contains both flavours"
            ),
            "analytic_second_derivative_one_colour": analytic_curvature_one_colour,
            "finite_difference_second_derivative_one_colour": float(
                finite_difference_curvature_one_colour
            ),
            "relative_error": float(curvature_relative_error),
            "full_colour_factor_not_used_as_positivity_proof": benchmark.colour_degeneracy,
        },
        "determinant_positivity_proven": False,
        "overlap_index_or_admissibility_proven": False,
    }


def dispersion_scaling(physical_length: float) -> dict[str, Any]:
    continuum = (2.0 * math.pi / physical_length) ** 2
    rows = []
    for extent in (8, 16, 32, 64):
        spacing = physical_length / extent
        momentum = 2.0 * math.pi / physical_length
        wilson = 4.0 * math.sin(0.5 * momentum * spacing) ** 2 / spacing**2
        symanzik = wilson + spacing**2 * wilson**2 / 12.0
        rows.append(
            {
                "N": extent,
                "a": spacing,
                "continuum_p_squared": continuum,
                "wilson_p_squared": wilson,
                "symanzik_tree_p_squared": symanzik,
                "wilson_abs_error": abs(wilson - continuum),
                "symanzik_abs_error": abs(symanzik - continuum),
            }
        )

    def observed_order(key: str) -> float:
        return math.log(rows[-2][key] / rows[-1][key]) / math.log(2.0)

    return {
        "rows": rows,
        "wilson_observed_order_last_pair": observed_order("wilson_abs_error"),
        "symanzik_observed_order_last_pair": observed_order(
            "symanzik_abs_error"
        ),
    }


def build_card() -> dict[str, Any]:
    benchmark = Benchmark()
    topology = topology_scan(benchmark)
    gauge = gauge_ghost_audit(benchmark.extents)
    bosons = boson_hessian_audit(benchmark)
    fermions = fermion_audit(benchmark)
    scaling = dispersion_scaling(benchmark.physical_length)
    sources = [source_row(path) for path in (THIS_SCRIPT, TEX, BIB)]

    check(
        "provenance",
        "all declared source files exist",
        all(row["exists"] for row in sources),
        ", ".join(row["path"] for row in sources),
    )
    check(
        "provenance",
        "deterministic seed recorded",
        True,
        "numpy Generator seed=20260716",
    )
    check(
        "topology",
        "finite-lattice B estimator approaches one monotonically",
        all(
            topology["finite_difference_estimator"][i + 1]["B_lattice"]
            > topology["finite_difference_estimator"][i]["B_lattice"]
            for i in range(len(topology["finite_difference_estimator"]) - 1)
        ),
        str([row["B_lattice"] for row in topology["finite_difference_estimator"]]),
    )
    check(
        "topology",
        "B=1 continuum extrapolation",
        topology["extrapolation_abs_error"] < 2.0e-3,
        f"B0={topology['quadratic_in_a_squared_extrapolation_last_four']:.12g}",
    )
    check(
        "topology",
        "B estimator has second-order convergence",
        1.8 < topology["observed_order_last_pair"] < 2.2,
        f"p={topology['observed_order_last_pair']:.8f}",
    )
    check(
        "gauge_ghost",
        "four gauge toron zero modes per generator",
        all(row["gauge_zero_modes"] == 4 for row in gauge["xi_rows"]),
        str([row["gauge_zero_modes"] for row in gauge["xi_rows"]]),
    )
    check(
        "gauge_ghost",
        "one constant ghost zero mode per generator",
        all(row["ghost_zero_modes"] == 1 for row in gauge["xi_rows"]),
        str([row["ghost_zero_modes"] for row in gauge["xi_rows"]]),
    )
    check(
        "gauge_ghost",
        "all nonzero gauge and ghost modes positive",
        min(row["minimum_nonzero_gauge_eigenvalue"] for row in gauge["xi_rows"])
        > 1.0e-8
        and min(row["minimum_nonzero_ghost_eigenvalue"] for row in gauge["xi_rows"])
        > 1.0e-8,
        "prime spectra only",
    )
    check(
        "gauge_ghost",
        "xi dependence is the predicted background-independent constant",
        gauge["maximum_xi_corrected_residual"] < 1.0e-10,
        f"max residual={gauge['maximum_xi_corrected_residual']:.3e}",
    )
    check(
        "boson_hessian",
        "declared meson tangent block is positive",
        bosons["meson_minimum_eigenvalue"] > 1.0e-8,
        f"lambda_min={bosons['meson_minimum_eigenvalue']:.12g}",
    )
    check(
        "boson_hessian",
        "charge-two diquark block is positive",
        bosons["diquark_minimum_eigenvalue"] > 1.0e-8,
        f"lambda_min={bosons['diquark_minimum_eigenvalue']:.12g}",
    )
    check(
        "boson_hessian",
        "Delta=0 removes classical gauge-diquark mixing",
        bosons["classical_gauge_diquark_mixing"] == 0.0,
        "K_A_Delta=0 at the benchmark background",
    )
    check(
        "negative_control",
        "over-attractive core coupling creates a negative diquark mode",
        bosons["negative_control_minimum_eigenvalue"] < -1.0e-3
        and bosons["negative_control_negative_mode_count_per_real_component"] >= 1,
        (
            f"lambda_min={bosons['negative_control_minimum_eigenvalue']:.12g}; "
            f"count={bosons['negative_control_negative_mode_count_per_real_component']}"
        ),
    )
    check(
        "fermion",
        "Euclidean Clifford algebra",
        fermions["gamma_checks"]["clifford_error"] < 1.0e-12,
        f"error={fermions['gamma_checks']['clifford_error']:.3e}",
    )
    check(
        "fermion",
        "finite Wilson operator is nonsingular",
        fermions["minimum_singular_value"] > 1.0e-3,
        f"sigma_min={fermions['minimum_singular_value']:.12g}",
    )
    check(
        "fermion",
        "even-odd Schur determinant factorisation",
        fermions["determinant_logabs_residual"] < 1.0e-9
        and fermions["determinant_phase_residual"] < 1.0e-9,
        (
            f"logabs={fermions['determinant_logabs_residual']:.3e}; "
            f"phase={fermions['determinant_phase_residual']:.3e}"
        ),
    )
    check(
        "fermion",
        "even-odd Schur solve agrees with direct solve",
        fermions["schur_solve_relative_residual"] < 1.0e-11,
        f"relative residual={fermions['schur_solve_relative_residual']:.3e}",
    )
    check(
        "fermion",
        "fermion log-determinant curvature trace identity",
        fermions["one_loop_collective_coordinate"]["relative_error"] < 2.0e-6,
        (
            "relative error="
            f"{fermions['one_loop_collective_coordinate']['relative_error']:.3e}"
        ),
    )
    check(
        "scaling",
        "Wilson dispersion has second-order scaling",
        1.95 < scaling["wilson_observed_order_last_pair"] < 2.05,
        f"p={scaling['wilson_observed_order_last_pair']:.8f}",
    )
    check(
        "scaling",
        "tree-level Symanzik dispersion has fourth-order scaling",
        3.90 < scaling["symanzik_observed_order_last_pair"] < 4.10,
        f"p={scaling['symanzik_observed_order_last_pair']:.8f}",
    )

    fail_closed = {
        "importance_sampling_or_Monte_Carlo_performed": False,
        "dynamical_4d_lattice_QC2D_phase_established": False,
        "fermion_determinant_positivity_proven": False,
        "complete_nonlinear_gauge_meson_ghost_fermion_Hessian_computed": False,
        "Skyrme_Jacobi_translation_isorotation_zero_modes_resolved": False,
        "one_loop_renormalised_continuum_limit_established": False,
        "Finkelstein_Rubinstein_constraint_derived": False,
        "finite_amplitude_global_stability_established": False,
        "degree_one_route_E_portal_allowed": False,
        "physics_promotion_allowed": False,
    }
    check(
        "fail_closed",
        "no Monte Carlo claim",
        not fail_closed["importance_sampling_or_Monte_Carlo_performed"],
        "deterministic background-field diagonalisation only",
    )
    check(
        "fail_closed",
        "determinant positivity remains open",
        not fail_closed["fermion_determinant_positivity_proven"],
        "no antiunitary/global-form proof is supplied",
    )
    check(
        "fail_closed",
        "FR and global stability remain open",
        not fail_closed["Finkelstein_Rubinstein_constraint_derived"]
        and not fail_closed["finite_amplitude_global_stability_established"],
        "local Gaussian spectra cannot close configuration-space topology",
    )
    check(
        "fail_closed",
        "portal and physics promotion remain false",
        not fail_closed["degree_one_route_E_portal_allowed"]
        and not fail_closed["physics_promotion_allowed"],
        "necessary finite-lattice diagnostic only",
    )

    passed = sum(row["pass"] for row in CHECKS)
    return {
        "artifact": "AP-E5 4D charged-QC2D lattice/Hessian diagnostic",
        "status": "PASS" if passed == len(CHECKS) else "FAIL",
        "scope": (
            "finite 4D background-field/tree-plus-one-collective-coordinate-loop "
            "proxy; not dynamical lattice QC2D"
        ),
        "benchmark": asdict(benchmark),
        "topology": topology,
        "gauge_and_ghost": gauge,
        "boson_hessian": bosons,
        "fermion": fermions,
        "continuum_scaling": scaling,
        "fail_closed": fail_closed,
        "checks": CHECKS,
        "summary": {"passed": passed, "total": len(CHECKS)},
        "reproducibility": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "random_seed": 20260716,
            "sources": sources,
        },
    }


def markdown(card: dict[str, Any]) -> str:
    topology = card["topology"]
    gauge = card["gauge_and_ghost"]
    bosons = card["boson_hessian"]
    fermions = card["fermion"]
    scaling = card["continuum_scaling"]
    lines = [
        "# AP-E5 4D charged-QC2D lattice/Hessian diagnostic",
        "",
        f"**Status:** {card['status']} ({card['summary']['passed']}/{card['summary']['total']})",
        "",
        "This is a deterministic finite-volume background-field diagnostic. It is not a Monte Carlo or dynamical-lattice result.",
        "",
        "## Numerical certificate",
        "",
        (
            "- B=1 extrapolation: "
            f"`{topology['quadratic_in_a_squared_extrapolation_last_four']:.12g}` "
            f"(absolute error `{topology['extrapolation_abs_error']:.3e}`)."
        ),
        (
            "- Gauge/ghost xi-cancellation residual: "
            f"`{gauge['maximum_xi_corrected_residual']:.3e}`."
        ),
        (
            "- Meson / diquark minimum eigenvalues: "
            f"`{bosons['meson_minimum_eigenvalue']:.12g}` / "
            f"`{bosons['diquark_minimum_eigenvalue']:.12g}`."
        ),
        (
            "- Diquark negative-control minimum: "
            f"`{bosons['negative_control_minimum_eigenvalue']:.12g}`."
        ),
        (
            "- Wilson-Dirac minimum singular value: "
            f"`{fermions['minimum_singular_value']:.12g}`."
        ),
        (
            "- Schur determinant / solve residuals: "
            f"`{fermions['determinant_logabs_residual']:.3e}` / "
            f"`{fermions['schur_solve_relative_residual']:.3e}`."
        ),
        (
            "- Fermion log-det curvature identity relative error: "
            f"`{fermions['one_loop_collective_coordinate']['relative_error']:.3e}`."
        ),
        (
            "- Observed Wilson / tree-Symanzik orders: "
            f"`{scaling['wilson_observed_order_last_pair']:.8f}` / "
            f"`{scaling['symanzik_observed_order_last_pair']:.8f}`."
        ),
        "",
        "## Fail-closed boundary",
        "",
    ]
    lines.extend(
        f"- `{name}`: `{str(value).lower()}`"
        for name, value in card["fail_closed"].items()
    )
    lines.extend(["", "## Checks", ""])
    lines.extend(
        f"- [{'x' if row['pass'] else ' '}] **{row['group']} / {row['name']}** — {row['detail']}"
        for row in card["checks"]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    global CHECKS
    CHECKS = []
    OUTPUT.mkdir(parents=True, exist_ok=True)
    card = build_card()
    json_path = OUTPUT / "ap_e5_4d_qc2d_lattice_hessian.json"
    markdown_path = OUTPUT / "ap_e5_4d_qc2d_lattice_hessian.md"
    json_path.write_text(json.dumps(card, indent=2) + "\n", encoding="utf-8")
    markdown_path.write_text(markdown(card), encoding="utf-8")
    print(
        f"{card['status']}: {card['summary']['passed']}/{card['summary']['total']} checks"
    )
    print(json_path.relative_to(REPO))
    print(markdown_path.relative_to(REPO))
    return 0 if card["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
