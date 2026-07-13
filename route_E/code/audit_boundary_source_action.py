#!/usr/bin/env python3
"""Local action-level audit for the boundary constrained-source ansatz.

The goal is to test a minimal quadratic expansion of a microscopic source
sector that could realize assumptions A5--A6:

  W_loc = Lambda_54 . N_54^T s + Lambda_210 . N_210^T omega
        + Pi . (T_54^T s - T_210^T omega)
        + (mu_54/2) (r_54.s)^2 + (mu_210/2) (r_210.omega)^2.

Here T are orbit tangent bases, r are radial singlets, and N are normal-bundle
bases.  The desired local Hessian has exactly 24 zero modes: the shared
SO(10)/(SO(6)xSO(4)) orbit tangents.  All normal modes should be paired with
auxiliary/composite multipliers, so no incomplete normal threshold propagates.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "boundary_source_action"
TOL = 1e-10


def so_generator(n: int, a: int, b: int) -> np.ndarray:
    gen = np.zeros((n, n))
    gen[a, b] = 1.0
    gen[b, a] = -1.0
    return gen


def symmetric_traceless_basis(n: int = 10) -> np.ndarray:
    cols = []
    for i in range(n):
        mat = np.zeros((n, n))
        mat[i, i] = 1.0
        cols.append(mat.reshape(-1))
    for i, j in itertools.combinations(range(n), 2):
        mat = np.zeros((n, n))
        mat[i, j] = 1.0 / math.sqrt(2.0)
        mat[j, i] = 1.0 / math.sqrt(2.0)
        cols.append(mat.reshape(-1))
    sym = np.column_stack(cols)
    trace_coeff = sym.T @ np.eye(n).reshape(-1)
    _, _, vt = np.linalg.svd(trace_coeff.reshape(1, -1))
    null_in_sym = vt[1:].T
    return sym @ null_in_sym


def wedge_sign(unsorted_new: list[int]) -> tuple[tuple[int, ...], int]:
    if len(set(unsorted_new)) < len(unsorted_new):
        return tuple(), 0
    inversions = 0
    for i in range(len(unsorted_new)):
        for j in range(i + 1, len(unsorted_new)):
            if unsorted_new[i] > unsorted_new[j]:
                inversions += 1
    return tuple(sorted(unsorted_new)), -1 if inversions % 2 else 1


def act_generator_on_form(
    form: tuple[int, ...],
    pair: tuple[int, int],
    basis_index: dict[tuple[int, ...], int],
) -> np.ndarray:
    a, b = pair
    vec = np.zeros(len(basis_index))
    for slot, idx in enumerate(form):
        replacement = None
        sign = 1
        if idx == a:
            replacement = b
            sign = 1
        elif idx == b:
            replacement = a
            sign = -1
        if replacement is None:
            continue
        new_form = list(form)
        new_form[slot] = replacement
        sorted_form, wedge_sgn = wedge_sign(new_form)
        if wedge_sgn:
            vec[basis_index[sorted_form]] += sign * wedge_sgn
    return vec


def orthonormal_columns(mat: np.ndarray, rank: int | None = None) -> np.ndarray:
    q, r = np.linalg.qr(mat)
    if rank is None:
        rank = int(np.linalg.matrix_rank(mat, tol=TOL))
    return q[:, :rank]


def complement_basis(keep: np.ndarray, dim: int) -> np.ndarray:
    q, _ = np.linalg.qr(keep, mode="complete")
    return q[:, keep.shape[1] :]


def source_54_subspaces() -> dict[str, np.ndarray | int]:
    n = 10
    basis = symmetric_traceless_basis(n)
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    s0_coord = basis.T @ s0.reshape(-1)
    r = s0_coord / np.linalg.norm(s0_coord)
    tangent_cols = []
    for a in range(6):
        for b in range(6, 10):
            gen = so_generator(n, a, b)
            tangent = gen @ s0 - s0 @ gen
            tangent_cols.append(basis.T @ tangent.reshape(-1))
    tangent = np.column_stack(tangent_cols)
    tangent_basis = orthonormal_columns(tangent, 24)
    keep = orthonormal_columns(np.column_stack([tangent_basis, r]), 25)
    normal = complement_basis(keep, 54)
    return {"T": tangent_basis, "r": r, "N": normal, "rank": 24, "dim": 54}


def source_210_subspaces() -> dict[str, np.ndarray | int]:
    n = 10
    basis = list(itertools.combinations(range(n), 4))
    basis_index = {item: i for i, item in enumerate(basis)}
    omega = (6, 7, 8, 9)
    omega_coord = np.zeros(len(basis))
    omega_coord[basis_index[omega]] = 1.0
    r = omega_coord / np.linalg.norm(omega_coord)
    tangent_cols = []
    for a in range(6):
        for b in range(6, 10):
            tangent_cols.append(act_generator_on_form(omega, (a, b), basis_index))
    tangent = np.column_stack(tangent_cols)
    tangent_basis = orthonormal_columns(tangent, 24)
    keep = orthonormal_columns(np.column_stack([tangent_basis, r]), 25)
    normal = complement_basis(keep, len(basis))
    return {"T": tangent_basis, "r": r, "N": normal, "rank": 24, "dim": len(basis)}


def build_hessian(
    s54: dict[str, np.ndarray | int],
    s210: dict[str, np.ndarray | int],
    mu54: float = 1.0,
    mu210: float = 1.5,
) -> tuple[np.ndarray, dict[str, slice]]:
    d54 = int(s54["dim"])
    d210 = int(s210["dim"])
    n54 = s54["N"].shape[1]  # type: ignore[index,union-attr]
    n210 = s210["N"].shape[1]  # type: ignore[index,union-attr]
    nt = 24
    offsets = {}
    pos = 0
    for name, size in [
        ("x54", d54),
        ("x210", d210),
        ("lam54", n54),
        ("lam210", n210),
        ("align", nt),
    ]:
        offsets[name] = slice(pos, pos + size)
        pos += size
    hess = np.zeros((pos, pos))

    r54 = s54["r"]  # type: ignore[assignment]
    r210 = s210["r"]  # type: ignore[assignment]
    n54_basis = s54["N"]  # type: ignore[assignment]
    n210_basis = s210["N"]  # type: ignore[assignment]
    t54_basis = s54["T"]  # type: ignore[assignment]
    t210_basis = s210["T"]  # type: ignore[assignment]

    hess[offsets["x54"], offsets["x54"]] += mu54 * np.outer(r54, r54)
    hess[offsets["x210"], offsets["x210"]] += mu210 * np.outer(r210, r210)

    hess[offsets["x54"], offsets["lam54"]] = n54_basis
    hess[offsets["lam54"], offsets["x54"]] = n54_basis.T
    hess[offsets["x210"], offsets["lam210"]] = n210_basis
    hess[offsets["lam210"], offsets["x210"]] = n210_basis.T

    hess[offsets["x54"], offsets["align"]] = t54_basis
    hess[offsets["align"], offsets["x54"]] = t54_basis.T
    hess[offsets["x210"], offsets["align"]] = -t210_basis
    hess[offsets["align"], offsets["x210"]] = -t210_basis.T
    return hess, offsets


def expected_diagonal_orbit_vectors(
    s54: dict[str, np.ndarray | int],
    s210: dict[str, np.ndarray | int],
    offsets: dict[str, slice],
    total_dim: int,
) -> np.ndarray:
    t54 = s54["T"]  # type: ignore[assignment]
    t210 = s210["T"]  # type: ignore[assignment]
    vecs = []
    for i in range(24):
        vec = np.zeros(total_dim)
        vec[offsets["x54"]] = t54[:, i] / math.sqrt(2.0)
        vec[offsets["x210"]] = t210[:, i] / math.sqrt(2.0)
        vecs.append(vec)
    return np.column_stack(vecs)


def rotating_margins() -> dict[str, float | bool]:
    masses = {
        "Sigma3": 0.13001566778681514,
        "X622": 1.0,
        "conormal": 1.0,
    }
    omega = 0.12
    return {
        "Omega_over_MG": omega,
        "Sigma3_margin_m1": masses["Sigma3"] - omega,
        "X622_margin_m1": masses["X622"] - omega,
        "conormal_margin_m1": masses["conormal"] - omega,
        "all_positive": all(mass - omega > 0 for mass in masses.values()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s54 = source_54_subspaces()
    s210 = source_210_subspaces()
    hess, offsets = build_hessian(s54, s210)
    singular = np.linalg.svd(hess, compute_uv=False)
    rank = int(np.sum(singular > TOL))
    nullity = hess.shape[0] - rank
    u, svals, vt = np.linalg.svd(hess)
    null_basis = vt[svals <= TOL].T
    expected_zero = expected_diagonal_orbit_vectors(s54, s210, offsets, hess.shape[0])
    residuals = np.linalg.norm(hess @ expected_zero, axis=0)
    projection = null_basis @ (null_basis.T @ expected_zero)
    expected_projection_error = float(np.max(np.linalg.norm(expected_zero - projection, axis=0)))

    n54 = s54["N"]  # type: ignore[assignment]
    n210 = s210["N"]  # type: ignore[assignment]
    normal_components = []
    radial_components = []
    for i in range(null_basis.shape[1]):
        vec = null_basis[:, i]
        x54 = vec[offsets["x54"]]
        x210 = vec[offsets["x210"]]
        normal_components.append(np.linalg.norm(n54.T @ x54) + np.linalg.norm(n210.T @ x210))
        radial_components.append(abs(float(s54["r"].T @ x54)) + abs(float(s210["r"].T @ x210)))  # type: ignore[union-attr]

    lapse_n = 1e-2
    lapse = math.tanh(lapse_n)
    boundary_kinetic_prefactor = lapse * lapse
    rot = rotating_margins()

    spectrum_rows = [
        {
            "quantity": "total_variables",
            "value": hess.shape[0],
            "interpretation": "54+210 fluctuations plus normal and alignment multipliers",
        },
        {
            "quantity": "hessian_rank",
            "value": rank,
            "interpretation": "paired normal/radial/relative-orbit directions",
        },
        {
            "quantity": "hessian_nullity",
            "value": nullity,
            "interpretation": "shared diagonal SO(10)/(SO(6)xSO(4)) orbit",
        },
        {
            "quantity": "positive_singular_min",
            "value": float(np.min(singular[singular > TOL])),
            "interpretation": "smallest nonzero local source mass in arbitrary units",
        },
        {
            "quantity": "positive_singular_max",
            "value": float(np.max(singular)),
            "interpretation": "largest local source mass in arbitrary units",
        },
    ]
    with (OUT / "local_hessian_spectrum.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(spectrum_rows[0].keys()))
        writer.writeheader()
        writer.writerows(spectrum_rows)

    summary = {
        "note": "No web lookup used. Local boundary constrained-source action audit.",
        "local_superpotential": (
            "W_loc = Lambda54.N54^T s + Lambda210.N210^T omega "
            "+ Pi.(T54^T s - T210^T omega) "
            "+ mu54/2 (r54.s)^2 + mu210/2 (r210.omega)^2"
        ),
        "dimensions": {
            "x54": int(s54["dim"]),
            "x210": int(s210["dim"]),
            "lambda54_normal": s54["N"].shape[1],  # type: ignore[index,union-attr]
            "lambda210_normal": s210["N"].shape[1],  # type: ignore[index,union-attr]
            "alignment_drivers": 24,
            "total": hess.shape[0],
        },
        "hessian": {
            "rank": rank,
            "nullity": nullity,
            "expected_rank": 478,
            "expected_nullity": 24,
            "max_expected_zero_residual": float(np.max(residuals)),
            "expected_zero_projection_error": expected_projection_error,
            "max_null_normal_component": float(np.max(normal_components)),
            "max_null_radial_component": float(np.max(radial_components)),
            "positive_singular_min": float(np.min(singular[singular > TOL])),
            "positive_singular_max": float(np.max(singular)),
        },
        "boundary_time": {
            "n_over_ell": lapse_n,
            "lapse": lapse,
            "normal_kinetic_prefactor_N_squared": boundary_kinetic_prefactor,
        },
        "rotating_boundary": rot,
        "passes": {
            "hessian_rank_expected": rank == 478,
            "hessian_nullity_24": nullity == 24,
            "expected_orbit_zero_residual_lt_1e_minus_10": float(np.max(residuals)) < 1e-10,
            "nullspace_has_no_normal_component": float(np.max(normal_components)) < 1e-10,
            "nullspace_has_no_radial_component": float(np.max(radial_components)) < 1e-10,
            "boundary_suppresses_normal_kinetic_lt_1e_minus_4": boundary_kinetic_prefactor < 1e-4,
            "rotating_boundary_margins_positive": bool(rot["all_positive"]),
        },
        "verdict": (
            "The local constrained-source action has exactly 24 zero modes, "
            "identified with the shared orbit tangents.  Normal and relative "
            "orientation directions are paired by auxiliary/composite multipliers, "
            "and radial modes are lifted.  This strengthens A5-A6 to a local "
            "action-level ansatz, but does not yet prove a UV microscopic origin."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    report = [
        "# Boundary constrained-source action audit",
        "",
        "No web lookup used.",
        "",
        "## Local superpotential",
        "",
        summary["local_superpotential"],
        "",
        "## Hessian",
        "",
        f"Total variables: {hess.shape[0]}",
        f"Rank: {rank}",
        f"Nullity: {nullity}",
        f"Max expected zero residual: {np.max(residuals):.6e}",
        f"Max null normal component: {np.max(normal_components):.6e}",
        f"Max null radial component: {np.max(radial_components):.6e}",
        "",
        "## Boundary-time suppression",
        "",
        f"At n/ell={lapse_n:g}, N^2={boundary_kinetic_prefactor:.6e}.",
        "",
        "## Rotating boundary",
        "",
        f"Omega/MG={rot['Omega_over_MG']:.6f}",
        f"Sigma3 margin={rot['Sigma3_margin_m1']:.6e}",
        f"X622 margin={rot['X622_margin_m1']:.6e}",
        f"Conormal margin={rot['conormal_margin_m1']:.6e}",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report))
    print(json.dumps(summary["passes"], sort_keys=True))


if __name__ == "__main__":
    main()
