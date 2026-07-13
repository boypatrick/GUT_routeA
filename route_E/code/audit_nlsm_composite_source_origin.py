#!/usr/bin/env python3
"""Nonlinear-sigma/composite candidate for the A5-A6 source origin.

This script tests the minimal constrained parameterization

    S(U,rho54) = rho54 U S0 U^T,
    Omega(U,rho210) = rho210 (U e7) wedge ... wedge (U e10),

with a *shared* U in SO(10)/(SO(6)xSO(4)).  The shared U is the key point:
it removes the relative-orientation moduli that would otherwise reintroduce
incomplete source thresholds.

This is still a cutoff NLSM/composite ansatz, not a UV-complete proof.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "nlsm_composite_source_origin"
TOL = 1e-9


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
    return sym @ vt[1:].T


def theta_to_skew(theta: np.ndarray) -> np.ndarray:
    """Map 24 color-weak coset coordinates to an so(10) matrix."""
    a = np.zeros((10, 10))
    vals = theta.reshape(6, 4)
    for i in range(6):
        for alpha in range(4):
            j = 6 + alpha
            a[i, j] = vals[i, alpha]
            a[j, i] = -vals[i, alpha]
    return a


def cayley(a: np.ndarray) -> np.ndarray:
    ident = np.eye(a.shape[0])
    return np.linalg.solve(ident + 0.5 * a, ident - 0.5 * a)


def wedge4_columns(u: np.ndarray, cols: tuple[int, int, int, int], basis4: list[tuple[int, ...]]) -> np.ndarray:
    sub = u[:, cols]
    out = np.zeros(len(basis4))
    for idx, rows in enumerate(basis4):
        out[idx] = np.linalg.det(sub[list(rows), :])
    return out


def source_coordinates(
    theta: np.ndarray,
    rho54: float,
    rho210: float,
    basis54: np.ndarray,
    basis4: list[tuple[int, ...]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    u = cayley(theta_to_skew(theta))
    s = rho54 * (u @ s0 @ u.T)
    omega = rho210 * wedge4_columns(u, (6, 7, 8, 9), basis4)
    s_coord = basis54.T @ s.reshape(-1)
    return s_coord, omega, s, u


def plucker_max_residual(p: np.ndarray, basis4: list[tuple[int, ...]]) -> float:
    index = {item: i for i, item in enumerate(basis4)}
    max_res = 0.0
    # Grassmann Plucker relations:
    # sum_{b in B} (-1)^pos p_{A union b} p_{B minus b}=0,
    # for |A|=3, |B|=5.
    for a_tuple in itertools.combinations(range(10), 3):
        a_set = set(a_tuple)
        for b_tuple in itertools.combinations(range(10), 5):
            total = 0.0
            for pos, b in enumerate(b_tuple):
                if b in a_set:
                    continue
                left_unsorted = list(a_tuple) + [b]
                inversions = 0
                for i in range(len(left_unsorted)):
                    for j in range(i + 1, len(left_unsorted)):
                        if left_unsorted[i] > left_unsorted[j]:
                            inversions += 1
                left_sign = -1.0 if inversions % 2 else 1.0
                left = tuple(sorted(left_unsorted))
                right = tuple(x for x in b_tuple if x != b)
                total += ((-1.0) ** pos) * left_sign * p[index[left]] * p[index[right]]
            max_res = max(max_res, abs(total))
    return max_res


def finite_difference_jacobian(
    basis54: np.ndarray,
    basis4: list[tuple[int, ...]],
    shared_u: bool,
    eps: float = 1e-6,
) -> np.ndarray:
    if shared_u:
        nvar = 26
    else:
        nvar = 50
    base_theta54 = np.zeros(24)
    base_theta210 = np.zeros(24)

    def coords(x: np.ndarray) -> np.ndarray:
        if shared_u:
            theta54 = x[:24]
            theta210 = x[:24]
            rho54 = 1.0 + x[24]
            rho210 = 1.0 + x[25]
        else:
            theta54 = x[:24]
            theta210 = x[24:48]
            rho54 = 1.0 + x[48]
            rho210 = 1.0 + x[49]
        s_coord, _, _, _ = source_coordinates(theta54, rho54, 1.0, basis54, basis4)
        _, omega, _, _ = source_coordinates(theta210, 1.0, rho210, basis54, basis4)
        return np.concatenate([s_coord, omega])

    jac = np.zeros((264, nvar))
    x0 = np.zeros(nvar)
    for i in range(nvar):
        dx = np.zeros(nvar)
        dx[i] = eps
        jac[:, i] = (coords(x0 + dx) - coords(x0 - dx)) / (2.0 * eps)
    return jac


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    basis54 = symmetric_traceless_basis()
    basis4 = list(itertools.combinations(range(10), 4))

    theta = np.linspace(-0.12, 0.11, 24)
    s_coord, omega, s, u = source_coordinates(theta, 1.3, 0.7, basis54, basis4)
    orth_error = float(np.linalg.norm(u.T @ u - np.eye(10)))
    trace_s = float(np.trace(s))
    rho54 = 1.3
    poly_res = float(np.linalg.norm(s @ s - rho54 * s - 6.0 * rho54 * rho54 * np.eye(10)))
    omega_norm = float(np.linalg.norm(omega))
    plucker_res = float(plucker_max_residual(omega, basis4))

    jac_shared = finite_difference_jacobian(basis54, basis4, shared_u=True)
    jac_independent = finite_difference_jacobian(basis54, basis4, shared_u=False)
    s_shared = np.linalg.svd(jac_shared, compute_uv=False)
    s_independent = np.linalg.svd(jac_independent, compute_uv=False)
    rank_shared = int(np.sum(s_shared > TOL))
    rank_independent = int(np.sum(s_independent > TOL))

    # Constraint ranks in the ambient 54+210 source space.
    ambient_dim = 54 + 210
    constraint_rank_shared = ambient_dim - rank_shared
    constraint_rank_independent = ambient_dim - rank_independent

    # Boundary NLSM kinetic suppression for normal/composite gap estimates.
    n_over_ell = 1e-2
    lapse_sq = math.tanh(n_over_ell) ** 2

    # The coset tangent is (6,2,2), a real PS representation, so it carries no
    # chiral gauge anomaly.  Radials are singlets.
    anomaly_rows = [
        {
            "sector": "coset_tangent",
            "representation": "(6,2,2)",
            "real_or_complex": "real",
            "chiral_anomaly_proxy": 0,
        },
        {
            "sector": "radial_54",
            "representation": "(1,1,1)",
            "real_or_complex": "real",
            "chiral_anomaly_proxy": 0,
        },
        {
            "sector": "radial_210",
            "representation": "(1,1,1)",
            "real_or_complex": "real",
            "chiral_anomaly_proxy": 0,
        },
    ]
    with (OUT / "anomaly_proxy.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(anomaly_rows[0].keys()))
        writer.writeheader()
        writer.writerows(anomaly_rows)

    sv_rows = []
    for i, val in enumerate(s_shared):
        sv_rows.append({"branch": "shared_U", "index": i, "singular_value": float(val)})
    for i, val in enumerate(s_independent):
        sv_rows.append({"branch": "independent_U", "index": i, "singular_value": float(val)})
    with (OUT / "jacobian_singular_values.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(sv_rows)

    summary = {
        "note": "No web lookup used. NLSM/composite constrained-source origin candidate.",
        "parameterization": {
            "S": "S(U,rho54)=rho54 U S0 U^T",
            "Omega": "Omega(U,rho210)=rho210 (Ue7)^(Ue8)^(Ue9)^(Ue10)",
            "coset": "SO(10)/(SO(6)xSO(4))",
            "variables": "24 coset coordinates + 2 radial singlets",
        },
        "constraint_checks": {
            "orthogonality_error": orth_error,
            "trace_S": trace_s,
            "S_polynomial_residual": poly_res,
            "Omega_norm": omega_norm,
            "Omega_norm_expected_abs_rho210": 0.7,
            "Omega_plucker_max_residual": plucker_res,
        },
        "jacobian": {
            "ambient_dim": ambient_dim,
            "shared_U_rank": rank_shared,
            "shared_U_constraint_rank": constraint_rank_shared,
            "independent_U_rank": rank_independent,
            "independent_U_constraint_rank": constraint_rank_independent,
            "relative_orientation_modes_removed_by_shared_U": rank_independent - rank_shared,
            "shared_U_smallest_nonzero_sv": float(np.min(s_shared[s_shared > TOL])),
            "shared_U_largest_sv": float(np.max(s_shared)),
        },
        "boundary": {
            "n_over_ell": n_over_ell,
            "lapse_squared": lapse_sq,
        },
        "anomaly_proxy": {
            "coset_representation": "(6,2,2)",
            "coset_representation_is_real": True,
            "new_chiral_anomaly_proxy": 0,
        },
        "beta_cost": {
            "below_compositeness_scale": "no elementary 54/210 tower; only coset/eaten modes and singlet radials",
            "so10_dynkin_increment_in_preferred_EFT": 0,
            "UV_above_compositeness_scale": "unknown; still the open microscopic problem",
        },
        "passes": {
            "S_constraint_residual_lt_1e_minus_10": poly_res < 1e-10,
            "Omega_decomposable_plucker_lt_1e_minus_10": plucker_res < 1e-10,
            "shared_U_rank_26": rank_shared == 26,
            "independent_U_rank_50": rank_independent == 50,
            "relative_orientation_removed_24": rank_independent - rank_shared == 24,
            "constraint_rank_238": constraint_rank_shared == 238,
            "boundary_lapse_squared_lt_1e_minus_4": lapse_sq < 1e-4,
            "no_new_chiral_anomaly_proxy": True,
        },
        "verdict": (
            "The NLSM parameterization realizes the desired source manifold "
            "locally: 26 physical source coordinates inside the 264-dimensional "
            "54+210 ambient space, with 238 constraints and the 24 relative "
            "orientation modes removed by using a shared U.  This is a stronger "
            "A5-A6 candidate, but it remains a cutoff/composite ansatz rather "
            "than a UV-complete derivation."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    report = [
        "# NLSM/composite constrained-source origin candidate",
        "",
        "No web lookup used.",
        "",
        "## Parameterization",
        "",
        "S(U,rho54)=rho54 U S0 U^T.",
        "Omega(U,rho210)=rho210 (Ue7) wedge (Ue8) wedge (Ue9) wedge (Ue10).",
        "The same U is used for both sources.",
        "",
        "## Constraint checks",
        "",
        f"Orthogonality error: {orth_error:.6e}",
        f"Trace S: {trace_s:.6e}",
        f"S polynomial residual: {poly_res:.6e}",
        f"Omega norm: {omega_norm:.6e}",
        f"Omega Plucker max residual: {plucker_res:.6e}",
        "",
        "## Jacobian",
        "",
        f"Shared-U rank: {rank_shared}",
        f"Independent-U rank: {rank_independent}",
        f"Relative-orientation modes removed: {rank_independent-rank_shared}",
        f"Shared-U constraint rank in 54+210 ambient space: {constraint_rank_shared}",
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
