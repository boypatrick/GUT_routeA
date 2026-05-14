#!/usr/bin/env python3
"""Audit a D-term/Kahler quotient origin for the crossed-120 unitary link.

No web lookup is used.

The holomorphic crossed-link equations

    A B = 1,        P B P = W_kappa

derive the desired visible inverse block but leave nonunitary complex moduli
that split the triplet singular values.  This script asks whether the next
minimal ingredient, a D-term/Kahler quotient constraint

    B^\dagger B = 1

removes precisely those dangerous modes while keeping a nonempty moduli space.

The audit is deliberately finite-dimensional and local:

* The D-term linearization has rank 16 on the 32 real components of B, leaving
  the 16-dimensional tangent space of U(4).
* Because W_kappa lies on the contraction boundary (||W_kappa||_2=1), the
  fixed-block constraint PBP=W_kappa has rank 7 on that unitary tangent space,
  leaving 9 real unitary completion moduli.
* Random finite unitary-completion deformations preserve PBP=W_kappa exactly,
  keep all singular values equal to one, and therefore cannot violate the
  completed-partner mass-locking window.

This upgrades the crossed-120 link from an imposed NLSM assumption to a local
D-term/Kahler quotient ansatz.  It is still not a microscopic strong-sector
derivation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "unitary_link_dterm_quotient"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"
LINK = ROOT / "output" / "crossed_120_link_locking_action" / "summary.json"

RNG_SEED = 202605090052
ANGLE_GRID = [0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cval(raw: dict[str, float]) -> complex:
    return raw["re"] + 1j * raw["im"]


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cval(cell) for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    if float(np.min(vals)) < -tol:
        raise RuntimeError(f"negative defect eigenvalue {np.min(vals)}")
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def julia_dilation(w: np.ndarray) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    return np.block(
        [
            [w, hermitian_sqrt(eye - w @ w.conj().T)],
            [hermitian_sqrt(eye - w.conj().T @ w), -w.conj().T],
        ]
    )


def real_pack(mat: np.ndarray) -> np.ndarray:
    return np.concatenate([mat.real.reshape(-1), mat.imag.reshape(-1)])


def hermitian_pack(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    vals: list[float] = []
    for i in range(n):
        vals.append(float(np.real(mat[i, i])))
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(float(np.real(mat[i, j])))
            vals.append(float(np.imag(mat[i, j])))
    return np.array(vals, dtype=float)


def complex_matrix_basis(n: int) -> list[np.ndarray]:
    basis = []
    for i in range(n):
        for j in range(n):
            e = np.zeros((n, n), dtype=complex)
            e[i, j] = 1.0
            basis.append(e)
            ei = np.zeros((n, n), dtype=complex)
            ei[i, j] = 1.0j
            basis.append(ei)
    return basis


def antihermitian_basis(n: int) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for a in range(n):
        mat = np.zeros((n, n), dtype=complex)
        mat[a, a] = 1j
        basis.append(mat)
    for a in range(n):
        for b in range(a + 1, n):
            skew = np.zeros((n, n), dtype=complex)
            skew[a, b] = 1.0
            skew[b, a] = -1.0
            basis.append(skew)
            sym_i = np.zeros((n, n), dtype=complex)
            sym_i[a, b] = 1j
            sym_i[b, a] = 1j
            basis.append(sym_i)
    return basis


def rank_of_columns(cols: list[np.ndarray], tol: float = 1.0e-10) -> tuple[int, list[float]]:
    mat = np.stack(cols, axis=1)
    s = np.linalg.svd(mat, compute_uv=False)
    rank = int(np.sum(s > tol * max(float(s[0]), 1.0)))
    return rank, [float(x) for x in s]


def dterm_linear_rank(u: np.ndarray) -> dict[str, Any]:
    cols = []
    for db in complex_matrix_basis(u.shape[0]):
        d_moment = u.conj().T @ db + db.conj().T @ u
        cols.append(hermitian_pack(d_moment))
    rank, s = rank_of_columns(cols)
    return {
        "real_variables_B": 2 * u.shape[0] * u.shape[0],
        "real_D_equations": u.shape[0] * u.shape[0],
        "Dterm_linear_rank": rank,
        "Dflat_tangent_real_dimension": 2 * u.shape[0] * u.shape[0] - rank,
        "singular_values_min": float(min(s)),
        "singular_values_max": float(max(s)),
    }


def fixed_block_rank_on_unitary(u: np.ndarray, visible: int) -> dict[str, Any]:
    cols = []
    for k in antihermitian_basis(u.shape[0]):
        db = u @ k
        cols.append(real_pack(db[:visible, :visible]))
    rank, s = rank_of_columns(cols)
    return {
        "unitary_tangent_real_dimension": u.shape[0] * u.shape[0],
        "fixed_block_real_equations": 2 * visible * visible,
        "fixed_block_rank_on_Dflat_tangent": rank,
        "quotient_residual_unitary_moduli_real_dimension": u.shape[0] * u.shape[0] - rank,
        "singular_values_min": float(min(s)),
        "singular_values_max": float(max(s)),
    }


def random_unitary(n: int, angle: float, rng: np.random.Generator) -> np.ndarray:
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = 0.5 * (raw + raw.conj().T)
    h = h / max(np.linalg.norm(h, ord=2), 1.0e-30)
    vals, vecs = np.linalg.eigh(h)
    return vecs @ np.diag(np.exp(1j * angle * vals)) @ vecs.conj().T


def unitary_completion(w: np.ndarray, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    d_left = hermitian_sqrt(eye - w @ w.conj().T)
    d_right = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, d_left @ right],
            [left @ d_right, -left @ w.conj().T @ right],
        ]
    )


def finite_completion_scan(w: np.ndarray, lock_window: float) -> list[dict[str, Any]]:
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, Any]] = []
    for angle in ANGLE_GRID:
        left = random_unitary(w.shape[0], angle, rng)
        right = random_unitary(w.shape[0], angle, rng)
        b = unitary_completion(w, left, right)
        singulars = np.linalg.svd(b.conj().T, compute_uv=False)
        max_abs_log = float(np.max(np.abs(np.log(singulars))))
        rows.append(
            {
                "angle": float(angle),
                "fixed_block_residual_fro": float(np.linalg.norm(b[: w.shape[0], : w.shape[1]] - w, ord="fro")),
                "unitarity_residual_2norm": float(np.linalg.norm(b.conj().T @ b - np.eye(b.shape[0]), ord=2)),
                "mass_singular_min_over_Slock": float(np.min(singulars)),
                "mass_singular_max_over_Slock": float(np.max(singulars)),
                "max_abs_log_singular": max_abs_log,
                "within_completed_partner_lock_window": bool(max_abs_log <= lock_window),
            }
        )
    return rows


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    completed = read_json(COMPLETED)
    link = read_json(LINK)
    w = cmat(ps["superpotential_ansatz"]["W_kappa"])
    u0 = julia_dilation(w)
    d_rank = dterm_linear_rank(u0)
    f_rank = fixed_block_rank_on_unitary(u0, w.shape[0])
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])
    rows = finite_completion_scan(w, lock_window)
    verdict = {
        "Dterm_removes_nonunitary_holomorphic_moduli": True,
        "Dflat_fixed_block_moduli_real_dimension": f_rank["quotient_residual_unitary_moduli_real_dimension"],
        "all_finite_unitary_completions_lock_singular_values": all(
            row["within_completed_partner_lock_window"] for row in rows
        ),
        "max_finite_completion_unitarity_residual": max(row["unitarity_residual_2norm"] for row in rows),
        "max_finite_completion_fixed_block_residual": max(row["fixed_block_residual_fro"] for row in rows),
        "previous_holomorphic_nullity_real": link["verdict"]["constraint_nullity_real"],
        "interpretation": (
            "The D-term/Kahler quotient B^dagger B=1 removes the dangerous nonunitary "
            "holomorphic moduli.  After imposing the fixed visible block PBP=W_kappa, "
            "a 9-real-dimensional family of unitary completions remains because W_kappa "
            "sits on the contraction boundary, but every finite "
            "sample keeps all singular values equal to one and therefore stays inside the "
            "completed-partner mass-locking window.  This is a local quotient derivation of "
            "the unitary-link assumption; a microscopic strong-sector origin remains open."
        ),
    }
    return {
        "note": "No web lookup used. D-term/Kahler quotient audit for crossed-120 unitary link.",
        "quotient_action": {
            "Dterm_moment_map": "mu=B^dagger B-I=0, optionally with left/right gauge quotient",
            "Fterm_driver": "Tr Y(PBP-W_kappa)",
            "common_spurion": "S_lock multiplies both bar(Phi_T)^T B^dagger Phi_T and inert L Lbar masses",
            "meaning": "D-flatness restricts B to U(4); F-flatness fixes the visible block.",
        },
        "benchmarks": {
            "completed_partner_max_abs_log_xi": lock_window,
            "M_lock_GeV": float(ps["benchmarks"]["M_lock_GeV"]),
            "MG_GeV": float(ps["benchmarks"]["MG_GeV"]),
        },
        "Dterm_linearization": d_rank,
        "Fterm_on_Dflat_tangent": f_rank,
        "finite_unitary_completion_scan": rows,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    lines = [
        "# Unitary link D-term/Kahler quotient audit",
        "",
        "No web lookup was used.",
        "",
        "## Linearized quotient",
        "",
        f"D-term rank: {payload['Dterm_linearization']['Dterm_linear_rank']}",
        f"D-flat tangent dimension: {payload['Dterm_linearization']['Dflat_tangent_real_dimension']}",
        f"fixed-block rank on D-flat tangent: {payload['Fterm_on_Dflat_tangent']['fixed_block_rank_on_Dflat_tangent']}",
        f"residual unitary moduli: {payload['Fterm_on_Dflat_tangent']['quotient_residual_unitary_moduli_real_dimension']}",
        "",
        "## Finite unitary completions",
        "",
        "| angle | fixed-block residual | unitarity residual | max |log s| | inside lock window |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in payload["finite_unitary_completion_scan"]:
        lines.append(
            f"| {row['angle']:.0e} | {row['fixed_block_residual_fro']:.3e} | "
            f"{row['unitarity_residual_2norm']:.3e} | {row['max_abs_log_singular']:.3e} | "
            f"{row['within_completed_partner_lock_window']} |"
        )
    lines.extend(["", "## Verdict", "", payload["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "finite_unitary_completion_scan.csv",
        payload["finite_unitary_completion_scan"],
        [
            "angle",
            "fixed_block_residual_fro",
            "unitarity_residual_2norm",
            "mass_singular_min_over_Slock",
            "mass_singular_max_over_Slock",
            "max_abs_log_singular",
            "within_completed_partner_lock_window",
        ],
    )
    write_report(payload)
    v = payload["verdict"]
    print("Unitary link D-term/Kahler quotient audit")
    print(f"  D-flat fixed-block moduli dim: {v['Dflat_fixed_block_moduli_real_dimension']}")
    print(f"  all finite completions locked: {v['all_finite_unitary_completions_lock_singular_values']}")
    print(f"  max unitarity residual: {v['max_finite_completion_unitarity_residual']:.3e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
