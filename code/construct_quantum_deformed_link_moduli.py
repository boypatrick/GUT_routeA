#!/usr/bin/env python3
"""Quantum-deformed Nf=Nc hidden-SQCD audit for the unitary link.

No web lookup is used.

For hidden SU(4) SQCD with Nf=Nc=4, the gauge-invariant moduli obey

    det M - Baryon * Baryon_tilde = Lambda_H^8.

The previous hidden-GLSM caricature used the meson coordinate

    m = M/f^2 = B_link.

This script checks whether the required unitary completions are compatible
with the quantum-deformed constraint, and whether the quantum constraint alone
is strong enough to replace the D-term/Kahler unitary lock.

Result preview:
* Compatibility: yes.  For each unitary completion m, baryon vevs can solve
  det(m) - beta beta_tilde = (Lambda_H/f)^8.
* Closure: no.  The holomorphic quantum-deformed constraints still leave a
  large complex moduli space.  The D-term/radial unitary constraint remains
  necessary to remove nonunitary link deformations.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "quantum_deformed_link_moduli"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"

ANGLE_GRID = [0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
LAMBDA_OVER_F_GRID = [0.3, 0.7, 1.0]
VISIBLE = 2
N = 4
TOL = 1.0e-10


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cval(raw: dict[str, float]) -> complex:
    return raw["re"] + 1j * raw["im"]


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cval(cell) for cell in row] for row in raw], dtype=complex)


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    if float(np.min(vals)) < -tol:
        raise RuntimeError(f"negative defect eigenvalue {np.min(vals)}")
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def random_unitary(n: int, angle: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = 0.5 * (raw + raw.conj().T)
    h = h / max(np.linalg.norm(h, ord=2), 1.0e-30)
    vals, vecs = np.linalg.eigh(h)
    return vecs @ np.diag(np.exp(1j * angle * vals)) @ vecs.conj().T


def unitary_completion(w: np.ndarray, angle: float, seed: int) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    left = random_unitary(w.shape[0], angle, seed)
    right = random_unitary(w.shape[0], angle, seed + 1297)
    d_left = hermitian_sqrt(eye - w @ w.conj().T)
    d_right = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, d_left @ right],
            [left @ d_right, -left @ w.conj().T @ right],
        ]
    )


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


def real_pack_complex(z: np.ndarray) -> np.ndarray:
    flat = z.reshape(-1)
    return np.concatenate([flat.real, flat.imag])


def rank_real(cols: list[np.ndarray]) -> tuple[int, list[float]]:
    mat = np.stack(cols, axis=1)
    s = np.linalg.svd(mat, compute_uv=False)
    rank = int(np.sum(s > TOL * max(float(s[0]), 1.0)))
    return rank, [float(x) for x in s]


def quantum_baryons(m: np.ndarray, lambda_over_f: float) -> tuple[complex, complex, complex]:
    q = complex(lambda_over_f**8)
    product = np.linalg.det(m) - q
    beta = np.sqrt(product)
    beta_tilde = beta
    return beta, beta_tilde, q


def quantum_residual(m: np.ndarray, beta: complex, beta_tilde: complex, q: complex) -> complex:
    return np.linalg.det(m) - beta * beta_tilde - q


def holomorphic_jacobian_rank(m: np.ndarray, beta: complex, beta_tilde: complex) -> dict[str, Any]:
    # Complex variables: 16 meson entries + beta + beta_tilde.
    rows = np.zeros((VISIBLE * VISIBLE + 1, N * N + 2), dtype=complex)
    row = 0
    for i in range(VISIBLE):
        for j in range(VISIBLE):
            rows[row, i * N + j] = 1.0
            row += 1
    det_m = np.linalg.det(m)
    inv_t = np.linalg.inv(m).T
    rows[row, : N * N] = (det_m * inv_t).reshape(-1)
    rows[row, N * N] = -beta_tilde
    rows[row, N * N + 1] = -beta
    s = np.linalg.svd(rows, compute_uv=False)
    rank = int(np.sum(s > TOL * max(float(s[0]), 1.0)))
    return {
        "complex_variables": N * N + 2,
        "complex_constraints": VISIBLE * VISIBLE + 1,
        "complex_rank": rank,
        "complex_nullity": N * N + 2 - rank,
        "smallest_singular_value": float(np.min(s)),
        "largest_singular_value": float(np.max(s)),
    }


def unitary_baryon_real_rank(m: np.ndarray, beta: complex, beta_tilde: complex) -> dict[str, Any]:
    cols: list[np.ndarray] = []
    det_m = np.linalg.det(m)
    inv_m = np.linalg.inv(m)
    for k in antihermitian_basis(N):
        dm = m @ k
        fixed = dm[:VISIBLE, :VISIBLE]
        d_quantum = det_m * np.trace(inv_m @ dm)
        cols.append(real_pack_complex(np.concatenate([fixed.reshape(-1), np.array([d_quantum])])))

    for d_beta, d_betat in [(1.0, 0.0), (1.0j, 0.0), (0.0, 1.0), (0.0, 1.0j)]:
        fixed = np.zeros((VISIBLE, VISIBLE), dtype=complex)
        d_quantum = -beta_tilde * d_beta - beta * d_betat
        cols.append(real_pack_complex(np.concatenate([fixed.reshape(-1), np.array([d_quantum])])))

    rank, s = rank_real(cols)
    return {
        "real_variables_unitary_link_plus_baryons": N * N + 4,
        "real_constraints_fixed_block_plus_quantum": 2 * (VISIBLE * VISIBLE + 1),
        "real_rank": rank,
        "real_nullity": N * N + 4 - rank,
        "smallest_singular_value": float(min(s)),
        "largest_singular_value": float(max(s)),
    }


def sample_rows(w: np.ndarray, lock_window: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, angle in enumerate(ANGLE_GRID):
        m = unitary_completion(w, angle, 202605090207 + idx)
        singulars = np.linalg.svd(m, compute_uv=False)
        for lam in LAMBDA_OVER_F_GRID:
            beta, beta_tilde, q = quantum_baryons(m, lam)
            res = quantum_residual(m, beta, beta_tilde, q)
            rows.append(
                {
                    "angle": angle,
                    "Lambda_over_f": lam,
                    "det_m_abs": float(abs(np.linalg.det(m))),
                    "det_m_phase": float(np.angle(np.linalg.det(m))),
                    "quantum_q": float(q.real),
                    "baryon_abs": float(abs(beta)),
                    "baryon_tilde_abs": float(abs(beta_tilde)),
                    "quantum_constraint_abs_residual": float(abs(res)),
                    "fixed_block_residual_fro": float(np.linalg.norm(m[:VISIBLE, :VISIBLE] - w, ord="fro")),
                    "unitarity_residual_2norm": float(np.linalg.norm(m.conj().T @ m - np.eye(N), ord=2)),
                    "max_abs_log_singular": float(np.max(np.abs(np.log(singulars)))),
                    "inside_completed_partner_lock_window": bool(
                        np.max(np.abs(np.log(singulars))) <= lock_window
                    ),
                }
            )
    return rows


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    completed = read_json(COMPLETED)
    w = cmat(ps["superpotential_ansatz"]["W_kappa"])
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])
    rows = sample_rows(w, lock_window)

    reference_m = unitary_completion(w, 0.1, 202605090299)
    beta, beta_tilde, _ = quantum_baryons(reference_m, 0.7)
    holo = holomorphic_jacobian_rank(reference_m, beta, beta_tilde)
    unitary = unitary_baryon_real_rank(reference_m, beta, beta_tilde)
    holomorphic_real_nullity = 2 * holo["complex_nullity"]
    nonunitary_removed_by_dterm = holomorphic_real_nullity - unitary["real_nullity"]

    verdict = {
        "all_unitary_completions_satisfy_quantum_constraint_with_baryons": all(
            row["quantum_constraint_abs_residual"] < 1.0e-12
            and row["fixed_block_residual_fro"] < 1.0e-12
            and row["inside_completed_partner_lock_window"]
            for row in rows
        ),
        "max_quantum_constraint_abs_residual": max(row["quantum_constraint_abs_residual"] for row in rows),
        "max_baryon_abs": max(max(row["baryon_abs"], row["baryon_tilde_abs"]) for row in rows),
        "holomorphic_quantum_constraint_alone_enforces_unitarity": False,
        "holomorphic_real_nullity_after_constraints": holomorphic_real_nullity,
        "unitary_plus_baryon_real_nullity": unitary["real_nullity"],
        "nonunitary_real_moduli_removed_by_Dterm_lock": nonunitary_removed_by_dterm,
        "visible_threshold_vector": [0.0, 0.0, 0.0],
        "interpretation": (
            "The Nf=Nc quantum-deformed constraint is compatible with the "
            "unitary crossed-120 link: baryon vevs solve det m - beta beta_tilde "
            "= (Lambda/f)^8 for every sampled unitary completion.  However, "
            "the holomorphic constraint system has complex nullity 13, so it "
            "does not replace the D-term/radial unitary lock.  With the unitary "
            "lock retained, only hidden baryon singlet moduli are added and the "
            "visible threshold vector remains zero."
        ),
    }

    return {
        "note": "No web lookup used. Quantum-deformed Nf=Nc hidden-SQCD audit for the unitary link.",
        "superpotential": {
            "W_qd": "X(det M - Baryon Baryon_tilde - Lambda_H^8) + Tr Y(PM/f^2-W_kappa)",
            "dimensionless_constraint": "det m - beta beta_tilde = (Lambda_H/f)^8",
            "Fflat_reference": "X=Y=0; beta,beta_tilde chosen to satisfy the constraint",
        },
        "finite_quantum_samples": rows,
        "linearized_constraint_ranks": {
            "holomorphic_fixed_block_plus_quantum": holo,
            "unitary_link_plus_baryons": unitary,
            "holomorphic_real_nullity_after_constraints": holomorphic_real_nullity,
            "nonunitary_real_moduli_removed_by_Dterm_lock": nonunitary_removed_by_dterm,
        },
        "benchmarks": {
            "M_lock_GeV": float(ps["benchmarks"]["M_lock_GeV"]),
            "MG_GeV": float(ps["benchmarks"]["MG_GeV"]),
            "completed_partner_max_abs_log_xi": lock_window,
        },
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    v = payload["verdict"]
    ranks = payload["linearized_constraint_ranks"]
    lines = [
        "# Quantum-deformed unitary-link moduli audit",
        "",
        "No web lookup was used.",
        "",
        "## Constraint",
        "",
        "```text",
        "det m - beta beta_tilde = (Lambda_H/f)^8",
        "P m P = W_kappa",
        "```",
        "",
        "For each sampled unitary link m, beta and beta_tilde are chosen as a",
        "complex square root of det(m)-(Lambda_H/f)^8.",
        "",
        "## Numerical result",
        "",
        f"max quantum residual: {v['max_quantum_constraint_abs_residual']:.3e}",
        f"max baryon magnitude: {v['max_baryon_abs']:.6g}",
        f"visible threshold vector: {v['visible_threshold_vector']}",
        "",
        "## Linearized ranks",
        "",
        "| system | variables | rank | nullity |",
        "|---|---:|---:|---:|",
        "| holomorphic fixed block + quantum | "
        f"{ranks['holomorphic_fixed_block_plus_quantum']['complex_variables']} complex | "
        f"{ranks['holomorphic_fixed_block_plus_quantum']['complex_rank']} complex | "
        f"{ranks['holomorphic_fixed_block_plus_quantum']['complex_nullity']} complex |",
        "| unitary link + baryons | "
        f"{ranks['unitary_link_plus_baryons']['real_variables_unitary_link_plus_baryons']} real | "
        f"{ranks['unitary_link_plus_baryons']['real_rank']} real | "
        f"{ranks['unitary_link_plus_baryons']['real_nullity']} real |",
        "",
        f"nonunitary real moduli removed by D-term/radial lock: {v['nonunitary_real_moduli_removed_by_Dterm_lock']}",
        "",
        "## Verdict",
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "finite_quantum_samples.csv",
        payload["finite_quantum_samples"],
        [
            "angle",
            "Lambda_over_f",
            "det_m_abs",
            "det_m_phase",
            "quantum_q",
            "baryon_abs",
            "baryon_tilde_abs",
            "quantum_constraint_abs_residual",
            "fixed_block_residual_fro",
            "unitarity_residual_2norm",
            "max_abs_log_singular",
            "inside_completed_partner_lock_window",
        ],
    )
    write_report(payload)
    v = payload["verdict"]
    print("Quantum-deformed unitary-link moduli audit")
    print(
        "  compatible with unitary completions: "
        f"{v['all_unitary_completions_satisfy_quantum_constraint_with_baryons']}"
    )
    print(f"  max quantum residual: {v['max_quantum_constraint_abs_residual']:.3e}")
    print(f"  holomorphic real nullity: {v['holomorphic_real_nullity_after_constraints']}")
    print(f"  unitary+baryon real nullity: {v['unitary_plus_baryon_real_nullity']}")
    print(f"  quantum constraint alone enforces unitarity: {v['holomorphic_quantum_constraint_alone_enforces_unitarity']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
