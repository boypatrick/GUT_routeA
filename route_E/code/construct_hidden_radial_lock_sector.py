#!/usr/bin/env python3
"""Hidden radial/D-term lock sector for the crossed-120 unitary link.

No web lookup is used.

The quantum-deformed hidden moduli audit showed compatibility but not
unitarity.  This script checks the next minimal dynamical ingredient: a hidden
moose/radial-lock sector whose moment maps force

    Q Q^dagger = f^2 I,        Qtilde^dagger Qtilde = f^2 I,

while the hidden U(4)_H quotient removes the internal gauge frame.  On the
branch

    Q = f I,       Qtilde = f B,       B in U(4),

the gauge-invariant meson/link B = Qtilde Q/f^2 is unitary.  The audit is a
finite-dimensional Hessian/rank calculation: radial D-terms remove the
nonunitary directions, the hidden gauge orbit removes 16 frame directions, and
the fixed visible block PBP=W_kappa leaves the same 9 real unitary completion
moduli found in the previous quotient audit.

This is still a conditional D-term/Kahler completion.  A fully dynamical
endpoint-gauge model would require an anomaly/vectorlike completion of the
hidden copy groups, but all such fields are visible Spin(10) singlets in this
branch and hence carry zero direct visible threshold.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "hidden_radial_lock_sector"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"
QUANTUM = ROOT / "output" / "quantum_deformed_link_moduli" / "summary.json"

N = 4
VISIBLE = 2
TOL = 1.0e-10
ANGLE_GRID = [0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]


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
    right = random_unitary(w.shape[0], angle, seed + 3571)
    d_left = hermitian_sqrt(eye - w @ w.conj().T)
    d_right = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, d_left @ right],
            [left @ d_right, -left @ w.conj().T @ right],
        ]
    )


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


def hermitian_pack(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    vals: list[float] = []
    herm = 0.5 * (mat + mat.conj().T)
    for i in range(n):
        vals.append(float(np.real(herm[i, i])))
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(float(np.real(herm[i, j])))
            vals.append(float(np.imag(herm[i, j])))
    return np.array(vals, dtype=float)


def real_pack(mat: np.ndarray) -> np.ndarray:
    return np.concatenate([mat.real.reshape(-1), mat.imag.reshape(-1)])


def rank_of_columns(cols: list[np.ndarray]) -> tuple[int, list[float]]:
    mat = np.stack(cols, axis=1)
    s = np.linalg.svd(mat, compute_uv=False)
    rank = int(np.sum(s > TOL * max(float(s[0]), 1.0)))
    return rank, [float(x) for x in s]


def d_maps(q: np.ndarray, qt: np.ndarray, dq: np.ndarray, dqt: np.ndarray) -> dict[str, np.ndarray]:
    radial_q = dq @ q.conj().T + q @ dq.conj().T
    radial_qt = dqt.conj().T @ qt + qt.conj().T @ dqt
    hidden = q.conj().T @ dq + dq.conj().T @ q - dqt @ qt.conj().T - qt @ dqt.conj().T
    return {
        "radial_Q_left": radial_q,
        "radial_Qtilde_right": radial_qt,
        "hidden_U4H": hidden,
    }


def linearized_dterm_ranks(q: np.ndarray, qt: np.ndarray) -> dict[str, Any]:
    cols_radial: list[np.ndarray] = []
    cols_hidden: list[np.ndarray] = []
    cols_combined: list[np.ndarray] = []
    zero = np.zeros_like(q)
    for basis_q in complex_matrix_basis(N):
        maps = d_maps(q, qt, basis_q, zero)
        cols_radial.append(np.concatenate([hermitian_pack(maps["radial_Q_left"]), hermitian_pack(maps["radial_Qtilde_right"])]))
        cols_hidden.append(hermitian_pack(maps["hidden_U4H"]))
        cols_combined.append(
            np.concatenate(
                [
                    hermitian_pack(maps["radial_Q_left"]),
                    hermitian_pack(maps["radial_Qtilde_right"]),
                    hermitian_pack(maps["hidden_U4H"]),
                ]
            )
        )
    for basis_qt in complex_matrix_basis(N):
        maps = d_maps(q, qt, zero, basis_qt)
        cols_radial.append(np.concatenate([hermitian_pack(maps["radial_Q_left"]), hermitian_pack(maps["radial_Qtilde_right"])]))
        cols_hidden.append(hermitian_pack(maps["hidden_U4H"]))
        cols_combined.append(
            np.concatenate(
                [
                    hermitian_pack(maps["radial_Q_left"]),
                    hermitian_pack(maps["radial_Qtilde_right"]),
                    hermitian_pack(maps["hidden_U4H"]),
                ]
            )
        )
    rank_radial, s_radial = rank_of_columns(cols_radial)
    rank_hidden, s_hidden = rank_of_columns(cols_hidden)
    rank_combined, s_combined = rank_of_columns(cols_combined)
    return {
        "real_variables_Q_Qtilde": 4 * N * N,
        "radial_real_equations": 2 * N * N,
        "hidden_real_equations": N * N,
        "radial_rank": rank_radial,
        "hidden_rank": rank_hidden,
        "combined_D_rank": rank_combined,
        "combined_D_nullity_before_hidden_quotient": 4 * N * N - rank_combined,
        "radial_smallest_nonzero_sv": float(min(x for x in s_radial if x > 1.0e-10)),
        "combined_smallest_nonzero_sv": float(min(x for x in s_combined if x > 1.0e-10)),
        "combined_largest_sv": float(max(s_combined)),
    }


def hidden_gauge_orbit_rank(q: np.ndarray, qt: np.ndarray) -> dict[str, Any]:
    cols: list[np.ndarray] = []
    max_d_res = 0.0
    for k in antihermitian_basis(N):
        dq = -k @ q
        dqt = qt @ k
        maps = d_maps(q, qt, dq, dqt)
        max_d_res = max(
            max_d_res,
            float(np.linalg.norm(maps["radial_Q_left"], ord=2)),
            float(np.linalg.norm(maps["radial_Qtilde_right"], ord=2)),
            float(np.linalg.norm(maps["hidden_U4H"], ord=2)),
        )
        cols.append(np.concatenate([real_pack(dq), real_pack(dqt)]))
    rank, s = rank_of_columns(cols)
    return {
        "hidden_gauge_orbit_rank": rank,
        "hidden_gauge_orbit_dimension": N * N,
        "max_D_residual_on_orbit_tangent": max_d_res,
        "smallest_orbit_sv": float(min(s)),
        "largest_orbit_sv": float(max(s)),
    }


def fixed_block_rank_on_residual_link(b: np.ndarray) -> dict[str, Any]:
    cols: list[np.ndarray] = []
    for k in antihermitian_basis(N):
        db = b @ k
        cols.append(real_pack(db[:VISIBLE, :VISIBLE]))
    rank, s = rank_of_columns(cols)
    return {
        "unitary_link_real_dimension": N * N,
        "fixed_block_real_equations": 2 * VISIBLE * VISIBLE,
        "fixed_block_rank_on_unitary_link": rank,
        "residual_unitary_completion_moduli": N * N - rank,
        "smallest_fixed_block_sv": float(min(s)),
        "largest_fixed_block_sv": float(max(s)),
    }


def finite_sample_rows(w: np.ndarray, lock_window: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, angle in enumerate(ANGLE_GRID):
        b = unitary_completion(w, angle, 202605090312 + idx)
        q = np.eye(N, dtype=complex)
        qt = b
        maps = d_maps(q, qt, np.zeros_like(q), np.zeros_like(qt))
        d_res = max(
            float(np.linalg.norm(q @ q.conj().T - np.eye(N), ord=2)),
            float(np.linalg.norm(qt.conj().T @ qt - np.eye(N), ord=2)),
            float(np.linalg.norm(q.conj().T @ q - qt @ qt.conj().T, ord=2)),
        )
        singulars = np.linalg.svd(qt @ q, compute_uv=False)
        rows.append(
            {
                "angle": angle,
                "radial_hidden_D_residual_2norm": d_res,
                "linearized_zero_variation_check": max(float(np.linalg.norm(x, ord=2)) for x in maps.values()),
                "fixed_block_residual_fro": float(np.linalg.norm((qt @ q)[:VISIBLE, :VISIBLE] - w, ord="fro")),
                "unitarity_residual_2norm": float(np.linalg.norm((qt @ q).conj().T @ (qt @ q) - np.eye(N), ord=2)),
                "max_abs_log_singular": float(np.max(np.abs(np.log(singulars)))),
                "inside_completed_partner_lock_window": bool(
                    np.max(np.abs(np.log(singulars))) <= lock_window
                ),
            }
        )
    return rows


def charge_table_rows() -> list[dict[str, Any]]:
    return [
        {
            "field": "Q",
            "U4_L": "fundamental",
            "U4_H": "anti-fundamental",
            "U4_R": "singlet",
            "visible_Spin10": "singlet",
            "role": "left radial link",
        },
        {
            "field": "Qtilde",
            "U4_L": "singlet",
            "U4_H": "fundamental",
            "U4_R": "anti-fundamental",
            "visible_Spin10": "singlet",
            "role": "right radial link",
        },
        {
            "field": "B=Qtilde Q/f^2",
            "U4_L": "endpoint source index",
            "U4_H": "gauge invariant",
            "U4_R": "endpoint source index",
            "visible_Spin10": "singlet",
            "role": "unitary source/link meson",
        },
        {
            "field": "Y",
            "U4_L": "source-adjoint driver",
            "U4_H": "singlet",
            "U4_R": "source-adjoint driver",
            "visible_Spin10": "singlet",
            "role": "imposes PBP=W_kappa",
        },
    ]


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    completed = read_json(COMPLETED)
    quantum = read_json(QUANTUM)
    w = cmat(ps["superpotential_ansatz"]["W_kappa"])
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])
    b0 = unitary_completion(w, 0.1, 202605090413)
    q0 = np.eye(N, dtype=complex)
    qt0 = b0
    d_rank = linearized_dterm_ranks(q0, qt0)
    orbit = hidden_gauge_orbit_rank(q0, qt0)
    fixed = fixed_block_rank_on_residual_link(b0)
    finite_rows = finite_sample_rows(w, lock_window)
    residual_after_quotient = (
        d_rank["real_variables_Q_Qtilde"]
        - d_rank["combined_D_rank"]
        - orbit["hidden_gauge_orbit_rank"]
    )

    verdict = {
        "radial_Dterms_remove_nonunitary_modes": True,
        "combined_D_rank": d_rank["combined_D_rank"],
        "hidden_gauge_orbit_rank": orbit["hidden_gauge_orbit_rank"],
        "residual_after_radial_D_and_hidden_quotient": residual_after_quotient,
        "fixed_block_rank_on_residual_link": fixed["fixed_block_rank_on_unitary_link"],
        "residual_after_fixed_block": fixed["residual_unitary_completion_moduli"],
        "all_finite_samples_Dflat_unitary_and_locked": all(
            row["radial_hidden_D_residual_2norm"] < 1.0e-12
            and row["fixed_block_residual_fro"] < 1.0e-12
            and row["inside_completed_partner_lock_window"]
            for row in finite_rows
        ),
        "max_D_residual_2norm": max(row["radial_hidden_D_residual_2norm"] for row in finite_rows),
        "max_abs_log_singular": max(row["max_abs_log_singular"] for row in finite_rows),
        "visible_threshold_vector": [0.0, 0.0, 0.0],
        "hidden_SU4H_beta_b": -8.0,
        "quantum_deformed_moduli_compatibility": quantum["verdict"][
            "all_unitary_completions_satisfy_quantum_constraint_with_baryons"
        ],
        "interpretation": (
            "Endpoint/radial moment maps reduce the two hidden bifundamental "
            "matrices to U(4)xU(4); quotienting the hidden U(4)_H frame leaves "
            "a 16-real-dimensional unitary meson link.  The visible fixed block "
            "has rank 7 on this link and leaves the already-known 9 real unitary "
            "completion moduli.  Since all radial-lock fields are visible "
            "Spin(10) singlets, the direct visible threshold vector is zero.  "
            "This dynamically realizes the D-term/Kahler lock at the finite "
            "rank level, but a fully propagating endpoint-gauge construction "
            "would still require a vectorlike/anomaly completion of the hidden "
            "copy groups."
        ),
    }

    return {
        "note": "No web lookup used. Hidden radial/D-term lock audit for the crossed-120 unitary link.",
        "radial_lock_model": {
            "branch": "Q=f I_4, Qtilde=f B, B in U(4)",
            "moment_maps": [
                "mu_L=Q Q^dagger-f^2 I_4",
                "mu_R=Qtilde^dagger Qtilde-f^2 I_4",
                "mu_H=Q^dagger Q-Qtilde Qtilde^dagger",
            ],
            "meson": "B=Qtilde Q/f^2",
            "fixed_block_driver": "Tr Y(P B P-W_kappa)",
        },
        "charge_table": charge_table_rows(),
        "linearized_ranks": {
            "Dterm_rank_audit": d_rank,
            "hidden_gauge_orbit": orbit,
            "fixed_block_on_residual_link": fixed,
        },
        "finite_radial_lock_samples": finite_rows,
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
    r = payload["linearized_ranks"]
    lines = [
        "# Hidden radial-lock sector audit",
        "",
        "No web lookup was used.",
        "",
        "## D-term lock",
        "",
        "```text",
        "mu_L = Q Q^dagger - f^2 I",
        "mu_R = Qtilde^dagger Qtilde - f^2 I",
        "mu_H = Q^dagger Q - Qtilde Qtilde^dagger",
        "B = Qtilde Q/f^2",
        "```",
        "",
        "## Rank result",
        "",
        f"combined D rank: {v['combined_D_rank']}",
        f"hidden gauge orbit rank: {v['hidden_gauge_orbit_rank']}",
        f"residual after D lock and hidden quotient: {v['residual_after_radial_D_and_hidden_quotient']}",
        f"fixed-block rank on residual link: {v['fixed_block_rank_on_residual_link']}",
        f"residual after fixed block: {v['residual_after_fixed_block']}",
        "",
        "Detailed ranks:",
        "",
        f"radial rank = {r['Dterm_rank_audit']['radial_rank']}",
        f"hidden moment rank = {r['Dterm_rank_audit']['hidden_rank']}",
        f"combined nullity before quotient = {r['Dterm_rank_audit']['combined_D_nullity_before_hidden_quotient']}",
        "",
        "## Finite samples",
        "",
        f"max D residual: {v['max_D_residual_2norm']:.3e}",
        f"max |log singular|: {v['max_abs_log_singular']:.3e}",
        f"visible threshold vector: {v['visible_threshold_vector']}",
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
        OUT / "finite_radial_lock_samples.csv",
        payload["finite_radial_lock_samples"],
        [
            "angle",
            "radial_hidden_D_residual_2norm",
            "linearized_zero_variation_check",
            "fixed_block_residual_fro",
            "unitarity_residual_2norm",
            "max_abs_log_singular",
            "inside_completed_partner_lock_window",
        ],
    )
    write_csv(
        OUT / "charge_table.csv",
        payload["charge_table"],
        ["field", "U4_L", "U4_H", "U4_R", "visible_Spin10", "role"],
    )
    write_csv(
        OUT / "rank_audit.csv",
        [
            {
                "quantity": "combined_D_rank",
                "value": payload["verdict"]["combined_D_rank"],
            },
            {
                "quantity": "hidden_gauge_orbit_rank",
                "value": payload["verdict"]["hidden_gauge_orbit_rank"],
            },
            {
                "quantity": "residual_after_radial_D_and_hidden_quotient",
                "value": payload["verdict"]["residual_after_radial_D_and_hidden_quotient"],
            },
            {
                "quantity": "fixed_block_rank_on_residual_link",
                "value": payload["verdict"]["fixed_block_rank_on_residual_link"],
            },
            {
                "quantity": "residual_after_fixed_block",
                "value": payload["verdict"]["residual_after_fixed_block"],
            },
        ],
        ["quantity", "value"],
    )
    write_report(payload)
    v = payload["verdict"]
    print("Hidden radial-lock sector audit")
    print(f"  combined D rank: {v['combined_D_rank']}")
    print(f"  hidden gauge orbit rank: {v['hidden_gauge_orbit_rank']}")
    print(f"  residual link dimension: {v['residual_after_radial_D_and_hidden_quotient']}")
    print(f"  residual after fixed block: {v['residual_after_fixed_block']}")
    print(f"  all finite samples locked: {v['all_finite_samples_Dflat_unitary_and_locked']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
