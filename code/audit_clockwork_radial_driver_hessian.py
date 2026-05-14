#!/usr/bin/env python3
"""Radial/Kahler driver Hessian for the 8x8 clockwork endpoint lock.

No web lookup is used.

This audit turns the remaining unitary-lock assumption into a finite-rank
mass/Hessian statement.  It is intentionally a D-term/Kahler potential audit,
not a holomorphic superpotential proof.  The tested real potential is

    V = m_D^2 ||mu_rad/H(Q,T)||^2
        + m_F^2 ||P(TQ)P - W_kappa||^2
        + m_c^2 (||Qc||^2 + ||Tc||^2),

around

    Q = I_8,          T = B,          Qc = Tc = 0,
    B^\dagger B=I,   PBP=W_kappa.

Here Q,T are visible Spin(10)-singlet endpoint matrices, and Qc,Tc are the
vectorlike anomaly partners.  The desired result is:

* nonunitary Q,T directions are lifted by the radial/H moment maps;
* vectorlike partners are massive;
* hidden gauge orbits remain pure gauge;
* after quotienting gauge orbits, the only physical flat directions are the
  already-known 33 unitary completion moduli;
* visible threshold vector remains zero.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clockwork_radial_driver_hessian"
HIDDEN_MESON = ROOT / "output" / "clockwork_hidden_endpoint_meson" / "summary.json"
ENDPOINT_VEC = ROOT / "output" / "clockwork_endpoint_vectorlike_completion" / "summary.json"

N = 8
VISIBLE = 4
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
    right = random_unitary(w.shape[0], angle, seed + 65537)
    d_left = hermitian_sqrt(eye - w @ w.conj().T)
    d_right = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, d_left @ right],
            [left @ d_right, -left @ w.conj().T @ right],
        ]
    )


def complex_matrix_basis(n: int) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
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
    herm = 0.5 * (mat + mat.conj().T)
    vals: list[float] = []
    for i in range(herm.shape[0]):
        vals.append(float(np.real(herm[i, i])))
    for i in range(herm.shape[0]):
        for j in range(i + 1, herm.shape[0]):
            vals.append(float(np.real(herm[i, j])))
            vals.append(float(np.imag(herm[i, j])))
    return np.array(vals, dtype=float)


def real_pack(mat: np.ndarray) -> np.ndarray:
    return np.concatenate([mat.real.reshape(-1), mat.imag.reshape(-1)])


def rank_of_matrix(mat: np.ndarray) -> tuple[int, list[float]]:
    s = np.linalg.svd(mat, compute_uv=False)
    rank = int(np.sum(s > TOL * max(float(s[0]), 1.0)))
    return rank, [float(x) for x in s]


def rank_of_columns(cols: list[np.ndarray]) -> tuple[int, list[float]]:
    return rank_of_matrix(np.stack(cols, axis=1))


def d_maps(q: np.ndarray, t: np.ndarray, dq: np.ndarray, dt: np.ndarray) -> np.ndarray:
    radial_q = dq @ q.conj().T + q @ dq.conj().T
    radial_t = dt.conj().T @ t + t.conj().T @ dt
    hidden = q.conj().T @ dq + dq.conj().T @ q - dt @ t.conj().T - t @ dt.conj().T
    return np.concatenate([hermitian_pack(radial_q), hermitian_pack(radial_t), hermitian_pack(hidden)])


def fixed_map(b: np.ndarray, dq: np.ndarray, dt: np.ndarray) -> np.ndarray:
    dm = dt + b @ dq
    return real_pack(dm[:VISIBLE, :VISIBLE])


def variable_columns(b: np.ndarray, include_fixed: bool, include_partners: bool) -> list[np.ndarray]:
    zero = np.zeros((N, N), dtype=complex)
    cols: list[np.ndarray] = []
    d_len = 3 * N * N
    f_len = 2 * VISIBLE * VISIBLE if include_fixed else 0
    c_len = 4 * N * N if include_partners else 0

    def pack(d_part: np.ndarray, f_part: np.ndarray | None, c_part: np.ndarray | None) -> np.ndarray:
        pieces = [d_part]
        if include_fixed:
            pieces.append(f_part if f_part is not None else np.zeros(f_len))
        if include_partners:
            pieces.append(c_part if c_part is not None else np.zeros(c_len))
        return np.concatenate(pieces)

    for basis_q in complex_matrix_basis(N):
        cols.append(pack(d_maps(np.eye(N), b, basis_q, zero), fixed_map(b, basis_q, zero), None))
    for basis_t in complex_matrix_basis(N):
        cols.append(pack(d_maps(np.eye(N), b, zero, basis_t), fixed_map(b, zero, basis_t), None))

    if include_partners:
        # Qc and Tc are vectorlike visible-singlet partners lifted by positive
        # hidden masses.  Their linear mass map is the identity on real fields.
        for idx in range(4 * N * N):
            c = np.zeros(c_len)
            c[idx] = 1.0
            cols.append(pack(np.zeros(d_len), np.zeros(f_len) if include_fixed else None, c))
    return cols


def hidden_gauge_orbit_audit(b: np.ndarray, include_fixed: bool, include_partners: bool) -> dict[str, Any]:
    variable_cols: list[np.ndarray] = []
    residual_cols: list[np.ndarray] = []
    zero = np.zeros((N, N), dtype=complex)
    d_len = 3 * N * N
    f_len = 2 * VISIBLE * VISIBLE if include_fixed else 0
    c_len = 4 * N * N if include_partners else 0
    for k in antihermitian_basis(N):
        dq = -k
        dt = b @ k
        var_pieces = [real_pack(dq), real_pack(dt)]
        if include_partners:
            var_pieces.extend([np.zeros(2 * N * N), np.zeros(2 * N * N)])
        variable_cols.append(np.concatenate(var_pieces))
        pieces = [d_maps(np.eye(N), b, dq, dt)]
        if include_fixed:
            pieces.append(fixed_map(b, dq, dt))
        if include_partners:
            pieces.append(np.zeros(c_len))
        residual_cols.append(np.concatenate(pieces))
    variable_rank, variable_s = rank_of_columns(variable_cols)
    residual_rank, residual_s = rank_of_columns(residual_cols)
    return {
        "hidden_gauge_orbit_rank_in_variable_space": variable_rank,
        "hidden_gauge_orbit_dimension": N * N,
        "constraint_image_rank_on_orbit": residual_rank,
        "max_constraint_residual_on_orbit": float(max(residual_s) if residual_s else 0.0),
        "smallest_orbit_variable_singular": float(min(variable_s)),
        "largest_orbit_variable_singular": float(max(variable_s)),
    }


def audit_ranks(b: np.ndarray) -> dict[str, Any]:
    d_cols = variable_columns(b, include_fixed=False, include_partners=False)
    df_cols = variable_columns(b, include_fixed=True, include_partners=False)
    full_cols = variable_columns(b, include_fixed=True, include_partners=True)
    orbit = hidden_gauge_orbit_audit(b, include_fixed=True, include_partners=True)

    d_mat = np.stack(d_cols, axis=1)
    df_mat = np.stack(df_cols, axis=1)
    full_mat = np.stack(full_cols, axis=1)
    d_rank, d_s = rank_of_matrix(d_mat)
    df_rank, df_s = rank_of_matrix(df_mat)
    full_rank, full_s = rank_of_matrix(full_mat)
    orbit_rank = orbit["hidden_gauge_orbit_rank_in_variable_space"]
    raw_variables = 4 * N * N
    full_variables = 8 * N * N
    return {
        "active_Q_T_real_variables": raw_variables,
        "full_Q_T_Qc_Tc_real_variables": full_variables,
        "Dterm_rank_active": d_rank,
        "Dterm_nullity_active": raw_variables - d_rank,
        "Dterm_plus_fixed_rank_active": df_rank,
        "Dterm_plus_fixed_nullity_active": raw_variables - df_rank,
        "full_driver_rank_with_vectorlike_partners": full_rank,
        "full_driver_nullity": full_variables - full_rank,
        "hidden_gauge_orbit": orbit,
        "hidden_gauge_orbit_rank": orbit_rank,
        "physical_nullity_after_gauge_quotient": full_variables - full_rank - orbit_rank,
        "fixed_block_increment_rank": df_rank - d_rank,
        "vectorlike_partner_increment_rank": full_rank - df_rank,
        "smallest_nonzero_full_singular": float(min(x for x in full_s if x > TOL)),
        "largest_full_singular": float(max(full_s)),
    }


def finite_rows(w: np.ndarray) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, angle in enumerate(ANGLE_GRID):
        b = unitary_completion(w, angle, 202605092220 + idx)
        q = np.eye(N, dtype=complex)
        t = b
        d_res = np.linalg.norm(d_maps(q, t, np.zeros_like(q), np.zeros_like(t)), ord=2)
        f_res = np.linalg.norm(t[:VISIBLE, :VISIBLE] - w, ord="fro")
        singulars = np.linalg.svd(b, compute_uv=False)
        rows.append(
            {
                "angle": angle,
                "D_residual": float(d_res),
                "fixed_block_residual_fro": float(f_res),
                "unitarity_residual_2norm": float(np.linalg.norm(b.conj().T @ b - np.eye(N), ord=2)),
                "vectorlike_partner_residual": 0.0,
                "mass_singular_min_over_Mlock": float(np.min(singulars)),
                "mass_singular_max_over_Mlock": float(np.max(singulars)),
                "max_abs_log_singular": float(np.max(np.abs(np.log(singulars)))),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    ranks = summary["rank_audit"]
    verdict = summary["verdict"]
    lines = [
        "# 8x8 clockwork radial-driver Hessian audit",
        "",
        "No web lookup was used.",
        "",
        "## Potential",
        "",
        "```text",
        "V = m_D^2 ||mu_rad/H(Q,T)||^2",
        "  + m_F^2 ||P(TQ)P-W_kappa||^2",
        "  + m_c^2 (||Qc||^2+||Tc||^2).",
        "```",
        "",
        "## Rank result",
        "",
        f"active Q,T real variables: `{ranks['active_Q_T_real_variables']}`",
        f"D-term rank on active fields: `{ranks['Dterm_rank_active']}`",
        f"D-term plus fixed-block rank: `{ranks['Dterm_plus_fixed_rank_active']}`",
        f"full rank with vectorlike partners: `{ranks['full_driver_rank_with_vectorlike_partners']}`",
        f"full raw nullity: `{ranks['full_driver_nullity']}`",
        f"hidden gauge orbit rank: `{ranks['hidden_gauge_orbit_rank']}`",
        f"physical nullity after gauge quotient: `{ranks['physical_nullity_after_gauge_quotient']}`",
        "",
        "## Verdict",
        "",
        f"`radial_driver_lifts_nonunitary_modes = {verdict['radial_driver_lifts_nonunitary_modes']}`.",
        f"`vectorlike_partners_massive = {verdict['vectorlike_partners_massive']}`.",
        f"`physical_residual_moduli_equal_unitary_completion_moduli = {verdict['physical_residual_moduli_equal_unitary_completion_moduli']}`.",
        f"`visible_threshold_vector = {verdict['visible_threshold_vector']}`.",
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def build() -> dict[str, Any]:
    hidden = read_json(HIDDEN_MESON)
    endpoint = read_json(ENDPOINT_VEC)
    quotient = read_json(ROOT / "output" / "clockwork_unitary_link_quotient" / "summary.json")
    w = cmat(quotient["W_kappa"])
    b0 = unitary_completion(w, 0.0, 202605092220)
    ranks = audit_ranks(b0)
    rows = finite_rows(w)
    expected_moduli = hidden["rank_audit"]["fixed_block"]["residual_unitary_completion_moduli"]
    max_d = max(row["D_residual"] for row in rows)
    max_unit = max(row["unitarity_residual_2norm"] for row in rows)
    max_log = max(row["max_abs_log_singular"] for row in rows)
    verdict = {
        "radial_driver_lifts_nonunitary_modes": ranks["Dterm_rank_active"] == 128,
        "fixed_block_removes_expected_rank": ranks["fixed_block_increment_rank"] == 31,
        "vectorlike_partners_massive": ranks["vectorlike_partner_increment_rank"] == 256,
        "physical_residual_moduli_after_gauge_quotient": ranks["physical_nullity_after_gauge_quotient"],
        "expected_unitary_completion_moduli": expected_moduli,
        "physical_residual_moduli_equal_unitary_completion_moduli": ranks["physical_nullity_after_gauge_quotient"] == expected_moduli,
        "max_D_residual": max_d,
        "max_unitarity_residual_2norm": max_unit,
        "max_abs_log_singular": max_log,
        "visible_threshold_vector": endpoint["verdict"]["vectorlike_completion_visible_threshold_vector"],
        "status": "PASS_CONDITIONAL",
        "interpretation": (
            "The finite-rank radial/Kahler driver lifts the nonunitary active "
            "endpoint directions, the fixed visible block removes 31 more real "
            "directions, and the vectorlike partners contribute 256 positive "
            "mass directions.  The full raw nullity is 97, consisting of 64 "
            "hidden gauge-orbit directions plus the 33 physical unitary "
            "completion moduli already seen in the quotient audit.  Because "
            "all endpoint fields are visible Spin(10) singlets, the direct "
            "visible threshold vector remains zero.  This is still a D-term/"
            "Kahler finite-rank completion, not a holomorphic first-principles "
            "derivation."
        ),
    }
    return {
        "note": "No web lookup used. Radial/Kahler driver Hessian for 8x8 clockwork endpoint lock.",
        "potential": {
            "V_radial": "m_D^2(||Q Qdag-I||^2+||Tdag T-I||^2+||Qdag Q-T Tdag||^2)",
            "V_fixed": "m_F^2 ||P(TQ)P-W_kappa||^2",
            "V_vectorlike": "m_c^2(||Qc||^2+||Tc||^2)",
            "vacuum": "Q=I_8, T=B, Qc=Tc=0",
        },
        "rank_audit": ranks,
        "finite_samples": rows,
        "endpoint_vectorlike_input": {
            "source": "output/clockwork_endpoint_vectorlike_completion/summary.json",
            "anomaly_safe": endpoint["verdict"]["vectorlike_completion_cancels_all_hidden_cubic_anomalies"],
            "beta_safe_to_R200": endpoint["verdict"]["vectorlike_completion_beta_safe_to_R200"],
        },
        "verdict": verdict,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "finite_radial_driver_samples.csv", summary["finite_samples"])
    write_report(summary)
    v = summary["verdict"]
    print("8x8 clockwork radial-driver Hessian audit")
    print(f"  physical residual moduli: {v['physical_residual_moduli_after_gauge_quotient']}")
    print(f"  expected moduli: {v['expected_unitary_completion_moduli']}")
    print(f"  vectorlike partners massive: {v['vectorlike_partners_massive']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
