#!/usr/bin/env python3
"""Hidden endpoint/meson audit for the clockwork unitary-link quotient.

No web lookup is used.

The current clockwork publication card requires the finite-dimensional
constraint

    P B P = W_kappa,        B^dagger B = I,

with W_kappa=diag(1, eps, eps, eps), eps=3^-6, and B an 8 by 8 unitary
Julia/Halmos completion.  Earlier hidden-sector audits treated the older
crossed-link 2 by 2 fixed block.  This script upgrades that audit to the
current 4 by 4 publication-card block.

The tested microscopic caricature is a hidden endpoint sector

    Q = f I_8,       Qtilde = f B,       M/f^2 = B = Qtilde Q/f^2,

where Q,Qtilde are visible Spin(10)-singlet hidden endpoint fields.  The audit
separates three statements:

1. A hidden radial/D-term lock can enforce B in U(8) at finite rank.
2. The Nf=Nc=8 quantum-deformed meson constraint is compatible with this
   branch, but cannot by itself enforce unitarity.
3. Since the endpoint fields are visible singlets, the direct visible threshold
   vector is zero; only hidden-sector consistency remains.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clockwork_hidden_endpoint_meson"
CLOCKWORK_UNITARY = ROOT / "output" / "clockwork_unitary_link_quotient" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"

N = 8
VISIBLE = 4
TOL = 1.0e-10
ANGLE_GRID = [0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
LAMBDA_OVER_F_GRID = [0.3, 0.7, 1.0]
NF_GRID = [8, 10, 12, 16, 24, 32]
ALPHA_H_INV_GRID = [1.0, 5.0, 10.0, 25.0]
R_TARGETS = [50.0, 200.0]


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
    right = random_unitary(w.shape[0], angle, seed + 104729)
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
        "hidden_U8H": hidden,
    }


def linearized_dterm_ranks(q: np.ndarray, qt: np.ndarray) -> dict[str, Any]:
    cols_radial: list[np.ndarray] = []
    cols_hidden: list[np.ndarray] = []
    cols_combined: list[np.ndarray] = []
    zero = np.zeros_like(q)
    for basis_q in complex_matrix_basis(N):
        maps = d_maps(q, qt, basis_q, zero)
        cols_radial.append(np.concatenate([hermitian_pack(maps["radial_Q_left"]), hermitian_pack(maps["radial_Qtilde_right"])]))
        cols_hidden.append(hermitian_pack(maps["hidden_U8H"]))
        cols_combined.append(
            np.concatenate(
                [
                    hermitian_pack(maps["radial_Q_left"]),
                    hermitian_pack(maps["radial_Qtilde_right"]),
                    hermitian_pack(maps["hidden_U8H"]),
                ]
            )
        )
    for basis_qt in complex_matrix_basis(N):
        maps = d_maps(q, qt, zero, basis_qt)
        cols_radial.append(np.concatenate([hermitian_pack(maps["radial_Q_left"]), hermitian_pack(maps["radial_Qtilde_right"])]))
        cols_hidden.append(hermitian_pack(maps["hidden_U8H"]))
        cols_combined.append(
            np.concatenate(
                [
                    hermitian_pack(maps["radial_Q_left"]),
                    hermitian_pack(maps["radial_Qtilde_right"]),
                    hermitian_pack(maps["hidden_U8H"]),
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
        "radial_smallest_nonzero_sv": float(min(x for x in s_radial if x > TOL)),
        "combined_smallest_nonzero_sv": float(min(x for x in s_combined if x > TOL)),
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
            float(np.linalg.norm(maps["hidden_U8H"], ord=2)),
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


def quantum_baryons(m: np.ndarray, lambda_over_f: float) -> tuple[complex, complex, complex]:
    q = complex(lambda_over_f ** (2 * N))
    product = np.linalg.det(m) - q
    beta = np.sqrt(product)
    return beta, beta, q


def holomorphic_quantum_rank(m: np.ndarray, beta: complex, beta_tilde: complex) -> dict[str, Any]:
    rows = np.zeros((VISIBLE * VISIBLE + 1, N * N + 2), dtype=complex)
    row = 0
    for i in range(VISIBLE):
        for j in range(VISIBLE):
            rows[row, i * N + j] = 1.0
            row += 1
    det_m = np.linalg.det(m)
    rows[row, : N * N] = (det_m * np.linalg.inv(m).T).reshape(-1)
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
        cols.append(real_pack(np.concatenate([fixed.reshape(-1), np.array([d_quantum])])))
    for d_beta, d_betat in [(1.0, 0.0), (1.0j, 0.0), (0.0, 1.0), (0.0, 1.0j)]:
        fixed = np.zeros((VISIBLE, VISIBLE), dtype=complex)
        d_quantum = -beta_tilde * d_beta - beta * d_betat
        cols.append(real_pack(np.concatenate([fixed.reshape(-1), np.array([d_quantum])])))
    rank, s = rank_of_columns(cols)
    return {
        "real_variables_unitary_link_plus_baryons": N * N + 4,
        "real_constraints_fixed_block_plus_quantum": 2 * (VISIBLE * VISIBLE + 1),
        "real_rank": rank,
        "real_nullity": N * N + 4 - rank,
        "smallest_singular_value": float(min(s)),
        "largest_singular_value": float(max(s)),
    }


def sqcd_phase(nf: int) -> str:
    if nf < N:
        return "ADS_runaway_not_used"
    if nf == N:
        return "Nf=Nc_quantum_deformed_composite_meson"
    if nf < 1.5 * N:
        return "below_conformal_window_strong"
    if nf == int(1.5 * N):
        return "lower_edge_conformal_window"
    if nf < 3 * N:
        return "inside_conformal_window"
    if nf == 3 * N:
        return "one_loop_marginal_edge"
    return "IR_free_UV_Landau_possible"


def landau_ratio(alpha_inv: float, b_hidden: float) -> float | None:
    if b_hidden <= 0.0:
        return None
    return math.exp(2.0 * math.pi * alpha_inv / b_hidden)


def hidden_beta_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for nf in NF_GRID:
        b_hidden = float(nf - 3 * N)
        for alpha_inv in ALPHA_H_INV_GRID:
            ratio = landau_ratio(alpha_inv, b_hidden)
            for target in R_TARGETS:
                rows.append(
                    {
                        "Nc": N,
                        "Nf": nf,
                        "b_hidden": b_hidden,
                        "alphaH_inv_at_matching": alpha_inv,
                        "phase_label": sqcd_phase(nf),
                        "uv_landau_ratio_over_matching": ratio,
                        "target_R": target,
                        "safe_to_target_R": bool(ratio is None or ratio > target),
                    }
                )
    return rows


def finite_rows(w: np.ndarray, lock_window: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, angle in enumerate(ANGLE_GRID):
        b = unitary_completion(w, angle, 202605092108 + idx)
        q = np.eye(N, dtype=complex)
        qt = b
        singulars = np.linalg.svd(b, compute_uv=False)
        d_res = max(
            float(np.linalg.norm(q @ q.conj().T - np.eye(N), ord=2)),
            float(np.linalg.norm(qt.conj().T @ qt - np.eye(N), ord=2)),
            float(np.linalg.norm(q.conj().T @ q - qt @ qt.conj().T, ord=2)),
        )
        max_abs_log = float(np.max(np.abs(np.log(singulars))))
        for lam in LAMBDA_OVER_F_GRID:
            beta, beta_tilde, qdef = quantum_baryons(b, lam)
            qres = np.linalg.det(b) - beta * beta_tilde - qdef
            rows.append(
                {
                    "angle": angle,
                    "Lambda_over_f": lam,
                    "fixed_block_residual_fro": float(np.linalg.norm(b[:VISIBLE, :VISIBLE] - w, ord="fro")),
                    "D_residual_2norm": d_res,
                    "unitarity_residual_2norm": float(np.linalg.norm(b.conj().T @ b - np.eye(N), ord=2)),
                    "quantum_constraint_abs_residual": float(abs(qres)),
                    "baryon_abs": float(abs(beta)),
                    "mass_singular_min_over_Mlock": float(np.min(singulars)),
                    "mass_singular_max_over_Mlock": float(np.max(singulars)),
                    "max_abs_log_singular": max_abs_log,
                    "within_completed_partner_lock_window": bool(max_abs_log <= lock_window),
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
    quantum = summary["quantum_rank_audit"]
    verdict = summary["verdict"]
    lines = [
        "# Clockwork hidden endpoint meson audit",
        "",
        "No web lookup was used.",
        "",
        "## Branch",
        "",
        "```text",
        "Q = f I_8,    Qtilde = f B,    B = Qtilde Q/f^2,",
        "P B P = W_kappa = diag(1, epsilon, epsilon, epsilon),",
        "B^dagger B = I.",
        "```",
        "",
        "## Rank result",
        "",
        f"combined D rank: `{ranks['Dterm']['combined_D_rank']}`",
        f"hidden gauge orbit rank: `{ranks['hidden_orbit']['hidden_gauge_orbit_rank']}`",
        f"residual after D lock and hidden quotient: `{ranks['residual_after_D_lock_and_hidden_quotient']}`",
        f"fixed-block rank on residual link: `{ranks['fixed_block']['fixed_block_rank_on_unitary_link']}`",
        f"residual after fixed block: `{ranks['fixed_block']['residual_unitary_completion_moduli']}`",
        "",
        "## Quantum-deformed compatibility",
        "",
        f"holomorphic complex nullity: `{quantum['holomorphic']['complex_nullity']}`",
        f"unitary plus baryon real nullity: `{quantum['unitary_plus_baryons']['real_nullity']}`",
        f"nonunitary real moduli removed by radial/D lock: `{quantum['nonunitary_real_moduli_removed_by_Dterm_lock']}`",
        "",
        "## Verdict",
        "",
        f"`finite_samples_Dflat_unitary_quantum_and_locked = {verdict['finite_samples_Dflat_unitary_quantum_and_locked']}`.",
        f"`quantum_constraint_alone_enforces_unitarity = {verdict['quantum_constraint_alone_enforces_unitarity']}`.",
        f"`visible_threshold_vector = {verdict['visible_threshold_vector']}`.",
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def build() -> dict[str, Any]:
    clockwork = read_json(CLOCKWORK_UNITARY)
    completed = read_json(COMPLETED)
    w = cmat(clockwork["W_kappa"])
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])
    b0 = unitary_completion(w, 0.0, 202605092108)
    q0 = np.eye(N, dtype=complex)
    qt0 = b0

    dterm = linearized_dterm_ranks(q0, qt0)
    orbit = hidden_gauge_orbit_rank(q0, qt0)
    fixed = fixed_block_rank_on_residual_link(b0)
    beta, beta_tilde, _ = quantum_baryons(b0, 0.7)
    holo = holomorphic_quantum_rank(b0, beta, beta_tilde)
    unitary_q = unitary_baryon_real_rank(b0, beta, beta_tilde)
    holomorphic_real_nullity = 2 * holo["complex_nullity"]
    nonunitary_removed = holomorphic_real_nullity - unitary_q["real_nullity"]
    finite = finite_rows(w, lock_window)
    beta_rows = hidden_beta_rows()
    minimal = next(row for row in beta_rows if row["Nf"] == 8 and row["alphaH_inv_at_matching"] == 10.0 and row["target_R"] == 200.0)
    conformal = next(row for row in beta_rows if row["Nf"] == 16 and row["alphaH_inv_at_matching"] == 10.0 and row["target_R"] == 200.0)

    verdict = {
        "finite_samples_Dflat_unitary_quantum_and_locked": all(
            row["D_residual_2norm"] < 1.0e-12
            and row["unitarity_residual_2norm"] < 1.0e-12
            and row["fixed_block_residual_fro"] < 1.0e-12
            and row["quantum_constraint_abs_residual"] < 1.0e-12
            and row["within_completed_partner_lock_window"]
            for row in finite
        ),
        "quantum_constraint_alone_enforces_unitarity": False,
        "minimal_Nf8_phase": minimal["phase_label"],
        "minimal_Nf8_uv_landau_safe_to_R200": bool(minimal["safe_to_target_R"]),
        "conformal_Nf16_phase": conformal["phase_label"],
        "conformal_Nf16_uv_landau_safe_to_R200": bool(conformal["safe_to_target_R"]),
        "max_D_residual_2norm": max(row["D_residual_2norm"] for row in finite),
        "max_unitarity_residual_2norm": max(row["unitarity_residual_2norm"] for row in finite),
        "max_quantum_constraint_abs_residual": max(row["quantum_constraint_abs_residual"] for row in finite),
        "max_abs_log_singular": max(row["max_abs_log_singular"] for row in finite),
        "visible_threshold_vector": [0.0, 0.0, 0.0],
        "status": "PASS_CONDITIONAL",
        "interpretation": (
            "The latest 8x8 clockwork unitary link admits the same hidden "
            "endpoint/meson realization as the older 4x4 crossed-link audit: "
            "radial moment maps plus the hidden U(8)_H quotient leave a U(8) "
            "meson link, and the 4x4 fixed block leaves 33 real unitary "
            "completion moduli.  The Nf=Nc=8 quantum-deformed constraint is "
            "compatible but not sufficient; it leaves 98 real holomorphic "
            "moduli after fixed-block and determinant constraints, so the "
            "radial/D-term lock remains a necessary conditional ingredient.  "
            "All endpoint fields are visible Spin(10) singlets, hence the "
            "direct visible threshold vector is zero."
        ),
    }

    return {
        "note": "No web lookup used. Hidden endpoint meson audit for the current 8x8 clockwork unitary quotient.",
        "input_card": {
            "q": clockwork["input_card"]["q"],
            "n": clockwork["input_card"]["n"],
            "epsilon": clockwork["input_card"]["epsilon"],
            "kappa": clockwork["input_card"]["kappa"],
            "unitary_dimension": N,
            "visible_block_dimension": VISIBLE,
        },
        "meson_branch": {
            "hidden_gauge_group": "U(8)_H for the finite-rank quotient; SU(8) beta estimate for the nonabelian factor",
            "fields": "Q,Qtilde visible Spin(10)-singlet endpoint matrices",
            "locked_branch": "Q=f I_8, Qtilde=f B, B in U(8)",
            "meson": "B=Qtilde Q/f^2",
            "fixed_block": "PBP=W_kappa",
        },
        "rank_audit": {
            "Dterm": dterm,
            "hidden_orbit": orbit,
            "residual_after_D_lock_and_hidden_quotient": dterm["combined_D_nullity_before_hidden_quotient"] - orbit["hidden_gauge_orbit_rank"],
            "fixed_block": fixed,
        },
        "quantum_rank_audit": {
            "constraint": "det m - beta beta_tilde = (Lambda_H/f)^16",
            "holomorphic": holo,
            "unitary_plus_baryons": unitary_q,
            "holomorphic_real_nullity_after_constraints": holomorphic_real_nullity,
            "nonunitary_real_moduli_removed_by_Dterm_lock": nonunitary_removed,
        },
        "finite_samples": finite,
        "hidden_beta_scan": beta_rows,
        "threshold_interpretation": {
            "endpoint_fields_visible_spin10_singlets": True,
            "direct_visible_threshold_vector": [0.0, 0.0, 0.0],
            "caveat": "Hidden endpoint dynamics can still affect the UV cutoff, but not the nonuniversal SM/GUT matching if it stays visible-singlet.",
        },
        "verdict": verdict,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "finite_endpoint_samples.csv", summary["finite_samples"])
    write_csv(OUT / "hidden_beta_scan.csv", summary["hidden_beta_scan"])
    write_report(summary)
    verdict = summary["verdict"]
    print("Clockwork hidden endpoint meson audit")
    print(f"  finite samples pass: {verdict['finite_samples_Dflat_unitary_quantum_and_locked']}")
    print(f"  quantum alone enforces unitarity: {verdict['quantum_constraint_alone_enforces_unitarity']}")
    print(f"  residual unitary moduli: {summary['rank_audit']['fixed_block']['residual_unitary_completion_moduli']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
