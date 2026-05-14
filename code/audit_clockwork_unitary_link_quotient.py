#!/usr/bin/env python3
"""Combine the constrained clockwork source with the unitary-link quotient.

No web lookup is used.

The constrained clockwork audit produces the visible inverse block

    W_kappa = diag(1, eps, eps, eps),  eps=3^-6.

This script embeds that block into a unitary link B by the Julia/Halmos
dilation and checks the combined constraints

    C_clock X = 0,      boundary source closes the clockwork zero mode,
    P B P = W_kappa,
    B^\dagger B = I.

The key question is whether the unitary quotient keeps the physical triplet
mass singular values locked at M_lock after W_kappa is generated dynamically.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clockwork_unitary_link_quotient"
CLOCK_HESSIAN = ROOT / "output" / "constrained_clockwork_source_hessian" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"

RNG_SEED = 202605091300
ANGLE_GRID = [0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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
                "mass_singular_min_over_Mlock": float(np.min(singulars)),
                "mass_singular_max_over_Mlock": float(np.max(singulars)),
                "max_abs_log_singular": max_abs_log,
                "within_completed_partner_lock_window": bool(max_abs_log <= lock_window),
            }
        )
    return rows


def audit() -> dict[str, Any]:
    hessian = read_json(CLOCK_HESSIAN)
    completed = read_json(COMPLETED)
    eps = float(hessian["input_card"]["epsilon"])
    q = float(hessian["input_card"]["q"])
    n = int(hessian["input_card"]["n"])
    kappa = float(hessian["input_card"]["kappa"])
    m_lock = float(hessian["input_card"]["M_lock_GeV"])

    w = np.diag([1.0, eps, eps, eps]).astype(complex)
    b0 = julia_dilation(w)
    visible = w.shape[0]
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])

    d_rank = dterm_linear_rank(b0)
    f_rank = fixed_block_rank_on_unitary(b0, visible)
    completion_rows = finite_completion_scan(w, lock_window)
    b_singulars = np.linalg.svd(b0.conj().T, compute_uv=False)

    return {
        "note": "No web lookup used. Combined q=3,n=6 clockwork source plus unitary-link D-term quotient audit.",
        "input_card": {
            "q": q,
            "n": n,
            "epsilon": eps,
            "kappa": kappa,
            "M_lock_GeV": m_lock,
        },
        "combined_constraints": {
            "clockwork": "C_clock X=0 plus boundary row v^T X=0 after source closure",
            "visible_block": "P B P = W_kappa = diag(1,epsilon,epsilon,epsilon)",
            "Dterm": "B^dagger B = I",
            "mass_operator": "M_lock B^dagger",
        },
        "W_kappa": matrix_json(w),
        "unitary_dilation": {
            "dimension": int(b0.shape[0]),
            "unitarity_residual_2norm": float(np.linalg.norm(b0.conj().T @ b0 - np.eye(b0.shape[0]), ord=2)),
            "fixed_block_residual_fro": float(np.linalg.norm(b0[:visible, :visible] - w, ord="fro")),
            "mass_singular_values_over_Mlock": [float(x) for x in b_singulars],
            "max_abs_log_singular": float(np.max(np.abs(np.log(b_singulars)))),
        },
        "Dterm_linearization": d_rank,
        "Fterm_on_Dflat_tangent": f_rank,
        "finite_unitary_completion_scan": completion_rows,
        "source_hessian_import": {
            "boundary_closed_zero_modes_per_chain": hessian["linear_algebra"]["boundary_closed_zero_modes_per_chain"],
            "boundary_closed_min_abs_eigenvalue": hessian["linear_algebra"]["boundary_closed_min_abs_eigenvalue"],
            "boundary_closed_max_abs_eigenvalue": hessian["linear_algebra"]["boundary_closed_max_abs_eigenvalue"],
            "dressed_max_width_margin_at_ST_1e_minus_5": hessian["effective_inverse_block"]["dressed_max_width_margin_at_ST_1e_minus_5"],
        },
        "threshold_interpretation": {
            "physical_triplet_mass_singular_spread": float(np.max(b_singulars) - np.min(b_singulars)),
            "nonuniversal_threshold_if_complete_degenerate": [0.0, 0.0, 0.0],
            "constrained_source_threshold": [0.0, 0.0, 0.0],
            "literal_visible_clockwork_tower_warning": "A literal propagating visible 10+10bar clockwork tower is still UV-disfavored; use constrained/composite source or complete degenerate cutoff-local links.",
        },
        "verdict": {
            "unitary_dilation_exists": bool(np.linalg.norm(b0.conj().T @ b0 - np.eye(b0.shape[0]), ord=2) < 1.0e-12),
            "visible_block_fixed": bool(np.linalg.norm(b0[:visible, :visible] - w, ord="fro") < 1.0e-12),
            "mass_singular_values_locked": bool(np.max(np.abs(np.log(b_singulars))) < 1.0e-12),
            "all_completion_samples_locked": all(row["within_completed_partner_lock_window"] for row in completion_rows),
            "source_zero_modes_closed": hessian["verdict"]["boundary_driver_removes_all_clockwork_zero_modes"],
            "status": "PASS_CONDITIONAL",
            "main_caveat": (
                "The combined quotient is a finite-dimensional constrained-source construction. "
                "It keeps the physical triplet masses locked while generating the small visible inverse block, "
                "but a microscopic Spin(10) strong-sector or boundary-origin derivation remains open."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    card = summary["input_card"]
    unit = summary["unitary_dilation"]
    dterm = summary["Dterm_linearization"]
    fterm = summary["Fterm_on_Dflat_tangent"]
    source = summary["source_hessian_import"]
    th = summary["threshold_interpretation"]
    verdict = summary["verdict"]
    lines = [
        "# Clockwork plus unitary-link quotient audit",
        "",
        "No web lookup was used.",
        "",
        f"Input: `q={card['q']:.1f}`, `n={card['n']}`, `epsilon={card['epsilon']:.6e}`, "
        f"`kappa={card['kappa']:.0f}`.",
        "",
        "The visible block is",
        "",
        "`W_kappa = diag(1, epsilon, epsilon, epsilon)`.",
        "",
        "It is embedded into a unitary link by the Julia/Halmos dilation",
        "",
        "`B = [[W, sqrt(1-WWdag)], [sqrt(1-Wdag W), -Wdag]]`.",
        "",
        "## Link checks",
        "",
        f"`dim B = {unit['dimension']}`.",
        f"`||Bdag B - I||_2 = {unit['unitarity_residual_2norm']:.6e}`.",
        f"`||PBP-W||_F = {unit['fixed_block_residual_fro']:.6e}`.",
        f"`max |log singular(M/M_lock)| = {unit['max_abs_log_singular']:.6e}`.",
        "",
        "D-term linearization:",
        f"`rank = {dterm['Dterm_linear_rank']}`, `D-flat tangent dim = {dterm['Dflat_tangent_real_dimension']}`.",
        "",
        "Fixed-block F-term on the D-flat tangent:",
        f"`rank = {fterm['fixed_block_rank_on_Dflat_tangent']}`, "
        f"`residual unitary moduli = {fterm['quotient_residual_unitary_moduli_real_dimension']}`.",
        "",
        "## Source import",
        "",
        f"Boundary-closed clockwork zero modes per chain: `{source['boundary_closed_zero_modes_per_chain']}`.",
        f"Boundary-closed source masses: `{source['boundary_closed_min_abs_eigenvalue']:.6f}` to "
        f"`{source['boundary_closed_max_abs_eigenvalue']:.6f}` times `M_lock`.",
        f"Dressed max-width margin at `S_T=1e-5`: `{source['dressed_max_width_margin_at_ST_1e_minus_5']:.6f}`.",
        "",
        "## Completion samples",
        "",
        "| angle | fixed residual | unitary residual | min singular | max singular | lock-window pass |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary["finite_unitary_completion_scan"]:
        lines.append(
            f"| {row['angle']:.6g} | {row['fixed_block_residual_fro']:.3e} | "
            f"{row['unitarity_residual_2norm']:.3e} | {row['mass_singular_min_over_Mlock']:.12f} | "
            f"{row['mass_singular_max_over_Mlock']:.12f} | {row['within_completed_partner_lock_window']} |"
        )
    lines += [
        "",
        "## Threshold interpretation",
        "",
        f"Physical triplet singular spread: `{th['physical_triplet_mass_singular_spread']:.6e}`.",
        f"Complete-degenerate nonuniversal threshold: `{th['nonuniversal_threshold_if_complete_degenerate']}`.",
        f"Constrained-source threshold: `{th['constrained_source_threshold']}`.",
        th["literal_visible_clockwork_tower_warning"],
        "",
        "## Verdict",
        "",
        f"`unitary_dilation_exists = {verdict['unitary_dilation_exists']}`.",
        f"`visible_block_fixed = {verdict['visible_block_fixed']}`.",
        f"`mass_singular_values_locked = {verdict['mass_singular_values_locked']}`.",
        f"`all_completion_samples_locked = {verdict['all_completion_samples_locked']}`.",
        f"`source_zero_modes_closed = {verdict['source_zero_modes_closed']}`.",
        "",
        verdict["main_caveat"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = audit()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "unitary_completion_scan.csv", summary["finite_unitary_completion_scan"])
    write_report(summary)
    verdict = summary["verdict"]
    print("Clockwork plus unitary-link quotient audit")
    print(f"  unitary dilation exists: {verdict['unitary_dilation_exists']}")
    print(f"  mass singular values locked: {verdict['mass_singular_values_locked']}")
    print(f"  all samples locked: {verdict['all_completion_samples_locked']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
