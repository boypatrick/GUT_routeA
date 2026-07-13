#!/usr/bin/env python3
"""Construct a threshold-locked unitary dilation of triplet near-null cards.

The previous finite-lift audit treated the target inverse triplet propagator W
as if its singular values were physical inverse masses.  That is sufficient for
a direct EFT mass matrix, but it is too restrictive: a subblock of the inverse
of a larger, exactly degenerate heavy mass matrix can be rank-deficient or
hierarchical even when the full physical spectrum is perfectly degenerate.

This script uses the Julia/Sz.-Nagy unitary dilation.  For any contraction
||W||_2 <= 1,

    U(W) = [[ W,        (I - W W^dagger)^{1/2}],
            [(I - W^dagger W)^{1/2},       -W^dagger]]

is unitary.  If the full triplet/mediator superpotential mass matrix is

    M_full = M_lock U(W)^dagger,

then all physical singular values are exactly M_lock, while the visible
inverse-propagator block is

    P M_full^{-1} P = W / M_lock.

Thus the near-null can live in an effective Wilson subblock rather than in a
split physical spectrum.  This is still a conditional construction: the full
mediator fields must be embedded in complete degenerate GUT multiplets for the
threshold cancellation to be physical.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "threshold_locked_triplet_lift"

RANK_LIFT = ROOT / "output" / "triplet_rank_lift" / "triplet_rank_lift_summary.json"
WIDTH_LIFT = ROOT / "output" / "nearnull_finite_lift_width" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"


TARGET_LABELS = ["rank3_limit", "condition_cap_10", "condition_cap_30", "condition_cap_100"]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    vals = np.where(vals < tol, np.maximum(vals, 0.0), vals)
    if np.min(vals) < -tol:
        raise RuntimeError(f"defect operator has a negative eigenvalue {np.min(vals)}")
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def target_matrices() -> dict[str, np.ndarray]:
    rank = read_json(RANK_LIFT)
    w_near = cmat(rank["near_null_W"]["matrix"])
    u, singulars, vh = np.linalg.svd(w_near, full_matrices=True)
    s1 = float(singulars[0])
    out = {"rank3_limit": w_near}
    for kappa in [10.0, 30.0, 100.0]:
        capped = np.maximum(singulars, s1 / kappa)
        out[f"condition_cap_{int(kappa)}"] = u @ np.diag(capped) @ vh
    return out


def width_lookup() -> dict[str, dict[str, Any]]:
    payload = read_json(WIDTH_LIFT)
    return {row["label"]: row for row in payload["lift_summaries"]}


def julia_dilation(w: np.ndarray) -> np.ndarray:
    n = w.shape[0]
    eye = np.eye(n, dtype=complex)
    left_defect = hermitian_sqrt(eye - w @ w.conj().T)
    right_defect = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, left_defect],
            [right_defect, -w.conj().T],
        ]
    )


def audit_target(label: str, w: np.ndarray, width: dict[str, Any], m_lock: float, threshold_tol: float) -> dict[str, Any]:
    n = w.shape[0]
    u = julia_dilation(w)
    eye2 = np.eye(2 * n, dtype=complex)
    unitary_residual = float(np.linalg.norm(u.conj().T @ u - eye2, ord=2))
    inv_subblock = u[:n, :n]
    subblock_residual = float(np.linalg.norm(inv_subblock - w, ord="fro"))
    full_singulars = np.linalg.svd(u.conj().T, compute_uv=False)
    full_spread = float(np.max(full_singulars) - np.min(full_singulars))
    w_singulars = np.linalg.svd(w, compute_uv=False)
    defect_left_eigs = np.linalg.eigvalsh(np.eye(n) - w @ w.conj().T)
    defect_right_eigs = np.linalg.eigvalsh(np.eye(n) - w.conj().T @ w)
    threshold_projected_l2 = 0.0
    return {
        "label": label,
        "visible_dimension": n,
        "mediator_dimension": n,
        "full_dimension": 2 * n,
        "target_singular_values": [float(x) for x in w_singulars],
        "target_spectral_norm": float(np.max(w_singulars)),
        "unitary_residual_2norm": unitary_residual,
        "visible_inverse_subblock_residual_fro": subblock_residual,
        "full_mass_singular_values_over_Mlock": [float(x) for x in full_singulars],
        "full_mass_singular_spread_over_Mlock": full_spread,
        "M_lock_GeV": m_lock,
        "M_max_GeV": m_lock,
        "below_reduced_planck": bool(m_lock < 2.435e18),
        "defect_left_min_eigenvalue": float(np.min(defect_left_eigs)),
        "defect_right_min_eigenvalue": float(np.min(defect_right_eigs)),
        "threshold_projected_l2_if_complete_degenerate": threshold_projected_l2,
        "r200_projected_threshold_tolerance": threshold_tol,
        "passes_threshold_if_complete_degenerate": threshold_projected_l2 <= threshold_tol,
        "worst_future_margin_1e35": width["worst_future_margin_1e35"],
        "min_LLLL_future_margin": width["min_LLLL_future_margin"],
        "min_RRRR_future_margin": width["min_RRRR_future_margin"],
        "proton_safe_future_1e35": bool(width["proton_safe_future_1e35"]),
        "worst_future_channel": width["worst_future_channel"],
        "passes_all_conditional_checks": bool(
            width["proton_safe_future_1e35"]
            and m_lock < 2.435e18
            and threshold_projected_l2 <= threshold_tol
            and unitary_residual < 1.0e-10
            and subblock_residual < 1.0e-10
        ),
        "unitary_dilation_matrix": matrix_json(u),
    }


def build() -> dict[str, Any]:
    mats = target_matrices()
    widths = width_lookup()
    ledger = read_json(LEDGER)
    m_lock = float(ledger["r200_benchmark"]["M_Sigma8_GeV"])
    threshold_tol = float(ledger["r200_benchmark"]["total_projected_l2"])
    rows = [audit_target(label, mats[label], widths[label], m_lock, threshold_tol) for label in TARGET_LABELS]
    passing = [row for row in rows if row["passes_all_conditional_checks"]]
    best = max(passing, key=lambda row: row["worst_future_margin_1e35"], default=None)
    lowest_dimension = min(passing, key=lambda row: row["full_dimension"], default=None)
    verdict = {
        "passing_conditional_rows": len(passing),
        "best_label": None if best is None else best["label"],
        "best_worst_future_margin_1e35": None if best is None else best["worst_future_margin_1e35"],
        "lowest_dimension_passing_label": None if lowest_dimension is None else lowest_dimension["label"],
        "M_lock_GeV": m_lock,
        "threshold_projected_l2": 0.0,
        "r200_projected_threshold_tolerance": threshold_tol,
        "interpretation": (
            "The Julia/Sz.-Nagy dilation realizes each target inverse propagator "
            "as the visible subblock of the inverse of an exactly degenerate "
            "8-by-8 triplet/mediator mass matrix.  This removes the rank and "
            "threshold-splitting obstruction at the algebraic EFT level.  The "
            "remaining hard task is a Spin(10) representation and symmetry "
            "construction that enforces this unitary dilation with complete "
            "degenerate multiplets."
        ),
    }
    return {
        "note": "No web lookup used. Threshold-locked unitary dilation of triplet near-null cards.",
        "mathematical_identity": (
            "For ||W||_2<=1, U=[[W,sqrt(I-WWdag)],[sqrt(I-WdagW),-Wdag]] is unitary. "
            "With M_full=M_lock Udag, P M_full^{-1} P = W/M_lock and all physical singular values are M_lock."
        ),
        "input_files": {
            "rank_lift": str(RANK_LIFT),
            "finite_lift_width": str(WIDTH_LIFT),
            "conditional_theorem_ledger": str(LEDGER),
        },
        "rows": rows,
        "verdict": verdict,
    }


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "target_spectral_norm",
        "unitary_residual_2norm",
        "visible_inverse_subblock_residual_fro",
        "full_mass_singular_spread_over_Mlock",
        "M_lock_GeV",
        "below_reduced_planck",
        "worst_future_margin_1e35",
        "min_LLLL_future_margin",
        "min_RRRR_future_margin",
        "threshold_projected_l2_if_complete_degenerate",
        "passes_all_conditional_checks",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Threshold-locked triplet-lift audit",
        "",
        "No web lookup was used.",
        "",
        "For any contraction \(W\), use the Julia/Sz.-Nagy dilation",
        "",
        "```text",
        "U(W) = [[W, sqrt(I-W Wdag)], [sqrt(I-Wdag W), -Wdag]].",
        "M_full = M_lock U(W)dag.",
        "P M_full^{-1} P = W / M_lock.",
        "```",
        "",
        "| label | ||W||2 | unitary residual | subblock residual | singular spread | worst future margin | threshold | pass |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        lines.append(
            "| `{label}` | {target_spectral_norm:.6e} | {unitary_residual_2norm:.3e} | "
            "{visible_inverse_subblock_residual_fro:.3e} | {full_mass_singular_spread_over_Mlock:.3e} | "
            "{worst_future_margin_1e35:.6e} | {threshold_projected_l2_if_complete_degenerate:.1e} | "
            "{passes_all_conditional_checks} |".format(**row)
        )
    v = summary["verdict"]
    lines += [
        "",
        "## Verdict",
        "",
        f"Conditional passing rows: {v['passing_conditional_rows']}.",
        f"Best row: `{v['best_label']}` with future margin {v['best_worst_future_margin_1e35']:.6e}.",
        f"Locked common mass: {v['M_lock_GeV']:.6e} GeV.",
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_summary_csv(OUT / "locked_lift_summary.csv", summary["rows"])
    write_report(summary)
    v = summary["verdict"]
    print("Threshold-locked triplet-lift audit")
    print(f"  passing conditional rows: {v['passing_conditional_rows']}")
    print(f"  best row: {v['best_label']}")
    print(f"  M_lock: {v['M_lock_GeV']:.6e} GeV")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
