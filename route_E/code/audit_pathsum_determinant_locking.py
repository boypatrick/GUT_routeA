#!/usr/bin/env python3
"""Path-sum determinant-locking audit.

No web lookup is used.  This script formalizes the graph hint from the
thread: gauge thresholds see closed determinants, while dimension-five
operators see open propagator paths.  For a block triplet mass matrix

    M = [[A, U],
         [V, G0 + V A^{-1} U]],

the Schur complement with respect to A is exactly G0, so

    det M = det A det G0,

but the open A-to-A inverse propagator is shifted:

    (M^{-1})_AA = A^{-1} + A^{-1} U G0^{-1} V A^{-1}.

This is a finite-dimensional, action-level version of the determinant-locked
10'_G bridge: it changes Wilson tensors without creating a one-loop
non-universal gauge threshold.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from audit_gtr_mediator_spectrum import B_TRIPLET_PAIR, projected_l2


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "pathsum_determinant_locking"


def threshold_norm_from_det_ratio(det_ratio: float) -> float:
    delta = B_TRIPLET_PAIR * math.log(1.0 / abs(det_ratio)) / (2.0 * math.pi)
    return projected_l2(delta)


def make_benchmark(n_a: int, n_g: int, eps: float) -> dict[str, np.ndarray]:
    a = np.diag(np.linspace(1.0, 1.0 + 0.25 * (n_a - 1), n_a))
    g0 = np.diag(np.linspace(1.0, 1.0 + 0.20 * (n_g - 1), n_g))
    u_raw = np.fromfunction(lambda i, j: ((i + 1.0) * (j + 2.0)) % 5.0 + 1.0, (n_a, n_g), dtype=float)
    v_raw = np.fromfunction(lambda i, j: ((i + 2.0) * (j + 1.0)) % 7.0 + 1.0, (n_g, n_a), dtype=float)
    u = eps * u_raw / np.linalg.norm(u_raw)
    v = eps * v_raw / np.linalg.norm(v_raw)
    return {"A": a, "G0": g0, "U": u, "V": v}


def block_matrix(a: np.ndarray, u: np.ndarray, v: np.ndarray, g: np.ndarray) -> np.ndarray:
    return np.block([[a, u], [v, g]])


def cauchy_pathsum_metrics(n_a: int, n_g: int, eps: float) -> dict[str, Any]:
    bench = make_benchmark(n_a, n_g, eps)
    a, g0, u, v = bench["A"], bench["G0"], bench["U"], bench["V"]
    a_inv = np.linalg.inv(a)
    g0_inv = np.linalg.inv(g0)

    g_locked = g0 + v @ a_inv @ u
    m_locked = block_matrix(a, u, v, g_locked)
    m_naive = block_matrix(a, u, v, g0)
    target_det = float(np.linalg.det(a) * np.linalg.det(g0))

    inv_locked = np.linalg.inv(m_locked)
    inv_naive = np.linalg.inv(m_naive)
    inv_aa_locked = inv_locked[:n_a, :n_a]
    inv_aa_naive = inv_naive[:n_a, :n_a]
    identity_inv_aa = a_inv + a_inv @ u @ g0_inv @ v @ a_inv

    y_left = np.linspace(1.0, 0.35, n_a)
    y_right = np.linspace(0.45, 1.15, n_a)
    c0 = float(y_left @ a_inv @ y_right)
    c_locked = float(y_left @ inv_aa_locked @ y_right)
    c_naive = float(y_left @ inv_aa_naive @ y_right)

    det_ratio_locked = float(np.linalg.det(m_locked) / target_det)
    det_ratio_naive = float(np.linalg.det(m_naive) / target_det)
    return {
        "n_A": n_a,
        "n_G": n_g,
        "epsilon": eps,
        "det_ratio_locked": det_ratio_locked,
        "det_ratio_naive": det_ratio_naive,
        "abs_log_det_locked": abs(math.log(abs(det_ratio_locked))),
        "abs_log_det_naive": abs(math.log(abs(det_ratio_naive))),
        "projected_threshold_locked_l2": threshold_norm_from_det_ratio(det_ratio_locked),
        "projected_threshold_naive_l2": threshold_norm_from_det_ratio(det_ratio_naive),
        "inverse_identity_error_locked": float(np.linalg.norm(inv_aa_locked - identity_inv_aa)),
        "inverse_shift_locked_fro": float(np.linalg.norm(inv_aa_locked - a_inv)),
        "inverse_shift_naive_fro": float(np.linalg.norm(inv_aa_naive - a_inv)),
        "wilson_proxy_baseline": c0,
        "wilson_proxy_locked": c_locked,
        "wilson_proxy_naive": c_naive,
        "wilson_proxy_locked_relative_shift": float((c_locked - c0) / c0),
        "wilson_proxy_naive_relative_shift": float((c_naive - c0) / c0),
        "min_singular_locked": float(np.linalg.svd(m_locked, compute_uv=False)[-1]),
        "min_singular_naive": float(np.linalg.svd(m_naive, compute_uv=False)[-1]),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Path-sum determinant-locking audit",
        "",
        "No web lookup was used.",
        "",
        "## Theorem",
        "",
        "For `M=[[A,U],[V,G0+V A^{-1} U]]`,",
        "",
        "```text",
        "det M = det A det G0",
        "(M^{-1})_AA = A^{-1}+A^{-1} U G0^{-1} V A^{-1}.",
        "```",
        "",
        "The determinant is closed-path silent, but the open propagator is shifted.",
        "",
        "## Numerical rows",
        "",
        "| nA | nG | eps | |log det locked| | threshold locked | inverse shift | Wilson shift |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["rows"]:
        lines.append(
            f"| {row['n_A']} | {row['n_G']} | {row['epsilon']:.2f} | "
            f"{row['abs_log_det_locked']:.3e} | {row['projected_threshold_locked_l2']:.3e} | "
            f"{row['inverse_shift_locked_fro']:.3e} | {row['wilson_proxy_locked_relative_shift']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", payload["verdict"]])
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for n_a, n_g in [(1, 1), (2, 1), (3, 2)]:
        for eps in [0.05, 0.10, 0.20, 0.30, 0.50]:
            rows.append(cauchy_pathsum_metrics(n_a, n_g, eps))

    write_csv(OUT / "pathsum_determinant_locking_scan.csv", rows)
    max_locked_threshold = max(row["projected_threshold_locked_l2"] for row in rows)
    max_identity_error = max(row["inverse_identity_error_locked"] for row in rows)
    max_wilson_shift = max(abs(row["wilson_proxy_locked_relative_shift"]) for row in rows)
    benchmark = next(row for row in rows if row["n_A"] == 3 and row["n_G"] == 2 and abs(row["epsilon"] - 0.30) < 1.0e-12)

    verdict = (
        "The path-sum construction verifies the image-inspired mechanism in finite-dimensional algebra. "
        f"Across the tested blocks the largest locked projected threshold is {max_locked_threshold:.3e} "
        f"and the Schur-complement inverse identity error is {max_identity_error:.3e}.  Nevertheless "
        f"the Wilson proxy changes by up to {max_wilson_shift:.3e}.  In the 3x2 benchmark at epsilon=0.30, "
        f"the determinant-locked threshold is {benchmark['projected_threshold_locked_l2']:.3e} while the "
        f"open-path Wilson proxy shifts by {benchmark['wilson_proxy_locked_relative_shift']:.3e}.  This "
        "makes determinant-silent open-path deformation a precise next principle for the 10'_G bridge: "
        "promote the Schur-complement equation G=G0+V A^{-1}U to F-flatness and then feed the shifted "
        "inverse block into the dimension-five Wilson scan."
    )
    payload = {
        "note": "No web lookup used. Schur-complement/path-sum determinant-locking audit.",
        "theorem": {
            "matrix": "M=[[A,U],[V,G0+V A^{-1}U]]",
            "determinant": "det M=det A det G0",
            "inverse_AA": "(M^{-1})_AA=A^{-1}+A^{-1}U G0^{-1} V A^{-1}",
        },
        "rows": rows,
        "diagnostics": {
            "max_locked_projected_threshold_l2": max_locked_threshold,
            "max_inverse_identity_error_locked": max_identity_error,
            "max_abs_wilson_proxy_locked_relative_shift": max_wilson_shift,
            "benchmark_3x2_eps_0p30": benchmark,
        },
        "verdict": verdict,
    }
    (OUT / "pathsum_determinant_locking_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "pathsum_determinant_locking_report.md").write_text(make_report(payload), encoding="utf-8")
    print("Path-sum determinant-locking audit")
    print(f"  max locked projected threshold: {max_locked_threshold:.3e}")
    print(f"  max inverse identity error: {max_identity_error:.3e}")
    print(f"  3x2 eps=0.30 Wilson shift: {benchmark['wilson_proxy_locked_relative_shift']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
