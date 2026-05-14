#!/usr/bin/env python3
"""Polynomial F-flat realization of Schur-complement determinant locking.

No web lookup is used.  The previous path-sum audit used the algebraic
condition

    G = G0 + V A^{-1} U.

This script replaces the inverse by polynomial F-term equations.  Introduce
auxiliary singlet matrices B and C plus driver matrices X, Z, Y, with

    W_lock = Tr X (A B - I)
           + Tr Z (C - B U)
           + Tr Y (G - G0 - V C).

All interaction terms are at most cubic.  The F-driver equations imply

    B = A^{-1},   C = B U,   G = G0 + V C,

so the Schur-complement determinant identity follows without writing a
non-polynomial inverse in the superpotential.  With X=Y=Z=0, all non-driver
F-terms vanish on this constraint locus.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from audit_pathsum_determinant_locking import cauchy_pathsum_metrics, make_benchmark


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "polynomial_schur_locking"


def fro(mat: np.ndarray) -> float:
    return float(np.linalg.norm(mat))


def block_matrix(a: np.ndarray, u: np.ndarray, v: np.ndarray, g: np.ndarray) -> np.ndarray:
    return np.block([[a, u], [v, g]])


def f_flat_solution(n_a: int, n_g: int, eps: float) -> dict[str, Any]:
    bench = make_benchmark(n_a, n_g, eps)
    a, g0, u, v = bench["A"], bench["G0"], bench["U"], bench["V"]
    b = np.linalg.inv(a)
    c = b @ u
    g = g0 + v @ c

    # Driver fields vanish on the chosen branch.
    x = np.zeros((n_a, n_a))
    z = np.zeros((n_g, n_a))
    y = np.zeros((n_g, n_g))

    f_x = a @ b - np.eye(n_a)
    f_z = c - b @ u
    f_y = g - g0 - v @ c

    # Non-driver F-terms are all proportional to X, Z, or Y and hence vanish
    # on the X=Z=Y=0 branch.  We keep explicit residuals to catch shape errors.
    f_a = x @ b.T
    f_b = a.T @ x.T - u @ z
    f_u = -b.T @ z.T
    f_c = z.T - v.T @ y
    f_v = -y @ c.T
    f_g = y

    m = block_matrix(a, u, v, g)
    det_target = float(np.linalg.det(a) * np.linalg.det(g0))
    inv_aa = np.linalg.inv(m)[:n_a, :n_a]
    inv_formula = b + b @ u @ np.linalg.inv(g0) @ v @ b
    path_metrics = cauchy_pathsum_metrics(n_a, n_g, eps)

    return {
        "n_A": n_a,
        "n_G": n_g,
        "epsilon": eps,
        "F_X_AB_minus_I": fro(f_x),
        "F_Z_C_minus_BU": fro(f_z),
        "F_Y_G_minus_G0_minus_VC": fro(f_y),
        "F_non_driver_max": max(fro(f_a), fro(f_b), fro(f_u), fro(f_c), fro(f_v), fro(f_g)),
        "det_ratio": float(np.linalg.det(m) / det_target),
        "abs_log_det_ratio": float(abs(np.log(abs(np.linalg.det(m) / det_target)))),
        "inverse_AA_identity_error": fro(inv_aa - inv_formula),
        "inverse_AA_shift_fro": fro(inv_aa - b),
        "pathsum_locked_threshold_l2": path_metrics["projected_threshold_locked_l2"],
        "pathsum_wilson_proxy_shift": path_metrics["wilson_proxy_locked_relative_shift"],
        "min_singular_full_matrix": float(np.linalg.svd(m, compute_uv=False)[-1]),
        "max_abs_B_entry": float(np.max(np.abs(b))),
        "max_abs_C_entry": float(np.max(np.abs(c))),
        "max_abs_G_shift_entry": float(np.max(np.abs(g - g0))),
    }


def term_ledger() -> list[dict[str, Any]]:
    return [
        {
            "term": "Tr X A B",
            "degree": 3,
            "role": "enforce B as right inverse of A",
            "renormalizable": True,
        },
        {
            "term": "-Tr X",
            "degree": 1,
            "role": "dimensionful identity source in F_X",
            "renormalizable": True,
        },
        {
            "term": "Tr Z C",
            "degree": 2,
            "role": "define intermediate path variable C",
            "renormalizable": True,
        },
        {
            "term": "-Tr Z B U",
            "degree": 3,
            "role": "polynomial replacement of C=B U",
            "renormalizable": True,
        },
        {
            "term": "Tr Y G - Tr Y G0",
            "degree": 2,
            "role": "source sterile block G around G0",
            "renormalizable": True,
        },
        {
            "term": "-Tr Y V C",
            "degree": 3,
            "role": "polynomial replacement of G=G0+V C",
            "renormalizable": True,
        },
    ]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Polynomial Schur-locking F-flatness audit",
        "",
        "No web lookup was used.",
        "",
        "## Superpotential",
        "",
        "```text",
        "W_lock = Tr X(A B - I) + Tr Z(C - B U) + Tr Y(G - G0 - V C).",
        "```",
        "",
        "All non-source interaction terms have degree <= 3.",
        "",
        "## Numerical F-flatness",
        "",
        "| nA | nG | eps | max driver F | max non-driver F | |log det| | inverse error | Wilson shift |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["rows"]:
        max_driver = max(
            row["F_X_AB_minus_I"],
            row["F_Z_C_minus_BU"],
            row["F_Y_G_minus_G0_minus_VC"],
        )
        lines.append(
            f"| {row['n_A']} | {row['n_G']} | {row['epsilon']:.2f} | "
            f"{max_driver:.3e} | {row['F_non_driver_max']:.3e} | "
            f"{row['abs_log_det_ratio']:.3e} | {row['inverse_AA_identity_error']:.3e} | "
            f"{row['pathsum_wilson_proxy_shift']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", payload["verdict"]])
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for n_a, n_g in [(1, 1), (2, 1), (3, 2)]:
        for eps in [0.10, 0.30, 0.50]:
            rows.append(f_flat_solution(n_a, n_g, eps))
    write_csv(OUT / "polynomial_schur_locking_scan.csv", rows)
    ledger = term_ledger()

    max_driver_f = max(
        max(row["F_X_AB_minus_I"], row["F_Z_C_minus_BU"], row["F_Y_G_minus_G0_minus_VC"])
        for row in rows
    )
    max_non_driver_f = max(row["F_non_driver_max"] for row in rows)
    max_log_det = max(row["abs_log_det_ratio"] for row in rows)
    max_inverse_error = max(row["inverse_AA_identity_error"] for row in rows)
    benchmark = next(row for row in rows if row["n_A"] == 3 and row["n_G"] == 2 and abs(row["epsilon"] - 0.30) < 1.0e-12)

    verdict = (
        "The inverse in the path-sum condition can be eliminated with a fully polynomial F-flat system. "
        "The auxiliary matrices B and C implement B=A^{-1} and C=BU through driver equations, while "
        "G=G0+VC gives the determinant-locked Schur complement.  Across the scan the largest driver "
        f"F-residual is {max_driver_f:.3e}, the largest non-driver F-residual on the X=Y=Z=0 branch is "
        f"{max_non_driver_f:.3e}, and the largest inverse identity error is {max_inverse_error:.3e}. "
        f"The 3+2, epsilon=0.30 point keeps |log det|={benchmark['abs_log_det_ratio']:.3e} while shifting "
        f"the Wilson proxy by {benchmark['pathsum_wilson_proxy_shift']:.3e}.  The next test is not algebraic "
        "but representational: assign Spin(10)/grading charges to A,U,V,G,B,C,X,Y,Z so that these singlet "
        "drivers do not reintroduce doublet mixing or incomplete thresholds."
    )
    payload = {
        "note": "No web lookup used. Polynomial F-flat realization of Schur-complement determinant locking.",
        "superpotential": "W_lock = Tr X(A B - I) + Tr Z(C - B U) + Tr Y(G - G0 - V C)",
        "term_ledger": ledger,
        "rows": rows,
        "diagnostics": {
            "max_driver_F_residual": max_driver_f,
            "max_non_driver_F_residual": max_non_driver_f,
            "max_abs_log_det_ratio": max_log_det,
            "max_inverse_identity_error": max_inverse_error,
            "benchmark_3x2_eps_0p30": benchmark,
        },
        "verdict": verdict,
    }
    (OUT / "polynomial_schur_locking_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "polynomial_schur_locking_report.md").write_text(make_report(payload), encoding="utf-8")
    print("Polynomial Schur-locking audit")
    print(f"  max driver F residual: {max_driver_f:.3e}")
    print(f"  max non-driver F residual: {max_non_driver_f:.3e}")
    print(f"  3x2 eps=0.30 Wilson shift: {benchmark['pathsum_wilson_proxy_shift']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
