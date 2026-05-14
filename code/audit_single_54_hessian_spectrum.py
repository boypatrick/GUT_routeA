#!/usr/bin/env python3
"""Single-54_H Hessian, threshold, and proton/RGE replay audit.

No web lookup is used.  This is the follow-up to the F54 orbit-stabilizer
audit.  It expands the renormalizable one-field breaking superpotential

    W54 = m/2 Tr S^2 + lambda/3 Tr S^3

around S = v S0 with S0=diag(-2^6,3^4), v=-m/lambda, and diagonalizes the
full 54-dimensional traceless-symmetric Hessian.  It then translates the
physical non-Goldstone masses into one-loop MSSM threshold vectors and feeds
them into the existing corrected two-loop/proton scan.

The key convention is

    kappa_54 = M_phys[(20',1,1) and (1,3,3)] / M_G = 5 |m| / M_G.

The radial singlet has mass |m| = kappa_54 M_G/5 and no SM gauge threshold.
The (6,2,2) modes are exact chiral Goldstones of Spin(10)->Pati-Salam and
are eaten by massive vector multiplets at M_G.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_corrected_downstream as down  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "single_54_hessian"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
CORRECTED_SCAN_CSV = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_corrected_scan.csv"

S0 = np.diag([-2.0] * 6 + [3.0] * 4)
N = 10
M = 1.0
LAM = 1.0
V = -M / LAM

# Chiral superfield one-loop MSSM b-vectors with GUT-normalized hypercharge.
B_20P = np.array([16.0 / 5.0, 0.0, 8.0], dtype=float)
B_133 = np.array([18.0 / 5.0, 6.0, 0.0], dtype=float)
B_NON_GOLDSTONE = B_20P + B_133
B_GOLDSTONE_622 = np.array([26.0 / 5.0, 6.0, 4.0], dtype=float)
B_FULL_54 = np.array([12.0, 12.0, 12.0], dtype=float)

KAPPA_GRID = [
    0.5,
    0.75,
    0.9,
    0.95,
    0.99,
    1.0,
    1.01,
    1.05,
    1.1,
    1.25,
    1.5,
    2.0,
    3.0,
    5.0,
]


def symmetric_traceless_basis() -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    # Off-diagonal symmetric directions.
    for i in range(N):
        for j in range(i + 1, N):
            mat = np.zeros((N, N), dtype=float)
            mat[i, j] = 1.0
            mat[j, i] = 1.0
            basis.append(mat)
    # Diagonal traceless directions e_i - e_10.
    for i in range(N - 1):
        mat = np.zeros((N, N), dtype=float)
        mat[i, i] = 1.0
        mat[N - 1, N - 1] = -1.0
        basis.append(mat)
    assert len(basis) == 54
    return basis


def coords_from_matrix(mat: np.ndarray) -> np.ndarray:
    coords: list[float] = []
    for i in range(N):
        for j in range(i + 1, N):
            coords.append(float(mat[i, j]))
    diag = np.diag(mat)
    # In the chosen diagonal basis, mat_ii = c_i for i<9 and
    # mat_10,10 = -sum_i c_i.  The map preserves tracelessness.
    coords.extend(float(diag[i]) for i in range(N - 1))
    return np.array(coords, dtype=float)


def hessian_action(phi: np.ndarray) -> np.ndarray:
    """Linearized F-term on traceless-symmetric fluctuations."""
    background = V * S0
    trace_term = 2.0 * np.trace(background @ phi) / N
    return M * phi + LAM * (background @ phi + phi @ background - trace_term * np.eye(N))


def hessian_matrix() -> np.ndarray:
    basis = symmetric_traceless_basis()
    mat = np.zeros((len(basis), len(basis)), dtype=float)
    for col, phi in enumerate(basis):
        mat[:, col] = coords_from_matrix(hessian_action(phi))
    return mat


def classify_eigenvalue(value: float) -> str:
    if abs(value) < 1.0e-9:
        return "Goldstone_(6,2,2)"
    if abs(value - 5.0) < 1.0e-9:
        return "PS_(20prime,1,1)"
    if abs(value + 5.0) < 1.0e-9:
        return "PS_(1,3,3)"
    if abs(value + 1.0) < 1.0e-9:
        return "radial_singlet"
    return "unclassified"


def spectrum_summary() -> dict[str, Any]:
    mat = hessian_matrix()
    eigvals = np.linalg.eigvals(mat)
    eigvals_real = sorted(float(np.real_if_close(value).real) for value in eigvals)
    counts = Counter(classify_eigenvalue(value) for value in eigvals_real)
    return {
        "eigenvalues": eigvals_real,
        "classification_counts": dict(counts),
        "unique_eigenvalues": {
            key: count
            for key, count in sorted(
                Counter(round(value, 10) for value in eigvals_real).items(),
                key=lambda item: item[0],
            )
        },
        "max_imag_abs": float(max(abs(np.imag(value)) for value in eigvals)),
        "matrix_rank": int(np.linalg.matrix_rank(mat)),
        "zero_modes": int(counts["Goldstone_(6,2,2)"]),
        "passes": counts
        == Counter(
            {
                "PS_(20prime,1,1)": 20,
                "PS_(1,3,3)": 9,
                "Goldstone_(6,2,2)": 24,
                "radial_singlet": 1,
            }
        ),
    }


def load_vacuum_summary() -> dict[str, Any]:
    return json.loads(VACUUM.read_text(encoding="utf-8"))


def mediator_delta_for_r200() -> np.ndarray:
    payload = load_vacuum_summary()
    return np.array(
        payload["benchmark_cards"]["200.0"]["mediator_heavy_threshold"]["delta_full"],
        dtype=float,
    )


def iter_corrected_rows() -> list[dict[str, str]]:
    with CORRECTED_SCAN_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def delta_54(kappa_54: float) -> np.ndarray:
    return B_NON_GOLDSTONE * math.log(1.0 / kappa_54) / (2.0 * math.pi)


def replay_kappa_grid() -> list[dict[str, Any]]:
    source_rows = iter_corrected_rows()
    med_delta = mediator_delta_for_r200()
    rows: list[dict[str, Any]] = []
    for kappa in KAPPA_GRID:
        d54 = delta_54(kappa)
        scan = down.scan_cached_with_delta(source_rows, med_delta + d54)
        best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
        rows.append(
            {
                "kappa_54": kappa,
                "M_nonGoldstone_over_MG": kappa,
                "M_radial_over_MG": kappa / 5.0,
                "log_MG_over_M54": math.log(1.0 / kappa),
                "delta54": d54.tolist(),
                "projected_delta54": (base.PROJECTOR @ d54).tolist(),
                "projected_l2_delta54": float(np.linalg.norm(base.PROJECTOR @ d54)),
                "total_projected_l2_with_mediator": float(np.linalg.norm(base.PROJECTOR @ (med_delta + d54))),
                "safe_points": int(scan["safe_points"]),
                "safe_single_scale_factor2_points": int(scan["safe_single_scale_factor2_points"]),
                "best_is_safe": scan["best_safe"] is not None,
                "best": compact_best(best),
            }
        )
    return rows


def compact_best(row: dict[str, float | bool]) -> dict[str, Any]:
    keys = [
        "score",
        "tan_beta",
        "MSUSY_GeV",
        "MG_GeV",
        "alphaG_inv",
        "lambda_T",
        "lambda_S",
        "chi",
        "kappa_3",
        "kappa_8",
        "log_HC",
        "log_Sigma3",
        "log_Sigma8",
        "M_HC_GeV",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "residual_l2_after_mediator",
    ]
    out: dict[str, Any] = {}
    for key in keys:
        if key in row:
            value = row[key]
            out[key] = bool(value) if isinstance(value, bool) else float(value)
    return out


def threshold_tolerances() -> dict[str, Any]:
    unit_projected_l2 = float(np.linalg.norm(base.PROJECTOR @ B_NON_GOLDSTONE) / (2.0 * math.pi))
    budgets = [5.0e-4, 1.0e-3, 1.0e-2, 5.0e-2, 1.0e-1]
    return {
        "B_20prime": B_20P.tolist(),
        "B_133": B_133.tolist(),
        "B_non_goldstone": B_NON_GOLDSTONE.tolist(),
        "B_goldstone_622": B_GOLDSTONE_622.tolist(),
        "B_full_54_check": (B_NON_GOLDSTONE + B_GOLDSTONE_622).tolist(),
        "B_full_54_expected": B_FULL_54.tolist(),
        "unit_projected_l2_per_abs_log_kappa": unit_projected_l2,
        "log_kappa_allowed_by_budget": {
            f"{budget:.1e}": budget / unit_projected_l2 for budget in budgets
        },
        "kappa_range_allowed_by_budget": {
            f"{budget:.1e}": [
                math.exp(-budget / unit_projected_l2),
                math.exp(budget / unit_projected_l2),
            ]
            for budget in budgets
        },
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "kappa_54",
        "M_nonGoldstone_over_MG",
        "M_radial_over_MG",
        "log_MG_over_M54",
        "projected_l2_delta54",
        "total_projected_l2_with_mediator",
        "safe_points",
        "safe_single_scale_factor2_points",
        "best_is_safe",
        "alphaG_inv",
        "M_Sigma3_GeV",
        "M_Sigma8_GeV",
        "tau_dim6_years",
        "tau_dim5_target_filter_years",
        "triplet_filter_required",
        "log_HC",
        "log_Sigma3",
        "log_Sigma8",
    ]
    with (OUT / "single_54_threshold_replay.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            best = row["best"]
            writer.writerow(
                {
                    "kappa_54": row["kappa_54"],
                    "M_nonGoldstone_over_MG": row["M_nonGoldstone_over_MG"],
                    "M_radial_over_MG": row["M_radial_over_MG"],
                    "log_MG_over_M54": row["log_MG_over_M54"],
                    "projected_l2_delta54": row["projected_l2_delta54"],
                    "total_projected_l2_with_mediator": row["total_projected_l2_with_mediator"],
                    "safe_points": row["safe_points"],
                    "safe_single_scale_factor2_points": row["safe_single_scale_factor2_points"],
                    "best_is_safe": row["best_is_safe"],
                    "alphaG_inv": best["alphaG_inv"],
                    "M_Sigma3_GeV": best["M_Sigma3_GeV"],
                    "M_Sigma8_GeV": best["M_Sigma8_GeV"],
                    "tau_dim6_years": best["tau_dim6_years"],
                    "tau_dim5_target_filter_years": best["tau_dim5_target_filter_years"],
                    "triplet_filter_required": best["triplet_filter_required"],
                    "log_HC": best["log_HC"],
                    "log_Sigma3": best["log_Sigma3"],
                    "log_Sigma8": best["log_Sigma8"],
                }
            )


def build_payload() -> dict[str, Any]:
    spectrum = spectrum_summary()
    replay = replay_kappa_grid()
    tolerances = threshold_tolerances()
    row_k1 = next(row for row in replay if abs(row["kappa_54"] - 1.0) < 1.0e-12)
    return {
        "note": "No web lookup used. Single-54_H Hessian and threshold replay audit.",
        "superpotential": "W54 = m/2 Tr S^2 + lambda/3 Tr S^3, S=v S0, v=-m/lambda",
        "S0": [-2.0] * 6 + [3.0] * 4,
        "hessian_formula": (
            "delta F = m phi + lambda[v(S0 phi + phi S0) "
            "- 2 v Tr(S0 phi) I/10]"
        ),
        "spectrum": spectrum,
        "fragment_summary": {
            "(20prime,1,1)": {"multiplicity": 20, "holomorphic_eigenvalue": "+5 m", "mass_abs_over_abs_m": 5.0, "b_vector": B_20P.tolist()},
            "(1,3,3)": {"multiplicity": 9, "holomorphic_eigenvalue": "-5 m", "mass_abs_over_abs_m": 5.0, "b_vector": B_133.tolist()},
            "(6,2,2)": {"multiplicity": 24, "holomorphic_eigenvalue": "0", "role": "eaten Goldstone", "b_vector_if_not_eaten": B_GOLDSTONE_622.tolist()},
            "(1,1,1)": {"multiplicity": 1, "holomorphic_eigenvalue": "-m", "mass_abs_over_abs_m": 1.0, "b_vector": [0.0, 0.0, 0.0]},
        },
        "threshold_convention": {
            "kappa_54": "M[(20prime,1,1) and (1,3,3)] / M_G = 5 |m| / M_G",
            "delta54": "B_non_goldstone log(1/kappa_54)/(2 pi)",
            "radial_singlet": "M_radial/M_G = kappa_54/5; no SM one-loop threshold",
        },
        "threshold_tolerances": tolerances,
        "kappa_replay_rows": replay,
        "verdict": {
            "hessian_passes": bool(spectrum["passes"]),
            "kappa54_equal_1_branch_safe": bool(row_k1["best_is_safe"]),
            "kappa54_equal_1_delta54_l2": row_k1["projected_l2_delta54"],
            "nontrivial_element": (
                "Choose the cubic mass parameter so that 5|m|=M_G. Then all "
                "non-Goldstone non-singlet 54_H fragments sit at M_G, the radial "
                "singlet lies at M_G/5 but is SM-neutral, and the 24 zero modes are "
                "eaten. The F54 order parameter adds no one-loop SM threshold while "
                "remaining a one-54_H, low-Dynkin-index realization."
            ),
            "if_split": (
                "If kappa_54 differs from one, the active threshold vector is "
                "(34/5,6,8) log(1/kappa_54)/(2 pi). It can be replayed through "
                "the existing two-loop/proton scan using the output kappa grid."
            ),
        },
    }


def write_report(payload: dict[str, Any]) -> None:
    spectrum = payload["spectrum"]
    tol = payload["threshold_tolerances"]
    lines: list[str] = []
    lines.append("# Single-54_H Hessian and threshold replay audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Hessian spectrum")
    lines.append("")
    lines.append("For `W54 = m/2 Tr S^2 + lambda/3 Tr S^3` at `S=v S0`, `v=-m/lambda`,")
    lines.append("the full 54-dimensional traceless-symmetric Hessian has:")
    lines.append("")
    lines.append("| fragment | multiplicity | holomorphic eigenvalue | role |")
    lines.append("|---|---:|---:|---|")
    lines.append("| (20',1,1) | 20 | +5 m | physical non-singlet |")
    lines.append("| (1,3,3) | 9 | -5 m | physical non-singlet |")
    lines.append("| (6,2,2) | 24 | 0 | eaten Goldstone |")
    lines.append("| (1,1,1) | 1 | -m | radial singlet |")
    lines.append("")
    lines.append(f"Matrix rank: {spectrum['matrix_rank']}; zero modes: {spectrum['zero_modes']}.")
    lines.append("")
    lines.append("## Threshold vector")
    lines.append("")
    lines.append("The active non-Goldstone one-loop vector is")
    lines.append("")
    lines.append("```text")
    lines.append("B_54_phys = B_(20',1,1)+B_(1,3,3) = (34/5, 6, 8)")
    lines.append("Delta_54 = B_54_phys log(1/kappa_54)/(2 pi)")
    lines.append("```")
    lines.append("")
    lines.append(
        f"`||P Delta_54||_2 = {tol['unit_projected_l2_per_abs_log_kappa']:.9e} |log kappa_54|`."
    )
    lines.append("")
    lines.append("## RGE/proton replay")
    lines.append("")
    lines.append("| kappa_54 | ||P Delta54|| | safe points | alphaG^-1 | M_Sigma3 [GeV] | tau_d6 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for row in payload["kappa_replay_rows"]:
        best = row["best"]
        if row["kappa_54"] in {0.5, 0.9, 1.0, 1.1, 2.0, 5.0}:
            lines.append(
                f"| {row['kappa_54']:.2f} | {row['projected_l2_delta54']:.6e} | "
                f"{row['safe_points']} | {best['alphaG_inv']:.6f} | "
                f"{best['M_Sigma3_GeV']:.6e} | {best['tau_dim6_years']:.6e} |"
            )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"]["nontrivial_element"])
    lines.append(payload["verdict"]["if_split"])
    lines.append("")
    (OUT / "single_54_hessian_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    write_csv(payload["kappa_replay_rows"])
    (OUT / "single_54_hessian_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Single-54_H Hessian audit complete")
    print(f"  spectrum_passes={payload['spectrum']['passes']}")
    print(f"  spectrum_counts={payload['spectrum']['classification_counts']}")
    print(
        "  unit ||P Delta54|| per |log kappa|="
        f"{payload['threshold_tolerances']['unit_projected_l2_per_abs_log_kappa']:.9e}"
    )
    print(
        f"  kappa54=1 branch safe={payload['verdict']['kappa54_equal_1_branch_safe']}"
    )
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
