#!/usr/bin/env python3
"""Projector superpotential audit for the sterile 10'_G completion.

No web lookup is used.  The previous audit selected a sterile Spin(10) 10'_G
as the least costly completion of the mediator-only G_tr source.  This script
checks whether the same 10'_G can be made triplet-active but doublet-inert.

The key object is the 54_H order parameter in the SO(10) vector 10:

    S0 = diag(-2,...,-2, +3,+3,+3,+3).

It defines exact vector-space projectors

    P_D = (3 I - S0)/5,       P_L = (S0 + 2 I)/5,

onto the Pati-Salam (6,1,1) color-triplet sector and the (1,2,2) doublet
sector.  The audit then compares three bridge patterns:

1. one-sided triangular bridge: threshold silent but cannot correct the
   physical 120 inverse propagator by itself;
2. generic two-sided bridge: corrects the inverse propagator but generates a
   non-universal triplet determinant threshold;
3. determinant-locked two-sided bridge: adds a projector mass counterterm so
   det M_triplet stays fixed while the inverse propagator changes.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_corrected_downstream as corrected  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
from audit_gtr_mediator_spectrum import B_DOUBLE_PAIR, B_TRIPLET_PAIR, projected_l2  # noqa: E402


OUT = ROOT / "output" / "ten_g_projector"

XI_GRID = [-0.10, -0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05, 0.10]
EPS_GRID = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80]
TOLERANCES = [1.0e-3, 1.0e-2, 5.4e-2]


def s54_vector() -> np.ndarray:
    return np.diag(np.array([-2.0] * 6 + [3.0] * 4, dtype=float))


def projector_audit() -> dict[str, Any]:
    eye = np.eye(10)
    s = s54_vector()
    p_d = (3.0 * eye - s) / 5.0
    p_l = (s + 2.0 * eye) / 5.0
    return {
        "S54_vector_diagonal": np.diag(s).tolist(),
        "P_D_diagonal": np.diag(p_d).tolist(),
        "P_L_diagonal": np.diag(p_l).tolist(),
        "P_D_from_singlet_plus_54": "P_D = (3/5) I - (1/5) S54",
        "P_L_from_singlet_plus_54": "P_L = (2/5) I + (1/5) S54",
        "idempotence_error_P_D": float(np.linalg.norm(p_d @ p_d - p_d)),
        "idempotence_error_P_L": float(np.linalg.norm(p_l @ p_l - p_l)),
        "orthogonality_error": float(np.linalg.norm(p_d @ p_l)),
        "completeness_error": float(np.linalg.norm(p_d + p_l - eye)),
        "rank_P_D": int(round(np.trace(p_d))),
        "rank_P_L": int(round(np.trace(p_l))),
    }


def delta_for_split_mass(xi: float) -> np.ndarray | None:
    """Mass term M 10_G 10bar_G + xi M 10_G S54 10bar_G."""

    m_d = 1.0 - 2.0 * xi
    m_l = 1.0 + 3.0 * xi
    if m_d <= 0.0 or m_l <= 0.0:
        return None
    return (B_TRIPLET_PAIR * math.log(1.0 / m_d) + B_DOUBLE_PAIR * math.log(1.0 / m_l)) / (2.0 * math.pi)


def replay_delta(delta: np.ndarray, source_rows: list[dict[str, str]]) -> dict[str, Any]:
    scan = corrected.scan_cached_with_delta(source_rows, delta)
    best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
    compact = corrected.compact_best(best)
    assert compact is not None
    return {
        "safe_points": int(scan["safe_points"]),
        "best_is_safe": scan["best_safe"] is not None,
        "alphaG_inv": compact["alphaG_inv"],
        "M_Sigma3_GeV": compact["M_Sigma3_GeV"],
        "tau_dim6_years": compact["tau_dim6_years"],
        "tau_dim5_target_filter_years": compact["tau_dim5_target_filter_years"],
    }


def mass_split_scan(source_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    rows = []
    for xi in XI_GRID:
        m_d = 1.0 - 2.0 * xi
        m_l = 1.0 + 3.0 * xi
        delta = delta_for_split_mass(xi)
        if delta is None:
            rows.append(
                {
                    "xi": xi,
                    "M_D_over_M": m_d,
                    "M_L_over_M": m_l,
                    "valid_positive_masses": False,
                }
            )
            continue
        replay = replay_delta(delta, source_rows)
        rows.append(
            {
                "xi": xi,
                "M_D_over_M": m_d,
                "M_L_over_M": m_l,
                "threshold_delta": delta.tolist(),
                "projected_delta_l2": projected_l2(delta),
                **replay,
                "valid_positive_masses": True,
            }
        )
    return rows


def max_abs_xi_for_tolerance(tol: float) -> float:
    """Binary-search the symmetric interval around xi=0 allowed by P Delta."""

    def ok(x: float) -> bool:
        dpos = delta_for_split_mass(x)
        dneg = delta_for_split_mass(-x)
        if dpos is None or dneg is None:
            return False
        return max(projected_l2(dpos), projected_l2(dneg)) <= tol

    lo, hi = 0.0, 0.30
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def inverse_metrics(mat: np.ndarray) -> dict[str, float]:
    inv = np.linalg.inv(mat)
    det = float(np.linalg.det(mat))
    singulars = np.linalg.svd(mat, compute_uv=False)
    return {
        "determinant": det,
        "singular_min": float(singulars[-1]),
        "singular_max": float(singulars[0]),
        "inverse_AA": float(inv[0, 0]),
        "inverse_AG": float(inv[0, 1]),
        "inverse_GA": float(inv[1, 0]),
        "inverse_GG": float(inv[1, 1]),
    }


def bridge_matrix(kind: str, eps: float) -> np.ndarray:
    if kind == "triangular_one_sided":
        return np.array([[1.0, eps], [0.0, 1.0]], dtype=float)
    if kind == "two_sided_equal":
        return np.array([[1.0, eps], [eps, 1.0]], dtype=float)
    if kind == "determinant_locked_equal":
        return np.array([[1.0, eps], [eps, 1.0 + eps * eps]], dtype=float)
    raise ValueError(kind)


def bridge_delta_from_det(det: float) -> np.ndarray:
    return B_TRIPLET_PAIR * math.log(1.0 / abs(det)) / (2.0 * math.pi)


def bridge_scan(source_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    rows = []
    for kind in ["triangular_one_sided", "two_sided_equal", "determinant_locked_equal"]:
        for eps in EPS_GRID:
            mat = bridge_matrix(kind, eps)
            metrics = inverse_metrics(mat)
            delta = bridge_delta_from_det(metrics["determinant"])
            replay = replay_delta(delta, source_rows)
            rows.append(
                {
                    "bridge_kind": kind,
                    "epsilon": eps,
                    **metrics,
                    "threshold_delta": delta.tolist(),
                    "projected_delta_l2": projected_l2(delta),
                    **replay,
                    "changes_inverse_AA": abs(metrics["inverse_AA"] - 1.0) > 1.0e-12,
                    "threshold_silent": projected_l2(delta) < 1.0e-12,
                }
            )
    return rows


def threshold_tolerances() -> list[dict[str, float]]:
    rows = []
    for tol in TOLERANCES:
        xi = max_abs_xi_for_tolerance(tol)
        # Two-sided equal has det=1-eps^2 and P Delta = pb log(1/(1-eps^2))/(2pi).
        pb = projected_l2(B_TRIPLET_PAIR)
        max_log = 2.0 * math.pi * tol / pb
        eps = math.sqrt(max(0.0, 1.0 - math.exp(-max_log)))
        rows.append(
            {
                "projected_delta_l2_tolerance": tol,
                "max_abs_xi_54_mass_split": xi,
                "max_equal_two_sided_epsilon_without_det_lock": eps,
                "max_log_det_shift": max_log,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for row in rows for key in row.keys()})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            if "threshold_delta" in out and isinstance(out["threshold_delta"], list):
                out["threshold_delta"] = " ".join(f"{x:.12e}" for x in out["threshold_delta"])
            writer.writerow(out)


def make_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Sterile 10G projector superpotential audit",
        "",
        "No web lookup was used.",
        "",
        "## Vector projectors",
        "",
        "`P_D=(3I-S54)/5`, `P_L=(S54+2I)/5`.",
        "",
        f"`rank(P_D)={payload['projector']['rank_P_D']}`, `rank(P_L)={payload['projector']['rank_P_L']}`.",
        f"`||P_D^2-P_D||={payload['projector']['idempotence_error_P_D']:.3e}`.",
        "",
        "## Threshold tolerances",
        "",
        "| tolerance | max |xi| in 54 mass split | max equal bridge epsilon without determinant lock |",
        "|---:|---:|---:|",
    ]
    for row in payload["threshold_tolerances"]:
        lines.append(
            f"| {row['projected_delta_l2_tolerance']:.3e} | "
            f"{row['max_abs_xi_54_mass_split']:.6e} | "
            f"{row['max_equal_two_sided_epsilon_without_det_lock']:.6e} |"
        )
    lines.extend(
        [
            "",
            "## Verdict",
            "",
            payload["verdict"],
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_rows = corrected.iter_corrected_rows()
    projector = projector_audit()
    split_rows = mass_split_scan(source_rows)
    bridge_rows = bridge_scan(source_rows)
    tolerances = threshold_tolerances()

    write_csv(OUT / "ten_g_mass_split_scan.csv", split_rows)
    write_csv(OUT / "ten_g_bridge_scan.csv", bridge_rows)
    write_csv(OUT / "ten_g_threshold_tolerances.csv", tolerances)

    tri_03 = next(row for row in bridge_rows if row["bridge_kind"] == "triangular_one_sided" and abs(row["epsilon"] - 0.30) < 1e-12)
    two_03 = next(row for row in bridge_rows if row["bridge_kind"] == "two_sided_equal" and abs(row["epsilon"] - 0.30) < 1e-12)
    det_03 = next(row for row in bridge_rows if row["bridge_kind"] == "determinant_locked_equal" and abs(row["epsilon"] - 0.30) < 1e-12)

    verdict = (
        "A 54_H insertion gives exact triplet/doublet projectors in the vector 10.  "
        "However the same 54_H insertion must not appear as an unconstrained 10'_G mass term: "
        f"keeping ||P Delta||<1e-2 requires |xi|<{tolerances[1]['max_abs_xi_54_mass_split']:.4f}. "
        "The clean superpotential uses a universal 10'_G mass plus projected triplet bridges. "
        "A one-sided triangular bridge is determinant-silent but cannot by itself change the physical "
        "120 inverse propagator.  A generic two-sided bridge changes the propagator but produces a "
        f"non-universal threshold; at epsilon=0.3 it gives ||P Delta||={two_03['projected_delta_l2']:.4e}. "
        "The nontrivial viable structure is a determinant-locked two-sided bridge with triplet mass "
        "counterterm m_G=1+epsilon_L epsilon_R.  At epsilon=0.3 it changes the AA inverse to "
        f"{det_03['inverse_AA']:.3f}, keeps det M_D={det_03['determinant']:.3f}, and has "
        f"||P Delta||={det_03['projected_delta_l2']:.1e}.  This should be promoted to a small "
        "F-term driver enforcing det M_triplet = M_G^2 while preserving the doublet mass."
    )

    payload = {
        "note": "No web lookup used. Projector superpotential audit for sterile Spin(10) 10'_G.",
        "projector": projector,
        "superpotential_skeleton": {
            "universal_mass": "W0 = M_G X_10 Xbar_10",
            "triplet_projector_bridge": "Wbr = eps_R A_10 P_D Xbar_10 + eps_L X_10 P_D Abar_10",
            "determinant_lock": "Wlock imposes m_G^D = 1 + eps_L eps_R, while m_G^L = 1",
            "warning": "Do not use a generic X_10 S54 Xbar_10 mass unless its coefficient xi is below the threshold tolerance.",
        },
        "mass_split_rows": split_rows,
        "bridge_rows": bridge_rows,
        "threshold_tolerances": tolerances,
        "benchmark_bridge_epsilon_0p3": {
            "triangular_one_sided": tri_03,
            "two_sided_equal": two_03,
            "determinant_locked_equal": det_03,
        },
        "verdict": verdict,
    }
    (OUT / "ten_g_projector_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "ten_g_projector_report.md").write_text(make_report(payload), encoding="utf-8")

    print("10'_G projector superpotential audit")
    print(f"  projector idempotence error: {projector['idempotence_error_P_D']:.3e}")
    print(f"  |xi| max for ||P Delta||<1e-2: {tolerances[1]['max_abs_xi_54_mass_split']:.6e}")
    print(f"  two-sided eps=0.3 projected l2: {two_03['projected_delta_l2']:.6e}")
    print(f"  det-locked eps=0.3 projected l2: {det_03['projected_delta_l2']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
