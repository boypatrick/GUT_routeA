#!/usr/bin/env python3
"""Mediator-only G_tr spectrum and threshold audit.

No web lookup is used.  The previous G-parity/Yukawa audit showed that the
extra Z2^G symmetry cannot act on the physical 120_H if the antisymmetric
Clebsch flavor branch is retained.  The clean option is therefore a
mediator-only G_tr copy.  This script asks what that copy may be:

1. a full Spin(10) 120'_H multiplet;
2. a minimal Spin(10) 10'_G multiplet whose 5+5bar contains the triplet pair;
3. a post-Spin(10) Pati-Salam / SM vectorlike triplet pair.

The audit is representation-theoretic and numerical.  Complete SU(5) or
Spin(10) multiplets give a universal one-loop threshold, while a triplet-only
remnant has beta vector

    b_D+Dbar = (2/5, 0, 1),

and therefore produces a non-universal matching vector unless its mass is
locked close to M_G.
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


OUT = ROOT / "output" / "gtr_mediator_spectrum"

KAPPA_GRID = [0.01, 0.03, 0.10, 0.30, 1.0, 3.0, 10.0, 30.0, 100.0]
PROJECTED_TOLERANCES = [1.0e-3, 1.0e-2, 5.4e-2]

# N=1 chiral superfield one-loop beta contributions in MSSM normalization.
B_TRIPLET_PAIR = np.array([2.0 / 5.0, 0.0, 1.0], dtype=float)
B_DOUBLE_PAIR = np.array([3.0 / 5.0, 1.0, 0.0], dtype=float)
B_FIVE_PAIR = B_TRIPLET_PAIR + B_DOUBLE_PAIR
B_TEN_PAIR_SU5 = np.array([3.0, 3.0, 3.0], dtype=float)
B_FORTYFIVE_PAIR_SU5 = np.array([24.0, 24.0, 24.0], dtype=float)
B_FULL_120 = B_FIVE_PAIR + B_TEN_PAIR_SU5 + B_FORTYFIVE_PAIR_SU5


def threshold_delta(beta_vector: np.ndarray, kappa: float) -> np.ndarray:
    """Return Delta_i=b_i log(M_G/M)/(2pi), with kappa=M/M_G."""

    return beta_vector * math.log(1.0 / kappa) / (2.0 * math.pi)


def projected_l2(vec: np.ndarray) -> float:
    return float(np.linalg.norm(base.PROJECTOR @ vec))


def beta_rows() -> list[dict[str, Any]]:
    return [
        {
            "name": "D+Dbar triplet pair",
            "sm_content": "(3,1,-1/3)+(bar3,1,+1/3)",
            "beta_vector": B_TRIPLET_PAIR.tolist(),
            "projected_beta_l2": projected_l2(B_TRIPLET_PAIR),
        },
        {
            "name": "L+Lbar doublet pair",
            "sm_content": "(1,2,+1/2)+(1,2,-1/2)",
            "beta_vector": B_DOUBLE_PAIR.tolist(),
            "projected_beta_l2": projected_l2(B_DOUBLE_PAIR),
        },
        {
            "name": "5+5bar common pair",
            "sm_content": "D+Dbar plus L+Lbar",
            "beta_vector": B_FIVE_PAIR.tolist(),
            "projected_beta_l2": projected_l2(B_FIVE_PAIR),
        },
        {
            "name": "full 120_H",
            "sm_content": "5+5bar+10+10bar+45+45bar",
            "beta_vector": B_FULL_120.tolist(),
            "projected_beta_l2": projected_l2(B_FULL_120),
        },
    ]


def branch_definitions() -> list[dict[str, Any]]:
    return [
        {
            "branch": "spin10_full_120prime_common",
            "interpretation": "A distinct full 120'_H chiral multiplet, Z2^G odd, with all fragments at a common mass.",
            "action_level": "Spin(10)-covariant",
            "beta_vector": B_FULL_120,
            "so10_dynkin_index_increment": 28.0,
            "main_caveat": "Threshold-silent only if complete and nearly degenerate; worsens UV perturbativity by T(120)=28.",
        },
        {
            "branch": "spin10_minimal_10prime_common",
            "interpretation": "A sterile 10'_G chiral multiplet supplies a 5+5bar triplet/doublet pair; G_tr is a mediator label, not the physical 120_H.",
            "action_level": "Spin(10)-covariant after relabelling the regulator source",
            "beta_vector": B_FIVE_PAIR,
            "so10_dynkin_index_increment": 1.0,
            "main_caveat": "Needs a projector/selection rule so only its triplet pair participates in the d=5 regulator.",
        },
        {
            "branch": "su5_completed_5pair_ps_eft",
            "interpretation": "Post-Spin(10) EFT keeps an SU(5)-complete 5+5bar pair with common mass.",
            "action_level": "not a standalone Spin(10) invariant unless embedded in 10'_G",
            "beta_vector": B_FIVE_PAIR,
            "so10_dynkin_index_increment": None,
            "main_caveat": "Acceptable as conditional EFT, but its Spin(10)-breaking origin must be stated.",
        },
        {
            "branch": "ps_triplet_only_pair",
            "interpretation": "Only the vectorlike color-triplet pair remains below the Spin(10) threshold.",
            "action_level": "post-Spin(10) EFT remnant",
            "beta_vector": B_TRIPLET_PAIR,
            "so10_dynkin_index_increment": None,
            "main_caveat": "Non-universal one-loop threshold unless M_D is locked close to M_G.",
        },
        {
            "branch": "split_10prime_triplet_only",
            "interpretation": "A 10'_G exists, but its doublet pair is at M_G while the triplet pair is split.",
            "action_level": "Spin(10)-origin with explicit doublet-triplet split",
            "beta_vector": B_TRIPLET_PAIR,
            "so10_dynkin_index_increment": 1.0,
            "main_caveat": "Same non-universal threshold as a triplet-only remnant; the split itself must be dynamically justified.",
        },
    ]


def replay_delta(delta: np.ndarray, source_rows: list[dict[str, str]]) -> dict[str, Any]:
    scan = corrected.scan_cached_with_delta(source_rows, delta)
    best = scan["best_safe"] if scan["best_safe"] is not None else scan["best"]
    compact = corrected.compact_best(best)
    assert compact is not None
    return {
        "safe_points": int(scan["safe_points"]),
        "safe_single_scale_factor2_points": int(scan["safe_single_scale_factor2_points"]),
        "best_is_safe": scan["best_safe"] is not None,
        "best_safe_or_best": compact,
    }


def scan_thresholds() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    source_rows = corrected.iter_corrected_rows()
    rows: list[dict[str, Any]] = []
    branches = branch_definitions()
    for branch in branches:
        beta = np.array(branch["beta_vector"], dtype=float)
        for kappa in KAPPA_GRID:
            delta = threshold_delta(beta, kappa)
            replay = replay_delta(delta, source_rows)
            best = replay["best_safe_or_best"]
            rows.append(
                {
                    "branch": branch["branch"],
                    "kappa_mass_over_MG": float(kappa),
                    "threshold_delta": delta.tolist(),
                    "projected_delta_l2": projected_l2(delta),
                    "universal_delta": float(np.mean(delta)),
                    "safe_points": replay["safe_points"],
                    "safe_single_scale_factor2_points": replay["safe_single_scale_factor2_points"],
                    "best_is_safe": replay["best_is_safe"],
                    "alphaG_inv": best.get("alphaG_inv"),
                    "M_Sigma3_GeV": best.get("M_Sigma3_GeV"),
                    "M_Sigma8_GeV": best.get("M_Sigma8_GeV"),
                    "tau_dim6_years": best.get("tau_dim6_years"),
                    "tau_dim5_target_filter_years": best.get("tau_dim5_target_filter_years"),
                    "triplet_filter_required": best.get("triplet_filter_required"),
                    "main_caveat": branch["main_caveat"],
                }
            )
    by_branch: dict[str, Any] = {}
    for branch in branches:
        brows = [row for row in rows if row["branch"] == branch["branch"]]
        by_branch[branch["branch"]] = {
            "min_projected_l2": min(row["projected_delta_l2"] for row in brows),
            "max_projected_l2": max(row["projected_delta_l2"] for row in brows),
            "safe_points_at_kappa_1": next(row["safe_points"] for row in brows if row["kappa_mass_over_MG"] == 1.0),
            "safe_points_min": min(row["safe_points"] for row in brows),
            "safe_points_max": max(row["safe_points"] for row in brows),
        }
    return rows, by_branch


def lock_windows() -> list[dict[str, float]]:
    pb = projected_l2(B_TRIPLET_PAIR)
    rows = []
    for tol in PROJECTED_TOLERANCES:
        max_abs_log = 2.0 * math.pi * tol / pb
        rows.append(
            {
                "projected_delta_l2_tolerance": float(tol),
                "projected_beta_l2_triplet_pair": pb,
                "max_abs_log_kappa": float(max_abs_log),
                "kappa_min": float(math.exp(-max_abs_log)),
                "kappa_max": float(math.exp(max_abs_log)),
                "fractional_half_width_approx": float(max_abs_log),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for row in rows for key in row.keys()})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_report(payload: dict[str, Any]) -> str:
    lines = [
        "# Mediator-only G_tr spectrum audit",
        "",
        "No web lookup was used.",
        "",
        "## One-loop beta ledger",
        "",
        "| field content | beta vector | projected beta l2 |",
        "|---|---:|---:|",
    ]
    for row in payload["beta_ledger"]:
        lines.append(
            f"| {row['name']} | `{row['beta_vector']}` | {row['projected_beta_l2']:.6e} |"
        )
    lines.extend(
        [
            "",
            "A complete 5+5bar or full 120_H is universal in the gauge-matching plane.",
            "A triplet-only pair has nonzero projected beta and must be mass-locked.",
            "",
            "## Triplet-only mass-lock windows",
            "",
            "| projected l2 tolerance | kappa min | kappa max | max |log kappa| |",
            "|---:|---:|---:|---:|",
        ]
    )
    for row in payload["triplet_only_lock_windows"]:
        lines.append(
            f"| {row['projected_delta_l2_tolerance']:.3e} | {row['kappa_min']:.6f} | "
            f"{row['kappa_max']:.6f} | {row['max_abs_log_kappa']:.6e} |"
        )
    lines.extend(
        [
            "",
            "## Branch verdict",
            "",
            payload["verdict"],
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    threshold_rows, branch_summary = scan_thresholds()
    beta_ledger = beta_rows()
    windows = lock_windows()
    branches = branch_definitions()

    csv_rows = []
    for row in threshold_rows:
        out = dict(row)
        out["threshold_delta"] = " ".join(f"{x:.12e}" for x in row["threshold_delta"])
        csv_rows.append(out)
    write_csv(OUT / "gtr_threshold_replay.csv", csv_rows)
    write_csv(OUT / "gtr_lock_windows.csv", windows)

    verdict = (
        "The action-level branch with the least new UV cost is not a full 120'_H, "
        "but a sterile Spin(10) 10'_G or an equivalent SU(5)-complete 5+5bar EFT "
        "pair whose doublet partner is degenerate with the triplet.  This keeps "
        "P Delta_Gtr=0 while increasing the SO(10) Dynkin index only by 1 if "
        "embedded as 10'_G.  A triplet-only remnant is usable only as a conditional "
        "post-Spin(10) EFT threshold: for ||P Delta||<1e-2 it requires "
        f"{windows[1]['kappa_min']:.3f}<M_D/M_G<{windows[1]['kappa_max']:.3f}; "
        "otherwise it must be included as a new non-universal threshold in every "
        "RGE/proton scan.  A full 120'_H is threshold-silent when degenerate but "
        "adds T(120)=28 to the UV beta function, so it should be a last-resort "
        "completion rather than the default regulator."
    )

    payload = {
        "note": "No web lookup used. Mediator-only G_tr spectrum and threshold audit.",
        "threshold_formula": "Delta_i = b_i log(M_G/M)/(2 pi), kappa=M/M_G",
        "beta_ledger": beta_ledger,
        "branch_definitions": [
            {**branch, "beta_vector": branch["beta_vector"].tolist()} for branch in branches
        ],
        "triplet_only_lock_windows": windows,
        "threshold_replay_summary_by_branch": branch_summary,
        "threshold_replay_rows": threshold_rows,
        "verdict": verdict,
    }
    (OUT / "gtr_mediator_spectrum_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "gtr_mediator_spectrum_report.md").write_text(make_report(payload), encoding="utf-8")

    print("G_tr mediator spectrum audit")
    print(f"  triplet projected beta l2: {projected_l2(B_TRIPLET_PAIR):.9e}")
    print(
        "  triplet-only window for ||P Delta||<1e-2: "
        f"{windows[1]['kappa_min']:.6f} < M_D/M_G < {windows[1]['kappa_max']:.6f}"
    )
    print("  preferred action-level branch: spin10_minimal_10prime_common")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
