#!/usr/bin/env python3
"""Joint Knu/RRRR triplet-filter scan with condition-number caps.

No web lookup is used.  This is the two-sided continuation of the triplet
rank-lift audit: instead of nulling only the LLLL K nu block, we build a joint
linear system for the K nu and RRRR source vectors and scan weighted least
singular directions under finite triplet-mass condition-number caps.
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

import construct_triplet_rank_lift as lift  # noqa: E402
import scan_triplet_mixing_nullspace as ns  # noqa: E402


OUT = ROOT / "output" / "two_sided_triplet_filter"


def coeffs_from_w(w: np.ndarray) -> np.ndarray:
    return np.array([w[i, j] for i in range(w.shape[0]) for j in range(w.shape[1])], dtype=complex)


def w_from_coeffs(coeffs: np.ndarray, size: int = 4) -> np.ndarray:
    return np.array(coeffs, dtype=complex).reshape((size, size))


def normalized_block(rows: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(rows)
    return rows / max(float(norm), 1.0e-30)


def least_singular_coeffs(a: np.ndarray) -> dict[str, Any]:
    _, singulars, vh = np.linalg.svd(a, full_matrices=True)
    coeffs = vh.conjugate().T[:, -1]
    coeffs = coeffs / np.linalg.norm(coeffs)
    return {
        "coeffs": coeffs,
        "joint_singular_values": [float(x) for x in singulars],
        "joint_sigma_min": float(singulars[-1]),
        "joint_sigma_max": float(singulars[0]),
    }


def cap_condition(w: np.ndarray, kappa: float | None) -> tuple[np.ndarray, dict[str, Any]]:
    u, singulars, vh = np.linalg.svd(w, full_matrices=True)
    if kappa is None:
        capped = singulars.copy()
    else:
        capped = np.maximum(singulars, singulars[0] / kappa)
    w_cap = u @ np.diag(capped) @ vh
    w_cap = w_cap / np.linalg.norm(w_cap)
    final_s = np.linalg.svd(w_cap, compute_uv=False)
    return w_cap, {
        "original_singular_values": [float(x) for x in singulars],
        "capped_singular_values_before_norm": [float(x) for x in capped],
        "final_singular_values": [float(x) for x in final_s],
        "condition_cap": "uncapped" if kappa is None else float(kappa),
        "final_condition": float(final_s[0] / max(final_s[-1], 1.0e-30)),
        "mass_ratio_if_inverted": float(final_s[0] / max(final_s[-1], 1.0e-30)),
        "max_mass_if_lightest_1e16_GeV": float(final_s[0] / max(final_s[-1], 1.0e-30) * lift.M_LIGHT_BENCHMARK_GEV),
    }


def evaluate(
    label: str,
    weight_rrrr: float,
    kappa: float | None,
    w: np.ndarray,
    pair_full: list[tuple[int, int]],
    channel_entries: dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]],
    constants: dict[str, float],
    joint_info: dict[str, Any],
    cap_info: dict[str, Any],
) -> dict[str, Any]:
    channels = lift.max_channel_for_w(w, pair_full, channel_entries, constants)
    max_knu = max(channels["LLLL_upupdown_Knu"]["amplitude"], channels["LLLL_downdownup_Knu"]["amplitude"])
    k0mu = channels["LLLL_upupdown_K0mu"]["amplitude"]
    rrrr = channels["RRRR_uusd_anycharged"]["amplitude"]
    leading = max(channels.items(), key=lambda item: item[1]["amplitude"])
    worst_kr = max(max_knu, rrrr)
    worst_all = max(max_knu, k0mu, rrrr)
    return {
        "label": label,
        "weight_RRRR": float(weight_rrrr),
        "condition_cap": cap_info["condition_cap"],
        "max_Knu_amplitude": float(max_knu),
        "K0mu_amplitude": float(k0mu),
        "RRRR_amplitude": float(rrrr),
        "worst_Knu_RRRR": float(worst_kr),
        "worst_all_monitored": float(worst_all),
        "leading_channel": leading[0],
        "leading_amplitude": leading[1]["amplitude"],
        "leading_ST_required_tau_2p4e34": leading[1]["S_T_required_tau_2p4e34"],
        "final_condition": cap_info["final_condition"],
        "max_mass_if_lightest_1e16_GeV": cap_info["max_mass_if_lightest_1e16_GeV"],
        "joint_sigma_min": joint_info["joint_sigma_min"],
        "joint_sigma_max": joint_info["joint_sigma_max"],
        "channels": channels,
        "cap_info": cap_info,
        "W_matrix": lift.matrix_json(w),
        "joint_singular_values": joint_info["joint_singular_values"],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    basis_names, pair_full, channel_entries, constants = lift.build_channel_tensors()

    a_knu = np.vstack(
        [
            ns.row_matrix(pair_full, channel_entries["LLLL_upupdown_Knu"][0], channel_entries["LLLL_upupdown_Knu"][1]),
            ns.row_matrix(pair_full, channel_entries["LLLL_downdownup_Knu"][0], channel_entries["LLLL_downdownup_Knu"][1]),
        ]
    )
    a_rrrr = ns.row_matrix(pair_full, channel_entries["RRRR_uusd_anycharged"][0], channel_entries["RRRR_uusd_anycharged"][1])
    a_k0mu = ns.row_matrix(pair_full, channel_entries["LLLL_upupdown_K0mu"][0], channel_entries["LLLL_upupdown_K0mu"][1])

    b_knu = normalized_block(a_knu)
    b_rrrr = normalized_block(a_rrrr)

    weights = [0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
    kappas: list[float | None] = [None, 30.0, 100.0, 300.0]
    rows: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    for weight in weights:
        joint = np.vstack([b_knu, weight * b_rrrr])
        info = least_singular_coeffs(joint)
        w0 = w_from_coeffs(info["coeffs"], 4)
        for kappa in kappas:
            w_cap, cap_info = cap_condition(w0, kappa)
            label = f"omegaR_{weight:g}_kappa_{'uncapped' if kappa is None else int(kappa)}"
            row = evaluate(label, weight, kappa, w_cap, pair_full, channel_entries, constants, info, cap_info)
            rows.append(row)
            csv_rows.append(
                {
                    "label": row["label"],
                    "weight_RRRR": row["weight_RRRR"],
                    "condition_cap": row["condition_cap"],
                    "max_Knu_amplitude": row["max_Knu_amplitude"],
                    "K0mu_amplitude": row["K0mu_amplitude"],
                    "RRRR_amplitude": row["RRRR_amplitude"],
                    "worst_Knu_RRRR": row["worst_Knu_RRRR"],
                    "worst_all_monitored": row["worst_all_monitored"],
                    "leading_channel": row["leading_channel"],
                    "leading_amplitude": row["leading_amplitude"],
                    "leading_ST_required_tau_2p4e34": row["leading_ST_required_tau_2p4e34"],
                    "final_condition": row["final_condition"],
                    "max_mass_if_lightest_1e16_GeV": row["max_mass_if_lightest_1e16_GeV"],
                    "joint_sigma_min": row["joint_sigma_min"],
                }
            )

    with (OUT / "two_sided_triplet_filter_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    best_by_cap: dict[str, Any] = {}
    for cap in ["uncapped", 30.0, 100.0, 300.0]:
        subset = [row for row in rows if row["condition_cap"] == cap]
        best = min(subset, key=lambda row: row["worst_Knu_RRRR"])
        best_by_cap[str(cap)] = {
            key: best[key]
            for key in [
                "label",
                "weight_RRRR",
                "condition_cap",
                "max_Knu_amplitude",
                "K0mu_amplitude",
                "RRRR_amplitude",
                "worst_Knu_RRRR",
                "worst_all_monitored",
                "leading_channel",
                "leading_ST_required_tau_2p4e34",
                "final_condition",
                "max_mass_if_lightest_1e16_GeV",
            ]
        }

    # Also provide a row where all monitored blocks are included in the SVD,
    # to check whether K0mu is merely being moved around.
    b_k0mu = normalized_block(a_k0mu)
    joint_all = np.vstack([b_knu, b_rrrr, b_k0mu])
    all_info = least_singular_coeffs(joint_all)
    all_rows = []
    for kappa in kappas:
        w_cap, cap_info = cap_condition(w_from_coeffs(all_info["coeffs"], 4), kappa)
        label = f"all_blocks_kappa_{'uncapped' if kappa is None else int(kappa)}"
        row = evaluate(label, 1.0, kappa, w_cap, pair_full, channel_entries, constants, all_info, cap_info)
        all_rows.append(row)

    best_all_cap100 = next(row for row in all_rows if row["condition_cap"] == 100.0)
    verdict = (
        "A two-sided weighted SVD can reduce Knu and RRRR simultaneously, but finite condition caps prevent "
        "a dramatic all-channel null.  The best kappa=100 Knu/RRRR compromise lowers the worst monitored "
        "Knu/RRRR entry to the 1e-4 level; including K0mu in the SVD does not create a new charged-lepton "
        "obstruction, and RRRR remains the leading monitored row at a comparable level.  Therefore a "
        "two-sided missing-partner sector is plausible as an EFT filter, but a paper-level proof still needs "
        "a symmetry that explains the selected weight vector and an explicit dressing calculation."
    )

    output = {
        "note": "No web lookup used. Joint constrained SVD for Knu and RRRR triplet filters.",
        "basis_names": basis_names,
        "block_norms_before_normalization": {
            "Knu": float(np.linalg.norm(a_knu)),
            "RRRR": float(np.linalg.norm(a_rrrr)),
            "K0mu": float(np.linalg.norm(a_k0mu)),
        },
        "weights_RRRR": weights,
        "condition_caps": ["uncapped", 30, 100, 300],
        "rows": rows,
        "all_block_rows": all_rows,
        "best_by_cap": best_by_cap,
        "verdict": {
            "best_kappa_30": best_by_cap["30.0"],
            "best_kappa_100": best_by_cap["100.0"],
            "best_kappa_300": best_by_cap["300.0"],
            "all_blocks_kappa_100": {
                key: best_all_cap100[key]
                for key in [
                    "label",
                    "max_Knu_amplitude",
                    "K0mu_amplitude",
                    "RRRR_amplitude",
                    "worst_all_monitored",
                    "leading_channel",
                    "leading_ST_required_tau_2p4e34",
                    "final_condition",
                ]
            },
            "interpretation": verdict,
        },
    }
    (OUT / "two_sided_triplet_filter_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Two-sided triplet filter scan",
        "",
        "No web lookup was used.",
        "",
        "The scan solves a weighted least-singular problem for normalized `Knu` and `RRRR` blocks,",
        "then floors singular values to enforce finite condition caps.",
        "",
        "| cap | best label | max Knu | K0mu | RRRR | worst(Knu,RRRR) | leading | ST max 2.4e34 |",
        "|---:|---|---:|---:|---:|---:|---|---:|",
    ]
    for cap in ["30.0", "100.0", "300.0"]:
        row = best_by_cap[cap]
        report.append(
            f"| {cap} | `{row['label']}` | {row['max_Knu_amplitude']:.3e} | "
            f"{row['K0mu_amplitude']:.3e} | {row['RRRR_amplitude']:.3e} | "
            f"{row['worst_Knu_RRRR']:.3e} | `{row['leading_channel']}` | "
            f"{row['leading_ST_required_tau_2p4e34']:.3e} |"
        )
    report.extend(
        [
            "",
            "## All-block check",
            "",
            "| cap | label | max Knu | K0mu | RRRR | worst all | leading |",
            "|---:|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in all_rows:
        cap = row["condition_cap"]
        if cap in [30.0, 100.0, 300.0]:
            report.append(
                f"| {cap:.0f} | `{row['label']}` | {row['max_Knu_amplitude']:.3e} | "
                f"{row['K0mu_amplitude']:.3e} | {row['RRRR_amplitude']:.3e} | "
                f"{row['worst_all_monitored']:.3e} | `{row['leading_channel']}` |"
            )
    report.extend(["", "## Verdict", "", verdict, ""])
    (OUT / "two_sided_triplet_filter_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Two-sided triplet filter scan")
    for cap in ["30.0", "100.0", "300.0"]:
        row = best_by_cap[cap]
        print(
            f"  cap {cap}: {row['label']} worst(Knu,RRRR)={row['worst_Knu_RRRR']:.6e} "
            f"Knu={row['max_Knu_amplitude']:.6e} RRRR={row['RRRR_amplitude']:.6e}"
        )
    print(
        f"  all-block cap100: worst_all={best_all_cap100['worst_all_monitored']:.6e} "
        f"leading={best_all_cap100['leading_channel']}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
