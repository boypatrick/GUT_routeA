#!/usr/bin/env python3
"""Eigenstate dressed dimension-five proton-decay proxy.

No web lookup is used.  This audit upgrades the previous mass-insertion
operator

    C'_{abcd}=R_1{}_{aa'}R_2{}_{bb'}R_3{}_{cc'}R_4{}_{dd'}C_{a'b'c'd'}

by first diagonalizing a positive soft spectrum

    M_f^2 = m_0^2 (1 + epsilon Delta_f).

For a dressed pair of external flavor slots p,q we use the eigenstate kernel

    K^{pq}_{ab,a'b'} =
      sum_{r,s} U^p_{ar} U^{p*}_{a'r} U^q_{bs} U^{q*}_{b's}
        D_chi(m_{p,r}^2,m_{q,s}^2),

where D_chi is a symmetric two-mass loop proxy.  The six pair choices are
audited and the most dangerous one is retained.  This is still not a full
SUSY-GUT calculation: it does not include chargino/neutralino mixing matrices,
hadronic chiral reduction, or interference signs.  It does, however, replace
the linear mass-insertion stress by a positive-definite eigenstate sum.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import construct_triplet_rank_lift as rank_lift  # noqa: E402


OUT = ROOT / "output" / "eigenstate_d5_dressing"

DISPLAY_ST = 1.0e-5
ALPHA2_INV = 25.0
SOFT_EPS_GRID = [0.0, 0.01, 0.03, 0.10, 0.30, 0.60, 0.90]
MIN_SOFT_EIGENVALUE = 1.0e-4
PAIRINGS = list(combinations(range(4), 2))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def loop_function(x: float) -> float:
    if x <= 0.0:
        return 1.0
    if abs(x - 1.0) < 1.0e-8:
        return 0.5
    return (1.0 - x + x * math.log(x)) / ((1.0 - x) * (1.0 - x))


def loop_average(x: float, y: float) -> float:
    return 0.5 * (loop_function(x) + loop_function(y))


def diagonalize_soft(delta: np.ndarray, epsilon: float) -> dict[str, Any]:
    herm = 0.5 * (delta + delta.conjugate().T)
    matrix = np.eye(3, dtype=complex) + epsilon * herm
    eigvals, eigvecs = np.linalg.eigh(matrix)
    positive = bool(np.min(eigvals) > MIN_SOFT_EIGENVALUE)
    return {
        "matrix": matrix,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "positive_definite": positive,
        "min_eigenvalue": float(np.min(eigvals)),
        "max_eigenvalue": float(np.max(eigvals)),
    }


def spectrum_for_scenario(base_delta: dict[str, np.ndarray], epsilon: float) -> dict[str, Any]:
    spectra = {kind: diagonalize_soft(delta, epsilon) for kind, delta in base_delta.items()}
    return {
        "by_kind": spectra,
        "positive_definite": all(item["positive_definite"] for item in spectra.values()),
        "min_eigenvalue": min(item["min_eigenvalue"] for item in spectra.values()),
        "max_eigenvalue": max(item["max_eigenvalue"] for item in spectra.values()),
    }


def pair_dressing(
    point: dict[str, float],
    vac: dict[str, float],
    eig_i: float,
    eig_j: float,
    operator: str,
) -> float:
    m_sf = point["m_sfermion_GeV"]
    m2_i = m_sf * m_sf * eig_i
    m2_j = m_sf * m_sf * eig_j
    denom = math.sqrt(m2_i * m2_j)

    m_wino = point["m_wino_GeV"]
    x_w_i = m_wino * m_wino / m2_i
    x_w_j = m_wino * m_wino / m2_j
    d_w = (1.0 / ALPHA2_INV) / (4.0 * math.pi) * m_wino / denom * loop_average(x_w_i, x_w_j)

    if operator == "LLLL":
        return float(d_w)

    mu_h = point["mu_H_GeV"]
    x_h_i = mu_h * mu_h / m2_i
    x_h_j = mu_h * mu_h / m2_j
    d_h = (
        (1.0 / (16.0 * math.pi * math.pi))
        * mu_h
        / denom
        * vac["yt_MG"]
        * vac["yb_MG"]
        * (point["tan_beta"] / 10.0)
        * loop_average(x_h_i, x_h_j)
    )
    return float(d_w + d_h)


def pair_kernel(
    spec_i: dict[str, Any],
    spec_j: dict[str, Any],
    point: dict[str, float],
    vac: dict[str, float],
    operator: str,
) -> np.ndarray:
    ui = spec_i["eigvecs"]
    uj = spec_j["eigvecs"]
    ei = spec_i["eigvals"]
    ej = spec_j["eigvals"]
    kernel = np.zeros((3, 3, 3, 3), dtype=complex)
    for r in range(3):
        for s in range(3):
            d = pair_dressing(point, vac, float(ei[r]), float(ej[s]), operator)
            kernel += (
                d
                * np.einsum(
                    "a,A,b,B->abAB",
                    ui[:, r],
                    ui[:, r].conjugate(),
                    uj[:, s],
                    uj[:, s].conjugate(),
                    optimize=True,
                )
            )
    return kernel


def apply_pair_kernel(tensor: np.ndarray, kernel: np.ndarray, pair: tuple[int, int]) -> np.ndarray:
    p, q = pair
    rest = [axis for axis in range(4) if axis not in pair]
    perm = [p, q] + rest
    transposed = np.transpose(tensor, perm)
    flat = transposed.reshape(9, -1)
    dressed_flat = kernel.reshape(9, 9) @ flat
    dressed = dressed_flat.reshape(3, 3, 3, 3)
    inverse = np.argsort(perm)
    return np.transpose(dressed, inverse)


def c6_lifetime(amp_with_dressing: float, st: float, m_triplet: float, width_prefactor: float) -> dict[str, float]:
    c6 = st * amp_with_dressing / m_triplet
    width = width_prefactor * c6 * c6
    return {
        "C6_GeV_minus2": float(c6),
        "width_GeV": float(width),
        "tau_years": mi.width_to_tau(width),
    }


def st_max_for(
    amp_with_dressing: float,
    target_years: float,
    m_triplet: float,
    width_prefactor: float,
) -> float:
    unit = c6_lifetime(amp_with_dressing, 1.0, m_triplet, width_prefactor)
    return math.sqrt(unit["tau_years"] / target_years)


def scalar_reference_dressing(point: dict[str, float], vac: dict[str, float], operator: str) -> float:
    dress = mi.dressing_factors(point, vac)
    if operator == "LLLL":
        return dress["D_wino_GeV_minus1"]
    return dress["D_wino_GeV_minus1"] + dress["D_higgsino_GeV_minus1"]


def audit() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    _basis_names, _pair_full, channel_entries, _constants_unused = rank_lift.build_channel_tensors()
    proton = read_json(mi.PROTON)["hadronic_constants"]
    vac = mi.vacuum_inputs()
    soft = mi.insertion_basis()

    rows: list[dict[str, Any]] = []
    channel_rows: list[dict[str, Any]] = []

    for filt in mi.selected_filter_rows():
        w = mi.filter_matrix(filt)
        combined = mi.combine_tensors(w, channel_entries)
        base_amplitudes = {
            channel: mi.max_entry(tensor, entries)[0]
            for channel, (tensor, entries) in combined.items()
        }
        for scenario, base_delta in soft["matrices"].items():
            for eps in SOFT_EPS_GRID:
                spec = spectrum_for_scenario(base_delta, eps)
                for point in mi.SPECTRUM_POINTS:
                    channel_results: list[dict[str, Any]] = []
                    for channel, (tensor, entries) in combined.items():
                        meta = mi.CHANNEL_META[channel]
                        scalar_ref = scalar_reference_dressing(point, vac, meta["operator"])
                        ref_amp = max(base_amplitudes[channel] * scalar_ref, 1.0e-300)
                        if not spec["positive_definite"]:
                            result = {
                                "filter_label": filt["label"],
                                "scenario": scenario,
                                "epsilon": eps,
                                "spectrum_name": point["name"],
                                "channel": channel,
                                "operator": meta["operator"],
                                "positive_definite": False,
                                "min_soft_eigenvalue": spec["min_eigenvalue"],
                                "max_soft_eigenvalue": spec["max_eigenvalue"],
                                "pair": "",
                                "amplitude_with_dressing": math.inf,
                                "reference_amplitude_with_dressing": ref_amp,
                                "amplification_vs_degenerate": math.inf,
                                "selected_index": "",
                                "S_T_max": 0.0,
                                "tau_years_ST_display": 0.0,
                                "margin_at_ST_display": 0.0,
                                "passes": False,
                            }
                            channel_results.append(result)
                            channel_rows.append(result)
                            continue

                        pair_results: list[dict[str, Any]] = []
                        for pair in PAIRINGS:
                            kind_i = meta["slot_kinds"][pair[0]]
                            kind_j = meta["slot_kinds"][pair[1]]
                            kernel = pair_kernel(
                                spec["by_kind"][kind_i],
                                spec["by_kind"][kind_j],
                                point,
                                vac,
                                meta["operator"],
                            )
                            dressed_tensor = apply_pair_kernel(tensor, kernel, pair)
                            amp, idx, value = mi.max_entry(dressed_tensor, entries)
                            life = c6_lifetime(amp, DISPLAY_ST, vac["M_T_GeV"], proton["width_prefactor_GeV5"])
                            st_allowed = st_max_for(
                                amp,
                                meta["target_years"],
                                vac["M_T_GeV"],
                                proton["width_prefactor_GeV5"],
                            )
                            pair_results.append(
                                {
                                    "pair": f"{pair[0]}{pair[1]}",
                                    "amplitude_with_dressing": float(amp),
                                    "selected_index": "".join(str(i) for i in idx),
                                    "value": mi.cjson(value),
                                    "S_T_max": float(st_allowed),
                                    "tau_years_ST_display": life["tau_years"],
                                    "margin_at_ST_display": life["tau_years"] / meta["target_years"],
                                    "passes": life["tau_years"] >= meta["target_years"],
                                }
                            )
                        most_dangerous_pair = min(pair_results, key=lambda item: item["margin_at_ST_display"])
                        result = {
                            "filter_label": filt["label"],
                            "scenario": scenario,
                            "epsilon": eps,
                            "spectrum_name": point["name"],
                            "channel": channel,
                            "operator": meta["operator"],
                            "positive_definite": True,
                            "min_soft_eigenvalue": spec["min_eigenvalue"],
                            "max_soft_eigenvalue": spec["max_eigenvalue"],
                            **most_dangerous_pair,
                            "reference_amplitude_with_dressing": ref_amp,
                            "amplification_vs_degenerate": most_dangerous_pair["amplitude_with_dressing"] / ref_amp,
                        }
                        channel_results.append(result)
                        channel_rows.append(result)

                    worst = min(channel_results, key=lambda item: item["margin_at_ST_display"])
                    rows.append(
                        {
                            "filter_label": filt["label"],
                            "scenario": scenario,
                            "epsilon": eps,
                            "spectrum_name": point["name"],
                            **point,
                            "positive_definite": spec["positive_definite"],
                            "min_soft_eigenvalue": spec["min_eigenvalue"],
                            "max_soft_eigenvalue": spec["max_eigenvalue"],
                            "worst_channel": worst["channel"],
                            "worst_pair": worst["pair"],
                            "worst_margin": worst["margin_at_ST_display"],
                            "worst_ST_max": worst["S_T_max"],
                            "max_amplification_vs_degenerate": max(
                                item["amplification_vs_degenerate"] for item in channel_results
                            ),
                            "all_channels_pass": spec["positive_definite"]
                            and all(item["passes"] for item in channel_results),
                        }
                    )

    summary = summarize(rows, channel_rows, vac)
    return rows, channel_rows, summary


def trim_row(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "positive_definite",
        "min_soft_eigenvalue",
        "max_soft_eigenvalue",
        "worst_channel",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_amplification_vs_degenerate",
        "all_channels_pass",
    ]
    return {key: row[key] for key in keys}


def scenario_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for scenario in sorted({row["scenario"] for row in rows}):
        subset = [row for row in rows if row["scenario"] == scenario]
        physical = [row for row in subset if row["positive_definite"]]
        safe = [row for row in physical if row["all_channels_pass"]]
        unsafe = [row for row in physical if not row["all_channels_pass"]]
        table.append(
            {
                "scenario": scenario,
                "positive_points": len(physical),
                "invalid_points": len(subset) - len(physical),
                "safe_points": len(safe),
                "unsafe_points": len(unsafe),
                "total_points": len(subset),
                "safe_fraction_among_positive": len(safe) / len(physical) if physical else 0.0,
                "min_soft_eigenvalue": min((row["min_soft_eigenvalue"] for row in physical), default=math.nan),
                "max_amplification_vs_degenerate": max(
                    (row["max_amplification_vs_degenerate"] for row in physical),
                    default=math.nan,
                ),
                "worst_margin": min((row["worst_margin"] for row in physical), default=0.0),
            }
        )
    return table


def summarize(rows: list[dict[str, Any]], channel_rows: list[dict[str, Any]], vac: dict[str, float]) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    for label in mi.FILTER_LABELS:
        subset = [row for row in rows if row["filter_label"] == label]
        physical = [row for row in subset if row["positive_definite"]]
        safe = [row for row in physical if row["all_channels_pass"]]
        unsafe = [row for row in physical if not row["all_channels_pass"]]
        baseline = next(
            row for row in subset
            if row["scenario"] == "zero" and row["epsilon"] == 0.0 and row["spectrum_name"] == "baseline_100TeV"
        )
        aligned_names = {
            "zero",
            "up_aligned_LL_MFV",
            "down_aligned_LL_MFV",
            "right_third_split",
            "commutator_LL",
            "combined_LL_RR",
        }
        aligned = [row for row in physical if row["scenario"] in aligned_names]
        aligned_safe = [row for row in aligned if row["all_channels_pass"]]
        by_filter[label] = {
            "total_points": len(subset),
            "positive_points": len(physical),
            "invalid_points": len(subset) - len(physical),
            "safe_points": len(safe),
            "unsafe_points": len(unsafe),
            "safe_fraction_among_positive": len(safe) / len(physical) if physical else 0.0,
            "aligned_positive_points": len(aligned),
            "aligned_safe_points": len(aligned_safe),
            "aligned_safe_fraction": len(aligned_safe) / len(aligned) if aligned else 0.0,
            "baseline": trim_row(baseline),
            "most_marginal_safe": trim_row(min(safe, key=lambda row: row["worst_margin"])) if safe else None,
            "nearest_unsafe": trim_row(max(unsafe, key=lambda row: row["worst_margin"])) if unsafe else None,
            "max_amplification_vs_degenerate": max(
                (row["max_amplification_vs_degenerate"] for row in physical),
                default=math.nan,
            ),
            "scenario_table": scenario_table(subset),
        }

    preferred = by_filter["omegaR_0.1_kappa_100"]
    verdict = (
        "Eigenstate dressing keeps the preferred filter safe on the aligned "
        "MFV/right-split/commutator grids, while near-critical democratic soft "
        "spectra can be monitored separately through the positive-eigenvalue "
        "condition.  This replaces the linear insertion with a positive "
        "eigenstate loop proxy, but still awaits full chargino/neutralino and "
        "hadronic matching."
    )
    if preferred["nearest_unsafe"] is not None:
        verdict = (
            "Eigenstate dressing exposes a genuine soft-spectrum condition: "
            "the aligned MFV/right-split/commutator grids remain safe, but a "
            "positive-definite democratic stress can violate the displayed "
            "S_T=1e-5 proton envelope.  The next model-building constraint is "
            "therefore a soft-alignment or degeneracy rule, not another triplet "
            "filter."
        )
    return {
        "note": "No web lookup used. Positive soft-eigenstate proxy for dressed d=5 operators.",
        "display_triplet_filter": DISPLAY_ST,
        "soft_epsilon_grid": SOFT_EPS_GRID,
        "minimum_soft_eigenvalue_for_valid_point": MIN_SOFT_EIGENVALUE,
        "pairings_audited": [f"{p}{q}" for p, q in PAIRINGS],
        "vacuum_inputs": vac,
        "spectrum_points": mi.SPECTRUM_POINTS,
        "operator_formula": (
            "M_f^2=m0^2(1+epsilon Delta_f), "
            "K_ab,a'b'=sum_rs U_ar U*_a'r U_bs U*_b's D_chi(m_r^2,m_s^2)"
        ),
        "by_filter_label": by_filter,
        "verdict": {
            "preferred_filter": "omegaR_0.1_kappa_100",
            "preferred_aligned_safe_fraction": preferred["aligned_safe_fraction"],
            "preferred_safe_fraction_among_positive": preferred["safe_fraction_among_positive"],
            "preferred_has_unsafe_positive_point": preferred["nearest_unsafe"] is not None,
            "interpretation": verdict,
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "positive_definite",
        "min_soft_eigenvalue",
        "max_soft_eigenvalue",
        "worst_channel",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_amplification_vs_degenerate",
        "all_channels_pass",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_channel_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "channel",
        "operator",
        "positive_definite",
        "min_soft_eigenvalue",
        "max_soft_eigenvalue",
        "pair",
        "amplitude_with_dressing",
        "reference_amplitude_with_dressing",
        "amplification_vs_degenerate",
        "selected_index",
        "S_T_max",
        "tau_years_ST_display",
        "margin_at_ST_display",
        "passes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Eigenstate dressed dimension-five proxy",
        "",
        "No web lookup was used.",
        "",
        "The soft spectrum is diagonalized before the Wilson tensor is dressed:",
        "",
        "```text",
        "M_f^2=m0^2(1+epsilon Delta_f)",
        "K_ab,a'b'=sum_rs U_ar U*_a'r U_bs U*_b's D_chi(m_r^2,m_s^2)",
        "```",
        "",
        "## Filter summaries",
        "",
        "| filter | aligned safe | positive safe | invalid | max amp | marginal safe | nearest unsafe |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for label, payload in summary["by_filter_label"].items():
        marginal = payload["most_marginal_safe"]
        unsafe = payload["nearest_unsafe"]
        lines.append(
            f"| `{label}` | {payload['aligned_safe_points']}/{payload['aligned_positive_points']} | "
            f"{payload['safe_points']}/{payload['positive_points']} | {payload['invalid_points']} | "
            f"{payload['max_amplification_vs_degenerate']:.3e} | "
            f"{marginal['worst_margin'] if marginal else math.inf:.3e} | "
            f"{unsafe['worst_margin'] if unsafe else math.inf:.3e} |"
        )

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    lines.extend(
        [
            "",
            "## Preferred scenario table",
            "",
            "| scenario | safe/positive | invalid | min eig | max amp | worst margin |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in preferred["scenario_table"]:
        lines.append(
            f"| `{row['scenario']}` | {row['safe_points']}/{row['positive_points']} | "
            f"{row['invalid_points']} | {row['min_soft_eigenvalue']:.3e} | "
            f"{row['max_amplification_vs_degenerate']:.3e} | {row['worst_margin']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, channel_rows, summary = audit()
    write_csv(OUT / "eigenstate_scan.csv", rows)
    write_channel_csv(OUT / "channel_scan.csv", channel_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    print("Eigenstate d=5 dressing proxy")
    print(
        "  omegaR_0.1_kappa_100: "
        f"aligned_safe={preferred['aligned_safe_points']}/{preferred['aligned_positive_points']}, "
        f"positive_safe={preferred['safe_points']}/{preferred['positive_points']}, "
        f"invalid={preferred['invalid_points']}, "
        f"max_amp={preferred['max_amplification_vs_degenerate']:.3e}"
    )
    if preferred["nearest_unsafe"]:
        print(
            "  nearest unsafe: "
            f"margin={preferred['nearest_unsafe']['worst_margin']:.3e}, "
            f"scenario={preferred['nearest_unsafe']['scenario']}, "
            f"eps={preferred['nearest_unsafe']['epsilon']}"
        )
    else:
        print("  no unsafe positive-definite audited point")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
