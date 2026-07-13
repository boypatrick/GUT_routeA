#!/usr/bin/env python3
"""MSSM-mixing dressed dimension-five proton-decay proxy.

No web lookup is used.  This script upgrades the positive soft-eigenstate
proxy by replacing the scalar wino/higgsino dressing factors with explicit
chargino and neutralino mass matrices.  It is still a controlled proxy, not a
publication-level SUSY proton-decay code: it uses coherent absolute sums,
minimal channel weights, and no lattice/chiral matrix elements beyond the
local hadronic prefactor already used by the project.

The new information is whether the previously safe triplet filter survives
when the electroweakino eigenstate content is resolved:

    X_C = [[M2, sqrt(2) mW sin beta],
           [sqrt(2) mW cos beta, mu]]

and the real neutralino matrix in the (B,W0,Hd0,Hu0) basis.  The dressed
pair kernel is

    K_ab,a'b' = sum_rs U_ar U*_a'r U_bs U*_b's
                sum_chi D_chi(m_r^2,m_s^2),

where D_chi includes chargino and neutralino component weights.
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

import audit_eigenstate_d5_dressing as eig  # noqa: E402
import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import construct_triplet_rank_lift as rank_lift  # noqa: E402


OUT = ROOT / "output" / "mssm_mixing_d5_dressing"

DISPLAY_ST = 1.0e-5
MW_GEV = 80.379
MZ_GEV = 91.1876
SIN2_THETA_W = 0.23122
COS2_THETA_W = 1.0 - SIN2_THETA_W
ALPHA2_INV = 25.0
ALPHA2 = 1.0 / ALPHA2_INV
ALPHA_Y = ALPHA2 * SIN2_THETA_W / COS2_THETA_W
MIN_SOFT_EIGENVALUE = eig.MIN_SOFT_EIGENVALUE

HYPERCHARGE = {
    "uL": 1.0 / 6.0,
    "dL": 1.0 / 6.0,
    "nuL": -1.0 / 2.0,
    "eL": -1.0 / 2.0,
    "uR": -2.0 / 3.0,
    "dR": 1.0 / 3.0,
    "eR": 1.0,
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def beta_angles(tan_beta: float) -> tuple[float, float]:
    cosb = 1.0 / math.sqrt(1.0 + tan_beta * tan_beta)
    sinb = tan_beta * cosb
    return sinb, cosb


def chargino_spectrum(m2: float, mu: float, tan_beta: float) -> list[dict[str, float]]:
    sinb, cosb = beta_angles(tan_beta)
    x = np.array(
        [
            [m2, math.sqrt(2.0) * MW_GEV * sinb],
            [math.sqrt(2.0) * MW_GEV * cosb, mu],
        ],
        dtype=float,
    )
    u, masses, vh = np.linalg.svd(x)
    v = vh.conjugate().T
    out: list[dict[str, float]] = []
    for n, mass in enumerate(masses):
        wino = 0.5 * (abs(u[0, n]) ** 2 + abs(v[0, n]) ** 2)
        higgsino = 1.0 - wino
        out.append(
            {
                "mass_GeV": float(abs(mass)),
                "wino_fraction": float(wino),
                "higgsino_fraction": float(higgsino),
            }
        )
    return out


def neutralino_spectrum(m1: float, m2: float, mu: float, tan_beta: float) -> list[dict[str, float]]:
    sinb, cosb = beta_angles(tan_beta)
    sw = math.sqrt(SIN2_THETA_W)
    cw = math.sqrt(COS2_THETA_W)
    mat = np.array(
        [
            [m1, 0.0, -MZ_GEV * sw * cosb, MZ_GEV * sw * sinb],
            [0.0, m2, MZ_GEV * cw * cosb, -MZ_GEV * cw * sinb],
            [-MZ_GEV * sw * cosb, MZ_GEV * cw * cosb, 0.0, -mu],
            [MZ_GEV * sw * sinb, -MZ_GEV * cw * sinb, -mu, 0.0],
        ],
        dtype=float,
    )
    eigvals, eigvecs = np.linalg.eigh(mat)
    order = np.argsort(np.abs(eigvals))
    out: list[dict[str, float]] = []
    for idx in order:
        vec = eigvecs[:, idx]
        out.append(
            {
                "mass_GeV": float(abs(eigvals[idx])),
                "bino_fraction": float(abs(vec[0]) ** 2),
                "wino_fraction": float(abs(vec[1]) ** 2),
                "higgsino_fraction": float(abs(vec[2]) ** 2 + abs(vec[3]) ** 2),
                "signed_mass_GeV": float(eigvals[idx]),
            }
        )
    return out


def electroweakino_spectrum(point: dict[str, float]) -> dict[str, Any]:
    m2 = point["m_wino_GeV"]
    # GUT-inspired low-scale proxy M1 ~= (5/3) tan^2 theta_W M2.
    m1 = (5.0 / 3.0) * SIN2_THETA_W / COS2_THETA_W * m2
    mu = point["mu_H_GeV"]
    tan_beta = point["tan_beta"]
    charginos = chargino_spectrum(m2, mu, tan_beta)
    neutralinos = neutralino_spectrum(m1, m2, mu, tan_beta)
    return {
        "M1_GeV": float(m1),
        "M2_GeV": float(m2),
        "mu_H_GeV": float(mu),
        "tan_beta": float(tan_beta),
        "charginos": charginos,
        "neutralinos": neutralinos,
        "lightest_chargino_GeV": min(row["mass_GeV"] for row in charginos),
        "lightest_neutralino_GeV": min(row["mass_GeV"] for row in neutralinos),
    }


def gaugino_pair_weight(kind_i: str, kind_j: str, neutralino: dict[str, float]) -> float:
    y_weight = abs(HYPERCHARGE[kind_i] * HYPERCHARGE[kind_j])
    return ALPHA2 * neutralino["wino_fraction"] + ALPHA_Y * y_weight * neutralino["bino_fraction"]


def charged_higgsino_weight(kind_i: str, kind_j: str) -> float:
    kinds = {kind_i, kind_j}
    if kinds in [{"uR", "dR"}, {"uR", "eR"}]:
        return 1.0
    if kind_i.endswith("R") and kind_j.endswith("R"):
        return 0.5
    return 0.0


def pair_dressing(
    point: dict[str, float],
    vac: dict[str, float],
    ewkino: dict[str, Any],
    eig_i: float,
    eig_j: float,
    kind_i: str,
    kind_j: str,
    operator: str,
) -> dict[str, float]:
    m_sf = point["m_sfermion_GeV"]
    m2_i = m_sf * m_sf * eig_i
    m2_j = m_sf * m_sf * eig_j
    denom = math.sqrt(m2_i * m2_j)

    d_chargino = 0.0
    d_neutralino = 0.0
    d_higgsino = 0.0

    for chi in ewkino["charginos"]:
        mass = chi["mass_GeV"]
        loop = eig.loop_average(mass * mass / m2_i, mass * mass / m2_j)
        if operator == "LLLL":
            d_chargino += (
                ALPHA2
                / (4.0 * math.pi)
                * chi["wino_fraction"]
                * mass
                / denom
                * loop
            )
        else:
            d_higgsino += (
                1.0
                / (16.0 * math.pi * math.pi)
                * chi["higgsino_fraction"]
                * charged_higgsino_weight(kind_i, kind_j)
                * vac["yt_MG"]
                * vac["yb_MG"]
                * (point["tan_beta"] / 10.0)
                * mass
                / denom
                * loop
            )

    for neu in ewkino["neutralinos"]:
        mass = neu["mass_GeV"]
        loop = eig.loop_average(mass * mass / m2_i, mass * mass / m2_j)
        if operator == "LLLL":
            d_neutralino += (
                gaugino_pair_weight(kind_i, kind_j, neu)
                / (4.0 * math.pi)
                * mass
                / denom
                * loop
            )
        else:
            y_weight = abs(HYPERCHARGE[kind_i] * HYPERCHARGE[kind_j])
            d_neutralino += (
                ALPHA_Y
                * y_weight
                * neu["bino_fraction"]
                / (4.0 * math.pi)
                * mass
                / denom
                * loop
            )
            d_higgsino += (
                1.0
                / (16.0 * math.pi * math.pi)
                * neu["higgsino_fraction"]
                * vac["yt_MG"]
                * vac["yb_MG"]
                * (point["tan_beta"] / 10.0)
                * mass
                / denom
                * loop
            )

    total = d_chargino + d_neutralino + d_higgsino
    return {
        "total": float(total),
        "chargino": float(d_chargino),
        "neutralino": float(d_neutralino),
        "higgsino": float(d_higgsino),
    }


def pair_kernel(
    spec_i: dict[str, Any],
    spec_j: dict[str, Any],
    point: dict[str, float],
    vac: dict[str, float],
    ewkino: dict[str, Any],
    kind_i: str,
    kind_j: str,
    operator: str,
) -> tuple[np.ndarray, dict[str, float]]:
    ui = spec_i["eigvecs"]
    uj = spec_j["eigvecs"]
    ei = spec_i["eigvals"]
    ej = spec_j["eigvals"]
    kernel = np.zeros((3, 3, 3, 3), dtype=complex)
    parts = {"total": 0.0, "chargino": 0.0, "neutralino": 0.0, "higgsino": 0.0}
    for r in range(3):
        for s in range(3):
            dress = pair_dressing(point, vac, ewkino, float(ei[r]), float(ej[s]), kind_i, kind_j, operator)
            for key in parts:
                parts[key] += dress[key] / 9.0
            kernel += (
                dress["total"]
                * np.einsum(
                    "a,A,b,B->abAB",
                    ui[:, r],
                    ui[:, r].conjugate(),
                    uj[:, s],
                    uj[:, s].conjugate(),
                    optimize=True,
                )
            )
    return kernel, parts


def audit() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    _basis_names, _pair_full, channel_entries, _constants_unused = rank_lift.build_channel_tensors()
    proton = read_json(mi.PROTON)["hadronic_constants"]
    vac = mi.vacuum_inputs()
    soft = mi.insertion_basis()
    eigen_summary = read_json(eig.OUT / "summary.json")

    rows: list[dict[str, Any]] = []
    channel_rows: list[dict[str, Any]] = []
    spectrum_cards: dict[str, Any] = {}

    for filt in mi.selected_filter_rows():
        w = mi.filter_matrix(filt)
        combined = mi.combine_tensors(w, channel_entries)
        for scenario, base_delta in soft["matrices"].items():
            for eps in eig.SOFT_EPS_GRID:
                spec = eig.spectrum_for_scenario(base_delta, eps)
                for point in mi.SPECTRUM_POINTS:
                    ewkino = electroweakino_spectrum(point)
                    spectrum_cards[point["name"]] = ewkino
                    channel_results: list[dict[str, Any]] = []
                    for channel, (tensor, entries) in combined.items():
                        meta = mi.CHANNEL_META[channel]
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
                                "amplification_vs_eigenstate_proxy": math.inf,
                                "selected_index": "",
                                "S_T_max": 0.0,
                                "tau_years_ST_display": 0.0,
                                "margin_at_ST_display": 0.0,
                                "passes": False,
                                "avg_chargino_part": 0.0,
                                "avg_neutralino_part": 0.0,
                                "avg_higgsino_part": 0.0,
                            }
                            channel_results.append(result)
                            channel_rows.append(result)
                            continue

                        pair_results: list[dict[str, Any]] = []
                        for pair in eig.PAIRINGS:
                            kind_i = meta["slot_kinds"][pair[0]]
                            kind_j = meta["slot_kinds"][pair[1]]
                            kernel, parts = pair_kernel(
                                spec["by_kind"][kind_i],
                                spec["by_kind"][kind_j],
                                point,
                                vac,
                                ewkino,
                                kind_i,
                                kind_j,
                                meta["operator"],
                            )
                            dressed_tensor = eig.apply_pair_kernel(tensor, kernel, pair)
                            amp, idx, value = mi.max_entry(dressed_tensor, entries)
                            life = eig.c6_lifetime(amp, DISPLAY_ST, vac["M_T_GeV"], proton["width_prefactor_GeV5"])
                            st_allowed = eig.st_max_for(
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
                                    "avg_chargino_part": parts["chargino"],
                                    "avg_neutralino_part": parts["neutralino"],
                                    "avg_higgsino_part": parts["higgsino"],
                                    "avg_total_part": parts["total"],
                                }
                            )
                        most_dangerous = min(pair_results, key=lambda item: item["margin_at_ST_display"])
                        # The previous eigenstate proxy is the closest same-grid
                        # reference.  It can differ by channel/pair, so we keep
                        # the comparison conservative and local to the same
                        # result row where possible.
                        reference = max(
                            (
                                row["amplitude_with_dressing"]
                                for row in read_reference_channel_rows()
                                if row["filter_label"] == filt["label"]
                                and row["scenario"] == scenario
                                and float(row["epsilon"]) == eps
                                and row["spectrum_name"] == point["name"]
                                and row["channel"] == channel
                            ),
                            default=1.0e-300,
                        )
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
                            **most_dangerous,
                            "amplification_vs_eigenstate_proxy": most_dangerous["amplitude_with_dressing"] / reference,
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
                            "M1_GeV": ewkino["M1_GeV"],
                            "lightest_chargino_GeV": ewkino["lightest_chargino_GeV"],
                            "lightest_neutralino_GeV": ewkino["lightest_neutralino_GeV"],
                            "positive_definite": spec["positive_definite"],
                            "min_soft_eigenvalue": spec["min_eigenvalue"],
                            "max_soft_eigenvalue": spec["max_eigenvalue"],
                            "worst_channel": worst["channel"],
                            "worst_pair": worst["pair"],
                            "worst_margin": worst["margin_at_ST_display"],
                            "worst_ST_max": worst["S_T_max"],
                            "max_amplification_vs_eigenstate_proxy": max(
                                item["amplification_vs_eigenstate_proxy"] for item in channel_results
                            ),
                            "all_channels_pass": spec["positive_definite"]
                            and all(item["passes"] for item in channel_results),
                        }
                    )

    summary = summarize(rows, channel_rows, vac, spectrum_cards, eigen_summary)
    return rows, channel_rows, summary


_REFERENCE_ROWS: list[dict[str, Any]] | None = None


def read_reference_channel_rows() -> list[dict[str, Any]]:
    global _REFERENCE_ROWS
    if _REFERENCE_ROWS is not None:
        return _REFERENCE_ROWS
    path = eig.OUT / "channel_scan.csv"
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "amplitude_with_dressing",
                "S_T_max",
                "tau_years_ST_display",
                "margin_at_ST_display",
                "amplification_vs_degenerate",
            ]:
                row[key] = float(row[key])
            rows.append(row)
    _REFERENCE_ROWS = rows
    return rows


def trim_row(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "M1_GeV",
        "mu_H_GeV",
        "tan_beta",
        "lightest_chargino_GeV",
        "lightest_neutralino_GeV",
        "positive_definite",
        "min_soft_eigenvalue",
        "worst_channel",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_amplification_vs_eigenstate_proxy",
        "all_channels_pass",
    ]
    return {key: row[key] for key in keys}


def scenario_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for scenario in sorted({row["scenario"] for row in rows}):
        subset = [row for row in rows if row["scenario"] == scenario]
        physical = [row for row in subset if row["positive_definite"]]
        safe = [row for row in physical if row["all_channels_pass"]]
        table.append(
            {
                "scenario": scenario,
                "positive_points": len(physical),
                "safe_points": len(safe),
                "total_points": len(subset),
                "safe_fraction_among_positive": len(safe) / len(physical) if physical else 0.0,
                "max_amplification_vs_eigenstate_proxy": max(
                    (row["max_amplification_vs_eigenstate_proxy"] for row in physical),
                    default=math.nan,
                ),
                "worst_margin": min((row["worst_margin"] for row in physical), default=0.0),
            }
        )
    return table


def summarize(
    rows: list[dict[str, Any]],
    channel_rows: list[dict[str, Any]],
    vac: dict[str, float],
    spectrum_cards: dict[str, Any],
    eigen_summary: dict[str, Any],
) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    aligned_names = {
        "zero",
        "up_aligned_LL_MFV",
        "down_aligned_LL_MFV",
        "right_third_split",
        "commutator_LL",
        "combined_LL_RR",
    }
    for label in mi.FILTER_LABELS:
        subset = [row for row in rows if row["filter_label"] == label]
        physical = [row for row in subset if row["positive_definite"]]
        safe = [row for row in physical if row["all_channels_pass"]]
        unsafe = [row for row in physical if not row["all_channels_pass"]]
        aligned = [row for row in physical if row["scenario"] in aligned_names]
        aligned_safe = [row for row in aligned if row["all_channels_pass"]]
        baseline = next(
            row for row in subset
            if row["scenario"] == "zero" and row["epsilon"] == 0.0 and row["spectrum_name"] == "baseline_100TeV"
        )
        by_filter[label] = {
            "total_points": len(subset),
            "positive_points": len(physical),
            "safe_points": len(safe),
            "unsafe_points": len(unsafe),
            "safe_fraction_among_positive": len(safe) / len(physical) if physical else 0.0,
            "aligned_positive_points": len(aligned),
            "aligned_safe_points": len(aligned_safe),
            "aligned_safe_fraction": len(aligned_safe) / len(aligned) if aligned else 0.0,
            "baseline": trim_row(baseline),
            "most_marginal_safe": trim_row(min(safe, key=lambda row: row["worst_margin"])) if safe else None,
            "nearest_unsafe": trim_row(max(unsafe, key=lambda row: row["worst_margin"])) if unsafe else None,
            "max_amplification_vs_eigenstate_proxy": max(
                (row["max_amplification_vs_eigenstate_proxy"] for row in physical),
                default=math.nan,
            ),
            "scenario_table": scenario_table(subset),
        }

    preferred = by_filter["omegaR_0.1_kappa_100"]
    verdict = (
        "Resolving chargino and neutralino mixing keeps the preferred filter "
        "safe over the audited positive soft grid.  The calculation is a "
        "coherent-absolute proxy, so it is deliberately conservative relative "
        "to interference-aware amplitudes, but it still lacks the final chiral "
        "and lattice reduction."
    )
    if preferred["nearest_unsafe"] is not None:
        verdict = (
            "Resolving chargino and neutralino mixing produces an unsafe "
            "positive soft point in the preferred filter.  The next model step "
            "must impose an electroweakino/soft-alignment condition or reduce "
            "the triplet filter target before claiming proton safety."
        )
    return {
        "note": "No web lookup used. MSSM chargino/neutralino mixing proxy for dressed d=5 operators.",
        "display_triplet_filter": DISPLAY_ST,
        "alpha2": ALPHA2,
        "alphaY_proxy": ALPHA_Y,
        "neutralino_M1_relation": "M1=(5/3) tan^2(theta_W) M2",
        "operator_formula": "chargino SVD + real-neutralino eigensystem, coherent absolute dressing sums",
        "vacuum_inputs": vac,
        "electroweakino_spectrum_cards": spectrum_cards,
        "reference_eigenstate_summary": {
            "preferred_safe": eigen_summary["by_filter_label"]["omegaR_0.1_kappa_100"]["safe_points"],
            "preferred_positive": eigen_summary["by_filter_label"]["omegaR_0.1_kappa_100"]["positive_points"],
        },
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
        "M1_GeV",
        "mu_H_GeV",
        "tan_beta",
        "lightest_chargino_GeV",
        "lightest_neutralino_GeV",
        "positive_definite",
        "min_soft_eigenvalue",
        "max_soft_eigenvalue",
        "worst_channel",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_amplification_vs_eigenstate_proxy",
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
        "pair",
        "amplitude_with_dressing",
        "amplification_vs_eigenstate_proxy",
        "selected_index",
        "S_T_max",
        "tau_years_ST_display",
        "margin_at_ST_display",
        "passes",
        "avg_chargino_part",
        "avg_neutralino_part",
        "avg_higgsino_part",
        "avg_total_part",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# MSSM-mixing dressed dimension-five proxy",
        "",
        "No web lookup was used.",
        "",
        "The electroweakino sector is diagonalized explicitly:",
        "",
        "```text",
        "X_C = [[M2, sqrt(2) mW sin beta], [sqrt(2) mW cos beta, mu]]",
        "M_N = real symmetric neutralino matrix in (B,W0,Hd0,Hu0)",
        "```",
        "",
        "## Filter summaries",
        "",
        "| filter | aligned safe | positive safe | max amp vs eigen proxy | marginal safe | nearest unsafe |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, payload in summary["by_filter_label"].items():
        marginal = payload["most_marginal_safe"]
        unsafe = payload["nearest_unsafe"]
        lines.append(
            f"| `{label}` | {payload['aligned_safe_points']}/{payload['aligned_positive_points']} | "
            f"{payload['safe_points']}/{payload['positive_points']} | "
            f"{payload['max_amplification_vs_eigenstate_proxy']:.3e} | "
            f"{marginal['worst_margin'] if marginal else math.inf:.3e} | "
            f"{unsafe['worst_margin'] if unsafe else math.inf:.3e} |"
        )

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    lines.extend(
        [
            "",
            "## Preferred scenario table",
            "",
            "| scenario | safe/positive | max amp vs eigen proxy | worst margin |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in preferred["scenario_table"]:
        lines.append(
            f"| `{row['scenario']}` | {row['safe_points']}/{row['positive_points']} | "
            f"{row['max_amplification_vs_eigenstate_proxy']:.3e} | {row['worst_margin']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, channel_rows, summary = audit()
    write_csv(OUT / "mssm_mixing_scan.csv", rows)
    write_channel_csv(OUT / "channel_scan.csv", channel_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    print("MSSM-mixing d=5 dressing proxy")
    print(
        "  omegaR_0.1_kappa_100: "
        f"aligned_safe={preferred['aligned_safe_points']}/{preferred['aligned_positive_points']}, "
        f"positive_safe={preferred['safe_points']}/{preferred['positive_points']}, "
        f"max_amp={preferred['max_amplification_vs_eigenstate_proxy']:.3e}"
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
