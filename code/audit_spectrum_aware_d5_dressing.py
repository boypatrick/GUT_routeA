#!/usr/bin/env python3
"""Spectrum-aware dressed dimension-five proton-decay proxy.

No web lookup is used.  This is a controlled upgrade of
audit_dressed_dimension5_channels.py.  Instead of using a single external
"relative dressing" multiplier, it scans a small SUSY-spectrum grid and
computes

    D_W = (alpha_2/4 pi) (m_wino/m_sf^2) F(m_wino^2/m_sf^2),
    D_H = (1/16 pi^2) (mu_H/m_sf^2) y_t y_b (tan beta/10) F(mu_H^2/m_sf^2),

with coherent RRRR dressing D_R = D_W + xi_H D_H and LLLL dressing
D_L = D_W.  The xi_H grid is a stress proxy for unresolved flavor rotations
and mass-insertion coherence.
The loop function is F(x)=(1-x+x log x)/(1-x)^2 with F(1)=1/2.

This is still a proxy: it does not replace a full mass-insertion
wino/higgsino dressing calculation.  Its job is to quantify the safe region
for the currently verified two-sided triplet filter.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
TWO_SIDED = ROOT / "output" / "two_sided_triplet_filter" / "two_sided_triplet_filter_summary.json"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "spectrum_aware_d5_dressing"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
DISPLAY_ST = 1.0e-5

CHANNELS = {
    "LLLL_upupdown_Knu": {
        "physical_proxy": "p_to_K+nu_bar_L_upupdown",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
        "dressing_class": "L",
    },
    "LLLL_downdownup_Knu": {
        "physical_proxy": "p_to_K+nu_bar_L_downdownup",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
        "dressing_class": "L",
    },
    "LLLL_upupdown_K0mu": {
        "physical_proxy": "p_to_K0mu+_L_proxy",
        "target_years": 1.0e34,
        "hadronic_relative": 1.0,
        "dressing_class": "L",
    },
    "RRRR_uusd_anycharged": {
        "physical_proxy": "p_to_K+nu_or_K0mu_R_proxy",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
        "dressing_class": "R",
    },
}

FILTER_LABELS = ["omegaR_0.1_kappa_100", "all_blocks_kappa_100"]
M_SFERMION_GRID = [1.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e5, 5.0e5]
M_WINO_GRID = [5.0e2, 1.0e3, 3.0e3]
MU_H_GRID = [5.0e2, 1.0e3, 3.0e3, 1.0e4]
TAN_BETA_GRID = [5.0, 10.0, 30.0, 50.0]
XI_H_GRID = [1.0, 3.0, 10.0, 30.0, 100.0]
ALPHA2_INV = 25.0


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def loop_function(x: float) -> float:
    if x <= 0.0:
        return 1.0
    if abs(x - 1.0) < 1.0e-8:
        return 0.5
    return (1.0 - x + x * math.log(x)) / ((1.0 - x) * (1.0 - x))


def width_to_tau(width_gev: float) -> float:
    if width_gev <= 0.0:
        return math.inf
    return HBAR_GEV_S / width_gev / SECONDS_PER_YEAR


def selected_filter_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = {row["label"]: row for row in payload["rows"]}
    rows.update({row["label"]: row for row in payload["all_block_rows"]})
    return [rows[label] for label in FILTER_LABELS]


def vacuum_params() -> dict[str, float]:
    payload = read_json(VACUUM)
    vac = payload["recommended_benchmark"]
    # The vacuum card records masses in the recommended benchmark, while the
    # Yukawa inputs are kept in the seed/rge replay cache.  Keep this reader
    # tolerant so regenerated cards do not silently fork the physics input.
    yuk = payload.get("seed", {}).get("best", {})
    return {
        "M_T_GeV": float(vac["M_HC_GeV"]),
        "yt_MG": float(vac.get("yt_MG", yuk["yt_MG"])),
        "yb_MG": float(vac.get("yb_MG", yuk["yb_MG"])),
        "ytau_MG": float(vac.get("ytau_MG", yuk["ytau_MG"])),
        "tan_beta_benchmark": float(vac["tan_beta"]),
        "MSUSY_GeV": float(vac["MSUSY_GeV"]),
    }


def dressing_factors(
    m_sfermion: float,
    m_wino: float,
    mu_h: float,
    tan_beta: float,
    yt: float,
    yb: float,
) -> dict[str, float]:
    x_w = (m_wino / m_sfermion) ** 2
    x_h = (mu_h / m_sfermion) ** 2
    d_w = (1.0 / ALPHA2_INV) / (4.0 * math.pi) * m_wino / (m_sfermion * m_sfermion) * loop_function(x_w)
    d_h = (
        (1.0 / (16.0 * math.pi * math.pi))
        * mu_h
        / (m_sfermion * m_sfermion)
        * yt
        * yb
        * (tan_beta / 10.0)
        * loop_function(x_h)
    )
    return {
        "D_wino_GeV_minus1": d_w,
        "D_higgsino_proxy_GeV_minus1": d_h,
        "D_R_coherent_GeV_minus1": d_w + d_h,
        "x_wino": x_w,
        "x_higgsino": x_h,
        "F_wino": loop_function(x_w),
        "F_higgsino": loop_function(x_h),
        "higgsino_to_wino": d_h / d_w if d_w > 0.0 else math.inf,
    }


def lifetime(
    amplitude: float,
    st: float,
    m_triplet: float,
    width_prefactor: float,
    dressing: float,
    hadronic_relative: float,
) -> dict[str, float]:
    c5 = st * amplitude / m_triplet
    c6 = c5 * dressing
    width = width_prefactor * hadronic_relative * c6 * c6
    return {
        "C5_GeV_minus1": c5,
        "C6_GeV_minus2": c6,
        "width_GeV": width,
        "tau_years": width_to_tau(width),
    }


def st_max_for(
    amplitude: float,
    target_years: float,
    m_triplet: float,
    width_prefactor: float,
    dressing: float,
    hadronic_relative: float,
) -> float:
    unit = lifetime(amplitude, 1.0, m_triplet, width_prefactor, dressing, hadronic_relative)
    return math.sqrt(unit["tau_years"] / target_years)


def audit_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    two_sided = read_json(TWO_SIDED)
    proton = read_json(PROTON)["hadronic_constants"]
    vac = vacuum_params()
    rows: list[dict[str, Any]] = []

    for filt in selected_filter_rows(two_sided):
        for m_sf in M_SFERMION_GRID:
            for m_wino in M_WINO_GRID:
                for mu_h in MU_H_GRID:
                    for tan_beta in TAN_BETA_GRID:
                        for xi_h in XI_H_GRID:
                            dress = dressing_factors(
                                m_sf,
                                m_wino,
                                mu_h,
                                tan_beta,
                                vac["yt_MG"],
                                vac["yb_MG"],
                            )
                            d_r_stressed = dress["D_wino_GeV_minus1"] + xi_h * dress["D_higgsino_proxy_GeV_minus1"]
                            channel_results: list[dict[str, Any]] = []
                            for key, channel in CHANNELS.items():
                                amplitude = float(filt["channels"][key]["amplitude"])
                                d_eff = (
                                    dress["D_wino_GeV_minus1"]
                                    if channel["dressing_class"] == "L"
                                    else d_r_stressed
                                )
                                life = lifetime(
                                    amplitude,
                                    DISPLAY_ST,
                                    vac["M_T_GeV"],
                                    proton["width_prefactor_GeV5"],
                                    d_eff,
                                    channel["hadronic_relative"],
                                )
                                st_allowed = st_max_for(
                                    amplitude,
                                    channel["target_years"],
                                    vac["M_T_GeV"],
                                    proton["width_prefactor_GeV5"],
                                    d_eff,
                                    channel["hadronic_relative"],
                                )
                                channel_results.append(
                                    {
                                        "channel_key": key,
                                        "physical_proxy": channel["physical_proxy"],
                                        "amplitude": amplitude,
                                        "dressing_class": channel["dressing_class"],
                                        "effective_dressing_GeV_minus1": d_eff,
                                        "target_years": channel["target_years"],
                                        "tau_years_ST_display": life["tau_years"],
                                        "margin_at_ST_display": life["tau_years"] / channel["target_years"],
                                        "S_T_max": st_allowed,
                                        "passes": life["tau_years"] >= channel["target_years"],
                                    }
                                )
                            worst = min(channel_results, key=lambda row: row["margin_at_ST_display"])
                            rows.append(
                                {
                                    "filter_label": filt["label"],
                                    "condition_cap": filt["condition_cap"],
                                    "m_sfermion_GeV": m_sf,
                                    "m_wino_GeV": m_wino,
                                    "mu_H_GeV": mu_h,
                                    "tan_beta": tan_beta,
                                    "xi_H": xi_h,
                                    **dress,
                                    "D_R_stressed_GeV_minus1": d_r_stressed,
                                    "higgsino_to_wino_stressed": (
                                        xi_h * dress["D_higgsino_proxy_GeV_minus1"] / dress["D_wino_GeV_minus1"]
                                        if dress["D_wino_GeV_minus1"] > 0.0
                                        else math.inf
                                    ),
                                    "worst_channel": worst["physical_proxy"],
                                    "worst_channel_key": worst["channel_key"],
                                    "worst_margin": worst["margin_at_ST_display"],
                                    "worst_tau_years_ST_display": worst["tau_years_ST_display"],
                                    "global_S_T_max": min(row["S_T_max"] for row in channel_results),
                                    "all_channels_pass": all(row["passes"] for row in channel_results),
                                    "channel_results": channel_results,
                                }
                            )
    summary = summarize(rows, vac)
    return rows, summary


def summarize(rows: list[dict[str, Any]], vac: dict[str, float]) -> dict[str, Any]:
    by_label: dict[str, dict[str, Any]] = {}
    for label in FILTER_LABELS:
        subset = [row for row in rows if row["filter_label"] == label]
        safe = [row for row in subset if row["all_channels_pass"]]
        unsafe = [row for row in subset if not row["all_channels_pass"]]
        baseline = next(
            row
            for row in subset
            if row["m_sfermion_GeV"] == 1.0e5
            and row["m_wino_GeV"] == 1.0e3
            and row["mu_H_GeV"] == 1.0e3
            and row["tan_beta"] == 10.0
            and row["xi_H"] == 1.0
        )
        unstressed = [row for row in subset if row["xi_H"] == 1.0]
        worst_safe = min(safe, key=lambda row: row["worst_margin"]) if safe else None
        best_unsafe = max(unsafe, key=lambda row: row["worst_margin"]) if unsafe else None
        by_label[label] = {
            "total_points": len(subset),
            "safe_points": len(safe),
            "unsafe_points": len(unsafe),
            "safe_fraction": len(safe) / len(subset),
            "unstressed_safe_points": sum(1 for row in unstressed if row["all_channels_pass"]),
            "unstressed_total_points": len(unstressed),
            "baseline": trim_row(baseline),
            "most_marginal_safe": trim_row(worst_safe) if worst_safe else None,
            "nearest_unsafe": trim_row(best_unsafe) if best_unsafe else None,
            "minimum_safe_m_sfermion_by_tan_beta": min_safe_mass_table(safe),
            "safe_by_xi_H": safe_by_xi_table(subset),
        }

    preferred = by_label["omegaR_0.1_kappa_100"]
    return {
        "note": "No web lookup used. Spectrum-aware wino/higgsino d=5 dressing proxy.",
        "formula": {
            "D_wino": "(alpha2/4pi) m_wino/m_sf^2 F(m_wino^2/m_sf^2)",
            "D_higgsino_proxy": "(1/16pi^2) mu_H/m_sf^2 y_t y_b (tan beta/10) F(mu_H^2/m_sf^2)",
            "D_L": "D_wino",
            "D_R": "D_wino + xi_H D_higgsino_proxy, coherent stress sign",
            "F": "(1-x+x log x)/(1-x)^2, F(1)=1/2",
        },
        "display_triplet_filter": DISPLAY_ST,
        "vacuum_inputs": vac,
        "grid": {
            "m_sfermion_GeV": M_SFERMION_GRID,
            "m_wino_GeV": M_WINO_GRID,
            "mu_H_GeV": MU_H_GRID,
            "tan_beta": TAN_BETA_GRID,
            "xi_H": XI_H_GRID,
            "alpha2_inv": ALPHA2_INV,
        },
        "by_filter_label": by_label,
        "verdict": {
            "preferred_filter": "omegaR_0.1_kappa_100",
            "preferred_safe_fraction": preferred["safe_fraction"],
            "preferred_baseline_global_S_T_max": preferred["baseline"]["global_S_T_max"],
            "preferred_baseline_margin": preferred["baseline"]["worst_margin"],
            "interpretation": (
                "The kappa=100 two-sided filter remains safe at S_T=1e-5 for "
                "the unstressed spectrum grid.  Unsafe points, when xi_H is made "
                "large, locate the unresolved flavor-coherence risk in the RRRR "
                "higgsino dressing.  "
                "This supports the conditional heavy-spectrum branch but does not "
                "replace a full mass-insertion dressing calculation."
            ),
        },
    }


def trim_row(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "xi_H",
        "higgsino_to_wino",
        "higgsino_to_wino_stressed",
        "worst_channel",
        "worst_margin",
        "worst_tau_years_ST_display",
        "global_S_T_max",
        "all_channels_pass",
    ]
    return {key: row[key] for key in keys}


def min_safe_mass_table(safe_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for tan_beta in TAN_BETA_GRID:
        for mu_h in MU_H_GRID:
            rows = [row for row in safe_rows if row["tan_beta"] == tan_beta and row["mu_H_GeV"] == mu_h]
            table.append(
                {
                    "tan_beta": tan_beta,
                    "mu_H_GeV": mu_h,
                    "minimum_safe_m_sfermion_GeV": min((row["m_sfermion_GeV"] for row in rows), default=math.inf),
                }
            )
    return table


def safe_by_xi_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for xi_h in XI_H_GRID:
        subset = [row for row in rows if row["xi_H"] == xi_h]
        safe = [row for row in subset if row["all_channels_pass"]]
        unsafe = [row for row in subset if not row["all_channels_pass"]]
        table.append(
            {
                "xi_H": xi_h,
                "safe_points": len(safe),
                "total_points": len(subset),
                "safe_fraction": len(safe) / len(subset),
                "worst_safe_margin": min((row["worst_margin"] for row in safe), default=math.inf),
                "nearest_unsafe_margin": max((row["worst_margin"] for row in unsafe), default=None),
            }
        )
    return table


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "condition_cap",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "xi_H",
        "D_wino_GeV_minus1",
        "D_higgsino_proxy_GeV_minus1",
        "D_R_coherent_GeV_minus1",
        "D_R_stressed_GeV_minus1",
        "higgsino_to_wino",
        "higgsino_to_wino_stressed",
        "worst_channel",
        "worst_margin",
        "worst_tau_years_ST_display",
        "global_S_T_max",
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
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "xi_H",
        "channel_key",
        "physical_proxy",
        "dressing_class",
        "amplitude",
        "effective_dressing_GeV_minus1",
        "target_years",
        "tau_years_ST_display",
        "margin_at_ST_display",
        "S_T_max",
        "passes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for channel in row["channel_results"]:
                writer.writerow(
                    {
                        "filter_label": row["filter_label"],
                        "m_sfermion_GeV": row["m_sfermion_GeV"],
                        "m_wino_GeV": row["m_wino_GeV"],
                        "mu_H_GeV": row["mu_H_GeV"],
                        "tan_beta": row["tan_beta"],
                        "xi_H": row["xi_H"],
                        **channel,
                    }
                )


def write_report(summary: dict[str, Any]) -> None:
    def fmt_optional(value: Any) -> str:
        if value is None:
            return "--"
        return f"{float(value):.3e}"

    lines: list[str] = []
    lines.append("# Spectrum-aware dimension-five dressing proxy")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("The scan uses")
    lines.append("")
    lines.append("```text")
    lines.append("D_W = (alpha2/4pi) m_wino/m_sf^2 F(m_wino^2/m_sf^2)")
    lines.append("D_H = (1/16pi^2) mu_H/m_sf^2 y_t y_b (tan beta/10) F(mu_H^2/m_sf^2)")
    lines.append("D_L = D_W,  D_R = D_W + xi_H D_H")
    lines.append("```")
    lines.append("")
    lines.append("## Filter summaries")
    lines.append("")
    lines.append("| filter | safe / total | xi_H=1 safe | baseline ST max | baseline margin | marginal safe m_sf | nearest unsafe m_sf |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for label, payload in summary["by_filter_label"].items():
        marginal = payload["most_marginal_safe"]
        unsafe = payload["nearest_unsafe"]
        lines.append(
            f"| `{label}` | {payload['safe_points']}/{payload['total_points']} | "
            f"{payload['unstressed_safe_points']}/{payload['unstressed_total_points']} | "
            f"{payload['baseline']['global_S_T_max']:.3e} | "
            f"{payload['baseline']['worst_margin']:.3e} | "
            f"{marginal['m_sfermion_GeV'] if marginal else math.inf:.3e} | "
            f"{unsafe['m_sfermion_GeV'] if unsafe else math.inf:.3e} |"
        )
    lines.append("")
    lines.append("## Preferred kappa=100 coherent-higgsino stress")
    lines.append("")
    lines.append("| xi_H | safe / total | worst safe margin | nearest unsafe margin |")
    lines.append("|---:|---:|---:|---:|")
    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    for row in preferred["safe_by_xi_H"]:
        lines.append(
            f"| {row['xi_H']:.0f} | {row['safe_points']}/{row['total_points']} | "
            f"{row['worst_safe_margin']:.3e} | {fmt_optional(row['nearest_unsafe_margin'])} |"
        )
    lines.append("")
    lines.append("## Preferred kappa=100 minimum safe sfermion mass")
    lines.append("")
    lines.append("| tan beta | mu_H [GeV] | min safe m_sf [GeV] |")
    lines.append("|---:|---:|---:|")
    for row in preferred["minimum_safe_m_sfermion_by_tan_beta"]:
        lines.append(
            f"| {row['tan_beta']:.0f} | {row['mu_H_GeV']:.0f} | "
            f"{row['minimum_safe_m_sfermion_GeV']:.3e} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(summary["verdict"]["interpretation"])
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = audit_rows()
    write_csv(OUT / "spectrum_scan.csv", rows)
    write_channel_csv(OUT / "channel_scan.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    write_report(summary)

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    print("Spectrum-aware d=5 dressing proxy")
    print(
        "  omegaR_0.1_kappa_100: "
        f"safe={preferred['safe_points']}/{preferred['total_points']}, "
        f"baseline STmax={preferred['baseline']['global_S_T_max']:.3e}, "
        f"baseline margin={preferred['baseline']['worst_margin']:.3e}"
    )
    if preferred["nearest_unsafe"]:
        unsafe = preferred["nearest_unsafe"]
        print(
            "  nearest unsafe: "
            f"m_sf={unsafe['m_sfermion_GeV']:.3e}, "
            f"mu={unsafe['mu_H_GeV']:.3e}, tanb={unsafe['tan_beta']:.0f}, "
            f"margin={unsafe['worst_margin']:.3e}"
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
