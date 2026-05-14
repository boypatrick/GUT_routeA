#!/usr/bin/env python3
"""Channel-specific dressed dimension-five proton-decay proxy.

No web lookup is used.  Earlier scripts built mass-basis d=5 Wilson tensors
and a two-sided triplet filter.  This audit takes the monitored channel
amplitudes from that filter and applies explicit dressing factors channel by
channel:

    C5_ch = S_T A_ch / M_T,
    C6_ch = C5_ch (alpha_2/4 pi) m_wino/m_sf^2 * d_ch,
    Gamma_ch = K_ch |C6_ch|^2.

The result is still a controlled proxy, not a final SUSY-GUT proton-decay
calculation.  It exposes which channels fail under stress dressings and records
the triplet-filter upper bound required per channel.
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
OUT = ROOT / "output" / "dressed_dimension5_channels"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
DISPLAY_ST = 1.0e-5


CHANNELS = {
    "LLLL_upupdown_Knu": {
        "physical_proxy": "p_to_K+nu_bar_L_upupdown",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
    },
    "LLLL_downdownup_Knu": {
        "physical_proxy": "p_to_K+nu_bar_L_downdownup",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
    },
    "LLLL_upupdown_K0mu": {
        "physical_proxy": "p_to_K0mu+_L_proxy",
        "target_years": 1.0e34,
        "hadronic_relative": 1.0,
    },
    "RRRR_uusd_anycharged": {
        "physical_proxy": "p_to_K+nu_or_K0mu_R_proxy",
        "target_years": 2.4e34,
        "hadronic_relative": 1.0,
    },
}


DRESSINGS = [
    {
        "name": "wino_100TeV",
        "m_wino_GeV": 1.0e3,
        "m_sfermion_GeV": 1.0e5,
        "alpha2_inv": 25.0,
        "relative_dressing": 1.0,
        "interpretation": "baseline local item-4 dressing",
    },
    {
        "name": "wino_10TeV",
        "m_wino_GeV": 1.0e3,
        "m_sfermion_GeV": 1.0e4,
        "alpha2_inv": 25.0,
        "relative_dressing": 1.0,
        "interpretation": "light-sfermion stress test",
    },
    {
        "name": "wino_100TeV_dressing_x10",
        "m_wino_GeV": 1.0e3,
        "m_sfermion_GeV": 1.0e5,
        "alpha2_inv": 25.0,
        "relative_dressing": 10.0,
        "interpretation": "constructive dressing or short-distance enhancement stress",
    },
    {
        "name": "wino_10TeV_dressing_x100",
        "m_wino_GeV": 1.0e3,
        "m_sfermion_GeV": 1.0e4,
        "alpha2_inv": 25.0,
        "relative_dressing": 100.0,
        "interpretation": "extreme light-spectrum and dressing-enhancement stress",
    },
]


FILTER_LABELS = [
    "omegaR_0.1_kappa_30",
    "omegaR_0.1_kappa_100",
    "omegaR_0.1_kappa_300",
    "all_blocks_kappa_100",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def width_to_tau(width_gev: float) -> float:
    if width_gev <= 0.0:
        return math.inf
    return HBAR_GEV_S / width_gev / SECONDS_PER_YEAR


def dressed_lifetime(
    amplitude: float,
    st: float,
    m_triplet_gev: float,
    width_prefactor_gev5: float,
    dressing: dict[str, float],
    hadronic_relative: float,
) -> dict[str, float]:
    loop = (1.0 / dressing["alpha2_inv"]) / (4.0 * math.pi)
    dress = (
        loop
        * dressing["m_wino_GeV"]
        / (dressing["m_sfermion_GeV"] * dressing["m_sfermion_GeV"])
        * dressing["relative_dressing"]
    )
    c5 = st * amplitude / m_triplet_gev
    c6 = c5 * dress
    width = width_prefactor_gev5 * hadronic_relative * c6 * c6
    return {
        "loop_alpha2_over_4pi": loop,
        "dressing_GeV_minus1": dress,
        "C5_GeV_minus1": c5,
        "C6_GeV_minus2": c6,
        "width_GeV": width,
        "tau_years": width_to_tau(width),
    }


def st_max(
    amplitude: float,
    target_years: float,
    m_triplet_gev: float,
    width_prefactor_gev5: float,
    dressing: dict[str, float],
    hadronic_relative: float,
) -> float:
    unit = dressed_lifetime(
        amplitude,
        1.0,
        m_triplet_gev,
        width_prefactor_gev5,
        dressing,
        hadronic_relative,
    )
    return math.sqrt(unit["tau_years"] / target_years)


def selected_filter_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = {row["label"]: row for row in payload["rows"]}
    rows.update({row["label"]: row for row in payload["all_block_rows"]})
    return [rows[label] for label in FILTER_LABELS]


def benchmark_triplet_mass() -> float:
    payload = read_json(VACUUM)
    return float(payload["recommended_benchmark"]["M_HC_GeV"])


def audit_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    two_sided = read_json(TWO_SIDED)
    constants = read_json(PROTON)["hadronic_constants"]
    m_triplet = benchmark_triplet_mass()
    rows: list[dict[str, Any]] = []
    for filt in selected_filter_rows(two_sided):
        for channel_key, channel in CHANNELS.items():
            amplitude = float(filt["channels"][channel_key]["amplitude"])
            for dressing in DRESSINGS:
                tau_display = dressed_lifetime(
                    amplitude,
                    DISPLAY_ST,
                    m_triplet,
                    constants["width_prefactor_GeV5"],
                    dressing,
                    channel["hadronic_relative"],
                )
                tau_unfiltered = dressed_lifetime(
                    amplitude,
                    1.0,
                    m_triplet,
                    constants["width_prefactor_GeV5"],
                    dressing,
                    channel["hadronic_relative"],
                )
                st_allowed = st_max(
                    amplitude,
                    channel["target_years"],
                    m_triplet,
                    constants["width_prefactor_GeV5"],
                    dressing,
                    channel["hadronic_relative"],
                )
                rows.append(
                    {
                        "filter_label": filt["label"],
                        "condition_cap": filt["condition_cap"],
                        "channel_key": channel_key,
                        "physical_proxy": channel["physical_proxy"],
                        "amplitude": amplitude,
                        "dressing": dressing["name"],
                        "dressing_interpretation": dressing["interpretation"],
                        "M_T_GeV": m_triplet,
                        "target_years": channel["target_years"],
                        "tau_years_ST_display": tau_display["tau_years"],
                        "tau_years_ST_1": tau_unfiltered["tau_years"],
                        "margin_at_ST_display": tau_display["tau_years"] / channel["target_years"],
                        "S_T_display": DISPLAY_ST,
                        "S_T_max": st_allowed,
                        "passes_display_ST": tau_display["tau_years"] >= channel["target_years"],
                        "passes_unfiltered": tau_unfiltered["tau_years"] >= channel["target_years"],
                        "C6_display_GeV_minus2": tau_display["C6_GeV_minus2"],
                    }
                )
    summary = summarize(rows, m_triplet)
    return rows, summary


def summarize(rows: list[dict[str, Any]], m_triplet: float) -> dict[str, Any]:
    display_rows = [
        row
        for row in rows
        if row["filter_label"] == "omegaR_0.1_kappa_100"
        and row["dressing"] == "wino_100TeV"
    ]
    worst_display = min(display_rows, key=lambda row: row["margin_at_ST_display"])
    stress_rows = [
        row
        for row in rows
        if row["filter_label"] == "omegaR_0.1_kappa_100"
        and row["dressing"] == "wino_10TeV_dressing_x100"
    ]
    worst_stress = min(stress_rows, key=lambda row: row["margin_at_ST_display"])
    all_block_rows = [
        row
        for row in rows
        if row["filter_label"] == "all_blocks_kappa_100"
        and row["dressing"] == "wino_100TeV"
    ]
    worst_all_block = min(all_block_rows, key=lambda row: row["margin_at_ST_display"])
    return {
        "note": "No web lookup used. Channel-specific dressed d=5 proxy; not a final chiral/lattice calculation.",
        "formula": {
            "C5": "S_T A_ch / M_T",
            "C6": "C5 (alpha2/4pi) m_wino/m_sf^2 * relative_dressing",
            "Gamma": "K_ch |C6|^2",
        },
        "display_triplet_filter": DISPLAY_ST,
        "M_T_GeV_from_R200_HC": m_triplet,
        "verdict": {
            "omegaR_0p1_kappa100_wino100TeV_passes_all_display_ST": all(
                row["passes_display_ST"] for row in display_rows
            ),
            "omegaR_0p1_kappa100_wino100TeV_worst": {
                key: worst_display[key]
                for key in [
                    "channel_key",
                    "physical_proxy",
                    "amplitude",
                    "S_T_max",
                    "tau_years_ST_display",
                    "margin_at_ST_display",
                ]
            },
            "all_blocks_kappa100_wino100TeV_worst": {
                key: worst_all_block[key]
                for key in [
                    "channel_key",
                    "physical_proxy",
                    "amplitude",
                    "S_T_max",
                    "tau_years_ST_display",
                    "margin_at_ST_display",
                ]
            },
            "extreme_stress_kappa100_worst": {
                key: worst_stress[key]
                for key in [
                    "channel_key",
                    "physical_proxy",
                    "amplitude",
                    "S_T_max",
                    "tau_years_ST_display",
                    "margin_at_ST_display",
                    "passes_display_ST",
                ]
            },
            "interpretation": (
                "With the current display filter S_T=1e-5, the kappa=100 two-sided "
                "filter passes the monitored Knu, K0mu, and RRRR proxy channels under "
                "the baseline 100 TeV wino dressing by a large margin.  The same branch "
                "is not robust against an extreme simultaneous 10 TeV sfermion and "
                "100x dressing-enhancement stress, so this remains a proxy until a "
                "spectrum-specific wino/higgsino calculation is supplied."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Channel-specific dressed dimension-five proxy")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("The audit uses")
    lines.append("")
    lines.append("```text")
    lines.append("C5_ch = S_T A_ch / M_T")
    lines.append("C6_ch = C5_ch (alpha2/4pi) m_wino/m_sf^2 * d_ch")
    lines.append("Gamma_ch = K_ch |C6_ch|^2")
    lines.append("```")
    lines.append("")
    lines.append(f"R=200 colored-triplet mass: `{summary['M_T_GeV_from_R200_HC']:.6e} GeV`.")
    lines.append(f"Displayed scalar filter: `S_T={summary['display_triplet_filter']:.1e}`.")
    lines.append("")
    lines.append("## Baseline kappa=100 two-sided filter")
    lines.append("")
    lines.append("| channel | amp | S_T max | tau(S_T=1e-5) [yr] | margin |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in rows:
        if row["filter_label"] == "omegaR_0.1_kappa_100" and row["dressing"] == "wino_100TeV":
            lines.append(
                f"| `{row['physical_proxy']}` | {row['amplitude']:.3e} | "
                f"{row['S_T_max']:.3e} | {row['tau_years_ST_display']:.3e} | "
                f"{row['margin_at_ST_display']:.3e} |"
            )
    lines.append("")
    lines.append("## Stress comparison")
    lines.append("")
    lines.append("| filter | dressing | worst channel | S_T max | tau(S_T=1e-5) [yr] | pass |")
    lines.append("|---|---|---|---:|---:|---|")
    for label in FILTER_LABELS:
        for dressing in ("wino_100TeV", "wino_10TeV", "wino_10TeV_dressing_x100"):
            subset = [row for row in rows if row["filter_label"] == label and row["dressing"] == dressing]
            worst = min(subset, key=lambda row: row["margin_at_ST_display"])
            lines.append(
                f"| `{label}` | `{dressing}` | `{worst['physical_proxy']}` | "
                f"{worst['S_T_max']:.3e} | {worst['tau_years_ST_display']:.3e} | "
                f"{worst['passes_display_ST']} |"
            )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(summary["verdict"]["interpretation"])
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = audit_rows()
    write_csv(OUT / "channel_dressing_rows.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(rows, summary)

    v = summary["verdict"]
    worst = v["omegaR_0p1_kappa100_wino100TeV_worst"]
    stress = v["extreme_stress_kappa100_worst"]
    print("Channel-specific dressed d=5 proxy")
    print(
        "  kappa=100 baseline worst: "
        f"{worst['physical_proxy']}, S_Tmax={worst['S_T_max']:.3e}, "
        f"tau(S_T=1e-5)={worst['tau_years_ST_display']:.3e} yr"
    )
    print(
        "  extreme stress worst: "
        f"{stress['physical_proxy']}, S_Tmax={stress['S_T_max']:.3e}, "
        f"pass={stress['passes_display_ST']}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
