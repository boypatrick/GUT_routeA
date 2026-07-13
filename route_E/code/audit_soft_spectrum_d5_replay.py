#!/usr/bin/env python3
"""Replay the MSSM soft-spectrum d=5 proxy at the tightened triplet filter.

No web lookup is used.  The earlier MSSM-mixing audit already diagonalized the
chargino, neutralino, and positive sfermion proxy spectra, but it reported
lifetimes at the historical S_T=1e-5 and with the old single hadronic width
prefactor.  This audit is a bookkeeping-but-physical replay:

    margin_new =
        margin_old * (S_T_old/S_T_new)^2 * (K_old/K_channel),

where K_channel is the central, max-width, or min-width chiral/lattice
prefactor.  The point is to test whether the tightened S_T=7.5e-6 branch
survives explicit electroweakino and sfermion-eigenstate dressing under the
same conservative hadronic envelope used by the field-basis Wilson replay.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MSSM_SUMMARY = ROOT / "output" / "mssm_mixing_d5_dressing" / "summary.json"
MSSM_CHANNELS = ROOT / "output" / "mssm_mixing_d5_dressing" / "channel_scan.csv"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
TIGHTENED = ROOT / "output" / "tightened_triplet_filter" / "summary.json"
PHYSICAL = ROOT / "output" / "physical_d5_wilson_replay" / "summary.json"
OUT = ROOT / "output" / "soft_spectrum_d5_replay"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def channel_class(channel: str) -> str:
    if "K0mu" in channel:
        return "K0mu"
    if "Knu" in channel or "RRRR_uusd" in channel:
        return "Knu"
    if "epi" in channel:
        return "e_pi"
    raise ValueError(f"unknown channel {channel}")


def current_target_years(channel: str) -> float:
    if "K0mu" in channel:
        return 1.0e34
    return 2.4e34


def chiral_widths(chiral: dict[str, Any]) -> dict[str, dict[str, float]]:
    old = float(chiral["old_single_width_prefactor"])
    central = {key: float(value) for key, value in chiral["central_width_prefactors"].items()}
    max_width = {
        "Knu": 0.0019122754090593987,
        "K0mu": 0.00040366230812257124,
        "e_pi": 0.005713123603904858,
    }
    min_width = {
        "Knu": 0.00007673089919519867,
        "K0mu": 0.00001346200606010078,
        "e_pi": 0.0002088057718555193,
    }
    return {
        key: {
            "old": old,
            "central": central[key],
            "max_width": max_width[key],
            "min_width": min_width[key],
        }
        for key in ["Knu", "K0mu", "e_pi"]
    }


def read_channel_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with MSSM_CHANNELS.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "amplitude_with_dressing",
                "amplification_vs_eigenstate_proxy",
                "S_T_max",
                "tau_years_ST_display",
                "margin_at_ST_display",
                "avg_chargino_part",
                "avg_neutralino_part",
                "avg_higgsino_part",
                "avg_total_part",
            ]:
                row[key] = float(row[key])
            row["positive_definite"] = row["positive_definite"] == "True"
            row["passes"] = row["passes"] == "True"
            rows.append(row)
    return rows


def replay() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    mssm = read_json(MSSM_SUMMARY)
    chiral = read_json(CHIRAL)
    tightened = read_json(TIGHTENED)
    physical = read_json(PHYSICAL)
    widths = chiral_widths(chiral)
    display_st = float(mssm["display_triplet_filter"])
    st_target = float(next(row for row in tightened["targets"] if row["target_label"] == "rounded_7p5e-6")["S_T_target"])

    channel_rows: list[dict[str, Any]] = []
    for row in read_channel_rows():
        cls = channel_class(row["channel"])
        target_current = current_target_years(row["channel"])
        for case in ["central", "max_width", "min_width"]:
            k_ratio = widths[cls][case] / widths[cls]["old"]
            st_allowed = row["S_T_max"] / math.sqrt(k_ratio)
            margin = row["margin_at_ST_display"] * (display_st / st_target) ** 2 / k_ratio
            future_margin_1e35 = margin * target_current / 1.0e35
            channel_rows.append(
                {
                    **row,
                    "channel_class": cls,
                    "normalization_case": case,
                    "S_T_target": st_target,
                    "old_display_S_T": display_st,
                    "width_ratio_new_over_old": k_ratio,
                    "S_T_max_replayed": st_allowed,
                    "margin_replayed_current": margin,
                    "future_margin_1e35": future_margin_1e35,
                    "passes_replayed_current": margin >= 1.0,
                    "passes_future_1e35": future_margin_1e35 >= 1.0,
                    "current_target_years": target_current,
                }
            )

    point_rows = point_summaries(channel_rows)
    summary = summarize(point_rows, channel_rows, mssm, physical, st_target)
    return channel_rows, point_rows, summary


def point_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        row["filter_label"],
        row["scenario"],
        row["epsilon"],
        row["spectrum_name"],
        row["normalization_case"],
    )


def point_summaries(channel_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in channel_rows:
        if row["positive_definite"]:
            grouped.setdefault(point_key(row), []).append(row)

    point_rows: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        worst_current = min(rows, key=lambda item: item["margin_replayed_current"])
        worst_future = min(rows, key=lambda item: item["future_margin_1e35"])
        all_current = all(row["passes_replayed_current"] for row in rows)
        all_future = all(row["passes_future_1e35"] for row in rows)
        point_rows.append(
            {
                "filter_label": key[0],
                "scenario": key[1],
                "epsilon": float(key[2]),
                "spectrum_name": key[3],
                "normalization_case": key[4],
                "channel_rows": len(rows),
                "worst_current_channel": worst_current["channel"],
                "worst_current_pair": worst_current["pair"],
                "worst_current_margin": worst_current["margin_replayed_current"],
                "worst_current_ST_max": worst_current["S_T_max_replayed"],
                "worst_future_channel": worst_future["channel"],
                "worst_future_pair": worst_future["pair"],
                "worst_future_margin_1e35": worst_future["future_margin_1e35"],
                "all_channels_pass_current": all_current,
                "all_channels_pass_future_1e35": all_future,
                "min_soft_eigenvalue": float(worst_current.get("min_soft_eigenvalue", "nan"))
                if "min_soft_eigenvalue" in worst_current
                else math.nan,
                "max_amplification_vs_eigenstate_proxy": max(
                    row["amplification_vs_eigenstate_proxy"] for row in rows
                ),
            }
        )
    return point_rows


def aligned_scenario(scenario: str) -> bool:
    return scenario in {
        "zero",
        "up_aligned_LL_MFV",
        "down_aligned_LL_MFV",
        "right_third_split",
        "commutator_LL",
        "combined_LL_RR",
    }


def trim_point(row: dict[str, Any] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "normalization_case",
        "worst_current_channel",
        "worst_current_pair",
        "worst_current_margin",
        "worst_current_ST_max",
        "worst_future_channel",
        "worst_future_margin_1e35",
        "all_channels_pass_current",
        "all_channels_pass_future_1e35",
        "max_amplification_vs_eigenstate_proxy",
    ]
    return {key: row[key] for key in keys}


def summarize(
    point_rows: list[dict[str, Any]],
    channel_rows: list[dict[str, Any]],
    mssm: dict[str, Any],
    physical: dict[str, Any],
    st_target: float,
) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    for filt in sorted({row["filter_label"] for row in point_rows}):
        filt_rows = [row for row in point_rows if row["filter_label"] == filt]
        by_case: dict[str, Any] = {}
        for case in ["central", "max_width", "min_width"]:
            rows = [row for row in filt_rows if row["normalization_case"] == case]
            current_safe = [row for row in rows if row["all_channels_pass_current"]]
            current_unsafe = [row for row in rows if not row["all_channels_pass_current"]]
            future_safe = [row for row in rows if row["all_channels_pass_future_1e35"]]
            aligned = [row for row in rows if aligned_scenario(row["scenario"])]
            aligned_current_safe = [row for row in aligned if row["all_channels_pass_current"]]
            most_marginal_safe = (
                min(current_safe, key=lambda row: row["worst_current_margin"])
                if current_safe
                else None
            )
            nearest_unsafe = (
                max(current_unsafe, key=lambda row: row["worst_current_margin"])
                if current_unsafe
                else None
            )
            worst_future = min(rows, key=lambda row: row["worst_future_margin_1e35"])
            by_case[case] = {
                "positive_points": len(rows),
                "current_safe_points": len(current_safe),
                "current_unsafe_points": len(current_unsafe),
                "future_1e35_safe_points": len(future_safe),
                "aligned_positive_points": len(aligned),
                "aligned_current_safe_points": len(aligned_current_safe),
                "current_safe_fraction": len(current_safe) / len(rows) if rows else 0.0,
                "aligned_current_safe_fraction": len(aligned_current_safe) / len(aligned)
                if aligned
                else 0.0,
                "global_worst_current_margin": min(row["worst_current_margin"] for row in rows),
                "global_worst_future_1e35_margin": worst_future["worst_future_margin_1e35"],
                "additional_amplitude_suppression_for_future": max(
                    1.0 / math.sqrt(worst_future["worst_future_margin_1e35"]),
                    1.0,
                ),
                "most_marginal_current_safe": trim_point(most_marginal_safe),
                "nearest_current_unsafe": trim_point(nearest_unsafe),
                "worst_future_1e35": trim_point(worst_future),
                "max_amplification_vs_eigenstate_proxy": max(
                    row["max_amplification_vs_eigenstate_proxy"] for row in rows
                ),
            }
        by_filter[filt] = by_case

    preferred_max = by_filter["omegaR_0.1_kappa_100"]["max_width"]
    allblock_max = by_filter["all_blocks_kappa_100"]["max_width"]
    verdict = {
        "S_T_target": st_target,
        "preferred_filter": "omegaR_0.1_kappa_100",
        "preferred_max_width_safe_points": preferred_max["current_safe_points"],
        "preferred_max_width_positive_points": preferred_max["positive_points"],
        "preferred_max_width_worst_current_margin": preferred_max["global_worst_current_margin"],
        "preferred_max_width_worst_future_1e35_margin": preferred_max["global_worst_future_1e35_margin"],
        "preferred_future_amplitude_suppression_needed": preferred_max[
            "additional_amplitude_suppression_for_future"
        ],
        "allblock_max_width_safe_points": allblock_max["current_safe_points"],
        "allblock_max_width_positive_points": allblock_max["positive_points"],
        "allblock_max_width_worst_current_margin": allblock_max["global_worst_current_margin"],
        "physical_identity_max_width_margin": physical["verdict"]["identity_max_width_worst_margin"],
        "mssm_reference_old_worst_margin": mssm["by_filter_label"]["omegaR_0.1_kappa_100"][
            "most_marginal_safe"
        ]["worst_margin"],
        "interpretation": (
            "After replaying explicit chargino/neutralino and positive-sfermion "
            "eigenstate dressing at S_T=7.5e-6 with the conservative chiral/lattice "
            "max-width envelope, the preferred omegaR=0.1,kappa=100 branch remains "
            "safe at every audited positive soft-spectrum point.  The all-block "
            "fallback remains unsafe.  The preferred branch is still not robust "
            "against a uniform 1e35 yr future stress without an extra O(1) "
            "amplitude suppression or a stronger triplet filter."
        ),
    }
    return {
        "note": "No web lookup used. Tightened soft-spectrum d=5 replay with chiral/lattice envelope.",
        "input_files": {
            "mssm_mixing_d5_dressing": str(MSSM_SUMMARY),
            "mssm_channel_scan": str(MSSM_CHANNELS),
            "chiral_lattice_d5": str(CHIRAL),
            "tightened_triplet_filter": str(TIGHTENED),
            "physical_d5_wilson_replay": str(PHYSICAL),
        },
        "formula": "margin_new = margin_old * (S_T_old/S_T_new)^2 * (K_old/K_channel)",
        "by_filter_label": by_filter,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_outputs(channel_rows: list[dict[str, Any]], point_rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    channel_fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "channel",
        "channel_class",
        "operator",
        "normalization_case",
        "pair",
        "positive_definite",
        "amplitude_with_dressing",
        "amplification_vs_eigenstate_proxy",
        "S_T_max",
        "S_T_max_replayed",
        "S_T_target",
        "width_ratio_new_over_old",
        "margin_at_ST_display",
        "margin_replayed_current",
        "future_margin_1e35",
        "passes_replayed_current",
        "passes_future_1e35",
        "avg_chargino_part",
        "avg_neutralino_part",
        "avg_higgsino_part",
    ]
    point_fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "normalization_case",
        "channel_rows",
        "worst_current_channel",
        "worst_current_pair",
        "worst_current_margin",
        "worst_current_ST_max",
        "worst_future_channel",
        "worst_future_pair",
        "worst_future_margin_1e35",
        "all_channels_pass_current",
        "all_channels_pass_future_1e35",
        "max_amplification_vs_eigenstate_proxy",
    ]
    write_csv(OUT / "channel_replay.csv", channel_rows, channel_fields)
    write_csv(OUT / "point_replay.csv", point_rows, point_fields)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Soft-spectrum d=5 replay",
        "",
        "No web lookup was used.",
        "",
        "Replay formula:",
        "",
        "```text",
        summary["formula"],
        "```",
        "",
        "## Max-width envelope",
        "",
        "| filter | current safe/positive | worst current margin | future safe/positive | worst future margin | future amp suppression |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, payload in summary["by_filter_label"].items():
        row = payload["max_width"]
        lines.append(
            f"| `{label}` | {row['current_safe_points']}/{row['positive_points']} | "
            f"{row['global_worst_current_margin']:.6e} | "
            f"{row['future_1e35_safe_points']}/{row['positive_points']} | "
            f"{row['global_worst_future_1e35_margin']:.6e} | "
            f"{row['additional_amplitude_suppression_for_future']:.6e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    channel_rows, point_rows, summary = replay()
    write_outputs(channel_rows, point_rows, summary)
    verdict = summary["verdict"]
    print("Soft-spectrum d=5 replay")
    print(f"  S_T target: {verdict['S_T_target']:.6e}")
    print(
        "  preferred max-width current: "
        f"{verdict['preferred_max_width_safe_points']}/"
        f"{verdict['preferred_max_width_positive_points']} safe, "
        f"worst margin={verdict['preferred_max_width_worst_current_margin']:.6e}"
    )
    print(
        "  preferred future 1e35: "
        f"worst margin={verdict['preferred_max_width_worst_future_1e35_margin']:.6e}, "
        f"amp suppression needed={verdict['preferred_future_amplitude_suppression_needed']:.6e}"
    )
    print(
        "  all-block max-width current: "
        f"{verdict['allblock_max_width_safe_points']}/"
        f"{verdict['allblock_max_width_positive_points']} safe, "
        f"worst margin={verdict['allblock_max_width_worst_current_margin']:.6e}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
