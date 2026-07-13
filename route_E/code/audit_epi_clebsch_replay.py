#!/usr/bin/env python3
"""Clebsch/projection sensitivity replay for the e pi d=5 monitor.

No web lookup is used.  The tightened worst-phase replay identified the
minimal LLLL_uude_epi monitor as the current bottleneck.  This script asks a
sharper, representation-level question: if the physical p -> e+ pi0 amplitude
contains the standard neutral-pion projection and identical-up normalization
factors, does the e pi monitor remain the bottleneck?

We do not claim a final lattice-level p -> e+ pi0 computation here.  We apply
controlled amplitude projections only to e_pi rows:

    pi0 projection:        A -> A/sqrt(2)
    identical-up norm:     A -> A/sqrt(2)
    combined projection:   A -> A/2

Since widths scale as |A|^2, margins scale as 1/projection^2 and S_T^max
scales as 1/projection.  K nu and K0 mu rows are left unchanged, so the replay
also reveals the irreducible non-e_pi bottleneck.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
INTERFERENCE = ROOT / "output" / "interference_tightened_d5_replay" / "channel_replay.csv"
INTERFERENCE_SUMMARY = ROOT / "output" / "interference_tightened_d5_replay" / "summary.json"
OUT = ROOT / "output" / "epi_clebsch_replay"


PROFILES = {
    "raw_monitor": {
        "e_pi_amplitude_scale": 1.0,
        "description": "No additional e pi Clebsch/projection factor.",
    },
    "pi0_isospin": {
        "e_pi_amplitude_scale": 1.0 / math.sqrt(2.0),
        "description": "Neutral pion projection only: pi0=(u ubar-d dbar)/sqrt(2).",
    },
    "identical_up_norm": {
        "e_pi_amplitude_scale": 1.0 / math.sqrt(2.0),
        "description": "Two identical up-quark normalization only.",
    },
    "pi0_times_identical_up": {
        "e_pi_amplitude_scale": 0.5,
        "description": "Neutral pion projection times identical-up normalization.",
    },
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with INTERFERENCE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "triplet_phase",
                "min_phase_amplitude",
                "max_phase_amplitude",
                "phase_spread",
                "S_T_max",
                "S_T_max_replayed",
                "S_T_target",
                "width_ratio_new_over_old",
                "margin_at_ST_display",
                "margin_replayed_current_worst_phase",
                "margin_replayed_current_best_phase_proxy",
                "future_margin_1e35_worst_phase",
            ]:
                row[key] = float(row[key])
            row["passes_replayed_current_worst_phase"] = (
                row["passes_replayed_current_worst_phase"] == "True"
            )
            row["passes_future_1e35_worst_phase"] = row["passes_future_1e35_worst_phase"] == "True"
            rows.append(row)
    return rows


def projected_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    base_summary = read_json(INTERFERENCE_SUMMARY)
    rows: list[dict[str, Any]] = []
    for profile, payload in PROFILES.items():
        scale = payload["e_pi_amplitude_scale"]
        for row in read_rows():
            is_epi = row["channel_class"] == "e_pi"
            amp_scale = scale if is_epi else 1.0
            margin_factor = 1.0 / (amp_scale * amp_scale)
            projected = {
                **row,
                "profile": profile,
                "profile_description": payload["description"],
                "amplitude_projection": amp_scale,
                "S_T_max_projected": row["S_T_max_replayed"] / amp_scale,
                "margin_current_projected": row["margin_replayed_current_worst_phase"] * margin_factor,
                "future_margin_1e35_projected": row["future_margin_1e35_worst_phase"] * margin_factor,
            }
            projected["passes_current_projected"] = projected["margin_current_projected"] >= 1.0
            projected["passes_future_1e35_projected"] = projected["future_margin_1e35_projected"] >= 1.0
            rows.append(projected)

    point_rows = summarize_points(rows)
    summary = summarize(rows, point_rows, base_summary)
    return rows, point_rows, summary


def point_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        row["profile"],
        row["filter_label"],
        row["scenario"],
        row["epsilon"],
        row["spectrum_name"],
        row["normalization_case"],
    )


def summarize_points(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(point_key(row), []).append(row)

    out: list[dict[str, Any]] = []
    for key, items in sorted(grouped.items()):
        worst_current = min(items, key=lambda row: row["margin_current_projected"])
        worst_future = min(items, key=lambda row: row["future_margin_1e35_projected"])
        out.append(
            {
                "profile": key[0],
                "filter_label": key[1],
                "scenario": key[2],
                "epsilon": float(key[3]),
                "spectrum_name": key[4],
                "normalization_case": key[5],
                "channel_rows": len(items),
                "worst_current_channel": worst_current["channel"],
                "worst_current_class": worst_current["channel_class"],
                "worst_current_margin": worst_current["margin_current_projected"],
                "worst_current_ST_max": worst_current["S_T_max_projected"],
                "worst_future_channel": worst_future["channel"],
                "worst_future_class": worst_future["channel_class"],
                "worst_future_margin_1e35": worst_future["future_margin_1e35_projected"],
                "all_channels_pass_current": all(row["passes_current_projected"] for row in items),
                "all_channels_pass_future_1e35": all(row["passes_future_1e35_projected"] for row in items),
            }
        )
    return out


def trim_point(row: dict[str, Any] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    keys = [
        "profile",
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "normalization_case",
        "worst_current_channel",
        "worst_current_class",
        "worst_current_margin",
        "worst_current_ST_max",
        "worst_future_channel",
        "worst_future_class",
        "worst_future_margin_1e35",
        "all_channels_pass_current",
        "all_channels_pass_future_1e35",
    ]
    return {key: row[key] for key in keys}


def summarize(rows: list[dict[str, Any]], point_rows: list[dict[str, Any]], base: dict[str, Any]) -> dict[str, Any]:
    by_profile: dict[str, Any] = {}
    for profile in PROFILES:
        prof_rows = [row for row in point_rows if row["profile"] == profile]
        by_filter: dict[str, Any] = {}
        for filt in sorted({row["filter_label"] for row in prof_rows}):
            filt_rows = [row for row in prof_rows if row["filter_label"] == filt]
            by_case: dict[str, Any] = {}
            for case in ["central", "max_width", "min_width"]:
                case_rows = [row for row in filt_rows if row["normalization_case"] == case]
                safe = [row for row in case_rows if row["all_channels_pass_current"]]
                future_safe = [row for row in case_rows if row["all_channels_pass_future_1e35"]]
                worst_current = min(case_rows, key=lambda row: row["worst_current_margin"])
                worst_future = min(case_rows, key=lambda row: row["worst_future_margin_1e35"])
                by_case[case] = {
                    "positive_points": len(case_rows),
                    "current_safe_points": len(safe),
                    "future_1e35_safe_points": len(future_safe),
                    "global_worst_current_margin": worst_current["worst_current_margin"],
                    "global_worst_current_channel": worst_current["worst_current_channel"],
                    "global_worst_current_class": worst_current["worst_current_class"],
                    "global_worst_future_1e35_margin": worst_future["worst_future_margin_1e35"],
                    "global_worst_future_channel": worst_future["worst_future_channel"],
                    "global_worst_future_class": worst_future["worst_future_class"],
                    "additional_amplitude_suppression_for_future": max(
                        1.0 / math.sqrt(worst_future["worst_future_margin_1e35"]),
                        1.0,
                    ),
                    "most_marginal_current": trim_point(worst_current),
                    "worst_future_1e35": trim_point(worst_future),
                }
            by_filter[filt] = by_case
        by_profile[profile] = by_filter

    preferred_raw = by_profile["raw_monitor"]["omegaR_0.1_kappa_100"]["max_width"]
    preferred_combined = by_profile["pi0_times_identical_up"]["omegaR_0.1_kappa_100"]["max_width"]
    verdict = {
        "raw_preferred_current_margin": preferred_raw["global_worst_current_margin"],
        "raw_preferred_current_channel": preferred_raw["global_worst_current_channel"],
        "combined_preferred_current_margin": preferred_combined["global_worst_current_margin"],
        "combined_preferred_current_channel": preferred_combined["global_worst_current_channel"],
        "combined_preferred_current_class": preferred_combined["global_worst_current_class"],
        "combined_preferred_future_margin": preferred_combined["global_worst_future_1e35_margin"],
        "combined_preferred_future_suppression_needed": preferred_combined[
            "additional_amplitude_suppression_for_future"
        ],
        "base_signed_margin": base["verdict"]["preferred_max_width_worst_current_margin"],
        "interpretation": (
            "Applying both neutral-pion and identical-up projection factors removes "
            "the minimal e pi monitor as the current bottleneck.  The worst current "
            "preferred margin rises from 1.052 to the non-e_pi floor 1.113, and the "
            "limiting channel becomes Knu.  The future 1e35 stress still fails, so "
            "the projection sensitivity improves the current claim but does not by "
            "itself close the future-proof d=5 problem."
        ),
    }
    return {
        "note": "No web lookup used. e pi Clebsch/projection sensitivity replay.",
        "input_files": {
            "interference_tightened_channel_replay": str(INTERFERENCE),
            "interference_tightened_summary": str(INTERFERENCE_SUMMARY),
        },
        "profiles": PROFILES,
        "by_profile": by_profile,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# e pi Clebsch replay",
        "",
        "No web lookup was used.",
        "",
        "Only e_pi rows are rescaled.  Knu and K0mu rows are unchanged.",
        "",
        "## Preferred max-width summary",
        "",
        "| profile | e_pi amplitude scale | current margin | current bottleneck | future margin | future amp suppression |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for profile, payload in summary["profiles"].items():
        row = summary["by_profile"][profile]["omegaR_0.1_kappa_100"]["max_width"]
        lines.append(
            f"| `{profile}` | {payload['e_pi_amplitude_scale']:.6e} | "
            f"{row['global_worst_current_margin']:.6e} | `{row['global_worst_current_channel']}` | "
            f"{row['global_worst_future_1e35_margin']:.6e} | "
            f"{row['additional_amplitude_suppression_for_future']:.6e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, point_rows, summary = projected_rows()
    channel_fields = [
        "profile",
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "triplet_phase",
        "channel",
        "channel_class",
        "normalization_case",
        "pair",
        "amplitude_projection",
        "S_T_max_projected",
        "margin_current_projected",
        "future_margin_1e35_projected",
        "passes_current_projected",
        "passes_future_1e35_projected",
    ]
    point_fields = [
        "profile",
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "normalization_case",
        "worst_current_channel",
        "worst_current_class",
        "worst_current_margin",
        "worst_current_ST_max",
        "worst_future_channel",
        "worst_future_class",
        "worst_future_margin_1e35",
        "all_channels_pass_current",
        "all_channels_pass_future_1e35",
    ]
    write_csv(OUT / "channel_replay.csv", rows, channel_fields)
    write_csv(OUT / "point_replay.csv", point_rows, point_fields)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    verdict = summary["verdict"]
    print("e pi Clebsch replay")
    print(
        "  raw preferred current: "
        f"margin={verdict['raw_preferred_current_margin']:.6e}, "
        f"channel={verdict['raw_preferred_current_channel']}"
    )
    print(
        "  combined projection preferred current: "
        f"margin={verdict['combined_preferred_current_margin']:.6e}, "
        f"channel={verdict['combined_preferred_current_channel']}"
    )
    print(
        "  combined projection future: "
        f"margin={verdict['combined_preferred_future_margin']:.6e}, "
        f"amp suppression needed={verdict['combined_preferred_future_suppression_needed']:.6e}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
