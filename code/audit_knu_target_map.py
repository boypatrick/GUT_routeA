#!/usr/bin/env python3
"""K nu target map for the tightened d=5 proton-decay bottleneck.

No web lookup is used.  The preceding e-pi Clebsch replay showed that the
minimal e-pi monitor is not the stable preferred-filter bottleneck after
neutral-pion and identical-up projections.  This script therefore asks a more
useful design question:

    if K nu is the robust bottleneck, how strong must the next triplet filter
    or K nu-specific amplitude suppression be?

The input rows already include the conservative max-width chiral/lattice
envelope, worst scanned phases, explicit positive soft-spectrum dressing, and
the e-pi projection sensitivity.  This script does not claim a full
channel-specific hadronic calculation.  It converts the existing row-level
margins into reproducible targets:

    margin -> allowed amplitude headroom sqrt(margin)
    margin -> required future suppression 1/sqrt(margin_future)
    margin -> S_T,max = S_T,target sqrt(margin)
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import median
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "output" / "epi_clebsch_replay" / "channel_replay.csv"
EPI_SUMMARY = ROOT / "output" / "epi_clebsch_replay" / "summary.json"
OUT = ROOT / "output" / "knu_target_map"

PREFERRED_PROFILE = "pi0_times_identical_up"
PREFERRED_FILTER = "omegaR_0.1_kappa_100"
CONSERVATIVE_CASE = "max_width"


FLOAT_FIELDS = [
    "epsilon",
    "triplet_phase",
    "amplitude_projection",
    "S_T_max_projected",
    "margin_current_projected",
    "future_margin_1e35_projected",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with INPUT.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in FLOAT_FIELDS:
                row[key] = float(row[key])
            row["passes_current_projected"] = row["passes_current_projected"] == "True"
            row["passes_future_1e35_projected"] = row["passes_future_1e35_projected"] == "True"
            row["S_T_target_inferred"] = row["S_T_max_projected"] / math.sqrt(
                row["margin_current_projected"]
            )
            rows.append(row)
    return rows


def trim(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
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
        "S_T_max_projected",
        "margin_current_projected",
        "future_margin_1e35_projected",
    ]
    return {key: row[key] for key in keys}


def target_row(row: dict[str, Any], S_T_target: float) -> dict[str, Any]:
    current_margin = row["margin_current_projected"]
    future_margin = row["future_margin_1e35_projected"]
    return {
        **trim(row),
        "current_amplitude_headroom": math.sqrt(current_margin),
        "future_amplitude_suppression_needed": max(1.0 / math.sqrt(future_margin), 1.0),
        "S_T_safe_max_current": S_T_target * math.sqrt(current_margin),
        "S_T_safe_max_future_1e35": S_T_target * math.sqrt(future_margin),
    }


def summarize_channel(items: list[dict[str, Any]], S_T_target: float) -> dict[str, Any]:
    worst_current = min(items, key=lambda row: row["margin_current_projected"])
    worst_future = min(items, key=lambda row: row["future_margin_1e35_projected"])
    return {
        "channel": worst_current["channel"],
        "channel_class": worst_current["channel_class"],
        "rows": len(items),
        "worst_current": target_row(worst_current, S_T_target),
        "worst_future_1e35": target_row(worst_future, S_T_target),
    }


def build_summary(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    targets = [row["S_T_target_inferred"] for row in rows]
    S_T_target = median(targets)
    max_target_deviation = max(abs(target - S_T_target) for target in targets)

    preferred = [
        row
        for row in rows
        if row["profile"] == PREFERRED_PROFILE
        and row["filter_label"] == PREFERRED_FILTER
        and row["normalization_case"] == CONSERVATIVE_CASE
    ]
    if not preferred:
        raise RuntimeError("Preferred conservative rows were not found.")

    channel_targets: list[dict[str, Any]] = []
    for channel in sorted({row["channel"] for row in preferred}):
        channel_rows = [row for row in preferred if row["channel"] == channel]
        summary = summarize_channel(channel_rows, S_T_target)
        current = summary["worst_current"]
        future = summary["worst_future_1e35"]
        channel_targets.append(
            {
                "channel": summary["channel"],
                "channel_class": summary["channel_class"],
                "rows": summary["rows"],
                "worst_current_margin": current["margin_current_projected"],
                "current_amplitude_headroom": current["current_amplitude_headroom"],
                "S_T_safe_max_current": current["S_T_safe_max_current"],
                "worst_future_margin_1e35": future["future_margin_1e35_projected"],
                "future_amplitude_suppression_needed": future[
                    "future_amplitude_suppression_needed"
                ],
                "S_T_safe_max_future_1e35": future["S_T_safe_max_future_1e35"],
                "worst_current_scenario": current["scenario"],
                "worst_current_spectrum": current["spectrum_name"],
                "worst_current_epsilon": current["epsilon"],
                "worst_current_phase": current["triplet_phase"],
                "worst_future_scenario": future["scenario"],
                "worst_future_spectrum": future["spectrum_name"],
                "worst_future_epsilon": future["epsilon"],
                "worst_future_phase": future["triplet_phase"],
            }
        )

    global_current = target_row(
        min(preferred, key=lambda row: row["margin_current_projected"]), S_T_target
    )
    global_future = target_row(
        min(preferred, key=lambda row: row["future_margin_1e35_projected"]), S_T_target
    )
    knu_rows = [row for row in preferred if row["channel_class"] == "Knu"]
    knu_current = target_row(
        min(knu_rows, key=lambda row: row["margin_current_projected"]), S_T_target
    )
    knu_future = target_row(
        min(knu_rows, key=lambda row: row["future_margin_1e35_projected"]), S_T_target
    )

    future_safe_rows = [row for row in preferred if row["future_margin_1e35_projected"] >= 1.0]
    summary = {
        "note": "No web lookup used. Knu bottleneck target map.",
        "input_files": {
            "epi_clebsch_channel_replay": str(INPUT),
            "epi_clebsch_summary": str(EPI_SUMMARY),
        },
        "preferred_profile": PREFERRED_PROFILE,
        "preferred_filter": PREFERRED_FILTER,
        "conservative_case": CONSERVATIVE_CASE,
        "S_T_target_inferred": S_T_target,
        "S_T_target_max_deviation": max_target_deviation,
        "preferred_rows": len(preferred),
        "preferred_future_safe_rows": len(future_safe_rows),
        "global_worst_current": global_current,
        "global_worst_future_1e35": global_future,
        "Knu_worst_current": knu_current,
        "Knu_worst_future_1e35": knu_future,
        "common_Knu_future_amplitude_suppression_needed": knu_future[
            "future_amplitude_suppression_needed"
        ],
        "common_Knu_future_amplitude_scale_max": 1.0
        / knu_future["future_amplitude_suppression_needed"],
        "triplet_filter_S_T_safe_max_future_1e35": global_future[
            "S_T_safe_max_future_1e35"
        ],
        "triplet_filter_tightening_factor_from_7p5e_minus_6": S_T_target
        / global_future["S_T_safe_max_future_1e35"],
        "verdict": {
            "current_bottleneck_channel": global_current["channel"],
            "current_bottleneck_class": global_current["channel_class"],
            "current_margin": global_current["margin_current_projected"],
            "future_bottleneck_channel": global_future["channel"],
            "future_bottleneck_class": global_future["channel_class"],
            "future_margin": global_future["future_margin_1e35_projected"],
            "future_ST_target_required": global_future["S_T_safe_max_future_1e35"],
            "future_amplitude_suppression_needed": global_future[
                "future_amplitude_suppression_needed"
            ],
            "interpretation": (
                "After the e-pi projection sensitivity replay, the preferred "
                "conservative bottleneck is LLLL_upupdown_Knu.  The current "
                "2.4e34 yr stress passes with margin 1.113, but a uniform "
                "1e35 yr stress requires either S_T <= 3.876e-6 or a common "
                "Knu amplitude suppression of about 1.935."
            ),
        },
    }
    return channel_targets, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "channel",
        "channel_class",
        "rows",
        "worst_current_margin",
        "current_amplitude_headroom",
        "S_T_safe_max_current",
        "worst_future_margin_1e35",
        "future_amplitude_suppression_needed",
        "S_T_safe_max_future_1e35",
        "worst_current_scenario",
        "worst_current_spectrum",
        "worst_current_epsilon",
        "worst_current_phase",
        "worst_future_scenario",
        "worst_future_spectrum",
        "worst_future_epsilon",
        "worst_future_phase",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(channel_targets: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    lines = [
        "# Knu target map",
        "",
        "No web lookup was used.",
        "",
        "The input is the combined e-pi projection profile from the previous replay.",
        "",
        "## Preferred conservative channel targets",
        "",
        "| channel | class | current margin | future margin | future amp suppression | future S_T max |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in channel_targets:
        lines.append(
            "| `{channel}` | `{channel_class}` | {worst_current_margin:.6e} | "
            "{worst_future_margin_1e35:.6e} | {future_amplitude_suppression_needed:.6e} | "
            "{S_T_safe_max_future_1e35:.6e} |".format(**row)
        )

    verdict = summary["verdict"]
    lines += [
        "",
        "## Verdict",
        "",
        (
            f"Current bottleneck: `{verdict['current_bottleneck_channel']}` "
            f"with margin {verdict['current_margin']:.6e}."
        ),
        (
            f"Future 1e35 bottleneck: `{verdict['future_bottleneck_channel']}` "
            f"with margin {verdict['future_margin']:.6e}."
        ),
        (
            f"Future-safe S_T target: {verdict['future_ST_target_required']:.6e}; "
            f"equivalent amplitude suppression: "
            f"{verdict['future_amplitude_suppression_needed']:.6e}."
        ),
        "",
        summary["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = read_rows()
    channel_targets, summary = build_summary(rows)
    write_csv(OUT / "channel_targets.csv", channel_targets)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(channel_targets, summary)
    print("Knu target map")
    print(
        "  current bottleneck: {channel}, margin={margin:.6e}".format(
            channel=summary["verdict"]["current_bottleneck_channel"],
            margin=summary["verdict"]["current_margin"],
        )
    )
    print(
        "  future S_T target: {target:.6e}, amplitude suppression={supp:.6e}".format(
            target=summary["verdict"]["future_ST_target_required"],
            supp=summary["verdict"]["future_amplitude_suppression_needed"],
        )
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
