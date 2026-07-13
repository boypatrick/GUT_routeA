#!/usr/bin/env python3
"""Compare the Knu target map with the older mass-basis Wilson tensor audit.

No web lookup is used.  The latest Knu target map says that, in the
conservative soft-spectrum/chiral replay, future 1e35 yr safety needs

    S_T <= 3.876330873474431e-6

or an equivalent Knu amplitude suppression.  The older mass-basis Wilson
tensor table used a common 100 TeV dressing normalization and therefore has a
different scale.  This script measures that normalization gap rather than
pretending the two audits are already the same calculation.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
WILSON = ROOT / "output" / "dimension5_wilson_tensors" / "dimension5_wilson_tensors.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "knu_wilson_normalization_gap"

S_T_CURRENT = 7.5e-6


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_knu(row: dict[str, Any]) -> bool:
    return "Knu" in row["channel_proxy"] or "RRRR_uusd" in row["channel_proxy"]


def build_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    wilson = read_json(WILSON)
    target = read_json(KNU_TARGET)
    target_st = float(target["verdict"]["future_ST_target_required"])
    target_margin = float(target["verdict"]["future_margin"])
    target_supp = float(target["verdict"]["future_amplitude_suppression_needed"])

    rows: list[dict[str, Any]] = []
    for row in wilson["rows"]:
        if not is_knu(row):
            continue
        st_req = float(row["S_T_required_tau_1e35"])
        current_headroom = st_req / S_T_CURRENT
        target_gap = st_req / target_st
        rows.append(
            {
                "hypothesis": row["hypothesis"],
                "operator": row["operator"],
                "channel_proxy": row["channel_proxy"],
                "selected_index": row["selected_index"],
                "amplitude": float(row["amplitude"]),
                "S_T_required_tau_1e35_common_dressing": st_req,
                "margin_at_S_T_7p5e_minus_6_common_dressing": current_headroom**2,
                "amplitude_headroom_at_S_T_7p5e_minus_6_common_dressing": current_headroom,
                "ratio_common_dressing_STmax_to_conservative_target": target_gap,
                "conservative_target_S_T_future": target_st,
                "conservative_target_future_margin": target_margin,
                "conservative_target_amplitude_suppression_needed": target_supp,
                "note": row["note"],
            }
        )

    rows.sort(key=lambda item: item["S_T_required_tau_1e35_common_dressing"])
    worst = rows[0]
    summary = {
        "note": "No web lookup used. Knu Wilson normalization-gap audit.",
        "input_files": {
            "dimension5_wilson_tensors": str(WILSON),
            "knu_target_map": str(KNU_TARGET),
        },
        "S_T_current": S_T_CURRENT,
        "conservative_target": {
            "future_ST_target_required": target_st,
            "future_margin_at_S_T_current": target_margin,
            "future_amplitude_suppression_needed": target_supp,
        },
        "rows": rows,
        "verdict": {
            "most_dangerous_common_dressing_channel": worst["channel_proxy"],
            "most_dangerous_common_dressing_hypothesis": worst["hypothesis"],
            "S_T_required_tau_1e35_common_dressing": worst[
                "S_T_required_tau_1e35_common_dressing"
            ],
            "margin_at_S_T_7p5e_minus_6_common_dressing": worst[
                "margin_at_S_T_7p5e_minus_6_common_dressing"
            ],
            "common_dressing_amplitude_headroom": worst[
                "amplitude_headroom_at_S_T_7p5e_minus_6_common_dressing"
            ],
            "normalization_gap_to_conservative_target": worst[
                "ratio_common_dressing_STmax_to_conservative_target"
            ],
            "interpretation": (
                "The older mass-basis Wilson tensor table has large Knu headroom "
                "under its common 100 TeV dressing normalization: even the most "
                "dangerous Knu row allows S_T about 24.6 times larger than the "
                "current 7.5e-6 filter.  However, this is about 47.6 times looser "
                "than the conservative soft-spectrum target map.  Therefore the "
                "remaining d=5 problem is a normalization/dressing/hadronic bridge, "
                "not merely a bare flavor-rotation problem."
            ),
        },
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "hypothesis",
        "operator",
        "channel_proxy",
        "selected_index",
        "amplitude",
        "S_T_required_tau_1e35_common_dressing",
        "margin_at_S_T_7p5e_minus_6_common_dressing",
        "amplitude_headroom_at_S_T_7p5e_minus_6_common_dressing",
        "ratio_common_dressing_STmax_to_conservative_target",
        "conservative_target_S_T_future",
        "conservative_target_future_margin",
        "conservative_target_amplitude_suppression_needed",
        "note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Knu Wilson normalization-gap audit",
        "",
        "No web lookup was used.",
        "",
        "This compares the older mass-basis Wilson tensor table, which uses a",
        "common 100 TeV dressing normalization, with the newer conservative",
        "soft-spectrum Knu target map.",
        "",
        "| hypothesis | operator | channel | S_T max at 1e35 | headroom at 7.5e-6 | gap to conservative target |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        lines.append(
            "| `{hypothesis}` | `{operator}` | `{channel_proxy}` | "
            "{S_T_required_tau_1e35_common_dressing:.6e} | "
            "{amplitude_headroom_at_S_T_7p5e_minus_6_common_dressing:.6e} | "
            "{ratio_common_dressing_STmax_to_conservative_target:.6e} |".format(**row)
        )

    verdict = summary["verdict"]
    lines += [
        "",
        "## Verdict",
        "",
        (
            f"Most dangerous common-dressing Knu row: "
            f"`{verdict['most_dangerous_common_dressing_channel']}` under "
            f"`{verdict['most_dangerous_common_dressing_hypothesis']}`."
        ),
        (
            f"It allows S_T = "
            f"{verdict['S_T_required_tau_1e35_common_dressing']:.6e} at 1e35 yr, "
            f"which is a lifetime margin "
            f"{verdict['margin_at_S_T_7p5e_minus_6_common_dressing']:.6e} "
            f"at S_T=7.5e-6."
        ),
        (
            f"The gap to the conservative target map is "
            f"{verdict['normalization_gap_to_conservative_target']:.6e} in amplitude."
        ),
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = build_rows()
    write_csv(OUT / "knu_wilson_normalization_gap.csv", rows)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(summary)
    verdict = summary["verdict"]
    print("Knu Wilson normalization-gap audit")
    print(
        "  most dangerous: {hyp} {ch}".format(
            hyp=verdict["most_dangerous_common_dressing_hypothesis"],
            ch=verdict["most_dangerous_common_dressing_channel"],
        )
    )
    print(
        "  S_T max common dressing: {st:.6e}".format(
            st=verdict["S_T_required_tau_1e35_common_dressing"]
        )
    )
    print(
        "  gap to conservative target: {gap:.6e}".format(
            gap=verdict["normalization_gap_to_conservative_target"]
        )
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
