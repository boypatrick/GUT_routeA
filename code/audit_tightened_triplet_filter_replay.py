#!/usr/bin/env python3
"""Replay the d=5 channel ledger at tightened triplet-filter targets.

No web lookup is used.  The chiral/lattice audit found that the displayed
S_T=1e-5 filter is not robust under the conservative hadronic-width envelope.
This script performs the algebraic replay

    tau(S_T) / tau_bound = margin(S_0) (S_0/S_T)^2,

with S_0=1e-5, and checks whether the preferred filter closes the envelope
when S_T is set to the exact max-width allowance or to a rounded conservative
target.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "chiral_replay.csv"
CHIRAL_SUMMARY = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
OUT = ROOT / "output" / "tightened_triplet_filter"

DISPLAY_ST = 1.0e-5
PREFERRED = "omegaR_0.1_kappa_100"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with CHIRAL.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "triplet_phase",
                "old_width_prefactor",
                "new_width_prefactor",
                "prefactor_ratio_new_over_old",
                "S_T_max",
                "S_T_max_replayed",
                "tau_years_ST_display",
                "tau_years_replayed",
                "margin_at_ST_display",
                "margin_replayed",
            ]:
                row[key] = float(row[key])
            row["passes_replayed"] = row["passes_replayed"] == "True"
            rows.append(row)
    return rows


def trimmed(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "triplet_phase",
        "channel",
        "channel_class",
        "operator",
        "pair",
        "normalization_case",
        "S_T_max_replayed",
        "margin_replayed",
    ]
    return {key: row[key] for key in keys}


def target_rows(chiral_summary: dict[str, Any]) -> list[dict[str, Any]]:
    exact = float(chiral_summary["verdict"]["max_width_global_S_T_max"])
    return [
        {
            "target_label": "display_1e-5",
            "S_T_target": DISPLAY_ST,
            "interpretation": "Legacy display target used before the chiral/lattice envelope.",
        },
        {
            "target_label": "exact_max_width_preferred",
            "S_T_target": exact,
            "interpretation": "Exact preferred-filter max-width envelope allowance.",
        },
        {
            "target_label": "rounded_7p5e-6",
            "S_T_target": 7.5e-6,
            "interpretation": "Rounded conservative target below the exact max-width allowance.",
        },
    ]


def audit() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    chiral_summary = read_json(CHIRAL_SUMMARY)
    rows = read_rows()
    targets = target_rows(chiral_summary)
    replay: list[dict[str, Any]] = []

    for target in targets:
        st = float(target["S_T_target"])
        scale = (DISPLAY_ST / st) ** 2
        for row in rows:
            margin = row["margin_replayed"] * scale
            replay.append(
                {
                    **row,
                    "target_label": target["target_label"],
                    "S_T_target": st,
                    "margin_at_target": margin,
                    "passes_target": st <= row["S_T_max_replayed"],
                }
            )

    summary = summarize(replay, targets)
    return replay, summary


def summarize(replay: list[dict[str, Any]], targets: list[dict[str, Any]]) -> dict[str, Any]:
    by_target: dict[str, Any] = {}
    for target in targets:
        label = target["target_label"]
        target_rows_all = [row for row in replay if row["target_label"] == label]
        by_filter: dict[str, Any] = {}
        for filt in sorted({row["filter_label"] for row in target_rows_all}):
            filt_rows = [row for row in target_rows_all if row["filter_label"] == filt]
            by_case: dict[str, Any] = {}
            for case in ["central", "min_width", "max_width"]:
                case_rows = [row for row in filt_rows if row["normalization_case"] == case]
                safe = [row for row in case_rows if row["passes_target"]]
                unsafe = [row for row in case_rows if not row["passes_target"]]
                by_case[case] = {
                    "channel_rows": len(case_rows),
                    "safe_channel_rows": len(safe),
                    "unsafe_channel_rows": len(unsafe),
                    "worst_margin_at_target": min(row["margin_at_target"] for row in case_rows),
                    "global_S_T_max": min(row["S_T_max_replayed"] for row in case_rows),
                    "worst_safe": trimmed(min(safe, key=lambda row: row["margin_at_target"])) if safe else None,
                    "nearest_unsafe": trimmed(max(unsafe, key=lambda row: row["margin_at_target"])) if unsafe else None,
                }
            by_filter[filt] = by_case
        by_target[label] = {
            "S_T_target": target["S_T_target"],
            "interpretation": target["interpretation"],
            "by_filter_label": by_filter,
        }

    preferred_exact = by_target["exact_max_width_preferred"]["by_filter_label"][PREFERRED]["max_width"]
    preferred_rounded = by_target["rounded_7p5e-6"]["by_filter_label"][PREFERRED]["max_width"]
    all_blocks_rounded = by_target["rounded_7p5e-6"]["by_filter_label"]["all_blocks_kappa_100"]["max_width"]
    verdict = {
        "preferred_exact_passes_max_width": preferred_exact["unsafe_channel_rows"] == 0,
        "preferred_rounded_passes_max_width": preferred_rounded["unsafe_channel_rows"] == 0,
        "all_blocks_rounded_passes_max_width": all_blocks_rounded["unsafe_channel_rows"] == 0,
        "preferred_exact_worst_margin": preferred_exact["worst_margin_at_target"],
        "preferred_rounded_worst_margin": preferred_rounded["worst_margin_at_target"],
        "all_blocks_rounded_worst_margin": all_blocks_rounded["worst_margin_at_target"],
        "interpretation": (
            "The preferred omegaR=0.1,kappa=100 branch closes the conservative "
            "chiral/lattice envelope once the triplet filter is tightened to "
            "7.5e-6.  The all-block fallback still fails the same max-width "
            "envelope, so the result is filter-specific rather than a generic "
            "d=5 proof."
        ),
    }
    return {
        "note": "No web lookup used. Algebraic S_T replay of chiral/lattice d=5 rows.",
        "display_S_T": DISPLAY_ST,
        "targets": targets,
        "by_target": by_target,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target_label",
        "S_T_target",
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "triplet_phase",
        "channel",
        "channel_class",
        "operator",
        "pair",
        "normalization_case",
        "S_T_max_replayed",
        "margin_replayed",
        "margin_at_target",
        "passes_target",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Tightened triplet-filter replay",
        "",
        "No web lookup was used.",
        "",
        "The replay uses",
        "",
        "```text",
        "tau(S_T)/tau_bound = margin(S_0) (S_0/S_T)^2,  S_0=1e-5.",
        "```",
        "",
        "| target | filter | case | global S_T max | worst margin | unsafe rows |",
        "|---|---|---|---:|---:|---:|",
    ]
    for target_label, target_payload in summary["by_target"].items():
        for filt, filt_payload in target_payload["by_filter_label"].items():
            for case in ["central", "min_width", "max_width"]:
                row = filt_payload[case]
                lines.append(
                    f"| `{target_label}` | `{filt}` | `{case}` | "
                    f"{row['global_S_T_max']:.6e} | "
                    f"{row['worst_margin_at_target']:.6e} | "
                    f"{row['unsafe_channel_rows']} |"
                )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = audit()
    write_csv(OUT / "tightened_replay.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    verdict = summary["verdict"]
    print("Tightened triplet-filter replay")
    print(
        "  preferred exact max-width: "
        f"pass={verdict['preferred_exact_passes_max_width']}, "
        f"worst_margin={verdict['preferred_exact_worst_margin']:.6f}"
    )
    print(
        "  preferred S_T=7.5e-6: "
        f"pass={verdict['preferred_rounded_passes_max_width']}, "
        f"worst_margin={verdict['preferred_rounded_worst_margin']:.6f}"
    )
    print(
        "  all-block S_T=7.5e-6: "
        f"pass={verdict['all_blocks_rounded_passes_max_width']}, "
        f"worst_margin={verdict['all_blocks_rounded_worst_margin']:.6f}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
