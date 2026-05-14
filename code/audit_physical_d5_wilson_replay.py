#!/usr/bin/env python3
"""First physical field-basis replay of d=5 Wilson tensors at tightened S_T.

No web lookup is used.  This audit reuses the existing mass-basis Wilson tensor
and triplet-mixing nullspace outputs, then applies the channel-dependent
chiral/lattice width envelope and the tightened S_T=7.5e-6 target.

It is not a full SUSY proton-decay calculation: the dressing is still the
common 100 TeV proxy used by the Wilson-tensor module.  The purpose is to
separate the field-basis flavor-rotation question from the later soft-spectrum
dressing question.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
NULLSPACE = ROOT / "output" / "triplet_mixing_nullspace" / "triplet_mixing_nullspace_summary.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
TIGHTENED = ROOT / "output" / "tightened_triplet_filter" / "summary.json"
DOWNSTREAM = ROOT / "output" / "tightened_downstream_proton" / "summary.json"
OUT = ROOT / "output" / "physical_d5_wilson_replay"


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


def chiral_widths(chiral: dict[str, Any]) -> dict[str, dict[str, float]]:
    old = float(chiral["old_single_width_prefactor"])
    central = {key: float(value) for key, value in chiral["central_width_prefactors"].items()}
    # The max/min widths follow from the local parameter box.  They are stable
    # constants already audited in output/chiral_lattice_d5/chiral_replay.csv;
    # keep them explicit here so this replay remains independent of row order.
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


def replay() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    nullspace = read_json(NULLSPACE)
    chiral = read_json(CHIRAL)
    tightened = read_json(TIGHTENED)
    downstream = read_json(DOWNSTREAM)
    widths = chiral_widths(chiral)
    st_target = float(next(t for t in tightened["targets"] if t["target_label"] == "rounded_7p5e-6")["S_T_target"])

    rows: list[dict[str, Any]] = []
    for audit in nullspace["audits"]:
        label = audit["label"]
        for channel, payload in audit["channels"].items():
            cls = channel_class(channel)
            old_st_max = float(payload["S_T_required_tau_2p4e34"])
            for case in ["central", "max_width", "min_width"]:
                ratio = widths[cls][case] / widths[cls]["old"]
                st_max = old_st_max / math.sqrt(ratio)
                margin = (st_max / st_target) ** 2
                rows.append(
                    {
                        "label": label,
                        "mode": audit["mode_info"]["mode"],
                        "channel": channel,
                        "channel_class": cls,
                        "normalization_case": case,
                        "amplitude": float(payload["amplitude"]),
                        "S_T_max_old_universal": old_st_max,
                        "width_ratio": ratio,
                        "S_T_max_replayed": st_max,
                        "S_T_target": st_target,
                        "margin_at_target": margin,
                        "passes_target": st_target <= st_max,
                        "identity_overlap_abs": float(audit["identity_overlap_abs"]),
                        "W_rank": int(audit["W_rank"]),
                        "W_condition": audit["W_condition"],
                    }
                )

    summary = summarize(rows, nullspace, st_target, downstream)
    return rows, summary


def summarize(rows: list[dict[str, Any]], nullspace: dict[str, Any], st_target: float, downstream: dict[str, Any]) -> dict[str, Any]:
    by_label: dict[str, Any] = {}
    for label in sorted({row["label"] for row in rows}):
        label_rows = [row for row in rows if row["label"] == label]
        by_case: dict[str, Any] = {}
        for case in ["central", "max_width", "min_width"]:
            case_rows = [row for row in label_rows if row["normalization_case"] == case]
            unsafe = [row for row in case_rows if not row["passes_target"]]
            worst = min(case_rows, key=lambda row: row["margin_at_target"])
            by_case[case] = {
                "channel_rows": len(case_rows),
                "unsafe_channel_rows": len(unsafe),
                "worst_margin_at_target": worst["margin_at_target"],
                "worst_row": trim(worst),
                "max_allowed_extra_amplitude_factor": math.sqrt(worst["margin_at_target"]),
            }
        by_label[label] = by_case

    identity = by_label["full_bipartite_mixing_identity"]["max_width"]
    near_null = by_label["full_bipartite_mixing_knu_null"]["max_width"]
    verdict = {
        "S_T_target": st_target,
        "identity_max_width_worst_margin": identity["worst_margin_at_target"],
        "identity_max_width_extra_amplitude_factor": identity["max_allowed_extra_amplitude_factor"],
        "near_null_max_width_worst_margin": near_null["worst_margin_at_target"],
        "near_null_extra_amplitude_factor": near_null["max_allowed_extra_amplitude_factor"],
        "downstream_R200_margin_2p4e34": downstream["verdict"]["R200_current_margin"],
        "interpretation": (
            "At the common 100 TeV dressing normalization, explicit field-basis CKM/PMNS "
            "Wilson tensors pass the tightened S_T=7.5e-6 target even under the conservative "
            "chiral/lattice width envelope.  The identity triplet-mixing branch still has "
            "O(10) amplitude headroom, while the Knu near-null branch has much larger flavor "
            "headroom but is rank deficient.  Therefore the next risk is not the bare flavor "
            "rotation; it is the full soft-spectrum dressing and UV origin of the required "
            "triplet filter."
        ),
    }
    return {
        "note": "No web lookup used. First physical field-basis d=5 Wilson replay at tightened S_T.",
        "input_files": {
            "triplet_mixing_nullspace": str(NULLSPACE),
            "chiral_lattice_d5": str(CHIRAL),
            "tightened_triplet_filter": str(TIGHTENED),
            "tightened_downstream_proton": str(DOWNSTREAM),
        },
        "nullspace_verdict": nullspace["verdict"],
        "by_label": by_label,
        "verdict": verdict,
    }


def trim(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "label",
        "channel",
        "channel_class",
        "normalization_case",
        "amplitude",
        "S_T_max_replayed",
        "S_T_target",
        "margin_at_target",
        "identity_overlap_abs",
        "W_rank",
        "W_condition",
    ]
    return {key: row[key] for key in keys}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "mode",
        "channel",
        "channel_class",
        "normalization_case",
        "amplitude",
        "S_T_max_old_universal",
        "width_ratio",
        "S_T_max_replayed",
        "S_T_target",
        "margin_at_target",
        "passes_target",
        "identity_overlap_abs",
        "W_rank",
        "W_condition",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Physical d=5 Wilson replay",
        "",
        "No web lookup was used.",
        "",
        "This first field-basis replay combines the triplet-mixing Wilson tensors,",
        "the channel-dependent chiral/lattice width envelope, and the tightened",
        "`S_T=7.5e-6` target.  Dressing is still the common 100 TeV proxy.",
        "",
        "| label | case | worst margin | extra amp factor | unsafe rows | worst channel |",
        "|---|---|---:|---:|---:|---|",
    ]
    for label, payload in summary["by_label"].items():
        for case in ["central", "max_width", "min_width"]:
            row = payload[case]
            worst = row["worst_row"]
            lines.append(
                f"| `{label}` | `{case}` | {row['worst_margin_at_target']:.6e} | "
                f"{row['max_allowed_extra_amplitude_factor']:.6e} | "
                f"{row['unsafe_channel_rows']} | `{worst['channel']}` |"
            )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = replay()
    write_csv(OUT / "physical_d5_wilson_replay.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    verdict = summary["verdict"]
    print("Physical d=5 Wilson replay")
    print(f"  S_T target: {verdict['S_T_target']:.6e}")
    print(f"  identity max-width margin: {verdict['identity_max_width_worst_margin']:.6e}")
    print(f"  identity extra amplitude factor: {verdict['identity_max_width_extra_amplitude_factor']:.6e}")
    print(f"  near-null max-width margin: {verdict['near_null_max_width_worst_margin']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
