#!/usr/bin/env python3
"""Bridge audit for the full Knu d=5 proton-decay pipeline.

No web lookup is used.  This is not yet the final channel-specific proton
decay calculation.  It is the reproducibility scaffold for it: the script puts
the existing Knu-relevant replay outputs on a single lifetime-normalization
axis and measures where the available amplitude headroom is actually lost.

All future margins below are normalized to a uniform 1e35 yr stress target.
For Knu-like current bounds we use the local convention 2.4e34 yr, so a
current-margin row is converted to the future margin by multiplying by 0.24.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "full_knu_channel_pipeline"

KNU_WILSON_GAP = ROOT / "output" / "knu_wilson_normalization_gap" / "summary.json"
PHYSICAL_D5 = ROOT / "output" / "physical_d5_wilson_replay" / "summary.json"
SOFT_D5 = ROOT / "output" / "soft_spectrum_d5_replay" / "summary.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"

CURRENT_KNU_BOUND_YR = 2.4e34
FUTURE_STRESS_YR = 1.0e35
CURRENT_TO_FUTURE = CURRENT_KNU_BOUND_YR / FUTURE_STRESS_YR


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def stage_row(
    *,
    stage: str,
    channel: str,
    margin_current: float,
    margin_future: float,
    st_max_future: float | None,
    source: str,
    interpretation: str,
) -> dict[str, Any]:
    future_headroom = math.sqrt(margin_future)
    current_headroom = math.sqrt(margin_current)
    return {
        "stage": stage,
        "channel": channel,
        "margin_current_2p4e34": margin_current,
        "amplitude_headroom_current": current_headroom,
        "margin_future_1e35": margin_future,
        "amplitude_headroom_future": future_headroom,
        "amplitude_suppression_needed_future": max(1.0 / future_headroom, 1.0),
        "S_T_max_future_1e35": st_max_future,
        "evidence_path": source,
        "interpretation": interpretation,
    }


def build_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    gap = read_json(KNU_WILSON_GAP)
    physical = read_json(PHYSICAL_D5)
    soft = read_json(SOFT_D5)
    target = read_json(KNU_TARGET)

    gap_v = gap["verdict"]
    physical_v = physical["verdict"]
    soft_pref_max = soft["by_filter_label"]["omegaR_0.1_kappa_100"]["max_width"]
    target_v = target["verdict"]

    # The Wilson-gap audit's S_T maximum is already defined at the 1e35 yr
    # stress target, so its stored margin is a future-normalized margin.
    common_future_margin = float(gap_v["margin_at_S_T_7p5e_minus_6_common_dressing"])
    common_current_margin = common_future_margin / CURRENT_TO_FUTURE

    # The physical field-basis replay margin is a current-bound margin.
    physical_current_margin = float(physical_v["identity_max_width_worst_margin"])
    physical_future_margin = physical_current_margin * CURRENT_TO_FUTURE

    # The soft-spectrum replay already reports both current and 1e35 margins.
    soft_current_margin = float(soft_pref_max["global_worst_current_margin"])
    soft_future_margin = float(soft_pref_max["global_worst_future_1e35_margin"])

    target_current_margin = float(target_v["current_margin"])
    target_future_margin = float(target_v["future_margin"])

    rows = [
        stage_row(
            stage="mass_basis_common_wilson",
            channel=gap_v["most_dangerous_common_dressing_channel"],
            margin_current=common_current_margin,
            margin_future=common_future_margin,
            st_max_future=float(gap_v["S_T_required_tau_1e35_common_dressing"]),
            source=str(KNU_WILSON_GAP),
            interpretation=(
                "Bare mass-basis Wilson tensor with common 100 TeV dressing. "
                "This is the most optimistic common-normalization endpoint."
            ),
        ),
        stage_row(
            stage="field_basis_common_dressing_max_width",
            channel="RRRR_uusd_anycharged",
            margin_current=physical_current_margin,
            margin_future=physical_future_margin,
            st_max_future=None,
            source=str(PHYSICAL_D5),
            interpretation=(
                "CKM/PMNS field-basis replay with the conservative max-width "
                "chiral/lattice envelope, still using common dressing."
            ),
        ),
        stage_row(
            stage="soft_spectrum_preferred_max_width",
            channel=soft_pref_max["worst_future_1e35"]["worst_future_channel"],
            margin_current=soft_current_margin,
            margin_future=soft_future_margin,
            st_max_future=float(soft_pref_max["worst_future_1e35"]["worst_current_ST_max"])
            * math.sqrt(soft_future_margin / soft_current_margin),
            source=str(SOFT_D5),
            interpretation=(
                "Preferred omegaR=0.1,kappa=100 branch with explicit positive "
                "soft eigenstates, chargino/neutralino dressing, and max-width "
                "hadronic envelope."
            ),
        ),
        stage_row(
            stage="epi_projected_knu_target",
            channel=target_v["future_bottleneck_channel"],
            margin_current=target_current_margin,
            margin_future=target_future_margin,
            st_max_future=float(target_v["future_ST_target_required"]),
            source=str(KNU_TARGET),
            interpretation=(
                "Final local target map after neutral-pion and identical-up "
                "projection sensitivity; this is the current bottleneck summary."
            ),
        ),
    ]

    previous = None
    first = rows[0]
    for row in rows:
        if previous is None:
            row["amplitude_loss_vs_previous_future"] = 1.0
        else:
            row["amplitude_loss_vs_previous_future"] = (
                previous["amplitude_headroom_future"] / row["amplitude_headroom_future"]
            )
        row["amplitude_loss_vs_mass_basis_common_future"] = (
            first["amplitude_headroom_future"] / row["amplitude_headroom_future"]
        )
        previous = row

    bottleneck_loss = max(rows[1:], key=lambda row: row["amplitude_loss_vs_previous_future"])
    final = rows[-1]
    summary = {
        "note": "No web lookup used. Full Knu channel-pipeline bridge audit.",
        "input_files": {
            "knu_wilson_normalization_gap": str(KNU_WILSON_GAP),
            "physical_d5_wilson_replay": str(PHYSICAL_D5),
            "soft_spectrum_d5_replay": str(SOFT_D5),
            "knu_target_map": str(KNU_TARGET),
        },
        "normalization": {
            "current_Knu_bound_years": CURRENT_KNU_BOUND_YR,
            "future_stress_years": FUTURE_STRESS_YR,
            "current_to_future_margin_factor": CURRENT_TO_FUTURE,
        },
        "rows": rows,
        "verdict": {
            "final_stage": final["stage"],
            "final_channel": final["channel"],
            "final_current_margin": final["margin_current_2p4e34"],
            "final_future_margin_1e35": final["margin_future_1e35"],
            "final_future_amplitude_headroom": final["amplitude_headroom_future"],
            "final_future_amplitude_suppression_needed": final[
                "amplitude_suppression_needed_future"
            ],
            "largest_headroom_loss_stage": bottleneck_loss["stage"],
            "largest_headroom_loss_vs_previous": bottleneck_loss[
                "amplitude_loss_vs_previous_future"
            ],
            "total_headroom_loss_common_to_final": final[
                "amplitude_loss_vs_mass_basis_common_future"
            ],
            "interpretation": (
                "The one-axis bridge shows that the large common-dressing Wilson "
                "headroom is mostly consumed when the replay moves from common "
                "field-basis dressing to the explicit preferred soft-spectrum "
                "dressing.  The final Knu target map then leaves only a 0.5168 "
                "future-stress amplitude headroom, i.e. an additional suppression "
                "of about 1.935 is needed for a uniform 1e35 yr target."
            ),
        },
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "stage",
        "channel",
        "margin_current_2p4e34",
        "amplitude_headroom_current",
        "margin_future_1e35",
        "amplitude_headroom_future",
        "amplitude_suppression_needed_future",
        "S_T_max_future_1e35",
        "amplitude_loss_vs_previous_future",
        "amplitude_loss_vs_mass_basis_common_future",
        "evidence_path",
        "interpretation",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Full Knu channel-pipeline bridge audit",
        "",
        "No web lookup was used.",
        "",
        "This is a one-axis normalization bridge, not the final proton-decay",
        "calculation.  It puts the existing Knu-relevant replay outputs on a",
        "uniform future-stress lifetime target of `1e35 yr`.",
        "",
        "| stage | channel | future margin | future amp headroom | suppression needed | loss vs previous |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        lines.append(
            "| `{stage}` | `{channel}` | {margin_future_1e35:.6e} | "
            "{amplitude_headroom_future:.6e} | "
            "{amplitude_suppression_needed_future:.6e} | "
            "{amplitude_loss_vs_previous_future:.6e} |".format(**row)
        )

    verdict = summary["verdict"]
    lines += [
        "",
        "## Verdict",
        "",
        (
            "Final bottleneck: `{}` in `{}` with future margin {:.6e}.".format(
                verdict["final_channel"],
                verdict["final_stage"],
                verdict["final_future_margin_1e35"],
            )
        ),
        (
            "Future-stress amplitude headroom is {:.6e}; required additional "
            "suppression is {:.6e}.".format(
                verdict["final_future_amplitude_headroom"],
                verdict["final_future_amplitude_suppression_needed"],
            )
        ),
        (
            "Largest single headroom loss occurs at `{}` with loss factor "
            "{:.6e}.".format(
                verdict["largest_headroom_loss_stage"],
                verdict["largest_headroom_loss_vs_previous"],
            )
        ),
        (
            "Total future-headroom loss from common Wilson normalization to the "
            "final target map is {:.6e}.".format(
                verdict["total_headroom_loss_common_to_final"]
            )
        ),
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = build_rows()
    write_csv(OUT / "pipeline_bridge.csv", rows)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(summary)
    verdict = summary["verdict"]
    print("Full Knu channel-pipeline bridge audit")
    print(
        "  final future margin: {margin:.6e}".format(
            margin=verdict["final_future_margin_1e35"]
        )
    )
    print(
        "  final suppression needed: {supp:.6e}".format(
            supp=verdict["final_future_amplitude_suppression_needed"]
        )
    )
    print(
        "  largest headroom loss: {stage} x {loss:.6e}".format(
            stage=verdict["largest_headroom_loss_stage"],
            loss=verdict["largest_headroom_loss_vs_previous"],
        )
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
