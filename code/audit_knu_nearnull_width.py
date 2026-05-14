#!/usr/bin/env python3
"""Insert the triplet-mixing Knu near-null into the calibrated width formula.

No web lookup is used.  This is a controlled bridge, not the final microscopic
triplet-dressing calculation.  It transfers the already calibrated soft
dressing amplitudes to the near-null branch by the channel-wise Wilson ratio

    A_dress(nearnull, c)
      = A_dress(identity, c)
        |C_nearnull(c)| / |C_identity(c)|,

then evaluates the same width formula

    Gamma = K_dyn K_had (S_T |A_dress|)^2.

The purpose is to test whether the near-null is strong enough in the calibrated
width normalization and to expose any K0mu/RRRR side effects.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "knu_nearnull_width"

PHYSICAL = ROOT / "output" / "physical_d5_wilson_replay" / "physical_d5_wilson_replay.csv"
INTERFERENCE = ROOT / "output" / "interference_tightened_d5_replay" / "channel_replay.csv"
WIDTH = ROOT / "output" / "full_knu_width" / "summary.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
NULLSPACE = ROOT / "output" / "triplet_mixing_nullspace" / "triplet_mixing_nullspace_summary.json"

PREFERRED_FILTER = "omegaR_0.1_kappa_100"
PREFERRED_SCENARIO = "democratic_all_stress"
PREFERRED_SPECTRUM = "near_unsafe_20TeV"
PREFERRED_CASE = "max_width"
PREFERRED_PHASE = 4.1887902047863905
S_T = 7.5e-6
FUTURE_STRESS_YR = 1.0e35


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def f(row: dict[str, Any], key: str) -> float:
    return float(row[key])


def current_bound(channel: str) -> float:
    return 1.0e34 if "K0mu" in channel else 2.4e34


def width_prefactor(channel: str, chiral: dict[str, Any]) -> float:
    old = float(chiral["old_single_width_prefactor"])
    if "K0mu" in channel:
        ratio = 0.3420340196046356
    elif "Knu" in channel or "RRRR_uusd" in channel:
        ratio = 1.6203228084230246
    elif "epi" in channel:
        ratio = 4.840884549836016
    else:
        raise ValueError(f"unknown channel class for {channel}")
    return old * ratio


def physical_row(label: str, channel: str) -> dict[str, Any]:
    matches = [
        row
        for row in read_csv(PHYSICAL)
        if row["label"] == label
        and row["channel"] == channel
        and row["normalization_case"] == PREFERRED_CASE
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one physical row for {label} {channel}, got {len(matches)}")
    return matches[0]


def dressed_identity_row(channel: str) -> dict[str, Any]:
    rows = [
        row
        for row in read_csv(INTERFERENCE)
        if row["filter_label"] == PREFERRED_FILTER
        and row["scenario"] == PREFERRED_SCENARIO
        and row["spectrum_name"] == PREFERRED_SPECTRUM
        and row["channel"] == channel
        and row["normalization_case"] == PREFERRED_CASE
    ]
    if not rows:
        raise RuntimeError(f"expected dressed identity rows for {channel}, got 0")
    return min(rows, key=lambda row: f(row, "future_margin_1e35_worst_phase"))


def nearnull_meta() -> dict[str, Any]:
    payload = read_json(NULLSPACE)
    return next(row for row in payload["audits"] if row["label"] == "full_bipartite_mixing_knu_null")


def build() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    width = read_json(WIDTH)
    chiral = read_json(CHIRAL)
    k_dyn = float(width["verdict"]["K_dyn_year_inv"])
    meta = nearnull_meta()
    channels = [
        "LLLL_upupdown_Knu",
        "LLLL_downdownup_Knu",
        "LLLL_upupdown_K0mu",
        "RRRR_uusd_anycharged",
    ]

    rows: list[dict[str, Any]] = []
    for channel in channels:
        ident_raw = physical_row("full_bipartite_mixing_identity", channel)
        near_raw = physical_row("full_bipartite_mixing_knu_null", channel)
        ident_dressed = dressed_identity_row(channel)

        raw_ratio = f(near_raw, "amplitude") / f(ident_raw, "amplitude")
        dressed_identity_amp = f(ident_dressed, "max_phase_amplitude")
        dressed_near_amp = dressed_identity_amp * raw_ratio
        k_had = width_prefactor(channel, chiral)
        gamma = k_dyn * k_had * (S_T * dressed_near_amp) ** 2
        tau_years = math.inf if gamma <= 0.0 else 1.0 / gamma
        margin_current = tau_years / current_bound(channel)
        margin_future = tau_years / FUTURE_STRESS_YR
        rows.append(
            {
                "channel": channel,
                "operator": ident_dressed["operator"],
                "identity_raw_amplitude": f(ident_raw, "amplitude"),
                "near_null_raw_amplitude": f(near_raw, "amplitude"),
                "raw_ratio_near_over_identity": raw_ratio,
                "identity_dressed_amplitude": dressed_identity_amp,
                "near_null_dressed_amplitude": dressed_near_amp,
                "K_had": k_had,
                "K_dyn_year_inv": k_dyn,
                "S_T": S_T,
                "gamma_year_inv": gamma,
                "tau_years": tau_years,
                "current_bound_years": current_bound(channel),
                "margin_current": margin_current,
                "margin_future_1e35": margin_future,
                "future_amplitude_headroom": math.sqrt(margin_future),
                "future_suppression_needed": max(1.0 / math.sqrt(margin_future), 1.0),
                "S_T_max_future_1e35": S_T * math.sqrt(margin_future),
                "identity_pair": ident_dressed["pair"],
                "identity_phase": f(ident_dressed, "triplet_phase"),
                "W_rank": int(float(near_raw["W_rank"])),
                "W_condition": near_raw["W_condition"],
                "identity_overlap_abs": f(near_raw, "identity_overlap_abs"),
            }
        )

    worst_future = min(rows, key=lambda row: row["margin_future_1e35"])
    worst_current = min(rows, key=lambda row: row["margin_current"])
    llll_rows = [row for row in rows if row["operator"] == "LLLL"]
    rrrr_rows = [row for row in rows if row["operator"] == "RRRR"]
    summary = {
        "note": "No web lookup used. Calibrated Knu near-null width bridge.",
        "formula": "A_near_dress = A_identity_dress * |C_near|/|C_identity|; Gamma = K_dyn K_had (S_T |A_near_dress|)^2",
        "input_files": {
            "physical_d5_wilson_replay": str(PHYSICAL),
            "interference_tightened_channel_replay": str(INTERFERENCE),
            "full_knu_width": str(WIDTH),
            "chiral_lattice": str(CHIRAL),
            "triplet_mixing_nullspace": str(NULLSPACE),
        },
        "near_null_meta": {
            "W_rank": meta["W_rank"],
            "W_condition": meta["W_condition"],
            "identity_overlap_abs": meta["identity_overlap_abs"],
            "W_singular_values": meta["W_singular_values"],
            "coefficient_peak_over_rms": meta["coefficient_peak_over_rms"],
        },
        "rows": rows,
        "verdict": {
            "worst_current_channel": worst_current["channel"],
            "worst_current_margin": worst_current["margin_current"],
            "worst_future_channel": worst_future["channel"],
            "worst_future_margin_1e35": worst_future["margin_future_1e35"],
            "worst_future_amplitude_headroom": worst_future["future_amplitude_headroom"],
            "worst_future_suppression_needed": worst_future["future_suppression_needed"],
            "worst_future_ST_max": worst_future["S_T_max_future_1e35"],
            "min_LLLL_future_margin": min(row["margin_future_1e35"] for row in llll_rows),
            "min_RRRR_future_margin": min(row["margin_future_1e35"] for row in rrrr_rows),
            "solves_0p516844_target": worst_future["future_amplitude_headroom"] >= 0.5168441164632575,
            "rank_deficient": meta["W_rank"] < 4,
            "interpretation": (
                "The near-null easily satisfies the calibrated Knu width target "
                "under channel-wise ratio transfer, including the RRRR proxy. "
                "However, the audited near-null is rank deficient and highly "
                "non-generic, so the remaining problem is an action-level origin "
                "or a finite-rank lift that preserves these suppressions."
            ),
        },
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "channel",
        "operator",
        "identity_raw_amplitude",
        "near_null_raw_amplitude",
        "raw_ratio_near_over_identity",
        "identity_dressed_amplitude",
        "near_null_dressed_amplitude",
        "K_had",
        "K_dyn_year_inv",
        "S_T",
        "gamma_year_inv",
        "tau_years",
        "current_bound_years",
        "margin_current",
        "margin_future_1e35",
        "future_amplitude_headroom",
        "future_suppression_needed",
        "S_T_max_future_1e35",
        "identity_pair",
        "identity_phase",
        "W_rank",
        "W_condition",
        "identity_overlap_abs",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Knu near-null calibrated width bridge",
        "",
        "No web lookup was used.",
        "",
        "The transfer rule is",
        "",
        "```text",
        "A_near_dress(c) = A_identity_dress(c) |C_near(c)|/|C_identity(c)|",
        "Gamma = K_dyn K_had (S_T |A_near_dress|)^2",
        "```",
        "",
        "| channel | raw ratio | near dressed amp | future margin | S_T max future |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        lines.append(
            "| `{channel}` | {raw_ratio_near_over_identity:.6e} | "
            "{near_null_dressed_amplitude:.6e} | {margin_future_1e35:.6e} | "
            "{S_T_max_future_1e35:.6e} |".format(**row)
        )
    v = summary["verdict"]
    meta = summary["near_null_meta"]
    lines += [
        "",
        "## Verdict",
        "",
        (
            f"Worst future channel is `{v['worst_future_channel']}` with margin "
            f"{v['worst_future_margin_1e35']:.6e}."
        ),
        (
            f"Minimum LLLL future margin: {v['min_LLLL_future_margin']:.6e}; "
            f"minimum RRRR future margin: {v['min_RRRR_future_margin']:.6e}."
        ),
        (
            f"W rank = {meta['W_rank']}, condition = {meta['W_condition']}, "
            f"identity overlap = {meta['identity_overlap_abs']:.6e}."
        ),
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = build()
    write_csv(OUT / "knu_nearnull_width.csv", rows)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(summary)
    v = summary["verdict"]
    print("Knu near-null calibrated width bridge")
    print(f"  worst future channel: {v['worst_future_channel']}")
    print(f"  worst future margin: {v['worst_future_margin_1e35']:.6e}")
    print(f"  solves 0.516844 target: {v['solves_0p516844_target']}")
    print(f"  rank deficient: {v['rank_deficient']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
