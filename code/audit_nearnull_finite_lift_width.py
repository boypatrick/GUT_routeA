#!/usr/bin/env python3
"""Audit finite lifts of the triplet near-null in the calibrated width formula.

This script is deliberately conservative.  It does not claim to be the final
microscopic colored-triplet sector.  It asks a sharper question:

  * if the rank-deficient near-null is replaced by the finite-rank cards from
    construct_triplet_rank_lift.py,
  * and the same calibrated channel-wise dressing transfer is used,

which cards remain proton-safe, which remain below the reduced Planck scale for
M_lightest = 1e16 GeV, and how large a naive one-unit beta threshold leakage
would be if the triplet mass spread were not protected by complete multiplets?
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "nearnull_finite_lift_width"

RANK_LIFT = ROOT / "output" / "triplet_rank_lift" / "triplet_rank_lift_summary.json"
PHYSICAL = ROOT / "output" / "physical_d5_wilson_replay" / "physical_d5_wilson_replay.csv"
INTERFERENCE = ROOT / "output" / "interference_tightened_d5_replay" / "channel_replay.csv"
WIDTH = ROOT / "output" / "full_knu_width" / "summary.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"

PREFERRED_FILTER = "omegaR_0.1_kappa_100"
PREFERRED_SCENARIO = "democratic_all_stress"
PREFERRED_SPECTRUM = "near_unsafe_20TeV"
PREFERRED_CASE = "max_width"
S_T = 7.5e-6
FUTURE_STRESS_YR = 1.0e35
M_LIGHTEST_TRIPLET_GEV = 1.0e16
REDUCED_PLANCK_GEV = 2.435e18

CHANNELS = [
    "LLLL_upupdown_Knu",
    "LLLL_downdownup_Knu",
    "LLLL_upupdown_K0mu",
    "RRRR_uusd_anycharged",
]


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


def physical_rows_by_channel() -> dict[str, dict[str, Any]]:
    rows = read_csv(PHYSICAL)
    out: dict[str, dict[str, Any]] = {}
    for channel in CHANNELS:
        matches = [
            row
            for row in rows
            if row["label"] == "full_bipartite_mixing_identity"
            and row["channel"] == channel
            and row["normalization_case"] == PREFERRED_CASE
        ]
        if len(matches) != 1:
            raise RuntimeError(f"expected one identity Wilson row for {channel}, got {len(matches)}")
        out[channel] = matches[0]
    return out


def dressed_identity_rows_by_channel() -> dict[str, dict[str, Any]]:
    rows = read_csv(INTERFERENCE)
    out: dict[str, dict[str, Any]] = {}
    for channel in CHANNELS:
        matches = [
            row
            for row in rows
            if row["filter_label"] == PREFERRED_FILTER
            and row["scenario"] == PREFERRED_SCENARIO
            and row["spectrum_name"] == PREFERRED_SPECTRUM
            and row["channel"] == channel
            and row["normalization_case"] == PREFERRED_CASE
        ]
        if not matches:
            raise RuntimeError(f"expected dressed identity rows for {channel}, got 0")
        out[channel] = min(matches, key=lambda row: f(row, "future_margin_1e35_worst_phase"))
    return out


def row_mass_numbers(row: dict[str, Any]) -> tuple[float, float, bool]:
    mass = row["mass_summary"]
    ratio = float(mass["max_finite_mass_ratio"])
    max_mass = float(mass["max_mass_if_lightest_1e16_GeV"])
    planck_safe = max_mass <= REDUCED_PLANCK_GEV
    return ratio, max_mass, planck_safe


def finite_lift_channel_width(
    row: dict[str, Any],
    channel: str,
    identity_raw: dict[str, Any],
    identity_dressed: dict[str, Any],
    k_dyn: float,
    chiral: dict[str, Any],
) -> dict[str, Any]:
    raw_amp = float(row["channels"][channel]["amplitude"])
    ident_raw_amp = f(identity_raw, "amplitude")
    raw_ratio = raw_amp / ident_raw_amp
    dressed_identity_amp = f(identity_dressed, "max_phase_amplitude")
    dressed_amp = dressed_identity_amp * raw_ratio
    k_had = width_prefactor(channel, chiral)
    gamma = k_dyn * k_had * (S_T * dressed_amp) ** 2
    tau_years = math.inf if gamma <= 0.0 else 1.0 / gamma
    margin_current = tau_years / current_bound(channel)
    margin_future = tau_years / FUTURE_STRESS_YR
    return {
        "label": row["label"],
        "kind": row["kind"],
        "epsilon": float(row["epsilon"]),
        "channel": channel,
        "operator": identity_dressed["operator"],
        "raw_amplitude": raw_amp,
        "identity_raw_amplitude": ident_raw_amp,
        "raw_ratio_over_identity": raw_ratio,
        "identity_dressed_amplitude": dressed_identity_amp,
        "finite_lift_dressed_amplitude": dressed_amp,
        "K_had": k_had,
        "K_dyn_year_inv": k_dyn,
        "S_T": S_T,
        "gamma_year_inv": gamma,
        "tau_years": tau_years,
        "current_bound_years": current_bound(channel),
        "margin_current": margin_current,
        "margin_future_1e35": margin_future,
        "future_amplitude_headroom": math.sqrt(margin_future),
        "S_T_max_future_1e35": S_T * math.sqrt(margin_future),
    }


def summarize_lift_row(row: dict[str, Any], channel_rows: list[dict[str, Any]], r200_tol: float) -> dict[str, Any]:
    max_ratio, max_mass, planck_safe = row_mass_numbers(row)
    worst_future = min(channel_rows, key=lambda item: item["margin_future_1e35"])
    worst_current = min(channel_rows, key=lambda item: item["margin_current"])
    llll_rows = [item for item in channel_rows if item["operator"] == "LLLL"]
    rrrr_rows = [item for item in channel_rows if item["operator"] == "RRRR"]
    # Conservative diagnostic only: a single split beta-vector of unit length.
    # Complete multiplet degeneracy or a locked triplet sector can cancel this.
    threshold_unit_beta_proxy = math.log(max_ratio) / (2.0 * math.pi) if max_ratio > 0 else math.inf
    return {
        "label": row["label"],
        "kind": row["kind"],
        "epsilon": float(row["epsilon"]),
        "max_mass_ratio": max_ratio,
        "max_mass_if_lightest_1e16_GeV": max_mass,
        "reduced_planck_GeV": REDUCED_PLANCK_GEV,
        "planck_safe_for_Mlight_1e16": planck_safe,
        "worst_current_channel": worst_current["channel"],
        "worst_current_margin": worst_current["margin_current"],
        "worst_future_channel": worst_future["channel"],
        "worst_future_margin_1e35": worst_future["margin_future_1e35"],
        "worst_future_ST_max": worst_future["S_T_max_future_1e35"],
        "min_LLLL_future_margin": min(item["margin_future_1e35"] for item in llll_rows),
        "min_RRRR_future_margin": min(item["margin_future_1e35"] for item in rrrr_rows),
        "proton_safe_future_1e35": worst_future["margin_future_1e35"] >= 1.0,
        "threshold_unit_beta_proxy": threshold_unit_beta_proxy,
        "r200_projected_threshold_tolerance": r200_tol,
        "passes_naive_threshold_proxy": threshold_unit_beta_proxy <= r200_tol,
        "channel_rows": channel_rows,
    }


def build() -> dict[str, Any]:
    rank = read_json(RANK_LIFT)
    width = read_json(WIDTH)
    chiral = read_json(CHIRAL)
    ledger = read_json(LEDGER)
    k_dyn = float(width["verdict"]["K_dyn_year_inv"])
    r200_tol = float(ledger["r200_benchmark"]["total_projected_l2"])
    identity_raw = physical_rows_by_channel()
    identity_dressed = dressed_identity_rows_by_channel()

    channel_rows: list[dict[str, Any]] = []
    lift_summaries: list[dict[str, Any]] = []
    for row in rank["rows"]:
        per_channel = [
            finite_lift_channel_width(
                row=row,
                channel=channel,
                identity_raw=identity_raw[channel],
                identity_dressed=identity_dressed[channel],
                k_dyn=k_dyn,
                chiral=chiral,
            )
            for channel in CHANNELS
        ]
        channel_rows.extend(per_channel)
        lift_summaries.append(summarize_lift_row(row, per_channel, r200_tol))

    planck_and_proton = [
        row
        for row in lift_summaries
        if row["planck_safe_for_Mlight_1e16"] and row["proton_safe_future_1e35"]
    ]
    planck_proton_threshold = [
        row for row in planck_and_proton if row["passes_naive_threshold_proxy"]
    ]
    best = max(planck_and_proton, key=lambda row: row["worst_future_margin_1e35"], default=None)
    lowest_mass = min(planck_and_proton, key=lambda row: row["max_mass_ratio"], default=None)
    verdict = {
        "planck_and_proton_safe_count": len(planck_and_proton),
        "planck_proton_and_naive_threshold_safe_count": len(planck_proton_threshold),
        "best_planck_proton_label": None if best is None else best["label"],
        "best_planck_proton_worst_future_margin": None if best is None else best["worst_future_margin_1e35"],
        "best_planck_proton_max_mass_GeV": None if best is None else best["max_mass_if_lightest_1e16_GeV"],
        "lowest_mass_planck_proton_label": None if lowest_mass is None else lowest_mass["label"],
        "lowest_mass_planck_proton_max_mass_ratio": None if lowest_mass is None else lowest_mass["max_mass_ratio"],
        "lowest_mass_planck_proton_worst_future_margin": None if lowest_mass is None else lowest_mass["worst_future_margin_1e35"],
        "r200_projected_threshold_tolerance": r200_tol,
        "interpretation": (
            "Condition-capped finite lifts can keep the calibrated d=5 widths "
            "and the M_lightest=1e16 GeV Planck benchmark safe.  However, no "
            "such row passes the deliberately conservative one-unit beta-vector "
            "threshold proxy.  The next microscopic requirement is therefore a "
            "threshold-degenerate complete-multiplet/locked triplet lift, not "
            "another arbitrary split mass matrix."
        ),
    }
    return {
        "note": "No web lookup used. Finite near-null lift inserted into the calibrated d=5 width bridge.",
        "formula": (
            "A_lift_dress(c) = A_identity_dress(c) |C_lift(c)|/|C_identity(c)|; "
            "Gamma = K_dyn K_had (S_T |A_lift_dress|)^2"
        ),
        "benchmarks": {
            "S_T": S_T,
            "future_stress_years": FUTURE_STRESS_YR,
            "M_lightest_triplet_GeV": M_LIGHTEST_TRIPLET_GEV,
            "reduced_planck_GeV": REDUCED_PLANCK_GEV,
            "planck_safe_mass_ratio": REDUCED_PLANCK_GEV / M_LIGHTEST_TRIPLET_GEV,
            "r200_projected_threshold_tolerance": r200_tol,
            "threshold_proxy": "ln(max_mass_ratio)/(2*pi) for a single uncancelled unit beta-vector",
        },
        "input_files": {
            "rank_lift": str(RANK_LIFT),
            "physical_d5_wilson_replay": str(PHYSICAL),
            "interference_tightened_channel_replay": str(INTERFERENCE),
            "full_knu_width": str(WIDTH),
            "chiral_lattice": str(CHIRAL),
            "conditional_theorem_ledger": str(LEDGER),
        },
        "lift_summaries": lift_summaries,
        "channel_rows": channel_rows,
        "verdict": verdict,
    }


def write_channel_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "kind",
        "epsilon",
        "channel",
        "operator",
        "raw_amplitude",
        "identity_raw_amplitude",
        "raw_ratio_over_identity",
        "identity_dressed_amplitude",
        "finite_lift_dressed_amplitude",
        "K_had",
        "K_dyn_year_inv",
        "S_T",
        "gamma_year_inv",
        "tau_years",
        "current_bound_years",
        "margin_current",
        "margin_future_1e35",
        "future_amplitude_headroom",
        "S_T_max_future_1e35",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "kind",
        "epsilon",
        "max_mass_ratio",
        "max_mass_if_lightest_1e16_GeV",
        "planck_safe_for_Mlight_1e16",
        "worst_future_channel",
        "worst_future_margin_1e35",
        "min_LLLL_future_margin",
        "min_RRRR_future_margin",
        "proton_safe_future_1e35",
        "threshold_unit_beta_proxy",
        "r200_projected_threshold_tolerance",
        "passes_naive_threshold_proxy",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    rows = summary["lift_summaries"]
    lines = [
        "# Finite near-null lift calibrated width audit",
        "",
        "No web lookup was used.",
        "",
        "The channel transfer rule is",
        "",
        "```text",
        "A_lift_dress(c) = A_identity_dress(c) |C_lift(c)|/|C_identity(c)|",
        "Gamma = K_dyn K_had (S_T |A_lift_dress|)^2",
        "```",
        "",
        "The threshold column is only a conservative diagnostic:",
        "",
        "```text",
        "Delta_proxy = ln(max_mass_ratio)/(2*pi)",
        "```",
        "",
        "| lift | Mmax/Mmin | Mmax [GeV] | worst future margin | Planck safe | proxy threshold |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        if row["kind"] not in {"condition_cap", "rank3_limit"}:
            continue
        lines.append(
            "| `{label}` | {max_mass_ratio:.6e} | {max_mass_if_lightest_1e16_GeV:.6e} | "
            "{worst_future_margin_1e35:.6e} | {planck_safe_for_Mlight_1e16} | "
            "{threshold_unit_beta_proxy:.6e} |".format(**row)
        )
    v = summary["verdict"]
    b = summary["benchmarks"]
    lines += [
        "",
        "## Verdict",
        "",
        f"Planck-safe and proton-safe finite lifts: {v['planck_and_proton_safe_count']}.",
        (
            "Planck-safe, proton-safe, and naive-threshold-safe finite lifts: "
            f"{v['planck_proton_and_naive_threshold_safe_count']}."
        ),
        (
            f"Reduced-Planck mass-ratio cap for M_lightest=1e16 GeV: "
            f"{b['planck_safe_mass_ratio']:.6e}."
        ),
        (
            f"Best Planck/proton row: `{v['best_planck_proton_label']}` "
            f"with worst future margin {v['best_planck_proton_worst_future_margin']:.6e} "
            f"and Mmax {v['best_planck_proton_max_mass_GeV']:.6e} GeV."
        ),
        (
            f"Lightest hierarchy Planck/proton row: `{v['lowest_mass_planck_proton_label']}` "
            f"with ratio {v['lowest_mass_planck_proton_max_mass_ratio']:.6e} "
            f"and worst future margin {v['lowest_mass_planck_proton_worst_future_margin']:.6e}."
        ),
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    write_channel_csv(OUT / "channel_widths.csv", summary["channel_rows"])
    write_summary_csv(OUT / "lift_summary.csv", summary["lift_summaries"])
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(summary)
    v = summary["verdict"]
    print("Finite near-null lift calibrated width audit")
    print(f"  Planck/proton safe count: {v['planck_and_proton_safe_count']}")
    print(
        "  Planck/proton/naive-threshold safe count: "
        f"{v['planck_proton_and_naive_threshold_safe_count']}"
    )
    print(f"  best Planck/proton row: {v['best_planck_proton_label']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
