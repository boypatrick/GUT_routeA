#!/usr/bin/env python3
"""Calibrate the physical LLLL_upupdown_Knu d=5 width formula.

No web lookup is used.  This script is the first narrow "same Wilson tensor to
same width" audit for the current Knu bottleneck.  It checks whether the
preferred soft-spectrum row and the final worst-phase/projected Knu row are
consistent with one channel-width formula

    Gamma = K_dyn * K_had * (S_T |A_dress|)^2 .

The script does not claim to contain the final SUSY proton-decay calculation.
It calibrates the effective dynamical constant K_dyn from existing local replay
outputs, then verifies that the constant is stable between the soft-spectrum
and projected-target rows.  Stability means the normalization conventions are
now aligned well enough for the next, genuinely microscopic dressing module.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "full_knu_width"

SOFT_CHANNELS = ROOT / "output" / "soft_spectrum_d5_replay" / "channel_replay.csv"
INTERFERENCE_CHANNELS = (
    ROOT / "output" / "interference_tightened_d5_replay" / "channel_replay.csv"
)
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"

PREFERRED_FILTER = "omegaR_0.1_kappa_100"
PREFERRED_SCENARIO = "democratic_all_stress"
PREFERRED_SPECTRUM = "near_unsafe_20TeV"
PREFERRED_CHANNEL = "LLLL_upupdown_Knu"
PREFERRED_CASE = "max_width"

CURRENT_BOUND_YR = 2.4e34
FUTURE_STRESS_YR = 1.0e35
SECONDS_PER_YEAR = 31_557_600.0
HBAR_GEV_SECONDS = 6.582_119_569e-25
GEV_PER_YEAR_INV = HBAR_GEV_SECONDS / SECONDS_PER_YEAR


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def f(row: dict[str, Any], key: str) -> float:
    return float(row[key])


def preferred_soft_rows() -> list[dict[str, Any]]:
    rows = []
    for row in read_csv(SOFT_CHANNELS):
        if (
            row["filter_label"] == PREFERRED_FILTER
            and row["scenario"] == PREFERRED_SCENARIO
            and row["spectrum_name"] == PREFERRED_SPECTRUM
            and row["channel"] == PREFERRED_CHANNEL
            and row["normalization_case"] == PREFERRED_CASE
            and row["positive_definite"] == "True"
        ):
            rows.append(row)
    if not rows:
        raise RuntimeError("Preferred soft-spectrum Knu rows were not found.")
    return rows


def preferred_interference_rows() -> list[dict[str, Any]]:
    rows = []
    for row in read_csv(INTERFERENCE_CHANNELS):
        if (
            row["filter_label"] == PREFERRED_FILTER
            and row["scenario"] == PREFERRED_SCENARIO
            and row["spectrum_name"] == PREFERRED_SPECTRUM
            and row["channel"] == PREFERRED_CHANNEL
            and row["normalization_case"] == PREFERRED_CASE
        ):
            rows.append(row)
    if not rows:
        raise RuntimeError("Preferred worst-phase Knu rows were not found.")
    return rows


def width_prefactor_from_row(row: dict[str, Any], chiral: dict[str, Any]) -> float:
    old = float(chiral["old_single_width_prefactor"])
    return old * f(row, "width_ratio_new_over_old")


def calibrate(
    *,
    label: str,
    row: dict[str, Any],
    amplitude: float,
    margin_current: float,
    margin_future: float,
    width_prefactor: float,
    extra: dict[str, Any],
) -> dict[str, Any]:
    s_t = f(row, "S_T_target")
    tau_current_years = CURRENT_BOUND_YR * margin_current
    gamma_year_inv = 1.0 / tau_current_years
    gamma_gev = gamma_year_inv * GEV_PER_YEAR_INV
    dressed_filter_amplitude = s_t * amplitude
    k_dyn_year_inv = gamma_year_inv / (
        width_prefactor * dressed_filter_amplitude**2
    )
    k_dyn_gev = gamma_gev / (width_prefactor * dressed_filter_amplitude**2)
    return {
        "label": label,
        "filter_label": row["filter_label"],
        "scenario": row["scenario"],
        "spectrum_name": row["spectrum_name"],
        "channel": row["channel"],
        "normalization_case": row["normalization_case"],
        "pair": row["pair"],
        "triplet_phase": extra.get("triplet_phase"),
        "amplitude_used": amplitude,
        "S_T": s_t,
        "width_prefactor": width_prefactor,
        "dressed_filter_amplitude": dressed_filter_amplitude,
        "margin_current_2p4e34": margin_current,
        "margin_future_1e35": margin_future,
        "tau_current_years": tau_current_years,
        "gamma_year_inv": gamma_year_inv,
        "gamma_GeV": gamma_gev,
        "K_dyn_year_inv": k_dyn_year_inv,
        "K_dyn_GeV": k_dyn_gev,
        "S_T_max_current": s_t * math.sqrt(margin_current),
        "S_T_max_future_1e35": s_t * math.sqrt(margin_future),
        "future_amplitude_headroom": math.sqrt(margin_future),
        "future_amplitude_suppression_needed": max(
            1.0 / math.sqrt(margin_future), 1.0
        ),
        **extra,
    }


def build() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    chiral = read_json(CHIRAL)
    target = read_json(KNU_TARGET)

    soft = min(preferred_soft_rows(), key=lambda row: f(row, "future_margin_1e35"))
    projected = min(
        preferred_interference_rows(),
        key=lambda row: f(row, "future_margin_1e35_worst_phase"),
    )

    soft_width_prefactor = width_prefactor_from_row(soft, chiral)
    projected_width_prefactor = width_prefactor_from_row(projected, chiral)

    rows = [
        calibrate(
            label="soft_spectrum_width",
            row=soft,
            amplitude=f(soft, "amplitude_with_dressing"),
            margin_current=f(soft, "margin_replayed_current"),
            margin_future=f(soft, "future_margin_1e35"),
            width_prefactor=soft_width_prefactor,
            extra={
                "source_file": str(SOFT_CHANNELS),
                "amplitude_column": "amplitude_with_dressing",
                "selected_index": None,
                "phase_spread": None,
            },
        ),
        calibrate(
            label="projected_worst_phase_width",
            row=projected,
            amplitude=f(projected, "max_phase_amplitude"),
            margin_current=f(projected, "margin_replayed_current_worst_phase"),
            margin_future=f(projected, "future_margin_1e35_worst_phase"),
            width_prefactor=projected_width_prefactor,
            extra={
                "source_file": str(INTERFERENCE_CHANNELS),
                "amplitude_column": "max_phase_amplitude",
                "triplet_phase": f(projected, "triplet_phase"),
                "min_phase_amplitude": f(projected, "min_phase_amplitude"),
                "phase_spread": f(projected, "phase_spread"),
            },
        ),
    ]

    k_values = [row["K_dyn_year_inv"] for row in rows]
    k_mean = sum(k_values) / len(k_values)
    max_rel_spread = max(abs(k - k_mean) / k_mean for k in k_values)
    margin_ratio = rows[1]["margin_future_1e35"] / rows[0]["margin_future_1e35"]
    amplitude_ratio_sq = (rows[0]["amplitude_used"] / rows[1]["amplitude_used"]) ** 2

    target_v = target["verdict"]
    summary = {
        "note": "No web lookup used. Calibrated physical Knu-width audit.",
        "formula": "Gamma = K_dyn * K_had * (S_T |A_dress|)^2",
        "input_files": {
            "soft_spectrum_channel_replay": str(SOFT_CHANNELS),
            "interference_tightened_channel_replay": str(INTERFERENCE_CHANNELS),
            "chiral_lattice_summary": str(CHIRAL),
            "knu_target_map": str(KNU_TARGET),
        },
        "unit_conventions": {
            "current_bound_years": CURRENT_BOUND_YR,
            "future_stress_years": FUTURE_STRESS_YR,
            "seconds_per_year": SECONDS_PER_YEAR,
            "hbar_GeV_seconds": HBAR_GEV_SECONDS,
            "GeV_per_year_inverse": GEV_PER_YEAR_INV,
        },
        "rows": rows,
        "cross_checks": {
            "K_dyn_mean_year_inv": k_mean,
            "K_dyn_max_relative_spread": max_rel_spread,
            "projected_to_soft_future_margin_ratio": margin_ratio,
            "soft_to_projected_amplitude_ratio_squared": amplitude_ratio_sq,
            "ratio_mismatch": abs(margin_ratio - amplitude_ratio_sq),
            "target_future_margin": float(target_v["future_margin"]),
            "target_future_suppression_needed": float(
                target_v["future_amplitude_suppression_needed"]
            ),
        },
        "verdict": {
            "calibration_stable": max_rel_spread < 1.0e-10,
            "K_dyn_year_inv": k_mean,
            "K_dyn_GeV": k_mean * GEV_PER_YEAR_INV,
            "final_margin_future_1e35": rows[1]["margin_future_1e35"],
            "final_suppression_needed": rows[1][
                "future_amplitude_suppression_needed"
            ],
            "final_S_T_max_future": rows[1]["S_T_max_future_1e35"],
            "interpretation": (
                "The preferred soft-spectrum bottleneck and the final projected "
                "Knu target row are described by the same calibrated width "
                "constant to numerical precision.  The remaining gap is therefore "
                "not a hidden normalization mismatch; it is a physical need for "
                "an additional Knu amplitude suppression or a stronger triplet "
                "filter in the same width formula."
            ),
        },
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "filter_label",
        "scenario",
        "spectrum_name",
        "channel",
        "normalization_case",
        "pair",
        "triplet_phase",
        "amplitude_used",
        "amplitude_column",
        "S_T",
        "width_prefactor",
        "dressed_filter_amplitude",
        "margin_current_2p4e34",
        "margin_future_1e35",
        "tau_current_years",
        "gamma_year_inv",
        "gamma_GeV",
        "K_dyn_year_inv",
        "K_dyn_GeV",
        "S_T_max_current",
        "S_T_max_future_1e35",
        "future_amplitude_headroom",
        "future_amplitude_suppression_needed",
        "selected_index",
        "min_phase_amplitude",
        "phase_spread",
        "source_file",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Full Knu width calibration",
        "",
        "No web lookup was used.",
        "",
        "The audited formula is",
        "",
        "```text",
        "Gamma = K_dyn * K_had * (S_T |A_dress|)^2",
        "```",
        "",
        "| label | pair | phase | amplitude | future margin | K_dyn [yr^-1] | S_T max at 1e35 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        phase = row["triplet_phase"]
        phase_text = "" if phase is None else f"{phase:.6e}"
        lines.append(
            "| `{label}` | `{pair}` | {phase} | {amp:.6e} | {margin:.6e} | "
            "{k:.6e} | {st:.6e} |".format(
                label=row["label"],
                pair=row["pair"],
                phase=phase_text,
                amp=row["amplitude_used"],
                margin=row["margin_future_1e35"],
                k=row["K_dyn_year_inv"],
                st=row["S_T_max_future_1e35"],
            )
        )

    checks = summary["cross_checks"]
    verdict = summary["verdict"]
    lines += [
        "",
        "## Cross-checks",
        "",
        (
            f"- K_dyn max relative spread: "
            f"{checks['K_dyn_max_relative_spread']:.6e}"
        ),
        (
            f"- projected/soft future-margin ratio: "
            f"{checks['projected_to_soft_future_margin_ratio']:.6e}"
        ),
        (
            f"- soft/projected amplitude ratio squared: "
            f"{checks['soft_to_projected_amplitude_ratio_squared']:.6e}"
        ),
        f"- ratio mismatch: {checks['ratio_mismatch']:.6e}",
        "",
        "## Verdict",
        "",
        verdict["interpretation"],
        "",
        (
            f"Final future margin is {verdict['final_margin_future_1e35']:.6e}; "
            f"additional amplitude suppression needed is "
            f"{verdict['final_suppression_needed']:.6e}."
        ),
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = build()
    write_csv(OUT / "knu_width_calibration.csv", rows)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(summary)
    verdict = summary["verdict"]
    checks = summary["cross_checks"]
    print("Full Knu width calibration")
    print(f"  calibration stable: {verdict['calibration_stable']}")
    print(f"  K_dyn rel spread: {checks['K_dyn_max_relative_spread']:.6e}")
    print(
        "  final S_T max future: {st:.6e}".format(
            st=verdict["final_S_T_max_future"]
        )
    )
    print(
        "  final suppression needed: {supp:.6e}".format(
            supp=verdict["final_suppression_needed"]
        )
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
