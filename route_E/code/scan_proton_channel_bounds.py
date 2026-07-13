#!/usr/bin/env python3
"""Replay proton-decay constraints with replaceable channel bounds.

The earlier item-4 scan used ``tau > 1e34 yr`` as a deliberately replaceable
benchmark.  This script keeps the local no-web convention and upgrades the
workflow: channel limits are explicit input numbers, every bound is propagated
by the exact lifetime scaling, and the 4pi-corrected downstream benchmarks are
checked without rerunning the expensive RGE scan.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTON_JSON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
DOWNSTREAM_JSON = ROOT / "output" / "corrected_downstream" / "corrected_downstream_summary.json"
OUT = ROOT / "output" / "proton_channel_bounds"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
LEGACY_TAU_TARGET_YEARS = 1.0e34
DISPLAY_ST = 1.0e-5

# These numbers are local inputs, not newly web-verified references.  The
# e+pi0 value is the stronger review-supplied target; the others are stress-test
# targets used to expose the required scaling before the bibliography pass.
DIM6_TARGETS = [
    {
        "channel": "p_to_e+pi0",
        "tau_bound_years": 2.4e34,
        "relative_width_factor": 1.0,
        "status": "review_supplied_local_input_verify_reference_before_submission",
    },
    {
        "channel": "p_to_mu+pi0",
        "tau_bound_years": 1.0e34,
        "relative_width_factor": 1.0,
        "status": "placeholder_channel_workflow_until_local_bound_is_supplied",
    },
]

DIM5_TARGETS = [
    {
        "channel": "p_to_K+nu_bar",
        "tau_bound_years": 1.0e34,
        "relative_width_factor": 1.0,
        "status": "legacy_channel_target",
    },
    {
        "channel": "p_to_K+nu_bar_strong",
        "tau_bound_years": 2.4e34,
        "relative_width_factor": 1.0,
        "status": "stress_test_same_strength_as_e+pi0",
    },
    {
        "channel": "p_to_K0mu+",
        "tau_bound_years": 1.0e34,
        "relative_width_factor": 1.0,
        "status": "placeholder_channel_workflow_until_local_bound_is_supplied",
    },
    {
        "channel": "future_1e35_stress",
        "tau_bound_years": 1.0e35,
        "relative_width_factor": 1.0,
        "status": "future_sensitivity_stress_test",
    },
]


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def mx_required(
    alpha_g_inv: float,
    tau_bound_years: float,
    width_prefactor_gev5: float,
    relative_width_factor: float,
) -> float:
    """Dimension-6 lower mass from Gamma = K r (g_G^2/M_X^2)^2."""
    g2 = 4.0 * math.pi / alpha_g_inv
    width_target = HBAR_GEV_S / (tau_bound_years * SECONDS_PER_YEAR)
    return ((width_prefactor_gev5 * relative_width_factor * g2 * g2) / width_target) ** 0.25


def st_max_from_lifetime(
    tau_at_display_st_years: float,
    tau_bound_years: float,
    relative_width_factor: float,
    display_st: float = DISPLAY_ST,
) -> float:
    """Triplet-filter upper bound using tau(S_T) proportional to S_T^{-2}."""
    return display_st * math.sqrt(tau_at_display_st_years / (tau_bound_years * relative_width_factor))


def corrected_benchmark_rows(payload: dict[str, object]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    base = payload["yukawa_4pi_audit"]["corrected_base_best"]
    fixed_r50 = payload["fixed_R50"]["best_safe_or_best"]
    rows.append({"name": "corrected_cache_base", **extract_benchmark(base)})
    rows.append({"name": "fixed_R50_mediator", **extract_benchmark(fixed_r50)})

    for row in payload["R_window"]["rows"]:
        r_value = float(row["R"])
        if r_value in (20.0, 50.0, 100.0, 200.0):
            rows.append({"name": f"R_window_R{r_value:g}", **extract_benchmark(row)})
    return rows


def extract_benchmark(row: dict[str, object]) -> dict[str, float]:
    mg = float(row.get("MG_GeV", 7.079458e15))
    return {
        "MG_GeV": mg,
        "M_X_GeV": float(row.get("M_X_GeV", mg)),
        "alphaG_inv": float(row["alphaG_inv"]),
        "tau_dim6_years": float(row["tau_dim6_years"]),
        "tau_dim5_ST1e-5_years": float(row["tau_dim5_target_filter_years"]),
    }


def channel_replay_rows(
    benchmarks: list[dict[str, float | str]],
    width_prefactor_gev5: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for bench in benchmarks:
        alpha_g_inv = float(bench["alphaG_inv"])
        m_x = float(bench["M_X_GeV"])
        for target in DIM6_TARGETS:
            tau_bound = float(target["tau_bound_years"])
            rel = float(target["relative_width_factor"])
            mx_min = mx_required(alpha_g_inv, tau_bound, width_prefactor_gev5, rel)
            tau = float(bench["tau_dim6_years"]) / rel
            rows.append(
                {
                    "benchmark": bench["name"],
                    "operator_dimension": 6,
                    "channel": target["channel"],
                    "target_status": target["status"],
                    "tau_bound_years": tau_bound,
                    "tau_pred_years": tau,
                    "margin_tau_over_bound": tau / tau_bound,
                    "M_X_GeV": m_x,
                    "M_X_min_GeV": mx_min,
                    "rho_X_min": max(1.0, mx_min / m_x),
                    "passes": tau >= tau_bound and m_x >= mx_min,
                    "S_T_max": "",
                }
            )

        for target in DIM5_TARGETS:
            tau_bound = float(target["tau_bound_years"])
            rel = float(target["relative_width_factor"])
            tau = float(bench["tau_dim5_ST1e-5_years"]) / rel
            rows.append(
                {
                    "benchmark": bench["name"],
                    "operator_dimension": 5,
                    "channel": target["channel"],
                    "target_status": target["status"],
                    "tau_bound_years": tau_bound,
                    "tau_pred_years": tau,
                    "margin_tau_over_bound": tau / tau_bound,
                    "M_X_GeV": "",
                    "M_X_min_GeV": "",
                    "rho_X_min": "",
                    "passes": tau >= tau_bound,
                    "S_T_max": st_max_from_lifetime(
                        float(bench["tau_dim5_ST1e-5_years"]), tau_bound, rel
                    ),
                }
            )
    return rows


def mx_scaling_table(width_prefactor_gev5: float) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for alpha_g_inv in (24.0, 30.0, 40.0):
        old = mx_required(alpha_g_inv, LEGACY_TAU_TARGET_YEARS, width_prefactor_gev5, 1.0)
        for target in DIM6_TARGETS:
            tau_bound = float(target["tau_bound_years"])
            rel = float(target["relative_width_factor"])
            new = mx_required(alpha_g_inv, tau_bound, width_prefactor_gev5, rel)
            rows.append(
                {
                    "channel": target["channel"],
                    "alphaG_inv": alpha_g_inv,
                    "old_tau_bound_years": LEGACY_TAU_TARGET_YEARS,
                    "new_tau_bound_years": tau_bound,
                    "M_X_min_old_GeV": old,
                    "M_X_min_new_GeV": new,
                    "scaling": new / old,
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, object]) -> None:
    lines: list[str] = []
    lines.append("# Proton channel-bound replay")
    lines.append("")
    lines.append("No web lookup was used.  The channel limits below are explicit local inputs,")
    lines.append("so final paper use still requires the bibliography/reference pass.")
    lines.append("")
    lines.append("## Scaling formulas")
    lines.append("")
    lines.append("For dimension-6 exchange,")
    lines.append("")
    lines.append("```text")
    lines.append("Gamma_6 = K r_ch (g_G^2/M_X^2)^2,")
    lines.append("M_X^min = [K r_ch g_G^4 tau_bound sec/yr / hbar]^(1/4).")
    lines.append("```")
    lines.append("")
    lines.append("For the dressed colored-triplet estimate at fixed spectrum,")
    lines.append("")
    lines.append("```text")
    lines.append("tau_5(S_T) = tau_5(S0) (S0/S_T)^2 / r_ch,")
    lines.append("S_T^max = S0 sqrt[tau_5(S0)/(r_ch tau_bound)].")
    lines.append("```")
    lines.append("")
    lines.append("## Dimension-6 mass rescaling")
    lines.append("")
    lines.append("| channel | alphaG^-1 | old M_X min [GeV] | new M_X min [GeV] | scaling |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in summary["mx_scaling"]:
        if row["channel"] == "p_to_e+pi0":
            lines.append(
                f"| `{row['channel']}` | {row['alphaG_inv']:.0f} | "
                f"{row['M_X_min_old_GeV']:.6e} | {row['M_X_min_new_GeV']:.6e} | "
                f"{row['scaling']:.6f} |"
            )
    lines.append("")
    lines.append("## Corrected downstream benchmark checks")
    lines.append("")
    lines.append("| benchmark | channel | dim | tau pred [yr] | tau bound [yr] | margin | extra bound | pass |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|")
    for row in summary["channel_replay"]:
        if row["benchmark"] not in ("fixed_R50_mediator", "R_window_R200"):
            continue
        if row["channel"] not in ("p_to_e+pi0", "p_to_K+nu_bar", "p_to_K+nu_bar_strong", "future_1e35_stress"):
            continue
        extra = row["M_X_min_GeV"] if row["operator_dimension"] == 6 else row["S_T_max"]
        lines.append(
            f"| `{row['benchmark']}` | `{row['channel']}` | {row['operator_dimension']} | "
            f"{row['tau_pred_years']:.6e} | {row['tau_bound_years']:.6e} | "
            f"{row['margin_tau_over_bound']:.3f} | {float(extra):.6e} | {row['passes']} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "The corrected branch easily passes the stronger dimension-6 e+pi0 target. "
        "The dimension-5 channel is the fragile one: at S_T=1e-5 it passes the "
        "1e34 and 2.4e34 stress tests, but a 1e35 yr future target would require "
        "S_T below the displayed 1e-5 design value."
    )
    (OUT / "proton_channel_bounds_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    proton = load_json(PROTON_JSON)
    downstream = load_json(DOWNSTREAM_JSON)
    width_prefactor = float(proton["hadronic_constants"]["width_prefactor_GeV5"])
    benchmarks = corrected_benchmark_rows(downstream)
    scaling = mx_scaling_table(width_prefactor)
    replay = channel_replay_rows(benchmarks, width_prefactor)

    summary = {
        "note": "No web lookup used; channel limits are explicit local inputs.",
        "input_files": {
            "proton_json": str(PROTON_JSON),
            "corrected_downstream_json": str(DOWNSTREAM_JSON),
        },
        "legacy_tau_target_years": LEGACY_TAU_TARGET_YEARS,
        "display_triplet_filter": DISPLAY_ST,
        "dimension6_targets": DIM6_TARGETS,
        "dimension5_targets": DIM5_TARGETS,
        "mx_scaling": scaling,
        "benchmarks": benchmarks,
        "channel_replay": replay,
        "all_displayed_current_targets_pass": all(
            bool(row["passes"])
            for row in replay
            if row["channel"] in ("p_to_e+pi0", "p_to_K+nu_bar", "p_to_K+nu_bar_strong")
        ),
    }
    write_csv(OUT / "proton_channel_bound_replay.csv", replay)
    write_csv(OUT / "proton_mx_bound_scaling.csv", scaling)
    (OUT / "proton_channel_bounds_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)

    r50 = next(row for row in replay if row["benchmark"] == "fixed_R50_mediator" and row["channel"] == "p_to_e+pi0")
    r50_k = next(row for row in replay if row["benchmark"] == "fixed_R50_mediator" and row["channel"] == "p_to_K+nu_bar_strong")
    future = next(row for row in replay if row["benchmark"] == "fixed_R50_mediator" and row["channel"] == "future_1e35_stress")
    print("Proton channel-bound replay")
    print(f"  e+pi0 stronger-target margin at R=50: {r50['margin_tau_over_bound']:.3f}")
    print(f"  K+nu strong-target S_T max at R=50: {r50_k['S_T_max']:.6e}")
    print(f"  Future 1e35 S_T max at R=50: {future['S_T_max']:.6e}")
    print(f"  all displayed current targets pass: {summary['all_displayed_current_targets_pass']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
