#!/usr/bin/env python3
"""Propagate the tightened triplet filter into downstream RGE/proton cards.

No web lookup is used.  This is an algebraic replay on top of the already
generated corrected downstream cache.  The previous cache displays
tau_d5(S_T=1e-5).  For the tightened d=5 proxy target the only required
rescaling is

    tau_d5(S_T) = tau_d5(1e-5) (1e-5/S_T)^2.

The scan does not claim a final physical d=5 calculation; it checks whether
the RGE/threshold benchmark cards remain viable once the hadronic-envelope
proxy demands S_T=7.5e-6.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DOWNSTREAM = ROOT / "output" / "corrected_downstream" / "corrected_downstream_summary.json"
TIGHTENED = ROOT / "output" / "tightened_triplet_filter" / "summary.json"
OUT = ROOT / "output" / "tightened_downstream_proton"

DISPLAY_ST = 1.0e-5
CURRENT_STRONG_D5_TARGET_YEARS = 2.4e34
FUTURE_D5_TARGET_YEARS = 1.0e35


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def tau_at_target(tau_at_display: float, st_target: float) -> float:
    return tau_at_display * (DISPLAY_ST / st_target) ** 2


def compact_benchmark(name: str, row: dict[str, Any], st_target: float) -> dict[str, Any]:
    tau_display = float(row.get("tau_dim5_target_filter_years", row.get("tau_dim5_years")))
    tau_tight = tau_at_target(tau_display, st_target)
    triplet_required = float(row.get("triplet_filter_required", row.get("filter_required", float("nan"))))
    return {
        "benchmark": name,
        "S_T_target": st_target,
        "MG_GeV": float(row.get("MG_GeV", 7.079458e15)),
        "alphaG_inv": float(row["alphaG_inv"]),
        "M_Sigma3_GeV": float(row["M_Sigma3_GeV"]),
        "M_Sigma8_GeV": float(row["M_Sigma8_GeV"]),
        "tau_dim6_years": float(row["tau_dim6_years"]),
        "tau_dim5_ST_1e_minus_5_years": tau_display,
        "tau_dim5_tightened_years": tau_tight,
        "margin_dim5_current_2p4e34": tau_tight / CURRENT_STRONG_D5_TARGET_YEARS,
        "margin_dim5_future_1e35": tau_tight / FUTURE_D5_TARGET_YEARS,
        "triplet_filter_required": triplet_required,
        "triplet_filter_required_over_target": triplet_required / st_target,
        "passes_current_dim5_target": tau_tight >= CURRENT_STRONG_D5_TARGET_YEARS,
        "passes_future_1e35_target": tau_tight >= FUTURE_D5_TARGET_YEARS,
        "passes_filter_requirement": st_target <= triplet_required,
    }


def collect_benchmarks(downstream: dict[str, Any], st_target: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base = downstream["yukawa_4pi_audit"]["corrected_base_best"]
    rows.append(compact_benchmark("corrected_base", base, st_target))

    fixed = downstream["fixed_R50"]["best_safe_or_best"]
    rows.append(compact_benchmark("fixed_R50", fixed, st_target))

    for r_value in ["20.0", "50.0", "100.0", "200.0"]:
        detail = downstream["R_window"]["detailed"][r_value]["scan_diagnostics"]["best_safe_or_best"]
        rows.append(compact_benchmark(f"R_window_R{float(r_value):g}", detail, st_target))

    # Goldstone-locked entries were explicitly scanned at S_T=1e-5.  Since
    # 7.5e-6 is less stringent, the locked epsilon bound is inherited; only the
    # displayed d=5 lifetime is rescaled.
    for r_value in ["50", "200"]:
        key = f"R={r_value},S_T=1.0e-05"
        locked = downstream["goldstone_locking"]["max_allowed_by_R_and_ST"][key]
        if locked is not None:
            rows.append(compact_benchmark(f"goldstone_locked_R{r_value}", locked, st_target))
            rows[-1]["epsilon_abs"] = float(locked["epsilon_abs"])
            rows[-1]["epsilon_over_lock_max"] = float(locked["epsilon_over_lock_max"])
    return rows


def summarize(rows: list[dict[str, Any]], st_target: float, tightened: dict[str, Any]) -> dict[str, Any]:
    selected = [row for row in rows if row["benchmark"] in {"R_window_R200", "goldstone_locked_R200", "fixed_R50"}]
    all_current = all(row["passes_current_dim5_target"] and row["passes_filter_requirement"] for row in selected)
    min_current_margin = min(row["margin_dim5_current_2p4e34"] for row in selected)
    min_filter_margin = min(row["triplet_filter_required_over_target"] for row in selected)
    return {
        "note": "No web lookup used. Tightened S_T replay for corrected downstream RGE/proton cards.",
        "input_files": {
            "corrected_downstream": str(DOWNSTREAM),
            "tightened_triplet_filter": str(TIGHTENED),
        },
        "display_S_T": DISPLAY_ST,
        "tightened_S_T": st_target,
        "lifetime_scaling": (DISPLAY_ST / st_target) ** 2,
        "benchmarks": rows,
        "selected_branch_checks": selected,
        "verdict": {
            "selected_current_dim5_targets_pass": all_current,
            "selected_min_current_margin": min_current_margin,
            "selected_min_filter_margin": min_filter_margin,
            "R200_tau_dim5_tightened_years": next(row for row in rows if row["benchmark"] == "R_window_R200")["tau_dim5_tightened_years"],
            "R200_current_margin": next(row for row in rows if row["benchmark"] == "R_window_R200")["margin_dim5_current_2p4e34"],
            "R200_future_1e35_margin": next(row for row in rows if row["benchmark"] == "R_window_R200")["margin_dim5_future_1e35"],
            "chiral_proxy_preferred_margin": tightened["verdict"]["preferred_rounded_worst_margin"],
            "interpretation": (
                "The corrected RGE/proton cards remain alive after replacing the display "
                "S_T=1e-5 by the tightened S_T=7.5e-6.  The R=200 branch exceeds the "
                "2.4e34 yr d=5 stress target, but it does not reach a future 1e35 yr "
                "stress target without further suppression."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Tightened downstream proton replay",
        "",
        "No web lookup was used.",
        "",
        "The replay uses",
        "",
        "```text",
        "tau_d5(S_T) = tau_d5(1e-5) (1e-5/S_T)^2.",
        "```",
        "",
        f"Tightened target: `S_T={summary['tightened_S_T']:.6e}`.",
        f"Lifetime scaling: `{summary['lifetime_scaling']:.6f}`.",
        "",
        "| benchmark | tau_d5 tight [yr] | margin 2.4e34 | margin 1e35 | filter req / target | pass current |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary["benchmarks"]:
        if row["benchmark"] in {"fixed_R50", "R_window_R200", "goldstone_locked_R200", "corrected_base"}:
            lines.append(
                f"| `{row['benchmark']}` | {row['tau_dim5_tightened_years']:.6e} | "
                f"{row['margin_dim5_current_2p4e34']:.6f} | "
                f"{row['margin_dim5_future_1e35']:.6f} | "
                f"{row['triplet_filter_required_over_target']:.6f} | "
                f"{row['passes_current_dim5_target']} |"
            )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    downstream = read_json(DOWNSTREAM)
    tightened = read_json(TIGHTENED)
    st_target = float(next(t for t in tightened["targets"] if t["target_label"] == "rounded_7p5e-6")["S_T_target"])
    rows = collect_benchmarks(downstream, st_target)
    summary = summarize(rows, st_target, tightened)
    write_csv(OUT / "benchmark_replay.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    verdict = summary["verdict"]
    print("Tightened downstream proton replay")
    print(f"  S_T target: {st_target:.6e}")
    print(f"  R200 tau_d5: {verdict['R200_tau_dim5_tightened_years']:.6e} yr")
    print(f"  R200 current margin: {verdict['R200_current_margin']:.6f}")
    print(f"  R200 future 1e35 margin: {verdict['R200_future_1e35_margin']:.6f}")
    print(f"  selected current targets pass: {verdict['selected_current_dim5_targets_pass']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
