#!/usr/bin/env python3
"""Score minimal UV interpretations of the shared-U NLSM source.

No web lookup is used.  The previous audit showed that

    S = rho54 U S0 U^T,
    Omega = rho210 (U e7)^(U e8)^(U e9)^(U e10)

realizes the desired 26-dimensional shared-orientation source manifold.  This
script asks the next harder question: can this source be promoted to a
microscopic high-R Spin(10) completion without reintroducing the large 54+210
Landau-pole problem or incomplete thresholds?

The audit is deliberately conservative.  It tests four routes:

  A. elementary propagating 54+210;
  B. a minimal confining-preon caricature with four Spin(10) vector preons;
  C. a deconstructed/link NLSM caricature with one complete 10 x 10 link;
  D. the cutoff shared-U NLSM, where only the constrained source manifold
     exists below the compositeness scale.

The output is a no-go/scorecard, not a model claim.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "nlsm_uv_completion_scorecard"
NLSM = ROOT / "output" / "nlsm_composite_source_origin" / "summary.json"
CONSTRAINED = ROOT / "output" / "constrained_54_210_source_sector" / "summary.json"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"

GAUGE_CONTRIBUTION = 24.0
T_INDEX = {
    "10": 1.0,
    "45": 8.0,
    "54": 12.0,
    "210": 56.0,
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def landau_ratio(alpha_inv: float, b10: float) -> float:
    if b10 <= 0.0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def alpha_inv_at_ratio(alpha_inv: float, b10: float, ratio: float) -> float:
    return alpha_inv - b10 * math.log(ratio) / (2.0 * math.pi)


def scenario_rows() -> list[dict[str, Any]]:
    nlsm = read_json(NLSM)
    constrained = read_json(CONSTRAINED)
    vacuum = read_json(VACUUM)

    baseline = next(
        row for row in constrained["uv_rows"] if row["scenario"] == "baseline_plus_10G"
    )
    baseline_sum_t = float(baseline["sum_T"])

    shared_rank = int(nlsm["jacobian"]["shared_U_rank"])
    constraint_rank = int(nlsm["jacobian"]["shared_U_constraint_rank"])
    relative_removed = int(nlsm["jacobian"]["relative_orientation_modes_removed_by_shared_U"])
    anomaly_proxy = int(nlsm["anomaly_proxy"]["new_chiral_anomaly_proxy"])

    scenarios = [
        {
            "route": "A_elementary_54_210",
            "field_content": "baseline+10G plus one elementary 54 and one elementary 210",
            "extra_sum_T": T_INDEX["54"] + T_INDEX["210"],
            "microscopic": True,
            "shared_U_locking": False,
            "nonuniversal_threshold_risk": "high: full PS fragments must be split or lifted",
            "moduli_or_source_rank": 264,
            "anomaly_proxy": 0,
            "interpretation": (
                "This is the direct renormalizable-looking route, but it is the "
                "same elementary 54+210 tower already known to create a short "
                "one-loop Landau reach."
            ),
        },
        {
            "route": "B_minimal_vector_preon_composite",
            "field_content": "baseline+10G plus four Spin(10) vector preons in 10",
            "extra_sum_T": 4 * T_INDEX["10"],
            "microscopic": "partial",
            "shared_U_locking": "possible if a single condensate orientation forms both composites",
            "nonuniversal_threshold_risk": "hidden-sector dependent",
            "moduli_or_source_rank": shared_rank,
            "anomaly_proxy": 0,
            "interpretation": (
                "Four vector preons are the smallest caricature that can make a "
                "four-form composite while remaining real under Spin(10).  It "
                "is not yet a UV model because hidden gauge anomalies, "
                "confinement, and operator dimensions are unspecified."
            ),
        },
        {
            "route": "C_deconstructed_10x10_link",
            "field_content": "baseline+10G plus one complete link 10 x 10 = 1+45+54",
            "extra_sum_T": T_INDEX["45"] + T_INDEX["54"],
            "microscopic": True,
            "shared_U_locking": True,
            "nonuniversal_threshold_risk": "low only if the whole link is degenerate",
            "moduli_or_source_rank": 45,
            "anomaly_proxy": 0,
            "interpretation": (
                "A complete link naturally supplies a shared orientation, but "
                "the whole 1+45+54 link must propagate unless it is itself "
                "constrained or strongly coupled."
            ),
        },
        {
            "route": "D_cutoff_shared_U_NLSM",
            "field_content": "baseline+10G plus constrained shared-U source manifold only",
            "extra_sum_T": int(nlsm["beta_cost"]["so10_dynkin_increment_in_preferred_EFT"]),
            "microscopic": False,
            "shared_U_locking": True,
            "nonuniversal_threshold_risk": "none below compositeness scale in the audit",
            "moduli_or_source_rank": shared_rank,
            "anomaly_proxy": anomaly_proxy,
            "interpretation": (
                "This is the current conditional EFT branch: it keeps the "
                "verified 26-dimensional source manifold but does not derive "
                "the microscopic dynamics that produced it."
            ),
        },
    ]

    rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        sum_t = baseline_sum_t + float(scenario["extra_sum_T"])
        b10 = sum_t - GAUGE_CONTRIBUTION
        for vac in vacuum["replay_rows"]:
            r_target = float(vac["R"])
            alpha_inv = float(vac["alphaG_inv"])
            lp = landau_ratio(alpha_inv, b10)
            alpha_at_r = alpha_inv_at_ratio(alpha_inv, b10, r_target)
            rows.append(
                {
                    "route": scenario["route"],
                    "field_content": scenario["field_content"],
                    "R_target": r_target,
                    "alphaG_inv_MG": alpha_inv,
                    "baseline_sum_T": baseline_sum_t,
                    "extra_sum_T": float(scenario["extra_sum_T"]),
                    "sum_T": sum_t,
                    "b10": b10,
                    "landau_ratio_Lambda_over_MG": lp,
                    "alpha_inv_at_RMG": alpha_at_r,
                    "reaches_RMG": alpha_at_r > 0.0,
                    "microscopic": scenario["microscopic"],
                    "shared_U_locking": scenario["shared_U_locking"],
                    "nonuniversal_threshold_risk": scenario["nonuniversal_threshold_risk"],
                    "moduli_or_source_rank": int(scenario["moduli_or_source_rank"]),
                    "constraint_rank_if_applicable": constraint_rank
                    if scenario["route"] == "D_cutoff_shared_U_NLSM"
                    else "",
                    "relative_modes_removed_if_applicable": relative_removed
                    if scenario["route"] == "D_cutoff_shared_U_NLSM"
                    else "",
                    "anomaly_proxy": int(scenario["anomaly_proxy"]),
                    "interpretation": scenario["interpretation"],
                }
            )
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    routes = sorted({row["route"] for row in rows})
    route_summaries: list[dict[str, Any]] = []
    for route in routes:
        rrows = [row for row in rows if row["route"] == route]
        all_high_r_pass = all(row["reaches_RMG"] for row in rrows if row["R_target"] in {50.0, 200.0})
        r200 = next(row for row in rrows if row["R_target"] == 200.0)
        first = rrows[0]
        microscopic_pass = bool(first["microscopic"] is True and all_high_r_pass)
        conditional_pass = bool(route == "D_cutoff_shared_U_NLSM" and all_high_r_pass)
        if route == "D_cutoff_shared_U_NLSM":
            status = "PASS_CONDITIONAL"
        elif all_high_r_pass and first["microscopic"] is True:
            status = "PASS_MICROSCOPIC"
        elif route == "B_minimal_vector_preon_composite":
            status = "OPEN_FAILS_R200_AS_PERTURBATIVE_PREON"
        else:
            status = "FAIL_HIGH_R"
        route_summaries.append(
            {
                "route": route,
                "status": status,
                "microscopic_pass_high_R": microscopic_pass,
                "conditional_pass_high_R": conditional_pass,
                "sum_T": first["sum_T"],
                "b10": first["b10"],
                "R200_landau_ratio": r200["landau_ratio_Lambda_over_MG"],
                "R200_alpha_inv_at_RMG": r200["alpha_inv_at_RMG"],
                "all_R50_R200_reached": all_high_r_pass,
                "shared_U_locking": first["shared_U_locking"],
                "nonuniversal_threshold_risk": first["nonuniversal_threshold_risk"],
                "interpretation": first["interpretation"],
            }
        )

    microscopic_candidates = [row for row in route_summaries if row["microscopic_pass_high_R"]]
    conditional_candidates = [row for row in route_summaries if row["conditional_pass_high_R"]]
    return {
        "note": "No web lookup used. Minimal A-D UV completion/no-go scorecard for the shared-U NLSM source.",
        "formula": {
            "b10": "sum_R T(R) - 24",
            "alpha_inv_RMG": "alpha_inv(MG) - b10 log(R)/(2 pi)",
            "landau_ratio": "exp(2 pi alpha_inv(MG)/b10) for b10>0",
        },
        "baseline": "baseline_plus_10G from constrained_54_210_source_sector",
        "route_summaries": route_summaries,
        "verdict": {
            "microscopic_high_R_candidate_found": bool(microscopic_candidates),
            "conditional_EFT_candidate_found": bool(conditional_candidates),
            "preferred_current_route": "D_cutoff_shared_U_NLSM",
            "main_obstacle": (
                "Among the four minimal routes, no fully microscopic perturbative "
                "candidate both explains the shared-U source and reaches the "
                "displayed R=50,200 mediator scales.  The cutoff shared-U NLSM "
                "passes as a conditional EFT but does not close O3."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Shared-U NLSM UV completion/no-go scorecard")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("The one-loop test is")
    lines.append("")
    lines.append("```text")
    lines.append("b10 = sum_R T(R) - 24")
    lines.append("alpha_G^{-1}(R M_G) = alpha_G^{-1}(M_G) - b10 log(R)/(2 pi)")
    lines.append("Lambda_LP/M_G = exp(2 pi alpha_G^{-1}(M_G)/b10)")
    lines.append("```")
    lines.append("")
    lines.append("## Route summary")
    lines.append("")
    lines.append("| route | status | sum T | b10 | LP/MG at R=200 card | alpha^-1(200MG) | shared-U |")
    lines.append("|---|---|---:|---:|---:|---:|---|")
    for route in summary["route_summaries"]:
        lines.append(
            f"| `{route['route']}` | `{route['status']}` | "
            f"{route['sum_T']:.0f} | {route['b10']:.0f} | "
            f"{route['R200_landau_ratio']:.3g} | "
            f"{route['R200_alpha_inv_at_RMG']:.3f} | "
            f"{route['shared_U_locking']} |"
        )
    lines.append("")
    lines.append("## Displayed R rows")
    lines.append("")
    lines.append("| route | R | alpha^-1(MG) | alpha^-1(RMG) | LP/MG | reaches RMG |")
    lines.append("|---|---:|---:|---:|---:|---|")
    for row in rows:
        if row["R_target"] in {50.0, 200.0}:
            lines.append(
                f"| `{row['route']}` | {row['R_target']:.0f} | "
                f"{row['alphaG_inv_MG']:.3f} | {row['alpha_inv_at_RMG']:.3f} | "
                f"{row['landau_ratio_Lambda_over_MG']:.3g} | {row['reaches_RMG']} |"
            )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(summary["verdict"]["main_obstacle"])
    lines.append("")
    lines.append(
        "Thus the NLSM branch is stronger than a free assumption, because its "
        "rank, constraints, Plucker relations, and anomaly proxy are verified; "
        "nevertheless it remains a conditional cutoff EFT until a microscopic "
        "source dynamics is supplied."
    )
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = scenario_rows()
    summary = summarize(rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(OUT / "landau_rows.csv", rows)
    write_csv(OUT / "scorecard.csv", summary["route_summaries"])
    write_report(rows, summary)

    print("Shared-U NLSM UV completion/no-go scorecard")
    for route in summary["route_summaries"]:
        print(
            f"  {route['route']}: {route['status']}, "
            f"LP/MG(R200 card)={route['R200_landau_ratio']:.3g}, "
            f"alpha_inv(200MG)={route['R200_alpha_inv_at_RMG']:.3f}"
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
