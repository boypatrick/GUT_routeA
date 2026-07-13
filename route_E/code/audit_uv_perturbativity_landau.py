#!/usr/bin/env python3
"""SO(10) UV perturbativity / Landau-pole audit.

No web lookup is used.  This script checks a point raised by the 2026-05-07
review: complete Spin(10) multiplets are silent in the two-dimensional
low-energy matching plane, but they are not silent in the SO(10) gauge beta
function above M_G.

For N=1 SUSY SO(10),

    d alpha_G^{-1} / d log mu = - b_10 / (2 pi),
    b_10 = sum_R T(R) - 3 C_2(SO(10)),   C_2(adj)=8.

Thus C_2 gives a -24 gauge contribution.  If b_10>0, the one-loop Landau
scale from M_G is

    Lambda_LP / M_G = exp(2 pi alpha_G^{-1}(M_G) / b_10).

The audit is deliberately scenario-based because the manuscript must decide
which source/link/driver fields are propagating chiral multiplets and which
are better treated as constrained spurions or auxiliary fields.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "uv_perturbativity_landau"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"

C2_SO10_ADJ = 8.0
GAUGE_CONTRIBUTION = 3.0 * C2_SO10_ADJ

T_INDEX = {
    "10": 1.0,
    "16": 2.0,
    "45": 8.0,
    "54": 12.0,
    "120": 28.0,
    "126": 35.0,
    "126bar": 35.0,
    "210": 56.0,
}


SCENARIOS: list[dict[str, Any]] = [
    {
        "name": "minimal_yukawa_higgs_plus_matter",
        "interpretation": "three 16 matter fields plus 10_H, 120_H, and 126bar_H only",
        "counts": {"16": 3, "10": 1, "120": 1, "126bar": 1},
        "status_if_fails": "would already endanger a conventional large-rep SO(10) fit",
    },
    {
        "name": "documented_complete_drive_pairs",
        "interpretation": "three 54 link/driver pairs plus one 210 link/driver pair active at M_G",
        "counts": {"54": 6, "210": 2},
        "status_if_fails": "drive/source completion cannot be a perturbative propagating sector up to R M_G",
    },
    {
        "name": "projector_45_mediator_triplet_only",
        "interpretation": "light Sigma plus two 45 mediators counted as propagating 45 fields",
        "counts": {"45": 3},
        "status_if_fails": "45 mediator sector alone would be inconsistent",
    },
    {
        "name": "wrel_source_blocks",
        "interpretation": "W_rel-completed 11-field 54 block plus 3-field 210 source/link/driver block",
        "counts": {"54": 11, "210": 3},
        "status_if_fails": "full source block must be spurionic/constrained or live above the perturbative cutoff",
    },
    {
        "name": "yukawa_plus_documented_drive_pairs",
        "interpretation": "minimal Yukawa/matter sector plus documented complete drive pairs",
        "counts": {"16": 3, "10": 1, "120": 1, "126bar": 1, "54": 6, "210": 2},
        "status_if_fails": "large-rep flavor plus complete drive pairs are mutually UV-tight",
    },
    {
        "name": "wrel_plus_projector_plus_yukawa",
        "interpretation": "W_rel source blocks, 45 projector sector, and minimal Yukawa/matter fields",
        "counts": {"16": 3, "10": 1, "120": 1, "126bar": 1, "45": 3, "54": 11, "210": 3},
        "status_if_fails": "the fully propagating completion is only a low-cutoff EFT",
    },
]


def sum_dynkin(counts: dict[str, int]) -> float:
    return float(sum(T_INDEX[rep] * count for rep, count in counts.items()))


def load_vacuum_rows() -> list[dict[str, Any]]:
    payload = json.loads(VACUUM.read_text(encoding="utf-8"))
    return payload["replay_rows"]


def landau_ratio(alpha_inv: float, b10: float) -> float:
    if b10 <= 0.0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def alpha_at_ratio(alpha_inv: float, b10: float, ratio: float) -> float:
    return alpha_inv - b10 * math.log(ratio) / (2.0 * math.pi)


def audit_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        sum_t = sum_dynkin(scenario["counts"])
        b10 = sum_t - GAUGE_CONTRIBUTION
        for vac in load_vacuum_rows():
            target_r = float(vac["R"])
            alpha_inv = float(vac["alphaG_inv"])
            ratio_lp = landau_ratio(alpha_inv, b10)
            alpha_at_mmed = alpha_at_ratio(alpha_inv, b10, target_r)
            rows.append(
                {
                    "scenario": scenario["name"],
                    "interpretation": scenario["interpretation"],
                    "counts": dict(scenario["counts"]),
                    "sum_T": sum_t,
                    "b10": b10,
                    "R_target": target_r,
                    "alphaG_inv_MG": alpha_inv,
                    "MG_GeV": float(vac["MG_GeV"]),
                    "Mmed_GeV": float(vac["MG_GeV"]) * target_r,
                    "alpha_inv_at_RMG": alpha_at_mmed,
                    "landau_ratio_from_MG": ratio_lp,
                    "landau_scale_GeV": float(vac["MG_GeV"]) * ratio_lp
                    if math.isfinite(ratio_lp)
                    else math.inf,
                    "reaches_RMG_perturbatively": alpha_at_mmed > 0.0,
                    "margin_log": math.log(ratio_lp / target_r)
                    if math.isfinite(ratio_lp)
                    else math.inf,
                    "status_if_fails": scenario["status_if_fails"],
                }
            )
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "note": "No web lookup used. One-loop SO(10) UV perturbativity audit.",
        "formula": {
            "b10": "sum_R T(R) - 24",
            "landau_ratio": "exp(2*pi*alphaG_inv/b10) for b10>0",
        },
        "dynkin_indices": T_INDEX,
        "gauge_contribution_3C2": GAUGE_CONTRIBUTION,
        "scenarios": [],
    }
    for scenario in SCENARIOS:
        srows = [row for row in rows if row["scenario"] == scenario["name"]]
        min_ratio = min(row["landau_ratio_from_MG"] for row in srows)
        max_target_pass = max(
            (row["R_target"] for row in srows if row["reaches_RMG_perturbatively"]),
            default=0.0,
        )
        out["scenarios"].append(
            {
                "name": scenario["name"],
                "interpretation": scenario["interpretation"],
                "counts": scenario["counts"],
                "sum_T": srows[0]["sum_T"],
                "b10": srows[0]["b10"],
                "min_landau_ratio_over_displayed_R": min_ratio,
                "max_displayed_R_reached_perturbatively": max_target_pass,
                "all_displayed_R_pass": all(row["reaches_RMG_perturbatively"] for row in srows),
            }
        )
    hard = next(row for row in out["scenarios"] if row["name"] == "documented_complete_drive_pairs")
    out["verdict"] = {
        "documented_drive_pairs_pass_R50_R200": hard["all_displayed_R_pass"],
        "main_blocker": (
            "If the complete 54/210 source/drive multiplets are propagating at M_G, "
            "the displayed R=50 and R=200 mediator interpretation is not UV-perturbative "
            "at one loop. They must be spurionic/constrained, lifted above the cutoff, "
            "or replaced by a smaller representation completion."
        ),
    }
    return out


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "scenario",
        "sum_T",
        "b10",
        "R_target",
        "alphaG_inv_MG",
        "alpha_inv_at_RMG",
        "landau_ratio_from_MG",
        "reaches_RMG_perturbatively",
        "margin_log",
        "interpretation",
    ]
    with (OUT / "uv_landau_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# SO(10) UV perturbativity / Landau-pole audit")
    lines.append("")
    lines.append("No web lookup was used.  The audit separates low-energy threshold")
    lines.append("universality from UV perturbativity above `M_G`.")
    lines.append("")
    lines.append("For N=1 SUSY SO(10),")
    lines.append("")
    lines.append("```text")
    lines.append("d alpha_G^{-1}/d log mu = - b10/(2 pi)")
    lines.append("b10 = sum_R T(R) - 3 C2(SO(10)) = sum_R T(R) - 24")
    lines.append("Lambda_LP/M_G = exp(2 pi alpha_G^{-1}(M_G)/b10), b10>0")
    lines.append("```")
    lines.append("")
    lines.append("The Dynkin indices used locally are")
    lines.append("")
    lines.append("| rep | T(R) |")
    lines.append("|---|---:|")
    for rep, index in T_INDEX.items():
        lines.append(f"| `{rep}` | {index:g} |")
    lines.append("")
    lines.append("## Scenario summary")
    lines.append("")
    lines.append("| scenario | sum T | b10 | min LP/MG | max displayed R reached | all pass |")
    lines.append("|---|---:|---:|---:|---:|---|")
    for scenario in summary["scenarios"]:
        lines.append(
            f"| `{scenario['name']}` | {scenario['sum_T']:.0f} | {scenario['b10']:.0f} | "
            f"{scenario['min_landau_ratio_over_displayed_R']:.3g} | "
            f"{scenario['max_displayed_R_reached_perturbatively']:.0f} | "
            f"{scenario['all_displayed_R_pass']} |"
        )
    lines.append("")
    lines.append("## Displayed R points")
    lines.append("")
    lines.append("| scenario | R | alphaG^-1(MG) | alpha^-1(R MG) | LP/MG | pass |")
    lines.append("|---|---:|---:|---:|---:|---|")
    for row in rows:
        if row["scenario"] in {
            "minimal_yukawa_higgs_plus_matter",
            "documented_complete_drive_pairs",
            "wrel_source_blocks",
            "wrel_plus_projector_plus_yukawa",
        }:
            lines.append(
                f"| `{row['scenario']}` | {row['R_target']:.0f} | "
                f"{row['alphaG_inv_MG']:.3f} | {row['alpha_inv_at_RMG']:.3f} | "
                f"{row['landau_ratio_from_MG']:.3g} | {row['reaches_RMG_perturbatively']} |"
            )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("Complete multiplets are projected-threshold silent but not beta-function")
    lines.append("silent.  If the documented three 54 link/driver pairs and one 210")
    lines.append("link/driver pair are propagating already at `M_G`, then `sum T=184`,")
    lines.append("`b10=160`, and the one-loop Landau ratio is only a few times `M_G`.")
    lines.append("This does not reach the displayed `R=50` or `R=200` mediator scales.")
    lines.append("")
    lines.append("Therefore the current high-R completion must be interpreted in one of")
    lines.append("three ways: source/link/driver multiplets are constrained spurions rather")
    lines.append("than propagating chiral fields; the perturbative cutoff is only a few")
    lines.append("times `M_G`; or the projector/source realization must be replaced by a")
    lines.append("smaller-representation completion.")
    (OUT / "uv_landau_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = audit_rows()
    summary = summarize(rows)
    (OUT / "uv_landau_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(rows, summary)
    documented = [
        row for row in rows if row["scenario"] == "documented_complete_drive_pairs" and row["R_target"] in {50.0, 200.0}
    ]
    print("SO(10) UV perturbativity / Landau audit")
    for row in documented:
        print(
            f"  R={row['R_target']:.0f}: b10={row['b10']:.0f}, "
            f"LP/MG={row['landau_ratio_from_MG']:.3g}, "
            f"alpha_inv(RMG)={row['alpha_inv_at_RMG']:.3f}, "
            f"pass={row['reaches_RMG_perturbatively']}"
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
