#!/usr/bin/env python3
"""Bridge the best flavor rows to the crossed-120 d=5 null projector.

No web lookup is used.  The previous closure pipeline showed that the finite
CP1 two-kernel/Clebsch profiles do not close full flavor plus dressed d=5
proton decay.  Separately, the component audits found a crossed 120A/120B
triplet projector that kills the monitored LLLL Knu and RRRR uusd tensors, but
only as a post-Spin(10)/constrained triplet-sector target.

This script asks the sharpened question:

    If the crossed-120 projector is accepted as the triplet sector of the
    source-consistent best operator row, what remains open?

It deliberately distinguishes:

  * ordinary finite-Clebsch profiles, which are action-level tied to the
    doublet fit but fail future Knu stress;
  * crossed-120 triplet projector, which closes d=5 with huge margins but is
    only a conditional post-Spin(10)/constrained projector.

The output is a closure bridge, not a new free-parameter fit.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clebsch_null_flavor_d5"

FULL_PIPELINE = ROOT / "output" / "full_flavor_d5_pipeline" / "summary.json"
TWO_KERNEL = ROOT / "output" / "two_kernel_flavor_then_d5" / "summary.json"
TRIPLET_120 = ROOT / "output" / "triplet_120_rrrr_mass_matrix" / "summary.json"
CROSSED_120 = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"

STRICT_CKM = 1.0e-3
LOOSE_CKM = 5.0e-2
LOOSE_MASS = 2.0e-1


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def suppression_needed(margin: float) -> float:
    if margin <= 0.0:
        return math.inf
    return max(1.0 / math.sqrt(margin), 1.0)


def find_score_row(full: dict[str, Any], label: str) -> dict[str, Any]:
    return next(row for row in full["scoreboard"] if row["label"] == label)


def finite_projector_margin(crossed: dict[str, Any], kappa: float) -> float:
    row = next(item for item in crossed["finite_leakage_replay"] if float(item["kappa"]) == float(kappa))
    return float(row["worst_margin_1e35_replayed"])


def candidate(
    *,
    label: str,
    source_label: str,
    ckm_score: float,
    mass_score: float | None,
    future_margin: float,
    source_consistent: bool,
    projector: str,
    field_only_allowed: bool,
    interpretation: str,
) -> dict[str, Any]:
    return {
        "label": label,
        "source_label": source_label,
        "projector": projector,
        "ckm_score": ckm_score,
        "mass_score": mass_score,
        "future_margin_1e35": future_margin,
        "d5_suppression_needed": suppression_needed(future_margin),
        "source_consistent_with_triplet_audit": source_consistent,
        "field_only_spin10_allowed": field_only_allowed,
        "passes_strict_ckm": ckm_score < STRICT_CKM,
        "passes_loose_ckm": ckm_score < LOOSE_CKM,
        "passes_loose_mass": bool(mass_score is not None and mass_score <= LOOSE_MASS),
        "passes_future_d5": future_margin >= 1.0,
        "passes_conditional_loose_closure": bool(
            source_consistent
            and ckm_score < LOOSE_CKM
            and mass_score is not None
            and mass_score <= LOOSE_MASS
            and future_margin >= 1.0
        ),
        "passes_conditional_strict_closure": bool(
            source_consistent
            and ckm_score < STRICT_CKM
            and mass_score is not None
            and mass_score <= LOOSE_MASS
            and future_margin >= 1.0
        ),
        "ckm_improvement_factor_to_strict": max(ckm_score / STRICT_CKM, 1.0),
        "interpretation": interpretation,
    }


def build() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    full = read_json(FULL_PIPELINE)
    two = read_json(TWO_KERNEL)
    triplet = read_json(TRIPLET_120)
    crossed = read_json(CROSSED_120)

    best_ckm = find_score_row(full, "operator_best_ckm")
    best_knu = find_score_row(full, "operator_best_knu")
    source_label = triplet["verdict"]["best_physical_row"]["label"]
    crossed_source_label = "d5_both_F_minus_0"

    kappa30_margin = finite_projector_margin(crossed, 30.0)
    kappa100_margin = finite_projector_margin(crossed, 100.0)
    crossed_field_only_allowed = bool(crossed["verdict"]["field_only_unbroken_spin10_projector_possible"])
    crossed_conditional_allowed = bool(crossed["verdict"]["ps_eft_or_constrained_projector_possible"])

    rows = [
        candidate(
            label="ordinary_best_ckm_profile",
            source_label=str(best_ckm["details"]["candidate_label"]),
            ckm_score=float(best_ckm["ckm_score"]),
            mass_score=float(best_ckm["mass_score"]),
            future_margin=float(best_ckm["future_margin_1e35"]),
            source_consistent=False,
            projector="finite_doublet_reused_Clebsch_profile",
            field_only_allowed=True,
            interpretation="Best CKM row from the ordinary operator scan; close to strict CKM but still future-d5 unsafe.",
        ),
        candidate(
            label="ordinary_source_consistent_best_knu",
            source_label=str(best_knu["details"]["candidate_label"]),
            ckm_score=float(best_knu["ckm_score"]),
            mass_score=float(best_knu["mass_score"]),
            future_margin=float(best_knu["future_margin_1e35"]),
            source_consistent=True,
            projector="finite_doublet_reused_Clebsch_profile",
            field_only_allowed=True,
            interpretation="The source-consistent row used by the 120 null audits; mass is good and CKM is loose-good, but finite profiles fail future d=5.",
        ),
        candidate(
            label="source_consistent_crossed120_kappa30",
            source_label=str(best_knu["details"]["candidate_label"]),
            ckm_score=float(best_knu["ckm_score"]),
            mass_score=float(best_knu["mass_score"]),
            future_margin=kappa30_margin,
            source_consistent=True,
            projector="crossed_120A_120B_finite_lift_kappa30",
            field_only_allowed=crossed_field_only_allowed,
            interpretation="Accept the crossed 120A/120B post-Spin(10) projector with the most conservative audited finite lift; d=5 closes, but strict CKM remains above 1e-3.",
        ),
        candidate(
            label="source_consistent_crossed120_kappa100",
            source_label=str(best_knu["details"]["candidate_label"]),
            ckm_score=float(best_knu["ckm_score"]),
            mass_score=float(best_knu["mass_score"]),
            future_margin=kappa100_margin,
            source_consistent=True,
            projector="crossed_120A_120B_finite_lift_kappa100",
            field_only_allowed=crossed_field_only_allowed,
            interpretation="Same crossed projector with the preferred kappa=100 lift; d=5 is extremely safe, but the field-only Spin(10) grading no-go remains.",
        ),
    ]

    strict = [row for row in rows if row["passes_conditional_strict_closure"]]
    loose = [row for row in rows if row["passes_conditional_loose_closure"]]
    source_projector = next(row for row in rows if row["label"] == "source_consistent_crossed120_kappa30")

    payload = {
        "note": "No web lookup used. Clebsch-null flavor plus d=5 bridge.",
        "formulae": {
            "LLLL_null": "If Y_QQ^T is pure antisymmetric 120-like, sym(Y_QQ^T)=0 and LLLL Knu vanishes.",
            "crossed_RRRR_null": "Choose 10-10 triplet source G_A and 10-5bar triplet source G_B-like so the RRRR uusd map is also nulled.",
            "finite_lift": "W_kappa = U diag(1,1/kappa,1/kappa,1/kappa) V^dagger; leakage is O(1/kappa).",
        },
        "input_files": {
            "full_pipeline": str(FULL_PIPELINE.relative_to(ROOT)),
            "two_kernel": str(TWO_KERNEL.relative_to(ROOT)),
            "triplet_120": str(TRIPLET_120.relative_to(ROOT)),
            "crossed_120": str(CROSSED_120.relative_to(ROOT)),
        },
        "source_consistency": {
            "triplet_120_source_label": crossed_source_label,
            "used_source_label": str(best_knu["details"]["candidate_label"]),
            "consistent": str(best_knu["details"]["candidate_label"]) == crossed_source_label,
            "component_physical_row": source_label,
        },
        "selection_rule": {
            "field_only_unbroken_spin10_projector_possible": crossed_field_only_allowed,
            "ps_eft_or_constrained_projector_possible": crossed_conditional_allowed,
            "interpretation": crossed["verdict"]["interpretation"],
        },
        "rows": rows,
        "verdict": {
            "strict_conditional_closure_count": len(strict),
            "loose_conditional_closure_count": len(loose),
            "d5_gap_closed_by_crossed_projector": bool(source_projector["passes_future_d5"]),
            "strict_ckm_still_open": not bool(source_projector["passes_strict_ckm"]),
            "source_consistent_ckm_score": source_projector["ckm_score"],
            "source_consistent_mass_score": source_projector["mass_score"],
            "source_consistent_kappa30_margin": source_projector["future_margin_1e35"],
            "ckm_factor_to_strict": source_projector["ckm_improvement_factor_to_strict"],
            "publication_complete": bool(len(strict) > 0 and crossed_conditional_allowed),
            "interpretation": (
                "The crossed-120 projector conditionally closes the dressed d=5 future-stress gap "
                "for the source-consistent flavor row, but it does not make the paper complete. "
                "The projector is forbidden by a field-only unbroken Spin(10) grading and must be "
                "realized as a post-Spin(10)/constrained triplet-sector operator; moreover the "
                "source-consistent CKM score remains above the strict 1e-3 target by the displayed factor."
            ),
        },
    }
    return payload, rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "label",
        "source_label",
        "projector",
        "ckm_score",
        "mass_score",
        "future_margin_1e35",
        "d5_suppression_needed",
        "source_consistent_with_triplet_audit",
        "field_only_spin10_allowed",
        "passes_strict_ckm",
        "passes_loose_ckm",
        "passes_loose_mass",
        "passes_future_d5",
        "passes_conditional_loose_closure",
        "passes_conditional_strict_closure",
        "ckm_improvement_factor_to_strict",
        "interpretation",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    v = payload["verdict"]
    lines = [
        "# Clebsch-null flavor plus d=5 bridge",
        "",
        "No web lookup was used.",
        "",
        "The key algebraic point is",
        "",
        "```text",
        "Y_QQ^T pure 120-like antisymmetric => sym(Y_QQ^T)=0 => LLLL Knu null.",
        "Crossed 120_A/120_B: 10-10 uses G_A, 10-5bar uses G_B-like.",
        "Finite lift W_kappa leaks as O(1/kappa).",
        "```",
        "",
        "| candidate | CKM score | mass score | future Knu margin | strict? | loose? | field-only Spin10? |",
        "|---|---:|---:|---:|---|---|---|",
    ]
    for row in payload["rows"]:
        lines.append(
            f"| `{row['label']}` | {row['ckm_score']:.6e} | "
            f"{row['mass_score']:.6e} | {row['future_margin_1e35']:.6e} | "
            f"{row['passes_conditional_strict_closure']} | "
            f"{row['passes_conditional_loose_closure']} | "
            f"{row['field_only_spin10_allowed']} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        f"Strict conditional closures: `{v['strict_conditional_closure_count']}`.",
        f"Loose conditional closures: `{v['loose_conditional_closure_count']}`.",
        f"Source-consistent crossed-projector Knu margin: `{v['source_consistent_kappa30_margin']:.6e}`.",
        f"CKM factor still needed to reach strict target: `{v['ckm_factor_to_strict']:.6e}`.",
        "",
        v["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload, rows = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "closure_rows.csv", rows)
    write_report(payload)
    v = payload["verdict"]
    print("Clebsch-null flavor + d=5 bridge")
    print(f"  strict closures: {v['strict_conditional_closure_count']}")
    print(f"  loose closures: {v['loose_conditional_closure_count']}")
    print(f"  d5 gap closed by crossed projector: {v['d5_gap_closed_by_crossed_projector']}")
    print(f"  CKM factor to strict: {v['ckm_factor_to_strict']:.6e}")
    print(f"  publication complete: {v['publication_complete']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
