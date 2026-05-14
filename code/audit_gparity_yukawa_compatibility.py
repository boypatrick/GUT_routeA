#!/usr/bin/env python3
"""Audit compatibility of the G120 triplet parity with the flavor sector.

No web lookup is used.  The triplet regulator is protected by an extra
Z2^G parity under which the G120-like triplet source is odd.  This script asks
whether that parity can act on the physical 120_H Yukawa multiplet, or whether
the protected triplet source must be a mediator-only copy.

The key mathematical test is a finite Z2 charge scan over family parities and
Higgs parities.  We require all symmetric 10_H and 126bar_H entries and all
three antisymmetric 120_H family entries to be allowed.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "gparity_yukawa_compatibility"

FLAVOR_CARD = ROOT / "output" / "flavor_fit" / "flavor_observable_audit.json"
SINGLE_120 = ROOT / "output" / "flavor_clebsch_scan" / "scan_clebsch_flavor_fit.json"
TWO_120 = ROOT / "output" / "flavor_clebsch_two120" / "scan_clebsch_flavor_two120.json"
SYMM_RELAX = ROOT / "output" / "flavor_symmetric_relaxation" / "scan_clebsch_flavor_symmetric_relaxation.json"
ALIGNMENT = ROOT / "output" / "z3_z2_alignment" / "z3_z2_alignment_summary.json"


def allowed_symmetric_entries(q16: tuple[int, int, int], qh: int) -> bool:
    return all((q16[i] + q16[j] + qh) % 2 == 0 for i in range(3) for j in range(i, 3))


def allowed_antisymmetric_entries(q16: tuple[int, int, int], qh: int) -> bool:
    return all((q16[i] + q16[j] + qh) % 2 == 0 for i in range(3) for j in range(i + 1, 3))


def scan_family_parities(force_120_odd: bool | None = None) -> list[dict[str, Any]]:
    rows = []
    for q16 in product([0, 1], repeat=3):
        for q10, q126, q120, q16bar in product([0, 1], repeat=4):
            if force_120_odd is not None and q120 != int(force_120_odd):
                continue
            sym10 = allowed_symmetric_entries(q16, q10)
            sym126 = allowed_symmetric_entries(q16, q126)
            anti120 = allowed_antisymmetric_entries(q16, q120)
            majorana = all((q16[i] + q16[j] + 2 * q16bar) % 2 == 0 for i in range(3) for j in range(i, 3))
            rows.append(
                {
                    "q16": q16,
                    "q10": q10,
                    "q126bar": q126,
                    "q120": q120,
                    "q16bar": q16bar,
                    "all_10_entries": sym10,
                    "all_126_entries": sym126,
                    "all_120_antisym_entries": anti120,
                    "all_majorana_entries": majorana,
                    "full_flavor_allowed": sym10 and sym126 and anti120 and majorana,
                }
            )
    return rows


def load_flavor_metrics() -> list[dict[str, Any]]:
    no120 = json.loads(FLAVOR_CARD.read_text(encoding="utf-8"))
    single = json.loads(SINGLE_120.read_text(encoding="utf-8"))["best"]
    two = json.loads(TWO_120.read_text(encoding="utf-8"))["best"]
    relax_payload = json.loads(SYMM_RELAX.read_text(encoding="utf-8"))
    relax = relax_payload["best"]

    def ckm_tuple_from_matrix(mat: list[list[float]]) -> tuple[float, float, float]:
        return float(mat[0][1]), float(mat[1][2]), float(mat[0][2])

    return [
        {
            "branch": "no_120_exact_Ominus4_card",
            "uses_physical_120": False,
            "Vus": no120["ckm_current_observables"]["abs_Vus"],
            "Vcb": no120["ckm_current_observables"]["abs_Vcb"],
            "Vub": no120["ckm_current_observables"]["abs_Vub"],
            "ckm_score": no120["ckm_magnitude_log_score"],
            "mass_score": sum(item["log_score"] for item in no120["mass_ratio_audit"].values()),
            "phenomenology_viable": False,
        },
        {
            "branch": "single_120_clebsch",
            "uses_physical_120": True,
            "Vus": ckm_tuple_from_matrix(single["CKM_abs"])[0],
            "Vcb": ckm_tuple_from_matrix(single["CKM_abs"])[1],
            "Vub": ckm_tuple_from_matrix(single["CKM_abs"])[2],
            "ckm_score": single["scores"]["ckm_magnitude_log_score"],
            "mass_score": single["scores"]["mass_log_score"],
            "phenomenology_viable": False,
        },
        {
            "branch": "two_120_clebsch",
            "uses_physical_120": True,
            "Vus": ckm_tuple_from_matrix(two["CKM_abs"])[0],
            "Vcb": ckm_tuple_from_matrix(two["CKM_abs"])[1],
            "Vub": ckm_tuple_from_matrix(two["CKM_abs"])[2],
            "ckm_score": two["scores"]["ckm_magnitude_log_score"],
            "mass_score": two["scores"]["mass_log_score"],
            "phenomenology_viable": False,
        },
        {
            "branch": "two_120_plus_symmetric_relaxation",
            "uses_physical_120": True,
            "Vus": ckm_tuple_from_matrix(relax["CKM_abs"])[0],
            "Vcb": ckm_tuple_from_matrix(relax["CKM_abs"])[1],
            "Vub": ckm_tuple_from_matrix(relax["CKM_abs"])[2],
            "ckm_score": relax["scores"]["ckm_magnitude_log_score"],
            "mass_score": relax["scores"]["mass_log_score"],
            "phenomenology_viable": relax_payload["phenomenology_viable"],
            "relative_symmetric_leakage": relax.get("epsilon_component_budget"),
        },
    ]


def operator_ledger() -> list[dict[str, Any]]:
    # Scenario A: the same physical 120_H is odd under Z2^G.
    # Scenario B: the physical 120_H is even, while the protected G_tr copy is odd.
    return [
        {
            "operator": "16_i 16_j 10_H",
            "role": "symmetric Yukawa H",
            "global_120_odd": "allowed",
            "mediator_only_Gtr": "allowed",
        },
        {
            "operator": "16_i 16_j 126bar_H",
            "role": "symmetric Yukawa F and/or Majorana channel",
            "global_120_odd": "allowed",
            "mediator_only_Gtr": "allowed",
        },
        {
            "operator": "16_i 16_j 120_H",
            "role": "antisymmetric Yukawa G needed by Clebsch flavor branch",
            "global_120_odd": "forbidden if all 16_i are even; impossible to allow all three antisymmetric entries with 120_H odd",
            "mediator_only_Gtr": "allowed",
        },
        {
            "operator": "T_Gtr barT_Gtr S1",
            "role": "desired Gtr-Gtr regulator lift",
            "global_120_odd": "allowed",
            "mediator_only_Gtr": "allowed",
        },
        {
            "operator": "single-Gtr dangerous entries",
            "role": "unused charged triplet leakage",
            "global_120_odd": "forbidden by G parity",
            "mediator_only_Gtr": "forbidden by G parity",
        },
        {
            "operator": "120_H physical - Gtr mediator mixing",
            "role": "would communicate parity to flavor sector",
            "global_120_odd": "not applicable: same field",
            "mediator_only_Gtr": "forbidden without an odd bridge spurion",
        },
    ]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scans_all = scan_family_parities()
    scans_odd = scan_family_parities(force_120_odd=True)
    scans_even = scan_family_parities(force_120_odd=False)
    full_all = [row for row in scans_all if row["full_flavor_allowed"]]
    full_odd = [row for row in scans_odd if row["full_flavor_allowed"]]
    full_even = [row for row in scans_even if row["full_flavor_allowed"]]
    metrics = load_flavor_metrics()
    alignment = json.loads(ALIGNMENT.read_text(encoding="utf-8"))
    ledger = operator_ledger()

    with (OUT / "gparity_family_parity_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "q16",
            "q10",
            "q126bar",
            "q120",
            "q16bar",
            "all_10_entries",
            "all_126_entries",
            "all_120_antisym_entries",
            "all_majorana_entries",
            "full_flavor_allowed",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in scans_all:
            out = dict(row)
            out["q16"] = " ".join(map(str, row["q16"]))
            writer.writerow(out)

    with (OUT / "gparity_operator_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ledger[0].keys()))
        writer.writeheader()
        writer.writerows(ledger)

    with (OUT / "gparity_flavor_branch_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = sorted({key for row in metrics for key in row.keys()})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics)

    no_go_proof = (
        "If q(120_H)=1 and all three antisymmetric 120 family entries are required, "
        "then q_i+q_j=1 for pairs 12,13,23.  Adding these three equations gives "
        "2(q_1+q_2+q_3)=3=1 mod 2, a contradiction.  Hence a Z2 that makes the physical "
        "120_H odd cannot preserve a full antisymmetric 120_H Yukawa tensor.  The scan "
        f"confirms full-flavor solutions with q(120_H)=1: {len(full_odd)}."
    )
    viable_statement = (
        "The viable selection-rule branch keeps the physical 120_H even and introduces a distinct "
        "mediator-only Gtr copy that is odd under Z2^G.  Then the previously found triplet regulator "
        "support is protected, while 16_i16_j120_H remains allowed.  Mixing between physical 120_H and "
        "Gtr is odd and therefore forbidden unless an explicit odd bridge spurion is introduced."
    )
    verdict = (
        "The G120 parity must not be applied to the full physical 120_H multiplet if the Clebsch flavor "
        "branch is retained.  Existing numerical flavor audits show that branches using physical 120_H "
        "directions are the only ones reaching small CKM angles, while the no-120 card has CKM score "
        f"{metrics[0]['ckm_score']:.6e}.  Therefore the symmetry-complete branch is a mediator-only "
        "Gtr parity, or equivalently a post-Spin(10) triplet-sector parity in the PS EFT.  The former is "
        "cleaner as an action-level bookkeeping device; the latter would need an explicit Spin(10)-breaking "
        "origin."
    )

    output = {
        "note": "No web lookup used. Compatibility audit of Z2^G triplet parity with Yukawa/Majorana sector.",
        "family_parity_scan": {
            "total_assignments": len(scans_all),
            "full_flavor_assignments": len(full_all),
            "full_flavor_assignments_with_physical_120_odd": len(full_odd),
            "full_flavor_assignments_with_physical_120_even": len(full_even),
            "example_even_solution": full_even[0] if full_even else None,
        },
        "no_go_proof": no_go_proof,
        "viable_branch": viable_statement,
        "operator_ledger": ledger,
        "flavor_branch_metrics": metrics,
        "triplet_alignment_reference": alignment["best_extra_grading"],
        "triplet_alignment_channels": alignment["channels"],
        "verdict": verdict,
    }
    (OUT / "gparity_yukawa_compatibility_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# G-parity Yukawa compatibility audit",
        "",
        "No web lookup was used.",
        "",
        "## Z2 family-parity theorem",
        "",
        no_go_proof,
        "",
        f"Full-flavor assignments with physical `120_H` odd: `{len(full_odd)}`.",
        f"Full-flavor assignments with physical `120_H` even: `{len(full_even)}`.",
        "",
        "## Flavor metrics",
        "",
        "| branch | uses 120 | Vus | Vcb | Vub | CKM score | mass score | viable |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics:
        report.append(
            f"| `{row['branch']}` | {int(row['uses_physical_120'])} | "
            f"{row['Vus']:.6e} | {row['Vcb']:.6e} | {row['Vub']:.6e} | "
            f"{row['ckm_score']:.6e} | {row['mass_score']:.6e} | {int(row['phenomenology_viable'])} |"
        )
    report.extend(
        [
            "",
            "## Operator ledger",
            "",
            "| operator | global physical 120 odd | mediator-only Gtr |",
            "|---|---|---|",
        ]
    )
    for row in ledger:
        report.append(f"| `{row['operator']}` | {row['global_120_odd']} | {row['mediator_only_Gtr']} |")
    report.extend(["", "## Verdict", "", verdict, ""])
    (OUT / "gparity_yukawa_compatibility_report.md").write_text("\n".join(report), encoding="utf-8")

    print("G-parity Yukawa compatibility audit")
    print(f"  full-flavor solutions with physical 120 odd: {len(full_odd)}")
    print(f"  full-flavor solutions with physical 120 even: {len(full_even)}")
    print(f"  no-120 CKM score: {metrics[0]['ckm_score']:.6e}")
    print(f"  best current viable branch uses 120: {metrics[-1]['uses_physical_120']} viable={metrics[-1]['phenomenology_viable']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
