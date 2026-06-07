#!/usr/bin/env python3
"""Route-C P5 anomaly and chiral-consistency ledger.

P5 compares the Route-C candidate set against anomaly and chiral-spectrum
requirements, using the P1 SM + nu^c family ledger as the low-energy baseline.

This is still not a full amplitude bootstrap.  It records which candidates are
anomaly-consistent, which need vectorlike lifting audits, and which fail the
Route-C minimal single-object filter for reasons other than anomalies.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P1_JSON = OUT / "external_state_charges.json"


def load_p1() -> dict:
    if not P1_JSON.exists():
        raise SystemExit(f"missing {P1_JSON}; run external_state_charge_table.py")
    return json.loads(P1_JSON.read_text(encoding="utf-8"))


def baseline_summary(p1: dict) -> dict:
    anomaly_passes = {
        name: row["passes"] for name, row in p1["checks"]["anomalies"].items()
    }
    return {
        "spectrum": "Q + L + u^c + d^c + nu^c + e^c",
        "hypercharge_normalization": p1["checks"]["hypercharge_normalization"],
        "all_sm_family_anomalies_pass": all(anomaly_passes.values())
        and p1["checks"]["witten_SU2"]["passes_even_doublet_count"],
        "anomaly_passes": anomaly_passes,
        "witten_SU2": p1["checks"]["witten_SU2"],
    }


def candidate_rows() -> list[dict]:
    return [
        {
            "candidate": "SU(5)",
            "category": "four-dimensional field-theory GUT candidate",
            "gauge_group_type": "simple",
            "family_realization": "10 + bar5 + 1",
            "contains_nu_c": "yes, but as an SU(5) singlet",
            "single_irreducible_family_object": False,
            "chiral_exotics_status": "none for one SM family plus singlet",
            "anomaly_checks": [
                {
                    "name": "SU(5)^3 perturbative anomaly",
                    "formula": "A(10)+A(bar5)+A(1)=1-1+0=0",
                    "passes": True,
                },
                {
                    "name": "global SU(5) anomaly",
                    "formula": "pi_4(SU(5))=0",
                    "passes": True,
                },
            ],
            "route_c_minimal_single_object_filter": False,
            "route_c_status": (
                "anomaly-consistent but fails the single irreducible family "
                "object criterion"
            ),
            "missing_audits": [
                "embedding comparison if the single-object criterion is relaxed",
                "proton and threshold audits for an SU(5)-specific branch",
            ],
        },
        {
            "candidate": "Pati-Salam",
            "category": "four-dimensional field-theory GUT-like candidate",
            "gauge_group_type": "product",
            "family_realization": "(4,2,1) + (bar4,1,2)",
            "contains_nu_c": "yes, inside the SU(2)_R doublet",
            "single_irreducible_family_object": False,
            "chiral_exotics_status": "none for the standard family pair",
            "anomaly_checks": [
                {
                    "name": "SU(4)_C^3 perturbative anomaly",
                    "formula": "2 A(4)+2 A(bar4)=2-2=0",
                    "passes": True,
                },
                {
                    "name": "SU(2)_L Witten anomaly",
                    "formula": "four SU(2)_L doublets from dim(4)=4, even",
                    "passes": True,
                },
                {
                    "name": "SU(2)_R Witten anomaly",
                    "formula": "four SU(2)_R doublets from dim(bar4)=4, even",
                    "passes": True,
                },
            ],
            "route_c_minimal_single_object_filter": False,
            "route_c_status": (
                "anomaly-consistent and physically natural, but fails the "
                "simple-group and single-object filters"
            ),
            "missing_audits": [
                "explicit product-group bootstrap branch if simplicity is relaxed",
                "broken-sector Ward and proton audits for Pati-Salam leptoquark currents",
            ],
        },
        {
            "candidate": "Spin(10)",
            "category": "four-dimensional field-theory GUT candidate",
            "gauge_group_type": "simple",
            "family_realization": "16 half-spinor",
            "contains_nu_c": "yes, inside the half-spinor",
            "single_irreducible_family_object": True,
            "chiral_exotics_status": "none in the standard 16",
            "anomaly_checks": [
                {
                    "name": "Spin(10)^3 perturbative anomaly",
                    "formula": "D5 has no independent cubic gauge-anomaly invariant; 16 is anomaly-safe",
                    "passes": True,
                },
                {
                    "name": "Spin(10) global anomaly",
                    "formula": "pi_4(Spin(10))=0 in the standard four-dimensional check",
                    "passes": True,
                },
            ],
            "route_c_minimal_single_object_filter": True,
            "route_c_status": (
                "standard minimal anomaly-safe single-object survivor under "
                "the stated Route-C pre-filter"
            ),
            "missing_audits": [
                "explicit broken-generator matrices for Ward checks",
                "Higgs/source sector for broken-vector completion",
                "full low-energy proton and threshold matching",
            ],
        },
        {
            "candidate": "E6",
            "category": "four-dimensional field-theory GUT candidate",
            "gauge_group_type": "simple",
            "family_realization": "27 -> 16 + 10 + 1 under Spin(10)",
            "contains_nu_c": "yes, through the embedded Spin(10) 16",
            "single_irreducible_family_object": True,
            "chiral_exotics_status": (
                "extra 10 + 1 must be lifted or otherwise accounted for"
            ),
            "anomaly_checks": [
                {
                    "name": "E6^3 perturbative anomaly",
                    "formula": "E6 has no independent cubic gauge-anomaly invariant for the 27",
                    "passes": True,
                },
                {
                    "name": "extra-state anomaly pairing after restriction",
                    "formula": "27 contains 16 plus vectorlike 10 and singlet under Spin(10)",
                    "passes": "conditional on a lifting/vectorlike audit",
                },
            ],
            "route_c_minimal_single_object_filter": False,
            "route_c_status": (
                "anomaly-consistent, but not minimal and requires an extra-state "
                "lifting audit"
            ),
            "missing_audits": [
                "mass/lifting mechanism for 10 + 1 extra states",
                "check that exotics do not reintroduce low-energy chiral matter",
                "E6-specific threshold and proton audits",
            ],
        },
        {
            "candidate": "Superstring-derived completions",
            "category": "UV-completion class, not a single four-dimensional group",
            "gauge_group_type": "compactification-dependent",
            "family_realization": (
                "model-dependent; may realize SU(5), Pati-Salam, Spin(10), E6, "
                "or other quiver/brane/heterotic branches"
            ),
            "contains_nu_c": "compactification-dependent",
            "single_irreducible_family_object": "not a universal requirement",
            "chiral_exotics_status": (
                "compactification-dependent; exotics must be absent, vectorlike, "
                "or lifted"
            ),
            "anomaly_checks": [
                {
                    "name": "higher-dimensional anomaly cancellation",
                    "formula": "Green-Schwarz/tadpole/modular consistency as appropriate to the string construction",
                    "passes": "construction-dependent",
                },
                {
                    "name": "four-dimensional anomaly cancellation",
                    "formula": "massless spectrum plus any 4D Green-Schwarz terms must cancel anomalies",
                    "passes": "construction-dependent",
                },
                {
                    "name": "worldsheet or compactification consistency",
                    "formula": "modular invariance, tadpole cancellation, K-theory/Freed-Witten-type constraints where applicable",
                    "passes": "construction-dependent",
                },
            ],
            "route_c_minimal_single_object_filter": "not directly comparable",
            "route_c_status": (
                "UV-completion class with stronger high-energy softness/tower "
                "structure, but requires compactification-specific spectrum and "
                "exotics audits"
            ),
            "missing_audits": [
                "explicit compactification and massless spectrum",
                "modular/tadpole/Green-Schwarz consistency data",
                "chiral index and exotics-lifting audit",
                "matching to the Route-A representation/index skeleton",
            ],
        },
    ]


def candidate_passes_anomaly(candidate: dict) -> bool | str:
    values = [check["passes"] for check in candidate["anomaly_checks"]]
    if all(value is True for value in values):
        return True
    if any(value is False for value in values):
        return False
    return "conditional"


def build_markdown(payload: dict) -> str:
    baseline = payload["baseline"]
    lines = [
        "# Route-C P5 Anomaly and Chiral-Consistency Ledger",
        "",
        "This ledger compares candidate GUT branches against anomaly and",
        "chiral-spectrum requirements.  It uses the P1 one-family",
        "`SM + nu^c` ledger as the low-energy baseline.",
        "",
        "## Baseline",
        "",
        "The P1 baseline spectrum is",
        "",
        "```text",
        baseline["spectrum"],
        "```",
        "",
        "and the exact checks pass:",
        "",
        "```text",
        f"Tr Y^2 / Tr T3L^2 = {baseline['hypercharge_normalization']['ratio']['fraction']}",
        f"all SM-family anomalies pass = {baseline['all_sm_family_anomalies_pass']}",
        "```",
        "",
        "## Candidate Summary",
        "",
        "| candidate | anomaly status | single object | exotics/lifting status | Route-C status |",
        "| --- | --- | --- | --- | --- |",
    ]
    for candidate in payload["candidates"]:
        lines.append(
            "| {name} | {anom} | {single} | {exotics} | {status} |".format(
                name=candidate["candidate"],
                anom=candidate["anomaly_status"],
                single=candidate["single_irreducible_family_object"],
                exotics=candidate["chiral_exotics_status"],
                status=candidate["route_c_status"],
            )
        )

    lines.extend(
        [
            "",
            "## Detailed Checks",
            "",
        ]
    )
    for candidate in payload["candidates"]:
        lines.extend(
            [
                f"### {candidate['candidate']}",
                "",
                f"- category: {candidate['category']}",
                f"- family realization: `{candidate['family_realization']}`",
                f"- contains `nu^c`: {candidate['contains_nu_c']}",
                f"- minimal single-object filter: {candidate['route_c_minimal_single_object_filter']}",
                "",
                "| check | formula | passes |",
                "| --- | --- | --- |",
            ]
        )
        for check in candidate["anomaly_checks"]:
            lines.append(
                "| {name} | `{formula}` | {passes} |".format(
                    name=check["name"],
                    formula=check["formula"],
                    passes=check["passes"],
                )
            )
        lines.extend(
            [
                "",
                "Missing audits:",
                "",
            ]
        )
        for audit in candidate["missing_audits"]:
            lines.append(f"- {audit}")
        lines.append("")

    lines.extend(
        [
            "## P5 Boundary",
            "",
            "Passing anomaly checks does not by itself select a GUT.  `SU(5)`,",
            "Pati-Salam, `Spin(10)`, and `E6` all have anomaly-consistent forms,",
            "but they satisfy different structural filters.  Superstring-derived",
            "completions are a UV-completion class rather than one fixed",
            "four-dimensional candidate; they require compactification-specific",
            "modular/tadpole/Green-Schwarz, spectrum, and exotics-lifting audits.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    p1 = load_p1()
    candidates = []
    for row in candidate_rows():
        candidate = dict(row)
        candidate["anomaly_status"] = candidate_passes_anomaly(candidate)
        candidates.append(candidate)

    payload = {
        "description": "Route-C P5 anomaly and chiral-consistency ledger",
        "baseline": baseline_summary(p1),
        "candidates": candidates,
        "summary": {
            "candidates_inspected": len(candidates),
            "field_theory_candidates": sum(
                1
                for row in candidates
                if row["category"].startswith("four-dimensional")
            ),
            "uv_completion_classes": sum(
                1 for row in candidates if "UV-completion" in row["category"]
            ),
            "route_c_minimal_single_object_survivors": [
                row["candidate"]
                for row in candidates
                if row["route_c_minimal_single_object_filter"] is True
            ],
            "baseline_all_anomalies_pass": baseline_summary(p1)[
                "all_sm_family_anomalies_pass"
            ],
        },
        "boundary": (
            "Anomaly consistency is necessary but not sufficient.  The ledger "
            "does not perform Ward, high-energy, proton, threshold, or compactification audits."
        ),
    }

    (OUT / "anomaly_chiral_consistency.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "anomaly_chiral_consistency.md").write_text(
        build_markdown(payload), encoding="utf-8"
    )

    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
