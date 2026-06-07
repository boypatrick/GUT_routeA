#!/usr/bin/env python3
"""Route-C candidate ledger.

This is not an amplitude bootstrap yet.  It is the first reproducible ledger
for the representation-level assumptions that will be fed into the amplitude
bootstrap.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"


CANDIDATES = [
    {
        "name": "SU(5)",
        "simple": True,
        "single_irrep_family": False,
        "visible_family": "10 + bar5 + 1",
        "contains_nu_c": "requires singlet 1",
        "chiral_sm_exotics": False,
        "minimal_single_object": False,
        "route_c_status": "fails single-object criterion",
    },
    {
        "name": "Pati-Salam",
        "simple": False,
        "single_irrep_family": False,
        "visible_family": "(4,2,1) + (bar4,1,2)",
        "contains_nu_c": True,
        "chiral_sm_exotics": False,
        "minimal_single_object": False,
        "route_c_status": "natural face, but not simple",
    },
    {
        "name": "Spin(10)",
        "simple": True,
        "single_irrep_family": True,
        "visible_family": "16",
        "contains_nu_c": True,
        "chiral_sm_exotics": False,
        "minimal_single_object": True,
        "route_c_status": "standard minimal survivor",
    },
    {
        "name": "E6",
        "simple": True,
        "single_irrep_family": True,
        "visible_family": "27 -> 16 + 10 + 1 under Spin(10)",
        "contains_nu_c": True,
        "chiral_sm_exotics": "requires lifting extra states",
        "minimal_single_object": False,
        "route_c_status": "allowed only with extra-state audit",
    },
]


def survives_minimal_single_object(candidate: dict) -> bool:
    return bool(
        candidate["simple"]
        and candidate["single_irrep_family"]
        and candidate["contains_nu_c"] is True
        and candidate["chiral_sm_exotics"] is False
        and candidate["minimal_single_object"]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for cand in CANDIDATES:
        row = dict(cand)
        row["survives_minimal_single_object_filter"] = (
            survives_minimal_single_object(cand)
        )
        rows.append(row)

    survivors = [r["name"] for r in rows if r["survives_minimal_single_object_filter"]]
    payload = {
        "description": "Route-C representation-level candidate pre-ledger",
        "filter": [
            "simple group",
            "one irreducible family object",
            "contains nu^c",
            "no chiral SM exotics",
            "minimal single-object realization",
        ],
        "survivors": survivors,
        "candidates": rows,
    }

    (OUT / "candidate_ledger.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    md = [
        "# Route-C Candidate Ledger",
        "",
        "This is a pre-amplitude representation ledger.  It does not prove the",
        "GUT; it records which candidates survive the minimal single-object",
        "filter before any pole-factorization or Ward-identity bootstrap.",
        "",
        "| candidate | visible family | status | survivor |",
        "| --- | --- | --- | --- |",
    ]
    for row in rows:
        md.append(
            "| {name} | `{visible_family}` | {route_c_status} | {survivor} |".format(
                name=row["name"],
                visible_family=row["visible_family"],
                route_c_status=row["route_c_status"],
                survivor="yes"
                if row["survives_minimal_single_object_filter"]
                else "no",
            )
        )
    md.extend(
        [
            "",
            "Survivors under the stated filter:",
            "",
            "```text",
            ", ".join(survivors) if survivors else "none",
            "```",
            "",
            "Next audit: attach symbolic pole residues and test factorization,",
            "Ward identities, anomaly constraints, and high-energy behavior.",
            "",
        ]
    )
    (OUT / "candidate_ledger.md").write_text("\n".join(md), encoding="utf-8")

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

