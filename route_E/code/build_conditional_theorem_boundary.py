#!/usr/bin/env python3
"""Write the theorem-boundary certificate for the current GUT draft.

This is a logical closure audit.  It does not introduce a new hidden-sector
mechanism.  Instead it records the strongest theorem statement supported by
the local no-web evidence ledger and separates conditional assumptions from
claims that have been tested and rejected.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "conditional_theorem_boundary"


def read_json(relpath: str) -> dict[str, Any]:
    return json.loads((ROOT / relpath).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_claim_rows() -> list[dict[str, str]]:
    path = ROOT / "output" / "conditional_theorem_ledger" / "claim_ledger.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def status_count(rows: list[dict[str, str]], status: str) -> int:
    return sum(row["status"] == status for row in rows)


def write_report(summary: dict[str, Any]) -> str:
    lines = [
        "# Conditional theorem boundary certificate",
        "",
        "No web lookup was used.  This certificate records the strongest theorem",
        "statement supported by the local evidence ledger.",
        "",
        "## Boundary verdict",
        "",
        f"- unconditional first-principles GUT derived: `{summary['boundary_verdict']['unconditional_first_principles_gut_derived']}`",
        f"- conditional Spin(10) EFT branch verified locally: `{summary['boundary_verdict']['conditional_spin10_eft_branch_verified_locally']}`",
        f"- zeta is a current conditional datum: `{summary['boundary_verdict']['majorana_zeta_is_current_conditional_datum']}`",
        "",
        "## Numerical keys",
        "",
    ]
    for key, value in summary["numerical_keys"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(
        [
            "",
            "## Theorem boundary",
            "",
            summary["theorem_boundary_statement"],
            "",
            "## Conditional assumptions retained",
            "",
        ]
    )
    for item in summary["conditional_assumptions_retained"]:
        lines.append(f"- {item}")
    lines.extend(["", "## No-go conclusions", ""])
    for item in summary["no_go_conclusions"]:
        lines.append(f"- {item}")
    lines.extend(["", "## Open publication tasks", ""])
    for item in summary["open_publication_tasks"]:
        lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    ledger_summary = read_json("output/conditional_theorem_ledger/summary.json")
    majorana_contact = read_json("output/majorana_contact_sensitivity/summary.json")
    majorana_rank = read_json("output/source_majorana_texture_rank/summary.json")
    majorana_hidden = read_json("output/majorana_hidden_quotient_origin/summary.json")
    majorana_monomial = read_json("output/majorana_monomial_clockwork_origin/summary.json")
    provenance = read_json("output/no_web_table_provenance_manifest/summary.json")
    rows = read_claim_rows()

    row_count = len(rows)
    no_go_count = status_count(rows, "NO_GO")
    open_count = status_count(rows, "OPEN")
    tuned_count = status_count(rows, "TUNED_FALLBACK")
    pass_conditional_count = status_count(rows, "PASS_CONDITIONAL")

    evidence_paths = [ROOT / row["evidence_path"] for row in rows]
    missing_evidence = [str(path.relative_to(ROOT)) for path in evidence_paths if not path.exists()]

    zeta_contact_fraction = majorana_rank["decomposition"]["contact_fraction"]
    zeta_scale_width = majorana_contact["real_scale_loose_interval"]["half_width"]
    zeta_phase_width = majorana_contact["phase_loose_interval"]["half_width"]
    best_small = majorana_monomial["best_small_denominator_row"]
    first_large = majorana_monomial["first_denominator_within_loose_phase_tolerance"]

    boundary_verdict = {
        "unconditional_first_principles_gut_derived": False,
        "conditional_spin10_eft_branch_verified_locally": (
            not missing_evidence
            and provenance["verdict"]["provenance_manifest_complete"]
            and ledger_summary["counts"].get("PASS_CONDITIONAL", 0) >= 50
        ),
        "majorana_zeta_is_current_conditional_datum": (
            zeta_contact_fraction > 0.1
            and zeta_scale_width < 1.0e-4
            and zeta_phase_width < 1.0e-4
            and not majorana_hidden["verdict"]["pure_quotient_predicts_zeta"]
            and not majorana_monomial["verdict"]["small_monomial_clockwork_predicts_phase"]
        ),
        "claim_ledger_has_missing_evidence": bool(missing_evidence),
        "no_web_provenance_complete": provenance["verdict"]["provenance_manifest_complete"],
    }

    summary = {
        "note": "No web lookup used. Logical theorem-boundary certificate.",
        "input_manifest": [
            {
                "label": "conditional_theorem_ledger",
                "path": "output/conditional_theorem_ledger/summary.json",
                "sha256": sha256(ROOT / "output" / "conditional_theorem_ledger" / "summary.json"),
            },
            {
                "label": "claim_ledger",
                "path": "output/conditional_theorem_ledger/claim_ledger.csv",
                "sha256": sha256(ROOT / "output" / "conditional_theorem_ledger" / "claim_ledger.csv"),
            },
        ],
        "ledger_counts": ledger_summary["counts"],
        "row_count": row_count,
        "missing_evidence": missing_evidence,
        "boundary_verdict": boundary_verdict,
        "numerical_keys": {
            "claim_rows": row_count,
            "no_go_count": no_go_count,
            "open_count": open_count,
            "tuned_fallback_count": tuned_count,
            "pass_conditional_count": pass_conditional_count,
            "zeta_contact_fraction": f"{zeta_contact_fraction:.6e}",
            "zeta_scale_loose_half_width": f"{zeta_scale_width:.6e}",
            "zeta_phase_loose_half_width_rad": f"{zeta_phase_width:.6e}",
            "monomial_best_small_N": best_small["denominator"],
            "monomial_best_small_phase_residual_rad": f"{best_small['phase_residual_rad']:.6e}",
            "monomial_first_acceptable_N": first_large["denominator"],
        },
        "theorem_boundary_statement": (
            "The local evidence supports the implication: conditional Spin(10) EFT "
            "+ CP1/O(2) family kinematics + constrained/composite 54/210 threshold "
            "sector + source-fixed Majorana contact zeta + clockwork-rescued local "
            "d=5 proxy implies the current internally checked benchmark branch.  "
            "It does not support the stronger implication PSLT/Another Physics "
            "alone implies a unique GUT or a predicted Majorana zeta."
        ),
        "conditional_assumptions_retained": [
            "Spin(10)/Pati-Salam EFT branch rather than PSLT-only microscopic derivation.",
            "CP1/O(2) protected three-family kinematic space.",
            "Constrained/composite 54/210 and auxiliary/conormal source sectors for thresholds.",
            "The PMNS-sensitive Majorana contact coefficient zeta as a source datum.",
            "Local no-web flavor/proton conventions pending external cited input refresh.",
        ],
        "no_go_conclusions": [
            "PSLT/Another Physics alone has not derived a unique GUT.",
            "Pure Veronese Majorana texture is insufficient for the source-consistent PMNS target.",
            "Minimal hidden D/product quotient does not determine zeta.",
            "Small monomial/clockwork root phases do not determine zeta within PMNS tolerance.",
        ],
        "open_publication_tasks": ledger_summary["open_blockers"],
        "verdict": (
            "Stop treating the current construction as an unconditional derivation.  "
            "It is a locally verified conditional EFT branch whose dominant remaining "
            "first-principles obstruction is the microscopic origin of zeta."
        ),
    }

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "report.md").write_text(write_report(summary), encoding="utf-8")
    print("Conditional theorem boundary written")
    print(f"  rows: {row_count}")
    print(f"  no-go: {no_go_count}")
    print(f"  open: {open_count}")
    print(f"  zeta datum: {boundary_verdict['majorana_zeta_is_current_conditional_datum']}")
    print(f"  conditional branch: {boundary_verdict['conditional_spin10_eft_branch_verified_locally']}")


if __name__ == "__main__":
    main()
