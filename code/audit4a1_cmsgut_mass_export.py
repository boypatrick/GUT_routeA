#!/usr/bin/env python3
"""Build Audit 4a.1 CMSGUT mass-export stage-2 schema.

This is the next gate after the vacuum-cubic validator.  It performs the
group-theoretic Goldstone count for the generic CMSGUT branch and exports the
schema for the doublet/triplet mass blocks needed before threshold and
dimension-five proton-decay audits can run.  It deliberately does not invent
the full Aulakh--Girdhar/BMSV mass matrices; those entries must be imported
from the literature convention table in a later audit.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit4a1"

SOURCES = {
    "audit4a1_mass_export_builder": ROOT / "code" / "audit4a1_cmsgut_mass_export.py",
    "audit4a1_vacuum_card": ROOT / "output" / "audit4a1" / "vacuum_branches.json",
    "audit4a_schema": ROOT / "output" / "audit4a" / "source_spectrum_schema.json",
    "audit0_card": ROOT / "output" / "audit0" / "invariant_card.json",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    encoded = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def manifest() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label, path in SOURCES.items():
        rows.append(
            {
                "label": label,
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path),
            }
        )
    return rows


def goldstone_sectors() -> list[dict[str, Any]]:
    return [
        {
            "sector_id": "SU4C_over_SU3C_U1BL_offdiagonal",
            "source_ps_block": "(15,1,1) subset of adjoint 45",
            "sm_representations": ["(3,1,2/3)", "(bar3,1,-2/3)"],
            "real_generator_count": 6,
            "interpretation": "color-B-L off-diagonal broken generators",
        },
        {
            "sector_id": "SU2R_charged",
            "source_ps_block": "(1,1,3) subset of adjoint 45",
            "sm_representations": ["(1,1,1)", "(1,1,-1)"],
            "real_generator_count": 2,
            "interpretation": "charged right-handed weak generators",
        },
        {
            "sector_id": "orthogonal_U1R_BL",
            "source_ps_block": "(15,1,1)+(1,1,3) Cartan mixture",
            "sm_representations": ["(1,1,0)_orthogonal_to_Y"],
            "real_generator_count": 1,
            "interpretation": "neutral Cartan direction orthogonal to hypercharge",
        },
        {
            "sector_id": "PS_broken_bidoublet",
            "source_ps_block": "(6,2,2) subset of adjoint 45",
            "sm_representations": [
                "(3,2,1/6)",
                "(bar3,2,-1/6)",
                "(3,2,-5/6)",
                "(bar3,2,5/6)",
            ],
            "real_generator_count": 24,
            "interpretation": "leptoquark-like bidoublet broken generators",
        },
    ]


def goldstone_count_card() -> dict[str, Any]:
    sectors = goldstone_sectors()
    unbroken = {
        "SU3C": 8,
        "SU2L": 3,
        "U1Y": 1,
    }
    spin10_dim = 45
    unbroken_dim = sum(unbroken.values())
    broken_from_difference = spin10_dim - unbroken_dim
    broken_from_sectors = sum(row["real_generator_count"] for row in sectors)
    return {
        "parent_group": "Spin(10)",
        "parent_dimension": spin10_dim,
        "unbroken_group": "SU(3)_C x SU(2)_L x U(1)_Y",
        "unbroken_dimension_by_factor": unbroken,
        "unbroken_dimension": unbroken_dim,
        "broken_dimension_from_difference": broken_from_difference,
        "broken_dimension_from_sector_sum": broken_from_sectors,
        "sectors": sectors,
        "pass_fail": broken_from_difference == 33 and broken_from_sectors == 33,
        "claim_boundary": "group-theoretic eaten-generator count only; no scalar Hessian or chiral Goldstone mass matrix is exported here",
    }


def matrix_placeholder(name: str, shape: tuple[int, int], row_basis: list[str], col_basis: list[str]) -> dict[str, Any]:
    rows, cols = shape
    entries = [[f"{name}_{i+1}{j+1}" for j in range(cols)] for i in range(rows)]
    return {
        "block_id": name,
        "shape": [rows, cols],
        "row_basis": row_basis,
        "column_basis": col_basis,
        "symbolic_entries": entries,
        "entry_status": "schema_placeholder_pending_literature_import",
        "normalization_convention": "Aulakh-Girdhar vacuum convention fixed by Audit 4a.1; matrix-entry signs and field normalizations still require dual-source import",
    }


def mass_block_schema() -> dict[str, Any]:
    doublet_basis = [
        "D_10_from_(1,2,2)",
        "D_126_from_(15,2,2)",
        "D_overline126_from_(15,2,2)",
        "D_210_from_(10+overline10_or_equivalent_doublet_slot)",
    ]
    triplet_basis = [
        "T_10",
        "T_126",
        "T_overline126",
        "T_210_slot_A",
        "T_210_slot_B",
    ]
    triplet_matrix = matrix_placeholder("M_T", (5, 5), triplet_basis, [b.replace("T_", "Tbar_") for b in triplet_basis])
    return {
        "doublet_mass_matrix": matrix_placeholder("M_D", (4, 4), doublet_basis, [b.replace("D_", "Dbar_") for b in doublet_basis]),
        "triplet_mass_matrix": triplet_matrix,
        "triplet_inverse_block": {
            "block_id": "M_T_inverse",
            "shape": triplet_matrix["shape"],
            "basis_inherited_from": "triplet_mass_matrix",
            "required_for": "Audit 2 dimension-five Wilson tensors",
            "entry_status": "blocked_until_triplet_mass_entries_imported_and_inverted",
        },
        "heavy_spectrum_export": {
            "status": "not_exported",
            "reason": "full symbolic mass matrices and SM beta-vector decomposition are not yet imported",
        },
    }


def stage_gates(goldstone: dict[str, Any]) -> dict[str, Any]:
    return {
        "goldstone_count_33_exported": bool(goldstone["pass_fail"]),
        "symbolic_mass_matrix_schema_exported": True,
        "doublet_mass_matrix_schema_exported": True,
        "triplet_mass_matrix_schema_exported": True,
        "triplet_inverse_block_schema_exported": True,
        "doublet_mass_matrix_entries_imported": False,
        "triplet_mass_matrix_entries_imported": False,
        "triplet_inverse_block_entries_exported": False,
        "heavy_spectrum_json_nonplaceholder": False,
        "threshold_beta_vectors_exported": False,
    }


def write_markdown(card: dict[str, Any]) -> None:
    goldstone = card["goldstone_count"]
    lines = [
        "# Audit 4a.1 CMSGUT Mass-Export Stage-2 Schema",
        "",
        "This artifact starts the mass-export layer after the vacuum-cubic",
        "validator.  It proves only the group-theoretic Goldstone count and",
        "exports the doublet/triplet block schema; it does not import the full",
        "CMSGUT mass matrices.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- upstream vacuum sha256: `{card['upstream_audit4a1_vacuum_sha256']}`",
        f"- upstream Audit 4a sha256: `{card['upstream_audit4a_sha256']}`",
        "",
        "## Goldstone Count",
        "",
        f"- Parent dimension: `{goldstone['parent_dimension']}`.",
        f"- Unbroken dimension: `{goldstone['unbroken_dimension']}`.",
        f"- Broken dimension: `{goldstone['broken_dimension_from_difference']}`.",
        f"- Sector sum: `{goldstone['broken_dimension_from_sector_sum']}`.",
        f"- Pass: `{goldstone['pass_fail']}`.",
        "",
        "| sector | SM reps | count |",
        "| --- | --- | --- |",
    ]
    for row in goldstone["sectors"]:
        reps = ", ".join(row["sm_representations"])
        lines.append(f"| `{row['sector_id']}` | {reps} | {row['real_generator_count']} |")
    lines += [
        "",
        "## Mass-Block Export Status",
        "",
        "| block | shape | status |",
        "| --- | --- | --- |",
    ]
    blocks = card["mass_block_schema"]
    for key in ["doublet_mass_matrix", "triplet_mass_matrix"]:
        block = blocks[key]
        lines.append(f"| `{key}` | `{block['shape']}` | {block['entry_status']} |")
    inv = blocks["triplet_inverse_block"]
    lines.append(f"| `triplet_inverse_block` | `{inv['shape']}` | {inv['entry_status']} |")
    lines += [
        "",
        "## Boundary",
        "",
        "Audit 4a.1 now exports the Goldstone-count gate and the mass-block schema.",
        "It still does not export CMSGUT mass entries, heavy beta vectors, a",
        "non-placeholder heavy spectrum, or a dimension-five Wilson input.",
        "",
    ]
    (OUT / "mass_export_schema.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit0 = read_json(SOURCES["audit0_card"])
    audit4a = read_json(SOURCES["audit4a_schema"])
    vacuum = read_json(SOURCES["audit4a1_vacuum_card"])

    goldstone = goldstone_count_card()
    card: dict[str, Any] = {
        "audit": "audit4a1_cmsgut_mass_export",
        "status": "stage-2 schema: Goldstone count plus doublet/triplet block contract; no full CMSGUT mass entries",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "upstream_audit0_sha256": audit0["card_sha256"],
        "upstream_audit4a_sha256": audit4a["card_sha256"],
        "upstream_audit4a1_vacuum_sha256": vacuum["card_sha256"],
        "vacuum_convention_inherited": vacuum["literature_convention_map"],
        "goldstone_count": goldstone,
        "mass_block_schema": mass_block_schema(),
        "stage_gates": stage_gates(goldstone),
        "publication_boundary": {
            "claim": "The generic Spin(10)->G_SM branch has 33 broken generators, and the Audit-2-relevant mass-block schema is fixed.",
            "not_claimed": [
                "scalar Hessian Goldstone eigenvectors",
                "Aulakh-Girdhar/BMSV mass-matrix entries",
                "doublet fine-tuning solution",
                "triplet inverse propagator entries",
                "non-placeholder heavy_spectrum.json",
                "threshold closure",
                "dimension-five proton safety",
            ],
        },
        "next_required_steps": [
            "Import doublet and triplet mass-matrix entries in the inherited Aulakh-Girdhar convention.",
            "Verify scalar/chiral Goldstone directions against the mass matrices, not only generator dimensions.",
            "Invert the triplet block and export the entries needed by Audit 2.",
            "Decompose all heavy states into SM beta-vector rows for Audit 3.",
        ],
    }
    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "mass_export_schema.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(card)

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {(OUT / 'mass_export_schema.md').relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")


if __name__ == "__main__":
    main()
