#!/usr/bin/env python3
"""Build Audit 4a source-sector / heavy-spectrum schema card.

Audit 4a is the field-theory source-sector interface needed before threshold
and proton-decay audits can become paper-grade.  This script records the schema,
existing local evidence, acceptance gates, and unresolved source-spectrum
blockers.  It does not claim a completed UV model or a complete heavy spectrum.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit4a"

SOURCES = {
    "audit4a_builder": ROOT / "code" / "audit4a_source_spectrum_schema.py",
    "audit0_card": ROOT / "output" / "audit0" / "invariant_card.json",
    "constrained_54_210": ROOT / "output" / "constrained_54_210_source_sector" / "summary.json",
    "combined_conormal_54_210": ROOT / "output" / "combined_conormal_54_210" / "summary.json",
    "drive_sector_spectrum": ROOT / "output" / "drive_sector_spectrum" / "drive_sector_spectrum_summary.json",
    "completed_120_partner_action": ROOT / "output" / "completed_120_partner_action" / "summary.json",
    "no_web_input_convention_ledger": ROOT / "output" / "no_web_input_convention_ledger" / "summary.json",
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


def heavy_state_schema() -> dict[str, Any]:
    required = [
        "state_id",
        "sector",
        "spin10_representation",
        "ps_representation",
        "sm_representation",
        "multiplicity",
        "mass_GeV",
        "mass_expression",
        "threshold_beta_vector_b1_b2_b3",
        "d5_role",
        "source_artifact",
        "status",
    ]
    optional = [
        "mixing_block_id",
        "triplet_inverse_propagator_entry",
        "doublet_triplet_partner",
        "is_complete_gut_multiplet",
        "projected_threshold_vector",
        "notes",
    ]
    return {
        "required_fields": required,
        "optional_fields": optional,
        "status_values": [
            "derived_mass",
            "benchmark_mass",
            "schema_placeholder",
            "excluded_or_auxiliary",
        ],
        "minimum_output_for_Audit3": [
            "all non-SM heavy fields with SM representation and mass",
            "one-loop beta vector for each split threshold",
            "common matching-scale convention",
        ],
        "minimum_output_for_Audit2": [
            "colored-triplet mass matrix or inverse propagator",
            "triplet Yukawa tensor map in the Audit 1 physical basis",
            "doublet-triplet selection or filtering rule",
        ],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    audit0 = read_json(SOURCES["audit0_card"])
    constrained = read_json(SOURCES["constrained_54_210"])
    combined = read_json(SOURCES["combined_conormal_54_210"])
    drive = read_json(SOURCES["drive_sector_spectrum"])
    partner = read_json(SOURCES["completed_120_partner_action"])
    no_web = read_json(SOURCES["no_web_input_convention_ledger"])

    drive_fields = drive.get("drive_sector_fields", [])
    drive_beta_l2 = [row.get("projected_beta_l2") for row in drive_fields if row.get("projected_beta_l2") is not None]

    card: dict[str, Any] = {
        "audit": "audit4a_source_spectrum_schema",
        "status": "source-sector/heavy-spectrum interface scaffold only; not a completed UV model",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "audit0_card_sha256": audit0["card_sha256"],
        "candidate_tracks": {
            "field_theory_constrained_54_210": {
                "role": "near-term Audit 4a track",
                "existing_local_status": "has constrained-orbit and conormal evidence, but not a full heavy spectrum",
                "primary_artifacts": [
                    "output/constrained_54_210_source_sector/summary.json",
                    "output/combined_conormal_54_210/summary.json",
                    "output/drive_sector_spectrum/drive_sector_spectrum_summary.json",
                ],
            },
            "cmsgut_like_210_126_10": {
                "role": "literature-compatible fallback track",
                "existing_local_status": "schema placeholder; requires explicit superpotential and spectrum import",
                "reason_to_keep": "can supply P_nu^c via a 126bar_H source and complete heavy multiplet bookkeeping",
            },
            "routeD_global_geometry": {
                "role": "optional parallel 4b track, not the near-term Audit 4a default",
                "existing_local_status": "local string-placement appendix only; no global compactification spectrum",
            },
        },
        "heavy_spectrum_schema": heavy_state_schema(),
        "existing_local_evidence": {
            "constrained_54_210_passes": constrained.get("passes"),
            "audit_54": constrained.get("audit_54"),
            "audit_210": constrained.get("audit_210"),
            "projector_checks": constrained.get("projector_checks"),
            "uv_rows": constrained.get("uv_rows"),
            "combined_conormal_passes": combined.get("passes"),
            "combined_conormal_verdict": combined.get("verdict"),
            "drive_sector_passes": drive.get("passes"),
            "drive_sector_field_count": len(drive_fields),
            "drive_projected_beta_l2_max": max(drive_beta_l2) if drive_beta_l2 else None,
            "completed_120_partner_benchmarks": partner.get("benchmarks"),
            "completed_120_partner_mass_locking_window": partner.get("mass_locking_window"),
        },
        "threshold_inputs_inherited_for_Audit3": {
            "rge_inputs_from_no_web_ledger": no_web.get("rge_inputs"),
            "audit0_dependency_order": audit0["deferred_audit_dependency_graph"]["order"],
            "required_before_Audit3": "a non-placeholder heavy_spectrum.json satisfying the schema above",
        },
        "source_sector_acceptance_gates": {
            "gauge_breaking": "vacuum little group must be the selected SM/Pati-Salam face",
            "no_extra_chiral_exotics": "no unpaired SM-charged zero modes beyond the protected three 16_i",
            "pseudo_goldstone_control": "all accidental charged flat directions must be lifted or explicitly declared auxiliary/composite",
            "landau_control": "visible/gut beta-function additions must not force a low Landau pole below the chosen UV trust ratio",
            "threshold_export": "every split heavy state must export a beta vector and mass convention for Audit 3",
            "proton_export": "colored triplet blocks must export mass inverse and Yukawa maps for Audit 2",
            "P_nu_source": "P_nu^c must be tied to a post-B-L source, e.g. a 126bar_H-like SM-singlet direction, or remain an explicit benchmark datum",
        },
        "known_blockers": [
            "No complete heavy_spectrum.json exists yet.",
            "M_star/v_R/source-scale convention remains unset for publication.",
            "Doublet-triplet splitting and colored-triplet inverse propagator are not exported in final Audit-2-ready form.",
            "The constrained 54/210 branch has local conormal evidence but not a complete microscopic UV completion.",
            "Route-D placement is local and cannot replace a global compactification spectrum.",
        ],
        "publication_boundary": {
            "not_claimed": [
                "complete Spin(10) source-sector UV action",
                "full heavy spectrum",
                "threshold closure",
                "d=5 proton safety",
                "global F-theory compactification",
            ],
            "claim": "Audit 4a schema fixes what a source-sector/heavy-spectrum artifact must provide.",
        },
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "source_spectrum_schema.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        "# Audit 4a Source-Sector / Heavy-Spectrum Schema",
        "",
        "Audit 4a has been opened as a schema and decision card, not as a completed UV source-sector model.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- audit0 sha256: `{audit0['card_sha256']}`",
        f"- source artifacts: {sum(1 for row in card['source_manifest'] if row['exists'])}/{len(card['source_manifest'])} present",
        "",
        "## Near-Term Track",
        "",
        "- Default: constrained 54/210 field-theory source-sector interface.",
        "- Fallback: CMSGUT-like 210+126+10 spectrum import if the constrained branch cannot export a complete heavy spectrum.",
        "- Optional parallel: Route-D global/string geometry, but not as a replacement for a heavy-spectrum artifact.",
        "",
        "## Required Heavy-Spectrum Fields",
        "",
        ", ".join(card["heavy_spectrum_schema"]["required_fields"]),
        "",
        "## Blockers",
        "",
    ]
    md_lines += [f"- {item}" for item in card["known_blockers"]]
    md_lines += [
        "",
        "## Boundary",
        "",
        "This file does not claim threshold closure, proton safety, a complete UV action, or a global compactification. It only defines the source/heavy-spectrum output contract required before Audit 3 and Audit 2 can be run as final audits.",
        "",
    ]
    md_path = OUT / "source_spectrum_schema.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {md_path.relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")


if __name__ == "__main__":
    main()
