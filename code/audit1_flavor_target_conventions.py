#!/usr/bin/env python3
"""Build Audit 1 target-table convention scaffold.

This is not a global flavor fit.  It fixes the observable list, legacy local
target anchors, basis/rotation conventions, delta-Z replay protocol, and output
contract required before a publication-grade CKM/PMNS/seesaw fit can be run.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit1"

SOURCES = {
    "audit1_builder": ROOT / "code" / "audit1_flavor_target_conventions.py",
    "audit0_card": ROOT / "output" / "audit0" / "invariant_card.json",
    "flavor_benchmark_card": ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json",
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


def sector_summary(flavor_card: dict[str, Any]) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for name, payload in flavor_card["yukawa_sectors"].items():
        rows[name] = {
            "normalization": payload.get("normalization"),
            "target": payload.get("target"),
            "normalized_singular_values": payload.get("normalized_singular_values"),
            "checks": payload.get("checks"),
        }
    return rows


def fit_observables_19() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name in ["y_u", "y_c", "y_t", "y_d", "y_s", "y_b"]:
        rows.append({"block": "quark_yukawas", "observable": name, "target_status": "REQUIRES_PUBLICATION_REFRESH"})
    for name in ["y_e", "y_mu", "y_tau"]:
        rows.append({"block": "charged_lepton_yukawas", "observable": name, "target_status": "REQUIRES_PUBLICATION_REFRESH"})
    for name in ["V_us", "V_cb", "V_ub", "J_CKM"]:
        rows.append({"block": "ckm", "observable": name, "target_status": "REQUIRES_PUBLICATION_REFRESH"})
    for name in ["Delta_m21_sq", "Delta_m31_sq", "sin2_theta12", "sin2_theta13", "sin2_theta23", "delta_CP"]:
        rows.append({"block": "neutrino_pmns", "observable": name, "target_status": "REQUIRES_PUBLICATION_REFRESH"})
    assert len(rows) == 19
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    audit0 = read_json(SOURCES["audit0_card"])
    flavor = read_json(SOURCES["flavor_benchmark_card"])
    no_web = read_json(SOURCES["no_web_input_convention_ledger"])

    local_flavor_targets = no_web["flavor_targets"]
    local_seesaw = no_web["seesaw_benchmark"]
    flavor_gates = no_web["flavor_gates"]
    seesaw_card = flavor["seesaw"]

    card: dict[str, Any] = {
        "audit": "audit1_flavor_target_conventions",
        "status": "target-table and fit-contract scaffold only; not a global flavor fit",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "audit0_card_sha256": audit0["card_sha256"],
        "observable_contract": {
            "fit_observable_count": 19,
            "fit_observables": fit_observables_19(),
            "pure_predictions_not_fit_targets": [
                "m_lightest",
                "Majorana_alpha21",
                "Majorana_alpha31",
                "m_beta_beta",
                "heavy_neutrino_masses",
            ],
            "publication_refresh_required": True,
            "publication_refresh_note": "This scaffold uses local legacy anchors only. A paper-grade fit must refresh quark/lepton/neutrino target values and uncertainties from cited PDG/nuFIT-like sources.",
        },
        "local_legacy_anchors": {
            "flavor_targets_from_no_web_ledger": local_flavor_targets,
            "flavor_gates_from_no_web_ledger": flavor_gates,
            "seesaw_benchmark_from_no_web_ledger": local_seesaw,
            "local_target_warning": "No web lookup is used here; values are local reproducibility anchors, not final publication targets.",
        },
        "basis_and_rotation_conventions": {
            "family_basis": audit0["scheme_labels"]["matrix_basis"],
            "flavor_basis_functions": flavor["basis"],
            "seesaw_basis_convention": seesaw_card["basis_convention"],
            "pmns_definition": seesaw_card["basis_convention"]["pmns_definition"],
            "charged_lepton_left_rotation": seesaw_card["basis_convention"]["charged_lepton_left_rotation"],
        },
        "benchmark_matrix_summary": {
            "flavor_card_all_checks_ok": flavor["all_checks_ok"],
            "seesaw_checks": seesaw_card["checks"],
            "light_mass_residual": seesaw_card["light_mass_residual"],
            "MR_condition_number": seesaw_card["MR_condition_number"],
            "heavy_neutrino_masses_GeV": seesaw_card["heavy_neutrino_masses_GeV"],
            "sector_summaries": sector_summary(flavor),
        },
        "delta_ZN_replay_protocol": audit0["delta_ZN_replay_protocol"],
        "required_fit_outputs_for_downstream_audits": {
            "for_Audit2_d5": [
                "U_uL",
                "U_dL",
                "U_eL",
                "U_nu",
                "physical Yukawa eigenvalues",
                "triplet-sector flavor tensors in physical basis",
            ],
            "for_Audit3_thresholds": [
                "matching scale choices",
                "tan_beta convention",
                "RGE scheme",
            ],
            "posterior_contract": [
                "best-fit parameter card",
                "pull table for all 19 observables",
                "covariance or scan cloud",
                "failure-gate report",
            ],
        },
        "parameter_policy": {
            "baseline_covariant_ansatz": "Use the CP1/O(2) covariant basis inherited from the Route-A note; do not add 120_H or a second spurion unless the failure gate is triggered.",
            "fallback_price": "Adding a 120_H antisymmetric channel or second spurion must be recorded as +6 real parameters per added channel.",
            "do_not_do": "Do not tune d=5 proton decay in Audit 1; Audit 1 only exports rotations/tensors consumed by Audit 2.",
        },
        "publication_boundary": {
            "not_claimed": [
                "current PDG/nuFIT target freshness",
                "global CKM/PMNS fit",
                "posterior uncertainty",
                "proton lifetime safety",
                "threshold closure",
            ],
            "claim": "Audit 1 scaffold fixes the target-table and fit-output contract.",
        },
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "target_table_conventions.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        "# Audit 1 Flavor / Seesaw Target-Table Convention",
        "",
        "Audit 1 has been opened as a convention scaffold, not a global fit.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- audit0 sha256: `{audit0['card_sha256']}`",
        f"- source artifacts: {sum(1 for row in card['source_manifest'] if row['exists'])}/{len(card['source_manifest'])} present",
        "",
        "## Fit Observable Contract",
        "",
        "- 19 fit observables: 6 quark Yukawas, 3 charged-lepton Yukawas, 4 CKM quantities, 6 PMNS/neutrino quantities.",
        "- Pure predictions: lightest mass, two Majorana phases, m_beta_beta, and heavy-neutrino spectrum.",
        "- Publication targets are intentionally marked `REQUIRES_PUBLICATION_REFRESH`.",
        "",
        "## Local Anchors",
        "",
        f"- legacy `V_us = {local_flavor_targets['Vus']}`",
        f"- legacy `V_cb = {local_flavor_targets['Vcb']}`",
        f"- legacy `V_ub = {local_flavor_targets['Vub']}`",
        f"- legacy `J_CKM = {local_flavor_targets['J']}`",
        f"- local seesaw `sin2_theta12 = {local_seesaw['sin2_theta12']}`",
        f"- local seesaw `sin2_theta23 = {local_seesaw['sin2_theta23']}`",
        f"- local seesaw `sin2_theta13 = {local_seesaw['sin2_theta13']}`",
        "",
        "## Boundary",
        "",
        "This file does not claim a completed flavor fit. It only fixes the target-table contract, legacy local anchors, basis conventions, and downstream outputs required by Audit 2 and Audit 3.",
        "",
    ]
    md_path = OUT / "target_table_conventions.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {md_path.relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")


if __name__ == "__main__":
    main()
