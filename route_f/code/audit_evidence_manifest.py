#!/usr/bin/env python3
"""Audit a P0 subset of Route-E dynamics evidence cited by the documents.

The canonical implementations live under ``route_E/code_dyn``.  The historical
root ``code`` entry points are thin compatibility delegates, so both the cited
path and the single implementation are recorded without treating duplicate
physics logic as independent evidence.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
OUTPUT = REPO / "route_f" / "output"

CORE = [
    "route_E/code/verify_route_e_first_principles.py",
    "route_E/code/verify_route_e_dependency_closure.py",
    "route_E/output/route_e_first_principles.json",
    "route_E/output/route_e_dependency_closure.json",
]

SOURCE_SCRIPTS = [
    "code/audit0_dyn_conventions_and_inventory.py",
    "code/audit1_dyn4a_seesaw_replay_zeta_posterior.py",
    "code/audit1_dyn4b_refreshed_card_unconditional_zeta.py",
    "code/audit1_dyn4c_kernel_dirac_refit.py",
    "code/audit2_dyn3_proton_d5_kill_criterion.py",
    "code/audit3_dyn2_thresholds_unification.py",
    "code/audit3_dyn2b_rescue_scan.py",
    "code/audit4a_dyn1a_vacuum_goldstone_audit.py",
    "code/audit4a_dyn1b_full_spectrum.py",
    "code/audit5_dyn5_messenger_one_loop.py",
    "code/audit7_dyn7_leptogenesis_argzeta.py",
    "code/audit8_dyn8_falsifiability_collection.py",
    "code/audit9_dyn9_nonsusy_intermediate.py",
    "code/audit9_dyn9b1_nonsusy_vacuum_thresholds.py",
    "code/audit9_dyn9b1b_210_quartic_descent.py",
    "code/audit9_dyn9b1c_eps_and_126_quartics.py",
    "code/audit9_dyn9b1d_lr_ratio_scan.py",
    "code/audit9_dyn9b2_nonsusy_flavor_refit.py",
    "code/audit9_dyn9b3_nonsusy_leptogenesis.py",
]

DYN_LEDGERS = [
    "output/audit8/dyn8_falsifiability_collection.json",
    "output/audit9/dyn9_nonsusy_intermediate.json",
    "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
    "output/audit9/dyn9b3_nonsusy_leptogenesis.json",
]

STRING_CARDS = [
    "route_d/code/verify_d3_instanton_majorana_pricing.py",
    "route_d/code/verify_d4_stueckelberg_protection.py",
    "route_d/code/verify_d5_susy_breaking_bridge.py",
    "route_d/output/d3_instanton_majorana_pricing.json",
    "route_d/output/d4_stueckelberg_protection.json",
    "route_d/output/d5_susy_breaking_bridge.json",
]

# New Route-F execution/validity evidence is reported separately from the
# historical 29-item cited-path subset so progress is visible without falsely
# changing the paper-citation closure denominator.
EXECUTION_GUARDS = [
    "route_E/code_dyn/route_e_paths.py",
    "route_E/code_dyn/run_route_e_dynamics.py",
    "route_E/code_dyn/dyn_claim_registry.json",
    "route_E/code_dyn/requirements.txt",
    "route_E/code_dyn/audit5_dyn5_model_validity.py",
    "route_E/code_dyn/audit7_dyn7_flavor_regime_gate.py",
    "output/audit1/dyn4a_seesaw_zeta_posterior.json",
    "output/audit5/dyn5_model_validity.json",
    "output/audit7/dyn7_flavor_regime_gate.json",
    "route_f/code/audit_blocker_promotion_gate.py",
    "route_f/output/blocker_promotion_gate.json",
    "route_f/output/blocker_promotion_gate.md",
]

DYNAMICS = SOURCE_SCRIPTS + DYN_LEDGERS + STRING_CARDS
REQUIRED_SOURCE_PATHS = (
    [path for path in CORE if path.endswith(".py")]
    + SOURCE_SCRIPTS
    + [path for path in STRING_CARDS if path.endswith(".py")]
)


def locate(relative: str) -> list[str]:
    candidates = [REPO / relative]
    if relative.startswith("code/") or relative.startswith("output/"):
        candidates.append(REPO / "route_E" / relative)
    if relative.startswith("code/"):
        candidates.append(REPO / "route_E" / "code_dyn" / Path(relative).name)
    return [str(path.relative_to(REPO)) for path in candidates if path.exists()]


def tracked(relative: str) -> bool:
    return subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", relative],
        cwd=REPO,
        capture_output=True,
        text=True,
    ).returncode == 0


def records(paths: list[str]) -> list[dict[str, object]]:
    return [
        {
            "cited_path": path,
            "found": bool(locate(path)),
            "exact_path_found": (REPO / path).exists(),
            "locations": locate(path),
            "tracked_locations": [location for location in locate(path) if tracked(location)],
        }
        for path in paths
    ]


def main() -> None:
    core_records = records(CORE)
    source_records = records(SOURCE_SCRIPTS)
    ledger_records = records(DYN_LEDGERS)
    string_records = records(STRING_CARDS)
    guard_records = records(EXECUTION_GUARDS)
    dynamics_records = records(DYNAMICS)
    core_present = sum(bool(item["found"]) for item in core_records)
    source_present = sum(bool(item["found"]) for item in source_records)
    source_exact = sum(bool(item["exact_path_found"]) for item in source_records)
    source_relocated_tracked = sum(
        bool(item["tracked_locations"]) for item in source_records
    )
    ledger_present = sum(bool(item["found"]) for item in ledger_records)
    string_present = sum(bool(item["found"]) for item in string_records)
    dynamics_present = sum(bool(item["found"]) for item in dynamics_records)
    guards_present = sum(bool(item["found"]) for item in guard_records)
    guards_tracked = sum(bool(item["tracked_locations"]) for item in guard_records)
    required_sources_exact = sum((REPO / path).exists() for path in REQUIRED_SOURCE_PATHS)
    required_sources_tracked = sum(tracked(path) for path in REQUIRED_SOURCE_PATHS)
    if core_present != len(core_records):
        status = "fail_missing_core_artifacts"
    elif dynamics_present != len(dynamics_records):
        status = "fail_incomplete_cited_dynamics_evidence"
    elif required_sources_exact != len(REQUIRED_SOURCE_PATHS):
        status = "fail_relocated_or_missing_cited_sources"
    elif required_sources_tracked != len(REQUIRED_SOURCE_PATHS):
        status = "fail_untracked_cited_sources"
    else:
        status = "pass_workspace_existence_subset"
    result = {
        "status": status,
        "scope": (
            "29-item dynamics P0 file-existence subset plus 4 core artifacts; "
            "not an exhaustive extraction of all cited paths and not a validation "
            "of scientific correctness"
        ),
        "interpretation": (
            "Root cited paths are compatibility delegates to one canonical "
            "route_E/code_dyn implementation. File existence and tracking satisfy "
            "only the workspace evidence gate; scientific promotion is governed "
            "separately by the claim registry and blocker promotion gate."
        ),
        "required_source_provenance": {
            "exact_paths_present": required_sources_exact,
            "git_tracked": required_sources_tracked,
            "expected": len(REQUIRED_SOURCE_PATHS),
            "paths": REQUIRED_SOURCE_PATHS,
        },
        "core": {
            "present": core_present,
            "expected": len(core_records),
            "records": core_records,
        },
        "source_scripts": {
            "present_or_relocated": source_present,
            "exact_cited_path_present": source_exact,
            "present_locations_git_tracked": source_relocated_tracked,
            "expected": len(source_records),
            "records": source_records,
        },
        "selected_dyn_ledgers": {
            "present": ledger_present,
            "expected": len(ledger_records),
            "records": ledger_records,
        },
        "string_cards": {
            "present": string_present,
            "expected": len(string_records),
            "records": string_records,
        },
        "execution_guard_evidence": {
            "present": guards_present,
            "git_tracked": guards_tracked,
            "expected": len(guard_records),
            "note": "supplemental 2026-07-14 evidence; excluded from the historical 29-item cited-path denominator",
            "records": guard_records,
        },
        "cited_dynamics": {
            "present_or_relocated": dynamics_present,
            "expected": len(dynamics_records),
            "missing": [item["cited_path"] for item in dynamics_records if not item["found"]],
            "records": dynamics_records,
        },
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "route_e_evidence_audit.json"
    md_path = OUTPUT / "route_e_evidence_audit.md"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    missing_lines = "\n".join(
        f"- `{path}`" for path in result["cited_dynamics"]["missing"]
    )
    md_content = (
        "# Route-E P0 Evidence-Subset Existence Audit\n\n"
        f"- status: **{result['status']}**\n"
        f"- core artifacts present: `{core_present}/{len(core_records)}`\n"
        f"- DYN source scripts present or relocated: "
        f"`{source_present}/{len(source_records)}`\n"
        f"- DYN source scripts at the exact cited paths: "
        f"`{source_exact}/{len(source_records)}`\n"
        f"- DYN source scripts with at least one Git-tracked found location: "
        f"`{source_relocated_tracked}/{len(source_records)}`\n"
        f"- all required source paths present/tracked: "
        f"`{required_sources_exact}/{len(REQUIRED_SOURCE_PATHS)}` present, "
        f"`{required_sources_tracked}/{len(REQUIRED_SOURCE_PATHS)}` tracked\n"
        f"- selected generated DYN ledgers present: "
        f"`{ledger_present}/{len(ledger_records)}`\n"
        f"- Route-E string-card scripts/ledgers present: "
        f"`{string_present}/{len(string_records)}`\n"
        f"- 2026-07-14 execution/guard evidence present/tracked: "
        f"`{guards_present}/{len(guard_records)}` present, "
        f"`{guards_tracked}/{len(guard_records)}` tracked "
        "(supplemental; not in the 29-item denominator)\n"
        f"- cited dynamics artifacts present: "
        f"`{dynamics_present}/{len(dynamics_records)}`\n"
        "- scope: 29-item dynamics subset plus 4 core artifacts; it is not "
        "every cited path, and a present file is not automatically correct\n"
        "- root DYN entry points delegate to one canonical "
        "`route_E/code_dyn` implementation; this prevents ledger overwrite drift\n\n"
        "## Missing cited dynamics artifacts\n\n"
        f"{missing_lines}\n"
    )
    md_path.write_text(md_content.rstrip() + "\n", encoding="utf-8")
    print(
        "route_e_evidence_audit: "
        f"{result['status']} ({dynamics_present}/{len(dynamics_records)} dynamics artifacts present)"
    )
    if not result["status"].startswith("pass_"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
