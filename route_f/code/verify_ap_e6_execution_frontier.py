#!/usr/bin/env python3
"""Integrate the three AP-E6 pre-portal execution lanes.

The source cards may be mechanically green while their physical closure
flags remain false.  This verifier therefore recomputes every lane closure
from an explicit conjunction and enforces the user rule that portal work may
start only after at least one complete pre-portal lane closes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()

SP4_CARD = OUTPUT / "ap_e6_sp4_eta_bordism_threshold.json"
B1_CARD = OUTPUT / "ap_e6_relaxed_b1_multigrid_hessian.json"
CALLIAS_CARD = OUTPUT / "ap_e6_same_soliton_yukawa_callias.json"

SOURCE_PATHS = [
    THIS_SCRIPT,
    ROUTE_F / "tex" / "ap_e6_sp4_eta_bordism_threshold.tex",
    ROUTE_F / "tex" / "ap_e6_sp4_eta_bordism_threshold.bib",
    ROUTE_F / "code" / "verify_ap_e6_sp4_eta_bordism_threshold.py",
    ROUTE_F / "tex" / "ap_e6_relaxed_b1_multigrid_hessian.tex",
    ROUTE_F / "tex" / "ap_e6_relaxed_b1_multigrid_hessian.bib",
    ROUTE_F / "code" / "scan_ap_e6_relaxed_b1_multigrid_hessian.py",
    ROUTE_F / "tex" / "ap_e6_same_soliton_yukawa_callias.tex",
    ROUTE_F / "tex" / "ap_e6_same_soliton_yukawa_callias.bib",
    ROUTE_F / "code" / "verify_ap_e6_same_soliton_yukawa_callias.py",
]

SP4_REQUIRED = (
    "sp4_euclidean_regulator_complete",
    "omega5_spin_bsp4_computed",
    "target_mapping_torus_eta_pair_selected",
    "heavy_threshold_eta_matched",
    "threshold_radiatively_protected",
    "strong_b1_phase_proven",
)

B1_REQUIRED = (
    "same_relaxed_b1_solution_computed",
    "multigrid_multivolume_converged",
    "sparse_projected_hessian_complete_in_declared_sector",
    "all_blocks_assembled_on_one_same_grid_and_boundary_condition",
    "discrete_stationarity_achieved",
    "lattice_topology_converged_to_B1",
    "physical_projected_hessian_at_stationary_discrete_solution",
    "projected_hessian_stability_gate_pass",
    "full_gauge_meson_ghost_fermion_hessian_complete",
    "dynamical_4d_importance_sampling_performed",
    "continuum_quantum_phase_proven",
)

CALLIAS_REQUIRED = (
    "same_physical_moduli_space_derived",
    "actual_yukawa_callias_operator_derived",
    "uniform_fredholm_gap_proven",
    "determinant_line_o2_derived",
    "physical_cpt_regulator_derived",
    "gauge_basic_wzw_descent_proven",
    "same_soliton_composition_closed",
)

CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    exists = path.is_file()
    return {
        "path": str(path.relative_to(REPO)),
        "exists": exists,
        "size_bytes": path.stat().st_size if exists else None,
        "sha256": sha256(path) if exists else None,
    }


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def mechanical_green(card: dict[str, Any]) -> bool:
    passed = card.get("checks_passed")
    total = card.get("checks_total")
    if passed is None or total is None:
        summary = card.get("summary", {})
        passed = summary.get("passed")
        total = summary.get("total")
    return bool(card.get("all_pass", passed == total) and passed == total and total)


def manifest_fresh(card: dict[str, Any]) -> tuple[bool, str]:
    rows = card.get("source_manifest", [])
    if not rows:
        return False, "source_manifest is absent or empty"
    stale: list[str] = []
    for row in rows:
        relative = row.get("path")
        recorded = row.get("sha256")
        if not relative or not recorded:
            stale.append(str(relative))
            continue
        path = REPO / relative
        if not path.is_file() or sha256(path) != recorded:
            stale.append(str(relative))
    return not stale, f"fresh={len(rows) - len(stale)}/{len(rows)}; stale={stale}"


def truth_table(card: dict[str, Any], names: tuple[str, ...]) -> dict[str, bool]:
    return {name: card.get(name) is True for name in names}


def markdown(result: dict[str, Any]) -> str:
    summary = result["summary"]
    portal = result["portal"]
    lines = [
        "# AP-E6 execution-frontier gate",
        "",
        f"- Mechanical checks: **{summary['passed']}/{summary['total']}**",
        f"- Sp(4) eta/threshold lane closed: **{result['route_closure']['sp4_eta_threshold']}**",
        f"- Relaxed B=1/full-Hessian lane closed: **{result['route_closure']['relaxed_b1_full_hessian']}**",
        f"- Same-soliton Callias/CPT/descent lane closed: **{result['route_closure']['same_soliton_callias_descent']}**",
        f"- Canonical same-background SHA-256: **{result['canonical_same_background_sha256']}**",
        f"- Any pre-portal route closed: **{result['any_preportal_route_closed']}**",
        f"- Portal start authorized: **{portal['start_authorized']}**",
        f"- Portal constructed: **{portal['constructed']}**",
        f"- Physics promotion allowed: **{result['physics_promotion_allowed']}**",
        "",
        "## Decision",
        "",
        portal["decision"],
        "",
        "## Recomputed closure inputs",
        "",
    ]
    for lane, flags in result["closure_inputs"].items():
        lines.append(f"### {lane}")
        lines.append("")
        lines.extend(f"- `{key}`: `{str(value).lower()}`" for key, value in flags.items())
        lines.append("")
    lines.extend(["## Checks", ""])
    lines.extend(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    source_manifest = [source_row(path) for path in SOURCE_PATHS]
    check(
        "provenance",
        "all AP-E6 source sets exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    sp4 = load(SP4_CARD)
    b1 = load(B1_CARD)
    callias = load(CALLIAS_CARD)
    cards = {"sp4": sp4, "b1": b1, "callias": callias}

    for name, card in cards.items():
        check(
            "source_cards",
            f"{name} mechanical card is green",
            mechanical_green(card),
            f"checks={card.get('checks_passed', card.get('summary', {}).get('passed'))}/"
            f"{card.get('checks_total', card.get('summary', {}).get('total'))}",
        )
        check(
            "source_cards",
            f"{name} card forbids premature physics promotion",
            card.get("physics_promotion_allowed") is False,
            f"physics_promotion_allowed={card.get('physics_promotion_allowed')}",
        )
        fresh, detail = manifest_fresh(card)
        check(
            "provenance",
            f"{name} source-card manifest is fresh",
            fresh,
            detail,
        )

    sp4_inputs = truth_table(sp4, SP4_REQUIRED)
    b1_inputs = truth_table(b1, B1_REQUIRED)
    callias_inputs = truth_table(callias, CALLIAS_REQUIRED)

    b1_checksum = (
        b1.get("coupled_bvp", {})
        .get("canonical_same_solution", {})
        .get("profile_sha256_r_F_s_little_endian_float64")
    )
    callias_checksum = callias.get("same_profile", {}).get(
        "profile_sha256_r_F_s_little_endian_float64"
    )
    check(
        "same_background",
        "Callias/CPT/WZW card consumes the canonical relaxed B=1 profile",
        isinstance(b1_checksum, str)
        and len(b1_checksum) == 64
        and callias_checksum == b1_checksum,
        f"B1={b1_checksum}; Callias={callias_checksum}",
    )

    sp4_closed = all(sp4_inputs.values())
    b1_closed = all(b1_inputs.values())
    callias_closed = all(callias_inputs.values())

    for name, card, closed, inputs in (
        ("sp4", sp4, sp4_closed, sp4_inputs),
        ("b1", b1, b1_closed, b1_inputs),
        ("callias", callias, callias_closed, callias_inputs),
    ):
        check(
            "closure_logic",
            f"{name} source lane flag equals the independent conjunction",
            card.get("lane_closed") is closed,
            f"source={card.get('lane_closed')}; recomputed={closed}; inputs={inputs}",
        )

    any_closed = sp4_closed or b1_closed or callias_closed
    portal_started = False
    check(
        "portal",
        "portal construction obeys the one-complete-lane prerequisite",
        not portal_started or any_closed,
        f"any_closed={any_closed}; constructed={portal_started}",
    )

    if any_closed:
        decision = (
            "At least one complete pre-portal lane is closed.  Portal work is now "
            "authorized, but no degree-one portal is constructed by this execution card."
        )
    else:
        decision = (
            "No complete pre-portal lane is closed.  In accordance with the ordered "
            "research rule, the degree-one Route-E portal is not started."
        )

    passed = sum(row["pass"] for row in CHECKS)
    total = len(CHECKS)
    result = {
        "schema_version": "ap_e6_execution_frontier_v1",
        "all_pass": passed == total,
        "checks_passed": passed,
        "checks_total": total,
        "summary": {"passed": passed, "total": total},
        "source_card_counts": {
            name: {
                "passed": card.get("checks_passed", card.get("summary", {}).get("passed")),
                "total": card.get("checks_total", card.get("summary", {}).get("total")),
            }
            for name, card in cards.items()
        },
        "closure_inputs": {
            "sp4_eta_threshold": sp4_inputs,
            "relaxed_b1_full_hessian": b1_inputs,
            "same_soliton_callias_descent": callias_inputs,
        },
        "route_closure": {
            "sp4_eta_threshold": sp4_closed,
            "relaxed_b1_full_hessian": b1_closed,
            "same_soliton_callias_descent": callias_closed,
        },
        "canonical_same_background_sha256": b1_checksum,
        "any_preportal_route_closed": any_closed,
        "portal": {
            "start_authorized": any_closed,
            "constructed": portal_started,
            "decision": decision,
        },
        "physics_promotion_allowed": False,
        "checks": CHECKS,
        "source_manifest": source_manifest,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e6_execution_frontier.json"
    md_path = OUTPUT / "ap_e6_execution_frontier.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(markdown(result), encoding="utf-8")

    print(f"AP-E6 execution frontier: {passed}/{total} checks pass")
    print(f"any_preportal_route_closed={str(any_closed).lower()}")
    print(f"portal_constructed={str(portal_started).lower()}")
    print("physics_promotion_allowed=false")
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
