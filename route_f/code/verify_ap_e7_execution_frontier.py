#!/usr/bin/env python3
"""Integrate the three AP-E7 pre-portal branches without promoting partial results.

The regulator, discrete-soliton, and determinant-line cards each contain
mathematically positive subresults.  A physical lane closes only when every
required same-model condition is true.  The FR and rank-three constructions
are alternatives, while the conditions inside either alternative are
conjunctive.  Portal work remains forbidden until at least one complete lane
closes.
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

APS_CARD = OUTPUT / "ap_e7_common_sp4_aps_pv_regulator.json"
LATTICE_CARD = OUTPUT / "ap_e7_discrete_rerelax_superhessian.json"
FAMILY_CARD = OUTPUT / "ap_e7_so3_fr_rank3_family.json"
E6_CALLIAS_CARD = OUTPUT / "ap_e6_same_soliton_yukawa_callias.json"

SOURCE_PATHS = [
    THIS_SCRIPT,
    ROUTE_F / "tex" / "ap_e7_common_sp4_aps_pv_regulator.tex",
    ROUTE_F / "tex" / "ap_e7_common_sp4_aps_pv_regulator.bib",
    ROUTE_F / "code" / "verify_ap_e7_common_sp4_aps_pv_regulator.py",
    ROUTE_F / "tex" / "ap_e7_discrete_rerelax_superhessian.tex",
    ROUTE_F / "tex" / "ap_e7_discrete_rerelax_superhessian.bib",
    ROUTE_F / "code" / "scan_ap_e7_discrete_rerelax_superhessian.py",
    ROUTE_F / "tex" / "ap_e7_so3_fr_rank3_family.tex",
    ROUTE_F / "tex" / "ap_e7_so3_fr_rank3_family.bib",
    ROUTE_F / "code" / "verify_ap_e7_so3_fr_rank3_family.py",
]

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


def counts(card: dict[str, Any]) -> tuple[int | None, int | None]:
    passed = card.get("checks_passed")
    total = card.get("checks_total")
    if passed is None or total is None:
        summary = card.get("summary", {})
        passed = summary.get("checks_passed", summary.get("passed"))
        total = summary.get("checks_total", summary.get("total"))
    return passed, total


def mechanically_green(card: dict[str, Any]) -> bool:
    passed, total = counts(card)
    all_pass = card.get("all_pass", passed == total)
    return bool(all_pass and isinstance(total, int) and total > 0 and passed == total)


def manifest_rows(card: dict[str, Any]) -> list[dict[str, Any]]:
    rows = card.get("source_manifest")
    if rows is None:
        rows = card.get("sources", [])
    return rows


def manifest_fresh(card: dict[str, Any]) -> tuple[bool, str]:
    rows = manifest_rows(card)
    if not rows:
        return False, "source manifest is absent or empty"
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


def explicit_false(card: dict[str, Any], name: str) -> bool:
    if name in card:
        return card[name] is False
    gates = card.get("gates", {})
    if name in gates:
        return gates[name] is False
    fail_closed = card.get("fail_closed", {})
    if name in fail_closed:
        return fail_closed[name] is False
    return False


def markdown(result: dict[str, Any]) -> str:
    lines = [
        "# AP-E7 execution-frontier gate",
        "",
        f"- Mechanical checks: **{result['checks_passed']}/{result['checks_total']}**",
        f"- Common APS/PV physical lane closed: **{result['route_closure']['common_aps_pv']}**",
        f"- Unanchored B=1/full-super-Hessian lane closed: **{result['route_closure']['unanchored_b1_superhessian']}**",
        f"- SO(3)/FR alternative closed: **{result['route_closure']['so3_fr_alternative']}**",
        f"- Rank-three alternative closed: **{result['route_closure']['rank_three_alternative']}**",
        f"- Complete determinant/family/descent lane closed: **{result['route_closure']['determinant_family_descent']}**",
        f"- Any pre-portal route closed: **{result['any_preportal_route_closed']}**",
        f"- Degree-one portal started: **{result['degree_one_portal_started']}**",
        f"- Physics promotion allowed: **{result['physics_promotion_allowed']}**",
        "",
        "## Decision",
        "",
        result["decision"],
        "",
        "## Recomputed closure inputs",
        "",
    ]
    for lane, values in result["closure_inputs"].items():
        lines.append(f"### {lane}")
        lines.append("")
        lines.extend(f"- `{name}`: `{str(value).lower()}`" for name, value in values.items())
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
        "all AP-E7 source sets exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    aps = load(APS_CARD)
    lattice = load(LATTICE_CARD)
    family = load(FAMILY_CARD)
    e6_callias = load(E6_CALLIAS_CARD)
    cards = {"aps": aps, "lattice": lattice, "family": family}

    for name, card in cards.items():
        passed, total = counts(card)
        check(
            "source_cards",
            f"{name} mechanical card is green",
            mechanically_green(card),
            f"checks={passed}/{total}; status={card.get('status')}",
        )
        fresh, detail = manifest_fresh(card)
        check("provenance", f"{name} source-card manifest is fresh", fresh, detail)

    check(
        "promotion_boundary",
        "APS card explicitly forbids promotion",
        explicit_false(aps, "physics_promotion_authorized"),
        f"value={aps.get('gates', {}).get('physics_promotion_authorized')}",
    )
    check(
        "promotion_boundary",
        "lattice card explicitly forbids promotion",
        explicit_false(lattice, "physics_promotion_allowed"),
        f"value={lattice.get('physics_promotion_allowed')}",
    )
    check(
        "promotion_boundary",
        "family card explicitly forbids promotion",
        explicit_false(family, "physics_promotion_allowed"),
        f"value={family.get('physics_promotion_allowed')}",
    )

    frozen_checksum = "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59"
    lattice_checksum = lattice.get("canonical_profile_checksum")
    family_checksum = family.get("same_profile", {}).get(
        "profile_sha256_r_F_s_little_endian_float64"
    )
    check(
        "same_background",
        "lattice and family branches consume the frozen AP-E6 profile",
        lattice_checksum == family_checksum == frozen_checksum,
        f"lattice={lattice_checksum}; family={family_checksum}",
    )

    aps_gates = aps.get("gates", {})
    aps_inputs = {
        "common_quadratic_aps_pv_regulator_specified": aps_gates.get(
            "common_quadratic_aps_pv_regulator_specified"
        )
        is True,
        "five_dimensional_aps_domain_specified": aps_gates.get(
            "five_dimensional_aps_domain_specified"
        )
        is True,
        "pure_unbroken_gauge_gravity_heavy_phase_computed": aps_gates.get(
            "pure_unbroken_gauge_gravity_heavy_phase_computed"
        )
        is True,
        "full_nonperturbative_sp4_measure_constructed": aps_gates.get(
            "full_nonperturbative_sp4_measure_constructed"
        )
        is True,
        "original_action_target_mass_family_lift_defined": aps_gates.get(
            "original_action_target_mass_family_lift_defined"
        )
        is True,
        "target_mapping_torus_eta_pair_selected": aps_gates.get(
            "target_mapping_torus_eta_pair_selected"
        )
        is True,
        "mixed_wzw_heavy_determinant_computed": aps_gates.get(
            "mixed_wzw_heavy_determinant_computed"
        )
        is True,
        "k_plus_two_all_scale_matched": aps_gates.get("k_plus_two_all_scale_matched")
        is True,
    }
    aps_closed = all(aps_inputs.values())
    check(
        "closure_logic",
        "APS/PV source lane equals the full physical conjunction",
        aps_gates.get("sp4_aps_pv_lane_closed") is aps_closed,
        f"source={aps_gates.get('sp4_aps_pv_lane_closed')}; recomputed={aps_closed}",
    )

    lattice_fail = lattice.get("fail_closed", {})
    lattice_inputs = {
        "genuine_unanchored_B1_lattice_stationary_point_found": lattice_fail.get(
            "genuine_unanchored_B1_lattice_stationary_point_found"
        )
        is True,
        "admissibility_independent_continuum_B1_solution_constructed": lattice_fail.get(
            "admissibility_independent_continuum_B1_solution_constructed"
        )
        is True,
        "physical_aggregate_superhessian_complete": lattice.get(
            "physical_aggregate_superhessian_complete"
        )
        is True,
        "interacting_gauge_meson_cross_blocks_computed": lattice_fail.get(
            "interacting_gauge_meson_cross_blocks_computed"
        )
        is True,
        "fermion_determinant_bosonic_second_variation_computed": lattice_fail.get(
            "fermion_determinant_bosonic_second_variation_computed"
        )
        is True,
        "projected_stability_continuum_extrapolated": lattice_fail.get(
            "projected_stability_continuum_extrapolated"
        )
        is True,
        "four_dimensional_dynamics_performed": lattice.get(
            "four_dimensional_dynamics_performed"
        )
        is True,
        "importance_sampling_or_hmc_performed": lattice_fail.get(
            "importance_sampling_or_hmc_performed"
        )
        is True,
    }
    lattice_closed = all(lattice_inputs.values())
    check(
        "closure_logic",
        "lattice source lane equals the unanchored physical conjunction",
        lattice.get("lane_closed") is lattice_closed,
        f"source={lattice.get('lane_closed')}; recomputed={lattice_closed}",
    )

    fr_inputs = {
        "so3_torsion_lines_classified": family.get("so3_torsion_lines_classified")
        is True,
        "actual_yukawa_mapping_torus_mod2_index_computed": family.get(
            "actual_yukawa_mapping_torus_mod2_index_computed"
        )
        is True,
        "actual_pfaffian_torsion_class_selected": family.get(
            "actual_pfaffian_torsion_class_selected"
        )
        is True,
        "microscopic_cpt_regulator_constructed": family.get(
            "microscopic_cpt_regulator_constructed"
        )
        is True,
        "same_soliton_family_uniform_fredholm_gap_proven": family.get(
            "same_soliton_family_uniform_fredholm_gap_proven"
        )
        is True,
    }
    rank_three_inputs = {
        "rank_three_topological_mass_family_constructed": family.get(
            "rank_three_topological_mass_family_constructed"
        )
        is True,
        "rank_three_c1_target_realized": family.get("rank_three_c1_target_realized")
        is True,
        "rank_three_uniform_matrix_gap_proven": family.get(
            "rank_three_uniform_matrix_gap_proven"
        )
        is True,
        "physical_rank_three_yukawa_embedding_derived": family.get(
            "physical_rank_three_yukawa_embedding_derived"
        )
        is True,
        "same_soliton_family_uniform_fredholm_gap_proven": family.get(
            "same_soliton_family_uniform_fredholm_gap_proven"
        )
        is True,
        "microscopic_cpt_regulator_constructed": family.get(
            "microscopic_cpt_regulator_constructed"
        )
        is True,
    }
    fr_closed = all(fr_inputs.values())
    rank_three_closed = all(rank_three_inputs.values())
    descent_proven = e6_callias.get("gauge_basic_wzw_descent_proven") is True
    family_closed = (fr_closed or rank_three_closed) and descent_proven
    check(
        "closure_logic",
        "family source lane agrees with alternatives plus prior descent gate",
        family.get("lane_closed") is family_closed,
        "source={}; FR={}; rank3={}; descent={}; recomputed={}".format(
            family.get("lane_closed"),
            fr_closed,
            rank_three_closed,
            descent_proven,
            family_closed,
        ),
    )

    any_closed = aps_closed or lattice_closed or family_closed
    portal_started = False
    physics_promotion_allowed = False
    check(
        "portal",
        "no portal is started while every complete pre-portal lane is open",
        (not portal_started) and (not any_closed),
        f"any_closed={any_closed}; portal_started={portal_started}",
    )
    check(
        "portal",
        "physics promotion remains false",
        not physics_promotion_allowed,
        f"physics_promotion_allowed={physics_promotion_allowed}",
    )

    passed = sum(row["pass"] for row in CHECKS)
    total = len(CHECKS)
    decision = (
        "The common quadratic regulator and pure heavy threshold are real advances; "
        "the finite-site topology theorem and rank-three bundle are also exact.  "
        "Nevertheless the original target mass-family lift, an unanchored discrete "
        "B=1 stationary point with its physical super-Hessian, and an actual same-"
        "soliton determinant/family/descent construction all remain open.  The "
        "degree-one Route-E portal is therefore not started."
    )
    result = {
        "schema_version": "ap_e7_execution_frontier_v1",
        "status": "mechanical_pass_physics_fail_closed"
        if passed == total
        else "mechanical_failure_physics_fail_closed",
        "all_pass": passed == total,
        "checks_passed": passed,
        "checks_total": total,
        "source_manifest": source_manifest,
        "source_card_counts": {
            name: {"passed": counts(card)[0], "total": counts(card)[1]}
            for name, card in cards.items()
        },
        "canonical_same_background_sha256": frozen_checksum,
        "closure_inputs": {
            "common_aps_pv": aps_inputs,
            "unanchored_b1_superhessian": lattice_inputs,
            "so3_fr_alternative": fr_inputs,
            "rank_three_alternative": rank_three_inputs,
            "prior_same_soliton_descent": {
                "gauge_basic_wzw_descent_proven": descent_proven
            },
        },
        "route_closure": {
            "common_aps_pv": aps_closed,
            "unanchored_b1_superhessian": lattice_closed,
            "so3_fr_alternative": fr_closed,
            "rank_three_alternative": rank_three_closed,
            "determinant_family_descent": family_closed,
        },
        "any_preportal_route_closed": any_closed,
        "degree_one_portal_started": portal_started,
        "degree_one_portal_constructed": False,
        "physics_promotion_allowed": physics_promotion_allowed,
        "decision": decision,
        "checks": CHECKS,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e7_execution_frontier.json"
    md_path = OUTPUT / "ap_e7_execution_frontier.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(markdown(result), encoding="utf-8")
    print(
        f"AP-E7 frontier: {passed}/{total}; status={result['status']}; "
        f"any_closed={any_closed}; portal_started={portal_started}"
    )
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
