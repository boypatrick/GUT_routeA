#!/usr/bin/env python3
"""Integrate the three 2026-07-16 AP-E completion-frontier audits.

This verifier is intentionally a *promotion gate*, not another microscopic
calculation.  It loads the finite-lattice/Hessian, APS/global-form, and
same-soliton SQM/Callias/descent cards, checks their positive mathematical or
numerical results, and then computes whether any lane is genuinely closed.
The degree-one Route-E portal is authorized only after at least one complete
pre-portal lane closes.  A green result may therefore (and presently does)
certify that portal construction remains forbidden.
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

LATTICE_CARD = OUTPUT / "ap_e5_4d_qc2d_lattice_hessian.json"
APS_CARD = OUTPUT / "ap_e3_aps_global_form_search.json"
SQM_CARD = OUTPUT / "ap_e4_same_soliton_callias_descent.json"

SOURCE_PATHS = [
    THIS_SCRIPT,
    ROUTE_F / "tex" / "ap_e5_4d_qc2d_lattice_hessian.tex",
    ROUTE_F / "tex" / "ap_e5_4d_qc2d_lattice_hessian.bib",
    ROUTE_F / "code" / "verify_ap_e5_4d_qc2d_lattice_hessian.py",
    ROUTE_F / "tex" / "ap_e3_aps_global_form_search.tex",
    ROUTE_F / "tex" / "ap_e3_aps_global_form_search.bib",
    ROUTE_F / "code" / "verify_ap_e3_aps_global_form_search.py",
    ROUTE_F / "tex" / "ap_e4_same_soliton_callias_descent.tex",
    ROUTE_F / "tex" / "ap_e4_same_soliton_callias_descent.bib",
    ROUTE_F / "code" / "verify_ap_e4_same_soliton_callias_descent.py",
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
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def markdown(result: dict[str, Any]) -> str:
    routes = result["route_closure"]
    portal = result["degree_one_portal"]
    lines = [
        "# AP-E completion-frontier promotion gate",
        "",
        f"- Mechanical checks: **{result['summary']['passed']}/{result['summary']['total']}**",
        f"- 4D lattice/full-Hessian lane closed: **{routes['four_dimensional_lattice_and_full_hessian']}**",
        f"- APS/semisimple all-scale lane closed: **{routes['aps_and_semisimple_uv']}**",
        f"- Same-soliton SQM/Callias/descent lane closed: **{routes['same_soliton_sqm_callias_descent']}**",
        f"- Any pre-portal lane closed: **{result['any_preportal_route_closed']}**",
        f"- Degree-one portal authorized: **{portal['authorized']}**",
        f"- Degree-one portal constructed: **{portal['constructed']}**",
        f"- Physics promotion allowed: **{result['physics_promotion_allowed']}**",
        "",
        "## Decision",
        "",
        portal["decision"],
        "",
        "## Portal theorem prerequisites",
        "",
    ]
    lines.extend(
        f"- `{name}`: `{str(value).lower()}`"
        for name, value in portal["prerequisites"].items()
    )
    lines.extend(["", "## Checks", ""])
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
        "all three audit source sets exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    lattice = load(LATTICE_CARD)
    aps = load(APS_CARD)
    sqm = load(SQM_CARD)

    lattice_summary = lattice["summary"]
    aps_passed = int(aps["checks_passed"])
    aps_total = int(aps["checks_total"])
    sqm_summary = sqm["summary"]
    check(
        "source_cards",
        "all three mechanical cards are green",
        lattice["status"] == "PASS"
        and lattice_summary["passed"] == lattice_summary["total"]
        and aps["all_pass"] is True
        and aps_passed == aps_total
        and sqm_summary["passed"] == sqm_summary["total"],
        (
            f"lattice={lattice_summary['passed']}/{lattice_summary['total']}; "
            f"APS={aps_passed}/{aps_total}; "
            f"SQM={sqm_summary['passed']}/{sqm_summary['total']}"
        ),
    )

    lattice_fail = lattice["fail_closed"]
    lattice_topology = lattice["topology"]
    lattice_gauge = lattice["gauge_and_ghost"]
    lattice_bosons = lattice["boson_hessian"]
    lattice_fermion = lattice["fermion"]
    lattice_scaling = lattice["continuum_scaling"]
    check(
        "lattice",
        "sampled B=1 background has a controlled continuum extrapolation",
        lattice_topology["extrapolation_abs_error"] < 2.0e-3,
        f"B0={lattice_topology['quadratic_in_a_squared_extrapolation_last_four']:.12g}",
    )
    check(
        "lattice",
        "gauge and ghost determinant xi dependence cancels up to the predicted constant",
        lattice_gauge["maximum_xi_corrected_residual"] < 1.0e-10,
        f"residual={lattice_gauge['maximum_xi_corrected_residual']:.3e}",
    )
    check(
        "lattice",
        "declared meson and diquark blocks are positive with a tachyonic control",
        lattice_bosons["meson_minimum_eigenvalue"] > 0.0
        and lattice_bosons["diquark_minimum_eigenvalue"] > 0.0
        and lattice_bosons["negative_control_minimum_eigenvalue"] < 0.0,
        (
            f"meson={lattice_bosons['meson_minimum_eigenvalue']:.9g}; "
            f"diquark={lattice_bosons['diquark_minimum_eigenvalue']:.9g}; "
            f"control={lattice_bosons['negative_control_minimum_eigenvalue']:.9g}"
        ),
    )
    check(
        "lattice",
        "Wilson-Dirac Schur and one-loop curvature controls agree",
        lattice_fermion["minimum_singular_value"] > 1.0e-3
        and lattice_fermion["determinant_logabs_residual"] < 1.0e-10
        and lattice_fermion["schur_solve_relative_residual"] < 1.0e-10
        and lattice_fermion["one_loop_collective_coordinate"]["relative_error"]
        < 2.0e-6,
        (
            f"sigma_min={lattice_fermion['minimum_singular_value']:.9g}; "
            f"Schur={lattice_fermion['schur_solve_relative_residual']:.3e}; "
            "curvature="
            f"{lattice_fermion['one_loop_collective_coordinate']['relative_error']:.3e}"
        ),
    )
    check(
        "lattice",
        "Wilson and tree-Symanzik scaling controls have orders two and four",
        1.95 < lattice_scaling["wilson_observed_order_last_pair"] < 2.05
        and 3.90 < lattice_scaling["symanzik_observed_order_last_pair"] < 4.10,
        (
            f"pW={lattice_scaling['wilson_observed_order_last_pair']:.8f}; "
            f"pS={lattice_scaling['symanzik_observed_order_last_pair']:.8f}"
        ),
    )

    lattice_required = (
        "importance_sampling_or_Monte_Carlo_performed",
        "dynamical_4d_lattice_QC2D_phase_established",
        "fermion_determinant_positivity_proven",
        "complete_nonlinear_gauge_meson_ghost_fermion_Hessian_computed",
        "one_loop_renormalised_continuum_limit_established",
        "finite_amplitude_global_stability_established",
    )
    lattice_lane_closed = all(lattice_fail[name] for name in lattice_required)
    check(
        "lattice_boundary",
        "finite diagnostic is not promoted to a dynamical 4D/full-Hessian result",
        not lattice_lane_closed and not any(lattice_fail[name] for name in lattice_required),
        "all six microscopic/full-closure flags remain false",
    )

    check(
        "aps",
        "both reduced bordism generators and reference product spectra are explicit",
        aps["reduced_bordism_generators_explicit"] is True
        and aps["product_dirac_spectra_verified"] is True
        and aps["plain_product_eta_phases"]["G3_S1_periodic_times_S3"] == 1
        and aps["plain_product_eta_phases"]["G2_T2_odd_times_S2"] == 1,
        f"plain phases={aps['plain_product_eta_phases']}",
    )
    regulator_phases = {
        (row["phase_on_G3"], row["phase_on_G2"])
        for row in aps["decorated_defect_regulator_table"]
    }
    check(
        "aps",
        "defect regulators realize all four torsion characters",
        regulator_phases == {(1, 1), (-1, 1), (1, -1), (-1, -1)},
        f"characters={sorted(regulator_phases)}",
    )
    check(
        "aps",
        "charged-QC2D does not yet select either torsion character",
        aps["epsilon_3_selected_by_charged_qc2d_uv"] is False
        and aps["epsilon_2_selected_by_charged_qc2d_uv"] is False
        and aps["aps_torsion_pair_uniquely_determined"] is False,
        "reference/defect calculations establish regulator dependence, not UV selection",
    )
    check(
        "semisimple",
        "the simply connected Sp(4) candidate passes the scanned algebraic screens",
        aps["best_semisimple_candidate"].startswith("Sp(4)")
        and aps["best_candidate_group_theory_pass"] is True
        and aps["best_candidate_tree_level_spectrum_pass"] is True
        and aps["best_candidate_one_loop_asymptotic_freedom_pass"] is True
        and aps["route_e_exactly_two_operator_group_theory_realized"] is True,
        (
            f"candidate={aps['best_semisimple_candidate']}; "
            f"b0={aps['one_loop_beta0']['Sp4']}"
        ),
    )
    aps_required = (
        "epsilon_3_selected_by_charged_qc2d_uv",
        "epsilon_2_selected_by_charged_qc2d_uv",
        "best_candidate_full_dai_freed_gauge_bordism_audit_complete",
        "best_candidate_threshold_decoupling_radiatively_protected",
        "best_candidate_heavy_threshold_eta_matched",
        "best_candidate_nonperturbative_phase_and_soliton_proven",
    )
    aps_lane_closed = all(aps[name] for name in aps_required)
    check(
        "aps_boundary",
        "Sp(4) remains a tree-level candidate rather than an all-scale UV completion",
        not aps_lane_closed and not any(aps[name] for name in aps_required),
        "APS selection, gauge bordism, protected thresholds, eta matching, and strong dynamics are open",
    )

    sqm_status = sqm["status"]
    callias = sqm["metrics"]["callias_pushforward"]
    cpt = sqm["metrics"]["cpt"]
    transgression = sqm["metrics"]["wzw_transgression"]
    check(
        "same_soliton",
        "conditional Callias boundary class gives rank one and determinant c1=+2",
        callias["required"] == {"rank": 1, "determinant_c1": 2}
        and sqm_status["conditional_callias_template_c1_plus_2_verified"] is True,
        f"template={callias['required']}",
    )
    check(
        "same_soliton",
        "Serre-dual anti-canonical CPT sector has three states while O(-2) has one",
        cpt["anti_canonical_minus_degree"] == -4
        and cpt["plus_kernel"] == {"positive": 3, "negative": 0}
        and cpt["minus_kernel"] == {"positive": 0, "negative": 3}
        and cpt["naive_inverse_kernel"] == {"positive": 0, "negative": 1},
        (
            f"plus={cpt['plus_kernel']}; minus={cpt['minus_kernel']}; "
            f"naive={cpt['naive_inverse_kernel']}"
        ),
    )
    check(
        "same_soliton",
        "spatial differential pushforward has degree two and conditional Chern number +2",
        transgression["input_degree"] == 5
        and transgression["fiber_dimension"] == 3
        and transgression["output_degree"] == 2
        and transgression["plus_sector"] == 2,
        (
            f"degree={transgression['input_degree']}-{transgression['fiber_dimension']}="
            f"{transgression['output_degree']}; c1={transgression['plus_sector']}"
        ),
    )
    check(
        "same_soliton",
        "gauge-basic descent conditions are explicit and all remain unsupplied",
        len(sqm["gauge_descent_conditions"]) == 5
        and not any(sqm["gauge_descent_conditions"].values())
        and sqm_status["wzw_gauge_basic_descent_same_model_proven"] is False,
        f"conditions={list(sqm['gauge_descent_conditions'])}",
    )

    sqm_required = (
        "same_soliton_cp1_moduli_derived",
        "worldline_n2_supersymmetry_same_model_derived",
        "physical_tangent_fermion_same_model_derived",
        "callias_operator_same_model_specified",
        "callias_fredholm_gap_same_model_verified",
        "callias_determinant_c1_same_model_computed",
        "cpt_anti_canonical_map_same_model_derived",
        "wzw_equivariant_refinement_same_model_constructed",
        "wzw_gauge_basic_descent_same_model_proven",
        "wzw_pullback_o2_same_model_proven",
        "ap_e3_ap_e4_composition_gate_closed",
    )
    sqm_lane_closed = all(sqm_status[name] for name in sqm_required)
    check(
        "same_soliton_boundary",
        "conditional template is not substituted for a same-mother-model derivation",
        not sqm_lane_closed and not any(sqm_status[name] for name in sqm_required),
        "all eleven same-model closure flags remain false",
    )

    route_closure = {
        "four_dimensional_lattice_and_full_hessian": lattice_lane_closed,
        "aps_and_semisimple_uv": aps_lane_closed,
        "same_soliton_sqm_callias_descent": sqm_lane_closed,
    }
    any_preportal_route_closed = any(route_closure.values())
    portal_prerequisites = {
        "at_least_one_preportal_lane_closed": any_preportal_route_closed,
        "same_moduli_space_identified": sqm_status["same_soliton_cp1_moduli_derived"],
        "uniform_gap_below_portal_scale": False,
        "gauge_basic_wzw_line_o_plus_2": sqm_status["wzw_pullback_o2_same_model_proven"],
        "physical_tangent_or_callias_line_matched": sqm_status[
            "callias_determinant_c1_same_model_computed"
        ],
        "physical_cpt_anti_canonical_map": sqm_status[
            "cpt_anti_canonical_map_same_model_derived"
        ],
        "anomaly_free_portal_operator": False,
        "degree_one_map_and_orientation_proven": False,
    }
    degree_one_portal = {
        "authorized": any_preportal_route_closed,
        "constructed": False,
        "prerequisites": portal_prerequisites,
        "decision": (
            "No pre-portal lane is closed.  The degree-one Route-E portal is "
            "therefore deliberately not constructed; only its theorem-level "
            "prerequisites are recorded for the next audit."
        ),
    }
    check(
        "portal",
        "no lane closure means no degree-one portal construction",
        not any_preportal_route_closed
        and degree_one_portal["authorized"] is False
        and degree_one_portal["constructed"] is False,
        f"routes={route_closure}",
    )
    check(
        "promotion",
        "all source cards and the integrated frontier remain non-promoting",
        lattice_fail["physics_promotion_allowed"] is False
        and aps["physics_promotion_allowed"] is False
        and sqm_status["physics_promotion_allowed"] is False,
        "physics_promotion_allowed=false",
    )

    passed = sum(row["pass"] for row in CHECKS)
    result = {
        "schema_version": "ap-e5-completion-frontier-v1",
        "summary": {"passed": passed, "total": len(CHECKS)},
        "all_pass": passed == len(CHECKS),
        "source_card_counts": {
            "lattice": lattice_summary,
            "aps": {"passed": aps_passed, "total": aps_total},
            "same_soliton": sqm_summary,
        },
        "route_closure": route_closure,
        "any_preportal_route_closed": any_preportal_route_closed,
        "degree_one_portal": degree_one_portal,
        "physics_promotion_allowed": False,
        "source_manifest": source_manifest,
        "checks": CHECKS,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e5_completion_frontier.json"
    md_path = OUTPUT / "ap_e5_completion_frontier.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(markdown(result), encoding="utf-8")
    print(f"AP-E completion frontier: {passed}/{len(CHECKS)} checks pass")
    print("any_preportal_route_closed=false")
    print("degree_one_portal_constructed=false")
    print("physics_promotion_allowed=false")
    return 0 if result["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
