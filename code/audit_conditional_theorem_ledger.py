#!/usr/bin/env python3
"""Build a compact evidence ledger for the conditional Spin(10) theorem.

This is a bookkeeping audit, not a new physics scan.  It reads the local
verification outputs produced by the earlier scans and writes a table that
separates proved/verified conditional claims from tuned fallbacks and open
problems.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "conditional_theorem_ledger"


def read_json(relpath: str) -> dict[str, Any]:
    path = ROOT / relpath
    with path.open() as handle:
        return json.load(handle)


def exists(relpath: str) -> bool:
    return (ROOT / relpath).exists()


def fmt(x: Any, digits: int = 6) -> str:
    if isinstance(x, float):
        return f"{x:.{digits}g}"
    return str(x)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    constrained = read_json("output/constrained_54_210_source_sector/summary.json")
    conormal = read_json("output/combined_conormal_54_210/summary.json")
    schur = read_json("output/schur_auxiliary_component_spectrum/summary.json")
    xi = read_json("output/combined_conormal_multiplier_threshold/summary.json")
    y4pi = read_json("output/yukawa_4pi_audit/yukawa_4pi_audit_summary.json")
    lock = read_json("output/polynomial_schur_locking/polynomial_schur_locking_summary.json")
    boundary_action = read_json("output/boundary_source_action/summary.json")
    nlsm_source = read_json("output/nlsm_composite_source_origin/summary.json")
    nlsm_uv_scorecard = read_json("output/nlsm_uv_completion_scorecard/summary.json")
    dressed_d5 = read_json("output/dressed_dimension5_channels/summary.json")
    spectrum_d5 = read_json("output/spectrum_aware_d5_dressing/summary.json")
    mass_insertion_d5 = read_json("output/mass_insertion_d5_dressing/summary.json")
    eigenstate_d5 = read_json("output/eigenstate_d5_dressing/summary.json")
    mssm_mixing_d5 = read_json("output/mssm_mixing_d5_dressing/summary.json")
    signed_interference_d5 = read_json("output/signed_interference_d5/summary.json")
    chiral_lattice_d5 = read_json("output/chiral_lattice_d5/summary.json")
    tightened_triplet_filter = read_json("output/tightened_triplet_filter/summary.json")
    tightened_downstream_proton = read_json("output/tightened_downstream_proton/summary.json")
    physical_d5_wilson = read_json("output/physical_d5_wilson_replay/summary.json")
    soft_spectrum_d5 = read_json("output/soft_spectrum_d5_replay/summary.json")
    interference_tightened_d5 = read_json("output/interference_tightened_d5_replay/summary.json")
    epi_clebsch = read_json("output/epi_clebsch_replay/summary.json")
    knu_target = read_json("output/knu_target_map/summary.json")
    knu_wilson_gap = read_json("output/knu_wilson_normalization_gap/summary.json")
    full_knu_pipeline = read_json("output/full_knu_channel_pipeline/summary.json")
    full_knu_width = read_json("output/full_knu_width/summary.json")
    knu_nearnull_width = read_json("output/knu_nearnull_width/summary.json")
    nearnull_finite_lift_width = read_json("output/nearnull_finite_lift_width/summary.json")
    threshold_locked_lift = read_json("output/threshold_locked_triplet_lift/summary.json")
    unitary_link_spin10 = read_json("output/unitary_link_spin10_embedding/summary.json")
    unitary_link_nlsm = read_json("output/unitary_link_nlsm_action/summary.json")
    locked_flavor_d5_card = read_json("output/locked_link_full_flavor_d5_card/summary.json")
    full_flavor_d5_pipeline = read_json("output/full_flavor_d5_pipeline/summary.json")
    shadow_contact_joint = read_json("output/shadow_contact_flavor_d5_joint/summary.json")
    two_direction_shadow_triplet = read_json("output/two_direction_shadow_triplet_joint/summary.json")
    two_kernel_flavor_d5 = read_json("output/two_kernel_flavor_then_d5/summary.json")
    clebsch_null_flavor_d5 = read_json("output/clebsch_null_flavor_d5/summary.json")
    source_ckm_crossed120 = read_json("output/source_consistent_ckm_crossed120/summary.json")
    crossed120_source_symmetry = read_json("output/crossed120_action_level_source_symmetry/summary.json")
    triplet_clebsch_nullspace = read_json("output/triplet_clebsch_nullspace/summary.json")
    triplet_120_rrrr = read_json("output/triplet_120_rrrr_mass_matrix/summary.json")
    crossed_120_projector = read_json("output/crossed_120_triplet_projector/summary.json")
    ps_crossed_120_source = read_json("output/ps_crossed_120_source_action/summary.json")
    completed_120_partner = read_json("output/completed_120_partner_action/summary.json")
    crossed_120_link_locking = read_json("output/crossed_120_link_locking_action/summary.json")
    unitary_link_dterm = read_json("output/unitary_link_dterm_quotient/summary.json")
    composite_unitary_link = read_json("output/composite_unitary_link_glsm/summary.json")
    quantum_deformed_link = read_json("output/quantum_deformed_link_moduli/summary.json")
    hidden_radial_lock = read_json("output/hidden_radial_lock_sector/summary.json")
    endpoint_vectorlike = read_json("output/endpoint_vectorlike_completion/summary.json")
    publication_closure_card = read_json("output/publication_closure_card/summary.json")
    publication_repro = read_json("output/publication_flavor_d5_reproducibility/summary.json")
    publication_channel_d5 = read_json("output/publication_channel_d5_tables/summary.json")
    publication_triplet_eigen = read_json("output/publication_triplet_eigenstate_card/summary.json")
    publication_dressed_c5 = read_json("output/publication_dressed_c5_from_eigenstate_card/summary.json")
    publication_dressed_c5_kappa = read_json("output/publication_dressed_c5_kappa_scan/summary.json")
    publication_dressed_c5_rescue = read_json("output/publication_dressed_c5_clockwork_rescue/summary.json")
    clockwork_publication_card = read_json("output/clockwork_rescued_publication_card/clockwork_rescued_publication_card.json")
    no_web_input_ledger = read_json("output/no_web_input_convention_ledger/summary.json")
    no_web_flavor_targets = read_json("output/no_web_flavor_target_provenance/summary.json")
    no_web_pmns_replay = read_json("output/no_web_pmns_benchmark_replay/summary.json")
    source_pmns_replay = read_json("output/source_consistent_pmns_replay/summary.json")
    source_majorana_rank = read_json("output/source_majorana_texture_rank/summary.json")
    majorana_contact_sensitivity = read_json("output/majorana_contact_sensitivity/summary.json")
    majorana_source_locking = read_json("output/majorana_source_locking_sector/summary.json")
    majorana_hidden_quotient = read_json("output/majorana_hidden_quotient_origin/summary.json")
    majorana_monomial_clockwork = read_json("output/majorana_monomial_clockwork_origin/summary.json")
    rank_one_clockwork = read_json("output/rank_one_clockwork_locking/summary.json")
    constrained_clockwork_hessian = read_json("output/constrained_clockwork_source_hessian/summary.json")
    clockwork_unitary_link = read_json("output/clockwork_unitary_link_quotient/summary.json")
    clockwork_hidden_endpoint = read_json("output/clockwork_hidden_endpoint_meson/summary.json")
    clockwork_endpoint_vectorlike = read_json("output/clockwork_endpoint_vectorlike_completion/summary.json")
    clockwork_radial_driver = read_json("output/clockwork_radial_driver_hessian/summary.json")
    clockwork_hidden_quotient = read_json("output/clockwork_hidden_gauge_quotient_origin/summary.json")

    safe_r200 = next(
        row
        for row in xi["rows"]
        if row["R"] == 200.0 and row["kappa_Xi_combined"] == 1.0
    )
    clean_r200 = safe_r200["best"]
    xi_window = xi["allowed_windows"]["PDelta_lt_1e_minus_2"]
    schur_window = next(
        row
        for row in schur["windows"]
        if row["branch"] == "triplet_only_bridge" and row["projected_tol"] == 0.01
    )
    elem_uv = next(
        row
        for row in constrained["uv_rows"]
        if row["scenario"] == "elementary_54_210_plus_10G"
    )
    constr_uv = next(
        row
        for row in constrained["uv_rows"]
        if row["scenario"] == "constrained_54_210_sources_plus_10G"
    )
    nlsm_uv_route_d = next(
        row
        for row in nlsm_uv_scorecard["route_summaries"]
        if row["route"] == "D_cutoff_shared_U_NLSM"
    )
    dressed_d5_worst = dressed_d5["verdict"]["omegaR_0p1_kappa100_wino100TeV_worst"]
    spectrum_d5_pref = spectrum_d5["by_filter_label"]["omegaR_0.1_kappa_100"]
    spectrum_d5_baseline = spectrum_d5_pref["baseline"]
    mass_insertion_pref = mass_insertion_d5["by_filter_label"]["omegaR_0.1_kappa_100"]
    eigenstate_pref = eigenstate_d5["by_filter_label"]["omegaR_0.1_kappa_100"]
    mssm_mixing_pref = mssm_mixing_d5["by_filter_label"]["omegaR_0.1_kappa_100"]
    signed_interference_pref = signed_interference_d5["by_filter_label"]["omegaR_0.1_kappa_100"]
    chiral_lattice_verdict = chiral_lattice_d5["verdict"]
    tightened_verdict = tightened_triplet_filter["verdict"]
    tightened_downstream_verdict = tightened_downstream_proton["verdict"]
    physical_d5_wilson_verdict = physical_d5_wilson["verdict"]
    soft_spectrum_d5_verdict = soft_spectrum_d5["verdict"]
    interference_tightened_d5_verdict = interference_tightened_d5["verdict"]
    epi_clebsch_verdict = epi_clebsch["verdict"]
    knu_target_verdict = knu_target["verdict"]
    knu_wilson_gap_verdict = knu_wilson_gap["verdict"]
    full_knu_pipeline_verdict = full_knu_pipeline["verdict"]
    full_knu_width_verdict = full_knu_width["verdict"]
    full_knu_width_checks = full_knu_width["cross_checks"]
    knu_nearnull_width_verdict = knu_nearnull_width["verdict"]
    finite_lift_verdict = nearnull_finite_lift_width["verdict"]
    threshold_locked_lift_verdict = threshold_locked_lift["verdict"]
    unitary_link_verdict = unitary_link_spin10["verdict"]
    unitary_link_nlsm_verdict = unitary_link_nlsm["verdict"]
    locked_flavor_d5_verdict = locked_flavor_d5_card["verdict"]
    shadow_contact_verdict = shadow_contact_joint["verdict"]
    two_direction_verdict = two_direction_shadow_triplet["verdict"]
    two_kernel_verdict = two_kernel_flavor_d5["verdict"]
    source_ckm_crossed120_verdict = source_ckm_crossed120["verdict"]
    crossed120_source_symmetry_verdict = crossed120_source_symmetry["verdict"]
    triplet_null_verdict = triplet_clebsch_nullspace["verdict"]
    triplet_120_verdict = triplet_120_rrrr["verdict"]
    crossed_120_verdict = crossed_120_projector["verdict"]
    ps_crossed_120_verdict = ps_crossed_120_source["verdict"]
    ps_crossed_120_flat = ps_crossed_120_source["flatness_and_hessian"]
    ps_crossed_120_literal = next(
        row
        for row in ps_crossed_120_source["threshold_interpretation_scan"]
        if row["branch"] == "two_literal_triplet_mediator_pairs"
    )
    completed_120_verdict = completed_120_partner["verdict"]
    completed_120_window = completed_120_partner["mass_locking_window"]
    crossed_120_link_verdict = crossed_120_link_locking["verdict"]
    crossed_120_link_geometry = crossed_120_link_locking["linearized_constraint_geometry"]
    unitary_link_dterm_verdict = unitary_link_dterm["verdict"]
    unitary_link_dterm_f = unitary_link_dterm["Fterm_on_Dflat_tangent"]
    composite_unitary_link_verdict = composite_unitary_link["verdict"]
    quantum_deformed_link_verdict = quantum_deformed_link["verdict"]
    hidden_radial_lock_verdict = hidden_radial_lock["verdict"]
    endpoint_vectorlike_verdict = endpoint_vectorlike["verdict"]
    clockwork_card = rank_one_clockwork["preferred_clockwork_card"]
    constrained_clockwork_verdict = constrained_clockwork_hessian["verdict"]
    constrained_clockwork_la = constrained_clockwork_hessian["linear_algebra"]
    clockwork_unitary_verdict = clockwork_unitary_link["verdict"]
    clockwork_unitary_unit = clockwork_unitary_link["unitary_dilation"]
    clockwork_unitary_f = clockwork_unitary_link["Fterm_on_Dflat_tangent"]
    clockwork_hidden_verdict = clockwork_hidden_endpoint["verdict"]
    clockwork_hidden_rank = clockwork_hidden_endpoint["rank_audit"]
    clockwork_hidden_quantum = clockwork_hidden_endpoint["quantum_rank_audit"]
    clockwork_endpoint_vectorlike_verdict = clockwork_endpoint_vectorlike["verdict"]
    clockwork_radial_verdict = clockwork_radial_driver["verdict"]
    clockwork_radial_rank = clockwork_radial_driver["rank_audit"]
    clockwork_hidden_quotient_verdict = clockwork_hidden_quotient["verdict"]
    clockwork_hidden_quotient_rank = clockwork_hidden_quotient["linearized_quotient_rank"]
    publication_dressed_c5_rescue_verdict = publication_dressed_c5_rescue["verdict"]
    publication_dressed_c5_rescue_exact = publication_dressed_c5_rescue["exact_kappa729_replay"]
    clockwork_publication_card_verdict = clockwork_publication_card["verdict"]
    clockwork_publication_card_d5 = clockwork_publication_card["d5_summary"]
    no_web_input_ledger_verdict = no_web_input_ledger["verdict"]
    no_web_d5_replay = no_web_input_ledger["d5_row_formula_verification"]
    no_web_flavor_verdict = no_web_flavor_targets["verdict"]
    no_web_flavor_scores = no_web_flavor_targets["scores"]
    no_web_pmns_verdict = no_web_pmns_replay["verdict"]
    no_web_pmns_residuals = no_web_pmns_replay["residuals"]
    source_pmns_verdict = source_pmns_replay["verdict"]
    source_pmns_residuals = source_pmns_replay["residuals"]
    source_majorana_rank_verdict = source_majorana_rank["verdict"]
    source_majorana_rank_decomp = source_majorana_rank["decomposition"]
    majorana_contact_verdict = majorana_contact_sensitivity["verdict"]
    majorana_contact_reps = majorana_contact_sensitivity["representative_rows"]
    majorana_contact_s_width = majorana_contact_sensitivity["real_scale_loose_interval"]["half_width"]
    majorana_contact_phase_width = majorana_contact_sensitivity["phase_loose_interval"]["half_width"]
    majorana_source_locking_verdict = majorana_source_locking["verdict"]
    majorana_hidden_quotient_verdict = majorana_hidden_quotient["verdict"]
    majorana_hidden_candidates = majorana_hidden_quotient["candidates"]
    majorana_monomial_verdict = majorana_monomial_clockwork["verdict"]
    majorana_monomial_best = majorana_monomial_clockwork["best_small_denominator_row"]
    majorana_monomial_first = majorana_monomial_clockwork["first_denominator_within_loose_phase_tolerance"]

    rows: list[dict[str, str]] = [
        {
            "claim": "PSLT/Another Physics alone fixes a unique GUT",
            "status": "NO_GO",
            "key_number": "not derivable from present axioms",
            "evidence_path": "paper/gut_framework.tex",
            "interpretation": "The result must be stated as a conditional Spin(10) EFT branch, not an unconditional first-principles derivation.",
        },
        {
            "claim": "Pati-Salam embedding has correct hypercharge and anomalies",
            "status": "PASS",
            "key_number": "k_Y=5/3",
            "evidence_path": "output/pati_salam/pati_salam_report.md",
            "interpretation": "The charge table and anomaly cancellation are already mechanically checked.",
        },
        {
            "claim": "CP1 O(2) gives a protected three-family kinematic space",
            "status": "PASS_AS_KINEMATICS",
            "key_number": "dim H0(CP1,O(2))=3",
            "evidence_path": "output/yukawa_o2/cp1_o2_yukawa_report.md",
            "interpretation": "This is a valid family-counting input, but it is not by itself a complete flavor fit.",
        },
        {
            "claim": "Type-I seesaw plus trace lift is algebraically consistent",
            "status": "PASS",
            "key_number": "trace-lift residual ~2.6e-16",
            "evidence_path": "output/seesaw/seesaw_item3_report.md",
            "interpretation": "The local seesaw matrix identities and trace-lift construction pass the numerical audits.",
        },
        {
            "claim": "Corrected 4pi two-loop branch remains viable",
            "status": "PASS",
            "key_number": (
                f"safe={y4pi['corrected_summary']['scan_diagnostics']['proton_filter_perturbative_safe_points']}, "
                f"M3={fmt(y4pi['corrected_summary']['best']['M_Sigma3_GeV'])} GeV"
            ),
            "evidence_path": "output/yukawa_4pi_audit/yukawa_4pi_audit_summary.json",
            "interpretation": "The branch survives after replacing the Yukawa drag by Tr(Ydag Y)/(4pi).",
        },
        {
            "claim": "Constrained 54/210 sources preserve the projector without UV blow-up",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"54 normals removed={constrained['audit_54']['normal_modes_removed_if_radial_retained']}, "
                f"210 normals removed={constrained['audit_210']['normal_modes_removed_if_radial_retained']}, "
                f"LP/MG={fmt(constr_uv['landau_ratio_Lambda_over_MG'])}"
            ),
            "evidence_path": "output/constrained_54_210_source_sector/summary.json",
            "interpretation": "Works only if F54 and Omega210 are constrained/composite order parameters rather than full elementary towers.",
        },
        {
            "claim": "Elementary propagating 54+210 source tower is UV-safe to R=50,200",
            "status": "FAIL",
            "key_number": f"LP/MG={fmt(elem_uv['landau_ratio_Lambda_over_MG'])}",
            "evidence_path": "output/constrained_54_210_source_sector/summary.json",
            "interpretation": "The Landau pole is below the mediator scales, so this is not the preferred branch.",
        },
        {
            "claim": "Combined conormal constraints leave only shared orbit modes",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"rank={conormal['combined_jacobian']['jacobian_rank']}, "
                f"nullity={conormal['combined_jacobian']['jacobian_nullity']}"
            ),
            "evidence_path": "output/combined_conormal_54_210/summary.json",
            "interpretation": "Normal-bundle multipliers pair the 240 normal directions; the 24 zero modes are gauge-orbit tangents.",
        },
        {
            "claim": "Propagating combined conormal multiplier is a natural threshold",
            "status": "TUNED_FALLBACK",
            "key_number": (
                f"{xi_window['kappa_min']:.6f}<kappa_Xi<{xi_window['kappa_max']:.6f}"
            ),
            "evidence_path": "output/combined_conormal_multiplier_threshold/summary.json",
            "interpretation": "If the conormal multiplier propagates, it must be mass-locked at the percent level for projected threshold <1e-2.",
        },
        {
            "claim": "Schur lock can be polynomial and threshold-silent",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"Fmax={fmt(lock['diagnostics']['max_driver_F_residual'])}, "
                f"projected_l2={fmt(schur['minimal_branch']['projected_threshold_l2'])}, "
                f"LP/MG={fmt(schur['minimal_branch']['landau_ratio_with_baseline'])}"
            ),
            "evidence_path": "output/polynomial_schur_locking/polynomial_schur_locking_summary.json",
            "interpretation": "The minimal singlet plus degenerate 10'_G branch is threshold-safe; triplet-only variants are fallbacks.",
        },
        {
            "claim": "Triplet-only Schur bridge is harmless without mass locking",
            "status": "TUNED_FALLBACK",
            "key_number": (
                f"{schur_window['kappa_min']:.6f}<kappa_D<{schur_window['kappa_max']:.6f}"
            ),
            "evidence_path": "output/schur_auxiliary_component_spectrum/summary.json",
            "interpretation": "The incomplete bridge is allowed only as a scanned EFT fallback.",
        },
        {
            "claim": "Local A5-A6 boundary source action leaves only shared orbit modes",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"rank={boundary_action['hessian']['rank']}, "
                f"nullity={boundary_action['hessian']['nullity']}, "
                f"Nnormal={boundary_action['hessian']['max_null_normal_component']:.2e}"
            ),
            "evidence_path": "output/boundary_source_action/summary.json",
            "interpretation": "This strengthens A5-A6 to a local action-level ansatz, while the UV microscopic origin remains open.",
        },
        {
            "claim": "Shared-U NLSM source removes relative 54/210 orientation modes",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"rank={nlsm_source['jacobian']['shared_U_rank']}, "
                f"constraints={nlsm_source['jacobian']['shared_U_constraint_rank']}, "
                f"Plucker={nlsm_source['constraint_checks']['Omega_plucker_max_residual']:.2e}"
            ),
            "evidence_path": "output/nlsm_composite_source_origin/summary.json",
            "interpretation": "This is a cutoff NLSM/composite realization of the source manifold; UV completion remains open.",
        },
        {
            "claim": "Minimal A-D routes give a microscopic high-R origin for A5-A6",
            "status": "NO_GO",
            "key_number": (
                f"micro={nlsm_uv_scorecard['verdict']['microscopic_high_R_candidate_found']}, "
                f"D_alpha200={nlsm_uv_route_d['R200_alpha_inv_at_RMG']:.3f}"
            ),
            "evidence_path": "output/nlsm_uv_completion_scorecard/summary.json",
            "interpretation": "Elementary, preon-caricature, and deconstructed-link routes do not provide a perturbative microscopic R=50,200 completion; only the cutoff NLSM branch passes conditionally.",
        },
        {
            "claim": "R=200 constrained branch is threshold/proton safe in the displayed scan",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"M3={fmt(clean_r200['M_Sigma3_GeV'])} GeV, "
                f"tau6={fmt(clean_r200['tau_dim6_years'])} yr, "
                f"safe={safe_r200['safe_points']}"
            ),
            "evidence_path": "output/combined_conormal_multiplier_threshold/summary.json",
            "interpretation": "This is the current numerical benchmark for the conditional theorem.",
        },
        {
            "claim": "Monitored dressed d=5 channels pass at S_T=1e-5 in the baseline spectrum",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"worst={dressed_d5_worst['physical_proxy']}, "
                f"STmax={dressed_d5_worst['S_T_max']:.3e}, "
                f"margin={dressed_d5_worst['margin_at_ST_display']:.3e}"
            ),
            "evidence_path": "output/dressed_dimension5_channels/summary.json",
            "interpretation": "The channel-specific wino-dressed proxy passes for the 100 TeV baseline, but an extreme 10 TeV plus 100x dressing stress fails.",
        },
        {
            "claim": "Spectrum-aware d=5 dressing scan preserves the unstressed branch",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"xi1={spectrum_d5_pref['unstressed_safe_points']}/{spectrum_d5_pref['unstressed_total_points']}, "
                f"allstress={spectrum_d5_pref['safe_points']}/{spectrum_d5_pref['total_points']}, "
                f"baseline_STmax={spectrum_d5_baseline['global_S_T_max']:.3e}"
            ),
            "evidence_path": "output/spectrum_aware_d5_dressing/summary.json",
            "interpretation": "The loop-function wino/higgsino proxy is safe across the xi_H=1 spectrum grid; xi_H stress localizes the remaining risk to coherent RRRR higgsino enhancement.",
        },
        {
            "claim": "Explicit mass-insertion d=5 proxy stays within the proton-safe envelope",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"safe={mass_insertion_pref['safe_points']}/{mass_insertion_pref['total_points']}, "
                f"xi_eff_max={mass_insertion_pref['max_effective_xi_H']:.3e}, "
                f"worst_margin={mass_insertion_pref['most_marginal_safe']['worst_margin']:.3e}"
            ),
            "evidence_path": "output/mass_insertion_d5_dressing/summary.json",
            "interpretation": "Replacing scalar xi_H by explicit Hermitian mass-insertion kernels keeps the preferred filter safe, including a democratic off-diagonal stress.",
        },
        {
            "claim": "Positive soft-eigenstate d=5 proxy stays within the proton-safe envelope",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"aligned={eigenstate_pref['aligned_safe_points']}/{eigenstate_pref['aligned_positive_points']}, "
                f"positive={eigenstate_pref['safe_points']}/{eigenstate_pref['positive_points']}, "
                f"worst_margin={eigenstate_pref['most_marginal_safe']['worst_margin']:.3e}"
            ),
            "evidence_path": "output/eigenstate_d5_dressing/summary.json",
            "interpretation": "Diagonalizing positive soft mass matrices before applying pairwise loop kernels keeps the preferred filter safe in the audited eigenstate proxy.",
        },
        {
            "claim": "MSSM chargino/neutralino mixing d=5 proxy preserves the preferred filter",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"aligned={mssm_mixing_pref['aligned_safe_points']}/{mssm_mixing_pref['aligned_positive_points']}, "
                f"positive={mssm_mixing_pref['safe_points']}/{mssm_mixing_pref['positive_points']}, "
                f"worst_margin={mssm_mixing_pref['most_marginal_safe']['worst_margin']:.3e}"
            ),
            "evidence_path": "output/mssm_mixing_d5_dressing/summary.json",
            "interpretation": "Resolving chargino and neutralino eigenstate content in a coherent-absolute proxy leaves the preferred triplet filter proton-safe.",
        },
        {
            "claim": "Signed CP-phase d=5 interference proxy preserves the preferred filter",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"aligned={signed_interference_pref['aligned_safe_points']}/{signed_interference_pref['aligned_positive_points']}, "
                f"positive={signed_interference_pref['safe_points']}/{signed_interference_pref['positive_points']}, "
                f"worst_margin={signed_interference_pref['most_marginal_safe']['worst_margin']:.3e}"
            ),
            "evidence_path": "output/signed_interference_d5/summary.json",
            "interpretation": "Scanning signed chargino/neutralino/higgsino phase combinations and minimal e-pi channel maps keeps the preferred filter just above the displayed proton bound.",
        },
        {
            "claim": "The d=5 filter target is robust under the conservative chiral/lattice envelope",
            "status": "TUNED_FALLBACK",
            "key_number": (
                f"central_STmax={chiral_lattice_verdict['central_global_S_T_max']:.3e}, "
                f"max_width_STmax={chiral_lattice_verdict['max_width_global_S_T_max']:.3e}, "
                f"max_width_margin={chiral_lattice_verdict['max_width_worst_margin']:.3e}"
            ),
            "evidence_path": "output/chiral_lattice_d5/summary.json",
            "interpretation": "The central chiral replay passes S_T=1e-5, but the conservative hadronic-width envelope requires tightening the filter to about 7.7e-6.",
        },
        {
            "claim": "Tightened triplet filter closes the conservative d=5 chiral/lattice envelope",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"S_T=7.5e-6, "
                f"pref_margin={tightened_verdict['preferred_rounded_worst_margin']:.3e}, "
                f"allblock_margin={tightened_verdict['all_blocks_rounded_worst_margin']:.3e}"
            ),
            "evidence_path": "output/tightened_triplet_filter/summary.json",
            "interpretation": "The preferred omegaR=0.1,kappa=100 filter passes all conservative max-width rows at S_T=7.5e-6; the all-block fallback still fails, so the closure is filter-specific.",
        },
        {
            "claim": "Tightened triplet filter remains compatible with downstream RGE/proton benchmark cards",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"R200_tau5={tightened_downstream_verdict['R200_tau_dim5_tightened_years']:.3e} yr, "
                f"margin_2p4e34={tightened_downstream_verdict['R200_current_margin']:.3e}, "
                f"future_margin={tightened_downstream_verdict['R200_future_1e35_margin']:.3e}"
            ),
            "evidence_path": "output/tightened_downstream_proton/summary.json",
            "interpretation": "Replacing S_T=1e-5 by S_T=7.5e-6 leaves the corrected R=200 threshold/proton card alive for the current 2.4e34 yr stress target, but not for a 1e35 yr future stress target.",
        },
        {
            "claim": "Field-basis CKM/PMNS d=5 Wilson replay passes the tightened filter at common dressing",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"identity_max_margin={physical_d5_wilson_verdict['identity_max_width_worst_margin']:.3e}, "
                f"identity_amp_headroom={physical_d5_wilson_verdict['identity_max_width_extra_amplitude_factor']:.3e}, "
                f"near_null_margin={physical_d5_wilson_verdict['near_null_max_width_worst_margin']:.3e}"
            ),
            "evidence_path": "output/physical_d5_wilson_replay/summary.json",
            "interpretation": "At common 100 TeV dressing normalization, physical flavor rotations do not close the tightened S_T=7.5e-6 margin; full soft-spectrum dressing remains open.",
        },
        {
            "claim": "Soft-spectrum d=5 Wilson replay passes the tightened filter on the audited positive grid",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"safe={soft_spectrum_d5_verdict['preferred_max_width_safe_points']}/"
                f"{soft_spectrum_d5_verdict['preferred_max_width_positive_points']}, "
                f"margin={soft_spectrum_d5_verdict['preferred_max_width_worst_current_margin']:.3e}, "
                f"future_margin={soft_spectrum_d5_verdict['preferred_max_width_worst_future_1e35_margin']:.3e}"
            ),
            "evidence_path": "output/soft_spectrum_d5_replay/summary.json",
            "interpretation": "The preferred filter survives explicit chargino/neutralino and positive-sfermion eigenstate dressing at the current stress target, but not a uniform 1e35 yr future stress.",
        },
        {
            "claim": "Worst-phase signed d=5 replay passes the tightened filter on the audited positive grid",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"safe={interference_tightened_d5_verdict['preferred_max_width_safe_points']}/"
                f"{interference_tightened_d5_verdict['preferred_max_width_positive_points']}, "
                f"margin={interference_tightened_d5_verdict['preferred_max_width_worst_current_margin']:.3e}, "
                f"future_margin={interference_tightened_d5_verdict['preferred_max_width_worst_future_1e35_margin']:.3e}"
            ),
            "evidence_path": "output/interference_tightened_d5_replay/summary.json",
            "interpretation": "The preferred filter survives the scanned worst CP phases at the current stress target, so the pass does not rely on destructive interference; the future-stress margin remains insufficient.",
        },
        {
            "claim": "Neutral-pion and identical-up projection removes the minimal e-pi monitor as the current d=5 bottleneck",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"raw_margin={epi_clebsch_verdict['raw_preferred_current_margin']:.3e}, "
                f"projected_margin={epi_clebsch_verdict['combined_preferred_current_margin']:.3e}, "
                f"future_margin={epi_clebsch_verdict['combined_preferred_future_margin']:.3e}"
            ),
            "evidence_path": "output/epi_clebsch_replay/summary.json",
            "interpretation": "Applying the standard e-pi projection sensitivity shifts the current bottleneck from the minimal e-pi monitor to Knu; the future 1e35 yr stress remains open.",
        },
        {
            "claim": "Knu target map localizes the remaining d=5 future-stress requirement",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"margin={knu_target_verdict['current_margin']:.3e}, "
                f"S_T_future={knu_target_verdict['future_ST_target_required']:.3e}, "
                f"amp_supp={knu_target_verdict['future_amplitude_suppression_needed']:.3e}"
            ),
            "evidence_path": "output/knu_target_map/summary.json",
            "interpretation": "The remaining proton-side target is sharply localized to LLLL uud Knu: future 1e35 yr safety needs either a stronger triplet filter or an equivalent Knu amplitude suppression.",
        },
        {
            "claim": "Knu Wilson normalization-gap audit separates flavor rotation from full dressing risk",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"headroom={knu_wilson_gap_verdict['common_dressing_amplitude_headroom']:.3e}, "
                f"gap={knu_wilson_gap_verdict['normalization_gap_to_conservative_target']:.3e}, "
                f"STmax={knu_wilson_gap_verdict['S_T_required_tau_1e35_common_dressing']:.3e}"
            ),
            "evidence_path": "output/knu_wilson_normalization_gap/summary.json",
            "interpretation": "The older mass-basis Wilson tensor table has large common-dressing headroom, so the remaining Knu problem is the bridge to full soft-spectrum dressing and hadronic normalization.",
        },
        {
            "claim": "Full Knu pipeline bridge localizes the d=5 headroom loss",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"future_margin={full_knu_pipeline_verdict['final_future_margin_1e35']:.3e}, "
                f"supp={full_knu_pipeline_verdict['final_future_amplitude_suppression_needed']:.3e}, "
                f"loss={full_knu_pipeline_verdict['largest_headroom_loss_vs_previous']:.3e}"
            ),
            "evidence_path": "output/full_knu_channel_pipeline/summary.json",
            "interpretation": "On a common 1e35 yr axis, the largest loss of Knu amplitude headroom occurs when moving from common field-basis dressing to the explicit preferred soft-spectrum replay.",
        },
        {
            "claim": "Calibrated Knu width formula aligns soft and projected target rows",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"Kspread={full_knu_width_checks['K_dyn_max_relative_spread']:.3e}, "
                f"mismatch={full_knu_width_checks['ratio_mismatch']:.3e}, "
                f"STmax={full_knu_width_verdict['final_S_T_max_future']:.3e}"
            ),
            "evidence_path": "output/full_knu_width/summary.json",
            "interpretation": "The preferred soft-spectrum bottleneck and the projected Knu target row share the same calibrated width constant, so the remaining d=5 gap is physical rather than a convention mismatch.",
        },
        {
            "claim": "Triplet near-null solves the calibrated Knu width target as a conditional branch",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"worst_future={knu_nearnull_width_verdict['worst_future_margin_1e35']:.3e}, "
                f"RRRR={knu_nearnull_width_verdict['min_RRRR_future_margin']:.3e}, "
                f"rank_def={knu_nearnull_width_verdict['rank_deficient']}"
            ),
            "evidence_path": "output/knu_nearnull_width/summary.json",
            "interpretation": "Under channel-wise transfer into the calibrated width formula, the Knu near-null has ample margin including RRRR, but it remains rank-deficient and needs an action-level origin.",
        },
        {
            "claim": "Finite triplet near-null lifts can be proton-safe and Planck-safe only with threshold locking",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"safe={finite_lift_verdict['planck_and_proton_safe_count']}, "
                f"thr_safe={finite_lift_verdict['planck_proton_and_naive_threshold_safe_count']}, "
                f"best={finite_lift_verdict['best_planck_proton_label']}"
            ),
            "evidence_path": "output/nearnull_finite_lift_width/summary.json",
            "interpretation": "Condition-capped finite lifts preserve calibrated d=5 margins below the reduced Planck benchmark, but no row passes a conservative uncancelled threshold-split proxy; the lift must be complete-multiplet/locked.",
        },
        {
            "claim": "Threshold-locked unitary dilation removes the triplet rank and split-threshold obstruction algebraically",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"pass={threshold_locked_lift_verdict['passing_conditional_rows']}, "
                f"best={threshold_locked_lift_verdict['best_label']}, "
                f"PDelta={threshold_locked_lift_verdict['threshold_projected_l2']:.1e}"
            ),
            "evidence_path": "output/threshold_locked_triplet_lift/summary.json",
            "interpretation": "The Julia/Sz.-Nagy dilation realizes the near-null as a visible inverse-propagator subblock of a degenerate 8x8 mass matrix; a Spin(10) representation and symmetry origin remain required.",
        },
        {
            "claim": "Perturbative elementary Spin(10) 10-copy unitary link closes the triplet sector",
            "status": "NO_GO",
            "key_number": (
                f"sym_fail={unitary_link_verdict['single_species_10_mass_fails_all_targets']}, "
                f"alpha200={unitary_link_verdict['elementary_spin10_two_species_alpha_inv_R200']:.3e}, "
                f"PDelta_trip={unitary_link_verdict['triplet_only_projected_threshold_l2']:.3e}"
            ),
            "evidence_path": "output/unitary_link_spin10_embedding/summary.json",
            "interpretation": "A single complete-10 species fails the symmetric mass requirement, two species lose R=200 perturbativity, and triplet-only remnants fail thresholds; only post-GUT complete-pair or constrained/composite branches remain.",
        },
        {
            "claim": "Constrained unitary-link NLSM gives a nonempty threshold-silent triplet-link origin",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"nonempty={unitary_link_nlsm_verdict['all_targets_nonempty']}, "
                f"moduli={unitary_link_nlsm_verdict['residual_moduli_real_dimension']}, "
                f"margin={unitary_link_nlsm_verdict['minimum_future_margin_1e35']:.3e}"
            ),
            "evidence_path": "output/unitary_link_nlsm_action/summary.json",
            "interpretation": "Treating the unitary block as a constrained NLSM/source manifold gives full-rank fixed-block constraints, 32 real residual moduli, and threshold-silent bookkeeping; it remains conditional rather than elementary holomorphic Spin(10).",
        },
        {
            "claim": "Locked-link flavor and d=5 reproducibility card localizes the remaining bottleneck",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"ckm={locked_flavor_d5_verdict['ckm_fit_completed']}, "
                f"d5cur={locked_flavor_d5_verdict['d5_current_bound_passes']}, "
                f"supp={locked_flavor_d5_verdict['d5_future_suppression_needed']:.3e}"
            ),
            "evidence_path": "output/locked_link_full_flavor_d5_card/summary.json",
            "interpretation": "The exact local Yukawa/seesaw matrices, locked unitary link, Wilson replays, and Knu width calibration are now in one manifest; CKM and future-stress d=5 safety remain open.",
        },
        {
            "claim": "Full flavor plus dressed d=5 pipeline is publication-complete",
            "status": "NO_GO",
            "key_number": (
                f"pub={full_flavor_d5_pipeline['closure']['publication_level_complete']}, "
                f"strict={full_flavor_d5_pipeline['closure']['operator_strict_closures']}, "
                f"loose={full_flavor_d5_pipeline['closure']['operator_loose_closures']}, "
                f"bestKnu={full_flavor_d5_pipeline['numerical_verdict']['best_operator_future_margin']:.3e}"
            ),
            "evidence_path": "output/full_flavor_d5_pipeline/summary.json",
            "interpretation": "Putting the reproducibility card, operator-level flavor scan, and dressed Knu target on one closure ledger shows that the current branch is not publication-complete: no strict or loose operator-level closure survives the future d=5 stress.",
        },
        {
            "claim": "Shadow/contact deformation has a joint CKM-Knu Pareto direction",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"joint={shadow_contact_verdict['joint_improver_count']}, "
                f"safeCKM={shadow_contact_verdict['future_safe_with_ckm_gain_count']}, "
                f"bestCKM={shadow_contact_verdict['best_ckm']['ckm_score']:.3e}"
            ),
            "evidence_path": "output/shadow_contact_flavor_d5_joint/summary.json",
            "interpretation": "A small transvectant/contact deformation can improve CKM magnitude and the Knu proxy together, but no sampled point also satisfies the uniform 1e35 yr Knu future-stress target with CKM gain.",
        },
        {
            "claim": "Minimal two-direction K plus triplet-only shadow scan closes flavor and d=5 together",
            "status": "NO_GO",
            "key_number": (
                f"closures={two_direction_verdict['local_closure_count']}, "
                f"contact={two_direction_verdict['contact_closure_count']}, "
                f"linear={two_direction_verdict['linear_cancel_closure_count']}"
            ),
            "evidence_path": "output/two_direction_shadow_triplet_joint/summary.json",
            "interpretation": "The conservative contact-K branch and the diagnostic minimum-norm linear-cancel branch both fail the simultaneous CKM<1e-3, mass-tolerance, and future Knu-margin targets in this reduced audit.",
        },
        {
            "claim": "Operator-level two-kernel finite-Clebsch ansatz closes flavor and d=5 together",
            "status": "NO_GO",
            "key_number": (
                f"strict={two_kernel_verdict['strict_closure_count']}, "
                f"loose={two_kernel_verdict['loose_closure_count']}, "
                f"bestCKM={two_kernel_verdict['best_ckm']['ckm_score']:.3e}, "
                f"bestKnu={two_kernel_verdict['best_knu']['future_margin']:.3e}"
            ),
            "evidence_path": "output/two_kernel_flavor_then_d5/summary.json",
            "interpretation": "Deriving triplet tensors from the same CP1 Veronese/contact operators improves CKM to the edge of the strict target, but all finite Clebsch profiles remain below the 1e35 yr Knu future margin.",
        },
        {
            "claim": "Crossed-120 Clebsch-null bridge closes the d=5 gap without closing strict flavor",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"loose={clebsch_null_flavor_d5['verdict']['loose_conditional_closure_count']}, "
                f"strict={clebsch_null_flavor_d5['verdict']['strict_conditional_closure_count']}, "
                f"Knu={clebsch_null_flavor_d5['verdict']['source_consistent_kappa30_margin']:.3e}, "
                f"ckmfac={clebsch_null_flavor_d5['verdict']['ckm_factor_to_strict']:.3e}"
            ),
            "evidence_path": "output/clebsch_null_flavor_d5/summary.json",
            "interpretation": "Accepting the post-Spin(10)/constrained crossed 120A/120B triplet projector gives huge d=5 future-stress margin for the source-consistent row, but the field-only Spin(10) projector is forbidden and strict CKM remains open.",
        },
        {
            "claim": "Source-consistent doublet-only CKM repair closes strict flavor with frozen crossed-120 d=5 projector",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"closures={source_ckm_crossed120_verdict['strict_conditional_closure_count']}, "
                f"bestCKM={source_ckm_crossed120_verdict['best_mass_valid_ckm']['ckm_score']:.3e}, "
                f"mass={source_ckm_crossed120_verdict['best_mass_valid_ckm']['mass_score']:.3e}, "
                f"Knu={source_ckm_crossed120_verdict['best_mass_valid_ckm']['future_margin']:.3e}"
            ),
            "evidence_path": "output/source_consistent_ckm_crossed120/summary.json",
            "interpretation": "Holding H,F,G_A,G_B and the crossed 120 triplet projector fixed, a local doublet-Higgs mixing refit reaches CKM<1e-3, mass<0.2, seesaw residual<1e-10, and huge d=5 margin.  This closes the local flavor+d=5 branch only conditionally; the crossed projector still needs an action-level post-Spin(10)/constrained origin.",
        },
        {
            "claim": "Post-Spin(10) source symmetry makes the crossed-120 flavor+d=5 closure action-level compatible",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"ledger={crossed120_source_symmetry_verdict['operator_ledger_passes']}, "
                f"branch={crossed120_source_symmetry_verdict['post_spin10_source_symmetry_branch_viable']}, "
                f"CKM={crossed120_source_symmetry_verdict['best_ckm_score']:.3e}, "
                f"Knu={crossed120_source_symmetry_verdict['best_d5_margin']:.3e}, "
                f"nullity={crossed120_source_symmetry_verdict['holomorphic_link_nullity_real']}"
            ),
            "evidence_path": "output/crossed120_action_level_source_symmetry/summary.json",
            "interpretation": "A PS fragment/source grading with an R-selection rule allows the required crossed triplet source, finite-lift reverse entry, physical doublet Yukawa/mixing sector, and inert threshold-completing partners while forbidding wrong triplet Yukawas, same-source bare masses, inert Yukawas, and triplet-doublet off-block mixing.  It remains conditional because the selection is post-breaking/spurion-assisted and the link unitarity still requires an NLSM/D-term or composite origin.",
        },
        {
            "claim": "Component triplet Clebsch nullspace contains an antisymmetric QQ Knu null",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"fixedQL={triplet_null_verdict['fixed_QL_all_Knu_null_exists']}, "
                f"rank1={triplet_null_verdict['rank_one_bilinear_null_found']}, "
                f"antiQQ={triplet_null_verdict['rank_one_null_is_antisymmetric_QQ_only']}"
            ),
            "evidence_path": "output/triplet_clebsch_nullspace/summary.json",
            "interpretation": "At component level, choosing QQT in the pure antisymmetric 120-like direction kills the LLLL Knu tensor because sym(Y_QQ)=0; this gives a concrete triplet-sector target, but RRRR channels and the Higgs triplet mass matrix still need derivation.",
        },
        {
            "claim": "Crossed 120A/120B triplet projector has a component-level LLLL+RRRR null",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"fixedR={triplet_120_verdict['fixed_anti_RRRR_exact_q_null']}, "
                f"joint={triplet_120_verdict['joint_rank_one_exact_LLLL_RRRR_null']}, "
                f"res={triplet_120_rrrr['linear_algebra']['joint_rank_one_residual']:.1e}"
            ),
            "evidence_path": "output/triplet_120_rrrr_mass_matrix/summary.json",
            "interpretation": "With QQT/UET in the antisymmetric G_A direction and QLT/UDT nearly pure orthogonal G_B, the monitored LLLL Knu and RRRR uusd maps vanish at component level; finite mass-matrix lift leakage remains to be derived and scanned.",
        },
        {
            "claim": "Field-only unbroken Spin(10) grading enforces the crossed 120 triplet projector",
            "status": "NO_GO",
            "key_number": (
                f"witnesses={crossed_120_projector['field_only_grading_no_go']['crossed_only_witness_count']}, "
                f"ZN<={crossed_120_projector['field_only_grading_no_go']['searched_ZN_up_to']}"
            ),
            "evidence_path": "output/crossed_120_triplet_projector/summary.json",
            "interpretation": "If both 16 16 120_A and 16 16 120_B are field-only invariant, the two 120 copies have identical charge, so a field-only unbroken Spin(10) grading cannot allow only the crossed mass entry.",
        },
        {
            "claim": "Post-breaking constrained crossed-120 projector is threshold-silent and d=5 safe in the replay",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"ps={crossed_120_verdict['ps_eft_or_constrained_projector_possible']}, "
                f"leak_safe={crossed_120_verdict['all_finite_leakage_rows_future_safe']}, "
                f"kappaPlanck={crossed_120_verdict['direct_split_planck_safe_kappas']}"
            ),
            "evidence_path": "output/crossed_120_triplet_projector/summary.json",
            "interpretation": "As a post-Spin(10)-breaking Pati-Salam/constrained triplet projector, the crossed 120 block has a degenerate unitary-dilation realization with zero projected threshold and very large replayed 1e35 yr d=5 margins.",
        },
        {
            "claim": "Post-breaking PS crossed-120 source action realizes the projector without doublet mixing",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"F0={ps_crossed_120_flat['F_norm_at_zero']:.1e}, "
                f"inv={ps_crossed_120_flat['triplet_inverse_subblock_residual_fro']:.1e}, "
                f"rhoD={ps_crossed_120_verdict['rho_max_doublet_preserved_1e_minus_6']:.1e}"
            ),
            "evidence_path": "output/ps_crossed_120_source_action/summary.json",
            "interpretation": "At the PS component-action level, a quadratic constrained source with a Julia-dilated triplet block is F-flat at the origin, realizes the crossed inverse block exactly, and leaves the doublet block diagonal up to the audited leakage tolerance.",
        },
        {
            "claim": "Literal triplet-only crossed-120 mediator completion is threshold silent",
            "status": "NO_GO",
            "key_number": (
                f"PDelta={ps_crossed_120_literal['projected_l2']:.3e}, "
                f"ratio={ps_crossed_120_verdict['literal_triplet_projected_l2_over_budget']:.1f}"
            ),
            "evidence_path": "output/ps_crossed_120_source_action/summary.json",
            "interpretation": "If the Julia defect fields are interpreted as literal propagating triplet-only PS mediator pairs at M_lock, their nonuniversal threshold overshoots the R=200 budget; they must be constrained/auxiliary or completed by inert doublet partners.",
        },
        {
            "claim": "Inert completed 5+5bar partner package repairs the crossed-120 threshold",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"silent={completed_120_verdict['exact_completed_pair_threshold_silent']}, "
                f"Z2={completed_120_verdict['inert_grading_decouples_doublets']}, "
                f"LP={completed_120_verdict['minimal_yukawa_plus_two_10G_R200_safe']}"
            ),
            "evidence_path": "output/completed_120_partner_action/summary.json",
            "interpretation": "Adding inert doublet partners to form two mass-locked 5+5bar packages cancels the one-loop nonuniversal threshold; a Z2 inert grading forbids their matter Yukawa and physical-doublet mixing operators.",
        },
        {
            "claim": "Inert completed 5+5bar partners can be freely split from the triplet mass",
            "status": "NO_GO",
            "key_number": (
                f"|logxi|<{completed_120_window['max_abs_log_xi']:.3e}, "
                f"window={completed_120_window['percent_window']:.3f}%"
            ),
            "evidence_path": "output/completed_120_partner_action/summary.json",
            "interpretation": "The completed package is threshold silent only when doublet partners are mass-locked to the triplets at the per-mille level for the R=200 budget; percent-level splitting already fails the projected threshold test.",
        },
        {
            "claim": "Common-spurion crossed-120 link action derives the inverse block and shared mass scale",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"inv={crossed_120_link_verdict['holomorphic_constraints_realize_inverse_block']}, "
                f"S={crossed_120_link_verdict['common_spurion_locks_triplet_and_doublet_scale']}, "
                f"rank={crossed_120_link_geometry['real_rank']}"
            ),
            "evidence_path": "output/crossed_120_link_locking_action/summary.json",
            "interpretation": "The holomorphic link constraints AB=1 and PBP=W_kappa derive the crossed inverse block, while a common S_lock spurion ties the triplet link and inert doublet mass scale; the branch remains conditional on a unitary/NLSM lock.",
        },
        {
            "claim": "Pure holomorphic crossed-120 link F-terms are sufficient for threshold-safe mass locking",
            "status": "NO_GO",
            "key_number": (
                f"nullity={crossed_120_link_verdict['constraint_nullity_real']}, "
                f"fail_amp={crossed_120_link_verdict['smallest_deformation_failing_lock_window']}"
            ),
            "evidence_path": "output/crossed_120_link_locking_action/summary.json",
            "interpretation": "The holomorphic constraints leave residual moduli; exact deformations preserving AB=1 and PBP=W_kappa can split the triplet singular values beyond the completed-partner mass-locking window, so a D-term/Kahler/NLSM or composite unitarity constraint is required.",
        },
        {
            "claim": "D-term/Kahler quotient removes dangerous crossed-120 link moduli",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"dim={unitary_link_dterm_verdict['Dflat_fixed_block_moduli_real_dimension']}, "
                f"locked={unitary_link_dterm_verdict['all_finite_unitary_completions_lock_singular_values']}, "
                f"rank={unitary_link_dterm_f['fixed_block_rank_on_Dflat_tangent']}"
            ),
            "evidence_path": "output/unitary_link_dterm_quotient/summary.json",
            "interpretation": "Imposing B^dagger B=1 before PBP=W_kappa removes the nonunitary holomorphic moduli; the remaining unitary completions preserve the visible block and keep all singular values locked, so the completed partner threshold remains silent in this local quotient ansatz.",
        },
        {
            "claim": "Hidden U(4) composite GLSM can realize the unitary crossed-120 link locally",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"Dmax={composite_unitary_link_verdict['max_D_residual_2norm']:.3e}, "
                f"logmax={composite_unitary_link_verdict['max_abs_log_singular']:.3e}, "
                f"Nf4={composite_unitary_link_verdict['minimal_Nf4_phase']}, "
                f"Nf8={composite_unitary_link_verdict['conformal_Nf8_phase']}"
            ),
            "evidence_path": "output/composite_unitary_link_glsm/summary.json",
            "interpretation": "A visible-singlet hidden U(4)_H meson branch B=Qtilde Q/f^2 reproduces the unitary link with zero projected visible threshold and no UV Landau obstruction for the asymptotically-free/conformal hidden choices; it remains a constrained composite-sector assumption rather than a completed nonperturbative proof.",
        },
        {
            "claim": "Nf=Nc quantum-deformed hidden moduli are compatible with the unitary link",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"res={quantum_deformed_link_verdict['max_quantum_constraint_abs_residual']:.3e}, "
                f"|Baryon|max={quantum_deformed_link_verdict['max_baryon_abs']:.3g}, "
                f"thr={quantum_deformed_link_verdict['visible_threshold_vector']}"
            ),
            "evidence_path": "output/quantum_deformed_link_moduli/summary.json",
            "interpretation": "For every sampled unitary completion, baryon vevs solve det m - beta beta_tilde=(Lambda/f)^8 at numerical precision; because these are hidden singlets, no visible threshold is induced.",
        },
        {
            "claim": "Nf=Nc quantum deformation alone enforces the unitary crossed-120 link",
            "status": "NO_GO",
            "key_number": (
                f"holo_null={quantum_deformed_link_verdict['holomorphic_real_nullity_after_constraints']}, "
                f"unitary_null={quantum_deformed_link_verdict['unitary_plus_baryon_real_nullity']}, "
                f"removed={quantum_deformed_link_verdict['nonunitary_real_moduli_removed_by_Dterm_lock']}"
            ),
            "evidence_path": "output/quantum_deformed_link_moduli/summary.json",
            "interpretation": "The quantum-deformed constraint has complex nullity 13 after the fixed visible block, so it is compatible with the D-term/radial unitary lock but cannot replace it.",
        },
        {
            "claim": "Hidden radial D-term moose dynamically realizes the crossed-120 unitary lock",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"Drank={hidden_radial_lock_verdict['combined_D_rank']}, "
                f"orbit={hidden_radial_lock_verdict['hidden_gauge_orbit_rank']}, "
                f"linkdim={hidden_radial_lock_verdict['residual_after_radial_D_and_hidden_quotient']}, "
                f"moduli={hidden_radial_lock_verdict['residual_after_fixed_block']}"
            ),
            "evidence_path": "output/hidden_radial_lock_sector/summary.json",
            "interpretation": "Endpoint/radial moment maps make Q and Qtilde unitary up to the hidden U(4)_H frame; quotienting leaves the 16-real-dimensional unitary link, and PBP=W_kappa leaves the expected 9 completion moduli with zero visible threshold.",
        },
        {
            "claim": "Raw endpoint gauging of the hidden radial link is anomaly-free",
            "status": "NO_GO",
            "key_number": f"raw_anom={endpoint_vectorlike_verdict['raw_endpoint_gauging_is_anomalous']}",
            "evidence_path": "output/endpoint_vectorlike_completion/summary.json",
            "interpretation": "Gauging U(4)_L x U(4)_H x U(4)_R with only Q and Qtilde leaves nonzero SU(4)_L and SU(4)_R cubic anomaly proxies, so endpoint gauging requires a vectorlike completion.",
        },
        {
            "claim": "Vectorlike endpoint completion is anomaly-free, beta-safe, and threshold-silent",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"anom={endpoint_vectorlike_verdict['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}, "
                f"beta={endpoint_vectorlike_verdict['vectorlike_completion_beta_safe_to_R200']}, "
                f"thr={endpoint_vectorlike_verdict['vectorlike_completion_visible_threshold_vector']}"
            ),
            "evidence_path": "output/endpoint_vectorlike_completion/summary.json",
            "interpretation": "Adding Qc and Qtilde_c makes every hidden endpoint factor vectorlike, gives b_L=-8, b_H=-4, b_R=-8 at one loop, and induces no visible Spin(10) threshold because all completion fields are visible singlets.",
        },
        {
            "claim": "Single-card source-consistent crossed-120 flavor+d=5 candidate is reproducible",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"local={publication_repro['verdict']['local_source_consistent_candidate_passes']}, "
                f"ckm={publication_repro['fresh_recompute']['scores']['ckm_score']:.3e}, "
                f"mass={publication_repro['fresh_recompute']['scores']['mass_score']:.3e}, "
                f"margin={publication_closure_card['crossed_projector']['future_margin_1e35']:.3e}"
            ),
            "evidence_path": "output/publication_flavor_d5_reproducibility/summary.json",
            "interpretation": "The best source-consistent crossed-120 row is now promoted to a single machine-readable card and recomputed from its own Yukawa matrices.  It passes the local strict CKM, mass, seesaw, post-Spin(10) source/link, and 1e35 yr d=5 margin gates.  It is still not publication-complete because the final channel-specific proton tables and cited input manifest must be regenerated from this exact card.",
        },
        {
            "claim": "Publication-card d=5 channel table scaffold is generated",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"rows={publication_channel_d5['verdict']['proxy_channel_rows']}, "
                f"exact={publication_channel_d5['verdict']['exact_channel_wilson_complete']}, "
                f"margin={publication_channel_d5['verdict']['worst_class_level_margin_1e35']:.3e}"
            ),
            "evidence_path": "output/publication_channel_d5_tables/summary.json",
            "interpretation": "The publication card now regenerates K+nu, K0mu, e+pi0, and mu+pi0 channel rows from its embedded source tensors and flavor rotations.  This is an honest proxy/class-level scaffold rather than the final d=5 theorem, because exact triplet mass-eigenstate tensors, SUSY dressing matrices, and chiral reduction operators remain absent.",
        },
        {
            "claim": "Publication-card triplet eigenstate card exports source-basis C5 tensors",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"cond={publication_triplet_eigen['verdict']['finite_block_condition_number']:.3g}, "
                f"tau0={publication_triplet_eigen['verdict']['finite_kappa_min_tau_years_unfiltered']:.3e}, "
                f"marginST={publication_triplet_eigen['verdict']['finite_kappa_min_margin_1e35_at_reference_filter']:.3e}"
            ),
            "evidence_path": "output/publication_triplet_eigenstate_card/summary.json",
            "interpretation": "The crossed triplet inverse block is now diagonalized as a finite source-basis matrix and contracted into explicit C5L/C5R proxy tensors.  With S_T=1e-5 the worst finite-kappa source-basis proxy margin is large, but the row remains conditional because final SUSY dressing and chiral reductions are not yet included.",
        },
        {
            "claim": "Finite kappa=30 publication-card eigenstate C5 replay is sufficient",
            "status": "NO_GO",
            "key_number": (
                f"central_unsafe={publication_dressed_c5['by_block_case']['finite']['central']['unsafe_1e35_rows']}, "
                f"max_unsafe={publication_dressed_c5['by_block_case']['finite']['max_width']['unsafe_1e35_rows']}, "
                f"STmax={publication_dressed_c5['by_block_case']['finite']['central']['global_S_T_max_1e35']:.3e}"
            ),
            "evidence_path": "output/publication_dressed_c5_from_eigenstate_card/summary.json",
            "interpretation": "The scalar leakage shortcut is removed: channel rows are dressed directly from the exported eigenstate C5 tensors.  The finite kappa=30 branch fails at S_T=1e-5 on the full soft/chiral stress grid, so it is a no-go for the publication branch unless replaced by a tighter or near-rank-one completion.",
        },
        {
            "claim": "Dressed C5 finite-block failure can be converted into a clockwork rank-one lock",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"kmin={publication_dressed_c5_kappa['verdict']['max_width_min_kappa_pass_1e35']:.0f}, "
                f"q={clockwork_card['q']:.1f}, n={clockwork_card['n']}, "
                f"kappa={clockwork_card['kappa']:.0f}, "
                f"margin={clockwork_card['max_width_margin_1e35_at_ST_1e_minus_5']:.3g}"
            ),
            "evidence_path": "output/rank_one_clockwork_locking/summary.json",
            "interpretation": "The near-rank-one requirement kappa>=300 can be generated by a finite graded link chain, for example q=3,n=6 gives kappa=729.  This is only a conditional closure: if implemented as propagating complete 10+10bar chains the estimated Landau-pole ratio is only about 12.6, so the preferred realization must be constrained, composite, or otherwise threshold-silent.",
        },
        {
            "claim": "Clockwork-rescued publication-card dressed C5 replay passes the full stress grid",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"kappa={publication_dressed_c5_rescue['clockwork_card']['kappa']:.0f}, "
                f"central_margin={publication_dressed_c5_rescue_verdict['central_margin_1e35_at_ST_1e_minus_5']:.3g}, "
                f"max_margin={publication_dressed_c5_rescue_verdict['max_width_margin_1e35_at_ST_1e_minus_5']:.3g}, "
                f"unsafe={publication_dressed_c5_rescue_exact['max_width_unsafe_1e35_rows']}"
            ),
            "evidence_path": "output/publication_dressed_c5_clockwork_rescue/summary.json",
            "interpretation": "Recomputing the exact dressed eigenstate-card channel grid at q=3,n=6, kappa=729 gives zero unsafe rows in the central, max-width, and min-width stress cases.  This closes the local d=5 proxy only conditionally on the hidden clockwork/Kahler/FI quotient and still needs the final paper-grade input manifest.",
        },
        {
            "claim": "Single clockwork-rescued flavor+d5 manifest is locally complete",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"flavor={clockwork_publication_card['gates']['flavor_local_pass']}, "
                f"d5={clockwork_publication_card['gates']['d5_all_channel_worst_rows_pass']}, "
                f"rows={clockwork_publication_card_d5['total_channel_rows']}, "
                f"worst={clockwork_publication_card_d5['global_worst']['worst_margin_1e35']:.3g}"
            ),
            "evidence_path": "output/clockwork_rescued_publication_card/clockwork_rescued_publication_card.json",
            "interpretation": "The source-consistent flavor observables, seesaw replay, hidden clockwork quotient, and kappa=729 dressed channel rows now live in one hash-locked local manifest.  This closes the cache-synchronization problem but remains conditional and not publication-final.",
        },
        {
            "claim": "No-web input convention ledger reproduces the dressed d=5 card arithmetic",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"constants={no_web_input_ledger['constant_row_count']}, "
                f"rows={no_web_d5_replay['rows_checked']}, "
                f"tau_relerr={no_web_d5_replay['max_tau_relative_error']:.1e}, "
                f"digest={no_web_input_ledger['constant_rows_sha256'][:12]}"
            ),
            "evidence_path": "output/no_web_input_convention_ledger/summary.json",
            "interpretation": "Every local flavor target, seesaw benchmark, RGE input, width prefactor, channel bound, soft spectrum point, and legacy proton constant used by the current card is now recorded in a hash-locked no-web convention ledger.  Replaying the lifetime formula over all 21168 dressed channel rows gives zero mismatches and zero displayed relative error.  This narrows the remaining publication d=5 blocker to replacing local conventions by a cited input table, not fixing an internal arithmetic inconsistency.",
        },
        {
            "claim": "Constrained clockwork source Hessian has no unintended light tower",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"raw_zero={constrained_clockwork_la['raw_zero_modes_per_chain']}, "
                f"closed_zero={constrained_clockwork_la['boundary_closed_zero_modes_per_chain']}, "
                f"min_closed={constrained_clockwork_la['boundary_closed_min_abs_eigenvalue']:.3g}, "
                f"condW={constrained_clockwork_hessian['effective_inverse_block']['condition_number']:.0f}"
            ),
            "evidence_path": "output/constrained_clockwork_source_hessian/summary.json",
            "interpretation": "The q=3,n=6 link Hessian has exactly one intended zero mode per chain; adding the boundary/source driver removes all clockwork zero modes and leaves masses between 1 and 3.925 M_lock.  This supports the constrained/composite source interpretation, while leaving the microscopic Spin(10) origin open.",
        },
        {
            "claim": "Clockwork W_kappa embeds into a unitary-link quotient with locked triplet masses",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"dimB={clockwork_unitary_unit['dimension']}, "
                f"unit={clockwork_unitary_unit['unitarity_residual_2norm']:.1e}, "
                f"spread={clockwork_unitary_link['threshold_interpretation']['physical_triplet_mass_singular_spread']:.1e}, "
                f"moduli={clockwork_unitary_f['quotient_residual_unitary_moduli_real_dimension']}"
            ),
            "evidence_path": "output/clockwork_unitary_link_quotient/summary.json",
            "interpretation": "The q=3,n=6 clockwork block W=diag(1,epsilon,epsilon,epsilon) admits an 8x8 Julia-Halmos unitary dilation with PBP=W and all physical triplet mass singular values equal to M_lock.  The generalized 4x4 block leaves 33 real unitary-completion moduli, but sampled completions keep the mass spectrum degenerate; microscopic boundary/strong-sector origin remains open.",
        },
        {
            "claim": "Clockwork unitary quotient has a hidden endpoint meson realization",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"D_rank={clockwork_hidden_rank['Dterm']['combined_D_rank']}, "
                f"fixed_rank={clockwork_hidden_rank['fixed_block']['fixed_block_rank_on_unitary_link']}, "
                f"moduli={clockwork_hidden_rank['fixed_block']['residual_unitary_completion_moduli']}, "
                f"qnull={clockwork_hidden_quantum['holomorphic']['complex_nullity']}"
            ),
            "evidence_path": "output/clockwork_hidden_endpoint_meson/summary.json",
            "interpretation": "The 8x8 clockwork link can be represented as a hidden visible-singlet meson with zero direct visible threshold vector.  The radial/D-term lock reproduces the 33 unitary completion moduli, while the Nf=Nc quantum-deformed determinant constraint is only compatible, not sufficient: it leaves 49 complex holomorphic flat directions and therefore does not derive Bdag B=1.",
        },
        {
            "claim": "8x8 endpoint gauging has a vectorlike anomaly/beta completion",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"raw_anom={clockwork_endpoint_vectorlike_verdict['raw_endpoint_gauging_is_anomalous']}, "
                f"vec_anom={clockwork_endpoint_vectorlike_verdict['vectorlike_completion_cancels_all_hidden_cubic_anomalies']}, "
                f"beta_safe={clockwork_endpoint_vectorlike_verdict['vectorlike_completion_beta_safe_to_R200']}, "
                f"thr={clockwork_endpoint_vectorlike_verdict['vectorlike_completion_visible_threshold_vector']}"
            ),
            "evidence_path": "output/clockwork_endpoint_vectorlike_completion/summary.json",
            "interpretation": "Raw SU(8)_L x SU(8)_H x SU(8)_R endpoint gauging is chiral, but adding Qc and Qtilde_c cancels all hidden cubic anomaly proxies.  The vectorlike-completed one-loop coefficients are b=(-16,-8,-16), beta-safe to R=200 in the tested range, and all endpoint fields are visible Spin(10) singlets.",
        },
        {
            "claim": "8x8 radial/Kahler driver lifts nonunitary endpoint modes",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"rank={clockwork_radial_rank['full_driver_rank_with_vectorlike_partners']}, "
                f"null={clockwork_radial_rank['full_driver_nullity']}, "
                f"gauge={clockwork_radial_rank['hidden_gauge_orbit_rank']}, "
                f"phys={clockwork_radial_verdict['physical_residual_moduli_after_gauge_quotient']}"
            ),
            "evidence_path": "output/clockwork_radial_driver_hessian/summary.json",
            "interpretation": "The finite-rank D-term/Kahler driver gives 415 massive real directions out of 512.  The 97 raw flat directions decompose into 64 hidden gauge orbits plus the expected 33 unitary-completion moduli, while vectorlike partners are massive and the visible threshold vector remains zero.  This closes the local Hessian, not the microscopic origin of the Kahler potential.",
        },
        {
            "claim": "8x8 radial equations derive from hidden gauge-quotient moment maps",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"mu_rank={clockwork_hidden_quotient_rank['moment_map_rank']}, "
                f"fixed_rank={clockwork_hidden_quotient_rank['moment_plus_fixed_rank']}, "
                f"TrFI={clockwork_hidden_quotient_verdict['U1_FI_trace_safe']}, "
                f"thr={clockwork_hidden_quotient_verdict['visible_threshold_vector']}"
            ),
            "evidence_path": "output/clockwork_hidden_gauge_quotient_origin/summary.json",
            "interpretation": "The radial equations QQdag=1, TdagT=1, and QdagQ=TTdag are the D-term moment maps of a vectorlike U(8)_L x U(8)_H x U(8)_R endpoint quotient.  The moment-map rank and fixed-block increment reproduce the radial Hessian, and U(1) charge traces vanish, but the existence/stabilization of the hidden Kahler/FI sector remains a conditional EFT assumption.",
        },
        {
            "claim": "Source-consistent Majorana matrix is Veronese plus contact, not Veronese-only",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"Vres={source_majorana_rank_decomp['veronese_only_relative_residual']:.3e}, "
                f"Cfrac={source_majorana_rank_decomp['contact_fraction']:.3e}, "
                f"lift={source_majorana_rank_decomp['veronese_plus_contact_relative_residual']:.1e}"
            ),
            "evidence_path": "output/source_majorana_texture_rank/report.md",
            "interpretation": (
                "The source-consistent inverse-seesaw M_R decomposes exactly into the "
                "five-complex-dimensional CP1/O(2) Veronese subspace plus the orthogonal "
                "trace/contact direction.  The contact fraction is O(0.13), so the "
                "Majorana sector is not close to a predictive Veronese-only texture; "
                "the trace/contact lift is a real conditional input."
            ),
        },
        {
            "claim": "Majorana contact coefficient must be dynamically source-locked",
            "status": "PASS_CONDITIONAL",
            "key_number": (
                f"s0_angle={majorana_contact_reps['veronese_only_s0']['angle_max_absdiff']:.3e}, "
                f"s_half={majorana_contact_s_width:.3e}, "
                f"phi_half={majorana_contact_phase_width:.3e}"
            ),
            "evidence_path": "output/majorana_contact_sensitivity/report.md",
            "interpretation": (
                "Scanning M_R(s)=M_*(M_V+s zeta K) and the contact phase ring shows "
                "that removing the contact spoils PMNS and splittings, while even a "
                "loose PMNS/dm window around the target requires the complex contact "
                "coefficient to be locked at the 1e-5 level.  A predictive Majorana "
                "sector therefore needs a source action for zeta rather than a free "
                "inverse-seesaw reconstruction."
            ),
        },
        {
            "claim": "Affine Majorana source-locking sector predicts zeta from first principles",
            "status": "TUNED_FALLBACK",
            "key_number": (
                f"Fnorm={majorana_source_locking['f_norm_at_vacuum']:.1e}, "
                f"Hrank={majorana_source_locking['hessian_rank']}/6, "
                f"minsv={min(majorana_source_locking['hessian_singular_values']):.3e}, "
                f"pred={majorana_source_locking_verdict['zeta_value_predicted']}"
            ),
            "evidence_path": "output/majorana_source_locking_sector/report.md",
            "interpretation": (
                "A minimal hidden-singlet superpotential W=A(Z-PQ)+B(P-p0)+C(Q-q0) "
                "locks the target contact coefficient with F-flatness, no flat Hessian "
                "mode, and zero visible threshold vector.  But p0 and q0 encode the "
                "target zeta, so this is a controlled tuned source-locking fallback, "
                "not a microscopic first-principles prediction of the PMNS Majorana "
                "texture."
            ),
        },
        {
            "claim": "Minimal hidden quotient predicts the Majorana contact coefficient",
            "status": "NO_GO",
            "key_number": (
                f"D={majorana_hidden_candidates['D_only_hidden_U1']['hessian_rank']}/4, "
                f"DP={majorana_hidden_candidates['D_plus_product_constraint']['hessian_rank']}/8, "
                f"pred={majorana_hidden_quotient_verdict['pure_quotient_predicts_zeta']}"
            ),
            "evidence_path": "output/majorana_hidden_quotient_origin/report.md",
            "interpretation": (
                "The D-only hidden quotient and the D+product-constraint extension leave "
                "flat directions and do not determine the magnitude or phase of zeta.  "
                "The affine source lock removes the flats only by inserting p0,q0.  "
                "Thus this minimal quotient route is ruled out as a first-principles "
                "origin of the PMNS-sensitive Majorana contact coefficient."
            ),
        },
        {
            "claim": "Small monomial/clockwork hidden sector predicts the Majorana contact coefficient",
            "status": "NO_GO",
            "key_number": (
                f"bestN={majorana_monomial_best['denominator']}, "
                f"res={majorana_monomial_best['phase_residual_rad']:.3e}, "
                f"firstN={majorana_monomial_first['denominator']}, "
                f"scale={majorana_monomial_verdict['scale_relative_tuning_if_direct']:.3e}"
            ),
            "evidence_path": "output/majorana_monomial_clockwork_origin/report.md",
            "interpretation": (
                "Enumerating small real-coefficient root-of-unity phases through "
                "denominator six misses the PMNS-sensitive zeta phase by thousands "
                "of loose tolerance widths; the first denominator within the loose "
                "phase window is 178.  A complex coefficient reintroduces the phase "
                "spurion, and the magnitude still requires a continuous hidden scale "
                "at the 1e-5 level.  Thus small monomial/clockwork is ruled out as a "
                "first-principles origin of zeta."
            ),
        },
        {
            "claim": "External cited flavor target refresh and PMNS completion are completed",
            "status": "OPEN",
            "key_number": (
                f"local={no_web_flavor_verdict['local_flavor_gates_pass']}, "
                f"ckm_diff={no_web_flavor_scores['ckm_score_absdiff']:.1e}, "
                f"mass_diff={no_web_flavor_scores['mass_score_absdiff']:.1e}, "
                f"pmns_bench={no_web_pmns_verdict['local_pmns_benchmark_replay_passes']}, "
                f"source_pmns={source_pmns_verdict['source_consistent_pmns_compatible']}, "
                f"pmns_pred={source_pmns_verdict['pmns_predictive_fit_done']}, "
                f"Cfrac={source_majorana_rank_decomp['contact_fraction']:.2e}"
            ),
            "evidence_path": "output/source_consistent_pmns_replay/report.md",
            "interpretation": (
                "The source-consistent CKM/mass/seesaw flavor card is locally reproducible "
                "from embedded Yukawa matrices and the no-web target rows.  A separate "
                f"exact CP1/O(2) PMNS benchmark replay passes locally with max angle residual "
                f"{no_web_pmns_residuals['max_pmns_angle_absdiff']:.1e}, and the source-consistent "
                f"publication card is PMNS-compatible by inverse-seesaw reconstruction with "
                f"max angle residual {source_pmns_residuals['max_pmns_angle_absdiff']:.1e}.  The remaining "
                "publication-level task is narrower: replace local targets by cited external "
                "CKM/PMNS/mass inputs and derive or fit the Majorana sector rather than "
                "reconstructing M_R from the chosen PMNS target."
            ),
        },
        {
            "claim": "External cited d=5 input refresh is completed",
            "status": "OPEN",
            "key_number": (
                f"local_rows={clockwork_publication_card_d5['total_channel_rows']}, "
                f"worst_margin={clockwork_publication_card_d5['global_worst']['worst_margin_1e35']:.3g}, "
                f"ledger={no_web_input_ledger_verdict['d5_rows_reproduce_formula']}"
            ),
            "evidence_path": "output/no_web_input_convention_ledger/report.md",
            "interpretation": "The exact clockwork-rescued dressed replay and local channel table package now pass, and the no-web ledger reproduces all 21168 displayed rows.  The open item is now strictly the external-input layer: replace local hadronic/chiral constants, present bounds, and flavor targets by a final cited publication input manifest, then rerun the same ledger replay.",
        },
        {
            "claim": "A5-A6 have a microscopic first-principles origin",
            "status": "OPEN",
            "key_number": "composite/auxiliary source dynamics not derived",
            "evidence_path": "roadmap.md",
            "interpretation": "This is the main remaining first-principles gap.",
        },
    ]

    missing = [row["evidence_path"] for row in rows if not exists(row["evidence_path"])]
    if missing:
        raise FileNotFoundError("Missing evidence files: " + ", ".join(missing))

    with (OUT / "claim_ledger.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["claim", "status", "key_number", "evidence_path", "interpretation"],
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "note": "No web lookup used. Conditional theorem evidence ledger.",
        "counts": {
            status: sum(1 for row in rows if row["status"] == status)
            for status in sorted({row["status"] for row in rows})
        },
        "r200_benchmark": {
            "R": safe_r200["R"],
            "kappa_Xi_combined": safe_r200["kappa_Xi_combined"],
            "safe_points": safe_r200["safe_points"],
            "M_Sigma3_GeV": clean_r200["M_Sigma3_GeV"],
            "M_Sigma8_GeV": clean_r200["M_Sigma8_GeV"],
            "alphaG_inv": clean_r200["alphaG_inv"],
            "tau_dim6_years": clean_r200["tau_dim6_years"],
            "triplet_filter_required": clean_r200["triplet_filter_required"],
            "total_projected_l2": safe_r200["total_projected_l2"],
        },
        "open_blockers": [
            row["claim"] for row in rows if row["status"] == "OPEN"
        ],
        "fallbacks": [
            row["claim"] for row in rows if row["status"] == "TUNED_FALLBACK"
        ],
        "verdict": (
            "The current result is a verified conditional Spin(10) EFT branch, "
            "not an unconditional first-principles GUT derivation.  The branch "
            "is internally consistent only with constrained/composite 54/210 "
            "sources and auxiliary/composite conormal multipliers; elementary "
            "propagating fallbacks require explicit threshold tuning."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    report_lines = [
        "# Conditional theorem evidence ledger",
        "",
        "No web lookup used.",
        "",
        f"Rows: {len(rows)}",
        "",
        "## Status counts",
        "",
    ]
    for status, count in summary["counts"].items():
        report_lines.append(f"- {status}: {count}")
    report_lines += [
        "",
        "## R=200 benchmark",
        "",
        f"- M_Sigma3 = {clean_r200['M_Sigma3_GeV']:.6e} GeV",
        f"- M_Sigma8 = {clean_r200['M_Sigma8_GeV']:.6e} GeV",
        f"- alphaG^-1 = {clean_r200['alphaG_inv']:.6f}",
        f"- tau_d6 = {clean_r200['tau_dim6_years']:.6e} yr",
        f"- safe points = {safe_r200['safe_points']}",
        f"- total projected threshold norm = {safe_r200['total_projected_l2']:.6e}",
        "",
        "## Open blockers",
        "",
    ]
    for item in summary["open_blockers"]:
        report_lines.append(f"- {item}")
    report_lines += ["", "## Tuned fallbacks", ""]
    for item in summary["fallbacks"]:
        report_lines.append(f"- {item}")
    report_lines += [
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report_lines))

    print(json.dumps(summary["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
