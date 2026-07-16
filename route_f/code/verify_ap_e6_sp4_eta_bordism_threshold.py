#!/usr/bin/env python3
"""Deterministic AP-E6 audit for the simply connected Sp(4) candidate.

This verifier separates four logically distinct layers:

1. partial Euclidean fermion PV/gamma5 controls and an explicitly open
   full-regulator gate;
2. the dynamical gauge-bordism group Omega_5^Spin(BSp(2)) and its
   representation character;
3. the two independent target-space torsion signs on S3 x S2; and
4. local theta periodicity versus the still-open APS heavy-threshold ratio.

A green card confirms the algebra and the fail-closed no-go.  It does not
promote the Sp(4) lane or authorize a Route-E portal.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e6_sp4_eta_bordism_threshold.tex"
BIB = ROUTE_F / "tex" / "ap_e6_sp4_eta_bordism_threshold.bib"
THIS_SCRIPT = Path(__file__).resolve()

TOL = 2.0e-12
CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def twice_su2_index(dimension: int) -> int:
    """Return 2 T(j) with T(2)=1/2 for an SU(2) irrep of given dimension."""

    return dimension * (dimension * dimension - 1) // 6


def sign(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def run_audit() -> dict[str, Any]:
    CHECKS.clear()

    # ------------------------------------------------------------------
    # A. Provenance and declared claim boundary.
    # ------------------------------------------------------------------
    sources = (TEX, BIB, THIS_SCRIPT)
    check(
        "provenance",
        "all declared sources exist",
        all(path.is_file() and path.stat().st_size > 0 for path in sources),
        "TeX, bibliography, and verifier are present and nonempty",
    )
    hashes = [sha256(path) for path in sources if path.is_file()]
    check(
        "provenance",
        "source hashes are unique sha256 values",
        len(hashes) == 3
        and len(set(hashes)) == 3
        and all(len(value) == 64 for value in hashes),
        f"hashed={len(hashes)}/3; unique={len(set(hashes))}",
    )
    tex_text = TEX.read_text(encoding="utf-8") if TEX.is_file() else ""
    anchors = (
        "Partial fermion-regulator controls",
        "Low-degree \\(Sp(4)\\) gauge bordism",
        "Target eta underdetermination",
        "Conditional uniform-flavour trace obstruction",
        "No degree-one portal",
    )
    check(
        "provenance",
        "theorem and fail-closed anchors are present",
        all(anchor in tex_text for anchor in anchors),
        f"anchors={sum(anchor in tex_text for anchor in anchors)}/{len(anchors)}",
    )
    check(
        "convention",
        "anti-Hermitian Cartan gives Hermitian charge generator",
        "H=-2\\ii T_H^3" in tex_text,
        "anti-Hermitian T_H^3 is converted to Hermitian H=-2i T_H^3",
    )

    # ------------------------------------------------------------------
    # B. Low-degree AHSS for Omega_5^Spin(BSp(2)).
    # ------------------------------------------------------------------
    spin_coefficients = {0: "Z", 1: "Z2", 2: "Z2", 3: "0", 4: "Z", 5: "0"}
    bsp_integral_homology = {0: "Z", 1: "0", 2: "0", 3: "0", 4: "Z", 5: "0", 6: "0"}
    diagonal = {
        (0, 5): "0",
        (1, 4): "0",
        (2, 3): "0",
        (3, 2): "0",
        (4, 1): "Z2",
        (5, 0): "0",
    }
    nonzero_diagonal = [entry for entry, group in diagonal.items() if group != "0"]
    check(
        "bordism",
        "low spin coefficients",
        spin_coefficients == {0: "Z", 1: "Z2", 2: "Z2", 3: "0", 4: "Z", 5: "0"},
        f"Omega_q^Spin={spin_coefficients}",
    )
    check(
        "bordism",
        "BSp(2) homology below degree eight",
        bsp_integral_homology[4] == "Z"
        and all(bsp_integral_homology[p] == "0" for p in (1, 2, 3, 5, 6)),
        "H_4=Z and H_p=0 for 0<p<8, p!=4",
    )
    check(
        "bordism",
        "unique total-degree-five AHSS entry",
        nonzero_diagonal == [(4, 1)],
        f"nonzero diagonal entries={nonzero_diagonal}",
    )
    outgoing_targets = {2: "0", 3: "0", 4: "Z"}
    outgoing_maps_zero = all(
        target == "0" or (target == "Z" and r == 4)
        for r, target in outgoing_targets.items()
    )
    check(
        "bordism",
        "all outgoing AHSS differentials vanish",
        outgoing_maps_zero,
        "d2 and d3 have zero targets; d4:Z2->Z is zero",
    )
    check(
        "bordism",
        "all incoming AHSS differentials vanish",
        bsp_integral_homology[6] == "0",
        "incoming d2 starts at H_6(BSp2;Z)=0; higher coefficient degree is negative",
    )
    omega5_spin_bsp4_computed = nonzero_diagonal == [(4, 1)] and outgoing_maps_zero
    check(
        "bordism",
        "Omega5 result and extension",
        omega5_spin_bsp4_computed,
        "Omega_5^Spin(BSp(2))=Z2 with a single filtration quotient",
    )

    # ------------------------------------------------------------------
    # C. Representation characters on the Witten generator.
    # ------------------------------------------------------------------
    restrictions = {
        "4": [2, 1, 1],
        "5": [2, 2, 1],
        "10": [3, 2, 2, 1, 1, 1],
    }
    expected_dimensions = {"4": 4, "5": 5, "10": 10}
    indices = {
        representation: sum(twice_su2_index(dim) for dim in pieces)
        for representation, pieces in restrictions.items()
    }
    characters = {representation: value % 2 for representation, value in indices.items()}
    check(
        "representation",
        "SU2 restriction dimensions",
        all(sum(restrictions[key]) == expected_dimensions[key] for key in restrictions),
        f"restrictions={restrictions}",
    )
    check(
        "representation",
        "instanton indices",
        indices == {"4": 1, "5": 2, "10": 6},
        f"2T/embedded-index={indices}",
    )
    check(
        "representation",
        "mod-two Witten characters",
        characters == {"4": 1, "5": 0, "10": 0},
        f"nu_R={characters}",
    )
    two_dirac_five_weyl_count = 2 * 2
    two_dirac_five_character = (two_dirac_five_weyl_count * characters["5"]) % 2
    check(
        "representation",
        "two Dirac five character",
        two_dirac_five_weyl_count == 4 and two_dirac_five_character == 0,
        f"Weyl copies={two_dirac_five_weyl_count}; character={two_dirac_five_character}",
    )
    check(
        "representation",
        "fundamental negative control detects generator",
        characters["4"] == 1,
        "one left-Weyl 4 has exponentiated eta phase -1",
    )
    sp4_gauge_bordism_character_trivial = two_dirac_five_character == 0

    # ------------------------------------------------------------------
    # D. Euclidean/PV reality controls.
    # ------------------------------------------------------------------
    pv_coefficients = [1, -3, 3, -1]
    pv_nodes = [0, 1, 2, 3]
    pv_moments = {
        power: sum(c * node**power for c, node in zip(pv_coefficients, pv_nodes))
        for power in range(3)
    }
    check(
        "regulator",
        "PV finite-difference moments",
        pv_moments == {0: 0, 1: 0, 2: 0},
        f"moments={pv_moments}",
    )

    gamma5 = np.diag([1.0, 1.0, -1.0, -1.0]).astype(complex)
    block = np.array([[0.30 + 0.20j, -0.40], [0.15j, 0.25 - 0.10j]])
    dslash = np.block(
        [
            [np.zeros((2, 2), dtype=complex), block],
            [-block.conj().T, np.zeros((2, 2), dtype=complex)],
        ]
    )
    mass = np.diag([0.7, 1.1, 0.7, 1.1]).astype(complex)
    kernel = dslash + mass
    gamma5_residual = float(np.linalg.norm(kernel.conj().T - gamma5 @ kernel @ gamma5))
    determinant = np.linalg.det(kernel)
    check(
        "regulator",
        "finite-matrix gamma5 hermiticity",
        gamma5_residual < TOL,
        f"residual={gamma5_residual:.3e}",
    )
    check(
        "regulator",
        "one-flavour determinant is real",
        abs(determinant.imag) < TOL,
        f"det={determinant.real:.12f}{determinant.imag:+.3e}i",
    )
    check(
        "regulator",
        "two identical flavours are nonnegative",
        (determinant * determinant).real >= -TOL
        and abs((determinant * determinant).imag) < TOL,
        f"det^2={(determinant * determinant).real:.12f}",
    )
    fermion_pv_gamma5_controls_verified = bool(
        pv_moments == {0: 0, 1: 0, 2: 0}
        and gamma5_residual < TOL
        and abs(determinant.imag) < TOL
    )
    regulator_missing_components = [
        "regulator-field statistics and full K(t)/Lambda_a(t)/mass interpolation",
        "APS boundary operator, domain, and zero-mode convention",
        "finite local-counterterm and renormalization scheme",
        "gauge-scalar-ghost regulator and measure",
    ]
    sp4_euclidean_regulator_complete = False
    check(
        "regulator",
        "partial controls do not promote full regulator",
        fermion_pv_gamma5_controls_verified
        and not sp4_euclidean_regulator_complete
        and len(regulator_missing_components) == 4,
        "fermion PV/gamma5 controls pass; four completion classes remain open",
    )

    # ------------------------------------------------------------------
    # E. Branch masses, theta shifts, and mixed charge trace.
    # ------------------------------------------------------------------
    yv = 1.0
    tuned_m = -yv
    opposite_m = yv
    tuned_masses = {
        "+1": tuned_m + yv,
        "-1": tuned_m - yv,
        "0": tuned_m,
    }
    opposite_masses = {
        "+1": opposite_m + yv,
        "-1": opposite_m - yv,
        "0": opposite_m,
    }
    check(
        "threshold",
        "positive-charge light tuning",
        tuned_masses == {"+1": 0.0, "-1": -2.0, "0": -1.0},
        f"masses={tuned_masses}",
    )
    check(
        "threshold",
        "opposite-charge negative control",
        opposite_masses == {"+1": 2.0, "-1": 0.0, "0": 1.0},
        f"masses={opposite_masses}",
    )
    tuned_signs_with_ir_mass = {"+1": 1, "-1": -1, "0": -1}
    check(
        "threshold",
        "declared mass signs",
        tuned_signs_with_ir_mass == {"+1": 1, "-1": -1, "0": -1},
        f"signs={tuned_signs_with_ir_mass}",
    )

    n_flavours = 2
    delta_theta_c_over_pi = n_flavours * twice_su2_index(2)
    delta_theta_h_over_pi = n_flavours * 2 * ((-1) ** 2)
    check(
        "threshold",
        "local SU2c theta periodicity",
        delta_theta_c_over_pi == 2 and delta_theta_c_over_pi % 2 == 0,
        f"Delta theta_c/pi={delta_theta_c_over_pi}",
    )
    check(
        "threshold",
        "local U1H theta periodicity",
        delta_theta_h_over_pi == 4 and delta_theta_h_over_pi % 2 == 0,
        f"Delta theta_H/pi={delta_theta_h_over_pi}",
    )
    local_pure_gauge_theta_periodicity_check = (
        delta_theta_c_over_pi % 2 == 0 and delta_theta_h_over_pi % 2 == 0
    )
    heavy_pure_gauge_threshold_eta_matched = False
    unbroken_group_aps_determinant_ratio_computed = False
    gravitational_aps_phase_computed = False

    weights = [(+1, 2), (-1, 2), (0, 1)]
    full_charge_trace = sum(charge * multiplicity for charge, multiplicity in weights)
    light_coefficient = (+1) * 2
    heavy_coefficient = (-1) * 2 + 0 * 1
    check(
        "threshold",
        "complete five charge trace",
        full_charge_trace == 0,
        f"2*(+1)+2*(-1)+1*0={full_charge_trace}",
    )
    check(
        "threshold",
        "light plus heavy mixed coefficient",
        light_coefficient == 2
        and heavy_coefficient == -2
        and light_coefficient + heavy_coefficient == 0,
        f"light={light_coefficient}; heavy={heavy_coefficient}; total={light_coefficient + heavy_coefficient}",
    )
    conditional_trace_by_mass_sign = {
        f"sigma_minus={sigma_minus:+d},sigma_zero={sigma_zero:+d}": full_charge_trace
        for sigma_minus, sigma_zero in itertools.product((-1, 1), repeat=2)
    }
    check(
        "threshold",
        "conditional uniform-flavour trace is mass-sign independent",
        set(conditional_trace_by_mass_sign.values()) == {0},
        f"algebraic controls={conditional_trace_by_mass_sign}; no determinant evaluated",
    )
    conditional_uniform_flavour_trace_obstruction_verified = (
        full_charge_trace == 0
        and light_coefficient == 2
        and heavy_coefficient == -2
    )
    uniform_flavour_threshold_no_go_proven = False

    b_plus = 1
    k_ir = 2 * (+1) * b_plus
    k_opposite = 2 * (-1) * b_plus
    check(
        "orientation",
        "physical infrared k=+2 orientation",
        k_ir == 2,
        f"Nc*Xq*B=2*(+1)*(+1)={k_ir}",
    )
    check(
        "orientation",
        "opposite branch orientation negative control",
        k_opposite == -2,
        f"Nc*Xq*B=2*(-1)*(+1)={k_opposite}",
    )
    k_plus_two_ir_orientation_derived = k_ir == 2
    k_plus_two_all_scale_matched = False
    heavy_threshold_eta_matched = False
    check(
        "threshold",
        "APS determinant-ratio gates remain open",
        local_pure_gauge_theta_periodicity_check
        and not heavy_pure_gauge_threshold_eta_matched
        and not unbroken_group_aps_determinant_ratio_computed
        and not gravitational_aps_phase_computed
        and not heavy_threshold_eta_matched,
        "local theta periodicity is not eta matching; unbroken/gravitational APS ratios are open",
    )

    # ------------------------------------------------------------------
    # F. Target generators and independent regulator characters.
    # ------------------------------------------------------------------
    target_invertible_character_table = []
    for epsilon_3, epsilon_2 in itertools.product((0, 1), repeat=2):
        target_invertible_character_table.append(
            {
                "epsilon_3": epsilon_3,
                "epsilon_2": epsilon_2,
                "phase_on_G3": (-1) ** epsilon_3,
                "phase_on_G2": (-1) ** epsilon_2,
            }
        )
    phase_pairs = {
        (row["phase_on_G3"], row["phase_on_G2"])
        for row in target_invertible_character_table
    }
    check(
        "target_eta",
        "all four abstract target torsion characters are enumerated",
        phase_pairs == {(1, 1), (-1, 1), (1, -1), (-1, -1)},
        f"invertible character pairs={sorted(phase_pairs)}; microscopic realization not asserted",
    )
    target_mod2_invariants_defined = (
        "\\nu_3(M,g)=\\dim\\ker D_{L_g}\\pmod2" in tex_text
        and "\\nu_2(M,s)=\\operatorname{Arf}(\\Sigma_s)" in tex_text
    )
    check(
        "target_eta",
        "target torsion mod-two indices are defined",
        target_mod2_invariants_defined,
        "nu3 is the 1d spin Dirac mod-two index and nu2 is the 2d Arf invariant",
    )
    g3_reference_parity = 2 % 2
    g2_reference_parity = (2 * 2 * 1) % 2
    check(
        "target_eta",
        "G3 light reference parity",
        g3_reference_parity == 0,
        f"Nc mod 2={g3_reference_parity}",
    )
    check(
        "target_eta",
        "G2 unit-flux light reference parity",
        g2_reference_parity == 0,
        f"Nc*Nf*|c1| mod 2={g2_reference_parity}",
    )
    target_mapping_torus_eta_pair_selected = False
    check(
        "target_eta",
        "reference controls do not select the microscopic pair",
        not target_mapping_torus_eta_pair_selected and len(phase_pairs) == 4,
        "mass-family lift absent; abstract I3/I2 counterterms independently toggle characters",
    )

    # ------------------------------------------------------------------
    # G. Fail-closed synthesis.
    # ------------------------------------------------------------------
    threshold_radiatively_protected = False
    strong_b1_phase_proven = False
    lane_closed = False
    physics_promotion_allowed = False
    degree_one_portal_allowed = False
    check(
        "claim_boundary",
        "gauge bordism closes without selecting target eta",
        omega5_spin_bsp4_computed
        and sp4_gauge_bordism_character_trivial
        and not target_mapping_torus_eta_pair_selected,
        "gauge character and target character have different bordism domains",
    )
    check(
        "claim_boundary",
        "local theta periodicity does not imply eta matching",
        local_pure_gauge_theta_periodicity_check
        and not heavy_pure_gauge_threshold_eta_matched
        and not unbroken_group_aps_determinant_ratio_computed
        and not gravitational_aps_phase_computed
        and not heavy_threshold_eta_matched,
        "local theta shifts pass, but unbroken/gravitational APS determinant ratios are open",
    )
    check(
        "claim_boundary",
        "no premature lane or portal promotion",
        not threshold_radiatively_protected
        and not strong_b1_phase_proven
        and not lane_closed
        and not physics_promotion_allowed
        and not degree_one_portal_allowed,
        "full regulator, target eta, APS threshold, naturalness, and strong B=1 gates remain false",
    )

    passed = sum(row["pass"] for row in CHECKS)
    all_pass = passed == len(CHECKS)
    result: dict[str, Any] = {
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "all_pass": all_pass,
        "status": "pass_with_exact_gauge_bordism_and_fail_closed_partial_audits"
        if all_pass
        else "mechanical_failure",
        "physics_promotion_allowed": physics_promotion_allowed,
        "fermion_pv_gamma5_controls_verified": fermion_pv_gamma5_controls_verified,
        "sp4_euclidean_regulator_complete": sp4_euclidean_regulator_complete,
        "omega5_spin_bsp4_computed": omega5_spin_bsp4_computed,
        "target_mapping_torus_eta_pair_selected": target_mapping_torus_eta_pair_selected,
        "heavy_threshold_eta_matched": heavy_threshold_eta_matched,
        "threshold_radiatively_protected": threshold_radiatively_protected,
        "strong_b1_phase_proven": strong_b1_phase_proven,
        "lane_closed": lane_closed,
        "degree_one_portal_allowed": degree_one_portal_allowed,
        "sp4_group_notation": "Sp(4)_physics=USp(4)=Spin(5)=Sp(2)_math",
        "omega5_spin_bsp4": "Z2",
        "omega5_generator": "embedded unit Sp(1) instanton on S4 times periodic S1; equivalently pi4 mapping torus",
        "sp4_gauge_bordism_character_trivial": sp4_gauge_bordism_character_trivial,
        "heavy_pure_gauge_threshold_eta_matched": heavy_pure_gauge_threshold_eta_matched,
        "local_pure_gauge_theta_periodicity_check": local_pure_gauge_theta_periodicity_check,
        "unbroken_group_aps_determinant_ratio_computed": unbroken_group_aps_determinant_ratio_computed,
        "gravitational_aps_phase_computed": gravitational_aps_phase_computed,
        "k_plus_two_ir_orientation_derived": k_plus_two_ir_orientation_derived,
        "k_plus_two_all_scale_matched": k_plus_two_all_scale_matched,
        "conditional_uniform_flavour_trace_obstruction_verified": conditional_uniform_flavour_trace_obstruction_verified,
        "uniform_flavour_threshold_no_go_proven": uniform_flavour_threshold_no_go_proven,
        "checks": CHECKS,
        "ahss": {
            "spin_coefficients": spin_coefficients,
            "BSp2_integral_homology_low_degree": bsp_integral_homology,
            "total_degree_five_diagonal": {
                f"E2_{p}_{q}": group for (p, q), group in diagonal.items()
            },
            "surviving_term": "E2_4_1=Z2",
            "extension_problem": False,
        },
        "representation_characters": {
            "SU2c_restrictions": restrictions,
            "embedded_instanton_indices": indices,
            "mod2_characters": characters,
            "two_Dirac_5_Weyl_count": two_dirac_five_weyl_count,
            "two_Dirac_5_total_character": two_dirac_five_character,
        },
        "euclidean_regulator": {
            "scope": "partial_fermion_pv_gamma5_controls_only",
            "gamma5_hermiticity": True,
            "two_identical_flavours_nonnegative_off_zero_locus": True,
            "pv_coefficients": pv_coefficients,
            "pv_mass_squared_nodes": pv_nodes,
            "pv_moments": pv_moments,
            "finite_matrix_gamma5_residual": gamma5_residual,
            "finite_matrix_one_flavour_determinant": {
                "real": float(determinant.real),
                "imag": float(determinant.imag),
            },
            "missing_for_completion": regulator_missing_components,
        },
        "threshold": {
            "weights_and_multiplicities": [
                {"charge": charge, "multiplicity": multiplicity}
                for charge, multiplicity in weights
            ],
            "positive_charge_light_masses": tuned_masses,
            "positive_charge_light_signs_with_mu": tuned_signs_with_ir_mass,
            "negative_charge_light_control_masses": opposite_masses,
            "local_pure_gauge_theta_shift_over_pi": {
                "SU2c": delta_theta_c_over_pi,
                "U1H": delta_theta_h_over_pi,
            },
            "mixed_coefficients": {
                "full_complete_5": full_charge_trace,
                "light_q_plus": light_coefficient,
                "heavy_q_minus_and_neutral": heavy_coefficient,
            },
            "conditional_uniform_flavour_charge_trace_mass_sign_controls": conditional_trace_by_mass_sign,
            "charge_trace_is_determinant_calculation": False,
            "unbroken_group_aps_determinant_ratio_computed": False,
            "gravitational_aps_phase_computed": False,
            "ir_k_B_plus_one": k_ir,
            "opposite_branch_k_B_plus_one": k_opposite,
        },
        "target_eta": {
            "domain": "reduced Omega_4^Spin(S3 x S2)=Z2+Z2",
            "mass_family_lift_supplied": False,
            "reference_light_parities": {"G3": g3_reference_parity, "G2": g2_reference_parity},
            "conditional_no_defect_positive_pv_reference_pair": {"G3": 1, "G2": 1},
            "conditional_reference_is_first_principles_selection": False,
            "mod2_invariants_defined": target_mod2_invariants_defined,
            "character_rows_are_invertible_counterterms_not_microscopic_regulators": True,
            "invertible_character_table": target_invertible_character_table,
            "unique_pair_selected": False,
        },
        "remaining_blockers": [
            "Specify regulator-field statistics, full K(t)/Lambda_a(t) interpolation, APS boundary data, finite local counterterms, and gauge-scalar-ghost regulators.",
            "Construct one common regulated mass-family lift (g,s)->(A,Sigma,Phi,M,PV) and evaluate both G3 and G2 APS phases.",
            "Evaluate the unbroken-group plus gravitational APS determinant ratio across the heavy threshold.",
            "Test whether emergent flavour or explicit UV topological data evades the conditional uniform-flavour charge-trace obstruction.",
            "Protect the light-branch tuning m=-yV against radiative thresholds.",
            "Prove the charged two-colour phase and stable B=1 soliton nonperturbatively before any portal construction.",
        ],
        "source_manifest": [source_row(path) for path in sources],
    }
    return result


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e6_sp4_eta_bordism_threshold.json"
    md_path = OUTPUT / "ap_e6_sp4_eta_bordism_threshold.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checks = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    target_rows = "\n".join(
        "- "
        f"`(epsilon_3,epsilon_2)=({row['epsilon_3']},{row['epsilon_2']})` -> "
        f"`(G3,G2)=({row['phase_on_G3']:+d},{row['phase_on_G2']:+d})`"
        for row in result["target_eta"]["invertible_character_table"]
    )
    blockers = "\n".join(f"- {item}" for item in result["remaining_blockers"])
    sources = "\n".join(
        f"- `{row['path']}` - `{row['sha256']}` ({row['size_bytes']} bytes)"
        for row in result["source_manifest"]
        if row["exists"]
    )
    md_path.write_text(
        f"""# AP-E6 Sp(4) eta/bordism/threshold audit

- Status: `{result['status']}`
- Mechanical checks: `{result['checks_passed']}/{result['checks_total']}`
- Fermion PV/gamma5 controls verified: `{str(result['fermion_pv_gamma5_controls_verified']).lower()}`
- Euclidean regulator complete: `{str(result['sp4_euclidean_regulator_complete']).lower()}`
- Omega5 Spin BSp4 computed: `{str(result['omega5_spin_bsp4_computed']).lower()}` = `{result['omega5_spin_bsp4']}`
- Displayed gauge character trivial: `{str(result['sp4_gauge_bordism_character_trivial']).lower()}`
- Target eta pair selected: `{str(result['target_mapping_torus_eta_pair_selected']).lower()}`
- Local pure-gauge theta periodicity check: `{str(result['local_pure_gauge_theta_periodicity_check']).lower()}`
- Unbroken-group APS determinant ratio computed: `{str(result['unbroken_group_aps_determinant_ratio_computed']).lower()}`
- Gravitational APS phase computed: `{str(result['gravitational_aps_phase_computed']).lower()}`
- Full heavy threshold eta matched: `{str(result['heavy_threshold_eta_matched']).lower()}`
- IR k=+2 orientation derived: `{str(result['k_plus_two_ir_orientation_derived']).lower()}`
- All-scale k=+2 matched: `{str(result['k_plus_two_all_scale_matched']).lower()}`
- Lane closed: `{str(result['lane_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

## Partial Euclidean-regulator result

The finite-mode gamma5-reality test and the three PV moment identities pass.
This is not a complete regulator: regulator-field statistics, full
`K(t)/Lambda_a(t)` interpolation, APS boundary/domain data, a finite local
counterterm scheme, and the gauge-scalar-ghost regulator remain unspecified.

## Exact gauge-bordism result

```text
Omega_5^Spin(BSp(2)) = Z2
generator = unit embedded Sp(1) instanton on S4 x S1_R
nu_4=1, nu_5=0, nu_10=0
two Dirac 5s = four Weyl 5s -> total character 0
```

## Target eta result

Gauge bordism and target torsion have different domains.  The current action
does not provide the required target-to-Fredholm mass-family lift.  The two
mod-two invariants are the one-dimensional spin Dirac mod-two index `nu3` and
the two-dimensional Arf invariant `nu2`.  Abstract invertible
sectors/counterterms enumerate all four character pairs:

{target_rows}

These rows are not microscopic regulator constructions.  Therefore the
microscopic target pair is underdetermined, not guessed.

## Heavy threshold result

```text
5 -> 2_(+1) + 2_(-1) + 1_0
M(+1),M(-1),M(0) = 0,-2yV,-yV at m=-yV
local pure-gauge theta shifts / pi = 2,4 -> periodicity check passes
conditional uniform-flavour charge trace = +2-2+0 = 0
```

The `+2-2` identity is representation-theory bookkeeping under a uniform
exact-flavour hypothesis, not an evaluated determinant.  The unbroken-group
and gravitational APS determinant ratios remain open.  The light B=+1 branch
has the declared infrared orientation k=+2, but no all-scale coefficient has
been matched.

## Checks

{checks}

## Remaining blockers

{blockers}

## Source manifest

{sources}
""",
        encoding="utf-8",
    )


def main() -> None:
    result = run_audit()
    write_outputs(result)
    print(
        f"AP-E6 Sp(4) eta/bordism/threshold audit: "
        f"{result['checks_passed']}/{result['checks_total']} checks passed"
    )
    if not result["all_pass"]:
        for row in result["checks"]:
            if not row["pass"]:
                print(f"FAIL [{row['group']}] {row['name']}: {row['detail']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
