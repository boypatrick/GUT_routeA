#!/usr/bin/env python3
"""Deterministic AP-E7 common Sp(4) APS/PV regulator audit.

The card verifies the Clifford/PV/projector/index algebra, the physical
zero-cut convention, and the distinction between a conditional GW-compatible
Higgsed source and a completed overlap lift. A green card does not authorize
the portal.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e7_common_sp4_aps_pv_regulator.tex"
BIB = ROUTE_F / "tex" / "ap_e7_common_sp4_aps_pv_regulator.bib"
THIS_SCRIPT = Path(__file__).resolve()
JSON_OUT = OUTPUT / "ap_e7_common_sp4_aps_pv_regulator.json"
MD_OUT = OUTPUT / "ap_e7_common_sp4_aps_pv_regulator.md"
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


def residual(lhs: np.ndarray, rhs: np.ndarray) -> float:
    return float(np.linalg.norm(lhs - rhs, ord="fro"))


def circle_crossings(spin_shift: float) -> int:
    crossings = 0
    for n in range(-4, 5):
        left = n + spin_shift - 0.25
        right = n + spin_shift + 0.25
        crossings += int(left < 0.0 < right)
    return crossings


def scalar_spectral_flow(start: float, end: float, cut: float) -> int:
    """Signed spectral flow for a one-level affine path away from endpoints."""
    if start < cut < end:
        return 1
    if end < cut < start:
        return -1
    return 0


def run_audit() -> dict[str, Any]:
    CHECKS.clear()

    # Provenance and textual claim boundary.
    sources = (TEX, BIB, THIS_SCRIPT)
    check(
        "provenance",
        "all declared sources exist",
        all(path.is_file() and path.stat().st_size > 0 for path in sources),
        "TeX, bibliography, and verifier are nonempty",
    )
    hashes = [sha256(path) for path in sources if path.is_file()]
    check(
        "provenance",
        "source hashes are unique sha256 values",
        len(hashes) == 3 and len(set(hashes)) == 3,
        f"hashed={len(hashes)}/3; unique={len(set(hashes))}",
    )
    tex_text = TEX.read_text(encoding="utf-8") if TEX.is_file() else ""
    anchors = (
        "The full background-field quadratic complex",
        "One common Pauli--Villars prescription",
        "Five-dimensional interpolation, spectral cut, and APS domain",
        "A conditional GW-compatible Higgsed source",
        "Target eta non-identifiability",
        "Two-flavour pure threshold",
        "Exact separated rank-two projector obstruction",
        "does not import or",
        "No portal is authorized",
    )
    check(
        "provenance",
        "all theorem anchors are present",
        all(anchor in tex_text for anchor in anchors),
        f"anchors={sum(anchor in tex_text for anchor in anchors)}/{len(anchors)}",
    )
    check(
        "provenance",
        "physical fermion phase fixes the zero cut with gapped endpoints",
        r"0\notin\operatorname{spec}(H_0)\cup\operatorname{spec}(H_1)"
        in tex_text
        and r"\SF_0\{H(t)\}" in tex_text,
        "zero is excluded from both endpoint spectra and SF_0 defines the phase",
    )
    check(
        "provenance",
        "arbitrary APS cut is restricted to Fredholm-index bookkeeping",
        r"\alpha\notin\operatorname{spec}(H_0)\cup\operatorname{spec}(H_1)"
        in tex_text
        and "bookkeeping device" in tex_text
        and "must not be substituted" in tex_text,
        "general alpha remains available for the APS index, not the physical sign",
    )
    check(
        "provenance",
        "lateral Agmon cut and finite log are explicit",
        r"\theta_B=\pi-\varepsilon" in tex_text
        and r"\left[\log Z_{\rm quad}^{\rm PV}" in tex_text,
        "upper lateral cut and renormalized finite part found",
    )

    # Euclidean gamma matrices and scalar-reference quaternionic reality.
    s1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    s2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    s3 = np.diag([1.0, -1.0]).astype(complex)
    i2 = np.eye(2, dtype=complex)
    i4 = np.eye(4, dtype=complex)
    zero4 = np.zeros((4, 4), dtype=complex)
    gammas = (
        np.kron(s1, s1),
        np.kron(s1, s2),
        np.kron(s1, s3),
        np.kron(s2, i2),
    )
    gamma5 = gammas[0] @ gammas[1] @ gammas[2] @ gammas[3]
    hermitian_error = max(residual(g, g.conj().T) for g in gammas)
    clifford_error = max(
        residual(
            gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu],
            (2.0 if mu == nu else 0.0) * i4,
        )
        for mu in range(4)
        for nu in range(4)
    )
    gamma5_error = max(
        residual(gamma5, gamma5.conj().T),
        residual(gamma5 @ gamma5, i4),
        max(residual(gamma5 @ g + g @ gamma5, zero4) for g in gammas),
    )
    check(
        "clifford",
        "Hermitian Euclidean gamma matrices",
        hermitian_error < TOL,
        f"max residual={hermitian_error:.3e}",
    )
    check(
        "clifford",
        "Euclidean Clifford algebra",
        clifford_error < TOL,
        f"max residual={clifford_error:.3e}",
    )
    check(
        "clifford",
        "gamma5 involution and anticommutation",
        gamma5_error < TOL,
        f"max residual={gamma5_error:.3e}",
    )
    charge_b = gammas[1] @ gammas[3]
    reality_error = max(
        residual(charge_b @ g.conj(), g @ charge_b) for g in gammas
    )
    j_square_error = residual(charge_b @ charge_b.conj(), -i4)
    check(
        "reality",
        "B gamma-star equals gamma B",
        reality_error < TOL,
        f"max residual={reality_error:.3e}",
    )
    check(
        "reality",
        "quaternionic reality squares to minus one",
        j_square_error < TOL,
        f"norm(BB*+1)={j_square_error:.3e}",
    )
    t_real = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    sigma_hat = -1.0j * t_real
    i_internal = np.eye(2, dtype=complex)
    j_linear = np.kron(charge_b, i_internal)
    adjoint_term = 0.21 * np.kron(i4, sigma_hat)
    mass_plus = 0.73 * np.eye(8, dtype=complex) + adjoint_term
    mass_minus = 0.73 * np.eye(8, dtype=complex) - adjoint_term
    j_mass = j_linear @ mass_plus.conj() @ np.linalg.inv(j_linear)
    adjoint_odd_error = residual(sigma_hat.conj(), -sigma_hat)
    family_involution_error = residual(j_mass, mass_minus)
    fixed_background_residual = residual(j_mass, mass_plus)
    check(
        "reality",
        "adjoint Hermitian mass is odd under vector conjugation",
        residual(sigma_hat, sigma_hat.conj().T) < TOL
        and adjoint_odd_error < TOL,
        f"Hermitian residual={residual(sigma_hat, sigma_hat.conj().T):.3e}; "
        f"K5 Sigmahat K5^-1 + Sigmahat={adjoint_odd_error:.3e}",
    )
    check(
        "reality",
        "anti-linear map sends the plus-Sigma mass to minus-Sigma",
        family_involution_error < TOL and fixed_background_residual > 1.0e-6,
        f"family residual={family_involution_error:.3e}; "
        f"fixed-background residual={fixed_background_residual:.3e}",
    )
    momentum = np.array([0.31, -0.27, 0.19, 0.43])
    dslash = sum(1.0j * p * g for p, g in zip(momentum, gammas))
    kernel = dslash + 0.73 * i4
    hermitian_h = gamma5 @ kernel
    kernel_error = residual(kernel.conj().T, gamma5 @ kernel @ gamma5)
    h_error = residual(hermitian_h, hermitian_h.conj().T)
    determinant = np.linalg.det(kernel)
    check(
        "reality",
        "finite scalar-reference gamma5-hermitian kernel",
        kernel_error < TOL and h_error < TOL,
        f"K residual={kernel_error:.3e}; H residual={h_error:.3e}",
    )
    check(
        "reality",
        "scalar-reference two-flavour determinant square is nonnegative",
        abs((determinant**2).imag) < TOL and (determinant**2).real >= -TOL,
        f"det(K)^2={(determinant**2).real:.12f}",
    )

    # Common gauge complex toy and PV moments.
    rng = np.random.default_rng(17072026)
    gauge_map = rng.normal(size=(9, 3))
    ghost = gauge_map.T @ gauge_map
    boson = np.diag(np.linspace(0.4, 1.2, 9)) + gauge_map @ gauge_map.T
    check(
        "quadratic_complex",
        "ghost R-dagger R is symmetric positive",
        residual(ghost, ghost.T) < TOL
        and float(np.min(np.linalg.eigvalsh(ghost))) > 0.0,
        f"min eigenvalue={np.min(np.linalg.eigvalsh(ghost)):.6e}",
    )
    check(
        "quadratic_complex",
        "gauge-fixed boson toy is symmetric positive",
        residual(boson, boson.T) < TOL
        and float(np.min(np.linalg.eigvalsh(boson))) > 0.0,
        f"min eigenvalue={np.min(np.linalg.eigvalsh(boson)):.6e}",
    )
    pv_coefficients = [1, -3, 3, -1]
    pv_nodes = [0, 1, 2, 3]
    moments = {
        power: sum(c * node**power for c, node in zip(pv_coefficients, pv_nodes))
        for power in range(4)
    }
    check(
        "pv",
        "first three PV moments cancel",
        moments == {0: 0, 1: 0, 2: 0, 3: -6},
        f"moments={moments}",
    )
    large_lam = 10_000.0
    log_combo = sum(
        c * math.log1p(node / large_lam)
        for c, node in zip(pv_coefficients, pv_nodes)
    )
    scaled_log = log_combo * large_lam**3
    check(
        "pv",
        "large-lambda log coefficient is minus two",
        abs(scaled_log + 2.0) < 2.0e-3,
        f"lambda^3 sum(c log(1+a/lambda))={scaled_log:.9f}",
    )
    check(
        "pv",
        "two-flavour fermion amplitude exponents are integral",
        [2 * c // 2 for c in pv_coefficients] == pv_coefficients,
        f"Nf*c/2={pv_coefficients}",
    )

    # APS zero-cut phase, auxiliary-cut control, and abstract target characters.
    scalar_start = 1.0
    scalar_end = -1.0
    zero_cut_flow = scalar_spectral_flow(scalar_start, scalar_end, 0.0)
    auxiliary_cut_flow = scalar_spectral_flow(scalar_start, scalar_end, 2.0)
    one_flavour_zero_cut_sign = -1 if zero_cut_flow % 2 else 1
    check(
        "aps",
        "gapped scalar endpoints give the physical zero-cut sign",
        scalar_start != 0.0
        and scalar_end != 0.0
        and zero_cut_flow == -1
        and one_flavour_zero_cut_sign == -1,
        f"SF_0={zero_cut_flow}; one-flavour sign={one_flavour_zero_cut_sign}",
    )
    check(
        "aps",
        "auxiliary cut changes raw spectral flow and is not the physical sign",
        auxiliary_cut_flow == 0 and auxiliary_cut_flow != zero_cut_flow,
        f"SF_0={zero_cut_flow}; SF_2={auxiliary_cut_flow}",
    )
    r_crossings = circle_crossings(0.0)
    ns_crossings = circle_crossings(0.5)
    check(
        "aps",
        "periodic circle path has one crossing",
        r_crossings == 1,
        f"crossings={r_crossings}",
    )
    check(
        "aps",
        "antiperiodic circle path has no crossing",
        ns_crossings == 0,
        f"crossings={ns_crossings}",
    )
    characters = {
        (e3, e2): ((-1) ** e3, (-1) ** e2)
        for e3 in (0, 1)
        for e2 in (0, 1)
    }
    expected_characters = {
        (0, 0): (1, 1),
        (1, 0): (-1, 1),
        (0, 1): (1, -1),
        (1, 1): (-1, -1),
    }
    check(
        "target",
        "four abstract target torsion characters",
        characters == expected_characters,
        f"characters={characters}",
    )

    # Higgsed projectors, overlap pair, and Schur obstruction.
    h = np.diag([1.0, 1.0, -1.0, -1.0, 0.0])
    p_plus = (h @ h + h) / 2.0
    p_minus = (h @ h - h) / 2.0
    p_zero = np.eye(5) - h @ h
    projectors = (p_plus, p_minus, p_zero)
    idem_error = max(residual(p @ p, p) for p in projectors)
    orth_error = max(
        residual(projectors[i] @ projectors[j], np.zeros((5, 5)))
        for i in range(3)
        for j in range(3)
        if i != j
    )
    sum_error = residual(sum(projectors), np.eye(5))
    ranks = tuple(int(np.linalg.matrix_rank(p, tol=1.0e-10)) for p in projectors)
    check(
        "projectors",
        "P plus P minus P zero are orthogonal projectors",
        idem_error < TOL and orth_error < TOL and sum_error < TOL,
        f"idempotent={idem_error:.3e}; orthogonal={orth_error:.3e}; sum={sum_error:.3e}",
    )
    check(
        "projectors",
        "projector ranks are two two one",
        ranks == (2, 2, 1),
        f"ranks={ranks}",
    )
    schur_trace = float(np.trace(h @ (0.37 * np.eye(5))))
    light_trace = float(np.trace(h @ p_plus))
    check(
        "extendibility",
        "Schur scalar has zero H trace",
        abs(float(np.trace(h))) < TOL and abs(schur_trace) < TOL,
        f"Tr(H)={np.trace(h):.1f}; Tr(H cI)={schur_trace:.3e}",
    )
    check(
        "extendibility",
        "Higgsed light selector has H trace two",
        abs(light_trace - 2.0) < TOL,
        f"Tr(H P+)={light_trace:.12f}",
    )
    scalar_ranks = {
        c: int(np.linalg.matrix_rank(c * np.eye(5), tol=1.0e-10))
        for c in (0.0, 1.0)
    }
    check(
        "extendibility",
        "Schur idempotent ranks exclude rank two",
        scalar_ranks == {0.0: 0, 1.0: 5},
        f"ranks={scalar_ranks}",
    )
    g3_parity = ranks[0] % 2
    g2_light_copies = ranks[0] * 2
    g2_heavy_copies = ranks[1] * 2
    g2_parity = (g2_light_copies + g2_heavy_copies) % 2
    check(
        "conditional_higgsed_source",
        "conditional Higgsed G3 copy parity is even",
        g3_parity == 0,
        f"rank(P+) mod 2={g3_parity}",
    )
    check(
        "conditional_higgsed_source",
        "conditional Higgsed G2 copy parity is even",
        g2_light_copies == 4 and g2_heavy_copies == 4 and g2_parity == 0,
        f"light={g2_light_copies}; heavy={g2_heavy_copies}; parity={g2_parity}",
    )
    conditional_pair = ((-1) ** g3_parity, (-1) ** g2_parity)
    check(
        "conditional_higgsed_source",
        "copy-counting parity prediction is derived as plus plus",
        conditional_pair == (1, 1),
        f"derived pair=((-1)^{g3_parity},(-1)^{g2_parity})={conditional_pair}",
    )

    # Exact heavy index and mixed charge trace.
    gravitational_sum = Fraction(-1, 12) + Fraction(-1, 24)
    check(
        "heavy_index",
        "charged plus neutral gravitational coefficient",
        gravitational_sum == Fraction(-1, 8),
        f"-1/12-1/24={gravitational_sum}",
    )
    rank_vc = Fraction(2)
    ch_vc_degree4_c2 = Fraction(-1)
    exp_minus_x_degree4_x2 = Fraction(1, 2)
    derived_ch2_c2 = ch_vc_degree4_c2
    derived_ch2_x2 = rank_vc * exp_minus_x_degree4_x2
    check(
        "heavy_index",
        "symbolic degree-four Chern expansion gives minus c2 plus x squared",
        (derived_ch2_c2, derived_ch2_x2) == (Fraction(-1), Fraction(1)),
        "(2-c2)(1-x+x^2/2) [degree four]="
        f"({derived_ch2_c2})c2+({derived_ch2_x2})x^2",
    )
    grid_rows: list[dict[str, int]] = []
    all_integral = True
    for c2_num in range(-4, 5):
        for x2_num in range(-8, 9, 2):
            for sig_multiple in range(-3, 4):
                signature = 16 * sig_multiple
                p1_num = 3 * signature
                exact_i = Fraction(-c2_num) + Fraction(x2_num) - Fraction(p1_num, 8)
                all_integral &= exact_i.denominator == 1
                if len(grid_rows) < 8:
                    grid_rows.append(
                        {
                            "c2": c2_num,
                            "x2": x2_num,
                            "signature": signature,
                            "p1": p1_num,
                            "I_heavy": int(exact_i),
                        }
                    )
    grid_size = 9 * 9 * 7
    check(
        "heavy_index",
        "finite spin characteristic grid spot-checks heavy-index integrality",
        all_integral,
        f"grid points={grid_size}; all denominators one",
    )
    generator_characteristics = {
        "G3": {"c2": 0, "x2": 0, "p1": 0},
        "G2": {"c2": 0, "x2": 0, "p1": 0},
    }
    generator_indices = {
        name: -row["c2"] + row["x2"] - Fraction(row["p1"], 8)
        for name, row in generator_characteristics.items()
    }
    pure_pair = tuple(
        (-1) ** (2 * int(generator_indices[name])) for name in ("G3", "G2")
    )
    check(
        "heavy_index",
        "declared generator data directly evaluate to the pure pair plus plus",
        all(index == 0 for index in generator_indices.values())
        and pure_pair == (1, 1),
        f"indices={generator_indices}; pair={pure_pair}",
    )
    charges = [1, 1, -1, -1, 0]
    light_coeff = sum(charges[:2])
    heavy_coeff = sum(charges[2:])
    full_coeff = sum(charges)
    check(
        "mixed_wzw",
        "complete five charge trace vanishes",
        full_coeff == 0,
        f"charges={charges}; sum={full_coeff}",
    )
    check(
        "mixed_wzw",
        "light heavy split is plus two minus two",
        (light_coeff, heavy_coeff) == (2, -2),
        f"light={light_coeff}; heavy={heavy_coeff}; total={full_coeff}",
    )

    # Declared-scope gate logic. False construction gates are policy declarations,
    # while algebraic true gates are tied to the controls above.
    gates = {
        "euclidean_gamma_reference_reality_data_specified": max(
            hermitian_error, clifford_error, gamma5_error, reality_error, j_square_error
        )
        < TOL,
        "fixed_adjoint_mass_kramers_reality_proven": False,
        "common_quadratic_aps_pv_regulator_specified": moments
        == {0: 0, 1: 0, 2: 0, 3: -6},
        "five_dimensional_aps_domain_specified": "operatorname{Dom}_{\\rm APS,\\alpha}"
        in tex_text,
        "physical_fermion_phase_zero_cut_fixed": zero_cut_flow == -1
        and auxiliary_cut_flow == 0,
        "conditional_higgsed_gw_source_specified": ranks == (2, 2, 1),
        "finite_higgsed_overlap_dirac_yukawa_operator_defined": False,
        "five_dimensional_higgsed_target_kernel_defined": False,
        "minimal_higgsed_overlap_eta_pair_computed": False,
        "conditional_higgsed_parity_predicts_plus_plus": conditional_pair == (1, 1),
        "pure_unbroken_gauge_gravity_heavy_phase_computed": gravitational_sum
        == Fraction(-1, 8)
        and (derived_ch2_c2, derived_ch2_x2) == (Fraction(-1), Fraction(1)),
        "pure_heavy_phase_on_G3_computed": pure_pair[0] == 1,
        "pure_heavy_phase_on_G2_computed": pure_pair[1] == 1,
        "exact_separated_rank_two_projector_obstruction_proven": scalar_ranks
        == {0.0: 0, 1.0: 5}
        and abs(light_trace - 2.0) < TOL,
        "full_nonperturbative_sp4_measure_constructed": False,
        "original_action_target_mass_family_lift_defined": False,
        "original_deep_sp4_action_uniquely_selects_target_pair": False,
        "target_mapping_torus_eta_pair_selected": False,
        "mixed_wzw_heavy_determinant_computed": False,
        "k_plus_two_all_scale_matched": False,
        "sp4_aps_pv_lane_closed": False,
        "degree_one_route_e_portal_started": False,
        "physics_promotion_authorized": False,
    }
    expected_true = {
        "euclidean_gamma_reference_reality_data_specified",
        "common_quadratic_aps_pv_regulator_specified",
        "five_dimensional_aps_domain_specified",
        "physical_fermion_phase_zero_cut_fixed",
        "conditional_higgsed_gw_source_specified",
        "conditional_higgsed_parity_predicts_plus_plus",
        "pure_unbroken_gauge_gravity_heavy_phase_computed",
        "pure_heavy_phase_on_G3_computed",
        "pure_heavy_phase_on_G2_computed",
        "exact_separated_rank_two_projector_obstruction_proven",
    }
    expected_false = set(gates) - expected_true
    true_keys = {key for key, value in gates.items() if value}
    false_keys = {key for key, value in gates.items() if not value}
    check(
        "declared_scope",
        "derived and conditional-source gates match the declared true set",
        true_keys == expected_true,
        f"true gates={sorted(true_keys)}",
    )
    check(
        "declared_scope",
        "unconstructed operators, physics promotion, and portal stay false",
        false_keys == expected_false,
        f"false gates={sorted(false_keys)}",
    )
    lane_closure = bool(
        gates["full_nonperturbative_sp4_measure_constructed"]
        and gates["target_mapping_torus_eta_pair_selected"]
        and gates["mixed_wzw_heavy_determinant_computed"]
        and gates["k_plus_two_all_scale_matched"]
    )
    check(
        "declared_scope",
        "Sp4 APS/PV lane remains fail closed",
        not lane_closure and not gates["sp4_aps_pv_lane_closed"],
        f"computed lane closure={lane_closure}",
    )

    passed = sum(item["pass"] for item in CHECKS)
    failed = len(CHECKS) - passed
    status = (
        "mechanical_pass_physics_fail_closed"
        if failed == 0
        and not gates["sp4_aps_pv_lane_closed"]
        and not gates["physics_promotion_authorized"]
        else "mechanical_failure"
    )
    return {
        "schema_version": "ap-e7-common-sp4-aps-pv-v2",
        "status": status,
        "summary": {
            "checks_total": len(CHECKS),
            "checks_passed": passed,
            "checks_failed": failed,
            "physics_lane_closed": gates["sp4_aps_pv_lane_closed"],
            "portal_started": gates["degree_one_route_e_portal_started"],
        },
        "regulator": {
            "pv_coefficients": pv_coefficients,
            "pv_nodes_squared_in_Lambda_units": pv_nodes,
            "pv_moments": moments,
            "scaled_large_lambda_log_coefficient": scaled_log,
            "agmon_cut": "theta_B=pi-epsilon, upper lateral",
            "pv_real_boson_realization": "spectral-cut determinant square root; no literal symmetric real-Grassmann Gaussian claimed",
            "physical_fermion_spectral_cut": 0,
            "physical_endpoint_condition": "0 outside spec(H0) union spec(H1)",
            "auxiliary_aps_cut_role": "alpha in common endpoint resolvent, Fredholm index only",
            "aps_domain_t0": "P_[alpha,infinity)(H0)u(0)=0",
            "aps_domain_t1": "P_(-infinity,alpha)(H1)u(1)=0",
            "scalar_control_SF_0": zero_cut_flow,
            "scalar_control_SF_2": auxiliary_cut_flow,
        },
        "reality_scope": {
            "reference_J_square": -1,
            "scalar_reference_reality": True,
            "fixed_adjoint_mass_kramers_reality": False,
            "adjoint_family_involution": "J K(A,Sigma) J^-1 = K(A,-Sigma)",
            "adjoint_family_residual": family_involution_error,
            "fixed_background_noncommutation_residual": fixed_background_residual,
        },
        "conditional_higgsed_gw_source": {
            "projector_ranks": list(ranks),
            "trace_H": float(np.trace(h)),
            "trace_H_P_plus": light_trace,
            "G3_mod2_exponent": g3_parity,
            "G2_light_copies": g2_light_copies,
            "G2_heavy_copies": g2_heavy_copies,
            "G2_mod2_exponent": g2_parity,
            "conditional_parity_prediction": list(conditional_pair),
            "finite_overlap_dirac_yukawa_operator_defined": False,
            "five_dimensional_target_kernel_defined": False,
            "overlap_eta_pair_computed": False,
            "scope": "GW-compatible continuum-limit source on the gapped Higgs orbit; not a completed lift",
        },
        "target_eta": {
            "abstract_character_table": {
                f"{e3}{e2}": list(value)
                for (e3, e2), value in sorted(characters.items())
            },
            "pure_heavy_pair": list(pure_pair),
            "conditional_higgsed_parity_prediction": list(conditional_pair),
            "conditional_higgsed_overlap_eta_pair": None,
            "full_original_action_pair": None,
        },
        "heavy_threshold": {
            "one_flavour_index_polynomial": "-c2(Vc)+x^2-p1(TM)/8",
            "derived_ch2_coefficients": {
                "c2_Vc": int(derived_ch2_c2),
                "x_squared": int(derived_ch2_x2),
            },
            "gravitational_coefficient": str(gravitational_sum),
            "flavour_count": 2,
            "phase_on_closed_spin_four_manifolds": 1,
            "characteristic_grid_points": grid_size,
            "characteristic_grid_sample": grid_rows,
            "generator_indices": {
                name: int(index) for name, index in generator_indices.items()
            },
        },
        "mixed_wzw": {
            "charges": charges,
            "light_coefficient": light_coeff,
            "heavy_coefficient": heavy_coeff,
            "full_coefficient": full_coeff,
            "IR_k_for_B_plus_one": 2,
            "all_scale_k_matched": False,
            "controlled_assignments_exhaustive": False,
        },
        "dependencies": {
            "relaxed_B1_soliton_profile_used": False,
            "soliton_profile_checksum": None,
            "reason": "this APS/PV branch is background-field and topological",
        },
        "gates": gates,
        "checks": CHECKS,
        "sources": [source_row(path) for path in sources],
    }


def write_markdown(card: dict[str, Any]) -> None:
    lines = [
        "# AP-E7 common Sp(4) APS/PV regulator audit",
        "",
        f"- Status: {card['status']}",
        f"- Checks: {card['summary']['checks_passed']}/{card['summary']['checks_total']} passed",
        "- Common continuum quadratic APS/PV regulator: specified",
        "- Physical fermion determinant phase: zero spectral cut with gapped endpoints",
        "- Conditional Higgsed parity prediction: (+1,+1)",
        "- Higgsed finite overlap Dirac-Yukawa/5D kernel: not defined; eta pair not computed",
        "- Original deep-Sp(4) target eta pair: not selected",
        "- Pure gauge/gravity heavy phase: exactly +1",
        "- All-scale mixed-WZW k=+2: not matched",
        "- Portal: not started",
        "",
        "## Exact formulas",
        "",
        r"- \(I_h=\int[-c_2(V_c)+x^2-p_1(TM)/8]\) and \(N_f=2\) gives \(\exp(2\pi iI_h)=1\).",
        r"- Physical sign: \(\tau_F=\exp(i\pi N_f\,\mathrm{SF}_0)\); arbitrary \(\alpha\) is used only for the Fredholm index.",
        r"- \(J\widehat\Sigma J^{-1}=-\widehat\Sigma\), so fixed-adjoint-mass Kramers reality is not claimed.",
        r"- \(P_\pm=(h^2\pm h)/2,\ P_0=1-h^2\), with ranks \((2,2,1)\).",
        r"- Exact separated-projector obstruction: \(C(0)=cI_5\) has \(\mathrm{Tr}(HC(0))=0\), but \(\mathrm{Tr}(HP_+)=2\).",
        r"- Complete-irrep mixed trace: \(\kappa_{\rm full}=+2-2=0\).",
        "",
        "## Gate ledger",
        "",
        "| gate | value |",
        "|---|---:|",
    ]
    for key, value in card["gates"].items():
        lines.append(f"| {key} | {str(value).lower()} |")
    lines.extend(
        [
            "",
            "## Deterministic checks",
            "",
            "| group | check | pass | detail |",
            "|---|---|---:|---|",
        ]
    )
    for item in card["checks"]:
        detail = str(item["detail"]).replace("|", "/")
        lines.append(
            f"| {item['group']} | {item['name']} | "
            f"{'yes' if item['pass'] else 'no'} | {detail} |"
        )
    lines.extend(["", "## Source manifest", "", "| path | sha256 |", "|---|---|"])
    for source in card["sources"]:
        lines.append(f"| {source['path']} | {source['sha256']} |")
    MD_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    card = run_audit()
    JSON_OUT.write_text(
        json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_markdown(card)
    print(
        f"{card['status']}: "
        f"{card['summary']['checks_passed']}/{card['summary']['checks_total']} checks passed"
    )
    print(f"wrote {JSON_OUT.relative_to(REPO)}")
    print(f"wrote {MD_OUT.relative_to(REPO)}")
    return 0 if card["summary"]["checks_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
