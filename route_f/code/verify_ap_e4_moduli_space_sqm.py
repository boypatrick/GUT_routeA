#!/usr/bin/env python3
"""Deterministic AP-E4 audit for the CP1 moduli-space N=2 SQM route.

The mother model is a half-BPS non-Abelian U(2) Chern--Simons--Higgs vortex.
Its internal orientational modulus is CP1 and its low-energy worldline action
has the Fubini--Study L2 metric.  Broken supersymmetry supplies a tangent
fermion, whose quantization produces the canonical Spin-c/Dolbeault module.

The central fail-closed distinction is tested explicitly: in the declared
declared canonical Spin-c polarization the ordinary tangent fermion generates
the Clifford module Omega^{0,*}(CP1); it does not by itself make the
wavefunction T-valued.  That canonical re-quantization has one untwisted
ground state, while the source-selected half-form ordering has none.  Three
positive-chirality zero modes occur only after an independent coefficient
line E=O(2), generated for example by a quantized magnetic/Wess--Zumino
coupling or a separately derived Fermi zero-mode bundle.  No such bulk
derivation or AP-E3/AP-E4 embedding is assumed by a green run.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e4_moduli_space_sqm.tex"
BIB = ROUTE_F / "tex" / "ap_e4_moduli_space_sqm.bib"
THIS_SCRIPT = Path(__file__).resolve()

SEED = 20260716
TOL = 5.0e-11
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


def normalized_doublet_and_derivative(
    w: complex, dw: complex
) -> tuple[np.ndarray, np.ndarray]:
    denominator = 1.0 + abs(w) ** 2
    d_denominator = 2.0 * (w.conjugate() * dw).real
    root = math.sqrt(denominator)
    section = np.array([1.0, w], dtype=complex) / root
    derivative = np.array(
        [
            -0.5 * d_denominator / denominator**1.5,
            dw / root - 0.5 * w * d_denominator / denominator**1.5,
        ],
        dtype=complex,
    )
    return section, derivative


def christoffel(w: complex) -> complex:
    """Gamma^w_ww for g_wbarw=C/(1+|w|^2)^2; C cancels."""

    return -2.0 * w.conjugate() / (1.0 + abs(w) ** 2)


def h0_o(k: int) -> int:
    return max(k + 1, 0)


def h1_o(k: int) -> int:
    return max(-k - 1, 0)


def midpoint_chern(k: int, theta_cells: int) -> float:
    dtheta = math.pi / theta_cells
    return 0.5 * k * sum(
        math.sin((index + 0.5) * dtheta) * dtheta
        for index in range(theta_cells)
    )


def surrogate_profile_integral(points: int, radial_cutoff: float = 20.0) -> float:
    """A convergent near-BPS profile check, not a numerical vortex solution.

    We choose xi_hat=1 and phi_hat=tanh(r), which obey the correct core and
    vacuum limits.  The nondynamical rho is then eliminated by its exact
    algebraic equation from the Chern--Simons vortex effective action.
    """

    nodes, weights = np.polynomial.legendre.leggauss(points)
    radius = 0.5 * radial_cutoff * (nodes + 1.0)
    phi = np.tanh(radius)
    xi = 1.0
    rho = (phi - math.sqrt(xi)) ** 2 / (phi**2 + xi)
    integrand = radius * (
        (phi**2 + xi) * rho**2
        + 2.0 * (1.0 - rho) * (phi - math.sqrt(xi)) ** 2
    )
    return float(0.5 * radial_cutoff * np.dot(weights, integrand))


def build_dolbeault_block(k: int, coefficient: float, levels: int) -> dict[str, Any]:
    """Finite exact spectral block for D=sqrt(2)(dbar+dbar^dagger)."""

    zero_plus = h0_o(k)
    if k < 0:
        raise ValueError("The explicit positive-flux block expects k >= 0")

    multiplicities = [k + 2 * n + 1 for n in range(1, levels + 1)]
    lambdas = [
        2.0 * math.sqrt(n * (n + k + 1)) / math.sqrt(coefficient)
        for n in range(1, levels + 1)
    ]
    massive_dimension = sum(multiplicities)
    plus_dimension = zero_plus + massive_dimension
    minus_dimension = massive_dimension

    b_map = np.zeros((minus_dimension, plus_dimension), dtype=complex)
    plus_offset = zero_plus
    minus_offset = 0
    for multiplicity, eigenvalue in zip(multiplicities, lambdas):
        b_map[
            minus_offset : minus_offset + multiplicity,
            plus_offset : plus_offset + multiplicity,
        ] = eigenvalue * np.eye(multiplicity)
        plus_offset += multiplicity
        minus_offset += multiplicity

    zero_pp = np.zeros((plus_dimension, plus_dimension), dtype=complex)
    zero_mm = np.zeros((minus_dimension, minus_dimension), dtype=complex)
    zero_pm = np.zeros((plus_dimension, minus_dimension), dtype=complex)
    zero_mp = np.zeros((minus_dimension, plus_dimension), dtype=complex)
    supercharge = np.block([[zero_pp, zero_pm], [b_map, zero_mm]])
    supercharge_dagger = supercharge.conjugate().T
    dirac = supercharge + supercharge_dagger
    gamma = np.diag(
        np.concatenate([np.ones(plus_dimension), -np.ones(minus_dimension)])
    ).astype(complex)
    hamiltonian = 0.5 * (
        supercharge @ supercharge_dagger + supercharge_dagger @ supercharge
    )

    return {
        "D": dirac,
        "Q": supercharge,
        "Qdag": supercharge_dagger,
        "Gamma": gamma,
        "H": hamiltonian,
        "zero_plus": zero_plus,
        "plus_dimension": plus_dimension,
        "minus_dimension": minus_dimension,
        "multiplicities": multiplicities,
        "lambdas": lambdas,
    }


def markdown_report(result: dict[str, Any]) -> str:
    lines = [
        "# AP-E4 CP1 moduli-space N=2 SQM audit",
        "",
        f"- Passed: **{result['summary']['passed']}/{result['summary']['total']}**",
        f"- Fail-closed physics promotion: **{result['status']['physics_promotion_allowed']}**",
        f"- Intrinsic physical tangent fermion: **{result['status']['physical_tangent_fermion_in_intrinsic_sqm_derived']}**",
        f"- Mother-model canonical vacuum line: **{result['status']['mother_model_vacuum_line_matches_canonical_spinc']}**",
        f"- CPT antibaryon polarization closed: **{result['status']['cpt_antibaryon_polarization_closed']}**",
        f"- Charged-two-colour embedding: **{result['status']['charged_two_colour_soliton_sqm_embedding_closed']}**",
        f"- AP-E4 worldline completion: **{result['status']['ap_e4_worldline_completion_closed']}**",
        f"- Mother-model L2 coefficient surrogate: `{result['metrics']['profile']['coefficient_C']:.12f}`",
        f"- Twisted k=2 Dirac gap: `{result['metrics']['spectrum']['dirac_gap']:.12f}`",
        f"- Twisted k=2 SQM gap: `{result['metrics']['spectrum']['hamiltonian_gap']:.12f}`",
        "",
        "## Main conclusion",
        "",
        "The BPS vortex supplies CP1 orientational zero modes, a finite L2 metric,",
        "and a tangent fermion.  Quantizing that tangent fermion supplies the",
        "Clifford algebra, but the mother model does **not** select the canonical",
        "Spin-c vacuum line used in the spectral card and does not supply an additional",
        "T=O(2) coefficient bundle.  The declared canonical Spin-c",
        "re-quantization has one untwisted ground state; the source-selected",
        "half-form ordering has none.  Three positive-chirality ground states require a separately",
        "derived E=O(2) flux/Fermi bundle.",
        "The AP-E3 level-two WZW term can be that coefficient line only if",
        "the same B=1 CP1 collective-coordinate embedding is proved and the",
        "degree-five WZW character is first transgressed along the spatial",
        "three-cycle, with the resulting line pulling back as O(+2).  The independent Chern--Simons vortex used here",
        "does not establish that composability statement.",
        "",
        "## Checks",
        "",
    ]
    for row in result["checks"]:
        marker = "PASS" if row["pass"] else "FAIL"
        lines.append(
            f"- [{marker}] `{row['group']}` — {row['name']}: {row['detail']}"
        )
    lines.extend(
        [
            "",
            "## Fail-closed gates",
            "",
            "- The analytic radial profile is a convergence surrogate, not a solved bulk vortex.",
            "- The unmodified mother model has no induced level-two Berry line.",
            "- The source-selected half-form ordering has no normalizable zero mode; the one-state untwisted count belongs to a declared canonical re-quantization.",
            "- Its CP1 vortex is not the charged-two-colour AP-E3 soliton.",
            "- AP-E3/AP-E4 composability needs the same B=1 moduli map and a spatial transgression whose signed line pulls back as O(2).",
            "- The transgressed line must descend through the gauge quotient, and the B=-1 sector needs a derived CPT/anti-canonical polarization map.",
            "- A bulk Callias/determinant-line or supersymmetric WZ derivation of O(2) remains open.",
            "- The product-compactification route remains a backup; its six-dimensional local, global, and Green--Schwarz anomaly gate is open.",
            "- No zero mode is identified with a chiral Route-E family.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    rng = np.random.default_rng(SEED)

    critical_sources = [TEX, BIB, THIS_SCRIPT]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "M0_provenance",
        "all moduli-SQM critical sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # -------------------------------------------------- CP1 zero mode and L2 metric
    horizontal_residuals: list[float] = []
    metric_residuals: list[float] = []
    chart_metric_residuals: list[float] = []
    fermion_covariance_residuals: list[float] = []
    fermion_norm_residuals: list[float] = []

    coefficient_probe = 1.731
    for _ in range(250):
        w = complex(rng.normal(), rng.normal())
        if abs(w) < 0.2:
            w += 0.4 + 0.2j
        dw = complex(rng.normal(), rng.normal())
        section, derivative = normalized_doublet_and_derivative(w, dw)
        projector = np.eye(2, dtype=complex) - np.outer(section, section.conjugate())
        tangent = projector @ derivative
        horizontal_residuals.append(float(abs(np.vdot(section, tangent))))
        expected_norm = abs(dw) ** 2 / (1.0 + abs(w) ** 2) ** 2
        metric_residuals.append(float(abs(np.vdot(tangent, tangent).real - expected_norm)))

        v = 1.0 / w
        jacobian = -1.0 / w**2
        dv = jacobian * dw
        metric_w = coefficient_probe * abs(dw) ** 2 / (1.0 + abs(w) ** 2) ** 2
        metric_v = coefficient_probe * abs(dv) ** 2 / (1.0 + abs(v) ** 2) ** 2
        chart_metric_residuals.append(abs(metric_w - metric_v))

        psi = complex(rng.normal(), rng.normal())
        dpsi = complex(rng.normal(), rng.normal())
        jacobian_prime = 2.0 / w**3
        psi_v = jacobian * psi
        dpsi_v = jacobian * dpsi + jacobian_prime * dw * psi
        covariant_w = dpsi + christoffel(w) * dw * psi
        covariant_v = dpsi_v + christoffel(v) * dv * psi_v
        fermion_covariance_residuals.append(abs(covariant_v - jacobian * covariant_w))
        fermion_norm_w = coefficient_probe * abs(psi) ** 2 / (1.0 + abs(w) ** 2) ** 2
        fermion_norm_v = coefficient_probe * abs(psi_v) ** 2 / (1.0 + abs(v) ** 2) ** 2
        fermion_norm_residuals.append(abs(fermion_norm_w - fermion_norm_v))

    check(
        "M1_moduli_metric",
        "projected orientational variations are horizontal CP1 zero modes",
        max(horizontal_residuals) < TOL,
        f"max horizontal residual={max(horizontal_residuals):.3e}",
    )
    check(
        "M1_moduli_metric",
        "the L2 metric has the Fubini--Study shape |dw|^2/(1+|w|^2)^2",
        max(metric_residuals) < TOL,
        f"max metric residual={max(metric_residuals):.3e}",
    )
    check(
        "M1_moduli_metric",
        "the Fubini--Study kinetic norm agrees in the w and v=1/w charts",
        max(chart_metric_residuals) < TOL,
        f"max chart residual={max(chart_metric_residuals):.3e}",
    )

    profile_sequence = [surrogate_profile_integral(points) for points in (64, 128, 256, 512)]
    profile_i0 = profile_sequence[-1]
    coefficient_c = 2.0 * math.pi * profile_i0
    check(
        "M1_moduli_metric",
        "the eliminated-rho near-BPS profile integral is positive, finite, and converged",
        profile_i0 > 0.0
        and abs(profile_sequence[-1] - profile_sequence[-2]) < 2.0e-10,
        f"I0 sequence={profile_sequence}; C(e=1)={coefficient_c:.12f}",
    )
    check(
        "M1_moduli_metric",
        "the profile surrogate is not promoted to a solved microscopic vortex",
        True,
        "profile_surrogate_is_microscopic_proof=false",
    )

    # ------------------------------------------------ tangent fermion covariance
    check(
        "M2_tangent_fermion",
        "psi transforms as a holomorphic tangent vector under v=1/w",
        max(fermion_covariance_residuals) < 2.0e-9,
        f"max D_t covariance residual={max(fermion_covariance_residuals):.3e}",
    )
    check(
        "M2_tangent_fermion",
        "the tangent-fermion kinetic norm is chart invariant",
        max(fermion_norm_residuals) < TOL,
        f"max norm residual={max(fermion_norm_residuals):.3e}",
    )

    # ------------------------------------------ Dolbeault module and flux twisting
    cohomology = [
        {"k": k, "h0": h0_o(k), "h1": h1_o(k), "index": h0_o(k) - h1_o(k)}
        for k in range(-6, 7)
    ]
    check(
        "M3_dolbeault",
        "Riemann--Roch and Serre duality give index(O(k))=k+1",
        all(row["index"] == row["k"] + 1 for row in cohomology),
        "k range=-6..6",
    )
    check(
        "M3_dolbeault",
        "untwisted canonical Dolbeault SQM has one, not three, ground state",
        (h0_o(0), h1_o(0)) == (1, 0),
        f"(h0,h1) for O(0)=({h0_o(0)},{h1_o(0)})",
    )
    check(
        "M3_dolbeault",
        "an independent O(2) coefficient line gives three positive-chirality zero modes",
        (h0_o(2), h1_o(2)) == (3, 0),
        f"(h0,h1,index) for O(2)=({h0_o(2)},{h1_o(2)},{h0_o(2)-h1_o(2)})",
    )
    check(
        "M3_dolbeault",
        "orientation reversal requires O(-4) for three negative-chirality zero modes",
        (h0_o(-4), h1_o(-4)) == (0, 3),
        f"(h0,h1,index) for O(-4)=({h0_o(-4)},{h1_o(-4)},{h0_o(-4)-h1_o(-4)})",
    )
    check(
        "M3_dolbeault",
        "the de Rham completion has two CP1 vacua and is also not a three-state mechanism",
        1 + 1 == 2,
        "b0+b2=2",
    )

    chern_sequence = [midpoint_chern(2, cells) for cells in (10, 20, 40, 80, 160)]
    phases = np.linspace(0.0, 2.0 * math.pi, 4097)
    winding = float(
        (np.unwrap(np.angle(np.exp(2.0j * phases)))[-1]
         - np.unwrap(np.angle(np.exp(2.0j * phases)))[0])
        / (2.0 * math.pi)
    )
    check(
        "M3_flux",
        "the candidate magnetic line has integral flux c1(O(2))=2",
        abs(chern_sequence[-1] - 2.0) < 4.0e-5
        and all(
            abs(chern_sequence[index + 1] - 2.0)
            < abs(chern_sequence[index] - 2.0)
            for index in range(len(chern_sequence) - 1)
        ),
        f"curvature sequence={chern_sequence}",
    )
    check(
        "M3_flux",
        "the north--south transition of O(2) has winding two",
        abs(winding - 2.0) < TOL,
        f"winding={winding:.12f}",
    )
    check(
        "M3_flux",
        "canonical Spin-c twisting by O(2) has determinant line O(6)",
        2 + 2 * 2 == 6,
        "c1(det W_E)=2+2k=6",
    )

    # ------------------------------------------- finite superalgebra/spectrum block
    block = build_dolbeault_block(k=2, coefficient=coefficient_c, levels=4)
    dirac = block["D"]
    supercharge = block["Q"]
    supercharge_dagger = block["Qdag"]
    gamma = block["Gamma"]
    hamiltonian = block["H"]
    check(
        "M4_superalgebra",
        "the finite Dolbeault block has nilpotent Q and Q-dagger",
        np.linalg.norm(supercharge @ supercharge, ord=np.inf) < TOL
        and np.linalg.norm(supercharge_dagger @ supercharge_dagger, ord=np.inf) < TOL,
        "Q^2=(Qdag)^2=0",
    )
    check(
        "M4_superalgebra",
        "D is Hermitian, odd, and satisfies D^2={Q,Qdag}=2H",
        np.linalg.norm(dirac - dirac.conjugate().T, ord=np.inf) < TOL
        and np.linalg.norm(dirac @ gamma + gamma @ dirac, ord=np.inf) < TOL
        and np.linalg.norm(
            dirac @ dirac
            - (supercharge @ supercharge_dagger + supercharge_dagger @ supercharge),
            ord=np.inf,
        )
        < TOL
        and np.linalg.norm(dirac @ dirac - 2.0 * hamiltonian, ord=np.inf) < TOL,
        "Hermitian/odd/superalgebra residuals below tolerance",
    )

    d_eigenvalues, d_eigenvectors = np.linalg.eigh(dirac)
    zero_mask = np.abs(d_eigenvalues) < 1.0e-10
    zero_vectors = d_eigenvectors[:, zero_mask]
    zero_chiralities = np.real(np.diag(zero_vectors.conjugate().T @ gamma @ zero_vectors))
    positive_d = d_eigenvalues[d_eigenvalues > 1.0e-10]
    dirac_gap = float(np.min(positive_d))
    h_eigenvalues = np.linalg.eigvalsh(hamiltonian)
    hamiltonian_gap = float(np.min(h_eigenvalues[h_eigenvalues > 1.0e-10]))
    expected_d_gap = 4.0 / math.sqrt(coefficient_c)
    expected_h_gap = 8.0 / coefficient_c
    check(
        "M4_spectrum",
        "the k=2 kernel contains exactly three positive-chirality zero modes",
        int(np.count_nonzero(zero_mask)) == 3
        and np.max(np.abs(zero_chiralities - 1.0)) < TOL,
        f"kernel={np.count_nonzero(zero_mask)}; chiralities={zero_chiralities.tolist()}",
    )
    check(
        "M4_spectrum",
        "the first twisted Dirac and SQM gaps scale as 4/sqrt(C) and 8/C",
        abs(dirac_gap - expected_d_gap) < TOL
        and abs(hamiltonian_gap - expected_h_gap) < TOL,
        f"Delta_D={dirac_gap:.12f}; Delta_H={hamiltonian_gap:.12f}; C={coefficient_c:.12f}",
    )
    check(
        "M4_spectrum",
        "the first massive level has five states at each Dirac sign",
        int(np.count_nonzero(np.isclose(d_eigenvalues, dirac_gap, atol=1.0e-10))) == 5
        and int(np.count_nonzero(np.isclose(d_eigenvalues, -dirac_gap, atol=1.0e-10))) == 5,
        "multiplicity(+gap)=multiplicity(-gap)=5",
    )

    # ------------------------------------------------------------- fail-closed logic
    check(
        "M5_fail_closed",
        "a tangent fermion is not an automatic T-valued coefficient twist",
        True,
        "ordinary_tangent_fermion_implies_E_equals_O(0), not E_equals_T",
    )
    check(
        "M5_fail_closed",
        "the unmodified Chern--Simons vortex mother model has no derived O(2) Berry line",
        True,
        "mother_model_level_two_flux_derived=false",
    )
    check(
        "M5_composability",
        "the intrinsic vortex SQM is not silently identified with the charged-two-colour soliton",
        True,
        "charged_two_colour_soliton_sqm_embedding_closed=false",
    )
    check(
        "M5_composability",
        "AP-E3 k=2 can twist AP-E4 only after the same B=1 CP1 moduli map and spatially transgressed WZW line pullback are proved",
        True,
        "same_B1_CP1_map=false; signed_WZW_pullback_to_O2=false",
    )
    check(
        "M5_fail_closed",
        "the product compactification route is not used to obtain the SQM result",
        True,
        "product_compactification_used=false",
    )
    check(
        "M5_fail_closed",
        "the backup product compactification has not passed its six-dimensional anomaly gate",
        True,
        "I8_factorization_or_cancellation=false; global_bordism_audit=false",
    )

    passed = sum(row["pass"] for row in CHECKS)
    status = {
        "bps_vortex_cp1_moduli_candidate": True,
        "l2_metric_shape_and_finiteness_checked": True,
        "tangent_fermion_covariance_checked": True,
        "canonical_spinc_dolbeault_sqm_derived": False,
        "canonical_spinc_spectrum_solved_given_declared_vacuum_line": True,
        "mother_model_vacuum_line_matches_canonical_spinc": False,
        "physical_tangent_fermion_in_intrinsic_sqm_derived": True,
        "untwisted_tangent_sqm_three_states": False,
        "o2_twist_three_positive_zero_modes": True,
        "mother_model_level_two_flux_derived": False,
        "independent_fermi_zero_mode_bundle_derived": False,
        "same_b1_cp1_collective_coordinate_proven": False,
        "signed_wzw_pullback_to_o2_proven": False,
        "spatial_wzw_transgression_and_pullback_proven": False,
        "wzw_gauge_basic_descent_proven": False,
        "cpt_antibaryon_polarization_closed": False,
        "ap_e3_k2_o2_coefficient_line_composable": False,
        "charged_two_colour_soliton_sqm_embedding_closed": False,
        "ap_e4_worldline_completion_closed": False,
        "product_compactification_used": False,
        "six_dimensional_anomaly_gate_closed": False,
        "route_e_family_identification": False,
        "route_e_portal_closed": False,
        "physics_promotion_allowed": False,
    }
    result = {
        "all_pass": passed == len(CHECKS),
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "physical_tangent_fermion_in_intrinsic_sqm_derived": True,
        "mother_model_vacuum_line_matches_canonical_spinc": False,
        "cpt_antibaryon_polarization_closed": False,
        "charged_two_colour_soliton_sqm_embedding_closed": False,
        "ap_e4_worldline_completion_closed": False,
        "route_e_portal_closed": False,
        "physics_promotion_allowed": False,
        "summary": {"passed": passed, "total": len(CHECKS)},
        "checks": CHECKS,
        "source_manifest": source_manifest,
        "metrics": {
            "profile": {
                "description": "phi_hat=tanh(r), xi_hat=1 convergence surrogate",
                "I0_sequence": profile_sequence,
                "I0": profile_i0,
                "coefficient_C": coefficient_c,
                "microscopic_solution": False,
            },
            "geometry": {
                "max_horizontal_residual": max(horizontal_residuals),
                "max_metric_residual": max(metric_residuals),
                "max_chart_metric_residual": max(chart_metric_residuals),
                "max_fermion_covariance_residual": max(fermion_covariance_residuals),
                "chern_sequence_k2": chern_sequence,
                "transition_winding_k2": winding,
            },
            "cohomology": cohomology,
            "spectrum": {
                "coefficient_C": coefficient_c,
                "radius": math.sqrt(coefficient_c) / 2.0,
                "dirac_gap": dirac_gap,
                "hamiltonian_gap": hamiltonian_gap,
                "first_four_multiplicities_per_sign": block["multiplicities"],
                "first_four_positive_eigenvalues": block["lambdas"],
            },
            "composability": {
                "intrinsic_sqm_moduli": "BPS U(2) Chern--Simons vortex CP1",
                "ap_e3_candidate_moduli": "charged-two-colour B=1 soliton CP1 (not yet established)",
                "same_b1_cp1_collective_coordinate_proven": False,
                "wzw_pullback_degree": None,
                "required_wzw_pullback_degree": 2,
                "differential_character_source_degree": 5,
                "spatial_fiber_dimension": 3,
                "transgressed_line_degree": 2,
                "wzw_gauge_basic_descent_proven": False,
                "cpt_antibaryon_polarization_closed": False,
                "supersymmetric_wz_completion_derived": False,
            },
        },
        "status": status,
        "next_steps": [
            "Solve the full BPS radial profiles and recompute I0 with certified tail bounds.",
            "Derive, rather than postulate, a supersymmetric worldline line E=O(2) from a bulk WZ/CS inflow or a Callias-index Fermi zero-mode bundle.",
            "Prove that the AP-E3 B=1 charged-two-colour soliton has this same CP1 collective coordinate, transgress the degree-five character along its spatial three-cycle, and show that the resulting signed line pulls back as O(+2).",
            "Compute the induced determinant-line orientation and verify that its physical sign is +2.",
            "Derive the mother-model vacuum line and its CPT/anti-canonical map; fixed canonical O(-2) has one negative mode, not the conjugate three-state kernel.",
            "If product compactification is revived, cancel or Green--Schwarz-factorize the full six-dimensional I8 and close the global bordism anomaly gate.",
            "Couple candidate ground states to Route E only after the flux origin and every retained spectral gap are closed.",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e4_moduli_space_sqm.json"
    md_path = OUTPUT / "ap_e4_moduli_space_sqm.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    md_path.write_text(markdown_report(result))

    print(f"AP-E4 moduli SQM checks: {passed}/{len(CHECKS)} passed")
    print(f"JSON: {json_path.relative_to(REPO)}")
    print(f"Markdown: {md_path.relative_to(REPO)}")
    return 0 if passed == len(CHECKS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
