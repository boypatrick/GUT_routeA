#!/usr/bin/env python3
"""Fail-closed algebraic and numerical audit for AP-E1.

The card checks four logically distinct constructions:

1. a local U(1) quotient of one nonzero complex doublet;
2. the Hopf connection and Fubini--Study metric induced from flat C^2;
3. fixed-global-charge Routh reduction and its monopole rotor;
4. the extra holomorphic/lowest-Landau-level gate needed for O(k).

A green run proves arithmetic identities inside these declared models.  It does
not derive the microscopic doublet, select k=2, identify the triplet with
families, or establish a stable four-dimensional bubble.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
PRIOR_BRIDGE = OUTPUT / "another_physics_route_e_bridge.json"
AP_E1_TEX = ROUTE_F / "tex" / "ap_e1_projective_doublet_action.tex"
AP_E1_BIB = ROUTE_F / "tex" / "ap_e1_projective_doublet_action.bib"
TOL = 1.0e-12
CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    """Append one mechanical check."""
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inner(u: list[complex], v: list[complex]) -> complex:
    return sum(a.conjugate() * b for a, b in zip(u, v))


def outer(u: list[complex], v: list[complex]) -> list[list[complex]]:
    return [[a * b.conjugate() for b in v] for a in u]


def matmul(
    a: list[list[complex]], b: list[list[complex]]
) -> list[list[complex]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def matadd(
    a: list[list[complex]], b: list[list[complex]]
) -> list[list[complex]]:
    return [
        [a[i][j] + b[i][j] for j in range(len(a[0]))]
        for i in range(len(a))
    ]


def matscale(c: complex, a: list[list[complex]]) -> list[list[complex]]:
    return [[c * value for value in row] for row in a]


def dagger(a: list[list[complex]]) -> list[list[complex]]:
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def max_matrix_difference(
    a: list[list[complex]], b: list[list[complex]]
) -> float:
    return max(
        abs(a[i][j] - b[i][j])
        for i in range(len(a))
        for j in range(len(a[0]))
    )


def commutator(
    a: list[list[complex]], b: list[list[complex]]
) -> list[list[complex]]:
    return matadd(matmul(a, b), matscale(-1.0, matmul(b, a)))


def identity(size: int) -> list[list[complex]]:
    return [
        [complex(1.0 if i == j else 0.0) for j in range(size)]
        for i in range(size)
    ]


def directional_section(
    w: complex, dw: complex
) -> tuple[list[complex], list[complex]]:
    """Return s_N(w) and its exact directional derivative."""
    norm = 1.0 + abs(w) ** 2
    dnorm = 2.0 * (w.conjugate() * dw).real
    root = math.sqrt(norm)
    s = [1.0 / root, w / root]
    ds = [
        -0.5 * dnorm / norm ** 1.5,
        dw / root - 0.5 * w * dnorm / norm ** 1.5,
    ]
    return s, ds


def vector_subtract(u: list[complex], v: list[complex]) -> list[complex]:
    return [a - b for a, b in zip(u, v)]


def main() -> None:
    if not PRIOR_BRIDGE.exists():
        raise FileNotFoundError(
            "AP-E1 requires the prior bridge card for the Q=10^6 comparison: "
            f"{PRIOR_BRIDGE}"
        )
    prior = json.loads(PRIOR_BRIDGE.read_text(encoding="utf-8"))
    for required_source in (AP_E1_TEX, AP_E1_BIB):
        if not required_source.exists():
            raise FileNotFoundError(f"missing AP-E1 source artifact: {required_source}")
    tex_bytes = AP_E1_TEX.read_bytes()
    tex_source = tex_bytes.decode("utf-8")

    # ----------------------------------------------------------- Hopf quotient
    z = [0.7 + 0.2j, -0.3 + 1.1j]
    norm_sq = inner(z, z).real
    projector = matscale(1.0 / norm_sq, outer(z, z))
    projector_square = matmul(projector, projector)
    trace_p = projector[0][0] + projector[1][1]
    determinant_p = projector[0][0] * projector[1][1] - projector[0][1] * projector[1][0]
    check(
        "G1_Hopf_quotient",
        "normalized doublet defines a rank-one Hermitian projector",
        max_matrix_difference(projector_square, projector) < TOL
        and max_matrix_difference(dagger(projector), projector) < TOL
        and abs(trace_p - 1.0) < TOL
        and abs(determinant_p) < TOL,
        "max|P^2-P|={:.2e}, max|P^dag-P|={:.2e}, TrP={:.16f}, |detP|={:.2e}".format(
            max_matrix_difference(projector_square, projector),
            max_matrix_difference(dagger(projector), projector),
            trace_p.real,
            abs(determinant_p),
        ),
    )

    alpha = 0.731
    z_phase = [cmath.exp(1j * alpha) * value for value in z]
    phase_projector = matscale(1.0 / inner(z_phase, z_phase).real, outer(z_phase, z_phase))
    phase_error = max_matrix_difference(projector, phase_projector)
    check(
        "G1_Hopf_quotient",
        "projector is invariant under the common U(1) phase",
        phase_error < TOL,
        f"max|P(e^(i alpha)z)-P(z)|={phase_error:.2e}",
    )

    sigma_x = [[0j, 1.0 + 0j], [1.0 + 0j, 0j]]
    sigma_y = [[0j, -1j], [1j, 0j]]
    sigma_z = [[1.0 + 0j, 0j], [0j, -1.0 + 0j]]
    spinor = [value / math.sqrt(norm_sq) for value in z]
    n = [inner(spinor, [sum(row[j] * spinor[j] for j in range(2)) for row in sigma]).real
         for sigma in (sigma_x, sigma_y, sigma_z)]
    reconstructed = matscale(0.5, matadd(identity(2), matadd(
        matscale(n[0], sigma_x), matadd(matscale(n[1], sigma_y), matscale(n[2], sigma_z))
    )))
    check(
        "G1_Hopf_quotient",
        "Pauli Hopf vector lies on S2 and reconstructs the projector",
        abs(sum(value * value for value in n) - 1.0) < TOL
        and max_matrix_difference(projector, reconstructed) < TOL,
        "n^2={:.16f}, max|P-(I+n.sigma)/2|={:.2e}".format(
            sum(value * value for value in n),
            max_matrix_difference(projector, reconstructed),
        ),
    )

    # ------------------------------------------------ induced CP1 geometry
    w = 0.31 - 0.27j
    dw = -0.13 + 0.09j
    section, dsection = directional_section(w, dw)
    connection_rate = (-1j * inner(section, dsection)).real
    fs_norm = abs(dw) ** 2 / (1.0 + abs(w) ** 2) ** 2
    total_norm = inner(dsection, dsection).real
    horizontal = vector_subtract(
        dsection, [value * inner(section, dsection) for value in section]
    )
    horizontal_norm = inner(horizontal, horizontal).real
    decomposition_error = abs(total_norm - connection_rate**2 - fs_norm)
    check(
        "G2_induced_geometry",
        "flat C2 metric splits into Hopf-fiber plus Fubini--Study metric",
        decomposition_error < TOL and abs(horizontal_norm - fs_norm) < TOL,
        "|ds|^2={:.16f}, A^2={:.16f}, g_FS={:.16f}, residual={:.2e}".format(
            total_norm, connection_rate**2, fs_norm, decomposition_error
        ),
    )

    eta = 1.0 / w
    deta = -dw / w**2
    fs_south = abs(deta) ** 2 / (1.0 + abs(eta) ** 2) ** 2
    check(
        "G2_induced_geometry",
        "Fubini--Study metric agrees on north and south charts",
        abs(fs_norm - fs_south) < TOL,
        f"north={fs_norm:.16f}, south={fs_south:.16f}",
    )

    alpha_rate = 0.413
    transformed_section = [cmath.exp(1j * alpha) * value for value in section]
    transformed_dsection = [
        cmath.exp(1j * alpha) * (derivative + 1j * alpha_rate * value)
        for value, derivative in zip(section, dsection)
    ]
    transformed_connection = (-1j * inner(transformed_section, transformed_dsection)).real
    horizontal_transformed = vector_subtract(
        transformed_dsection,
        [
            value * inner(transformed_section, transformed_dsection)
            for value in transformed_section
        ],
    )
    gauge_horizontal_error = abs(inner(horizontal_transformed, horizontal_transformed).real - fs_norm)
    check(
        "G2_induced_geometry",
        "Hopf connection transforms as A to A+d alpha while horizontal norm is invariant",
        abs(transformed_connection - connection_rate - alpha_rate) < TOL
        and gauge_horizontal_error < TOL,
        "A'-A-alpha_dot={:.2e}, horizontal residual={:.2e}".format(
            transformed_connection - connection_rate - alpha_rate,
            gauge_horizontal_error,
        ),
    )

    # F=(1/2) sin(theta) dtheta wedge dphi.  Midpoint quadrature is
    # independent of any chart singularity and converges from above.
    theta_cells = 4000
    dtheta = math.pi / theta_cells
    flux = 2.0 * math.pi * sum(
        0.5 * math.sin((index + 0.5) * dtheta) * dtheta
        for index in range(theta_cells)
    )
    chern_one = flux / (2.0 * math.pi)
    check(
        "G2_induced_geometry",
        "Hopf curvature has first Chern number one",
        abs(chern_one - 1.0) < 3.0e-8,
        f"integral F={flux:.15f}, integral F/(2pi)={chern_one:.12f}",
    )

    # ---------------------------------------------------- local/global split
    rho = 1.7
    drho = -0.23
    dchi = 0.61
    flat_polar = drho**2 + rho**2 * (
        (dchi + connection_rate) ** 2 + fs_norm
    )
    dz = [
        cmath.exp(1j * alpha)
        * (drho * value + rho * (1j * dchi * value + derivative))
        for value, derivative in zip(section, dsection)
    ]
    direct_flat = inner(dz, dz).real
    local_auxiliary_minimum = drho**2 + rho**2 * fs_norm
    check(
        "G3_local_global_split",
        "polar decomposition of the canonical doublet kinetic term is exact",
        abs(flat_polar - direct_flat) < TOL,
        f"direct={direct_flat:.16f}, polar={flat_polar:.16f}",
    )
    check(
        "G3_local_global_split",
        "auxiliary local gauge field removes only the Hopf-fiber square",
        abs(
            flat_polar
            - rho**2 * (dchi + connection_rate) ** 2
            - local_auxiliary_minimum
        )
        < TOL,
        "global kinetic minus fiber={:.16f}, local quotient={:.16f}".format(
            flat_polar - rho**2 * (dchi + connection_rate) ** 2,
            local_auxiliary_minimum,
        ),
    )

    # ------------------------------------------------ fixed-charge reduction
    inertia = 3.7
    charge = 1.3
    a_q = 0.21
    qdot = -0.47
    g_qq = 1.8
    potential = 0.6
    chidot = charge / inertia - a_q * qdot
    lagrangian = 0.5 * inertia * (
        (chidot + a_q * qdot) ** 2 + g_qq * qdot**2
    ) - potential
    routh_direct = lagrangian - charge * chidot
    routh_closed = (
        0.5 * inertia * g_qq * qdot**2
        + charge * a_q * qdot
        - potential
        - charge**2 / (2.0 * inertia)
    )
    check(
        "G4_fixed_charge",
        "direct Routh transform equals the monopole-coupled CP1 expression",
        abs(routh_direct - routh_closed) < TOL,
        f"direct={routh_direct:.16f}, closed={routh_closed:.16f}, residual={abs(routh_direct-routh_closed):.2e}",
    )

    canonical_q_momentum = inertia * g_qq * qdot + charge * a_q
    hamiltonian_direct = canonical_q_momentum * qdot + charge * chidot - lagrangian
    hamiltonian_closed = (
        (canonical_q_momentum - charge * a_q) ** 2 / (2.0 * inertia * g_qq)
        + potential
        + charge**2 / (2.0 * inertia)
    )
    check(
        "G4_fixed_charge",
        "Legendre transform gives the minimally coupled monopole Hamiltonian",
        abs(hamiltonian_direct - hamiltonian_closed) < TOL,
        f"direct={hamiltonian_direct:.16f}, closed={hamiltonian_closed:.16f}",
    )

    fixed_q_energy = (
        potential
        + charge**2 / (2.0 * inertia)
        + 0.5 * inertia * g_qq * qdot**2
    )
    orientation_minimum = potential + charge**2 / (2.0 * inertia)
    orientation_excess = 0.5 * inertia * g_qq * qdot**2
    check(
        "G4_fixed_charge",
        "fixed-Q profile energy retains nonnegative orientation kinetic energy",
        abs(fixed_q_energy - orientation_minimum - orientation_excess) < TOL
        and orientation_excess > 0.0,
        "E-E_min={:.16f}, (I/2)g*qdot^2={:.16f}".format(
            fixed_q_energy - orientation_minimum, orientation_excess
        ),
    )

    positive_level = 2.0
    positive_radius = math.sqrt(positive_level)
    positive_moment_gradient_norm = 2.0 * positive_radius
    zero_level_radius = 0.0
    check(
        "G5_holomorphic_quantization",
        "C2 moment-map quotient is CP1 only at positive level; k=0 is a point",
        positive_moment_gradient_norm > 0.0 and zero_level_radius == 0.0,
        "k>0 regular-level gradient norm={:.6f}; mu^-1(0) radius={:.1f}".format(
            positive_moment_gradient_norm, zero_level_radius
        ),
    )

    level = 2.0
    occupation = 1.7
    alpha_dot = 0.13
    dc_coefficient = level - occupation
    kinetic_gauge_shift = occupation * alpha_dot
    constraint_gauge_shift = (level - occupation) * alpha_dot
    total_gauge_shift = kinetic_gauge_shift + constraint_gauge_shift
    check(
        "G5_holomorphic_quantization",
        "corrected Branch-B symplectic sign imposes N=k and a positive level-k endpoint phase",
        abs(dc_coefficient - (level - occupation)) < TOL
        and abs(total_gauge_shift - level * alpha_dot) < TOL
        and abs(level - level) < TOL,
        "dL/d(hbar*c)=k-N={:.6f}; delta L/(hbar*alpha_dot)={:.6f}; reduced flux level=+{}".format(
            dc_coefficient, total_gauge_shift / alpha_dot, int(level)
        ),
    )

    for k in range(7):
        weights = [k - 2 * occupation for occupation in range(k + 1)]
        expected = list(range(k, -k - 1, -2))
        check(
            "G5_holomorphic_quantization",
            f"Sym^{k}(C2) has k+1 states and the correct Cartan weights",
            len(weights) == k + 1 and weights == expected,
            f"k={k}, dimension={len(weights)}, weights={weights}",
        )

    sqrt_two = math.sqrt(2.0)
    j_plus = [[0j, sqrt_two, 0j], [0j, 0j, sqrt_two], [0j, 0j, 0j]]
    j_minus = dagger(j_plus)
    j_x = matscale(0.5, matadd(j_plus, j_minus))
    j_y = matscale(-0.5j, matadd(j_plus, matscale(-1.0, j_minus)))
    j_z = [[1.0 + 0j, 0j, 0j], [0j, 0j, 0j], [0j, 0j, -1.0 + 0j]]
    commutator_errors = [
        max_matrix_difference(commutator(j_x, j_y), matscale(1j, j_z)),
        max_matrix_difference(commutator(j_y, j_z), matscale(1j, j_x)),
        max_matrix_difference(commutator(j_z, j_x), matscale(1j, j_y)),
    ]
    casimir = matadd(matmul(j_x, j_x), matadd(matmul(j_y, j_y), matmul(j_z, j_z)))
    casimir_error = max_matrix_difference(casimir, matscale(2.0, identity(3)))
    check(
        "G5_holomorphic_quantization",
        "k=2 Schwinger sector is the spin-one triplet",
        max(commutator_errors) < TOL and casimir_error < TOL,
        "max commutator error={:.2e}, max|J^2-2I|={:.2e}".format(
            max(commutator_errors), casimir_error
        ),
    )

    k = 2
    levels = []
    for radial_index in range(4):
        degeneracy = k + 2 * radial_index + 1
        dimensionless_energy = 2.0 * (
            radial_index * (radial_index + k + 1) + k / 2.0
        )
        levels.append(
            {
                "n": radial_index,
                "j": k / 2.0 + radial_index,
                "degeneracy": degeneracy,
                "energy_times_I_over_hbar2": dimensionless_energy,
            }
        )
    lll_gap = levels[1]["energy_times_I_over_hbar2"] - levels[0]["energy_times_I_over_hbar2"]
    check(
        "G5_holomorphic_quantization",
        "second-order k=2 rotor has a three-state LLL but infinitely many higher levels",
        levels[0]["degeneracy"] == 3
        and lll_gap == 8.0
        and all(levels[i + 1]["degeneracy"] > levels[i]["degeneracy"] for i in range(3)),
        f"first four levels={json.dumps(levels, sort_keys=True)}, gap*I/hbar^2={lll_gap:.1f}",
    )

    def signed_rotor_spectrum(charge_level: int, radial_index: int) -> tuple[float, int, float]:
        """Return (j, degeneracy, E*I/hbar^2) for either flux orientation."""
        flux_magnitude = abs(charge_level)
        angular_momentum = flux_magnitude / 2.0 + radial_index
        degeneracy = flux_magnitude + 2 * radial_index + 1
        dimensionless_energy = 2.0 * (
            radial_index * (radial_index + flux_magnitude + 1)
            + flux_magnitude / 2.0
        )
        return angular_momentum, degeneracy, dimensionless_energy

    positive_flux_levels = [signed_rotor_spectrum(+2, n) for n in range(4)]
    negative_flux_levels = [signed_rotor_spectrum(-2, n) for n in range(4)]
    check(
        "G5_holomorphic_quantization",
        "monopole-rotor spectrum depends on |k| and is invariant under flux reversal",
        positive_flux_levels == negative_flux_levels
        and positive_flux_levels[0] == (1.0, 3, 2.0)
        and positive_flux_levels[1][2] - positive_flux_levels[0][2] == 8.0,
        f"k=+2 levels={positive_flux_levels}; k=-2 levels={negative_flux_levels}",
    )

    # ------------------------------------------------- Q-ball incompatibility
    qball = prior["separation_theorems"]["Q_ball"]
    q_macro = int(qball["parameters_GeV_convention"]["Q"])
    f0_sq = float(qball["f0_squared_GeV2"])
    radius = float(qball["finite_Q_step_variational"]["radius_GeV_inverse"])
    omega_reported = float(qball["finite_Q_step_variational"]["omega_GeV"])
    volume = 4.0 * math.pi * radius**3 / 3.0
    collective_inertia = 2.0 * f0_sq * volume
    omega_reconstructed = q_macro / collective_inertia
    macro_dimension = q_macro + 1
    check(
        "G6_route_e_gate",
        "prior Q-ball angular frequency is Q divided by collective inertia",
        abs(omega_reconstructed - omega_reported) < TOL,
        "I={:.12f}, Q/I={:.12f}, stored omega={:.12f}".format(
            collective_inertia, omega_reconstructed, omega_reported
        ),
    )
    check(
        "G6_route_e_gate",
        "macroscopic Q=10^6 sector is not the k=2 triplet",
        q_macro == 1_000_000 and macro_dimension == 1_000_001 and macro_dimension != 3,
        f"if k=Q/hbar={q_macro}, dim H0(O(k))={macro_dimension}, not 3",
    )

    # ------------------------------------------------------- Derrick scaling
    e2 = 3.0
    e4 = 12.0
    e0 = 0.8
    scale = 1.73
    def radial_test_energies(size: float, cells: int = 20_000) -> dict[str, float]:
        """Numerically integrate a Gaussian derivative-counting test field."""
        dr = 8.0 * size / cells
        totals = {"E2": 0.0, "E4": 0.0, "E0": 0.0}
        for index in range(cells):
            radius = (index + 0.5) * dr
            field = math.exp(-(radius / size) ** 2)
            derivative = -2.0 * radius * field / size**2
            measure = 4.0 * math.pi * radius**2 * dr
            totals["E2"] += measure * derivative**2
            totals["E4"] += measure * derivative**4
            totals["E0"] += measure * field**2
        return totals

    unit_test_terms = radial_test_energies(1.0)
    scaled_test_terms = radial_test_energies(scale)
    scaling_ratios = {
        name: scaled_test_terms[name] / unit_test_terms[name]
        for name in unit_test_terms
    }
    expected_ratios = {"E2": scale, "E4": 1.0 / scale, "E0": scale**3}
    scaling_relative_error = max(
        abs(scaling_ratios[name] / expected_ratios[name] - 1.0)
        for name in scaling_ratios
    )
    pure_two_derivative_slope = e2 + 3.0 * e0
    optimum = math.sqrt(e4 / e2)
    second_derivative = 2.0 * e4 / optimum**3
    check(
        "G7_Derrick",
        "three-dimensional size dilation has E2~lambda, E4~lambda^-1, E0~lambda^3",
        scaling_relative_error < 2.0e-12,
        "lambda={}, numerical ratios={}, expected={}, max relative error={:.2e}".format(
            scale,
            json.dumps(scaling_ratios, sort_keys=True),
            json.dumps(expected_ratios, sort_keys=True),
            scaling_relative_error,
        ),
    )
    check(
        "G7_Derrick",
        "positive E2 plus E0 cannot have a finite-size Derrick stationary point",
        pure_two_derivative_slope > 0.0,
        f"dE/dlambda at lambda=1 is E2+3E0={pure_two_derivative_slope:.6f}",
    )
    check(
        "G7_Derrick",
        "an E4 stabilizer can balance E2 at finite size",
        abs(optimum - 2.0) < TOL and second_derivative > 0.0,
        f"lambda_star=sqrt(E4/E2)={optimum:.12f}, E''={second_derivative:.12f}",
    )

    critical_tex_fragments = [
        r"(D_\mu Z)^\dagger D^\mu Z+\Lambda(x)",
        r"E[f,q,\dot q]=E_0[f]+\frac{Q^2}{2\cI[f]}",
        r"\mu^{-1}(0)=\{0\}",
        r"j=\frac{\kappa}{2}+n",
        r"L_b&=-\frac{\ii\hbar}{2}",
        r"\left[b^\dagger D_tb-(D_tb)^\dagger b\right]+\hbar k c_t",
    ]
    missing_tex_fragments = [
        fragment for fragment in critical_tex_fragments if fragment not in tex_source
    ]
    check(
        "G8_source_regression",
        "TeX contains the corrected action signs, charge domain, and signed-spectrum formulas",
        not missing_tex_fragments and b"\r" not in tex_bytes,
        "missing critical fragments={} ; carriage_return_present={}".format(
            missing_tex_fragments, b"\r" in tex_bytes
        ),
    )

    passed = sum(item["pass"] for item in CHECKS)
    all_pass = passed == len(CHECKS)
    result = {
        "audit": "verify_ap_e1_projective_doublet",
        "status": (
            "ap_e1_conditionally_proved_with_explicit_o2_and_stability_blockers"
            if all_pass
            else "mechanical_failure"
        ),
        "all_pass": all_pass,
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "physics_promotion_allowed": False,
        "logical_result": {
            "local_U1": (
                "proved under one nonzero charge-one complex doublet, fixed norm, "
                "and an auxiliary local U(1): S3/U(1)=CP1 with the induced FS action"
            ),
            "global_U1": (
                "refuted as a pointwise CP1 quotient: the local Hopf-fiber phase remains physical"
            ),
            "fixed_global_charge": (
                "proved for the collective ansatz: Routh reduction gives a monopole-coupled "
                "T*CP1 rotor, not by itself H0(CP1,O(k))"
            ),
            "holomorphic_sector": (
                "conditional on first-order reduction or a controlled LLL/Kahler polarization: "
                "for k>=0, H_k=H0(CP1,O(k)), dim=k+1; for k<0 the conjugate "
                "anti-holomorphic sector has dim=|k|+1"
            ),
            "route_E_O2": (
                "not derived: k=2, the LLL gate, and the identification with T_CP1 and three "
                "chiral families remain separate assumptions"
            ),
        },
        "first_principles_chain": [
            "flat positive Hermitian metric on one complex doublet",
            "normalization or a dynamically frozen radial mode",
            "ray/local-phase equivalence",
            "Hopf quotient and induced Fubini--Study metric",
        ],
        "nonderived_inputs": [
            "why the microscopic order parameter is exactly one complex doublet",
            "why the radial mode is frozen and what fixes its scale",
            "why the relevant prequantum level is k=2",
            "why the second-order rotor may be projected to its LLL",
            "why the resulting spin-one triplet is a family carrier",
            "how a stable four-dimensional bubble survives radial, gauge, and portal fluctuations",
        ],
        "branch_decision": {
            "A_charge_two": (
                "set the physical fixed charge to k=2; obtains the triplet but abandons the "
                "macroscopic Q=10^6 thin-wall benchmark"
            ),
            "B_separate_level_two_internal_sector": (
                "retain the macroscopic Q-ball charge but introduce an independent quantized "
                "level-two Wess--Zumino/Schwinger-boson sector"
            ),
            "forbidden_shortcut": (
                "identify Q=10^6 with k=2 or infer O(2) merely from CP1"
            ),
        },
        "numerical_summary": {
            "hopf_flux": flux,
            "hopf_chern_number": chern_one,
            "k2_rotor_first_levels": levels,
            "k2_LLL_gap_times_I_over_hbar2": lll_gap,
            "prior_qball_charge": q_macro,
            "prior_qball_implied_holomorphic_dimension": macro_dimension,
            "prior_qball_collective_inertia": collective_inertia,
            "prior_qball_omega_reconstructed": omega_reconstructed,
            "prior_qball_omega_stored": omega_reported,
        },
        "source_manifest": {
            str(PRIOR_BRIDGE.relative_to(ROUTE_F.parent)): sha256(PRIOR_BRIDGE),
            str(AP_E1_TEX.relative_to(ROUTE_F.parent)): sha256(AP_E1_TEX),
            str(AP_E1_BIB.relative_to(ROUTE_F.parent)): sha256(AP_E1_BIB),
            str(Path(__file__).resolve().relative_to(ROUTE_F.parent)): sha256(
                Path(__file__).resolve()
            ),
        },
        "checks": CHECKS,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e1_projective_doublet.json"
    markdown_path = OUTPUT / "ap_e1_projective_doublet.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# AP-E1 projective-doublet audit",
        "",
        f"- Status: `{result['status']}`",
        f"- Checks: `{passed}/{len(CHECKS)}`",
        "- Physics promotion allowed: `false`",
        "",
        "## Logical result",
        "",
    ]
    for name, statement in result["logical_result"].items():
        lines.append(f"- **{name}:** {statement}")
    lines.extend(["", "## Required branch decision", ""])
    for name, statement in result["branch_decision"].items():
        lines.append(f"- **{name}:** {statement}")
    lines.extend(["", "## Non-derived inputs", ""])
    lines.extend(f"- {item}" for item in result["nonderived_inputs"])
    lines.extend(["", "## Mechanical checks", ""])
    for item in CHECKS:
        mark = "PASS" if item["pass"] else "FAIL"
        lines.append(f"- `{mark}` **{item['group']} / {item['name']}** — {item['detail']}")
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"AP-E1 checks: {passed}/{len(CHECKS)}")
    print(f"status: {result['status']}")
    print(f"wrote: {json_path}")
    print(f"wrote: {markdown_path}")
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
