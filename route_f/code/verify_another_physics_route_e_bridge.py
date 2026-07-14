#!/usr/bin/env python3
"""Verify algebraic/numerical facts for the Another-Physics/Route-E bridge.

This card deliberately separates two logical levels:

1. ``separation_theorem`` entries are consequences of an explicitly stated
   mathematical model (CP1/O(2), the fixed Route-E pairing, a chosen Q-ball
   action, or the displayed oscillatory equation).
2. ``bridge_axiom`` entries identify objects belonging to different models.
   They are research hypotheses, not consequences of either source theory.

A green run certifies reproducible arithmetic and symbolic identities only.
It never promotes the bridge to a physical theorem.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
ROUTE_E_INVARIANT_CARD = REPO / "route_E" / "output" / "audit0" / "invariant_card.json"

TOL = 1.0e-12
CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    """Record a fail-closed mechanical check."""
    CHECKS.append(
        {
            "group": group,
            "name": name,
            "pass": bool(condition),
            "detail": detail,
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matmul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a: list[list[complex]]) -> list[list[complex]]:
    return [list(row) for row in zip(*a)]


def matadd(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def max_matrix_abs(a: list[list[complex]]) -> float:
    return max(abs(cell) for row in a for cell in row)


def max_matrix_difference(a: list[list[complex]], b: list[list[complex]]) -> float:
    return max(
        abs(a[i][j] - b[i][j])
        for i in range(len(a))
        for j in range(len(a[0]))
    )


def quadratic_form(vector: list[complex], matrix: list[list[complex]]) -> complex:
    return sum(
        vector[i] * matrix[i][j] * vector[j]
        for i in range(len(vector))
        for j in range(len(vector))
    )


def invariant_residual(generator: list[list[complex]], pairing: list[list[complex]]) -> float:
    """Return max |rho(X)^T K + K rho(X)|."""
    return max_matrix_abs(
        matadd(matmul(transpose(generator), pairing), matmul(pairing, generator))
    )


def wrap_phase(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def oscillatory_function(x: float) -> float:
    return math.sin(x * x) + math.cos(x * x)


def qball_step_energy(
    radius: float,
    charge: float,
    f0_sq: float,
    potential_f0: float,
    surface_tension: float,
) -> float:
    volume = 4.0 * math.pi * radius**3 / 3.0
    area = 4.0 * math.pi * radius**2
    return (
        charge**2 / (4.0 * f0_sq * volume)
        + potential_f0 * volume
        + surface_tension * area
    )


def main() -> None:
    # ------------------------------------------------------------------ O(2)
    degree = 2
    basis = ["X^2", "sqrt(2) X Y", "Y^2"]
    # diag(t,t^{-1}) acts on X^(2-k)Y^k with weight 2-2k.
    weights = [degree - 2 * k for k in range(degree + 1)]
    h0 = degree + 1
    check(
        "S1_CP1_O2",
        "H0(CP1,O(2)) has dimension three",
        h0 == 3 and len(basis) == 3,
        "h^0(O(d))=d+1 for d>=0; d=2 gives 3",
    )
    check(
        "S1_CP1_O2",
        "Cartan weights are (+2,0,-2)",
        weights == [2, 0, -2],
        f"weights={weights} in basis {basis}",
    )

    # ------------------------------------------ quadratic sections / divisors
    # Homogeneous section s=a X^2+b X Y+c Y^2 has Delta=b^2-4ac.
    regular = {"a": 1.0, "b": 0.0, "c": -1.0}  # X^2-Y^2
    degenerate = {"a": 1.0, "b": 0.0, "c": 0.0}  # X^2
    delta_regular = regular["b"] ** 2 - 4.0 * regular["a"] * regular["c"]
    delta_degenerate = (
        degenerate["b"] ** 2 - 4.0 * degenerate["a"] * degenerate["c"]
    )
    roots_regular = [
        (-regular["b"] + sign * cmath.sqrt(delta_regular)) / (2.0 * regular["a"])
        for sign in (1.0, -1.0)
    ]
    root_degenerate = -degenerate["b"] / (2.0 * degenerate["a"])
    regular_root_residual = max(
        abs(regular["a"] * z * z + regular["b"] * z + regular["c"])
        for z in roots_regular
    )
    degenerate_root_residual = abs(
        degenerate["a"] * root_degenerate**2
        + degenerate["b"] * root_degenerate
        + degenerate["c"]
    )
    check(
        "S2_quadratic_divisor",
        "nonzero discriminant gives two distinct projective zeros",
        delta_regular != 0.0
        and abs(roots_regular[0] - roots_regular[1]) > 1.0
        and regular_root_residual < TOL,
        f"Delta={delta_regular:.1f}, roots={roots_regular}, residual={regular_root_residual:.2e}",
    )
    check(
        "S2_quadratic_divisor",
        "zero discriminant gives one double zero",
        delta_degenerate == 0.0 and degenerate_root_residual < TOL,
        f"Delta={delta_degenerate:.1f}, double root={root_degenerate:.1f}",
    )

    # ---------------------------------------------------- Route-E K_tr pairing
    if not ROUTE_E_INVARIANT_CARD.exists():
        raise FileNotFoundError(
            "Route-E convention source is required for an independent K_tr comparison: "
            f"{ROUTE_E_INVARIANT_CARD}"
        )
    route_e_card = json.loads(ROUTE_E_INVARIANT_CARD.read_text(encoding="utf-8"))
    raw_k = route_e_card["routeA_invariant_anchors"]["k_tr"]["matrix"]
    loaded_k = [
        [complex(float(cell["re"]), float(cell["im"])) for cell in row]
        for row in raw_k
    ]
    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    local_k = [
        [0j, 0j, -inv_sqrt3],
        [0j, inv_sqrt3, 0j],
        [-inv_sqrt3, 0j, 0j],
    ]
    k_source_difference = max_matrix_difference(loaded_k, local_k)
    k_squared = matmul(loaded_k, loaded_k)
    one_third_identity = [
        [complex(1.0 / 3.0 if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]
    k_square_residual = max_matrix_difference(k_squared, one_third_identity)

    rt2 = math.sqrt(2.0)
    rho_h = [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -2.0]]
    rho_e = [[0.0, rt2, 0.0], [0.0, 0.0, rt2], [0.0, 0.0, 0.0]]
    rho_f = transpose(rho_e)
    invariance = {
        "H": invariant_residual(rho_h, loaded_k),
        "E": invariant_residual(rho_e, loaded_k),
        "F": invariant_residual(rho_f, loaded_k),
    }
    check(
        "S3_Ktr",
        "locally reconstructed pairing matches the Route-E convention card",
        k_source_difference < 1.0e-15,
        f"max matrix difference={k_source_difference:.2e}",
    )
    check(
        "S3_Ktr",
        "K_tr squared equals I/3",
        k_square_residual < 1.0e-15,
        f"max|K_tr^2-I/3|={k_square_residual:.2e}",
    )
    check(
        "S3_Ktr",
        "K_tr is invariant under the full sl2 generator triple",
        max(invariance.values()) < 1.0e-14,
        "residuals=" + json.dumps(invariance, sort_keys=True),
    )

    # Let x=(x_+,x_0,x_-) be Route E's spherical coordinates.  For the
    # standard projective matrix A=[[-b/2,-c],[a,b/2]], the corresponding
    # binary quadratic has
    #   a=x_-/sqrt(2), b=-x_0, c=x_+/sqrt(2).
    # Thus B_Kill(x,x)=2 Delta(p)=2 sqrt(3) K_tr(x,x).  Keeping the factor two
    # explicit avoids an impermissible rescaling of Lie-algebra generators.
    x_regular = [complex(0.7), complex(-0.4), complex(0.2)]
    a_regular = x_regular[2] / math.sqrt(2.0)
    b_regular = -x_regular[1]
    c_regular = x_regular[0] / math.sqrt(2.0)
    polynomial_discriminant = b_regular**2 - 4.0 * a_regular * c_regular
    killing_matrix = [
        [2.0 * math.sqrt(3.0) * cell for cell in row] for row in loaded_k
    ]
    killing_norm = quadratic_form(x_regular, killing_matrix)
    discriminant_killing_error = abs(2.0 * polynomial_discriminant - killing_norm)

    projective_matrix = [
        [-b_regular / 2.0, -c_regular],
        [a_regular, b_regular / 2.0],
    ]
    projective_matrix_squared = matmul(projective_matrix, projective_matrix)
    delta_over_four_identity = [
        [polynomial_discriminant / 4.0 if i == j else 0j for j in range(2)]
        for i in range(2)
    ]
    projective_square_error = max_matrix_difference(
        projective_matrix_squared, delta_over_four_identity
    )
    canonical_killing_from_matrix = 4.0 * sum(
        projective_matrix_squared[i][i] for i in range(2)
    )
    canonical_factor_two_error = abs(
        canonical_killing_from_matrix - 2.0 * polynomial_discriminant
    )

    x_null = [complex(1.0), complex(1.0), complex(0.5)]
    a_null = x_null[2] / math.sqrt(2.0)
    b_null = -x_null[1]
    c_null = x_null[0] / math.sqrt(2.0)
    null_discriminant = b_null**2 - 4.0 * a_null * c_null
    null_killing_norm = quadratic_form(x_null, killing_matrix)
    check(
        "S3_Ktr",
        "Route-E Killing norm equals twice the binary-quadratic discriminant",
        discriminant_killing_error < 1.0e-14
        and abs(polynomial_discriminant) > 1.0e-2
        and abs(null_discriminant) < 1.0e-14
        and abs(null_killing_norm) < 1.0e-14,
        (
            f"regular Delta={polynomial_discriminant.real:.6f}, "
            f"B(x,x)={killing_norm.real:.6f}, |B-2Delta|={discriminant_killing_error:.2e}; "
            f"null Delta={abs(null_discriminant):.2e}, "
            f"null B={abs(null_killing_norm):.2e}"
        ),
    )
    check(
        "S3_Ktr",
        "standard projective matrix gives A^2=Delta I/4 and canonical B=2 Delta",
        projective_square_error < 1.0e-14 and canonical_factor_two_error < 1.0e-14,
        (
            f"max|A^2-Delta I/4|={projective_square_error:.2e}; "
            f"|4 Tr(A^2)-2Delta|={canonical_factor_two_error:.2e}"
        ),
    )

    # A bounded sign label is available from the CP1 moment map; this is a
    # charge/orientation label, not a negative Hamiltonian energy.
    def moment_map(z0: complex, z1: complex) -> float:
        norm_sq = abs(z0) ** 2 + abs(z1) ** 2
        return (abs(z0) ** 2 - abs(z1) ** 2) / norm_sq

    moment_map_values = {
        "p_plus": moment_map(1.0 + 0j, 0j),
        "equator": moment_map(1.0 + 0j, 1.0 + 0j),
        "p_minus": moment_map(0j, 1.0 + 0j),
    }
    check(
        "S3_Ktr",
        "the normalized CP1 moment map labels the two fixed points by plus/minus one",
        moment_map_values == {"p_plus": 1.0, "equator": 0.0, "p_minus": -1.0},
        "mu=" + json.dumps(moment_map_values, sort_keys=True),
    )

    # ------------------------------------------------------- zeta phase bridge
    zeta_source = route_e_card["majorana_invariant_anchors"]["zeta_tar"]
    zeta = complex(float(zeta_source["re"]), float(zeta_source["im"]))
    declared_zeta = complex(0.1076472949, 0.0736514853)
    zeta_source_error = abs(zeta - declared_zeta)
    arg_zeta = cmath.phase(zeta)
    visibility_arg = math.cos(arg_zeta) ** 2
    # If zeta is lifted as a weight-two composite of a weight-one phase delta,
    # then delta=arg(zeta)/2 modulo pi; this is not the same bridge convention.
    visibility_half_arg = math.cos(arg_zeta / 2.0) ** 2
    check(
        "S4_zeta",
        "benchmark magnitude and argument are reproduced",
        zeta_source_error < 1.0e-15
        and abs(abs(zeta) - 0.13043190325293763) < 1.0e-15
        and abs(arg_zeta - 0.600038020318215) < 1.0e-15,
        (
            f"canonical-card error={zeta_source_error:.2e}; "
            f"|zeta|={abs(zeta):.16f}, arg(zeta)={arg_zeta:.15f}"
        ),
    )
    check(
        "S4_zeta",
        "both explicit cos-squared phase conventions are reproducible",
        abs(visibility_arg - 0.6811434402919294) < 1.0e-15
        and abs(visibility_half_arg - 0.9126570732133188) < 1.0e-15,
        f"cos^2(arg zeta)={visibility_arg:.16f}; cos^2(arg zeta/2)={visibility_half_arg:.16f}",
    )

    scale_modulus = 1.7
    scale_phase = 0.37
    y = scale_modulus * cmath.exp(1j * scale_phase)
    weight_two_factor = y * y / abs(y) ** 2
    zeta_scaled = weight_two_factor * zeta
    zeta_expected = cmath.exp(2j * scale_phase) * zeta
    covariance_error = abs(zeta_scaled - zeta_expected)
    modulus_error = abs(abs(zeta_scaled) - abs(zeta))
    phase_error = abs(wrap_phase(cmath.phase(zeta_scaled) - arg_zeta - 2.0 * scale_phase))
    check(
        "S4_zeta",
        "complex rescaling obeys weight-two covariance",
        covariance_error < 1.0e-15 and modulus_error < 1.0e-15 and phase_error < 1.0e-14,
        f"|error|={covariance_error:.2e}, modulus error={modulus_error:.2e}, phase error={phase_error:.2e}",
    )
    check(
        "S4_zeta",
        "weight-two covariance is not complex-phase invariance",
        abs(zeta_scaled - zeta) > 1.0e-2,
        f"|zeta'-zeta|={abs(zeta_scaled-zeta):.6f} for arg(y)={scale_phase:.2f}",
    )

    # ----------------------------------------------------- Q-ball benchmark
    # U(f)=m^2 f^2-lambda f^4+(eta/M^2)f^6, Phi=f(r)e^{i omega t}.
    mass = 1.0
    cutoff = 1.0
    quartic = 1.0
    sextic = 1.0
    charge = 1.0e6
    f0_sq = quartic * cutoff**2 / (2.0 * sextic)
    f0 = math.sqrt(f0_sq)
    omega0_sq = mass**2 - quartic**2 * cutoff**2 / (4.0 * sextic)
    omega0 = math.sqrt(omega0_sq)
    potential_f0 = (
        mass**2 * f0_sq
        - quartic * f0_sq**2
        + sextic * f0_sq**3 / cutoff**2
    )
    factorization_residual = 0.0
    for i in range(101):
        field = f0 * i / 100.0
        potential = (
            mass**2 * field**2
            - quartic * field**4
            + sextic * field**6 / cutoff**2
        )
        factored = sextic / cutoff**2 * field**2 * (field**2 - f0_sq) ** 2
        factorization_residual = max(
            factorization_residual,
            abs((potential - omega0_sq * field**2) - factored),
        )
    check(
        "S5_Qball",
        "chosen sextic potential has a nonempty Q-ball frequency window",
        0.0 < omega0_sq < mass**2 and factorization_residual < 1.0e-14,
        f"omega0^2={omega0_sq:.6f}; window {omega0:.12f}<omega<{mass:.1f}; factorization residual={factorization_residual:.2e}",
    )

    surface_tension = quartic**2 * cutoff**3 / (8.0 * sextic ** 1.5)
    # Independent Simpson integration of sigma=2 int_0^f0 sqrt(U_omega0) df.
    n_simpson = 10000
    step = f0 / n_simpson
    surface_sum = 0.0
    for i in range(n_simpson + 1):
        field = i * step
        u_eff = sextic / cutoff**2 * field**2 * (field**2 - f0_sq) ** 2
        weight = 1.0 if i in (0, n_simpson) else (4.0 if i % 2 else 2.0)
        surface_sum += weight * math.sqrt(max(0.0, u_eff))
    surface_tension_numeric = 2.0 * step * surface_sum / 3.0
    check(
        "S5_Qball",
        "analytic wall tension agrees with independent quadrature",
        abs(surface_tension_numeric - surface_tension) < 1.0e-10,
        f"sigma_analytic={surface_tension:.12f}, sigma_numeric={surface_tension_numeric:.12f}",
    )

    radius_leading = (
        3.0 * charge / (8.0 * math.pi * omega0 * f0_sq)
    ) ** (1.0 / 3.0)
    radius_leading_fm = radius_leading * 0.1973269804
    eq_leading = omega0 + 4.0 * math.pi * radius_leading**2 * surface_tension / charge
    check(
        "S5_Qball",
        "leading thin-wall benchmark reproduces the quoted radius and E/Q",
        abs(radius_leading_fm - 12.842415685843626) < 1.0e-10
        and abs(eq_leading - 0.872678753981896) < 1.0e-12
        and eq_leading < mass,
        f"R={radius_leading:.9f} GeV^-1={radius_leading_fm:.9f} fm; E/Q={eq_leading:.12f} GeV",
    )

    # Finite-Q cross-check: minimize the same step-profile energy after
    # eliminating omega through Q=2 omega f0^2 V.
    def derivative_step_energy(radius: float) -> float:
        volume_coefficient = 4.0 * math.pi / 3.0
        a_coeff = charge**2 / (4.0 * f0_sq * volume_coefficient)
        b_coeff = potential_f0 * volume_coefficient
        c_coeff = 4.0 * math.pi * surface_tension
        return (
            -3.0 * a_coeff / radius**4
            + 3.0 * b_coeff * radius**2
            + 2.0 * c_coeff * radius
        )

    lo = 0.5 * radius_leading
    hi = 1.5 * radius_leading
    if not (derivative_step_energy(lo) < 0.0 < derivative_step_energy(hi)):
        raise RuntimeError("failed to bracket the finite-Q step-profile minimum")
    for _ in range(160):
        mid = 0.5 * (lo + hi)
        if derivative_step_energy(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    radius_variational = 0.5 * (lo + hi)
    volume_variational = 4.0 * math.pi * radius_variational**3 / 3.0
    omega_variational = charge / (2.0 * f0_sq * volume_variational)
    eq_variational = qball_step_energy(
        radius_variational, charge, f0_sq, potential_f0, surface_tension
    ) / charge
    check(
        "S5_Qball",
        "finite-Q step-profile variational cross-check stays inside the window and below free-particle threshold",
        omega0 < omega_variational < mass
        and eq_variational < mass
        and abs(eq_variational - eq_leading) / eq_variational < 2.0e-5,
        f"R_var={radius_variational:.9f} GeV^-1, omega_var={omega_variational:.12f}, E/Q={eq_variational:.12f} GeV",
    )

    # --------------------------------------- sin(x^2)+cos(x^2) solution sets
    # sin A=sin B gives A=B+2pi n or A=pi-B+2pi n after shifting by pi/4.
    circle_residual = 0.0
    for n in (0, 1, 2, 10):
        radius = math.sqrt(math.pi / 2.0 + 2.0 * math.pi * n)
        for j in range(65):
            theta = 2.0 * math.pi * j / 64.0
            x = radius * math.cos(theta)
            y_coord = radius * math.sin(theta)
            circle_residual = max(
                circle_residual, abs(oscillatory_function(x) - oscillatory_function(y_coord))
            )
    hyperbola_residual = 0.0
    for n in (1, 3, 10):
        for j in range(-40, 41):
            y_coord = 3.0 * j / 40.0
            x = math.sqrt(y_coord**2 + 2.0 * math.pi * n)
            hyperbola_residual = max(
                hyperbola_residual, abs(oscillatory_function(x) - oscillatory_function(y_coord))
            )
    check(
        "S6_quadratic_phase",
        "circle family x^2+y^2=pi/2+2pi n solves f(x)=f(y)",
        circle_residual < 1.0e-12,
        f"maximum sampled residual={circle_residual:.2e}",
    )
    check(
        "S6_quadratic_phase",
        "hyperbola family x^2-y^2=2pi n solves f(x)=f(y)",
        hyperbola_residual < 1.0e-12,
        f"maximum sampled residual={hyperbola_residual:.2e}",
    )

    pitch_n = 1_000_000
    radius_n = math.sqrt(math.pi / 2.0 + 2.0 * math.pi * pitch_n)
    radius_next = math.sqrt(math.pi / 2.0 + 2.0 * math.pi * (pitch_n + 1))
    pitch = radius_next - radius_n
    pitch_exact_stable = 2.0 * math.pi / (radius_next + radius_n)
    pitch_asymptotic = math.pi / radius_n
    exact_pitch_error = abs(pitch - pitch_exact_stable)
    asymptotic_relative_error = abs(pitch / pitch_asymptotic - 1.0)
    check(
        "S6_quadratic_phase",
        "radial pitch obeys the exact rationalized identity and pi/r asymptotic",
        exact_pitch_error < 1.0e-12 and asymptotic_relative_error < 1.0e-6,
        f"n={pitch_n}, Delta r={pitch:.12e}, identity error={exact_pitch_error:.2e}, asymptotic relative error={asymptotic_relative_error:.2e}",
    )

    # ----------------------------------------------------- logical boundary
    bridge_axioms = [
        {
            "id": "BA1_phase_identification",
            "statement": (
                "Identify an Another-Physics visibility phase delta with either arg(zeta) "
                "or the weight-one lift arg(zeta)/2.  These choices predict different "
                "cos^2 values and are not selected by Route-E."
            ),
            "status": "not_derived",
        },
        {
            "id": "BA2_bubble_identification",
            "statement": (
                "Identify the proposed energy bubble with the charge-Q soliton of the "
                "displayed complex-scalar action.  Q-ball existence does not derive this map."
            ),
            "status": "not_derived",
        },
        {
            "id": "BA3_two_center_identification",
            "statement": (
                "Identify a quadratic-phase circle/hyperbola morphology with the two-zero "
                "divisor of a CP1/O(2) section.  Shared degree-two algebra is insufficient "
                "to identify their state spaces or dynamics."
            ),
            "status": "not_derived",
        },
        {
            "id": "BA4_UV_embedding",
            "statement": (
                "Embed the Q-ball scalar and its conserved U(1) charge in the Route-E "
                "Spin(10)/family sector while preserving gauge, anomaly, and threshold gates."
            ),
            "status": "open",
        },
    ]
    check(
        "S7_boundary",
        "every cross-theory identification remains an explicit bridge axiom",
        all(item["status"] != "derived" for item in bridge_axioms),
        f"{len(bridge_axioms)} bridge axioms retained as non-derived",
    )

    passed = sum(item["pass"] for item in CHECKS)
    all_pass = passed == len(CHECKS)
    result: dict[str, Any] = {
        "audit": "verify_another_physics_route_e_bridge",
        "status": (
            "mechanical_separation_checks_pass_no_physics_promotion"
            if all_pass
            else "mechanical_check_failure"
        ),
        "all_pass": all_pass,
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "physics_promotion_allowed": False,
        "logical_contract": {
            "separation_theorem": (
                "A result derived inside one explicitly stated mathematical model; it may "
                "rule out a proposed identification but cannot establish a cross-model map."
            ),
            "bridge_axiom": (
                "Additional cross-theory identification requiring action-level derivation "
                "and independent phenomenological tests."
            ),
            "green_run_meaning": "arithmetic/symbolic reproducibility only",
        },
        "source_manifest": {
            "route_E_invariant_card": {
                "path": str(ROUTE_E_INVARIANT_CARD.relative_to(REPO)),
                "sha256": sha256(ROUTE_E_INVARIANT_CARD),
            }
        },
        "separation_theorems": {
            "CP1_O2": {
                "basis": basis,
                "dimension": h0,
                "cartan_weights": weights,
                "scope": "conditional on the CP1/O(2) carrier",
            },
            "quadratic_divisor": {
                "formula": "Delta=b^2-4ac for a X^2+b X Y+c Y^2",
                "regular_example": {
                    "coefficients": regular,
                    "discriminant": delta_regular,
                    "affine_roots": [[z.real, z.imag] for z in roots_regular],
                    "two_distinct_zeros": True,
                },
                "degenerate_example": {
                    "coefficients": degenerate,
                    "discriminant": delta_degenerate,
                    "double_root": root_degenerate,
                },
            },
            "K_tr": {
                "matrix": [[cell.real for cell in row] for row in loaded_k],
                "source_comparison_max_abs": k_source_difference,
                "K_squared_minus_I_over_3_max_abs": k_square_residual,
                "sl2_invariance_residuals": invariance,
                "invariance_equation": "rho(A)^T K_tr+K_tr rho(A)=0 for A=H,E,F",
                "discriminant_contact_identity": {
                    "route_E_spherical_coordinates": "x=(x_+,x_0,x_-)",
                    "section_map": "p=(x_-/sqrt(2)) X^2-x_0 X Y+(x_+/sqrt(2)) Y^2",
                    "formula": "B_Kill(x,x)=2 Delta(p)=2 sqrt(3) K_tr(x,x)",
                    "regular_example_error": discriminant_killing_error,
                    "projective_matrix_square_error": projective_square_error,
                    "canonical_factor_two_error": canonical_factor_two_error,
                    "contact_null_example_discriminant_abs": abs(null_discriminant),
                    "contact_null_example_killing_norm_abs": abs(null_killing_norm),
                },
                "CP1_moment_map": {
                    "formula": "mu=([z0]^2-[z1]^2)/([z0]^2+[z1]^2), with brackets denoting modulus",
                    "fixed_point_values": moment_map_values,
                    "normalization_note": "the displayed generator and symmetric additive convention set the endpoints to +/-1; generator rescaling changes their magnitude",
                    "interpretation": "bounded charge/orientation polarity, not negative energy",
                },
            },
            "zeta": {
                "value": [zeta.real, zeta.imag],
                "loaded_from_route_E_invariant_card": True,
                "declared_benchmark_source_error": zeta_source_error,
                "abs": abs(zeta),
                "arg_rad": arg_zeta,
                "cos_squared_arg": visibility_arg,
                "cos_squared_half_arg_weight_one_lift": visibility_half_arg,
                "complex_rescaling": {
                    "y": [y.real, y.imag],
                    "arg_y": scale_phase,
                    "law": "zeta'=(y^2/|y|^2)zeta=exp(2 i arg(y))zeta",
                    "zeta_prime": [zeta_scaled.real, zeta_scaled.imag],
                    "covariance_error": covariance_error,
                    "modulus_error": modulus_error,
                    "phase_error_rad": phase_error,
                    "phase_invariant": False,
                },
            },
            "Q_ball": {
                "action_scope": "flat-space classical complex scalar; thin-wall/step-profile approximation",
                "potential": "U(f)=m^2 f^2-lambda f^4+(eta/M^2)f^6",
                "parameters_GeV_convention": {
                    "m": mass,
                    "M": cutoff,
                    "lambda": quartic,
                    "eta": sextic,
                    "Q": charge,
                },
                "f0_squared_GeV2": f0_sq,
                "omega0_squared_GeV2": omega0_sq,
                "omega0_GeV": omega0,
                "existence_window_GeV": [omega0, mass],
                "surface_tension_GeV3": surface_tension,
                "leading_thin_wall": {
                    "radius_GeV_inverse": radius_leading,
                    "radius_fm": radius_leading_fm,
                    "E_over_Q_GeV": eq_leading,
                },
                "finite_Q_step_variational": {
                    "radius_GeV_inverse": radius_variational,
                    "omega_GeV": omega_variational,
                    "E_over_Q_GeV": eq_variational,
                },
                "below_free_particle_threshold": eq_variational < mass,
                "not_an_exact_radial_ODE_solution": True,
            },
            "quadratic_phase_solution": {
                "equation": "sin(x^2)+cos(x^2)=sin(y^2)+cos(y^2)",
                "families": [
                    "x^2-y^2=2 pi n",
                    "x^2+y^2=pi/2+2 pi n",
                ],
                "sampled_max_residual": {
                    "hyperbola": hyperbola_residual,
                    "circle": circle_residual,
                },
                "pitch": {
                    "exact": "r_(n+1)-r_n=2 pi/(r_(n+1)+r_n)",
                    "asymptotic": "Delta r_n~pi/r_n",
                    "test_n": pitch_n,
                    "relative_error": asymptotic_relative_error,
                },
                "scope": "continuous interference level sets, not discrete atomic orbitals",
            },
        },
        "bridge_axioms": bridge_axioms,
        "checks": CHECKS,
        "open_physics_gates": [
            "derive rather than identify the phase map between zeta and a visibility observable",
            "supply a gauge/anomaly-consistent Route-E embedding of the Q-ball U(1)",
            "solve the radial Q-ball boundary-value problem and fluctuation spectrum",
            "derive a dynamical map, if any, between quadratic-phase level sets and CP1 divisors",
            "rerun Route-E flavor, threshold, proton, and cosmology gates after any embedding",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "another_physics_route_e_bridge.json"
    md_path = OUTPUT / "another_physics_route_e_bridge.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    check_lines = "\n".join(
        f"- [{'PASS' if item['pass'] else 'FAIL'}] `{item['group']}` — "
        f"{item['name']}: {item['detail']}"
        for item in CHECKS
    )
    bridge_lines = "\n".join(
        f"- **{item['id']}** (`{item['status']}`): {item['statement']}"
        for item in bridge_axioms
    )
    md_path.write_text(
        f"""# Another Physics / Route-E Bridge Verification

Status: **{result['status']}** — `{passed}/{len(CHECKS)}` checks pass.

This is a fail-closed algebraic/numerical card.  A green result means only
that the displayed identities are reproducible.  It does **not** identify the
two theories and `physics_promotion_allowed=false`.

## Separation theorems

1. On the assumed carrier, `H^0(CP1,O(2))` has basis
   `(X^2,sqrt(2)XY,Y^2)`, dimension `3`, and Cartan weights `(+2,0,-2)`.
2. A quadratic section `aX^2+bXY+cY^2` has discriminant
   `Delta=b^2-4ac`: `Delta != 0` gives two distinct projective zeros, while
   `Delta=0` gives a double zero.  This separates the regular and degenerate
   cases; it does not make every degree-two pattern the same physical object.
3. The pairing loaded from Route-E obeys
   `rho(A)^T K_tr+K_tr rho(A)=0` for `A=H,E,F` and `K_tr^2=I/3`.
   In Route-E spherical coordinates it also obeys the exact identity
   `B_Kill(x,x)=2 Delta(p)=2 sqrt(3) K_tr(x,x)`, while the standard
   projective matrix satisfies `A^2=Delta I/4`.  Thus the contact-null cone is
   the double-zero locus.  With the displayed generator normalization the CP1
   moment map assigns `+1/-1` to the two fixed points; their opposition, not
   that numerical magnitude, is invariant.  These are charge/orientation
   labels, not energies.
4. For `zeta={zeta.real:.10f}+{zeta.imag:.10f}i`,
   `arg(zeta)={arg_zeta:.15f}`,
   `cos^2(arg zeta)={visibility_arg:.15f}`, and the weight-one lift gives
   `cos^2(arg zeta/2)={visibility_half_arg:.15f}`.  Under complex `y`,
   `zeta'=(y^2/|y|^2)zeta`; this is weight-two covariance, not phase invariance.
5. For the explicitly chosen sextic Q-ball action, the existence window is
   `{omega0:.12f}<omega<1 GeV`.  At `Q=10^6`, the leading thin-wall card gives
   `R={radius_leading_fm:.9f} fm` and `E/Q={eq_leading:.12f} GeV`; the independent
   finite-Q step-profile minimization gives `E/Q={eq_variational:.12f} GeV < m`.
   Both are controlled approximations, not an exact radial-ODE solution.
6. The displayed oscillatory equation has continuous solution families
   `x^2-y^2=2 pi n` and `x^2+y^2=pi/2+2 pi n`.  For the circular family,
   `Delta r=2 pi/(r_(n+1)+r_n)~pi/r_n`; this proves a shrinking interference
   pitch, not atomic-orbital quantization.

## Bridge axioms retained as non-derived

{bridge_lines}

## Mechanical checks

{check_lines}

## Promotion boundary

- `physics_promotion_allowed=false`.
- A physical bridge requires an action-level phase map, a gauge/anomaly-safe
  embedding, an exact soliton/fluctuation analysis, and reruns of Route-E's
  flavor, threshold, proton, and cosmology gates.
""",
        encoding="utf-8",
    )

    print(
        f"another_physics_route_e_bridge: {passed}/{len(CHECKS)} checks pass; "
        "physics_promotion_allowed=false"
    )
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
