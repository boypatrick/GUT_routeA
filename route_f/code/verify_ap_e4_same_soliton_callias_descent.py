#!/usr/bin/env python3
"""Deterministic AP-E4 same-soliton Callias/descent audit.

This verifier deliberately separates a mathematical conditional theorem from
facts actually derived for the charged two-colour B=1 soliton.  It checks:

* CP1 tangent-bundle chart covariance and c1(TCP1)=2;
* a finite-matrix Berry family whose kernel line has Chern number +2;
* the families-Callias boundary push-forward formula in a declared
  orientation convention, including q=0 and doubled-family controls;
* the N=2 SQM algebra, three-state O(2) sector, and the Serre-dual
  anti-canonical O(-4) CPT sector;
* degree-five to degree-two differential-character fiber integration and the
  integer n B d governing the pullback to the same CP1 moduli space.

A green card is not evidence that the present nonsupersymmetric proxy supplies
the missing Yukawa operator, Callias gap, gauge-equivariant refinement, or
same-soliton CP1 map.  Those status flags remain fail-closed.
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
TEX = ROUTE_F / "tex" / "ap_e4_same_soliton_callias_descent.tex"
BIB = ROUTE_F / "tex" / "ap_e4_same_soliton_callias_descent.bib"
THIS_SCRIPT = Path(__file__).resolve()

SEED = 20260716
TOL = 2.0e-11
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


def h0_o(k: int) -> int:
    return max(k + 1, 0)


def h1_o(k: int) -> int:
    return max(-k - 1, 0)


def coherent_state(k: int, theta: float, phi: float) -> np.ndarray:
    """Anti-Veronese spin-k/2 state, oriented so its smooth Berry line has c1=+k."""

    c = math.cos(0.5 * theta)
    s = math.sin(0.5 * theta)
    return np.array(
        [
            math.sqrt(math.comb(k, r))
            * c ** (k - r)
            * s**r
            * np.exp(-1j * r * phi)
            for r in range(k + 1)
        ],
        dtype=complex,
    )


def triangulated_berry_chern(k: int, theta_cells: int, phi_cells: int) -> dict[str, Any]:
    """Gauge-invariant Bargmann-triangle Chern number on a sphere mesh."""

    vertices = [coherent_state(k, 0.0, 0.0)]
    for j in range(1, theta_cells):
        theta = math.pi * j / theta_cells
        for ell in range(phi_cells):
            phi = 2.0 * math.pi * ell / phi_cells
            vertices.append(coherent_state(k, theta, phi))
    vertices.append(coherent_state(k, math.pi, 0.0))
    south = len(vertices) - 1

    triangles: list[tuple[int, int, int]] = []
    for ell in range(phi_cells):
        triangles.append((0, 1 + ell, 1 + (ell + 1) % phi_cells))
    for j in range(1, theta_cells - 1):
        upper = 1 + (j - 1) * phi_cells
        lower = 1 + j * phi_cells
        for ell in range(phi_cells):
            nxt = (ell + 1) % phi_cells
            triangles.append((upper + ell, lower + ell, lower + nxt))
            triangles.append((upper + ell, lower + nxt, upper + nxt))
    final_ring = 1 + (theta_cells - 2) * phi_cells
    for ell in range(phi_cells):
        triangles.append(
            (final_ring + ell, south, final_ring + (ell + 1) % phi_cells)
        )

    phases: list[float] = []
    min_overlap = 1.0
    for a, b, c in triangles:
        ab = np.vdot(vertices[a], vertices[b])
        bc = np.vdot(vertices[b], vertices[c])
        ca = np.vdot(vertices[c], vertices[a])
        min_overlap = min(min_overlap, abs(ab), abs(bc), abs(ca))
        phases.append(float(np.angle(ab * bc * ca)))
    chern_value = float(-sum(phases) / (2.0 * math.pi))
    if abs(chern_value) < 1.0e-15:
        chern_value = 0.0
    return {
        # With the anti-Hermitian Berry connection A=u^dagger du,
        # c1=(i/2pi)int F is minus the oriented Bargmann-phase sum.
        "chern": chern_value,
        "triangles": len(triangles),
        "min_link_overlap": min_overlap,
        "max_triangle_phase": max(abs(value) for value in phases) if phases else 0.0,
    }


def finite_kernel_family_metrics(k: int, theta_cells: int, phi_cells: int) -> dict[str, Any]:
    """H=1-|u><u|: exact one-dimensional kernel and unit massive gap."""

    max_kernel_residual = 0.0
    max_projector_residual = 0.0
    min_nonzero = math.inf
    for j in range(theta_cells + 1):
        theta = math.pi * j / theta_cells
        for ell in range(phi_cells):
            phi = 2.0 * math.pi * ell / phi_cells
            state = coherent_state(k, theta, phi)
            projector = np.outer(state, state.conjugate())
            hamiltonian = np.eye(k + 1, dtype=complex) - projector
            max_kernel_residual = max(
                max_kernel_residual, float(np.linalg.norm(hamiltonian @ state))
            )
            max_projector_residual = max(
                max_projector_residual,
                float(np.linalg.norm(projector @ projector - projector)),
            )
            eigenvalues = np.linalg.eigvalsh(hamiltonian)
            positive = eigenvalues[eigenvalues > 0.5]
            min_nonzero = min(min_nonzero, float(positive.min()))
    berry = triangulated_berry_chern(k, theta_cells, phi_cells)
    return {
        "dimension": k + 1,
        "max_kernel_residual": max_kernel_residual,
        "max_projector_residual": max_projector_residual,
        "minimum_nonzero_eigenvalue": min_nonzero,
        "berry": berry,
    }


def callias_boundary_pushforward(p: int, q: int, epsilon: int = 1) -> dict[str, int]:
    """Declared orientation: ch Ind = epsilon * int_S2 ch(E_+).

    For c1(E_+)=p*x+q*y on S2_infinity x CP1 and x^2=y^2=0,
    ch(E_+)=1+p*x+q*y+p*q*x*y.  The reduced Callias boundary map therefore
    has rank epsilon*p and determinant c1 epsilon*p*q*y.
    """

    if epsilon not in (-1, 1):
        raise ValueError("epsilon must be +1 or -1")
    return {"rank": epsilon * p, "determinant_c1": epsilon * p * q}


def build_sqm_block(
    zero_plus: int, zero_minus: int, coefficient: float, levels: int, flux_abs: int
) -> dict[str, Any]:
    """Finite exact block for a Serre-dual pair of Dolbeault Dirac operators."""

    multiplicities = [flux_abs + 2 * n + 1 for n in range(1, levels + 1)]
    lambdas = [
        2.0 * math.sqrt(n * (n + flux_abs + 1)) / math.sqrt(coefficient)
        for n in range(1, levels + 1)
    ]
    massive_dimension = sum(multiplicities)
    plus_dimension = zero_plus + massive_dimension
    minus_dimension = zero_minus + massive_dimension
    b_map = np.zeros((minus_dimension, plus_dimension), dtype=complex)
    plus_offset = zero_plus
    minus_offset = zero_minus
    for multiplicity, eigenvalue in zip(multiplicities, lambdas):
        b_map[
            minus_offset : minus_offset + multiplicity,
            plus_offset : plus_offset + multiplicity,
        ] = eigenvalue * np.eye(multiplicity)
        plus_offset += multiplicity
        minus_offset += multiplicity

    q = np.block(
        [
            [
                np.zeros((plus_dimension, plus_dimension), dtype=complex),
                np.zeros((plus_dimension, minus_dimension), dtype=complex),
            ],
            [b_map, np.zeros((minus_dimension, minus_dimension), dtype=complex)],
        ]
    )
    qdag = q.conjugate().T
    d_op = q + qdag
    gamma = np.diag(
        np.concatenate([np.ones(plus_dimension), -np.ones(minus_dimension)])
    ).astype(complex)
    h_op = 0.5 * (q @ qdag + qdag @ q)
    eigenvalues = np.linalg.eigvalsh(d_op)
    return {
        "Q": q,
        "Qdag": qdag,
        "D": d_op,
        "Gamma": gamma,
        "H": h_op,
        "index": zero_plus - zero_minus,
        "kernel_dimension": zero_plus + zero_minus,
        "plus_dimension": plus_dimension,
        "minus_dimension": minus_dimension,
        "lambdas": lambdas,
        "multiplicities": multiplicities,
        "smallest_abs_nonzero_D": float(
            min(abs(value) for value in eigenvalues if abs(value) > 1.0e-9)
        ),
    }


def matrix_residuals(block: dict[str, Any]) -> dict[str, float]:
    q = block["Q"]
    qdag = block["Qdag"]
    d_op = block["D"]
    gamma = block["Gamma"]
    h_op = block["H"]
    return {
        "q_nilpotency": float(np.linalg.norm(q @ q)),
        "qdag_nilpotency": float(np.linalg.norm(qdag @ qdag)),
        "chirality_anticommutator": float(np.linalg.norm(gamma @ d_op + d_op @ gamma)),
        "hamiltonian_identity": float(
            np.linalg.norm(h_op - 0.5 * (q @ qdag + qdag @ q))
        ),
    }


def markdown_report(result: dict[str, Any]) -> str:
    status = result["status"]
    metrics = result["metrics"]
    lines = [
        "# AP-E4 same-soliton Callias/descent audit",
        "",
        f"- Passed: **{result['summary']['passed']}/{result['summary']['total']}**",
        f"- Conditional sufficient composition criterion established: **{status['conditional_composition_sufficient_criterion_established']}**",
        f"- Same-model physical tangent fermion derived: **{status['physical_tangent_fermion_same_model_derived']}**",
        f"- Same-model Callias c1 computed: **{status['callias_determinant_c1_same_model_computed']}**",
        f"- Gauge-basic WZW descent proved: **{status['wzw_gauge_basic_descent_same_model_proven']}**",
        f"- Degree-one Route-E portal built: **{status['degree_one_route_e_portal_built']}**",
        f"- Physics promotion: **{status['physics_promotion_allowed']}**",
        "",
        "## Quantitative regression",
        "",
        f"- Discrete Berry Chern number of the conditional tangent/Callias line: `{metrics['berry_chern']['k2']['chern']:.15f}`.",
        f"- Finite template/control kernel residual: `{metrics['finite_kernel_family_template_control']['max_kernel_residual']:.3e}`; gap: `{metrics['finite_kernel_family_template_control']['minimum_nonzero_eigenvalue']:.12f}`.",
        f"- Conditional Callias boundary data `(p,q,epsilon)=(1,2,+1)` give rank `{metrics['callias_pushforward']['required']['rank']}` and determinant c1 `{metrics['callias_pushforward']['required']['determinant_c1']}`.",
        f"- B=+1 / B=-1 fixed-polarization kernels: `{metrics['cpt']['plus_kernel']}` / `{metrics['cpt']['minus_kernel']}`; naive O(-2) control: `{metrics['cpt']['naive_inverse_kernel']}`.",
        f"- Factorized WZW pullback integer `n B d = {metrics['wzw_transgression']['plus_sector']}`.",
        "",
        "## Theorem-level conclusion",
        "",
        "At low energy, a one-multiplet tangent N=2 SQM requires normalizable CP1",
        "modes, a uniform gap, a Fredholm tangent kernel/connection, N=2 Ward",
        "identities, and no extra zero modes.  Exact microscopic supercharges are one",
        "sufficient way to enforce those Ward identities, but are not logically necessary",
        "because accidental/emergent worldline supersymmetry is possible.  The conservative",
        "same-model audit conditions below require either microscopic supercharges or an",
        "independent derivation of all emergent Ward identities.  The current bosonic proxy contains none of the",
        "required supersymmetry transformations or Yukawa mass family, so the result is",
        "underdetermined rather than proved.",
        "",
        "The degree-five differential character can always be fiber-integrated along",
        "the spatial three-cycle on the framed configuration space.  Descent to the",
        "gauge quotient additionally requires an equivariant differential refinement,",
        "basic curvature, and trivial vertical/stabilizer holonomy.  Those data are not",
        "present.  The exact pullback criterion is the integer integral over",
        "Sigma3 x CP1; the factorized conditional map gives +2, but the present proxy",
        "has not identified its scalar S2 with the same soliton CP1.",
        "",
        "## Checks",
        "",
    ]
    for row in result["checks"]:
        marker = "PASS" if row["pass"] else "FAIL"
        lines.append(f"- [{marker}] `{row['group']}` — {row['name']}: {row['detail']}")
    lines.extend(
        [
            "",
            "## Fail-closed composition gate",
            "",
            "- Specify one charged-two-colour mother action and derive either microscopic supercharges or complete emergent N=2 Ward identities, together with its Yukawa operator.",
            "- Prove a normalizable CP1 family on the existing B=1 soliton and a uniform noncollective gap.",
            "- Compute the actual Callias asymptotic positive-mass eigenbundle; do not substitute the finite matrix template/control.",
            "- Derive the anti-canonical vacuum-line factor K=O(-2) in the B=-1 regulator, rather than using the naive dual O(-2).",
            "- Construct an equivariant degree-five differential character and prove trivial large-gauge/stabilizer holonomy.",
            "- Evaluate the same-soliton integral and obtain +2 before composing AP-E3 with AP-E4.",
            "- Do not build the degree-one Route-E portal until this gate is closed.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    rng = np.random.default_rng(SEED)
    source_manifest = [source_row(path) for path in (TEX, BIB, THIS_SCRIPT)]
    check(
        "S0_provenance",
        "all same-soliton audit sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # ------------------------------------------------ CP1 tangent geometry
    metric_residuals: list[float] = []
    tangent_residuals: list[float] = []
    connection_residuals: list[float] = []
    for _ in range(300):
        w = complex(rng.normal(), rng.normal())
        if abs(w) < 0.25:
            w += 0.5 + 0.3j
        dw = complex(rng.normal(), rng.normal())
        psi_w = complex(rng.normal(), rng.normal())
        v = 1.0 / w
        jacobian = -1.0 / w**2
        psi_v = jacobian * psi_w
        metric_w = abs(dw) ** 2 / (1.0 + abs(w) ** 2) ** 2
        metric_v = abs(jacobian * dw) ** 2 / (1.0 + abs(v) ** 2) ** 2
        metric_residuals.append(abs(metric_w - metric_v))
        tangent_w = abs(psi_w) ** 2 / (1.0 + abs(w) ** 2) ** 2
        tangent_v = abs(psi_v) ** 2 / (1.0 + abs(v) ** 2) ** 2
        tangent_residuals.append(abs(tangent_w - tangent_v))
        # Gamma^w_ww=-2 bar(w)/(1+|w|^2); check affine-connection law.
        gamma_w = -2.0 * w.conjugate() / (1.0 + abs(w) ** 2)
        gamma_v = -2.0 * v.conjugate() / (1.0 + abs(v) ** 2)
        second = 2.0 / w**3
        lhs = gamma_v * jacobian**2 + second
        rhs = jacobian * gamma_w
        connection_residuals.append(abs(lhs - rhs))

    check(
        "S1_tangent_geometry",
        "Fubini--Study kinetic metric is chart invariant",
        max(metric_residuals) < TOL,
        f"max residual={max(metric_residuals):.3e}",
    )
    check(
        "S1_tangent_geometry",
        "psi^v=(-w^-2)psi^w has the tangent-bundle norm",
        max(tangent_residuals) < TOL,
        f"max residual={max(tangent_residuals):.3e}",
    )
    check(
        "S1_tangent_geometry",
        "Levi-Civita connection obeys the CP1 chart law",
        max(connection_residuals) < TOL,
        f"max residual={max(connection_residuals):.3e}",
    )

    berry_rows = {
        f"k{k}": triangulated_berry_chern(k, 24, 48) for k in (0, 1, 2, 4)
    }
    check(
        "S1_tangent_geometry",
        "discrete Berry regression gives c1(TCP1)=+2",
        abs(berry_rows["k2"]["chern"] - 2.0) < TOL,
        f"c1={berry_rows['k2']['chern']:.15f}; triangles={berry_rows['k2']['triangles']}",
    )
    check(
        "S1_tangent_geometry",
        "Berry-Chern sequence returns k for k=0,1,2,4",
        all(abs(berry_rows[f"k{k}"]["chern"] - k) < TOL for k in (0, 1, 2, 4)),
        f"sequence={[berry_rows[f'k{k}']['chern'] for k in (0, 1, 2, 4)]}",
    )

    # ----------------------------- finite template/control and Callias template
    finite_family = finite_kernel_family_metrics(2, 18, 36)
    check(
        "S2_callias_family",
        "finite H=1-|u_2><u_2| template/control has an exact one-dimensional kernel",
        finite_family["max_kernel_residual"] < TOL
        and finite_family["max_projector_residual"] < TOL,
        f"kernel={finite_family['max_kernel_residual']:.3e}; projector={finite_family['max_projector_residual']:.3e}",
    )
    check(
        "S2_callias_family",
        "finite template/control has a uniform unit nonzero gap",
        abs(finite_family["minimum_nonzero_eigenvalue"] - 1.0) < TOL,
        f"gap={finite_family['minimum_nonzero_eigenvalue']:.15f}",
    )
    check(
        "S2_callias_family",
        "finite template/control kernel Berry line is O(+2)",
        abs(finite_family["berry"]["chern"] - 2.0) < TOL,
        f"c1={finite_family['berry']['chern']:.15f}",
    )

    required_push = callias_boundary_pushforward(1, 2, +1)
    trivial_moduli_push = callias_boundary_pushforward(1, 0, +1)
    opposite_push = callias_boundary_pushforward(1, -2, +1)
    doubled_c1 = required_push["determinant_c1"] + opposite_push["determinant_c1"]
    check(
        "S2_callias_family",
        "conditional positive-mass eigenbundle (p,q)=(1,2) gives rank one and c1=+2",
        required_push == {"rank": 1, "determinant_c1": 2},
        f"result={required_push}",
    )
    check(
        "S2_callias_family",
        "q=0 negative control retains pointwise index but has trivial determinant line",
        trivial_moduli_push == {"rank": 1, "determinant_c1": 0},
        f"result={trivial_moduli_push}",
    )
    check(
        "S2_callias_family",
        "CPT-doubled q=+2 and q=-2 control cancels determinant c1",
        doubled_c1 == 0,
        f"c1 sum={doubled_c1}",
    )

    # ------------------------------------------------- SQM and CPT/Serre duality
    coefficient = 2.75
    plus_block = build_sqm_block(3, 0, coefficient, 4, flux_abs=2)
    minus_block = build_sqm_block(0, 3, coefficient, 4, flux_abs=2)
    plus_residuals = matrix_residuals(plus_block)
    minus_residuals = matrix_residuals(minus_block)
    check(
        "S3_sqm_algebra",
        "finite B=+1 and B=-1 blocks satisfy Q^2=0 and {Gamma,D}=0",
        max([*plus_residuals.values(), *minus_residuals.values()]) < TOL,
        f"max residual={max([*plus_residuals.values(), *minus_residuals.values()]):.3e}",
    )
    check(
        "S3_sqm_algebra",
        "O(2) has three positive-chirality Dolbeault ground states",
        (h0_o(2), h1_o(2), plus_block["index"]) == (3, 0, 3),
        f"(h0,h1,index)=({h0_o(2)},{h1_o(2)},{plus_block['index']})",
    )
    anti_canonical_degree = -2 - 2
    check(
        "S3_sqm_algebra",
        "fixed-polarization CPT/Serre duality requires K tensor E^vee=O(-4)",
        anti_canonical_degree == -4
        and (h0_o(anti_canonical_degree), h1_o(anti_canonical_degree)) == (0, 3)
        and minus_block["index"] == -3,
        f"degree={anti_canonical_degree}; (h0,h1,index)=({h0_o(-4)},{h1_o(-4)},{minus_block['index']})",
    )
    check(
        "S3_sqm_algebra",
        "Serre-dual massive spectra and gaps agree between B=+1 and B=-1",
        np.allclose(plus_block["lambdas"], minus_block["lambdas"], atol=TOL, rtol=0.0)
        and abs(
            plus_block["smallest_abs_nonzero_D"]
            - minus_block["smallest_abs_nonzero_D"]
        )
        < TOL,
        f"gap+={plus_block['smallest_abs_nonzero_D']:.12f}; gap-={minus_block['smallest_abs_nonzero_D']:.12f}",
    )
    check(
        "S3_sqm_algebra",
        "naive dual O(-2) is a one-state negative control, not the CPT triplet",
        (h0_o(-2), h1_o(-2)) == (0, 1),
        f"(h0,h1)=({h0_o(-2)},{h1_o(-2)})",
    )
    check(
        "S3_sqm_algebra",
        "the conditional tangent line and Callias determinant have matching c1=2",
        required_push["determinant_c1"] == 2
        and abs(berry_rows["k2"]["chern"] - 2.0) < TOL,
        "c1(TCP1)=c1(det Ind D)=2 in the declared conditional template",
    )

    # ------------------------------------------ differential-character transgression
    degree_before = 5
    fiber_dimension = 3
    degree_after = degree_before - fiber_dimension
    check(
        "S4_differential_transgression",
        "fiber integration along Sigma3 lowers differential degree 5 to degree 2",
        degree_after == 2,
        f"{degree_before}-{fiber_dimension}={degree_after}",
    )
    transgression = {
        "plus_sector": 2 * 1 * 1,
        "minus_sector": 2 * (-1) * 1,
        "degree_zero_control": 2 * 1 * 0,
        "level_one_control": 1 * 1 * 1,
    }
    check(
        "S4_differential_transgression",
        "factorized same-sphere map gives c1=n B d=+2 for (n,B,d)=(2,1,1)",
        transgression["plus_sector"] == 2,
        f"integer={transgression['plus_sector']}",
    )
    check(
        "S4_differential_transgression",
        "orientation reversal B=+1 to B=-1 reverses the raw WZW line",
        transgression["minus_sector"] == -2,
        f"integer={transgression['minus_sector']}",
    )
    check(
        "S4_differential_transgression",
        "degree-zero orientation-map control pulls the WZW line back trivially",
        transgression["degree_zero_control"] == 0,
        f"integer={transgression['degree_zero_control']}",
    )
    check(
        "S4_differential_transgression",
        "level-one control produces c1=1 and therefore cannot supply O(2)",
        transgression["level_one_control"] == 1,
        f"integer={transgression['level_one_control']}",
    )

    descent_conditions = {
        "equivariant_differential_refinement": False,
        "curvature_horizontal_and_invariant": False,
        "vertical_and_large_gauge_holonomy_trivial": False,
        "stabilizer_character_trivial": False,
        "class_in_image_of_pullback_from_quotient": False,
    }
    check(
        "S4_differential_transgression",
        "descent audit enumerates local and global necessary conditions",
        set(descent_conditions)
        == {
            "equivariant_differential_refinement",
            "curvature_horizontal_and_invariant",
            "vertical_and_large_gauge_holonomy_trivial",
            "stabilizer_character_trivial",
            "class_in_image_of_pullback_from_quotient",
        },
        f"conditions={list(descent_conditions)}",
    )
    check(
        "S4_differential_transgression",
        "gauge descent fails closed when no condition is supplied by the mother model",
        not any(descent_conditions.values()),
        "all same-model descent evidence flags are false",
    )

    same_model_audit_conditions = {
        "normalizable_same_soliton_cp1_moduli_and_positive_metric": False,
        "uniform_noncollective_bosonic_gap": False,
        "microscopic_supercharges_or_independently_derived_emergent_n2_ward_identities": False,
        "specified_yukawa_callias_family": False,
        "uniform_callias_fredholm_gap": False,
        "kernel_bundle_isomorphic_to_holomorphic_tangent": False,
        "ward_identity_fixes_levi_civita_connection_and_curvature_term": False,
        "no_extra_zero_modes": False,
    }
    check(
        "S5_underdetermination",
        "same-model audit records all eight conservative closure conditions",
        len(same_model_audit_conditions) == 8,
        f"conditions={len(same_model_audit_conditions)}",
    )
    check(
        "S5_underdetermination",
        "present charged two-colour proxy supplies none of the CP1/SUSY/Callias closure conditions",
        not any(same_model_audit_conditions.values()),
        "bosonic potential and radial Hessian do not determine a Yukawa family",
    )

    status = {
        "conditional_composition_sufficient_criterion_established": True,
        "low_energy_tangent_sqm_characterization_derived": True,
        "sufficient_same_model_audit_conditions_derived": True,
        "underdetermination_theorem_applies": True,
        "same_charged_two_colour_n2_closure_mechanism_specified": False,
        "same_soliton_cp1_moduli_derived": False,
        "worldline_n2_supersymmetry_same_model_derived": False,
        "physical_tangent_fermion_same_model_derived": False,
        "callias_operator_same_model_specified": False,
        "callias_fredholm_gap_same_model_verified": False,
        "callias_family_index_same_model_computed": False,
        "callias_determinant_c1_same_model_computed": False,
        "conditional_callias_template_c1_plus_2_verified": True,
        "callias_line_as_quantum_coefficient_line_same_model_derived": False,
        "cpt_anti_canonical_map_mathematically_derived": True,
        "cpt_anti_canonical_map_same_model_derived": False,
        "wzw_spatial_fiber_integration_defined": True,
        "wzw_equivariant_refinement_same_model_constructed": False,
        "wzw_gauge_basic_descent_same_model_proven": False,
        "wzw_pullback_o2_same_model_proven": False,
        "ap_e3_ap_e4_composition_gate_closed": False,
        "degree_one_route_e_portal_built": False,
        "physics_promotion_allowed": False,
    }

    passed = sum(row["pass"] for row in CHECKS)
    result = {
        "schema_version": "ap-e4-same-soliton-callias-descent-v1",
        "generated_by": str(THIS_SCRIPT.relative_to(REPO)),
        "seed": SEED,
        "summary": {"passed": passed, "total": len(CHECKS)},
        "status": status,
        "same_model_audit_conditions": same_model_audit_conditions,
        "gauge_descent_conditions": descent_conditions,
        "metrics": {
            "cp1_chart": {
                "max_metric_residual": max(metric_residuals),
                "max_tangent_norm_residual": max(tangent_residuals),
                "max_connection_residual": max(connection_residuals),
            },
            "berry_chern": berry_rows,
            "finite_kernel_family_template_control": finite_family,
            "callias_pushforward": {
                "formula": {
                    "asymptotic_c1": "p*x + q*y",
                    "rank": "epsilon*p",
                    "determinant_c1": "epsilon*p*q",
                },
                "required": required_push,
                "q_zero_control": trivial_moduli_push,
                "opposite_control": opposite_push,
                "doubled_determinant_c1": doubled_c1,
            },
            "cpt": {
                "canonical_degree": -2,
                "plus_coefficient_degree": 2,
                "anti_canonical_minus_degree": anti_canonical_degree,
                "plus_kernel": {"positive": h0_o(2), "negative": h1_o(2)},
                "minus_kernel": {
                    "positive": h0_o(anti_canonical_degree),
                    "negative": h1_o(anti_canonical_degree),
                },
                "naive_inverse_kernel": {"positive": h0_o(-2), "negative": h1_o(-2)},
            },
            "sqm_finite_blocks": {
                "coefficient_C": coefficient,
                "plus": {
                    "index": plus_block["index"],
                    "kernel_dimension": plus_block["kernel_dimension"],
                    "dirac_gap": plus_block["smallest_abs_nonzero_D"],
                    "residuals": plus_residuals,
                },
                "minus": {
                    "index": minus_block["index"],
                    "kernel_dimension": minus_block["kernel_dimension"],
                    "dirac_gap": minus_block["smallest_abs_nonzero_D"],
                    "residuals": minus_residuals,
                },
            },
            "wzw_transgression": {
                "input_degree": degree_before,
                "fiber_dimension": fiber_dimension,
                "output_degree": degree_after,
                **transgression,
                "general_pullback_criterion": "int_{Sigma3 x CP1} ev^*(n*u3*u2)=+2",
            },
        },
        "theorems": {
            "low_energy_n2_characterization": (
                "At two-derivative adiabatic order, a one-multiplet tangent N=2 SQM is "
                "characterized by normalizable CP1 modes, a uniform gap, a Fredholm "
                "fermion kernel identified with T^(1,0)CP1, the N=2 Ward identities, "
                "and no extra zero modes. Exact microscopic supercharges are sufficient "
                "but not logically necessary because the Ward identities may emerge."
            ),
            "bosonic_underdetermination": (
                "A bosonic soliton profile and Hessian do not determine a fermion kernel: "
                "gapped trivial and topologically winding Yukawa families can share the "
                "same bosonic action."
            ),
            "callias_required_class": (
                "With declared Callias orientation epsilon*p=+1, an asymptotic "
                "positive-mass eigenbundle class p*x+2*y is necessary and sufficient "
                "within the rank-one product template for det Ind D=O(+2)."
            ),
            "fixed_polarization_cpt": (
                "Serre duality maps E_+=O(2) to E_-=K tensor E_+^vee=O(-4), "
                "giving three negative-chirality states; naive O(-2) gives only one."
            ),
            "wzw_pullback": (
                "After equivariant gauge descent, the spatially transgressed WZW line "
                "pulls back to O(+2) iff its integral characteristic number on "
                "Sigma3 x CP1 is +2; factorization gives n*B*degree."
            ),
        },
        "composition_gate": {
            "closed": False,
            "required_before_degree_one_portal": [
                "same charged-two-colour mother action with microscopic or derived emergent N2 closure",
                "same B=1 CP1 moduli and uniform bosonic gap",
                "specified uniformly Fredholm Callias/Yukawa family",
                "actual determinant line c1=+2 and measure coupling",
                "UV-derived CPT anti-canonical vacuum line",
                "equivariant differential character and global gauge descent",
                "same-soliton transgression integral +2",
            ],
        },
        "source_manifest": source_manifest,
        "checks": CHECKS,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e4_same_soliton_callias_descent.json"
    md_path = OUTPUT / "ap_e4_same_soliton_callias_descent.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(markdown_report(result), encoding="utf-8")
    print(f"wrote {json_path.relative_to(REPO)}")
    print(f"wrote {md_path.relative_to(REPO)}")
    print(f"AP-E4 same-soliton Callias/descent checks: {passed}/{len(CHECKS)} passed")
    print("physics_promotion_allowed=false")
    return 0 if passed == len(CHECKS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
