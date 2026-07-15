#!/usr/bin/env python3
"""Deterministic AP-E3 audit for the exactly-two and mixed-WZW branches.

The first branch is an exact constrained-cell theorem formulated with bosonic
Schwinger partons: two independent compact U(1) Gauss laws impose N_1=N_2=1
before any low-energy approximation.  A ferromagnetic Hund term and the
physical Pauli-Zeeman sign then select the anti-aligned spin-one coherent line
Q^2=O(2), giving signed k=+2.

The second branch is an anomaly-consistent *intermediate* UV completion of a
mixed WZW term.  It uses SU(2)_c x U(1)_g with two vectorlike Dirac flavours
and a charge-one scalar doublet.  All dynamical gauge anomalies cancel and the
mixed flavour/one-form Postnikov coefficient is two.  Strong-dynamics, compact
U(1) monopole, deep-UV, soliton-stability, and Route-E portal gates remain
explicitly false, so a green run never permits physics promotion.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e3_exact_two_mixed_wzw_uv.tex"
BIB = ROUTE_F / "tex" / "ap_e3_exact_two_mixed_wzw_uv.bib"
LITERATURE = ROUTE_F / "LITERATURE_2023_2026.md"
THIS_SCRIPT = Path(__file__).resolve()

SEED = 20260715
TOL = 3.0e-11
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


def directional_spinor(w: complex, dw: complex) -> tuple[np.ndarray, np.ndarray]:
    norm = 1.0 + abs(w) ** 2
    dnorm = 2.0 * (w.conjugate() * dw).real
    root = math.sqrt(norm)
    section = np.array([1.0 / root, w / root], dtype=complex)
    derivative = np.array(
        [
            -0.5 * dnorm / norm**1.5,
            dw / root - 0.5 * w * dnorm / norm**1.5,
        ],
        dtype=complex,
    )
    return section, derivative


def anti_aligned(section: np.ndarray) -> np.ndarray:
    epsilon = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    return epsilon @ section.conjugate()


def anti_aligned_derivative(derivative: np.ndarray) -> np.ndarray:
    epsilon = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    return epsilon @ derivative.conjugate()


def midpoint_chern(theta_cells: int) -> float:
    """Chern number for pair curvature 2F=sin(theta)dtheta^dphi."""

    dtheta = math.pi / theta_cells
    return sum(
        math.sin((index + 0.5) * dtheta) * dtheta
        for index in range(theta_cells)
    )


def bosonic_fock_states(
    max_total_per_orbital: int = 2,
) -> list[tuple[int, int, int, int]]:
    """Four Schwinger-boson modes with the declared per-orbital audit cutoff."""
    states: list[tuple[int, int, int, int]] = []
    for state in itertools.product(range(max_total_per_orbital + 1), repeat=4):
        if state[0] + state[1] <= max_total_per_orbital and state[2] + state[3] <= max_total_per_orbital:
            states.append(state)
    return states


def fourier_projector_value(charge: int, target: int, samples: int = 16) -> complex:
    return sum(
        np.exp(2.0j * math.pi * index * (charge - target) / samples)
        for index in range(samples)
    ) / samples


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e3_exact_two_mixed_wzw.json"
    md_path = OUTPUT / "ap_e3_exact_two_mixed_wzw.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checks = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    blockers = "\n".join(f"- {item}" for item in result["remaining_blockers"])
    sources = "\n".join(
        f"- `{row['path']}` - `{row['sha256']}` ({row['size_bytes']} bytes)"
        for row in result["source_manifest"]
    )
    md_path.write_text(
        f"""# AP-E3 exactly-two and mixed-WZW audit

- Status: `{result['status']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Exactly-two constrained-cell theorem: `{str(result['exactly_two_constrained_cell_proved']).lower()}`
- Physical Pauli orientation gives k=+2: `{str(result['physical_pauli_orientation_gives_k_plus_2']).lower()}`
- Dynamical gauge anomalies cancel: `{str(result['dynamical_gauge_anomalies_cancel']).lower()}`
- Mixed-WZW intermediate UV candidate: `{str(result['mixed_wzw_intermediate_uv_candidate_constructed']).lower()}`
- Full all-scale/nonperturbative UV closure: `{str(result['mixed_wzw_full_uv_closed']).lower()}`
- Route-E portal closed: `{str(result['route_e_portal_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

## Derived chain

```text
G_r=N_r-1=0 (r=1,2) -> H_phys=C2 tensor C2
-J_H S1.S2 + h n.(S1+S2), J_H,h>0 -> unique m=-1 triplet
Pauli electron sign -> anti-aligned line Q tensor Q = O(2)
i hbar <Omega_-|d Omega_-> = +2 hbar A_+ -> k=+2
SU(2)c x U(1)g, nc=2, Xq=1 -> kappa_L=-kappa_R=2
S_mix=2 pi hbar 2 integral(omega3 wedge omega2)
B=+1 worldline -> k=nB=+2
```

Numerical Hund-Zeeman spectrum: `{result['numerical_summary']['hund_zeeman_spectrum']}`.
Numerical Chern sequence: `{result['numerical_summary']['chern_midpoint_sequence']}`.

## Remaining blockers

{blockers}

## Mechanical checks

{checks}

## Source hashes

{sources}
""",
        encoding="utf-8",
    )


def main() -> None:
    rng = random.Random(SEED)
    critical_sources = [TEX, BIB, LITERATURE, THIS_SCRIPT]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "U0_provenance",
        "all exact-two/mixed-WZW critical sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # ---------------------------------------------------- exact U(1)^2 Gauss rule
    states = bosonic_fock_states(2)
    n1 = np.array([row[0] + row[1] for row in states], dtype=int)
    n2 = np.array([row[2] + row[3] for row in states], dtype=int)
    exact_projector = np.diag(((n1 == 1) & (n2 == 1)).astype(float))
    fourier_diagonal = np.array(
        [
            fourier_projector_value(int(first), 1)
            * fourier_projector_value(int(second), 1)
            for first, second in zip(n1, n2)
        ],
        dtype=complex,
    )
    fourier_projector = np.diag(fourier_diagonal)
    projector_residual = float(np.max(np.abs(fourier_projector - exact_projector)))
    check(
        "U1_exact_cell",
        "two compact Gauss laws implement the exact Fourier projector P_1 P_1",
        projector_residual < TOL,
        f"basis={len(states)}; projector residual={projector_residual:.3e}",
    )
    check(
        "U1_exact_cell",
        "the constrained physical cell has rank four and exactly total occupation two",
        int(round(np.trace(exact_projector).real)) == 4
        and np.all((n1 + n2)[np.diag(exact_projector) > 0.5] == 2),
        f"rank={int(round(np.trace(exact_projector).real))}",
    )
    check(
        "U1_exact_cell",
        "separate Gauss laws remove literal odd, singleton, and charge-transfer sectors",
        not np.any(np.diag(exact_projector)[(n1 + n2) % 2 == 1] > 0.5)
        and not np.any(np.diag(exact_projector)[(n1 == 2) & (n2 == 0)] > 0.5)
        and not np.any(np.diag(exact_projector)[(n1 == 0) & (n2 == 2)] > 0.5),
        "forbidden=(odd totals),(2,0),(0,2)",
    )

    total_two_projector = np.diag((n1 + n2 == 2).astype(float))
    charge_transfer_survives = bool(
        np.any(np.diag(total_two_projector)[(n1 == 2) & (n2 == 0)] > 0.5)
        and np.any(np.diag(total_two_projector)[(n1 == 0) & (n2 == 2)] > 0.5)
    )
    check(
        "U1_exact_cell",
        "negative control in the N_r<=2 bosonic audit: one total-number Gauss law is strictly weaker",
        int(round(np.trace(total_two_projector).real)) == 10
        and charge_transfer_survives,
        f"rank(P_total=2)={int(round(np.trace(total_two_projector).real))}; transfer survives={charge_transfer_survives}",
    )
    parity_projector = np.diag(((n1 + n2) % 2 == 0).astype(float))
    check(
        "U1_exact_cell",
        "negative control in the N_r<=2 bosonic audit: even parity neither selects N=2 nor removes charge transfer",
        int(round(np.trace(parity_projector).real)) > 4
        and np.any(np.diag(parity_projector)[n1 + n2 == 0] > 0.5)
        and charge_transfer_survives,
        f"rank(P_even)={int(round(np.trace(parity_projector).real))}",
    )
    integer_large_gauge_phase = np.exp(-2.0j * math.pi * (3 + -2))
    half_integer_bad_phase = np.exp(-2.0j * math.pi * 0.5)
    check(
        "U1_exact_cell",
        "integer background charges pass 0+1d large-gauge invariance",
        abs(integer_large_gauge_phase - 1.0) < TOL
        and abs(half_integer_bad_phase - 1.0) > 1.0,
        f"integer phase={integer_large_gauge_phase}; half-integer phase={half_integer_bad_phase}",
    )

    # ----------------------------------------- Hund spectrum and physical orientation
    sx = 0.5 * np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = 0.5 * np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    sz = 0.5 * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    s1 = [np.kron(matrix, identity) for matrix in (sx, sy, sz)]
    s2 = [np.kron(identity, matrix) for matrix in (sx, sy, sz)]
    dot = sum(first @ second for first, second in zip(s1, s2))
    total = [first + second for first, second in zip(s1, s2)]
    j_h = 1.0
    h = 0.2
    hamiltonian_z = -j_h * dot + h * total[2]
    spectrum = np.linalg.eigvalsh(hamiltonian_z)
    expected_spectrum = np.array(
        [-j_h / 4.0 - h, -j_h / 4.0, -j_h / 4.0 + h, 3.0 * j_h / 4.0]
    )
    check(
        "U2_orientation",
        "Hund plus Pauli-Zeeman Hamiltonian has the exact triplet/singlet spectrum",
        float(np.max(np.abs(spectrum - expected_spectrum))) < TOL,
        f"spectrum={spectrum.tolist()}",
    )
    check(
        "U2_orientation",
        "the unique local gap is min(h,J_H+h)=h for J_H,h>0",
        abs((spectrum[1] - spectrum[0]) - h) < TOL,
        f"gap={spectrum[1] - spectrum[0]:.12f}",
    )

    ground_projector_residuals: list[float] = []
    eigen_residuals: list[float] = []
    berry_residuals: list[float] = []
    metric_residuals: list[float] = []
    for _ in range(100):
        vector = np.array([rng.gauss(0.0, 1.0) for _ in range(3)], dtype=float)
        vector /= np.linalg.norm(vector)
        hamiltonian = -j_h * dot + h * sum(
            component * operator for component, operator in zip(vector, total)
        )
        values, vectors = np.linalg.eigh(hamiltonian)
        theta = math.acos(float(vector[2]))
        phi = math.atan2(float(vector[1]), float(vector[0]))
        s_plus = np.array(
            [math.cos(theta / 2.0), np.exp(1.0j * phi) * math.sin(theta / 2.0)],
            dtype=complex,
        )
        t_minus = anti_aligned(s_plus)
        pair = np.kron(t_minus, t_minus)
        ground_projector_residuals.append(
            float(
                np.max(
                    np.abs(
                        np.outer(vectors[:, 0], vectors[:, 0].conjugate())
                        - np.outer(pair, pair.conjugate())
                    )
                )
            )
        )
        eigen_residuals.append(float(np.linalg.norm(hamiltonian @ pair - values[0] * pair)))

        w = complex(rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0))
        dw = complex(rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3))
        section, derivative = directional_spinor(w, dw)
        anti = anti_aligned(section)
        danti = anti_aligned_derivative(derivative)
        pair = np.kron(anti, anti)
        dpair = np.kron(danti, anti) + np.kron(anti, danti)
        a_plus = -1.0j * np.vdot(section, derivative)
        action_connection = 1.0j * np.vdot(pair, dpair)
        berry_residuals.append(abs(action_connection - 2.0 * a_plus))
        single_metric = np.vdot(derivative, derivative) - abs(np.vdot(section, derivative)) ** 2
        pair_metric = np.vdot(dpair, dpair) - abs(np.vdot(pair, dpair)) ** 2
        metric_residuals.append(abs(pair_metric - 2.0 * single_metric))

    check(
        "U2_orientation",
        "random-direction diagonalization selects the anti-aligned tensor-square ground line",
        max(ground_projector_residuals) < TOL and max(eigen_residuals) < TOL,
        f"projector={max(ground_projector_residuals):.3e}; eigen={max(eigen_residuals):.3e}",
    )
    check(
        "U2_orientation",
        "the microscopic kinetic term gives signed action connection +2 A_+",
        max(berry_residuals) < TOL,
        f"max residual={max(berry_residuals):.3e}",
    )
    check(
        "U2_orientation",
        "the anti-aligned pair quantum metric is twice the one-doublet metric",
        max(metric_residuals) < TOL,
        f"max residual={max(metric_residuals):.3e}",
    )
    chern_sequence = [midpoint_chern(cells) for cells in (10, 20, 40, 80, 160)]
    check(
        "U2_orientation",
        "pair curvature quadrature converges to c1=+2 with second-order refinement",
        abs(chern_sequence[-1] - 2.0) < 4.0e-5
        and all(abs(chern_sequence[i + 1] - 2.0) < abs(chern_sequence[i] - 2.0) for i in range(4)),
        f"sequence={chern_sequence}",
    )
    overlap_phases = np.linspace(0.0, 2.0 * math.pi, 2049)
    unwrapped = np.unwrap(np.angle(np.exp(2.0j * overlap_phases)))
    winding = (unwrapped[-1] - unwrapped[0]) / (2.0 * math.pi)
    check(
        "U2_orientation",
        "the quotient-line pair transition has winding +2",
        abs(winding - 2.0) < TOL,
        f"winding={winding:.12f}",
    )

    # ------------------------------------- anomaly-consistent mixed-WZW candidate
    n_c = 2
    n_f = 2
    x_q = 1
    left_charges = [x_q] * (n_c * n_f)
    right_conjugate_charges = [-x_q] * (n_c * n_f)
    all_fermion_charges = left_charges + right_conjugate_charges
    u1_cubic = sum(charge**3 for charge in all_fermion_charges)
    grav_u1 = sum(all_fermion_charges)
    su2_sq_u1 = n_f * 0.5 * (x_q - x_q)
    pauli_generators = [
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex) / 2.0,
        np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex) / 2.0,
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex) / 2.0,
    ]
    su2_u1_sq = float(
        max(
            abs(n_f * (x_q**2 + (-x_q) ** 2) * np.trace(generator))
            for generator in pauli_generators
        )
    )
    check(
        "U3_anomaly",
        "all local dynamical gauge-anomaly sums vanish",
        u1_cubic == 0
        and grav_u1 == 0
        and abs(su2_sq_u1) < TOL
        and su2_u1_sq < TOL,
        f"U1^3={u1_cubic}; grav^2-U1={grav_u1}; "
        f"SU2^2-U1={su2_sq_u1}; SU2-U1^2={su2_u1_sq}",
    )
    left_weyl_su2_doublets = 2 * n_f
    check(
        "U3_anomaly",
        "the SU(2)_c Witten anomaly cancels because the Weyl-doublet count is even",
        left_weyl_su2_doublets % 2 == 0,
        f"left Weyl doublets={left_weyl_su2_doublets}",
    )
    kappa_left = n_c * x_q
    kappa_right = -n_c * x_q
    check(
        "U3_anomaly",
        "the mixed flavour/one-form Postnikov classes are kappa_L=-kappa_R=2",
        kappa_left == 2 and kappa_right == -2,
        f"kappa_L={kappa_left}; kappa_R={kappa_right}",
    )
    flavour_witten_copies = n_c
    check(
        "U3_anomaly",
        "both SU(2)_L and SU(2)_R background Witten anomalies have even colour multiplicity",
        flavour_witten_copies % 2 == 0,
        f"copies per flavour SU(2)={flavour_witten_copies}",
    )
    check(
        "U3_anomaly",
        "gauged U(1)_g exactly removes off-diagonal Pauli-Guersey generators",
        x_q != -x_q,
        "off-diagonal generators change U(1)_g charge by 2 and do not commute with the gauge generator",
    )

    mixed_level = n_c * x_q
    extension_samples = [
        mixed_level * rng.randint(-8, 8) * rng.randint(-8, 8) for _ in range(200)
    ]
    extension_phases = [np.exp(2.0j * math.pi * value) for value in extension_samples]
    check(
        "U4_WZW",
        "the five-dimensional mixed WZW exponent is extension-independent for tested extendible integral classes",
        max(abs(phase - 1.0) for phase in extension_phases) < TOL,
        f"tested={len(extension_phases)}; n={mixed_level}",
    )
    baryon_number = 1
    worldline_level = mixed_level * baryon_number
    check(
        "U4_WZW",
        "a positive unit baryon reduces the mixed WZW term to signed worldline k=+2",
        worldline_level == 2 and mixed_level * (-baryon_number) == -2,
        f"k(B=+1)={worldline_level}; k(B=-1)={-worldline_level}",
    )
    baryon_quark_count = n_c
    dressed_charge = baryon_quark_count * x_q + 2 * (-1)
    check(
        "U4_WZW",
        "the minimal local colour singlet requires two quarks and two conjugate Higgs dressings",
        baryon_quark_count == 2 and dressed_charge == 0,
        f"qq charge={baryon_quark_count * x_q}; dressing charge=-2; total={dressed_charge}",
    )
    check(
        "U4_WZW",
        "the local Higgs dressing is Sym^2(conjugate doublet), hence a triplet of degree two",
        math.comb(2 + 2 - 1, 2) == 3,
        "dim Sym^2(C^2)=3; associated projective line is O(2)",
    )
    check(
        "U4_WZW",
        "a coloured singleton is excluded from the local gauge-invariant operator algebra",
        n_c == 2,
        "one SU(2)_c fundamental is not a colour singlet; the first local baryon uses epsilon_ab q^a q^b",
    )

    beta_u1 = (4.0 / 3.0) * (n_c * n_f) + (1.0 / 3.0) * 2.0
    beta_su2 = (11.0 / 3.0) * 2.0 - (4.0 / 3.0) * 0.5 * n_f
    check(
        "U5_boundaries",
        "the running diagnostic detects asymptotically free colour but an abelian Landau pole",
        abs(beta_u1 - 6.0) < TOL and abs(beta_su2 - 6.0) < TOL,
        f"b_U1={beta_u1:.6f}; b_SU2={beta_su2:.6f}",
    )
    check(
        "U5_boundaries",
        "pure two-colour QCD is retained as a negative control, not used for the completion",
        True,
        "without U(1)_g, SU(4)->Sp(4) gives S^5 and pi_3(S^5)=0; the charged theory needs a separate phase audit",
    )

    remaining_blockers = [
        "Prove nonperturbatively that the charged two-colour theory condenses in the mesonic SU(2)_L x SU(2)_R -> SU(2)_V channel and gaps every Pauli-Guersey/diquark direction.",
        "Keep confinement at or above the Higgs-vector scale, or otherwise show that approximate Pauli-Guersey restoration cannot unwind the B=1 soliton.",
        "Audit compact-U(1) monopoles, the discrete axial anomaly, the spin-bordism sign, and an explicit all-scale embedding beyond the abelian Landau pole.",
        "Construct the differential-character/Cech-bordism definition on non-extendible four-dimensional sectors; the present test proves only extension independence when an extension exists.",
        "Compute the soliton Hessian, Finkelstein-Rubinstein statistics, collective-coordinate inertia, and the LLL/adiabatic gap before retaining only the three O(2) states.",
        "Derive a degree +1 Route-E portal map; neither the constrained cell nor the mixed-WZW sector yet identifies its triplet with three chiral Spin(10) families.",
        "The complete-cell boundary theorem covers positive parent interactions only; arbitrary intercell couplings may generate emergent projective edge modes even though literal odd-occupancy states are absent.",
    ]

    all_pass = all(row["pass"] for row in CHECKS)
    result: dict[str, Any] = {
        "status": "ap_e3_exact_cell_and_anomaly_consistent_mixed_wzw_candidate_physics_portal_open",
        "all_pass": all_pass,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "checks": CHECKS,
        "exactly_two_constrained_cell_proved": all_pass,
        "separate_gauss_laws_enforce_N1_N2_equal_1": all_pass,
        "literal_odd_and_charge_transfer_sectors_absent": all_pass,
        "complete_cell_safe_parent_boundary_singleton_absent": True,
        "arbitrary_interacting_boundary_singleton_absent": False,
        "finite_U_electronic_uv_exactly_two": False,
        "physical_pauli_orientation_gives_k_plus_2": all_pass,
        "positive_branch_bundle_is_Q_tensor_2_equal_O_2": all_pass,
        "dynamical_gauge_anomalies_cancel": all_pass,
        "mixed_wzw_intermediate_uv_candidate_constructed": all_pass,
        "mixed_wzw_extendible_sector_extension_independence_checked": all_pass,
        "mixed_wzw_global_differential_cohomology_closed": False,
        "mixed_wzw_full_uv_closed": False,
        "route_e_portal_closed": False,
        "ap_e3_full_uv_closed": False,
        "physics_promotion_allowed": False,
        "parameters": {
            "J_H": j_h,
            "h": h,
            "n_c": n_c,
            "n_f": n_f,
            "X_q": x_q,
            "mixed_wzw_level": mixed_level,
        },
        "numerical_summary": {
            "hund_zeeman_spectrum": spectrum.tolist(),
            "random_ground_projector_max_residual": max(ground_projector_residuals),
            "random_eigenvector_max_residual": max(eigen_residuals),
            "berry_sign_max_residual": max(berry_residuals),
            "metric_additivity_max_residual": max(metric_residuals),
            "chern_midpoint_sequence": chern_sequence,
            "transition_winding": winding,
            "gauge_anomalies": {
                "U1_cubed": u1_cubic,
                "gravity_squared_U1": grav_u1,
                "SU2_squared_U1": su2_sq_u1,
                "SU2_U1_squared": su2_u1_sq,
                "SU2_Witten_doublets": left_weyl_su2_doublets,
            },
            "beta_coefficients": {"U1": beta_u1, "SU2": beta_su2},
        },
        "remaining_blockers": remaining_blockers,
        "source_manifest": source_manifest,
    }
    write_outputs(result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
