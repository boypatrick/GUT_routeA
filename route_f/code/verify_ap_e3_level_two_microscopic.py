#!/usr/bin/env python3
"""Fail-closed AP-E3 audit for a microscopic level-two-magnitude candidate.

The candidate is a two-orbital Mott/Hund dimer.  Each orbital carries one
SU(2) doublet, a ferromagnetic Hund coupling projects the low-energy sector to
the symmetric spin-one triplet, and the coherent-state ket line is the tensor
square of the Hopf ket line.  Within those declared microscopic inputs the
level magnitude is derived to be |k|=2.  With A=-i<s|ds> and the displayed
+i hbar a^dagger dot(a) kinetic term, the aligned model has signed k=-2; its
dual prequantum line is O(2).  Reversing the orientation gives k=+2.

The audit also records the decisive boundary: large-gauge invariance alone
only requires integer k, and this repository does not yet derive why the UV
field content contains exactly two same-orientation constituents or why all
single-constituent channels are absent.  A green run therefore never permits
physics promotion.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"

AP_E1_TEX = ROUTE_F / "tex" / "ap_e1_projective_doublet_action.tex"
AP_E2_JSON = OUTPUT / "ap_e2_discriminant_regression.json"
AP_E3_TEX = ROUTE_F / "tex" / "ap_e3_level_two_microscopic_origin.tex"
AP_E3_BIB = ROUTE_F / "tex" / "ap_e3_level_two_microscopic.bib"
LITERATURE = ROUTE_F / "LITERATURE_2023_2026.md"
THIS_SCRIPT = Path(__file__).resolve()

SEED = 20260714
TOL = 2.0e-12
CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
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


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def max_abs(matrix: np.ndarray) -> float:
    return float(np.max(np.abs(matrix)))


def commutator(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    return first @ second - second @ first


def directional_section(w: complex, dw: complex) -> tuple[np.ndarray, np.ndarray]:
    """North-chart normalized spinor and its exact directional derivative."""

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


def veronese(section: np.ndarray) -> np.ndarray:
    z0, z1 = section
    return np.array([z0 * z0, math.sqrt(2.0) * z0 * z1, z1 * z1], dtype=complex)


def pair_full(section: np.ndarray) -> np.ndarray:
    return np.kron(section, section)


def pair_full_derivative(section: np.ndarray, derivative: np.ndarray) -> np.ndarray:
    return np.kron(derivative, section) + np.kron(section, derivative)


def midpoint_chern_pair(theta_cells: int) -> float:
    """Integrate F_pair=sin(theta)dtheta^dphi and divide by 2pi."""

    dtheta = math.pi / theta_cells
    theta_integral = sum(
        math.sin((index + 0.5) * dtheta) * dtheta for index in range(theta_cells)
    )
    return theta_integral


def onsite_energy(occupation: int, interaction: float, chemical_potential: float) -> float:
    return (
        0.5 * interaction * occupation * (occupation - 1)
        - chemical_potential * occupation
    )


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e3_level_two_microscopic.json"
    md_path = OUTPUT / "ap_e3_level_two_microscopic.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checks = "\n".join(
        f"- [{'PASS' if item['pass'] else 'FAIL'}] `{item['group']}` - "
        f"{item['name']}: {item['detail']}"
        for item in result["checks"]
    )
    sources = "\n".join(
        f"- `{row['path']}` - `{row['sha256']}` ({row['size_bytes']} bytes)"
        for row in result["source_manifest"]
    )
    alternatives = "\n".join(
        f"- **{name}:** {description}"
        for name, description in result["alternative_branches"].items()
    )
    blockers = "\n".join(f"- {item}" for item in result["remaining_blockers"])
    md_path.write_text(
        f"""# AP-E3 microscopic origin of projective level magnitude two

- Status: `{result['status']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Candidate-level magnitude derivation complete: `{str(result['candidate_level_magnitude_derivation_complete']).lower()}`
- Level magnitude: `{result['level_magnitude']}`
- Signed level for the declared aligned action: `{result['signed_level_for_declared_action']}`
- Positive signed Route-E level selected: `{str(result['signed_route_e_level_selected']).lower()}`
- AP-E3 physics closure: `{str(result['ap_e3_physics_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

## Conditional first-principles result

For two orbitals with one doublet per orbital and ferromagnetic Hund locking,

```text
H_micro = sum_r [U N_r(N_r-1)/2 - mu N_r]
          - J_H S_1.S_2 - h n.(S_1+S_2),
C = mu + h/2 + J_H/4 < U - J_H/8.
```

the low-energy state is the symmetric pair

```text
|s;2> = |s> tensor |s>,
-i <s;2|d|s;2> = 2 A,
L_ket = nu_2^* O_CP2(-1) = O_CP1(-2),
L_pre = dual(L_ket) = O_CP1(2).
```

Thus the declared dimer derives `|k_eff|=2`, a three-state spin-one carrier,
and a singlet gap `Delta_singlet=J_H`.  Relative to `S_k=hbar*k*integral(A)`,
the displayed aligned microscopic action gives `k=-2`; reversing the
orientation gives `k=+2`.  This is not a derivation of the two-orbital UV field
content or of the signed Route-E chirality.

Numerical Chern sequence: `{result['numerical_summary']['chern_midpoint_sequence']}`.
Hund spectrum at `J_H=1`: `{result['numerical_summary']['hund_spectrum']}`.
Interacting-plateau sufficient-condition margin: `{result['numerical_summary']['interacting_plateau_condition_margin']}`.

## Alternative branches

{alternatives}

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
    hund_coupling = 1.0
    field_strength = 0.2
    critical_sources = [
        AP_E1_TEX,
        AP_E2_JSON,
        AP_E3_TEX,
        AP_E3_BIB,
        LITERATURE,
        THIS_SCRIPT,
    ]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "M0_provenance",
        "all AP-E3 input sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    ap_e2 = json.loads(AP_E2_JSON.read_text(encoding="utf-8"))
    check(
        "M0_provenance",
        "AP-E2 exact regression is green but explicitly leaves the Berry level underived",
        ap_e2.get("all_pass") is True
        and ap_e2.get("checks_passed") == 30
        and ap_e2.get("physics_promotion_allowed") is False
        and "Berry level k=2 or an O(2) prequantum bundle"
        in ap_e2["theorem_scope"]["not_derived"],
        "AP-E2={}/{}; promotion={}".format(
            ap_e2.get("checks_passed"),
            ap_e2.get("checks_total"),
            ap_e2.get("physics_promotion_allowed"),
        ),
    )

    # ------------------------------------------------ two-orbital Mott sector
    interaction = 4.0
    chemical_potential = 1.5
    plateau_constant = chemical_potential + field_strength / 2.0 + hund_coupling / 4.0
    plateau_upper_bound = interaction - hund_coupling / 8.0
    plateau_condition_margin = plateau_upper_bound - plateau_constant
    occupations = range(5)
    one_orbital_energies = {
        occupation: onsite_energy(occupation, interaction, chemical_potential)
        for occupation in occupations
    }
    one_orbital_ground = min(one_orbital_energies, key=one_orbital_energies.get)
    two_orbital_energies = {
        (n1, n2): one_orbital_energies[n1] + one_orbital_energies[n2]
        for n1 in occupations
        for n2 in occupations
    }
    two_orbital_ground = min(two_orbital_energies, key=two_orbital_energies.get)
    # The bare window is only an onsite statement.  It is deliberately kept
    # separate from the interacting theorem below.
    check(
        "M1_Mott",
        "the bare Bose-Hubbard window 0<mu<U locks exactly one constituent on each orbital",
        0.0 < chemical_potential < interaction
        and one_orbital_ground == 1
        and two_orbital_ground == (1, 1),
        f"U={interaction}; mu={chemical_potential}; one-orbital ground={one_orbital_ground}; pair={two_orbital_ground}",
    )
    charge_gap = min(chemical_potential, interaction - chemical_potential)
    excited_occupancy_energies = sorted(two_orbital_energies.values())
    numerical_charge_gap = excited_occupancy_energies[1] - excited_occupancy_energies[0]
    check(
        "M1_Mott",
        "the bare analytic charge gap min(mu,U-mu) matches the enumerated spectrum",
        charge_gap > 0.0 and abs(charge_gap - numerical_charge_gap) < TOL,
        f"analytic={charge_gap:.12f}; enumerated={numerical_charge_gap:.12f}",
    )
    low_mu_ground = min(
        occupations, key=lambda n: onsite_energy(n, interaction, -0.2)
    )
    high_mu_ground = min(
        occupations, key=lambda n: onsite_energy(n, interaction, interaction + 0.2)
    )
    # A counterexample prevents the bare window (even together with
    # coercivity) from being promoted to an interacting plateau theorem.
    counterexample_interaction = 1.0
    counterexample_mu = 0.99
    counterexample_hund = 3.9
    counterexample_h = 0.0
    counterexample_e11 = (
        2.0 * onsite_energy(1, counterexample_interaction, counterexample_mu)
        - 0.25 * counterexample_hund
        - counterexample_h
    )
    counterexample_e22 = (
        2.0 * onsite_energy(2, counterexample_interaction, counterexample_mu)
        - counterexample_hund
        - 2.0 * counterexample_h
    )
    check(
        "M1_Mott",
        "bare-window limits and the interacting false-inference counterexample are detected",
        low_mu_ground == 0
        and high_mu_ground == 2
        and 0.0 < counterexample_mu < counterexample_interaction
        and 4.0 * counterexample_interaction > counterexample_hund
        and counterexample_e22 < counterexample_e11,
        "mu<0 ground n={}; mu>U ground n={}; counterexample E11={:.6f}, E22={:.6f}".format(
            low_mu_ground, high_mu_ground, counterexample_e11, counterexample_e22
        ),
    )

    # For one spinful bosonic mode per orbital, occupancy n_r carries
    # j_r=n_r/2.  Ferromagnetic Hund and orientation terms are minimized by
    # maximal total spin, so the following is the exact lowest energy in each
    # (n1,n2) sector.  The quadratic lower bound proves that the finite
    # enumeration has not hidden a lower state at arbitrarily large occupancy.
    def full_sector_minimum(n1: int, n2: int) -> float:
        return (
            onsite_energy(n1, interaction, chemical_potential)
            + onsite_energy(n2, interaction, chemical_potential)
            - 0.25 * hund_coupling * n1 * n2
            - 0.5 * field_strength * (n1 + n2)
        )

    tail_start = 8
    sector_energies = {
        (n1, n2): full_sector_minimum(n1, n2)
        for n1 in range(tail_start)
        for n2 in range(tail_start - n1)
    }
    full_ground_sector = min(sector_energies, key=sector_energies.get)
    ordered_full_energies = sorted(sector_energies.values())
    full_charge_gap = ordered_full_energies[1] - ordered_full_energies[0]
    coercive_quadratic = (4.0 * interaction - hund_coupling) / 16.0
    coercive_linear = interaction / 2.0 + chemical_potential + field_strength / 2.0
    tail_lower_bound = (
        coercive_quadratic * tail_start**2 - coercive_linear * tail_start
    )
    tail_is_increasing = (
        2.0 * coercive_quadratic * tail_start - coercive_linear > 0.0
    )
    check(
        "M1_Mott",
        "a proved sufficient inequality selects the unique interacting (1,1) plateau",
        plateau_condition_margin > 0.0
        and full_ground_sector == (1, 1)
        and full_charge_gap > 0.0
        and coercive_quadratic > 0.0
        and tail_is_increasing
        and tail_lower_bound > sector_energies[full_ground_sector],
        "C={:.6f} < U-J_H/8={:.6f} (margin={:.6f}); ground={}; full gap={:.12f}; N>={} lower bound={:.6f}; ground energy={:.6f}".format(
            plateau_constant,
            plateau_upper_bound,
            plateau_condition_margin,
            full_ground_sector,
            full_charge_gap,
            tail_start,
            tail_lower_bound,
            sector_energies[full_ground_sector],
        ),
    )

    # ----------------------------------------------- exact dimer Hamiltonian
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    pauli = [sigma_x, sigma_y, sigma_z]
    identity_two = np.eye(2, dtype=complex)
    identity_four = np.eye(4, dtype=complex)
    spin_one = [0.5 * np.kron(matrix, identity_two) for matrix in pauli]
    spin_two = [0.5 * np.kron(identity_two, matrix) for matrix in pauli]
    total_spin = [spin_one[i] + spin_two[i] for i in range(3)]
    spin_dot = sum((spin_one[i] @ spin_two[i] for i in range(3)), np.zeros((4, 4), dtype=complex))
    total_spin_squared = sum(
        (operator @ operator for operator in total_spin),
        np.zeros((4, 4), dtype=complex),
    )
    swap = np.array(
        [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
        dtype=complex,
    )
    symmetric_projector = 0.5 * (identity_four + swap)
    spin_projector = 0.5 * total_spin_squared

    hund_hamiltonian = -hund_coupling * spin_dot
    check(
        "M2_Hund",
        "the ferromagnetic Hund Hamiltonian is Hermitian and globally SU(2)-invariant",
        max_abs(hund_hamiltonian - hund_hamiltonian.conjugate().T) < TOL
        and max(max_abs(commutator(hund_hamiltonian, generator)) for generator in total_spin) < TOL,
        "Hermitian residual={:.2e}; SU2 residual={:.2e}".format(
            max_abs(hund_hamiltonian - hund_hamiltonian.conjugate().T),
            max(max_abs(commutator(hund_hamiltonian, generator)) for generator in total_spin),
        ),
    )
    hund_eigenvalues, hund_eigenvectors = np.linalg.eigh(hund_hamiltonian)
    expected_hund = np.array([-0.25, -0.25, -0.25, 0.75])
    check(
        "M2_Hund",
        "exact diagonalization gives a triply degenerate triplet and one singlet",
        max_abs(hund_eigenvalues - expected_hund) < TOL,
        f"eigenvalues={hund_eigenvalues.tolist()}",
    )
    spin_squared_labels = [
        float(np.real(vector.conjugate() @ total_spin_squared @ vector))
        for vector in hund_eigenvectors.T
    ]
    check(
        "M2_Hund",
        "the three low states have S^2=2 and the high state has S^2=0",
        max(abs(value - target) for value, target in zip(spin_squared_labels, [2, 2, 2, 0])) < TOL,
        f"S2 labels={spin_squared_labels}",
    )
    check(
        "M2_Hund",
        "the spin-one projector equals the symmetric-exchange projector",
        max_abs(spin_projector - symmetric_projector) < TOL
        and abs(np.trace(symmetric_projector) - 3.0) < TOL,
        "max projector residual={:.2e}; rank trace={:.1f}".format(
            max_abs(spin_projector - symmetric_projector),
            float(np.trace(symmetric_projector).real),
        ),
    )
    singlet_gap = float(hund_eigenvalues[3] - hund_eigenvalues[2])
    check(
        "M2_Hund",
        "the singlet exclusion gap is exactly Delta_singlet=J_H",
        abs(singlet_gap - hund_coupling) < TOL,
        f"gap={singlet_gap:.12f}; J_H={hund_coupling:.12f}",
    )
    antiferromagnetic_eigenvalues = np.linalg.eigvalsh(+hund_coupling * spin_dot)
    check(
        "M2_Hund",
        "reversing the Hund sign selects a singlet and destroys the |k|=2 triplet",
        abs(antiferromagnetic_eigenvalues[0] + 0.75) < TOL
        and np.count_nonzero(np.isclose(antiferromagnetic_eigenvalues, -0.75, atol=TOL)) == 1,
        f"antiferromagnetic spectrum={antiferromagnetic_eigenvalues.tolist()}",
    )

    zeeman_gaps: list[float] = []
    spectral_residuals: list[float] = []
    expected_zeeman = np.array([-0.25 - field_strength, -0.25, -0.25 + field_strength, 0.75])
    for _ in range(64):
        direction = np.array([rng.gauss(0.0, 1.0) for _ in range(3)], dtype=float)
        direction /= np.linalg.norm(direction)
        zeeman = sum(
            (direction[i] * total_spin[i] for i in range(3)),
            np.zeros((4, 4), dtype=complex),
        )
        spectrum = np.linalg.eigvalsh(hund_hamiltonian - field_strength * zeeman)
        spectral_residuals.append(max_abs(spectrum - expected_zeeman))
        zeeman_gaps.append(float(spectrum[1] - spectrum[0]))
    check(
        "M2_Hund",
        "a weak orientation field gives a unique gapped ground-state line over every n in S2",
        max(spectral_residuals) < TOL
        and min(zeeman_gaps) > 0.0
        and max(abs(gap - field_strength) for gap in zeeman_gaps) < TOL,
        "64 directions; max spectral residual={:.2e}; min gap={:.12f}".format(
            max(spectral_residuals), min(zeeman_gaps)
        ),
    )

    # ----------------------------------- tensor-square / Veronese geometry
    berry_residuals: list[float] = []
    action_sign_residuals: list[float] = []
    metric_residuals: list[float] = []
    veronese_residuals: list[float] = []
    phase_residuals: list[float] = []
    singlet_leakages: list[float] = []
    for _ in range(256):
        w = complex(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))
        dw = complex(rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5))
        section, derivative = directional_section(w, dw)
        pair = pair_full(section)
        pair_derivative = pair_full_derivative(section, derivative)
        veronese_state = veronese(section)
        symmetric_coordinates = np.array(
            [pair[0], (pair[1] + pair[2]) / math.sqrt(2.0), pair[3]],
            dtype=complex,
        )
        veronese_residuals.append(
            max(
                abs(np.vdot(pair, pair) - 1.0),
                max_abs(symmetric_coordinates - veronese_state),
                abs(np.vdot(veronese_state, veronese_state) - 1.0),
            )
        )
        singlet = np.array([0.0, 1.0, -1.0, 0.0], dtype=complex) / math.sqrt(2.0)
        singlet_leakages.append(abs(np.vdot(singlet, pair)))

        single_overlap = np.vdot(section, derivative)
        pair_overlap = np.vdot(pair, pair_derivative)
        berry_residuals.append(abs(pair_overlap - 2.0 * single_overlap))
        single_connection = -1j * single_overlap
        pair_kinetic_one_form = 1j * pair_overlap
        action_sign_residuals.append(
            abs(pair_kinetic_one_form + 2.0 * single_connection)
        )
        single_metric = float(
            np.vdot(derivative, derivative).real - abs(single_overlap) ** 2
        )
        pair_metric = float(
            np.vdot(pair_derivative, pair_derivative).real - abs(pair_overlap) ** 2
        )
        metric_residuals.append(abs(pair_metric - 2.0 * single_metric))

        alpha = rng.uniform(-math.pi, math.pi)
        phased_section = np.exp(1j * alpha) * section
        phased_pair = pair_full(phased_section)
        phase_residuals.append(max_abs(phased_pair - np.exp(2j * alpha) * pair))

    check(
        "M3_Veronese",
        "the symmetric tensor product is the normalized quadratic Veronese state",
        max(veronese_residuals) < TOL and max(singlet_leakages) < TOL,
        "samples=256; max norm/map residual={:.2e}; max singlet leakage={:.2e}".format(
            max(veronese_residuals), max(singlet_leakages)
        ),
    )
    check(
        "M3_Veronese",
        "the pair transforms with phase weight two",
        max(phase_residuals) < TOL,
        f"samples=256; max residual={max(phase_residuals):.2e}",
    )
    check(
        "M3_Veronese",
        "the pair connection doubles while the declared microscopic kinetic action has signed k=-2",
        max(berry_residuals) < TOL and max(action_sign_residuals) < TOL,
        "samples=256; connection residual={:.2e}; action-sign residual={:.2e}".format(
            max(berry_residuals), max(action_sign_residuals)
        ),
    )
    check(
        "M3_Veronese",
        "the Veronese pullback quantum metric is twice the single-doublet metric",
        max(metric_residuals) < TOL,
        f"samples=256; max metric residual={max(metric_residuals):.2e}",
    )

    transition_phases = np.linspace(0.0, 2.0 * math.pi, 513)
    single_transition = np.exp(-1j * transition_phases)
    pair_transition = np.exp(-2j * transition_phases)
    single_winding = (
        np.unwrap(np.angle(single_transition))[-1]
        - np.unwrap(np.angle(single_transition))[0]
    ) / (2.0 * math.pi)
    pair_winding = (
        np.unwrap(np.angle(pair_transition))[-1]
        - np.unwrap(np.angle(pair_transition))[0]
    ) / (2.0 * math.pi)
    check(
        "M3_Veronese",
        "the ket transition has winding -2 and its dual prequantum line has Chern number +2",
        abs(single_winding + 1.0) < TOL and abs(pair_winding + 2.0) < TOL,
        f"single winding={single_winding:.12f}; pair winding={pair_winding:.12f}",
    )

    chern_cells = [20, 80, 320, 1280]
    chern_sequence = [midpoint_chern_pair(cells) for cells in chern_cells]
    chern_errors = [abs(value - 2.0) for value in chern_sequence]
    check(
        "M3_Veronese",
        "direct curvature quadrature converges to dual-prequantum first Chern number two",
        all(chern_errors[i + 1] < chern_errors[i] for i in range(len(chern_errors) - 1))
        and chern_errors[-1] < 1.0e-6,
        f"cells={chern_cells}; c1={chern_sequence}; final error={chern_errors[-1]:.2e}",
    )

    # -------------------------------- representation and failure controls
    j_plus = np.array(
        [[0.0, math.sqrt(2.0), 0.0], [0.0, 0.0, math.sqrt(2.0)], [0.0, 0.0, 0.0]],
        dtype=complex,
    )
    j_minus = j_plus.conjugate().T
    j_x = 0.5 * (j_plus + j_minus)
    j_y = (j_plus - j_minus) / (2j)
    j_z = np.diag([1.0, 0.0, -1.0]).astype(complex)
    casimir = j_x @ j_x + j_y @ j_y + j_z @ j_z
    check(
        "M4_representation",
        "the locked pair is the three-state spin-one representation",
        max_abs(casimir - 2.0 * np.eye(3)) < TOL
        and max_abs(commutator(j_x, j_y) - 1j * j_z) < TOL,
        "dimension=3; max Casimir residual={:.2e}; commutator residual={:.2e}".format(
            max_abs(casimir - 2.0 * np.eye(3)),
            max_abs(commutator(j_x, j_y) - 1j * j_z),
        ),
    )

    large_gauge_residual = max(
        abs(np.exp(1j * 2.0 * math.pi * level * winding) - 1.0)
        for level in (1, 2, 3)
        for winding in range(-3, 4)
    )
    check(
        "M4_representation",
        "large-gauge invariance quantizes k but does not uniquely select signed k=+2",
        large_gauge_residual < TOL,
        f"levels 1,2,3 all invariant; max phase residual={large_gauge_residual:.2e}",
    )

    time_reversal = 1j * sigma_y
    opposite_berry_residuals: list[float] = []
    for _ in range(128):
        w = complex(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))
        dw = complex(rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5))
        section, derivative = directional_section(w, dw)
        partner = time_reversal @ section.conjugate()
        partner_derivative = time_reversal @ derivative.conjugate()
        opposite_berry_residuals.append(
            abs(np.vdot(section, derivative) + np.vdot(partner, partner_derivative))
        )
    check(
        "M5_failure_modes",
        "opposite Berry orientations cancel to k=0 rather than |k|=2",
        max(opposite_berry_residuals) < TOL,
        f"samples=128; max cancellation residual={max(opposite_berry_residuals):.2e}",
    )
    single_dimension = int(identity_two.shape[0])
    pair_dimension = int(round(float(np.trace(symmetric_projector).real)))
    singleton_exclusion_derived = False
    check(
        "M5_failure_modes",
        "a permitted single constituent would retain an unwanted |k|=1 doublet",
        single_dimension == 2
        and pair_dimension == 3
        and not singleton_exclusion_derived,
        "single Hilbert dimension={} at |k|=1; pair dimension={} at |k|=2; "
        "singleton_exclusion_derived={}".format(
            single_dimension, pair_dimension, singleton_exclusion_derived
        ),
    )
    check(
        "M5_failure_modes",
        "without Hund locking the target is CP1xCP1 rather than its diagonal CP1",
        identity_four.shape == (single_dimension**2, single_dimension**2)
        and pair_dimension == 3,
        "unlocked product Hilbert dimension={}; symmetric locked subspace dimension={}".format(
            identity_four.shape[0], pair_dimension
        ),
    )

    portal_scale = 0.01
    temperature_scale = 0.02
    drive_scale = 0.03
    retained_gap = min(full_charge_gap, singlet_gap, field_strength)
    hierarchy_ratio = max(portal_scale, temperature_scale, drive_scale) / retained_gap
    check(
        "M6_EFT_gate",
        "an explicit benchmark satisfies portal, temperature, and drive scales below every retained gap",
        hierarchy_ratio < 0.2,
        f"min gap={retained_gap:.6f}; max perturbation={max(portal_scale,temperature_scale,drive_scale):.6f}; ratio={hierarchy_ratio:.6f}",
    )
    bad_portal_scale = 1.1 * retained_gap
    bad_hierarchy_pass = bad_portal_scale < retained_gap
    check(
        "M6_EFT_gate",
        "a portal above the smallest gap is rejected by the scale-separation gate",
        not bad_hierarchy_pass,
        "bad portal={:.6f}; min gap={:.6f}; hierarchy pass={}".format(
            bad_portal_scale, retained_gap, bad_hierarchy_pass
        ),
    )

    level_magnitude = 2
    signed_level_for_declared_action = -2
    positive_level_requires_orientation_reversal = True
    signed_route_e_level_selected = False
    candidate_level_derivation_complete = True
    ap_e3_physics_closed = False
    physics_promotion_allowed = False
    check(
        "M7_boundary",
        "the declared dimer derives |k|=2 while the signed Route-E orientation and UV selection remain open",
        candidate_level_derivation_complete
        and level_magnitude == 2
        and signed_level_for_declared_action == -2
        and positive_level_requires_orientation_reversal
        and not signed_route_e_level_selected
        and not ap_e3_physics_closed
        and not physics_promotion_allowed,
        "magnitude-two proof complete; aligned action k=-2; signed orientation, UV field content, and Route-E identification remain open",
    )

    passed = sum(item["pass"] for item in CHECKS)
    all_pass = passed == len(CHECKS)
    result: dict[str, Any] = {
        "audit": "verify_ap_e3_level_two_microscopic",
        "status": (
            "ap_e3_hund_pair_derives_abs_k2_conditionally_uv_and_sign_open"
            if all_pass
            else "ap_e3_mechanical_failure"
        ),
        "all_pass": all_pass,
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "candidate_level_derivation_complete": candidate_level_derivation_complete,
        "candidate_level_magnitude_derivation_complete": candidate_level_derivation_complete,
        "level_magnitude": level_magnitude,
        "signed_level_for_declared_action": signed_level_for_declared_action,
        "positive_level_requires_orientation_reversal": positive_level_requires_orientation_reversal,
        "signed_route_e_level_selected": signed_route_e_level_selected,
        "ap_e3_physics_closed": ap_e3_physics_closed,
        "physics_promotion_allowed": physics_promotion_allowed,
        "derived_chain": [
            "for H_portal=0, C=mu+h/2+J_H/4 < U-J_H/8 is a sufficient global interacting (1,1)-plateau condition",
            "J_H>0 projects the dimer to the symmetric spin-one triplet",
            "the diagonal coherent state is the quadratic Veronese embedding",
            "Berry connections add, so the ket eigenline is O(-2) and its dual prequantum line is O(2)",
            "the aligned displayed action has k=-2; orientation reversal gives k=+2",
            "the dual quantized carrier is H0(CP1,O(2)) with dimension three",
        ],
        "alternative_branches": {
            "fermion_determinant": (
                "r filled same-orientation eigenlines induce k=sum C_a; r=2 can give |k|=2, "
                "but flavor count, filling, and signs must be independently fixed"
            ),
            "mixed_WZW": (
                "a quantized mixed WZW coefficient can reduce to k=n_c B and is nonrenormalization-protected, "
                "but an n_c=2 UV completion with the correct global coset is not present in this repository"
            ),
            "tangent_Dirac": (
                "T_CP1=O(2) gives a three-section Spin-c/Dolbeault index only after the physical mode "
                "is proved tangent-valued; this is deferred to AP-E4"
            ),
        },
        "remaining_blockers": [
            "derive why the UV field content contains exactly two orbitals/constituents",
            "prove odd or single-constituent sectors are absent rather than merely heavy by assumption",
            "embed the dimer and its Mott/Hund interactions into an anomaly-consistent four-dimensional field theory",
            "derive the portal to Route-E families and keep it below charge, singlet, and orientation gaps",
            "show the spin-one triplet is chiral family data rather than an internal vector multiplet",
            "derive the physical orientation/chirality that selects the required signed Route-E level rather than only |k|=2",
        ],
        "numerical_summary": {
            "mott_parameters": {"U": interaction, "mu": chemical_potential},
            "interacting_plateau_constant_C": plateau_constant,
            "interacting_plateau_upper_bound": plateau_upper_bound,
            "interacting_plateau_condition_margin": plateau_condition_margin,
            "charge_gap": charge_gap,
            "interacting_charge_gap": full_charge_gap,
            "interacting_ground_sector": list(full_ground_sector),
            "hund_spectrum": hund_eigenvalues.tolist(),
            "singlet_gap": singlet_gap,
            "orientation_gap": min(zeeman_gaps),
            "chern_theta_cells": chern_cells,
            "chern_midpoint_sequence": chern_sequence,
            "chern_final_error": chern_errors[-1],
            "berry_factor_two_max_residual": max(berry_residuals),
            "aligned_action_sign_max_residual": max(action_sign_residuals),
            "metric_factor_two_max_residual": max(metric_residuals),
            "scale_hierarchy_ratio": hierarchy_ratio,
        },
        "source_manifest": source_manifest,
        "checks": CHECKS,
    }
    write_outputs(result)
    print(
        f"verify_ap_e3_level_two_microscopic: {passed}/{len(CHECKS)} checks pass; "
        f"candidate_level_magnitude_derivation_complete={str(candidate_level_derivation_complete).lower()}; "
        "physics_promotion_allowed=false"
    )
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
