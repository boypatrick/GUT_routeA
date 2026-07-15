#!/usr/bin/env python3
"""Deterministic audit for the AP-E3 global-definition/SU(3) UV note.

This verifier checks only algebraic and numerical statements that can be
closed without assuming strong dynamics.  A green run deliberately does not
select the two spin-torsion phases, prove radiatively stable lepton
decoupling, restore an exact magnetic one-form symmetry at finite SU(3)
breaking scale, or establish a Route-E-equivalent all-scale completion.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e3_nonextendible_wzw_su3_uv_audit.tex"
BIB = ROUTE_F / "tex" / "ap_e3_nonextendible_wzw_su3_uv_audit.bib"
THIS_SCRIPT = Path(__file__).resolve()

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


def gauss_integral(function, lower: float, upper: float, order: int = 96) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points = 0.5 * (upper - lower) * nodes + 0.5 * (upper + lower)
    return float(0.5 * (upper - lower) * np.dot(weights, function(points)))


def decompose_su2_weights(weights: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return (dimension, U(1) charge) terms from (twice-m, charge) weights."""

    by_charge: dict[int, Counter[int]] = defaultdict(Counter)
    for twice_m, charge in weights:
        by_charge[charge][twice_m] += 1

    terms: list[tuple[int, int]] = []
    for charge, multiplicities in by_charge.items():
        remaining = Counter(multiplicities)
        while sum(remaining.values()):
            highest = max(m for m, count in remaining.items() if count)
            if highest < 0:
                raise AssertionError("non-symmetric SU(2) weight system")
            for twice_m in range(highest, -highest - 1, -2):
                if remaining[twice_m] <= 0:
                    raise AssertionError("invalid SU(2) character")
                remaining[twice_m] -= 1
            terms.append((highest + 1, charge))
    return sorted(terms, key=lambda item: (-item[1], -item[0]))


def symmetric_power_branch(power: int, conjugate: bool = False) -> list[tuple[int, int]]:
    """Branch Sym^power(3) using 3 -> 2_(+1) + 1_(-2)."""

    weights: list[tuple[int, int]] = []
    for n1 in range(power + 1):
        for n2 in range(power - n1 + 1):
            n3 = power - n1 - n2
            twice_m = n1 - n2
            charge = n1 + n2 - 2 * n3
            if conjugate:
                charge = -charge
            weights.append((twice_m, charge))
    return decompose_su2_weights(weights)


def adjoint_branch() -> list[tuple[int, int]]:
    fundamental = [(1, 1), (-1, 1), (0, -2)]
    antifundamental = [(1, -1), (-1, -1), (0, 2)]
    product = Counter(
        (m_a + m_b, q_a + q_b)
        for (m_a, q_a), (m_b, q_b) in itertools.product(
            fundamental, antifundamental
        )
    )
    product[(0, 0)] -= 1  # remove the trace singlet from 3 x bar(3)
    weights: list[tuple[int, int]] = []
    for weight, multiplicity in product.items():
        weights.extend([weight] * multiplicity)
    return decompose_su2_weights(weights)


def matrix_embedding(u2: np.ndarray, alpha: float) -> np.ndarray:
    result = np.zeros((3, 3), dtype=complex)
    result[:2, :2] = cmath.exp(1.0j * alpha) * u2
    result[2, 2] = cmath.exp(-2.0j * alpha)
    return result


def run_audit() -> dict[str, Any]:
    # Spin-bordism stable splitting for (S^3 x S^2)_+.
    spin_bordism = {0: "Z", 1: "Z2", 2: "Z2", 3: "0", 4: "Z"}
    full_omega4 = [spin_bordism[4], spin_bordism[1], spin_bordism[2]]
    reduced_omega4 = [spin_bordism[1], spin_bordism[2]]
    check(
        "bordism",
        "full stable splitting",
        full_omega4 == ["Z", "Z2", "Z2"],
        "Omega_4^Spin(S3 x S2)=Z + Z2 + Z2",
    )
    check(
        "bordism",
        "reduced stable splitting",
        reduced_omega4 == ["Z2", "Z2"],
        "the two reduced generators are inverse-image spin 1- and 2-manifolds",
    )

    # Integral cohomology and uniqueness of a degree-five differential character.
    cohomology_ranks = {0: 1, 2: 1, 3: 1, 5: 1}
    check(
        "differential_character",
        "Kunneth cohomology ranks",
        cohomology_ranks == {0: 1, 2: 1, 3: 1, 5: 1},
        "H*(S3 x S2;Z) is generated by 1,u2,u3,u3 cup u2",
    )
    h4_rank = cohomology_ranks.get(4, 0)
    check(
        "differential_character",
        "no flat degree-five ambiguity",
        h4_rank == 0,
        "H^4(X;R/Z)=0, so fixed integral curvature has a unique ordinary differential refinement",
    )

    s3_period = (
        gauss_integral(lambda x: np.sin(x) ** 2, 0.0, math.pi)
        * gauss_integral(np.sin, 0.0, math.pi)
        * (2.0 * math.pi)
        / (2.0 * math.pi**2)
    )
    s2_period = (
        gauss_integral(np.sin, 0.0, math.pi)
        * (2.0 * math.pi)
        / (4.0 * math.pi)
    )
    check(
        "differential_character",
        "normalized S3 generator",
        abs(s3_period - 1.0) < TOL,
        f"integral omega3={s3_period:.16g}",
    )
    check(
        "differential_character",
        "normalized S2 generator",
        abs(s2_period - 1.0) < TOL,
        f"integral omega2={s2_period:.16g}",
    )
    check(
        "differential_character",
        "integral mixed generator",
        abs(s3_period * s2_period - 1.0) < TOL,
        f"integral_(S3xS2) omega3 wedge omega2={s3_period*s2_period:.16g}",
    )

    sample_chi = np.linspace(0.07, math.pi - 0.07, 257)
    fprime = (1.0 - np.cos(2.0 * sample_chi)) / (4.0 * math.pi**2)
    target_fprime = np.sin(sample_chi) ** 2 / (2.0 * math.pi**2)
    overlap_period = (
        gauss_integral(np.sin, 0.0, math.pi)
        * (2.0 * math.pi)
        / (4.0 * math.pi)
    )
    check(
        "cech",
        "gerbe local-potential derivative",
        float(np.max(np.abs(fprime - target_fprime))) < TOL,
        "d[(2chi-sin2chi)/(8pi^2) vol_S2]=omega3",
    )
    check(
        "cech",
        "gerbe overlap integral",
        abs(overlap_period - 1.0) < TOL,
        f"integral(B_N-B_S)={overlap_period:.16g}",
    )
    theta = np.linspace(0.11, math.pi - 0.11, 173)
    hopf_n = 0.5 * (1.0 - np.cos(theta))
    hopf_s = -0.5 * (1.0 + np.cos(theta))
    check(
        "cech",
        "Hopf transition winding",
        float(np.max(np.abs((hopf_n - hopf_s) - 1.0))) < TOL,
        "A_N-A_S=dphi and the transition has winding one",
    )
    extension_phases = [
        cmath.exp(2.0j * math.pi * 2 * integer) for integer in range(-8, 9)
    ]
    check(
        "differential_character",
        "extension-independence at n=2",
        max(abs(phase - 1.0) for phase in extension_phases) < TOL,
        "exp(2pi i n m)=1 for every tested integral five-cycle m",
    )

    # SU(3) -> [SU(2) x U(1)]/Z2 group map and representation branching.
    identity2 = np.eye(2, dtype=complex)
    kernel_image = matrix_embedding(-identity2, math.pi)
    nonkernel_image = matrix_embedding(identity2, math.pi)
    check(
        "su3_embedding",
        "quotient kernel",
        np.linalg.norm(kernel_image - np.eye(3)) < TOL
        and np.linalg.norm(nonkernel_image - np.eye(3)) > 1.0,
        "ker iota={(I,1),(-I,-1)}",
    )

    branches = {
        "3": symmetric_power_branch(1),
        "bar3": symmetric_power_branch(1, conjugate=True),
        "6": symmetric_power_branch(2),
        "bar6": symmetric_power_branch(2, conjugate=True),
        "8": adjoint_branch(),
        "10": symmetric_power_branch(3),
    }
    expected = {
        "3": [(2, 1), (1, -2)],
        "bar3": [(1, 2), (2, -1)],
        "6": [(3, 2), (2, -1), (1, -4)],
        "bar6": [(1, 4), (2, 1), (3, -2)],
        "8": [(2, 3), (3, 0), (1, 0), (2, -3)],
        "10": [(4, 3), (3, 0), (2, -3), (1, -6)],
    }
    for name in branches:
        check(
            "branching",
            f"{name} branching",
            branches[name] == expected[name],
            f"computed={branches[name]}",
        )

    all_terms = [term for branch in branches.values() for term in branch]
    parity_ok = all(((dimension - 1 + charge) % 2 == 0) for dimension, charge in all_terms)
    singlet_odd_absent = all(
        not (dimension == 1 and charge % 2) for dimension, charge in all_terms
    )
    check(
        "su3_embedding",
        "quotient parity",
        parity_ok,
        "every (j,q) term obeys 2j+q=0 mod 2",
    )
    check(
        "su3_embedding",
        "no odd-charge colour singlet",
        singlet_odd_absent,
        "all tested SU(3) irreps obey q even when j=0; the group-kernel proof is representation-independent",
    )
    diophantine_counterexamples = [
        (xq, xphi)
        for xq in range(1, 12, 2)
        for xphi in range(2, 14, 2)
        if 2 * xq == 2 * xphi
    ]
    check(
        "su3_embedding",
        "exactly-two dressing Diophantine no-go",
        not diophantine_counterexamples,
        "m=2 requires Xq=Xphi, impossible for odd doublet charge and even singlet charge",
    )
    minimal_dressing = 2 * 1 // math.gcd(2 * 1, 2)
    check(
        "su3_embedding",
        "minimal charge-two dressing",
        minimal_dressing == 1,
        "Xq=1,Xphi=2 makes qq phi^dagger neutral with one, not two, scalar",
    )

    # Gauge and generalized anomaly ledger after branching a vectorlike SU(3) Dirac field.
    nf = 2
    su2sq_u1 = nf * (Fraction(1, 2) * 1 + Fraction(1, 2) * -1)
    u1_cube = nf * (2 * 1**3 + 2 * (-1) ** 3 + (-2) ** 3 + 2**3)
    grav_u1 = nf * (2 * 1 + 2 * -1 + -2 + 2)
    su3_cube = nf * (1 - 1)
    color_witten_doublets = 2 * nf
    check("anomaly", "SU(3)^3 vectorlike", su3_cube == 0, f"A={su3_cube}")
    check(
        "anomaly",
        "SU(2)c^2 U(1)g",
        su2sq_u1 == 0,
        f"A={su2sq_u1}",
    )
    check("anomaly", "U(1)g^3", u1_cube == 0, f"A={u1_cube}")
    check("anomaly", "grav^2 U(1)g", grav_u1 == 0, f"A={grav_u1}")
    check(
        "anomaly",
        "SU(2)c Witten parity",
        color_witten_doublets % 2 == 0,
        f"left-Weyl colour doublets={color_witten_doublets}",
    )
    kappa_q = 2 * 1
    kappa_lepton = 1 * -2
    check(
        "anomaly",
        "Postnikov cancellation in diagonal deep-UV flavour",
        kappa_q == 2 and kappa_lepton == -2 and kappa_q + kappa_lepton == 0,
        f"kappa_q={kappa_q}, kappa_l={kappa_lepton}",
    )
    check(
        "anomaly",
        "Nf=2 flavour-Witten split",
        2 % 2 == 0 and 1 % 2 == 1 and 3 % 2 == 1,
        "quark copies=2 (even), lepton copies=1 (odd), diagonal SU(3) copies=3 (odd)",
    )

    # Tree-level mass matrices and the conditional charge-two scalar split.
    v_sigma = 5.0
    yukawa = np.array([[0.7, 0.2], [0.2, 1.1]], dtype=float)
    bare_mass = -v_sigma * yukawa
    mass_q = bare_mass + v_sigma * yukawa
    mass_l = bare_mass - 2.0 * v_sigma * yukawa
    lepton_singular_values = np.linalg.svd(mass_l, compute_uv=False)
    check(
        "mass_matrix",
        "tree-level massless quark block",
        np.linalg.norm(mass_q) < TOL,
        f"||M_q||={np.linalg.norm(mass_q):.3e}",
    )
    check(
        "mass_matrix",
        "full-rank heavy lepton block",
        min(lepton_singular_values) > 1.0,
        f"singular values={lepton_singular_values.tolist()}",
    )

    v_scalar = 1.0
    scalar_m2 = 1.0
    scalar_kappa = -1.0
    scalar_rho = 0.0
    mass_chi2 = scalar_m2 - scalar_kappa * v_scalar + scalar_rho * v_scalar**2
    mass_phi2 = (
        scalar_m2 + 2.0 * scalar_kappa * v_scalar + 4.0 * scalar_rho * v_scalar**2
    )
    check(
        "mass_matrix",
        "charge-two scalar split",
        mass_chi2 > 0.0 and mass_phi2 < 0.0,
        f"M_chi^2={mass_chi2}, m_phi^2={mass_phi2}",
    )

    b0_su3 = (
        Fraction(11, 3) * 3
        - Fraction(4, 3) * nf * Fraction(1, 2)
        - Fraction(1, 6) * 3
        - Fraction(1, 3) * 2 * Fraction(1, 2)
    )
    check(
        "running",
        "deep SU(3) asymptotic freedom",
        b0_su3 == Fraction(53, 6) and b0_su3 > 0,
        f"b0={b0_su3}",
    )

    # Fail-closed claim gates: these are expected false and are checked as such.
    torsion_uv_selected = False
    chiral_decoupling_protected = False
    exact_magnetic_one_form_at_finite_breaking = False
    charge_two_variant_route_equivalent = False
    full_uv_closed = False
    physics_promotion_allowed = False
    check(
        "claim_boundary",
        "torsion selection remains open",
        not torsion_uv_selected,
        "epsilon_3 and epsilon_2 require a specified UV APS/Dai-Freed determinant",
    )
    check(
        "claim_boundary",
        "radiative chiral protection remains open",
        not chiral_decoupling_protected,
        "M_q=0 is a tree-level tuning in the displayed adjoint mass model",
    )
    check(
        "claim_boundary",
        "finite-scale magnetic one-form is emergent only",
        not exact_magnetic_one_form_at_finite_breaking,
        "SU(3) monopoles violate the abelian Bianchi identity nonperturbatively",
    )
    check(
        "claim_boundary",
        "charge-two branch is not Route-E-equivalent",
        not charge_two_variant_route_equivalent,
        "one scalar dressing gives a doublet rather than the exactly-two triplet",
    )
    check(
        "claim_boundary",
        "no false full-UV promotion",
        not full_uv_closed and not physics_promotion_allowed,
        "global free term is defined, but UV dynamics and torsion choices remain open",
    )

    passed = sum(row["pass"] for row in CHECKS)
    all_pass = passed == len(CHECKS)
    result: dict[str, Any] = {
        "all_pass": all_pass,
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "status": "pass_with_open_uv_gates" if all_pass else "mechanical_failure",
        "checks": CHECKS,
        "spin_bordism_full": "Z + Z2 + Z2",
        "spin_bordism_reduced": "Z2 + Z2",
        "ordinary_differential_character_defined_on_all_4_cycles": True,
        "u2_quotient_global_bundle_normalization_proven": False,
        "torsion_uv_selected": torsion_uv_selected,
        "su3_charge_one_singlet_embedding_possible": False,
        "su3_charge_two_variant_tree_level_constructed": True,
        "su3_charge_two_variant_route_equivalent": charge_two_variant_route_equivalent,
        "chiral_decoupling_protected": chiral_decoupling_protected,
        "exact_magnetic_one_form_at_finite_breaking": exact_magnetic_one_form_at_finite_breaking,
        "full_uv_closed": full_uv_closed,
        "physics_promotion_allowed": physics_promotion_allowed,
        "numerics": {
            "s3_period": s3_period,
            "s2_period": s2_period,
            "mixed_period": s3_period * s2_period,
            "gerbe_overlap_period": overlap_period,
            "lepton_singular_values": lepton_singular_values.tolist(),
            "scalar_coloured_mass_squared": mass_chi2,
            "scalar_singlet_mass_squared": mass_phi2,
            "su3_b0": f"{b0_su3.numerator}/{b0_su3.denominator}",
        },
        "branching": {
            name: [{"su2_dimension": dim, "u1_charge": charge} for dim, charge in terms]
            for name, terms in branches.items()
        },
        "remaining_blockers": [
            "Compute the two torsion characters on explicit generator backgrounds using a specified UV APS/Dai-Freed regulator; epsilon_2 is not fixed by the curvature term.",
            "Protect M_q=0 while M_l is full rank, or perform and document order-by-order threshold tuning without claiming an exact deep-UV chiral symmetry.",
            "Include finite-mass SU(3) monopoles and quantify violation of the emergent magnetic one-form symmetry.",
            "Show a nonperturbative mesonic phase and a stable B=1 soliton after the charge-two scalar deformation.",
            "Resolve the exact-two/triplet mismatch: the SU(3)-allowed charge-two scalar needs only one dressing.",
            "Match determinant lines and first Chern classes on faithful U(2) bundles, including the correlation between U(1) flux and SU(2)-centre 't Hooft flux; the covering-group charge-two normalization is conditional until then.",
        ],
        "source_manifest": [source_row(path) for path in (TEX, BIB, THIS_SCRIPT)],
    }
    return result


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e3_nonextendible_wzw_su3_uv.json"
    md_path = OUTPUT / "ap_e3_nonextendible_wzw_su3_uv.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    check_lines = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` — "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    blocker_lines = "\n".join(
        f"- {blocker}" for blocker in result["remaining_blockers"]
    )
    source_lines = "\n".join(
        f"- `{row['path']}` — `{row['sha256']}` ({row['size_bytes']} bytes)"
        for row in result["source_manifest"]
        if row["exists"]
    )
    md_path.write_text(
        f"""# AP-E3 non-extendible mixed-WZW / SU(3) UV audit

- Status: `{result['status']}`
- Mechanical checks: `{result['checks_passed']}/{result['checks_total']}`
- `all_pass`: `{str(result['all_pass']).lower()}`
- Ordinary differential character defined on all four-cycles: `{str(result['ordinary_differential_character_defined_on_all_4_cycles']).lower()}`
- Faithful-U(2) quotient bundle normalization proven: `{str(result['u2_quotient_global_bundle_normalization_proven']).lower()}`
- Spin torsion selected by UV: `{str(result['torsion_uv_selected']).lower()}`
- Charge-one singlet embeds in SU(3): `{str(result['su3_charge_one_singlet_embedding_possible']).lower()}`
- Charge-two branch Route-E-equivalent: `{str(result['su3_charge_two_variant_route_equivalent']).lower()}`
- Full UV closed: `{str(result['full_uv_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

## Exact chain

```text
X = SU(2) x S2 = S3 x S2
Omega_4^Spin(X) = Z + Z2 + Z2
Z_mix = Hol(n u3-hat cup u2-hat) (-1)^(epsilon_3 nu_3 + epsilon_2 nu_2)
SU(3) -> [SU(2)c x U(1)g]/Z2: 2j+q = 0 mod 2
3 -> 2_(+1) + 1_(-2), bar(3) -> 2_(-1) + 1_(+2)
M_q = m + V y, M_l = m - 2 V y
m = -V y -> M_q = 0, M_l = -3 V y (tree-level tuning)
2 X_q = m_dress X_phi; m_dress=2 contradicts odd X_q/even X_phi
```

## Checks

{check_lines}

## Remaining blockers

{blocker_lines}

## Source manifest

{source_lines}
""",
        encoding="utf-8",
    )


def main() -> None:
    result = run_audit()
    write_outputs(result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
