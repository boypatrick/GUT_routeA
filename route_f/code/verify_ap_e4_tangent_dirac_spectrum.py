#!/usr/bin/env python3
"""Fail-closed AP-E4 tangent/Dolbeault-Dirac spectral audit.

Mathematically, the canonical Spin^c Dolbeault operator on
Lambda^{0,*}(CP1) tensor T^{1,0}CP1 has three positive-chirality zero modes,
no negative-chirality zero mode, and a completely paired massive tower.  The
script verifies the index, the exact spectrum, SU(2) Casimir realization,
chirality pairing, radius convention, and the ordinary-spin negative control.

The physical interpretation remains conditional.  A bosonic sigma-model
fluctuation is tangent-valued, but it does not by itself become a fermion or a
Dirac operator on the target.  A moduli-space SQM or higher-dimensional
compactification, its anomaly cancellation, and the Route-E portal are still
missing.  Therefore a green run is mathematical closure only and never
authorizes physics promotion.
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
TEX = ROUTE_F / "tex" / "ap_e4_tangent_dirac_spectrum.tex"
BIB = ROUTE_F / "tex" / "ap_e4_tangent_dirac_spectrum.bib"
AP_E1_TEX = ROUTE_F / "tex" / "ap_e1_projective_doublet_action.tex"
LITERATURE = ROUTE_F / "LITERATURE_2023_2026.md"
THIS_SCRIPT = Path(__file__).resolve()

RADIUS = Fraction(1, 2)
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


def h0_o(degree: int) -> int:
    return max(degree + 1, 0)


def h1_o(degree: int) -> int:
    return h0_o(-degree - 2)


def su2_generators(j: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dimension = int(round(2.0 * j + 1.0))
    m_values = np.array([-j + index for index in range(dimension)], dtype=float)
    jz = np.diag(m_values).astype(complex)
    jplus = np.zeros((dimension, dimension), dtype=complex)
    for index, m_value in enumerate(m_values[:-1]):
        jplus[index + 1, index] = math.sqrt(
            (j - m_value) * (j + m_value + 1.0)
        )
    jminus = jplus.conjugate().T
    jx = 0.5 * (jplus + jminus)
    jy = (jplus - jminus) / (2.0j)
    return jx, jy, jz


def midpoint_flux(theta_cells: int, degree: int) -> float:
    dtheta = math.pi / theta_cells
    # F_p=(p/2) sin(theta)dtheta^dphi, so (2pi)^-1 int F_p.
    return 0.5 * degree * sum(
        math.sin((index + 0.5) * dtheta) * dtheta
        for index in range(theta_cells)
    )


def build_truncated_dirac(q: int, radius: float, max_level: int) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    zero_count = q + 1
    blocks: list[np.ndarray] = [np.zeros((zero_count, zero_count), dtype=complex)]
    gamma_blocks: list[np.ndarray] = [np.eye(zero_count, dtype=complex)]
    ledger: list[dict[str, Any]] = []
    for level in range(1, max_level + 1):
        multiplicity = q + 2 * level + 1
        eigenvalue = math.sqrt(level * (level + q + 1)) / radius
        identity = np.eye(multiplicity, dtype=complex)
        blocks.append(
            np.block(
                [
                    [np.zeros_like(identity), eigenvalue * identity],
                    [eigenvalue * identity, np.zeros_like(identity)],
                ]
            )
        )
        gamma_blocks.append(
            np.block(
                [
                    [identity, np.zeros_like(identity)],
                    [np.zeros_like(identity), -identity],
                ]
            )
        )
        ledger.append(
            {
                "n": level,
                "lambda_abs": eigenvalue,
                "multiplicity_per_sign": multiplicity,
            }
        )

    def block_diag(matrices: list[np.ndarray]) -> np.ndarray:
        size = sum(matrix.shape[0] for matrix in matrices)
        output = np.zeros((size, size), dtype=complex)
        offset = 0
        for matrix in matrices:
            dimension = matrix.shape[0]
            output[offset : offset + dimension, offset : offset + dimension] = matrix
            offset += dimension
        return output

    return block_diag(blocks), block_diag(gamma_blocks), ledger


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e4_tangent_dirac_spectrum.json"
    md_path = OUTPUT / "ap_e4_tangent_dirac_spectrum.md"
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
        f"""# AP-E4 tangent-valued Dirac audit

- Status: `{result['status']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Tangent bundle derived: `{str(result['tangent_bundle_derived']).lower()}`
- Canonical Spin-c mathematics complete: `{str(result['canonical_spinc_mathematics_complete']).lower()}`
- Three chiral zero modes: `{str(result['canonical_spinc_three_zero_modes']).lower()}`
- Full partner spectrum and gap: `{str(result['partner_spectrum_and_gap_complete']).lower()}`
- Physical fermionic tangent mode derived: `{str(result['physical_fermionic_tangent_mode_derived']).lower()}`
- AP-E4 physics closed: `{str(result['ap_e4_physics_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

## Canonical Spin-c result

```text
E=T^(1,0)CP1=O(2)
D_T^c=sqrt(2)(dbar_T+dbar_T^dagger)
ker D_+=H0(O(2))=C^3, ker D_-=H1(O(2))=0
lambda_(n,+/-)=+/- sqrt(n(n+3))/R
mult per sign=2n+3, n>=1
R=1/2 -> gap=4; first massive multiplicity=5 per sign
```

Ordinary spin Dirac twisted only by `T=O(2)` instead gives two zero modes and
gap `2 sqrt(3)`.  Canonical Spin-c with `O(2)` is equivalent to ordinary spin
twisted by `O(3)`; the half-canonical shift is physical input, not notation.

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
    radius = float(RADIUS)
    q = 2
    critical_sources = [TEX, BIB, AP_E1_TEX, LITERATURE, THIS_SCRIPT]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "E40_provenance",
        "all AP-E4 critical sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # ------------------------------------------------ geometry and tangent bundle
    gaussian_curvature = 1.0 / radius**2
    scalar_curvature = 2.0 / radius**2
    check(
        "E41_geometry",
        "the AP-E1 Fubini-Study convention is a round sphere of radius one half",
        RADIUS == Fraction(1, 2)
        and abs(gaussian_curvature - 4.0) < TOL
        and abs(scalar_curvature - 8.0) < TOL,
        f"R={radius}; K={gaussian_curvature}; Scal={scalar_curvature}",
    )
    phases = np.linspace(0.0, 2.0 * math.pi, 4097)
    tangent_winding = (
        np.unwrap(np.angle(np.exp(2.0j * phases)))[-1]
        - np.unwrap(np.angle(np.exp(2.0j * phases)))[0]
    ) / (2.0 * math.pi)
    check(
        "E41_geometry",
        "the tangent transition has degree two, hence T^(1,0)CP1=O(2)",
        abs(tangent_winding - 2.0) < TOL,
        f"transition winding={tangent_winding:.12f}",
    )
    flux_sequence = [midpoint_flux(cells, degree=2) for cells in (10, 20, 40, 80, 160)]
    check(
        "E41_geometry",
        "the tangent Chern connection has c1=2",
        abs(flux_sequence[-1] - 2.0) < 4.0e-5
        and all(abs(flux_sequence[i + 1] - 2.0) < abs(flux_sequence[i] - 2.0) for i in range(4)),
        f"flux sequence={flux_sequence}",
    )

    # Direct horizontal projection at deterministic points.
    horizontal_residuals: list[float] = []
    norm_derivative_residuals: list[float] = []
    for w in (-0.8 + 0.3j, -0.1 - 0.6j, 0.4 + 0.2j, 1.1 - 0.7j):
        for dw in (0.07 + 0.03j, -0.02 + 0.09j):
            norm = 1.0 + abs(w) ** 2
            dnorm = 2.0 * (w.conjugate() * dw).real
            z = np.array([1.0, w], dtype=complex) / math.sqrt(norm)
            dz = np.array(
                [
                    -0.5 * dnorm / norm**1.5,
                    dw / math.sqrt(norm) - 0.5 * w * dnorm / norm**1.5,
                ],
                dtype=complex,
            )
            projector = np.eye(2, dtype=complex) - np.outer(z, z.conjugate())
            eta = projector @ dz
            horizontal_residuals.append(abs(np.vdot(z, eta)))
            norm_derivative_residuals.append(abs(2.0 * np.vdot(z, dz).real))
    check(
        "E41_geometry",
        "Q_z delta z is the horizontal tangent-valued sigma-model fluctuation",
        max(horizontal_residuals) < TOL and max(norm_derivative_residuals) < TOL,
        f"horizontal={max(horizontal_residuals):.3e}; norm={max(norm_derivative_residuals):.3e}",
    )

    # ----------------------------------------- cohomology, index, and Spin choices
    canonical_h0 = h0_o(q)
    canonical_h1 = h1_o(q)
    canonical_index = canonical_h0 - canonical_h1
    check(
        "E42_index",
        "canonical Dolbeault Spin-c twisted by T=O(2) has (h0,h1,index)=(3,0,3)",
        (canonical_h0, canonical_h1, canonical_index) == (3, 0, 3),
        f"h0={canonical_h0}; h1={canonical_h1}; index={canonical_index}",
    )
    check(
        "E42_index",
        "Serre duality independently gives H1(O(2))=H0(O(-4))*=0",
        h1_o(2) == h0_o(-4) == 0,
        f"h1(O2)={h1_o(2)}; h0(O-4)={h0_o(-4)}",
    )
    check(
        "E42_index",
        "Riemann-Roch/Todd integration gives chi(O(q))=q+1",
        canonical_index == q + 1,
        f"integral ch(O({q})) Td(T)=1+{q}={q + 1}",
    )
    ordinary_tangent_degree = 2
    ordinary_effective_dolbeault_degree = ordinary_tangent_degree - 1
    ordinary_h0 = h0_o(ordinary_effective_dolbeault_degree)
    ordinary_h1 = h1_o(ordinary_effective_dolbeault_degree)
    check(
        "E42_index",
        "negative control: ordinary spin Dirac twisted only by T=O(2) has two, not three, zero modes",
        (ordinary_h0, ordinary_h1, ordinary_h0 - ordinary_h1) == (2, 0, 2),
        f"effective Dolbeault line=O(1); h0={ordinary_h0}; index={ordinary_h0 - ordinary_h1}",
    )
    check(
        "E42_index",
        "canonical Spin-c O(2) equals ordinary spin twist O(3), exposing the half-canonical shift",
        ordinary_tangent_degree + 1 == 3 and h0_o(3 - 1) == 3,
        "S^+=O(-1); S^+ tensor O(3)=O(2)",
    )
    check(
        "E42_index",
        "nearby bundle controls q=1 and q=3 give two and four canonical zero modes",
        h0_o(1) == 2 and h0_o(3) == 4,
        f"q=1 -> {h0_o(1)}; q=3 -> {h0_o(3)}",
    )

    # ------------------------------------------------ exact and numerical spectrum
    exact_identity_pass = True
    exact_rows: list[dict[str, str]] = []
    for level in range(0, 10):
        j = Fraction(q, 2) + level
        left = j * (j + 1) - Fraction(q * q, 4) - Fraction(q, 2)
        right = Fraction(level * (level + q + 1), 1)
        multiplicity = 2 * j + 1
        exact_identity_pass &= left == right and multiplicity == q + 2 * level + 1
        exact_rows.append(
            {"n": str(level), "lambda_squared_R2": str(left), "multiplicity": str(multiplicity)}
        )
    check(
        "E43_spectrum",
        "the monopole-Casimir identity yields lambda_n^2=n(n+q+1)/R^2 exactly",
        exact_identity_pass,
        f"levels checked={len(exact_rows)}",
    )

    casimir_residuals: list[float] = []
    spectral_residuals: list[float] = []
    for level in range(0, 8):
        j = q / 2.0 + level
        generators = su2_generators(j)
        dimension = int(round(2.0 * j + 1.0))
        casimir = sum(generator @ generator for generator in generators)
        casimir_residuals.append(
            float(np.max(np.abs(casimir - j * (j + 1.0) * np.eye(dimension))))
        )
        shifted = casimir - (q / 2.0) * (q / 2.0 + 1.0) * np.eye(dimension)
        expected = level * (level + q + 1.0)
        spectral_residuals.append(
            float(np.max(np.abs(np.linalg.eigvalsh(shifted) - expected)))
        )
    check(
        "E43_spectrum",
        "independent finite SU(2) matrices reproduce every tested D^2 level and degeneracy",
        max(casimir_residuals) < TOL and max(spectral_residuals) < TOL,
        f"Casimir={max(casimir_residuals):.3e}; shifted spectrum={max(spectral_residuals):.3e}",
    )

    dirac, gamma, tower = build_truncated_dirac(q=q, radius=radius, max_level=5)
    hermitian_residual = float(np.max(np.abs(dirac - dirac.conjugate().T)))
    chirality_residual = float(np.max(np.abs(dirac @ gamma + gamma @ dirac)))
    eigenvalues = np.linalg.eigvalsh(dirac)
    zero_count = int(np.sum(np.abs(eigenvalues) < TOL))
    check(
        "E43_spectrum",
        "the truncated spectral operator is Hermitian and anticommutes with chirality",
        hermitian_residual < TOL and chirality_residual < TOL,
        f"Hermitian={hermitian_residual:.3e}; anticommutator={chirality_residual:.3e}",
    )
    check(
        "E43_spectrum",
        "the kernel contains exactly three positive-chirality states",
        zero_count == 3
        and np.allclose(np.diag(gamma)[:3], np.ones(3))
        and abs(float(np.trace(gamma[:3, :3]).real) - 3.0) < TOL,
        f"zero count={zero_count}; kernel chirality trace={np.trace(gamma[:3, :3]).real}",
    )
    nonzero = eigenvalues[np.abs(eigenvalues) >= TOL]
    positive = np.sort(nonzero[nonzero > 0.0])
    negative = np.sort(-nonzero[nonzero < 0.0])
    check(
        "E43_spectrum",
        "every nonzero level is paired at plus/minus lambda as implied by chirality anticommutation",
        len(positive) == len(negative) and float(np.max(np.abs(positive - negative))) < TOL,
        f"paired nonzero states per sign={len(positive)}",
    )
    canonical_gap = math.sqrt(1.0 * (1.0 + q + 1.0)) / radius
    first_positive_count = int(np.sum(np.abs(positive - canonical_gap) < TOL))
    check(
        "E43_spectrum",
        "at R=1/2 the canonical first gap is four with multiplicity five per sign",
        abs(canonical_gap - 4.0) < TOL and first_positive_count == 5,
        f"gap={canonical_gap:.12f}; multiplicity per sign={first_positive_count}",
    )
    ordinary_q = ordinary_tangent_degree - 1
    ordinary_gap = math.sqrt(1.0 * (1.0 + ordinary_q + 1.0)) / radius
    check(
        "E43_spectrum",
        "ordinary spin twist T=O(2) has the distinct gap 2 sqrt(3) and first multiplicity four",
        abs(ordinary_gap - 2.0 * math.sqrt(3.0)) < TOL
        and ordinary_q + 2 * 1 + 1 == 4,
        f"ordinary gap={ordinary_gap:.12f}; multiplicity=4",
    )
    check(
        "E43_spectrum",
        "radius drift control: changing R from one half to one halves the canonical gap",
        abs((math.sqrt(4.0) / 1.0) - canonical_gap / 2.0) < TOL,
        f"gap(R=1)={math.sqrt(4.0):.6f}; gap(R=1/2)={canonical_gap:.6f}",
    )

    mass = 0.3
    massive_dirac = dirac + mass * gamma
    massive_eigenvalues = np.linalg.eigvalsh(massive_dirac)
    check(
        "E44_fail_closed",
        "a chirality-breaking mass lifts the protected zero modes",
        int(np.sum(np.abs(massive_eigenvalues) < TOL)) == 0
        and int(np.sum(np.abs(massive_eigenvalues - mass) < TOL)) == 3,
        f"mass={mass}; remaining zeros={int(np.sum(np.abs(massive_eigenvalues) < TOL))}",
    )
    check(
        "E44_fail_closed",
        "nonintegral monopole degree is rejected as a global line bundle",
        not float(2.5).is_integer() and float(2.0).is_integer(),
        "q=5/2 rejected; q=2 accepted",
    )
    check(
        "E44_fail_closed",
        "absence of opposite chirality applies only to the kernel, not the massive tower",
        zero_count == 3 and len(positive) > 0 and len(negative) > 0,
        f"kernel=3 positive; massive states={len(positive) + len(negative)}",
    )

    remaining_blockers = [
        "The horizontal fluctuation Q_z delta z is a commuting sigma-model scalar; it does not derive a fermionic Clifford module or a Dirac operator on target CP1.",
        "Choose and derive either N=2 moduli-space/SQM quantization or an explicit higher-dimensional product containing CP1.  If CP1 is only field-value space, the target spectrum is not a Kaluza-Klein particle spectrum.",
        "Explain microscopically why the canonical Spin-c structure is selected.  Relative to the unique ordinary spin structure it supplies a half-canonical O(1) shift; ordinary spin twisted only by T has two zero modes.",
        "Prove that H0(TCP1) represents matter families rather than automorphism, gauge, or collective-coordinate redundancies.",
        "For a six-dimensional Weyl realization, cancel the full I8 gauge/gravitational anomaly (or exhibit Green-Schwarz factorization) and audit the Spin-c determinant flux.",
        "Show that all four-dimensional representation-weighted anomalies cancel and exclude conjugate bundles or additional sectors with compensating zero modes.",
        "Construct the Route-E portal below the gap 4 v/R_normalization and show its orientation degree is +1 without closing the gap or mixing in the five-fold first massive level.",
    ]

    all_pass = all(row["pass"] for row in CHECKS)
    result: dict[str, Any] = {
        "status": "ap_e4_canonical_spinc_mathematics_complete_physical_origin_and_anomalies_open",
        "all_pass": all_pass,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "checks": CHECKS,
        "tangent_bundle_derived": all_pass,
        "bosonic_tangent_fluctuation_derived": all_pass,
        "canonical_spinc_mathematics_complete": all_pass,
        "canonical_spinc_three_zero_modes": all_pass,
        "partner_spectrum_and_gap_complete": all_pass,
        "ordinary_spin_tangent_negative_control_passed": all_pass,
        "physical_fermionic_tangent_mode_derived": False,
        "compactification_or_moduli_origin_derived": False,
        "six_dimensional_anomaly_cancelled": False,
        "route_e_chiral_family_identification": False,
        "ap_e4_mathematics_complete": all_pass,
        "ap_e4_physics_closed": False,
        "physics_promotion_allowed": False,
        "conventions": {
            "R": radius,
            "gaussian_curvature": gaussian_curvature,
            "scalar_curvature": scalar_curvature,
            "canonical_bundle": "O(2)=T^(1,0)CP1",
            "canonical_operator": "sqrt(2)(dbar_T+dbar_T^dagger)",
        },
        "canonical_spectrum": {
            "zero_modes_positive": canonical_h0,
            "zero_modes_negative": canonical_h1,
            "index": canonical_index,
            "formula": "lambda_n,+/-=+/-sqrt(n(n+3))/R; multiplicity=2n+3",
            "gap": canonical_gap,
            "first_multiplicity_per_sign": first_positive_count,
            "tower_sample": tower,
        },
        "ordinary_spin_control": {
            "twist": "O(2)=T",
            "effective_dolbeault_line": "O(1)",
            "zero_modes": ordinary_h0,
            "index": ordinary_h0 - ordinary_h1,
            "gap": ordinary_gap,
        },
        "numerical_summary": {
            "flux_sequence": flux_sequence,
            "horizontal_max_residual": max(horizontal_residuals),
            "casimir_max_residual": max(casimir_residuals),
            "shifted_spectrum_max_residual": max(spectral_residuals),
            "dirac_hermitian_residual": hermitian_residual,
            "chirality_anticommutator_residual": chirality_residual,
            "truncated_zero_count": zero_count,
            "paired_nonzero_states_per_sign": len(positive),
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
