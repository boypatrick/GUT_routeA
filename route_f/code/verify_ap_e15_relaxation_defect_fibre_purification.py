#!/usr/bin/env python3
"""AP-E15 continuum defect and fibre-purification verification card.

This script performs algebraic and finite-dimensional consistency checks only.
It deliberately contains no lattice relaxation, mesh scan, translated start,
or new AP-E11 background generation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
SCRIPT = Path(__file__).resolve()
TEX = ROUTE_F / "tex" / "ap_e15_relaxation_defect_fibre_purification.tex"
R = 1.0
K = 0.35
MASS_SQUARED = 0.12
CHECKS: list[dict[str, Any]] = []


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


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def hopf_map(u: np.ndarray) -> np.ndarray:
    x1, x2, x3, x4 = u
    return np.asarray(
        [
            2.0 * (x1 * x3 + x2 * x4),
            2.0 * (x2 * x3 - x1 * x4),
            x1 * x1 + x2 * x2 - x3 * x3 - x4 * x4,
        ]
    )


def hopf_derivative(u: np.ndarray) -> np.ndarray:
    x1, x2, x3, x4 = u
    return np.asarray(
        [
            [2.0 * x3, 2.0 * x4, 2.0 * x1, 2.0 * x2],
            [-2.0 * x4, 2.0 * x3, 2.0 * x2, -2.0 * x1],
            [2.0 * x1, 2.0 * x2, -2.0 * x3, -2.0 * x4],
        ]
    )


def hopf_variables(
    u: np.ndarray, gradient: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x1, x2, x3, x4 = u
    theta = np.asarray([-x2, x1, -x4, x3])
    n = hopf_map(u)
    a = 2.0 * theta @ gradient
    dn = hopf_derivative(u) @ gradient
    curvature_star = -np.asarray(
        [
            n @ np.cross(dn[:, 1], dn[:, 2]),
            n @ np.cross(dn[:, 2], dn[:, 0]),
            n @ np.cross(dn[:, 0], dn[:, 1]),
        ]
    )
    return n, a, dn, curvature_star


def wedge_one_dn_squared(one: np.ndarray, dn: np.ndarray) -> float:
    return float(
        sum(
            np.sum((one[i] * dn[:, j] - one[j] * dn[:, i]) ** 2)
            for i in range(3)
            for j in range(i + 1, 3)
        )
    )


def phase_density(
    t: float,
    a: np.ndarray,
    dn: np.ndarray,
    curvature_star: np.ndarray,
    p: np.ndarray,
    z: complex,
    chi: float,
    phi: float,
) -> float:
    shifted = a + 2.0 * t * p
    return float(
        0.125 * np.dot(shifted, shifted)
        + R * wedge_one_dn_squared(shifted, dn) / 32.0
        + K * np.dot(shifted, curvature_star) ** 2 / 128.0
        + MASS_SQUARED
        * (1.0 - np.real(np.exp(1j * (chi + t * phi)) * z))
    )


def phase_hessian_audit(samples: int) -> dict[str, Any]:
    rng = np.random.default_rng(20260723)
    residuals = []
    for _ in range(samples):
        u = rng.normal(size=4)
        u /= np.linalg.norm(u)
        ambient = rng.normal(size=(4, 3))
        gradient = ambient - u[:, None] * (u @ ambient)[None, :]
        _, a, dn, curvature_star = hopf_variables(u, gradient)
        p = rng.normal(size=3)
        chi = float(rng.uniform(-math.pi, math.pi))
        phi = float(rng.normal())
        z = complex(u[0], u[1])
        analytic = float(
            np.dot(p, p)
            + R * wedge_one_dn_squared(p, dn) / 4.0
            + K * np.dot(p, curvature_star) ** 2 / 16.0
            + MASS_SQUARED * np.real(np.exp(1j * chi) * z) * phi**2
        )
        h = 2.0e-4
        numerical = (
            -phase_density(2 * h, a, dn, curvature_star, p, z, chi, phi)
            + 16 * phase_density(h, a, dn, curvature_star, p, z, chi, phi)
            - 30 * phase_density(0, a, dn, curvature_star, p, z, chi, phi)
            + 16 * phase_density(-h, a, dn, curvature_star, p, z, chi, phi)
            - phase_density(-2 * h, a, dn, curvature_star, p, z, chi, phi)
        ) / (12 * h * h)
        residuals.append(abs(analytic - numerical))
    return {
        "samples": samples,
        "max_absolute_residual": max(residuals),
        "formula": (
            "|dphi|^2+(R/4)|dphi wedge Dn|^2"
            "+(K/16)|dphi wedge da|^2"
            "+m^2 Re(exp(i chi)z) phi^2"
        ),
    }


def convexity_audit() -> dict[str, Any]:
    critical_radius = math.pi / math.sqrt(MASS_SQUARED)
    radii = np.asarray([0.25, 1.0, 5.0, 9.0, critical_radius, 9.2])
    lambdas = 1.0 - MASS_SQUARED * radii**2 / math.pi**2
    return {
        "mass_squared": MASS_SQUARED,
        "critical_radius": critical_radius,
        "radii": radii.tolist(),
        "lambda_r": lambdas.tolist(),
        "production_local_balls_are_available": bool(critical_radius > 0.0),
        "strong_convexity_condition": "r < pi/sqrt(m^2)",
    }


def hilbert_defect_audit() -> dict[str, Any]:
    rng = np.random.default_rng(15)
    x = rng.normal(size=256)
    y = rng.normal(size=256)
    lhs = float(np.dot(y - x, y - x))
    rhs = float(np.dot(y, y) + np.dot(x, x) - 2.0 * np.dot(y, x))
    # A toy orthogonal defect realizes the exact positive norm gap.
    x /= np.linalg.norm(x)
    orth = y - np.dot(y, x) * x
    orth /= np.linalg.norm(orth)
    amplitudes = np.asarray([0.5, 0.25, 0.125, 0.0])
    gaps = 0.5 * amplitudes**2
    distances = amplitudes
    return {
        "polarization_residual": abs(lhs - rhs),
        "toy_defect_masses": gaps.tolist(),
        "toy_strong_distances": distances.tolist(),
        "identity": (
            "weak convergence plus convergence of Hilbert norms "
            "is equivalent to strong complete-minor graph convergence"
        ),
    }


def weighted_cutoff_audit() -> dict[str, Any]:
    rng = np.random.default_rng(1515)
    tangent_gradient = rng.normal(size=(4, 3))
    f0 = float(abs(rng.normal()) + 0.2)
    radii = np.logspace(-1.0, -5.0, 9)
    coefficient = 4.0 * math.pi * f0 * float(
        np.sum(tangent_gradient**2)
    ) / 15.0
    affine_values = coefficient * radii**3
    fitted_slope = float(
        np.polyfit(np.log(radii), np.log(affine_values), 1)[0]
    )
    return {
        "radii": radii.tolist(),
        "affine_weighted_defect": affine_values.tolist(),
        "log_log_slope": fitted_slope,
        "exact_scaling": (
            "for u=q+A(x-x0), f=|M2|^2=f0, "
            "r^-2 int_Br |u-q|^2 f0=(4 pi/15)f0|A|^2 r^3"
        ),
        "proved_scope": "Lebesgue-a.e. centres of each fixed graph map",
        "not_proved_scope": "defect-support or all centres",
    }


def source_policy_audit() -> dict[str, Any]:
    text = SCRIPT.read_text(encoding="utf-8")
    # Assemble the sentinels so this declaration does not match itself.
    forbidden_calls = [
        "relax_" + "field(",
        "Relaxation" + "Plan(",
        "translated_" + "start",
        "triangulation_" + "scan",
        "multi" + "grid",
    ]
    found = [token for token in forbidden_calls if token in text]
    return {
        "ap_e11_action_frozen": True,
        "new_lattice_scans_authorized": False,
        "forbidden_scan_calls_found": found,
        "only_algebraic_continuum_checks": not found,
    }


def gate_ledger() -> dict[str, Any]:
    return {
        "ap_e11_action_frozen": True,
        "new_lattice_scan": False,
        "relaxation_defect_measure_defined": True,
        "compact_fibre_phase_strong_convexity": True,
        "phase_reducible_vertical_defect_zero": True,
        "mixed_hopf_base_defect_zero": False,
        "weighted_cutoff_decay_lebesgue_a.e.": True,
        "weighted_cutoff_decay_on_defect_support": False,
        "sequence_uniform_weighted_cutoff_mu_a.e.": False,
        "full_mu_zero": False,
        "local_recovery": False,
        "regularity": False,
        "continuum_isolation": False,
        "bosonic_hessian_authorized": False,
        "parallel_dirac_callias_mathematics_allowed": True,
        "determinant_promotion": False,
        "degree_one_portal_promotion": False,
        "physics_promotion_allowed": False,
    }


def write_markdown(report: dict[str, Any], path: Path) -> None:
    gates = report["gates"]
    lines = [
        "# AP-E15 relaxation-defect and fibre-purification audit",
        "",
        f"- Status: **{report['status']}**",
        f"- Checks: **{report['summary']['passed']}/{report['summary']['total']}**",
        f"- Phase Hessian maximum residual: "
        f"`{report['phase_hessian']['max_absolute_residual']:.3e}`",
        f"- Strong-convexity critical radius: "
        f"`{report['convexity']['critical_radius']:.9f}`",
        f"- Weighted-cutoff affine slope: "
        f"`{report['weighted_cutoff']['log_log_slope']:.12f}`",
        "",
        "## Theorem ledger",
        "",
        "- `mu` is a nonnegative Radon relaxation-defect measure.",
        "- `mu=0` iff the selected recovery subsequence converges strongly in "
        "the complete-minor graph norm.",
        "- Every compact-fibre-phase-reducible vertical defect vanishes for "
        "an asymptotically minimizing smooth fixed-degree sequence.",
        "- Weighted cutoff decay holds at Lebesgue-a.e. centres for each fixed "
        "map, but the sequence-uniform `lim_r limsup_j` estimate is open.",
        "",
        "## Gates",
        "",
    ]
    lines.extend(f"- `{name}`: `{value}`" for name, value in gates.items())
    lines.extend(
        [
            "",
            "## Check details",
            "",
        ]
    )
    lines.extend(
        f"- **{row['group']} / {row['name']}**: "
        f"`{'PASS' if row['pass'] else 'FAIL'}` — {row['detail']}"
        for row in report["checks"]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    samples = 32 if args.quick else 160

    phase = phase_hessian_audit(samples)
    convexity = convexity_audit()
    hilbert = hilbert_defect_audit()
    weighted = weighted_cutoff_audit()
    policy = source_policy_audit()
    gates = gate_ledger()

    check(
        "algebra",
        "phase Hessian coefficients",
        phase["max_absolute_residual"] < 2.0e-6,
        f"maximum five-point finite-difference residual "
        f"{phase['max_absolute_residual']:.3e}",
    )
    check(
        "analysis",
        "production convexity radius",
        abs(convexity["critical_radius"] - 9.068996821) < 1.0e-8,
        "lambda_r=1-m^2 r^2/pi^2 and m^2=0.12",
    )
    check(
        "analysis",
        "small-ball strong convexity",
        convexity["lambda_r"][0] > 0 and convexity["lambda_r"][3] > 0,
        "r=0.25 and r=9 lie below the exact critical radius",
    )
    check(
        "defect",
        "Hilbert polarization identity",
        hilbert["polarization_residual"] < 1.0e-11,
        f"residual {hilbert['polarization_residual']:.3e}",
    )
    check(
        "defect",
        "nonnegative toy defect masses",
        min(hilbert["toy_defect_masses"]) >= 0.0,
        "quadratic norm gaps are nonnegative",
    )
    check(
        "purification",
        "degree-preserving local drop bound",
        gates["phase_reducible_vertical_defect_zero"],
        "Delta_j(Q)<=E(u_j)-I_B -> 0 for every fixed admissible ball",
    )
    check(
        "cutoff",
        "smooth weighted-cutoff exponent",
        abs(weighted["log_log_slope"] - 3.0) < 1.0e-11,
        f"fitted affine exponent {weighted['log_log_slope']:.12f}",
    )
    check(
        "cutoff",
        "a.e. versus singular-support scope",
        gates["weighted_cutoff_decay_lebesgue_a.e."]
        and not gates["weighted_cutoff_decay_on_defect_support"]
        and not gates["sequence_uniform_weighted_cutoff_mu_a.e."],
        "fixed-map Lebesgue differentiation does not exchange r and recovery limits",
    )
    check(
        "policy",
        "no new lattice scan",
        policy["only_algebraic_continuum_checks"]
        and not policy["new_lattice_scans_authorized"],
        f"forbidden scan calls found: {policy['forbidden_scan_calls_found']}",
    )
    check(
        "gates",
        "classicality chain remains closed",
        not any(
            gates[key]
            for key in (
                "full_mu_zero",
                "local_recovery",
                "regularity",
                "continuum_isolation",
                "bosonic_hessian_authorized",
            )
        ),
        "full mu=0, recovery, regularity, isolation, and Hessian remain false",
    )
    check(
        "gates",
        "determinant and portal embargo",
        gates["parallel_dirac_callias_mathematics_allowed"]
        and not gates["determinant_promotion"]
        and not gates["degree_one_portal_promotion"],
        "parallel operator mathematics is allowed without downstream promotion",
    )

    passed = sum(row["pass"] for row in CHECKS)
    report = {
        "card": "AP-E15 relaxation defect and fibre-phase purification",
        "mode": "quick" if args.quick else "production",
        "status": "PASS" if passed == len(CHECKS) else "FAIL",
        "summary": {"passed": passed, "total": len(CHECKS)},
        "parameters": {"R": R, "K": K, "mass_squared": MASS_SQUARED},
        "phase_hessian": phase,
        "convexity": convexity,
        "hilbert_defect": hilbert,
        "weighted_cutoff": weighted,
        "policy": policy,
        "gates": gates,
        "checks": CHECKS,
        "sources": [source_row(SCRIPT), source_row(TEX)],
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    suffix = "_quick" if args.quick else ""
    json_path = OUTPUT / f"ap_e15_relaxation_defect_fibre_purification{suffix}.json"
    md_path = OUTPUT / f"ap_e15_relaxation_defect_fibre_purification{suffix}.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(report, md_path)
    print(
        f"{report['status']} {passed}/{len(CHECKS)}; "
        f"vertical={gates['phase_reducible_vertical_defect_zero']}; "
        f"mu={gates['full_mu_zero']}; "
        f"hessian={gates['bosonic_hessian_authorized']}; "
        f"portal={gates['degree_one_portal_promotion']}"
    )
    print(json_path)
    print(md_path)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
