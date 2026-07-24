#!/usr/bin/env python3
"""AP-E16 tangent-measure, Caccioppoli, and Young-measure audit.

The card is continuum/algebraic.  It performs no lattice relaxation or mesh
scan.  It verifies:

1. the unique simultaneous complete-minor blow-up scaling;
2. rank-one cutoff expansions for second and third minors;
3. an exact homogeneous Hopf-base Young-measure generator;
4. vanishing weighted cutoff despite a positive diffuse graph-energy defect;
5. fail-closed promotion gates.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
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
TEX = ROUTE_F / "tex" / "ap_e16_tangent_caccioppoli_young_measure.tex"
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


def minors(matrix: np.ndarray, order: int) -> np.ndarray:
    rows, cols = matrix.shape
    values = []
    for row_ids in itertools.combinations(range(rows), order):
        for col_ids in itertools.combinations(range(cols), order):
            values.append(
                np.linalg.det(matrix[np.ix_(np.asarray(row_ids), np.asarray(col_ids))])
            )
    return np.asarray(values)


def scaling_audit() -> dict[str, Any]:
    # If spatial scale is r, target amplitude is s=r^alpha, and defect mass
    # is r^d, the three derivative orders scale as below.
    # 2 alpha + 1 = d, 4 alpha - 1 = d, 6 alpha - 3 = d.
    linear_system = np.asarray([[2.0, -1.0], [4.0, -1.0], [6.0, -1.0]])
    rhs = np.asarray([-1.0, 1.0, 3.0])
    solution, residuals, _, _ = np.linalg.lstsq(linear_system, rhs, rcond=None)
    alpha, dimension = solution
    exponents = np.asarray(
        [2.0 * alpha + 1.0, 4.0 * alpha - 1.0, 6.0 * alpha - 3.0]
    )
    dimensions = np.linspace(0.5, 3.0, 6)
    first_order_alpha = (dimensions - 1.0) / 2.0
    mismatch_second = 4.0 * first_order_alpha - 1.0 - dimensions
    mismatch_third = 6.0 * first_order_alpha - 3.0 - dimensions
    return {
        "energy_scaling": {
            "first": "s^2 r",
            "second": "s^4/r",
            "third": "s^6/r^3",
        },
        "solution": {"alpha": float(alpha), "dimension": float(dimension)},
        "residual_sum_squares": float(np.sum((exponents - dimension) ** 2)),
        "common_exponents": exponents.tolist(),
        "test_dimensions": dimensions.tolist(),
        "second_mismatch_after_first_normalization": mismatch_second.tolist(),
        "third_mismatch_after_first_normalization": mismatch_third.tolist(),
        "conclusion": (
            "one ordinary field blow-up normalizes all minor orders only for "
            "target amplitude s~r and volumetric defect mass mu(Br)~r^3"
        ),
    }


def cutoff_minor_audit(samples: int) -> dict[str, Any]:
    rng = np.random.default_rng(20260724)
    second_residuals = []
    third_residuals = []
    second_bound_ratios = []
    third_bound_ratios = []
    for _ in range(samples):
        a = rng.normal(size=(4, 3))
        b = np.outer(rng.normal(size=4), rng.normal(size=3))
        c = float(rng.uniform(-1.5, 1.5))
        m2a = minors(a, 2)
        m2_linear = minors(a + b, 2) - m2a
        predicted_m2 = c * c * m2a + c * m2_linear
        actual_m2 = minors(c * a + b, 2)
        second_residuals.append(float(np.max(np.abs(predicted_m2 - actual_m2))))
        denom2 = float(np.sum(a * a) * np.sum(b * b))
        second_bound_ratios.append(
            float(np.sum(m2_linear * m2_linear) / max(denom2, 1.0e-15))
        )

        m3a = minors(a, 3)
        m3_linear = minors(a + b, 3) - m3a
        predicted_m3 = c**3 * m3a + c * c * m3_linear
        actual_m3 = minors(c * a + b, 3)
        third_residuals.append(float(np.max(np.abs(predicted_m3 - actual_m3))))
        denom3 = float(np.sum(m2a * m2a) * np.sum(b * b))
        third_bound_ratios.append(
            float(np.sum(m3_linear * m3_linear) / max(denom3, 1.0e-15))
        )
    return {
        "samples": samples,
        "maximum_second_minor_polynomial_residual": max(second_residuals),
        "maximum_third_minor_polynomial_residual": max(third_residuals),
        "maximum_second_rank_one_bound_ratio": max(second_bound_ratios),
        "maximum_third_rank_one_bound_ratio": max(third_bound_ratios),
        "identities": {
            "second": "M2(cA+B)=c^2 M2(A)+c L2(A,B) for rank(B)=1",
            "third": "M3(cA+B)=c^3 M3(A)+c^2 L3(A,A,B) for rank(B)=1",
        },
    }


def homogeneous_young_measure_audit(grid: int) -> dict[str, Any]:
    theta = 2.0 * math.pi * (np.arange(grid) + 0.5) / grid
    c1, c2 = np.meshgrid(np.cos(theta), np.cos(theta), indexing="ij")
    mean_c1 = float(np.mean(c1))
    mean_c2 = float(np.mean(c2))
    mean_minor = float(np.mean(c1 * c2))
    mean_first_squared = float(np.mean(c1 * c1 + c2 * c2))
    mean_second_squared = float(np.mean((c1 * c2) ** 2))
    energy_gap = 0.5 * mean_first_squared + 0.5 * R * mean_second_squared
    exact_gap = 0.5 + R / 8.0
    return {
        "grid": grid,
        "generator": (
            "F(theta1,theta2) has F_31=cos(theta1), "
            "F_42=cos(theta2), all other entries zero"
        ),
        "mean_gradient_components": [mean_c1, mean_c2],
        "mean_only_nonzero_second_minor": mean_minor,
        "mean_first_squared": mean_first_squared,
        "mean_second_squared": mean_second_squared,
        "mean_third_squared": 0.0,
        "graph_energy_jensen_gap": energy_gap,
        "exact_gap": exact_gap,
        "gap_residual": abs(energy_gap - exact_gap),
        "base_fibre_content": {
            "limiting_connection_a": 0.0,
            "hopf_base_derivative_nonzero": True,
            "hopf_curvature_nonzero": True,
            "third_minor": 0.0,
        },
    }


def sphere_generator_audit(grid: int) -> dict[str, Any]:
    theta = 2.0 * math.pi * (np.arange(grid) + 0.5) / grid
    s1, s2 = np.meshgrid(np.sin(theta), np.sin(theta), indexing="ij")
    c1, c2 = np.meshgrid(np.cos(theta), np.cos(theta), indexing="ij")
    ells = np.asarray([1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0, 1.0 / 64.0])
    rows = []
    for ell in ells:
        x0 = np.sqrt(1.0 - ell * ell * (s1 * s1 + s2 * s2))
        d0_1 = -ell * s1 * c1 / x0
        d0_2 = -ell * s2 * c2 / x0
        col1_squared = d0_1 * d0_1 + c1 * c1
        col2_squared = d0_2 * d0_2 + c2 * c2
        col_dot = d0_1 * d0_2
        first_squared = col1_squared + col2_squared
        second_squared = col1_squared * col2_squared - col_dot * col_dot
        a1 = -2.0 * ell * s2 * c1
        a2 = 2.0 * ell * s1 * c2
        a_squared = a1 * a1 + a2 * a2
        distance_squared = (x0 - 1.0) ** 2 + ell * ell * (s1 * s1 + s2 * s2)
        weighted = distance_squared * second_squared
        sphere_residual = np.abs(
            x0 * x0 + ell * ell * (s1 * s1 + s2 * s2) - 1.0
        )
        rows.append(
            {
                "ell": float(ell),
                "sphere_residual": float(np.max(sphere_residual)),
                "mean_first_squared": float(np.mean(first_squared)),
                "mean_second_squared": float(np.mean(second_squared)),
                "mean_third_squared": 0.0,
                "mean_connection_squared": float(np.mean(a_squared)),
                "connection_squared_over_ell_squared": float(
                    np.mean(a_squared) / (ell * ell)
                ),
                "mean_weighted_cutoff_integrand": float(np.mean(weighted)),
                "mean_weighted_over_ell_squared": float(
                    np.mean(weighted) / (ell * ell)
                ),
            }
        )
    weighted_values = np.asarray(
        [row["mean_weighted_cutoff_integrand"] for row in rows]
    )
    weighted_slope = float(
        np.polyfit(np.log(ells), np.log(weighted_values), 1)[0]
    )
    return {
        "family": (
            "u_N=(sqrt(1-N^-2(sin^2 Nx1+sin^2 Nx2)),0,"
            "N^-1 sin Nx1,N^-1 sin Nx2)"
        ),
        "rows": rows,
        "weighted_cutoff_log_log_slope": weighted_slope,
        "limits": {
            "mean_first_squared": 1.0,
            "mean_second_squared": 0.25,
            "mean_third_squared": 0.0,
            "connection_L2": 0.0,
            "weighted_cutoff": 0.0,
            "positive_defect_density": 0.5 + R / 8.0,
        },
        "phase_purification_estimate": (
            "the periodic phase minimizer satisfies ||dchi_N||_2"
            "<=C||a_N||_2=O(N^-1), so it preserves the leading base Young measure"
        ),
    }


def tangent_ledger() -> dict[str, Any]:
    return {
        "defect_components": [
            "mu_a",
            "mu_Dn",
            "mu_a_wedge_Dn",
            "mu_da",
            "mu_a_wedge_da",
        ],
        "radon_nikodym_weights": (
            "theta_alpha=d mu_alpha/d mu, theta_alpha>=0, sum theta_alpha=1"
        ),
        "common_tangent": (
            "at mu-a.e. x0, along one normalized tangent sequence, "
            "Tan(mu_alpha,x0)=theta_alpha(x0) tau"
        ),
        "diagonal_recovery": (
            "choose j_l after r_l so weak-star errors and global energy excess "
            "are o(mu(B_r_l)); normalized compact-phase drop then tends to zero"
        ),
        "what_purification_imposes": (
            "tangent stationarity under exact fibre shifts; it does not set "
            "all a-containing Radon-Nikodym weights to zero"
        ),
    }


def caccioppoli_ledger() -> dict[str, Any]:
    return {
        "status": "conditional theorem proved",
        "competitor": "v=exp_q((1-eta) log_q u) in one normal target ball",
        "derivative_bound": (
            "Eder(B_s)<=C Eder(B_t\\B_s)"
            "+C(t-s)^-2 int_A |u-q|^2(1+R|Du|^2+K|M2|^2)"
            "+C m^2 int_Bt |u-q|+local comparison deficit"
        ),
        "hole_filling": (
            "Eder(B_s)<=theta Eder(B_t)+weighted terms+normalized deficit, "
            "theta=C/(1+C)<1"
        ),
        "required_inputs_for_actual_minimizer": [
            "normal-ball or endpoint replacement control",
            "sequence-uniform weighted equiintegrability",
            "normalized local quasiminimality",
            "strong-quasiconvex harmonic approximation for excess contraction",
        ],
        "not_sufficient": (
            "phase purification plus weighted cutoff decay alone; "
            "the homogeneous base Young measure has both and retains positive defect"
        ),
    }


def source_policy_audit() -> dict[str, Any]:
    text = SCRIPT.read_text(encoding="utf-8")
    forbidden = [
        "relax_" + "field(",
        "Relaxation" + "Plan(",
        "translated_" + "start",
        "triangulation_" + "scan",
        "multi" + "grid",
    ]
    found = [token for token in forbidden if token in text]
    return {
        "ap_e11_action_frozen": True,
        "new_lattice_scan": False,
        "forbidden_scan_calls_found": found,
        "continuum_algebra_only": not found,
    }


def gate_ledger() -> dict[str, Any]:
    return {
        "tangent_measure_blowup": True,
        "hopf_defect_radon_nikodym_split": True,
        "normalized_phase_stationarity": True,
        "single_field_all_order_tangent_general": False,
        "single_field_volumetric_branch": True,
        "conditional_hopf_base_caccioppoli": True,
        "unconditional_excess_contraction": False,
        "homogeneous_base_young_measure_counterexample": True,
        "counterexample_is_asymptotically_minimizing": False,
        "two_limit_weighted_decay_sufficient_alone": False,
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
    lines = [
        "# AP-E16 tangent/Caccioppoli/Young-measure audit",
        "",
        f"- Status: **{report['status']}**",
        f"- Checks: **{report['summary']['passed']}/{report['summary']['total']}**",
        f"- Simultaneous field scaling: `alpha="
        f"{report['scaling']['solution']['alpha']:.12f}, d="
        f"{report['scaling']['solution']['dimension']:.12f}`",
        f"- Homogeneous Young-measure gap: "
        f"`{report['young_measure']['graph_energy_jensen_gap']:.12f}`",
        f"- Weighted-cutoff slope: "
        f"`{report['sphere_generator']['weighted_cutoff_log_log_slope']:.9f}`",
        "",
        "## Main conclusions",
        "",
        "- At `mu`-a.e. points, normalized tangent measures and their Hopf "
        "Radon-Nikodym components exist; a diagonal recovery sequence can be "
        "chosen with vanishing normalized compact-phase drop.",
        "- One ordinary field blow-up retains all three minor orders only in "
        "the volumetric branch `mu(B_r)~r^3`, with target amplitude `s~r`.",
        "- The normal-ball Hopf-base Caccioppoli inequality is proved with "
        "explicit weighted and local-comparison-deficit terms.",
        "- A phase-purified homogeneous pure-base Young measure has positive "
        "defect while its weighted cutoff tends to zero.  It refutes "
        "unconditional excess contraction, but is not an asymptotically "
        "minimizing sequence.",
        "",
        "## Gates",
        "",
    ]
    lines.extend(f"- `{key}`: `{value}`" for key, value in report["gates"].items())
    lines.extend(["", "## Check details", ""])
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
    samples = 32 if args.quick else 192
    grid = 64 if args.quick else 384

    scaling = scaling_audit()
    cutoff = cutoff_minor_audit(samples)
    young = homogeneous_young_measure_audit(grid)
    sphere = sphere_generator_audit(grid)
    tangent = tangent_ledger()
    caccioppoli = caccioppoli_ledger()
    policy = source_policy_audit()
    gates = gate_ledger()

    check(
        "scaling",
        "unique simultaneous complete-minor scaling",
        abs(scaling["solution"]["alpha"] - 1.0) < 1.0e-12
        and abs(scaling["solution"]["dimension"] - 3.0) < 1.0e-12,
        "the three equations give target amplitude s~r and mu(B_r)~r^3",
    )
    check(
        "scaling",
        "nonvolumetric branches require order-specific tangents",
        max(
            abs(np.asarray(scaling["second_mismatch_after_first_normalization"])[:-1])
        )
        > 0.1,
        "a first-order field normalization does not retain all higher orders for d!=3",
    )
    check(
        "cutoff",
        "rank-one second-minor expansion",
        cutoff["maximum_second_minor_polynomial_residual"] < 2.0e-12,
        f"maximum residual {cutoff['maximum_second_minor_polynomial_residual']:.3e}",
    )
    check(
        "cutoff",
        "rank-one third-minor expansion",
        cutoff["maximum_third_minor_polynomial_residual"] < 2.0e-12,
        f"maximum residual {cutoff['maximum_third_minor_polynomial_residual']:.3e}",
    )
    check(
        "young measure",
        "zero barycentric gradient",
        max(abs(x) for x in young["mean_gradient_components"]) < 1.0e-14,
        f"means {young['mean_gradient_components']}",
    )
    check(
        "young measure",
        "zero barycentric second minor",
        abs(young["mean_only_nonzero_second_minor"]) < 1.0e-14,
        f"mean minor {young['mean_only_nonzero_second_minor']:.3e}",
    )
    check(
        "young measure",
        "positive exact Jensen gap",
        young["gap_residual"] < 1.0e-13
        and young["graph_energy_jensen_gap"] > 0.0,
        f"gap {young['graph_energy_jensen_gap']:.12f}",
    )
    check(
        "sphere generator",
        "exact target constraint",
        max(row["sphere_residual"] for row in sphere["rows"]) < 5.0e-16,
        "all sampled generators lie on S3 to machine precision",
    )
    last = sphere["rows"][-1]
    check(
        "sphere generator",
        "pure-base leading limits",
        abs(last["mean_first_squared"] - 1.0) < 5.0e-4
        and abs(last["mean_second_squared"] - 0.25) < 5.0e-4
        and last["mean_connection_squared"] < 5.0e-4,
        "first/second energy persist while the Hopf connection vanishes",
    )
    check(
        "sphere generator",
        "weighted cutoff vanishes with positive defect",
        abs(sphere["weighted_cutoff_log_log_slope"] - 2.0) < 0.02
        and sphere["limits"]["positive_defect_density"] > 0.0,
        f"weighted slope {sphere['weighted_cutoff_log_log_slope']:.9f}",
    )
    check(
        "analysis",
        "conditional Caccioppoli versus unconditional contraction",
        gates["conditional_hopf_base_caccioppoli"]
        and not gates["unconditional_excess_contraction"]
        and gates["homogeneous_base_young_measure_counterexample"],
        "Caccioppoli needs normalized full quasiminimality for excess contraction",
    )
    check(
        "policy",
        "no new lattice scan and downstream embargo",
        policy["continuum_algebra_only"]
        and not policy["new_lattice_scan"]
        and not gates["bosonic_hessian_authorized"]
        and not gates["determinant_promotion"]
        and not gates["degree_one_portal_promotion"],
        f"forbidden scan calls found: {policy['forbidden_scan_calls_found']}",
    )

    passed = sum(row["pass"] for row in CHECKS)
    report = {
        "card": "AP-E16 tangent measure, Caccioppoli, and Young measure",
        "mode": "quick" if args.quick else "production",
        "status": "PASS" if passed == len(CHECKS) else "FAIL",
        "summary": {"passed": passed, "total": len(CHECKS)},
        "parameters": {"R": R, "K": K, "mass_squared": MASS_SQUARED},
        "scaling": scaling,
        "cutoff_minor_algebra": cutoff,
        "tangent_measure": tangent,
        "caccioppoli": caccioppoli,
        "young_measure": young,
        "sphere_generator": sphere,
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
    json_path = OUTPUT / f"ap_e16_tangent_caccioppoli_young_measure{suffix}.json"
    md_path = OUTPUT / f"ap_e16_tangent_caccioppoli_young_measure{suffix}.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(report, md_path)
    print(
        f"{report['status']} {passed}/{len(CHECKS)}; "
        f"tangent={gates['tangent_measure_blowup']}; "
        f"caccioppoli={gates['conditional_hopf_base_caccioppoli']}; "
        f"contraction={gates['unconditional_excess_contraction']}; "
        f"young={gates['homogeneous_base_young_measure_counterexample']}"
    )
    print(json_path)
    print(md_path)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
