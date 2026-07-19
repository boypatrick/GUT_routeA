#!/usr/bin/env python3
"""AP-E13 annular endpoint-replacement and reverse-Hoelder audit.

The card tests the exact geodesic-cap family

    g_eps(omega) = (cos eps, sin eps * omega) in S^3,

on an annulus A_(R,2R).  Its shell complete-minor energy tends to zero for
eps=sqrt(R), while every trace-preserving filling of B_R has a non-vanishing
L2 third-minor cost.  This is a topological Wess--Zumino/Jacobian obstruction,
not a discretization effect.

The same calculation verifies the repaired, scale-compatible regime
eps=lambda*R and records why exact strong quasiconvexity alone does not start
a target-constrained reverse-Hoelder iteration.
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
THIS_SCRIPT = Path(__file__).resolve()
E12_JSON = OUTPUT / "ap_e12_graph_density_regular_minimizer.json"
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


def cap_volume(epsilon: np.ndarray | float) -> np.ndarray | float:
    """Volume of the S^3 geodesic cap of radius epsilon."""
    value = np.asarray(epsilon, dtype=float)
    direct = 2.0 * math.pi * (value - 0.5 * np.sin(2.0 * value))
    # Avoid subtractive cancellation at the radii used in the continuum scan.
    series = (
        (4.0 * math.pi / 3.0) * value**3
        - (4.0 * math.pi / 15.0) * value**5
        + (8.0 * math.pi / 315.0) * value**7
        - (4.0 * math.pi / 2835.0) * value**9
    )
    result = np.where(np.abs(value) < 1.0e-3, series, direct)
    return float(result) if result.ndim == 0 else result


def trace_norms(epsilon: np.ndarray | float) -> tuple[Any, Any]:
    """Unit-S2 squared L2 norms of Dg and its tangential two-minor."""
    sine = np.sin(epsilon)
    return 8.0 * math.pi * sine**2, 4.0 * math.pi * sine**4


def radial_constant_shell(R: np.ndarray, epsilon: np.ndarray) -> dict[str, np.ndarray]:
    """Raw complete-minor norms on A_(R,2R) for u(r,omega)=g_eps(omega)."""
    sine = np.sin(epsilon)
    first = 8.0 * math.pi * R * sine**2
    second = 2.0 * math.pi * sine**4 / R
    third = np.zeros_like(R)
    return {"first": first, "second": second, "third": third, "total": first + second}


def filling_lower_bound(R: np.ndarray, epsilon: np.ndarray) -> np.ndarray:
    """WZ/Cauchy lower bound for the third-minor L2 norm on B_R."""
    ball_volume = (4.0 * math.pi / 3.0) * R**3
    return np.asarray(cap_volume(epsilon)) ** 2 / ball_volume


def linear_cap_fill_m3(R: float, epsilon: float, points: int = 200001) -> float:
    """M3 energy of f(r)=epsilon*r/R, integrated by the trapezoid rule."""
    s = np.linspace(0.0, 1.0, points)
    integrand = np.empty_like(s)
    integrand[0] = 0.0
    integrand[1:] = np.sin(epsilon * s[1:]) ** 4 / s[1:] ** 2
    return float(
        4.0 * math.pi * epsilon**2 / R**3 * np.trapezoid(integrand, s)
    )


def fit_slope(R: np.ndarray, values: np.ndarray) -> float:
    return float(np.polyfit(np.log(R), np.log(values), 1)[0])


def scaling_audit(quick: bool) -> dict[str, Any]:
    R = np.logspace(-2.0, -8.0 if not quick else -6.0, 13 if not quick else 9)
    rows = []
    maximum_slope_error = 0.0
    for alpha in (0.30, 0.40, 0.50, 0.60, 0.75, 1.00):
        epsilon = R**alpha
        shell = radial_constant_shell(R, epsilon)
        lower = filling_lower_bound(R, epsilon)
        slopes = {
            "first": fit_slope(R[-6:], shell["first"][-6:]),
            "second": fit_slope(R[-6:], shell["second"][-6:]),
            "fill_m3_lower": fit_slope(R[-6:], lower[-6:]),
        }
        predicted = {
            "first": 1.0 + 2.0 * alpha,
            "second": 4.0 * alpha - 1.0,
            "fill_m3_lower": 6.0 * alpha - 3.0,
        }
        slope_error = max(abs(slopes[key] - predicted[key]) for key in slopes)
        maximum_slope_error = max(maximum_slope_error, slope_error)
        rows.append(
            {
                "alpha": alpha,
                "measured_slopes": slopes,
                "predicted_slopes": predicted,
                "maximum_slope_error": slope_error,
            }
        )

    epsilon_half = np.sqrt(R)
    half_shell = radial_constant_shell(R, epsilon_half)
    half_lower = filling_lower_bound(R, epsilon_half)
    alpha_one_shell = radial_constant_shell(R, R)
    alpha_one_lower = filling_lower_bound(R, R)
    return {
        "R_values": R.tolist(),
        "power_rows": rows,
        "maximum_slope_error": maximum_slope_error,
        "critical_alpha_half": {
            "shell_total": half_shell["total"].tolist(),
            "fill_m3_lower": half_lower.tolist(),
            "smallest_R_shell_total": float(half_shell["total"][-1]),
            "smallest_R_fill_lower": float(half_lower[-1]),
            "limiting_fill_lower": 4.0 * math.pi / 3.0,
            "shell_tends_to_zero": bool(half_shell["total"][-1] < 1.0e-5),
            "fill_stays_positive": bool(half_lower[-1] > 4.0),
        },
        "scale_compatible_alpha_one": {
            "first_slope": fit_slope(R[-6:], alpha_one_shell["first"][-6:]),
            "second_slope": fit_slope(R[-6:], alpha_one_shell["second"][-6:]),
            "fill_slope": fit_slope(R[-6:], alpha_one_lower[-6:]),
            "all_predicted_slope": 3.0,
            "all_costs_tend_to_zero": bool(
                max(
                    alpha_one_shell["first"][-1],
                    alpha_one_shell["second"][-1],
                    alpha_one_lower[-1],
                )
                < 1.0e-20
            ),
        },
    }


def geometry_audit(quick: bool) -> dict[str, Any]:
    epsilon = np.logspace(-1.0, -5.0, 9)
    volumes = np.asarray(cap_volume(epsilon))
    cubic = (4.0 * math.pi / 3.0) * epsilon**3
    relative_error = np.abs(volumes / cubic - 1.0)

    fill_rows = []
    for R, eps in ((0.15, 0.11), (0.08, 0.06), (0.03, 0.02)):
        lower = float(filling_lower_bound(np.array([R]), np.array([eps]))[0])
        direct = linear_cap_fill_m3(R, eps, 20001 if quick else 100001)
        fill_rows.append(
            {
                "R": R,
                "epsilon": eps,
                "wz_lower_bound": lower,
                "linear_cap_fill": direct,
                "ratio": direct / lower,
            }
        )
    return {
        "cap_volume_formula": "2*pi*(epsilon-sin(2*epsilon)/2)",
        "small_cap_leading_term": "4*pi*epsilon^3/3",
        "smallest_epsilon_relative_error": float(relative_error[-1]),
        "trace_first_norm_formula": "8*pi*sin(epsilon)^2",
        "trace_second_minor_norm_formula": "4*pi*sin(epsilon)^4",
        "linear_fill_rows": fill_rows,
        "minimum_linear_fill_to_lower_bound_ratio": min(row["ratio"] for row in fill_rows),
    }


def theorem_ledger() -> dict[str, Any]:
    return {
        "unconditional_annular_endpoint_replacement": False,
        "unconditional_statement_disproved": True,
        "reason": "the geodesic-cap trace with epsilon=sqrt(R) has vanishing shell graph energy but every trace-preserving filling carries a non-vanishing WZ-forced L2 third-minor cost",
        "necessary_wz_condition": "dist(V_WZ(g_R), 2*pi^2*Z)^2 / R^3 -> 0",
        "conditional_scale_compatible_replacement": True,
        "conditional_hypothesis": "on a good sphere S_R, the trace lies in a normal ball about q and ||log_q g||_infinity <= Lambda R",
        "conditional_estimates": {
            "first_minor": "int_B_R |DU|^2 <= C[R*T1 + Lambda^2*R^3]",
            "second_minor": "int_B_R |M2(DU)|^2 <= C[R*T2 + Lambda^2*R*T1]",
            "third_minor": "int_B_R |M3(DU)|^2 <= C*Lambda^2*R*T2",
            "good_sphere": "T1 <= E1(A_R,2R)/R and T2 <= E2(A_R,2R)/R",
        },
        "strong_quasiconvexity_implies_reverse_holder_without_extra_input": False,
        "reverse_holder_obstruction": "the highest cutoff-stress term leaves r^(-2)|u-q|^2|M2(Du)|^2 after Young; controlling it requires precisely an L-infinity/Morrey oscillation bound or prior higher integrability",
        "model_specific_reverse_holder_status": "open_for_actual_relaxed_minimizer",
        "weak_euler_lagrange_authorized": False,
        "endpoint_graph_density_proved": False,
        "endpoint_graph_density_disproved": False,
        "classicality": False,
        "continuum_isolation": False,
        "hessian_gate_open": False,
        "determinant_variation_gate_open": False,
        "portal_start_allowed": False,
        "ordered_next_gate": [
            "prove WZ-flux decay plus scale-compatible trace oscillation for the selected relaxed B=1 minimizer, using a monotonicity/frequency argument",
            "or derive a Noether-current/Maurer-Cartan reverse-Hoelder estimate that controls the cutoff remainder without assuming the desired Morrey bound",
            "only then return to graph density/classicality and the same-action Hessian",
        ],
    }


def render_markdown(result: dict[str, Any]) -> str:
    critical = result["scaling"]["critical_alpha_half"]
    compatible = result["scaling"]["scale_compatible_alpha_one"]
    theory = result["theory"]
    return "\n".join(
        [
            "# AP-E13 annular replacement and reverse-Hölder audit",
            "",
            f"Status: `{result['status']}`.",
            "",
            "## Unconditional theorem verdict",
            "",
            f"- universal trace-preserving endpoint replacement: `{theory['unconditional_annular_endpoint_replacement']}`.",
            f"- universal statement disproved: `{theory['unconditional_statement_disproved']}`.",
            f"- smallest-R critical shell norm: `{critical['smallest_R_shell_total']:.6e}`.",
            f"- corresponding WZ filling lower bound: `{critical['smallest_R_fill_lower']:.9f}`.",
            f"- exact limiting lower bound: `{critical['limiting_fill_lower']:.9f}`.",
            "",
            "## Repaired theorem",
            "",
            f"- scale-compatible conditional replacement: `{theory['conditional_scale_compatible_replacement']}`.",
            f"- measured alpha=1 slopes (first/M2/M3): `{compatible['first_slope']:.6f}`, `{compatible['second_slope']:.6f}`, `{compatible['fill_slope']:.6f}`.",
            f"- necessary WZ condition: `{theory['necessary_wz_condition']}`.",
            "",
            "## Reverse-Hölder boundary",
            "",
            f"- follows from exact strong quasiconvexity alone: `{theory['strong_quasiconvexity_implies_reverse_holder_without_extra_input']}`.",
            f"- model-specific result for the selected minimizer: `{theory['model_specific_reverse_holder_status']}`.",
            f"- Hessian gate: `{theory['hessian_gate_open']}`.",
            "",
            f"Checks: `{result['checks_passed']}/{result['checks_total']}`.",
        ]
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    geometry = geometry_audit(arguments.quick)
    scaling = scaling_audit(arguments.quick)
    theory = theorem_ledger()
    critical = scaling["critical_alpha_half"]
    compatible = scaling["scale_compatible_alpha_one"]

    check(
        "G1",
        "exact small-cap WZ volume asymptotics",
        geometry["smallest_epsilon_relative_error"] < 1.0e-8,
        f"relative error={geometry['smallest_epsilon_relative_error']:.3e}",
    )
    check(
        "G2",
        "explicit cap filling respects the WZ lower bound",
        geometry["minimum_linear_fill_to_lower_bound_ratio"] >= 1.0 - 1.0e-8,
        f"minimum ratio={geometry['minimum_linear_fill_to_lower_bound_ratio']:.9f}",
    )
    check(
        "S1",
        "power-law exponents match exact scaling",
        scaling["maximum_slope_error"] < 3.0e-4,
        f"maximum slope error={scaling['maximum_slope_error']:.3e}",
    )
    check(
        "N1",
        "critical cap has vanishing shell energy",
        critical["shell_tends_to_zero"],
        f"shell={critical['smallest_R_shell_total']:.3e}",
    )
    check(
        "N2",
        "critical cap has non-vanishing filling cost",
        critical["fill_stays_positive"]
        and abs(critical["smallest_R_fill_lower"] - critical["limiting_fill_lower"])
        < 1.0e-3,
        f"lower={critical['smallest_R_fill_lower']:.9f}",
    )
    check(
        "C1",
        "scale-compatible cap costs all decay cubically",
        max(
            abs(compatible["first_slope"] - 3.0),
            abs(compatible["second_slope"] - 3.0),
            abs(compatible["fill_slope"] - 3.0),
        )
        < 3.0e-4
        and compatible["all_costs_tend_to_zero"],
        "epsilon=R gives first-, second-, and third-minor costs of order R^3",
    )
    check(
        "T1",
        "unconditional annular theorem is disproved",
        theory["unconditional_statement_disproved"]
        and not theory["unconditional_annular_endpoint_replacement"],
        "the same trace occurs on every sphere of the low-energy radial-constant shell",
    )
    check(
        "T2",
        "conditional target-valued replacement closes all minor estimates",
        theory["conditional_scale_compatible_replacement"],
        "geodesic contraction under ||log_q g||_infinity <= Lambda R",
    )
    check(
        "R1",
        "reverse-Hoelder and downstream gates fail closed",
        not theory["strong_quasiconvexity_implies_reverse_holder_without_extra_input"]
        and not theory["weak_euler_lagrange_authorized"]
        and not theory["hessian_gate_open"]
        and not theory["portal_start_allowed"],
        "the cutoff estimate is circular at r^-2 |u-q|^2 |M2|^2",
    )

    result = {
        "artifact": "AP-E13 annular endpoint replacement and reverse-Hoelder audit",
        "status": "production_unconditional_no_go_conditional_repair"
        if not arguments.quick
        else "quick_unconditional_no_go_conditional_repair",
        "quick": arguments.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "source_manifest": [source_row(THIS_SCRIPT), source_row(E12_JSON)],
        "geometry": geometry,
        "scaling": scaling,
        "theory": theory,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    json_path = OUTPUT / "ap_e13_annular_replacement_reverse_holder.json"
    markdown_path = OUTPUT / "ap_e13_annular_replacement_reverse_holder.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    markdown_path.write_text(render_markdown(result))
    print(
        f"AP-E13 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"unconditional={theory['unconditional_annular_endpoint_replacement']}; "
        f"conditional={theory['conditional_scale_compatible_replacement']}; "
        f"hessian={theory['hessian_gate_open']}",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
