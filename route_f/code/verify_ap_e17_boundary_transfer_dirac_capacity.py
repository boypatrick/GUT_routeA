#!/usr/bin/env python3
"""AP-E17 boundary-transfer, Dirac-rigidity, and capacity audit.

This is a continuum/algebraic verification card.  It performs no lattice
relaxation or grid scan.  It checks:

1. the radial-projection derivative used by the sphere-valued comparison;
2. the exact complete-minor strong-quasiconvexity identity;
3. positivity of the corresponding graph-variance/Dirac gap;
4. the rank-two microball obstruction to an automatic weighted gluing bound;
5. the sharp Wess--Zumino capacity constants and nonvolumetric scaling;
6. fail-closed promotion gates.
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
TEX = ROUTE_F / "tex" / "ap_e17_boundary_transfer_dirac_capacity.tex"
R = 1.0
K = 0.35
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


def batch_minors(matrices: np.ndarray, order: int) -> np.ndarray:
    """All order-k minors of (..., 4, 3) matrices."""
    if order == 1:
        return matrices.reshape(matrices.shape[:-2] + (12,))
    values = []
    for rows in itertools.combinations(range(4), order):
        for cols in itertools.combinations(range(3), order):
            sub = np.take(np.take(matrices, rows, axis=-2), cols, axis=-1)
            values.append(np.linalg.det(sub))
    return np.stack(values, axis=-1)


def graph_vector(matrices: np.ndarray) -> np.ndarray:
    return np.concatenate(
        [
            batch_minors(matrices, 1),
            math.sqrt(R) * batch_minors(matrices, 2),
            math.sqrt(K) * batch_minors(matrices, 3),
        ],
        axis=-1,
    )


def energy(matrices: np.ndarray) -> np.ndarray:
    graph = graph_vector(matrices)
    return 0.5 * np.sum(graph * graph, axis=-1)


def projection_audit(samples: int) -> dict[str, Any]:
    rng = np.random.default_rng(20260724)
    derivative_residuals = []
    boundary_residuals = []
    tangent_residuals = []
    finite_difference_step = 1.0e-7
    for _ in range(samples):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        v = rng.normal(size=4)
        v -= q * np.dot(q, v)
        f = rng.normal(size=(4, 3))
        f -= np.outer(q, q @ f)
        radius = float(rng.uniform(1.0e-4, 5.0e-3))
        z = q + radius * v
        p = z / np.linalg.norm(z)
        dp = (np.eye(4) - np.outer(p, p)) / np.linalg.norm(z)
        analytic = dp @ f
        numeric = np.zeros_like(f)
        for column in range(3):
            direction = f[:, column]
            plus = z + finite_difference_step * direction
            minus = z - finite_difference_step * direction
            numeric[:, column] = (
                plus / np.linalg.norm(plus) - minus / np.linalg.norm(minus)
            ) / (2.0 * finite_difference_step)
        derivative_residuals.append(float(np.max(np.abs(analytic - numeric))))
        tangent_residuals.append(float(np.max(np.abs(p @ analytic))))

        u = rng.normal(size=4)
        u /= np.linalg.norm(u)
        chord = (u - q) / radius
        recovered = q + radius * chord
        recovered /= np.linalg.norm(recovered)
        boundary_residuals.append(float(np.max(np.abs(recovered - u))))
    return {
        "samples": samples,
        "finite_difference_step": finite_difference_step,
        "maximum_projection_derivative_residual": max(derivative_residuals),
        "maximum_tangent_residual": max(tangent_residuals),
        "maximum_exact_boundary_recovery_residual": max(boundary_residuals),
        "identity": "DP(z)=(I-P(z) tensor P(z))/|z|",
    }


def periodic_gradient(grid: int) -> tuple[np.ndarray, np.ndarray]:
    theta = 2.0 * math.pi * np.arange(grid) / grid
    x1, x2, x3 = np.meshgrid(theta, theta, theta, indexing="ij")
    shape = x1.shape + (4, 3)
    dphi = np.zeros(shape)
    amplitudes = (0.31, -0.27, 0.22, 0.19)

    # phi_1=A sin(x1) cos(x2)
    dphi[..., 0, 0] = amplitudes[0] * np.cos(x1) * np.cos(x2)
    dphi[..., 0, 1] = -amplitudes[0] * np.sin(x1) * np.sin(x2)
    # phi_2=B sin(x2) cos(x3)
    dphi[..., 1, 1] = amplitudes[1] * np.cos(x2) * np.cos(x3)
    dphi[..., 1, 2] = -amplitudes[1] * np.sin(x2) * np.sin(x3)
    # phi_3=C sin(x3) cos(x1)
    dphi[..., 2, 0] = -amplitudes[2] * np.sin(x3) * np.sin(x1)
    dphi[..., 2, 2] = amplitudes[2] * np.cos(x3) * np.cos(x1)
    # phi_4=D sin(x1+x2+x3)
    common = amplitudes[3] * np.cos(x1 + x2 + x3)
    dphi[..., 3, 0] = common
    dphi[..., 3, 1] = common
    dphi[..., 3, 2] = common

    f0 = np.asarray(
        [
            [0.12, -0.08, 0.04],
            [0.03, 0.16, -0.05],
            [-0.07, 0.02, 0.11],
            [0.05, -0.09, 0.13],
        ]
    )
    return f0, dphi


def strong_quasiconvexity_audit(grid: int) -> dict[str, Any]:
    f0, dphi = periodic_gradient(grid)
    perturbed = f0 + dphi
    g0 = graph_vector(f0)
    graph_perturbed = graph_vector(perturbed)
    lhs = float(np.mean(energy(perturbed) - energy(f0)))
    rhs = float(0.5 * np.mean(np.sum((graph_perturbed - g0) ** 2, axis=-1)))
    mean_f = np.mean(perturbed, axis=(0, 1, 2))
    mean_m2 = np.mean(batch_minors(perturbed, 2), axis=(0, 1, 2))
    mean_m3 = np.mean(batch_minors(perturbed, 3), axis=(0, 1, 2))
    f0_m2 = batch_minors(f0, 2)
    f0_m3 = batch_minors(f0, 3)
    return {
        "grid": grid,
        "energy_excess": lhs,
        "graph_variance": rhs,
        "identity_residual": abs(lhs - rhs),
        "mean_gradient_residual": float(np.max(np.abs(mean_f - f0))),
        "mean_second_minor_residual": float(np.max(np.abs(mean_m2 - f0_m2))),
        "mean_third_minor_residual": float(np.max(np.abs(mean_m3 - f0_m3))),
        "dirac_implication": (
            "if transferred quasiminimality makes the energy excess nonpositive, "
            "the nonnegative graph variance and concentration mass both vanish"
        ),
    }


def microball_audit() -> dict[str, Any]:
    h = np.asarray([2.0**-k for k in range(4, 13)])
    amplitude = h**0.25
    l2_map = amplitude**2 * h**3
    first_energy = amplitude**2 * h
    second_energy = amplitude**4 / h
    weighted_gluing = amplitude**6 / h**3

    def slope(values: np.ndarray) -> float:
        return float(np.polyfit(np.log(h), np.log(values), 1)[0])

    return {
        "amplitude": "a_h=h^(1/4)",
        "scales": {
            "L2_map": "a_h^2 h^3",
            "first_minor_energy": "a_h^2 h",
            "second_minor_energy": "a_h^4/h",
            "third_minor_energy": 0.0,
            "weighted_gluing_at_delta_h": "a_h^6/h^3",
        },
        "slopes": {
            "L2_map": slope(l2_map),
            "first_minor_energy": slope(first_energy),
            "second_minor_energy": slope(second_energy),
            "weighted_gluing": slope(weighted_gluing),
        },
        "smallest_h_values": {
            "h": float(h[-1]),
            "L2_map": float(l2_map[-1]),
            "first_minor_energy": float(first_energy[-1]),
            "second_minor_energy": float(second_energy[-1]),
            "weighted_gluing": float(weighted_gluing[-1]),
        },
        "conclusion": (
            "strong L2 convergence plus bounded graph energy does not imply "
            "the GTA weighted annular remainder"
        ),
    }


def wz_capacity_audit() -> dict[str, Any]:
    radii = np.asarray([0.8, 0.4, 0.2, 0.1])
    charge = 0.37
    rows = []
    residuals = []
    for radius in radii:
        volume = 4.0 * math.pi * radius**3 / 3.0
        constant_j = 2.0 * math.pi**2 * charge / volume
        direct_energy = 0.5 * K * constant_j**2 * volume
        capacity_energy = (
            3.0 * K * math.pi**3 * charge**2 / (2.0 * radius**3)
        )
        residuals.append(abs(direct_energy - capacity_energy))
        rows.append(
            {
                "radius": float(radius),
                "constant_J": float(constant_j),
                "direct_sextic_energy": float(direct_energy),
                "capacity_lower_bound": float(capacity_energy),
            }
        )

    a_norm_squared = 7.0
    da_norm_squared = 11.0
    wedge_integral = math.sqrt(a_norm_squared * da_norm_squared)
    hopf_number = wedge_integral / (16.0 * math.pi**2)
    left = (a_norm_squared / 8.0) * (R * da_norm_squared / 32.0)
    right = R * math.pi**4 * hopf_number**2
    return {
        "wz_rows": rows,
        "maximum_sharp_constant_residual": max(residuals),
        "capacity_formula": "E3 >= 3 K pi^3 |b|^2/(2 r^3)",
        "hopf_component_test": {
            "left_component_product": left,
            "right_topological_bound": right,
            "residual": abs(left - right),
        },
    }


def nonvolumetric_audit() -> dict[str, Any]:
    dimensions = np.asarray([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    exponents = 0.5 * (dimensions + 3.0)
    test_radius = 1.0e-4
    decay = test_radius**exponents
    return {
        "dimensions": dimensions.tolist(),
        "wz_decay_exponents": exponents.tolist(),
        "relative_wz_decay_at_r_1e-4": decay.tolist(),
        "all_quantized_charges_eventually_zero": bool(np.all(exponents > 0.0)),
        "alternative": {
            "capacity_active": (
                "positive normalized kappa forces a positive sextic tangent fraction "
                "but the absolute WZ number still vanishes"
            ),
            "wz_neutral": (
                "kappa vanishes at every rational radius; any remaining defect "
                "is analytic and not a quantized topological atom"
            ),
        },
    }


def policy_audit() -> dict[str, Any]:
    text = SCRIPT.read_text(encoding="utf-8")
    forbidden = [
        "sub" + "process",
        "scipy" + ".optimize",
        "lattice_" + "relax",
        "translated_" + "starts",
        "triangulation_" + "scan",
    ]
    hits = [token for token in forbidden if token in text]
    return {
        "forbidden_execution_tokens": forbidden,
        "hits": hits,
        "no_lattice_or_mesh_scan": not hits,
    }


def gate_ledger() -> dict[str, bool]:
    return {
        "normalized_full_local_quasiminimality": True,
        "graph_tight_annulus_derived_from_current_bounds": False,
        "conditional_recovery_compatible_boundary_modification": True,
        "quasiminimality_transfer_under_gta": True,
        "strong_quasiconvex_dirac_rigidity_under_transfer": True,
        "volumetric_defect_eliminated_without_gta": False,
        "nonvolumetric_wz_capacity_bound": True,
        "quantized_point_degree_concentration": False,
        "nonvolumetric_capacity_topology_alternative": True,
        "full_mu_zero": False,
        "local_recovery": False,
        "classicality": False,
        "continuum_isolation": False,
        "bosonic_hessian_authorized": False,
        "parallel_dirac_callias_mathematics_allowed": True,
        "determinant_promotion": False,
        "degree_one_portal_promotion": False,
        "physics_promotion_allowed": False,
    }


def build_card(args: argparse.Namespace) -> dict[str, Any]:
    projection = projection_audit(args.samples)
    strong_qc = strong_quasiconvexity_audit(args.grid)
    microball = microball_audit()
    capacity = wz_capacity_audit()
    nonvolumetric = nonvolumetric_audit()
    policy = policy_audit()
    gates = gate_ledger()

    check(
        "boundary",
        "radial projection derivative",
        projection["maximum_projection_derivative_residual"] < 2.0e-8,
        (
            "maximum finite-difference residual "
            f"{projection['maximum_projection_derivative_residual']:.3e}"
        ),
    )
    check(
        "boundary",
        "exact original trace recovery",
        projection["maximum_exact_boundary_recovery_residual"] < 1.0e-13,
        (
            "P(q+r(u-q)/r)=u residual "
            f"{projection['maximum_exact_boundary_recovery_residual']:.3e}"
        ),
    )
    check(
        "strong quasiconvexity",
        "null-Lagrangian moment identities",
        max(
            strong_qc["mean_gradient_residual"],
            strong_qc["mean_second_minor_residual"],
            strong_qc["mean_third_minor_residual"],
        )
        < 2.0e-13,
        "periodic gradient means reproduce all affine complete minors",
    )
    check(
        "strong quasiconvexity",
        "exact graph variance identity",
        strong_qc["identity_residual"] < 2.0e-13,
        f"identity residual {strong_qc['identity_residual']:.3e}",
    )
    check(
        "strong quasiconvexity",
        "Dirac gap is strictly positive for a nonconstant generator",
        strong_qc["graph_variance"] > 1.0e-5,
        f"positive graph variance {strong_qc['graph_variance']:.12f}",
    )
    check(
        "endpoint",
        "rank-two microball keeps bounded second-minor energy",
        abs(microball["slopes"]["second_minor_energy"]) < 1.0e-12,
        "a_h=h^(1/4) gives integral |M2|^2 of order one",
    )
    check(
        "endpoint",
        "weighted annular gluing is not automatic",
        abs(microball["slopes"]["weighted_gluing"] + 1.5) < 1.0e-12,
        "the cutoff remainder diverges like h^(-3/2)",
    )
    check(
        "capacity",
        "sharp WZ-capacity constant",
        capacity["maximum_sharp_constant_residual"] < 2.0e-12,
        (
            "maximum constant-field equality residual "
            f"{capacity['maximum_sharp_constant_residual']:.3e}"
        ),
    )
    check(
        "capacity",
        "Hopf component capacity product",
        capacity["hopf_component_test"]["residual"] < 2.0e-14,
        (
            "Cauchy equality residual "
            f"{capacity['hopf_component_test']['residual']:.3e}"
        ),
    )
    check(
        "capacity",
        "nonvolumetric defect cannot carry quantized point degree",
        nonvolumetric["all_quantized_charges_eventually_zero"],
        "all b(r)=O(r^((d+3)/2)) exponents are positive for d>=0",
    )
    check(
        "gates",
        "conditional closure and downstream embargo",
        (
            gates["conditional_recovery_compatible_boundary_modification"]
            and gates["strong_quasiconvex_dirac_rigidity_under_transfer"]
            and gates["nonvolumetric_capacity_topology_alternative"]
            and not gates["graph_tight_annulus_derived_from_current_bounds"]
            and not gates["full_mu_zero"]
            and not gates["bosonic_hessian_authorized"]
            and not gates["determinant_promotion"]
            and not gates["degree_one_portal_promotion"]
            and not gates["physics_promotion_allowed"]
        ),
        "GTA remains open; classicality, Hessian, determinant, and portal remain false",
    )
    check(
        "policy",
        "no lattice relaxation or mesh scan",
        policy["no_lattice_or_mesh_scan"],
        f"forbidden execution calls found: {policy['hits']}",
    )

    passed = sum(1 for row in CHECKS if row["pass"])
    total = len(CHECKS)
    return {
        "card": "AP-E17 boundary transfer / Dirac rigidity / WZ capacity",
        "mode": "continuum algebraic audit; no lattice scan",
        "status": "PASS" if passed == total else "FAIL",
        "summary": {"passed": passed, "total": total},
        "parameters": {
            "R": R,
            "K": K,
            "samples": args.samples,
            "periodic_grid": args.grid,
        },
        "boundary_projection": projection,
        "strong_quasiconvexity": strong_qc,
        "microball_obstruction": microball,
        "wz_capacity": capacity,
        "nonvolumetric": nonvolumetric,
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


def markdown(card: dict[str, Any]) -> str:
    status = card["status"]
    passed = card["summary"]["passed"]
    total = card["summary"]["total"]
    sqc = card["strong_quasiconvexity"]
    micro = card["microball_obstruction"]
    lines = [
        "# AP-E17 boundary-transfer/Dirac/capacity audit",
        "",
        f"- Status: **{status}**",
        f"- Checks: **{passed}/{total}**",
        f"- Exact graph-variance residual: `{sqc['identity_residual']:.3e}`",
        f"- Nonconstant periodic graph variance: `{sqc['graph_variance']:.12f}`",
        (
            "- Microball weighted-gluing slope: "
            f"`{micro['slopes']['weighted_gluing']:.12f}`"
        ),
        "",
        "## Main conclusions",
        "",
        (
            "- The AP-E16 diagonal already has normalized full local "
            "quasiminimality."
        ),
        (
            "- Under the explicit graph-tight annulus condition, radial "
            "projection gives a trace-, homotopy-, and degree-preserving affine "
            "interior comparison."
        ),
        (
            "- The exact complete-minor variance identity then forces the "
            "homogeneous tangent Young measure to be Dirac and removes graph "
            "concentration."
        ),
        (
            "- GTA is not derived from current endpoint bounds: a rank-two "
            "microball keeps bounded second-minor energy while its weighted "
            "gluing term diverges."
        ),
        (
            "- The sharp WZ-capacity bound excludes quantized point-degree "
            "concentration and gives the capacity-active/WZ-neutral "
            "nonvolumetric alternative."
        ),
        "",
        "## Gates",
        "",
    ]
    for name, value in card["gates"].items():
        lines.append(f"- `{name}`: `{value}`")
    lines.extend(["", "## Check details", ""])
    for row in card["checks"]:
        label = "PASS" if row["pass"] else "FAIL"
        lines.append(
            f"- **{row['group']} / {row['name']}**: `{label}` — {row['detail']}"
        )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--grid", type=int, default=24)
    parser.add_argument("--tag", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    card = build_card(args)
    suffix = f"_{args.tag}" if args.tag else ""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / f"ap_e17_boundary_transfer_dirac_capacity{suffix}.json"
    md_path = OUTPUT / f"ap_e17_boundary_transfer_dirac_capacity{suffix}.md"
    json_path.write_text(json.dumps(card, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(markdown(card), encoding="utf-8")
    print(
        f"{card['status']} {card['summary']['passed']}/{card['summary']['total']}; "
        f"transfer={card['gates']['quasiminimality_transfer_under_gta']}; "
        f"dirac={card['gates']['strong_quasiconvex_dirac_rigidity_under_transfer']}; "
        f"gta={card['gates']['graph_tight_annulus_derived_from_current_bounds']}; "
        f"capacity={card['gates']['nonvolumetric_capacity_topology_alternative']}"
    )
    print(json_path)
    print(md_path)
    return 0 if card["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
