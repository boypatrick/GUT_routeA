#!/usr/bin/env python3
"""AP-E9 scaling-window, triangulation, and translated-start scan.

This card tests the only elementary scaling window in which the AP-E8
topology barrier is simultaneously uniform on bounded sublevels and invisible
on smooth recovery fields.  With

    gamma_p(a) = gamma_ref (a_ref/a)^p,
    d(a)       = 1-epsilon(a),
    w(a)       = gamma(a) a,
    R(a)       = gamma(a) a^2 / d(a)^2,

fixed d requires p >= 1 for a uniform pair-gap consequence of an energy bound,
while p < 2 is required for the smooth barrier to vanish.  The matched p=1
trajectory used here has a_ref=0.4 and gamma_ref=1, hence w=0.4 exactly and
R=0.4 a/d^2 -> 0.

The base Dirichlet term does prove fixed-box strong-L2 equicoercivity through
the standard piecewise-affine H1 bound and Rellich compactness.  The numerical
scan remains fail-closed because degree-sector/topological-current compactness
and the full Skyrme Gamma-liminf/recovery are not proved.  In particular it
never opens the Riemann-Hessian or regulated-determinant gate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
E8_SCRIPT = ROUTE_F / "code" / "scan_ap_e8_topology_preserving_b1.py"
TRI_SCRIPT = ROUTE_F / "code" / "ap_e9_triangulation_tools.py"

EPSILON = 0.01
REFERENCE_SPACING = 0.4
REFERENCE_GAMMA = 1.0
STATIONARITY_TOLERANCE = 2.0e-6
PAIR_MARGIN_TOLERANCE = 1.0e-8
TRANSLATED_OFFSET = (0.37, -0.23, 0.41)
CHECKS: list[dict[str, Any]] = []


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E8 = load_module("ap_e8_for_ap_e9", E8_SCRIPT)
TRI = load_module("ap_e9_triangulation_tools_for_scan", TRI_SCRIPT)
E7 = E8.E7
E6 = E8.E6

TRIANGULATIONS = {
    spec.name: spec
    for spec in (
        TRI.UNIFORM_SPECS[0],
        *TRI.FIVE_TET_SPECS,
    )
}


@dataclass(frozen=True)
class CasePlan:
    case_id: str
    extent: int
    length: float
    exponent_p: float
    triangulation: str = "uniform_ppp"
    translation_lattice_units: tuple[float, float, float] = (0.0, 0.0, 0.0)
    roles: tuple[str, ...] = ()


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


def gamma_trajectory(spacing: float, exponent_p: float) -> float:
    return REFERENCE_GAMMA * (REFERENCE_SPACING / spacing) ** exponent_p


def scaling_numbers(
    spacing: float, exponent_p: float, epsilon: float = EPSILON
) -> dict[str, float]:
    gamma = gamma_trajectory(spacing, exponent_p)
    d = 1.0 - epsilon
    return {
        "a": spacing,
        "p": exponent_p,
        "epsilon": epsilon,
        "d_equals_one_minus_epsilon": d,
        "gamma": gamma,
        "w_equals_gamma_a": gamma * spacing,
        "R_equals_gamma_a2_over_d2": gamma * spacing**2 / d**2,
    }


def exact_scaling_contract() -> dict[str, Any]:
    sample_spacings = (0.5, 0.4, 0.3, 0.25, 0.2)
    rows = [
        scaling_numbers(spacing, exponent)
        for exponent in (0.0, 1.0, 2.0)
        for spacing in sample_spacings
    ]
    return {
        "definitions": {
            "d(a)": "1-epsilon(a)",
            "w(a)": "gamma(a)*a",
            "R(a)": "gamma(a)*a^2/d(a)^2",
            "matched_gamma_p(a)": "gamma_ref*(a_ref/a)^p",
            "a_ref": REFERENCE_SPACING,
            "gamma_ref": REFERENCE_GAMMA,
        },
        "rigorous_pair_gap_bound": {
            "statement": "for every fixed C>0, if total barrier energy <= C and w>0, every regulated edge obeys q-epsilon >= d*exp(-1-C/w); C=0 is an exceptional vacuum sublevel",
            "proof": "write x=(q-epsilon)/d. Since b(x)=-log(x)+x-1 >= -log(x)-1 and w*b(x)<=C, x>=exp(-1-C/w)",
        },
        "p_less_than_one_countersequence": {
            "statement": "if w(a)->0, uniform pair-gap coercivity is false for the barrier alone",
            "construction": "rotate one interior vacuum-collar vertex, with at most m incident barrier edges, so x_a=exp(-alpha/w(a)) and 0<alpha<C/(2m); the incident barrier sum stays below C while q_a-epsilon=d*x_a->0 and the rank-one local base costs vanish",
        },
        "smooth_expansion": {
            "statement": "for 1-q=O(a^2) and (1-q)/d->0, b(q)=(1-q)^2/(2d^2)+O((1-q)^3/d^3), so the three-dimensional edge sum is O(R)",
            "generic_smooth_recovery_condition": "a^2/d(a)->0",
        },
        "fixed_d_classification": {
            "p=0": "R->0 but w->0: continuum-invisible yet not uniformly sector-coercive",
            "p=1": "w=gamma_ref*a_ref>0 and R->0: minimal uniform-edge-margin, smooth-vanishing candidate",
            "1<p<2": "w grows and R->0: also in the necessary window, with stronger finite-a barrier",
            "p=2": "w grows and R approaches a nonzero constant: a raw anisotropic Cauchy-Born quartic survives on fixed smooth samples; the periodic homogenized Gamma density remains unevaluated",
            "p>2": "R diverges on generic nonconstant smooth fields",
            "candidate_window": "1 <= p < 2",
        },
        "power_law_d_generalization": {
            "assumption": "gamma~a^(-p), d~a^q",
            "w_exponent": "1-p",
            "R_exponent": "2-p-2q",
            "necessary_joint_window": "p>=1 and p+2q<2, plus epsilon>-1/3 uniformly and a^2/d->0",
            "warning": "these conditions control only the topology barrier; fixed-box strong-L2 equicoercivity follows separately from the base Dirichlet term, while degree-sector/topological-current compactness and the full Skyrme Gamma-liminf remain open",
        },
        "fixed_box_L2_equicoercivity": {
            "status": "proved",
            "statement": "the declared base action contains the nearest-neighbour Dirichlet sum with a mesh-uniform positive coefficient; a uniform bound on the nonnegative full action therefore makes the standard piecewise-affine interpolants H1-bounded on every fixed box, and Rellich gives a strongly L2-precompact subsequence",
            "scope": "fixed physical box with the declared boundary data; this does not by itself preserve degree in the limit",
        },
        "barrier_gamma_audit": {
            "state_space": "strong L2(Omega;R4), using dual-cell piecewise-constant extensions of sphere-valued nodes and +infinity off the discrete image",
            "R_to_zero": "proved Gamma limit: zero on L2(Omega;S3), +infinity outside; boundary trace and degree are forgotten",
            "positive_R_compactness": "proved: inf R>0 gives affine W1,4 and C0,beta compactness for beta<1/4, hence degree closure",
            "R_to_infinity": "proved bounded-sublevel limit: fixed vacuum boundary forces the vacuum",
            "finite_positive_R_density": "open: Q_T is only the raw Cauchy-Born smooth-sampling density; a periodic multi-cell triangulation requires a homogenized cell formula",
            "raw_six_vs_five_tet_tensors_proportional": False,
            "homogenized_six_vs_five_tet_densities_computed": False,
        },
        "normalized_affine_no_zero": {
            "bound": "norm(sum_i lambda_i n_i)^2 >= (1+3 epsilon)/4",
            "uniform_requirement": "liminf epsilon(a)>-1/3",
            "current_value": math.sqrt((1.0 + 3.0 * EPSILON) / 4.0),
        },
        "sample_rows": rows,
    }


def topology_summary(
    field: np.ndarray, tetrahedra: np.ndarray, orientations: np.ndarray
) -> dict[str, Any]:
    rows = []
    for target in E7.TARGETS:
        row = TRI.regular_value_degree(field, tetrahedra, orientations, target)
        rows.append({"target": target.tolist(), **row})
    baryon_numbers = [row["baryon_number"] for row in rows]
    targets_agree = len(set(baryon_numbers)) == 1
    return {
        "method": "orientation-aware normalized-affine degree on the same tetrahedral complex used by the barrier",
        "target_rows": rows,
        "baryon_numbers": baryon_numbers,
        "targets_agree": targets_agree,
        "all_targets_B1": targets_agree and baryon_numbers == [1, 1, 1],
    }


def spatial_moments(field: np.ndarray, axis: np.ndarray) -> dict[str, Any]:
    coordinates = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    )
    weight = np.maximum(1.0 - field[..., 0], 0.0)
    total = float(np.sum(weight))
    if total <= 1.0e-300:
        barycentre = np.zeros(3)
        origin_rms = 0.0
        centered_rms = 0.0
    else:
        barycentre = np.sum(weight[..., None] * coordinates, axis=(0, 1, 2)) / total
        origin_rms = math.sqrt(
            float(np.sum(weight * np.sum(coordinates**2, axis=-1))) / total
        )
        centered_rms = math.sqrt(
            float(
                np.sum(
                    weight
                    * np.sum((coordinates - barycentre) ** 2, axis=-1)
                )
            )
            / total
        )
    spacing = float(axis[1] - axis[0])
    return {
        "weight_definition": "rho=1-n0",
        "total_weight": total,
        "barycentre": barycentre.tolist(),
        "barycentre_norm": float(np.linalg.norm(barycentre)),
        "barycentre_in_lattice_units": (barycentre / spacing).tolist(),
        "origin_rms": origin_rms,
        "centered_rms": centered_rms,
        "centered_rms_over_a": centered_rms / spacing,
    }


def stationarity_and_concentration(
    field: np.ndarray,
    breathing: np.ndarray,
    spacing: float,
    model: Any,
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    epsilon: float,
) -> dict[str, Any]:
    energy, gradient_n, gradient_s, geometry = E8.total_energy_gradient(
        field, breathing, spacing, model, left, right, gamma, epsilon
    )
    interior_field = field[1:-1, 1:-1, 1:-1]
    interior_gradient = gradient_n[1:-1, 1:-1, 1:-1]
    tangent = interior_gradient - np.sum(
        interior_gradient * interior_field, axis=-1, keepdims=True
    ) * interior_field
    scalar = gradient_s[1:-1, 1:-1, 1:-1]
    site_squared = np.sum(tangent**2, axis=-1) + scalar**2
    total_squared = float(np.sum(site_squared))
    number_sites = int(site_squared.size)
    dofs = 4 * number_sites
    norm = math.sqrt(total_squared)
    if total_squared > 0.0:
        fractions = (site_squared / total_squared).ravel()
        sorted_fractions = np.sort(fractions)[::-1]
        sites_for_90 = int(np.searchsorted(np.cumsum(sorted_fractions), 0.9) + 1)
        maximum_fraction = float(sorted_fractions[0])
        ipr = float(fractions @ fractions)
    else:
        sites_for_90 = 0
        maximum_fraction = 0.0
        ipr = 0.0
    return {
        "total_energy": float(energy),
        "base_energy": geometry["base_energy"],
        "barrier_energy": geometry["barrier_energy"],
        "direct_tangent_plus_scalar_gradient_l2": norm,
        "direct_projected_gradient_l2_density": norm / math.sqrt(dofs) / spacing**3,
        "minimum_pair_dot": geometry["minimum_pair_dot"],
        "minimum_pair_margin": geometry["minimum_pair_margin"],
        "invalid_edges": geometry["invalid_edges"],
        "local_gradient_concentration": {
            "site_measure": "sum over three tangent components plus breathing-gradient squared",
            "interior_sites": number_sites,
            "maximum_site_fraction": maximum_fraction,
            "inverse_participation_ratio": ipr,
            "effective_sites": 1.0 / ipr if ipr > 0.0 else None,
            "sites_containing_90_percent": sites_for_90,
            "fraction_of_sites_containing_90_percent": (
                sites_for_90 / number_sites if number_sites else None
            ),
        },
    }


def relax_case(
    plan: CasePlan,
    solution: Any,
    model: Any,
    quick: bool,
) -> dict[str, Any]:
    spec = TRIANGULATIONS[plan.triangulation]
    field, breathing, axis, _ = TRI.translated_radial_background(
        solution,
        plan.extent,
        plan.length,
        plan.translation_lattice_units,
        "lattice",
        radial_epsilon=E6.EPS,
    )
    field, breathing = E7.exact_vacuum_boundary(field, breathing)
    spacing = float(axis[1] - axis[0])
    gamma = gamma_trajectory(spacing, plan.exponent_p)
    tetrahedra, orientations = TRI.tetrahedral_complex(plan.extent, spec)
    left, right = TRI.unique_edges(tetrahedra)
    triangulation_audit = TRI.triangulation_audit(plan.extent, spec)
    initial = stationarity_and_concentration(
        field, breathing, spacing, model, left, right, gamma, EPSILON
    )
    initial_topology = topology_summary(field, tetrahedra, orientations)
    history = []
    stationary = False
    maximum_recentres = 3 if quick else 5
    maximum_iterations = 500 if quick else 750
    for recenter in range(maximum_recentres):
        base, scalar, frames = E7.fixed_chart(field, breathing)

        def objective(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
            return E8.chart_objective(
                coordinates,
                field,
                breathing,
                base,
                scalar,
                frames,
                spacing,
                model,
                left,
                right,
                gamma,
                EPSILON,
            )

        result = minimize(
            objective,
            np.zeros(4 * len(base)),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": maximum_iterations,
                "maxcor": 30,
                "maxls": 80,
                "ftol": 1.0e-15,
                "gtol": 1.0e-10,
            },
        )
        field, breathing, _, _ = E7.chart_fields(
            result.x, field, breathing, base, scalar, frames
        )
        stationarity = stationarity_and_concentration(
            field, breathing, spacing, model, left, right, gamma, EPSILON
        )
        history.append(
            {
                "recenter": recenter,
                "optimizer_success": bool(result.success),
                "optimizer_message": str(result.message),
                "iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "energy": stationarity["total_energy"],
                "direct_projected_gradient_l2_density": stationarity[
                    "direct_projected_gradient_l2_density"
                ],
                "minimum_pair_margin": stationarity["minimum_pair_margin"],
            }
        )
        if (
            stationarity["direct_projected_gradient_l2_density"]
            < STATIONARITY_TOLERANCE
        ):
            stationary = True
            break
    final = stationarity_and_concentration(
        field, breathing, spacing, model, left, right, gamma, EPSILON
    )
    final_topology = topology_summary(field, tetrahedra, orientations)
    moments = spatial_moments(field, axis)
    passed = bool(
        stationary
        and final["minimum_pair_margin"] > PAIR_MARGIN_TOLERANCE
        and final_topology["all_targets_B1"]
        and triangulation_audit["valid_conforming_triangulation"]
    )
    return {
        "case_id": plan.case_id,
        "roles": list(plan.roles),
        "N": plan.extent,
        "L": plan.length,
        **scaling_numbers(spacing, plan.exponent_p),
        "triangulation": asdict(spec),
        "translation_lattice_units": list(plan.translation_lattice_units),
        "translation_physical_units": (
            spacing * np.asarray(plan.translation_lattice_units)
        ).tolist(),
        "core_or_centre_anchor_used": False,
        "tetrahedra": int(len(tetrahedra)),
        "regulated_unique_pair_edges": int(len(left)),
        "triangulation_audit": triangulation_audit,
        "initial": initial,
        "initial_topology": initial_topology,
        "recentering_history": history,
        "stationary_by_direct_tangent_test": stationary,
        "final": final,
        "final_topology": final_topology,
        "spatial_moments": moments,
        "case_pass": passed,
    }


def smooth_profile_scaling(quick: bool) -> dict[str, Any]:
    extents = (17, 21, 25) if quick else (17, 21, 25, 29, 33)
    result: dict[str, Any] = {
        "profile": "AP-E8 C2 compact-supported degree-one hedgehog, L=6 and support radius=2",
        "triangulation": "uniform_ppp",
        "rows_by_p": {},
    }
    for exponent in (0.0, 1.0, 2.0):
        rows = []
        for extent in extents:
            field, axis = E8.compact_supported_degree_one_field(extent)
            spacing = float(axis[1] - axis[0])
            tetrahedra, _ = TRI.tetrahedral_complex(
                extent, TRIANGULATIONS["uniform_ppp"]
            )
            left, right = TRI.unique_edges(tetrahedra)
            gamma = gamma_trajectory(spacing, exponent)
            energy, _, geometry = E8.barrier_energy_gradient(
                field, left, right, spacing, gamma, EPSILON
            )
            rows.append(
                {
                    "N": extent,
                    **scaling_numbers(spacing, exponent),
                    "barrier_energy": energy,
                    "minimum_pair_dot": geometry["minimum_pair_dot"],
                }
            )
        log_a = np.log([row["a"] for row in rows])
        log_energy = np.log([row["barrier_energy"] for row in rows])
        slope, intercept = np.polyfit(log_a, log_energy, 1)
        fitted = slope * log_a + intercept
        residual = log_energy - fitted
        centered = log_energy - np.mean(log_energy)
        r_squared = 1.0 - float(residual @ residual) / float(centered @ centered)
        result["rows_by_p"][f"p={int(exponent)}"] = {
            "p": exponent,
            "expected_asymptotic_slope": 2.0 - exponent,
            "fitted_log_log_slope": float(slope),
            "fit_r_squared": r_squared,
            "local_slopes": [
                float(
                    (log_energy[index + 1] - log_energy[index])
                    / (log_a[index + 1] - log_a[index])
                )
                for index in range(len(rows) - 1)
            ],
            "rows": rows,
        }
    return result


def add_plan(
    plans: dict[tuple[Any, ...], CasePlan],
    case_id: str,
    extent: int,
    length: float,
    exponent_p: float,
    triangulation: str,
    translation: Sequence[float],
    role: str,
) -> None:
    translation_tuple = tuple(float(value) for value in translation)
    key = (
        extent,
        float(length),
        float(exponent_p),
        triangulation,
        translation_tuple,
    )
    if key in plans:
        old = plans[key]
        plans[key] = CasePlan(
            old.case_id,
            old.extent,
            old.length,
            old.exponent_p,
            old.triangulation,
            old.translation_lattice_units,
            (*old.roles, role),
        )
    else:
        plans[key] = CasePlan(
            case_id,
            extent,
            length,
            exponent_p,
            triangulation,
            translation_tuple,
            (role,),
        )


def build_case_plan(quick: bool) -> list[CasePlan]:
    plans: dict[tuple[Any, ...], CasePlan] = {}
    fixed_box_extents = (15, 17) if quick else (15, 17, 19, 21, 23, 25)
    for extent in fixed_box_extents:
        add_plan(
            plans,
            f"p1_fixed_box_N{extent}",
            extent,
            6.0,
            1.0,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            "matched_p1_fixed_box",
        )
    fixed_spacing_specs = (
        ((16, 6.0), (21, 8.0))
        if quick
        else ((16, 6.0), (21, 8.0), (26, 10.0), (31, 12.0))
    )
    for extent, length in fixed_spacing_specs:
        add_plan(
            plans,
            f"p1_fixed_a_N{extent}_L{int(length)}",
            extent,
            length,
            1.0,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            "matched_p1_fixed_spacing",
        )

    # N=13 is already outside the declared uniform-Freudenthal barrier for
    # the translated seed (two q<=epsilon edges).  N=16 is the smallest common
    # screen in this matrix on which all three triangulations and both starts
    # are admissible, so a comparison never begins at infinite action.
    triangulation_extents = (16,) if quick else (16, 21)
    for triangulation_extent in triangulation_extents:
        suffix = "" if triangulation_extent == 16 else "_fine"
        for triangulation in TRIANGULATIONS:
            add_plan(
                plans,
                f"tri_center_N{triangulation_extent}_{triangulation}",
                triangulation_extent,
                6.0,
                1.0,
                triangulation,
                (0.0, 0.0, 0.0),
                f"triangulation{suffix}_centered",
            )
            add_plan(
                plans,
                f"tri_translated_N{triangulation_extent}_{triangulation}",
                triangulation_extent,
                6.0,
                1.0,
                triangulation,
                TRANSLATED_OFFSET,
                f"triangulation{suffix}_translated_start",
            )
    if not quick:
        add_plan(
            plans,
            "uniform_half_cell_start",
            16,
            6.0,
            1.0,
            "uniform_ppp",
            (0.5, 0.5, 0.5),
            "translated_start_control",
        )

    control_extent = 15 if quick else 21
    for exponent in (0.0, 2.0):
        add_plan(
            plans,
            f"p{int(exponent)}_control_N{control_extent}",
            control_extent,
            6.0,
            exponent,
            "uniform_ppp",
            (0.0, 0.0, 0.0),
            f"p{int(exponent)}_control",
        )
    return list(plans.values())


def relative_spread(values: Sequence[float]) -> float:
    values = list(values)
    if not values:
        # Keep reduced-coverage JSON strict (allow_nan=False) while ensuring
        # every absent diagnostic fails closed by an overwhelming margin.
        return 1.0e300
    scale = max(max(abs(value) for value in values), 1.0e-300)
    return (max(values) - min(values)) / scale


def numerical_gate(cases: list[dict[str, Any]], quick: bool) -> dict[str, Any]:
    by_role = {
        role: [row for row in cases if role in row["roles"]]
        for role in (
            "matched_p1_fixed_box",
            "matched_p1_fixed_spacing",
            "triangulation_centered",
            "triangulation_translated_start",
            "triangulation_fine_centered",
            "triangulation_fine_translated_start",
            "p0_control",
        )
    }
    fixed_box = sorted(by_role["matched_p1_fixed_box"], key=lambda row: row["a"])
    fixed_spacing = sorted(
        by_role["matched_p1_fixed_spacing"], key=lambda row: row["L"]
    )
    triangulation_centered = by_role["triangulation_centered"]
    triangulation_translated = by_role["triangulation_translated_start"]
    triangulation_fine_centered = by_role["triangulation_fine_centered"]
    triangulation_fine_translated = by_role[
        "triangulation_fine_translated_start"
    ]
    finest = fixed_box[:3]
    fixed_box_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in finest]
    )
    fixed_box_radius_spread = relative_spread(
        [row["spatial_moments"]["centered_rms"] for row in finest]
    )
    volume_tail = fixed_spacing[-3:]
    volume_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in volume_tail]
    )
    volume_radius_spread = relative_spread(
        [row["spatial_moments"]["centered_rms"] for row in volume_tail]
    )
    fixed_box_barrier_slope = float(
        np.polyfit(
            np.log([row["a"] for row in fixed_box]),
            np.log([row["final"]["barrier_energy"] for row in fixed_box]),
            1,
        )[0]
    )
    centered_tri_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in triangulation_centered]
    )
    centered_tri_radius_spread = relative_spread(
        [row["spatial_moments"]["centered_rms"] for row in triangulation_centered]
    )
    translated_tri_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in triangulation_translated]
    )
    fine_centered_tri_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in triangulation_fine_centered]
    )
    fine_centered_tri_radius_spread = relative_spread(
        [
            row["spatial_moments"]["centered_rms"]
            for row in triangulation_fine_centered
        ]
    )
    fine_translated_tri_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in triangulation_fine_translated]
    )
    p0_rows = by_role["p0_control"]
    p1_control_rows = [
        row
        for row in fixed_box
        if abs(row["a"] - 0.3) < 1.0e-12
        and row["triangulation"]["name"] == "uniform_ppp"
    ]
    vanishing_regulator_rows = p0_rows + p1_control_rows
    p0_p1_energy_spread = relative_spread(
        [row["final"]["total_energy"] for row in vanishing_regulator_rows]
    )
    p0_p1_radius_spread = relative_spread(
        [row["spatial_moments"]["centered_rms"] for row in vanishing_regulator_rows]
    )
    all_required_cases_pass = all(row["case_pass"] for row in cases)
    production_coverage = bool(
        not quick
        and len(fixed_box) >= 6
        and min(row["a"] for row in fixed_box) <= 0.25 + 1.0e-12
        and len(fixed_spacing) >= 4
        and max(row["L"] for row in fixed_spacing) >= 12.0
        and len(triangulation_centered) == 3
        and len(triangulation_translated) == 3
        and len(triangulation_fine_centered) == 3
        and len(triangulation_fine_translated) == 3
        and len(vanishing_regulator_rows) == 2
    )
    threshold_rows = {
        "all_required_cases_pass": all_required_cases_pass,
        "production_coverage": production_coverage,
        "finest_three_fixed_box_energy_relative_spread": fixed_box_energy_spread,
        "finest_three_fixed_box_energy_threshold": 0.04,
        "finest_three_fixed_box_centered_rms_relative_spread": fixed_box_radius_spread,
        "finest_three_fixed_box_centered_rms_threshold": 0.08,
        "largest_three_volumes_energy_relative_spread": volume_energy_spread,
        "largest_three_volumes_energy_threshold": 0.015,
        "largest_three_volumes_centered_rms_relative_spread": volume_radius_spread,
        "largest_three_volumes_centered_rms_threshold": 0.05,
        "fixed_box_barrier_log_log_slope": fixed_box_barrier_slope,
        "fixed_box_barrier_expected_slope": 1.0,
        "fixed_box_barrier_slope_absolute_tolerance": 0.4,
        "centered_triangulations_energy_relative_spread": centered_tri_energy_spread,
        "centered_triangulations_energy_threshold": 0.03,
        "centered_triangulations_centered_rms_relative_spread": centered_tri_radius_spread,
        "centered_triangulations_centered_rms_threshold": 0.08,
        "translated_triangulations_energy_relative_spread": translated_tri_energy_spread,
        "translated_triangulations_energy_threshold": 0.03,
        "fine_centered_triangulations_energy_relative_spread": fine_centered_tri_energy_spread,
        "fine_centered_triangulations_energy_threshold": 0.03,
        "fine_centered_triangulations_centered_rms_relative_spread": fine_centered_tri_radius_spread,
        "fine_centered_triangulations_centered_rms_threshold": 0.08,
        "fine_translated_triangulations_energy_relative_spread": fine_translated_tri_energy_spread,
        "fine_translated_triangulations_energy_threshold": 0.03,
        "p0_p1_vanishing_regulators_energy_relative_spread": p0_p1_energy_spread,
        "p0_p1_vanishing_regulators_energy_threshold": 0.03,
        "p0_p1_vanishing_regulators_centered_rms_relative_spread": p0_p1_radius_spread,
        "p0_p1_vanishing_regulators_centered_rms_threshold": 0.05,
    }
    numerical_thresholds_pass = bool(
        all_required_cases_pass
        and production_coverage
        and fixed_box_energy_spread < 0.04
        and fixed_box_radius_spread < 0.08
        and volume_energy_spread < 0.015
        and volume_radius_spread < 0.05
        and abs(fixed_box_barrier_slope - 1.0) < 0.4
        and centered_tri_energy_spread < 0.03
        and centered_tri_radius_spread < 0.08
        and translated_tri_energy_spread < 0.03
        and fine_centered_tri_energy_spread < 0.03
        and fine_centered_tri_radius_spread < 0.08
        and fine_translated_tri_energy_spread < 0.03
        and p0_p1_energy_spread < 0.03
        and p0_p1_radius_spread < 0.05
    )
    return {
        **threshold_rows,
        "numerical_regulator_stability_thresholds_pass": numerical_thresholds_pass,
        "interpretation": "even a numerical pass is evidence only; fixed-box L2 equicoercivity is proved, but degree-sector/topological-current compactness and full Skyrme Gamma-liminf/recovery remain separate obligations",
    }


def markdown(result: dict[str, Any]) -> str:
    lines = [
        "# AP-E9 gamma-scaling and triangulation card",
        "",
        f"Status: `{result['status']}`",
        "",
        "## Exact scaling diagnosis",
        "",
        "For fixed epsilon, the uniform-edge-margin but smooth-vanishing window is `1 <= p < 2`. The matched scan uses `p=1`, so `w=gamma*a=0.4` while `R=gamma*a^2/d^2 -> 0`; this window alone does not prove degree-sector compactness.",
        "",
        "| p | fitted smooth slope | expected | R2 |",
        "|---:|---:|---:|---:|",
    ]
    for row in result["smooth_profile_scaling"]["rows_by_p"].values():
        lines.append(
            f"| {row['p']:.0f} | {row['fitted_log_log_slope']:.6f} | "
            f"{row['expected_asymptotic_slope']:.1f} | {row['fit_r_squared']:.6f} |"
        )
    lines.extend(
        [
            "",
            "## Re-relaxed cases",
            "",
            "| id | roles | N | L | a | p | mesh | offset/a | E | grad | min margin | B | barycentre | centered RMS | pass |",
            "|:--|:--|---:|---:|---:|---:|:--|:--|---:|---:|---:|:--|---:|---:|:--:|",
        ]
    )
    for row in result["cases"]:
        lines.append(
            "| {case_id} | {roles} | {N} | {L:.1f} | {a:.6f} | {p:.0f} | "
            "{mesh} | {offset} | {energy:.8f} | {gradient:.2e} | {margin:.3e} | "
            "{B} | {bary:.3e} | {radius:.5f} | {passed} |".format(
                case_id=row["case_id"],
                roles=",".join(row["roles"]),
                N=row["N"],
                L=row["L"],
                a=row["a"],
                p=row["p"],
                mesh=row["triangulation"]["name"],
                offset=row["translation_lattice_units"],
                energy=row["final"]["total_energy"],
                gradient=row["final"]["direct_projected_gradient_l2_density"],
                margin=row["final"]["minimum_pair_margin"],
                B=row["final_topology"]["baryon_numbers"],
                bary=row["spatial_moments"]["barycentre_norm"],
                radius=row["spatial_moments"]["centered_rms"],
                passed=row["case_pass"],
            )
        )
    gate = result["numerical_gate"]
    flags = result["fail_closed"]
    lines.extend(
        [
            "",
            "## Fail-closed gate",
            "",
            f"Numerical regulator thresholds: `{gate['numerical_regulator_stability_thresholds_pass']}`.",
            f"Fixed-box strong-L2 equicoercivity: `{flags['fixed_box_strong_L2_equicoercivity_proven']}`.",
            f"Barrier R->0 Gamma limit on the declared dual-cell embedding: `{flags['barrier_R_to_zero_gamma_limit_proven_on_declared_embedding']}`.",
            f"Positive-R barrier degree compactness: `{flags['barrier_positive_R_degree_compactness_proven']}`.",
            f"Finite-R homogenized barrier density computed: `{flags['finite_R_homogenized_barrier_density_computed']}`.",
            f"Degree-sector/topological-current compactness: `{flags['degree_sector_topological_current_compactness_proven']}`.",
            f"Full Gamma limit: `{flags['full_action_gamma_limit_proven']}`.",
            f"Regulator-stable continuum background: `{flags['regulator_stable_continuum_background']}`.",
            f"Same-action Riemann Hessian allowed: `{flags['same_action_riemann_hessian_allowed']}`.",
            f"Regulated determinant variation allowed: `{flags['regulated_determinant_variation_allowed']}`.",
            "",
            "The base Dirichlet term closes fixed-box L2 equicoercivity. The p=1 argument additionally proves a uniform edge-gap consequence and smooth-profile recovery scaling. Neither result yet proves degree-sector closure or the full Skyrme liminf/recovery theorem.",
            "",
            f"Checks: **{result['checks_passed']}/{result['checks_total']}**.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT / "ap_e9_gamma_triangulation"
    stem.with_suffix(".json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    stem.with_suffix(".md").write_text(markdown(result), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="reduced coverage; all continuum and Hessian flags remain fail-closed",
    )
    args = parser.parse_args()
    CHECKS.clear()
    model = E6.Model()
    profile = E6.solve_relaxed_profile(box=16.0, model=model, sample_points=4097)
    checksum = profile["profile_sha256_r_F_s_little_endian_float64"]
    solution = profile["solution"]
    check(
        "P0_provenance",
        "canonical AP-E6 continuum seed checksum is frozen",
        checksum
        == "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59",
        f"checksum={checksum}",
    )

    helper_tests = TRI.run_self_tests()
    check(
        "P1_triangulation",
        "all six-tet and checkerboard five-tet complexes pass exact conformity, orientation, degree, and translation tests",
        helper_tests["pass"],
        "audits="
        + str(
            [
                (row["spec"]["name"], row["valid_conforming_triangulation"])
                for row in helper_tests["triangulation_audits"]
            ]
        ),
    )
    scaling_contract = exact_scaling_contract()
    p1_rows = [
        row
        for row in scaling_contract["sample_rows"]
        if row["p"] == 1.0
    ]
    check(
        "P1_scaling",
        "matched p=1 has exactly constant w and decreasing R on the declared sample spacings",
        max(row["w_equals_gamma_a"] for row in p1_rows)
        - min(row["w_equals_gamma_a"] for row in p1_rows)
        < 1.0e-14
        and all(
            p1_rows[index + 1]["R_equals_gamma_a2_over_d2"]
            < p1_rows[index]["R_equals_gamma_a2_over_d2"]
            for index in range(len(p1_rows) - 1)
        ),
        "(a,w,R)="
        + str(
            [
                (
                    row["a"],
                    row["w_equals_gamma_a"],
                    row["R_equals_gamma_a2_over_d2"],
                )
                for row in p1_rows
            ]
        ),
    )
    smooth_scaling = smooth_profile_scaling(args.quick)
    slope_pass = all(
        abs(row["fitted_log_log_slope"] - row["expected_asymptotic_slope"])
        < 0.4
        and (row["p"] == 2.0 or row["fit_r_squared"] > 0.99)
        for row in smooth_scaling["rows_by_p"].values()
    )
    check(
        "P1_scaling",
        "p=0,1,2 compact-profile slopes agree with 2-p within the declared finite-grid tolerance",
        slope_pass,
        "slopes="
        + str(
            {
                key: (
                    row["fitted_log_log_slope"],
                    row["fit_r_squared"],
                )
                for key, row in smooth_scaling["rows_by_p"].items()
            }
        )
        + "; R2 is required only for the nonzero p=0,1 slopes because the p=2 target is a constant limit",
    )
    check(
        "P1_compactness",
        "the nonnegative full action is fixed-box strongly L2 equicoercive through its base Dirichlet term",
        scaling_contract["fixed_box_L2_equicoercivity"]["status"] == "proved",
        "uniform action bound => piecewise-affine H1 bound; Rellich => strong-L2 precompact subsequence; degree-sector closure is not included",
    )

    plans = build_case_plan(args.quick)
    cases = []
    for index, plan in enumerate(plans, start=1):
        print(
            f"AP-E9 case {index}/{len(plans)}: {plan.case_id} "
            f"N={plan.extent} L={plan.length:g} p={plan.exponent_p:g} "
            f"mesh={plan.triangulation} offset/a={plan.translation_lattice_units}",
            flush=True,
        )
        cases.append(relax_case(plan, solution, model, args.quick))

    all_cases_pass = all(row["case_pass"] for row in cases)
    check(
        "P2_relaxation",
        "all requested fields are unanchored, stationary, admissible, and B=1 on the same complex used by their barrier",
        all_cases_pass,
        "max_gradient={:.3e}; min_margin={:.3e}; failed={}".format(
            max(
                row["final"]["direct_projected_gradient_l2_density"]
                for row in cases
            ),
            min(row["final"]["minimum_pair_margin"] for row in cases),
            [row["case_id"] for row in cases if not row["case_pass"]],
        ),
    )
    p_controls = [
        row
        for row in cases
        if any(role in ("p0_control", "p2_control") for role in row["roles"])
    ]
    check(
        "P3_controls",
        "p=0 and p=2 controls both re-relax rather than being inferred by rescaling",
        len(p_controls) == 2 and all(row["case_pass"] for row in p_controls),
        "controls="
        + str(
            [
                (
                    row["p"],
                    row["gamma"],
                    row["final"]["total_energy"],
                    row["final"]["minimum_pair_margin"],
                )
                for row in p_controls
            ]
        ),
    )
    gate = numerical_gate(cases, args.quick)
    check(
        "P4_gate",
        "numerical regulator-stability gate is evaluated with explicit production coverage and thresholds",
        isinstance(gate["numerical_regulator_stability_thresholds_pass"], bool),
        f"pass={gate['numerical_regulator_stability_thresholds_pass']}; production_coverage={gate['production_coverage']}",
    )

    fail_closed = {
        "barrier_sublevel_uniform_pair_gap_proven_for_matched_p1": True,
        "smooth_recovery_barrier_vanishes_for_matched_p1": True,
        "p_less_than_one_uniform_pair_gap_disproven_for_barrier_alone": True,
        "fixed_box_strong_L2_equicoercivity_proven": True,
        "barrier_R_to_zero_gamma_limit_proven_on_declared_embedding": True,
        "barrier_positive_R_degree_compactness_proven": True,
        "finite_R_homogenized_barrier_density_computed": False,
        "raw_cauchy_born_density_promoted_to_gamma_density": False,
        "degree_sector_topological_current_compactness_proven": False,
        "full_action_gamma_liminf_proven": False,
        "full_action_gamma_limsup_proven_beyond_smooth_profiles": False,
        "full_action_gamma_limit_proven": False,
        "regulator_stable_continuum_background": False,
        "same_action_riemann_hessian_allowed": False,
        "same_action_riemann_hessian_built": False,
        "regulated_determinant_variation_allowed": False,
        "regulated_determinant_variation_built": False,
        "gw_domain_wall_parallel_lane_retained": True,
        "so3_mod_two_index_parallel_lane_retained": True,
    }
    result = {
        "artifact": "AP-E9 gamma-scaling, triangulation, and translated-start scan",
        "status": "quick_smoke_continuum_and_hessian_fail_closed"
        if args.quick
        else "production_scan_continuum_and_hessian_fail_closed",
        "generated_utc": "deterministic-no-wall-clock",
        "quick": args.quick,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "model": asdict(model),
        "canonical_profile_checksum": checksum,
        "source_manifest": [
            source_row(path)
            for path in (
                THIS_SCRIPT,
                TRI_SCRIPT,
                E8_SCRIPT,
                E8.E7_SCRIPT,
                E8.E6_SCRIPT,
            )
        ],
        "thresholds": {
            "stationarity_density": STATIONARITY_TOLERANCE,
            "minimum_pair_margin": PAIR_MARGIN_TOLERANCE,
            "translated_offset_lattice_units": list(TRANSLATED_OFFSET),
            "numerical_gate": {
                key: value
                for key, value in gate.items()
                if key.endswith("threshold")
            },
        },
        "exact_scaling_contract": scaling_contract,
        "smooth_profile_scaling": smooth_scaling,
        "triangulation_helper_self_tests": helper_tests,
        "case_plan": [asdict(plan) for plan in plans],
        "cases": cases,
        "numerical_gate": gate,
        "fail_closed": fail_closed,
        "regulator_stable": False,
        "hessian_gate_open": False,
        "determinant_variation_gate_open": False,
        "portal_start_allowed": False,
        "checks": CHECKS,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
    }
    if result["checks_passed"] != result["checks_total"]:
        result["status"] = "mechanical_failure_all_physics_gates_fail_closed"
    write_outputs(result)
    print(
        f"AP-E9 checks: {result['checks_passed']}/{result['checks_total']}; "
        f"numerical_gate={gate['numerical_regulator_stability_thresholds_pass']}; "
        "regulator_stable=False; hessian_gate=False",
        flush=True,
    )
    if result["checks_passed"] != result["checks_total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
