#!/usr/bin/env python3
"""Parent-resolved two-site thresholds for the hierarchical P54PQ vacuum.

This verifier does not assign broken-phase eigenvectors to a parent by their
largest component.  It constructs the exact Pati--Salam Casimir projectors on
all 328 canonical real scalar coordinates and evaluates the basis-independent
spectral traces Tr(P_parent t_i^2 log M).  The same construction is applied to
the 45 gauge generators, including the -21 vector coefficient.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("JAX_ENABLE_X64", "True")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

import jax
import numpy as np
from scipy.optimize import least_squares


REPO = Path(__file__).resolve().parents[2]
P1_SCRIPT = REPO / "route_f/code/verify_p54_p1_hessian_spectrum.py"
P2_SCRIPT = REPO / "route_f/code/verify_p54_p2_running_thresholds.py"
HIERARCHICAL_JSON = REPO / "route_f/output/p54_hierarchical_p1_search.json"
SEARCH_SCRIPT = REPO / "route_f/code/search_p54_hierarchical_p1.py"
OUTPUT_JSON = REPO / "route_f/output/p54_p2_two_site_matching.json"
OUTPUT_MD = REPO / "route_f/output/p54_p2_two_site_matching.md"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def so_generator(a: int, b: int, n: int = 10) -> np.ndarray:
    result = np.zeros((n, n))
    result[a, b] = 1.0 / math.sqrt(2.0)
    result[b, a] = -1.0 / math.sqrt(2.0)
    return result


def ps_generators(p1) -> dict[str, list[np.ndarray]]:
    su4 = [so_generator(a, b) for a in range(6) for b in range(a + 1, 6)]
    su2l = p1.sm_generators()["SU2"]
    j = {(a, b): so_generator(a, b) for a in range(6, 10) for b in range(a + 1, 10)}
    su2r = [
        (j[6, 7] + j[8, 9]) / math.sqrt(2.0),
        (j[6, 8] - j[7, 9]) / math.sqrt(2.0),
        (j[6, 9] + j[7, 8]) / math.sqrt(2.0),
    ]
    assert max(np.linalg.norm(a @ b - b @ a) for a in su2l for b in su2r) < 1e-13
    return {"SU4": su4, "SU2L": su2l, "SU2R": su2r}


def casimir_operators(p1, generators: dict[str, list[np.ndarray]]) -> dict[str, np.ndarray]:
    result = {}
    for group, group_generators in generators.items():
        casimir = np.zeros((p1.N_REAL, p1.N_REAL))
        for generator in group_generators:
            representation = p1.representation_matrix(generator)
            casimir -= representation @ representation
        result[group] = (casimir + casimir.T) / 2.0
    return result


def sm_index_operators(p1) -> list[np.ndarray]:
    reps = p1.sm_representation_matrices()
    c3 = -sum(matrix @ matrix for matrix in reps["SU3"])
    c2 = -sum(matrix @ matrix for matrix in reps["SU2"])
    y2 = -reps["Y"][0] @ reps["Y"][0]
    return [(c3 + c3.T) / 16.0, (c2 + c2.T) / 6.0, 3.0 * (y2 + y2.T) / 10.0]


def ps_index_operators(casimir: dict[str, np.ndarray]) -> list[np.ndarray]:
    return [casimir["SU4"] / 15.0, casimir["SU2L"] / 3.0, casimir["SU2R"] / 3.0]


def su4_label(casimir: float) -> str:
    labels = {0.0: "1", 2.5: "6", 4.0: "15", 4.5: "10-pair", 6.0: "20prime"}
    value, label = min(labels.items(), key=lambda item: abs(item[0] - casimir))
    return label if abs(value - casimir) < 2e-6 else f"C4={casimir:.7g}"


def su2_dimension(casimir: float) -> int:
    j = max(0.0, (-1.0 + math.sqrt(max(1.0, 1.0 + 4.0 * casimir))) / 2.0)
    return int(round(2.0 * j + 1.0))


def parent_projectors(p1, casimir: dict[str, np.ndarray]) -> list[dict[str, Any]]:
    blocks = [
        ("Phi54", p1.SL_PHI),
        ("Sigma126", slice(p1.SL_SIGMA_RE.start, p1.SL_SIGMA_IM.stop)),
        ("phi10", slice(p1.SL_H_RE.start, p1.SL_H_IM.stop)),
        ("S", p1.SL_S),
    ]
    combined = casimir["SU4"] + math.sqrt(2.0) * casimir["SU2L"] + math.pi * casimir["SU2R"]
    parents = []
    for field, block_slice in blocks:
        indices = np.arange(block_slice.start, block_slice.stop)
        values, vectors = np.linalg.eigh(combined[np.ix_(indices, indices)])
        start = 0
        while start < len(values):
            stop = start + 1
            while stop < len(values) and abs(values[stop] - values[start]) < 2e-7:
                stop += 1
            local = vectors[:, start:stop]
            basis = np.zeros((p1.N_REAL, stop - start))
            basis[indices, :] = local
            dim = stop - start
            c4 = float(np.trace(basis.T @ casimir["SU4"] @ basis) / dim)
            cl = float(np.trace(basis.T @ casimir["SU2L"] @ basis) / dim)
            cr = float(np.trace(basis.T @ casimir["SU2R"] @ basis) / dim)
            projector = basis @ basis.T
            label = f"{field}:({su4_label(c4)},{su2_dimension(cl)},{su2_dimension(cr)})"
            parents.append(
                {
                    "label": label,
                    "field": field,
                    "dimension": dim,
                    "C4": c4,
                    "C2L": cl,
                    "C2R": cr,
                    "projector": projector,
                }
            )
            start = stop
    return parents


def spectral_log(eigenvalues: np.ndarray, eigenvectors: np.ndarray, scale: float, tolerance: float) -> np.ndarray:
    mask = eigenvalues > tolerance
    basis = eigenvectors[:, mask]
    logs = np.log(np.sqrt(eigenvalues[mask]) / scale)
    return (basis * logs) @ basis.T


def trace_vector(index_operators: list[np.ndarray], projector: np.ndarray, log_operator: np.ndarray) -> np.ndarray:
    return np.array([float(np.trace(index @ projector @ log_operator)) for index in index_operators])


def adjoint_representation(basis: list[np.ndarray], generator: np.ndarray) -> np.ndarray:
    result = np.zeros((len(basis), len(basis)))
    for column, element in enumerate(basis):
        variation = generator @ element - element @ generator
        result[:, column] = [np.trace(out.T @ variation) for out in basis]
    return result


def adjoint_index_operators(groups: dict[str, list[np.ndarray]], basis: list[np.ndarray], dimensions: dict[str, int], hypercharge: bool = False) -> list[np.ndarray]:
    operators = []
    for group, group_generators in groups.items():
        c = np.zeros((len(basis), len(basis)))
        for generator in group_generators:
            representation = adjoint_representation(basis, generator)
            c -= representation @ representation
        divisor = dimensions[group]
        if hypercharge and group == "Y":
            divisor = 5.0 / 3.0
        operators.append((c + c.T) / (2.0 * divisor))
    return operators


def vector_thresholds(p1, x0: np.ndarray, ratio: float, g_u: float) -> dict[str, Any]:
    basis = [so_generator(a, b) for a in range(10) for b in range(a + 1, 10)]
    mass2 = p1.gauge_orbit(x0).T @ p1.gauge_orbit(x0)
    values, vectors = np.linalg.eigh((mass2 + mass2.T) / 2.0)
    sm = p1.sm_generators()
    ps = ps_generators(p1)
    sm_ops = adjoint_index_operators(
        {"SU3": sm["SU3"], "SU2": sm["SU2"], "Y": sm["Y"]},
        basis,
        {"SU3": 8, "SU2": 3, "Y": 1},
        hypercharge=True,
    )
    ps_ops = adjoint_index_operators(ps, basis, {"SU4": 15, "SU2L": 3, "SU2R": 3})
    # The exact spectrum is {r^2 x 8, 5 r^2 x 1, 1 x 12,
    # (1+r^2) x 12}; the gap remains below 0.8 throughout the declared
    # hierarchy search r<0.4.
    low = (values > 2e-10) & (values < 0.8)
    high = values >= 0.8

    def logop(mask: np.ndarray, scale: float) -> np.ndarray:
        selected = vectors[:, mask]
        logs = np.log(g_u * np.sqrt(values[mask]) / scale)
        return (selected * logs) @ selected.T

    # Matching surfaces are tied to the representative gauge-boson masses,
    # mu_I=g*sigma and mu_U=g*omega.  This removes an otherwise spurious
    # universal log(g) from the vector threshold ledger.
    lambda_i = -21.0 * np.array([np.trace(op @ logop(low, g_u * ratio)) for op in sm_ops])
    lambda_u = -21.0 * np.array([np.trace(op @ logop(high, g_u)) for op in ps_ops])
    nuisance_rows = []
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and abs(values[stop] - values[start]) < 2e-9:
            stop += 1
        if values[start] > 2e-10:
            block = vectors[:, start:stop]
            projector = block @ block.T
            if values[start] < 0.8:
                coefficient_i = -21.0 * np.array([np.trace(op @ projector) for op in sm_ops])
                coefficient_u = np.zeros(3)
                site = "MI"
            else:
                coefficient_i = np.zeros(3)
                coefficient_u = -21.0 * np.array([np.trace(op @ projector) for op in ps_ops])
                site = "MU"
            nuisance_rows.append(
                {
                    "name": f"vector_m2_{float(np.mean(values[start:stop])):.9g}",
                    "site": site,
                    "multiplicity": stop - start,
                    "d_lambda_I_d_log_mass": coefficient_i.tolist(),
                    "d_lambda_U_d_log_mass": coefficient_u.tolist(),
                }
            )
        start = stop
    return {
        "mass2_groups_over_g2omega2": p1.group_eigenvalues(values, abs_tol=2e-10),
        "intermediate_count": int(np.sum(low)),
        "gut_count": int(np.sum(high)),
        "lambda_vector_MI_3_2_1": lambda_i.tolist(),
        "lambda_vector_MU_4_L_R": lambda_u.tolist(),
        "log_mass_nuisance_rows": nuisance_rows,
    }


def solve_thresholded(
    p2,
    lambda_i: np.ndarray,
    lambda_u: np.ndarray,
    initial: np.ndarray,
    inputs: dict[str, float] | None = None,
) -> dict[str, Any]:
    actual_inputs = p2.CURRENT_INPUT if inputs is None else inputs
    alpha0 = p2.sm_inverse_couplings(actual_inputs)

    def residual(z: np.ndarray) -> np.ndarray:
        if np.min(z) <= 0:
            return np.array([1e3, 1e3])
        at_i = p2.run_two_loop(alpha0, z[0], p2.A_SM, p2.B_SM)
        ps_i = p2.sm_to_ps(at_i, additive_shift=lambda_i / (12.0 * math.pi))
        ps_u = p2.run_two_loop(ps_i, z[1], p2.A_PS, p2.B_PS)
        implied = p2.so10_implied(ps_u) + lambda_u / (12.0 * math.pi)
        return np.array([implied[0] - implied[1], implied[0] - implied[2]])

    solved = least_squares(residual, initial, xtol=2e-13, ftol=2e-13, gtol=2e-13, max_nfev=160)
    if not solved.success or np.linalg.norm(solved.fun) > 2e-7:
        raise RuntimeError(f"threshold solve failed: {solved.fun}")
    l_sm, l_ps = solved.x
    mi = p2.MZ * math.exp(l_sm)
    mu = mi * math.exp(l_ps)
    at_i = p2.run_two_loop(alpha0, l_sm, p2.A_SM, p2.B_SM)
    ps_i = p2.sm_to_ps(at_i, additive_shift=lambda_i / (12.0 * math.pi))
    ps_u = p2.run_two_loop(ps_i, l_ps, p2.A_PS, p2.B_PS)
    implied = p2.so10_implied(ps_u) + lambda_u / (12.0 * math.pi)
    return {
        "MI_GeV": mi,
        "MU_GeV": mu,
        "MI_over_MU": mi / mu,
        "alpha_U_inverse": float(np.mean(implied)),
        "matching_residual": solved.fun.tolist(),
        "solver_coordinates": solved.x.tolist(),
    }


def scalar_nuisance_rows(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    tolerance: float,
    p_intermediate: np.ndarray,
    p_gut: np.ndarray,
    sm_ops: list[np.ndarray],
    ps_ops: list[np.ndarray],
) -> list[dict[str, Any]]:
    rows = []
    start = 0
    while start < len(eigenvalues):
        stop = start + 1
        while stop < len(eigenvalues) and abs(eigenvalues[stop] - eigenvalues[start]) <= max(
            tolerance, 2e-7 * max(1.0, abs(eigenvalues[start]))
        ):
            stop += 1
        if eigenvalues[start] > tolerance:
            block = eigenvectors[:, start:stop]
            projector = block @ block.T
            coefficient_i = np.array(
                [np.trace(op @ p_intermediate @ projector) for op in sm_ops]
            )
            coefficient_u = np.array(
                [np.trace(op @ p_gut @ projector) for op in ps_ops]
            )
            rows.append(
                {
                    "name": f"scalar_m2_{float(np.mean(eigenvalues[start:stop])):.9g}",
                    "multiplicity": stop - start,
                    "d_lambda_I_d_log_mass": coefficient_i.tolist(),
                    "d_lambda_U_d_log_mass": coefficient_u.tolist(),
                }
            )
        start = stop
    return rows


def covariance_report(
    p2,
    base_solution: dict[str, Any],
    lambda_i: np.ndarray,
    lambda_u: np.ndarray,
    nuisance_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    def coordinates(solution: dict[str, Any]) -> np.ndarray:
        return np.array(
            [
                math.log10(solution["MI_GeV"]),
                math.log10(solution["MU_GeV"]),
                solution["alpha_U_inverse"],
            ]
        )

    initial = np.array(base_solution["solver_coordinates"])
    input_keys = ["alpha_em_inverse", "sin2_theta_w", "alpha_s"]
    input_steps = np.array([2e-4, 5e-7, 2e-5])
    input_point = np.array([p2.CURRENT_INPUT[key] for key in input_keys])
    input_columns = []
    for index, step in enumerate(input_steps):
        plus = input_point.copy(); plus[index] += step
        minus = input_point.copy(); minus[index] -= step
        plus_inputs = {key: float(value) for key, value in zip(input_keys, plus)}
        minus_inputs = {key: float(value) for key, value in zip(input_keys, minus)}
        y_plus = coordinates(solve_thresholded(p2, lambda_i, lambda_u, initial, plus_inputs))
        y_minus = coordinates(solve_thresholded(p2, lambda_i, lambda_u, initial, minus_inputs))
        input_columns.append((y_plus - y_minus) / (2.0 * step))
    input_jacobian = np.column_stack(input_columns)
    input_sigma = np.array([p2.CURRENT_SIGMA[key] for key in input_keys])
    covariance_exp = input_jacobian @ np.diag(input_sigma**2) @ input_jacobian.T

    lambda_point = np.concatenate([lambda_i, lambda_u])
    lambda_columns = []
    step = 2e-5
    for index in range(6):
        plus = lambda_point.copy(); plus[index] += step
        minus = lambda_point.copy(); minus[index] -= step
        y_plus = coordinates(solve_thresholded(p2, plus[:3], plus[3:], initial))
        y_minus = coordinates(solve_thresholded(p2, minus[:3], minus[3:], initial))
        lambda_columns.append((y_plus - y_minus) / (2.0 * step))
    lambda_jacobian = np.column_stack(lambda_columns)

    nuisance_matrix = np.array(
        [row["d_lambda_I_d_log_mass"] + row["d_lambda_U_d_log_mass"] for row in nuisance_rows],
        dtype=float,
    ).T
    output_nuisance_jacobian = lambda_jacobian @ nuisance_matrix
    sigma_log_mass = math.log(1.10)
    covariance_threshold = (
        output_nuisance_jacobian
        @ (np.eye(len(nuisance_rows)) * sigma_log_mass**2)
        @ output_nuisance_jacobian.T
    )
    covariance_total = covariance_exp + covariance_threshold
    return {
        "output_coordinates": ["log10_MI_GeV", "log10_MU_GeV", "alpha_U_inverse"],
        "mass_nuisance_convention": "independent 10 percent log-mass error per exactly degenerate scalar/vector block",
        "mass_nuisance_count": len(nuisance_rows),
        "sigma_log_mass": sigma_log_mass,
        "experimental_input_jacobian": input_jacobian.tolist(),
        "lambda_I_U_to_output_jacobian": lambda_jacobian.tolist(),
        "log_mass_to_lambda_I_U_matrix": nuisance_matrix.tolist(),
        "experimental_covariance": covariance_exp.tolist(),
        "threshold_covariance": covariance_threshold.tolist(),
        "total_covariance": covariance_total.tolist(),
        "experimental_sigma": np.sqrt(np.diag(covariance_exp)).tolist(),
        "threshold_sigma": np.sqrt(np.diag(covariance_threshold)).tolist(),
        "total_sigma": np.sqrt(np.diag(covariance_total)).tolist(),
        "threshold_covariance_eigenvalues": np.linalg.eigvalsh(covariance_threshold).tolist(),
        "nuisance_rows": nuisance_rows,
    }


def run(target_ratio: float | None = None) -> dict[str, Any]:
    started = time.time()
    p1 = load_module("p54_p1", P1_SCRIPT)
    p2 = load_module("p54_p2", P2_SCRIPT)
    search = load_module("p54_hierarchy_search", SEARCH_SCRIPT)
    hierarchy = json.loads(HIERARCHICAL_JSON.read_text(encoding="utf-8"))
    ratio = float(hierarchy["vevs"]["sigma"] if target_ratio is None else target_ratio)
    vs = float(hierarchy["vevs"]["vs"])
    p = search.parameters_at(p1, ratio, vs)
    x0 = p1.vacuum_vector(1.0, ratio, vs)
    potential = p1.potential_factory(p)
    hessian0 = np.asarray(jax.jit(p1.hessian(potential))(p1.anp.asarray(x0)), dtype=float)
    hessian0 = (hessian0 + hessian0.T) / 2.0
    symmetry, _, complement = search.symmetry_complement(p1, x0)
    projector10 = np.zeros((p1.N_REAL, p1.N_REAL))
    projector10[p1.SL_H_RE, p1.SL_H_RE] = np.eye(p1.SL_H_RE.stop - p1.SL_H_RE.start)
    projector10[p1.SL_H_IM, p1.SL_H_IM] = np.eye(p1.SL_H_IM.stop - p1.SL_H_IM.start)
    tuning = search.tune_first_instability(hessian0, projector10, complement)
    if not tuning.get("tunable", False):
        raise RuntimeError(f"single-doublet tuning failed at ratio {ratio}: {tuning}")
    p["xi02"] = float(tuning["xi02"])
    hessian = hessian0 - float(tuning["xi02"]) * projector10
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    tolerance = float(tuning["zero_tolerance"])
    extra_zero = p1.classify_extra_zero_modes(
        eigenvalues, eigenvectors, tolerance, symmetry
    )
    scalar_irreps = p1.classify_full_spectrum(
        eigenvalues, eigenvectors, tolerance, p1.sm_representation_matrices()
    )
    positive_doublets = sorted(
        float(row["m2"])
        for row in scalar_irreps["rows"]
        if row["SU3"] == "1"
        and int(row["SU2_dimension"]) == 2
        and abs(float(row["abs_hypercharge"]) - 0.5) < 1e-5
        and float(row["m2"]) > tolerance
    )
    weight10 = float(extra_zero["field_weights"]["phi10"])
    retuning = {
        "tree_xi02": float(tuning["xi02"]),
        "light_field_weight_10": weight10,
        "light_field_weight_126": float(extra_zero["field_weights"]["Sigma126"]),
        "renormalized_condition": "m_h,pole^2=hT[H_D+Pi_D(0;mu)-delta_xi02*P10]h=0",
        "retuning_formula": "delta_xi02=hT Pi_D(0;mu) h / light_field_weight_10",
        "delta_xi02_per_unit_projected_Pi_over_omega2": 1.0 / weight10,
        "nearest_heavy_doublet_gap_m2_over_omega2": positive_doublets[0],
        "single_doublet_sufficient_condition": "norm_Q(Pi_D-delta_xi02*P10)Q < nearest_heavy_doublet_gap",
        "condition_reimposed_as_renormalization_condition": True,
        "finite_self_energy_numerically_predicted": False,
        "prediction_boundary": "requires field-dependent CW Hessian, gauge/ghost prescription, and renormalized heavy Yukawa matrices",
    }

    ps = ps_generators(p1)
    casimir = casimir_operators(p1, ps)
    parents = parent_projectors(p1, casimir)
    parent_census = [
        {key: row[key] for key in ("label", "field", "dimension", "C4", "C2L", "C2R")}
        for row in parents
    ]
    sigma_vector = np.zeros(p1.N_REAL)
    sigma_vector[p1.SL_SIGMA_RE] = x0[p1.SL_SIGMA_RE]
    sigma_vector[p1.SL_SIGMA_IM] = x0[p1.SL_SIGMA_IM]
    sigma_norm2 = float(sigma_vector @ sigma_vector)
    sigma_weights = {
        row["label"]: float(sigma_vector @ row["projector"] @ sigma_vector / sigma_norm2)
        for row in parents
        if row["field"] == "Sigma126"
    }
    breaking_parent = max(sigma_weights, key=sigma_weights.get)
    p_intermediate = next(row["projector"] for row in parents if row["label"] == breaking_parent)
    p_gut = np.eye(p1.N_REAL) - p_intermediate

    baseline = p2.solve_two_stage(p2.CURRENT_INPUT)
    sm_ops = sm_index_operators(p1)
    ps_ops = ps_index_operators(casimir)
    g_u = math.sqrt(4.0 * math.pi / baseline["alpha_U_inverse"])
    gauge_iterations = []
    thresholded = baseline
    for iteration in range(8):
        log_i = spectral_log(eigenvalues, eigenvectors, g_u * ratio, tolerance)
        log_u = spectral_log(eigenvalues, eigenvectors, g_u, tolerance)
        scalar_i = trace_vector(sm_ops, p_intermediate, log_i)
        scalar_u = trace_vector(ps_ops, p_gut, log_u)
        vector = vector_thresholds(p1, x0, ratio, g_u)
        lambda_i = scalar_i + np.array(vector["lambda_vector_MI_3_2_1"])
        lambda_u = scalar_u + np.array(vector["lambda_vector_MU_4_L_R"])
        thresholded = solve_thresholded(
            p2, lambda_i, lambda_u, np.array(thresholded["solver_coordinates"])
        )
        updated_g = math.sqrt(4.0 * math.pi / thresholded["alpha_U_inverse"])
        gauge_iterations.append(
            {"iteration": iteration, "g_input": g_u, "g_output": updated_g}
        )
        if abs(updated_g / g_u - 1.0) < 2e-10:
            g_u = updated_g
            break
        g_u = updated_g

    scalar_nuisances = scalar_nuisance_rows(
        eigenvalues, eigenvectors, tolerance, p_intermediate, p_gut, sm_ops, ps_ops
    )
    nuisance_rows = scalar_nuisances + vector["log_mass_nuisance_rows"]
    covariance = covariance_report(p2, thresholded, lambda_i, lambda_u, nuisance_rows)
    consistency_factor = thresholded["MI_over_MU"] / ratio

    parent_dim = sum(row["dimension"] for row in parent_census)
    ps_decomposition_ok = (
        parent_dim == p1.N_REAL
        and all(not row["label"].split(":", 1)[1].startswith("(C4=") for row in parent_census)
    )
    checks = [
        {"name": "PS parent projectors cover 328 real coordinates", "pass": parent_dim == 328},
        {"name": "fixed-point scalar ledger has 35 SM sectors", "pass": len(scalar_irreps["rows"]) == 35 and scalar_irreps["real_dimension_sum"] == 328},
        {"name": "fixed-point Hessian has 38 zeros and no tachyon", "pass": int(np.sum(np.abs(eigenvalues) < tolerance)) == 38 and int(np.sum(eigenvalues < -tolerance)) == 0},
        {"name": "four nonsymmetry zeros are one doublet", "pass": extra_zero["rank"] == 4},
        {"name": "all PS Casimirs have recognized labels", "pass": ps_decomposition_ok},
        {"name": "126 VEV lies in a unique PS parent", "pass": sigma_weights[breaking_parent] > 1.0 - 2e-10},
        {"name": "vector site census is 9 plus 24", "pass": vector["intermediate_count"] == 9 and vector["gut_count"] == 24},
        {"name": "two-site unification residual closes", "pass": max(abs(x) for x in thresholded["matching_residual"]) < 2e-7},
        {"name": "one replay keeps the hierarchy within a factor two", "pass": 0.5 < consistency_factor < 2.0},
        {"name": "gU threshold iteration converges", "pass": abs(gauge_iterations[-1]["g_output"] / gauge_iterations[-1]["g_input"] - 1.0) < 2e-8},
        {"name": "two-site threshold covariance is positive semidefinite", "pass": min(covariance["threshold_covariance_eigenvalues"]) > -1e-11},
    ]
    return {
        "schema": "route-f-p54-p2-two-site-spectral-matching-v1",
        "date": "2026-08-30",
        "status": {
            "parent_resolved_spectral_trace_built": True,
            "two_site_scalar_and_vector_matching_built": True,
            "single_replay_hierarchy_consistent_within_factor_two": 0.5 < consistency_factor < 2.0,
            "fixed_point_iteration_closed": abs(math.log(consistency_factor)) < 2e-3,
            "interpretation": "extended-survival site projector; no nearest-mass assignment",
        },
        "hierarchical_vacuum_ratio": ratio,
        "doublet_tuning": tuning,
        "loop_corrected_doublet_condition": retuning,
        "fixed_point_scalar_spectrum": {
            "eigenvalue_groups": p1.group_eigenvalues(eigenvalues, abs_tol=tolerance),
            "sm_irrep_ledger": scalar_irreps,
        },
        "parent_census": parent_census,
        "sigma_parent_weights": sigma_weights,
        "intermediate_breaking_parent": breaking_parent,
        "thresholds": {
            "matching_scale_convention": "mu_I=gU*sigma, mu_U=gU*omega",
            "lambda_scalar_MI_3_2_1": scalar_i.tolist(),
            "lambda_scalar_MU_4_L_R": scalar_u.tolist(),
            "lambda_total_MI_3_2_1": lambda_i.tolist(),
            "lambda_total_MU_4_L_R": lambda_u.tolist(),
        },
        "vector": vector,
        "gauge_coupling_iteration": gauge_iterations,
        "no_threshold_baseline": baseline,
        "thresholded_solution": thresholded,
        "hierarchy_replay_factor": consistency_factor,
        "covariance": covariance,
        "checks": checks,
        "summary": {"passed": sum(row["pass"] for row in checks), "total": len(checks), "all_pass": all(row["pass"] for row in checks)},
        "runtime_seconds": time.time() - started,
        "sources": [
            {"path": str(P1_SCRIPT.relative_to(REPO)), "sha256": sha256(P1_SCRIPT)},
            {"path": str(P2_SCRIPT.relative_to(REPO)), "sha256": sha256(P2_SCRIPT)},
            {"path": str(HIERARCHICAL_JSON.relative_to(REPO)), "sha256": sha256(HIERARCHICAL_JSON)},
            {"path": str(SEARCH_SCRIPT.relative_to(REPO)), "sha256": sha256(SEARCH_SCRIPT)},
            {"path": str(Path(__file__).resolve().relative_to(REPO)), "sha256": sha256(Path(__file__).resolve())},
        ],
    }


def markdown(report: dict[str, Any]) -> str:
    thresholded = report["thresholded_solution"]
    lines = [
        "# Hierarchical P54PQ two-site matching",
        "",
        f"Status: **{report['summary']['passed']}/{report['summary']['total']} checks passed**.",
        "",
        "The calculation uses exact PS Casimir projectors and spectral traces, not a nearest-mass assignment of broken-phase states.",
        "",
        f"- PS-breaking parent: `{report['intermediate_breaking_parent']}`;",
        f"- scalar `lambda_I(3,2,1) = {report['thresholds']['lambda_scalar_MI_3_2_1']}`;",
        f"- scalar `lambda_U(4,L,R) = {report['thresholds']['lambda_scalar_MU_4_L_R']}`;",
        f"- total `lambda_I = {report['thresholds']['lambda_total_MI_3_2_1']}`;",
        f"- total `lambda_U = {report['thresholds']['lambda_total_MU_4_L_R']}`;",
        f"- thresholded `MI={thresholded['MI_GeV']:.6e} GeV`, `MU={thresholded['MU_GeV']:.6e} GeV`, `alphaU^-1={thresholded['alpha_U_inverse']:.7f}`;",
        f"- resulting `MI/MU={thresholded['MI_over_MU']:.9g}` versus the input vacuum ratio `{report['hierarchical_vacuum_ratio']:.9g}` (replay factor `{report['hierarchy_replay_factor']:.6g}`).",
        f"- covariance sigma in `(log10 MI, log10 MU, alphaU^-1)`: experiment `{report['covariance']['experimental_sigma']}`, thresholds `{report['covariance']['threshold_sigma']}`, total `{report['covariance']['total_sigma']}`.",
        "",
        "## PS parent census",
        "",
        "| parent | real dimension | C4 | C2L | C2R |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in report["parent_census"]:
        lines.append(f"| {row['label']} | {row['dimension']} | {row['C4']:.6g} | {row['C2L']:.6g} | {row['C2R']:.6g} |")
    lines.extend(["", "## Checks", "", "| check | result |", "|---|---|"])
    for row in report["checks"]:
        lines.append(f"| {row['name']} | {'PASS' if row['pass'] else 'FAIL'} |")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    target = os.environ.get("ROUTE_F_TARGET_RATIO")
    report = run(None if target is None else float(target))
    OUTPUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUTPUT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    print(json.dumps(report["status"], sort_keys=True))
    print(json.dumps(report["thresholded_solution"], sort_keys=True))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
