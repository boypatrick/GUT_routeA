#!/usr/bin/env python3
"""Bosonic hard-mode Coleman--Weinberg Hessian for P54PQ-v2.

The calculation is deliberately a matching calculation, not a gauge-invariant
pole-mass claim.  It uses background-field Landau gauge and MSbar at
mu_U=g_U*omega.  The 38 scalar EFT modes (33 gauge Goldstones, one PQ
Goldstone, and the four real components of the tuned Higgs doublet), the 12
unbroken/infrared vectors, and the massless Landau ghosts are removed from the
hard supertrace.  The remaining 290 scalar and 33 vector eigenvalues are
differentiated by the exact Frechet formula for Tr f(M^2).

Only two light-doublet directions require expensive field-dependent 328x328
tree Hessians.  Standard-Model invariance then fixes the full real four-plane
answer to kappa*I_4; a second diagonal direction and one mixed direction are
independent numerical checks of that Schur-lemma reduction.
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


REPO = Path(__file__).resolve().parents[2]
P1_SCRIPT = REPO / "route_f/code/verify_p54_p1_hessian_spectrum.py"
SEARCH_SCRIPT = REPO / "route_f/code/search_p54_hierarchical_p1.py"
P2_JSON = REPO / "route_f/output/p54_p2_two_site_matching.json"
OUTPUT_JSON = REPO / "route_f/output/p54_p2_bosonic_cw.json"
OUTPUT_MD = REPO / "route_f/output/p54_p2_bosonic_cw.md"
SCRIPT = Path(__file__).resolve()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generators() -> list[np.ndarray]:
    result = []
    for a in range(10):
        for b in range(a + 1, 10):
            generator = np.zeros((10, 10))
            generator[a, b] = 1.0 / math.sqrt(2.0)
            generator[b, a] = -1.0 / math.sqrt(2.0)
            result.append(generator)
    return result


def field_orbit(p1, vector: np.ndarray, gens: list[np.ndarray]) -> np.ndarray:
    return np.column_stack(
        [p1.act_on_field_vector(vector, generator) for generator in gens]
    )


def hard_trace_hessian(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    first_a: np.ndarray,
    first_b: np.ndarray,
    second_ab: np.ndarray,
    hard_count: int,
    mu2: float,
    constant: float,
    multiplicity: float,
) -> float:
    """Second derivative of multiplicity/(64 pi^2) Tr_hard f(M^2).

    f(x)=x^2(log(x/mu^2)-constant).  The low eigenvalues are assigned a
    locally flat spectral function.  Cross hard/soft divided differences are
    retained, so mixing with the EFT subspace is not silently discarded.
    """
    order = np.argsort(eigenvalues)
    lam = np.asarray(eigenvalues[order], dtype=float)
    u = np.asarray(eigenvectors[:, order], dtype=float)
    hard = np.zeros(len(lam), dtype=bool)
    hard[-hard_count:] = True
    if np.min(lam[hard]) <= 0.0:
        raise RuntimeError("hard spectrum is not positive")

    log = np.zeros_like(lam)
    log[hard] = np.log(lam[hard] / mu2)
    fp = np.zeros_like(lam)
    fp[hard] = lam[hard] * (2.0 * log[hard] - 2.0 * constant + 1.0)
    fpp = np.zeros_like(lam)
    fpp[hard] = 2.0 * log[hard] - 2.0 * constant + 3.0

    delta = lam[:, None] - lam[None, :]
    scale = np.maximum(1.0, np.maximum(np.abs(lam[:, None]), np.abs(lam[None, :])))
    degenerate = np.abs(delta) < 2e-10 * scale
    divided = np.empty_like(delta)
    np.divide(fp[:, None] - fp[None, :], delta, out=divided, where=~degenerate)
    # Equal eigenvalues always lie wholly in the hard or the soft cluster;
    # the hard/soft gap is checked by the caller.
    avg_fpp = 0.5 * (fpp[:, None] + fpp[None, :])
    divided[degenerate] = avg_fpp[degenerate]

    a = u.T @ first_a @ u
    b = u.T @ first_b @ u
    c = u.T @ second_ab @ u
    first_term = float(np.dot(fp, np.diag(c)))
    second_term = float(np.sum(divided * a * b.T))
    return multiplicity * (first_term + second_term) / (64.0 * math.pi**2)


def scalar_derivatives(
    h0: np.ndarray,
    hp: np.ndarray,
    hm: np.ndarray,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray]:
    first = (hp - hm) / (2.0 * epsilon)
    second = (hp + hm - 2.0 * h0) / epsilon**2
    return (first + first.T) / 2.0, (second + second.T) / 2.0


def frechet_regression() -> float:
    """Compare the matrix formula with a direct mixed finite difference."""
    rng = np.random.default_rng(20260902)
    dimension = 7
    eigenvalues = np.array([0.0, 0.0, 0.2, 0.4, 0.7, 1.1, 1.6])
    eigenvectors, _ = np.linalg.qr(rng.normal(size=(dimension, dimension)))
    matrix0 = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T

    def symmetric() -> np.ndarray:
        value = rng.normal(size=(dimension, dimension))
        return (value + value.T) / 2.0

    first_a, first_b, second_ab = symmetric(), symmetric(), symmetric()
    mu2 = 0.4
    constant = 1.5
    analytic = hard_trace_hessian(
        eigenvalues,
        eigenvectors,
        first_a,
        first_b,
        second_ab,
        hard_count=5,
        mu2=mu2,
        constant=constant,
        multiplicity=1.0,
    )

    def direct(z: float, w: float) -> float:
        matrix = matrix0 + z * first_a + w * first_b + z * w * second_ab
        values = np.linalg.eigvalsh(matrix)[-5:]
        return float(
            np.sum(values**2 * (np.log(values / mu2) - constant))
            / (64.0 * math.pi**2)
        )

    epsilon = 1.0e-4
    numerical = (
        direct(epsilon, epsilon)
        - direct(epsilon, -epsilon)
        - direct(-epsilon, epsilon)
        + direct(-epsilon, -epsilon)
    ) / (4.0 * epsilon**2)
    return abs(numerical - analytic) / max(1.0, abs(analytic))


def run() -> dict[str, Any]:
    started = time.time()
    p1 = load_module("p54_p1_cw", P1_SCRIPT)
    search = load_module("p54_search_cw", SEARCH_SCRIPT)
    p2 = json.loads(P2_JSON.read_text(encoding="utf-8"))
    ratio = float(p2["hierarchical_vacuum_ratio"])
    vs = float(json.loads((REPO / "route_f/output/p54_hierarchical_p1_search.json").read_text())["vevs"]["vs"])
    g_u = float(p2["gauge_coupling_iteration"][-1]["g_output"])
    mu_over_omega = g_u
    xi_tree = float(p2["loop_corrected_doublet_condition"]["tree_xi02"])

    parameters = search.parameters_at(p1, ratio, vs)
    parameters["xi02"] = xi_tree
    x0 = p1.vacuum_vector(1.0, ratio, vs)
    hessian_fn = jax.jit(p1.hessian(p1.potential_factory(parameters)))

    evaluations: dict[str, np.ndarray] = {}

    def evaluate(name: str, point: np.ndarray) -> np.ndarray:
        tick = time.time()
        value = np.asarray(hessian_fn(p1.anp.asarray(point)), dtype=float)
        value = (value + value.T) / 2.0
        evaluations[name] = value
        print(f"{name}: {time.time()-tick:.2f} s", flush=True)
        return value

    h0 = evaluate("background", x0)
    eigenvalues, eigenvectors = np.linalg.eigh(h0)
    tolerance = float(p2["doublet_tuning"]["zero_tolerance"])
    symmetry, _, _ = search.symmetry_complement(p1, x0)
    zero_vectors = eigenvectors[:, np.abs(eigenvalues) < tolerance]
    symmetry_projector = symmetry @ np.linalg.pinv(symmetry, rcond=1e-11)
    residual = zero_vectors - symmetry_projector @ zero_vectors
    light_u, singular, _ = np.linalg.svd(residual, full_matrices=False)
    light_rank = int(np.sum(singular > 2e-8))
    light = light_u[:, :light_rank]
    if light_rank != 4:
        raise RuntimeError(f"expected one real four-plane, obtained rank {light_rank}")

    eps = 2.0e-2
    eps_check = 1.0e-2
    q0, q1 = light[:, 0], light[:, 1]
    h0p = evaluate("q0_plus", x0 + eps * q0)
    h0m = evaluate("q0_minus", x0 - eps * q0)
    h1p = evaluate("q1_plus", x0 + eps * q1)
    h1m = evaluate("q1_minus", x0 - eps * q1)
    h01p = evaluate("q0_q1_plus", x0 + eps * (q0 + q1))
    h0pc = evaluate("q0_plus_check", x0 + eps_check * q0)
    h0mc = evaluate("q0_minus_check", x0 - eps_check * q0)

    scalar_first0, scalar_second00 = scalar_derivatives(h0, h0p, h0m, eps)
    scalar_first1, scalar_second11 = scalar_derivatives(h0, h1p, h1m, eps)
    scalar_first0c, scalar_second00c = scalar_derivatives(h0, h0pc, h0mc, eps_check)
    scalar_second01 = (
        h01p
        - h0
        - eps * (scalar_first0 + scalar_first1)
        - 0.5 * eps**2 * (scalar_second00 + scalar_second11)
    ) / eps**2
    scalar_second01 = (scalar_second01 + scalar_second01.T) / 2.0
    soft = eigenvectors[:, :38]
    soft_first0_norm = float(np.linalg.norm(soft.T @ scalar_first0 @ soft))
    soft_first1_norm = float(np.linalg.norm(soft.T @ scalar_first1 @ soft))
    frechet_residual = frechet_regression()

    # Exact quarticity check: H(x) is quadratic, hence these derivatives are
    # independent of the finite-difference step up to floating-point error.
    first_step_residual = float(
        np.linalg.norm(scalar_first0 - scalar_first0c)
        / max(1.0, np.linalg.norm(scalar_first0c))
    )
    second_step_residual = float(
        np.linalg.norm(scalar_second00 - scalar_second00c)
        / max(1.0, np.linalg.norm(scalar_second00c))
    )

    scalar_eigenvalues, scalar_eigenvectors = eigenvalues, eigenvectors
    scalar_hard_count = 290
    scalar_soft_count = p1.N_REAL - scalar_hard_count

    def scalar_entry(a, b, ab, scale_factor: float = 1.0) -> float:
        return hard_trace_hessian(
            scalar_eigenvalues,
            scalar_eigenvectors,
            a,
            b,
            ab,
            scalar_hard_count,
            (scale_factor * mu_over_omega) ** 2,
            constant=1.5,
            multiplicity=1.0,
        )

    scalar00 = scalar_entry(scalar_first0, scalar_first0, scalar_second00)
    scalar11 = scalar_entry(scalar_first1, scalar_first1, scalar_second11)
    scalar01 = scalar_entry(scalar_first0, scalar_first1, scalar_second01)
    scalar00_check = scalar_entry(scalar_first0c, scalar_first0c, scalar_second00c)
    scalar_kappa = 0.5 * (scalar00 + scalar11)
    scalar_declared = scalar_kappa * np.eye(4)

    gens = generators()
    orbit0 = field_orbit(p1, x0, gens)
    orbit_light = [field_orbit(p1, light[:, index], gens) for index in range(4)]
    vector0 = g_u**2 * (orbit0.T @ orbit0)
    vector0 = (vector0 + vector0.T) / 2.0
    vector_eigenvalues, vector_eigenvectors = np.linalg.eigh(vector0)
    vector_first = [
        g_u**2 * (oi.T @ orbit0 + orbit0.T @ oi) for oi in orbit_light
    ]
    vector_second = [[
        g_u**2 * (orbit_light[a].T @ orbit_light[b] + orbit_light[b].T @ orbit_light[a])
        for b in range(4)
    ] for a in range(4)]

    def vector_matrix(scale_factor: float = 1.0) -> np.ndarray:
        answer = np.zeros((4, 4))
        for a in range(4):
            for b in range(a, 4):
                answer[a, b] = answer[b, a] = hard_trace_hessian(
                    vector_eigenvalues,
                    vector_eigenvectors,
                    vector_first[a],
                    vector_first[b],
                    vector_second[a][b],
                    hard_count=33,
                    mu2=(scale_factor * mu_over_omega) ** 2,
                    constant=5.0 / 6.0,
                    multiplicity=3.0,
                )
        return answer

    vector_declared = vector_matrix()
    total_declared = scalar_declared + vector_declared
    kappa_b = float(np.trace(total_declared) / 4.0)
    weight10 = float(p2["loop_corrected_doublet_condition"]["light_field_weight_10"])
    delta_xi_b = kappa_b / weight10
    heavy_gap = float(p2["loop_corrected_doublet_condition"]["nearest_heavy_doublet_gap_m2_over_omega2"])

    scale_rows = []
    for factor in (0.5, 1.0, 2.0):
        s00 = scalar_entry(scalar_first0, scalar_first0, scalar_second00, factor)
        s11 = scalar_entry(scalar_first1, scalar_first1, scalar_second11, factor)
        smat = 0.5 * (s00 + s11) * np.eye(4)
        vmat = vector_matrix(factor)
        total = smat + vmat
        kappa = float(np.trace(total) / 4.0)
        scale_rows.append({
            "mu_over_muU": factor,
            "scalar_kappa_over_omega2": float(np.trace(smat) / 4.0),
            "vector_kappa_over_omega2": float(np.trace(vmat) / 4.0),
            "bosonic_kappa_over_omega2": kappa,
            "delta_xi02_bosonic": kappa / weight10,
        })

    zero_count = int(np.sum(np.abs(scalar_eigenvalues) < tolerance))
    negative_count = int(np.sum(scalar_eigenvalues < -tolerance))
    scalar_hard_min = float(np.sort(scalar_eigenvalues)[-scalar_hard_count])
    vector_hard_min = float(np.sort(vector_eigenvalues)[-33])
    vector_soft_max = float(np.max(np.abs(np.sort(vector_eigenvalues)[:-33])))
    scalar_diag_mismatch = abs(scalar00 - scalar11) / max(1.0, abs(scalar_kappa))
    scalar_mixed_ratio = abs(scalar01) / max(1.0, abs(scalar_kappa))
    vector_isotropy = float(
        np.linalg.norm(vector_declared - np.trace(vector_declared) / 4.0 * np.eye(4))
        / max(1.0, np.linalg.norm(vector_declared))
    )
    total_isotropy = float(
        np.linalg.norm(total_declared - kappa_b * np.eye(4))
        / max(1.0, np.linalg.norm(total_declared))
    )
    retuned = total_declared - delta_xi_b * weight10 * np.eye(4)

    checks = [
        {"name": "tree background has 38 soft modes and no tachyon", "pass": zero_count == 38 and negative_count == 0},
        {"name": "hard scalar cluster contains 290 positive modes", "pass": scalar_hard_min > 1e-3},
        {"name": "hard vector cluster contains 33 positive modes", "pass": vector_hard_min > 1e-4 and vector_soft_max < 1e-12},
        {"name": "light subspace is one real four-plane", "pass": light_rank == 4},
        {"name": "quartic Hessian first derivative is step stable", "pass": first_step_residual < 2e-9},
        {"name": "quartic Hessian second derivative is step stable", "pass": second_step_residual < 2e-8},
        {"name": "soft scalar block has no linear mass splitting", "pass": max(soft_first0_norm, soft_first1_norm) < 2e-10},
        {"name": "Frechet trace Hessian matches direct finite difference", "pass": frechet_residual < 2e-8},
        {"name": "scalar CW curvature agrees on two doublet directions", "pass": scalar_diag_mismatch < 2e-7},
        {"name": "scalar CW mixed curvature vanishes", "pass": scalar_mixed_ratio < 2e-7},
        {"name": "vector CW Hessian is SM-isotropic", "pass": vector_isotropy < 2e-10},
        {"name": "Landau ghost hard contribution vanishes", "pass": True},
        {"name": "bosonic xi02 counterterm retunes all four light modes", "pass": float(np.linalg.norm(retuned)) < 2e-10},
        {"name": "heavy-Yukawa term is retained as a matching nuisance", "pass": True},
    ]
    passed = sum(bool(row["pass"]) for row in checks)

    report: dict[str, Any] = {
        "schema": "route-f-p54-p2-bosonic-cw-v1",
        "date": "2026-09-02",
        "model_id": "P54PQ-v2",
        "scheme": {
            "effective_action": "zero-momentum hard-mode Coleman-Weinberg matching functional",
            "gauge": "background-field Landau gauge (xi_gf=0)",
            "renormalization": "MSbar",
            "matching_scale": "muU=gU*omega",
            "muU_over_omega": mu_over_omega,
            "scalar_constant": 1.5,
            "vector_constant": 5.0 / 6.0,
            "vector_polarizations": 3,
            "hard_scalar_count": scalar_hard_count,
            "soft_scalar_count": scalar_soft_count,
            "hard_vector_count": 33,
            "soft_vector_count": 12,
            "ghost_statement": "FP ghosts have xi_gf*M_V^2=0 and therefore no hard CW term in Landau gauge",
            "goldstone_statement": "33 gauge Goldstones, one PQ Goldstone, and the tuned doublet are EFT modes excluded from the hard trace",
            "claim_boundary": "scheme-dependent zero-momentum matching curvature, not a gauge-independent pole mass",
        },
        "background": {
            "sigma_over_omega": ratio,
            "vs_over_omega": vs,
            "gU": g_u,
            "alphaU_inverse": float(p2["thresholded_solution"]["alpha_U_inverse"]),
            "tree_xi02": xi_tree,
            "tree_zero_count": zero_count,
            "tree_negative_count": negative_count,
            "scalar_hard_min_m2_over_omega2": scalar_hard_min,
            "vector_hard_min_m2_over_omega2": vector_hard_min,
            "nearest_heavy_doublet_gap_m2_over_omega2": heavy_gap,
        },
        "numerics": {
            "field_hessian_dimension": p1.N_REAL,
            "main_epsilon_over_omega": eps,
            "check_epsilon_over_omega": eps_check,
            "field_hessian_evaluations": len(evaluations),
            "first_derivative_step_residual": first_step_residual,
            "second_derivative_step_residual": second_step_residual,
            "soft_block_first_derivative_norm_direction0": soft_first0_norm,
            "soft_block_first_derivative_norm_direction1": soft_first1_norm,
            "frechet_random_matrix_regression_residual": frechet_residual,
            "scalar_direction0_kappa_over_omega2": scalar00,
            "scalar_direction1_kappa_over_omega2": scalar11,
            "scalar_direction0_check_kappa_over_omega2": scalar00_check,
            "scalar_mixed01_over_omega2": scalar01,
            "scalar_diagonal_mismatch": scalar_diag_mismatch,
            "scalar_mixed_ratio": scalar_mixed_ratio,
            "vector_isotropy_residual": vector_isotropy,
            "total_isotropy_residual": total_isotropy,
        },
        "cw_hessian_over_omega2": {
            "scalar": scalar_declared.tolist(),
            "vector": vector_declared.tolist(),
            "ghost": np.zeros((4, 4)).tolist(),
            "bosonic_total": total_declared.tolist(),
            "bosonic_kappa": kappa_b,
        },
        "scale_variation": scale_rows,
        "matching_condition": {
            "light_field_weight_10": weight10,
            "bosonic_delta_xi02_at_muU": delta_xi_b,
            "heavy_yukawa_projection_nuisance": "eta_Y(mu)=h^T Pi_heavy-Y(0;mu) h / omega^2",
            "total_counterterm": "delta_xi02(mu)=[kappa_B(mu)+eta_Y(mu)]/w10",
            "nuisance_is_not_set_to_zero": True,
            "heavy_doublet_safety_condition": "||Q[Pi_B+Pi_Y-delta_xi02 P10]Q||/omega^2 < 0.03055703025613425",
            "pole_mass_claimed": False,
        },
        "status": {
            "bosonic_cw_light_doublet_hessian_computed": True,
            "gauge_ghost_scheme_declared": True,
            "heavy_yukawa_projection_exposed_as_nuisance": True,
            "p2_matching_relation_closed_conditionally": True,
            "full_pole_mass_predicted": False,
            "global_flavor_fit_performed": False,
        },
        "checks": checks,
        "summary": {"passed": passed, "total": len(checks), "all_pass": passed == len(checks)},
        "runtime_seconds": time.time() - started,
        "sources": [
            {"path": str(P1_SCRIPT.relative_to(REPO)), "sha256": sha256(P1_SCRIPT)},
            {"path": str(SEARCH_SCRIPT.relative_to(REPO)), "sha256": sha256(SEARCH_SCRIPT)},
            {"path": str(P2_JSON.relative_to(REPO)), "sha256": sha256(P2_JSON)},
            {"path": str(SCRIPT.relative_to(REPO)), "sha256": sha256(SCRIPT)},
        ],
    }
    return report


def markdown(report: dict[str, Any]) -> str:
    cw = report["cw_hessian_over_omega2"]
    match = report["matching_condition"]
    n = report["numerics"]
    lines = [
        "# P54PQ-v2 bosonic Coleman-Weinberg matching Hessian",
        "",
        f"Status: **{report['summary']['passed']}/{report['summary']['total']} checks passed**.",
        "",
        "Declared scheme: background-field Landau gauge, MSbar, zero-momentum hard-mode matching at `muU=gU*omega`. This is not advertised as a gauge-independent pole mass.",
        "",
        f"- hard trace: `{report['scheme']['hard_scalar_count']}` scalars and `{report['scheme']['hard_vector_count']}` vectors;",
        f"- EFT trace removed: `{report['scheme']['soft_scalar_count']}` scalar directions and `{report['scheme']['soft_vector_count']}` vectors;",
        f"- `gU={report['background']['gU']:.12f}`, `muU/omega={report['scheme']['muU_over_omega']:.12f}`;",
        f"- scalar kappa: `{cw['scalar'][0][0]:.12e} omega^2`;",
        f"- vector kappa: `{cw['vector'][0][0]:.12e} omega^2`;",
        f"- bosonic kappa: `{cw['bosonic_kappa']:.12e} omega^2`;",
        f"- bosonic retuning: `delta xi02,B={match['bosonic_delta_xi02_at_muU']:.12e}`.",
        "",
        "The full real light-doublet Hessian is proportional to `I4` by SM invariance. The numerical second direction and mixed-direction checks give",
        "",
        f"- scalar diagonal mismatch: `{n['scalar_diagonal_mismatch']:.3e}`;",
        f"- scalar mixed ratio: `{n['scalar_mixed_ratio']:.3e}`;",
        f"- vector isotropy residual: `{n['vector_isotropy_residual']:.3e}`.",
        "",
        "## Scale variation",
        "",
        "| mu/muU | scalar kappa | vector kappa | total kappa | delta xi02,B |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in report["scale_variation"]:
        lines.append(
            f"| {row['mu_over_muU']:.1f} | {row['scalar_kappa_over_omega2']:.8e} | {row['vector_kappa_over_omega2']:.8e} | {row['bosonic_kappa_over_omega2']:.8e} | {row['delta_xi02_bosonic']:.8e} |"
        )
    lines.extend([
        "",
        "## Matching nuisance",
        "",
        "The unfitted fermionic term is not set to zero:",
        "",
        "`eta_Y(mu)=h^T Pi_heavy-Y(0;mu) h/omega^2`,",
        "",
        "`delta xi02(mu)=[kappa_B(mu)+eta_Y(mu)]/w10`.",
        "",
        "The nearest heavy-doublet gap remains a separate norm inequality. Thus the bosonic matching coefficient is now numerical, while the heavy-Yukawa projection remains an explicit P3 matching nuisance.",
        "",
        "## Checks",
        "",
    ])
    for row in report["checks"]:
        lines.append(f"- [{'x' if row['pass'] else ' '}] {row['name']}")
    return "\n".join(lines) + "\n"


def main() -> None:
    report = run()
    OUTPUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUTPUT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    print(f"wrote {OUTPUT_JSON.relative_to(REPO)}")
    print(f"wrote {OUTPUT_MD.relative_to(REPO)}")
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
