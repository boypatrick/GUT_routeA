#!/usr/bin/env python3
"""Search the original P54PQ-v2 coupling neighbourhood at hierarchical sigma/omega.

The script is intentionally a search/diagnostic rather than a promotion
verifier.  It imports the exact P1 tensor action, derives the three radial
mass parameters at the requested VEVs, and computes one complete 328-real
Hessian at xi02=0.  Because xi02 enters as -xi02 phi^dagger phi, subsequent
doublet tuning is an exact linear rank-20 update and requires no re-fit.
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
from scipy.linalg import null_space
from scipy.optimize import brentq


REPO = Path(__file__).resolve().parents[2]
P1_SCRIPT = REPO / "route_f/code/verify_p54_p1_hessian_spectrum.py"
P2_JSON = REPO / "route_f/output/p54_p2_running_thresholds_cosmology.json"
OUTPUT = REPO / "route_f/output/p54_hierarchical_p1_search.json"
OUTPUT_MD = REPO / "route_f/output/p54_hierarchical_p1_search.md"


def load_p1():
    spec = importlib.util.spec_from_file_location("p54_p1", P1_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load P1 verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameters_at(p1, sigma: float, vs: float) -> dict[str, float | complex]:
    p, _ = p1.benchmark()
    p = dict(p)
    omega = 1.0
    p["xi02"] = 0.0
    radial_a = 36.0 / 25.0 * p["a"] + 42.0 / 125.0 * p["b"]
    d = p["alpha"] - p["beta"]
    p["mu2"] = p["c"] * omega / 5.0 + (5.0 / 3.0) * radial_a * omega**2 + d * sigma**2 + p["chi3"] * vs**2
    p["nu2"] = p["lambda0"] * sigma**2 + 12.0 / 5.0 * d * omega**2 + 120.0 * p["chi2"] * vs**2
    p["mus2"] = p["chi1"] * vs**2 + 120.0 * p["chi2"] * sigma**2 + 12.0 / 5.0 * p["chi3"] * omega**2
    return p


def symmetry_complement(p1, x0: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    orbit = p1.gauge_orbit(x0)
    qvec = p1.pq_direction(x0)
    qphysical = qvec - orbit @ np.linalg.pinv(orbit, rcond=1e-11) @ qvec
    symmetry = np.column_stack([orbit, qphysical])
    u, singular, _ = np.linalg.svd(symmetry, full_matrices=False)
    rank = int(np.sum(singular > 2e-10))
    sym_basis = u[:, :rank]
    complement = null_space(sym_basis.T, rcond=2e-11)
    return symmetry, sym_basis, complement


def tune_first_instability(
    h0: np.ndarray, projector10: np.ndarray, complement: np.ndarray
) -> dict[str, Any]:
    h0p = complement.T @ h0 @ complement
    pp = complement.T @ projector10 @ complement
    h0p = (h0p + h0p.T) / 2.0
    pp = (pp + pp.T) / 2.0
    eig0 = np.linalg.eigvalsh(h0p)

    def minimum(xi: float) -> float:
        return float(np.linalg.eigvalsh(h0p - xi * pp)[0])

    if eig0[0] <= 1e-9:
        return {
            "tunable": False,
            "reason": "xi02=0 physical Hessian is already nonpositive",
            "minimum_at_xi02_zero": float(eig0[0]),
            "negative_count_at_xi02_zero": int(np.sum(eig0 < -1e-8)),
        }
    high = 0.05
    while minimum(high) > 0 and high < 10.0:
        high *= 2.0
    if high >= 10.0 and minimum(high) > 0:
        return {
            "tunable": False,
            "reason": "no physical crossing below xi02=10",
            "minimum_at_xi02_zero": float(eig0[0]),
        }
    xi = brentq(minimum, 0.0, high, xtol=2e-13, rtol=2e-13)
    tuned = h0p - xi * pp
    eig = np.linalg.eigvalsh(tuned)
    tolerance = max(2e-8, 3e-7 * max(1.0, float(np.max(np.abs(eig)))))
    return {
        "tunable": True,
        "xi02": float(xi),
        "minimum_at_xi02_zero": float(eig0[0]),
        "zero_tolerance": tolerance,
        "physical_zero_count": int(np.sum(np.abs(eig) < tolerance)),
        "physical_negative_count": int(np.sum(eig < -tolerance)),
        "nearest_positive_m2": float(eig[np.where(eig > tolerance)[0][0]]),
        "maximum_m2": float(eig[-1]),
    }


def run() -> dict[str, Any]:
    started = time.time()
    p1 = load_p1()
    p2 = json.loads(P2_JSON.read_text(encoding="utf-8"))
    sigma = float(p2["current_no_threshold_baseline"]["MI_over_MU"])
    vs = 0.25
    p = parameters_at(p1, sigma, vs)
    x0 = p1.vacuum_vector(1.0, sigma, vs)
    potential = p1.potential_factory(p)
    gradient = np.asarray(jax.jit(p1.grad(potential))(p1.anp.asarray(x0)), dtype=float)
    h_start = time.time()
    h0 = np.asarray(jax.jit(p1.hessian(potential))(p1.anp.asarray(x0)), dtype=float)
    hessian_seconds = time.time() - h_start
    h0 = (h0 + h0.T) / 2.0

    symmetry, _, complement = symmetry_complement(p1, x0)
    projector10 = np.zeros((p1.N_REAL, p1.N_REAL))
    projector10[p1.SL_H_RE, p1.SL_H_RE] = np.eye(p1.SL_H_RE.stop - p1.SL_H_RE.start)
    projector10[p1.SL_H_IM, p1.SL_H_IM] = np.eye(p1.SL_H_IM.stop - p1.SL_H_IM.start)
    tuning = tune_first_instability(h0, projector10, complement)

    classification: dict[str, Any] = {}
    scalar_spectrum: dict[str, Any] = {}
    if tuning.get("tunable"):
        xi = float(tuning["xi02"])
        h = h0 - xi * projector10
        eig, eigvec = np.linalg.eigh(h)
        extra = p1.classify_extra_zero_modes(
            eig,
            eigvec,
            float(tuning["zero_tolerance"]),
            symmetry,
        )
        classification = {
            "extra_zero_modes": extra,
            "full_negative_count": int(np.sum(eig < -float(tuning["zero_tolerance"]))),
            "full_zero_count": int(np.sum(np.abs(eig) < float(tuning["zero_tolerance"]))),
            "first_crossing_is_one_doublet": bool(
                extra["rank"] == 4
                and max(abs(value) for value in extra["casimir_eigenvalues"]["SU3"]) < 3e-6
                and max(abs(value - 0.75) for value in extra["casimir_eigenvalues"]["SU2"]) < 3e-6
                and max(abs(value - 0.25) for value in extra["casimir_eigenvalues"]["Y"]) < 3e-6
            ),
        }
        scalar_spectrum = {
            "eigenvalue_groups": p1.group_eigenvalues(
                eig, abs_tol=float(tuning["zero_tolerance"])
            ),
            "sm_irrep_ledger": p1.classify_full_spectrum(
                eig,
                eigvec,
                float(tuning["zero_tolerance"]),
                p1.sm_representation_matrices(),
            ),
        }

    vector_m2 = np.linalg.eigvalsh(p1.gauge_orbit(x0).T @ p1.gauge_orbit(x0))
    vector_spectrum = p1.group_eigenvalues(vector_m2, abs_tol=2e-10)

    points = p1.stationary_points(p)
    target = min(
        points,
        key=lambda row: (row["omega"] - 1.0) ** 2
        + (abs(row["sigma"]) - sigma) ** 2
        + (abs(row["vs"]) - vs) ** 2,
    )
    radial_global = bool(
        target["V"] <= points[0]["V"] + 2e-9
        and min(target["radial_hessian_eigenvalues"]) > 1e-8
    )
    accepted = bool(
        np.linalg.norm(gradient) < 2e-8
        and radial_global
        and tuning.get("tunable", False)
        and classification.get("full_negative_count") == 0
        and classification.get("full_zero_count") == 38
        and classification.get("first_crossing_is_one_doublet")
    )
    return {
        "schema": "route-f-p54-hierarchical-p1-search-v1",
        "date": "2026-08-30",
        "status": {
            "accepted_hierarchical_p1": accepted,
            "search_scope": "original_P1_couplings_with_radial_masses_rederived",
            "next_action": (
                "parent-resolve spectrum and rerun P2"
                if accepted
                else "identify unstable irrep and vary only the responsible invariant couplings"
            ),
        },
        "vevs": {"omega": 1.0, "sigma": sigma, "vs": vs},
        "derived_mass_parameters_at_xi02_zero": {
            key: float(p[key]) for key in ("mu2", "nu2", "mus2")
        },
        "gradient_norm": float(np.linalg.norm(gradient)),
        "symmetry_rank": int(np.linalg.matrix_rank(symmetry, tol=2e-10)),
        "physical_complement_dimension": int(complement.shape[1]),
        "tuning": tuning,
        "classification": classification,
        "scalar_spectrum": scalar_spectrum,
        "vector_spectrum_over_g2": vector_spectrum,
        "radial": {
            "root_count": len(points),
            "target": target,
            "minimum": points[0],
            "target_is_radial_global": radial_global,
        },
        "parameters": {
            key: ([value.real, value.imag] if isinstance(value, complex) else value)
            for key, value in p.items()
        },
        "runtime": {
            "hessian_seconds": hessian_seconds,
            "total_seconds": time.time() - started,
        },
        "sources": [
            {"path": str(P1_SCRIPT.relative_to(REPO)), "sha256": sha256(P1_SCRIPT)},
            {"path": str(P2_JSON.relative_to(REPO)), "sha256": sha256(P2_JSON)},
            {"path": str(Path(__file__).resolve().relative_to(REPO)), "sha256": sha256(Path(__file__).resolve())},
        ],
    }


def markdown(report: dict[str, Any]) -> str:
    tuning = report["tuning"]
    lines = [
        "# Hierarchical P54PQ P1 search",
        "",
        f"Status: **{'ACCEPTED' if report['status']['accepted_hierarchical_p1'] else 'NOT ACCEPTED'}** within the declared search scope.",
        "",
        "This is a hierarchy-restoration diagnostic: the dimensionless quartic/cubic couplings are unchanged from P1, while the three radial quadratic masses are re-derived at the RGE-selected vacuum.",
        "",
        f"- `sigma/omega = {report['vevs']['sigma']:.12g}`",
        f"- `xi02 = {tuning.get('xi02', float('nan')):.12g}`",
        f"- full gradient norm: `{report['gradient_norm']:.3e}`",
        f"- physical zeros / negatives: `{tuning.get('physical_zero_count', 'n/a')}` / `{tuning.get('physical_negative_count', 'n/a')}`",
        f"- first crossing is one SM doublet: `{report['classification'].get('first_crossing_is_one_doublet', False)}`",
        f"- target is the lowest enumerated radial branch: `{report['radial']['target_is_radial_global']}`",
        "",
        "## Scalar eigenvalue groups",
        "",
        "| `m^2/omega^2` | real multiplicity |",
        "|---:|---:|",
    ]
    for row in report.get("scalar_spectrum", {}).get("eigenvalue_groups", []):
        lines.append(f"| {row['m2']:.9e} | {row['multiplicity']} |")
    lines.extend([
        "",
        "## Complete SM-Casimir ledger",
        "",
        "| `m^2/omega^2` | SU(3) | SU(2) dim | `|Y|` | real mult. |",
        "|---:|---|---:|---:|---:|",
    ])
    for row in report.get("scalar_spectrum", {}).get("sm_irrep_ledger", {}).get("rows", []):
        lines.append(
            f"| {row['m2']:.9e} | {row['SU3']} | {row['SU2_dimension']} | "
            f"{row['abs_hypercharge']:.6g} | {row['real_multiplicity']} |"
        )
    lines.extend([
        "",
        "## Vector eigenvalue groups",
        "",
        "The entries are `M_V^2/(g^2 omega^2)`.",
        "",
        "| value | multiplicity |",
        "|---:|---:|",
    ])
    for row in report["vector_spectrum_over_g2"]:
        lines.append(f"| {row['m2']:.9e} | {row['multiplicity']} |")
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "Acceptance proves stationarity, radial ordering, Hessian semipositivity and the single-doublet crossing for this benchmark. It does not by itself assign every SM row to a Pati--Salam parent or complete two-site one-loop matching.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    report = run()
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUTPUT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["status"], sort_keys=True))
    print(json.dumps(report["tuning"], sort_keys=True))


if __name__ == "__main__":
    main()
