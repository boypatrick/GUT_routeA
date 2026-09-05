#!/usr/bin/env python3
"""Reconstruct the P2 light-doublet overlaps needed by the P3 Yukawa map.

The P1/P2 ledgers only retained the total 10_H and 126_H weights.  P3 needs
the four separate coefficients in

    c = (phi_hol, phi_anti, Sigma_hol, Sigma_anti).

This script recomputes the fixed-point 328-real Hessian, removes the gauge and
PQ zero modes, diagonalises hypercharge on the remaining real four-plane, and
separates holomorphic from antiholomorphic canonical coordinates.  Only
absolute values are exported. The tree diagnostic uses a phase-uniform
operator-norm bound. These are not loop-corrected eigenvectors, and the
opposite conjugations in h phi* and f Sigma must be retained in the map.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
P1_SCRIPT = ROUTE_F / "code" / "verify_p54_p1_hessian_spectrum.py"
SEARCH_SCRIPT = ROUTE_F / "code" / "search_p54_hierarchical_p1.py"
P2_JSON = ROUTE_F / "output" / "p54_p2_two_site_matching.json"
HIERARCHY_JSON = ROUTE_F / "output" / "p54_hierarchical_p1_search.json"
OUTPUT_JSON = ROUTE_F / "output" / "p54_p3_light_doublet_overlaps.json"
OUTPUT_MD = ROUTE_F / "output" / "p54_p3_light_doublet_overlaps.md"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> dict[str, Any]:
    started = time.time()
    p1 = load_module("p54_p1_p3_overlap", P1_SCRIPT)
    search = load_module("p54_search_p3_overlap", SEARCH_SCRIPT)
    p2 = json.loads(P2_JSON.read_text(encoding="utf-8"))
    hierarchy = json.loads(HIERARCHY_JSON.read_text(encoding="utf-8"))

    sigma = 0.1265
    vs = float(hierarchy["vevs"]["vs"])
    parameters = search.parameters_at(p1, sigma, vs)
    parameters["xi02"] = float(p2["doublet_tuning"]["xi02"])
    x0 = p1.vacuum_vector(1.0, sigma, vs)
    hessian_fn = p1.jax.jit(p1.hessian(p1.potential_factory(parameters)))
    hessian = np.asarray(hessian_fn(p1.anp.asarray(x0)), dtype=float)
    hessian = (hessian + hessian.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    tolerance = float(p2["doublet_tuning"]["zero_tolerance"])

    symmetry, _, _ = search.symmetry_complement(p1, x0)
    zero = eigenvectors[:, np.abs(eigenvalues) < tolerance]
    symmetry_projector = symmetry @ np.linalg.pinv(symmetry, rcond=1.0e-11)
    residual = zero - symmetry_projector @ zero
    light_u, singular, _ = np.linalg.svd(residual, full_matrices=False)
    light_rank = int(np.sum(singular > 2.0e-8))
    light = light_u[:, :light_rank]
    if light_rank != 4:
        raise RuntimeError(f"expected one real Higgs four-plane, got {light_rank}")

    hypercharge = p1.representation_matrix(p1.sm_generators()["Y"][0])
    restricted_y = light.T @ hypercharge @ light
    y_values, y_vectors = np.linalg.eig(restricted_y.astype(complex))

    # p1 uses R_Y=-i Y: eigenvalue -i/2 has physical charge +1/2.
    # Hol/anti here refers to the stored scalar, not to the field appearing
    # in the Yukawa invariant: h couples phi*, whereas f couples Sigma.
    sector_rows: list[dict[str, float]] = []
    for index, value in enumerate(y_values):
        if abs(value.imag + 0.5) > 2.0e-8 or abs(value.real) > 2.0e-8:
            continue
        vector = light @ y_vectors[:, index]
        vector /= np.linalg.norm(vector)
        h_re = vector[p1.SL_H_RE]
        h_im = vector[p1.SL_H_IM]
        sigma_re = vector[p1.SL_SIGMA_RE]
        sigma_im = vector[p1.SL_SIGMA_IM]
        h_hol = (h_re + 1j * h_im) / math.sqrt(2.0)
        h_anti = (h_re - 1j * h_im) / math.sqrt(2.0)
        sigma_hol = (sigma_re + 1j * sigma_im) / math.sqrt(2.0)
        sigma_anti = (sigma_re - 1j * sigma_im) / math.sqrt(2.0)
        sector_rows.append(
            {
                "H10u_sq": float(np.vdot(h_hol, h_hol).real),
                "eps_H10d_star_sq": float(np.vdot(h_anti, h_anti).real),
                "H126u_sq": float(np.vdot(sigma_hol, sigma_hol).real),
                "eps_H126d_star_sq": float(np.vdot(sigma_anti, sigma_anti).real),
            }
        )
    if len(sector_rows) != 2:
        raise RuntimeError(f"expected two Y=-1/2 weak components, got {len(sector_rows)}")

    mean_squares = {
        key: float(np.mean([row[key] for row in sector_rows]))
        for key in sector_rows[0]
    }
    spread = {
        key: float(np.ptp([row[key] for row in sector_rows]))
        for key in sector_rows[0]
    }
    coefficients = {key.removesuffix("_sq"): math.sqrt(value) for key, value in mean_squares.items()}

    # In the standard two-matrix sum-rule notation:
    # alpha_{1,2} multiply the down-type VEV and beta_{1,2} the up-type VEV.
    alpha1 = coefficients["H10u"]
    beta1 = coefficients["eps_H10d_star"]
    alpha2 = coefficients["eps_H126d_star"]
    beta2 = coefficients["H126u"]
    r_abs = beta1 / alpha1
    s_abs = alpha1 * beta2 / (alpha2 * beta1)
    weight10 = mean_squares["H10u_sq"] + mean_squares["eps_H10d_star_sq"]
    weight126 = mean_squares["H126u_sq"] + mean_squares["eps_H126d_star_sq"]
    expected_weight10 = float(
        json.loads((ROUTE_F / "output" / "p54_p2_running_thresholds_cosmology.json").read_text())[
            "doublet_retuning"
        ]["light_field_weight_10"]
    )

    checks = [
        {
            "name": "one complex SM doublet is recovered",
            "pass": light_rank == 4 and len(sector_rows) == 2,
        },
        {
            "name": "two weak components have identical sector weights",
            "pass": max(spread.values()) < 2.0e-13,
        },
        {
            "name": "four overlap weights are normalized",
            "pass": abs(sum(mean_squares.values()) - 1.0) < 2.0e-13,
        },
        {
            "name": "10_H weight agrees with the P2 ledger",
            "pass": abs(weight10 - expected_weight10) < 2.0e-12,
        },
        {
            "name": "126_H weight is the complement",
            "pass": abs(weight126 - (1.0 - expected_weight10)) < 2.0e-12,
        },
    ]
    passed = sum(bool(row["pass"]) for row in checks)
    report: dict[str, Any] = {
        "schema": "route-f-p54-p3-light-overlaps-v2",
        "date": "2026-09-05",
        "model_id": "P54PQ-v2",
        "basis": {
            "light_doublet_definition": "geometric copy order: phi_hol, phi_anti, Sigma_hol, Sigma_anti",
            "hypercharge_eigenspace": "R_Y eigenvalue -i/2, physical Y=+1/2",
            "legacy_key_warning": "coefficient labels H10u etc are historical geometric keys; use geometric_sector_magnitudes and corrected sum_rule_overlap_map",
            "yukawa_conjugation": "h couples phi*, f couples Sigma; up overlaps are phi_anti and Sigma_hol",
            "loop_order": "tree Hessian only; no loop-corrected light eigenvector is computed",
            "phase_policy": "absolute overlaps only; the downstream no-go bound is uniform in all relative phases",
        },
        "coefficient_magnitudes": {
            "abs_c1_H10u": coefficients["H10u"],
            "abs_c2_eps_H10d_star": coefficients["eps_H10d_star"],
            "abs_c3_H126u": coefficients["H126u"],
            "abs_c4_eps_H126d_star": coefficients["eps_H126d_star"],
        },
        "coefficient_squares": mean_squares,
        "geometric_sector_magnitudes": {
            "phi_hol": coefficients["H10u"],
            "phi_anti": coefficients["eps_H10d_star"],
            "Sigma_hol": coefficients["H126u"],
            "Sigma_anti": coefficients["eps_H126d_star"],
        },
        "weak_component_spread": spread,
        "sum_rule_overlap_map": {
            "abs_alpha1_down_10": alpha1,
            "abs_beta1_up_10": beta1,
            "abs_alpha2_down_126": alpha2,
            "abs_beta2_up_126": beta2,
            "abs_r_beta1_over_alpha1": r_abs,
            "abs_s_alpha1beta2_over_alpha2beta1": s_abs,
        },
        "field_weights": {"10_H": weight10, "126_H": weight126},
        "checks": checks,
        "summary": {"passed": passed, "total": len(checks), "all_pass": passed == len(checks)},
        "runtime_seconds": time.time() - started,
        "sources": [
            {"path": str(P1_SCRIPT.relative_to(REPO)), "sha256": sha256(P1_SCRIPT)},
            {"path": str(SEARCH_SCRIPT.relative_to(REPO)), "sha256": sha256(SEARCH_SCRIPT)},
            {"path": str(P2_JSON.relative_to(REPO)), "sha256": sha256(P2_JSON)},
            {"path": str(HIERARCHY_JSON.relative_to(REPO)), "sha256": sha256(HIERARCHY_JSON)},
            {"path": str(Path(__file__).resolve().relative_to(REPO)), "sha256": sha256(Path(__file__))},
        ],
    }
    return report


def markdown(report: dict[str, Any]) -> str:
    c = report["coefficient_magnitudes"]
    m = report["sum_rule_overlap_map"]
    lines = [
        "# P54 P3 light-doublet overlap reconstruction",
        "",
        "The tree zero mode at the P2 background gives (legacy geometric c labels)",
        "",
        "```text",
        f"|c1| = {c['abs_c1_H10u']:.12g}",
        f"|c2| = {c['abs_c2_eps_H10d_star']:.12g}",
        f"|c3| = {c['abs_c3_H126u']:.12g}",
        f"|c4| = {c['abs_c4_eps_H126d_star']:.12g}",
        f"|r|  = {m['abs_r_beta1_over_alpha1']:.12g}",
        f"|s|  = {m['abs_s_alpha1beta2_over_alpha2beta1']:.12g}",
        "```",
        "",
        "Corrected 2026-09-05: R_Y=-iY, so the selected mode has Y=+1/2. Since h couples phi* and f couples Sigma, alpha1=c1, beta1=c2, alpha2=c4, beta2=c3. Only magnitudes are exported. Full loop eigenvectors remain uncomputed.",
        "",
        f"Checks: **{report['summary']['passed']}/{report['summary']['total']}**.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    report = run()
    OUTPUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUTPUT_MD.write_text(markdown(report), encoding="utf-8")
    print(f"P54 P3 light overlaps: {report['summary']['passed']}/{report['summary']['total']}")
    print(f"wrote {OUTPUT_JSON.relative_to(REPO)}")
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
