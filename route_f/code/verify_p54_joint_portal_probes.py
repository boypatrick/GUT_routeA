#!/usr/bin/env python3
"""Two declared second-stage, joint scalar-portal engineering probes.

These are NOT an optimizer or a continuation of an experimental fit.
After the first six nonuniform cards failed, they jointly weaken the
existing Phi--Sigma and Sigma--S interactions and test two coarse radii.
Only cached unit-invariant Hessians and exact polynomial assembly are
permitted; this script fails before running if a required cache is absent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

import verify_p54_common_renormalization as common
import verify_p54_invariant_decomposition as decomposition
import verify_p54_nonuniform_reselection as candidates


RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_joint_portal_probes"


def require_cached_units(p1, parameters, radii, cache_dir):
    """Preflight all unit actions so InvariantJets cannot trigger new H."""
    x = p1.vacuum_vector(*radii)
    paths = []
    for name in decomposition.DEGREES:
        if name in ("xi1", "xi2"):
            continue
        unit = {key: 0. for key in parameters}
        unit[name] = 1.
        store = common.scalar.ActionJets(p1, unit, cache_dir)
        key = hashlib.sha256(store.keybase + np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        hits = [directory / (key + ".npz") for directory in (store.read_cache, store.write_cache)
                if (directory / (key + ".npz")).exists()]
        if not hits:
            raise FileNotFoundError(f"Zero-new-Hessian policy: missing cached unit invariant {name}")
        paths.append(hits[0])
    return paths


def run(cache_dir):
    oldpath = RF / "output/p54_full_doublet_cw.json"
    firststage = RF / "output/p54_nonuniform_reselection.json"
    ledgerpath = RF / "output/p54_invariant_decomposition.json"
    old = json.loads(oldpath.read_text())
    prior = json.loads(firststage.read_text())
    if prior["accepted_candidate_labels"]:
        raise ValueError("These second-stage cards are conditional on the first-stage rejection")
    p1 = common.yuk.module("joint_probe_p1", RF / "code/verify_p54_p1_hessian_spectrum.py")
    cw = common.yuk.module("joint_probe_cw", RF / "code/verify_p54_p2_bosonic_cw.py")
    search = common.yuk.module("joint_probe_search", RF / "code/search_p54_hierarchical_p1.py")
    original = old["tree_parameters"]
    reference_radii = np.array([old["vacuum"][key] for key in ("omega", "sigma", "vs")])
    mu = old["scheme"]["mu_over_omega"]
    rad = np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    metric = rad.T @ rad
    cached_paths = require_cached_units(p1, original, reference_radii, cache_dir)
    inv = decomposition.InvariantJets(p1, original, reference_radii, cache_dir)
    checks = []

    def check(name, value, expected=0., tolerance=3e-8):
        value, expected = np.asarray(value), np.asarray(expected)
        error = float(np.linalg.norm(value - expected) / max(1., np.linalg.norm(value), np.linalg.norm(expected)))
        checks.append({"name": name, "residual": error, "tolerance": tolerance,
                       "pass": bool(np.isfinite(error) and error < tolerance)})

    check("no_new_unit_Hessians", inv.evaluated)
    check("all_27_nonzero_unit_Hessians_cached", inv.hits, 27)
    rows = []
    # Predeclared bounded design: TWO cards, no fitted data, no adaptation
    # after seeing either loop result. lambda0=.2 is an absolute value.
    modified = dict(original)
    modified["lambda0"] = .2
    for key in ("lambda2", "lambda4", "lambda4p", "chi4"):
        modified[key] *= .01
    for key in ("alpha", "beta", "chi2"):
        modified[key] *= .1
    design = {
        "lambda0_absolute_value": .2,
        "multipliers": {"lambda2": .01, "lambda4": .01, "lambda4p": .01,
                        "alpha": .1, "beta": .1, "chi2": .1, "chi4": .01},
        "unchanged_independent_coefficients": [
            key for key in original if key not in
            {"lambda0", "lambda2", "lambda4", "lambda4p", "alpha", "beta", "chi2", "chi4",
             "mu2", "nu2", "mus2", "xi02"}],
        "dependent_reanchor_and_tune": ["mu2", "nu2", "mus2", "xi02"],
        "fixed_omega": 1., "fixed_vs": .25,
        "declared_sigma_values_before_candidate_evaluation": [.3, .5],
        "no_observed_mass_or_flavor_input": True,
    }
    for sigma in (.3, .5):
        label = f"joint_portals_sigma_{sigma:.1f}"
        r = np.array([1., sigma, .25])
        h, t, second = inv.assemble(modified, r)
        row = candidates.evaluate_candidate(
            p1, cw, search, modified, r, h, t, second, mu, rad, checks, label,
            gauge_coupling=mu,
        )
        row["design"] = design
        if "radial_hard_scalar_kinetic" in row:
            dz = np.asarray(row["radial_hard_scalar_kinetic"])
            eigenvalues, eigenvectors = eigh(dz, metric)
            leading = eigenvectors[:, -1]
            row["leading_kinetic_direction_radial_coordinates"] = leading.tolist()
            row["leading_kinetic_canonical_coordinate_fractions"] = (
                metric.diagonal() * leading**2).tolist()
            check(label + "_leading_kinetic_metric_norm", leading @ metric @ leading, 1.)
            check(label + "_leading_kinetic_eigenvalue", eigenvalues[-1], row["rho_Z"])
        # The candidate evaluator must leave the input independent card
        # untouched, regardless of its dependent stationary mass shifts.
        check(label + "_input_parameter_card_not_mutated",
              [modified[key] for key in original],
              [(.2 if key == "lambda0" else original[key] * design["multipliers"].get(key, 1.))
               for key in original])
        rows.append(row)
        print(label, {key: row.get(key) for key in
              ("failure_reason", "minimum_hard_radial_mass2", "rho_H", "rho_Z",
               "accepted_for_next_radial_audit")}, flush=True)
    sources = [Path(__file__), Path(decomposition.__file__), Path(candidates.__file__),
               Path(common.__file__), Path(common.scalar.__file__), Path(common.yuk.__file__),
               Path(p1.__file__), Path(cw.__file__), Path(search.__file__),
               oldpath, firststage, ledgerpath]
    return {
        "schema": "p54-two-coarse-joint-portal-engineering-probes-v1", "date": "2026-09-14",
        "scope": "declared second stage of two existing-action coupling/VEV cards; no optimization or data fit",
        "motivation": "first six nonuniform cards failed; decreasing transverse Sigma quartics alone can transfer the leading kinetic correction toward the Phi direction",
        "design": design,
        "reference_parameters": original, "reference_radii": reference_radii.tolist(),
        "engineering_margin": .3, "candidates": rows,
        "accepted_for_next_radial_audit": [row["label"] for row in rows if row["accepted_for_next_radial_audit"]],
        "selected_physical_benchmark": None, "default_parameters_changed": False,
        "physical_fit_enabled": False, "complete_mixed_BFB": False,
        "complete_gauged_pole_kernel": False, "complete_matching": False,
        "independent_new_full_Hessian_crosscheck_performed": False,
        "cache": {"unit_hits": inv.hits, "new_unit_Hessians": inv.evaluated,
                  "new_full_Hessians": 0, "preflighted_unit_count": len(cached_paths)},
        "source_sha256": {str(path.relative_to(RF)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in sources},
        "checks": checks,
        "summary": {"passed": sum(c["pass"] for c in checks), "total": len(checks),
                    "all_pass": all(c["pass"] for c in checks)},
    }


def markdown(report):
    lines = ["# Two coarse joint-portal engineering probes", "",
             "Explicit second-stage design, not an optimizer or fit. No default benchmark changes.", "",
             f"Checks: {report['summary']['passed']}/{report['summary']['total']}; new Hessians: 0.", "",
             "| sigma/omega | min hard radial mass2 | rho_H | rho_Z | Next radial audit only |",
             "|---:|---:|---:|---:|:---:|"]
    for row in report["candidates"]:
        if "rho_H" in row:
            lines.append(f"| {row['radii'][1]} | {row['minimum_hard_radial_mass2']:.9g} | {row['rho_H']:.9g} | {row['rho_Z']:.9g} | {row['accepted_for_next_radial_audit']} |")
        else:
            lines.append(f"| {row['radii'][1]} | tree gate rejected | - | - | False |")
    lines += ["", "The quartic parameter lambda0 is set to 0.2 (not multiplied by 0.2).",
              "lambda2/lambda4/lambda4p/chi4 are multiplied by 0.01; alpha/beta/chi2 by 0.1.",
              "All other independent scalar couplings are retained; three stationary masses and xi02 are solved anew.", "",
              "No additional full Hessian evaluation is used: cache preflight plus exact invariant continuation only.",
              "Passing a radial gate would still leave full mixed BFB, global vacua, nonradial stability, gauge/ghost/fermion momentum, loop doublet tuning, Wilson/box, lower matching and finite seesaw open."]
    for row in report["candidates"]:
        if row.get("failure_reason"):
            lines.append(f"- {row['label']}: {row['failure_reason']}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    arguments = parser.parse_args()
    report = run(arguments.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    OUT.with_suffix(".md").write_text(markdown(report))
    print(markdown(report))
    if not report["summary"]["all_pass"]:
        print("failed", [check for check in report["checks"] if not check["pass"]])
        raise SystemExit(1)
