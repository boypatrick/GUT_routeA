#!/usr/bin/env python3
"""Cache-only P54 kinetic re-attribution and transverse-gap obstruction.

The two rational projectors have exact integer algebra certificates. Their
identification with the tensor action is checked numerically against all 29
unit-invariant Hessians; it is NOT an interval/symbolic contraction proof.
The ensuing all-tau inequality is an analytic theorem conditional on those
displayed action identities and an admissible positive hard spectrum.
No new Hessian evaluation, parameter search, CW fit or benchmark mutation.
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
import verify_p54_nonuniform_reselection as nonuniform

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_transverse_gap_obstruction"
TRANSVERSE = ("lambda2", "lambda4", "lambda4p")
LOOP2 = 32 * np.pi**2


def run(cache_dir):
    oldpath = RF / "output/p54_full_doublet_cw.json"
    ledgerpath = RF / "output/p54_invariant_decomposition.json"
    old = json.loads(oldpath.read_text())
    ledger = json.loads(ledgerpath.read_text())
    if not ledger["summary"]["all_pass"]:
        raise ValueError("The prerequisite invariant reconstruction failed")
    p, r0 = old["tree_parameters"], np.array(ledger["reference_radii"])
    p1 = common.yuk.module("gap_obstruction_p1", RF / "code/verify_p54_p1_hessian_spectrum.py")
    search = common.yuk.module("gap_obstruction_search", RF / "code/search_p54_hierarchical_p1.py")
    rad = np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    metric = rad.T @ rad
    checks = []

    def check(name, actual, expected=0., tol=2e-10):
        actual, expected = np.asarray(actual), np.asarray(expected)
        residual = float(np.linalg.norm(actual - expected) /
                         max(1., np.linalg.norm(actual), np.linalg.norm(expected)))
        checks.append(dict(name=name, residual=residual, tolerance=tol,
                           passed=bool(np.isfinite(residual) and residual < tol)))

    # Refuse to compile/evaluate a Hessian even if a cache entry is missing.
    original = common.scalar.ActionJets.hessian
    cache_sources = {}

    def cached_only(store, x, label):
        key = hashlib.sha256(store.keybase + np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        paths = [d / (key + ".npz") for d in (store.read_cache, store.write_cache)]
        present = next((path for path in paths if path.exists()), None)
        if present is None:
            raise RuntimeError(label + ": missing cache; new Hessians are forbidden")
        cache_sources[key] = hashlib.sha256(present.read_bytes()).hexdigest()
        return original(store, x, label)

    common.scalar.ActionJets.hessian = cached_only
    try:
        inv = decomposition.InvariantJets(p1, p, r0, cache_dir)
    finally:
        common.scalar.ActionJets.hessian = original
    check("zero_new_Hessian_evaluations", inv.evaluated)
    h0, t0, _ = inv.assemble(p, r0)
    lam0, q0 = np.linalg.eigh(h0)
    p10 = np.zeros_like(h0)
    p10[306:326, 306:326] = np.eye(20)

    # A dyadic reconstruction supplies EXACT rational projectors. It does
    # not, by itself, prove their identification with the exact tensor action.
    blocks = []
    for rank, coeffs in ((50, (4., 4., 16.)), (36, (6., 6., 8.))):
        kappa = sum(c * p[name] for c, name in zip(coeffs, TRANSVERSE))
        selected = q0[:, abs(lam0 - kappa * r0[1]**2) < 1e-8]
        check(f"rank_{rank}_reference_cluster_dimension", selected.shape[1], rank)
        numerical = selected @ selected.T
        integer = np.rint(8 * numerical[54:306, 54:306]).astype(np.int64)
        projector = np.zeros_like(h0)
        projector[54:306, 54:306] = integer / 8
        check(f"rank_{rank}_rational_reconstruction", projector, numerical)
        # These boolean certificates use only exact integer arithmetic.
        exact = dict(symmetric=bool(np.array_equal(integer, integer.T)),
                     idempotent=bool(np.array_equal(integer @ integer, 8 * integer)),
                     trace_gives_rank=bool(np.trace(integer) == 8 * rank))
        check(f"rank_{rank}_exact_integer_projector_certificate", int(not all(exact.values())))
        check(f"rank_{rank}_annihilates_tuned_H10", p10 @ projector)
        expected = {name: 0. for name in inv.names}
        expected.update(nu2=-.5, lambda0=.5 * r0[1]**2,
                        alpha=6 / 5 * r0[0]**2,
                        beta=-6 / 5 * r0[0]**2, chi2=60 * r0[2]**2)
        for c, name in zip(coeffs, TRANSVERSE):
            expected[name] = c * r0[1]**2
        invariant_identities = []
        for a, name in enumerate(inv.names):
            actual = inv.units[a] @ projector
            predicted = expected[name] * projector
            check(f"rank_{rank}_unit_action_{name}", actual, predicted)
            invariant_identities.append(dict(invariant=name,
                expected_unit_Hessian_eigenvalue=expected[name],
                action_residual=float(np.linalg.norm(actual - predicted))))
        indices = np.argwhere(integer != 0)
        blocks.append(dict(rank=rank, kappa=kappa, transverse_coefficients=coeffs,
            projector=projector, integer=integer,
            public=dict(rank=rank, kappa=kappa, transverse_coefficients=coeffs,
                rational_projector={"sector": "canonical real Sigma coordinates 54:306",
                    "denominator": 8, "shape": [252, 252],
                    "nonzero_integer_entries": [[int(i), int(j), int(integer[i, j])] for i, j in indices]},
                exact_integer_certificate=exact, unit_action_identities=invariant_identities)))
    check("two_projectors_exact_integer_orthogonality",
          int(np.any(blocks[0]["integer"] @ blocks[1]["integer"])))
    transverse_h = sum(p[name] * inv.units[inv.names.index(name)] for name in TRANSVERSE)
    transverse_t = sum(p[name] * inv.component(inv.names.index(name), r0)[1] for name in TRANSVERSE)
    w, sigma, vs = r0
    b0 = np.array([12 / 5 * (p["alpha"] - p["beta"]) * w,
                   p["lambda0"] * sigma, 120 * p["chi2"] * vs])
    for block in blocks:
        rank, kappa, projector = block["rank"], block["kappa"], block["projector"]
        check(f"rank_{rank}_stationary_constant_mass_cancellation", (h0 - transverse_h) @ projector)
        check(f"rank_{rank}_transverse_mass_slope", transverse_h @ projector,
              kappa * sigma**2 * projector)
        for j in range(3):
            check(f"rank_{rank}_fixed_vertex_constant_{j}",
                  projector @ (t0[j] - transverse_t[j]) @ projector, b0[j] * projector)
            slope = 2 * kappa * sigma if j == 1 else 0.
            check(f"rank_{rank}_fixed_vertex_slope_{j}",
                  projector @ transverse_t[j] @ projector, slope * projector)

    def gram_lower(tau):
        result = np.zeros((3, 3))
        for block in blocks:
            b = b0 + np.array([0., 2 * block["kappa"] * tau * sigma, 0.])
            result += block["rank"] * np.outer(b, b) / (
                192 * np.pi**2 * block["kappa"] * tau * sigma**2)
        return result

    factor = sum(block["rank"] / block["kappa"] for block in blocks) / (192 * np.pi**2 * sigma**2)
    phi_constant = factor * b0[0]**2 / metric[0, 0]
    norm2 = float(b0 @ np.linalg.solve(metric, b0))
    fixed_witness = np.linalg.solve(metric, b0) / np.sqrt(norm2)
    stronger_constant = factor * norm2
    check("positive_witness_cross_term_for_all_tau", min(0., fixed_witness[1]))
    check("fixed_witness_canonical_unit_norm", fixed_witness @ metric @ fixed_witness, 1.)

    rows = []
    # Existing engineering cards only: this is re-attribution, not a scan.
    designs = (("reference", 1., 1., sigma), ("T10", .1, 1., sigma),
               ("T03", .03, 1., sigma), ("T03_X_S30", .03, .03, .3))
    for label, tau, chi4_factor, new_sigma in designs:
        pars = dict(p)
        for name in TRANSVERSE:
            pars[name] *= tau
        pars["chi4"] *= chi4_factor
        r = np.array([w, new_sigma, vs])
        h, t, _ = inv.assemble(pars, r)
        dp = -np.linalg.solve(common.radial_mass_map(r), p1.radial_gradient(r, pars))
        for name, value in zip(common.MASS_KEYS, dp):
            pars[name] += float(value)
        h += common.mass_operator(dp) + pars["xi02"] * p10
        pars["xi02"] = 0.
        _, sym, complement = search.symmetry_complement(p1, rad @ r)
        check(label + "_tree_stationarity", p1.radial_gradient(r, pars))
        check(label + "_tree_Ward_identity", h @ sym)
        tuning = search.tune_first_instability(h, p10, complement)
        if not tuning["tunable"]:
            raise ValueError(label + " has no admissible first crossing")
        pars["xi02"] = tuning["xi02"]
        h -= pars["xi02"] * p10
        lam, q = np.linalg.eigh(h)
        check(label + "_38_tree_soft_modes", int(np.sum(abs(lam) < tuning["zero_tolerance"])), 38)
        check(label + "_positive_hard_spectrum", int(lam[38] <= tuning["zero_tolerance"]))
        weight = nonuniform.bubble_slope_matrix(lam, 38)
        te = np.array([q.T @ v @ q for v in t])
        dz = np.einsum("ij,rij,sij->rs", weight, te, te) / LOOP2
        values, vectors = eigh(dz, metric)
        lead = vectors[:, -1]
        lead *= np.sign(lead[np.argmax(abs(lead))])
        total_vertex = np.einsum("r,rij->ij", lead, te)
        allocations = []
        for a, name in enumerate(inv.names):
            _, ta, _ = inv.component(a, r)
            vertex = q.T @ np.einsum("r,rij->ij", lead, pars[name] * ta) @ q
            allocation = float(np.sum(weight * vertex * total_vertex) / LOOP2)
            self_term = float(np.sum(weight * vertex**2) / LOOP2)
            allocations.append(dict(invariant=name, signed_allocation=allocation, self_term=self_term))
        allocations.sort(key=lambda row: abs(row["signed_allocation"]), reverse=True)
        check(label + "_all_signed_allocations_recombine", sum(a["signed_allocation"] for a in allocations), values[-1])
        row = dict(label=label, tau=tau, radii=r.tolist(), tree_parameters=pars,
            tree_min_hard_mass2=float(lam[38]), kinetic_matrix=dz.tolist(),
            kinetic_eigenvalues=values.tolist(), leading_radial_vector=lead.tolist(),
            leading_canonical_fractions=(np.diag(metric) * lead**2).tolist(),
            kinetic_signed_allocations=allocations,
            allocation_is_frozen_at_this_cards_own_mass_basis=True,
            hard_soft_leading_contribution=float(2 * np.sum(weight[:38, 38:] * total_vertex[:38, 38:]**2) / LOOP2))
        if new_sigma == sigma and chi4_factor == 1.:
            lower = gram_lower(tau)
            check(label + "_full_kinetic_dominates_two_projector_Gram",
                  min(0., eigh(dz - lower, metric, eigvals_only=True).min()))
            check(label + "_analytic_fixed_witness_inverse_tau_bound",
                  min(0., fixed_witness @ lower @ fixed_witness - stronger_constant / tau))
            row.update(two_projector_Gram=lower.tolist(),
                       two_projector_Gram_max_eigenvalue=float(eigh(lower, metric, eigvals_only=True)[-1]),
                       analytic_phi_bound=phi_constant / tau,
                       analytic_fixed_witness_bound=stronger_constant / tau)
        rows.append(row)
    for block in blocks:
        block.pop("projector"); block.pop("integer")
    sources = [Path(__file__), Path(p1.__file__), Path(search.__file__),
               Path(common.__file__), Path(common.scalar.__file__), Path(common.yuk.__file__),
               Path(decomposition.__file__), Path(nonuniform.__file__), oldpath, ledgerpath]
    return dict(schema="p54-transverse-gap-obstruction-v1", date="2026-09-14",
        scope="real radial scalar hard kinetic kernel; no pole/stability or physical fit certification",
        projectors=[block["public"] for block in blocks],
        reference_radii=r0.tolist(), radial_metric=metric.tolist(),
        analytic_bound=dict(fixed_vertex_b0=b0.tolist(), fixed_canonical_witness=fixed_witness.tolist(),
            phi_inverse_tau_constant=phi_constant, stronger_inverse_tau_constant=stronger_constant,
            margin=.3, tau_below_which_phi_bound_fails_margin=phi_constant / .3,
            tau_below_which_stronger_bound_fails_margin=stronger_constant / .3,
            domain="Every tau>0 on the stated transverse-only ray for which the retained hard spectrum is positive; the two projectors must be retained hard modes.",
            exact_integer_projector_algebra=True,
            tensor_action_identification="Numerically verified for every unit invariant, not independently interval/symbolically certified.",
            theorem_status="Analytic implication conditional on the displayed finite action-block identities; no unconditional model-wide exclusion."),
        cards=rows, new_Hessian_evaluations=inv.evaluated, unit_cache_hits=inv.hits,
        unit_cache_content_sha256=cache_sources,
        source_sha256={str(path.relative_to(RF)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        physical_fit_enabled=False, default_parameters_changed=False, complete_gauged_pole_kernel=False,
        checks=checks, summary=dict(passed=sum(c["passed"] for c in checks), total=len(checks),
                                  all_pass=all(c["passed"] for c in checks)))


def markdown(report):
    bound = report["analytic_bound"]
    lines = ["# P54 transverse-gap obstruction and new-basis kinetic attribution", "",
        "No new Hessian evaluations. These are existing cards, not a new parameter scan.", "",
        "Exact integer projector algebra is certified. Identification with the tensor action is numerically checked, not a symbolic/interval proof; the analytic all-tau implication is conditional on those identities.", "",
        f"Checks: {report['summary']['passed']}/{report['summary']['total']}.", "",
        "| Card | hard mass-squared gap | largest kinetic eigenvalue | canonical Phi / Sigma / S fractions |",
        "|---|---:|---:|---|"]
    for row in report["cards"]:
        lines.append(f"| {row['label']} | {row['tree_min_hard_mass2']:.9g} | {row['kinetic_eigenvalues'][-1]:.9g} | "
                     + " / ".join(f"{f:.6f}" for f in row["leading_canonical_fractions"]) + " |")
    lines += ["", "## Structural mechanism", "",
        "The stationary mass condition cancels the static portal contribution in two pure-Sigma spectral branches of real ranks 50 and 36; it does not cancel fixed-Lagrangian radial vertices.", "",
        "Their transverse coefficients are kappa_50=4 lambda2+4 lambda4+16 lambda4p and kappa_36=6 lambda2+6 lambda4+8 lambda4p. On the declared ray they are 12 tau and 11.6 tau, while T_w remains 0.36.", "",
        f"Consequently rho_Z >= {bound['phi_inverse_tau_constant']:.12g}/tau on the Phi witness, or >= {bound['stronger_inverse_tau_constant']:.12g}/tau on the displayed fixed canonical witness.",
        f"The latter excludes tau < {bound['tau_below_which_stronger_bound_fails_margin']:.12g} from the 0.3 engineering gate wherever the positive hard kernel is admissible. This is not a convergence theorem or physical pole exclusion.", "",
        "Further transverse-only weakening cannot uniformly fix this kinetic diagnostic. A coordinated portal/radial stiffness change, or an EFT retaining the shrinking branches, is the next distinct question; neither is established here."]
    for row in report["cards"]:
        lines += ["", "## " + row["label"] + " own-basis signed allocations", "",
                  "| Invariant | signed leading allocation | self term |", "|---|---:|---:|"]
        for item in row["kinetic_signed_allocations"][:10]:
            lines.append(f"| {item['invariant']} | {item['signed_allocation']:.9g} | {item['self_term']:.9g} |")
    lines += ["", "The JSON retains complete parameter cards, all allocations, both explicit rational projectors, every unit-action identity residual, and source/cache hashes.",
              "Gauge/ghost/fermion momentum, full Wilson/box, lower matching and finite seesaw remain incomplete. No default benchmark or physical fit is promoted."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    result = run(parser.parse_args().cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    OUT.with_suffix(".md").write_text(markdown(result))
    print(markdown(result))
    if not result["summary"]["all_pass"]:
        print([row for row in result["checks"] if not row["passed"]])
        raise SystemExit(1)
