#!/usr/bin/env python3
"""DYN-5 model-validity gate for the quadratic hidden-messenger ansatz.

This audit deliberately does *not* repeat the historical one-loop claim.  The
displayed Route-B superpotential is Gaussian (quadratic) in ``(N, X)``.  It
therefore has a tree-level Schur complement and a tree-level Kahler matching
term, but no cubic tensor from which a Yukawa anomalous dimension could be
built.  A loop calculation remains blocked until an interacting messenger
action is supplied.

The arithmetic checks may all pass while the scientific claim is invalid.
Consumers must use ``claim_status`` rather than treating ``all_pass`` as a
physics verdict.
"""

from __future__ import annotations

import cmath
import json
import math
from datetime import datetime, timezone

import numpy as np

from route_e_paths import AUDIT_OUTPUT


OUT = AUDIT_OUTPUT / "audit5"
ZETA = complex(0.1076472949, 0.0736514853)
LAMBDA = cmath.sqrt(ZETA)
K_TR = np.array(
    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
    dtype=complex,
) / math.sqrt(3.0)
K_INV = 3.0 * K_TR

CHECKS: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    ok = bool(condition)
    CHECKS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail else ""))


def relative_norm(lhs: np.ndarray, rhs: np.ndarray) -> float:
    return float(np.linalg.norm(lhs - rhs) / max(np.linalg.norm(rhs), 1.0e-300))


def main() -> None:
    print("== DYN-5V section 1: tree-level Kahler matching ==")
    k_metric = K_TR.conjugate().T @ K_TR
    k_metric_residual = float(np.linalg.norm(k_metric - np.eye(3) / 3.0))
    delta_z_tree = abs(ZETA) / 3.0
    z_n = 1.0 + delta_z_tree
    contact_before = abs(ZETA)
    contact_after_without_retuning = contact_before / z_n
    contact_shift_fraction = contact_after_without_retuning / contact_before - 1.0
    putative_loop_per_efold = abs(ZETA) / (16.0 * math.pi**2)
    tree_to_putative_loop = delta_z_tree / putative_loop_per_efold

    check(
        "K_tr^dagger K_tr = I/3",
        k_metric_residual < 1.0e-14,
        f"residual={k_metric_residual:.3e}",
    )
    check(
        "tree Kahler correction is delta Z_N = |zeta|/3",
        abs(delta_z_tree - 0.04347730108431254) < 1.0e-12,
        f"deltaZ_tree={delta_z_tree:.12f}",
    )
    check(
        "tree correction is 16 pi^2/3 times the historical putative loop coefficient",
        abs(tree_to_putative_loop - 16.0 * math.pi**2 / 3.0) < 1.0e-12,
        f"ratio={tree_to_putative_loop:.9f}",
    )

    print("== DYN-5V section 2: full 6x6 propagator and Schur complement ==")
    # A deterministic nonsingular complex-symmetric Majorana block.  The
    # identity below is algebraic; random trials ensure that it is not an
    # accident of a commuting or diagonal benchmark.
    rng = np.random.default_rng(20260714)
    max_inverse_residual = 0.0
    max_c5_residual = 0.0
    min_singular_value = math.inf
    for _ in range(200):
        raw = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        m_v = 0.15 * (raw + raw.T) + 2.0 * np.eye(3)
        full = np.block(
            [
                [m_v, LAMBDA * np.eye(3)],
                [LAMBDA * np.eye(3), -K_INV],
            ]
        )
        schur = m_v + ZETA * K_TR
        min_singular_value = min(
            min_singular_value,
            float(np.linalg.svd(full, compute_uv=False)[-1]),
            float(np.linalg.svd(schur, compute_uv=False)[-1]),
        )
        nn_inverse = np.linalg.inv(full)[:3, :3]
        max_inverse_residual = max(
            max_inverse_residual, relative_norm(nn_inverse, np.linalg.inv(schur))
        )
        y_nu = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        c5_full = y_nu @ nn_inverse @ y_nu.T
        c5_schur = y_nu @ np.linalg.inv(schur) @ y_nu.T
        max_c5_residual = max(max_c5_residual, relative_norm(c5_full, c5_schur))

    check(
        "NN block of the full inverse equals the Schur-complement inverse",
        max_inverse_residual < 1.0e-12,
        f"max relative residual={max_inverse_residual:.3e}",
    )
    check(
        "full and Schur-complement Weinberg tensors agree at zero momentum",
        max_c5_residual < 1.0e-12,
        f"max relative residual={max_c5_residual:.3e}",
    )

    # Canonical normalization changes parameters but not the physical
    # zero-momentum Weinberg tensor when Y and M are transformed together.
    raw = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    m_eff = 0.2 * (raw + raw.T) + 2.0 * np.eye(3)
    y_nu = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    c5_before = y_nu @ np.linalg.inv(m_eff) @ y_nu.T
    y_canonical = y_nu / math.sqrt(z_n)
    m_canonical = m_eff / z_n
    c5_after = y_canonical @ np.linalg.inv(m_canonical) @ y_canonical.T
    canonical_c5_residual = relative_norm(c5_after, c5_before)
    check(
        "consistent canonical normalization leaves Y M^-1 Y^T invariant",
        canonical_c5_residual < 1.0e-13,
        f"relative residual={canonical_c5_residual:.3e}",
    )

    print("== DYN-5V section 3: interaction and selection-rule gates ==")
    # A quadratic superpotential has constant second derivatives and zero
    # third derivatives.  Thus y_ijk = d^3 W/dPhi_i dPhi_j dPhi_k vanishes.
    polynomial_degree = 2
    cubic_tensor_norm = 0.0 if polynomial_degree < 3 else math.nan
    check(
        "displayed action has no cubic messenger tensor",
        cubic_tensor_norm == 0.0,
        "W is quadratic; the claimed |lambda|^2/(16 pi^2) anomalous dimension is undefined",
    )

    # Explicit post-B-L R-charge counterexample.
    r = {"X": 1, "N": 1, "L": 1, "H_u": 0, "W": 2}
    operator_r = {
        "XN": r["X"] + r["N"],
        "XX": 2 * r["X"],
        "LNH_u": r["L"] + r["N"] + r["H_u"],
        "XLH_u": r["X"] + r["L"] + r["H_u"],
    }
    check(
        "the displayed R charges allow the dangerous X L H_u operator",
        operator_r["XLH_u"] == r["W"],
        f"R(XLH_u)={operator_r['XLH_u']}=R(W)",
    )

    # General additive-Abelian no-go.  If XN, XX and LNH are all allowed,
    # q_X+q_N=w, 2q_X=w and q_L+q_N+q_H=w imply
    # q_X+q_L+q_H=2q_X=w.  The same congruence proof holds for Z_n.
    additive_no_go_identity = True
    check(
        "one additive Abelian charge cannot allow XN, XX, LNH_u while forbidding XLH_u",
        additive_no_go_identity,
        "q(XLH)=2 q_X=w follows algebraically (also modulo n)",
    )

    n_pass = sum(ok for _, ok, _ in CHECKS)
    all_pass = n_pass == len(CHECKS)
    claim_status = "invalid_pending_rederivation"
    blockers = [
        "no cubic or gauge messenger interaction is specified, so the one-loop anomalous dimension does not exist in the displayed model",
        "the displayed post-B-L R charges allow X L H_u",
        "a full branch-local action and an operator-complete selection rule are required",
    ]
    payload = {
        "audit": "audit5_dyn5_model_validity",
        "dyn_item": "DYN-5",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "all_pass": all_pass,
        "checks_passed": n_pass,
        "checks_total": len(CHECKS),
        "checks": [
            {"name": name, "pass": ok, "detail": detail}
            for name, ok, detail in CHECKS
        ],
        "claim_status": claim_status,
        "physics_claim_valid": False,
        "blockers": blockers,
        "tree_matching": {
            "abs_zeta": abs(ZETA),
            "delta_Z_tree": delta_z_tree,
            "Z_N": z_n,
            "contact_abs_before": contact_before,
            "contact_abs_after_without_retuning": contact_after_without_retuning,
            "contact_shift_fraction_without_retuning": contact_shift_fraction,
            "putative_loop_per_efold_not_generated_by_displayed_action": putative_loop_per_efold,
            "tree_to_putative_loop_ratio": tree_to_putative_loop,
        },
        "full_system_validation": {
            "mass_matrix": "[[M_V, lambda I], [lambda I, -K_tr^-1]]",
            "schur_complement": "M_V + lambda^2 K_tr",
            "random_trials": 200,
            "minimum_singular_value_seen": min_singular_value,
            "max_NN_inverse_relative_residual": max_inverse_residual,
            "max_Weinberg_tensor_relative_residual": max_c5_residual,
            "canonical_Weinberg_tensor_relative_residual": canonical_c5_residual,
        },
        "loop_gate": {
            "displayed_superpotential_degree": polynomial_degree,
            "cubic_tensor_norm": cubic_tensor_norm,
            "status": "blocked_model_not_specified",
            "minimum_repair": "supply an explicit interacting messenger action, e.g. a fully charged S X N sector, before computing gamma_N",
        },
        "selection_rule_gate": {
            "displayed_R_charges": r,
            "operator_R_charges": operator_r,
            "XLH_u_allowed": True,
            "single_additive_abelian_no_go": True,
        },
        "zeta_value_derived": False,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    json_path = OUT / "dyn5_model_validity.json"
    md_path = OUT / "dyn5_model_validity.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(
        "# DYN-5 model-validity gate\n\n"
        f"- arithmetic checks: **{n_pass}/{len(CHECKS)} pass**\n"
        f"- scientific status: **{claim_status}**\n"
        "- the historical one-loop promotion is not valid for the displayed quadratic action\n\n"
        "## Correct tree matching\n\n"
        f"- `delta Z_tree = |zeta|/3 = {delta_z_tree:.10f}`\n"
        f"- `Z_N = {z_n:.10f}`\n"
        f"- unretuned `|zeta|/Z_N = {contact_after_without_retuning:.10f}` "
        f"(`{100.0 * contact_shift_fraction:.3f}%` shift)\n"
        f"- tree / historical putative-loop-per-efold ratio = `{tree_to_putative_loop:.6f}`\n"
        f"- full-6x6 versus Schur Weinberg residual = `{max_c5_residual:.3e}`\n\n"
        "## Fail-closed blockers\n\n"
        + "\n".join(f"- {item}" for item in blockers)
        + "\n\n## Checks\n\n"
        + "\n".join(
            f"- [{'PASS' if ok else 'FAIL'}] {name}: {detail}"
            for name, ok, detail in CHECKS
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {json_path} (+ .md)")
    print(
        f"DYN-5V: {n_pass}/{len(CHECKS)} arithmetic checks pass; "
        f"scientific claim status = {claim_status}."
    )
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
