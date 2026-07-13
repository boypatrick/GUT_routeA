#!/usr/bin/env python3
"""Derive and audit a crossed-120 link-locking component action.

No web lookup is used.

The previous completed-partner audit showed that a propagating crossed-120
repair is threshold safe only if the inert doublet partners are mass-locked to
the crossed triplet package at the per-mille level.  This script tests whether
a single holomorphic link/spurion action can enforce both requirements:

    visible inverse triplet block = W_kappa,
    inert doublet mass             = M_lock.

The answer is intentionally nuanced.

* The holomorphic constraints

      A B = 1,       P B P = W_kappa

  do derive the crossed inverse block.
* A common spurion S_lock derives a common overall mass for the triplet link
  and inert doublet partners.
* However, holomorphic F-terms alone do not force A to be unitary or
  degenerate.  The constraint space has residual complex moduli; generic exact
  holomorphic deformations keep A B=1 and PBP=W_kappa but split the singular
  values.  A D-term/Kahler/NLSM unitary constraint is still required for the
  threshold-silent branch.

Thus this is a conditional component-action derivation, not an unconditional
microscopic Spin(10) proof.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "crossed_120_link_locking_action"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"

RNG_SEED = 202605090048
DEFORMATION_AMPLITUDES = [0.0, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 3.0e-3, 1.0e-2]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cval(raw: dict[str, float]) -> complex:
    return raw["re"] + 1j * raw["im"]


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cval(cell) for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    if float(np.min(vals)) < -tol:
        raise RuntimeError(f"negative defect eigenvalue {np.min(vals)}")
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def julia_dilation(w: np.ndarray) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    return np.block(
        [
            [w, hermitian_sqrt(eye - w @ w.conj().T)],
            [hermitian_sqrt(eye - w.conj().T @ w), -w.conj().T],
        ]
    )


def constraint_vector(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> np.ndarray:
    c1 = a @ b - np.eye(a.shape[0], dtype=complex)
    c2 = b[: w.shape[0], : w.shape[1]] - w
    return np.concatenate([c1.reshape(-1), c2.reshape(-1)])


def real_pack(values: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(values), np.imag(values)])


def constraint_jacobian(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    cols = []
    for block in ["A", "B"]:
        for i in range(n):
            for j in range(n):
                for phase in [1.0, 1.0j]:
                    da = np.zeros_like(a)
                    db = np.zeros_like(b)
                    if block == "A":
                        da[i, j] = phase
                    else:
                        db[i, j] = phase
                    dc1 = da @ b + a @ db
                    dc2 = db[: w.shape[0], : w.shape[1]]
                    cols.append(real_pack(np.concatenate([dc1.reshape(-1), dc2.reshape(-1)])))
    return np.stack(cols, axis=1)


def random_delta_b(n: int, visible: int) -> np.ndarray:
    rng = np.random.default_rng(RNG_SEED)
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    raw[:visible, :visible] = 0.0
    return raw / np.linalg.norm(raw, ord="fro")


def deformation_scan(a0: np.ndarray, b0: np.ndarray, w: np.ndarray, max_log_window: float) -> list[dict[str, Any]]:
    n = a0.shape[0]
    db = random_delta_b(n, w.shape[0])
    rows = []
    for amp in DEFORMATION_AMPLITUDES:
        b = b0 + amp * db
        a = np.linalg.inv(b)
        residual = float(np.linalg.norm(constraint_vector(a, b, w)))
        unitary_residual = float(np.linalg.norm(a.conj().T @ a - np.eye(n), ord=2))
        singulars = np.linalg.svd(a, compute_uv=False)
        logs = np.log(singulars)
        max_abs_log = float(np.max(np.abs(logs)))
        rows.append(
            {
                "deformation_amplitude": float(amp),
                "holomorphic_constraint_residual": residual,
                "unitary_residual_2norm": unitary_residual,
                "mass_singular_min_over_Slock": float(np.min(singulars)),
                "mass_singular_max_over_Slock": float(np.max(singulars)),
                "max_abs_log_singular": max_abs_log,
                "within_completed_partner_lock_window": bool(max_abs_log <= max_log_window),
            }
        )
    return rows


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    completed = read_json(COMPLETED)
    w = cmat(ps["superpotential_ansatz"]["W_kappa"])
    u = julia_dilation(w)

    # Triplet mass is A0, and B0=A0^{-1}=U gives P B0 P = W.
    a0 = u.conj().T
    b0 = u
    max_log_window = float(completed["mass_locking_window"]["max_abs_log_xi"])

    constraint_residual = float(np.linalg.norm(constraint_vector(a0, b0, w)))
    unitary_residual = float(np.linalg.norm(a0.conj().T @ a0 - np.eye(a0.shape[0]), ord=2))
    jac = constraint_jacobian(a0, b0, w)
    singulars_j = np.linalg.svd(jac, compute_uv=False)
    rank = int(np.linalg.matrix_rank(jac, tol=1.0e-10))
    rows = deformation_scan(a0, b0, w, max_log_window)
    pure_holomorphic_fail_rows = [
        row
        for row in rows
        if row["holomorphic_constraint_residual"] < 1.0e-10
        and not row["within_completed_partner_lock_window"]
    ]

    spurion_rules = [
        {"operator": "S_lock * bar(Phi_T) A Phi_T", "status": "allowed"},
        {"operator": "S_lock * L_inert Lbar_inert", "status": "allowed"},
        {"operator": "independent S_T or S_D mass spurion", "status": "absent by single-spurion assumption"},
        {"operator": "16_i 16_j L_inert", "status": "forbidden by Z2_inert"},
        {"operator": "H_phys Lbar_inert", "status": "forbidden by Z2_inert"},
    ]
    verdict = {
        "holomorphic_constraints_realize_inverse_block": constraint_residual < 1.0e-12,
        "common_spurion_locks_triplet_and_doublet_scale": True,
        "nlsm_unitarity_required": len(pure_holomorphic_fail_rows) > 0,
        "pure_holomorphic_constraints_sufficient_for_threshold_lock": len(pure_holomorphic_fail_rows) == 0,
        "constraint_rank_real": rank,
        "constraint_nullity_real": int(jac.shape[1] - rank),
        "expected_constraint_rank_real": 40,
        "expected_nullity_real": 24,
        "smallest_deformation_failing_lock_window": None
        if not pure_holomorphic_fail_rows
        else pure_holomorphic_fail_rows[0]["deformation_amplitude"],
        "interpretation": (
            "The link action derives the crossed inverse block and a common S_lock can tie the "
            "triplet and inert-doublet mass scales together.  But the holomorphic equations "
            "AB=1 and PBP=W leave 24 real moduli; generic exact deformations along those "
            "moduli split the triplet singular values beyond the completed-partner mass-locking "
            "window.  Therefore the conditional branch still needs an explicit NLSM/D-term or "
            "composite unitarity constraint on the link."
        ),
    }
    return {
        "note": "No web lookup used. Crossed 120 link-locking action audit.",
        "superpotential_ansatz": {
            "formula": (
                "W = Tr X(AB-I) + Tr Y(PBP-W_kappa) + Z(S_lock-M_lock) "
                "+ S_lock[bar(Phi_T)^T A Phi_T + L_inert Lbar_inert]"
            ),
            "A0": matrix_json(a0),
            "B0": matrix_json(b0),
            "W_kappa": matrix_json(w),
            "unitarity_status": "A0 is unitary, but unitarity is not implied by the holomorphic F-terms.",
        },
        "benchmarks": {
            "M_lock_GeV": float(ps["benchmarks"]["M_lock_GeV"]),
            "MG_GeV": float(ps["benchmarks"]["MG_GeV"]),
            "completed_partner_max_abs_log_xi": max_log_window,
            "completed_partner_xi_min": float(completed["mass_locking_window"]["xi_min"]),
            "completed_partner_xi_max": float(completed["mass_locking_window"]["xi_max"]),
        },
        "flatness_at_reference": {
            "constraint_residual": constraint_residual,
            "unitarity_residual_2norm": unitary_residual,
            "driver_fields_zero_make_nonconstraint_F_terms_zero": True,
            "zero_field_matter_F_norm": 0.0,
        },
        "linearized_constraint_geometry": {
            "real_variables": int(jac.shape[1]),
            "real_equations": int(jac.shape[0]),
            "real_rank": rank,
            "real_nullity": int(jac.shape[1] - rank),
            "singular_values_min": float(np.min(singulars_j)),
            "singular_values_max": float(np.max(singulars_j)),
        },
        "holomorphic_moduli_deformation_scan": rows,
        "selection_rules": spurion_rules,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    lines = [
        "# Crossed 120 link-locking action audit",
        "",
        "No web lookup was used.",
        "",
        "## Superpotential ansatz",
        "",
        "```text",
        payload["superpotential_ansatz"]["formula"],
        "```",
        "",
        "## Constraint geometry",
        "",
        f"real rank: {payload['linearized_constraint_geometry']['real_rank']}",
        f"real nullity: {payload['linearized_constraint_geometry']['real_nullity']}",
        "",
        "## Holomorphic moduli deformation scan",
        "",
        "| amp | constraint residual | unitary residual | max |log s| | inside lock window |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in payload["holomorphic_moduli_deformation_scan"]:
        lines.append(
            f"| {row['deformation_amplitude']:.0e} | {row['holomorphic_constraint_residual']:.3e} | "
            f"{row['unitary_residual_2norm']:.3e} | {row['max_abs_log_singular']:.3e} | "
            f"{row['within_completed_partner_lock_window']} |"
        )
    lines.extend(["", "## Verdict", "", payload["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "holomorphic_moduli_deformation_scan.csv",
        payload["holomorphic_moduli_deformation_scan"],
        [
            "deformation_amplitude",
            "holomorphic_constraint_residual",
            "unitary_residual_2norm",
            "mass_singular_min_over_Slock",
            "mass_singular_max_over_Slock",
            "max_abs_log_singular",
            "within_completed_partner_lock_window",
        ],
    )
    write_report(payload)
    verdict = payload["verdict"]
    print("Crossed 120 link-locking action audit")
    print(f"  inverse block derived: {verdict['holomorphic_constraints_realize_inverse_block']}")
    print(f"  common spurion locks scale: {verdict['common_spurion_locks_triplet_and_doublet_scale']}")
    print(f"  pure holomorphic sufficient: {verdict['pure_holomorphic_constraints_sufficient_for_threshold_lock']}")
    print(f"  NLSM/D-term unitarity required: {verdict['nlsm_unitarity_required']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
