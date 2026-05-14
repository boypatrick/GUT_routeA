#!/usr/bin/env python3
"""Endpoint-source-copy Hessian and relative-alignment audit.

No web lookup is used.  The previous combined charge table forced the single
54 source S_54 to split into endpoint-labelled copies S_AA, S_AB, S_AC.  This
script puts those copies back into the component Hessian and checks two things:

1. If each endpoint source has a Spin(10)-scalar source-link-driver mass, the
   source-copy spectrum is a sum of complete 54 multiplets and gives no
   projected gauge threshold.
2. Without an extra relative-alignment driver, S_AB and S_AC may rotate away
   from S_AA inside the 54 orbit.  Such relative-orientation moduli deform the
   AB/AC F54 operators, generate off-block component mixing, and must either be
   scanned as a threshold parameter or set to zero by a compatible driver.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_dynamic_source_mixing as dyn  # noqa: E402
import audit_untruncated_spin10_invariants as inv  # noqa: E402
import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "endpoint_source_copies"


def sym_inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(a * b))


def normalized_traceless(mat: np.ndarray) -> np.ndarray:
    sym = 0.5 * (mat + mat.T)
    sym -= np.eye(sym.shape[0]) * float(np.trace(sym)) / sym.shape[0]
    norm = math.sqrt(max(sym_inner(sym, sym), 1.0e-300))
    return sym / norm


def aligned_s54_tensor() -> np.ndarray:
    return np.diag([-2.0 / 3.0] * 6 + [1.0] * 4)


def orthogonalize_to_s0(mat: np.ndarray, s0: np.ndarray) -> np.ndarray:
    out = normalized_traceless(mat)
    out = out - sym_inner(out, s0) / sym_inner(s0, s0) * s0
    return normalized_traceless(out)


def deformation_tensors() -> dict[str, np.ndarray]:
    s0 = aligned_s54_tensor()

    color_diag = np.zeros((10, 10), dtype=float)
    color_diag[0, 0] = 1.0
    color_diag[1, 1] = -1.0

    weak_diag = np.zeros((10, 10), dtype=float)
    weak_diag[6, 6] = 1.0
    weak_diag[7, 7] = -1.0

    color_weak = np.zeros((10, 10), dtype=float)
    color_weak[0, 6] = color_weak[6, 0] = 1.0

    mixed_two_plane = np.zeros((10, 10), dtype=float)
    mixed_two_plane[0, 6] = mixed_two_plane[6, 0] = 1.0
    mixed_two_plane[1, 7] = mixed_two_plane[7, 1] = -1.0

    return {
        "color_diag_split": orthogonalize_to_s0(color_diag, s0),
        "weak_diag_split": orthogonalize_to_s0(weak_diag, s0),
        "color_weak_offblock": orthogonalize_to_s0(color_weak, s0),
        "mixed_two_plane_offblock": orthogonalize_to_s0(mixed_two_plane, s0),
    }


def operator_from_54_tensor(s: np.ndarray) -> np.ndarray:
    """Map a symmetric traceless 54 tensor to its operator on 45 two-forms.

    The bilinear is B(A,B)=1/2 Tr(A^T S B+B^T S A).  For the aligned diagonal
    tensor diag(-2/3^6,1^4), this exactly reproduces inv.f54_operator().
    """

    mats = inv.BASIS_45
    out = np.zeros((45, 45), dtype=float)
    for i, a in enumerate(mats):
        for j, b in enumerate(mats):
            out[i, j] = 0.5 * (np.trace(a.T @ s @ b) + np.trace(b.T @ s @ a))
    return 0.5 * (out + out.T)


def operator_checks() -> dict[str, Any]:
    s0 = aligned_s54_tensor()
    op0 = operator_from_54_tensor(s0)
    ref = inv.f54_operator()
    projectors = inv.block_projectors()
    tensors = deformation_tensors()
    mode_checks: dict[str, Any] = {}
    for name, tensor in tensors.items():
        op = operator_from_54_tensor(tensor)
        mode_checks[name] = {
            "tensor_norm": math.sqrt(sym_inner(tensor, tensor)),
            "inner_with_aligned_s54": sym_inner(tensor, s0),
            "operator_frobenius": float(np.linalg.norm(op)),
            "broad_PS_offblock_max_abs": inv.off_block_norms(op, projectors)["max_abs"],
            "broad_PS_offblock_frobenius": inv.off_block_norms(op, projectors)["frobenius"],
        }
    return {
        "aligned_operator_max_abs_error": float(np.max(np.abs(op0 - ref))),
        "aligned_operator_frobenius_error": float(np.linalg.norm(op0 - ref)),
        "mode_checks": mode_checks,
    }


def endpoint_source_matrix(source_mass: float, ell: float, kappa: float = 1.0) -> np.ndarray:
    # Basis [S_e, L_e, D_e].
    return np.array(
        [
            [source_mass, 0.0, -kappa * ell],
            [0.0, 0.0, kappa],
            [-kappa * ell, kappa, 0.0],
        ],
        dtype=float,
    )


def endpoint_source_threshold() -> dict[str, Any]:
    pars = dyn.card_ells()
    edges = [
        ("AA", float(pars["ells_54"][0])),
        ("AB", float(pars["ells_54"][1])),
        ("AC", float(pars["ells_54"][2])),
    ]
    beta_total = sum((item["beta"] for item in dyn.FRAGMENTS_54.values()), np.zeros(3))
    rows: dict[str, Any] = {}
    total_delta = np.zeros(3, dtype=float)
    for edge, ell in edges:
        matrix = endpoint_source_matrix(1.0, ell)
        evals = np.linalg.eigvalsh(matrix)
        log_sum = float(np.sum(np.log(1.0 / np.maximum(np.abs(evals), 1.0e-14))))
        # Every eigenvalue is repeated over a complete 54, so only the complete
        # 54 beta vector appears.
        delta = beta_total * log_sum / (2.0 * math.pi)
        rows[edge] = {
            "ell": ell,
            "matrix_over_MG": matrix.tolist(),
            "eigenvalues_over_MG": evals.tolist(),
            "log_sum": log_sum,
            "delta": delta.tolist(),
            "projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ delta)),
        }
        total_delta += delta
    return {
        "edges": rows,
        "complete_54_beta_sum": beta_total.tolist(),
        "complete_54_projected_beta_l2": float(np.linalg.norm(base.PROJECTOR @ beta_total)),
        "total_delta": total_delta.tolist(),
        "total_projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ total_delta)),
        "theorem_pass": float(np.linalg.norm(base.PROJECTOR @ total_delta)) < 1.0e-12,
    }


def insert_block(big: np.ndarray, i: int, j: int, op: np.ndarray) -> None:
    si = 45 * i
    sj = 45 * j
    big[si : si + 45, sj : sj + 45] += op
    if i != j:
        big[sj : sj + 45, si : si + 45] += op.T


def hessian_with_endpoint_orientations(
    card: dict[str, Any],
    eps_ab: float,
    eps_ac: float,
    tensor_ab: np.ndarray,
    tensor_ac: np.ndarray,
) -> np.ndarray:
    ops = route_a.baseline_operators(card)
    pars = card["vev_parameters_and_couplings"]
    eta = float(pars["eta"])
    i45 = np.eye(45)
    f0 = inv.f54_operator()
    fab = f0 + eps_ab * operator_from_54_tensor(tensor_ab)
    fac = f0 + eps_ac * operator_from_54_tensor(tensor_ac)

    h = np.zeros((135, 135), dtype=float)
    # Keep the verified AA/color blocks exactly as in the R=200 component card;
    # only the endpoint-copy relative orientation in AB/AC is varied.
    h += np.array(ops["H0"], copy=True)
    h[0:45, 45:90] = 0.0
    h[45:90, 0:45] = 0.0
    h[0:45, 90:135] = 0.0
    h[90:135, 0:45] = 0.0
    insert_block(h, 0, 1, eta * (fab - 2.0 * i45))
    insert_block(h, 0, 2, eta * (fac + (4.0 / 3.0) * i45))
    return 0.5 * (h + h.T)


def evaluate_hessian(h: np.ndarray, ops: dict[str, Any], grid: dict[str, np.ndarray], prefactor: float) -> tuple[dict[str, Any], dict[str, Any]]:
    metrics = route_a.spectrum_metrics(h, ops)
    raw = np.array(metrics["heavy_threshold_delta"], dtype=float)
    calibrated = (
        np.array(ops["card_heavy_delta"], dtype=float)
        + raw
        - np.array(ops["baseline_raw_heavy_delta"], dtype=float)
    )
    metrics["heavy_threshold_delta"] = calibrated.tolist()
    metrics["heavy_threshold_projected_l2"] = float(np.linalg.norm(base.PROJECTOR @ calibrated))
    rge = route_a.fixed_spectrum_rge_scan(metrics, grid, prefactor)
    return metrics, rge


def scenario_rows() -> list[dict[str, Any]]:
    card = comp.load_card()
    ops = route_a.baseline_operators(card)
    baseline_raw = route_a.spectrum_metrics(ops["H0"], ops)
    ops["baseline_raw_heavy_delta"] = baseline_raw["heavy_threshold_delta"]
    ops["card_heavy_delta"] = card["mediator_heavy_threshold"]["delta_full"]
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    tensors = deformation_tensors()

    specs: list[tuple[str, str, str, float, float]] = [
        ("aligned_endpoint_copies", "color_weak_offblock", "color_weak_offblock", 0.0, 0.0),
        ("offblock_common_1e-12", "color_weak_offblock", "color_weak_offblock", 1.0e-12, 1.0e-12),
        ("offblock_common_1e-10", "color_weak_offblock", "color_weak_offblock", 1.0e-10, 1.0e-10),
        ("offblock_common_1e-8", "color_weak_offblock", "color_weak_offblock", 1.0e-8, 1.0e-8),
        ("offblock_common_1e-6", "color_weak_offblock", "color_weak_offblock", 1.0e-6, 1.0e-6),
        ("offblock_antisym_1e-6", "color_weak_offblock", "color_weak_offblock", 1.0e-6, -1.0e-6),
        ("diag_common_1e-4", "color_diag_split", "color_diag_split", 1.0e-4, 1.0e-4),
        ("diag_common_1e-3", "color_diag_split", "color_diag_split", 1.0e-3, 1.0e-3),
        ("mixed_two_plane_1e-6", "mixed_two_plane_offblock", "mixed_two_plane_offblock", 1.0e-6, 1.0e-6),
    ]

    rows: list[dict[str, Any]] = []
    base_k3 = route_a.REPS[0][2]
    base_k8 = route_a.REPS[3][2]
    for name, mode_ab, mode_ac, eps_ab, eps_ac in specs:
        h = hessian_with_endpoint_orientations(
            card, eps_ab, eps_ac, tensors[mode_ab], tensors[mode_ac]
        )
        metrics, rge = evaluate_hessian(h, ops, grid, prefactor)
        safe = (
            metrics["goldstone_max_abs_kappa"] < 1.0e-6
            and metrics["x622_unwanted_projected_l2"] < 5.022738709841473e-4
            and metrics["fine_rep_offblock_max_abs"] < 1.0e-10
            and rge["safe_count_residual_lt_1e_3"] > 0
        )
        rows.append(
            {
                "scenario": name,
                "mode_ab": mode_ab,
                "mode_ac": mode_ac,
                "epsilon_ab": eps_ab,
                "epsilon_ac": eps_ac,
                "fine_rep_offblock_max_abs": float(metrics["fine_rep_offblock_max_abs"]),
                "fine_rep_offblock_frobenius": float(metrics["fine_rep_offblock_frobenius"]),
                "kappa3_eff": float(metrics["kappa3_eff"]),
                "kappa8_eff": float(metrics["kappa8_eff"]),
                "sigma3_rel_shift": abs(float(metrics["kappa3_eff"]) - base_k3) / base_k3,
                "sigma8_rel_shift": abs(float(metrics["kappa8_eff"]) - base_k8) / base_k8,
                "x622_unwanted_projected_l2": float(metrics["x622_unwanted_projected_l2"]),
                "goldstone_max_abs_kappa": float(metrics["goldstone_max_abs_kappa"]),
                "heavy_threshold_projected_l2": float(metrics["heavy_threshold_projected_l2"]),
                "rge_safe_count": int(rge["safe_count_residual_lt_1e_3"]),
                "rge_best_residual_l2": float(rge["best_residual_l2"]),
                "best_alphaG_inv": float(rge["best_alphaG_inv"]),
                "tau_d6": float(rge["best_tau_dim6_years"]),
                "tau_d5": float(rge["best_tau_dim5_ST_1e_5_years"]),
                "safe_without_relative_driver": bool(safe),
            }
        )
    return rows


def relative_alignment_driver_card() -> dict[str, Any]:
    return {
        "superpotential": [
            "W_rel = <Y_AB^54, S_AB^54 - v_AB^54 S_AA^54>",
            "      + <Y_AC^54, S_AC^54 - v_AC^54 S_AA^54>",
        ],
        "charges": {
            "Y_AB^54": {"R": 2, "U1M": 1, "Z3E": 1},
            "Y_AC^54": {"R": 2, "U1M": -1, "Z3E": 1},
            "v_AB^54": {"R": 0, "U1M": -1, "Z3E": 2},
            "v_AC^54": {"R": 0, "U1M": 1, "Z3E": 2},
        },
        "allowed_terms": {
            "Y_AB S_AB": [2, 0, 0],
            "Y_AB v_AB S_AA": [2, 0, 0],
            "Y_AC S_AC": [2, 0, 0],
            "Y_AC v_AC S_AA": [2, 0, 0],
        },
        "forbidden_cross_terms": [
            "Y_AB S_AC",
            "Y_AC S_AB",
            "Y_AB S_AA without v_AB spurion",
            "Y_AC S_AA without v_AC spurion",
        ],
        "effect": "Sets the relative 54 orientation moduli to zero, so AB and AC use the same F54 operator as AA.",
    }


def summarize(rows: list[dict[str, Any]], endpoint: dict[str, Any]) -> dict[str, Any]:
    aligned = next(row for row in rows if row["scenario"] == "aligned_endpoint_copies")
    unsafe_due_to_offblock = [
        row
        for row in rows
        if row["fine_rep_offblock_max_abs"] >= 1.0e-10
        or not row["safe_without_relative_driver"]
    ]
    first_unsafe = min(
        [row for row in rows if row["scenario"] != "aligned_endpoint_copies"],
        key=lambda row: max(abs(row["epsilon_ab"]), abs(row["epsilon_ac"]))
        if not row["safe_without_relative_driver"]
        else float("inf"),
    )
    return {
        "endpoint_source_projected_threshold_l2": endpoint["total_projected_delta_l2"],
        "aligned_safe_without_relative_driver": aligned["safe_without_relative_driver"],
        "aligned_offblock_max_abs": aligned["fine_rep_offblock_max_abs"],
        "number_of_tested_misalignment_scenarios": len(rows) - 1,
        "number_unsafe_or_offblock": len(unsafe_due_to_offblock),
        "first_tested_unsafe": first_unsafe,
        "interpretation": (
            "Endpoint source copies are harmless as complete multiplets when aligned, "
            "but the charge table leaves relative 54 orientation moduli.  Generic "
            "off-block misalignment immediately breaks the PS block structure, so "
            "the minimal relative-alignment driver should be added."
        ),
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "scenario",
        "mode_ab",
        "mode_ac",
        "epsilon_ab",
        "epsilon_ac",
        "fine_rep_offblock_max_abs",
        "fine_rep_offblock_frobenius",
        "sigma3_rel_shift",
        "sigma8_rel_shift",
        "x622_unwanted_projected_l2",
        "goldstone_max_abs_kappa",
        "heavy_threshold_projected_l2",
        "rge_safe_count",
        "rge_best_residual_l2",
        "best_alphaG_inv",
        "tau_d6",
        "tau_d5",
        "safe_without_relative_driver",
    ]
    with (OUT / "endpoint_source_misalignment_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fields})


def write_report(payload: dict[str, Any]) -> None:
    rows = payload["misalignment_scan"]
    endpoint = payload["endpoint_source_threshold"]
    summary = payload["summary"]
    driver = payload["relative_alignment_driver"]
    lines: list[str] = []
    lines.append("# Endpoint source-copy Hessian audit")
    lines.append("")
    lines.append("No web lookup was used.  This pass puts `S_AA^54`, `S_AB^54`,")
    lines.append("and `S_AC^54` back into the component Hessian.")
    lines.append("")
    lines.append("## Complete-multiplet source-copy theorem")
    lines.append("")
    lines.append("For each endpoint `e`, use the local quadratic block")
    lines.append("")
    lines.append("```text")
    lines.append("W_e = 1/2 m_e ||S_e||^2 + kappa_e <D_e, L_e - ell_e S_e>.")
    lines.append("```")
    lines.append("")
    lines.append("If `m_e` and `kappa_e` are Spin(10)-scalar, this 3x3 matrix is")
    lines.append("identically repeated on every component of a full `54`; each eigenvalue")
    lines.append("is therefore a complete `54` multiplet.")
    lines.append("")
    lines.append("```text")
    lines.append(f"P beta_54 l2        = {endpoint['complete_54_projected_beta_l2']:.3e}")
    lines.append(f"P Delta_source l2   = {endpoint['total_projected_delta_l2']:.3e}")
    lines.append(f"theorem pass        = {endpoint['theorem_pass']}")
    lines.append("```")
    lines.append("")
    lines.append("## Relative-orientation modulus")
    lines.append("")
    lines.append("The charge table aligns each link to its own source copy, but does not")
    lines.append("force `S_AB^54` and `S_AC^54` to point in the same `54` direction as")
    lines.append("`S_AA^54`.  A relative orientation")
    lines.append("")
    lines.append("```text")
    lines.append("S_AB = S0 + epsilon_AB T_AB,")
    lines.append("S_AC = S0 + epsilon_AC T_AC,")
    lines.append("<T,S0>=0")
    lines.append("```")
    lines.append("")
    lines.append("changes the AB/AC operators from `F54` to `F54+epsilon F_T`.")
    lines.append("")
    lines.append("## Misalignment scan")
    lines.append("")
    lines.append("| scenario | eps_AB | eps_AC | offblock max | safe points | residual | safe without driver |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| `{row['scenario']}` | {row['epsilon_ab']:.1e} | {row['epsilon_ac']:.1e} | "
            f"{row['fine_rep_offblock_max_abs']:.3e} | {row['rge_safe_count']} | "
            f"{row['rge_best_residual_l2']:.3e} | {row['safe_without_relative_driver']} |"
        )
    lines.append("")
    lines.append("The RGE residual can remain numerically small for tiny misalignment, but")
    lines.append("nonzero off-block entries mean the Pati--Salam fragment threshold basis is")
    lines.append("no longer the one used in the paper.  This is a structural modulus, not")
    lines.append("just a harmless small correction.")
    lines.append("")
    lines.append("## Minimal relative-alignment driver")
    lines.append("")
    lines.append("Add")
    lines.append("")
    lines.append("```text")
    for row in driver["superpotential"]:
        lines.append(row)
    lines.append("```")
    lines.append("")
    lines.append("with charges")
    lines.append("")
    lines.append("```text")
    for name, charge in driver["charges"].items():
        lines.append(f"{name}: R={charge['R']}, U1M={charge['U1M']}, Z3E={charge['Z3E']}")
    lines.append("```")
    lines.append("")
    lines.append("This sets the relative orientation moduli to zero.  With the driver, the")
    lines.append("aligned endpoint-copy row is the physical result:")
    lines.append("")
    lines.append("```text")
    lines.append(f"aligned offblock max = {summary['aligned_offblock_max_abs']:.3e}")
    lines.append(f"aligned safe         = {summary['aligned_safe_without_relative_driver']}")
    lines.append(f"source P Delta l2    = {summary['endpoint_source_projected_threshold_l2']:.3e}")
    lines.append("```")
    (OUT / "endpoint_source_copy_hessian_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    op = operator_checks()
    endpoint = endpoint_source_threshold()
    rows = scenario_rows()
    payload = {
        "note": "No web lookup used. Endpoint source-copy Hessian and relative-orientation audit.",
        "operator_checks": op,
        "endpoint_source_threshold": endpoint,
        "misalignment_scan": rows,
        "relative_alignment_driver": relative_alignment_driver_card(),
        "summary": summarize(rows, endpoint),
        "passes": {
            "aligned_54_operator_reproduces_F54": op["aligned_operator_max_abs_error"] < 1.0e-12,
            "endpoint_source_copies_complete_multiplets": endpoint["theorem_pass"],
            "aligned_endpoint_hessian_safe": next(
                row for row in rows if row["scenario"] == "aligned_endpoint_copies"
            )["safe_without_relative_driver"],
            "relative_driver_required": any(
                row["fine_rep_offblock_max_abs"] >= 1.0e-10
                for row in rows
                if row["scenario"] != "aligned_endpoint_copies"
            ),
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "endpoint_source_copy_hessian_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(payload)

    aligned = next(row for row in rows if row["scenario"] == "aligned_endpoint_copies")
    first = payload["summary"]["first_tested_unsafe"]
    print("Endpoint source-copy Hessian audit")
    print(f"  aligned F54 operator error: {op['aligned_operator_max_abs_error']:.3e}")
    print(f"  source-copy P Delta l2: {endpoint['total_projected_delta_l2']:.3e}")
    print(f"  aligned offblock max: {aligned['fine_rep_offblock_max_abs']:.3e}")
    print(f"  aligned safe points: {aligned['rge_safe_count']}")
    print(
        "  first tested unsafe/misaligned: "
        f"{first['scenario']} offblock={first['fine_rep_offblock_max_abs']:.3e}"
    )
    print(f"  relative driver required: {payload['passes']['relative_driver_required']}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
