#!/usr/bin/env python3
"""Construct a post-breaking PS source action for the crossed 120 projector.

No web lookup is used.

The crossed 120A/120B projector cannot be enforced by field-only unbroken
Spin(10) charges.  This script treats it honestly as a post-Spin(10)-breaking
Pati-Salam component source action and checks the minimal algebraic claims:

1. The triplet block realizes the desired visible inverse propagator

       W_T = |120_A><120_B| + epsilon |120_B><120_A|

   as a Julia/Sz.-Nagy unitary dilation.
2. The electroweak doublet block is left diagonal and unmixed.
3. The zero-field configuration is F-flat and D-flat for the quadratic
   source action.
4. If the triplet source is interpreted as a constrained/nonpropagating
   source, it is threshold silent.  If it is instead interpreted as literal
   propagating triplet-only PS mediator matter, the induced nonuniversal
   threshold is quantified and compared to the R=200 budget.

This is therefore a component-action audit, not a UV Spin(10) derivation.
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

from audit_gtr_mediator_spectrum import B_DOUBLE_PAIR, B_FIVE_PAIR, B_TRIPLET_PAIR, projected_l2  # noqa: E402


OUT = ROOT / "output" / "ps_crossed_120_source_action"
CROSSED = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"
LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"
DOWNSTREAM = ROOT / "output" / "tightened_downstream_proton" / "summary.json"

KAPPA = 100.0
OFFBLOCK_RHOS = [0.0, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def load_mg() -> float:
    payload = read_json(DOWNSTREAM)
    row = next(item for item in payload["benchmarks"] if item["benchmark"] == "R_window_R200")
    return float(row["MG_GeV"])


def target_w(kappa: float) -> np.ndarray:
    eps = 1.0 / kappa
    return np.array([[0.0, 1.0], [eps, 0.0]], dtype=complex)


def full_dimensionless_hessian(rho: float, w: np.ndarray) -> np.ndarray:
    """Return dimensionless mass Hessian in triplet+doublet blocks.

    The first four entries are the triplet visible+mediator fields of the
    Julia dilation; entries 4,5 are the 120_A/120_B doublet sources.  The
    off-block matrix is a deterministic unit-norm probe of PS-breaking leakage.
    """

    triplet = julia_dilation(w).conj().T
    doublet = np.eye(2, dtype=complex)
    bridge = np.zeros((4, 2), dtype=complex)
    bridge[0, 0] = 1.0 / math.sqrt(2.0)
    bridge[1, 1] = 1.0 / math.sqrt(2.0)
    return np.block(
        [
            [triplet, rho * bridge],
            [rho * bridge.conj().T, doublet],
        ]
    )


def zero_field_flatness(mat: np.ndarray) -> dict[str, float]:
    fields = np.zeros(mat.shape[0], dtype=complex)
    f_terms = mat @ fields
    return {
        "F_norm_at_zero": float(np.linalg.norm(f_terms)),
        "D_norm_at_zero_proxy": 0.0,
        "tadpole_norm": 0.0,
    }


def offblock_rows(w: np.ndarray, leak_row: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target_inv_t = w
    target_inv_d = np.eye(2, dtype=complex)
    llll0 = float(leak_row["LLLL_margin_1e35_replayed"])
    rrrr0 = float(leak_row["RRRR_margin_1e35_replayed"])
    llll_leak0 = float(leak_row["LLLL_leakage_ratio"])
    rrrr_leak0 = float(leak_row["RRRR_leakage_ratio"])
    for rho in OFFBLOCK_RHOS:
        mat = full_dimensionless_hessian(rho, w)
        inv = np.linalg.inv(mat)
        triplet_inv = inv[:2, :2]
        doublet_inv = inv[4:, 4:]
        triplet_residual = float(np.linalg.norm(triplet_inv - target_inv_t, ord="fro"))
        doublet_residual = float(np.linalg.norm(doublet_inv - target_inv_d, ord="fro"))
        offblock_inverse_norm = float(np.linalg.norm(inv[:4, 4:], ord="fro"))
        spectral = np.linalg.svd(mat, compute_uv=False)

        # A conservative proxy: treat triplet subblock deformation as an
        # additional linear leakage on top of the finite-lift leakage.
        llll_leak = llll_leak0 + triplet_residual
        rrrr_leak = rrrr_leak0 + triplet_residual
        llll_margin = llll0 * (llll_leak0 / max(llll_leak, 1.0e-300)) ** 2
        rrrr_margin = rrrr0 * (rrrr_leak0 / max(rrrr_leak, 1.0e-300)) ** 2
        rows.append(
            {
                "rho": float(rho),
                "triplet_inverse_residual_fro": triplet_residual,
                "doublet_inverse_residual_fro": doublet_residual,
                "offblock_inverse_norm_fro": offblock_inverse_norm,
                "mass_singular_min_over_Mlock": float(np.min(spectral)),
                "mass_singular_max_over_Mlock": float(np.max(spectral)),
                "mass_singular_spread_over_Mlock": float(np.max(spectral) - np.min(spectral)),
                "proxy_LLLL_margin_1e35": float(llll_margin),
                "proxy_RRRR_margin_1e35": float(rrrr_margin),
                "proxy_worst_margin_1e35": float(min(llll_margin, rrrr_margin)),
                "doublet_preserved_to_1e_minus_6": bool(doublet_residual < 1.0e-6),
                "proton_proxy_future_safe": bool(min(llll_margin, rrrr_margin) > 1.0),
            }
        )
    return rows


def threshold_rows(mg: float, m_lock: float, tol: float) -> list[dict[str, Any]]:
    kappa_mass = m_lock / mg
    log_factor = math.log(mg / m_lock) / (2.0 * math.pi)
    branches = [
        {
            "branch": "constrained_nonpropagating_source",
            "beta_vector": np.array([0.0, 0.0, 0.0], dtype=float),
            "interpretation": "Source is auxiliary/constrained; no propagating incomplete PS matter.",
        },
        {
            "branch": "two_literal_triplet_mediator_pairs",
            "beta_vector": 2.0 * B_TRIPLET_PAIR,
            "interpretation": "Julia defect fields are literal triplet pairs without doublet partners.",
        },
        {
            "branch": "two_complete_5_plus_5bar_mediator_pairs",
            "beta_vector": 2.0 * B_FIVE_PAIR,
            "interpretation": "Triplet mediators are completed by inert doublet partners at the same mass.",
        },
        {
            "branch": "two_missing_doublet_partners_only",
            "beta_vector": 2.0 * B_DOUBLE_PAIR,
            "interpretation": "Reference size of the inert doublet partners needed to complete the pair.",
        },
    ]
    rows = []
    for item in branches:
        beta = item["beta_vector"]
        delta = beta * log_factor
        p_l2 = projected_l2(delta)
        rows.append(
            {
                "branch": item["branch"],
                "kappa_Mlock_over_MG": float(kappa_mass),
                "log_MG_over_Mlock": float(math.log(mg / m_lock)),
                "beta_vector": [float(x) for x in beta],
                "delta_vector": [float(x) for x in delta],
                "projected_l2": float(p_l2),
                "within_R200_budget": bool(p_l2 <= tol),
                "interpretation": item["interpretation"],
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def build() -> dict[str, Any]:
    crossed = read_json(CROSSED)
    ledger = read_json(LEDGER)
    mg = load_mg()
    m_lock = float(ledger["r200_benchmark"]["M_Sigma8_GeV"])
    tol = float(ledger["r200_benchmark"]["total_projected_l2"])
    leak_row = next(row for row in crossed["finite_leakage_replay"] if row["kappa"] == KAPPA)

    w = target_w(KAPPA)
    u = julia_dilation(w)
    triplet_mass = u.conj().T
    flat = zero_field_flatness(triplet_mass)
    triplet_inverse_residual = float(np.linalg.norm(np.linalg.inv(triplet_mass)[:2, :2] - w, ord="fro"))
    triplet_unitary_residual = float(np.linalg.norm(u.conj().T @ u - np.eye(4), ord=2))
    singulars = np.linalg.svd(triplet_mass, compute_uv=False)
    off_rows = offblock_rows(w, leak_row)
    th_rows = threshold_rows(mg, m_lock, tol)
    exact_row = off_rows[0]
    constrained_threshold = next(row for row in th_rows if row["branch"] == "constrained_nonpropagating_source")
    literal_triplet_threshold = next(row for row in th_rows if row["branch"] == "two_literal_triplet_mediator_pairs")
    complete_threshold = next(row for row in th_rows if row["branch"] == "two_complete_5_plus_5bar_mediator_pairs")

    verdict = {
        "ps_component_action_realizes_projector": bool(
            flat["F_norm_at_zero"] == 0.0
            and triplet_inverse_residual < 1.0e-12
            and triplet_unitary_residual < 1.0e-12
            and exact_row["doublet_inverse_residual_fro"] < 1.0e-12
            and exact_row["proxy_worst_margin_1e35"] > 1.0
        ),
        "constrained_source_threshold_silent": bool(constrained_threshold["within_R200_budget"]),
        "literal_triplet_only_threshold_safe": bool(literal_triplet_threshold["within_R200_budget"]),
        "complete_5_pair_threshold_safe": bool(complete_threshold["within_R200_budget"]),
        "rho_max_doublet_preserved_1e_minus_6": max(
            row["rho"] for row in off_rows if row["doublet_preserved_to_1e_minus_6"]
        ),
        "rho_max_proton_proxy_future_safe": max(row["rho"] for row in off_rows if row["proton_proxy_future_safe"]),
        "literal_triplet_projected_l2_over_budget": float(
            literal_triplet_threshold["projected_l2"] / max(tol, 1.0e-300)
        ),
        "interpretation": (
            "A post-breaking PS component source action can realize the crossed 120 projector "
            "with zero F-terms at the origin, exact doublet decoupling, and a threshold-silent "
            "constrained-source interpretation.  The same algebra is not automatically a "
            "healthy literal propagating PS triplet-mediator completion: two triplet-only "
            "mediator pairs at M_lock overshoot the R=200 projected-threshold budget unless "
            "they are completed by inert doublet partners or treated as constrained/auxiliary sources."
        ),
    }

    return {
        "note": "No web lookup used. Post-breaking PS crossed 120 source action audit.",
        "ps_fragment_labels": {
            "120_decomposition_under_PS": [
                "(1,2,2)",
                "(15,2,2)",
                "(6,1,1)",
                "(10,1,1)",
                "(10bar,1,1)",
            ],
            "doublet_fragments": ["(1,2,2)", "(15,2,2)"],
            "triplet_fragments": ["color-triplet components inside (6,1,1)+(10,1,1)+(10bar,1,1)"],
            "source_copies": ["120_A", "120_B"],
        },
        "superpotential_ansatz": {
            "triplet": "W_T = M_lock * bar(Phi_T)^T U_T(W_kappa)^dagger Phi_T",
            "doublet": "W_D = M_lock * bar(D)^T I_2 D",
            "offblock": "W_mix = rho M_lock * bar(Phi_T)^T C D + h.c. for leakage audit only",
            "W_kappa": matrix_json(w),
            "Julia_identity": "U=[[W,sqrt(I-WWdag)],[sqrt(I-WdagW),-Wdag]]; P(Udag)^(-1)P=W.",
        },
        "benchmarks": {
            "kappa": KAPPA,
            "M_lock_GeV": m_lock,
            "MG_GeV": mg,
            "R200_projected_threshold_budget": tol,
        },
        "flatness_and_hessian": {
            **flat,
            "triplet_inverse_subblock_residual_fro": triplet_inverse_residual,
            "triplet_unitary_residual_2norm": triplet_unitary_residual,
            "triplet_mass_singular_values_over_Mlock": [float(x) for x in singulars],
            "triplet_mass_singular_spread_over_Mlock": float(np.max(singulars) - np.min(singulars)),
        },
        "offblock_mixing_scan": off_rows,
        "threshold_interpretation_scan": th_rows,
        "finite_leakage_input": leak_row,
        "verdict": verdict,
    }


def write_report(payload: dict[str, Any]) -> None:
    rows = payload["offblock_mixing_scan"]
    th = payload["threshold_interpretation_scan"]
    verdict = payload["verdict"]
    lines = [
        "# PS crossed 120 source action audit",
        "",
        "No web lookup was used.",
        "",
        "## Source action",
        "",
        "```text",
        payload["superpotential_ansatz"]["triplet"],
        payload["superpotential_ansatz"]["doublet"],
        "```",
        "",
        "The triplet block uses the Julia/Sz.-Nagy dilation, while the doublet block is diagonal.",
        "",
        "## Flatness and Hessian",
        "",
        f"F-norm at zero: {payload['flatness_and_hessian']['F_norm_at_zero']:.1e}",
        f"Triplet inverse residual: {payload['flatness_and_hessian']['triplet_inverse_subblock_residual_fro']:.3e}",
        f"Triplet unitary residual: {payload['flatness_and_hessian']['triplet_unitary_residual_2norm']:.3e}",
        "",
        "## Off-block mixing scan",
        "",
        "| rho | triplet residual | doublet residual | inverse offblock | worst 1e35 margin | doublet preserved | proton proxy safe |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['rho']:.0e} | {row['triplet_inverse_residual_fro']:.3e} | "
            f"{row['doublet_inverse_residual_fro']:.3e} | {row['offblock_inverse_norm_fro']:.3e} | "
            f"{row['proxy_worst_margin_1e35']:.3e} | {row['doublet_preserved_to_1e_minus_6']} | "
            f"{row['proton_proxy_future_safe']} |"
        )
    lines.extend(
        [
            "",
            "## Threshold interpretations",
            "",
            "| branch | projected l2 | within R=200 budget | beta vector |",
            "|---|---:|---:|---|",
        ]
    )
    for row in th:
        lines.append(
            f"| {row['branch']} | {row['projected_l2']:.6e} | "
            f"{row['within_R200_budget']} | `{row['beta_vector']}` |"
        )
    lines.extend(["", "## Verdict", "", verdict["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "offblock_mixing_scan.csv",
        payload["offblock_mixing_scan"],
        [
            "rho",
            "triplet_inverse_residual_fro",
            "doublet_inverse_residual_fro",
            "offblock_inverse_norm_fro",
            "mass_singular_min_over_Mlock",
            "mass_singular_max_over_Mlock",
            "mass_singular_spread_over_Mlock",
            "proxy_LLLL_margin_1e35",
            "proxy_RRRR_margin_1e35",
            "proxy_worst_margin_1e35",
            "doublet_preserved_to_1e_minus_6",
            "proton_proxy_future_safe",
        ],
    )
    write_csv(
        OUT / "threshold_interpretation_scan.csv",
        payload["threshold_interpretation_scan"],
        [
            "branch",
            "kappa_Mlock_over_MG",
            "log_MG_over_Mlock",
            "projected_l2",
            "within_R200_budget",
            "beta_vector",
            "delta_vector",
            "interpretation",
        ],
    )
    write_report(payload)

    verdict = payload["verdict"]
    print("PS crossed 120 source action audit")
    print(f"  component action realizes projector: {verdict['ps_component_action_realizes_projector']}")
    print(f"  constrained source threshold silent: {verdict['constrained_source_threshold_silent']}")
    print(f"  literal triplet-only threshold safe: {verdict['literal_triplet_only_threshold_safe']}")
    print(f"  complete 5+5bar threshold safe: {verdict['complete_5_pair_threshold_safe']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
