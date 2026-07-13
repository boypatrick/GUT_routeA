#!/usr/bin/env python3
"""Composite hidden-GLSM audit for the crossed-120 unitary link.

No web lookup is used.

The previous local quotient audit imposed

    B^dagger B = I,        P B P = W_kappa,

and showed that all remaining link moduli are unitary and threshold-silent.
This script tests the next microscopic caricature: realize the unitary link as
a hidden-sector meson

    B = Qtilde Q / f^2

on a locked-radial D-flat branch of a hidden U(4)_H gauge theory.  The hidden
fields are visible Spin(10) singlets and carry only link/source copy indices,
so their direct one-loop SM/GUT threshold vector is zero in the preferred
branch.

This is not a full nonperturbative derivation.  The minimal Nf=Nc=4 SQCD-like
choice is a constrained composite-meson branch, while the Nf=8 option is the
cleaner conformal-window route.  The audit records which hidden flavor choices
avoid a UV Landau obstruction and whether the finite meson samples reproduce
the already-verified unitary completion card.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "composite_unitary_link_glsm"
PS_SOURCE = ROOT / "output" / "ps_crossed_120_source_action" / "summary.json"
COMPLETED = ROOT / "output" / "completed_120_partner_action" / "summary.json"
DTERM = ROOT / "output" / "unitary_link_dterm_quotient" / "summary.json"

N_C = 4
ALPHA_H_INV_GRID = [1.0, 5.0, 10.0, 25.0]
N_F_GRID = [4, 5, 6, 8, 12, 16]
R_TARGETS = [50.0, 200.0]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cval(raw: dict[str, float]) -> complex:
    return raw["re"] + 1j * raw["im"]


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cval(cell) for cell in row] for row in raw], dtype=complex)


def hermitian_sqrt(mat: np.ndarray, tol: float = 1.0e-13) -> np.ndarray:
    herm = 0.5 * (mat + mat.conj().T)
    vals, vecs = np.linalg.eigh(herm)
    if float(np.min(vals)) < -tol:
        raise RuntimeError(f"negative defect eigenvalue {np.min(vals)}")
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def random_unitary(n: int, angle: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = 0.5 * (raw + raw.conj().T)
    h = h / max(np.linalg.norm(h, ord=2), 1.0e-30)
    vals, vecs = np.linalg.eigh(h)
    return vecs @ np.diag(np.exp(1j * angle * vals)) @ vecs.conj().T


def unitary_completion(w: np.ndarray, angle: float, seed: int) -> np.ndarray:
    eye = np.eye(w.shape[0], dtype=complex)
    left = random_unitary(w.shape[0], angle, seed)
    right = random_unitary(w.shape[0], angle, seed + 7919)
    d_left = hermitian_sqrt(eye - w @ w.conj().T)
    d_right = hermitian_sqrt(eye - w.conj().T @ w)
    return np.block(
        [
            [w, d_left @ right],
            [left @ d_right, -left @ w.conj().T @ right],
        ]
    )


def sqcd_phase(nf: int) -> str:
    if nf < N_C:
        return "ADS_runaway_not_used"
    if nf == N_C:
        return "Nf=Nc_quantum_deformed_composite_meson"
    if nf < 1.5 * N_C:
        return "below_conformal_window_strong"
    if nf == int(1.5 * N_C):
        return "lower_edge_conformal_window"
    if nf < 3 * N_C:
        return "inside_conformal_window"
    if nf == 3 * N_C:
        return "one_loop_marginal_edge"
    return "IR_free_UV_Landau_possible"


def landau_ratio(alpha_inv: float, b_hidden: float) -> float | None:
    # Convention: d alpha^{-1}/d log mu = -b/(2 pi).
    # Positive b drives alpha^{-1} down in the UV and can hit a Landau pole.
    if b_hidden <= 0.0:
        return None
    return math.exp(2.0 * math.pi * alpha_inv / b_hidden)


def hidden_beta_scan() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for nf in N_F_GRID:
        b_hidden = float(nf - 3 * N_C)
        for alpha_inv in ALPHA_H_INV_GRID:
            ratio = landau_ratio(alpha_inv, b_hidden)
            for r_target in R_TARGETS:
                if ratio is None:
                    safe_to_r = True
                    ratio_print = None
                else:
                    safe_to_r = bool(ratio > r_target)
                    ratio_print = ratio
                rows.append(
                    {
                        "Nc": N_C,
                        "Nf": nf,
                        "b_hidden": b_hidden,
                        "alphaH_inv_at_matching": alpha_inv,
                        "phase_label": sqcd_phase(nf),
                        "uv_landau_ratio_over_matching": ratio_print,
                        "target_R": r_target,
                        "safe_to_target_R": safe_to_r,
                    }
                )
    return rows


def meson_sample_rows(w: np.ndarray, lock_window: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, angle in enumerate([0.0, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]):
        b = unitary_completion(w, angle, 202605090157 + idx)
        f = 1.0
        q = f * np.eye(b.shape[0], dtype=complex)
        qtilde = f * b
        meson = qtilde @ q / (f * f)
        d_left = q @ q.conj().T - qtilde.conj().T @ qtilde
        d_right = q.conj().T @ q - qtilde @ qtilde.conj().T
        singulars = np.linalg.svd(meson, compute_uv=False)
        max_abs_log = float(np.max(np.abs(np.log(singulars))))
        rows.append(
            {
                "angle": angle,
                "fixed_block_residual_fro": float(np.linalg.norm(meson[: w.shape[0], : w.shape[1]] - w, ord="fro")),
                "meson_reconstruction_residual_fro": float(np.linalg.norm(meson - b, ord="fro")),
                "left_D_residual_2norm": float(np.linalg.norm(d_left, ord=2)),
                "right_D_residual_2norm": float(np.linalg.norm(d_right, ord=2)),
                "unitarity_residual_2norm": float(np.linalg.norm(meson.conj().T @ meson - np.eye(meson.shape[0]), ord=2)),
                "mass_singular_min_over_Slock": float(np.min(singulars)),
                "mass_singular_max_over_Slock": float(np.max(singulars)),
                "max_abs_log_singular": max_abs_log,
                "within_completed_partner_lock_window": bool(max_abs_log <= lock_window),
            }
        )
    return rows


def build() -> dict[str, Any]:
    ps = read_json(PS_SOURCE)
    completed = read_json(COMPLETED)
    dterm = read_json(DTERM)
    w = cmat(ps["superpotential_ansatz"]["W_kappa"])
    lock_window = float(completed["mass_locking_window"]["max_abs_log_xi"])
    beta_rows = hidden_beta_scan()
    sample_rows = meson_sample_rows(w, lock_window)

    preferred_minimal = next(
        row
        for row in beta_rows
        if row["Nf"] == 4 and row["alphaH_inv_at_matching"] == 10.0 and row["target_R"] == 200.0
    )
    preferred_conformal = next(
        row
        for row in beta_rows
        if row["Nf"] == 8 and row["alphaH_inv_at_matching"] == 10.0 and row["target_R"] == 200.0
    )

    all_samples_ok = all(
        row["within_completed_partner_lock_window"]
        and row["fixed_block_residual_fro"] < 1.0e-12
        and row["left_D_residual_2norm"] < 1.0e-12
        and row["right_D_residual_2norm"] < 1.0e-12
        for row in sample_rows
    )
    max_d_res = max(max(row["left_D_residual_2norm"], row["right_D_residual_2norm"]) for row in sample_rows)
    max_unitarity = max(row["unitarity_residual_2norm"] for row in sample_rows)
    max_log = max(row["max_abs_log_singular"] for row in sample_rows)

    verdict = {
        "all_finite_meson_samples_dflat_and_locked": all_samples_ok,
        "max_D_residual_2norm": max_d_res,
        "max_unitarity_residual_2norm": max_unitarity,
        "max_abs_log_singular": max_log,
        "visible_threshold_vector": [0.0, 0.0, 0.0],
        "minimal_Nf4_uv_landau_safe_to_R200": bool(preferred_minimal["safe_to_target_R"]),
        "minimal_Nf4_phase": preferred_minimal["phase_label"],
        "conformal_Nf8_uv_landau_safe_to_R200": bool(preferred_conformal["safe_to_target_R"]),
        "conformal_Nf8_phase": preferred_conformal["phase_label"],
        "interpretation": (
            "A locked-radial hidden U(4)_H branch can realize the link as the "
            "meson B=Qtilde Q/f^2 with D-flat residuals at numerical zero and "
            "with PBP=W_kappa inherited from the unitary completion.  Because "
            "Q and Qtilde are visible Spin(10) singlets in this branch, the "
            "direct visible threshold vector is zero.  The minimal Nf=Nc=4 "
            "route is a constrained composite-meson assumption, not a weakly "
            "coupled UV proof; an Nf=8 hidden conformal-window route is a more "
            "plausible microscopic direction but adds hidden spectators."
        ),
    }

    return {
        "note": "No web lookup used. Hidden U(4)_H composite GLSM audit for the unitary link.",
        "meson_model": {
            "hidden_gauge_group": "U(4)_H, with the SU(4) beta estimate used for the nonabelian factor",
            "fields": "Q, Qtilde as 4x4 hidden fundamental/antifundamental matrices with visible copy indices",
            "locked_radial_branch": "Q=f I_4, Qtilde=f B, B in U(4)",
            "meson_coordinate": "B=Qtilde Q/f^2",
            "hidden_Dflat_condition": "Q Q^dagger - Qtilde^dagger Qtilde = 0 and Q^dagger Q - Qtilde Qtilde^dagger = 0 on the locked branch",
            "fixed_block_driver": "Tr Y(P B P-W_kappa)",
        },
        "dimension_count": {
            "raw_Q_Qtilde_real_components": 64,
            "hidden_Dterm_rank": 16,
            "hidden_gauge_orbit_rank": 16,
            "radial_lock_rank": 16,
            "unitary_link_real_dimension": dterm["Fterm_on_Dflat_tangent"]["unitary_tangent_real_dimension"],
            "fixed_block_rank_on_unitary_tangent": dterm["Fterm_on_Dflat_tangent"]["fixed_block_rank_on_Dflat_tangent"],
            "residual_unitary_moduli_real_dimension": dterm["verdict"]["Dflat_fixed_block_moduli_real_dimension"],
        },
        "benchmarks": {
            "MG_GeV": float(ps["benchmarks"]["MG_GeV"]),
            "M_lock_GeV": float(ps["benchmarks"]["M_lock_GeV"]),
            "completed_partner_max_abs_log_xi": lock_window,
        },
        "finite_meson_samples": sample_rows,
        "hidden_beta_scan": beta_rows,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(payload: dict[str, Any]) -> None:
    v = payload["verdict"]
    lines = [
        "# Composite unitary-link hidden GLSM audit",
        "",
        "No web lookup was used.",
        "",
        "## Meson branch",
        "",
        "The tested branch is",
        "",
        "```text",
        "Q = f I_4,    Qtilde = f B,    B in U(4),    B = Qtilde Q/f^2.",
        "```",
        "",
        "The hidden D-flatness residuals are computed as",
        "",
        "```text",
        "||Q Q^dagger - Qtilde^dagger Qtilde||_2,",
        "||Q^dagger Q - Qtilde Qtilde^dagger||_2.",
        "```",
        "",
        "## Finite meson samples",
        "",
        "| angle | fixed block | D residual | unitary residual | max abs log s | locked |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["finite_meson_samples"]:
        d_res = max(row["left_D_residual_2norm"], row["right_D_residual_2norm"])
        lines.append(
            f"| {row['angle']:.0e} | {row['fixed_block_residual_fro']:.3e} | "
            f"{d_res:.3e} | {row['unitarity_residual_2norm']:.3e} | "
            f"{row['max_abs_log_singular']:.3e} | {row['within_completed_partner_lock_window']} |"
        )
    lines.extend(
        [
            "",
            "## Hidden beta scan",
            "",
            "With the convention d alpha_H^{-1}/d log mu = -b_H/(2 pi),",
            "",
            "```text",
            "b_H = N_f - 3 N_c,    N_c=4.",
            "```",
            "",
            "The preferred minimal row Nf=4 is asymptotically free but is a",
            "quantum-deformed/constrained meson branch.  The Nf=8 row lies in",
            "the SQCD conformal-window range and is the cleaner microscopic",
            "direction if hidden spectators are allowed.",
            "",
            "## Verdict",
            "",
            f"max D residual: {v['max_D_residual_2norm']:.3e}",
            f"max unitary residual: {v['max_unitarity_residual_2norm']:.3e}",
            f"max |log singular|: {v['max_abs_log_singular']:.3e}",
            f"visible threshold vector: {v['visible_threshold_vector']}",
            "",
            v["interpretation"],
            "",
        ]
    )
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "finite_meson_samples.csv",
        payload["finite_meson_samples"],
        [
            "angle",
            "fixed_block_residual_fro",
            "meson_reconstruction_residual_fro",
            "left_D_residual_2norm",
            "right_D_residual_2norm",
            "unitarity_residual_2norm",
            "mass_singular_min_over_Slock",
            "mass_singular_max_over_Slock",
            "max_abs_log_singular",
            "within_completed_partner_lock_window",
        ],
    )
    write_csv(
        OUT / "hidden_beta_scan.csv",
        payload["hidden_beta_scan"],
        [
            "Nc",
            "Nf",
            "b_hidden",
            "alphaH_inv_at_matching",
            "phase_label",
            "uv_landau_ratio_over_matching",
            "target_R",
            "safe_to_target_R",
        ],
    )
    write_report(payload)
    v = payload["verdict"]
    print("Composite unitary-link hidden GLSM audit")
    print(f"  all finite meson samples D-flat and locked: {v['all_finite_meson_samples_dflat_and_locked']}")
    print(f"  max D residual: {v['max_D_residual_2norm']:.3e}")
    print(f"  max |log singular|: {v['max_abs_log_singular']:.3e}")
    print(f"  Nf=4 phase: {v['minimal_Nf4_phase']}")
    print(f"  Nf=8 phase: {v['conformal_Nf8_phase']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
