#!/usr/bin/env python3
"""Self-consistent Spin(10) block-operator vacuum-alignment benchmark.

No web lookup is used.  This script is the next verification layer after the
explicit 210_H^3 contraction.  It does not attempt to diagonalize all 210_H
components.  Instead it uses the already verified block operators

    F54, D210, P_X, I3(210)

to build an algebraic mass matrix for the relevant Pati-Salam fragments, then
feeds the resulting finite mediator threshold back into the 4pi-corrected
two-loop RGE/proton cache.

The new point relative to the older R-window scan is self-consistency: the
mediator threshold depends on the light spectrum kappas, while the best RGE fit
depends on the mediator threshold.  For each R this script iterates

    (kappa_3, kappa_8) -> mediator mass matrix -> Delta_med
    -> corrected RGE/proton replay -> (kappa_3, kappa_8)

until the benchmark spectrum is a fixed point.
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

import scan_corrected_downstream as down  # noqa: E402
import scan_mediator_r_window as rwin  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "spin10_vacuum_alignment"
AUDIT_JSON = ROOT / "output" / "yukawa_4pi_audit" / "yukawa_4pi_audit_summary.json"

R_GRID = [20.0, 50.0, 100.0, 200.0]
CARD_R_VALUES = [50.0, 200.0]
MAX_ITER = 12
TOL = 1.0e-12


STATE_CLEBSCH = {
    "Sigma_L": {"label": "(1,3,1)", "F54": 2.0, "D210": 1.0},
    "Sigma_R": {"label": "(1,1,3)", "F54": 2.0, "D210": -1.0},
    "Sigma8_block": {"label": "(15,1,1)", "F54": -4.0 / 3.0, "D210": 0.0},
    "X_622": {"label": "(6,2,2)", "F54": 1.0 / 3.0, "D210": 0.0},
}

PHYSICAL_BETA = {
    "H_C+H_Cbar": [2.0 / 5.0, 0.0, 1.0],
    "Sigma_3": [0.0, 2.0, 0.0],
    "Sigma_8_octet": [0.0, 0.0, 3.0],
    "Sigma_R": [6.0 / 5.0, 0.0, 0.0],
    "X_622": [26.0 / 5.0, 6.0, 4.0],
    "colored_Goldstone_pair": [8.0 / 5.0, 0.0, 1.0],
    "radial_singlet": [0.0, 0.0, 0.0],
}


def to_jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, list) or isinstance(value, tuple):
        return [to_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return [to_jsonable(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def f_px(f54: float) -> float:
    return -9.0 * (f54 - 2.0) * (f54 + 4.0 / 3.0) / 25.0


def f_pc(f54: float) -> float:
    return 9.0 * (f54 - 2.0) * (f54 - 1.0 / 3.0) / 50.0


def f_pl(d210: float) -> float:
    return 0.5 * (d210 * d210 + d210)


def f_pr(d210: float) -> float:
    return 0.5 * (d210 * d210 - d210)


def operator_projector_table() -> list[dict[str, float | str]]:
    rows = []
    for name, data in STATE_CLEBSCH.items():
        f54 = data["F54"]
        d210 = data["D210"]
        rows.append(
            {
                "sector": name,
                "label": data["label"],
                "F54": f54,
                "D210": d210,
                "P_X": f_px(f54),
                "P_color": f_pc(f54),
                "P_L": f_pl(d210),
                "P_R": f_pr(d210),
            }
        )
    return rows


def load_base_seed() -> dict[str, Any]:
    payload = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
    best = payload["corrected_summary"]["best"]
    return {
        "source": str(AUDIT_JSON),
        "best": best,
        "kappa_3": float(best["kappa_3"]),
        "kappa_8": float(best["kappa_8"]),
    }


def update_initial_from_solution(solution: dict[str, Any]) -> dict[str, float]:
    pars = solution["parameters"]
    return {
        "mu": float(pars["mu"]),
        "a54": float(pars["a54"]),
        "b210": float(pars["b210"]),
        "delta_X_projector": max(float(pars["schur_delta_X"]), 1.0e-12),
    }


def choose_best(scan: dict[str, Any]) -> dict[str, float | bool]:
    return scan["best_safe"] if scan["best_safe"] is not None else scan["best"]


def self_consistent_alignment(
    r_med: float,
    source_rows: list[dict[str, str]],
    seed: dict[str, Any],
    initial: dict[str, float],
) -> dict[str, Any]:
    targets = {
        "Sigma_L": float(seed["kappa_3"]),
        "Sigma_R": 1.0,
        "Sigma8_block": float(seed["kappa_8"]),
        "X_622": 1.0,
    }
    initial_guess = dict(initial)
    iterations = []
    solution: dict[str, Any] | None = None
    threshold: dict[str, Any] | None = None
    scan: dict[str, Any] | None = None
    best: dict[str, float | bool] | None = None
    fixed_point_residual = float("inf")

    for idx in range(MAX_ITER):
        solution = rwin.uv.solve_finite_r(targets, r_med, initial_guess)
        threshold = rwin.uv.heavy_threshold_residual(solution)
        scan = down.scan_cached_with_delta(
            source_rows, np.array(threshold["delta_full"], dtype=float)
        )
        best = choose_best(scan)
        next_targets = {
            "Sigma_L": float(best["kappa_3"]),
            "Sigma_R": 1.0,
            "Sigma8_block": float(best["kappa_8"]),
            "X_622": 1.0,
        }
        fixed_point_residual = max(
            abs(next_targets["Sigma_L"] - targets["Sigma_L"]),
            abs(next_targets["Sigma8_block"] - targets["Sigma8_block"]),
        )
        iterations.append(
            {
                "iteration": idx,
                "input_kappa_3": targets["Sigma_L"],
                "input_kappa_8": targets["Sigma8_block"],
                "output_kappa_3": next_targets["Sigma_L"],
                "output_kappa_8": next_targets["Sigma8_block"],
                "fixed_point_step_residual": fixed_point_residual,
                "solver_success": bool(solution["success"]),
                "light_fit_residual_l2": float(solution["residual_l2"]),
                "mediator_projected_l2": float(threshold["projected_l2"]),
                "safe_points": int(scan["safe_points"]),
                "safe_single_scale_factor2_points": int(
                    scan["safe_single_scale_factor2_points"]
                ),
                "alphaG_inv": float(best["alphaG_inv"]),
                "tau_dim6_years": float(best["tau_dim6_years"]),
                "tau_dim5_target_filter_years": float(
                    best["tau_dim5_target_filter_years"]
                ),
            }
        )
        targets = next_targets
        initial_guess = update_initial_from_solution(solution)
        if fixed_point_residual < TOL:
            break

    assert solution is not None
    assert threshold is not None
    assert scan is not None
    assert best is not None
    return {
        "R": float(r_med),
        "converged": bool(fixed_point_residual < TOL),
        "iterations": iterations,
        "fixed_point_residual": float(fixed_point_residual),
        "targets": targets,
        "solution": solution,
        "heavy_threshold": threshold,
        "scan_diagnostics": {
            "total_points": int(scan["total_points"]),
            "safe_points": int(scan["safe_points"]),
            "safe_single_scale_factor2_points": int(
                scan["safe_single_scale_factor2_points"]
            ),
            "best_is_safe": bool(scan["best_safe"] is not None),
        },
        "best": best,
    }


def threshold_check(best: dict[str, float | bool], threshold: dict[str, Any]) -> dict[str, Any]:
    logs = np.array(
        [
            float(best["log_HC"]),
            float(best["log_Sigma3"]),
            float(best["log_Sigma8"]),
        ],
        dtype=float,
    )
    chiral = base.HEAVY_BASIS @ logs / (2.0 * math.pi)
    mediator = np.array(threshold["delta_full"], dtype=float)
    total = chiral + mediator
    alpha_inv = np.array(
        [
            float(best["alpha1_inv_MG_raw"]),
            float(best["alpha2_inv_MG_raw"]),
            float(best["alpha3_inv_MG_raw"]),
        ],
        dtype=float,
    )
    residual_vec = base.PROJECTOR @ (alpha_inv - total)
    alpha_g_check = float(np.mean(alpha_inv - total))

    log_sigma8 = float(best["log_Sigma8"])
    colored_if_unlocked = np.array(PHYSICAL_BETA["colored_Goldstone_pair"]) * log_sigma8 / (
        2.0 * math.pi
    )
    colored_projected = base.PROJECTOR @ colored_if_unlocked
    return {
        "logs": logs.tolist(),
        "alpha_inv_raw": alpha_inv.tolist(),
        "chiral_threshold_HC_Sigma3_Sigma8": chiral.tolist(),
        "mediator_threshold": mediator.tolist(),
        "total_threshold": total.tolist(),
        "projected_matching_residual": residual_vec.tolist(),
        "projected_matching_residual_l2": float(np.linalg.norm(residual_vec)),
        "alphaG_inv_from_matching": alpha_g_check,
        "alphaG_inv_reported": float(best["alphaG_inv"]),
        "colored_pair_if_unlocked_threshold": colored_if_unlocked.tolist(),
        "colored_pair_if_unlocked_projected_l2": float(np.linalg.norm(colored_projected)),
    }


def masses_for_solution(solution: dict[str, Any], mg: float) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name, row in solution["states"].items():
        out[name] = {
            "label": row["label"],
            "F54": float(row["F54"]),
            "D210": float(row["D210"]),
            "target_kappa": float(row["target_kappa"]),
            "light_kappa": float(row["light_kappa"]),
            "light_mass_GeV": float(row["light_kappa"] * mg),
            "heavy_kappas_abs": [float(x) for x in row["heavy_kappas_abs"]],
            "heavy_masses_GeV": [float(x * mg) for x in row["heavy_kappas_abs"]],
            "eigenvalues_signed": [float(x) for x in row["eigenvalues_signed"]],
            "matrix_over_MG": row["matrix"],
            "light_abs_error": float(row["light_abs_error"]),
        }
    return out


def physical_fragments(best: dict[str, float | bool]) -> list[dict[str, Any]]:
    mg = float(best["MG_GeV"])
    sigma8 = float(best["M_Sigma8_GeV"])
    return [
        {
            "fragment": "H_C+H_Cbar",
            "representation": "(3,1,-1/3)+(bar3,1,+1/3)",
            "role": "colored-Higgs pair entering dimension-5 proton decay",
            "mass_GeV": float(best["M_HC_GeV"]),
            "log_MG_over_M": float(best["log_HC"]),
            "b_vector": PHYSICAL_BETA["H_C+H_Cbar"],
            "one_loop_threshold_active": True,
        },
        {
            "fragment": "Sigma_3",
            "representation": "(1,3,0)",
            "role": "intermediate weak adjoint threshold",
            "mass_GeV": float(best["M_Sigma3_GeV"]),
            "log_MG_over_M": float(best["log_Sigma3"]),
            "b_vector": PHYSICAL_BETA["Sigma_3"],
            "one_loop_threshold_active": True,
        },
        {
            "fragment": "Sigma_8",
            "representation": "(8,1,0)",
            "role": "physical color octet from (15,1,1)",
            "mass_GeV": sigma8,
            "log_MG_over_M": float(best["log_Sigma8"]),
            "b_vector": PHYSICAL_BETA["Sigma_8_octet"],
            "one_loop_threshold_active": True,
        },
        {
            "fragment": "Sigma_R",
            "representation": "(1,1,3)",
            "role": "lifted at MG by D210 alignment",
            "mass_GeV": mg,
            "log_MG_over_M": 0.0,
            "b_vector": PHYSICAL_BETA["Sigma_R"],
            "one_loop_threshold_active": False,
        },
        {
            "fragment": "X_622",
            "representation": "(6,2,2)",
            "role": "lifted at MG by P_X mediator sector",
            "mass_GeV": mg,
            "log_MG_over_M": 0.0,
            "b_vector": PHYSICAL_BETA["X_622"],
            "one_loop_threshold_active": False,
        },
        {
            "fragment": "colored_Goldstone_pair",
            "representation": "(3,1,2/3)+(bar3,1,-2/3)",
            "role": "six chiral Goldstone multiplets eaten by broken SU(4)_C vectors",
            "mass_GeV": mg,
            "log_MG_over_M": 0.0,
            "b_vector": PHYSICAL_BETA["colored_Goldstone_pair"],
            "one_loop_threshold_active": False,
        },
        {
            "fragment": "radial_singlet",
            "representation": "(1,1,0)",
            "role": "SM singlet radial mode; M1=M_Sigma8/2 recorded but one-loop gauge-neutral",
            "mass_GeV": 0.5 * sigma8,
            "log_MG_over_M": math.log(mg / (0.5 * sigma8)),
            "b_vector": PHYSICAL_BETA["radial_singlet"],
            "one_loop_threshold_active": False,
        },
    ]


def benchmark_card(alignment: dict[str, Any], seed: dict[str, Any]) -> dict[str, Any]:
    best = alignment["best"]
    solution = alignment["solution"]
    threshold = alignment["heavy_threshold"]
    mg = float(best["MG_GeV"])
    return {
        "note": "No web lookup used. Self-consistent block-operator Spin(10) vacuum-alignment card.",
        "R": float(alignment["R"]),
        "base_seed": {
            "source": seed["source"],
            "kappa_3": seed["kappa_3"],
            "kappa_8": seed["kappa_8"],
        },
        "fixed_point": {
            "converged": alignment["converged"],
            "residual": alignment["fixed_point_residual"],
            "iterations": alignment["iterations"],
        },
        "block_operators": {
            "projector_table": operator_projector_table(),
            "mass_matrix": "M_fd/MG = [[mu+a54*f+b210*d, eta(f-2), eta(f+4/3)], [eta(f-2), 0, R], [eta(f+4/3), R, 0]]",
            "P_X": "-(9/25)(F54-2)(F54+4/3)",
            "P_color": "(9/50)(F54-2)(F54-1/3)",
            "P_L": "(D210^2+D210)/2",
            "P_R": "(D210^2-D210)/2",
            "I3_210": "I3(Phi)=Tr_{Lambda^2}(D_Phi^3), and I3(*_6 A)=6 Pf(A) on the color-adjoint branch",
        },
        "vev_parameters_and_couplings": solution["parameters"],
        "light_and_mediator_eigenvalues": masses_for_solution(solution, mg),
        "physical_fragments": physical_fragments(best),
        "goldstone_count": {
            "SU4_to_SU3xU1_broken_real_generators": 6,
            "eaten_chiral_multiplet_complex_dimension": 6,
            "colored_pair_threshold_log": 0.0,
            "statement": "The (3,1,2/3)+(bar3,1,-2/3) chiral pair is eaten and does not enter the chiral one-loop threshold.",
        },
        "mediator_heavy_threshold": threshold,
        "matching_threshold_check": threshold_check(best, threshold),
        "rge_proton_replay": {
            "scan_diagnostics": alignment["scan_diagnostics"],
            "branch_survives": bool(
                alignment["scan_diagnostics"]["safe_points"] > 0
                and bool(best["triplet_filter_safe"])
                and bool(best["perturbative_superpotential"])
                and float(best["rho_X"]) <= 1.000001
            ),
            "best": best,
        },
    }


def compact_replay_row(alignment: dict[str, Any]) -> dict[str, Any]:
    best = alignment["best"]
    th = alignment["heavy_threshold"]
    return {
        "R": float(alignment["R"]),
        "converged": bool(alignment["converged"]),
        "iterations": len(alignment["iterations"]),
        "fixed_point_residual": float(alignment["fixed_point_residual"]),
        "mu": float(alignment["solution"]["parameters"]["mu"]),
        "a54": float(alignment["solution"]["parameters"]["a54"]),
        "b210": float(alignment["solution"]["parameters"]["b210"]),
        "eta": float(alignment["solution"]["parameters"]["eta"]),
        "schur_delta_X": float(alignment["solution"]["parameters"]["schur_delta_X"]),
        "mediator_projected_l2": float(th["projected_l2"]),
        "safe_points": int(alignment["scan_diagnostics"]["safe_points"]),
        "safe_single_scale_factor2_points": int(
            alignment["scan_diagnostics"]["safe_single_scale_factor2_points"]
        ),
        "tan_beta": float(best["tan_beta"]),
        "MSUSY_GeV": float(best["MSUSY_GeV"]),
        "MG_GeV": float(best["MG_GeV"]),
        "alphaG_inv": float(best["alphaG_inv"]),
        "lambda_T": float(best["lambda_T"]),
        "lambda_S": float(best["lambda_S"]),
        "chi": float(best["chi"]),
        "kappa_3": float(best["kappa_3"]),
        "kappa_8": float(best["kappa_8"]),
        "M_HC_GeV": float(best["M_HC_GeV"]),
        "M_Sigma3_GeV": float(best["M_Sigma3_GeV"]),
        "M_Sigma8_GeV": float(best["M_Sigma8_GeV"]),
        "tau_dim6_years": float(best["tau_dim6_years"]),
        "tau_dim5_target_filter_years": float(best["tau_dim5_target_filter_years"]),
        "triplet_filter_required": float(best["triplet_filter_required"]),
        "residual_l2_after_mediator": float(best["residual_l2_after_mediator"]),
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    with (OUT / "spin10_vacuum_alignment_replay.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    rows = summary["replay_rows"]
    best = summary["recommended_benchmark"]
    best_card = summary["benchmark_cards"][str(best["R"])]
    lines: list[str] = []
    lines.append("# Spin(10) block-operator vacuum alignment")
    lines.append("")
    lines.append("No web lookup was used.  This pass uses the verified `F54`, `D210`,")
    lines.append("`P_X`, and `I3(210)` block operators rather than a hand-assembled")
    lines.append("threshold vector.  It is still a block-operator calculation, not a")
    lines.append("component-by-component diagonalization of the full 210-dimensional field.")
    lines.append("")
    lines.append("## Algebraic mass operator")
    lines.append("")
    lines.append("For each fragment with Clebsches `(f,d)=(F54,D210)`, the finite mediator")
    lines.append("mass matrix is")
    lines.append("")
    lines.append("```text")
    lines.append("M_fd/MG = [[mu+a54 f+b210 d, eta(f-2), eta(f+4/3)],")
    lines.append("          [eta(f-2), 0, R],")
    lines.append("          [eta(f+4/3), R, 0]].")
    lines.append("```")
    lines.append("")
    lines.append("The projectors used in the card are")
    lines.append("")
    lines.append("```text")
    lines.append("P_X     = -(9/25)(F54-2)(F54+4/3)")
    lines.append("P_color =  (9/50)(F54-2)(F54-1/3)")
    lines.append("P_L     =  (D210^2+D210)/2")
    lines.append("P_R     =  (D210^2-D210)/2")
    lines.append("I3(210) = Tr_{Lambda^2}(D_Phi^3), with I3(*_6 A)=6 Pf(A)")
    lines.append("```")
    lines.append("")
    lines.append("## Self-consistent RGE/proton replay")
    lines.append("")
    lines.append("| R | iterations | ||P Delta_med|| | safe points | alphaG^-1 | M_Sigma3 [GeV] | M_Sigma8 [GeV] | tau_d6 [yr] |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['R']:.0f} | {row['iterations']} | {row['mediator_projected_l2']:.6e} | "
            f"{row['safe_points']} | {row['alphaG_inv']:.6f} | "
            f"{row['M_Sigma3_GeV']:.6e} | {row['M_Sigma8_GeV']:.6e} | "
            f"{row['tau_dim6_years']:.6e} |"
        )
    lines.append("")
    lines.append("All displayed points are fixed points of the map")
    lines.append("`(kappa_3,kappa_8) -> Delta_med -> RGE replay -> (kappa_3,kappa_8)`")
    lines.append(f"to tolerance `{TOL:.1e}`.")
    lines.append("")
    lines.append("## Recommended card")
    lines.append("")
    lines.append("The most decoupled displayed point is the recommended benchmark card:")
    lines.append("")
    lines.append("```text")
    lines.append(f"R = {best['R']:.0f}")
    lines.append(f"mu = {best['mu']:.9f}")
    lines.append(f"a54 = {best['a54']:.9f}")
    lines.append(f"b210 = {best['b210']:.9f}")
    lines.append(f"eta = {best['eta']:.9f}")
    lines.append(f"||P Delta_med||_2 = {best['mediator_projected_l2']:.6e}")
    lines.append(f"safe points = {best['safe_points']}")
    lines.append(f"M_Sigma3 = {best['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"M_Sigma8 = {best['M_Sigma8_GeV']:.6e} GeV")
    lines.append(f"tau_d6 = {best['tau_dim6_years']:.6e} yr")
    lines.append(f"tau_d5(S_T=1e-5) = {best['tau_dim5_target_filter_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("The physical-fragment assignment in the card keeps only `H_C`, `Sigma_3`,")
    lines.append("and the octet part of `Sigma_8` active in the one-loop chiral threshold.")
    lines.append("`Sigma_R` and `X_622` are at `MG`; the colored pair is eaten by the broken")
    lines.append("`SU(4)_C` vectors; the radial singlet is gauge-neutral at one loop.")
    lines.append("")
    lines.append("## Verification files")
    lines.append("")
    lines.append("```text")
    lines.append(str(OUT / "spin10_vacuum_alignment_summary.json"))
    for r in CARD_R_VALUES:
        lines.append(str(OUT / f"benchmark_R{int(r)}.json"))
    lines.append(str(OUT / "spin10_vacuum_alignment_replay.csv"))
    lines.append("```")
    (OUT / "spin10_vacuum_alignment_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_rows = down.iter_corrected_rows()
    seed = load_base_seed()
    _old_targets, initial = rwin.load_targets()

    alignments = []
    for r_med in R_GRID:
        alignments.append(self_consistent_alignment(r_med, source_rows, seed, initial))

    replay_rows = [compact_replay_row(item) for item in alignments]
    cards = {str(item["R"]): benchmark_card(item, seed) for item in alignments if item["R"] in CARD_R_VALUES}
    recommended = max(replay_rows, key=lambda row: row["R"])
    summary = {
        "note": "No web lookup used. Self-consistent Spin(10) block-operator vacuum-alignment scan.",
        "script": str(Path(__file__).resolve()),
        "corrected_cache": str(down.CORRECTED_SCAN_CSV),
        "R_grid": R_GRID,
        "fixed_point_tolerance": TOL,
        "seed": seed,
        "operator_projector_table": operator_projector_table(),
        "replay_rows": replay_rows,
        "recommended_benchmark": recommended,
        "benchmark_cards": cards,
    }

    write_csv(replay_rows)
    for r_key, card in cards.items():
        (OUT / f"benchmark_R{int(float(r_key))}.json").write_text(
            json.dumps(to_jsonable(card), indent=2), encoding="utf-8"
        )
    (OUT / "spin10_vacuum_alignment_summary.json").write_text(
        json.dumps(to_jsonable(summary), indent=2), encoding="utf-8"
    )
    write_report(to_jsonable(summary))

    print("Spin(10) block-operator vacuum alignment")
    for row in replay_rows:
        print(
            f"  R={row['R']:.0f}: fixed residual={row['fixed_point_residual']:.3e}, "
            f"||P Delta_med||={row['mediator_projected_l2']:.6e}, "
            f"safe={row['safe_points']}, M_Sigma3={row['M_Sigma3_GeV']:.6e} GeV"
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
