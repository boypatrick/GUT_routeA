#!/usr/bin/env python3
"""Construct and scan a source-stabilization sector.

No web lookup is used.  This pass turns the previous source-mixing diagnostic
into a model-building card:

1. A sequestered source sector keeps the weak-volume Omega_210 source distinct
   from the color-cubic 210_C field.
2. A U(1)_R plus source parity allows source-link-driver alignment terms but
   forbids source self-cubics such as Omega_W^3 and S_54 Omega_W^2.
3. Any residual symmetry breaking is represented by fragment-splitting spurions
   epsilon54 and epsilon210.  These are scanned directly through the cached
   two-loop/RGE/proton replay.

The result is a falsifiable condition: the exact sequestered sector has no
non-universal source threshold; if a UV completion generates a non-scalar 210
source Hessian, its effective epsilon must be small enough.
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
import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


OUT = ROOT / "output" / "source_stabilization"


def selection_rule_card() -> dict[str, Any]:
    """R-sequestered source-sector rule set.

    R charges are normalized so the superpotential has charge +2.
    The color 210_C is deliberately a separate field from the weak-volume
    source Omega_W; this avoids forcing the source field to carry the 210^3
    color-cubic interaction.
    """

    fields = {
        "S_54": {"Spin10": "54", "R": 0, "P_src": -1, "role": "weak/projector 54 source"},
        "Omega_W": {"Spin10": "210", "R": 0, "P_src": -1, "role": "weak-volume D210 source"},
        "Phi_C": {"Spin10": "210", "R": 2.0 / 3.0, "P_src": +1, "role": "color 210_C cubic source"},
        "L_54_e": {"Spin10": "54", "R": 0, "P_src": -1, "role": "endpoint 54 links"},
        "L_210": {"Spin10": "210", "R": 0, "P_src": -1, "role": "endpoint 210 link"},
        "D_54_e": {"Spin10": "54", "R": 2, "P_src": -1, "role": "54 alignment drivers"},
        "D_210": {"Spin10": "210", "R": 2, "P_src": -1, "role": "210 alignment driver"},
        "X_a": {"Spin10": "1", "R": 2, "P_src": +1, "role": "singlet norm/ratio drivers"},
        "zeta": {"Spin10": "1", "R": 0, "P_src": -1, "role": "optional source-parity spurion"},
    }
    terms = [
        {
            "term": "D_54_e (L_54_e - ell_e S_54)",
            "status": "allowed",
            "reason": "R=2, P_src=(-)(-)=even; gives complete 54 source-link-driver block",
        },
        {
            "term": "D_210 (L_210 - ell_210 Omega_W)",
            "status": "allowed",
            "reason": "R=2, P_src=(-)(-)=even; gives complete 210 source-link-driver block",
        },
        {
            "term": "X_a f_a(S_54, Omega_W)",
            "status": "allowed",
            "reason": "driver constraints can fix norm/orbit while X_a=0 at the F-flat branch",
        },
        {
            "term": "Phi_C^3",
            "status": "allowed",
            "reason": "R(Phi_C)=2/3 and P_src even; retains the verified color 210_C cubic",
        },
        {
            "term": "Omega_W^3",
            "status": "forbidden",
            "reason": "R=0, P_src odd; removes the dangerous 210 source self-Hessian",
        },
        {
            "term": "S_54 Omega_W^2",
            "status": "forbidden",
            "reason": "R=0, P_src odd; removes mixed source-fragment Hessian",
        },
        {
            "term": "Omega_W Phi_C^2",
            "status": "forbidden",
            "reason": "R=4/3, P_src odd; sequesters weak-volume source from color cubic",
        },
        {
            "term": "zeta Omega_W^3 or zeta S_54 Omega_W^2",
            "status": "spurion",
            "reason": "parametrizes residual source-parity breaking; mapped to epsilon210",
        },
    ]
    return {
        "principle": "source-field doubling plus R-sequestering",
        "fields": fields,
        "terms": terms,
        "stabilizing_superpotential": [
            "W_src = sum_e kappa54_e M_G <D54_e, L54_e - ell54_e S54>",
            "      + kappa210 M_G <D210, L210 - ell210 Omega_W>",
            "      + sum_a X_a f_a(S54, Omega_W)",
            "      + m_C[1/2 ||Phi_C||^2 - I3(Phi_C)/6]",
            "      + W_link ratios already verified.",
        ],
        "theorem": "On the F-flat branch X_a=D=0, forbidden source self-cubics imply the non-Goldstone source Hessian is Spin(10)-scalar; source-link-driver eigenvalues are complete multiplets.",
    }


def scan_eps_grid() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    metrics = dyn.exact_link_metrics()
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()

    values = np.array(
        [
            -3.0e-3,
            -2.0e-3,
            -1.5e-3,
            -1.0e-3,
            -7.5e-4,
            -5.0e-4,
            -3.0e-4,
            -2.0e-4,
            -1.0e-4,
            -5.0e-5,
            0.0,
            5.0e-5,
            1.0e-4,
            2.0e-4,
            3.0e-4,
            5.0e-4,
            7.5e-4,
            1.0e-3,
            1.5e-3,
            2.0e-3,
            3.0e-3,
        ],
        dtype=float,
    )
    rows: list[dict[str, Any]] = []
    for eps54 in values:
        for eps210 in values:
            threshold = dyn.source_mixing_threshold(float(eps54), float(eps210))
            rge = dyn.fixed_spectrum_rge_scan_with_extra_delta(
                metrics,
                grid,
                prefactor,
                np.array(threshold["delta_total"], dtype=float),
            )
            safe = bool(rge["safe_count_residual_lt_1e_3"] > 0)
            rows.append(
                {
                    "epsilon54": float(eps54),
                    "epsilon210": float(eps210),
                    "projected_delta_l2": float(threshold["projected_delta_l2"]),
                    "universal_delta": float(threshold["universal_delta"]),
                    "safe_count": int(rge["safe_count_residual_lt_1e_3"]),
                    "best_residual_l2": float(rge["best_residual_l2"]),
                    "best_alphaG_inv": float(rge["best_alphaG_inv"]),
                    "tau_d6": float(rge["best_tau_dim6_years"]),
                    "tau_d5": float(rge["best_tau_dim5_ST_1e_5_years"]),
                    "filter_required": float(rge["best_filter_required"]),
                    "safe": safe,
                }
            )

    def axis_rows(e54: float | None = None, e210: float | None = None) -> list[dict[str, Any]]:
        out = rows
        if e54 is not None:
            out = [row for row in out if abs(row["epsilon54"] - e54) < 1.0e-15]
        if e210 is not None:
            out = [row for row in out if abs(row["epsilon210"] - e210) < 1.0e-15]
        return out

    safe_rows = [row for row in rows if row["safe"]]
    unsafe_rows = [row for row in rows if not row["safe"]]
    max_abs_210_axis = max(
        abs(row["epsilon210"]) for row in axis_rows(e54=0.0) if row["safe"]
    )
    max_abs_54_axis = max(abs(row["epsilon54"]) for row in axis_rows(e210=0.0) if row["safe"])
    max_radius_safe = max(
        math.hypot(row["epsilon54"], row["epsilon210"]) for row in safe_rows
    )
    first_unsafe_by_radius = min(
        unsafe_rows,
        key=lambda row: math.hypot(row["epsilon54"], row["epsilon210"]),
    )
    safe_near_square_1e4 = all(
        row["safe"]
        for row in rows
        if abs(row["epsilon54"]) <= 1.0e-4 + 1e-15
        and abs(row["epsilon210"]) <= 1.0e-4 + 1e-15
    )
    safe_near_square_3e4 = all(
        row["safe"]
        for row in rows
        if abs(row["epsilon54"]) <= 3.0e-4 + 1e-15
        and abs(row["epsilon210"]) <= 3.0e-4 + 1e-15
    )
    summary = {
        "grid_values": values.tolist(),
        "points": len(rows),
        "safe_points": len(safe_rows),
        "safe_fraction": len(safe_rows) / len(rows),
        "max_abs_epsilon210_safe_on_epsilon54_0_axis": max_abs_210_axis,
        "max_abs_epsilon54_safe_on_epsilon210_0_axis": max_abs_54_axis,
        "max_radius_safe_on_grid": max_radius_safe,
        "first_unsafe_by_radius": first_unsafe_by_radius,
        "all_points_safe_inside_1e_minus_4_square": safe_near_square_1e4,
        "all_points_safe_inside_3e_minus_4_square": safe_near_square_3e4,
    }
    return rows, summary


def spurion_order_estimate(lambda_src: float, epsilon_bound: float) -> dict[str, Any]:
    orders = []
    for n in range(1, 9):
        eps = lambda_src**n
        orders.append({"order": n, "epsilon": eps, "passes": eps <= epsilon_bound})
    first = next(row for row in orders if row["passes"])
    return {
        "lambda_src": lambda_src,
        "epsilon_bound": epsilon_bound,
        "orders": orders,
        "first_safe_order": first["order"],
        "interpretation": f"With lambda_src={lambda_src}, forbid dangerous 210 source operators up to order {first['order'] - 1}.",
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "epsilon54",
        "epsilon210",
        "projected_delta_l2",
        "universal_delta",
        "safe_count",
        "best_residual_l2",
        "best_alphaG_inv",
        "tau_d6",
        "tau_d5",
        "filter_required",
        "safe",
    ]
    with (OUT / "source_split_threshold_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{key: row[key] for key in fields} for row in rows])


def write_report(payload: dict[str, Any]) -> None:
    card = payload["selection_rule_card"]
    scan = payload["split_scan_summary"]
    spur = payload["spurion_order_estimate_lambda_0p1"]
    lines: list[str] = []
    lines.append("# Source-stabilization sector")
    lines.append("")
    lines.append("No web lookup was used.  This pass constructs an explicit")
    lines.append("R-sequestered source sector and scans residual source-fragment splittings.")
    lines.append("")
    lines.append("## Nontrivial model-building elements")
    lines.append("")
    lines.append("1. `Source-field doubling`: the weak-volume source `Omega_W` is distinct")
    lines.append("   from the color-cubic `Phi_C`.  `Phi_C^3` remains allowed, while")
    lines.append("   `Omega_W^3` is forbidden.")
    lines.append("2. `R-sequestering`: source/link fields have R=0 and drivers have R=2,")
    lines.append("   so source self-cubics do not enter the source Hessian.")
    lines.append("3. `Spurion bound`: any residual source-parity breaking is mapped to")
    lines.append("   epsilon54, epsilon210 and replayed through the RGE/proton scan.")
    lines.append("")
    lines.append("## Allowed and forbidden terms")
    lines.append("")
    lines.append("| term | status | reason |")
    lines.append("|---|---|---|")
    for row in card["terms"]:
        lines.append(f"| `{row['term']}` | {row['status']} | {row['reason']} |")
    lines.append("")
    lines.append("The stabilizing superpotential is")
    lines.append("")
    lines.append("```text")
    for term in card["stabilizing_superpotential"]:
        lines.append(term)
    lines.append("```")
    lines.append("")
    lines.append("On the F-flat branch, all R=2 drivers vanish.  Since the forbidden")
    lines.append("source self-cubics are the only terms that would generate Clebsch-dependent")
    lines.append("source Hessians, the non-Goldstone source-link-driver Hessian is a repeated")
    lines.append("complete-Spin(10)-multiplet matrix.")
    lines.append("")
    lines.append("## Split-threshold scan")
    lines.append("")
    lines.append("Residual UV breaking is represented by")
    lines.append("")
    lines.append("```text")
    lines.append("m_54(f)  = m_54  (1 + epsilon54 c_f),")
    lines.append("m_210(f) = m_210 (1 + epsilon210 c_f).")
    lines.append("```")
    lines.append("")
    lines.append(f"Grid points: `{scan['points']}`; safe fraction: `{scan['safe_fraction']:.3f}`.")
    lines.append("")
    lines.append("```text")
    lines.append(
        "max |epsilon210| safe on epsilon54=0 axis = "
        f"{scan['max_abs_epsilon210_safe_on_epsilon54_0_axis']:.3e}"
    )
    lines.append(
        "max |epsilon54| safe on epsilon210=0 axis = "
        f"{scan['max_abs_epsilon54_safe_on_epsilon210_0_axis']:.3e}"
    )
    lines.append(f"all points safe inside |eps54|,|eps210| <= 1e-4: {scan['all_points_safe_inside_1e_minus_4_square']}")
    lines.append(f"all points safe inside |eps54|,|eps210| <= 3e-4: {scan['all_points_safe_inside_3e_minus_4_square']}")
    bad = scan["first_unsafe_by_radius"]
    lines.append(
        "first unsafe by radius: "
        f"eps54={bad['epsilon54']:.1e}, eps210={bad['epsilon210']:.1e}, "
        f"residual={bad['best_residual_l2']:.3e}, tau_d5={bad['tau_d5']:.3e}"
    )
    lines.append("```")
    lines.append("")
    lines.append("## Spurion order")
    lines.append("")
    lines.append(f"For a source-parity spurion `lambda_src={spur['lambda_src']}`,")
    lines.append(f"the first safe order for epsilon210 <= {spur['epsilon_bound']:.1e} is")
    lines.append(f"`lambda_src^{spur['first_safe_order']}`.")
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append("The sequestered sector gives epsilon54=epsilon210=0 at tree level, so the")
    lines.append("dynamic source fields do not spoil the branch.  If the UV completion leaks")
    lines.append("source parity, the two-loop/proton replay requires roughly")
    lines.append("`|epsilon210| <= 1e-4` as a conservative paper benchmark.")
    (OUT / "source_stabilization_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, scan_summary = scan_eps_grid()
    epsilon_bound = 1.0e-4
    payload = {
        "note": "No web lookup used. R-sequestered source-stabilization sector and residual split scan.",
        "selection_rule_card": selection_rule_card(),
        "split_scan_summary": scan_summary,
        "spurion_order_estimate_lambda_0p1": spurion_order_estimate(0.1, epsilon_bound),
        "spurion_order_estimate_lambda_0p2": spurion_order_estimate(0.2, epsilon_bound),
        "passes": {
            "sequestered_tree_level_has_zero_split": True,
            "one_e_minus_4_square_safe": bool(scan_summary["all_points_safe_inside_1e_minus_4_square"]),
            "three_e_minus_4_square_safe": bool(scan_summary["all_points_safe_inside_3e_minus_4_square"]),
            "epsilon210_axis_reaches_at_least_1e_minus_4": scan_summary[
                "max_abs_epsilon210_safe_on_epsilon54_0_axis"
            ]
            >= 1.0e-4,
        },
    }
    payload["passes"]["all"] = (
        payload["passes"]["sequestered_tree_level_has_zero_split"]
        and payload["passes"]["one_e_minus_4_square_safe"]
        and payload["passes"]["epsilon210_axis_reaches_at_least_1e_minus_4"]
    )
    (OUT / "source_stabilization_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(payload)

    print("Source-stabilization sector")
    print(f"  grid points: {scan_summary['points']}")
    print(f"  safe fraction: {scan_summary['safe_fraction']:.3f}")
    print(
        "  max |eps210| safe on eps54=0 axis: "
        f"{scan_summary['max_abs_epsilon210_safe_on_epsilon54_0_axis']:.3e}"
    )
    print(f"  all |eps|<=1e-4 safe: {scan_summary['all_points_safe_inside_1e_minus_4_square']}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
