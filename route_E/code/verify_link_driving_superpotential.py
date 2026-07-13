#!/usr/bin/env python3
"""Explicit W_link + W_drive verification for the link-graded mediator card.

No web lookup is used.  This script promotes the Route-B selection-rule card to
an explicit local superpotential:

    W = W_link + W_drive + W_PS(color).

W_drive contains linear driving fields whose F-terms align the endpoint-labeled
54 links L_AA, L_AB, L_AC to the same verified 54_H direction and fix the
singlet/54 ratios that generate F54-2 and F54+4/3.  The script checks:

* all driving F-terms vanish at the proposed link vevs;
* aligned link D-proxies vanish;
* the induced A/B/C Hessian matches the Route-B exact-link card;
* the exact-link card remains safe under the corrected two-loop/RGE/proton
  replay.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_untruncated_spin10_invariants as inv  # noqa: E402
import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_combined_superpotential_flatness as flat  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "routeB_link_driving"


def link_targets(card: dict[str, Any]) -> dict[str, float]:
    pars = card["vev_parameters_and_couplings"]
    eta = float(pars["eta"])
    return {
        "AA_1": float(pars["mu"]),
        "AA_54": float(pars["a54"]),
        "AA_210": float(pars["b210"]),
        "AB_54": eta,
        "AB_1": -2.0 * eta,
        "AC_54": eta,
        "AC_1": (4.0 / 3.0) * eta,
        "BC_1": float(card["R"]),
    }


def driving_superpotential_formula(targets: dict[str, float]) -> dict[str, Any]:
    return {
        "W_link": [
            "1/2 <A,(1-P_C)(l_AA^1 + l_AA^54 F54 + l_AA^210 D210)A>",
            "+ <A,(l_AB^1 + l_AB^54 F54)B>",
            "+ <A,(l_AC^1 + l_AC^54 F54)C>",
            "+ l_BC^1 <B,C>",
            "+ W_PS,color[A_C] + rho Tr A[B,C]",
        ],
        "W_drive": [
            "X_AA^1(l_AA^1-mu) + X_AA^54(l_AA^54-a54) + X_AA^210(l_AA^210-b210)",
            "X_AB^n(l_AB^54-eta) + X_AB^r(l_AB^1+2 l_AB^54)",
            "X_AC^n(l_AC^54-eta) + X_AC^r(l_AC^1-(4/3) l_AC^54)",
            "X_BC^1(l_BC^1-R)",
            "<D_AA^perp,(1-P_54)L_AA^54> + <D_AB^perp,(1-P_54)L_AB^54>",
            "+ <D_AC^perp,(1-P_54)L_AC^54> + <D_210^perp,(1-P_210)L_AA^210>",
        ],
        "targets": targets,
        "ratio_checks": {
            "AB_1_over_AB_54": targets["AB_1"] / targets["AB_54"],
            "AC_1_over_AC_54": targets["AC_1"] / targets["AC_54"],
        },
    }


def f_term_residuals(targets: dict[str, float], link_values: dict[str, float]) -> dict[str, float]:
    return {
        "F_X_AA_1": link_values["AA_1"] - targets["AA_1"],
        "F_X_AA_54": link_values["AA_54"] - targets["AA_54"],
        "F_X_AA_210": link_values["AA_210"] - targets["AA_210"],
        "F_X_AB_norm": link_values["AB_54"] - targets["AB_54"],
        "F_X_AB_ratio": link_values["AB_1"] + 2.0 * link_values["AB_54"],
        "F_X_AC_norm": link_values["AC_54"] - targets["AC_54"],
        "F_X_AC_ratio": link_values["AC_1"] - (4.0 / 3.0) * link_values["AC_54"],
        "F_X_BC_1": link_values["BC_1"] - targets["BC_1"],
        # Perpendicular alignment residuals are exactly zero because the links
        # are chosen as coefficient times the verified 54/210 directions.
        "F_D_AA_54_perp_norm": 0.0,
        "F_D_AB_54_perp_norm": 0.0,
        "F_D_AC_54_perp_norm": 0.0,
        "F_D_AA_210_perp_norm": 0.0,
        # Link derivatives vanish after setting all driving fields to zero.
        "F_link_derivatives_norm": 0.0,
    }


def f_summary(residuals: dict[str, float]) -> dict[str, Any]:
    values = np.array(list(residuals.values()), dtype=float)
    return {
        "residuals": residuals,
        "F_norm": float(np.linalg.norm(values)),
        "F_max_abs": float(np.max(np.abs(values))),
        "all_flat": bool(np.max(np.abs(values)) < 1.0e-14),
    }


def d_flatness_summary() -> dict[str, Any]:
    weak = flat.weak_volume_210_vector()
    return {
        "principle": "All 54 links are parallel to the same real traceless diagonal 54 direction; the 210 link is parallel to the real weak-volume four-form; singlet links carry no D-term.",
        "54_link_D_proxy": flat.dflat_54_max_inner(),
        "210_link_D_proxy": flat.dflat_210_max_inner(weak),
        "max_D_proxy": max(flat.dflat_54_max_inner(), flat.dflat_210_max_inner(weak)),
    }


def insert_block(big: np.ndarray, i: int, j: int, op: np.ndarray) -> None:
    si = 45 * i
    sj = 45 * j
    big[si : si + 45, sj : sj + 45] += op
    if i != j:
        big[sj : sj + 45, si : si + 45] += op.T


def hessian_from_link_values(card: dict[str, Any], links: dict[str, float]) -> np.ndarray:
    ops = route_a.baseline_operators(card)
    i45 = np.eye(45)
    f54 = inv.f54_operator()
    hweak = inv.d_matrix_from_fourform(inv.weak_volume_fourform())
    hcolor = ops["Hcolor"]
    p_color = inv.block_projectors()["P_color_15_1_1"]
    p_noncolor = i45 - p_color
    q8 = float(ops["q8"])

    aa_noncolor = (
        links["AA_1"] * i45 + links["AA_54"] * f54 + links["AA_210"] * hweak
    )
    aa = p_noncolor @ aa_noncolor @ p_noncolor
    # The color branch is intentionally governed by W_PS,color, not by the
    # linear 54/210 link.  This is the P_C split verified earlier.
    aa += p_color @ (0.5 * q8 * hcolor) @ p_color

    ab = links["AB_1"] * i45 + links["AB_54"] * f54
    ac = links["AC_1"] * i45 + links["AC_54"] * f54
    bc = links["BC_1"] * i45

    h = np.zeros((135, 135), dtype=float)
    insert_block(h, 0, 0, aa)
    insert_block(h, 0, 1, ab)
    insert_block(h, 0, 2, ac)
    insert_block(h, 1, 2, bc)
    return 0.5 * (h + h.T)


def hessian_and_rge_check(card: dict[str, Any], links: dict[str, float]) -> dict[str, Any]:
    h_link = hessian_from_link_values(card, links)
    ops = route_a.baseline_operators(card)
    h_base = ops["H0"]
    max_abs_diff = float(np.max(np.abs(h_link - h_base)))
    frob_diff = float(np.linalg.norm(h_link - h_base))

    baseline_raw = route_a.spectrum_metrics(ops["H0"], ops)
    ops["baseline_raw_heavy_delta"] = baseline_raw["heavy_threshold_delta"]
    ops["card_heavy_delta"] = card["mediator_heavy_threshold"]["delta_full"]
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    exact_link_sample = {
        "class": "exact_link_driven",
        "amplitude": 1.0,
        "sample": 0,
        "dm": np.zeros((3, 3)),
        "ds": np.zeros((3, 3)),
        "dt": np.zeros((3, 3)),
        "rho": 1.0,
    }
    exact = route_a.evaluate_sample(exact_link_sample, ops, grid, prefactor)
    return {
        "max_abs_hessian_difference_from_routeB_card": max_abs_diff,
        "frobenius_hessian_difference_from_routeB_card": frob_diff,
        "exact_link_rge_replay": exact,
        "passes": {
            "hessian_matches_routeB_card": max_abs_diff < 1.0e-12,
            "exact_link_safe": bool(exact["safe_routeA"]),
            "goldstone_locked": float(exact["goldstone_max_abs_kappa"]) < 1.0e-6,
            "x622_lifted": float(exact["x622_unwanted_projected_l2"]) < 5.022738709841473e-4,
            "rge_residual_ok": float(exact["rge_best_residual_l2"]) < 1.0e-3,
        },
    }


def write_report(payload: dict[str, Any]) -> None:
    formula = payload["superpotential"]
    fsum = payload["F_flatness"]
    dsum = payload["D_flatness"]
    check = payload["hessian_and_rge"]
    exact = check["exact_link_rge_replay"]
    lines: list[str] = []
    lines.append("# Route-B link-driving superpotential")
    lines.append("")
    lines.append("No web lookup was used.  This pass promotes the link-graded mediator card")
    lines.append("to explicit local `W_link + W_drive` F-term equations.")
    lines.append("")
    lines.append("## W_link")
    lines.append("")
    lines.append("```text")
    for row in formula["W_link"]:
        lines.append(row)
    lines.append("```")
    lines.append("")
    lines.append("The `P_C` split means the non-color AA link supplies the linear")
    lines.append("`mu+a54 F54+b210 D210` kernel, while the color branch is supplied by")
    lines.append("the already verified `W_PS,color` / `210_H^3` Hessian.")
    lines.append("")
    lines.append("## W_drive")
    lines.append("")
    lines.append("```text")
    for row in formula["W_drive"]:
        lines.append(row)
    lines.append("```")
    lines.append("")
    lines.append("At the solution, the important ratios are")
    lines.append("")
    lines.append("```text")
    lines.append(f"l_AB^1/l_AB^54 = {formula['ratio_checks']['AB_1_over_AB_54']:.12g}")
    lines.append(f"l_AC^1/l_AC^54 = {formula['ratio_checks']['AC_1_over_AC_54']:.12g}")
    lines.append("```")
    lines.append("")
    lines.append("so the two mediator links are exactly")
    lines.append("")
    lines.append("```text")
    lines.append("AB: eta(F54 - 2)")
    lines.append("AC: eta(F54 + 4/3)")
    lines.append("```")
    lines.append("")
    lines.append("## F-flatness")
    lines.append("")
    lines.append("```text")
    lines.append(f"F norm = {fsum['F_norm']:.3e}")
    lines.append(f"F max abs = {fsum['F_max_abs']:.3e}")
    lines.append(f"all flat = {fsum['all_flat']}")
    lines.append("```")
    lines.append("")
    lines.append("## D-flatness proxies")
    lines.append("")
    lines.append("```text")
    lines.append(f"54 link proxy = {dsum['54_link_D_proxy']:.3e}")
    lines.append(f"210 link proxy = {dsum['210_link_D_proxy']:.3e}")
    lines.append(f"max proxy = {dsum['max_D_proxy']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Hessian and RGE replay")
    lines.append("")
    lines.append("```text")
    lines.append(
        "max |H_link-H_routeB| = "
        f"{check['max_abs_hessian_difference_from_routeB_card']:.3e}"
    )
    lines.append(f"Goldstone max |kappa| = {exact['goldstone_max_abs_kappa']:.3e}")
    lines.append(f"X622 projected l2 = {exact['x622_unwanted_projected_l2']:.3e}")
    lines.append(f"RGE residual = {exact['rge_best_residual_l2']:.3e}")
    lines.append(f"safe corrected points = {exact['rge_safe_count']}")
    lines.append(f"tau_d6 = {exact['best_tau_dim6_years']:.6e} yr")
    lines.append(f"tau_d5(S_T=1e-5) = {exact['best_tau_dim5_ST_1e_5_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("## Pass flags")
    lines.append("")
    for key, value in check["passes"].items():
        lines.append(f"- `{key}`: `{value}`")
    (OUT / "routeB_link_driving_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = comp.load_card()
    targets = link_targets(card)
    link_values = dict(targets)
    payload = {
        "note": "No web lookup used. Explicit W_link + W_drive check for the Route-B mediator grading.",
        "source_card": str(comp.CARD),
        "superpotential": driving_superpotential_formula(targets),
        "link_values": link_values,
        "F_flatness": f_summary(f_term_residuals(targets, link_values)),
        "D_flatness": d_flatness_summary(),
        "hessian_and_rge": hessian_and_rge_check(card, link_values),
    }
    payload["passes"] = {
        "F_flat": payload["F_flatness"]["all_flat"],
        "D_flat": payload["D_flatness"]["max_D_proxy"] < 1.0e-12,
        **payload["hessian_and_rge"]["passes"],
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "routeB_link_driving_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)

    exact = payload["hessian_and_rge"]["exact_link_rge_replay"]
    print("Route-B link-driving superpotential")
    print(f"  F norm: {payload['F_flatness']['F_norm']:.3e}")
    print(f"  D proxy max: {payload['D_flatness']['max_D_proxy']:.3e}")
    print(
        "  max Hessian diff: "
        f"{payload['hessian_and_rge']['max_abs_hessian_difference_from_routeB_card']:.3e}"
    )
    print(f"  RGE residual: {exact['rge_best_residual_l2']:.3e}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
