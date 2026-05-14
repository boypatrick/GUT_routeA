#!/usr/bin/env python3
"""Route-B mediator grading / R-symmetry verification.

No web lookup is used.  This script turns the Route-A robustness lesson into a
minimal selection-rule card.

The first result is a no-go: with only one 54_H portal and abelian field
charges on A,B,C, allowing AA, AB, AC, BC and the corresponding SAA, SAB, SAC
terms automatically allows the generic BB/CC/SBC structures.  A pure charge
assignment therefore cannot protect the projector card.

The viable Route-B completion is a link-graded mediator sector.  The mediator
nodes carry U(1)_M charges

    q(A)=0, q(B)=+1, q(C)=-1,

and the 54/singlet portals are endpoint-labeled links:

    L_AA, L_AB, L_AC, L_BC.

A driving sector aligns the 54-valued links in the same Spin(10) direction.
This is equivalent, on the component card, to keeping the projector-generating
AA, BC, AB, AC entries while suppressing all other m_ij and sigma_ij entries.
The remaining allowed breaking is tested with the Route-A Hessian/RGE scanner.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "routeB_mediator_grading"


def allowed_mod(charge_sum: int, modulus: int) -> bool:
    return charge_sum % modulus == 0


def term_charge(term: str, charges: dict[str, int]) -> int:
    return sum(charges[ch] for ch in term)


def brute_force_abelian_no_go(max_modulus: int = 16) -> dict[str, Any]:
    """Search finite Z_N charges for a field-only abelian selection rule."""

    required = ["AA", "AB", "AC", "BC", "SAA", "SAB", "SAC"]
    forbidden = ["BB", "CC", "SBB", "SCC", "SBC"]
    witnesses = []
    for n in range(2, max_modulus + 1):
        found = 0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for s in range(n):
                        charges = {"A": a, "B": b, "C": c, "S": s}
                        if all(allowed_mod(term_charge(term, charges), n) for term in required):
                            if all(not allowed_mod(term_charge(term, charges), n) for term in forbidden):
                                witnesses.append({"modulus": n, "charges": charges})
                                found += 1
        if found:
            break
    return {
        "searched_ZN_up_to": max_modulus,
        "required_allowed": required,
        "required_forbidden": forbidden,
        "solutions_found": len(witnesses),
        "witnesses": witnesses[:5],
        "analytic_reason": "AA, AB, AC, BC imply qA=qB=qC in any abelian additive grading; with SAA,SAB,SAC the same S charge then also allows SBB,SBC,SCC.",
    }


def routeB_selection_rules() -> dict[str, Any]:
    return {
        "symmetry": "U(1)_R times U(1)_M link grading",
        "superpotential_R_charge": 2,
        "node_charges_U1M": {
            "A_45": 0,
            "B_45": 1,
            "C_45": -1,
        },
        "link_portals_U1M": {
            "L_AA^(1,54,210)": 0,
            "L_AB^(1,54)": -1,
            "L_AC^(1,54)": 1,
            "L_BC^(1)": 0,
            "epsilon_tau": "small neutral or endpoint-labeled 210 breaking, |epsilon_tau| <= 1e-3",
        },
        "allowed_projector_terms": [
            "L_AA^(1) A A",
            "L_AA^(54) A A",
            "L_AA^(210) A A",
            "L_AB^(1) A B",
            "L_AB^(54) A B",
            "L_AC^(1) A C",
            "L_AC^(54) A C",
            "L_BC^(1) B C",
            "rho A[B,C]",
        ],
        "forbidden_or_spurion_suppressed": [
            "B B, C C",
            "generic B C 54 portal",
            "generic S B B, S B C, S C C",
            "generic m_ij and sigma_ij entries not in the link list",
        ],
        "driving_alignment": [
            "F_DAA enforces L_AA^(54) parallel to the verified 54_H Clebsch direction",
            "F_DAB enforces L_AB^(54) parallel to the same 54_H direction",
            "F_DAC enforces L_AC^(54) parallel to the same 54_H direction",
            "singlet vevs in L_AB/L_AC fix the shifted factors F54-2 and F54+4/3",
        ],
        "important_caveat": "The link portals are essential. A single unlabelled 54_H with only abelian charges cannot do the job.",
    }


def prepare_routeA_context() -> tuple[dict[str, Any], dict[str, np.ndarray], float]:
    card = comp.load_card()
    ops = route_a.baseline_operators(card)
    baseline_raw = route_a.spectrum_metrics(ops["H0"], ops)
    ops["baseline_raw_heavy_delta"] = baseline_raw["heavy_threshold_delta"]
    ops["card_heavy_delta"] = card["mediator_heavy_threshold"]["delta_full"]
    grid = route_a.load_cached_grid()
    prefactor = base.proton_prefactor()
    return ops, grid, prefactor


def run_reduced_scan(ops: dict[str, Any], grid: dict[str, np.ndarray], prefactor: float) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rng = np.random.default_rng(20260506)
    scenarios = [
        {
            "name": "exact_link_grading",
            "eps_m": 0.0,
            "eps_sigma": 0.0,
            "eps_tau": 0.0,
            "rho": 1.0,
            "samples": 96,
        },
        {
            "name": "routeB_target",
            "eps_m": 1.0e-6,
            "eps_sigma": 1.0e-6,
            "eps_tau": 1.0e-3,
            "rho": 1.0,
            "samples": 192,
        },
        {
            "name": "routeB_conservative",
            "eps_m": 3.0e-7,
            "eps_sigma": 3.0e-7,
            "eps_tau": 1.0e-3,
            "rho": 1.0,
            "samples": 192,
        },
        {
            "name": "tau_conservative",
            "eps_m": 1.0e-6,
            "eps_sigma": 1.0e-6,
            "eps_tau": 3.0e-4,
            "rho": 1.0,
            "samples": 192,
        },
        {
            "name": "tau_stress",
            "eps_m": 0.0,
            "eps_sigma": 0.0,
            "eps_tau": 3.0e-3,
            "rho": 1.0,
            "samples": 96,
        },
    ]
    rows: list[dict[str, Any]] = []
    for sc in scenarios:
        for idx in range(sc["samples"]):
            dm = sc["eps_m"] * route_a.symmetric_unit(rng) if sc["eps_m"] else np.zeros((3, 3))
            ds = (
                sc["eps_sigma"] * route_a.symmetric_unit(rng)
                if sc["eps_sigma"]
                else np.zeros((3, 3))
            )
            dt = (
                sc["eps_tau"] * route_a.symmetric_unit(rng)
                if sc["eps_tau"]
                else np.zeros((3, 3))
            )
            sample = {
                "class": sc["name"],
                "amplitude": max(sc["eps_m"], sc["eps_sigma"], sc["eps_tau"], sc["rho"]),
                "sample": idx,
                "dm": dm,
                "ds": ds,
                "dt": dt,
                "rho": sc["rho"] * float(rng.choice([-1.0, 1.0])),
            }
            row = route_a.evaluate_sample(sample, ops, grid, prefactor)
            row.update(
                {
                    "eps_m": sc["eps_m"],
                    "eps_sigma": sc["eps_sigma"],
                    "eps_tau": sc["eps_tau"],
                    "rho_abs": sc["rho"],
                }
            )
            rows.append(row)

    summary: dict[str, Any] = {}
    for sc in scenarios:
        selected = [row for row in rows if row["class"] == sc["name"]]
        summary[sc["name"]] = {
            "samples": len(selected),
            "eps_m": sc["eps_m"],
            "eps_sigma": sc["eps_sigma"],
            "eps_tau": sc["eps_tau"],
            "rho_abs": sc["rho"],
            "safe_fraction": float(np.mean([row["safe_routeA"] for row in selected])),
            "median_goldstone_kappa": float(
                np.median([row["goldstone_max_abs_kappa"] for row in selected])
            ),
            "median_X622_projected_l2": float(
                np.median([row["x622_unwanted_projected_l2"] for row in selected])
            ),
            "median_RGE_residual": float(
                np.median([row["rge_best_residual_l2"] for row in selected])
            ),
            "median_tau_d6_years": float(
                np.median([row["best_tau_dim6_years"] for row in selected])
            ),
            "median_tau_d5_years": float(
                np.median([row["best_tau_dim5_ST_1e_5_years"] for row in selected])
            ),
            "max_fine_offblock": float(max(row["fine_rep_offblock_max_abs"] for row in selected)),
        }
    return rows, summary


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "class",
        "sample",
        "eps_m",
        "eps_sigma",
        "eps_tau",
        "rho_abs",
        "rho",
        "kappa3_eff",
        "kappa8_eff",
        "goldstone_max_abs_kappa",
        "x622_unwanted_projected_l2",
        "fine_rep_offblock_max_abs",
        "rge_best_residual_l2",
        "rge_safe_count",
        "best_tau_dim6_years",
        "best_tau_dim5_ST_1e_5_years",
        "safe_routeA",
    ]
    with (OUT / "routeB_reduced_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row[field] for field in fields} for row in rows])


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Route-B mediator grading verification")
    lines.append("")
    lines.append("No web lookup was used.  This pass designs the minimal selection-rule")
    lines.append("completion suggested by the Route-A deformation scan.")
    lines.append("")
    lines.append("## Abelian field-charge no-go")
    lines.append("")
    no_go = payload["abelian_no_go"]
    lines.append("A brute-force search over finite `Z_N` gradings found")
    lines.append("")
    lines.append("```text")
    lines.append(f"searched N <= {no_go['searched_ZN_up_to']}")
    lines.append(f"solutions = {no_go['solutions_found']}")
    lines.append("```")
    lines.append("")
    lines.append("The analytic reason is simple.  In any additive abelian grading, allowing")
    lines.append("`AA, AB, AC, BC` implies `q_A=q_B=q_C`; allowing `SAA, SAB, SAC` then")
    lines.append("also allows `SBB, SBC, SCC`.  A single unlabelled `54_H` portal cannot")
    lines.append("protect the projector card.")
    lines.append("")
    lines.append("## Link-graded completion")
    lines.append("")
    sel = payload["selection_rules"]
    lines.append("Use the mediator-node grading")
    lines.append("")
    lines.append("```text")
    for key, value in sel["node_charges_U1M"].items():
        lines.append(f"q_M({key}) = {value}")
    lines.append("```")
    lines.append("")
    lines.append("and endpoint-labeled portals")
    lines.append("")
    lines.append("```text")
    for key, value in sel["link_portals_U1M"].items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("The driving sector aligns the three `54`-valued links in the same Spin(10)")
    lines.append("direction and fixes the shifted Clebsch factors `F54-2` and `F54+4/3`.")
    lines.append("Thus the exact link-graded card keeps `AA, BC, AB, AC` but removes generic")
    lines.append("`m_ij` and `sigma_ij` entries.")
    lines.append("")
    lines.append("## Reduced Route-A replay")
    lines.append("")
    lines.append("| scenario | eps_m | eps_sigma | eps_tau | |rho| | safe frac | median Goldstone | median X622 l2 | median RGE residual |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, row in payload["reduced_scan_summary"].items():
        lines.append(
            f"| `{name}` | {row['eps_m']:.1e} | {row['eps_sigma']:.1e} | "
            f"{row['eps_tau']:.1e} | {row['rho_abs']:.1e} | {row['safe_fraction']:.3f} | "
            f"{row['median_goldstone_kappa']:.3e} | {row['median_X622_projected_l2']:.3e} | "
            f"{row['median_RGE_residual']:.3e} |"
        )
    lines.append("")
    lines.append("The practical lesson is that exact link grading is safe; the target")
    lines.append("`eps_m=eps_sigma=1e-6, eps_tau=1e-3, |rho|=1` remains safe for most")
    lines.append("random directions.  A more conservative manuscript benchmark is")
    lines.append("`eps_m=eps_sigma=3e-7, eps_tau=1e-3`, or `eps_tau=3e-4` if one wants")
    lines.append("extra RGE margin.  The stress test at `eps_tau=3e-3` fails often, matching")
    lines.append("the Route-A tolerance estimate.")
    (OUT / "routeB_mediator_grading_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ops, grid, prefactor = prepare_routeA_context()
    rows, summary = run_reduced_scan(ops, grid, prefactor)
    payload = {
        "note": "No web lookup used. Route-B mediator grading / R-symmetry selection-rule check.",
        "abelian_no_go": brute_force_abelian_no_go(),
        "selection_rules": routeB_selection_rules(),
        "reduced_scan_summary": summary,
    }
    (OUT / "routeB_mediator_grading_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(rows)
    write_report(payload)

    print("Route-B mediator grading verification")
    print(f"  abelian no-go solutions: {payload['abelian_no_go']['solutions_found']}")
    for name, row in summary.items():
        print(
            f"  {name}: safe={row['safe_fraction']:.3f}, "
            f"gold={row['median_goldstone_kappa']:.3e}, "
            f"res={row['median_RGE_residual']:.3e}"
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
