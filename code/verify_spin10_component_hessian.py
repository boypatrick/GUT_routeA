#!/usr/bin/env python3
"""Component-level Hessian check for the R=200 Spin(10) benchmark card.

No web lookup is used.  This script promotes the block-operator
``benchmark_R200.json`` into an explicit direct-sum Hessian over the relevant
component sectors:

* weak triplets Sigma_L and Sigma_R;
* the mixed X_(6,2,2) block;
* the color (15,1,1) branch split by the 210_H^3 cubic into
  octet + colored Goldstone pair + radial singlet.

It is not yet a full 210-component vacuum solver.  It is the next verifiable
component Hessian card: every block is an explicit finite mediator mass matrix,
and the color branch uses the already verified 210_H^3 Hessian pattern
(-1)^1, 0^6, 2^8.
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

import compute_210_cubic_matching as cubic210  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402


CARD = ROOT / "output" / "spin10_vacuum_alignment" / "benchmark_R200.json"
OUT = ROOT / "output" / "spin10_component_hessian"

PHYSICAL_BETA = {
    "H_C+H_Cbar": np.array([2.0 / 5.0, 0.0, 1.0], dtype=float),
    "Sigma_3": np.array([0.0, 2.0, 0.0], dtype=float),
    "Sigma_8_octet": np.array([0.0, 0.0, 3.0], dtype=float),
    "Sigma_R": np.array([6.0 / 5.0, 0.0, 0.0], dtype=float),
    "X_622": np.array([26.0 / 5.0, 6.0, 4.0], dtype=float),
    "colored_Goldstone_pair": np.array([8.0 / 5.0, 0.0, 1.0], dtype=float),
    "radial_singlet": np.array([0.0, 0.0, 0.0], dtype=float),
}


def mass_matrix(q: float, u: float, v: float, r_med: float) -> np.ndarray:
    return np.array([[q, u, v], [u, 0.0, r_med], [v, r_med, 0.0]], dtype=float)


def pick_light(evals: np.ndarray, target: float) -> float:
    return float(evals[np.argmin(np.abs(evals - target))])


def projector_l2(vec: np.ndarray) -> float:
    return float(np.linalg.norm(base.PROJECTOR @ vec))


def load_card() -> dict[str, Any]:
    return json.loads(CARD.read_text(encoding="utf-8"))


def fragment_by_name(card: dict[str, Any], fragment: str) -> dict[str, Any]:
    for row in card["physical_fragments"]:
        if row["fragment"] == fragment:
            return row
    raise KeyError(fragment)


def component_blocks(card: dict[str, Any]) -> list[dict[str, Any]]:
    r_med = float(card["R"])
    mg = float(card["rge_proton_replay"]["best"]["MG_GeV"])
    k8 = float(card["rge_proton_replay"]["best"]["kappa_8"])
    states = card["light_and_mediator_eigenvalues"]
    sigma8_matrix = np.array(states["Sigma8_block"]["matrix_over_MG"], dtype=float)
    q_octet = float(sigma8_matrix[0, 0])
    u_color = float(sigma8_matrix[0, 1])
    v_color = float(sigma8_matrix[0, 2])

    rows: list[dict[str, Any]] = []

    simple_specs = [
        ("Sigma_L", 3, float(states["Sigma_L"]["target_kappa"]), "physical weak-adjoint threshold"),
        ("Sigma_R", 3, 1.0, "lifted at MG"),
        ("X_622", 24, 1.0, "lifted at MG by P_X mediator sector"),
    ]
    for name, mult, target, role in simple_specs:
        matrix = np.array(states[name]["matrix_over_MG"], dtype=float)
        evals = np.linalg.eigvalsh(matrix)
        light = pick_light(evals, target)
        rows.append(
            {
                "block": name,
                "multiplicity": mult,
                "component_hessian_source": "F54/D210/P_X finite-mediator matrix",
                "role": role,
                "matrix_over_MG": matrix.tolist(),
                "target_signed_kappa": target,
                "light_signed_kappa": light,
                "light_abs_kappa": abs(light),
                "light_mass_GeV": abs(light) * mg,
                "benchmark_mass_GeV": target * mg,
                "relative_mass_error": abs(abs(light) - target) / max(abs(target), 1.0e-300),
                "heavy_signed_kappas": [float(x) for x in evals if abs(x - light) > 1.0e-7],
            }
        )

    color_specs = [
        ("Sigma_8_octet", 8, 2.0, k8, "physical color-octet threshold"),
        ("colored_Goldstone_pair", 6, 0.0, 0.0, "eaten chiral Goldstone pair"),
        ("radial_singlet", 1, -1.0, -0.5 * k8, "gauge-neutral radial singlet"),
    ]
    for name, mult, hessian_eigenvalue, target, role in color_specs:
        q = 0.5 * q_octet * hessian_eigenvalue
        matrix = mass_matrix(q, u_color, v_color, r_med)
        evals = np.linalg.eigvalsh(matrix)
        light = pick_light(evals, target)
        benchmark_mass = abs(target) * mg
        if name == "radial_singlet":
            benchmark_mass = 0.5 * float(fragment_by_name(card, "Sigma_8")["mass_GeV"])
        if benchmark_mass > 0.0:
            relative_mass_error = abs(abs(light) * mg - benchmark_mass) / benchmark_mass
        else:
            # For exact Goldstone targets the meaningful error is dimensionless
            # absolute light kappa, not a relative mass error against zero.
            relative_mass_error = abs(light)
        rows.append(
            {
                "block": name,
                "multiplicity": mult,
                "component_hessian_source": "210_H^3 color Hessian eigenvalue times q_octet/2 plus finite mediator mixing",
                "role": role,
                "I3_hessian_eigenvalue": hessian_eigenvalue,
                "top_left_q_over_MG": q,
                "matrix_over_MG": matrix.tolist(),
                "target_signed_kappa": target,
                "light_signed_kappa": light,
                "light_abs_kappa": abs(light),
                "light_mass_GeV": abs(light) * mg,
                "benchmark_mass_GeV": benchmark_mass,
                "relative_mass_error": relative_mass_error,
                "heavy_signed_kappas": [float(x) for x in evals if abs(x - light) > 1.0e-7],
            }
        )
    return rows


def threshold_audit(card: dict[str, Any], blocks: list[dict[str, Any]]) -> dict[str, Any]:
    mg = float(card["rge_proton_replay"]["best"]["MG_GeV"])
    logs = np.array(
        [
            float(card["rge_proton_replay"]["best"]["log_HC"]),
            float(card["rge_proton_replay"]["best"]["log_Sigma3"]),
            float(card["rge_proton_replay"]["best"]["log_Sigma8"]),
        ],
        dtype=float,
    )
    chiral = base.HEAVY_BASIS @ logs / (2.0 * math.pi)
    mediator = np.array(card["mediator_heavy_threshold"]["delta_full"], dtype=float)
    total = chiral + mediator
    alpha_inv = np.array(card["matching_threshold_check"]["alpha_inv_raw"], dtype=float)
    residual = base.PROJECTOR @ (alpha_inv - total)

    sigma8_log = float(card["rge_proton_replay"]["best"]["log_Sigma8"])
    colored_if_unlocked = PHYSICAL_BETA["colored_Goldstone_pair"] * sigma8_log / (
        2.0 * math.pi
    )
    x622_log_if_unlifted = math.log(mg / float(fragment_by_name(card, "Sigma_8")["mass_GeV"]))
    x622_if_unlifted = PHYSICAL_BETA["X_622"] * x622_log_if_unlifted / (2.0 * math.pi)
    active_component_threshold = (
        PHYSICAL_BETA["Sigma_3"] * float(card["rge_proton_replay"]["best"]["log_Sigma3"])
        + PHYSICAL_BETA["Sigma_8_octet"] * sigma8_log
    ) / (2.0 * math.pi)
    radial_threshold = PHYSICAL_BETA["radial_singlet"] * float(
        next(row for row in blocks if row["block"] == "radial_singlet")["light_abs_kappa"]
    )
    return {
        "active_component_threshold_without_HC": active_component_threshold.tolist(),
        "HC_Sigma3_Sigma8_chiral_threshold": chiral.tolist(),
        "mediator_threshold": mediator.tolist(),
        "total_threshold": total.tolist(),
        "projected_matching_residual": residual.tolist(),
        "projected_matching_residual_l2": float(np.linalg.norm(residual)),
        "colored_pair_if_unlocked_at_Msigma8_projected_l2": projector_l2(colored_if_unlocked),
        "X_622_if_unlifted_at_Msigma8_projected_l2": projector_l2(x622_if_unlifted),
        "radial_singlet_threshold_l2": projector_l2(radial_threshold),
    }


def cubic_hessian_audit() -> dict[str, Any]:
    canonical = cubic210.canonical_formula_check()
    ratio = float(canonical["canonical_block_rows"][0]["ratio_I3_over_pfaffian"])
    hessian = cubic210.superpotential_hessian_check(ratio)
    return {
        "canonical_formula": canonical,
        "cubic_ratio_I3_over_pfaffian": ratio,
        "hessian_match": hessian,
        "interpreted_component_multiplicities": {
            "radial_singlet": 1,
            "colored_Goldstone_pair": 6,
            "Sigma_8_octet": 8,
        },
    }


def write_csv(blocks: list[dict[str, Any]]) -> None:
    fields = [
        "block",
        "multiplicity",
        "role",
        "target_signed_kappa",
        "light_signed_kappa",
        "light_abs_kappa",
        "light_mass_GeV",
        "benchmark_mass_GeV",
        "relative_mass_error",
    ]
    with (OUT / "spin10_component_hessian_blocks.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in blocks:
            writer.writerow({key: row.get(key, "") for key in fields})


def write_report(payload: dict[str, Any]) -> None:
    blocks = payload["component_blocks"]
    audit = payload["threshold_audit"]
    card = payload["benchmark_card"]
    lines: list[str] = []
    lines.append("# Spin(10) component Hessian check")
    lines.append("")
    lines.append("No web lookup was used.  This is a component-Hessian upgrade of")
    lines.append("`benchmark_R200.json`.  The calculation is still restricted to the")
    lines.append("relevant subspace, but every displayed component is represented by an")
    lines.append("explicit finite-mediator mass matrix.")
    lines.append("")
    lines.append("## Component superpotential Hessian")
    lines.append("")
    lines.append("For non-color blocks the Hessian is the finite-mediator matrix")
    lines.append("")
    lines.append("```text")
    lines.append("M_fd/MG = [[mu+a54 f+b210 d, eta(f-2), eta(f+4/3)],")
    lines.append("          [eta(f-2), 0, R],")
    lines.append("          [eta(f+4/3), R, 0]].")
    lines.append("```")
    lines.append("")
    lines.append("For the color `(15,1,1)` branch the explicit `210_H^3` Hessian supplies")
    lines.append("the component eigenvalues")
    lines.append("")
    lines.append("```text")
    lines.append("radial singlet: -1")
    lines.append("colored Goldstone pair: 0^6")
    lines.append("color octet: +2^8")
    lines.append("```")
    lines.append("")
    hess = payload["cubic_210_audit"]["hessian_match"]
    lines.append("The direct numerical 210 cubic check gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"gradient norm = {hess['gradient_norm_at_vacuum']:.3e}")
    lines.append(f"clusters = {hess['eigenvalue_clusters']}")
    lines.append(f"max deviation = {hess['max_abs_deviation_from_expected']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Block-by-block match")
    lines.append("")
    lines.append("| block | multiplicity | light kappa | benchmark mass [GeV] | component mass [GeV] | relative error |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in blocks:
        lines.append(
            f"| `{row['block']}` | {row['multiplicity']} | {row['light_signed_kappa']:.9g} | "
            f"{row['benchmark_mass_GeV']:.6e} | {row['light_mass_GeV']:.6e} | "
            f"{row['relative_mass_error']:.3e} |"
        )
    lines.append("")
    lines.append("The crucial checks are:")
    lines.append("")
    lines.append("```text")
    x622 = next(row for row in blocks if row["block"] == "X_622")
    gold = next(row for row in blocks if row["block"] == "colored_Goldstone_pair")
    rad = next(row for row in blocks if row["block"] == "radial_singlet")
    lines.append(f"X_622 light kappa = {x622['light_signed_kappa']:.12g} = 1")
    lines.append(f"colored Goldstone light kappa = {gold['light_signed_kappa']:.3e}")
    lines.append(f"radial mass / (M_Sigma8/2) - 1 = {rad['relative_mass_error']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Threshold audit")
    lines.append("")
    lines.append("The component Hessian keeps only `H_C`, `Sigma_3`, and the octet part of")
    lines.append("`Sigma_8` active in the one-loop chiral threshold.  The matching residual is")
    lines.append("")
    lines.append("```text")
    lines.append(f"||P(alpha^-1 - Delta_med - Delta_chiral)|| = {audit['projected_matching_residual_l2']:.3e}")
    lines.append(f"colored pair if unlocked at M_Sigma8: {audit['colored_pair_if_unlocked_at_Msigma8_projected_l2']:.6e}")
    lines.append(f"X_622 if unlifted at M_Sigma8: {audit['X_622_if_unlifted_at_Msigma8_projected_l2']:.6e}")
    lines.append(f"radial singlet one-loop threshold: {audit['radial_singlet_threshold_l2']:.1e}")
    lines.append("```")
    lines.append("")
    lines.append("Thus the component Hessian reproduces the R=200 card on the relevant")
    lines.append("subspace.  The remaining paper-level task is now sharper: derive this")
    lines.append("component Hessian from a single written `W_54+W_210+W_med+W_PS`")
    lines.append("superpotential and check full F-flatness/D-flatness, instead of treating")
    lines.append("the color cubic and projector sector as aligned blocks.")
    lines.append("")
    lines.append("## Source card")
    lines.append("")
    lines.append("```text")
    lines.append(str(CARD))
    lines.append(f"R = {card['R']:.0f}")
    lines.append(f"mu = {card['vev_parameters_and_couplings']['mu']:.9f}")
    lines.append(f"a54 = {card['vev_parameters_and_couplings']['a54']:.9f}")
    lines.append(f"b210 = {card['vev_parameters_and_couplings']['b210']:.9f}")
    lines.append(f"eta = {card['vev_parameters_and_couplings']['eta']:.9f}")
    lines.append("```")
    (OUT / "spin10_component_hessian_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = load_card()
    blocks = component_blocks(card)
    payload = {
        "note": "No web lookup used. Component Hessian check of benchmark_R200 on the relevant Spin(10)/Pati-Salam subspace.",
        "benchmark_card": {
            "source": str(CARD),
            "R": float(card["R"]),
            "vev_parameters_and_couplings": card["vev_parameters_and_couplings"],
        },
        "cubic_210_audit": cubic_hessian_audit(),
        "component_blocks": blocks,
        "threshold_audit": threshold_audit(card, blocks),
        "passes": {
            "X_622_at_MG": abs(next(row for row in blocks if row["block"] == "X_622")["light_signed_kappa"] - 1.0)
            < 1.0e-12,
            "colored_goldstone_zero": abs(
                next(row for row in blocks if row["block"] == "colored_Goldstone_pair")["light_signed_kappa"]
            )
            < 1.0e-12,
            "radial_matches_Msigma8_over_2": next(
                row for row in blocks if row["block"] == "radial_singlet"
            )["relative_mass_error"]
            < 1.0e-7,
            "matching_residual_roundoff": threshold_audit(card, blocks)[
                "projected_matching_residual_l2"
            ]
            < 1.0e-12,
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "spin10_component_hessian_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(blocks)
    write_report(payload)

    print("Spin(10) component Hessian check")
    for row in blocks:
        print(
            f"  {row['block']}: mult={row['multiplicity']}, "
            f"light={row['light_signed_kappa']:.12g}, rel_err={row['relative_mass_error']:.3e}"
        )
    print(f"  matching residual: {payload['threshold_audit']['projected_matching_residual_l2']:.3e}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
