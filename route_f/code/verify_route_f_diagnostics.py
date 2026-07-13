#!/usr/bin/env python3
"""Diagnostic arithmetic card supporting the Route-F theory audit."""

from __future__ import annotations

import cmath
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"


def matmul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def max_abs_matrix(a: list[list[complex]]) -> float:
    return max(abs(x) for row in a for x in row)


def tangent_sections(genus: int) -> int:
    """Return h^0(T) from the standard curve cases used in the hand proof."""
    degree_tangent = 2 - 2 * genus
    if genus == 0:
        h1 = 0  # Serre dual h^0(K^2)=h^0(O(-4)).
        return degree_tangent + 1 - genus + h1
    if genus == 1:
        return 1  # T is the trivial line bundle.
    if degree_tangent < 0:
        return 0
    raise AssertionError("unhandled genus case")


def ad_invariance_residual_1d(
    structure_constant: list[list[list[float]]], bilinear: list[list[float]]
) -> float:
    """Compute max |f_ij^l B_lk + f_ik^l B_jl| in one dimension."""
    residual = 0.0
    for i in range(1):
        for j in range(1):
            for k in range(1):
                value = sum(
                    structure_constant[i][j][ell] * bilinear[ell][k]
                    + structure_constant[i][k][ell] * bilinear[j][ell]
                    for ell in range(1)
                )
                residual = max(residual, abs(value))
    return residual


def main() -> None:
    zeta = complex(0.1076472949, 0.0736514853)
    lam = cmath.sqrt(zeta)
    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    ktr = [
        [0j, 0j, -inv_sqrt3],
        [0j, inv_sqrt3, 0j],
        [-inv_sqrt3, 0j, 0j],
    ]
    ktr2 = matmul(ktr, ktr)
    target_ktr2 = [
        [1.0 / 3.0 if i == j else 0.0 for j in range(3)] for i in range(3)
    ]
    ktr2_residual = max_abs_matrix(
        [[ktr2[i][j] - target_ktr2[i][j] for j in range(3)] for i in range(3)]
    )

    genus_ladder = {
        "g=0": tangent_sections(0),
        "g=1": tangent_sections(1),
        "g=2": tangent_sections(2),
        "g=3": tangent_sections(3),
    }

    # One-dimensional abelian Lie algebra: all structure constants vanish.
    # B=[1] is symmetric, non-degenerate, and ad-invariance is vacuous.
    abelian_structure_constant = [[[0.0]]]
    abelian_bilinear = [[1.0]]
    abelian_form_det = abelian_bilinear[0][0]
    abelian_ad_invariance_residual = ad_invariance_residual_1d(
        abelian_structure_constant, abelian_bilinear
    )

    tree_delta_z = abs(zeta) / 3.0
    loop_delta_z = abs(zeta) / (16.0 * math.pi**2)
    delta_z_ratio = tree_delta_z / loop_delta_z
    z_n = 1.0 + tree_delta_z
    canonically_normalized_contact_abs = abs(zeta) / z_n

    reconstruction_error = abs(lam * lam - zeta)
    unitarity_em_ratio_g07 = math.sqrt(8.0 * math.pi) / 0.7

    results = {
        "status": "pass_diagnostic_arithmetic_card",
        "checks": {
            "genus_ladder": genus_ladder,
            "abelian_counterexample": {
                "dimension": 1,
                "bracket": "zero",
                "bilinear_form": [[1]],
                "determinant": abelian_form_det,
                "ad_invariance_residual": abelian_ad_invariance_residual,
                "nondegenerate_ad_invariant_form_exists": True,
                "implication": "H3 as ad-invariance alone does not exclude g=1",
            },
            "ktr": {
                "max_abs_K2_minus_I_over_3": ktr2_residual,
                "lambda_squared_minus_zeta": reconstruction_error,
            },
            "messenger_normalization": {
                "abs_zeta": abs(zeta),
                "abs_lambda_squared": abs(lam) ** 2,
                "tree_delta_Z_abs_zeta_over_3": tree_delta_z,
                "loop_delta_Z_abs_zeta_over_16pi2": loop_delta_z,
                "tree_to_loop_ratio": delta_z_ratio,
                "Z_N": z_n,
                "canonically_normalized_contact_abs_without_retuning": (
                    canonically_normalized_contact_abs
                ),
            },
            "longitudinal_unitarity": {
                "coupling_g": 0.7,
                "E_over_M_at_abs_Re_a0_equal_half": unitarity_em_ratio_g07,
            },
            "proton_lifetime_mass_sensitivity": {
                "tau_factor_for_Mx_times_2": 2.0**4,
                "tau_factor_for_Mx_times_10": 10.0**4,
            },
        },
    }

    assert ktr2_residual < 1e-15
    assert reconstruction_error < 1e-15
    assert abelian_form_det != 0.0
    assert abelian_ad_invariance_residual == 0.0
    assert genus_ladder == {"g=0": 3, "g=1": 1, "g=2": 0, "g=3": 0}
    assert tree_delta_z > 0.04
    assert delta_z_ratio > 50.0

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "route_f_diagnostics.json"
    md_path = OUTPUT / "route_f_diagnostics.md"
    json_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(
        "# Route-F Diagnostic Arithmetic Card\n\n"
        f"- status: **{results['status']}**\n"
        f"- genus ladder: `{genus_ladder}`\n"
        "- one-dimensional abelian counterexample: `det B = 1`, "
        "ad-invariance residual `0`\n"
        f"- `max|K_tr^2-I/3| = {ktr2_residual:.3e}`\n"
        f"- `|lambda^2-zeta| = {reconstruction_error:.3e}`\n"
        f"- tree Kahler `delta Z = {tree_delta_z:.10f}`\n"
        f"- quoted loop scale `delta Z = {loop_delta_z:.10e}`\n"
        f"- tree/loop ratio: `{delta_z_ratio:.3f}`\n"
        f"- `Z_N = {z_n:.10f}` and unretuned canonical contact magnitude "
        f"`{canonically_normalized_contact_abs:.10f}`\n"
        f"- uncancelled longitudinal unitarity scale at `g=0.7`: "
        f"`E/M = {unitarity_em_ratio_g07:.3f}`\n"
        "- proton-lifetime scaling: `Mx x2 -> tau x16`; "
        "`Mx x10 -> tau x10000`\n",
        encoding="utf-8",
    )
    print(f"route_f_diagnostics: pass_diagnostic_arithmetic_card ({json_path})")


if __name__ == "__main__":
    main()
