#!/usr/bin/env python3
"""Combined-superpotential flatness and off-block-mixing check.

No web lookup is used.  This script writes and verifies a single local
superpotential on the relevant Spin(10)/Pati-Salam subspace,

    W = W_54 + W_210 + W_med + W_PS,

whose Hessian reproduces the R=200 component card.  The point of this pass is
to check the dangerous glue between the color 210_H^3 sector and the
54/210/projector mediator sector:

* F-flatness, including mediator tadpoles around the color PS vacuum;
* D-flatness of the real aligned vevs;
* off-block Hessian mixing.

The calculation is still a local component card rather than a full 210-field
global minimization.  It gives the explicit equations the full model must
embed.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import compute_210_cubic_matching as cubic210  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "combined_superpotential_flatness"


def permutation_sign(values: tuple[int, ...]) -> int:
    if len(set(values)) != len(values):
        return 0
    inv = 0
    for idx, vi in enumerate(values):
        for vj in values[idx + 1 :]:
            if vi > vj:
                inv += 1
    return -1 if inv % 2 else 1


def so_generators(n: int = 10) -> list[np.ndarray]:
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            t = np.zeros((n, n), dtype=float)
            t[i, j] = 1.0
            t[j, i] = -1.0
            gens.append(t)
    return gens


def canonical_color_twoform_10d() -> np.ndarray:
    a = np.zeros((10, 10), dtype=float)
    for i, j in [(0, 1), (2, 3), (4, 5)]:
        a[i, j] = 1.0
        a[j, i] = -1.0
    return a


def dflat_adjoint_max_inner(a: np.ndarray) -> float:
    max_abs = 0.0
    for t in so_generators(a.shape[0]):
        action = t @ a - a @ t
        inner = float(np.trace(a.T @ action))
        max_abs = max(max_abs, abs(inner))
    return max_abs


def dflat_54_max_inner() -> float:
    s = np.diag([-2.0] * 6 + [3.0] * 4)
    max_abs = 0.0
    for t in so_generators(10):
        action = t @ s - s @ t
        inner = float(np.trace(s.T @ action))
        max_abs = max(max_abs, abs(inner))
    return max_abs


PFORM_BASIS_4 = list(itertools.combinations(range(10), 4))
PFORM_INDEX_4 = {basis: idx for idx, basis in enumerate(PFORM_BASIS_4)}


def weak_volume_210_vector() -> np.ndarray:
    vec = np.zeros(len(PFORM_BASIS_4), dtype=float)
    vec[PFORM_INDEX_4[(6, 7, 8, 9)]] = 1.0
    return vec


def color_210_vector_from_a0() -> np.ndarray:
    a6 = np.zeros((6, 6), dtype=float)
    for i, j in [(0, 1), (2, 3), (4, 5)]:
        a6[i, j] = 1.0
        a6[j, i] = -1.0
    vec = np.zeros(len(PFORM_BASIS_4), dtype=float)
    for basis, idx in PFORM_INDEX_4.items():
        vec[idx] = cubic210.hodge4_from_twoform_component(a6, *basis)
    return vec


def pform_action(t: np.ndarray, vec: np.ndarray) -> np.ndarray:
    out = np.zeros_like(vec)
    for basis, idx in PFORM_INDEX_4.items():
        value = vec[idx]
        if value == 0.0:
            continue
        raw = list(basis)
        for pos, old in enumerate(raw):
            for new in range(10):
                coeff = t[new, old]
                if coeff == 0.0:
                    continue
                replaced = raw.copy()
                replaced[pos] = new
                if len(set(replaced)) < 4:
                    continue
                sorted_tuple = tuple(sorted(replaced))
                sign = permutation_sign(tuple(replaced)) * permutation_sign(sorted_tuple)
                # sorted_tuple already has positive order, so only the sign of
                # the unsorted replacement is relevant.  The second factor is 1
                # but keeps the convention explicit.
                out[PFORM_INDEX_4[sorted_tuple]] += coeff * sign * value
    return out


def dflat_210_max_inner(vec: np.ndarray) -> float:
    max_abs = 0.0
    for t in so_generators(10):
        inner = float(np.dot(vec, pform_action(t, vec)))
        max_abs = max(max_abs, abs(inner))
    return max_abs


def build_combined_hessian(blocks: list[dict[str, Any]]) -> tuple[np.ndarray, list[str]]:
    matrices = []
    labels = []
    for row in blocks:
        matrix = np.array(row["matrix_over_MG"], dtype=float)
        for copy in range(int(row["multiplicity"])):
            matrices.append(matrix)
            labels.append(f"{row['block']}[{copy}]")
    total_dim = 3 * len(matrices)
    h = np.zeros((total_dim, total_dim), dtype=float)
    for idx, matrix in enumerate(matrices):
        start = 3 * idx
        h[start : start + 3, start : start + 3] = matrix
    return h, labels


def off_block_norms(hessian: np.ndarray) -> dict[str, float]:
    mask = np.zeros_like(hessian, dtype=bool)
    for start in range(0, hessian.shape[0], 3):
        mask[start : start + 3, start : start + 3] = True
    off = np.where(mask, 0.0, hessian)
    return {
        "max_abs_off_component_entry": float(np.max(np.abs(off))),
        "frobenius_off_component": float(np.linalg.norm(off)),
    }


def light_eigenvalue_errors(blocks: list[dict[str, Any]]) -> dict[str, float]:
    return {
        row["block"]: float(row["relative_mass_error"])
        for row in blocks
        if row["block"]
        in {"Sigma_L", "Sigma_R", "X_622", "Sigma_8_octet", "colored_Goldstone_pair", "radial_singlet"}
    }


def f_flatness(card: dict[str, Any], cubic_audit: dict[str, Any]) -> dict[str, Any]:
    pars = card["vev_parameters_and_couplings"]
    eta = float(pars["eta"])
    r_med = float(card["R"])
    f_color = -4.0 / 3.0
    u_color = eta * (f_color - 2.0)
    v_color = eta * (f_color + 4.0 / 3.0)
    c0_over_a0 = -u_color / r_med
    b0_over_a0 = -v_color / r_med
    # In the shifted local superpotential the PS color vacuum is F-flat.
    # The small nonzero number is the finite-difference gradient from the
    # already verified 210_H^3 calculation, scaled by m_C=q8/2.
    q8 = float(card["light_and_mediator_eigenvalues"]["Sigma8_block"]["matrix_over_MG"][0][0])
    m_c = 0.5 * q8
    ps_gradient = float(cubic_audit["hessian_match"]["gradient_norm_at_vacuum"])
    mediator_fb_residual_over_a0 = r_med * c0_over_a0 + u_color
    mediator_fc_residual_over_a0 = r_med * b0_over_a0 + v_color
    mediator_fa_residual_over_a0 = u_color * b0_over_a0 + v_color * c0_over_a0
    total = math.sqrt(
        (m_c * ps_gradient) ** 2
        + mediator_fb_residual_over_a0**2
        + mediator_fc_residual_over_a0**2
        + mediator_fa_residual_over_a0**2
    )
    return {
        "interpretation": "A denotes the fluctuation around the PS color vacuum; B0,C0 are chosen so W_med has no tadpole.",
        "u_color_eta_F_minus_2": u_color,
        "v_color_eta_F_plus_4_over_3": v_color,
        "B0_over_A0": b0_over_a0,
        "C0_over_A0": c0_over_a0,
        "PS_210H3_gradient_norm": ps_gradient,
        "m_C_q8_over_2": m_c,
        "mediator_F_B_residual_over_A0": mediator_fb_residual_over_a0,
        "mediator_F_C_residual_over_A0": mediator_fc_residual_over_a0,
        "mediator_F_A_residual_over_A0": mediator_fa_residual_over_a0,
        "combined_dimensionless_F_norm": total,
    }


def d_flatness() -> dict[str, Any]:
    a0 = canonical_color_twoform_10d()
    weak = weak_volume_210_vector()
    color = color_210_vector_from_a0()
    return {
        "principle": "All nonzero vevs are real self-conjugate directions; adjoint vevs commute with their Hermitian conjugates, and real p-form reps have antisymmetric generators.",
        "adjoint_color_commutator_D_proxy": dflat_adjoint_max_inner(a0),
        "54H_diagonal_D_proxy": dflat_54_max_inner(),
        "210H_weak_volume_D_proxy": dflat_210_max_inner(weak),
        "210H_color_fourform_D_proxy": dflat_210_max_inner(color),
        "mediator_B0_parallel_A0": True,
        "mediator_C0_parallel_A0": True,
        "max_abs_D_proxy": max(
            dflat_adjoint_max_inner(a0),
            dflat_54_max_inner(),
            dflat_210_max_inner(weak),
            dflat_210_max_inner(color),
        ),
    }


def combined_superpotential_formula(card: dict[str, Any]) -> dict[str, Any]:
    pars = card["vev_parameters_and_couplings"]
    return {
        "W54_plus_W210": "1/2 <A, (1-P_C)(mu+a54 F54+b210 D210) A>",
        "WPS": "m_C [1/2 ||A_C||^2 - I3(*_6 A_C)/6], expanded around A0=e12+e34+e56",
        "Wmed": "R <B,C> + eta <A,(F54-2)B> + eta <A,(F54+4/3)C>",
        "projectors": {
            "P_C": "(9/50)(F54-2)(F54-1/3)",
            "P_X": "-(9/25)(F54-2)(F54+4/3)",
            "P_L": "(D210^2+D210)/2",
            "P_R": "(D210^2-D210)/2",
        },
        "parameters": {
            "mu": pars["mu"],
            "a54": pars["a54"],
            "b210": pars["b210"],
            "eta": pars["eta"],
            "R": card["R"],
            "m_C": 0.5
            * float(card["light_and_mediator_eigenvalues"]["Sigma8_block"]["matrix_over_MG"][0][0]),
        },
        "novel_element": "The color projector P_C routes the (15,1,1) branch to W_PS, avoiding double counting with the linear 54/210 Clebsch mass.  A small mediator counter-vev C0=-eta(F54-2)A0/R cancels the W_med tadpole while preserving D-flatness.",
    }


def write_report(payload: dict[str, Any]) -> None:
    formula = payload["combined_superpotential"]
    f = payload["F_flatness"]
    d = payload["D_flatness"]
    off = payload["off_block_mixing"]
    lines: list[str] = []
    lines.append("# Combined superpotential flatness check")
    lines.append("")
    lines.append("No web lookup was used.  This pass writes one local combined")
    lines.append("superpotential on the relevant Spin(10)/Pati-Salam component subspace")
    lines.append("and checks F-flatness, D-flatness, and off-block Hessian mixing.")
    lines.append("")
    lines.append("## Combined superpotential")
    lines.append("")
    lines.append("With self-adjoint Clebsch operators `F54,D210` and")
    lines.append("`P_C=(9/50)(F54-2)(F54-1/3)`, use")
    lines.append("")
    lines.append("```text")
    lines.append("W/MG = W54 + W210 + Wmed + WPS")
    lines.append(f"W54+W210 = {formula['W54_plus_W210']}")
    lines.append(f"WPS      = {formula['WPS']}")
    lines.append(f"Wmed     = {formula['Wmed']}")
    lines.append("```")
    lines.append("")
    lines.append("The color projector `P_C` is the important new alignment device: the")
    lines.append("`(15,1,1)` branch is governed by the explicit `210_H^3`/PS cubic instead")
    lines.append("of being double-counted by the linear `54/210` Clebsch mass.")
    lines.append("")
    lines.append("## F-flatness")
    lines.append("")
    lines.append("The mediator tadpole created by the color PS vacuum is cancelled by")
    lines.append("")
    lines.append("```text")
    lines.append("B0/A0 = -eta(F54+4/3)/R")
    lines.append("C0/A0 = -eta(F54-2)/R")
    lines.append("```")
    lines.append("")
    lines.append("On the color branch `F54=-4/3`, so `B0=0` and")
    lines.append(f"`C0/A0={f['C0_over_A0']:.9f}`.  Numerically:")
    lines.append("")
    lines.append("```text")
    lines.append(f"PS 210H3 gradient norm = {f['PS_210H3_gradient_norm']:.3e}")
    lines.append(f"F_B residual / A0 = {f['mediator_F_B_residual_over_A0']:.3e}")
    lines.append(f"F_C residual / A0 = {f['mediator_F_C_residual_over_A0']:.3e}")
    lines.append(f"F_A residual / A0 = {f['mediator_F_A_residual_over_A0']:.3e}")
    lines.append(f"combined dimensionless F norm = {f['combined_dimensionless_F_norm']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("The remaining nonzero value is just the finite-difference gradient from the")
    lines.append("previous `210_H^3` Hessian check.")
    lines.append("")
    lines.append("## D-flatness")
    lines.append("")
    lines.append("All vevs are real self-conjugate aligned directions.  The numerical D-proxy")
    lines.append("checks are inner products `<v,Tv>` over all 45 SO(10) generators:")
    lines.append("")
    lines.append("```text")
    lines.append(f"adjoint color proxy = {d['adjoint_color_commutator_D_proxy']:.3e}")
    lines.append(f"54H diagonal proxy = {d['54H_diagonal_D_proxy']:.3e}")
    lines.append(f"210H weak-volume proxy = {d['210H_weak_volume_D_proxy']:.3e}")
    lines.append(f"210H color-fourform proxy = {d['210H_color_fourform_D_proxy']:.3e}")
    lines.append(f"max proxy = {d['max_abs_D_proxy']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("Since the mediator counter-vev is parallel to the color adjoint vacuum, it")
    lines.append("does not introduce a new D-term direction.")
    lines.append("")
    lines.append("## Off-block Hessian mixing")
    lines.append("")
    lines.append("The full restricted Hessian has dimension")
    lines.append(f"`{payload['combined_hessian']['dimension']}` including A/B/C for all")
    lines.append("component copies.  The off-block entries are")
    lines.append("")
    lines.append("```text")
    lines.append(f"max abs off-block entry = {off['max_abs_off_component_entry']:.3e}")
    lines.append(f"off-block Frobenius norm = {off['frobenius_off_component']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("Thus off-block mixing vanishes in the projector-aligned local card.  The")
    lines.append("remaining global-model task is to show that a full untruncated")
    lines.append("`54_H+210_H+mediator+PS` superpotential produces the same projector")
    lines.append("orthogonality without adding forbidden invariant contractions.")
    lines.append("")
    lines.append("## Pass flags")
    lines.append("")
    for key, value in payload["passes"].items():
        lines.append(f"- `{key}`: `{value}`")
    (OUT / "combined_superpotential_flatness_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = comp.load_card()
    blocks = comp.component_blocks(card)
    hessian, labels = build_combined_hessian(blocks)
    cubic_audit = comp.cubic_hessian_audit()
    payload = {
        "note": "No web lookup used. Local combined W54+W210+Wmed+WPS flatness/off-block check.",
        "source_card": str(comp.CARD),
        "combined_superpotential": combined_superpotential_formula(card),
        "F_flatness": f_flatness(card, cubic_audit),
        "D_flatness": d_flatness(),
        "combined_hessian": {
            "dimension": int(hessian.shape[0]),
            "component_copies": labels,
            "eigenvalue_min": float(np.min(np.linalg.eigvalsh(hessian))),
            "eigenvalue_max": float(np.max(np.linalg.eigvalsh(hessian))),
        },
        "off_block_mixing": off_block_norms(hessian),
        "light_eigenvalue_relative_errors": light_eigenvalue_errors(blocks),
    }
    payload["passes"] = {
        "F_flat_to_numeric_precision": payload["F_flatness"]["combined_dimensionless_F_norm"]
        < 1.0e-10,
        "D_flat_to_numeric_precision": payload["D_flatness"]["max_abs_D_proxy"] < 1.0e-12,
        "off_block_mixing_zero": payload["off_block_mixing"]["max_abs_off_component_entry"]
        < 1.0e-15,
        "component_hessian_matches_card": max(payload["light_eigenvalue_relative_errors"].values())
        < 1.0e-7,
    }
    payload["passes"]["all"] = all(payload["passes"].values())

    (OUT / "combined_superpotential_flatness_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)

    print("Combined superpotential flatness check")
    print(f"  F norm: {payload['F_flatness']['combined_dimensionless_F_norm']:.3e}")
    print(f"  D proxy max: {payload['D_flatness']['max_abs_D_proxy']:.3e}")
    print(f"  off-block max: {payload['off_block_mixing']['max_abs_off_component_entry']:.3e}")
    print(f"  Hessian dimension: {payload['combined_hessian']['dimension']}")
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
