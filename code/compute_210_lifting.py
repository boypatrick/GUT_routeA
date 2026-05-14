#!/usr/bin/env python3
"""Direct Spin(10) 45-45-210 Clebsch and lifting-sector test.

The 210 of Spin(10) is represented as an SO(10) four-form Phi_ijkl.  The
Pati-Salam preserving D-parity-odd singlet is the SO(4) volume form on the
weak block.  Contracting it with two adjoint two-forms gives the Hodge-star
operator on SO(4):

    W_210 = (zeta/8) Phi_ijkl A_ij A_kl.

Thus the weak adjoint splits into self-dual and anti-self-dual triplets, while
the color adjoint and mixed (6,2,2) blocks are annihilated.  This script builds
the full 45-by-45 matrix from the four-form, diagonalizes it, and tests a
minimal 54+210 spectral-projector lifting sector against the existing
threshold/proton-safe point.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "output" / "yukawa_superpotential_rge" / "yukawa_superpotential_rge_summary.json"
OUT = ROOT / "output" / "clebsch_210"

PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0
B_SIGMA_R = np.array([6.0 / 5.0, 0.0, 0.0], dtype=float)
B_X_622 = np.array([26.0 / 5.0, 6.0, 4.0], dtype=float)

# Normalized Clebsch operators on the blocks relevant for threshold matching.
C54 = {
    "Sigma_L_(1,3,1)": 2.0,
    "Sigma_R_(1,1,3)": 2.0,
    "Sigma8_(15,1,1)": -4.0 / 3.0,
    "X_(6,2,2)": 1.0 / 3.0,
}
C210 = {
    "Sigma_L_(1,3,1)": 1.0,
    "Sigma_R_(1,1,3)": -1.0,
    "Sigma8_(15,1,1)": 0.0,
    "X_(6,2,2)": 0.0,
}


def pair_basis() -> list[tuple[int, int]]:
    return [(i, j) for i in range(10) for j in range(i + 1, 10)]


def permutation_sign(values: tuple[int, ...]) -> int:
    inv = 0
    for i, vi in enumerate(values):
        for vj in values[i + 1 :]:
            if vi > vj:
                inv += 1
    return -1 if inv % 2 else 1


def phi_weak_volume(indices: tuple[int, int, int, int]) -> float:
    """SO(4) volume form on zero-based indices 6,7,8,9.

    These are the one-based SO(10) labels 7,8,9,10 used in the paper.
    """

    if len(set(indices)) != 4 or set(indices) != {6, 7, 8, 9}:
        return 0.0
    return float(permutation_sign(indices))


def matrix_210() -> tuple[list[tuple[int, int]], np.ndarray]:
    basis = pair_basis()
    h = np.zeros((len(basis), len(basis)), dtype=float)
    for a, (i, j) in enumerate(basis):
        for b, (k, l) in enumerate(basis):
            h[a, b] = phi_weak_volume((i, j, k, l))
    return basis, h


def classify_pair(i: int, j: int) -> str:
    if i < 6 and j < 6:
        return "Sigma8_(15,1,1)"
    if i >= 6 and j >= 6:
        return "weak_SO4"
    return "X_(6,2,2)"


def direct_210_spectrum() -> dict[str, object]:
    basis, h = matrix_210()
    evals, evecs = np.linalg.eigh(h)
    rounded = np.round(evals, 12)
    unique = {}
    for val in rounded:
        unique[str(float(val))] = int(np.sum(np.isclose(rounded, val)))

    rows = []
    for n, val in enumerate(evals):
        support = [
            {
                "pair_1based": [basis[i][0] + 1, basis[i][1] + 1],
                "coefficient": float(evecs[i, n]),
            }
            for i in np.where(np.abs(evecs[:, n]) > 1.0e-9)[0]
        ]
        if abs(val) > 1.0e-9:
            rows.append({"eigenvalue": float(val), "support": support})

    block_zero_counts = {"Sigma8_(15,1,1)": 0, "X_(6,2,2)": 0}
    for i, j in basis:
        block = classify_pair(i, j)
        if block in block_zero_counts:
            block_zero_counts[block] += 1

    return {
        "invariant": "W_210=(zeta/8) Phi_ijkl A_ij A_kl",
        "vev": "Phi_{7,8,9,10}=+1 on the weak SO(4) block, using one-based SO(10) labels",
        "eigenvalue_multiplicities": unique,
        "nonzero_eigenvectors": rows,
        "block_clebsch": {
            "Sigma_L_(1,3,1)": 1.0,
            "Sigma_R_(1,1,3)": -1.0,
            "Sigma8_(15,1,1)": 0.0,
            "X_(6,2,2)": 0.0,
        },
        "zero_block_counts": block_zero_counts,
    }


def solve_coefficients(k3: float, k8: float, kx: float) -> dict[str, float]:
    """Solve k_r = mu + a C54_r + b C210_r from k3,k8,kx."""

    a54 = (kx - k8) / (C54["X_(6,2,2)"] - C54["Sigma8_(15,1,1)"])
    mu = kx - a54 * C54["X_(6,2,2)"]
    b210 = k3 - mu - a54 * C54["Sigma_L_(1,3,1)"]
    k_r = mu + a54 * C54["Sigma_R_(1,1,3)"] + b210 * C210["Sigma_R_(1,1,3)"]
    return {
        "mu": mu,
        "a54": a54,
        "b210": b210,
        "kappa_Sigma3": k3,
        "kappa_Sigma8": k8,
        "kappa_X_622": kx,
        "kappa_SigmaR": k_r,
    }


def extra_residual(coeffs: dict[str, float]) -> dict[str, object]:
    kr = coeffs["kappa_SigmaR"]
    kx = coeffs["kappa_X_622"]
    if kr <= 0.0 or kx <= 0.0:
        return {"valid": False, "residual_l2": float("inf")}
    log_r = -math.log(kr)
    log_x = -math.log(kx)
    vec = PROJECTOR @ (B_SIGMA_R * log_r + B_X_622 * log_x) / (2.0 * math.pi)
    return {
        "valid": True,
        "log_SigmaR": log_r,
        "log_X_622": log_x,
        "residual_vector": vec.tolist(),
        "residual_l2": float(np.linalg.norm(vec)),
    }


def scan_kx(k3: float, k8: float) -> dict[str, object]:
    lower = (2.0 * k8 + k3) / 4.0 + 1.0e-9
    upper = 3.0
    best = None
    for kx in np.linspace(lower, upper, 300_000):
        coeffs = solve_coefficients(k3, k8, float(kx))
        resid = extra_residual(coeffs)
        if not resid["valid"]:
            continue
        item = {"kappa_X_trial": float(kx), "coefficients": coeffs, "extra_threshold": resid}
        if best is None or resid["residual_l2"] < best["extra_threshold"]["residual_l2"]:
            best = item
    if best is None:
        raise RuntimeError("No positive 54+210 point found")
    return best


def scenario_payload(best: dict[str, float]) -> dict[str, object]:
    mg = float(best["MG_GeV"])
    k3 = float(best["M_Sigma3_GeV"]) / mg
    k8 = float(best["M_Sigma8_GeV"]) / mg

    kx_for_r_at_mg = (1.0 + 2.0 * k8 + k3) / 4.0
    r_at_mg = solve_coefficients(k3, k8, kx_for_r_at_mg)
    r_at_mg_resid = extra_residual(r_at_mg)
    delta_x_lift = 1.0 - r_at_mg["kappa_X_622"]

    lifted = dict(r_at_mg)
    lifted["kappa_X_622_before_lift"] = lifted["kappa_X_622"]
    lifted["delta_X_projector"] = delta_x_lift
    lifted["kappa_X_622"] = 1.0
    lifted["kappa_SigmaR"] = 1.0
    lifted_resid = extra_residual(lifted)

    kx_at_mg = solve_coefficients(k3, k8, 1.0)
    kx_at_mg_resid = extra_residual(kx_at_mg)

    no_lift_best = scan_kx(k3, k8)

    def masses(coeffs: dict[str, float]) -> dict[str, float]:
        return {
            "M_Sigma3_GeV": coeffs["kappa_Sigma3"] * mg,
            "M_Sigma8_GeV": coeffs["kappa_Sigma8"] * mg,
            "M_SigmaR_GeV": coeffs["kappa_SigmaR"] * mg,
            "M_X_622_GeV": coeffs["kappa_X_622"] * mg,
        }

    return {
        "target": {
            "MG_GeV": mg,
            "kappa_Sigma3": k3,
            "kappa_Sigma8": k8,
            "base_residual_l2": best["residual_l2"],
            "tau_dim6_years": best["tau_dim6_years"],
            "tau_dim5_target_filter_years": best["tau_dim5_target_filter_years"],
        },
        "best_54_plus_210_without_projector_lift": no_lift_best,
        "X_at_MG_without_R_lift": {
            "coefficients": kx_at_mg,
            "masses": masses(kx_at_mg),
            "extra_threshold": kx_at_mg_resid,
        },
        "R_at_MG_before_X_lift": {
            "coefficients": r_at_mg,
            "masses": masses(r_at_mg),
            "extra_threshold": r_at_mg_resid,
        },
        "projector_lifted_solution": {
            "coefficients": lifted,
            "masses": masses(lifted),
            "extra_threshold": lifted_resid,
            "projectors": {
                "P_R": "(D_210^2-D_210)/2",
                "P_X": "-(9/25)(F_54-2)(F_54+4/3)",
            },
            "interpretation": "Use the 210 D-parity operator D_210 and the 54 Clebsch operator F_54 as spectral projectors.  The listed delta_X_projector raises the (6,2,2) fragment to MG while leaving Sigma_L and Sigma8 untouched.",
        },
    }


def write_report(payload: dict[str, object]) -> None:
    spec = payload["direct_210"]
    scen = payload["scenarios"]
    target = scen["target"]
    lifted = scen["projector_lifted_solution"]
    before = scen["R_at_MG_before_X_lift"]
    no_lift = scen["best_54_plus_210_without_projector_lift"]

    lines: list[str] = []
    lines.append("# Direct 45-45-210 Clebsch and lifting-sector calculation")
    lines.append("")
    lines.append("No web lookup was used.  This is a direct SO(10) four-form computation.")
    lines.append("")
    lines.append("## Direct 210 invariant")
    lines.append("")
    lines.append("The 210 is represented by a four-form `Phi_ijkl`.  For the Pati-Salam")
    lines.append("preserving D-parity direction `Phi_{7,8,9,10}=+1`")
    lines.append("(one-based SO(10) labels),")
    lines.append("")
    lines.append("```text")
    lines.append("W_210 = (zeta/8) Phi_ijkl A_ij A_kl")
    lines.append("```")
    lines.append("")
    lines.append("acts as the SO(4) Hodge star on weak two-forms.  The direct 45-by-45")
    lines.append("matrix has eigenvalue multiplicities")
    lines.append("")
    lines.append("```text")
    for eig, mult in sorted(spec["eigenvalue_multiplicities"].items(), key=lambda kv: float(kv[0])):
        lines.append(f"{eig}: {mult}")
    lines.append("```")
    lines.append("")
    lines.append("Thus, after normalization,")
    lines.append("")
    lines.append("```text")
    lines.append("D_210(Sigma_L) = +1")
    lines.append("D_210(Sigma_R) = -1")
    lines.append("D_210(Sigma8)  =  0")
    lines.append("D_210(X_622)   =  0")
    lines.append("```")
    lines.append("")
    lines.append("The sign of the first two lines is an orientation convention; flipping")
    lines.append("`Phi_{7,8,9,10}` exchanges the labels.")
    lines.append("")
    lines.append("## 54+210 linear spectrum")
    lines.append("")
    lines.append("Use")
    lines.append("")
    lines.append("```text")
    lines.append("k_r = mu + a54 F_54(r) + b210 D_210(r)")
    lines.append("```")
    lines.append("")
    lines.append("with `F_54=(2,2,-4/3,1/3)` and `D_210=(+1,-1,0,0)` on")
    lines.append("`(Sigma_L,Sigma_R,Sigma8,X_622)`.  Fixing `Sigma_R` at `MG` gives")
    lines.append("")
    lines.append("```text")
    coeff = before["coefficients"]
    lines.append(f"mu   = {coeff['mu']:.9f}")
    lines.append(f"a54  = {coeff['a54']:.9f}")
    lines.append(f"b210 = {coeff['b210']:.9f}")
    lines.append(f"kappa_Sigma3 = {coeff['kappa_Sigma3']:.9f}")
    lines.append(f"kappa_Sigma8 = {coeff['kappa_Sigma8']:.9f}")
    lines.append(f"kappa_SigmaR = {coeff['kappa_SigmaR']:.9f}")
    lines.append(f"kappa_X_622  = {coeff['kappa_X_622']:.9f}")
    lines.append("```")
    lines.append("")
    lines.append("The unwanted `(6,2,2)` fragment is then lifted by the spectral projector")
    lines.append("")
    lines.append("```text")
    lines.append("P_X = -(9/25)(F_54-2)(F_54+4/3)")
    lines.append("delta_X = 1-kappa_X_622")
    lines.append("```")
    lines.append("")
    lines.append("which is one on `X_622` and zero on `Sigma_L`, `Sigma_R`, and `Sigma8`.")
    lines.append("")
    lines.append("## Numerical verification")
    lines.append("")
    masses = lifted["masses"]
    coeff = lifted["coefficients"]
    lines.append("```text")
    lines.append(f"MG = {target['MG_GeV']:.6e} GeV")
    lines.append(f"M_Sigma3 = {masses['M_Sigma3_GeV']:.6e} GeV")
    lines.append(f"M_Sigma8 = {masses['M_Sigma8_GeV']:.6e} GeV")
    lines.append(f"M_SigmaR = {masses['M_SigmaR_GeV']:.6e} GeV")
    lines.append(f"M_X_622  = {masses['M_X_622_GeV']:.6e} GeV")
    lines.append(f"delta_X_projector = {coeff['delta_X_projector']:.9f}")
    lines.append(f"extra residual after projector lift = {lifted['extra_threshold']['residual_l2']:.3e}")
    lines.append(f"base threshold residual = {target['base_residual_l2']:.3e}")
    lines.append(f"tau_dim6 = {target['tau_dim6_years']:.6e} yr")
    lines.append(f"tau_dim5(S_T=1e-5) = {target['tau_dim5_target_filter_years']:.6e} yr")
    lines.append("```")
    lines.append("")
    lines.append("Without the projector lift, the best one-parameter 54+210 choice leaves")
    lines.append(f"an extra-threshold residual `{no_lift['extra_threshold']['residual_l2']:.6e}`.")
    lines.append("So the direct 210 contraction fixes the chirality problem, but an explicit")
    lines.append("lifting/projector sector is still needed for exact matching.")
    (OUT / "clebsch_210_lifting_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source = json.loads(IN.read_text(encoding="utf-8"))
    best = source["best"]
    direct = direct_210_spectrum()
    scenarios = scenario_payload(best)
    payload = {
        "input_summary": str(IN),
        "direct_210": direct,
        "clebsch_operators": {"F54": C54, "D210": C210},
        "scenarios": scenarios,
    }
    (OUT / "clebsch_210_lifting_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)
    lifted = scenarios["projector_lifted_solution"]
    print("Direct 45-45-210 Clebsch and lifting-sector calculation")
    print("  210 eigenvalue multiplicities:", direct["eigenvalue_multiplicities"])
    coeff = lifted["coefficients"]
    print(f"  mu={coeff['mu']:.9f} a54={coeff['a54']:.9f} b210={coeff['b210']:.9f}")
    print(f"  kappa Sigma3={coeff['kappa_Sigma3']:.9f}")
    print(f"  kappa Sigma8={coeff['kappa_Sigma8']:.9f}")
    print(f"  kappa SigmaR={coeff['kappa_SigmaR']:.9f}")
    print(f"  kappa X_622={coeff['kappa_X_622']:.9f}")
    print(f"  delta_X_projector={coeff['delta_X_projector']:.9f}")
    print(f"  extra residual after lift={lifted['extra_threshold']['residual_l2']:.3e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
