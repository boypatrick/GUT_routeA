#!/usr/bin/env python3
"""Direct Spin(10) 45-45-54 Clebsch calculation.

Spin(10) is represented through SO(10) vector indices.  The adjoint 45 is an
antisymmetric tensor A_ij, and the 54 is a symmetric traceless tensor S_ij.  The
unique quadratic invariant relevant for masses is

    W_54 = (eta/2) S_ij A_ik A_jk.

For a diagonal Pati-Salam preserving 54 vev

    S = diag(-2,-2,-2,-2,-2,-2, 3,3,3,3),

the Clebsch eigenvalue on the plane generator A_ab is s_a+s_b.  This script
computes the full 45 Hessian, decomposes the entries by SO(6)xSO(4) blocks, and
tests whether a single 54 filter can explain the required Sigma_3 threshold.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "output" / "yukawa_superpotential_rge" / "yukawa_superpotential_rge_summary.json"
OUT = ROOT / "output" / "clebsch_54"

HEAVY_BASIS_BASE = np.array(
    [
        [2.0 / 5.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [1.0, 0.0, 3.0],
    ],
    dtype=float,
)
PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0

# Extra chiral fragments of a 45 if the same 54 mass term is applied literally.
B_SIGMA_R = np.array([6.0 / 5.0, 0.0, 0.0], dtype=float)
B_X_622 = np.array([26.0 / 5.0, 6.0, 4.0], dtype=float)


def pair_basis() -> list[tuple[int, int]]:
    return [(i, j) for i in range(10) for j in range(i + 1, 10)]


def classify_pair(i: int, j: int) -> str:
    if i < 6 and j < 6:
        return "(15,1,1)"
    if i >= 6 and j >= 6:
        # SO(4) adjoint = (3,1)+(1,3).  The 54 cannot split these two.
        return "(1,3,1)+(1,1,3)"
    return "(6,2,2)"


def clebsch_spectrum() -> dict[str, object]:
    s = np.array([-2.0] * 6 + [3.0] * 4, dtype=float)
    rows = []
    by_block: dict[str, list[float]] = defaultdict(list)
    for i, j in pair_basis():
        c = float(s[i] + s[j])
        block = classify_pair(i, j)
        rows.append({"i": i + 1, "j": j + 1, "block": block, "clebsch_raw": c})
        by_block[block].append(c)
    summary = {}
    for block, vals in by_block.items():
        summary[block] = {
            "multiplicity": len(vals),
            "raw_values": sorted(Counter(vals).items()),
            "normalized_to_weak_2": sorted(Counter([v / 3.0 for v in vals]).items()),
        }
    return {"diag_54": s.tolist(), "basis_rows": rows, "block_summary": summary}


def solve_54_parameters(best: dict[str, float]) -> dict[str, float]:
    mg = float(best["MG_GeV"])
    k3 = float(best["M_Sigma3_GeV"]) / mg
    k8 = float(best["M_Sigma8_GeV"]) / mg
    # Use normalized Clebsches C_3=+2 and C_8=-4/3.
    lam = 0.4 * k3 + 0.6 * k8
    chi = (k3 / lam - 1.0) / 2.0
    kx = lam * (1.0 + chi / 3.0)  # normalized C_(6,2,2)=+1/3
    return {
        "lambda_54": lam,
        "chi_54": chi,
        "C_Sigma3": 2.0,
        "C_Sigma8": -4.0 / 3.0,
        "C_X_622": 1.0 / 3.0,
        "kappa_Sigma3": k3,
        "kappa_Sigma8": k8,
        "kappa_SigmaR": k3,
        "kappa_X_622": kx,
        "M_SigmaR_GeV": k3 * mg,
        "M_X_622_GeV": kx * mg,
        "fits_target_Sigma3": lam * (1.0 + 2.0 * chi),
        "fits_target_Sigma8": lam * (1.0 - 4.0 * chi / 3.0),
    }


def threshold_with_literal_54_extras(best: dict[str, float], pars: dict[str, float]) -> dict[str, object]:
    alpha_inv = np.array(
        [best["alpha1_inv_MG"], best["alpha2_inv_MG"], best["alpha3_inv_MG"]],
        dtype=float,
    )
    logs_base = np.array([best["log_HC"], best["log_Sigma3"], best["log_Sigma8"]], dtype=float)
    delta_base = HEAVY_BASIS_BASE @ logs_base / (2.0 * math.pi)
    log_sigma_r = -math.log(pars["kappa_SigmaR"])
    log_x = -math.log(pars["kappa_X_622"])
    delta_extra = (B_SIGMA_R * log_sigma_r + B_X_622 * log_x) / (2.0 * math.pi)
    residual_base = PROJECTOR @ (alpha_inv - delta_base)
    residual_with_extra = PROJECTOR @ (alpha_inv - delta_base - delta_extra)
    return {
        "logs_base": logs_base.tolist(),
        "delta_base": delta_base.tolist(),
        "extra_fragments": {
            "Sigma_R_(1,1,3)": {
                "mass_GeV": pars["M_SigmaR_GeV"],
                "log_MG_over_M": log_sigma_r,
                "b": B_SIGMA_R.tolist(),
            },
            "X_(6,2,2)": {
                "mass_GeV": pars["M_X_622_GeV"],
                "log_MG_over_M": log_x,
                "b": B_X_622.tolist(),
            },
        },
        "delta_extra": delta_extra.tolist(),
        "base_residual_l2": float(np.linalg.norm(residual_base)),
        "literal_54_extra_residual_l2": float(np.linalg.norm(residual_with_extra)),
        "literal_54_extra_residual_vector": residual_with_extra.tolist(),
    }


def write_report(payload: dict[str, object]) -> None:
    spec = payload["clebsch"]
    pars = payload["fit_54_parameters"]
    thresh = payload["threshold_literal_54_extras"]
    lines: list[str] = []
    lines.append("# Direct 45-45-54 Clebsch calculation")
    lines.append("")
    lines.append("No web lookup was used.  This is a direct SO(10)-index computation.")
    lines.append("")
    lines.append("## Invariant")
    lines.append("")
    lines.append("Let `A_ij=-A_ji` be the adjoint `45` and let `S_ij=S_ji`, `Tr S=0`, be")
    lines.append("the `54`.  The Spin(10)-invariant quadratic contraction is")
    lines.append("")
    lines.append("```text")
    lines.append("W_54 = (eta/2) S_ij A_ik A_jk.")
    lines.append("```")
    lines.append("")
    lines.append("For the Pati-Salam preserving direction")
    lines.append("")
    lines.append("```text")
    lines.append("S = diag(-2,-2,-2,-2,-2,-2, 3,3,3,3),")
    lines.append("```")
    lines.append("")
    lines.append("a plane generator `A_ab` has Clebsch eigenvalue `s_a+s_b`.")
    lines.append("")
    lines.append("## Direct spectrum")
    lines.append("")
    lines.append("| SO(6)xSO(4) block | multiplicity | raw Clebsch | normalized weak=2 |")
    lines.append("|---|---:|---:|---:|")
    for block, row in spec["block_summary"].items():
        raw = row["raw_values"][0][0]
        norm = row["normalized_to_weak_2"][0][0]
        lines.append(f"| `{block}` | {row['multiplicity']} | {raw:.6g} | {norm:.6g} |")
    lines.append("")
    lines.append("Consequently a single `54_H` predicts")
    lines.append("")
    lines.append("```text")
    lines.append("C_(15,1,1) = -4/3,")
    lines.append("C_(1,3,1) = C_(1,1,3) = +2,")
    lines.append("C_(6,2,2) = +1/3,")
    lines.append("```")
    lines.append("")
    lines.append("after normalizing the weak adjoint Clebsch to `+2`.")
    lines.append("")
    lines.append("## Fit to the required Sigma3/Sigma8 masses")
    lines.append("")
    lines.append("Using")
    lines.append("")
    lines.append("```text")
    lines.append("M_Sigma3/MG = lambda_54 (1 + 2 chi_54)")
    lines.append("M_Sigma8/MG = lambda_54 (1 - 4 chi_54/3)")
    lines.append("```")
    lines.append("")
    lines.append("gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"lambda_54 = {pars['lambda_54']:.9f}")
    lines.append(f"chi_54 = {pars['chi_54']:.9f}")
    lines.append(f"M_SigmaR = {pars['M_SigmaR_GeV']:.6e} GeV")
    lines.append(f"M_(6,2,2) = {pars['M_X_622_GeV']:.6e} GeV")
    lines.append("```")
    lines.append("")
    lines.append("The target `Sigma3` and `Sigma8` masses are reproduced, but the literal")
    lines.append("single-54 contraction also leaves a degenerate `Sigma_R` and a moderately")
    lines.append("light `(6,2,2)` fragment.")
    lines.append("")
    lines.append("## Threshold check with literal 54 extras")
    lines.append("")
    lines.append("If those extra `45` fragments are kept at the same 54-Hessian masses, their")
    lines.append("additional one-loop threshold vector spoils the previous exact matching:")
    lines.append("")
    lines.append("```text")
    lines.append(f"base residual_l2 = {thresh['base_residual_l2']:.3e}")
    lines.append(f"literal 54 extra residual_l2 = {thresh['literal_54_extra_residual_l2']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("Therefore a lone `54_H` is not enough.  It is useful as a direct tensor")
    lines.append("calculation and as a no-go theorem: either a `210_H`-type contraction must")
    lines.append("split/lift the unwanted fragments, or an additional sector must remove")
    lines.append("`Sigma_R` and `(6,2,2)` while preserving the light `Sigma3` threshold.")
    (OUT / "clebsch_54_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads(IN.read_text(encoding="utf-8"))
    best = summary["best"]
    clebsch = clebsch_spectrum()
    pars = solve_54_parameters(best)
    threshold = threshold_with_literal_54_extras(best, pars)
    payload = {
        "input_summary": str(IN),
        "invariant": "W_54=(eta/2) S_ij A_ik A_jk",
        "clebsch": clebsch,
        "fit_54_parameters": pars,
        "threshold_literal_54_extras": threshold,
    }
    (OUT / "clebsch_54_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)
    print("Direct 45-45-54 Clebsch calculation")
    for block, row in clebsch["block_summary"].items():
        raw = row["raw_values"][0][0]
        norm = row["normalized_to_weak_2"][0][0]
        print(f"  {block}: mult={row['multiplicity']} raw={raw:.6g} norm={norm:.6g}")
    print(f"  lambda_54={pars['lambda_54']:.9f}")
    print(f"  chi_54={pars['chi_54']:.9f}")
    print(f"  literal 54 extra residual={threshold['literal_54_extra_residual_l2']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
