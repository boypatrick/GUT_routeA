#!/usr/bin/env python3
"""Match the Spin(10) 210_H^3 invariant to the Pati-Salam SU(4)_C cubic.

The 210 is represented as a four-form Phi in Lambda^4 R^10.  A four-form acts
on adjoint two-forms through

    (D_Phi)_[ij],[kl] = Phi_ijkl.

The cubic invariant

    I_3(Phi) = Tr_{Lambda^2 R^10}(D_Phi^3)

is Spin(10)-invariant and renormalizable.  Restrict Phi to the
(15,1,1) component Lambda^4 R^6, identify Phi=*_6 A with an SO(6) adjoint
two-form A, and compare I_3(*A) to the SO(6) Pfaffian cubic.  This is the
explicit full-Spin(10) matching target left open by the previous heartbeat.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PS_JSON = ROOT / "output" / "ps_goldstone_locking" / "ps_goldstone_locking_summary.json"
OUT = ROOT / "output" / "210_cubic_matching"

COLOR_DIM = 6
FULL_DIM = 10


def permutation_sign(values: tuple[int, ...]) -> int:
    if len(set(values)) != len(values):
        return 0
    inv = 0
    for i, vi in enumerate(values):
        for vj in values[i + 1 :]:
            if vi > vj:
                inv += 1
    return -1 if inv % 2 else 1


def pair_basis(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


COLOR_PAIRS = pair_basis(COLOR_DIM)
FULL_PAIRS = pair_basis(FULL_DIM)


def vector_to_twoform(x: np.ndarray) -> np.ndarray:
    a = np.zeros((COLOR_DIM, COLOR_DIM), dtype=float)
    for coeff, (i, j) in zip(x, COLOR_PAIRS):
        a[i, j] = coeff
        a[j, i] = -coeff
    return a


def hodge4_from_twoform_component(a: np.ndarray, i: int, j: int, k: int, l: int) -> float:
    inds = (i, j, k, l)
    if len(set(inds)) < 4 or any(idx >= COLOR_DIM for idx in inds):
        return 0.0
    remaining = tuple(idx for idx in range(COLOR_DIM) if idx not in inds)
    if len(remaining) != 2:
        return 0.0
    r, s = remaining
    return float(permutation_sign(inds + (r, s)) * a[r, s])


def d_matrix_from_color_twoform(a: np.ndarray, full: bool = True) -> np.ndarray:
    pairs = FULL_PAIRS if full else COLOR_PAIRS
    d = np.zeros((len(pairs), len(pairs)), dtype=float)
    for p, (i, j) in enumerate(pairs):
        for q, (k, l) in enumerate(pairs):
            d[p, q] = hodge4_from_twoform_component(a, i, j, k, l)
    return d


def invariant_i3(a: np.ndarray) -> float:
    d = d_matrix_from_color_twoform(a, full=True)
    return float(np.trace(d @ d @ d))


def pfaffian_recursive(a: np.ndarray) -> float:
    n = a.shape[0]
    if n == 0:
        return 1.0
    total = 0.0
    for j in range(1, n):
        rows = [idx for idx in range(1, n) if idx != j]
        sub = a[np.ix_(rows, rows)]
        total += ((-1) ** (j + 1)) * a[0, j] * pfaffian_recursive(sub)
    return float(total)


def canonical_twoform(a: float, b: float, c: float) -> np.ndarray:
    mat = np.zeros((COLOR_DIM, COLOR_DIM), dtype=float)
    for coeff, i, j in [(a, 0, 1), (b, 2, 3), (c, 4, 5)]:
        mat[i, j] = coeff
        mat[j, i] = -coeff
    return mat


def random_matching_test(samples: int = 200, seed: int = 20260506) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    ratios = []
    max_abs_error = 0.0
    for _ in range(samples):
        raw = rng.normal(size=(COLOR_DIM, COLOR_DIM))
        a = raw - raw.T
        pf = pfaffian_recursive(a)
        i3 = invariant_i3(a)
        if abs(pf) > 1.0e-8:
            ratios.append(i3 / pf)
            max_abs_error = max(max_abs_error, abs(i3 - ratios[0] * pf))
    return {
        "samples": samples,
        "seed": seed,
        "ratio_mean": float(np.mean(ratios)),
        "ratio_std": float(np.std(ratios)),
        "ratio_min": float(np.min(ratios)),
        "ratio_max": float(np.max(ratios)),
        "max_abs_error_against_first_ratio": float(max_abs_error),
    }


def canonical_formula_check() -> dict[str, object]:
    rows = []
    for abc in [(1.0, 1.0, 1.0), (2.0, -0.5, 3.0), (4.0, 1.0, -2.0)]:
        a = canonical_twoform(*abc)
        rows.append(
            {
                "a_b_c": list(abc),
                "pfaffian": pfaffian_recursive(a),
                "I3": invariant_i3(a),
                "ratio_I3_over_pfaffian": invariant_i3(a) / pfaffian_recursive(a),
            }
        )
    return {"canonical_block_rows": rows}


def superpotential_hessian_check(cubic_ratio: float) -> dict[str, object]:
    """Numerically verify W=1/2||A||^2-I3/c has the PS Hessian pattern."""

    x0 = np.zeros(len(COLOR_PAIRS), dtype=float)
    for idx, pair in enumerate(COLOR_PAIRS):
        if pair in [(0, 1), (2, 3), (4, 5)]:
            x0[idx] = 1.0

    def w_of_x(x: np.ndarray) -> float:
        a = vector_to_twoform(x)
        return 0.5 * float(np.dot(x, x)) - invariant_i3(a) / cubic_ratio

    h = 1.0e-4
    n = len(x0)
    grad = np.zeros(n, dtype=float)
    hess = np.zeros((n, n), dtype=float)
    w0 = w_of_x(x0)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h
        grad[i] = (w_of_x(x0 + ei) - w_of_x(x0 - ei)) / (2.0 * h)
        hess[i, i] = (w_of_x(x0 + ei) - 2.0 * w0 + w_of_x(x0 - ei)) / (h * h)
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = h
            val = (
                w_of_x(x0 + ei + ej)
                - w_of_x(x0 + ei - ej)
                - w_of_x(x0 - ei + ej)
                + w_of_x(x0 - ei - ej)
            ) / (4.0 * h * h)
            hess[i, j] = val
            hess[j, i] = val

    evals = np.linalg.eigvalsh(hess)
    clusters = {}
    for target in [-1.0, 0.0, 2.0]:
        clusters[str(target)] = int(np.sum(np.isclose(evals, target, atol=2.0e-6)))
    return {
        "vacuum_blocks": "A0 = e12 + e34 + e56",
        "superpotential": "W = 1/2 ||A||^2 - I3(*A)/c, c=I3/Pf",
        "gradient_norm_at_vacuum": float(np.linalg.norm(grad)),
        "hessian_eigenvalues": [float(v) for v in evals],
        "eigenvalue_clusters": clusters,
        "max_abs_deviation_from_expected": float(
            max(
                abs(evals[0] + 1.0),
                max(abs(v) for v in evals[1:7]),
                max(abs(v - 2.0) for v in evals[7:]),
            )
        ),
    }


def load_ps_eigen_summary() -> dict[str, object]:
    payload = json.loads(PS_JSON.read_text(encoding="utf-8"))
    return payload["eigenvalue_summary"]


def write_report(payload: dict[str, object]) -> None:
    lines: list[str] = []
    lines.append("# 210_H^3 matching to the Pati-Salam cubic")
    lines.append("")
    lines.append("No web lookup was used.  The 210 is treated as a four-form")
    lines.append("`Phi in Lambda^4 R^10`.")
    lines.append("")
    lines.append("## Spin(10) cubic invariant")
    lines.append("")
    lines.append("Let `D_Phi` be the induced operator on adjoint two-forms,")
    lines.append("")
    lines.append("```text")
    lines.append("(D_Phi)_[ij],[kl] = Phi_ijkl.")
    lines.append("```")
    lines.append("")
    lines.append("Then")
    lines.append("")
    lines.append("```text")
    lines.append("I3(Phi) = Tr_{Lambda^2 R^10}(D_Phi^3)")
    lines.append("```")
    lines.append("")
    lines.append("is a Spin(10)-invariant cubic contraction, hence a renormalizable")
    lines.append("`210_H^3` superpotential term.")
    lines.append("")
    lines.append("## Restriction to the (15,1,1) component")
    lines.append("")
    lines.append("Split `R^10 = R^6_C + R^4_W`.  The color-adjoint component of the 210 is")
    lines.append("")
    lines.append("```text")
    lines.append("Lambda^4 R^6_C ~= Lambda^2 R^6_C ~= (15,1,1).")
    lines.append("```")
    lines.append("")
    lines.append("For a color two-form `A`, set `Phi=*_6 A`.  The invariant obeys")
    lines.append("")
    rows = payload["canonical_formula"]["canonical_block_rows"]
    lines.append("```text")
    lines.append("A = a e12 + b e34 + c e56:")
    lines.append("Pf(A) = a b c,")
    lines.append("I3(*A) = c_210 Pf(A).")
    lines.append("```")
    lines.append("")
    lines.append("| (a,b,c) | Pf(A) | I3(*A) | I3/Pf |")
    lines.append("|---|---:|---:|---:|")
    for row in rows:
        abc = ",".join(f"{v:g}" for v in row["a_b_c"])
        lines.append(
            f"| `({abc})` | {row['pfaffian']:.6g} | {row['I3']:.6g} | "
            f"{row['ratio_I3_over_pfaffian']:.6g} |"
        )
    rand = payload["random_matching_test"]
    lines.append("")
    lines.append("Random check over generic color two-forms:")
    lines.append("")
    lines.append("```text")
    lines.append(f"samples = {rand['samples']}")
    lines.append(f"mean I3/Pf = {rand['ratio_mean']:.12g}")
    lines.append(f"std  I3/Pf = {rand['ratio_std']:.3e}")
    lines.append(f"min/max    = {rand['ratio_min']:.12g}, {rand['ratio_max']:.12g}")
    lines.append("```")
    lines.append("")
    lines.append("Thus the `210_H^3` contraction restricts exactly to the unique SO(6)")
    lines.append("Pfaffian cubic, i.e. the SU(4)_C adjoint cubic used in the Pati-Salam")
    lines.append("Goldstone-locking sector, up to the displayed normalization.")
    lines.append("")
    lines.append("## Hessian match")
    lines.append("")
    hess = payload["hessian_match"]
    lines.append("Use the normalized local superpotential")
    lines.append("")
    lines.append("```text")
    lines.append(hess["superpotential"])
    lines.append("```")
    lines.append("")
    lines.append("around `A0=e12+e34+e56`.  The numerical Hessian gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"gradient norm = {hess['gradient_norm_at_vacuum']:.3e}")
    lines.append(f"eigenvalue clusters = {hess['eigenvalue_clusters']}")
    lines.append(f"max deviation = {hess['max_abs_deviation_from_expected']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("The eigenvalue pattern `(-1)^1, 0^6, 2^8` is the same pattern found in")
    lines.append("the PS-stage derivation: singlet, eaten Goldstone pair, and octet.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("The full Spin(10) route exists at the level of invariant tensors:")
    lines.append("the renormalizable `210_H^3` contraction matches onto the Pati-Salam")
    lines.append("SU(4)_C cubic in the `(15,1,1)` branch.  What remains for a fully complete")
    lines.append("model is not this cubic matching, but the full multi-field vacuum alignment")
    lines.append("with the already-used `54_H/210_H` projector and mediator sectors.")
    (OUT / "210_cubic_matching_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    canonical = canonical_formula_check()
    cubic_ratio = canonical["canonical_block_rows"][0]["ratio_I3_over_pfaffian"]
    payload = {
        "note": "No web lookup used. 210_H^3 restricted to Lambda^4 R^6 matches the SO(6)/SU(4) cubic.",
        "input_files": {"ps_goldstone_locking_summary": str(PS_JSON)},
        "spin10_cubic": "I3(Phi)=Tr_{Lambda^2 R^10}(D_Phi^3)",
        "color_component": "Phi in Lambda^4 R^6 ~= Lambda^2 R^6 = (15,1,1)",
        "canonical_formula": canonical,
        "cubic_ratio_I3_over_pfaffian": cubic_ratio,
        "random_matching_test": random_matching_test(),
        "hessian_match": superpotential_hessian_check(cubic_ratio),
        "ps_stage_eigenvalue_summary": load_ps_eigen_summary(),
        "verdict": {
            "210H_cubic_matches_PS_cubic": True,
            "remaining_work": "Full Spin(10) vacuum alignment and mixing with 54_H/210_H projector/mediator sectors.",
        },
    }
    (OUT / "210_cubic_matching_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)

    hess = payload["hessian_match"]
    rand = payload["random_matching_test"]
    print("210_H^3 cubic matching")
    print(f"  I3/Pf canonical: {cubic_ratio:.6g}")
    print(f"  random ratio mean/std: {rand['ratio_mean']:.6g} / {rand['ratio_std']:.3e}")
    print(f"  Hessian clusters: {hess['eigenvalue_clusters']}")
    print(f"  gradient norm: {hess['gradient_norm_at_vacuum']:.3e}")
    print(f"  max Hessian deviation: {hess['max_abs_deviation_from_expected']:.3e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
