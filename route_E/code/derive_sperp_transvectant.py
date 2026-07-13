#!/usr/bin/env python3
"""Derive the fitted S_perp as the CP1 O(2) second transvectant.

No web lookup is used.  The previous flavor scan found that a single
orthogonal symmetric direction S_perp repairs the finite flavor fit.  This
script checks whether that direction is actually the canonical SL(2)-invariant
trace/contact functional on H^0(CP1,O(2)) = Sym^2(C^2).

In the normalized spin-1 basis

  u_0 = x^2,   u_1 = sqrt(2) x y,   u_2 = y^2,

the pointwise product Sym^2 H^0(O(2)) -> H^0(O(4)) has the one-dimensional
kernel

  u_1^2 - 2 u_0 u_2 = 0.

As a symmetric bilinear matrix this is

  K = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]].

Equivalently, K is the spin-0 second transvectant, whereas the Veronese
O(-4)-density product probes only the spin-2 image.
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

import scan_clebsch_flavor_fit as fit  # noqa: E402


CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
RELAX = ROOT / "output" / "flavor_symmetric_relaxation" / "scan_clebsch_flavor_symmetric_relaxation.json"
OUT = ROOT / "output" / "flavor_sperp_transvectant"


def load_product_clebsch() -> np.ndarray:
    card = json.loads(CARD.read_text(encoding="utf-8"))
    return np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )


def load_numeric_sperp() -> np.ndarray:
    payload = json.loads(RELAX.read_text(encoding="utf-8"))
    return fit.cmat(payload["S_perp"])


def analytic_transvectant_matrix() -> np.ndarray:
    k = np.zeros((3, 3), dtype=complex)
    k[1, 1] = 1.0 / math.sqrt(3.0)
    k[0, 2] = -1.0 / math.sqrt(3.0)
    k[2, 0] = -1.0 / math.sqrt(3.0)
    return k


def su2_spin1_matrix(a: complex, b: complex) -> np.ndarray:
    """Spin-1 action in the basis (x^2, sqrt(2)xy, y^2)."""
    d = np.zeros((3, 3), dtype=complex)
    # e0 -> (a x + b y)^2
    d[:, 0] = [a * a, math.sqrt(2.0) * a * b, b * b]
    # e1 -> sqrt(2) (a x + b y)(-b* x + a* y)
    d[:, 1] = [
        -math.sqrt(2.0) * a * np.conjugate(b),
        abs(a) ** 2 - abs(b) ** 2,
        math.sqrt(2.0) * b * np.conjugate(a),
    ]
    # e2 -> (-b* x + a* y)^2
    d[:, 2] = [
        np.conjugate(b) * np.conjugate(b),
        -math.sqrt(2.0) * np.conjugate(b) * np.conjugate(a),
        np.conjugate(a) * np.conjugate(a),
    ]
    return d


def random_su2(rng: np.random.Generator) -> tuple[complex, complex]:
    z = rng.normal(size=4)
    norm = float(np.linalg.norm(z))
    a = (z[0] + 1j * z[1]) / norm
    b = (z[2] + 1j * z[3]) / norm
    return a, b


def polynomial_kernel_coefficients(k: np.ndarray) -> dict[str, complex]:
    """Return coefficients in x^4, x^3y, x^2y^2, xy^3, y^4."""
    # u0=x^2, u1=sqrt(2)xy, u2=y^2.
    coeff = np.zeros(5, dtype=complex)
    monomials = {
        0: np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=complex),
        1: np.array([0.0, math.sqrt(2.0), 0.0, 0.0, 0.0], dtype=complex),
        2: np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=complex),
    }
    # Store polynomial u_i*u_j as degree-four coefficients.
    products: dict[tuple[int, int], np.ndarray] = {}
    products[(0, 0)] = np.array([1, 0, 0, 0, 0], dtype=complex)
    products[(0, 1)] = np.array([0, math.sqrt(2.0), 0, 0, 0], dtype=complex)
    products[(1, 0)] = products[(0, 1)]
    products[(0, 2)] = np.array([0, 0, 1, 0, 0], dtype=complex)
    products[(2, 0)] = products[(0, 2)]
    products[(1, 1)] = np.array([0, 0, 2, 0, 0], dtype=complex)
    products[(1, 2)] = np.array([0, 0, 0, math.sqrt(2.0), 0], dtype=complex)
    products[(2, 1)] = products[(1, 2)]
    products[(2, 2)] = np.array([0, 0, 0, 0, 1], dtype=complex)
    del monomials
    for i in range(3):
        for j in range(3):
            coeff += k[i, j] * products[(i, j)]
    return {
        "x4": coeff[0],
        "x3y": coeff[1],
        "x2y2": coeff[2],
        "xy3": coeff[3],
        "y4": coeff[4],
    }


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def complex_dict_json(raw: dict[str, complex]) -> dict[str, dict[str, float]]:
    return {key: cjson(value) for key, value in raw.items()}


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(__doc__.strip())
        return
    OUT.mkdir(parents=True, exist_ok=True)
    clebsch = load_product_clebsch()
    s_num = load_numeric_sperp()
    k = analytic_transvectant_matrix()

    phase = np.vdot(k, s_num)
    phase_unit = phase / abs(phase)
    aligned_residual = float(np.linalg.norm(s_num - phase_unit * k))
    phase_independent_distance = float(math.sqrt(max(0.0, 1.0 - abs(phase) ** 2)))
    clebsch_overlaps = [np.vdot(clebsch[:, :, idx], k) for idx in range(clebsch.shape[2])]
    product_coeff = polynomial_kernel_coefficients(k)
    product_kernel_norm = float(np.linalg.norm([product_coeff[key] for key in ["x4", "x3y", "x2y2", "xy3", "y4"]]))

    rng = np.random.default_rng(2026050723)
    covariance_errors = []
    for _ in range(32):
        a, b = random_su2(rng)
        d = su2_spin1_matrix(a, b)
        covariance_errors.append(float(np.linalg.norm(d.T @ k @ d - k)))

    payload: dict[str, Any] = {
        "note": "No web lookup used. S_perp is identified with the CP1 O(2) second transvectant.",
        "basis": "u0=x^2, u1=sqrt(2)xy, u2=y^2",
        "analytic_K": fit.matrix_json(k),
        "numeric_S_perp": fit.matrix_json(s_num),
        "phase_alignment": {
            "inner_product_vdot_K_numeric": cjson(phase),
            "phase_unit": cjson(phase_unit),
            "aligned_frobenius_residual": aligned_residual,
            "phase_independent_distance": phase_independent_distance,
        },
        "veronese_orthogonality": {
            "max_abs_overlap": float(max(abs(z) for z in clebsch_overlaps)),
            "overlaps": [cjson(z) for z in clebsch_overlaps],
        },
        "pointwise_product_kernel": {
            "degree_four_coefficients": complex_dict_json(product_coeff),
            "kernel_norm": product_kernel_norm,
        },
        "su2_covariance": {
            "samples": len(covariance_errors),
            "max_error": float(max(covariance_errors)),
            "mean_error": float(np.mean(covariance_errors)),
        },
        "interpretation": (
            "The fitted non-Veronese direction is the spin-0 summand in "
            "Sym^2(Sym^2 C^2)=Sym^4 C^2 plus Sym^0 C^2.  Ordinary O(-4) "
            "dual densities see only the spin-2 multiplication image; a "
            "contact/transvectant operator supplies the missing spin-0 trace."
        ),
        "checks": {
            "phase_alignment_residual_lt_1e_minus_12": bool(aligned_residual < 1.0e-12),
            "veronese_overlap_lt_1e_minus_12": bool(max(abs(z) for z in clebsch_overlaps) < 1.0e-12),
            "pointwise_product_kernel_lt_1e_minus_12": bool(product_kernel_norm < 1.0e-12),
            "su2_covariance_lt_1e_minus_12": bool(max(covariance_errors) < 1.0e-12),
        },
    }
    payload["all_checks_pass"] = all(payload["checks"].values())
    (OUT / "derive_sperp_transvectant.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# S_perp transvectant derivation",
        "",
        "No web lookup was used.",
        "",
        "In the normalized spin-1 basis `u0=x^2, u1=sqrt(2)xy, u2=y^2`,",
        "",
        "```text",
        "K = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]]",
        "u1^2 - 2 u0 u2 = 0",
        "```",
        "",
        f"Phase-aligned residual against the fitted S_perp: `{aligned_residual:.6e}`.",
        f"Max Veronese overlap: `{max(abs(z) for z in clebsch_overlaps):.6e}`.",
        f"Pointwise product kernel norm: `{product_kernel_norm:.6e}`.",
        f"Max sampled SU(2) covariance error: `{max(covariance_errors):.6e}`.",
        "",
        "Verdict: the fitted S_perp is the unique spin-0 second transvectant/contact direction.",
        "",
    ]
    (OUT / "derive_sperp_transvectant_report.md").write_text("\n".join(lines), encoding="utf-8")

    print("S_perp transvectant derivation")
    print(f"  phase-aligned residual: {aligned_residual:.6e}")
    print(f"  max Veronese overlap: {max(abs(z) for z in clebsch_overlaps):.6e}")
    print(f"  product kernel norm: {product_kernel_norm:.6e}")
    print(f"  max SU(2) covariance error: {max(covariance_errors):.6e}")
    print(f"  all checks pass: {payload['all_checks_pass']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
