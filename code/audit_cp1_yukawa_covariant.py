#!/usr/bin/env python3
"""Covariant audit for the CP1 O(2) Yukawa geometry.

No web lookup is used.  The audit separates two mathematically well-defined
ways to repair the original toy overlap:

1. Superpotential route:
   s_i s_j is a section of O(4), so the Higgs functional must be an O(-4)
   dual density.  Numerically, this is represented by dual coefficients to the
   orthonormal O(4) basis.

2. Hermitian Toeplitz route:
   T_ij = int s_i^* h s_j dmu_FS for a scalar real h.  This is coordinate
   invariant and naturally belongs to the kinetic/Kahler or wavefunction
   sector rather than directly to a holomorphic superpotential.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_cp1_o2_yukawa as legacy  # noqa: E402


OUT = ROOT / "output" / "cp1_yukawa_covariant"
LEGACY_JSON = ROOT / "output" / "yukawa_o2" / "cp1_o2_yukawa_scan.json"


@dataclass(frozen=True)
class Target:
    name: str
    small: float
    mid: float


TARGETS = [
    Target("up", 1.0e-5, 7.0e-3),
    Target("down", 1.0e-3, 2.0e-2),
    Target("charged_lepton", 3.0e-4, 6.0e-2),
    Target("neutrino_dirac", 1.0e-2, 2.0e-1),
]


def o4_basis(grid: legacy.CP1O2Grid) -> np.ndarray:
    u = grid.u
    phi = grid.phi
    c = np.sqrt((1.0 + u) / 2.0)
    s = np.sqrt((1.0 - u) / 2.0)
    rows = []
    for m in range(5):
        coeff = math.sqrt(5.0 * math.comb(4, m))
        rows.append(coeff * np.exp(1j * m * phi) * (s**m) * (c ** (4 - m)))
    return np.stack(rows, axis=0)


def singular_ratios(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.linalg.svd(matrix, compute_uv=False)
    values = np.sort(values)[::-1]
    return values, values / values[0]


def log_score(norm: np.ndarray, target: Target) -> float:
    return float((math.log10(max(norm[2], 1.0e-16) / target.small)) ** 2 + (math.log10(max(norm[1], 1.0e-16) / target.mid)) ** 2)


def profile_real(grid: legacy.CP1O2Grid, params: legacy.KernelParams) -> np.ndarray:
    c1 = legacy.unit_vector(params.theta1, params.phi1)
    c2 = legacy.unit_vector(params.theta2, params.phi2)
    dot1 = c1 @ grid.points
    dot2 = c2 @ grid.points
    k1 = np.exp(params.kappa1 * (dot1 - 1.0))
    k2 = np.exp(params.kappa2 * (dot2 - 1.0))
    q2 = 0.5 * (3.0 * grid.nz**2 - 1.0)
    h = k1 + params.amp2 * k2 + params.eps + params.quad * q2
    # Shift by a tiny floor if the quadrupole makes isolated points negative.
    min_h = float(np.min(h))
    if min_h < 1.0e-10:
        h = h + (1.0e-10 - min_h)
    return h


def toeplitz_matrix(grid: legacy.CP1O2Grid, params: legacy.KernelParams) -> np.ndarray:
    h = profile_real(grid, params)
    weighted = grid.psi.conjugate() * (grid.weight * h)
    mat = weighted @ grid.psi.T
    return 0.5 * (mat + mat.conjugate().T)


def product_map(grid: legacy.CP1O2Grid, chi: np.ndarray) -> dict[str, Any]:
    # C_ij^m = int chi_m^* psi_i psi_j dmu.
    coeff = np.zeros((3, 3, 5), dtype=complex)
    max_reconstruction_error = 0.0
    for i in range(3):
        for j in range(3):
            q = grid.psi[i] * grid.psi[j]
            for m in range(5):
                coeff[i, j, m] = np.sum(grid.weight * np.conjugate(chi[m]) * q)
            reconstructed = np.sum(coeff[i, j, :, None] * chi, axis=0)
            max_reconstruction_error = max(max_reconstruction_error, float(np.max(np.abs(q - reconstructed))))
    relation = grid.psi[1] * grid.psi[1] - 2.0 * grid.psi[0] * grid.psi[2]
    return {
        "coefficients": coeff,
        "max_product_reconstruction_error": max_reconstruction_error,
        "veronese_relation_max_abs": float(np.max(np.abs(relation))),
    }


def legacy_dual_completion_check(grid: legacy.CP1O2Grid, chi: np.ndarray, coeff: np.ndarray) -> dict[str, Any]:
    if not LEGACY_JSON.exists():
        return {"available": False}
    payload = json.loads(LEGACY_JSON.read_text(encoding="utf-8"))
    checks = {}
    for name, item in payload["results"].items():
        params = legacy.KernelParams(**{key: float(value) for key, value in item["params"].items() if key != "center_separation"})
        h = grid.higgs_profile(params)
        direct = grid.y_matrix(params)
        dual_a = np.array([np.sum(grid.weight * chi[m] * h) for m in range(5)], dtype=complex)
        dual = np.einsum("ijm,m->ij", coeff, dual_a)
        _, norm = singular_ratios(direct)
        checks[name] = {
            "direct_vs_dual_max_abs": float(np.max(np.abs(direct - dual))),
            "veronese_matrix_relation_abs": float(abs(direct[1, 1] - 2.0 * direct[0, 2])),
            "normalized_singular_values": [float(x) for x in norm],
        }
    return {"available": True, "checks": checks}


def scan_toeplitz(grid: legacy.CP1O2Grid, samples: int = 7000) -> dict[str, Any]:
    rng = np.random.default_rng(20260507)
    keep = 40
    best: dict[str, list[tuple[float, legacy.KernelParams, np.ndarray]]] = {target.name: [] for target in TARGETS}
    for _ in range(samples):
        params = legacy.sample_params(rng)
        matrix = toeplitz_matrix(grid, params)
        _values, norm = singular_ratios(matrix)
        for target in TARGETS:
            score = log_score(norm, target)
            bucket = best[target.name]
            if len(bucket) < keep:
                bucket.append((score, params, norm))
                bucket.sort(key=lambda item: item[0])
            elif score < bucket[-1][0]:
                bucket[-1] = (score, params, norm)
                bucket.sort(key=lambda item: item[0])
    out: dict[str, Any] = {}
    for target in TARGETS:
        score, params, norm = best[target.name][0]
        matrix = toeplitz_matrix(grid, params)
        eigvals = np.linalg.eigvalsh(matrix)
        out[target.name] = {
            "target": asdict(target),
            "score": float(score),
            "params": legacy.params_to_json(params),
            "normalized_singular_values": [float(x) for x in norm],
            "eigenvalues": [float(x) for x in eigvals],
            "condition_number": float(norm[0] / max(norm[2], 1.0e-300)),
            "matrix": legacy.matrix_to_json(matrix / np.linalg.svd(matrix, compute_uv=False)[0]),
        }
    return {"samples": samples, "results": out}


def write_toeplitz_csv(scan: dict[str, Any]) -> None:
    rows = []
    for name, item in scan["results"].items():
        norm = item["normalized_singular_values"]
        rows.append(
            {
                "sector": name,
                "target_small": item["target"]["small"],
                "target_mid": item["target"]["mid"],
                "found_small": norm[2],
                "found_mid": norm[1],
                "score": item["score"],
                "condition_number": item["condition_number"],
                "center_separation": item["params"]["center_separation"],
                "kappa1": item["params"]["kappa1"],
                "kappa2": item["params"]["kappa2"],
                "amp2": item["params"]["amp2"],
                "eps": item["params"]["eps"],
                "quad": item["params"]["quad"],
            }
        )
    with (OUT / "toeplitz_yukawa_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# CP1 O(2) covariant Yukawa audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Bundle checks")
    lines.append("")
    bundle = payload["bundle_checks"]
    lines.append(f"O(2) Gram max error: `{bundle['gram_o2_max_error']:.3e}`")
    lines.append(f"O(4) Gram max error: `{bundle['gram_o4_max_error']:.3e}`")
    lines.append(f"Product reconstruction max error: `{bundle['product_reconstruction_max_error']:.3e}`")
    lines.append(f"Veronese relation max error: `{bundle['veronese_relation_max_abs']:.3e}`")
    lines.append("")
    lines.append("## Dual O(-4) completion of legacy matrices")
    lines.append("")
    for name, check in payload["dual_completion"]["checks"].items():
        norm = check["normalized_singular_values"]
        lines.append(
            f"- `{name}`: direct-dual error `{check['direct_vs_dual_max_abs']:.3e}`, "
            f"relation error `{check['veronese_matrix_relation_abs']:.3e}`, "
            f"singulars `[{norm[2]:.3e}, {norm[1]:.3e}, 1]`"
        )
    lines.append("")
    lines.append("## Hermitian Toeplitz scan")
    lines.append("")
    lines.append("| sector | target small | target mid | found small | found mid | score |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for name, item in payload["toeplitz_scan"]["results"].items():
        target = item["target"]
        norm = item["normalized_singular_values"]
        lines.append(
            f"| {name} | {target['small']:.1e} | {target['mid']:.1e} | "
            f"{norm[2]:.3e} | {norm[1]:.3e} | {item['score']:.3f} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "cp1_yukawa_covariant_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    legacy.require_verified_fields()
    OUT.mkdir(parents=True, exist_ok=True)
    grid = legacy.CP1O2Grid(n_theta=54, n_phi=108)
    chi = o4_basis(grid)
    gram_o2 = grid.gram()
    gram_o4 = (chi.conjugate() * grid.weight) @ chi.T
    pmap = product_map(grid, chi)
    coeff = pmap["coefficients"]
    dual = legacy_dual_completion_check(grid, chi, coeff)
    toeplitz = scan_toeplitz(grid)
    payload = {
        "note": "No web lookup used. Coordinate-covariant audit for CP1/O(2) Yukawa geometry.",
        "bundle_checks": {
            "gram_o2_max_error": float(np.max(np.abs(gram_o2 - np.eye(3)))),
            "gram_o4_max_error": float(np.max(np.abs(gram_o4 - np.eye(5)))),
            "product_reconstruction_max_error": pmap["max_product_reconstruction_error"],
            "veronese_relation_max_abs": pmap["veronese_relation_max_abs"],
        },
        "dual_completion": dual,
        "toeplitz_scan": toeplitz,
        "verdict": (
            "The legacy holomorphic matrices can be made coordinate-covariant only by "
            "declaring the Higgs functional to be an O(-4)-valued dual density.  This "
            "preserves the same Veronese relation and hierarchy scan.  A scalar "
            "Hermitian Toeplitz kernel is manifestly invariant and can generate "
            "hierarchies, but it belongs naturally to a Kahler/wavefunction texture "
            "unless a separate holomorphic superpotential map is supplied."
        ),
    }
    write_toeplitz_csv(toeplitz)
    (OUT / "cp1_yukawa_covariant_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(payload)
    print("CP1 O(2) covariant Yukawa audit complete")
    print(f"  O2 gram error={payload['bundle_checks']['gram_o2_max_error']:.3e}")
    print(f"  O4 gram error={payload['bundle_checks']['gram_o4_max_error']:.3e}")
    print(f"  product reconstruction={payload['bundle_checks']['product_reconstruction_max_error']:.3e}")
    for name, item in toeplitz["results"].items():
        norm = item["normalized_singular_values"]
        print(f"  Toeplitz {name}: small={norm[2]:.3e} mid={norm[1]:.3e} score={item['score']:.3f}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
