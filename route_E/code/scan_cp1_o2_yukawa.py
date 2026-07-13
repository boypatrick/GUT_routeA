#!/usr/bin/env python3
"""CP1 O(2) Yukawa overlap texture scanner.

This script implements the item-2 toy model:

    Y_ij^(a) = lambda_a int_CP1 psi_i psi_j h_a dmu_FS,

where psi_i are the three normalized holomorphic sections of O(2) on CP1.

The Higgs profile h_a is modeled as a two-center PSLT/Ed kernel with a small
curvature-leakage term.  The scan asks whether hierarchical singular values can
appear from geometry without inserting a hand-made Yukawa matrix.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIELD_TABLE = ROOT / "output" / "pati_salam" / "pati_salam_decomposition.csv"
OUT = ROOT / "output" / "yukawa_o2"


@dataclass(frozen=True)
class SectorTarget:
    name: str
    left_field: str
    right_field: str
    higgs_label: str
    target_small_over_large: float
    target_mid_over_large: float


@dataclass(frozen=True)
class KernelParams:
    theta1: float
    phi1: float
    kappa1: float
    theta2: float
    phi2: float
    kappa2: float
    amp2: float
    phase2: float
    eps: float
    phase_eps: float
    quad: float
    phase_quad: float


def require_verified_fields() -> None:
    required = {"Q", "u^c", "d^c", "L", "e^c", "nu^c"}
    with FIELD_TABLE.open(newline="", encoding="utf-8") as handle:
        seen = {row["name"] for row in csv.DictReader(handle)}
    missing = sorted(required - seen)
    if missing:
        raise RuntimeError(f"Missing Pati-Salam fields from item 1 output: {missing}")


def unit_vector(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ],
        dtype=float,
    )


def angular_separation(a_theta: float, a_phi: float, b_theta: float, b_phi: float) -> float:
    a = unit_vector(a_theta, a_phi)
    b = unit_vector(b_theta, b_phi)
    return float(math.acos(np.clip(np.dot(a, b), -1.0, 1.0)))


class CP1O2Grid:
    def __init__(self, n_theta: int = 54, n_phi: int = 108) -> None:
        # Gauss-Legendre in u=cos(theta), uniform in phi.  The normalized measure is
        # dmu = du dphi/(4*pi), so the phi weight is 1/(2*n_phi).
        u, wu = np.polynomial.legendre.leggauss(n_theta)
        phi = (2.0 * np.pi / n_phi) * np.arange(n_phi)
        uu, pp = np.meshgrid(u, phi, indexing="ij")
        ww = np.repeat(wu[:, None], n_phi, axis=1) / (2.0 * n_phi)

        theta = np.arccos(uu)
        half_c = np.sqrt((1.0 + uu) / 2.0)
        half_s = np.sqrt((1.0 - uu) / 2.0)

        self.u = uu.reshape(-1)
        self.phi = pp.reshape(-1)
        self.theta = theta.reshape(-1)
        self.weight = ww.reshape(-1)
        self.nx = (np.sin(theta) * np.cos(pp)).reshape(-1)
        self.ny = (np.sin(theta) * np.sin(pp)).reshape(-1)
        self.nz = uu.reshape(-1)
        self.points = np.stack([self.nx, self.ny, self.nz], axis=0)

        # Normalized O(2) sections:
        # psi_m = sqrt(3*C(2,m)) exp(i*m*phi) sin(theta/2)^m cos(theta/2)^(2-m).
        psi0 = math.sqrt(3.0) * half_c**2
        psi1 = math.sqrt(6.0) * np.exp(1j * pp) * half_s * half_c
        psi2 = math.sqrt(3.0) * np.exp(2j * pp) * half_s**2
        self.psi = np.stack([psi0.reshape(-1), psi1.reshape(-1), psi2.reshape(-1)], axis=0)

        pair = []
        for i in range(3):
            for j in range(3):
                pair.append(self.weight * self.psi[i] * self.psi[j])
        self.pair = np.stack(pair, axis=0)

    def gram(self) -> np.ndarray:
        return (self.psi.conjugate() * self.weight) @ self.psi.T

    def y_matrix(self, params: KernelParams) -> np.ndarray:
        h = self.higgs_profile(params)
        y = (self.pair @ h).reshape(3, 3)
        return 0.5 * (y + y.T)

    def higgs_profile(self, params: KernelParams) -> np.ndarray:
        c1 = unit_vector(params.theta1, params.phi1)
        c2 = unit_vector(params.theta2, params.phi2)
        dot1 = c1 @ self.points
        dot2 = c2 @ self.points

        # exp(kappa*(n.c-1)) has peak 1 and avoids numerical overflow.
        k1 = np.exp(params.kappa1 * (dot1 - 1.0))
        k2 = np.exp(params.kappa2 * (dot2 - 1.0))
        q2 = 0.5 * (3.0 * self.nz**2 - 1.0)

        return (
            k1
            + params.amp2 * np.exp(1j * params.phase2) * k2
            + params.eps * np.exp(1j * params.phase_eps)
            + params.quad * np.exp(1j * params.phase_quad) * q2
        )


def singular_data(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.linalg.svd(y, compute_uv=False)
    values = np.sort(values)[::-1]
    normalized = values / values[0]
    return values, normalized


def log_score(norm_values: np.ndarray, target_small: float, target_mid: float) -> float:
    small = max(float(norm_values[2]), 1e-16)
    mid = max(float(norm_values[1]), 1e-16)
    return float((math.log10(small / target_small)) ** 2 + (math.log10(mid / target_mid)) ** 2)


def sample_params(rng: np.random.Generator) -> KernelParams:
    theta1 = rng.uniform(0.08, math.pi - 0.08)
    phi1 = rng.uniform(0.0, 2.0 * math.pi)

    # Draw the second center independently, but reject near-identical centers.
    for _ in range(100):
        theta2 = rng.uniform(0.08, math.pi - 0.08)
        phi2 = rng.uniform(0.0, 2.0 * math.pi)
        if angular_separation(theta1, phi1, theta2, phi2) > 0.25:
            break

    return KernelParams(
        theta1=theta1,
        phi1=phi1,
        kappa1=float(10 ** rng.uniform(math.log10(8.0), math.log10(85.0))),
        theta2=theta2,
        phi2=phi2,
        kappa2=float(10 ** rng.uniform(math.log10(8.0), math.log10(85.0))),
        amp2=float(10 ** rng.uniform(-2.2, 0.2)),
        phase2=rng.uniform(0.0, 2.0 * math.pi),
        eps=float(10 ** rng.uniform(-6.0, -2.2)),
        phase_eps=rng.uniform(0.0, 2.0 * math.pi),
        quad=float(10 ** rng.uniform(-6.0, -2.0)),
        phase_quad=rng.uniform(0.0, 2.0 * math.pi),
    )


def perturb_params(params: KernelParams, rng: np.random.Generator) -> KernelParams:
    def wrap_phi(x: float) -> float:
        return float(x % (2.0 * math.pi))

    def clip_theta(x: float) -> float:
        return float(np.clip(x, 0.03, math.pi - 0.03))

    def log_perturb(x: float, sigma: float, lo: float, hi: float) -> float:
        return float(np.clip(x * 10 ** rng.normal(0.0, sigma), lo, hi))

    return KernelParams(
        theta1=clip_theta(params.theta1 + rng.normal(0.0, 0.018)),
        phi1=wrap_phi(params.phi1 + rng.normal(0.0, 0.018)),
        kappa1=log_perturb(params.kappa1, 0.025, 4.0, 120.0),
        theta2=clip_theta(params.theta2 + rng.normal(0.0, 0.018)),
        phi2=wrap_phi(params.phi2 + rng.normal(0.0, 0.018)),
        kappa2=log_perturb(params.kappa2, 0.025, 4.0, 120.0),
        amp2=log_perturb(params.amp2, 0.04, 1e-4, 3.0),
        phase2=wrap_phi(params.phase2 + rng.normal(0.0, 0.035)),
        eps=log_perturb(params.eps, 0.05, 1e-8, 1e-1),
        phase_eps=wrap_phi(params.phase_eps + rng.normal(0.0, 0.04)),
        quad=log_perturb(params.quad, 0.05, 1e-8, 1e-1),
        phase_quad=wrap_phi(params.phase_quad + rng.normal(0.0, 0.04)),
    )


def stability_summary(grid: CP1O2Grid, params: KernelParams, target: SectorTarget, rng: np.random.Generator, n: int = 180) -> dict[str, object]:
    rows = []
    scores = []
    for _ in range(n):
        p = perturb_params(params, rng)
        _, norm = singular_data(grid.y_matrix(p))
        score = log_score(norm, target.target_small_over_large, target.target_mid_over_large)
        rows.append([float(norm[2]), float(norm[1]), score])
        scores.append(score)
    arr = np.array(rows)
    return {
        "perturbations": n,
        "small_median": float(np.median(arr[:, 0])),
        "small_log10_iqr": float(np.quantile(np.log10(arr[:, 0]), 0.75) - np.quantile(np.log10(arr[:, 0]), 0.25)),
        "mid_median": float(np.median(arr[:, 1])),
        "mid_log10_iqr": float(np.quantile(np.log10(arr[:, 1]), 0.75) - np.quantile(np.log10(arr[:, 1]), 0.25)),
        "score_median": float(np.median(scores)),
        "score_p90": float(np.quantile(scores, 0.90)),
    }


def params_to_json(params: KernelParams) -> dict[str, float]:
    data = asdict(params)
    data["center_separation"] = angular_separation(params.theta1, params.phi1, params.theta2, params.phi2)
    return data


def rank_staircase(grid: CP1O2Grid, params: KernelParams) -> dict[str, list[float]]:
    single = KernelParams(
        theta1=params.theta1,
        phi1=params.phi1,
        kappa1=params.kappa1,
        theta2=params.theta2,
        phi2=params.phi2,
        kappa2=params.kappa2,
        amp2=0.0,
        phase2=params.phase2,
        eps=0.0,
        phase_eps=params.phase_eps,
        quad=0.0,
        phase_quad=params.phase_quad,
    )
    double = KernelParams(
        theta1=params.theta1,
        phi1=params.phi1,
        kappa1=params.kappa1,
        theta2=params.theta2,
        phi2=params.phi2,
        kappa2=params.kappa2,
        amp2=params.amp2,
        phase2=params.phase2,
        eps=0.0,
        phase_eps=params.phase_eps,
        quad=0.0,
        phase_quad=params.phase_quad,
    )
    out = {}
    for label, p in [("single_center", single), ("two_centers", double), ("full_kernel", params)]:
        _, norm = singular_data(grid.y_matrix(p))
        out[label] = [float(norm[2]), float(norm[1]), float(norm[0])]
    return out


def matrix_to_json(y: np.ndarray) -> list[list[dict[str, float]]]:
    return [[{"re": float(z.real), "im": float(z.imag)} for z in row] for row in y]


def matrix_to_markdown(y: np.ndarray) -> str:
    lines = []
    for row in y:
        pieces = []
        for z in row:
            pieces.append(f"{z.real:+.4e}{z.imag:+.4e}i")
        lines.append("[ " + ", ".join(pieces) + " ]")
    return "\n".join(lines)


def scan_sector(grid: CP1O2Grid, target: SectorTarget, rng: np.random.Generator, samples: int, keep: int = 36) -> dict[str, object]:
    best: list[tuple[float, KernelParams, np.ndarray, np.ndarray]] = []
    for _ in range(samples):
        params = sample_params(rng)
        y = grid.y_matrix(params)
        singular_values, norm = singular_data(y)
        score = log_score(norm, target.target_small_over_large, target.target_mid_over_large)
        if len(best) < keep:
            best.append((score, params, singular_values, y))
            best.sort(key=lambda item: item[0])
        elif score < best[-1][0]:
            best[-1] = (score, params, singular_values, y)
            best.sort(key=lambda item: item[0])

    stability_seeds = {
        "up": 1101,
        "down": 1102,
        "charged_lepton": 1103,
        "neutrino_dirac": 1104,
    }
    stability_rng = np.random.default_rng(stability_seeds[target.name])
    reranked = []
    for score, params, singular_values, y in best[:12]:
        stability = stability_summary(grid, params, target, stability_rng)
        combined = score + 0.65 * float(stability["small_log10_iqr"]) + 0.35 * float(stability["mid_log10_iqr"])
        reranked.append((combined, score, params, singular_values, y, stability))
    reranked.sort(key=lambda item: item[0])

    combined, score, params, singular_values, y, stability = reranked[0]
    normalized = singular_values / singular_values[0]
    return {
        "target": asdict(target),
        "score": float(score),
        "combined_score": float(combined),
        "params": params_to_json(params),
        "singular_values": [float(x) for x in singular_values],
        "normalized_singular_values": [float(x) for x in normalized],
        "matrix": matrix_to_json(y / singular_values[0]),
        "matrix_markdown": matrix_to_markdown(y / singular_values[0]),
        "stability": stability,
        "rank_staircase": rank_staircase(grid, params),
    }


def write_report(results: dict[str, object], gram: np.ndarray, targets: Iterable[SectorTarget], samples: int) -> None:
    lines: list[str] = []
    lines.append("# CP1 O(2) Yukawa texture scan")
    lines.append("")
    lines.append("Family basis: normalized holomorphic sections of `O(2)` on `CP^1`:")
    lines.append("")
    lines.append("```text")
    lines.append("psi_0 = sqrt(3) cos^2(theta/2)")
    lines.append("psi_1 = sqrt(6) exp(i phi) sin(theta/2) cos(theta/2)")
    lines.append("psi_2 = sqrt(3) exp(2 i phi) sin^2(theta/2)")
    lines.append("dmu_FS = sin(theta) dtheta dphi/(4 pi)")
    lines.append("```")
    lines.append("")
    lines.append("Yukawa ansatz:")
    lines.append("")
    lines.append("```text")
    lines.append("Y_ij^(a) = lambda_a int_CP1 psi_i psi_j h_a(theta,phi; Ed) dmu_FS")
    lines.append("h_a = K_1 + A exp(i chi) K_2 + eps exp(i chi_eps) + q exp(i chi_q) P_2(cos theta)")
    lines.append("K_A = exp[kappa_A (n dot n_A - 1)]")
    lines.append("```")
    lines.append("")
    lines.append(f"Scan samples per sector: `{samples}`.")
    lines.append("")
    lines.append("Orthonormality check `int psi_i^* psi_j dmu_FS`:")
    lines.append("")
    for row in gram:
        lines.append("`" + "  ".join(f"{z.real:+.8f}{z.imag:+.1e}i" for z in row) + "`")
    lines.append("")
    lines.append("## Best sector candidates")
    lines.append("")
    for target in targets:
        item = results[target.name]
        norm = item["normalized_singular_values"]
        stability = item["stability"]
        params = item["params"]
        lines.append(f"### {target.name}")
        lines.append("")
        lines.append(f"Coupling: `{target.left_field} {target.right_field} {target.higgs_label}`.")
        lines.append(
            "Target normalized singulars: "
            f"`[{target.target_small_over_large:.3e}, {target.target_mid_over_large:.3e}, 1]`."
        )
        lines.append(
            "Found normalized singulars: "
            f"`[{norm[2]:.3e}, {norm[1]:.3e}, {norm[0]:.3e}]`."
        )
        lines.append(
            "Stability medians under perturbations: "
            f"small `{stability['small_median']:.3e}`, mid `{stability['mid_median']:.3e}`; "
            f"log10-IQRs `{stability['small_log10_iqr']:.3f}`, `{stability['mid_log10_iqr']:.3f}`."
        )
        lines.append(
            "Geometry: "
            f"separation `{params['center_separation']:.3f}`, "
            f"kappas `({params['kappa1']:.2f}, {params['kappa2']:.2f})`, "
            f"amp2 `{params['amp2']:.3e}`, eps `{params['eps']:.3e}`, quad `{params['quad']:.3e}`."
        )
        lines.append("Rank staircase `[small, mid, large]`:")
        for label, values in item["rank_staircase"].items():
            lines.append(f"- `{label}`: `[{values[0]:.3e}, {values[1]:.3e}, {values[2]:.3e}]`")
        lines.append("")
        lines.append("Normalized matrix:")
        lines.append("")
        lines.append("```text")
        lines.append(item["matrix_markdown"])
        lines.append("```")
        lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append("- A single infinitely localized center gives an outer product `v v^T`, hence rank one.")
    lines.append("- A second center raises the generic rank to two; finite-width and curvature-leakage terms generate the smallest singular value.")
    lines.append("- Sector-dependent Higgs kernels give misaligned symmetric matrices while keeping the same topological three-family basis.")
    lines.append("- The scan is a hierarchy existence proof for the toy geometry, not yet a full CKM/PMNS fit.")
    (OUT / "cp1_o2_yukawa_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_candidate_csv(results: dict[str, object], targets: Iterable[SectorTarget]) -> None:
    rows = []
    for target in targets:
        item = results[target.name]
        params = item["params"]
        norm = item["normalized_singular_values"]
        stability = item["stability"]
        rows.append(
            {
                "sector": target.name,
                "left_field": target.left_field,
                "right_field": target.right_field,
                "higgs_label": target.higgs_label,
                "target_small": target.target_small_over_large,
                "target_mid": target.target_mid_over_large,
                "found_small": norm[2],
                "found_mid": norm[1],
                "score": item["score"],
                "combined_score": item["combined_score"],
                "center_separation": params["center_separation"],
                "kappa1": params["kappa1"],
                "kappa2": params["kappa2"],
                "amp2": params["amp2"],
                "eps": params["eps"],
                "quad": params["quad"],
                "small_log10_iqr": stability["small_log10_iqr"],
                "mid_log10_iqr": stability["mid_log10_iqr"],
            }
        )
    with (OUT / "cp1_o2_yukawa_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    require_verified_fields()
    OUT.mkdir(parents=True, exist_ok=True)

    samples = 9000
    rng = np.random.default_rng(20260505)
    grid = CP1O2Grid(n_theta=54, n_phi=108)
    gram = grid.gram()
    gram_error = float(np.max(np.abs(gram - np.eye(3))))

    targets = [
        SectorTarget("up", "Q", "u^c", "H_u", 1.0e-5, 7.0e-3),
        SectorTarget("down", "Q", "d^c", "H_d", 1.0e-3, 2.0e-2),
        SectorTarget("charged_lepton", "L", "e^c", "H_d/e", 3.0e-4, 6.0e-2),
        SectorTarget("neutrino_dirac", "L", "nu^c", "H_u/nu", 1.0e-2, 2.0e-1),
    ]

    results: dict[str, object] = {}
    for target in targets:
        results[target.name] = scan_sector(grid, target, rng, samples=samples)

    payload = {
        "model": "CP1 O(2) two-center Higgs-kernel Yukawa scan",
        "samples_per_sector": samples,
        "grid": {
            "n_theta": 54,
            "n_phi": 108,
            "points": 54 * 108,
            "gram_max_abs_error": gram_error,
        },
        "basis": [
            "sqrt(3) cos^2(theta/2)",
            "sqrt(6) exp(i phi) sin(theta/2) cos(theta/2)",
            "sqrt(3) exp(2 i phi) sin^2(theta/2)",
        ],
        "results": results,
    }
    (OUT / "cp1_o2_yukawa_scan.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_candidate_csv(results, targets)
    write_report(results, gram, targets, samples)

    print("CP1 O(2) Yukawa scan")
    print(f"  field table: {FIELD_TABLE}")
    print(f"  grid points: {54 * 108}")
    print(f"  Gram max abs error: {gram_error:.3e}")
    print(f"  samples per sector: {samples}")
    for target in targets:
        item = results[target.name]
        norm = item["normalized_singular_values"]
        stability = item["stability"]
        print(
            "  {name:15s} target=({ts:.1e},{tm:.1e}) found=({fs:.3e},{fm:.3e}) "
            "iqr=({iqs:.3f},{iqm:.3f}) score={score:.3f}".format(
                name=target.name,
                ts=target.target_small_over_large,
                tm=target.target_mid_over_large,
                fs=norm[2],
                fm=norm[1],
                iqs=stability["small_log10_iqr"],
                iqm=stability["mid_log10_iqr"],
                score=item["score"],
            )
        )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
