#!/usr/bin/env python3
"""Audit hidden-quotient attempts to derive the Majorana contact coefficient.

No web lookup is used.  The previous affine source-locking sector can lock the
PMNS-sensitive contact coefficient zeta, but it inserts p0 and q0 as target
constants.  This script tests whether the obvious hidden quotient/clockwork
upgrade can remove those constants.

Three local mechanisms are compared at the same target point:

1. U(1) D-flat quotient for P,Q with charges (+1,-1).
2. D-flat quotient plus a holomorphic product constraint W=X(PQ-Z).
3. The affine source fallback W=A(Z-PQ)+B(P-p0)+C(Q-q0).

The audit computes real Hessian ranks for the first two potentials and imports
the affine Hessian from the previous source-locking audit.  The point is not to
claim a new model, but to distinguish real prediction from hidden spurions.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "majorana_hidden_quotient_origin"
SOURCE_LOCK = ROOT / "output" / "majorana_source_locking_sector" / "summary.json"
SENSITIVITY = ROOT / "output" / "majorana_contact_sensitivity" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cnum(raw: dict[str, float]) -> complex:
    return complex(raw["re"], raw["im"])


def real_hessian(fn: Callable[[np.ndarray], float], x0: np.ndarray, step: float = 1.0e-5) -> np.ndarray:
    n = len(x0)
    h = np.zeros((n, n), dtype=float)
    f0 = fn(x0)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = step
        h[i, i] = (fn(x0 + ei) - 2.0 * f0 + fn(x0 - ei)) / (step * step)
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = step
            h[i, j] = (
                fn(x0 + ei + ej)
                - fn(x0 + ei - ej)
                - fn(x0 - ei + ej)
                + fn(x0 - ei - ej)
            ) / (4.0 * step * step)
            h[j, i] = h[i, j]
    return h


def rank_from_hessian(h: np.ndarray, tol: float = 1.0e-8) -> tuple[int, list[float]]:
    eig = np.linalg.eigvalsh(0.5 * (h + h.T))
    rank = int(np.sum(np.abs(eig) > tol))
    return rank, [float(x) for x in eig]


def unpack_complex(vec: np.ndarray) -> list[complex]:
    return [vec[2 * i] + 1j * vec[2 * i + 1] for i in range(len(vec) // 2)]


def d_only_potential(vec: np.ndarray) -> float:
    p, q = unpack_complex(vec)
    d = abs(p) ** 2 - abs(q) ** 2
    return 0.5 * d * d


def product_potential(vec: np.ndarray) -> float:
    p, q, z, x = unpack_complex(vec)
    d = abs(p) ** 2 - abs(q) ** 2
    f_x = p * q - z
    f_p = x * q
    f_q = x * p
    f_z = -x
    return abs(f_x) ** 2 + abs(f_p) ** 2 + abs(f_q) ** 2 + abs(f_z) ** 2 + 0.5 * d * d


def build() -> dict[str, Any]:
    source = read_json(SOURCE_LOCK)
    sensitivity = read_json(SENSITIVITY)
    target = source["target"]
    zeta = cnum(target["zeta"])
    p0 = cnum(target["p0"])
    q0 = cnum(target["q0"])
    d_x0 = np.array([p0.real, p0.imag, q0.real, q0.imag], dtype=float)
    product_x0 = np.array(
        [p0.real, p0.imag, q0.real, q0.imag, zeta.real, zeta.imag, 0.0, 0.0],
        dtype=float,
    )
    d_h = real_hessian(d_only_potential, d_x0)
    product_h = real_hessian(product_potential, product_x0)
    d_rank, d_eigs = rank_from_hessian(d_h)
    product_rank, product_eigs = rank_from_hessian(product_h)
    d_flat = len(d_x0) - d_rank
    product_flat = len(product_x0) - product_rank
    loose_scale = sensitivity["real_scale_loose_interval"]["half_width"]
    loose_phase = sensitivity["phase_loose_interval"]["half_width"]
    summary = {
        "note": "No web lookup used. Hidden quotient origin audit for Majorana contact zeta.",
        "input_manifest": [
            {
                "label": "majorana_source_locking_sector",
                "path": str(SOURCE_LOCK.relative_to(ROOT)),
                "sha256": sha256(SOURCE_LOCK),
            },
            {
                "label": "majorana_contact_sensitivity",
                "path": str(SENSITIVITY.relative_to(ROOT)),
                "sha256": sha256(SENSITIVITY),
            },
        ],
        "target": {
            "zeta_abs": float(abs(zeta)),
            "zeta_arg_rad": float(math.atan2(zeta.imag, zeta.real)),
            "p0_abs": float(abs(p0)),
            "q0_abs": float(abs(q0)),
            "p0q0_minus_zeta_abs": float(abs(p0 * q0 - zeta)),
        },
        "pmns_locking_requirement": {
            "loose_scale_half_width": float(loose_scale),
            "loose_phase_half_width_rad": float(loose_phase),
        },
        "candidates": {
            "D_only_hidden_U1": {
                "potential": "V=1/2(|P|^2-|Q|^2)^2",
                "real_dimension": len(d_x0),
                "hessian_rank": d_rank,
                "flat_real_directions_before_gauge_quotient": d_flat,
                "hessian_eigenvalues": d_eigs,
                "predicts_abs_zeta": False,
                "predicts_arg_zeta": False,
                "interpretation": (
                    "The D-term fixes only the radial difference.  The radial sum "
                    "and gauge-invariant product phase/magnitude remain moduli."
                ),
            },
            "D_plus_product_constraint": {
                "potential": "V=|PQ-Z|^2+|XQ|^2+|XP|^2+|X|^2+1/2D^2",
                "superpotential": "W=X(PQ-Z)",
                "real_dimension": len(product_x0),
                "hessian_rank": product_rank,
                "flat_real_directions_before_gauge_quotient": product_flat,
                "hessian_eigenvalues": product_eigs,
                "predicts_abs_zeta": False,
                "predicts_arg_zeta": False,
                "interpretation": (
                    "The product constraint transfers the modulus into Z: it enforces "
                    "Z=PQ but does not choose the numerical complex value of PQ."
                ),
            },
            "affine_source_fallback": {
                "superpotential": source["superpotential"],
                "f_norm_at_vacuum": source["f_norm_at_vacuum"],
                "hessian_rank": source["hessian_rank"],
                "hessian_singular_value_floor": float(min(source["hessian_singular_values"])),
                "visible_threshold_vector": source["visible_threshold_vector"],
                "predicts_abs_zeta": False,
                "predicts_arg_zeta": False,
                "interpretation": (
                    "The affine fallback removes local moduli but only because p0 and q0 "
                    "are inserted target constants."
                ),
            },
        },
        "verdict": {
            "pure_quotient_predicts_zeta": False,
            "affine_fallback_locks_but_does_not_predict": True,
            "next_needed_mechanism": "A hidden dynamics that fixes p0 and q0, or a declared conditional EFT datum.",
            "interpretation": (
                "The minimal hidden quotient route does not predict the Majorana contact "
                "coefficient.  D-flatness and the product constraint leave moduli, while "
                "the affine sector removes them by inserting p0 and q0.  Therefore the "
                "first-principles bottleneck is now the microscopic origin of the hidden "
                "source constants, not the algebraic source lock."
            ),
        },
    }
    return summary


def report(summary: dict[str, Any]) -> str:
    d = summary["candidates"]["D_only_hidden_U1"]
    p = summary["candidates"]["D_plus_product_constraint"]
    a = summary["candidates"]["affine_source_fallback"]
    t = summary["target"]
    req = summary["pmns_locking_requirement"]
    return "\n".join(
        [
            "# Majorana Hidden-Quotient Origin Audit",
            "",
            "No web lookup was used.  This audit checks whether a minimal hidden quotient can replace the affine source constants p0,q0.",
            "",
            "## Target and required precision",
            "",
            f"- |zeta|: `{t['zeta_abs']:.6e}`",
            f"- arg(zeta): `{t['zeta_arg_rad']:.6e}` rad",
            f"- loose scale half-width: `{req['loose_scale_half_width']:.6e}`",
            f"- loose phase half-width: `{req['loose_phase_half_width_rad']:.6e}` rad",
            "",
            "## Candidate A: D-only hidden U(1)",
            "",
            f"- Hessian rank: `{d['hessian_rank']}/{d['real_dimension']}`",
            f"- flat real directions before gauge quotient: `{d['flat_real_directions_before_gauge_quotient']}`",
            f"- eigenvalues: `{', '.join(f'{x:.6e}' for x in d['hessian_eigenvalues'])}`",
            f"- predicts |zeta|: `{d['predicts_abs_zeta']}`",
            f"- predicts arg(zeta): `{d['predicts_arg_zeta']}`",
            "",
            "## Candidate B: D plus product constraint",
            "",
            f"- Hessian rank: `{p['hessian_rank']}/{p['real_dimension']}`",
            f"- flat real directions before gauge quotient: `{p['flat_real_directions_before_gauge_quotient']}`",
            f"- eigenvalues: `{', '.join(f'{x:.6e}' for x in p['hessian_eigenvalues'])}`",
            f"- predicts |zeta|: `{p['predicts_abs_zeta']}`",
            f"- predicts arg(zeta): `{p['predicts_arg_zeta']}`",
            "",
            "## Candidate C: affine source fallback",
            "",
            f"- F-term norm: `{a['f_norm_at_vacuum']:.6e}`",
            f"- Hessian rank: `{a['hessian_rank']}/6`",
            f"- singular-value floor: `{a['hessian_singular_value_floor']:.6e}`",
            f"- visible threshold vector: `{a['visible_threshold_vector']}`",
            f"- predicts |zeta|: `{a['predicts_abs_zeta']}`",
            f"- predicts arg(zeta): `{a['predicts_arg_zeta']}`",
            "",
            "## Verdict",
            "",
            summary["verdict"]["interpretation"],
            "",
        ]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    d = summary["candidates"]["D_only_hidden_U1"]
    p = summary["candidates"]["D_plus_product_constraint"]
    print("Majorana hidden-quotient origin audit")
    print(f"  D-only rank: {d['hessian_rank']}/{d['real_dimension']}")
    print(f"  D+product rank: {p['hessian_rank']}/{p['real_dimension']}")
    print(f"  pure quotient predicts zeta: {summary['verdict']['pure_quotient_predicts_zeta']}")


if __name__ == "__main__":
    main()
