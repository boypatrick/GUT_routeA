#!/usr/bin/env python3
"""Audit a minimal holomorphic source-locking sector for the Majorana contact.

No web lookup is used.  This is an algebraic stress test, not a microscopic
derivation.  The previous scan showed that the complex contact coefficient
zeta must be fixed at roughly the 1e-5 level.  Here we ask a narrower question:

Can a renormalizable holomorphic driving sector lock zeta without producing
flat directions or visible gauge thresholds?

We use the affine hidden-source superpotential

    W = A (Z - P Q) + B (P - p0) + C (Q - q0),

where P,Q,Z,A,B,C are gauge-singlet chiral fields and p0 q0 = zeta_target.
At the SUSY vacuum

    Z=zeta_target, P=p0, Q=q0, A=B=C=0,

all F terms vanish.  The Hessian is then inspected to determine whether the
locking sector has moduli.  Because p0 and q0 are inserted target constants,
this sector is a tuned/conditional source lock rather than a first-principles
prediction of zeta.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "majorana_source_locking_sector"
RANK_SUMMARY = ROOT / "output" / "source_majorana_texture_rank" / "summary.json"
SENSITIVITY = ROOT / "output" / "majorana_contact_sensitivity" / "summary.json"


FIELDS = ["Z", "P", "Q", "A", "B", "C"]


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


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(x) for x in row] for row in m]


def f_terms(phi: dict[str, complex], p0: complex, q0: complex) -> dict[str, complex]:
    z, p, q, a, b, c = (phi[name] for name in FIELDS)
    return {
        "F_Z": a,
        "F_P": -a * q + b,
        "F_Q": -a * p + c,
        "F_A": z - p * q,
        "F_B": p - p0,
        "F_C": q - q0,
    }


def hessian_at_vacuum(p0: complex, q0: complex) -> np.ndarray:
    # Field order: Z, P, Q, A, B, C.
    h = np.zeros((6, 6), dtype=complex)
    idx = {name: i for i, name in enumerate(FIELDS)}

    def set_sym(a: str, b: str, value: complex) -> None:
        h[idx[a], idx[b]] = value
        h[idx[b], idx[a]] = value

    set_sym("Z", "A", 1.0)
    set_sym("P", "A", -q0)
    set_sym("Q", "A", -p0)
    set_sym("P", "B", 1.0)
    set_sym("Q", "C", 1.0)
    return h


def build() -> dict[str, Any]:
    rank = read_json(RANK_SUMMARY)
    sensitivity = read_json(SENSITIVITY)
    zeta = cnum(rank["decomposition"]["contact_coefficient"])
    rho = abs(zeta)
    phase = math.atan2(zeta.imag, zeta.real)
    # A balanced hidden factorization p0 q0 = zeta avoids a large hierarchy in
    # the source fields.  This is a convention, not an extra prediction.
    p0 = math.sqrt(rho)
    q0 = math.sqrt(rho) * np.exp(1j * phase)
    vacuum = {
        "Z": zeta,
        "P": p0,
        "Q": q0,
        "A": 0.0 + 0.0j,
        "B": 0.0 + 0.0j,
        "C": 0.0 + 0.0j,
    }
    f = f_terms(vacuum, p0, q0)
    f_norm = float(math.sqrt(sum(abs(value) ** 2 for value in f.values())))
    h = hessian_at_vacuum(p0, q0)
    singulars = np.linalg.svd(h, compute_uv=False)
    rank_h = int(np.linalg.matrix_rank(h, tol=1.0e-12))
    mass_sq = singulars**2
    loose_s_half = sensitivity["real_scale_loose_interval"]["half_width"]
    loose_phi_half = sensitivity["phase_loose_interval"]["half_width"]
    tight_s_half = sensitivity["real_scale_tight_interval"]["half_width"]
    tight_phi_half = sensitivity["phase_tight_interval"]["half_width"]
    summary = {
        "note": "No web lookup used. Minimal affine Majorana source-locking sector audit.",
        "input_manifest": [
            {
                "label": "source_majorana_texture_rank",
                "path": str(RANK_SUMMARY.relative_to(ROOT)),
                "sha256": sha256(RANK_SUMMARY),
            },
            {
                "label": "majorana_contact_sensitivity",
                "path": str(SENSITIVITY.relative_to(ROOT)),
                "sha256": sha256(SENSITIVITY),
            },
        ],
        "superpotential": "W = A (Z - P Q) + B (P - p0) + C (Q - q0)",
        "field_order": FIELDS,
        "target": {
            "zeta": cjson(zeta),
            "abs_zeta": float(rho),
            "arg_zeta_rad": float(phase),
            "p0": cjson(p0),
            "q0": cjson(q0),
            "p0_times_q0_minus_zeta_abs": float(abs(p0 * q0 - zeta)),
        },
        "f_terms_at_vacuum": {key: cjson(value) for key, value in f.items()},
        "f_norm_at_vacuum": f_norm,
        "hessian": matrix_json(h),
        "hessian_singular_values": [float(x) for x in singulars],
        "hessian_rank": rank_h,
        "scalar_mass_squared_eigenvalues_canonical_susy": [float(x) for x in mass_sq],
        "source_locking_requirements_from_pmns_scan": {
            "loose_scale_half_width": float(loose_s_half),
            "loose_phase_half_width_rad": float(loose_phi_half),
            "tight_scale_half_width": float(tight_s_half),
            "tight_phase_half_width_rad": float(tight_phi_half),
        },
        "visible_threshold_vector": [0.0, 0.0, 0.0],
        "verdict": {
            "F_flat": f_norm < 1.0e-12,
            "hessian_full_rank": rank_h == len(FIELDS),
            "visible_threshold_safe_if_hidden_singlets": True,
            "zeta_value_predicted": False,
            "interpretation": (
                "The affine source sector can lock the reconstructed complex zeta "
                "without flat directions and without visible gauge thresholds if all "
                "new fields are hidden singlets.  However, p0 and q0 encode the target "
                "zeta, so this is a controlled conditional source lock, not a "
                "first-principles derivation of the PMNS Majorana texture."
            ),
        },
    }
    return summary


def report(summary: dict[str, Any]) -> str:
    target = summary["target"]
    req = summary["source_locking_requirements_from_pmns_scan"]
    v = summary["verdict"]
    return "\n".join(
        [
            "# Majorana Source-Locking Sector Audit",
            "",
            "No web lookup was used.  This audit checks a minimal holomorphic driving sector for the contact coefficient.",
            "",
            "## Superpotential",
            "",
            "`W = A (Z - P Q) + B (P - p0) + C (Q - q0)`",
            "",
            "At the vacuum, `Z=zeta`, `P=p0`, `Q=q0`, `A=B=C=0`.",
            "",
            "## Target",
            "",
            f"- |zeta|: `{target['abs_zeta']:.6e}`",
            f"- arg(zeta): `{target['arg_zeta_rad']:.6e}` rad",
            f"- |p0 q0 - zeta|: `{target['p0_times_q0_minus_zeta_abs']:.6e}`",
            "",
            "## F-flatness and Hessian",
            "",
            f"- F-term norm at vacuum: `{summary['f_norm_at_vacuum']:.6e}`",
            f"- Hessian rank: `{summary['hessian_rank']}/{len(FIELDS)}`",
            f"- Hessian singular values: `{', '.join(f'{x:.6e}' for x in summary['hessian_singular_values'])}`",
            f"- visible threshold vector: `{summary['visible_threshold_vector']}`",
            "",
            "## Required locking precision inherited from PMNS sensitivity",
            "",
            f"- loose scale half-width: `{req['loose_scale_half_width']:.6e}`",
            f"- loose phase half-width: `{req['loose_phase_half_width_rad']:.6e}` rad",
            f"- tight scale half-width: `{req['tight_scale_half_width']:.6e}`",
            f"- tight phase half-width: `{req['tight_phase_half_width_rad']:.6e}` rad",
            "",
            "## Verdict",
            "",
            f"- F-flat: `{v['F_flat']}`",
            f"- Hessian full rank: `{v['hessian_full_rank']}`",
            f"- zeta value predicted: `{v['zeta_value_predicted']}`",
            "",
            v["interpretation"],
            "",
        ]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    print("Majorana source-locking sector audit")
    print(f"  F norm: {summary['f_norm_at_vacuum']:.6e}")
    print(f"  Hessian rank: {summary['hessian_rank']}/{len(FIELDS)}")
    print(f"  min singular: {min(summary['hessian_singular_values']):.6e}")
    print(f"  zeta predicted: {summary['verdict']['zeta_value_predicted']}")


if __name__ == "__main__":
    main()
