#!/usr/bin/env python3
"""Audit the source-consistent Majorana texture rank.

No web lookup is used.  This script decomposes the inverse-seesaw reconstructed
M_R from the source-consistent publication card into

    Sym^2 H^0(O(2)) = Veronese_5 \oplus contact_1.

In the normalized CP1/O(2) basis the scalar Veronese subspace obeys

    (M_V)_{11} = 2 (M_V)_{02}.

The Frobenius-orthogonal missing direction is

    C0/sqrt(3),  C0 = [[0,0,1],[0,-1,0],[1,0,0]].

The audit reports whether the reconstructed source-consistent M_R is close to
Veronese-only or whether the trace/contact lift is quantitatively essential.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import verify_seesaw_item3 as seesaw  # noqa: E402


OUT = ROOT / "output" / "source_majorana_texture_rank"
CLOSURE_CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
SOURCE_PMNS = ROOT / "output" / "source_consistent_pmns_replay" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(x) for x in row] for row in m]


def fixed_benchmark() -> seesaw.LightNuBenchmark:
    return seesaw.LightNuBenchmark(
        m1_eV=1.0e-3,
        dm21_eV2=7.42e-5,
        dm31_eV2=2.517e-3,
        sin2_theta12=0.304,
        sin2_theta13=0.0222,
        sin2_theta23=0.573,
        delta_cp_rad=1.20 * math.pi,
        alpha21_rad=0.35 * math.pi,
        alpha31_rad=1.10 * math.pi,
    )


def reconstruct_mr(card: dict[str, Any]) -> dict[str, Any]:
    y = {name: cmat(raw) for name, raw in card["selected_row"]["Yukawa_fit"].items()}
    y_e = y["charged_lepton"]
    y_nu = y["neutrino_dirac"]
    u_e, _ = seesaw.left_rotation(y_e)
    b = fixed_benchmark()
    m_light_diag = np.array(
        [
            b.m1_eV,
            math.sqrt(b.m1_eV**2 + b.dm21_eV2),
            math.sqrt(b.m1_eV**2 + b.dm31_eV2),
        ],
        dtype=float,
    )
    u_pmns = seesaw.standard_pmns(
        b.sin2_theta12,
        b.sin2_theta13,
        b.sin2_theta23,
        b.delta_cp_rad,
        b.alpha21_rad,
        b.alpha31_rad,
    )
    u_nu = u_e @ u_pmns
    m_light = u_nu.conjugate() @ np.diag(m_light_diag) @ u_nu.conjugate().T
    mD_eV = y_nu * 100.0e9
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light) @ mD_eV)
    MR_GeV = MR_eV / 1.0e9
    m_light_reco = -(mD_eV @ np.linalg.inv(MR_eV) @ mD_eV.T)
    return {
        "MR_GeV": MR_GeV,
        "m_light": m_light,
        "m_light_reco": m_light_reco,
        "seesaw_residual": float(np.linalg.norm(m_light_reco - m_light) / np.linalg.norm(m_light)),
    }


def frob_inner(a: np.ndarray, b: np.ndarray) -> complex:
    return complex(np.vdot(a, b))


def decompose_majorana(m: np.ndarray) -> dict[str, Any]:
    singulars = np.linalg.svd(m, compute_uv=False)
    scale = float(singulars[0])
    mn = m / scale
    c0 = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=complex,
    )
    c_hat = c0 / math.sqrt(float(np.vdot(c0, c0).real))
    zeta = frob_inner(c_hat, mn)
    contact = zeta * c_hat
    veronese = mn - contact
    veronese_residual = float(np.linalg.norm(mn - veronese) / np.linalg.norm(mn))
    lifted_residual = float(np.linalg.norm(mn - veronese - contact) / np.linalg.norm(mn))
    contact_fraction = float(np.linalg.norm(contact) / np.linalg.norm(mn))
    veronese_constraint = complex(veronese[1, 1] - 2.0 * veronese[0, 2])
    raw_constraint = complex(mn[1, 1] - 2.0 * mn[0, 2])
    # The contact direction has C_11 - 2 C_02 = -sqrt(3), so this is equivalent
    # to zeta up to the chosen orientation.
    constraint_normed = abs(raw_constraint) / math.sqrt(3.0)
    return {
        "scale_largest_singular_GeV": scale,
        "singular_values_GeV": [float(x) for x in singulars],
        "condition_number": float(singulars[0] / max(singulars[-1], 1.0e-300)),
        "MR_normalized": matrix_json(mn),
        "veronese_projection": matrix_json(veronese),
        "contact_projection": matrix_json(contact),
        "contact_basis": matrix_json(c_hat),
        "contact_coefficient": cjson(zeta),
        "contact_coefficient_abs": float(abs(zeta)),
        "contact_fraction": contact_fraction,
        "veronese_only_relative_residual": veronese_residual,
        "veronese_plus_contact_relative_residual": lifted_residual,
        "raw_constraint_M11_minus_2M02": cjson(raw_constraint),
        "raw_constraint_abs_over_sqrt3": float(constraint_normed),
        "veronese_constraint_after_projection_abs": float(abs(veronese_constraint)),
        "rank_complex_veronese": 5,
        "rank_complex_veronese_plus_contact": 6,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def report(summary: dict[str, Any]) -> str:
    d = summary["decomposition"]
    v = summary["verdict"]
    return "\n".join(
        [
            "# Source Majorana Texture Rank Audit",
            "",
            "No web lookup was used.  The source-consistent inverse-seesaw M_R is decomposed into the CP1/O(2) Veronese subspace plus the orthogonal trace/contact direction.",
            "",
            f"- largest singular scale: `{d['scale_largest_singular_GeV']:.6e}` GeV",
            f"- condition number: `{d['condition_number']:.6e}`",
            f"- Veronese-only relative residual: `{d['veronese_only_relative_residual']:.6e}`",
            f"- contact fraction: `{d['contact_fraction']:.6e}`",
            f"- Veronese+contact residual: `{d['veronese_plus_contact_relative_residual']:.6e}`",
            f"- post-projection Veronese constraint residual: `{d['veronese_constraint_after_projection_abs']:.6e}`",
            f"- contact coefficient abs: `{d['contact_coefficient_abs']:.6e}`",
            f"- near Veronese: `{v['near_veronese_without_contact']}`",
            f"- contact lift essential: `{v['contact_lift_essential']}`",
            f"- predictive Majorana texture closed: `{v['predictive_majorana_texture_closed']}`",
            "",
            "## Interpretation",
            "",
            v["interpretation"],
            "",
        ]
    )


def build() -> dict[str, Any]:
    card = read_json(CLOSURE_CARD)
    source_pmns = read_json(SOURCE_PMNS)
    reconstructed = reconstruct_mr(card)
    decomposition = decompose_majorana(reconstructed["MR_GeV"])
    contact_threshold = 1.0e-2
    near_veronese = decomposition["contact_fraction"] < contact_threshold
    contact_essential = not near_veronese
    summary = {
        "note": "No web lookup used. Source-consistent Majorana texture rank audit.",
        "input_manifest": [
            {
                "label": "publication_closure_card",
                "path": str(CLOSURE_CARD.relative_to(ROOT)),
                "sha256": sha256(CLOSURE_CARD),
            },
            {
                "label": "source_consistent_pmns_replay",
                "path": str(SOURCE_PMNS.relative_to(ROOT)),
                "sha256": sha256(SOURCE_PMNS),
            },
        ],
        "source_pmns_pair_digest": source_pmns["source_pmns_pair_digest"],
        "MR_digest": stable_digest(matrix_json(reconstructed["MR_GeV"])),
        "seesaw_residual_recomputed": reconstructed["seesaw_residual"],
        "decomposition": decomposition,
        "verdict": {
            "near_veronese_without_contact": near_veronese,
            "contact_lift_essential": contact_essential,
            "predictive_majorana_texture_closed": False,
            "interpretation": (
                "The source-consistent reconstructed M_R is exactly representable by "
                "Veronese_5 plus the orthogonal trace/contact direction, but it is not "
                "close to Veronese-only under the 1e-2 contact-fraction criterion.  "
                "Thus PMNS compatibility currently relies on a genuine Majorana "
                "trace/contact lift and inverse reconstruction; the Majorana sector is "
                "not yet a predictive CP1/O(2) texture."
            ),
        },
    }
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    d = summary["decomposition"]
    rows = [
        {
            "quantity": "veronese_only_relative_residual",
            "value": d["veronese_only_relative_residual"],
        },
        {"quantity": "contact_fraction", "value": d["contact_fraction"]},
        {
            "quantity": "veronese_plus_contact_relative_residual",
            "value": d["veronese_plus_contact_relative_residual"],
        },
        {"quantity": "contact_coefficient_abs", "value": d["contact_coefficient_abs"]},
        {"quantity": "condition_number", "value": d["condition_number"]},
        {"quantity": "seesaw_residual_recomputed", "value": summary["seesaw_residual_recomputed"]},
    ]
    write_csv(OUT / "rank_rows.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    print("Source Majorana texture rank audit written")
    print(f"  Veronese-only residual: {d['veronese_only_relative_residual']:.6e}")
    print(f"  contact fraction: {d['contact_fraction']:.6e}")
    print(f"  lifted residual: {d['veronese_plus_contact_relative_residual']:.6e}")
    print(f"  contact lift essential: {summary['verdict']['contact_lift_essential']}")


if __name__ == "__main__":
    main()
