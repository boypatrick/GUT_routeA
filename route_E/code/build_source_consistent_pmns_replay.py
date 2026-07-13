#!/usr/bin/env python3
"""Replay PMNS compatibility for the source-consistent publication card.

No web lookup is used.  This script answers a narrow question raised by the
heartbeat roadmap: after the source-consistent crossed-120 CKM repair, do the
published Y_e and Y_nu matrices still support the same PMNS convention?

The answer is deliberately phrased as compatibility, not prediction.  We use
the same inverse type-I seesaw construction as the local flavor scans: for the
given Y_e and Y_nu, construct M_R so that the fixed local PMNS benchmark is
reproduced.  Passing this audit means the source-consistent card is internally
PMNS-compatible; it does not mean the PMNS angles have been derived without
choosing M_R.
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

import scan_clebsch_flavor_fit as fit  # noqa: E402
import verify_seesaw_item3 as seesaw  # noqa: E402


OUT = ROOT / "output" / "source_consistent_pmns_replay"
CLOSURE_CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
NO_WEB_PMNS_BENCH = ROOT / "output" / "no_web_pmns_benchmark_replay" / "summary.json"


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


def mixing_angles(u: np.ndarray) -> dict[str, float]:
    uabs2 = np.abs(u) ** 2
    s13 = float(uabs2[0, 2])
    denom = max(1.0 - s13, 1e-30)
    return {
        "sin2_theta12": float(uabs2[0, 1] / denom),
        "sin2_theta13": s13,
        "sin2_theta23": float(uabs2[1, 2] / denom),
    }


def jarlskog(u: np.ndarray) -> float:
    return float(np.imag(u[0, 0] * u[1, 1] * np.conjugate(u[0, 1]) * np.conjugate(u[1, 0])))


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


def build() -> dict[str, Any]:
    card = read_json(CLOSURE_CARD)
    bench_replay = read_json(NO_WEB_PMNS_BENCH)
    y = {name: cmat(raw) for name, raw in card["selected_row"]["Yukawa_fit"].items()}
    y_e = y["charged_lepton"]
    y_nu = y["neutrino_dirac"]
    u_e, e_singulars = seesaw.left_rotation(y_e)
    b = fixed_benchmark()
    m_light_diag = np.array(
        [
            b.m1_eV,
            math.sqrt(b.m1_eV**2 + b.dm21_eV2),
            math.sqrt(b.m1_eV**2 + b.dm31_eV2),
        ],
        dtype=float,
    )
    u_pmns_target = seesaw.standard_pmns(
        b.sin2_theta12,
        b.sin2_theta13,
        b.sin2_theta23,
        b.delta_cp_rad,
        b.alpha21_rad,
        b.alpha31_rad,
    )
    u_nu_target = u_e @ u_pmns_target
    u_pmns_replayed = u_e.conjugate().T @ u_nu_target
    m_light_target = u_nu_target.conjugate() @ np.diag(m_light_diag) @ u_nu_target.conjugate().T

    mD_largest_GeV = 100.0
    mD_eV = y_nu * (mD_largest_GeV * 1.0e9)
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light_target) @ mD_eV)
    MR_GeV = MR_eV / 1.0e9
    m_light_reco = -(mD_eV @ np.linalg.inv(MR_eV) @ mD_eV.T)
    seesaw_residual = float(np.linalg.norm(m_light_reco - m_light_target) / np.linalg.norm(m_light_target))
    theta_norm = float(np.linalg.norm(mD_eV @ np.linalg.inv(MR_eV), ord=2))
    heavy = np.sort(np.linalg.svd(MR_GeV, compute_uv=False))

    diag_target = u_nu_target.T @ m_light_target @ u_nu_target
    diag_reco = u_nu_target.T @ m_light_reco @ u_nu_target
    offdiag_target = diag_target - np.diag(np.diag(diag_target))
    offdiag_reco = diag_reco - np.diag(np.diag(diag_reco))
    mass_reco = np.real(np.diag(diag_reco))

    angles_target = {
        "sin2_theta12": b.sin2_theta12,
        "sin2_theta13": b.sin2_theta13,
        "sin2_theta23": b.sin2_theta23,
    }
    angles_replayed = mixing_angles(u_pmns_replayed)
    rows: list[dict[str, Any]] = []
    for key in ["sin2_theta12", "sin2_theta13", "sin2_theta23"]:
        rows.append(
            {
                "block": "PMNS_angles",
                "observable": key,
                "target": angles_target[key],
                "observed": angles_replayed[key],
                "difference": angles_replayed[key] - angles_target[key],
                "abs_difference": abs(angles_replayed[key] - angles_target[key]),
                "passes_1e_minus_10": abs(angles_replayed[key] - angles_target[key]) < 1.0e-10,
                "interpretation": "PMNS target replayed by inverse-seesaw construction",
            }
        )
    rows.append(
        {
            "block": "PMNS_CP",
            "observable": "J_PMNS",
            "target": jarlskog(u_pmns_target),
            "observed": jarlskog(u_pmns_replayed),
            "difference": jarlskog(u_pmns_replayed) - jarlskog(u_pmns_target),
            "abs_difference": abs(jarlskog(u_pmns_replayed) - jarlskog(u_pmns_target)),
            "passes_1e_minus_10": abs(jarlskog(u_pmns_replayed) - jarlskog(u_pmns_target)) < 1.0e-10,
            "interpretation": "CP invariant replayed by fixed PMNS target",
        }
    )
    for idx, (target, observed) in enumerate(zip(m_light_diag, mass_reco), start=1):
        rows.append(
            {
                "block": "light_masses",
                "observable": f"m{idx}_eV",
                "target": float(target),
                "observed": float(observed),
                "difference": float(observed - target),
                "abs_difference": float(abs(observed - target)),
                "passes_1e_minus_10": float(abs(observed - target)) < 1.0e-10,
                "interpretation": "diagonalized in the imposed U_nu target basis",
            }
        )
    dm_target = {
        "dm21_eV2": m_light_diag[1] ** 2 - m_light_diag[0] ** 2,
        "dm31_eV2": m_light_diag[2] ** 2 - m_light_diag[0] ** 2,
    }
    dm_reco = {
        "dm21_eV2": mass_reco[1] ** 2 - mass_reco[0] ** 2,
        "dm31_eV2": mass_reco[2] ** 2 - mass_reco[0] ** 2,
    }
    for key in ["dm21_eV2", "dm31_eV2"]:
        rows.append(
            {
                "block": "mass_splittings",
                "observable": key,
                "target": float(dm_target[key]),
                "observed": float(dm_reco[key]),
                "difference": float(dm_reco[key] - dm_target[key]),
                "abs_difference": float(abs(dm_reco[key] - dm_target[key])),
                "passes_1e_minus_10": float(abs(dm_reco[key] - dm_target[key])) < 1.0e-10,
                "interpretation": "mass splitting replayed after inverse-seesaw reconstruction",
            }
        )

    residuals = {
        "max_pmns_angle_absdiff": max(
            row["abs_difference"] for row in rows if row["block"] == "PMNS_angles"
        ),
        "max_light_mass_absdiff_eV": max(
            row["abs_difference"] for row in rows if row["block"] == "light_masses"
        ),
        "max_mass_splitting_absdiff_eV2": max(
            row["abs_difference"] for row in rows if row["block"] == "mass_splittings"
        ),
        "seesaw_matrix_residual": seesaw_residual,
        "target_basis_offdiag_norm": float(np.linalg.norm(offdiag_target)),
        "reco_basis_offdiag_norm": float(np.linalg.norm(offdiag_reco)),
        "theta_norm": theta_norm,
    }
    pass_compat = (
        residuals["max_pmns_angle_absdiff"] < 1.0e-10
        and residuals["max_light_mass_absdiff_eV"] < 1.0e-10
        and residuals["max_mass_splitting_absdiff_eV2"] < 1.0e-10
        and residuals["seesaw_matrix_residual"] < 1.0e-10
    )
    matrix_manifest = [
        {
            "object": "publication_Y_e",
            "digest": stable_digest(card["selected_row"]["Yukawa_fit"]["charged_lepton"]),
        },
        {
            "object": "publication_Y_nu",
            "digest": stable_digest(card["selected_row"]["Yukawa_fit"]["neutrino_dirac"]),
        },
        {
            "object": "source_consistent_PMNS_pair",
            "digest": stable_digest(
                {
                    "U_PMNS_target": matrix_json(u_pmns_target),
                    "U_PMNS_replayed": matrix_json(u_pmns_replayed),
                }
            ),
        },
        {
            "object": "source_consistent_MR_GeV",
            "digest": stable_digest(matrix_json(MR_GeV)),
        },
    ]
    return {
        "note": "No web lookup used. Source-consistent PMNS compatibility replay.",
        "input_manifest": [
            {
                "label": "publication_closure_card",
                "path": str(CLOSURE_CARD.relative_to(ROOT)),
                "sha256": sha256(CLOSURE_CARD),
            },
            {
                "label": "no_web_pmns_benchmark_replay",
                "path": str(NO_WEB_PMNS_BENCH.relative_to(ROOT)),
                "sha256": sha256(NO_WEB_PMNS_BENCH),
            },
        ],
        "observable_row_count": len(rows),
        "matrix_row_count": len(matrix_manifest),
        "observable_rows": rows,
        "matrix_manifest": matrix_manifest,
        "residuals": residuals,
        "heavy_neutrino_masses_GeV": [float(x) for x in heavy],
        "MR_condition_number": float(heavy[-1] / max(heavy[0], 1.0e-30)),
        "Y_e_singular_values_normalized": [float(x) for x in e_singulars],
        "benchmark_pmns_pair_digest": bench_replay["pmns_pair_digest"],
        "source_pmns_pair_digest": next(
            row["digest"] for row in matrix_manifest if row["object"] == "source_consistent_PMNS_pair"
        ),
        "verdict": {
            "source_consistent_pmns_compatible": pass_compat,
            "pmns_predictive_fit_done": False,
            "external_pmns_target_refresh_done": False,
            "pmns_completion_done": False,
            "interpretation": (
                "The source-consistent crossed-120 publication card is PMNS-compatible: "
                "using its Y_e and Y_nu, the inverse type-I seesaw reconstructs an M_R "
                "for the fixed local PMNS benchmark with sub-1e-10 residuals.  This is "
                "not yet a predictive PMNS fit, because the right-handed Majorana matrix "
                "is reconstructed from the chosen PMNS target rather than derived from a "
                "separate Majorana texture or external global-fit manifest."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def report(summary: dict[str, Any]) -> str:
    v = summary["verdict"]
    r = summary["residuals"]
    return "\n".join(
        [
            "# Source-consistent PMNS replay",
            "",
            "No web lookup was used.  This replay checks PMNS compatibility for the source-consistent publication card.",
            "",
            f"- observable rows: `{summary['observable_row_count']}`",
            f"- matrix rows: `{summary['matrix_row_count']}`",
            f"- source PMNS pair digest: `{summary['source_pmns_pair_digest']}`",
            f"- max PMNS angle residual: `{r['max_pmns_angle_absdiff']:.3e}`",
            f"- max light-mass residual: `{r['max_light_mass_absdiff_eV']:.3e}` eV",
            f"- max mass-splitting residual: `{r['max_mass_splitting_absdiff_eV2']:.3e}` eV^2",
            f"- seesaw matrix residual: `{r['seesaw_matrix_residual']:.3e}`",
            f"- source-consistent PMNS compatible: `{v['source_consistent_pmns_compatible']}`",
            f"- PMNS predictive fit done: `{v['pmns_predictive_fit_done']}`",
            "",
            "## Interpretation",
            "",
            v["interpretation"],
            "",
        ]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    write_csv(OUT / "observable_rows.csv", summary["observable_rows"])
    write_csv(OUT / "matrix_manifest.csv", summary["matrix_manifest"])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    print("Source-consistent PMNS replay written")
    print(f"  observable rows: {summary['observable_row_count']}")
    print(f"  source-consistent PMNS compatible: {summary['verdict']['source_consistent_pmns_compatible']}")
    print(f"  PMNS predictive fit done: {summary['verdict']['pmns_predictive_fit_done']}")


if __name__ == "__main__":
    main()
