#!/usr/bin/env python3
"""Build a no-web PMNS benchmark replay ledger.

This audit is intentionally narrow.  It verifies that the existing exact
CP1/O(2) seesaw benchmark reconstructs the local PMNS angles and light-neutrino
mass splittings from the stored matrices, while keeping the publication-level
statement open until the source-consistent flavor card is integrated with final
cited external PMNS targets.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "no_web_pmns_benchmark_replay"
BENCHMARK_CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
FLAVOR_TARGETS = ROOT / "output" / "no_web_flavor_target_provenance" / "summary.json"


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


def complex_matrix(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def mixing_angles(u: np.ndarray) -> dict[str, float]:
    uabs2 = np.abs(u) ** 2
    s13 = float(uabs2[0, 2])
    denom13 = max(1.0 - s13, 1e-30)
    return {
        "sin2_theta12": float(uabs2[0, 1] / denom13),
        "sin2_theta13": s13,
        "sin2_theta23": float(uabs2[1, 2] / denom13),
    }


def jarlskog(u: np.ndarray) -> float:
    return float(np.imag(u[0, 0] * u[1, 1] * np.conjugate(u[0, 1]) * np.conjugate(u[1, 0])))


def jarlskog_denominator(angles: dict[str, float]) -> float:
    s12 = math.sqrt(angles["sin2_theta12"])
    s13 = math.sqrt(angles["sin2_theta13"])
    s23 = math.sqrt(angles["sin2_theta23"])
    c12 = math.sqrt(max(1.0 - angles["sin2_theta12"], 0.0))
    c13 = math.sqrt(max(1.0 - angles["sin2_theta13"], 0.0))
    c23 = math.sqrt(max(1.0 - angles["sin2_theta23"], 0.0))
    return c12 * c23 * c13 * c13 * s12 * s23 * s13


def observable_rows(seesaw: dict[str, Any]) -> list[dict[str, Any]]:
    benchmark = seesaw["benchmark"]
    u_target = complex_matrix(seesaw["U_PMNS_target"])
    u_reco = complex_matrix(seesaw["U_PMNS_reconstructed"])
    angles_target = {
        key: float(benchmark[key])
        for key in ["sin2_theta12", "sin2_theta13", "sin2_theta23"]
    }
    angles_reco = mixing_angles(u_reco)
    rows: list[dict[str, Any]] = []
    for key in ["sin2_theta12", "sin2_theta13", "sin2_theta23"]:
        target = angles_target[key]
        observed = angles_reco[key]
        rows.append(
            {
                "block": "PMNS_angles",
                "observable": key,
                "target": target,
                "observed": observed,
                "difference": observed - target,
                "abs_difference": abs(observed - target),
                "passes_1e_minus_10": abs(observed - target) < 1e-10,
                "external_refresh_needed": True,
            }
        )

    j_target = jarlskog(u_target)
    j_reco = jarlskog(u_reco)
    den_target = jarlskog_denominator(angles_target)
    den_reco = jarlskog_denominator(angles_reco)
    rows.extend(
        [
            {
                "block": "PMNS_CP",
                "observable": "J_PMNS",
                "target": j_target,
                "observed": j_reco,
                "difference": j_reco - j_target,
                "abs_difference": abs(j_reco - j_target),
                "passes_1e_minus_10": abs(j_reco - j_target) < 1e-10,
                "external_refresh_needed": True,
            },
            {
                "block": "PMNS_CP",
                "observable": "sin_delta_from_J",
                "target": j_target / den_target,
                "observed": j_reco / den_reco,
                "difference": j_reco / den_reco - j_target / den_target,
                "abs_difference": abs(j_reco / den_reco - j_target / den_target),
                "passes_1e_minus_10": abs(j_reco / den_reco - j_target / den_target) < 1e-10,
                "external_refresh_needed": True,
            },
        ]
    )

    target_masses = [float(x) for x in seesaw["target_light_masses_eV"]]
    reco_masses = [float(x) for x in seesaw["reconstructed_light_masses_eV"]]
    for idx, (target, observed) in enumerate(zip(target_masses, reco_masses), start=1):
        rows.append(
            {
                "block": "light_masses",
                "observable": f"m{idx}_eV",
                "target": target,
                "observed": observed,
                "difference": observed - target,
                "abs_difference": abs(observed - target),
                "passes_1e_minus_10": abs(observed - target) < 1e-10,
                "external_refresh_needed": True,
            }
        )

    dm_target = {
        "dm21_eV2": target_masses[1] ** 2 - target_masses[0] ** 2,
        "dm31_eV2": target_masses[2] ** 2 - target_masses[0] ** 2,
    }
    dm_reco = {
        "dm21_eV2": reco_masses[1] ** 2 - reco_masses[0] ** 2,
        "dm31_eV2": reco_masses[2] ** 2 - reco_masses[0] ** 2,
    }
    for key in ["dm21_eV2", "dm31_eV2"]:
        rows.append(
            {
                "block": "mass_splittings",
                "observable": key,
                "target": dm_target[key],
                "observed": dm_reco[key],
                "difference": dm_reco[key] - dm_target[key],
                "abs_difference": abs(dm_reco[key] - dm_target[key]),
                "passes_1e_minus_10": abs(dm_reco[key] - dm_target[key]) < 1e-10,
                "external_refresh_needed": True,
            }
        )

    for key in ["delta_cp_rad", "alpha21_rad", "alpha31_rad"]:
        rows.append(
            {
                "block": "phase_inputs",
                "observable": key,
                "target": float(benchmark[key]),
                "observed": float(benchmark[key]),
                "difference": 0.0,
                "abs_difference": 0.0,
                "passes_1e_minus_10": True,
                "external_refresh_needed": True,
            }
        )
    return rows


def matrix_rows(seesaw: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for key in ["U_PMNS_target", "U_PMNS_reconstructed", "m_light_target_eV", "m_light_reconstructed_eV"]:
        rows.append(
            {
                "object": key,
                "digest": stable_digest(seesaw[key]),
                "source_card": str(BENCHMARK_CARD.relative_to(ROOT)),
                "source_card_sha256": sha256(BENCHMARK_CARD),
            }
        )
    rows.append(
        {
            "object": "pmns_pair",
            "digest": stable_digest(
                {
                    "U_PMNS_target": seesaw["U_PMNS_target"],
                    "U_PMNS_reconstructed": seesaw["U_PMNS_reconstructed"],
                }
            ),
            "source_card": str(BENCHMARK_CARD.relative_to(ROOT)),
            "source_card_sha256": sha256(BENCHMARK_CARD),
        }
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def report(summary: dict[str, Any]) -> str:
    v = summary["verdict"]
    return "\n".join(
        [
            "# No-web PMNS benchmark replay",
            "",
            "No web lookup was used.  This audit verifies the local exact CP1/O(2) PMNS benchmark replay and keeps publication-level PMNS completion separate.",
            "",
            f"- observable rows: `{summary['observable_row_count']}`",
            f"- matrix rows: `{summary['matrix_row_count']}`",
            f"- PMNS pair digest: `{summary['pmns_pair_digest']}`",
            f"- max PMNS angle residual: `{summary['residuals']['max_pmns_angle_absdiff']:.3e}`",
            f"- max light-mass residual: `{summary['residuals']['max_light_mass_absdiff_eV']:.3e}` eV",
            f"- max mass-splitting residual: `{summary['residuals']['max_mass_splitting_absdiff_eV2']:.3e}` eV^2",
            f"- local PMNS benchmark replay passes: `{v['local_pmns_benchmark_replay_passes']}`",
            f"- source-consistent publication card integrated: `{v['source_consistent_publication_card_integrated']}`",
            f"- external PMNS target refresh done: `{v['external_pmns_target_refresh_done']}`",
            "",
            "## Interpretation",
            "",
            v["interpretation"],
            "",
        ]
    )


def build() -> dict[str, Any]:
    card = read_json(BENCHMARK_CARD)
    flavor_targets = read_json(FLAVOR_TARGETS)
    seesaw = card["seesaw"]
    rows = observable_rows(seesaw)
    matrices = matrix_rows(seesaw)
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
        "seesaw_matrix_residual": seesaw["seesaw_matrix_residual"],
        "light_mass_residual": seesaw["light_mass_residual"],
    }
    local_pass = (
        residuals["max_pmns_angle_absdiff"] < 1e-10
        and residuals["max_light_mass_absdiff_eV"] < 1e-10
        and residuals["max_mass_splitting_absdiff_eV2"] < 1e-10
        and seesaw["checks"]["pmns_angle_max_residual_lt_1e-10"]
        and seesaw["checks"]["light_mass_residual_lt_1e-10"]
    )
    summary = {
        "note": "No web lookup used. PMNS benchmark replay is separated from publication-level PMNS completion.",
        "input_manifest": [
            {
                "label": "flavor_benchmark_card",
                "path": str(BENCHMARK_CARD.relative_to(ROOT)),
                "sha256": sha256(BENCHMARK_CARD),
            },
            {
                "label": "no_web_flavor_target_provenance",
                "path": str(FLAVOR_TARGETS.relative_to(ROOT)),
                "sha256": sha256(FLAVOR_TARGETS),
            },
        ],
        "observable_row_count": len(rows),
        "matrix_row_count": len(matrices),
        "pmns_pair_digest": next(row["digest"] for row in matrices if row["object"] == "pmns_pair"),
        "observable_rows": rows,
        "matrix_rows": matrices,
        "residuals": residuals,
        "linked_flavor_target_verdict": flavor_targets["verdict"],
        "verdict": {
            "local_pmns_benchmark_replay_passes": local_pass,
            "source_consistent_publication_card_integrated": False,
            "external_pmns_target_refresh_done": False,
            "pmns_completion_done": False,
            "interpretation": (
                "The exact CP1/O(2) seesaw benchmark reconstructs the local PMNS angles, "
                "light masses, and mass splittings at the 1e-10 gate.  This is a useful "
                "no-web benchmark replay, but it is not yet the final source-consistent "
                "publication flavor fit: the crossed-120 source-consistent card still "
                "needs the same PMNS observable table after final cited target refresh."
            ),
        },
    }
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    write_csv(OUT / "observable_rows.csv", summary["observable_rows"])
    write_csv(OUT / "matrix_rows.csv", summary["matrix_rows"])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    print("No-web PMNS benchmark replay written")
    print(f"  observable rows: {summary['observable_row_count']}")
    print(f"  matrix rows: {summary['matrix_row_count']}")
    print(f"  local PMNS benchmark replay passes: {summary['verdict']['local_pmns_benchmark_replay_passes']}")


if __name__ == "__main__":
    main()
