#!/usr/bin/env python3
"""Audit the publication closure card from its own machine-readable data.

No web lookup is used.  The point is to prevent a common failure mode in this
project: a later local improvement lives in one output directory while the
older global pipeline continues to report the old obstruction.  This audit
recomputes the local flavor observables directly from the publication card and
then separates local conditional closure from publication-final completeness.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402


CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
SUMMARY = ROOT / "output" / "publication_closure_card" / "summary.json"
OUT = ROOT / "output" / "publication_flavor_d5_reproducibility"

STRICT_CKM = 1.0e-3
LOOSE_MASS = 2.0e-1
SEESAW_RESIDUAL_MAX = 1.0e-10
TOL = 1.0e-10


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def recompute(y: dict[str, np.ndarray]) -> dict[str, Any]:
    v_ckm, ckm, _ = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    replay = fit.seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    return {
        "CKM_abs": fit.matrix_abs_json(v_ckm),
        "CKM_observables": ckm,
        "mass_ratios": {key: {"small": vals[0], "mid": vals[1]} for key, vals in masses.items()},
        "scores": {
            "ckm_score": float(ckm_score),
            "mass_score": float(mass_score),
        },
        "seesaw_replay": replay,
    }


def absdiff(a: float, b: float) -> float:
    return float(abs(float(a) - float(b)))


def build() -> dict[str, Any]:
    card = read_json(CARD)
    summary = read_json(SUMMARY)
    y = {name: cmat(raw) for name, raw in card["selected_row"]["Yukawa_fit"].items()}
    fresh = recompute(y)
    stored = card["recomputed_observables"]

    diffs = {
        "ckm_score": absdiff(fresh["scores"]["ckm_score"], stored["scores"]["ckm_score"]),
        "mass_score": absdiff(fresh["scores"]["mass_score"], stored["scores"]["mass_score"]),
        "seesaw_residual": absdiff(
            fresh["seesaw_replay"]["seesaw_matrix_residual"],
            stored["seesaw_replay"]["seesaw_matrix_residual"],
        ),
        "Vus": absdiff(fresh["CKM_observables"]["Vus"], stored["CKM_observables"]["Vus"]),
        "Vcb": absdiff(fresh["CKM_observables"]["Vcb"], stored["CKM_observables"]["Vcb"]),
        "Vub": absdiff(fresh["CKM_observables"]["Vub"], stored["CKM_observables"]["Vub"]),
    }
    reproducible = all(value < TOL for value in diffs.values())

    future_margin = float(card["crossed_projector"]["future_margin_1e35"])
    gates = {
        "card_recomputes_within_tolerance": reproducible,
        "strict_ckm": fresh["scores"]["ckm_score"] < STRICT_CKM,
        "loose_mass": fresh["scores"]["mass_score"] <= LOOSE_MASS,
        "seesaw_residual": fresh["seesaw_replay"]["seesaw_matrix_residual"] < SEESAW_RESIDUAL_MAX,
        "future_1e35_d5_margin": future_margin >= 1.0,
        "post_spin10_source_symmetry": bool(card["gates"]["post_spin10_source_symmetry"]),
        "dterm_unitary_lock": bool(card["gates"]["dterm_unitary_lock"]),
        "hidden_radial_lock": bool(card["gates"]["hidden_radial_lock"]),
        "endpoint_vectorlike_safe": bool(card["gates"]["endpoint_vectorlike_safe"]),
    }
    local_pass = all(gates.values())
    missing_publication_items = [
        "Replace internal flavor targets by a final cited input table.",
        "Regenerate all channel-specific d=5 widths from this exact card.",
        "Emit final paper tables for K+nu, e+pi0, mu+pi0, and K0mu+.",
        "Lock TeX tables to the generated JSON/CSV artifact hashes.",
    ]

    return {
        "note": "No web lookup used. Reproducibility audit for the publication closure card.",
        "input_card": str(CARD.relative_to(ROOT)),
        "fresh_recompute": fresh,
        "stored_summary_verdict": summary["verdict"],
        "absolute_differences": diffs,
        "gates": gates,
        "missing_publication_items": missing_publication_items,
        "verdict": {
            "local_source_consistent_candidate_reproducible": reproducible,
            "local_source_consistent_candidate_passes": local_pass,
            "publication_level_complete": False,
            "interpretation": (
                "The selected source-consistent crossed-120 card recomputes from its own "
                "Yukawa matrices and passes the local strict CKM, mass, seesaw, source/link, "
                "and 1e35 yr d=5 margin gates.  It remains publication-incomplete until the "
                "full channel-specific proton tables and final input manifest are regenerated "
                "from this exact card."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = [
        "# Publication flavor+d5 reproducibility audit",
        "",
        "No web lookup was used.",
        "",
        "## Gates",
        "",
        "| gate | pass |",
        "|---|---:|",
    ]
    for key, val in payload["gates"].items():
        lines.append(f"| `{key}` | `{val}` |")
    lines += [
        "",
        "## Recompute diffs",
        "",
        "| quantity | absolute diff |",
        "|---|---:|",
    ]
    for key, val in payload["absolute_differences"].items():
        lines.append(f"| `{key}` | {val:.6e} |")
    lines += [
        "",
        "## Missing publication items",
        "",
    ]
    for item in payload["missing_publication_items"]:
        lines.append(f"- {item}")
    lines += [
        "",
        "## Verdict",
        "",
        f"Local candidate reproducible: `{payload['verdict']['local_source_consistent_candidate_reproducible']}`.",
        f"Local candidate passes: `{payload['verdict']['local_source_consistent_candidate_passes']}`.",
        f"Publication complete: `{payload['verdict']['publication_level_complete']}`.",
        "",
        payload["verdict"]["interpretation"],
    ]
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "gate_status.csv",
        [{"gate": key, "passes": val} for key, val in payload["gates"].items()],
        ["gate", "passes"],
    )
    write_csv(
        OUT / "recompute_differences.csv",
        [{"quantity": key, "absolute_diff": val} for key, val in payload["absolute_differences"].items()],
        ["quantity", "absolute_diff"],
    )
    write_report(payload)
    v = payload["verdict"]
    fresh = payload["fresh_recompute"]
    print("Publication flavor+d5 reproducibility audit")
    print(f"  local reproducible: {v['local_source_consistent_candidate_reproducible']}")
    print(f"  local passes: {v['local_source_consistent_candidate_passes']}")
    print(f"  publication complete: {v['publication_level_complete']}")
    print(f"  CKM score: {fresh['scores']['ckm_score']:.6e}")
    print(f"  mass score: {fresh['scores']['mass_score']:.6e}")
    print(f"  seesaw residual: {fresh['seesaw_replay']['seesaw_matrix_residual']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
