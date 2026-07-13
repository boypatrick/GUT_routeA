#!/usr/bin/env python3
"""Build a no-web provenance ledger for the local flavor target layer.

The current branch has a reproducible local flavor/seesaw card, but it is not
a publication-final flavor fit.  This script separates those two statements:

* local matrix reproducibility and target scoring are audited here;
* external CKM/PMNS/mass target replacement remains explicitly open.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_publication_flavor_d5_reproducibility as repro  # noqa: E402
import scan_clebsch_flavor_fit as fit  # noqa: E402


OUT = ROOT / "output" / "no_web_flavor_target_provenance"
CLOSURE_CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
FLAVOR_REPRO = ROOT / "output" / "publication_flavor_d5_reproducibility" / "summary.json"
NO_WEB_LEDGER = ROOT / "output" / "no_web_input_convention_ledger" / "summary.json"


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


def log10_score(observed: float, target: float) -> float:
    return math.log10(observed / target) ** 2


def target_rows(flavor: dict[str, Any]) -> list[dict[str, Any]]:
    obs = flavor["fresh_recompute"]
    rows: list[dict[str, Any]] = []

    for key in ["Vus", "Vcb", "Vub", "J"]:
        observed = float(obs["CKM_observables"][key])
        target = float(fit.TARGETS[key])
        rows.append(
            {
                "block": "CKM",
                "observable": key,
                "target": target,
                "observed": observed,
                "ratio_observed_over_target": observed / target,
                "log10_score_contribution": log10_score(observed, target),
                "included_in_strict_ckm_gate": key in {"Vus", "Vcb", "Vub"},
                "external_refresh_needed": True,
            }
        )

    for sector, vals in obs["mass_ratios"].items():
        for label in ["small", "mid"]:
            key = f"{sector}_{label}"
            observed = float(vals[label])
            target = float(fit.TARGETS[key])
            rows.append(
                {
                    "block": "mass_ratio",
                    "observable": key,
                    "target": target,
                    "observed": observed,
                    "ratio_observed_over_target": observed / target,
                    "log10_score_contribution": log10_score(observed, target),
                    "included_in_strict_ckm_gate": False,
                    "external_refresh_needed": True,
                }
            )

    return rows


def seesaw_rows(no_web: dict[str, Any], flavor: dict[str, Any]) -> list[dict[str, Any]]:
    replay = flavor["fresh_recompute"]["seesaw_replay"]
    rows = [
        {
            "block": "seesaw_replay",
            "observable": "seesaw_matrix_residual",
            "target": repro.SEESAW_RESIDUAL_MAX,
            "observed": replay["seesaw_matrix_residual"],
            "ratio_observed_over_target": replay["seesaw_matrix_residual"] / repro.SEESAW_RESIDUAL_MAX,
            "passes": replay["seesaw_matrix_residual"] < repro.SEESAW_RESIDUAL_MAX,
            "external_refresh_needed": False,
        },
        {
            "block": "seesaw_replay",
            "observable": "theta_norm",
            "target": "",
            "observed": replay["theta_norm"],
            "ratio_observed_over_target": "",
            "passes": True,
            "external_refresh_needed": False,
        },
        {
            "block": "seesaw_replay",
            "observable": "MR_condition_number",
            "target": "",
            "observed": replay["MR_condition_number"],
            "ratio_observed_over_target": "",
            "passes": True,
            "external_refresh_needed": False,
        },
    ]
    for idx, mass in enumerate(replay["heavy_neutrino_masses_GeV"], start=1):
        rows.append(
            {
                "block": "seesaw_replay",
                "observable": f"M_R{idx}_GeV",
                "target": "",
                "observed": mass,
                "ratio_observed_over_target": "",
                "passes": True,
                "external_refresh_needed": False,
            }
        )
    for key, value in no_web["seesaw_benchmark"].items():
        rows.append(
            {
                "block": "fixed_light_neutrino_benchmark",
                "observable": key,
                "target": value,
                "observed": value,
                "ratio_observed_over_target": 1.0,
                "passes": True,
                "external_refresh_needed": True,
            }
        )
    return rows


def matrix_manifest(closure: dict[str, Any]) -> list[dict[str, Any]]:
    matrices = closure["selected_row"]["Yukawa_fit"]
    rows: list[dict[str, Any]] = []
    for name in ["up", "down", "charged_lepton", "neutrino_dirac"]:
        raw = matrices[name]
        rows.append(
            {
                "matrix": name,
                "shape": "3x3 complex",
                "entry_count": 9,
                "digest": stable_digest(raw),
                "source_card": str(CLOSURE_CARD.relative_to(ROOT)),
                "source_card_sha256": sha256(CLOSURE_CARD),
            }
        )
    rows.append(
        {
            "matrix": "all_yukawa_matrices",
            "shape": "4 x 3x3 complex",
            "entry_count": 36,
            "digest": stable_digest(matrices),
            "source_card": str(CLOSURE_CARD.relative_to(ROOT)),
            "source_card_sha256": sha256(CLOSURE_CARD),
        }
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build() -> dict[str, Any]:
    closure = read_json(CLOSURE_CARD)
    flavor = read_json(FLAVOR_REPRO)
    no_web = read_json(NO_WEB_LEDGER)
    targets_match = no_web["flavor_targets"] == fit.TARGETS
    rows = target_rows(flavor)
    seesaw = seesaw_rows(no_web, flavor)
    matrices = matrix_manifest(closure)
    ckm_score_from_rows = sum(row["log10_score_contribution"] for row in rows if row["observable"] in {"Vus", "Vcb", "Vub"})
    mass_score_from_rows = sum(row["log10_score_contribution"] for row in rows if row["block"] == "mass_ratio")
    gates = flavor["gates"]
    summary = {
        "note": "No web lookup used. Flavor target provenance separates local reproducibility from external target refresh.",
        "input_manifest": [
            {
                "label": "publication_closure_card",
                "path": str(CLOSURE_CARD.relative_to(ROOT)),
                "sha256": sha256(CLOSURE_CARD),
            },
            {
                "label": "publication_flavor_reproducibility",
                "path": str(FLAVOR_REPRO.relative_to(ROOT)),
                "sha256": sha256(FLAVOR_REPRO),
            },
            {
                "label": "no_web_input_convention_ledger",
                "path": str(NO_WEB_LEDGER.relative_to(ROOT)),
                "sha256": sha256(NO_WEB_LEDGER),
            },
        ],
        "target_rows": rows,
        "seesaw_rows": seesaw,
        "matrix_manifest": matrices,
        "target_row_count": len(rows),
        "seesaw_row_count": len(seesaw),
        "matrix_row_count": len(matrices),
        "scores": {
            "ckm_score_from_rows": ckm_score_from_rows,
            "ckm_score_stored": flavor["fresh_recompute"]["scores"]["ckm_score"],
            "ckm_score_absdiff": abs(ckm_score_from_rows - flavor["fresh_recompute"]["scores"]["ckm_score"]),
            "mass_score_from_rows": mass_score_from_rows,
            "mass_score_stored": flavor["fresh_recompute"]["scores"]["mass_score"],
            "mass_score_absdiff": abs(mass_score_from_rows - flavor["fresh_recompute"]["scores"]["mass_score"]),
        },
        "gates": gates,
        "verdict": {
            "local_targets_match_no_web_ledger": targets_match,
            "local_card_recomputes_within_tolerance": flavor["verdict"]["local_source_consistent_candidate_reproducible"],
            "local_flavor_gates_pass": bool(gates["strict_ckm"] and gates["loose_mass"] and gates["seesaw_residual"]),
            "score_rows_reproduce_stored_scores": (
                abs(ckm_score_from_rows - flavor["fresh_recompute"]["scores"]["ckm_score"]) < repro.TOL
                and abs(mass_score_from_rows - flavor["fresh_recompute"]["scores"]["mass_score"]) < repro.TOL
            ),
            "external_target_refresh_done": False,
            "pmns_completion_done": False,
            "interpretation": (
                "The local source-consistent flavor card is reproducible from its "
                "embedded Yukawa matrices and the no-web target table.  The open "
                "publication item is narrower: replace local targets by cited "
                "external CKM/PMNS/mass inputs and complete a true PMNS/global "
                "flavor fit."
            ),
        },
    }
    return summary


def report(summary: dict[str, Any]) -> str:
    v = summary["verdict"]
    s = summary["scores"]
    lines = [
        "# No-web flavor target provenance",
        "",
        "No web lookup was used.  This ledger separates local flavor reproducibility from the external target-refresh task.",
        "",
        f"- target rows: `{summary['target_row_count']}`",
        f"- seesaw rows: `{summary['seesaw_row_count']}`",
        f"- matrix rows: `{summary['matrix_row_count']}`",
        f"- CKM score diff: `{s['ckm_score_absdiff']:.3e}`",
        f"- mass score diff: `{s['mass_score_absdiff']:.3e}`",
        f"- local targets match no-web ledger: `{v['local_targets_match_no_web_ledger']}`",
        f"- local card recomputes: `{v['local_card_recomputes_within_tolerance']}`",
        f"- local flavor gates pass: `{v['local_flavor_gates_pass']}`",
        f"- external target refresh done: `{v['external_target_refresh_done']}`",
        f"- PMNS completion done: `{v['pmns_completion_done']}`",
        "",
        "## Interpretation",
        "",
        v["interpretation"],
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    write_csv(OUT / "target_rows.csv", summary["target_rows"])
    write_csv(OUT / "seesaw_rows.csv", summary["seesaw_rows"])
    write_csv(OUT / "matrix_manifest.csv", summary["matrix_manifest"])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    print("No-web flavor target provenance written")
    print(f"  target rows: {summary['target_row_count']}")
    print(f"  score rows reproduce stored scores: {summary['verdict']['score_rows_reproduce_stored_scores']}")
    print(f"  local flavor gates pass: {summary['verdict']['local_flavor_gates_pass']}")


if __name__ == "__main__":
    main()
