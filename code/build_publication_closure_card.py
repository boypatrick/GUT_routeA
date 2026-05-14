#!/usr/bin/env python3
"""Build a single-source publication closure candidate card.

No web lookup is used.  This is not another fit.  It promotes the best
source-consistent crossed-120 local closure into one machine-readable card and
records every upstream artifact needed to reproduce the conditional claim.

The card is deliberately honest about scope:

* local flavor+d=5 closure is checked from the crossed-120 source row;
* the crossed projector is post-Spin(10) / source-sector conditional;
* publication-level completeness still requires a final channel-by-channel
  d=5 proton-decay table and a single refreshed phenomenology input table.
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


OUT = ROOT / "output" / "publication_closure_card"

SOURCE_CKM = ROOT / "output" / "source_consistent_ckm_crossed120" / "summary.json"
CROSSED_PROJECTOR = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"
CROSSED_SOURCE_SYMM = ROOT / "output" / "crossed120_action_level_source_symmetry" / "summary.json"
UNITARY_DTERM = ROOT / "output" / "unitary_link_dterm_quotient" / "summary.json"
HIDDEN_RADIAL = ROOT / "output" / "hidden_radial_lock_sector" / "summary.json"
ENDPOINT_VECTORLIKE = ROOT / "output" / "endpoint_vectorlike_completion" / "summary.json"
LOCKED_LINK_CARD = ROOT / "output" / "locked_link_full_flavor_d5_card" / "locked_link_full_flavor_d5_card.json"
FULL_FLAVOR_OLD = ROOT / "output" / "full_flavor_d5_pipeline" / "summary.json"
THEOREM_LEDGER = ROOT / "output" / "conditional_theorem_ledger" / "summary.json"
TWO_KERNEL = ROOT / "output" / "two_kernel_flavor_then_d5" / "summary.json"

STRICT_CKM = 1.0e-3
LOOSE_MASS = 2.0e-1
SEESAW_RESIDUAL_MAX = 1.0e-10


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path) if path.exists() else None,
            }
        )
    return rows


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def recompute_flavor(y: dict[str, np.ndarray]) -> dict[str, Any]:
    v_ckm, ckm, _ = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    jarlskog_score = fit.log10_ratio(ckm["J"], fit.TARGETS["J"]) ** 2
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
            "jarlskog_score": float(jarlskog_score),
            "mass_score": float(mass_score),
        },
        "seesaw_replay": replay,
    }


def find_row(rows: list[dict[str, Any]], label: str) -> dict[str, Any]:
    return next(row for row in rows if row["label"] == label)


def finite_margin_kappa(payload: dict[str, Any], kappa: float = 30.0) -> float:
    row = next(item for item in payload["finite_leakage_replay"] if float(item["kappa"]) == float(kappa))
    return float(row["worst_margin_1e35_replayed"])


def diff(a: float, b: float) -> float:
    return float(abs(a - b))


def build() -> dict[str, Any]:
    source = read_json(SOURCE_CKM)
    crossed = read_json(CROSSED_PROJECTOR)
    source_symm = read_json(CROSSED_SOURCE_SYMM)
    dterm = read_json(UNITARY_DTERM)
    radial = read_json(HIDDEN_RADIAL)
    endpoint = read_json(ENDPOINT_VECTORLIKE)
    old_pipeline = read_json(FULL_FLAVOR_OLD)
    ledger = read_json(THEOREM_LEDGER)
    two_kernel = read_json(TWO_KERNEL)

    best = source["verdict"]["best_closure"]
    row = find_row(source["rows"], best["label"])
    y = {name: cmat(raw) for name, raw in row["Yukawa_fit"].items()}
    recomputed = recompute_flavor(y)
    margin = finite_margin_kappa(crossed, 30.0)

    source_checks = {
        "ckm_score_abs_diff": diff(recomputed["scores"]["ckm_score"], float(row["ckm_score"])),
        "mass_score_abs_diff": diff(recomputed["scores"]["mass_score"], float(row["mass_score"])),
        "seesaw_residual_abs_diff": diff(
            float(recomputed["seesaw_replay"]["seesaw_matrix_residual"]),
            float(row["seesaw_replay"]["seesaw_matrix_residual"]),
        ),
        "d5_margin_abs_diff": diff(margin, float(row["future_margin_1e35_crossed120_kappa30"])),
    }

    gates = {
        "strict_ckm": recomputed["scores"]["ckm_score"] < STRICT_CKM,
        "loose_mass": recomputed["scores"]["mass_score"] <= LOOSE_MASS,
        "seesaw_residual": recomputed["seesaw_replay"]["seesaw_matrix_residual"] < SEESAW_RESIDUAL_MAX,
        "future_1e35_d5_margin": margin >= 1.0,
        "post_spin10_source_symmetry": bool(
            source_symm["verdict"]["post_spin10_source_symmetry_branch_viable"]
        ),
        "dterm_unitary_lock": bool(dterm["verdict"]["Dterm_removes_nonunitary_holomorphic_moduli"]),
        "hidden_radial_lock": bool(radial["verdict"]["radial_Dterms_remove_nonunitary_modes"]),
        "endpoint_vectorlike_safe": bool(endpoint["verdict"]["vectorlike_completion_beta_safe_to_R200"]),
    }
    local_candidate_passes = all(gates.values())

    payload = {
        "note": "No web lookup used. Publication closure candidate card for the source-consistent crossed-120 branch.",
        "scope": {
            "claim": "conditional local flavor+d=5 closure with crossed 120 source projector",
            "not_claimed": [
                "field-only unbroken Spin(10) derivation of the crossed projector",
                "final literature-refresh input table",
                "complete channel-by-channel d=5 proton-decay paper table",
                "unconditional first-principles GUT theorem",
            ],
        },
        "formulae": {
            "CKM": "V_CKM = U_u^dagger U_d from left rotations of Y_u Y_u^dagger and Y_d Y_d^dagger",
            "mass_ratios": "sector ratios are singular values normalized to the largest singular value",
            "type_I_seesaw": "M_R = -m_D^T m_nu^{-1} m_D, replayed against the fixed PMNS benchmark",
            "crossed_d5_margin": "finite-lift crossed-120 kappa=30 worst Knu margin against 1e35 yr stress",
        },
        "input_manifest": manifest(
            [
                SOURCE_CKM,
                CROSSED_PROJECTOR,
                CROSSED_SOURCE_SYMM,
                UNITARY_DTERM,
                HIDDEN_RADIAL,
                ENDPOINT_VECTORLIKE,
                LOCKED_LINK_CARD,
                FULL_FLAVOR_OLD,
                THEOREM_LEDGER,
                TWO_KERNEL,
            ]
        ),
        "selected_row": {
            "label": row["label"],
            "mode": row["mode"],
            "mixing_coefficients": row["mixing_coefficients"],
            "doublet_deformation": row["doublet_deformation"],
            "Yukawa_fit": {name: matrix_json(mat) for name, mat in y.items()},
        },
        "recomputed_observables": recomputed,
        "source_row_observables": {
            "CKM_observables": row["CKM_observables"],
            "mass_ratios": row["mass_ratios"],
            "seesaw_replay": row["seesaw_replay"],
            "ckm_score": row["ckm_score"],
            "mass_score": row["mass_score"],
            "future_margin_1e35_crossed120_kappa30": row["future_margin_1e35_crossed120_kappa30"],
        },
        "crossed_projector": {
            "kappa": 30.0,
            "future_margin_1e35": margin,
            "field_only_unbroken_spin10_projector_possible": crossed["verdict"][
                "field_only_unbroken_spin10_projector_possible"
            ],
            "ps_eft_or_constrained_projector_possible": crossed["verdict"]["ps_eft_or_constrained_projector_possible"],
        },
        "triplet_tensor_source": {
            "source": "output/two_kernel_flavor_then_d5/summary.json::matrix_cards.d5_both_F_minus_0",
            "interpretation": (
                "The source-consistent H,F,G_A,G_B tensors are included so channel-specific "
                "d=5 proton tables can be rebuilt from the publication card instead of an "
                "implicit older scan artifact.  The physical crossed projector remains a "
                "post-Spin(10) source-sector assumption."
            ),
            "model_data": two_kernel["matrix_cards"]["d5_both_F_minus_0"]["model_data"],
            "source_best_triplet_profile": two_kernel["matrix_cards"]["d5_both_F_minus_0"]["best_triplet_profile"],
        },
        "source_and_link_status": {
            "source_symmetry": source_symm["verdict"],
            "dterm_unitary_link": dterm["verdict"],
            "hidden_radial_lock": radial["verdict"],
            "endpoint_vectorlike": endpoint["verdict"],
        },
        "reproducibility_checks": source_checks,
        "gates": gates,
        "legacy_pipeline_context": {
            "old_publication_level_complete": old_pipeline["closure"]["publication_level_complete"],
            "old_reason": old_pipeline["numerical_verdict"]["interpretation"],
        },
        "ledger_status_counts": ledger["counts"],
        "verdict": {
            "local_source_consistent_candidate_passes": local_candidate_passes,
            "publication_level_complete": False,
            "interpretation": (
                "The source-consistent crossed-120 row closes the local strict CKM, mass, seesaw, "
                "and 1e35 yr d=5 stress gates when paired with the post-Spin(10) source projector "
                "and unitary/radial/endpoint completion audits.  This supersedes the older "
                "full_flavor_d5_pipeline bottleneck for the local branch, but it is still not a "
                "publication-final proof because the full channel table and final unified "
                "phenomenology input manifest are not yet rebuilt around this card."
            ),
        },
    }
    return payload


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    obs = payload["recomputed_observables"]
    checks = payload["reproducibility_checks"]
    gates = payload["gates"]
    lines = [
        "# Publication closure candidate card",
        "",
        "No web lookup was used.",
        "",
        "## Selected branch",
        "",
        f"Row: `{payload['selected_row']['label']}`.",
        f"Mode: `{payload['selected_row']['mode']}`.",
        "",
        "## Recomputed local gates",
        "",
        "| gate | pass |",
        "|---|---:|",
    ]
    for key, val in gates.items():
        lines.append(f"| `{key}` | `{val}` |")
    lines += [
        "",
        "## Numerical observables",
        "",
        f"CKM score: `{obs['scores']['ckm_score']:.6e}`.",
        f"Mass score: `{obs['scores']['mass_score']:.6e}`.",
        f"Seesaw residual: `{obs['seesaw_replay']['seesaw_matrix_residual']:.6e}`.",
        f"Future d=5 margin: `{payload['crossed_projector']['future_margin_1e35']:.6e}`.",
        "",
        "## Reproducibility diffs",
        "",
        "| quantity | absolute diff |",
        "|---|---:|",
    ]
    for key, val in checks.items():
        lines.append(f"| `{key}` | {val:.6e} |")
    lines += [
        "",
        "## Verdict",
        "",
        f"Local source-consistent candidate passes: `{payload['verdict']['local_source_consistent_candidate_passes']}`.",
        f"Publication-level complete: `{payload['verdict']['publication_level_complete']}`.",
        "",
        payload["verdict"]["interpretation"],
        "",
        "Machine-readable outputs:",
        "",
        "- `output/publication_closure_card/summary.json`",
        "- `output/publication_closure_card/publication_closure_card.json`",
        "- `output/publication_closure_card/input_manifest.csv`",
    ]
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "publication_closure_card.json").write_text(
        json.dumps(
            {
                "selected_row": payload["selected_row"],
                "recomputed_observables": payload["recomputed_observables"],
                "crossed_projector": payload["crossed_projector"],
                "triplet_tensor_source": payload["triplet_tensor_source"],
                "source_and_link_status": payload["source_and_link_status"],
                "gates": payload["gates"],
                "verdict": payload["verdict"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_csv(OUT / "input_manifest.csv", payload["input_manifest"], ["path", "exists", "size_bytes", "sha256"])
    write_report(payload)
    v = payload["verdict"]
    obs = payload["recomputed_observables"]
    print("Publication closure candidate card")
    print(f"  local candidate passes: {v['local_source_consistent_candidate_passes']}")
    print(f"  publication complete: {v['publication_level_complete']}")
    print(f"  CKM score: {obs['scores']['ckm_score']:.6e}")
    print(f"  mass score: {obs['scores']['mass_score']:.6e}")
    print(f"  future d5 margin: {payload['crossed_projector']['future_margin_1e35']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
