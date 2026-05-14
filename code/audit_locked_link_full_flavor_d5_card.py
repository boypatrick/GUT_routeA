#!/usr/bin/env python3
"""Locked-link full-flavor and d=5 reproducibility card.

No web lookup is used.  This is an aggregator audit: it does not rerun the
flavor or proton scans, but it puts the exact local flavor benchmark, the
threshold-locked unitary-link branch, the mass/field-basis d=5 Wilson pipeline,
and the calibrated Knu width formula into one machine-readable card.

The card deliberately separates three statuses:

  * exact reproducible inputs that already exist;
  * conditional EFT checks that pass for the locked-link branch;
  * publication-level tasks that remain open.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "locked_link_full_flavor_d5_card"

FLAVOR_CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
FLAVOR_AUDIT = ROOT / "output" / "flavor_fit" / "flavor_observable_audit.json"
WILSON_TENSORS = ROOT / "output" / "dimension5_wilson_tensors" / "dimension5_wilson_tensors.json"
PHYSICAL_D5 = ROOT / "output" / "physical_d5_wilson_replay" / "summary.json"
SOFT_D5 = ROOT / "output" / "soft_spectrum_d5_replay" / "summary.json"
KNU_PIPELINE = ROOT / "output" / "full_knu_channel_pipeline" / "summary.json"
KNU_WIDTH = ROOT / "output" / "full_knu_width" / "summary.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
LOCKED_LINK = ROOT / "output" / "threshold_locked_triplet_lift" / "summary.json"
UNITARY_NLSM = ROOT / "output" / "unitary_link_nlsm_action" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict[str, Any]]:
    rows = []
    for path in paths:
        rel = path.relative_to(ROOT)
        rows.append(
            {
                "path": str(rel),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path) if path.exists() else None,
            }
        )
    return rows


def singulars_from_sector(sector: dict[str, Any]) -> list[float]:
    return [float(x) for x in sector["normalized_singular_values"]]


def top_wilson_rows(wilson: dict[str, Any], n: int = 8) -> list[dict[str, Any]]:
    rows = sorted(wilson["rows"], key=lambda row: float(row["S_T_required_tau_1e35"]))
    keep = []
    for row in rows[:n]:
        keep.append(
            {
                "hypothesis": row["hypothesis"],
                "operator": row["operator"],
                "channel_proxy": row["channel_proxy"],
                "selected_index": row["selected_index"],
                "amplitude": float(row["amplitude"]),
                "S_T_required_tau_1e35": float(row["S_T_required_tau_1e35"]),
                "S_T_required_tau_2p4e34": float(row["S_T_required_tau_2p4e34"]),
                "tau_years_ST_1": float(row["tau_years_ST_1"]),
            }
        )
    return keep


def locked_best_rows(locked: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in locked["rows"]:
        rows.append(
            {
                "label": row["label"],
                "target_spectral_norm": float(row["target_spectral_norm"]),
                "min_defect_left": float(row["defect_left_min_eigenvalue"]),
                "M_lock_GeV": float(row["M_lock_GeV"]),
                "singular_spread": float(row["full_mass_singular_spread_over_Mlock"]),
                "threshold_projected_l2": float(row["threshold_projected_l2_if_complete_degenerate"]),
                "min_LLLL_future_margin": float(row["min_LLLL_future_margin"]),
                "min_RRRR_future_margin": float(row["min_RRRR_future_margin"]),
                "passes": bool(row["passes_all_conditional_checks"]),
            }
        )
    return rows


def nlsm_geometry_rows(nlsm: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "label": row["label"],
            "nonempty": bool(row["nonempty_by_contraction_theorem"]),
            "rank": int(row["linearized_fixed_block_rank"]),
            "residual_moduli_real_dimension": int(row["residual_moduli_real_dimension"]),
            "strict_contraction_margin": float(row["strict_contraction_margin"]),
            "worst_future_margin_1e35": float(row["worst_future_margin_1e35"]),
            "threshold_projected_l2_constrained": float(row["threshold_projected_l2_constrained"]),
        }
        for row in nlsm["rows"]
    ]


def build() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    flavor = read_json(FLAVOR_CARD)
    flavor_audit = read_json(FLAVOR_AUDIT)
    wilson = read_json(WILSON_TENSORS)
    physical = read_json(PHYSICAL_D5)
    soft = read_json(SOFT_D5)
    pipeline = read_json(KNU_PIPELINE)
    width = read_json(KNU_WIDTH)
    target = read_json(KNU_TARGET)
    locked = read_json(LOCKED_LINK)
    nlsm = read_json(UNITARY_NLSM)

    sectors = flavor["yukawa_sectors"]
    seesaw = flavor["seesaw"]
    pverdict = pipeline["verdict"]
    wverdict = width["verdict"]
    tverdict = target["verdict"]
    nverdict = nlsm["verdict"]

    exact_yukawa_matrices = {
        name: sectors[name]["Y_matrix_normalized"]
        for name in sorted(sectors)
    }
    singular_summary = {
        name: singulars_from_sector(sectors[name])
        for name in sorted(sectors)
    }

    exact_inputs_ok = bool(flavor["all_checks_ok"])
    pmns_ok = bool(seesaw["all_checks_ok"])
    ckm_completed = not bool(flavor_audit["checks"]["ckm_needs_additional_fit"])
    locked_threshold_silent = bool(
        nverdict["all_threshold_silent_in_constrained_bookkeeping"]
    )
    d5_current_pass = float(tverdict["current_margin"]) >= 1.0
    d5_future_pass = float(tverdict["future_margin"]) >= 1.0
    width_calibration_ok = bool(wverdict["calibration_stable"])

    component_rows = [
        {
            "component": "exact_CP1_O2_yukawa_and_seesaw_card",
            "status": "PASS",
            "key_number": "all_checks_ok={}".format(exact_inputs_ok),
            "evidence_path": str(FLAVOR_CARD.relative_to(ROOT)),
            "interpretation": "Exact normalized Yukawa matrices and seesaw reconstruction are machine-readable.",
        },
        {
            "component": "joint_CKM_PMNS_flavor_fit",
            "status": "OPEN" if not ckm_completed else "PASS",
            "key_number": "CKM_log_score={:.6e}".format(float(flavor_audit["ckm_magnitude_log_score"])),
            "evidence_path": str(FLAVOR_AUDIT.relative_to(ROOT)),
            "interpretation": "PMNS benchmark is reconstructed, but CKM remains the explicit flavor-fit obstruction.",
        },
        {
            "component": "locked_unitary_link_threshold",
            "status": "PASS_CONDITIONAL" if locked_threshold_silent else "FAIL",
            "key_number": "NLSM_margin={:.6e}".format(float(nverdict["minimum_future_margin_1e35"])),
            "evidence_path": str(UNITARY_NLSM.relative_to(ROOT)),
            "interpretation": "The constrained unitary link is nonempty, full-rank, and threshold-silent in conditional bookkeeping.",
        },
        {
            "component": "mass_and_field_basis_d5_replay",
            "status": "PASS_CONDITIONAL",
            "key_number": "pipeline_final_margin={:.6e}".format(float(pverdict["final_future_margin_1e35"])),
            "evidence_path": str(KNU_PIPELINE.relative_to(ROOT)),
            "interpretation": "Existing replays put Wilson tensors, field rotations, and soft dressing on one Knu axis.",
        },
        {
            "component": "calibrated_Knu_width_formula",
            "status": "PASS_CONDITIONAL" if width_calibration_ok else "FAIL",
            "key_number": "Kspread={:.6e}".format(float(width["cross_checks"]["K_dyn_max_relative_spread"])),
            "evidence_path": str(KNU_WIDTH.relative_to(ROOT)),
            "interpretation": "The same width formula describes the soft-spectrum and projected-target Knu rows.",
        },
        {
            "component": "future_stress_d5_proton_safety",
            "status": "OPEN" if not d5_future_pass else "PASS_CONDITIONAL",
            "key_number": "future_suppression_needed={:.6e}".format(float(tverdict["future_amplitude_suppression_needed"])),
            "evidence_path": str(KNU_TARGET.relative_to(ROOT)),
            "interpretation": "The current bound passes, but the uniform 1e35 yr stress still needs stronger suppression.",
        },
    ]

    card = {
        "note": "No web lookup used. Locked-link full flavor and d=5 reproducibility card.",
        "formulae": {
            "unitary_link": "L in U(8), P L P = W_target",
            "d5_LLLL": "C_5L^{abcd}=sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl} U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld}",
            "d5_RRRR": "C_5R^{abcd}=sum_AB (Y_UE^A)_{ij} W_AB (Y_UD^B)_{kl} U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}",
            "width": width["formula"],
        },
        "input_manifest": manifest(
            [
                FLAVOR_CARD,
                FLAVOR_AUDIT,
                WILSON_TENSORS,
                PHYSICAL_D5,
                SOFT_D5,
                KNU_PIPELINE,
                KNU_WIDTH,
                KNU_TARGET,
                LOCKED_LINK,
                UNITARY_NLSM,
            ]
        ),
        "exact_flavor_inputs": {
            "basis": flavor["basis"],
            "bundle_checks": flavor["bundle_checks"],
            "Y_matrices_normalized": exact_yukawa_matrices,
            "singular_values_normalized": singular_summary,
            "PMNS_definition": seesaw["basis_convention"]["pmns_definition"],
            "seesaw_checks": seesaw["checks"],
            "heavy_neutrino_masses_GeV": seesaw["heavy_neutrino_masses_GeV"],
            "MR_GeV": seesaw["MR_GeV"],
            "mD_GeV": seesaw["mD_GeV"],
        },
        "flavor_fit_status": {
            "PMNS_reconstructed": pmns_ok,
            "CKM_completed": ckm_completed,
            "CKM_log_score": float(flavor_audit["ckm_magnitude_log_score"]),
            "CKM_current_observables": flavor_audit["ckm_current_observables"],
            "CKM_target_observables": flavor_audit["ckm_target_observables"],
            "mass_ratio_audit": flavor_audit["mass_ratio_audit"],
            "interpretation": flavor_audit["verdict"],
        },
        "locked_link_status": {
            "mathematical_identity": locked["mathematical_identity"],
            "locked_rows": locked_best_rows(locked),
            "nlsm_geometry": nlsm_geometry_rows(nlsm),
            "verdict": nverdict,
        },
        "d5_status": {
            "operator_definitions": wilson["operator_definitions"],
            "normalization": wilson["normalization"],
            "most_dangerous_mass_basis_rows": top_wilson_rows(wilson),
            "physical_field_basis_verdict": physical["verdict"],
            "soft_spectrum_verdict": soft["verdict"],
            "pipeline_verdict": pverdict,
            "width_verdict": wverdict,
            "target_verdict": tverdict,
        },
        "component_status": component_rows,
        "verdict": {
            "exact_inputs_available": exact_inputs_ok,
            "pmns_reconstructed": pmns_ok,
            "ckm_fit_completed": ckm_completed,
            "locked_link_threshold_silent": locked_threshold_silent,
            "d5_current_bound_passes": d5_current_pass,
            "d5_future_1e35_passes": d5_future_pass,
            "d5_future_suppression_needed": float(tverdict["future_amplitude_suppression_needed"]),
            "width_calibration_stable": width_calibration_ok,
            "publication_level_complete": bool(
                exact_inputs_ok
                and pmns_ok
                and ckm_completed
                and locked_threshold_silent
                and d5_current_pass
                and d5_future_pass
                and width_calibration_ok
            ),
            "interpretation": (
                "This card makes the current branch reproducible at the matrix/card level. "
                "It does not close the paper: CKM remains unfitted and the conservative Knu "
                "future-stress target still needs an additional amplitude suppression of "
                "{:.3f} or an equivalent stronger triplet filter."
            ).format(float(tverdict["future_amplitude_suppression_needed"])),
        },
    }
    return card, component_rows


def write_component_csv(rows: list[dict[str, Any]]) -> None:
    with (OUT / "component_status.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["component", "status", "key_number", "evidence_path", "interpretation"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_manifest_csv(rows: list[dict[str, Any]]) -> None:
    with (OUT / "input_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "exists", "size_bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def write_report(card: dict[str, Any]) -> None:
    verdict = card["verdict"]
    flavor = card["flavor_fit_status"]
    d5 = card["d5_status"]
    nlsm = card["locked_link_status"]["verdict"]
    lines = [
        "# Locked-link full flavor and d=5 reproducibility card",
        "",
        "No web lookup was used.",
        "",
        "This card aggregates the exact local Yukawa/seesaw matrices, the locked",
        "unitary-link branch, and the Knu-focused d=5 proton-decay replay into one",
        "manifest.  It is a reproducibility scaffold, not a publication-level",
        "phenomenology closure.",
        "",
        "## Status",
        "",
        "| component | status | key number |",
        "|---|---|---:|",
    ]
    for row in card["component_status"]:
        lines.append(
            f"| `{row['component']}` | `{row['status']}` | `{row['key_number']}` |"
        )
    lines += [
        "",
        "## Key Numbers",
        "",
        f"- CKM magnitude log-score: `{flavor['CKM_log_score']:.6e}`.",
        f"- Constrained-link minimum future margin: `{nlsm['minimum_future_margin_1e35']:.6e}`.",
        (
            "- Knu final future margin: "
            f"`{d5['target_verdict']['future_margin']:.6e}`; "
            "suppression needed: "
            f"`{d5['target_verdict']['future_amplitude_suppression_needed']:.6e}`."
        ),
        f"- Knu width calibration stable: `{d5['width_verdict']['calibration_stable']}`.",
        "",
        "## Verdict",
        "",
        verdict["interpretation"],
        "",
        "Machine-readable files:",
        "",
        "- `output/locked_link_full_flavor_d5_card/locked_link_full_flavor_d5_card.json`",
        "- `output/locked_link_full_flavor_d5_card/component_status.csv`",
        "- `output/locked_link_full_flavor_d5_card/input_manifest.csv`",
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card, component_rows = build()
    (OUT / "locked_link_full_flavor_d5_card.json").write_text(
        json.dumps(card, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "note": card["note"],
                "component_status": card["component_status"],
                "verdict": card["verdict"],
                "input_manifest": card["input_manifest"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_component_csv(component_rows)
    write_manifest_csv(card["input_manifest"])
    write_report(card)
    verdict = card["verdict"]
    print("Locked-link full flavor + d=5 reproducibility card")
    print(f"  exact inputs available: {verdict['exact_inputs_available']}")
    print(f"  CKM fit completed: {verdict['ckm_fit_completed']}")
    print(f"  d5 current bound passes: {verdict['d5_current_bound_passes']}")
    print(f"  d5 future 1e35 passes: {verdict['d5_future_1e35_passes']}")
    print(f"  suppression needed: {verdict['d5_future_suppression_needed']:.6e}")
    print(f"  publication complete: {verdict['publication_level_complete']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
