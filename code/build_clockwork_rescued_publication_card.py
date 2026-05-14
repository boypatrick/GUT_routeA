#!/usr/bin/env python3
"""Build one local publication card for the clockwork-rescued flavor+d5 branch.

No web lookup is used.

This is a packaging and reproducibility audit, not a new physical assumption.
It merges the source-consistent flavor card with the exact kappa=729 dressed
C5 replay and emits a single manifest with hashes, flavor observables, seesaw
data, channel worst rows, and explicit remaining caveats.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "clockwork_rescued_publication_card"

INPUTS = {
    "publication_closure_card": ROOT / "output" / "publication_closure_card" / "publication_closure_card.json",
    "flavor_reproducibility": ROOT / "output" / "publication_flavor_d5_reproducibility" / "summary.json",
    "triplet_eigenstate_card": ROOT / "output" / "publication_triplet_eigenstate_card" / "triplet_eigenstate_card.json",
    "clockwork_d5_rescue": ROOT / "output" / "publication_dressed_c5_clockwork_rescue" / "summary.json",
    "clockwork_d5_rows": ROOT / "output" / "publication_dressed_c5_clockwork_rescue" / "dressed_channel_rows_kappa729.csv",
    "hidden_quotient": ROOT / "output" / "clockwork_hidden_gauge_quotient_origin" / "summary.json",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def input_manifest() -> list[dict[str, Any]]:
    rows = []
    for label, path in INPUTS.items():
        rows.append(
            {
                "label": label,
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path) if path.exists() else None,
            }
        )
    return rows


def read_channel_rows() -> list[dict[str, Any]]:
    path = INPUTS["clockwork_d5_rows"]
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "amplitude_with_dressing",
                "M_triplet_GeV",
                "S_T_display",
                "width_prefactor_GeV5",
                "C6_GeV_minus2_at_ST_display",
                "tau_years_at_ST_display",
                "present_bound_years",
                "margin_present_at_ST_display",
                "margin_1e35_at_ST_display",
                "S_T_max_present",
                "S_T_max_1e35",
                "avg_chargino_part",
                "avg_neutralino_part",
                "avg_higgsino_part",
                "avg_total_part",
                "kappa",
            ]:
                row[key] = float(row[key])
            for key in ["passes_present_at_ST_display", "passes_1e35_at_ST_display"]:
                row[key] = row[key] == "True"
            rows.append(row)
    return rows


def worst_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["normalization_case"], row["operator"], row["channel"])
        grouped.setdefault(key, []).append(row)

    out = []
    for key, items in sorted(grouped.items()):
        worst_1e35 = min(items, key=lambda row: row["margin_1e35_at_ST_display"])
        worst_present = min(items, key=lambda row: row["margin_present_at_ST_display"])
        out.append(
            {
                "normalization_case": key[0],
                "operator": key[1],
                "channel": key[2],
                "rows": len(items),
                "unsafe_1e35_rows": sum(not row["passes_1e35_at_ST_display"] for row in items),
                "unsafe_present_rows": sum(not row["passes_present_at_ST_display"] for row in items),
                "worst_margin_1e35": worst_1e35["margin_1e35_at_ST_display"],
                "worst_tau_years_1e35_row": worst_1e35["tau_years_at_ST_display"],
                "S_T_max_1e35": worst_1e35["S_T_max_1e35"],
                "worst_margin_present": worst_present["margin_present_at_ST_display"],
                "S_T_max_present": worst_present["S_T_max_present"],
                "worst_scenario": worst_1e35["scenario"],
                "worst_epsilon": worst_1e35["epsilon"],
                "worst_spectrum": worst_1e35["spectrum_name"],
                "worst_pair": worst_1e35["pair"],
                "worst_selected_index": worst_1e35["selected_index"],
                "worst_amplitude_with_dressing": worst_1e35["amplitude_with_dressing"],
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def flavor_rows(flavor: dict[str, Any]) -> list[dict[str, Any]]:
    fresh = flavor["fresh_recompute"]
    rows = [
        {"sector": "CKM", "observable": key, "value": val}
        for key, val in fresh["CKM_observables"].items()
    ]
    rows.extend(
        [
            {"sector": "score", "observable": key, "value": val}
            for key, val in fresh["scores"].items()
        ]
    )
    for sector, vals in fresh["mass_ratios"].items():
        rows.append({"sector": sector, "observable": "small", "value": vals["small"]})
        rows.append({"sector": sector, "observable": "mid", "value": vals["mid"]})
    seesaw = fresh["seesaw_replay"]
    for key in ["seesaw_matrix_residual", "theta_norm", "MR_condition_number"]:
        rows.append({"sector": "seesaw", "observable": key, "value": seesaw[key]})
    for idx, mass in enumerate(seesaw["heavy_neutrino_masses_GeV"], start=1):
        rows.append({"sector": "seesaw", "observable": f"M_R{idx}_GeV", "value": mass})
    return rows


def build() -> dict[str, Any]:
    closure = read_json(INPUTS["publication_closure_card"])
    flavor = read_json(INPUTS["flavor_reproducibility"])
    rescue = read_json(INPUTS["clockwork_d5_rescue"])
    hidden = read_json(INPUTS["hidden_quotient"])
    all_rows = read_channel_rows()
    channel_worst = worst_rows(all_rows)

    flavor_gate = bool(flavor["verdict"]["local_source_consistent_candidate_passes"])
    d5_gate = bool(rescue["verdict"]["exact_dressed_grid_passes"])
    hidden_gate = bool(rescue["verdict"]["hidden_clockwork_branch_passes"])
    all_channel_rows_pass = all(row["unsafe_1e35_rows"] == 0 for row in channel_worst)
    worst_global = min(channel_worst, key=lambda row: row["worst_margin_1e35"])

    remaining = [
        "External literature-refresh and final cited input table are not done because this run used no web lookup.",
        "The d=5 calculation is still the local dressed proxy, not a full independent lattice/chiral publication review.",
        "The hidden clockwork/Kahler/FI quotient is a conditional EFT assumption rather than a microscopic first-principles derivation.",
        "Full CKM/PMNS target replacement by a final cited global-fit table remains open.",
    ]
    return {
        "note": "No web lookup used. Single local clockwork-rescued flavor+d5 publication manifest.",
        "input_manifest": input_manifest(),
        "clockwork_card": rescue["clockwork_card"],
        "flavor_summary": {
            "CKM_observables": flavor["fresh_recompute"]["CKM_observables"],
            "scores": flavor["fresh_recompute"]["scores"],
            "mass_ratios": flavor["fresh_recompute"]["mass_ratios"],
            "seesaw_replay": flavor["fresh_recompute"]["seesaw_replay"],
            "gates": flavor["gates"],
        },
        "d5_summary": {
            "total_channel_rows": len(all_rows),
            "channel_worst_table_rows": len(channel_worst),
            "all_channel_rows_pass_1e35": all_channel_rows_pass,
            "global_worst": worst_global,
            "exact_kappa729_replay": rescue["exact_kappa729_replay"],
        },
        "hidden_clockwork_summary": {
            "direct_clockwork_checks": rescue["direct_clockwork_checks"],
            "moment_maps_match_radial_lock": hidden["verdict"]["moment_maps_match_radial_lock"],
            "visible_threshold_vector": hidden["verdict"]["visible_threshold_vector"],
        },
        "gates": {
            "flavor_card_reproducible": flavor["verdict"]["local_source_consistent_candidate_reproducible"],
            "flavor_local_pass": flavor_gate,
            "d5_exact_dressed_grid_pass": d5_gate,
            "d5_all_channel_worst_rows_pass": all_channel_rows_pass,
            "hidden_clockwork_branch_pass": hidden_gate,
        },
        "scope": {
            "claim": "single local conditional flavor+d5 manifest for the clockwork-rescued crossed-120 branch",
            "not_claimed": remaining,
        },
        "verdict": {
            "local_clockwork_rescued_card_complete": bool(flavor_gate and d5_gate and hidden_gate and all_channel_rows_pass),
            "publication_level_complete": False,
            "remaining_open_items": [
                "Full CKM/PMNS/flavor fit is completed",
                "External cited d=5 input refresh is completed",
                "A5-A6 have a microscopic first-principles origin",
            ],
            "interpretation": (
                "The local flavor, seesaw, hidden-clockwork, and dressed d=5 proton rows now live in one "
                "hash-locked manifest.  This closes the local cache-synchronization problem and gives a "
                "conditional publication-card candidate, but it is not a final paper-grade proof until the "
                "external input table, full CKM/PMNS targets, and hadronic/chiral convention ledger are finalized."
            ),
        },
    }


def write_report(card: dict[str, Any]) -> None:
    f = card["flavor_summary"]
    d5 = card["d5_summary"]
    worst = d5["global_worst"]
    lines = [
        "# Clockwork-rescued publication card",
        "",
        "No web lookup was used.",
        "",
        "## Flavor",
        "",
        f"- CKM score: {f['scores']['ckm_score']:.6e}",
        f"- mass score: {f['scores']['mass_score']:.6e}",
        f"- seesaw residual: {f['seesaw_replay']['seesaw_matrix_residual']:.6e}",
        f"- heavy neutrino masses GeV: {', '.join(f'{x:.6e}' for x in f['seesaw_replay']['heavy_neutrino_masses_GeV'])}",
        "",
        "## D5 Proton",
        "",
        f"- total dressed channel rows: {d5['total_channel_rows']}",
        f"- channel/operator/normalization worst rows: {d5['channel_worst_table_rows']}",
        f"- all rows pass 1e35 yr stress: `{d5['all_channel_rows_pass_1e35']}`",
        f"- global worst: {worst['normalization_case']} {worst['operator']} {worst['channel']}",
        f"- global worst margin at S_T=1e-5: {worst['worst_margin_1e35']:.6e}",
        "",
        "## Gates",
        "",
        "| gate | pass |",
        "|---|---:|",
    ]
    for key, val in card["gates"].items():
        lines.append(f"| `{key}` | `{val}` |")
    lines += [
        "",
        "## Remaining Scope Limits",
        "",
    ]
    for item in card["scope"]["not_claimed"]:
        lines.append(f"- {item}")
    lines += [
        "",
        "## Verdict",
        "",
        f"Local card complete: `{card['verdict']['local_clockwork_rescued_card_complete']}`.",
        f"Publication complete: `{card['verdict']['publication_level_complete']}`.",
        "",
        card["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = build()
    rows = read_channel_rows()
    channel_worst = worst_rows(rows)
    frows = flavor_rows(read_json(INPUTS["flavor_reproducibility"]))
    (OUT / "clockwork_rescued_publication_card.json").write_text(
        json.dumps(card, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(OUT / "channel_worst_rows.csv", channel_worst)
    write_csv(OUT / "flavor_observables.csv", frows)
    write_csv(OUT / "input_manifest.csv", card["input_manifest"])
    write_report(card)
    print("Clockwork-rescued publication card")
    print(f"  local card complete: {card['verdict']['local_clockwork_rescued_card_complete']}")
    print(f"  publication complete: {card['verdict']['publication_level_complete']}")
    print(f"  worst d5 margin: {card['d5_summary']['global_worst']['worst_margin_1e35']:.6g}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
