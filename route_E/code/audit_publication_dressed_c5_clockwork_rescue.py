#!/usr/bin/env python3
"""Replay the dressed publication C5 audit at the clockwork value kappa=729.

No web lookup is used.

The earlier eigenstate-card replay showed a real failure for the finite
``kappa=30`` crossed-triplet inverse block.  The subsequent clockwork/hidden
endpoint audits produced a conditional source mechanism with

    q=3, n=6, kappa=q^n=729,

and an 8x8 unitary dilation whose visible threshold vector is zero.  This
script closes the loop by recomputing the same dressed-channel stress grid
directly at kappa=729, rather than only extrapolating from the coarse kappa
scan.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import scan_publication_dressed_c5_kappa as kscan


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "publication_dressed_c5_clockwork_rescue"

INPUTS = {
    "finite_kappa_failure": ROOT / "output" / "publication_dressed_c5_from_eigenstate_card" / "summary.json",
    "kappa_scan": ROOT / "output" / "publication_dressed_c5_kappa_scan" / "summary.json",
    "rank_one_clockwork": ROOT / "output" / "rank_one_clockwork_locking" / "summary.json",
    "constrained_clockwork": ROOT / "output" / "constrained_clockwork_source_hessian" / "summary.json",
    "unitary_link": ROOT / "output" / "clockwork_unitary_link_quotient" / "summary.json",
    "endpoint_meson": ROOT / "output" / "clockwork_hidden_endpoint_meson" / "summary.json",
    "endpoint_vectorlike": ROOT / "output" / "clockwork_endpoint_vectorlike_completion" / "summary.json",
    "radial_driver": ROOT / "output" / "clockwork_radial_driver_hessian" / "summary.json",
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


def manifest() -> list[dict[str, Any]]:
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


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def audit() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    failure = read_json(INPUTS["finite_kappa_failure"])
    kappa_scan = read_json(INPUTS["kappa_scan"])
    rank = read_json(INPUTS["rank_one_clockwork"])
    constrained = read_json(INPUTS["constrained_clockwork"])
    unitary = read_json(INPUTS["unitary_link"])
    meson = read_json(INPUTS["endpoint_meson"])
    vectorlike = read_json(INPUTS["endpoint_vectorlike"])
    radial = read_json(INPUTS["radial_driver"])
    hidden = read_json(INPUTS["hidden_quotient"])

    kappa = float(constrained["input_card"]["kappa"])
    q = float(constrained["input_card"]["q"])
    n = int(constrained["input_card"]["n"])
    required = float(kappa_scan["verdict"]["max_width_min_kappa_pass_1e35"])

    base = kscan.dressed.read_json(kscan.dressed.CARD)
    rows = kscan.rows_for_kappa(base, kappa)
    exact = kscan.summarize_kappa(rows, kappa)

    old_finite = failure["by_block_case"]["finite"]
    direct_clockwork_checks = {
        "q_n_equals_kappa": abs(q**n - kappa) < 1.0e-12,
        "kappa_exceeds_required": kappa >= required,
        "constrained_clockwork_condition_matches": constrained["effective_inverse_block"]["matches_required_kappa"],
        "constrained_clockwork_has_no_closed_zero_modes": (
            constrained["linear_algebra"]["boundary_closed_zero_modes_per_chain"] == 0
        ),
        "unitary_link_locks_physical_masses": unitary["verdict"]["mass_singular_values_locked"],
        "unitary_link_threshold_zero": (
            unitary["threshold_interpretation"]["constrained_source_threshold"] == [0.0, 0.0, 0.0]
        ),
        "meson_threshold_zero": meson["verdict"]["visible_threshold_vector"] == [0.0, 0.0, 0.0],
        "endpoint_vectorlike_safe": vectorlike["verdict"]["vectorlike_completion_cancels_all_hidden_cubic_anomalies"],
        "radial_driver_lifts_nonunitary_modes": radial["verdict"]["radial_driver_lifts_nonunitary_modes"],
        "hidden_quotient_moment_maps_match": hidden["verdict"]["moment_maps_match_radial_lock"],
        "hidden_quotient_threshold_zero": hidden["verdict"]["visible_threshold_vector"] == [0.0, 0.0, 0.0],
    }

    exact_passes = (
        exact["central_unsafe_1e35_rows"] == 0
        and exact["max_width_unsafe_1e35_rows"] == 0
    )
    hidden_branch_passes = all(direct_clockwork_checks.values())
    summary = {
        "note": (
            "No web lookup used. Exact dressed-channel replay at the q=3,n=6 "
            "clockwork value kappa=729."
        ),
        "input_manifest": manifest(),
        "clockwork_card": {
            "q": q,
            "n": n,
            "kappa": kappa,
            "epsilon": float(constrained["input_card"]["epsilon"]),
            "required_kappa_from_previous_scan": required,
            "M_lock_GeV": float(constrained["input_card"]["M_lock_GeV"]),
        },
        "old_finite_kappa30_failure": {
            "central_unsafe_1e35_rows": old_finite["central"]["unsafe_1e35_rows"],
            "central_global_STmax_1e35": old_finite["central"]["global_S_T_max_1e35"],
            "max_width_unsafe_1e35_rows": old_finite["max_width"]["unsafe_1e35_rows"],
            "max_width_global_STmax_1e35": old_finite["max_width"]["global_S_T_max_1e35"],
        },
        "exact_kappa729_replay": exact,
        "direct_clockwork_checks": direct_clockwork_checks,
        "verdict": {
            "status": "PASS_CONDITIONAL" if exact_passes and hidden_branch_passes else "OPEN",
            "exact_dressed_grid_passes": exact_passes,
            "hidden_clockwork_branch_passes": hidden_branch_passes,
            "central_margin_1e35_at_ST_1e_minus_5": exact["central_worst_margin_1e35"],
            "max_width_margin_1e35_at_ST_1e_minus_5": exact["max_width_worst_margin_1e35"],
            "min_width_margin_1e35_at_ST_1e_minus_5": exact["min_width_worst_margin_1e35"],
            "remaining_caveat": (
                "This closes the local publication-card dressed proxy at the "
                "clockwork-rescued inverse block.  It is still conditional on "
                "the hidden Kahler/FI clockwork sector and does not by itself "
                "replace a final paper-grade hadronic/chiral input manifest."
            ),
        },
    }
    return rows, summary


def write_report(summary: dict[str, Any]) -> None:
    exact = summary["exact_kappa729_replay"]
    old = summary["old_finite_kappa30_failure"]
    card = summary["clockwork_card"]
    lines = [
        "# Clockwork-rescued publication dressed C5 replay",
        "",
        "No web lookup was used.",
        "",
        "## Input",
        "",
        f"- q = {card['q']:.6g}",
        f"- n = {card['n']}",
        f"- kappa = q^n = {card['kappa']:.6g}",
        f"- required kappa from previous max-width scan = {card['required_kappa_from_previous_scan']:.6g}",
        f"- reference filter S_T = {kscan.dressed.DISPLAY_ST:.1e}",
        "",
        "## Old finite-block failure",
        "",
        "| block | central unsafe | central S_T max | max-width unsafe | max-width S_T max |",
        "|---|---:|---:|---:|---:|",
        (
            f"| kappa=30 finite | {old['central_unsafe_1e35_rows']} | "
            f"{old['central_global_STmax_1e35']:.6e} | "
            f"{old['max_width_unsafe_1e35_rows']} | "
            f"{old['max_width_global_STmax_1e35']:.6e} |"
        ),
        "",
        "## Exact kappa=729 replay",
        "",
        "| case | rows | unsafe 1e35 | unsafe present | global S_T max 1e35 | worst margin at S_T=1e-5 | worst channel |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for case in ["central", "max_width", "min_width"]:
        lines.append(
            f"| {case} | {exact[f'{case}_rows']} | "
            f"{exact[f'{case}_unsafe_1e35_rows']} | "
            f"{exact[f'{case}_unsafe_present_rows']} | "
            f"{exact[f'{case}_global_STmax_1e35']:.6e} | "
            f"{exact[f'{case}_worst_margin_1e35']:.6e} | "
            f"{exact[f'{case}_worst_channel']} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        f"Status: `{summary['verdict']['status']}`.",
        "",
        summary["verdict"]["remaining_caveat"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = audit()
    write_csv(OUT / "dressed_channel_rows_kappa729.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    print("Clockwork-rescued publication dressed C5 replay")
    print(f"  status: {summary['verdict']['status']}")
    print(f"  central margin: {summary['verdict']['central_margin_1e35_at_ST_1e_minus_5']:.6g}")
    print(f"  max-width margin: {summary['verdict']['max_width_margin_1e35_at_ST_1e_minus_5']:.6g}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
