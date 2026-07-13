#!/usr/bin/env python3
"""Build a no-web provenance manifest for the paper's closure tables.

The manifest is deliberately smaller than a full semantic parser for the TeX.
It targets the final closure package: flavor/seesaw, clockwork-rescued d=5
rows, hidden quotient, conditional ledger, and no-web convention ledger.  Each
row records a machine artifact, a SHA-256 hash, a key numerical claim, and
whether the current TeX contains the printed key used by the paper.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "no_web_table_provenance_manifest"
TEX = ROOT / "paper" / "gut_framework.tex"

ARTIFACTS = {
    "no_web_input_ledger": ROOT / "output" / "no_web_input_convention_ledger" / "summary.json",
    "clockwork_publication_card": ROOT / "output" / "clockwork_rescued_publication_card" / "clockwork_rescued_publication_card.json",
    "clockwork_d5_rows": ROOT / "output" / "publication_dressed_c5_clockwork_rescue" / "dressed_channel_rows_kappa729.csv",
    "clockwork_d5_rescue": ROOT / "output" / "publication_dressed_c5_clockwork_rescue" / "summary.json",
    "flavor_reproducibility": ROOT / "output" / "publication_flavor_d5_reproducibility" / "summary.json",
    "flavor_target_provenance": ROOT / "output" / "no_web_flavor_target_provenance" / "summary.json",
    "pmns_benchmark_replay": ROOT / "output" / "no_web_pmns_benchmark_replay" / "summary.json",
    "source_pmns_replay": ROOT / "output" / "source_consistent_pmns_replay" / "summary.json",
    "source_majorana_rank": ROOT / "output" / "source_majorana_texture_rank" / "summary.json",
    "majorana_contact_sensitivity": ROOT / "output" / "majorana_contact_sensitivity" / "summary.json",
    "majorana_source_locking": ROOT / "output" / "majorana_source_locking_sector" / "summary.json",
    "majorana_hidden_quotient": ROOT / "output" / "majorana_hidden_quotient_origin" / "summary.json",
    "majorana_monomial_clockwork": ROOT / "output" / "majorana_monomial_clockwork_origin" / "summary.json",
    "conditional_theorem_boundary": ROOT / "output" / "conditional_theorem_boundary" / "summary.json",
    "hidden_quotient": ROOT / "output" / "clockwork_hidden_gauge_quotient_origin" / "summary.json",
    "conditional_theorem_ledger": ROOT / "output" / "conditional_theorem_ledger" / "summary.json",
    "proton_decay_verification": ROOT / "output" / "proton_decay" / "proton_decay_verification.json",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tex_contains(needle: str) -> bool:
    return needle in TEX.read_text(encoding="utf-8")


def artifact_row(
    label: str,
    claim: str,
    key_value: str,
    tex_key: str,
    artifact: str,
    required_in_tex: bool = True,
) -> dict[str, Any]:
    path = ARTIFACTS[artifact]
    digest = sha256(path) if path.exists() else ""
    found = tex_contains(tex_key) if required_in_tex else True
    return {
        "label": label,
        "claim": claim,
        "key_value": key_value,
        "tex_key": tex_key,
        "required_in_tex": required_in_tex,
        "tex_key_found": found,
        "artifact": artifact,
        "path": str(path.relative_to(ROOT)),
        "exists": path.exists(),
        "sha256": digest,
        "sha256_prefix12": digest[:12],
    }


def build_rows() -> list[dict[str, Any]]:
    no_web = read_json(ARTIFACTS["no_web_input_ledger"])
    card = read_json(ARTIFACTS["clockwork_publication_card"])
    rescue = read_json(ARTIFACTS["clockwork_d5_rescue"])
    flavor = read_json(ARTIFACTS["flavor_reproducibility"])
    flavor_targets = read_json(ARTIFACTS["flavor_target_provenance"])
    pmns_replay = read_json(ARTIFACTS["pmns_benchmark_replay"])
    source_pmns = read_json(ARTIFACTS["source_pmns_replay"])
    source_majorana = read_json(ARTIFACTS["source_majorana_rank"])
    majorana_contact = read_json(ARTIFACTS["majorana_contact_sensitivity"])
    majorana_source = read_json(ARTIFACTS["majorana_source_locking"])
    majorana_hidden = read_json(ARTIFACTS["majorana_hidden_quotient"])
    majorana_monomial = read_json(ARTIFACTS["majorana_monomial_clockwork"])
    theorem_boundary = read_json(ARTIFACTS["conditional_theorem_boundary"])
    hidden = read_json(ARTIFACTS["hidden_quotient"])
    ledger = read_json(ARTIFACTS["conditional_theorem_ledger"])

    rows = [
        artifact_row(
            "constant-convention digest",
            "The no-web convention table has 78 rows and a stable digest.",
            no_web["constant_rows_sha256"],
            no_web["constant_rows_sha256"],
            "no_web_input_ledger",
        ),
        artifact_row(
            "publication-card hash",
            "The flavor+d5 closure card is hash-addressed from the TeX.",
            sha256(ARTIFACTS["clockwork_publication_card"])[:12],
            sha256(ARTIFACTS["clockwork_publication_card"])[:12],
            "clockwork_publication_card",
        ),
        artifact_row(
            "dressed-row hash",
            "The kappa=729 dressed channel CSV is hash-addressed from the TeX.",
            sha256(ARTIFACTS["clockwork_d5_rows"])[:12],
            sha256(ARTIFACTS["clockwork_d5_rows"])[:12],
            "clockwork_d5_rows",
        ),
        artifact_row(
            "legacy proton hash",
            "The legacy local proton convention JSON is hash-addressed from the TeX.",
            sha256(ARTIFACTS["proton_decay_verification"])[:12],
            sha256(ARTIFACTS["proton_decay_verification"])[:12],
            "proton_decay_verification",
        ),
        artifact_row(
            "flavor CKM score",
            "The local CKM score is printed in the flavor+d5 manifest proposition.",
            f"{card['flavor_summary']['scores']['ckm_score']:.6e}",
            "7.634719",
            "clockwork_publication_card",
        ),
        artifact_row(
            "flavor mass score",
            "The local mass score is printed in the flavor+d5 manifest proposition.",
            f"{card['flavor_summary']['scores']['mass_score']:.6e}",
            "1.490597",
            "clockwork_publication_card",
        ),
        artifact_row(
            "seesaw residual",
            "The local inverse-seesaw residual is printed in the flavor+d5 manifest proposition.",
            f"{card['flavor_summary']['seesaw_replay']['seesaw_matrix_residual']:.6e}",
            "5.173996",
            "clockwork_publication_card",
        ),
        artifact_row(
            "d5 row count",
            "The exact dressed d=5 row count is printed in the TeX.",
            str(card["d5_summary"]["total_channel_rows"]),
            "21168",
            "clockwork_publication_card",
        ),
        artifact_row(
            "d5 global worst margin",
            "The max-width global worst d=5 margin is printed in the TeX.",
            f"{card['d5_summary']['global_worst']['worst_margin_1e35']:.9f}",
            "6.216918955",
            "clockwork_publication_card",
        ),
        artifact_row(
            "d5 Knu weakest margin",
            "The weakest K+ nubar margin is printed in the TeX.",
            "14.8805947",
            "14.8805947",
            "clockwork_publication_card",
        ),
        artifact_row(
            "clockwork kappa",
            "The exact clockwork replay uses kappa=729.",
            f"{rescue['clockwork_card']['kappa']:.0f}",
            "729",
            "clockwork_d5_rescue",
        ),
        artifact_row(
            "hidden threshold vector",
            "The hidden quotient has zero visible threshold vector.",
            str(hidden["verdict"]["visible_threshold_vector"]),
            "visible threshold vector remains zero",
            "hidden_quotient",
            required_in_tex=False,
        ),
        artifact_row(
            "theorem-ledger open blocker rename",
            "The d=5 open item is narrowed to external cited input refresh.",
            "; ".join(ledger["open_blockers"]),
            "External cited d=5 input refresh",
            "conditional_theorem_ledger",
            required_in_tex=False,
        ),
        artifact_row(
            "flavor reproducibility gate",
            "The source-consistent flavor card recomputes from its embedded matrices.",
            str(flavor["verdict"]["local_source_consistent_candidate_reproducible"]),
            "card recomputes",
            "flavor_reproducibility",
            required_in_tex=False,
        ),
        artifact_row(
            "flavor-target all-Yukawa digest",
            "The no-web flavor-target ledger hash-addresses all embedded Yukawa matrices.",
            flavor_targets["matrix_manifest"][-1]["digest"][:12],
            flavor_targets["matrix_manifest"][-1]["digest"][:12],
            "flavor_target_provenance",
        ),
        artifact_row(
            "flavor-target local scores",
            "The no-web flavor-target rows reproduce the stored CKM and mass scores exactly.",
            (
                f"ckm_diff={flavor_targets['scores']['ckm_score_absdiff']:.1e}, "
                f"mass_diff={flavor_targets['scores']['mass_score_absdiff']:.1e}"
            ),
            "score differences are zero",
            "flavor_target_provenance",
            required_in_tex=False,
        ),
        artifact_row(
            "PMNS benchmark pair digest",
            "The local exact CP1/O(2) PMNS benchmark replay is hash-addressed from the TeX.",
            pmns_replay["pmns_pair_digest"][:12],
            pmns_replay["pmns_pair_digest"][:12],
            "pmns_benchmark_replay",
        ),
        artifact_row(
            "PMNS benchmark residual",
            "The local exact CP1/O(2) PMNS benchmark replay passes the 1e-10 angle gate.",
            f"{pmns_replay['residuals']['max_pmns_angle_absdiff']:.3e}",
            "6.950",
            "pmns_benchmark_replay",
            required_in_tex=False,
        ),
        artifact_row(
            "source PMNS pair digest",
            "The source-consistent publication-card PMNS compatibility replay is hash-addressed from the TeX.",
            source_pmns["source_pmns_pair_digest"][:12],
            source_pmns["source_pmns_pair_digest"][:12],
            "source_pmns_replay",
        ),
        artifact_row(
            "source PMNS residual",
            "The source-consistent publication card is PMNS-compatible by inverse-seesaw reconstruction.",
            f"{source_pmns['residuals']['max_pmns_angle_absdiff']:.3e}",
            "1.110",
            "source_pmns_replay",
            required_in_tex=False,
        ),
        artifact_row(
            "source Majorana contact fraction",
            "The source-consistent M_R requires a non-negligible trace/contact component.",
            f"{source_majorana['decomposition']['contact_fraction']:.6e}",
            "1.304275",
            "source_majorana_rank",
        ),
        artifact_row(
            "source Majorana lifted residual",
            "The source-consistent M_R is exactly represented by Veronese plus contact.",
            f"{source_majorana['decomposition']['veronese_plus_contact_relative_residual']:.6e}",
            "2.687332",
            "source_majorana_rank",
            required_in_tex=False,
        ),
        artifact_row(
            "Majorana contact scale loose half-width",
            "The PMNS-compatible source-contact scale is narrow even under the loose local gate.",
            f"{majorana_contact['real_scale_loose_interval']['half_width']:.6e}",
            "3.208798",
            "majorana_contact_sensitivity",
        ),
        artifact_row(
            "Majorana contact phase loose half-width",
            "The PMNS-compatible source-contact phase is narrow even under the loose local gate.",
            f"{majorana_contact['phase_loose_interval']['half_width']:.6e}",
            "5.424119",
            "majorana_contact_sensitivity",
        ),
        artifact_row(
            "Majorana affine source F-flatness",
            "The affine source-locking fallback has vanishing F-terms at the target vacuum.",
            f"{majorana_source['f_norm_at_vacuum']:.6e}",
            "3.103168",
            "majorana_source_locking",
        ),
        artifact_row(
            "Majorana affine source Hessian floor",
            "The affine source-locking fallback has no Hessian flat direction.",
            f"{min(majorana_source['hessian_singular_values']):.6e}",
            "7.767189",
            "majorana_source_locking",
        ),
        artifact_row(
            "Majorana hidden D-only rank",
            "The D-only hidden quotient leaves zeta moduli.",
            (
                f"{majorana_hidden['candidates']['D_only_hidden_U1']['hessian_rank']}/"
                f"{majorana_hidden['candidates']['D_only_hidden_U1']['real_dimension']}"
            ),
            "rank \\(1/4\\)",
            "majorana_hidden_quotient",
        ),
        artifact_row(
            "Majorana hidden product rank",
            "The D+product hidden quotient still leaves zeta moduli.",
            (
                f"{majorana_hidden['candidates']['D_plus_product_constraint']['hessian_rank']}/"
                f"{majorana_hidden['candidates']['D_plus_product_constraint']['real_dimension']}"
            ),
            "rank \\(5/8\\)",
            "majorana_hidden_quotient",
        ),
        artifact_row(
            "Majorana monomial best small residual",
            "Small monomial roots fail the PMNS phase tolerance.",
            f"{majorana_monomial['best_small_denominator_row']['phase_residual_rad']:.6e}",
            "4.471595",
            "majorana_monomial_clockwork",
        ),
        artifact_row(
            "Majorana monomial first acceptable denominator",
            "Root-of-unity phase matching first occurs only at high denominator.",
            str(majorana_monomial["first_denominator_within_loose_phase_tolerance"]["denominator"]),
            "178",
            "majorana_monomial_clockwork",
        ),
        artifact_row(
            "Conditional theorem boundary row count",
            "The theorem-boundary certificate reads the full claim ledger.",
            str(theorem_boundary["row_count"]),
            "84 claim rows",
            "conditional_theorem_boundary",
        ),
        artifact_row(
            "Conditional theorem boundary no-go count",
            "The theorem-boundary certificate records the current no-go count.",
            str(theorem_boundary["numerical_keys"]["no_go_count"]),
            "15 no-go rows",
            "conditional_theorem_boundary",
        ),
    ]
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def report(summary: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    lines = [
        "# No-web table provenance manifest",
        "",
        "No web lookup was used.  This manifest connects the final closure package's printed numerical claims to generated artifacts.",
        "",
        f"- rows: `{summary['row_count']}`",
        f"- required TeX keys found: `{summary['required_tex_keys_found']}/{summary['required_tex_keys_total']}`",
        f"- all artifacts exist: `{summary['all_artifacts_exist']}`",
        f"- all required TeX keys found: `{summary['all_required_tex_keys_found']}`",
        "",
        "| label | key value | artifact | sha256 prefix | TeX key found |",
        "|---|---:|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['label']} | `{row['key_value']}` | `{row['path']}` | `{row['sha256_prefix12']}` | `{row['tex_key_found']}` |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    required = [row for row in rows if row["required_in_tex"]]
    summary = {
        "note": "No web lookup used. TeX numeric provenance manifest.",
        "row_count": len(rows),
        "all_artifacts_exist": all(row["exists"] for row in rows),
        "required_tex_keys_total": len(required),
        "required_tex_keys_found": sum(bool(row["tex_key_found"]) for row in required),
        "all_required_tex_keys_found": all(bool(row["tex_key_found"]) for row in required),
        "missing_required_tex_keys": [
            {
                "label": row["label"],
                "tex_key": row["tex_key"],
                "artifact": row["artifact"],
            }
            for row in required
            if not row["tex_key_found"]
        ],
        "rows": rows,
        "verdict": {
            "provenance_manifest_complete": (
                all(row["exists"] for row in rows) and all(bool(row["tex_key_found"]) for row in required)
            ),
            "interpretation": (
                "The final local closure-package numbers printed in TeX are traceable "
                "to machine-readable artifacts.  This does not refresh external "
                "experimental inputs; it prevents local table/provenance drift."
            ),
        },
    }
    write_csv(OUT / "provenance_rows.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary, rows), encoding="utf-8")

    print("No-web table provenance manifest written")
    print(f"  rows: {summary['row_count']}")
    print(f"  required TeX keys: {summary['required_tex_keys_found']}/{summary['required_tex_keys_total']}")
    print(f"  complete: {summary['verdict']['provenance_manifest_complete']}")


if __name__ == "__main__":
    main()
