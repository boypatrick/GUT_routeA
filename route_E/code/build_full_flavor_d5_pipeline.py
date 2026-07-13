#!/usr/bin/env python3
"""Build the full-flavor plus d=5 proton-decay closure ledger.

No web lookup is used.  This script is intentionally a closure audit rather
than another toy scan: it reads the exact flavor/seesaw card, the
operator-level two-kernel flavor->triplet scan, and the Knu dressed d=5
pipeline, then asks whether the present local evidence is publication-complete.

Mathematical normalization used below:

    C5L^{abcd} = sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl}
                 U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld},

    C6_ch = D_ch C5_ch,
    Gamma_ch = K_ch |C6_ch|^2,
    margin_ch = tau_ch / tau_target.

If a row has margin mu < 1, the required common amplitude suppression is
mu^{-1/2}.  That conversion is the key bridge between flavor/Wilson output
and the proton lifetime target.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "full_flavor_d5_pipeline"

LOCKED_SUMMARY = ROOT / "output" / "locked_link_full_flavor_d5_card" / "summary.json"
LOCKED_CARD = ROOT / "output" / "locked_link_full_flavor_d5_card" / "locked_link_full_flavor_d5_card.json"
TWO_KERNEL = ROOT / "output" / "two_kernel_flavor_then_d5" / "summary.json"
FLAVOR_AUDIT = ROOT / "output" / "flavor_fit" / "flavor_observable_audit.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"


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


def suppression_needed(margin: float) -> float:
    if margin <= 0.0:
        return math.inf
    return max(1.0 / math.sqrt(margin), 1.0)


def score_row(
    *,
    label: str,
    source: str,
    ckm_score: float,
    mass_score: float | None,
    future_margin: float,
    strict_ckm: float,
    loose_ckm: float,
    loose_mass: float,
    details: dict[str, Any],
) -> dict[str, Any]:
    return {
        "label": label,
        "source": source,
        "ckm_score": ckm_score,
        "mass_score": mass_score,
        "future_margin_1e35": future_margin,
        "ckm_over_strict": ckm_score / strict_ckm if strict_ckm > 0 else math.inf,
        "ckm_over_loose": ckm_score / loose_ckm if loose_ckm > 0 else math.inf,
        "mass_over_loose": (mass_score / loose_mass if mass_score is not None and loose_mass > 0 else None),
        "d5_amplitude_suppression_needed": suppression_needed(future_margin),
        "passes_strict_ckm": ckm_score < strict_ckm,
        "passes_loose_ckm": ckm_score < loose_ckm,
        "passes_loose_mass": (mass_score is not None and mass_score <= loose_mass),
        "passes_future_d5": future_margin >= 1.0,
        "details": details,
    }


def build() -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    locked = read_json(LOCKED_SUMMARY)
    locked_card = read_json(LOCKED_CARD)
    two = read_json(TWO_KERNEL)
    flavor = read_json(FLAVOR_AUDIT)
    knu = read_json(KNU_TARGET)

    two_cfg = two["config"]
    strict_ckm = float(two_cfg["strict_ckm"])
    loose_ckm = float(two_cfg["loose_ckm"])
    loose_mass = float(two_cfg["loose_mass"])

    exact_verdict = locked["verdict"]
    knu_verdict = knu["verdict"]
    two_verdict = two["verdict"]

    rows = [
        score_row(
            label="exact_card_current",
            source=str(FLAVOR_AUDIT.relative_to(ROOT)),
            ckm_score=float(flavor["ckm_magnitude_log_score"]),
            mass_score=None,
            future_margin=float(knu_verdict["future_margin"]),
            strict_ckm=strict_ckm,
            loose_ckm=loose_ckm,
            loose_mass=loose_mass,
            details={
                "Vus": flavor["ckm_current_observables"]["abs_Vus"],
                "Vcb": flavor["ckm_current_observables"]["abs_Vcb"],
                "Vub": flavor["ckm_current_observables"]["abs_Vub"],
                "J": flavor["ckm_current_observables"]["jarlskog"],
            },
        )
    ]
    for key in ["best_ckm", "best_mass", "best_knu", "best_balanced"]:
        item = two_verdict[key]
        rows.append(
            score_row(
                label=f"operator_{key}",
                source=str(TWO_KERNEL.relative_to(ROOT)),
                ckm_score=float(item["ckm_score"]),
                mass_score=float(item["mass_score"]),
                future_margin=float(item["future_margin"]),
                strict_ckm=strict_ckm,
                loose_ckm=loose_ckm,
                loose_mass=loose_mass,
                details={
                    "candidate_label": item["label"],
                    "profile_hint": item["profile_hint"],
                    "best_triplet_profile": item["best_triplet_profile"],
                    "Vus": item["Vus"],
                    "Vcb": item["Vcb"],
                    "Vub": item["Vub"],
                    "J": item["J_abs"],
                    "seesaw_residual": item["seesaw_residual"],
                },
            )
        )

    # Deduplicate labels that point to the same optimizer row while preserving
    # their role labels in the CSV/report.
    best_operator_future_margin = max(row["future_margin_1e35"] for row in rows if row["label"].startswith("operator_"))
    best_operator_ckm = min(row["ckm_score"] for row in rows if row["label"].startswith("operator_"))

    component_rows = [
        {
            "component": "exact_reproducibility_card",
            "status": "PASS" if exact_verdict["exact_inputs_available"] and exact_verdict["pmns_reconstructed"] else "FAIL",
            "metric": "exact_inputs={}, pmns={}".format(
                exact_verdict["exact_inputs_available"],
                exact_verdict["pmns_reconstructed"],
            ),
            "evidence": str(LOCKED_CARD.relative_to(ROOT)),
        },
        {
            "component": "baseline_full_flavor",
            "status": "OPEN",
            "metric": "CKM_log_score={:.6e}".format(float(flavor["ckm_magnitude_log_score"])),
            "evidence": str(FLAVOR_AUDIT.relative_to(ROOT)),
        },
        {
            "component": "operator_level_two_kernel_closure",
            "status": "NO_CLOSURE",
            "metric": "strict={}, loose={}, best_ckm={:.6e}, best_margin={:.6e}".format(
                int(two_verdict["strict_closure_count"]),
                int(two_verdict["loose_closure_count"]),
                best_operator_ckm,
                best_operator_future_margin,
            ),
            "evidence": str(TWO_KERNEL.relative_to(ROOT)),
        },
        {
            "component": "current_Knu_bound",
            "status": "PASS_CONDITIONAL" if exact_verdict["d5_current_bound_passes"] else "FAIL",
            "metric": "current_margin={:.6e}".format(float(knu_verdict["current_margin"])),
            "evidence": str(KNU_TARGET.relative_to(ROOT)),
        },
        {
            "component": "future_1e35_Knu_stress",
            "status": "OPEN",
            "metric": "future_margin={:.6e}, suppression_needed={:.6e}".format(
                float(knu_verdict["future_margin"]),
                float(knu_verdict["future_amplitude_suppression_needed"]),
            ),
            "evidence": str(KNU_TARGET.relative_to(ROOT)),
        },
    ]

    closure = {
        "exact_inputs_available": bool(exact_verdict["exact_inputs_available"]),
        "pmns_reconstructed": bool(exact_verdict["pmns_reconstructed"]),
        "baseline_ckm_completed": bool(exact_verdict["ckm_fit_completed"]),
        "operator_strict_closures": int(two_verdict["strict_closure_count"]),
        "operator_loose_closures": int(two_verdict["loose_closure_count"]),
        "d5_current_bound_passes": bool(exact_verdict["d5_current_bound_passes"]),
        "d5_future_1e35_passes": bool(exact_verdict["d5_future_1e35_passes"]),
        "publication_level_complete": False,
    }
    closure["publication_level_complete"] = bool(
        closure["exact_inputs_available"]
        and closure["pmns_reconstructed"]
        and closure["baseline_ckm_completed"]
        and closure["operator_strict_closures"] > 0
        and closure["d5_current_bound_passes"]
        and closure["d5_future_1e35_passes"]
    )

    payload = {
        "note": "No web lookup used. Full flavor plus dressed d=5 closure pipeline.",
        "formulae": {
            "C5L": "sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl} U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld}",
            "C5R": "sum_AB (Y_UE^A)_{ij} W_AB (Y_UD^B)_{kl} U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}",
            "dressing": "C6_ch = D_ch C5_ch",
            "width": "Gamma_ch = K_ch |C6_ch|^2",
            "suppression_rule": "required common amplitude suppression = margin^{-1/2} for margin<1",
        },
        "input_manifest": manifest([LOCKED_SUMMARY, LOCKED_CARD, TWO_KERNEL, FLAVOR_AUDIT, KNU_TARGET]),
        "scoreboard": rows,
        "component_status": component_rows,
        "closure": closure,
        "numerical_verdict": {
            "baseline_ckm_score": float(flavor["ckm_magnitude_log_score"]),
            "best_operator_ckm_score": best_operator_ckm,
            "best_operator_future_margin": best_operator_future_margin,
            "best_operator_suppression_needed": suppression_needed(best_operator_future_margin),
            "knu_future_margin": float(knu_verdict["future_margin"]),
            "knu_future_suppression_needed": float(knu_verdict["future_amplitude_suppression_needed"]),
            "interpretation": (
                "The reproducible card exists and the current Knu bound passes, "
                "but neither the baseline card nor the operator-level two-kernel "
                "scan closes full flavor and future-stress d=5 proton safety. "
                "The best operator-level row improves CKM to about the strict "
                "threshold but still has future Knu margin below one."
            ),
        },
    }
    return payload, rows, component_rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    verdict = payload["numerical_verdict"]
    closure = payload["closure"]
    lines = [
        "# Full flavor plus dressed d=5 closure pipeline",
        "",
        "No web lookup was used.",
        "",
        "## Mathematical bridge",
        "",
        "The channel Wilson coefficient is built in the same mass basis as the",
        "flavor fit:",
        "",
        "```text",
        "C5L^{abcd} = sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl}",
        "             U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld}",
        "C6_ch = D_ch C5_ch",
        "Gamma_ch = K_ch |C6_ch|^2",
        "```",
        "",
        "Thus a lifetime margin `mu=tau/tau_target` below one requires a common",
        "amplitude suppression `mu^{-1/2}`.",
        "",
        "## Component Status",
        "",
        "| component | status | metric |",
        "|---|---|---:|",
    ]
    for row in payload["component_status"]:
        lines.append(f"| `{row['component']}` | `{row['status']}` | `{row['metric']}` |")
    lines += [
        "",
        "## Candidate Scoreboard",
        "",
        "| label | CKM score | mass score | future margin | suppression needed |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["scoreboard"]:
        mass = "NA" if row["mass_score"] is None else f"{row['mass_score']:.6e}"
        lines.append(
            f"| `{row['label']}` | {row['ckm_score']:.6e} | {mass} | "
            f"{row['future_margin_1e35']:.6e} | "
            f"{row['d5_amplitude_suppression_needed']:.6e} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        f"Publication-level closure: `{closure['publication_level_complete']}`.",
        f"Best operator CKM score: `{verdict['best_operator_ckm_score']:.6e}`.",
        f"Best operator future Knu margin: `{verdict['best_operator_future_margin']:.6e}`.",
        (
            "Required suppression at best operator point: "
            f"`{verdict['best_operator_suppression_needed']:.6e}`."
        ),
        "",
        verdict["interpretation"],
        "",
        "Machine-readable outputs:",
        "",
        "- `output/full_flavor_d5_pipeline/summary.json`",
        "- `output/full_flavor_d5_pipeline/candidate_scoreboard.csv`",
        "- `output/full_flavor_d5_pipeline/component_status.csv`",
        "- `output/full_flavor_d5_pipeline/input_manifest.csv`",
    ]
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload, score_rows, component_rows = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        OUT / "candidate_scoreboard.csv",
        score_rows,
        [
            "label",
            "source",
            "ckm_score",
            "mass_score",
            "future_margin_1e35",
            "ckm_over_strict",
            "ckm_over_loose",
            "mass_over_loose",
            "d5_amplitude_suppression_needed",
            "passes_strict_ckm",
            "passes_loose_ckm",
            "passes_loose_mass",
            "passes_future_d5",
        ],
    )
    write_csv(OUT / "component_status.csv", component_rows, ["component", "status", "metric", "evidence"])
    write_csv(
        OUT / "input_manifest.csv",
        payload["input_manifest"],
        ["path", "exists", "size_bytes", "sha256"],
    )
    write_report(payload)
    verdict = payload["numerical_verdict"]
    closure = payload["closure"]
    print("Full flavor + dressed d=5 closure pipeline")
    print(f"  publication complete: {closure['publication_level_complete']}")
    print(f"  operator strict closures: {closure['operator_strict_closures']}")
    print(f"  operator loose closures: {closure['operator_loose_closures']}")
    print(f"  best operator CKM score: {verdict['best_operator_ckm_score']:.6e}")
    print(f"  best operator future margin: {verdict['best_operator_future_margin']:.6e}")
    print(f"  best operator suppression needed: {verdict['best_operator_suppression_needed']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
