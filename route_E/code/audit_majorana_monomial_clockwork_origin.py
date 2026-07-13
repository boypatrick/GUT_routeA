#!/usr/bin/env python3
"""Audit monomial/clockwork origins for the Majorana contact coefficient.

No web lookup is used.  The question is whether replacing the affine constants
p0,q0 by a monomial hidden sector can predict the PMNS-sensitive complex
coefficient zeta rather than simply re-encoding it.

The conservative audit is algebraic:

* A real-coefficient monomial source can quantize phases as roots of unity.
* A complex monomial coefficient can fit the target phase, but then the phase is
  a spurion.
* A real scale Lambda can fit |zeta|, but then the required tolerance is inherited
  from the PMNS contact-sensitivity window.

We enumerate small monomial/clockwork orders N<=6, and also compute the first
root-of-unity denominator N that can approximate the target phase within the
PMNS loose phase tolerance.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "majorana_monomial_clockwork_origin"
SOURCE_LOCK = ROOT / "output" / "majorana_source_locking_sector" / "summary.json"
SENSITIVITY = ROOT / "output" / "majorana_contact_sensitivity" / "summary.json"
HIDDEN_QUOTIENT = ROOT / "output" / "majorana_hidden_quotient_origin" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cnum(raw: dict[str, float]) -> complex:
    return complex(raw["re"], raw["im"])


def wrap_angle(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def best_root_residual(phase: float, denominator: int) -> dict[str, Any]:
    k = round(denominator * phase / (2.0 * math.pi))
    candidates = []
    for kk in {k - 1, k, k + 1}:
        kk_mod = kk % denominator
        angle = 2.0 * math.pi * kk_mod / denominator
        residual = abs(wrap_angle(phase - angle))
        candidates.append((residual, kk_mod, angle))
    residual, kk, angle = min(candidates, key=lambda item: item[0])
    return {
        "denominator": denominator,
        "k": int(kk),
        "phase_rad": float(angle),
        "phase_residual_rad": float(residual),
    }


def first_denominator_within(phase: float, tolerance: float, max_denominator: int) -> dict[str, Any]:
    best = None
    for denom in range(1, max_denominator + 1):
        row = best_root_residual(phase, denom)
        if best is None or row["phase_residual_rad"] < best["phase_residual_rad"]:
            best = dict(row)
        if row["phase_residual_rad"] <= tolerance:
            out = dict(row)
            out["searched_to"] = max_denominator
            out["best_seen"] = dict(best)
            return out
    assert best is not None
    return {
        "denominator": None,
        "k": None,
        "phase_rad": None,
        "phase_residual_rad": None,
        "searched_to": max_denominator,
        "best_seen": best,
    }


def build() -> dict[str, Any]:
    source = read_json(SOURCE_LOCK)
    sensitivity = read_json(SENSITIVITY)
    hidden = read_json(HIDDEN_QUOTIENT)
    zeta = cnum(source["target"]["zeta"])
    phase = math.atan2(zeta.imag, zeta.real)
    abs_zeta = abs(zeta)
    scale_half = sensitivity["real_scale_loose_interval"]["half_width"]
    phase_half = sensitivity["phase_loose_interval"]["half_width"]
    rows = []
    for denom in range(1, 7):
        row = best_root_residual(phase, denom)
        row["passes_loose_phase"] = row["phase_residual_rad"] <= phase_half
        row["residual_over_tolerance"] = row["phase_residual_rad"] / phase_half
        # If |zeta| is generated as Lambda^power in units of a reference scale,
        # fractional source-scale stability must be at least this good.
        row["continuous_scale_input_required"] = True
        row["scale_relative_tolerance_if_zeta_direct"] = scale_half
        row["scale_relative_tolerance_if_power_denominator"] = scale_half / denom
        rows.append(row)
    best_small = min(rows, key=lambda row: row["phase_residual_rad"])
    first = first_denominator_within(phase, phase_half, 200_000)
    phase_spurion_needed_small = not any(row["passes_loose_phase"] for row in rows)
    continuous_scale_tuning = scale_half
    summary = {
        "note": "No web lookup used. Monomial/clockwork origin audit for Majorana zeta.",
        "input_manifest": [
            {
                "label": "majorana_source_locking_sector",
                "path": str(SOURCE_LOCK.relative_to(ROOT)),
                "sha256": sha256(SOURCE_LOCK),
            },
            {
                "label": "majorana_contact_sensitivity",
                "path": str(SENSITIVITY.relative_to(ROOT)),
                "sha256": sha256(SENSITIVITY),
            },
            {
                "label": "majorana_hidden_quotient_origin",
                "path": str(HIDDEN_QUOTIENT.relative_to(ROOT)),
                "sha256": sha256(HIDDEN_QUOTIENT),
            },
        ],
        "target": {
            "abs_zeta": float(abs_zeta),
            "arg_zeta_rad": float(phase),
            "loose_scale_half_width": float(scale_half),
            "loose_phase_half_width_rad": float(phase_half),
        },
        "small_denominator_rows": rows,
        "best_small_denominator_row": best_small,
        "first_denominator_within_loose_phase_tolerance": first,
        "hidden_quotient_baseline": {
            "D_only_rank": hidden["candidates"]["D_only_hidden_U1"]["hessian_rank"],
            "D_only_dimension": hidden["candidates"]["D_only_hidden_U1"]["real_dimension"],
            "D_plus_product_rank": hidden["candidates"]["D_plus_product_constraint"]["hessian_rank"],
            "D_plus_product_dimension": hidden["candidates"]["D_plus_product_constraint"]["real_dimension"],
        },
        "verdict": {
            "small_monomial_clockwork_predicts_phase": not phase_spurion_needed_small,
            "small_monomial_clockwork_predicts_scale": False,
            "complex_spurion_needed_for_small_orders": phase_spurion_needed_small,
            "continuous_scale_tuning_required": True,
            "scale_relative_tuning_if_direct": float(continuous_scale_tuning),
            "interpretation": (
                "Small monomial/clockwork orders do not predict the target contact "
                "phase within the PMNS tolerance.  A real-coefficient root-of-unity "
                "mechanism would need a much larger denominator, while any complex "
                "coefficient simply reintroduces the phase as a spurion.  The magnitude "
                "still requires a continuous hidden scale tuned to the contact-sensitivity "
                "window.  Thus the monomial clockwork route is not a first-principles "
                "origin for zeta in the small-order regime."
            ),
        },
    }
    return summary


def report(summary: dict[str, Any]) -> str:
    target = summary["target"]
    best = summary["best_small_denominator_row"]
    first = summary["first_denominator_within_loose_phase_tolerance"]
    lines = [
        "# Majorana Monomial/Clockwork Origin Audit",
        "",
        "No web lookup was used.  This audit tests whether small monomial hidden sectors can predict the PMNS-sensitive contact coefficient.",
        "",
        "## Target and tolerance",
        "",
        f"- |zeta|: `{target['abs_zeta']:.6e}`",
        f"- arg(zeta): `{target['arg_zeta_rad']:.6e}` rad",
        f"- loose scale half-width: `{target['loose_scale_half_width']:.6e}`",
        f"- loose phase half-width: `{target['loose_phase_half_width_rad']:.6e}` rad",
        "",
        "## Small root-of-unity scan",
        "",
        "| N | k | phase residual | residual/tolerance | passes |",
        "|---:|---:|---:|---:|---|",
    ]
    for row in summary["small_denominator_rows"]:
        lines.append(
            f"| {row['denominator']} | {row['k']} | `{row['phase_residual_rad']:.6e}` | "
            f"`{row['residual_over_tolerance']:.6e}` | `{row['passes_loose_phase']}` |"
        )
    lines += [
        "",
        "## Best small order",
        "",
        f"- best N<=6: `{best['denominator']}`",
        f"- best k: `{best['k']}`",
        f"- residual: `{best['phase_residual_rad']:.6e}` rad",
        f"- residual/tolerance: `{best['residual_over_tolerance']:.6e}`",
        "",
        "## First denominator within tolerance",
        "",
        f"- denominator: `{first['denominator']}`",
        f"- k: `{first['k']}`",
        f"- residual: `{first['phase_residual_rad']:.6e}` rad",
        "",
        "## Magnitude",
        "",
        f"- direct zeta scale relative tolerance: `{target['loose_scale_half_width']:.6e}`",
        "- a real scale Lambda can fit this, but Lambda is then a continuous input rather than a prediction.",
        "",
        "## Verdict",
        "",
        summary["verdict"]["interpretation"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")
    best = summary["best_small_denominator_row"]
    first = summary["first_denominator_within_loose_phase_tolerance"]
    print("Majorana monomial/clockwork origin audit")
    print(f"  best N<=6: {best['denominator']} residual {best['phase_residual_rad']:.6e}")
    print(f"  first N within phase tolerance: {first['denominator']}")
    print(f"  predicts small-order phase: {summary['verdict']['small_monomial_clockwork_predicts_phase']}")


if __name__ == "__main__":
    main()
