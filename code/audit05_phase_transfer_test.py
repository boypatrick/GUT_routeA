#!/usr/bin/env python3
"""Audit 0.5 phase-transfer diagnostic.

This cheap diagnostic tests whether the benchmark Majorana contact phase could
be inherited from the visible Veronese spurion invariants already recorded in
Audit 0.  It compares arg(zeta) against a deliberately small table:

    arg I, 2 arg I, arg J, arg I + arg J.

For real single-coefficient monomials, the sign of the coefficient makes the
phase meaningful modulo pi.  The exception kept here is two_arg_I, viewed as
lambda=c I followed by zeta=lambda^2, where c^2>0 and the phase is tested
modulo 2 pi.

The script also solves the two-real-coefficient diagnostic

    zeta = c1 I + c2 J,  c1,c2 real,

which generally fits exactly and is therefore a naturalness diagnostic rather
than a prediction.  It is a falsification-oriented side test, not a
hidden-sector derivation and not a global flavor fit.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit05"

SOURCES = {
    "audit05_builder": ROOT / "code" / "audit05_phase_transfer_test.py",
    "audit0_card": ROOT / "output" / "audit0" / "invariant_card.json",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    encoded = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def manifest() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label, path in SOURCES.items():
        rows.append(
            {
                "label": label,
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path),
            }
        )
    return rows


def complex_from_pair(pair: dict[str, float]) -> complex:
    return complex(float(pair["re"]), float(pair["im"]))


def principal_phase(x: float) -> float:
    return float(x % (2.0 * math.pi))


def signed_angular_delta(candidate: float, target: float) -> float:
    return float(math.atan2(math.sin(candidate - target), math.cos(candidate - target)))


def signed_delta_mod_pi(candidate: float, target: float) -> float:
    delta = signed_angular_delta(candidate, target)
    if delta > math.pi / 2.0:
        delta -= math.pi
    elif delta < -math.pi / 2.0:
        delta += math.pi
    return float(delta)


def phase_row(
    label: str,
    phase: float,
    target: float,
    loose: float,
    tight: float,
    distance_mode: str,
) -> dict[str, Any]:
    if distance_mode == "mod_pi":
        delta = signed_delta_mod_pi(phase, target)
    elif distance_mode == "mod_2pi":
        delta = signed_angular_delta(phase, target)
    else:
        raise ValueError(f"unknown distance mode: {distance_mode}")
    return {
        "candidate": label,
        "distance_mode": distance_mode,
        "phase_rad": principal_phase(phase),
        "target_arg_zeta_rad": principal_phase(target),
        "signed_delta_rad": delta,
        "abs_delta_rad": abs(delta),
        "hits_loose_window": abs(delta) <= loose,
        "hits_tight_window": abs(delta) <= tight,
    }


def solve_real_two_term(i_inv: complex, j_inv: complex, zeta: complex) -> dict[str, Any]:
    det = i_inv.real * j_inv.imag - j_inv.real * i_inv.imag
    if abs(det) < 1.0e-15:
        return {
            "status": "singular_real_system",
            "determinant": det,
            "interpretation": "I and J are real-linearly dependent; no stable two-real-coefficient diagnostic is available.",
        }
    c1 = (zeta.real * j_inv.imag - j_inv.real * zeta.imag) / det
    c2 = (i_inv.real * zeta.imag - zeta.real * i_inv.imag) / det
    term_i = c1 * i_inv
    term_j = c2 * j_inv
    residual = term_i + term_j - zeta
    zabs = abs(zeta)
    max_term = max(abs(term_i), abs(term_j))
    sum_terms = abs(term_i) + abs(term_j)
    return {
        "status": "solved",
        "determinant": det,
        "c1_real": c1,
        "c2_real": c2,
        "term_c1I": {"re": term_i.real, "im": term_i.imag, "abs": abs(term_i)},
        "term_c2J": {"re": term_j.real, "im": term_j.imag, "abs": abs(term_j)},
        "residual_abs": abs(residual),
        "diagnostic_definitions": {
            "naturalness_index_max_term_over_abs_zeta": "max(|c1 I|, |c2 J|)/|zeta|",
            "cancellation_index_sum_terms_over_abs_zeta": "(|c1 I|+|c2 J|)/|zeta|",
            "residual_abs": "|c1 I + c2 J - zeta|",
            "coefficient_normalization_caveat": "The raw sizes of c1 and c2 are not invariant under rescaling I or J. Interpret c1,c2 only together with the quartic invariant convention, basis, and term magnitudes recorded in this card.",
        },
        "naturalness_index_max_term_over_abs_zeta": max_term / zabs,
        "cancellation_index_sum_terms_over_abs_zeta": sum_terms / zabs,
        "interpretation": (
            "A two-real-coefficient combination can fit the phase by construction. "
            "The cancellation index is mild at O(1) for this card, so the class is "
            "not excluded as a possibility, but it is an unconstrained fit rather "
            "than a phase prediction."
        ),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    audit0 = read_json(SOURCES["audit0_card"])
    mv_inv = audit0["majorana_invariant_anchors"]["M_V_quartic_invariants"]
    i_inv = complex_from_pair(mv_inv["I"])
    j_inv = complex_from_pair(mv_inv["J"])
    zeta = complex_from_pair(audit0["majorana_invariant_anchors"]["zeta_tar"])
    target = float(audit0["majorana_invariant_anchors"]["zeta_arg_rad"])
    loose = float(audit0["legacy_local_conventions_index"]["contact_sensitivity"]["phase_loose_half_width_rad"])
    tight = float(audit0["legacy_local_conventions_index"]["contact_sensitivity"]["phase_tight_half_width_rad"])

    phase_i = math.atan2(i_inv.imag, i_inv.real)
    phase_j = math.atan2(j_inv.imag, j_inv.real)
    candidates = [
        phase_row("arg_I", phase_i, target, loose, tight, "mod_pi"),
        phase_row("two_arg_I", 2.0 * phase_i, target, loose, tight, "mod_2pi"),
        phase_row("arg_J", phase_j, target, loose, tight, "mod_pi"),
        phase_row("arg_I_plus_arg_J", phase_i + phase_j, target, loose, tight, "mod_pi"),
    ]
    two_term_fit = solve_real_two_term(i_inv, j_inv, zeta)

    loose_hits = [r for r in candidates if r["hits_loose_window"]]
    tight_hits = [r for r in candidates if r["hits_tight_window"]]
    closest = min(candidates, key=lambda r: r["abs_delta_rad"])

    card: dict[str, Any] = {
        "audit": "audit05_phase_transfer_test",
        "status": "finite phase-transfer diagnostic only; not a hidden-sector derivation",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "audit0_card_sha256": audit0["card_sha256"],
        "inputs": {
            "I": mv_inv["I"],
            "J": mv_inv["J"],
            "zeta": audit0["majorana_invariant_anchors"]["zeta_tar"],
            "arg_zeta_rad": target,
            "phase_loose_half_width_rad": loose,
            "phase_tight_half_width_rad": tight,
            "scheme_labels": audit0["scheme_labels"],
            "quartic_invariant_scheme": mv_inv["quartic_coefficients_convention"],
            "quartic_basis": mv_inv["basis"],
            "normalization_caveat": "I and J are binary-quartic invariants in the recorded coefficient convention and normalized CP1/O(2) basis. Raw fitted c1,c2 values are convention-dependent.",
        },
        "diagnostic_definitions": {
            "mod_pi_rows": "For arg I, arg J, and arg I+arg J, the signed phase distance is minimized modulo pi because a real coefficient may be positive or negative.",
            "mod_2pi_rows": "For 2 arg I, the signed phase distance is measured modulo 2 pi because the diagnostic is lambda=cI followed by zeta=lambda^2 with c^2>0.",
            "single_real_parameter_no_hit": "No row in the finite table lands inside the selected phase window under the row-specific distance convention.",
            "two_real_parameter_fit_boundary": "The two-term real fit zeta=c1 I+c2 J is generally solvable and is reported only with cost indices.",
        },
        "candidate_table": candidates,
        "real_two_term_fit_zeta_equals_c1I_plus_c2J": two_term_fit,
        "closest_candidate": closest,
        "loose_hit_count": len(loose_hits),
        "tight_hit_count": len(tight_hits),
        "verdict": {
            "simple_visible_spurion_phase_transfer_passes_loose": bool(loose_hits),
            "simple_visible_spurion_phase_transfer_passes_tight": bool(tight_hits),
            "interpretation": (
                "No candidate in {arg I, 2 arg I, arg J, arg I+arg J} lands inside "
                "the selected window under the stated mod-pi/mod-2pi distance rules; "
                "this falsifies only the simplest zero-parameter and single-real-parameter "
                "visible-spurion phase-transfer table for the current local benchmark."
                if not loose_hits
                else "At least one finite phase-transfer candidate lands in the loose window; this is a mechanism hint, not a derivation."
            ),
        },
        "publication_boundary": {
            "not_claimed": [
                "derivation of zeta",
                "global hidden-sector phase model",
                "flavor fit",
                "proton or threshold closure",
            ],
            "claim": "Audit 0.5 only tests a finite visible-spurion phase-transfer table against the local benchmark.",
        },
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "phase_transfer_test.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        "# Audit 0.5 Phase-Transfer Diagnostic",
        "",
        "This diagnostic compares `arg(zeta)` against the finite visible-spurion table `{arg I, 2 arg I, arg J, arg I + arg J}`.",
        "",
        "Real single-coefficient monomials are tested modulo pi because the real coefficient may have either sign. The `2 arg I` row is tested modulo 2 pi because it is interpreted as `lambda = c I`, followed by `zeta = lambda^2` with `c^2 > 0`.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- audit0 sha256: `{audit0['card_sha256']}`",
        f"- source artifacts: {sum(1 for item in card['source_manifest'] if item['exists'])}/{len(card['source_manifest'])} present",
        "",
        "## Result",
        "",
        f"- loose hit count: {len(loose_hits)}",
        f"- tight hit count: {len(tight_hits)}",
        f"- closest candidate: `{closest['candidate']}`",
        f"- closest absolute phase error: `{closest['abs_delta_rad']:.12e}` rad",
        f"- loose window: `{loose:.12e}` rad",
        f"- tight window: `{tight:.12e}` rad",
        "",
        "## Two-Term Diagnostic",
        "",
        "The real two-term ansatz `zeta = c1 I + c2 J` is solved exactly when the real 2x2 system is nonsingular. This is a fit-cost diagnostic, not a prediction.",
        "",
        "Raw `c1,c2` sizes are convention-dependent because they scale with the chosen normalization of `I,J`; interpret them together with the recorded basis, invariant convention, and term magnitudes.",
        "",
        f"- status: `{two_term_fit['status']}`",
    ]
    if two_term_fit["status"] == "solved":
        md_lines += [
            f"- c1: `{two_term_fit['c1_real']:.12e}`",
            f"- c2: `{two_term_fit['c2_real']:.12e}`",
            f"- residual abs: `{two_term_fit['residual_abs']:.12e}`",
            f"- naturalness index max(|c1 I|, |c2 J|)/|zeta|: `{two_term_fit['naturalness_index_max_term_over_abs_zeta']:.12e}`",
            f"- cancellation index (|c1 I|+|c2 J|)/|zeta|: `{two_term_fit['cancellation_index_sum_terms_over_abs_zeta']:.12e}`",
        ]
    md_lines += [
        "",
        "## Candidate Table",
        "",
        "| candidate | distance | phase rad | signed delta rad | loose | tight |",
        "|---|---|---:|---:|---|---|",
    ]
    for item in candidates:
        md_lines.append(
            f"| `{item['candidate']}` | `{item['distance_mode']}` | {item['phase_rad']:.12e} | {item['signed_delta_rad']:.12e} | {item['hits_loose_window']} | {item['hits_tight_window']} |"
        )
    md_lines += [
        "",
        "## Boundary",
        "",
        card["publication_boundary"]["claim"],
        "",
    ]

    md_path = OUT / "phase_transfer_test.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {md_path.relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")
    print(f"loose_hit_count = {len(loose_hits)}")
    print(f"closest = {closest['candidate']} ({closest['abs_delta_rad']:.12e} rad)")


if __name__ == "__main__":
    main()
