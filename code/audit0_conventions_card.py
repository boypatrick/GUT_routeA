#!/usr/bin/env python3
"""Build Audit 0 convention/invariant card for deferred companion audits.

Audit 0 is the single source of truth for follow-up audits.  It does not fit
flavor, compute thresholds, or evaluate proton lifetimes.  It records the
basis, invariant anchors, precision policy, delta-Z replay convention, and
audit dependency graph used by the deferred package.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit0"

SOURCES = {
    "audit0_builder": ROOT / "code" / "audit0_conventions_card.py",
    "routeB_hidden_zeta": ROOT
    / "output"
    / "hidden_zeta_origin"
    / "routeB_hidden_zeta_verification.json",
    "source_majorana_texture_rank": ROOT
    / "output"
    / "source_majorana_texture_rank"
    / "summary.json",
    "d5_half_spin": ROOT / "output" / "d5_half_spin" / "d5_half_spin_summary.json",
    "majorana_contact_sensitivity": ROOT
    / "output"
    / "majorana_contact_sensitivity"
    / "summary.json",
    "no_web_input_convention_ledger": ROOT
    / "output"
    / "no_web_input_convention_ledger"
    / "summary.json",
    "flavor_benchmark_card": ROOT
    / "output"
    / "flavor_benchmark"
    / "flavor_benchmark_card.json",
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
    rows = []
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


def cpair(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def from_pair(item: dict[str, float]) -> complex:
    return complex(float(item["re"]), float(item["im"]))


def cmatrix(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[from_pair(cell) for cell in row] for row in raw], dtype=complex)


def continued_fraction_convergents(x: float, n_terms: int = 8) -> list[dict[str, float]]:
    terms: list[int] = []
    y = x
    for _ in range(n_terms):
        a = math.floor(y)
        terms.append(a)
        frac = y - a
        if abs(frac) < 1.0e-18:
            break
        y = 1.0 / frac

    p0, p1 = 0, 1
    q0, q1 = 1, 0
    rows: list[dict[str, float]] = []
    for a in terms:
        p = a * p1 + p0
        q = a * q1 + q0
        phase_error = abs(2.0 * math.pi * p / q - 2.0 * math.pi * x)
        rows.append({"p": p, "q": q, "phase_error_rad": phase_error})
        p0, p1 = p1, p
        q0, q1 = q1, q
    return rows


def veronese_quartic_invariants(mv: np.ndarray) -> dict[str, Any]:
    """Map a Veronese symmetric matrix to binary-quartic invariants.

    Convention:
      V = Sym^2 C^2 has basis (x^2, sqrt(2)xy, y^2).
      A Veronese matrix M with M_11 = 2 M_02 maps to

        a x^4 + 4 b x^3 y + 6 c x^2 y^2 + 4 d x y^3 + e y^4

      via a=M00, b=M01/sqrt(2), c=M02, d=M12/sqrt(2), e=M22.
      The classical binary quartic invariants are

        I = a e - 4 b d + 3 c^2,
        J = a c e + 2 b c d - a d^2 - b^2 e - c^3.

    These are internal convention anchors, not publication-level modular data.
    """

    rt2 = math.sqrt(2.0)
    a = mv[0, 0]
    b = mv[0, 1] / rt2
    c = mv[0, 2]
    d = mv[1, 2] / rt2
    e = mv[2, 2]
    i_inv = a * e - 4.0 * b * d + 3.0 * c * c
    j_inv = a * c * e + 2.0 * b * c * d - a * d * d - b * b * e - c**3
    discr = i_inv**3 - 27.0 * j_inv**2
    j_like = None if abs(discr) == 0 else i_inv**3 / discr
    constraint = mv[1, 1] - 2.0 * mv[0, 2]
    return {
        "quartic_coefficients_convention": "a x^4 + 4 b x^3 y + 6 c x^2 y^2 + 4 d x y^3 + e y^4",
        "basis": "(x^2, sqrt(2) x y, y^2)",
        "coefficients": {
            "a": cpair(a),
            "b": cpair(b),
            "c": cpair(c),
            "d": cpair(d),
            "e": cpair(e),
        },
        "I": cpair(i_inv),
        "J": cpair(j_inv),
        "discriminant_like_I3_minus_27J2": cpair(discr),
        "j_like_I3_over_discriminant": None if j_like is None else cpair(j_like),
        "veronese_constraint_M11_minus_2M02": cpair(constraint),
        "veronese_constraint_abs": float(abs(constraint)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    zeta_payload = read_json(SOURCES["routeB_hidden_zeta"])
    source_payload = read_json(SOURCES["source_majorana_texture_rank"])
    d5_payload = read_json(SOURCES["d5_half_spin"])
    sensitivity_payload = read_json(SOURCES["majorana_contact_sensitivity"])
    no_web_payload = read_json(SOURCES["no_web_input_convention_ledger"])
    flavor_payload = read_json(SOURCES["flavor_benchmark_card"])

    zeta = complex(*zeta_payload["zeta"])
    sqrt_zeta = complex(*zeta_payload["sqrt_zeta"])
    k_tr = (1.0 / math.sqrt(3.0)) * np.array(
        [[0, 0, -1], [0, 1, 0], [-1, 0, 0]],
        dtype=complex,
    )

    decomp = source_payload["decomposition"]
    mv = cmatrix(decomp["veronese_projection"])
    mr_normalized = cmatrix(decomp["MR_normalized"])

    x_phase = float(np.angle(zeta) / (2.0 * math.pi))
    cf_rows = continued_fraction_convergents(x_phase, n_terms=8)
    relevant_cf = [row for row in cf_rows if row["q"] in {178, 733, 911, 1644}]

    card: dict[str, Any] = {
        "audit": "audit0_conventions_card",
        "status": "foundation card for deferred companion audits; not a flavor/proton/threshold closure",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "scheme_labels": {
            "matrix_basis": "CP1/O(2) normalized family basis psi=(x^2, sqrt(2)xy, y^2)",
            "quartic_invariant_scheme": "binary quartic coefficients (a,b,c,d,e) with combinatorial weights 1,4,6,4,1",
            "rge_scheme": "UNSET_FOR_PUBLICATION: follow-up audits must explicitly choose DRbar/MSbar and scale thresholds",
            "tan_beta_convention": "UNSET_FOR_PUBLICATION: required by Audit 1 and Audit 3",
            "M_star_definition": "benchmark/source-sector scale; must be tied to v_R or source UV in Audit 4a",
            "precision_policy": "long decimals are reproducibility anchors; physical claims should be limited to a few significant figures until UV replay",
        },
        "routeA_invariant_anchors": {
            "d5_half_spin_counts": d5_payload["field_counts"],
            "d5_expected_counts_match": d5_payload["field_counts_match"],
            "hypercharge_trace": {
                "TrY2": d5_payload["trace_Y2"],
                "TrT3L2": d5_payload["trace_T3L2"],
                "kY": d5_payload["kY"],
            },
            "k_tr": {
                "matrix": [[cpair(z) for z in row] for row in k_tr],
                "frobenius_norm": float(np.linalg.norm(k_tr)),
                "symmetric_error": float(np.linalg.norm(k_tr - k_tr.T)),
                "inverse_check_norm_K_times_3K_minus_I": float(np.linalg.norm(k_tr @ (3.0 * k_tr) - np.eye(3))),
                "eigenvalues": [cpair(z) for z in np.linalg.eigvals(k_tr)],
                "meaning": "unique second-transvectant/Killing-form line up to basis, sign, and normalization",
            },
        },
        "majorana_invariant_anchors": {
            "zeta_tar": cpair(zeta),
            "zeta_abs": float(abs(zeta)),
            "zeta_arg_rad": float(np.angle(zeta)),
            "sqrt_zeta": cpair(sqrt_zeta),
            "sqrt_zeta_abs": float(abs(sqrt_zeta)),
            "lambda_loop_prefactor_abs_lambda_sq_over_16pi2": float(abs(sqrt_zeta) ** 2 / (16.0 * math.pi**2)),
            "instanton_action_minus_log_abs_zeta": float(-math.log(abs(zeta))),
            "source_rank_contact_coefficient": source_payload["decomposition"]["contact_coefficient"],
            "source_rank_contact_fraction": float(decomp["contact_fraction"]),
            "veronese_only_relative_residual": float(decomp["veronese_only_relative_residual"]),
            "veronese_plus_contact_relative_residual": float(decomp["veronese_plus_contact_relative_residual"]),
            "MR_singular_values_GeV": decomp["singular_values_GeV"],
            "MR_condition_number": float(decomp["condition_number"]),
            "MR_normalized_digest": stable_digest([[cpair(z) for z in row] for row in mr_normalized]),
            "M_V_quartic_invariants": veronese_quartic_invariants(mv),
        },
        "z178_diagnostic": {
            "phase_fraction_arg_zeta_over_2pi": x_phase,
            "z178_order": zeta_payload["z178_order"],
            "z178_power": zeta_payload["z178_power"],
            "z178_phase_error_rad": zeta_payload["phase_error"],
            "z178_complex_error": zeta_payload["zeta_178_error"],
            "continued_fraction_convergents_selected": relevant_cf,
            "policy": "Diophantine side diagnostic only; not a phase prediction or benchmark premise",
        },
        "delta_ZN_replay_protocol": {
            "estimate": "delta Z_N ~ |lambda|^2/(16 pi^2) log(M_*/M_X)",
            "coefficient_without_log": float(abs(sqrt_zeta) ** 2 / (16.0 * math.pi**2)),
            "Y_nuD_replay": "Y_nuD -> Y_nuD (I - 1/2 deltaZ_N)",
            "MR_replay": "M_R -> M_R - 1/2 (deltaZ_N^T M_R + M_R deltaZ_N)",
            "status": "must be applied in Audit 1 if a finite messenger interval is specified",
        },
        "legacy_local_conventions_index": {
            "no_web_constant_row_count": no_web_payload["constant_row_count"],
            "no_web_constant_rows_sha256": no_web_payload["constant_rows_sha256"],
            "flavor_card_all_checks_ok": flavor_payload.get("all_checks_ok"),
            "flavor_basis": flavor_payload["basis"],
            "contact_sensitivity": {
                "real_scale_loose_half_width": sensitivity_payload["real_scale_loose_interval"]["half_width"],
                "phase_loose_half_width_rad": sensitivity_payload["phase_loose_interval"]["half_width"],
                "real_scale_tight_half_width": sensitivity_payload["real_scale_tight_interval"]["half_width"],
                "phase_tight_half_width_rad": sensitivity_payload["phase_tight_interval"]["half_width"],
            },
        },
        "deferred_audit_dependency_graph": {
            "order": "0 -> (1 || 4a) -> 3 -> 2; 4b optional parallel",
            "nodes": {
                "0": "convention/invariant card",
                "1": "full flavor and seesaw fit",
                "4a": "field-theory source-sector UV/heavy-spectrum interface",
                "4b": "optional Route-D global/string geometry track",
                "3": "threshold and unification audit using source spectrum",
                "2": "d=5 proton Wilson/lifetime audit consuming flavor+source+threshold outputs",
            },
            "fallback_prices": {
                "flavor_fit_fail": "add 120_H antisymmetric channel or second spurion; record +6 real-parameter cost per added channel",
                "source_uv_fail": "switch breaking chain or demote source-sector origin to benchmark datum",
                "threshold_fail": "change breaking chain, open intermediate scale, or compute Route-D flux threshold; record threshold naturalness cost",
                "d5_fail": "raise triplet masses within threshold band, use flux doublet-triplet mechanism, or move to non-SUSY branch with rewritten beta/superpotential language",
            },
        },
        "publication_boundary": {
            "not_claimed_by_audit0": [
                "full CKM/PMNS fit",
                "d=5 proton lifetime bound",
                "threshold closure",
                "source-sector UV completion",
                "global F-theory compactification",
                "ten-digit physical prediction of zeta",
            ],
            "claim": "Audit 0 fixes conventions and invariant anchors only.",
        },
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "invariant_card.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Audit 0 Convention / Invariant Card",
        "",
        "Audit 0 fixes the common convention layer for the deferred companion audits. It is not a flavor fit, threshold closure, proton lifetime calculation, or UV completion.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- source artifacts: {sum(1 for row in card['source_manifest'] if row['exists'])}/{len(card['source_manifest'])} present",
        "",
        "## Core Anchors",
        "",
        f"- `zeta = {zeta.real:.10f} + {zeta.imag:.10f} i`",
        f"- `|zeta| = {abs(zeta):.17g}`",
        f"- `arg(zeta) = {np.angle(zeta):.15g} rad`",
        f"- `|sqrt(zeta)| = {abs(sqrt_zeta):.16g}`",
        f"- `|lambda|^2/(16 pi^2) = {card['delta_ZN_replay_protocol']['coefficient_without_log']:.9e}`",
        f"- `K_tr inverse check = {card['routeA_invariant_anchors']['k_tr']['inverse_check_norm_K_times_3K_minus_I']:.3e}`",
        f"- `M_R contact fraction = {decomp['contact_fraction']:.9e}`",
        f"- `Veronese + contact residual = {decomp['veronese_plus_contact_relative_residual']:.3e}`",
        "",
        "## Dependency Graph",
        "",
        "`0 -> (1 || 4a) -> 3 -> 2`, with `4b` optional and parallel.",
        "",
        "| audit | role |",
        "|---|---|",
        "| 0 | convention/invariant card |",
        "| 1 | full flavor and seesaw fit |",
        "| 4a | field-theory source-sector UV/heavy-spectrum interface |",
        "| 4b | optional Route-D global/string geometry track |",
        "| 3 | threshold and unification audit using source spectrum |",
        "| 2 | d=5 proton Wilson/lifetime audit consuming all previous outputs |",
        "",
        "## Precision Policy",
        "",
        "Long decimals are retained as reproducibility anchors. They are not ten-digit physical predictions; follow-up audits should treat only the first few significant figures as physically meaningful unless a concrete UV replay justifies more.",
        "",
        "## Next Required Inputs",
        "",
        "- Audit 1 must set the final CKM/PMNS/mass target table and fit convention.",
        "- Audit 4a must provide a heavy-spectrum schema before Audit 3 can be paper-grade.",
        "- Audit 2 must not be run as a final proton-lifetime audit until flavor rotations, source spectrum, and threshold-fixed scales are supplied.",
        "",
    ]
    md_path = OUT / "invariant_card.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {md_path.relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")


if __name__ == "__main__":
    main()
