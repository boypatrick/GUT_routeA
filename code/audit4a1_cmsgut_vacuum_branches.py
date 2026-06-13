#!/usr/bin/env python3
"""Build Audit 4a.1 CMSGUT vacuum-branch stage-1 scaffold.

This is deliberately a stage-1 contract, not a CMSGUT spectrum calculation.
It fixes the Pati-Salam singlet-vev dictionary, D-flatness convention,
literature special-point acceptance gates, and the convention-diff file that
must be closed before Audit 4a can export a non-placeholder heavy spectrum.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit4a1"

SOURCES = {
    "audit4a1_builder": ROOT / "code" / "audit4a1_cmsgut_vacuum_branches.py",
    "audit4a_schema": ROOT / "output" / "audit4a" / "source_spectrum_schema.json",
    "audit0_card": ROOT / "output" / "audit0" / "invariant_card.json",
}


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def ps_singlet_vevs() -> list[dict[str, Any]]:
    return [
        {
            "symbol": "p",
            "spin10_field": "210_H",
            "ps_block": "(1,1,1)",
            "sm_status": "SM singlet",
            "role": "Pati-Salam singlet vev",
            "stage1_status": "registered",
        },
        {
            "symbol": "a",
            "spin10_field": "210_H",
            "ps_block": "(15,1,1)",
            "sm_status": "SM singlet in the SU(4)_C adjoint",
            "role": "B-L-like adjoint vev",
            "stage1_status": "registered",
        },
        {
            "symbol": "omega",
            "spin10_field": "210_H",
            "ps_block": "(15,1,3)",
            "sm_status": "SM singlet in the SU(2)_R triplet-aligned block",
            "role": "right-triplet-aligned adjoint vev",
            "stage1_status": "registered",
        },
        {
            "symbol": "sigma",
            "spin10_field": "126_H",
            "ps_block": "(overline{10},1,3)",
            "sm_status": "SM singlet after right-handed alignment",
            "role": "B-L breaking source direction",
            "stage1_status": "registered",
        },
        {
            "symbol": "bar_sigma",
            "spin10_field": "overline{126}_H",
            "ps_block": "(10,1,3)",
            "sm_status": "SM singlet conjugate source direction",
            "role": "post-B-L P_nu^c source and D-flat conjugate",
            "stage1_status": "registered",
        },
    ]


def frac(value: int | tuple[int, int]) -> Fraction:
    if isinstance(value, tuple):
        return Fraction(value[0], value[1])
    return Fraction(value, 1)


def fstr(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def cubic_residual(x: Fraction, xi: Fraction) -> Fraction:
    return 8 * x**3 - 15 * x**2 + 14 * x - 3 + xi * (1 - x) ** 2


def vevs_for_x(x: Fraction) -> dict[str, Fraction]:
    return {
        "omega_tilde": -x,
        "a_tilde": (x**2 + 2 * x - 1) / (1 - x),
        "p_tilde": x * (5 * x**2 - 1) / (1 - x) ** 2,
        "eta_sigma_bar_sigma_tilde": 2 * x * (1 - 3 * x) * (1 + x**2) / (1 - x) ** 2,
    }


def special_point_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = [
        {
            "point_label": "SU5_x_half",
            "source_name": "SU(5)",
            "x": frac((1, 2)),
            "xi": frac(-5),
            "expected_vev_relation": "p_tilde=a_tilde=-omega_tilde=1/2",
        },
        {
            "point_label": "SU5_x_minus_one",
            "source_name": "SU(5)",
            "x": frac(-1),
            "xi": frac(10),
            "expected_vev_relation": "p_tilde=a_tilde=-omega_tilde=-1",
        },
        {
            "point_label": "GLR_x_zero",
            "source_name": "G_LR",
            "x": frac(0),
            "xi": frac(3),
            "expected_vev_relation": "p_tilde=omega_tilde=eta_sigma_bar_sigma_tilde=0, a_tilde=-1",
            "warning": "Primary source labels xi=3 as G_LR, not full Pati-Salam SU(4)_C x SU(2)_L x SU(2)_R.",
        },
        {
            "point_label": "flipped_SU5_x_third",
            "source_name": "flipped SU(5) x U(1)",
            "x": frac((1, 3)),
            "xi": frac((-2, 3)),
            "expected_vev_relation": "p_tilde=a_tilde=omega_tilde=-1/3 and eta_sigma_bar_sigma_tilde=0",
        },
    ]
    out: list[dict[str, Any]] = []
    for row in rows:
        x = row["x"]
        xi = row["xi"]
        residual = cubic_residual(x, xi)
        out.append(
            {
                **{k: v for k, v in row.items() if k not in {"x", "xi"}},
                "x_value": fstr(x),
                "xi_value": fstr(xi),
                "polynomial_residual": fstr(residual),
                "pass_fail": residual == 0,
                "dimensionless_vevs": {k: fstr(v) for k, v in vevs_for_x(x).items()},
            }
        )
    return out


def convention_items() -> list[dict[str, Any]]:
    return [
        {
            "item": "normalization_of_p_a_omega",
            "why_it_matters": "F-term cubic coefficients and mass matrices are convention-sensitive.",
            "status": "unresolved; must be mapped to the chosen BMSV/Aulakh-style table",
        },
        {
            "item": "sigma_bar_sigma_phase",
            "why_it_matters": "D-flatness fixes magnitudes, while the relative phase is a gauge/source convention.",
            "status": "stage-1 convention: record |sigma|=|bar_sigma| and allow sigma=bar_sigma after gauge fixing",
        },
        {
            "item": "branch_variable_x_and_xi",
            "why_it_matters": "The literature cubic P3(x;xi)=0 has multiple equivalent definitions.",
            "status": "fixed to Aulakh-Girdhar 2005: x=-lambda*omega/m, xi=lambda*M/(eta*m)",
        },
        {
            "item": "cubic_polynomial_coefficients",
            "why_it_matters": "Acceptance depends on reproducing enhanced-symmetry limits, not on a custom cubic.",
            "status": "fixed to 8*x^3 - 15*x^2 + 14*x - 3 = -xi*(1-x)^2",
        },
        {
            "item": "special_point_values",
            "why_it_matters": "SU(5), flipped SU(5), and left-right/Pati-Salam-related limits are the stage-1 validation anchors.",
            "status": "source-named SU(5), flipped SU(5), and G_LR points validated; no full-Pati-Salam point claimed",
        },
        {
            "item": "BMSV_vs_Aulakh_symbol_map",
            "why_it_matters": "Historical spectrum tables differ by signs and normalizations.",
            "status": "conventions-diff file created; unresolved until dual-source comparison is filled",
        },
    ]


def special_point_contract() -> list[dict[str, Any]]:
    rows = special_point_rows()
    rows.append(
        {
            "point_label": "full_Pati_Salam_requested_anchor",
            "source_name": "full Pati-Salam",
            "x_value": None,
            "xi_value": None,
            "polynomial_residual": None,
            "pass_fail": None,
            "stage1_status": "not claimed: primary source names xi=3 as G_LR rather than full Pati-Salam",
        }
    )
    return rows


def stage_gates() -> dict[str, Any]:
    return {
        "ps_singlet_vev_registry_complete": True,
        "d_flatness_condition_recorded": True,
        "cubic_coefficients_fixed_from_literature": True,
        "source_named_special_points_reproduced": True,
        "full_pati_salam_special_point_claimed": False,
        "symbolic_mass_matrices_exported": False,
        "goldstone_count_33_exported": False,
        "doublet_mass_matrix_exported": False,
        "triplet_inverse_block_exported": False,
        "heavy_spectrum_json_nonplaceholder": False,
    }


def write_conventions_diff(card: dict[str, Any]) -> None:
    lines = [
        "# Audit 4a.1 CMSGUT Conventions-Diff",
        "",
        "This file is intentionally unresolved.  It records the sign and",
        "normalization items that must be closed before the CMSGUT-like branch can",
        "claim a literature-compatible vacuum cubic or heavy spectrum.",
        "",
        "## Pati-Salam Singlet VEV Dictionary",
        "",
        "| symbol | field | PS block | role |",
        "| --- | --- | --- | --- |",
    ]
    for row in card["ps_singlet_vevs"]:
        lines.append(f"| `{row['symbol']}` | `{row['spin10_field']}` | `{row['ps_block']}` | {row['role']} |")

    lines += [
        "",
        "## D-Flatness Convention",
        "",
        "- Required condition: `|sigma| = |bar_sigma|`.",
        "- Working gauge convention allowed for stage 1: `sigma = bar_sigma` after",
        "  fixing the conjugate source phase.",
        "- This is only a convention statement; it is not a solved F/D-flat vacuum.",
        "",
        "## Unresolved Convention Items",
        "",
    ]
    for item in card["conventions_to_resolve"]:
        lines.append(f"- `{item['item']}`: {item['status']} ({item['why_it_matters']})")

    lines += [
        "",
        "## Literature Special-Point Acceptance",
        "",
        "Stage 1 reproduces the source-named enhanced-symmetry special points",
        "of the chosen Aulakh-Girdhar convention.  The source labels `xi=3` as",
        "`G_LR`; this scaffold therefore does not claim a full Pati-Salam",
        "special point.",
        "",
        "| point | source name | x | xi | residual | pass |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for item in card["literature_special_point_contract"]:
        lines.append(
            f"| `{item['point_label']}` | {item.get('source_name')} | "
            f"{item.get('x_value')} | {item.get('xi_value')} | "
            f"{item.get('polynomial_residual')} | {item.get('pass_fail')} |"
        )

    lines += [
        "",
        "## Boundary",
        "",
        "No threshold, proton-decay, or full heavy-spectrum conclusion may use this",
        "stage-1 scaffold until the cubic coefficients, special-point values,",
        "Goldstone count, and triplet inverse block are exported by a later audit.",
        "",
    ]
    (OUT / "conventions_diff.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    audit4a = read_json(SOURCES["audit4a_schema"])
    audit0 = read_json(SOURCES["audit0_card"])

    card: dict[str, Any] = {
        "audit": "audit4a1_cmsgut_vacuum_branches",
        "status": "stage-1 CMSGUT vacuum-branch scaffold only; no heavy spectrum is claimed",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "upstream_audit0_sha256": audit0.get("card_sha256"),
        "upstream_audit4a_sha256": audit4a.get("card_sha256"),
        "branch": {
            "name": "CMSGUT-like 210+126+overline126+10",
            "reason": "literature-compatible source-sector fallback with a 126bar_H post-B-L source for P_nu^c",
            "language": "Pati-Salam blocks, not raw SO(10) tensor components",
        },
        "superpotential_data_required_later": [
            "m_210",
            "m_126",
            "m_10",
            "lambda_210",
            "eta_210_126_126bar",
            "10-sector Yukawa/mixing couplings if doublet-triplet blocks are added",
        ],
        "ps_singlet_vevs": ps_singlet_vevs(),
        "d_flatness_contract": {
            "condition": "|sigma| = |bar_sigma|",
            "stage1_gauge_choice": "sigma = bar_sigma may be used after fixing the conjugate source phase",
            "claim_boundary": "recorded convention only; no solved F/D-flat vacuum claimed",
        },
        "literature_convention_map": {
            "primary_source": {
                "key": "AulakhGirdhar2005",
                "arxiv": "hep-ph/0405074",
                "journal": "Nucl. Phys. B711:275-313,2005",
                "scope": "Pati-Salam decomposition, vev dictionary, F-term equations, cubic, spectra/couplings",
            },
            "crosscheck_source": {
                "key": "AulakhBajcMelfoSenjanovicVissani2004",
                "arxiv": "hep-ph/0306242",
                "journal": "Phys. Lett. B588:196-202,2004",
                "scope": "minimal SUSY SO(10) setup and compact cubic/vev summary",
            },
            "branch_variable": "x = -lambda*omega/m",
            "coupling_ratio": "xi = lambda*M/(eta*m)",
            "cubic_equation": "8*x^3 - 15*x^2 + 14*x - 3 = -xi*(1-x)^2",
            "dimensionless_vevs": {
                "omega_tilde": "-x",
                "a_tilde": "(x^2+2*x-1)/(1-x)",
                "p_tilde": "x*(5*x^2-1)/(1-x)^2",
                "eta_sigma_bar_sigma_tilde": "2*x*(1-3*x)*(1+x^2)/(1-x)^2",
            },
            "abmsv_2003_sign_note": "The 2003 short note writes omega with the opposite sign; p and a formulas match after omega_ABMSV=-omega_AG.",
        },
        "vacuum_cubic_contract": {
            "dimensionless_branch_variable": "x",
            "dimensionless_coupling_ratio": "xi",
            "literature_template": "8*x^3 - 15*x^2 + 14*x - 3 + xi*(1-x)^2 = 0",
            "coefficients_status": "fixed from Aulakh-Girdhar 2005 / ABMSV 2004",
            "why_not_filled_here": "now filled for vacuum branch only; mass matrices remain future work",
        },
        "literature_special_point_contract": special_point_contract(),
        "conventions_to_resolve": convention_items(),
        "stage_gates": stage_gates(),
        "blocked_downstream_audits": {
            "Audit3_thresholds": "blocked until symbolic or numerical heavy_spectrum.json exports masses and beta vectors",
            "Audit2_d5_proton": "blocked until triplet inverse-propagator block and physical flavor rotations are exported",
        },
        "publication_boundary": {
            "claim": "The PS singlet-vev dictionary, Aulakh-Girdhar cubic convention, and source-named special-point validator for CMSGUT stage 1 are fixed.",
            "not_claimed": [
                "F/D-flat solution",
                "full Pati-Salam special point",
                "Goldstone count",
                "complete heavy spectrum",
                "doublet-triplet fine tuning",
                "colored-triplet inverse block",
                "threshold closure",
                "proton safety",
            ],
        },
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "vacuum_branches.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report_lines = [
        "# Audit 4a.1 CMSGUT Vacuum-Branch Stage-1 Scaffold",
        "",
        "This artifact fixes the Pati-Salam singlet-vev dictionary and the",
        "literature-special-point acceptance contract for the CMSGUT-like branch.",
        "It does not export a heavy spectrum.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- upstream Audit 4a sha256: `{card['upstream_audit4a_sha256']}`",
        f"- upstream Audit 0 sha256: `{card['upstream_audit0_sha256']}`",
        "",
        "## Registered Singlet VEVs",
        "",
        "| symbol | field | Pati-Salam block | role |",
        "| --- | --- | --- | --- |",
    ]
    for row in card["ps_singlet_vevs"]:
        report_lines.append(f"| `{row['symbol']}` | `{row['spin10_field']}` | `{row['ps_block']}` | {row['role']} |")

    report_lines += [
        "",
        "## Cubic Convention",
        "",
        "- Primary convention: Aulakh-Girdhar 2005 (`hep-ph/0405074`).",
        "- `x = -lambda*omega/m`.",
        "- `xi = lambda*M/(eta*m)`.",
        "- Cubic: `8*x^3 - 15*x^2 + 14*x - 3 = -xi*(1-x)^2`.",
        "- Cross-check: ABMSV 2004 (`hep-ph/0306242`) matches after the noted",
        "  omega-sign convention change.",
        "",
        "## Source-Named Special Points",
        "",
        "| point | source name | x | xi | residual | pass |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in card["literature_special_point_contract"]:
        report_lines.append(
            f"| `{row['point_label']}` | {row.get('source_name')} | "
            f"{row.get('x_value')} | {row.get('xi_value')} | "
            f"{row.get('polynomial_residual')} | {row.get('pass_fail')} |"
        )

    report_lines += [
        "",
        "## Stage Gates",
        "",
        "| gate | value |",
        "| --- | --- |",
    ]
    for key, value in card["stage_gates"].items():
        report_lines.append(f"| `{key}` | `{value}` |")

    report_lines += [
        "",
        "## Downstream Boundary",
        "",
        "- Audit 3 remains blocked until a non-placeholder heavy spectrum exports",
        "  masses and beta vectors.",
        "- Audit 2 remains blocked until physical flavor rotations and the colored",
        "  triplet inverse block are exported.",
        "",
    ]
    (OUT / "vacuum_branch_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    write_conventions_diff(card)

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {(OUT / 'vacuum_branch_report.md').relative_to(ROOT)}")
    print(f"wrote {(OUT / 'conventions_diff.md').relative_to(ROOT)}")
    print(f"digest {card['card_sha256']}")


if __name__ == "__main__":
    main()
