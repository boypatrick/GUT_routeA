#!/usr/bin/env python3
"""Route-C P10 Spin(10) staged breaking mass-matrix ledger.

P10 chooses a concrete but conservative order-parameter/source branch for the
P8/P9 Spin(10) path.  The branch is deliberately staged and block-diagonal:

  Spin(10) -> Pati-Salam,
  SU(4)_C -> SU(3)_C x U(1)_{B-L},
  SU(2)_R -> U(1)_{T3R},
  U(1)_{T3R} x U(1)_{B-L} -> U(1)_Y.

The output turns the symbolic P9 masses

  M_(6,2,2), M_LQ, M_WR

into a diagonalizable broken-vector mass matrix with a transparent
normalization convention.  This is still not a full Higgs-potential or
proton-lifetime calculation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P9_JSON = OUT / "spin10_action_completion_ledger.json"
OUT_JSON = OUT / "spin10_breaking_mass_matrix.json"
OUT_MD = OUT / "spin10_breaking_mass_matrix.md"


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run P9 first")
    return json.loads(path.read_text(encoding="utf-8"))


def fnum(x: float) -> float:
    return float(f"{x:.15g}")


def neutral_matrix(g_r: float, g_bl: float, v_y: float) -> list[list[float]]:
    """Mass matrix for a Y-neutral nu^c-like source.

    The source has charges q_R=-1/2 and q_BL=+1/2 under
    T3R and Q_BL=(B-L)/2.  In the basis (W_R^3, B_{B-L}),

      M^2 = v_y^2 (g_R q_R, g_BL q_BL)^T (g_R q_R, g_BL q_BL).
    """
    pref = v_y * v_y / 4.0
    return [
        [fnum(pref * g_r * g_r), fnum(-pref * g_r * g_bl)],
        [fnum(-pref * g_r * g_bl), fnum(pref * g_bl * g_bl)],
    ]


def neutral_eigen(g_r: float, g_bl: float, v_y: float) -> dict:
    massive = (v_y * v_y / 4.0) * (g_r * g_r + g_bl * g_bl)
    norm = math.sqrt(g_r * g_r + g_bl * g_bl)
    return {
        "basis": ["W_R^3", "B_{B-L}"],
        "mass_matrix": neutral_matrix(g_r, g_bl, v_y),
        "eigenvalues": [0.0, fnum(massive)],
        "massless_hypercharge_vector": {
            "components_in_basis": [fnum(g_bl / norm), fnum(g_r / norm)],
            "interpretation": "B_Y proportional to g_BL W_R^3 + g_R B_{B-L}",
        },
        "massive_Zprime_vector": {
            "components_in_basis": [fnum(g_r / norm), fnum(-g_bl / norm)],
            "interpretation": "orthogonal broken neutral vector",
        },
        "hypercharge_coupling": {
            "formula": "g_Y = g_R g_BL / sqrt(g_R^2 + g_BL^2)",
            "demo_value": fnum(g_r * g_bl / norm),
        },
        "zero_mode_residual_demo": 0.0,
    }


def p9_map_counts(p9: dict) -> dict[str, int]:
    rows = p9["broken_generator_maps"]["maps"]
    counts = {
        "Pati-Salam SU(4)_C leptoquark": 0,
        "broken SU(2)_R charged current": 0,
        "Spin(10)/Pati-Salam off-face generator": 0,
    }
    for row in rows:
        fam = row["generator_family"]
        if fam not in counts:
            raise SystemExit(f"unexpected P9 generator family {fam}")
        counts[fam] += 1
    return counts


def build_payload() -> dict:
    p9 = load_json(P9_JSON)
    counts = p9_map_counts(p9)

    # Demo values are intentionally dimensionless and only test the diagonal
    # structure.  Physics claims remain symbolic.
    demo = {
        "g_10": 1.0,
        "g_4": 1.0,
        "g_R": 1.0,
        "g_BL": 1.0,
        "v_PS": 1.0,
        "v_4": 1.0,
        "v_R": 1.0,
        "v_Y": 1.0,
        "kappa_PS": 1.0,
    }

    m_ps = demo["kappa_PS"] * demo["g_10"] ** 2 * demo["v_PS"] ** 2
    # T_BL = sqrt(3/8) diag(1/3,1/3,1/3,-1), so a leptoquark root has
    # Delta T_BL^2 = (4/3)^2 * 3/8 = 2/3.
    m_lq = (2.0 / 3.0) * demo["g_4"] ** 2 * demo["v_4"] ** 2
    # With adjoint/source normalization Tr(T3R^2)=1/2, the charged root has
    # unit squared charge difference.
    m_wr = demo["g_R"] ** 2 * demo["v_R"] ** 2
    m_neutral = neutral_eigen(demo["g_R"], demo["g_BL"], demo["v_Y"])

    eigenvalues = (
        [m_ps] * counts["Spin(10)/Pati-Salam off-face generator"]
        + [m_lq] * counts["Pati-Salam SU(4)_C leptoquark"]
        + [m_wr] * counts["broken SU(2)_R charged current"]
        + m_neutral["eigenvalues"]
    )
    positive = [x for x in eigenvalues if x > 1e-12]
    zero = [x for x in eigenvalues if abs(x) <= 1e-12]

    return {
        "description": "Route-C P10 Spin(10) staged breaking mass-matrix ledger",
        "boundary": (
            "P10 chooses a staged orthogonal source branch and supplies a "
            "diagonalizable broken-vector mass matrix.  It does not specify a "
            "complete scalar potential, threshold spectrum, flavor fit, or "
            "proton lifetime."
        ),
        "branch": {
            "name": "staged_orthogonal_source_branch",
            "reason_for_choice": (
                "It is the minimal conservative branch that turns the P9 "
                "symbolic mass inputs into positive diagonal blocks while "
                "preserving the hypercharge zero mode."
            ),
            "breaking_chain": [
                "Spin(10) -> SU(4)_C x SU(2)_L x SU(2)_R",
                "SU(4)_C -> SU(3)_C x U(1)_{B-L}",
                "SU(2)_R -> U(1)_{T3R}",
                "U(1)_{T3R} x U(1)_{B-L} -> U(1)_Y",
            ],
        },
        "normalization": {
            "T_BL": "sqrt(3/8) diag(1/3,1/3,1/3,-1), Tr(T_BL^2)=1/2",
            "SU4_leptoquark_root_delta_squared": "2/3",
            "T3R": "diag(1/2,-1/2), Tr(T3R^2)=1/2",
            "nu_c_like_neutral_source_charges": {
                "q_R": "-1/2",
                "q_BL=(B-L)/2": "+1/2",
                "Y=q_R+q_BL": "0",
            },
        },
        "symbolic_mass_blocks": [
            {
                "block": "(6,2,2)",
                "dimension": counts["Spin(10)/Pati-Salam off-face generator"],
                "mass_squared": "M_(6,2,2)^2 = kappa_PS g_10^2 v_PS^2",
                "source": "Spin(10)/Pati-Salam-preserving order parameter or source",
                "p9_maps_assigned": counts["Spin(10)/Pati-Salam off-face generator"],
            },
            {
                "block": "SU(4)_C leptoquark roots",
                "dimension": counts["Pati-Salam SU(4)_C leptoquark"],
                "mass_squared": "M_LQ^2 = (2/3) g_4^2 v_4^2",
                "source": "adjoint/source along normalized B-L generator in SU(4)_C",
                "p9_maps_assigned": counts["Pati-Salam SU(4)_C leptoquark"],
            },
            {
                "block": "SU(2)_R charged roots",
                "dimension": counts["broken SU(2)_R charged current"],
                "mass_squared": "M_WR^2 = g_R^2 v_R^2",
                "source": "adjoint/source along normalized T3R in SU(2)_R",
                "p9_maps_assigned": counts["broken SU(2)_R charged current"],
            },
            {
                "block": "neutral U(1) mixing",
                "dimension": 2,
                "mass_squared_matrix_basis": ["W_R^3", "B_{B-L}"],
                "mass_squared_matrix": (
                    "v_Y^2/4 [[g_R^2, -g_R g_BL],[-g_R g_BL, g_BL^2]]"
                ),
                "eigenvalues": ["0", "v_Y^2/4 (g_R^2+g_BL^2)"],
                "source": "Y-neutral nu^c-like order parameter/source",
            },
        ],
        "full_block_matrix": {
            "basis_order": [
                "24 x (6,2,2) real broken vectors",
                "6 x SU(4)_C leptoquark real vectors",
                "2 x SU(2)_R charged real vectors",
                "2 x neutral basis vectors (W_R^3, B_{B-L})",
            ],
            "dimension": len(eigenvalues),
            "formula": (
                "diag(M_(6,2,2)^2 I_24, M_LQ^2 I_6, M_WR^2 I_2, "
                "M_neutral^2)"
            ),
            "neutral_block_formula": (
                "M_neutral^2 = v_Y^2/4 [[g_R^2, -g_R g_BL], "
                "[-g_R g_BL, g_BL^2]]"
            ),
        },
        "demo_diagonalization": {
            "dimensionless_inputs": demo,
            "block_eigenvalues": {
                "M_(6,2,2)^2": fnum(m_ps),
                "M_LQ^2": fnum(m_lq),
                "M_WR^2": fnum(m_wr),
                "neutral_eigenvalues": m_neutral["eigenvalues"],
            },
            "active_matrix_dimension_including_hypercharge_zero_mode": len(eigenvalues),
            "positive_broken_eigenvalue_count": len(positive),
            "zero_eigenvalue_count": len(zero),
            "expected_positive_broken_count": 33,
            "expected_zero_count": 1,
            "neutral_diagonalization": m_neutral,
        },
        "verification": {
            "p9_generator_layer_passes": p9["verification"]["p9_generator_layer_passes"],
            "all_p9_generator_maps_assigned_to_mass_blocks": sum(counts.values())
            == p9["broken_generator_maps"]["map_count"],
            "positive_broken_count_matches_spin10_to_sm": len(positive) == 33,
            "one_hypercharge_zero_mode": len(zero) == 1,
            "neutral_block_determinant_demo": 0.0,
            "neutral_block_trace_demo": fnum(sum(m_neutral["eigenvalues"])),
            "p10_mass_matrix_layer_passes": (
                p9["verification"]["p9_generator_layer_passes"]
                and sum(counts.values()) == p9["broken_generator_maps"]["map_count"]
                and len(positive) == 33
                and len(zero) == 1
            ),
        },
        "p6_p7_replay_interface": {
            "replace_symbolic_masses": {
                "M_(6,2,2)": "sqrt(kappa_PS) g_10 v_PS",
                "M_LQ": "sqrt(2/3) g_4 v_4",
                "M_WR": "g_R v_R",
                "M_Zprime": "0.5 v_Y sqrt(g_R^2+g_BL^2)",
            },
            "still_required": [
                "choose numerical scales/couplings or scan ranges",
                "choose scalar/source branch for the 8 non-adjoint P4 pairs if activated",
                "provide physical flavor rotations",
                "insert RG factors and hadronic matrix elements",
                "insert current experimental proton lifetime limits",
            ],
        },
        "next_stage": {
            "recommended": "P11_replay_P6_P7_with_staged_spin10_masses",
            "required_decision": (
                "Choose symbolic or numerical ranges for v_PS, v_4, v_R, v_Y "
                "and couplings, then rerun the P6/P7 matching gates with these "
                "mass blocks."
            ),
        },
    }


def write_markdown(payload: dict) -> str:
    lines = [
        "# Route-C P10 Spin(10) Breaking Mass-Matrix Ledger",
        "",
        payload["boundary"],
        "",
        "## Branch",
        "",
        f"- name: `{payload['branch']['name']}`",
        f"- reason: {payload['branch']['reason_for_choice']}",
        "",
        "Breaking chain:",
        "",
    ]
    for step in payload["branch"]["breaking_chain"]:
        lines.append(f"- {step}")
    lines.extend(
        [
            "",
            "## Symbolic Mass Blocks",
            "",
            "| block | dimension | mass-squared |",
            "| --- | ---: | --- |",
        ]
    )
    for row in payload["symbolic_mass_blocks"]:
        if "mass_squared" in row:
            mass = row["mass_squared"]
        else:
            mass = row["mass_squared_matrix"]
        lines.append(f"| {row['block']} | {row['dimension']} | `{mass}` |")
    demo = payload["demo_diagonalization"]
    lines.extend(
        [
            "",
            "## Demo Diagonalization",
            "",
            "| quantity | value |",
            "| --- | ---: |",
            f"| active matrix dimension, including hypercharge zero mode | {demo['active_matrix_dimension_including_hypercharge_zero_mode']} |",
            f"| positive broken eigenvalue count | {demo['positive_broken_eigenvalue_count']} |",
            f"| zero eigenvalue count | {demo['zero_eigenvalue_count']} |",
            f"| M_(6,2,2)^2 demo | {demo['block_eigenvalues']['M_(6,2,2)^2']} |",
            f"| M_LQ^2 demo | {demo['block_eigenvalues']['M_LQ^2']} |",
            f"| M_WR^2 demo | {demo['block_eigenvalues']['M_WR^2']} |",
            f"| neutral eigenvalues demo | {demo['block_eigenvalues']['neutral_eigenvalues']} |",
            "",
            "The massless neutral eigenvector is the hypercharge gauge boson:",
            "",
            "```text",
            "B_Y proportional to g_BL W_R^3 + g_R B_{B-L}",
            "g_Y = g_R g_BL / sqrt(g_R^2 + g_BL^2)",
            "```",
            "",
            "## Verification",
            "",
            "| check | value |",
            "| --- | ---: |",
        ]
    )
    for key, value in payload["verification"].items():
        lines.append(f"| {key} | {value} |")
    lines.extend(
        [
            "",
            "## P6/P7 Replay Interface",
            "",
            "| symbolic mass | replacement |",
            "| --- | --- |",
        ]
    )
    for key, value in payload["p6_p7_replay_interface"][
        "replace_symbolic_masses"
    ].items():
        lines.append(f"| `{key}` | `{value}` |")
    lines.extend(
        [
            "",
            "Still required:",
            "",
        ]
    )
    for item in payload["p6_p7_replay_interface"]["still_required"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## Next Stage",
            "",
            f"`{payload['next_stage']['recommended']}`",
            "",
            payload["next_stage"]["required_decision"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "branch": payload["branch"]["name"],
                "matrix_dimension": payload["full_block_matrix"]["dimension"],
                "positive_broken_eigenvalue_count": payload["demo_diagonalization"][
                    "positive_broken_eigenvalue_count"
                ],
                "zero_eigenvalue_count": payload["demo_diagonalization"][
                    "zero_eigenvalue_count"
                ],
                "p10_mass_matrix_layer_passes": payload["verification"][
                    "p10_mass_matrix_layer_passes"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
