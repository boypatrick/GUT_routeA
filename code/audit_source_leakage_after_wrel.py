#!/usr/bin/env python3
"""Source self-coupling and kinetic-leakage audit after W_rel.

No web lookup is used.  This script returns to the full source-stabilization
sector after adding

    W_rel = <Y_AB, S_AB - v_AB S_AA> + <Y_AC, S_AC - v_AC S_AA>.

It checks whether source self-couplings or source/link Kahler mixing can
regenerate non-scalar Pati-Salam fragment splittings.  The answer is:

* scalar source masses and scalar kinetic mixings remain complete multiplets;
* neutral source self-coupling leakage such as Omega^3 or S_AA Omega^2 is not
  removed by the Route-B endpoint grading alone once an R=2 leakage spurion is
  admitted;
* non-scalar Kahler metrics also canonicalize into fragment-dependent masses.

Both leakage mechanisms are therefore mapped to epsilon54/epsilon210-like
threshold parameters and replayed through the cached RGE/proton scan.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_drive_sector_spectrum as drive  # noqa: E402
import audit_dynamic_source_mixing as dyn  # noqa: E402
import scan_untruncated_invariant_deformations as route_a  # noqa: E402
import scan_yukawa_superpotential_rge as base  # noqa: E402
import verify_spin10_component_hessian as comp  # noqa: E402


OUT = ROOT / "output" / "source_leakage_after_wrel"
CARD = comp.load_card()
PARS = CARD["vev_parameters_and_couplings"]
ELL210 = float(PARS["b210"])
ELLS54 = [float(PARS["a54"]), float(PARS["eta"]), float(PARS["eta"])]
_RGE_CONTEXT: dict[str, Any] | None = None


SCAN_VALUES = np.array(
    [
        -3.0e-3,
        -2.0e-3,
        -1.5e-3,
        -1.0e-3,
        -7.5e-4,
        -5.0e-4,
        -3.0e-4,
        -2.0e-4,
        -1.0e-4,
        -5.0e-5,
        0.0,
        5.0e-5,
        1.0e-4,
        2.0e-4,
        3.0e-4,
        5.0e-4,
        7.5e-4,
        1.0e-3,
        1.5e-3,
        2.0e-3,
        3.0e-3,
    ],
    dtype=float,
)


def selection_rule_audit() -> list[dict[str, Any]]:
    """Human-readable selection-rule audit.

    The combined G_sel keeps W_rel compatible with the rest of the model, but
    source-only neutral operators can reappear if one admits an R=2 leakage
    spurion.  Those operators are exactly the ones represented by epsilons in
    the scan below.
    """

    return [
        {
            "term": "scalar source masses M_e S_e^2 and M_Omega Omega^2",
            "status": "allowed as R=2 mass spurions",
            "effect": "safe if Spin(10)-scalar; complete multiplet threshold only",
        },
        {
            "term": "Omega_AA^3",
            "status": "forbidden in strict R-sequestered sector; allowed if an R=2 neutral leakage spurion is admitted",
            "effect": "dangerous 210 source Hessian; mapped to epsilon210",
        },
        {
            "term": "S_AA Omega_AA^2",
            "status": "same neutral leakage class as Omega^3",
            "effect": "dangerous 210 PS-fragment splitting once S_AA gets aligned vev",
        },
        {
            "term": "S_AB Omega_AA^2 or S_AC Omega_AA^2",
            "status": "requires endpoint charge spurions",
            "effect": "higher-order leakage; mapped to epsilon210 if present",
        },
        {
            "term": "Phi_C^3",
            "status": "allowed",
            "effect": "safe because Phi_C is a distinct color source and is sequestered from Omega_AA",
        },
        {
            "term": "same-endpoint Kahler S_e^dagger L_e",
            "status": "allowed",
            "effect": "safe if Spin(10)-scalar; dangerous if fragment-dependent",
        },
        {
            "term": "cross-endpoint Kahler S_AA^dagger S_AB and S_AA^dagger S_AC",
            "status": "allowed only with v_AB^dagger/v_AC^dagger spurions",
            "effect": "safe if scalar and aligned; non-scalar part is kinetic epsilon54",
        },
    ]


def source210_matrix(source_mass: float = 1.0) -> np.ndarray:
    return drive.source_link_driver_matrix(source_mass, ELL210)


def relative_source_alignment_matrix(source_mass: float = 1.0) -> np.ndarray:
    """Fast W_rel-completed 54 Hessian with cached Spin(10) vev ratios."""

    matrix = np.zeros((11, 11), dtype=float)
    for endpoint, ell in enumerate(ELLS54):
        s_pos = endpoint
        l_pos = 3 + endpoint
        d_pos = 6 + endpoint
        matrix[s_pos, s_pos] = source_mass
        matrix[d_pos, l_pos] = 1.0
        matrix[l_pos, d_pos] = 1.0
        matrix[d_pos, s_pos] = -ell
        matrix[s_pos, d_pos] = -ell

    y_ab = 9
    y_ac = 10
    matrix[y_ab, 1] = matrix[1, y_ab] = 1.0
    matrix[y_ab, 0] = matrix[0, y_ab] = -1.0
    matrix[y_ac, 2] = matrix[2, y_ac] = 1.0
    matrix[y_ac, 0] = matrix[0, y_ac] = -1.0
    return 0.5 * (matrix + matrix.T)


def kinetic_pattern_54() -> np.ndarray:
    """Representative source/link Kahler leakage pattern for the 11x11 block."""

    pattern = np.zeros((11, 11), dtype=float)
    for endpoint in range(3):
        s_pos = endpoint
        l_pos = 3 + endpoint
        pattern[s_pos, s_pos] = 1.0
        pattern[l_pos, l_pos] = -0.5
        pattern[s_pos, l_pos] = pattern[l_pos, s_pos] = 0.25
    # Include the relative-source directions; this tests that W_rel does not
    # hide a kinetic relative-orientation leak.
    pattern[0, 1] = pattern[1, 0] = 0.2
    pattern[0, 2] = pattern[2, 0] = -0.2
    norm = float(np.linalg.norm(pattern, 2))
    return pattern / max(norm, 1.0e-300)


def kinetic_pattern_210() -> np.ndarray:
    pattern = np.zeros((3, 3), dtype=float)
    pattern[0, 0] = 1.0
    pattern[1, 1] = -0.5
    pattern[0, 1] = pattern[1, 0] = 0.25
    norm = float(np.linalg.norm(pattern, 2))
    return pattern / max(norm, 1.0e-300)


def canonical_log_sum(matrix: np.ndarray, kinetic_strength: float, split: float, pattern: np.ndarray) -> tuple[np.ndarray, float]:
    if abs(kinetic_strength) < 1.0e-300:
        return drive.log_sum_abs_eigenvalues(matrix)
    zmat = np.eye(matrix.shape[0]) + kinetic_strength * split * pattern
    evals_z, evecs_z = np.linalg.eigh(0.5 * (zmat + zmat.T))
    if float(np.min(evals_z)) <= 0.0:
        raise ValueError("non-positive kinetic matrix in leakage scan")
    invsqrt = evecs_z @ np.diag(1.0 / np.sqrt(evals_z)) @ evecs_z.T
    canon = invsqrt.T @ matrix @ invsqrt
    return drive.log_sum_abs_eigenvalues(canon)


def leakage_threshold(
    epsilon54: float = 0.0,
    epsilon210: float = 0.0,
    zeta54: float = 0.0,
    zeta210: float = 0.0,
) -> dict[str, Any]:
    delta54 = np.zeros(3, dtype=float)
    delta210 = np.zeros(3, dtype=float)
    rows54: dict[str, Any] = {}
    rows210: dict[str, Any] = {}
    kpat54 = kinetic_pattern_54()
    kpat210 = kinetic_pattern_210()

    for name, item in dyn.FRAGMENTS_54.items():
        split = float(item["split"])
        mass = 1.0 + epsilon54 * split
        matrix = relative_source_alignment_matrix(source_mass=mass)
        evals, log_sum = canonical_log_sum(matrix, zeta54, split, kpat54)
        beta = np.array(item["beta"], dtype=float)
        delta54 += beta * log_sum / (2.0 * math.pi)
        rows54[name] = {
            "split": split,
            "source_mass": mass,
            "log_sum": log_sum,
            "min_abs_eigenvalue": float(np.min(np.abs(evals))),
            "beta": beta.tolist(),
        }

    for name, item in dyn.FRAGMENTS_210.items():
        split = float(item["split"])
        mass = 1.0 + epsilon210 * split
        matrix = source210_matrix(source_mass=mass)
        evals, log_sum = canonical_log_sum(matrix, zeta210, split, kpat210)
        beta = np.array(item["beta"], dtype=float)
        delta210 += beta * log_sum / (2.0 * math.pi)
        rows210[name] = {
            "split": split,
            "source_mass": mass,
            "log_sum": log_sum,
            "min_abs_eigenvalue": float(np.min(np.abs(evals))),
            "beta": beta.tolist(),
        }

    total = delta54 + delta210
    eig54 = [row["log_sum"] for row in rows54.values()]
    eig210 = [row["log_sum"] for row in rows210.values()]
    return {
        "epsilon54": epsilon54,
        "epsilon210": epsilon210,
        "zeta54": zeta54,
        "zeta210": zeta210,
        "delta54": delta54.tolist(),
        "delta210": delta210.tolist(),
        "delta_total": total.tolist(),
        "projected_delta_l2": float(np.linalg.norm(base.PROJECTOR @ total)),
        "universal_delta": float(np.mean(total)),
        "log_sum_spread_54": float(max(eig54) - min(eig54)),
        "log_sum_spread_210": float(max(eig210) - min(eig210)),
        "fragments54": rows54,
        "fragments210": rows210,
    }


def replay_with_delta(delta: np.ndarray) -> dict[str, Any]:
    global _RGE_CONTEXT
    if _RGE_CONTEXT is None:
        _RGE_CONTEXT = {
            "metrics": dyn.exact_link_metrics(),
            "grid": route_a.load_cached_grid(),
            "prefactor": base.proton_prefactor(),
        }
    return dyn.fixed_spectrum_rge_scan_with_extra_delta(
        _RGE_CONTEXT["metrics"],
        _RGE_CONTEXT["grid"],
        _RGE_CONTEXT["prefactor"],
        delta,
    )


def scan_grid(kind: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for a in SCAN_VALUES:
        for b in SCAN_VALUES:
            if kind == "holomorphic":
                threshold = leakage_threshold(epsilon54=float(a), epsilon210=float(b))
                x_name, y_name = "epsilon54", "epsilon210"
            elif kind == "kinetic":
                threshold = leakage_threshold(zeta54=float(a), zeta210=float(b))
                x_name, y_name = "zeta54", "zeta210"
            else:
                raise ValueError(kind)
            rge = replay_with_delta(np.array(threshold["delta_total"], dtype=float))
            safe = bool(rge["safe_count_residual_lt_1e_3"] > 0)
            rows.append(
                {
                    "kind": kind,
                    x_name: float(a),
                    y_name: float(b),
                    "projected_delta_l2": threshold["projected_delta_l2"],
                    "universal_delta": threshold["universal_delta"],
                    "log_sum_spread_54": threshold["log_sum_spread_54"],
                    "log_sum_spread_210": threshold["log_sum_spread_210"],
                    "safe_count": int(rge["safe_count_residual_lt_1e_3"]),
                    "best_residual_l2": float(rge["best_residual_l2"]),
                    "best_alphaG_inv": float(rge["best_alphaG_inv"]),
                    "tau_d6": float(rge["best_tau_dim6_years"]),
                    "tau_d5": float(rge["best_tau_dim5_ST_1e_5_years"]),
                    "safe": safe,
                }
            )

    x_key, y_key = ("epsilon54", "epsilon210") if kind == "holomorphic" else ("zeta54", "zeta210")
    safe_rows = [row for row in rows if row["safe"]]
    unsafe_rows = [row for row in rows if not row["safe"]]

    def axis_rows(x: float | None = None, y: float | None = None) -> list[dict[str, Any]]:
        out = rows
        if x is not None:
            out = [row for row in out if abs(row[x_key] - x) < 1.0e-15]
        if y is not None:
            out = [row for row in out if abs(row[y_key] - y) < 1.0e-15]
        return out

    first_unsafe = (
        min(
            unsafe_rows,
            key=lambda row: math.hypot(float(row[x_key]), float(row[y_key])),
        )
        if unsafe_rows
        else None
    )
    summary = {
        "kind": kind,
        "points": len(rows),
        "safe_points": len(safe_rows),
        "safe_fraction": len(safe_rows) / len(rows),
        f"max_abs_{x_key}_safe_on_{y_key}_0_axis": max(
            (abs(row[x_key]) for row in axis_rows(y=0.0) if row["safe"]),
            default=0.0,
        ),
        f"max_abs_{y_key}_safe_on_{x_key}_0_axis": max(
            (abs(row[y_key]) for row in axis_rows(x=0.0) if row["safe"]),
            default=0.0,
        ),
        "all_points_safe_inside_1e_minus_4_square": all(
            row["safe"]
            for row in rows
            if abs(row[x_key]) <= 1.0e-4 + 1.0e-15 and abs(row[y_key]) <= 1.0e-4 + 1.0e-15
        ),
        "all_points_safe_inside_3e_minus_4_square": all(
            row["safe"]
            for row in rows
            if abs(row[x_key]) <= 3.0e-4 + 1.0e-15 and abs(row[y_key]) <= 3.0e-4 + 1.0e-15
        ),
        "first_unsafe_by_radius": first_unsafe,
    }
    return rows, summary


def spurion_order(lambda_src: float, epsilon_bound: float) -> dict[str, Any]:
    rows = []
    for order in range(1, 10):
        value = lambda_src**order
        rows.append({"order": order, "value": value, "passes": value <= epsilon_bound})
    first = next(row for row in rows if row["passes"])
    return {
        "lambda": lambda_src,
        "epsilon_bound": epsilon_bound,
        "orders": rows,
        "first_safe_order": first["order"],
    }


def write_csv(rows: list[dict[str, Any]], filename: str) -> None:
    fieldnames = list(rows[0].keys())
    with (OUT / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    hsum = payload["holomorphic_summary"]
    ksum = payload["kinetic_summary"]
    base0 = payload["baseline_wrel_source"]
    lines: list[str] = []
    lines.append("# Source leakage after W_rel")
    lines.append("")
    lines.append("No web lookup was used.  This audit checks whether source self-couplings")
    lines.append("or source/link Kahler mixing regenerate non-scalar PS-fragment splittings")
    lines.append("after the relative-alignment driver `W_rel` is included.")
    lines.append("")
    lines.append("## Selection-rule audit")
    lines.append("")
    lines.append("| term/class | status | effect |")
    lines.append("|---|---|---|")
    for row in payload["selection_rule_audit"]:
        lines.append(f"| `{row['term']}` | {row['status']} | {row['effect']} |")
    lines.append("")
    lines.append("The strict model keeps the source self-Hessian Spin(10)-scalar.  Any")
    lines.append("R=2 neutral source leakage spurion is therefore not a harmless detail; it")
    lines.append("is precisely an epsilon54/epsilon210 deformation.")
    lines.append("")
    lines.append("## W_rel-completed scalar baseline")
    lines.append("")
    lines.append("```text")
    lines.append(f"P Delta_source+rel l2 = {base0['projected_delta_l2']:.3e}")
    lines.append(f"universal Delta        = {base0['universal_delta']:.6e}")
    lines.append(f"log-sum spread 54      = {base0['log_sum_spread_54']:.3e}")
    lines.append(f"log-sum spread 210     = {base0['log_sum_spread_210']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Holomorphic source-Hessian leakage")
    lines.append("")
    lines.append("The scan uses")
    lines.append("")
    lines.append("```text")
    lines.append("m_54(f)  = m_54  [1 + epsilon54 c_f],")
    lines.append("m_210(f) = m_210 [1 + epsilon210 c_f],")
    lines.append("```")
    lines.append("")
    lines.append("but with the 54 matrix replaced by the W_rel-completed 11x11 block.")
    lines.append("")
    lines.append("```text")
    lines.append(f"points                         = {hsum['points']}")
    lines.append(f"safe fraction                  = {hsum['safe_fraction']:.3f}")
    lines.append(
        "max |epsilon54| safe at epsilon210=0 = "
        f"{hsum['max_abs_epsilon54_safe_on_epsilon210_0_axis']:.3e}"
    )
    lines.append(
        "max |epsilon210| safe at epsilon54=0 = "
        f"{hsum['max_abs_epsilon210_safe_on_epsilon54_0_axis']:.3e}"
    )
    lines.append(f"all safe inside 1e-4 square    = {hsum['all_points_safe_inside_1e_minus_4_square']}")
    lines.append(f"all safe inside 3e-4 square    = {hsum['all_points_safe_inside_3e_minus_4_square']}")
    bad = hsum["first_unsafe_by_radius"]
    if bad is None:
        lines.append("first unsafe by radius         = none on scanned grid")
    else:
        lines.append(
            "first unsafe by radius         = "
            f"eps54={bad['epsilon54']:.1e}, eps210={bad['epsilon210']:.1e}, "
            f"residual={bad['best_residual_l2']:.3e}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Kahler kinetic leakage")
    lines.append("")
    lines.append("A scalar Kahler metric is harmless.  The test deformation is")
    lines.append("")
    lines.append("```text")
    lines.append("Z_f = 1 + zeta c_f K_source-link,")
    lines.append("M_canon(f) = Z_f^{-1/2 T} M(f) Z_f^{-1/2}.")
    lines.append("```")
    lines.append("")
    lines.append("This represents the leading non-scalar source/link kinetic mixing.")
    lines.append("")
    lines.append("```text")
    lines.append(f"points                      = {ksum['points']}")
    lines.append(f"safe fraction               = {ksum['safe_fraction']:.3f}")
    lines.append(
        "max |zeta54| safe at zeta210=0 = "
        f"{ksum['max_abs_zeta54_safe_on_zeta210_0_axis']:.3e}"
    )
    lines.append(
        "max |zeta210| safe at zeta54=0 = "
        f"{ksum['max_abs_zeta210_safe_on_zeta54_0_axis']:.3e}"
    )
    lines.append(f"all safe inside 1e-4 square = {ksum['all_points_safe_inside_1e_minus_4_square']}")
    lines.append(f"all safe inside 3e-4 square = {ksum['all_points_safe_inside_3e_minus_4_square']}")
    bad = ksum["first_unsafe_by_radius"]
    if bad is None:
        lines.append("first unsafe by radius      = none on scanned grid")
    else:
        lines.append(
            "first unsafe by radius      = "
            f"zeta54={bad['zeta54']:.1e}, zeta210={bad['zeta210']:.1e}, "
            f"residual={bad['best_residual_l2']:.3e}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Spurion order")
    lines.append("")
    for row in payload["spurion_orders"]:
        lines.append(
            f"For lambda={row['lambda']}, epsilon <= {row['epsilon_bound']:.1e} "
            f"requires first dangerous order >= {row['first_safe_order']}."
        )
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append("Adding `W_rel` does not itself regenerate a non-scalar threshold.  The")
    lines.append("branch remains safe if source self-couplings and source/link Kahler terms")
    lines.append("are Spin(10)-scalar, or if their non-scalar leakage is kept at roughly")
    lines.append("`1e-4` or below.  Neutral source self-couplings such as `Omega^3` and")
    lines.append("`S_AA Omega^2` must therefore remain absent by the source-stabilization")
    lines.append("symmetry, or be counted explicitly as epsilon210 leakage.")
    (OUT / "source_leakage_after_wrel_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = leakage_threshold()
    hol_rows, hol_summary = scan_grid("holomorphic")
    kin_rows, kin_summary = scan_grid("kinetic")
    epsilon_bound = 1.0e-4
    payload = {
        "note": "No web lookup used. Source self-coupling and kinetic-leakage audit after W_rel.",
        "selection_rule_audit": selection_rule_audit(),
        "baseline_wrel_source": baseline,
        "holomorphic_scan": hol_summary,
        "kinetic_scan": kin_summary,
        "holomorphic_summary": hol_summary,
        "kinetic_summary": kin_summary,
        "spurion_orders": [
            spurion_order(0.1, epsilon_bound),
            spurion_order(0.2, epsilon_bound),
        ],
        "passes": {
            "baseline_projected_threshold_zero": baseline["projected_delta_l2"] < 1.0e-12,
            "holomorphic_1e_minus_4_square_safe": hol_summary[
                "all_points_safe_inside_1e_minus_4_square"
            ],
            "kinetic_1e_minus_4_square_safe": kin_summary[
                "all_points_safe_inside_1e_minus_4_square"
            ],
        },
    }
    payload["passes"]["all"] = all(payload["passes"].values())
    (OUT / "source_leakage_after_wrel_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_csv(hol_rows, "source_self_leakage_scan.csv")
    write_csv(kin_rows, "source_kinetic_leakage_scan.csv")
    write_report(payload)

    print("Source leakage after W_rel")
    print(f"  baseline projected threshold: {baseline['projected_delta_l2']:.3e}")
    print(f"  holomorphic safe fraction: {hol_summary['safe_fraction']:.3f}")
    print(f"  kinetic safe fraction: {kin_summary['safe_fraction']:.3f}")
    print(
        "  holomorphic max |eps210| axis: "
        f"{hol_summary['max_abs_epsilon210_safe_on_epsilon54_0_axis']:.3e}"
    )
    print(
        "  kinetic max |zeta210| axis: "
        f"{kin_summary['max_abs_zeta210_safe_on_zeta54_0_axis']:.3e}"
    )
    print(f"  all checks: {payload['passes']['all']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
