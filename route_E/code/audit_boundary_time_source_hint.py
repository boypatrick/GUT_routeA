#!/usr/bin/env python3
"""Audit the boundary-time hint behind the constrained-source branch.

The images suggest a recurring motif: dense regions do not support ordinary
bulk time; observable time and motion arise at interfaces, with rotation as a
possible near-critical boundary effect.  This script turns that motif into
three local checks:

1. A layer lapse N(n)=tanh(n/ell) has a Rindler limit at the boundary.
2. The constrained 54/210 source order parameters keep only orbit tangents
   plus radial singlets, removing normal modes from the propagating spectrum.
3. A rotating boundary obeys the superradiance stability inequality
   omega - m Omega > 0; near-critical Sigma_3 is possible without destabilizing
   the heavy modes only in a narrow Omega window.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "boundary_time_source_hint"


def upper_vec(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    return np.array([mat[i, j] for i in range(n) for j in range(i, n)])


def so_generator(n: int, a: int, b: int) -> np.ndarray:
    gen = np.zeros((n, n))
    gen[a, b] = 1.0
    gen[b, a] = -1.0
    return gen


def audit_54_orbit() -> dict[str, float | int]:
    n = 10
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    tangents = []
    for a, b in itertools.combinations(range(n), 2):
        gen = so_generator(n, a, b)
        delta = gen @ s0 - s0 @ gen
        tangents.append(upper_vec(delta))
    mat = np.vstack(tangents)
    rank = int(np.linalg.matrix_rank(mat, tol=1e-10))
    dim_sym_traceless = n * (n + 1) // 2 - 1
    return {
        "ambient_dim": dim_sym_traceless,
        "orbit_rank": rank,
        "radial_singlets_retained": 1,
        "normal_modes_removed": dim_sym_traceless - rank - 1,
        "stabilizer_dim": 45 - rank,
    }


def wedge_sign(sorted_old: tuple[int, ...], unsorted_new: list[int]) -> tuple[tuple[int, ...], int]:
    if len(set(unsorted_new)) < len(unsorted_new):
        return tuple(), 0
    inversions = 0
    for i in range(len(unsorted_new)):
        for j in range(i + 1, len(unsorted_new)):
            if unsorted_new[i] > unsorted_new[j]:
                inversions += 1
    return tuple(sorted(unsorted_new)), -1 if inversions % 2 else 1


def act_generator_on_form(
    form: tuple[int, ...],
    gen_pair: tuple[int, int],
    basis_index: dict[tuple[int, ...], int],
) -> np.ndarray:
    a, b = gen_pair
    vec = np.zeros(len(basis_index))
    for slot, idx in enumerate(form):
        replacement = None
        sign = 1
        if idx == a:
            replacement = b
            sign = 1
        elif idx == b:
            replacement = a
            sign = -1
        if replacement is None:
            continue
        new_form = list(form)
        new_form[slot] = replacement
        sorted_form, wedge_sgn = wedge_sign(form, new_form)
        if wedge_sgn:
            vec[basis_index[sorted_form]] += sign * wedge_sgn
    return vec


def audit_210_orbit() -> dict[str, float | int | dict[str, int]]:
    n = 10
    basis = list(itertools.combinations(range(n), 4))
    basis_index = {item: i for i, item in enumerate(basis)}
    omega = (6, 7, 8, 9)
    tangents = []
    nonzero_generators = 0
    for pair in itertools.combinations(range(n), 2):
        vec = act_generator_on_form(omega, pair, basis_index)
        if float(np.linalg.norm(vec)) > 1e-12:
            nonzero_generators += 1
        tangents.append(vec)
    mat = np.vstack(tangents)
    rank = int(np.linalg.matrix_rank(mat, tol=1e-10))
    dim_4form = math.comb(n, 4)
    return {
        "ambient_dim": dim_4form,
        "orbit_rank": rank,
        "nonzero_generators": nonzero_generators,
        "radial_singlets_retained": 1,
        "normal_modes_removed": dim_4form - rank - 1,
        "stabilizer_dim": 45 - rank,
        "D210_twoform_counts": {"-1": 3, "0": 39, "+1": 3},
    }


def audit_rindler_lapse() -> tuple[list[dict[str, float]], dict[str, float]]:
    ell = 1.0
    samples = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 5e-1]
    rows = []
    for n in samples:
        lapse = math.tanh(n / ell)
        rindler = n / ell
        rel_error = abs(lapse / rindler - 1.0)
        cubic_prediction = (n / ell) ** 2 / 3.0
        rows.append(
            {
                "n_over_ell": n / ell,
                "N_tanh": lapse,
                "N_Rindler_linear": rindler,
                "relative_error": rel_error,
                "cubic_error_estimate": cubic_prediction,
            }
        )
    small = [row for row in rows if row["n_over_ell"] <= 1e-2]
    summary = {
        "max_relative_error_n_le_1e_minus_2": max(row["relative_error"] for row in small),
        "rindler_acceleration_a": 1.0 / ell,
        "series": "tanh(n/ell)=n/ell-(n/ell)^3/3+O(n^5)",
    }
    return rows, summary


def audit_rotating_boundary() -> tuple[list[dict[str, float | str]], dict[str, float | str]]:
    modes = {
        "Sigma3_intermediate": 0.13001566778681514,
        "Sigma8_heavy": 0.5496247194531114,
        "X622_lifted": 1.0,
        "conormal_pair": 1.0,
    }
    omegas = [0.0, 0.03, 0.05, 0.08, 0.10, 0.12, 0.13, 0.15, 0.20]
    rows = []
    for name, omega0 in modes.items():
        for m in [1, 2]:
            critical = omega0 / m
            for Omega in omegas:
                eco = omega0 - m * Omega
                rows.append(
                    {
                        "mode": name,
                        "mass_over_MG": omega0,
                        "azimuthal_m": m,
                        "Omega_over_MG": Omega,
                        "E_corotating_over_MG": eco,
                        "status": "stable" if eco > 0 else "superradiant_or_tachyonic",
                        "critical_Omega_over_MG": critical,
                    }
                )
    sigma_m1 = [
        row
        for row in rows
        if row["mode"] == "Sigma3_intermediate" and row["azimuthal_m"] == 1
    ]
    stable = [row["Omega_over_MG"] for row in sigma_m1 if row["status"] == "stable"]
    unstable = [row["Omega_over_MG"] for row in sigma_m1 if row["status"] != "stable"]
    critical = modes["Sigma3_intermediate"]
    summary = {
        "Sigma3_m1_critical_Omega_over_MG": critical,
        "largest_sampled_stable_Omega_for_Sigma3_m1": max(stable),
        "smallest_sampled_unstable_Omega_for_Sigma3_m1": min(unstable),
        "critical_is_bracketed_by_scan": max(stable) < critical < min(unstable),
        "X622_m1_margin_at_Omega_0p12": modes["X622_lifted"] - 0.12,
        "conormal_m1_margin_at_Omega_0p12": modes["conormal_pair"] - 0.12,
        "interpretation": (
            "A rotating boundary can stress Sigma3 near criticality around "
            "Omega/MG ~= 0.13 for m=1 while keeping X622 and conormal modes "
            "stable; this is a candidate mechanism, not yet a derivation."
        ),
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lapse_rows, lapse_summary = audit_rindler_lapse()
    rot_rows, rot_summary = audit_rotating_boundary()
    source_54 = audit_54_orbit()
    source_210 = audit_210_orbit()

    write_csv(OUT / "rindler_lapse_scan.csv", lapse_rows)
    write_csv(OUT / "rotating_boundary_stability.csv", rot_rows)

    summary = {
        "note": "No web lookup used. Boundary-time/source hint audit.",
        "rindler_lapse": lapse_summary,
        "source_54": source_54,
        "source_210": source_210,
        "rotating_boundary": rot_summary,
        "passes": {
            "rindler_limit_error_lt_1e_minus_4": lapse_summary[
                "max_relative_error_n_le_1e_minus_2"
            ]
            < 1e-4,
            "source_54_orbit_rank_24": source_54["orbit_rank"] == 24,
            "source_210_orbit_rank_24": source_210["orbit_rank"] == 24,
            "normal_mode_counts_match_A5": (
                source_54["normal_modes_removed"] == 29
            and source_210["normal_modes_removed"] == 185
            ),
            "rotating_window_brackets_sigma3_m1_critical": rot_summary[
                "critical_is_bracketed_by_scan"
            ],
            "heavy_modes_stable_at_sigma3_near_critical": (
                rot_summary["X622_m1_margin_at_Omega_0p12"] > 0.8
                and rot_summary["conormal_m1_margin_at_Omega_0p12"] > 0.8
            ),
        },
        "verdict": (
            "The sketches can be formalized as a boundary-time constrained-source "
            "ansatz: lapse degeneracy gives a Rindler interface, hard-core source "
            "constraints remove 54/210 normal modes, and rotating-boundary "
            "stability can make Sigma3 near-critical without releasing X622. "
            "This is a concrete A5-A6 research route, not yet a microscopic proof."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    report = [
        "# Boundary-time constrained-source hint audit",
        "",
        "No web lookup used.",
        "",
        "## Rindler layer",
        "",
        "Use ds^2=-N(n)^2 dt^2+dn^2+dy^2+dz^2 with N(n)=tanh(n/ell).",
        "Near the interface n=0, N=n/ell-(n/ell)^3/3+O(n^5), giving the Rindler form.",
        f"Max relative error for n/ell <= 1e-2: {lapse_summary['max_relative_error_n_le_1e_minus_2']:.6e}.",
        "",
        "## Constrained sources",
        "",
        f"54 source: ambient={source_54['ambient_dim']}, orbit rank={source_54['orbit_rank']}, normal removed={source_54['normal_modes_removed']}.",
        f"210 source: ambient={source_210['ambient_dim']}, orbit rank={source_210['orbit_rank']}, normal removed={source_210['normal_modes_removed']}.",
        "",
        "## Rotating boundary",
        "",
        "Stability condition: E_co/MG = omega/MG - m Omega/MG > 0.",
        f"Sigma3 m=1 critical Omega/MG = {rot_summary['Sigma3_m1_critical_Omega_over_MG']:.6f}.",
        f"Largest sampled stable Omega/MG = {rot_summary['largest_sampled_stable_Omega_for_Sigma3_m1']:.6f}.",
        f"Smallest sampled unstable Omega/MG = {rot_summary['smallest_sampled_unstable_Omega_for_Sigma3_m1']:.6f}.",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report))
    print(json.dumps(summary["passes"], sort_keys=True))


if __name__ == "__main__":
    main()
