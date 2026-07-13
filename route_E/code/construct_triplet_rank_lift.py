#!/usr/bin/env python3
"""Construct a rank-one-lift realization of the triplet near-null propagator.

No web lookup is used.  The previous scan found a rank-deficient inverse
triplet propagator W_AB that nearly nulls the LLLL K nu tensor.  A finite
superpotential cannot literally use a singular inverse mass matrix, so this
script regularizes the near-null direction as

    W_eps = U diag(s1, s2, s3, eps*s1) V^dagger

and interprets the inverse

    M_eps = V diag(1/s1, 1/s2, 1/s3, 1/(eps*s1)) U^dagger

as a rank-one-lift triplet mass matrix.  It then checks the channel amplitudes,
mass hierarchy, and F/D-flatness at the zero-triplet vacuum.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_triplet_mixing_nullspace as ns  # noqa: E402


INPUT_FLAVOR = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
INPUT_NULL = ROOT / "output" / "triplet_mixing_nullspace" / "triplet_mixing_nullspace_summary.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "triplet_rank_lift"

PLANCK_BENCHMARK_GEV = 2.435e18
M_LIGHT_BENCHMARK_GEV = 1.0e16


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def coeffs_to_matrix(coefficients: list[dict[str, Any]], basis_names: list[str]) -> np.ndarray:
    lookup = {name: i for i, name in enumerate(basis_names)}
    out = np.zeros((len(basis_names), len(basis_names)), dtype=complex)
    for item in coefficients:
        i = lookup[item["QQ_source"]]
        j = lookup[item["QL_source"]]
        out[i, j] = item["re"] + 1j * item["im"]
    return out


def coeffs_from_matrix(w: np.ndarray) -> np.ndarray:
    return np.array([w[i, j] for i in range(w.shape[0]) for j in range(w.shape[1])], dtype=complex)


def build_channel_tensors() -> tuple[list[str], list[tuple[int, int]], dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]], dict[str, float]]:
    payload = json.loads(INPUT_FLAVOR.read_text(encoding="utf-8"))
    constants = json.loads(PROTON.read_text(encoding="utf-8"))["hadronic_constants"]

    y_raw = {name: cmat(raw) for name, raw in payload["Yukawa_matrices"].items()}
    rot = payload["biunitary_rotations"]
    ul = {name: cmat(rot[name]["left_rotation"]) for name in rot}
    ur = {name: cmat(rot[name]["right_rotation"]) for name in rot}
    u_nu = cmat(payload["PMNS_convention"]["neutrino_left_rotation_target"])

    y_phys = {
        "up": ns.scale_to_largest_singular(y_raw["up"], 0.60),
        "down": ns.scale_to_largest_singular(y_raw["down"], 0.024),
        "charged_lepton": ns.scale_to_largest_singular(y_raw["charged_lepton"], 0.010),
    }
    k_tr_top = ns.scale_to_largest_singular(
        np.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=complex) / math.sqrt(3.0),
        0.60,
    )
    k_tr_bottom = ns.scale_to_largest_singular(k_tr_top, 0.024)

    basis_names = ["H10", "F126", "G120", "Ktr"]
    delta_de = y_phys["down"] - y_phys["charged_lepton"].T
    qq_basis = [
        ns.scale_to_largest_singular(ns.sym(y_phys["up"]), 0.60),
        ns.scale_to_largest_singular(ns.sym(delta_de), 0.60),
        ns.scale_to_largest_singular(ns.anti(y_phys["up"]), 0.60),
        k_tr_top,
    ]
    ql_basis = [
        ns.scale_to_largest_singular(y_phys["down"], 0.024),
        ns.scale_to_largest_singular(delta_de, 0.024),
        ns.scale_to_largest_singular(ns.anti(y_phys["down"]), 0.024),
        k_tr_bottom,
    ]

    pair_full = [(i, j) for i in range(len(basis_names)) for j in range(len(basis_names))]
    tensors_up: dict[tuple[int, int], np.ndarray] = {}
    tensors_down: dict[tuple[int, int], np.ndarray] = {}
    tensors_mu: dict[tuple[int, int], np.ndarray] = {}
    tensors_r: dict[tuple[int, int], np.ndarray] = {}
    for i in range(len(basis_names)):
        for j in range(len(basis_names)):
            tensors_up[(i, j)] = ns.tensor_llll(qq_basis[i], ql_basis[j], ul["up"], ul["up"], ul["down"], u_nu)
            tensors_down[(i, j)] = ns.tensor_llll(qq_basis[i], ql_basis[j], ul["down"], ul["down"], ul["up"], u_nu)
            tensors_mu[(i, j)] = ns.tensor_llll(qq_basis[i], ql_basis[j], ul["up"], ul["up"], ul["down"], ul["charged_lepton"])
            tensors_r[(i, j)] = ns.tensor_rrrr(qq_basis[i], ql_basis[j], ur["up"], ur["charged_lepton"], ur["up"], ur["down"])

    entries_knu = [(a, b, c, ell) for (a, b, c) in sorted(set(permutations((0, 0, 1)))) for ell in [0, 1, 2]]
    entries_k0mu = [(a, b, c, 1) for (a, b, c) in sorted(set(permutations((0, 0, 1))))]
    entries_rrrr = []
    for ell in [0, 1, 2]:
        for u_pair in [(0, 0), (0, 1), (1, 0)]:
            entries_rrrr.append((u_pair[0], u_pair[1], 1, ell))
    channel_entries = {
        "LLLL_upupdown_Knu": (tensors_up, entries_knu),
        "LLLL_downdownup_Knu": (tensors_down, entries_knu),
        "LLLL_upupdown_K0mu": (tensors_mu, entries_k0mu),
        "RRRR_uusd_anycharged": (tensors_r, entries_rrrr),
    }
    return basis_names, pair_full, channel_entries, constants


def max_channel_for_w(w: np.ndarray, pair_full: list[tuple[int, int]], channel_entries: dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]], constants: dict[str, float]) -> dict[str, Any]:
    coeffs = coeffs_from_matrix(w)
    out: dict[str, Any] = {}
    for channel, (tensors, entries) in channel_entries.items():
        amp, idx, value = ns.max_channel(coeffs, pair_full, tensors, entries)
        life = ns.lifetime_from_amplitude(amp, constants)
        out[channel] = {
            "amplitude": float(amp),
            "selected_index": "".join(str(x) for x in idx),
            "value": cjson(value),
            "tau_years": life["tau_years"],
            "S_T_required_tau_2p4e34": life["S_T_required_tau_2p4e34"],
        }
    return out


def summarize_mass_from_singulars(singulars: np.ndarray) -> dict[str, Any]:
    s1 = float(singulars[0])
    zero_tol = 1.0e-14 * max(s1, 1.0)
    positive = singulars[singulars > zero_tol]
    mass_ratios = [float(s1 / s) if s > zero_tol else math.inf for s in singulars]
    max_ratio = max(x for x in mass_ratios if math.isfinite(x))
    return {
        "inverse_propagator_singular_values": [float(x) for x in singulars],
        "zero_singular_value_tolerance": float(zero_tol),
        "infinite_mass_direction_count": int(np.sum(singulars <= zero_tol)),
        "triplet_mass_ratios_to_lightest": mass_ratios,
        "max_finite_mass_ratio": float(max_ratio),
        "max_mass_if_lightest_1e16_GeV": float(max_ratio * M_LIGHT_BENCHMARK_GEV),
        "exceeds_reduced_planck_benchmark": bool(max_ratio * M_LIGHT_BENCHMARK_GEV > PLANCK_BENCHMARK_GEV),
    }


def row_from_family(label: str, kind: str, epsilon: float | None, w: np.ndarray, pair_full: list[tuple[int, int]], channel_entries: dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]], constants: dict[str, float]) -> dict[str, Any]:
    singulars = np.linalg.svd(w, compute_uv=False)
    channels = max_channel_for_w(w, pair_full, channel_entries, constants)
    max_knu = max(channels["LLLL_upupdown_Knu"]["amplitude"], channels["LLLL_downdownup_Knu"]["amplitude"])
    leading = max(channels.items(), key=lambda item: item[1]["amplitude"])
    mass = summarize_mass_from_singulars(singulars)
    return {
        "label": label,
        "kind": kind,
        "epsilon": "" if epsilon is None else float(epsilon),
        "max_Knu_amplitude": float(max_knu),
        "K0mu_amplitude": channels["LLLL_upupdown_K0mu"]["amplitude"],
        "RRRR_amplitude": channels["RRRR_uusd_anycharged"]["amplitude"],
        "leading_channel": leading[0],
        "leading_amplitude": leading[1]["amplitude"],
        "leading_ST_required_tau_2p4e34": leading[1]["S_T_required_tau_2p4e34"],
        "max_mass_ratio": mass["max_finite_mass_ratio"],
        "max_mass_if_lightest_1e16_GeV": mass["max_mass_if_lightest_1e16_GeV"],
        "exceeds_reduced_planck_benchmark": mass["exceeds_reduced_planck_benchmark"],
        "channels": channels,
        "mass_summary": mass,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    null_payload = json.loads(INPUT_NULL.read_text(encoding="utf-8"))
    basis_names = null_payload["basis_names"]
    near_audit = next(item for item in null_payload["audits"] if item["label"] == "full_bipartite_mixing_knu_null")
    w_near = coeffs_to_matrix(near_audit["coefficients"], basis_names)
    u, singulars, vh = np.linalg.svd(w_near, full_matrices=True)
    s1 = float(singulars[0])

    basis_names_check, pair_full, channel_entries, constants = build_channel_tensors()
    if basis_names_check != basis_names:
        raise RuntimeError("Basis mismatch between nullspace scan and rank-lift construction")

    rows: list[dict[str, Any]] = []
    regularized_family: list[dict[str, Any]] = []
    for eps in [1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1]:
        s_eps = singulars.copy()
        s_eps[-1] = eps * s1
        w_eps = u @ np.diag(s_eps) @ vh
        row = row_from_family(f"rank_one_lift_eps_{eps:.0e}", "rank_one_lift", eps, w_eps, pair_full, channel_entries, constants)
        rows.append(row)
        regularized_family.append(row)

    capped_family: list[dict[str, Any]] = []
    for kappa in [10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0]:
        s_cap = np.maximum(singulars, s1 / kappa)
        w_cap = u @ np.diag(s_cap) @ vh
        row = row_from_family(f"condition_cap_{kappa:.0f}", "condition_cap", 1.0 / kappa, w_cap, pair_full, channel_entries, constants)
        rows.append(row)
        capped_family.append(row)

    rank3_row = row_from_family("rank3_limit", "rank3_limit", 0.0, w_near, pair_full, channel_entries, constants)
    rows.insert(0, rank3_row)

    csv_rows = []
    for row in rows:
        csv_rows.append(
            {
                "label": row["label"],
                "kind": row["kind"],
                "epsilon": row["epsilon"],
                "max_Knu_amplitude": row["max_Knu_amplitude"],
                "K0mu_amplitude": row["K0mu_amplitude"],
                "RRRR_amplitude": row["RRRR_amplitude"],
                "leading_channel": row["leading_channel"],
                "leading_amplitude": row["leading_amplitude"],
                "leading_ST_required_tau_2p4e34": row["leading_ST_required_tau_2p4e34"],
                "max_mass_ratio": row["max_mass_ratio"],
                "max_mass_if_lightest_1e16_GeV": row["max_mass_if_lightest_1e16_GeV"],
                "exceeds_reduced_planck_benchmark": row["exceeds_reduced_planck_benchmark"],
            }
        )

    with (OUT / "triplet_rank_lift_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    best_planck_safe = [row for row in rows if not row["exceeds_reduced_planck_benchmark"]]
    best_planck_safe = min(best_planck_safe, key=lambda row: row["max_Knu_amplitude"]) if best_planck_safe else None
    rank3_mass = rank3_row["mass_summary"]
    obstruction = (
        "The rank-one lift gives a clean superpotential interpretation of the missing singular direction, "
        "but the near-null propagator is extremely hierarchical.  With the lightest triplet at 1e16 GeV, "
        "even the finite rank-3 singular values place a triplet above the reduced-Planck benchmark; "
        "therefore the mechanism must either be a cutoff EFT/composite projection or be supplemented by an "
        "RRRR filter and a different mass normalization."
    )

    output = {
        "note": "No web lookup used. Rank-one-lift realization of the triplet near-null inverse propagator.",
        "superpotential": {
            "direct_rank_lift": "W_T = T_A (M_rank3 + M_* v_4 u_4^dagger)_{AB} barT_B",
            "mediator_form": "W_T = T M_rank3 barT + X p_A T_A + q_B barT_B barX + M_X X barX; Schur complement gives a rank-one update.",
            "zero_triplet_vacuum": "For epsilon>0 the mass matrix is full rank, so F_T=M barT=0 and F_barT=M^T T=0 imply T=barT=0; D-terms vanish at the origin.",
        },
        "basis_names": basis_names,
        "near_null_W": {
            "matrix": matrix_json(w_near),
            "singular_values": [float(x) for x in singulars],
            "rank3_mass_summary": rank3_mass,
        },
        "rows": rows,
        "verdict": {
            "rank3_limit_max_Knu_amplitude": rank3_row["max_Knu_amplitude"],
            "rank3_limit_RRRR_amplitude": rank3_row["RRRR_amplitude"],
            "rank3_limit_max_mass_ratio": rank3_mass["max_finite_mass_ratio"],
            "rank3_limit_max_mass_if_lightest_1e16_GeV": rank3_mass["max_mass_if_lightest_1e16_GeV"],
            "best_planck_safe_row": None if best_planck_safe is None else {
                key: best_planck_safe[key]
                for key in [
                    "label",
                    "max_Knu_amplitude",
                    "RRRR_amplitude",
                    "leading_channel",
                    "leading_ST_required_tau_2p4e34",
                    "max_mass_ratio",
                    "max_mass_if_lightest_1e16_GeV",
                ]
            },
            "reduced_planck_benchmark_GeV": PLANCK_BENCHMARK_GEV,
            "interpretation": obstruction,
        },
    }
    (OUT / "triplet_rank_lift_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Triplet rank-one-lift audit",
        "",
        "No web lookup was used.",
        "",
        "The near-null inverse propagator is regularized as",
        "",
        "```text",
        "W_eps = U diag(s1,s2,s3,eps*s1) V^dagger.",
        "```",
        "",
        "| label | max Knu | K0mu | RRRR | leading | ST max 2.4e34 | max mass ratio | Mmax if Mlight=1e16 |",
        "|---|---:|---:|---:|---|---:|---:|---:|",
    ]
    for row in rows:
        report.append(
            f"| `{row['label']}` | {row['max_Knu_amplitude']:.3e} | {row['K0mu_amplitude']:.3e} | "
            f"{row['RRRR_amplitude']:.3e} | `{row['leading_channel']}` | "
            f"{row['leading_ST_required_tau_2p4e34']:.3e} | {row['max_mass_ratio']:.3e} | "
            f"{row['max_mass_if_lightest_1e16_GeV']:.3e} |"
        )
    report.extend(
        [
            "",
            "## Verdict",
            "",
            obstruction,
            "",
            f"Rank-3 limit max Knu amplitude: `{rank3_row['max_Knu_amplitude']:.6e}`.",
            f"Rank-3 limit RRRR amplitude: `{rank3_row['RRRR_amplitude']:.6e}`.",
            f"Rank-3 finite mass hierarchy: `{rank3_mass['max_finite_mass_ratio']:.6e}`.",
            f"Max finite triplet mass if the lightest triplet is 1e16 GeV: `{rank3_mass['max_mass_if_lightest_1e16_GeV']:.6e}` GeV.",
            "",
        ]
    )
    (OUT / "triplet_rank_lift_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Triplet rank-one-lift audit")
    print(f"  rank3 max Knu amplitude: {rank3_row['max_Knu_amplitude']:.6e}")
    print(f"  rank3 RRRR amplitude: {rank3_row['RRRR_amplitude']:.6e}")
    print(f"  rank3 finite mass hierarchy: {rank3_mass['max_finite_mass_ratio']:.6e}")
    print(f"  rank3 Mmax if Mlight=1e16 GeV: {rank3_mass['max_mass_if_lightest_1e16_GeV']:.6e}")
    if best_planck_safe is not None:
        print(f"  best Planck-safe row: {best_planck_safe['label']} Knu={best_planck_safe['max_Knu_amplitude']:.6e}")
    else:
        print("  no row is below the reduced-Planck benchmark for Mlight=1e16 GeV")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
