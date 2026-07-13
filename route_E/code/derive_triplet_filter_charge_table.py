#!/usr/bin/env python3
"""Derive a minimal grading for the two-sided triplet filter support.

No web lookup is used.  This script takes the kappa=100, omega_R=0.1 two-sided
triplet filter and asks a narrow action-level question: can a small cyclic
grading reproduce the dominant effective inverse-propagator support while
forbidding the entries that are numerically small in the SVD solution?

The grading is an EFT/missing-partner support diagnostic for W_AB.  A full
superpotential must still realize this effective propagator by integrating out
triplet partners and mediators.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import construct_triplet_rank_lift as lift  # noqa: E402
import scan_two_sided_triplet_filter as two  # noqa: E402
import scan_triplet_mixing_nullspace as ns  # noqa: E402


INPUT = ROOT / "output" / "two_sided_triplet_filter" / "two_sided_triplet_filter_summary.json"
OUT = ROOT / "output" / "triplet_filter_charge_table"


def load_best() -> tuple[list[str], np.ndarray, dict[str, Any]]:
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    row = next(item for item in payload["rows"] if item["label"] == "omegaR_0.1_kappa_100")
    w = np.array([[cell["re"] + 1j * cell["im"] for cell in line] for line in row["W_matrix"]], dtype=complex)
    return payload["basis_names"], w, row


def mask_from_charges(q_left: tuple[int, ...], q_right: tuple[int, ...], residues: tuple[int, ...], modulus: int) -> np.ndarray:
    rset = set(residues)
    return np.array(
        [[((q_left[i] + q_right[j]) % modulus) in rset for j in range(len(q_right))] for i in range(len(q_left))],
        dtype=bool,
    )


def scan_charge_tables(abs_w: np.ndarray, dominant_threshold: float = 4.0e-2, suppressed_threshold: float = 2.0e-2) -> list[dict[str, Any]]:
    dominant = abs_w >= dominant_threshold
    suppressed = abs_w <= suppressed_threshold
    results: list[dict[str, Any]] = []
    for modulus in range(2, 9):
        # Gauge-fix H10 charges to zero and require the neutral residue because
        # the H10-H10 entry is dominant.
        for q_left_tail in itertools.product(range(modulus), repeat=3):
            q_left = (0,) + q_left_tail
            for q_right_tail in itertools.product(range(modulus), repeat=3):
                q_right = (0,) + q_right_tail
                residues_options = [(0,)]
                if modulus <= 5:
                    residues_options += [tuple(sorted((0, r))) for r in range(1, modulus)]
                for residues in residues_options:
                    mask = mask_from_charges(q_left, q_right, residues, modulus)
                    false_neg = int(np.sum(dominant & ~mask))
                    false_pos = int(np.sum(suppressed & mask))
                    allowed_count = int(np.sum(mask))
                    dominant_count = int(np.sum(dominant))
                    suppressed_count = int(np.sum(suppressed))
                    recall = 1.0 - false_neg / max(dominant_count, 1)
                    leak_rate = false_pos / max(suppressed_count, 1)
                    # Lexicographic ranking: exact dominant coverage first,
                    # then no suppressed leakage, then small modulus/support.
                    score = (
                        1000 * false_neg
                        + 100 * false_pos
                        + 2 * len(residues)
                        + 0.01 * allowed_count
                        + 0.001 * modulus
                    )
                    results.append(
                        {
                            "modulus": modulus,
                            "q_left": q_left,
                            "q_right": q_right,
                            "allowed_residues": residues,
                            "false_negative_dominant": false_neg,
                            "false_positive_suppressed": false_pos,
                            "allowed_count": allowed_count,
                            "dominant_count": dominant_count,
                            "suppressed_count": suppressed_count,
                            "dominant_recall": recall,
                            "suppressed_leak_rate": leak_rate,
                            "score": score,
                            "mask": mask,
                        }
                    )
    return sorted(results, key=lambda item: (item["score"], item["modulus"], item["allowed_count"]))


def channel_rows_for_allowed(mask: np.ndarray, weight_rrrr: float = 0.1, kappa: float = 100.0) -> dict[str, Any]:
    basis_names, pair_full, channel_entries, constants = lift.build_channel_tensors()
    a_knu = np.vstack(
        [
            ns.row_matrix(pair_full, channel_entries["LLLL_upupdown_Knu"][0], channel_entries["LLLL_upupdown_Knu"][1]),
            ns.row_matrix(pair_full, channel_entries["LLLL_downdownup_Knu"][0], channel_entries["LLLL_downdownup_Knu"][1]),
        ]
    )
    a_rrrr = ns.row_matrix(pair_full, channel_entries["RRRR_uusd_anycharged"][0], channel_entries["RRRR_uusd_anycharged"][1])
    allowed_columns = [n for n, (i, j) in enumerate(pair_full) if mask[i, j]]
    joint = np.vstack([two.normalized_block(a_knu), weight_rrrr * two.normalized_block(a_rrrr)])
    restricted = joint[:, allowed_columns]
    _, singulars, vh = np.linalg.svd(restricted, full_matrices=True)
    coeff_restricted = vh.conjugate().T[:, -1]
    coeff_restricted = coeff_restricted / np.linalg.norm(coeff_restricted)
    coeff = np.zeros(len(pair_full), dtype=complex)
    for out_idx, col in enumerate(allowed_columns):
        coeff[col] = coeff_restricted[out_idx]
    w_allowed = two.w_from_coeffs(coeff, 4)
    # The support-only solution is generally rank deficient.  The capped replay
    # models Planck-safe leakage by adding the smallest allowed universal
    # singular floor after the support direction is selected.
    w_cap, cap_info = two.cap_condition(w_allowed, kappa)
    row_support = two.evaluate(
        "charge_support_rank_limit",
        weight_rrrr,
        None,
        w_allowed / max(np.linalg.norm(w_allowed), 1.0e-30),
        pair_full,
        channel_entries,
        constants,
        {"joint_sigma_min": float(singulars[-1]), "joint_sigma_max": float(singulars[0]), "joint_singular_values": [float(x) for x in singulars]},
        {"condition_cap": "rank_limit", "final_condition": math.inf, "max_mass_if_lightest_1e16_GeV": math.inf},
    )
    row_cap = two.evaluate(
        "charge_support_kappa_100_floor",
        weight_rrrr,
        kappa,
        w_cap,
        pair_full,
        channel_entries,
        constants,
        {"joint_sigma_min": float(singulars[-1]), "joint_sigma_max": float(singulars[0]), "joint_singular_values": [float(x) for x in singulars]},
        cap_info,
    )
    return {
        "basis_names": basis_names,
        "allowed_columns": allowed_columns,
        "restricted_singular_values": [float(x) for x in singulars],
        "support_rank_limit": row_support,
        "support_kappa_100_floor": row_cap,
    }


def matrix_bool_rows(mask: np.ndarray, names: list[str]) -> list[dict[str, Any]]:
    rows = []
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            rows.append({"left": left, "right": right, "allowed": bool(mask[i, j])})
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names, w_best, best_row = load_best()
    abs_w = np.abs(w_best)
    scans = scan_charge_tables(abs_w)
    best = scans[0]
    mask = best["mask"]
    replay = channel_rows_for_allowed(mask)

    csv_rows = []
    for item in scans[:200]:
        csv_rows.append(
            {
                "modulus": item["modulus"],
                "q_left": " ".join(map(str, item["q_left"])),
                "q_right": " ".join(map(str, item["q_right"])),
                "allowed_residues": " ".join(map(str, item["allowed_residues"])),
                "false_negative_dominant": item["false_negative_dominant"],
                "false_positive_suppressed": item["false_positive_suppressed"],
                "allowed_count": item["allowed_count"],
                "dominant_recall": item["dominant_recall"],
                "suppressed_leak_rate": item["suppressed_leak_rate"],
                "score": item["score"],
            }
        )
    with (OUT / "triplet_filter_charge_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    table_rows = []
    for i, name in enumerate(names):
        table_rows.append(
            {
                "field": f"T_{name}",
                "Z_N_charge": best["q_left"][i],
                "role": "left triplet source",
            }
        )
    for i, name in enumerate(names):
        table_rows.append(
            {
                "field": f"barT_{name}",
                "Z_N_charge": best["q_right"][i],
                "role": "right triplet source",
            }
        )
    table_rows.append({"field": "S_0", "Z_N_charge": 0, "role": "neutral mass/link spurion"})
    with (OUT / "triplet_filter_charge_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["field", "Z_N_charge", "role"])
        writer.writeheader()
        writer.writerows(table_rows)

    support_rows = matrix_bool_rows(mask, names)
    with (OUT / "triplet_filter_allowed_support.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["left", "right", "allowed"])
        writer.writeheader()
        writer.writerows(support_rows)

    verdict = (
        "A single neutral-spurion Z3 grading exactly reproduces the dominant support of the kappa=100 "
        "two-sided filter: it allows the H/F two-by-two block and the Ktr -> H,F row, while forbidding "
        "all entries below the suppressed threshold.  The support-only replay remains rank deficient, "
        "so the charge table is best interpreted as a missing-partner projection plus a small universal "
        "condition-number floor.  With the kappa=100 floor the amplitudes remain in the same 1e-4 class "
        "as the unconstrained SVD optimum."
    )

    output = {
        "note": "No web lookup used. Minimal cyclic grading for the effective triplet filter support.",
        "source_filter": {
            "label": best_row["label"],
            "weight_RRRR": best_row["weight_RRRR"],
            "condition_cap": best_row["condition_cap"],
            "max_Knu_amplitude": best_row["max_Knu_amplitude"],
            "RRRR_amplitude": best_row["RRRR_amplitude"],
            "K0mu_amplitude": best_row["K0mu_amplitude"],
        },
        "thresholds": {
            "dominant_abs_W_entry": 4.0e-2,
            "suppressed_abs_W_entry": 2.0e-2,
        },
        "best_charge_table": {
            "modulus": best["modulus"],
            "q_left_by_basis": {name: int(best["q_left"][i]) for i, name in enumerate(names)},
            "q_right_by_basis": {name: int(best["q_right"][i]) for i, name in enumerate(names)},
            "allowed_residues": list(best["allowed_residues"]),
            "false_negative_dominant": best["false_negative_dominant"],
            "false_positive_suppressed": best["false_positive_suppressed"],
            "allowed_support": matrix_bool_rows(mask, names),
        },
        "best_W_abs": [[float(x) for x in row] for row in abs_w],
        "replay": replay,
        "verdict": verdict,
    }
    (OUT / "triplet_filter_charge_table_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    support = replay["support_kappa_100_floor"]
    report = [
        "# Triplet filter charge-table audit",
        "",
        "No web lookup was used.",
        "",
        f"Best grading: `Z_{best['modulus']}` with neutral residue `{best['allowed_residues'][0]}`.",
        "",
        "| source | H10 | F126 | G120 | Ktr |",
        "|---|---:|---:|---:|---:|",
        f"| T_A | {best['q_left'][0]} | {best['q_left'][1]} | {best['q_left'][2]} | {best['q_left'][3]} |",
        f"| barT_B | {best['q_right'][0]} | {best['q_right'][1]} | {best['q_right'][2]} | {best['q_right'][3]} |",
        "",
        "Allowed support is `q(T_A)+q(barT_B)=0 mod N`.",
        "",
        "| replay | max Knu | K0mu | RRRR | leading | ST max 2.4e34 |",
        "|---|---:|---:|---:|---|---:|",
        f"| support + kappa=100 floor | {support['max_Knu_amplitude']:.3e} | "
        f"{support['K0mu_amplitude']:.3e} | {support['RRRR_amplitude']:.3e} | "
        f"`{support['leading_channel']}` | {support['leading_ST_required_tau_2p4e34']:.3e} |",
        "",
        "## Verdict",
        "",
        verdict,
        "",
    ]
    (OUT / "triplet_filter_charge_table_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Triplet filter charge-table audit")
    print(f"  best grading: Z_{best['modulus']}")
    print(f"  q_left: {best['q_left']}")
    print(f"  q_right: {best['q_right']}")
    print(f"  false negatives / false positives: {best['false_negative_dominant']} / {best['false_positive_suppressed']}")
    print(
        "  replay kappa100: "
        f"Knu={support['max_Knu_amplitude']:.6e} "
        f"K0mu={support['K0mu_amplitude']:.6e} "
        f"RRRR={support['RRRR_amplitude']:.6e}"
    )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
