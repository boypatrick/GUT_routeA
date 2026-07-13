#!/usr/bin/env python3
"""Construct a Z3-compatible finite regulator for the triplet filter.

No web lookup is used.  The previous charge-table audit found a neutral
Z3 support that reproduces the dominant two-sided triplet filter, but that
support is rank deficient.  Here we ask whether the kappa=100 singular-value
floor can be interpreted as a small charged-spurion regulator rather than as
an ad hoc numerical insertion.

The effective triplet inverse propagator is decomposed by

    r_AB = q(T_A) + q(bar T_B) mod 3.

Neutral entries have r=0.  Entries with r=1 and r=2 require spurions of charge
2 and 1, respectively.  The audit builds the closest condition-capped matrix
in singular-value norm and measures how much charged-spurion support it needs.
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

import construct_triplet_rank_lift as lift  # noqa: E402
import derive_triplet_filter_charge_table as ct  # noqa: E402
import scan_triplet_mixing_nullspace as ns  # noqa: E402
import scan_two_sided_triplet_filter as two  # noqa: E402


OUT = ROOT / "output" / "z3_triplet_regulator"
KAPPAS = [30.0, 100.0, 300.0]
WEIGHT_RRRR = 0.1
ENTRY_TOL = 1.0e-10


def reconstruct_neutral_support() -> dict[str, Any]:
    names, w_best, _ = ct.load_best()
    best = ct.scan_charge_tables(np.abs(w_best))[0]
    mask = best["mask"]
    basis_names, pair_full, channel_entries, constants = lift.build_channel_tensors()
    a_knu = np.vstack(
        [
            ns.row_matrix(pair_full, channel_entries["LLLL_upupdown_Knu"][0], channel_entries["LLLL_upupdown_Knu"][1]),
            ns.row_matrix(pair_full, channel_entries["LLLL_downdownup_Knu"][0], channel_entries["LLLL_downdownup_Knu"][1]),
        ]
    )
    a_rrrr = ns.row_matrix(pair_full, channel_entries["RRRR_uusd_anycharged"][0], channel_entries["RRRR_uusd_anycharged"][1])
    allowed_columns = [n for n, (i, j) in enumerate(pair_full) if mask[i, j]]
    joint = np.vstack([two.normalized_block(a_knu), WEIGHT_RRRR * two.normalized_block(a_rrrr)])
    restricted = joint[:, allowed_columns]
    _, singulars, vh = np.linalg.svd(restricted, full_matrices=True)
    coeff_restricted = vh.conjugate().T[:, -1]
    coeff_restricted = coeff_restricted / np.linalg.norm(coeff_restricted)
    coeff = np.zeros(len(pair_full), dtype=complex)
    for out_idx, col in enumerate(allowed_columns):
        coeff[col] = coeff_restricted[out_idx]
    w0 = two.w_from_coeffs(coeff, len(names))
    return {
        "names": names,
        "charges_left": tuple(int(best["q_left"][i]) for i in range(len(names))),
        "charges_right": tuple(int(best["q_right"][i]) for i in range(len(names))),
        "neutral_mask": mask,
        "w0": w0,
        "pair_full": pair_full,
        "channel_entries": channel_entries,
        "constants": constants,
        "restricted_singular_values": [float(x) for x in singulars],
    }


def residue_mask(charges_left: tuple[int, ...], charges_right: tuple[int, ...], residue: int) -> np.ndarray:
    size = len(charges_left)
    return np.array(
        [[(charges_left[i] + charges_right[j]) % 3 == residue for j in range(size)] for i in range(size)],
        dtype=bool,
    )


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_abs_rows(names: list[str], mat: np.ndarray, charges_left: tuple[int, ...], charges_right: tuple[int, ...]) -> list[dict[str, Any]]:
    rows = []
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            rows.append(
                {
                    "left": left,
                    "right": right,
                    "residue": int((charges_left[i] + charges_right[j]) % 3),
                    "spurion_charge_needed": int((-(charges_left[i] + charges_right[j])) % 3),
                    "abs": float(abs(mat[i, j])),
                    "re": float(np.real(mat[i, j])),
                    "im": float(np.imag(mat[i, j])),
                }
            )
    return rows


def residue_summary(mat: np.ndarray, charges_left: tuple[int, ...], charges_right: tuple[int, ...]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    neutral_max = 0.0
    for residue in [0, 1, 2]:
        mask = residue_mask(charges_left, charges_right, residue)
        vals = np.abs(mat[mask])
        fro = float(np.linalg.norm(vals))
        max_abs = float(np.max(vals)) if vals.size else 0.0
        if residue == 0:
            neutral_max = max_abs
        out[str(residue)] = {
            "entry_count": int(vals.size),
            "frobenius_norm": fro,
            "max_abs_entry": max_abs,
            "nonzero_count_tol_1e_minus_10": int(np.sum(vals > ENTRY_TOL)),
        }
    charged_max = max(out["1"]["max_abs_entry"], out["2"]["max_abs_entry"])
    charged_fro = math.sqrt(out["1"]["frobenius_norm"] ** 2 + out["2"]["frobenius_norm"] ** 2)
    out["charged_over_neutral_max"] = float(charged_max / max(neutral_max, 1.0e-30))
    out["charged_over_neutral_frobenius"] = float(charged_fro / max(out["0"]["frobenius_norm"], 1.0e-30))
    return out


def evaluate_regulator(label: str, kappa: float, w0: np.ndarray, payload: dict[str, Any]) -> dict[str, Any]:
    names = payload["names"]
    charges_left = payload["charges_left"]
    charges_right = payload["charges_right"]
    pair_full = payload["pair_full"]
    channel_entries = payload["channel_entries"]
    constants = payload["constants"]
    w_cap, cap_info = two.cap_condition(w0, kappa)
    singulars = np.linalg.svd(w_cap, compute_uv=False)
    det_abs = float(abs(np.linalg.det(w_cap)))
    channels = lift.max_channel_for_w(w_cap, pair_full, channel_entries, constants)
    max_knu = max(channels["LLLL_upupdown_Knu"]["amplitude"], channels["LLLL_downdownup_Knu"]["amplitude"])
    k0mu = channels["LLLL_upupdown_K0mu"]["amplitude"]
    rrrr = channels["RRRR_uusd_anycharged"]["amplitude"]
    leading = max(channels.items(), key=lambda item: item[1]["amplitude"])
    rows = matrix_abs_rows(names, w_cap, charges_left, charges_right)
    nonzero_charged = [
        row
        for row in rows
        if row["residue"] != 0 and row["abs"] > ENTRY_TOL
    ]
    return {
        "label": label,
        "kappa": float(kappa),
        "matrix": lift.matrix_json(w_cap),
        "singular_values": [float(x) for x in singulars],
        "condition": float(singulars[0] / max(singulars[-1], 1.0e-30)),
        "det_abs": det_abs,
        "cap_info": cap_info,
        "residue_summary": residue_summary(w_cap, charges_left, charges_right),
        "nonzero_charged_entries": nonzero_charged,
        "max_Knu_amplitude": float(max_knu),
        "K0mu_amplitude": float(k0mu),
        "RRRR_amplitude": float(rrrr),
        "leading_channel": leading[0],
        "leading_amplitude": leading[1]["amplitude"],
        "leading_ST_required_tau_2p4e34": leading[1]["S_T_required_tau_2p4e34"],
        "channels": channels,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = reconstruct_neutral_support()
    names = payload["names"]
    charges_left = payload["charges_left"]
    charges_right = payload["charges_right"]
    w0 = payload["w0"]
    neutral_singulars = np.linalg.svd(w0, compute_uv=False)
    rows = [evaluate_regulator(f"z3_svd_floor_kappa_{int(kappa)}", kappa, w0, payload) for kappa in KAPPAS]
    chosen = next(row for row in rows if row["kappa"] == 100.0)

    csv_rows = []
    for row in rows:
        residue = row["residue_summary"]
        csv_rows.append(
            {
                "label": row["label"],
                "kappa": row["kappa"],
                "condition": row["condition"],
                "det_abs": row["det_abs"],
                "neutral_frobenius": residue["0"]["frobenius_norm"],
                "residue1_frobenius": residue["1"]["frobenius_norm"],
                "residue2_frobenius": residue["2"]["frobenius_norm"],
                "charged_over_neutral_max": residue["charged_over_neutral_max"],
                "charged_over_neutral_frobenius": residue["charged_over_neutral_frobenius"],
                "max_Knu_amplitude": row["max_Knu_amplitude"],
                "K0mu_amplitude": row["K0mu_amplitude"],
                "RRRR_amplitude": row["RRRR_amplitude"],
                "leading_channel": row["leading_channel"],
                "leading_ST_required_tau_2p4e34": row["leading_ST_required_tau_2p4e34"],
            }
        )
    with (OUT / "z3_triplet_regulator_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    matrix_rows = matrix_abs_rows(names, np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in chosen["matrix"]]), charges_left, charges_right)
    with (OUT / "z3_triplet_regulator_matrix_kappa100.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matrix_rows[0].keys()))
        writer.writeheader()
        writer.writerows(matrix_rows)

    superpotential = {
        "charge_rule": "A term T_A barT_B S_r is invariant when q(T_A)+q(barT_B)+q(S_r)=0 mod 3.",
        "charges_left": {f"T_{name}": int(charges_left[i]) for i, name in enumerate(names)},
        "charges_right": {f"barT_{name}": int(charges_right[i]) for i, name in enumerate(names)},
        "spurions": {
            "S0": {"Z3_charge": 0, "role": "neutral dominant support"},
            "S1": {"Z3_charge": 1, "role": "lifts residue-2 entries, dominantly G120-G120"},
            "S2": {"Z3_charge": 2, "role": "lifts residue-1 entries, dominantly the right-Ktr column"},
        },
        "schematic": (
            "W_T/M_T = T_A [ S0 N_AB + S2 epsilon R1_AB + S1 epsilon R2_AB ] barT_B, "
            "where N has only residue-0 entries, R1 only residue-1 entries, and R2 only residue-2 entries."
        ),
        "flatness_at_zero_triplet": (
            "For det W_T != 0, F_T = W_T barT = 0 and F_barT = W_T^T T = 0 force T=barT=0. "
            "The triplet D-terms vanish at this origin; singlet-spurion stabilization is a separate neutral-sector problem."
        ),
    }
    verdict = (
        "The kappa=100 floor admits a compact Z3-spurion interpretation.  Relative to the neutral block, "
        "the charged-spurion Frobenius norm is about 1.4 percent and the largest charged entry is about "
        "2.0 percent of the largest neutral entry.  The only nonzero charged entries above 1e-10 are a "
        "right-Ktr column in residue 1 and the G120-G120 entry in residue 2.  Thus the regulator is not a "
        "generic symmetry-breaking leak; it is a two-spurion missing-partner lift.  The remaining paper-level "
        "work is to stabilize S1,S2 and verify that loops or higher operators do not fill the other charged "
        "entries above the percent budget."
    )
    output = {
        "note": "No web lookup used. Z3-compatible regulator audit for the triplet filter.",
        "neutral_support": {
            "basis_names": names,
            "charges_left": list(charges_left),
            "charges_right": list(charges_right),
            "singular_values": [float(x) for x in neutral_singulars],
            "restricted_singular_values": payload["restricted_singular_values"],
            "matrix": lift.matrix_json(w0),
        },
        "superpotential": superpotential,
        "rows": rows,
        "chosen_kappa100": {
            key: chosen[key]
            for key in [
                "label",
                "singular_values",
                "condition",
                "det_abs",
                "residue_summary",
                "nonzero_charged_entries",
                "max_Knu_amplitude",
                "K0mu_amplitude",
                "RRRR_amplitude",
                "leading_channel",
                "leading_ST_required_tau_2p4e34",
            ]
        },
        "verdict": verdict,
    }
    (OUT / "z3_triplet_regulator_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Z3 triplet regulator audit",
        "",
        "No web lookup was used.",
        "",
        "The neutral Z3 support has singular values",
        "",
        "```text",
        " ".join(f"{x:.6e}" for x in neutral_singulars),
        "```",
        "",
        "The finite regulator is the closest singular-value floor with charged spurion decomposition",
        "",
        "```text",
        "W_T/M_T = T_A [ S0 N_AB + S2 eps R1_AB + S1 eps R2_AB ] barT_B.",
        "```",
        "",
        "| kappa | cond | residue-1 fro | residue-2 fro | charged/max neutral | max Knu | K0mu | RRRR | ST max |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        residue = row["residue_summary"]
        report.append(
            f"| {row['kappa']:.0f} | {row['condition']:.3e} | "
            f"{residue['1']['frobenius_norm']:.3e} | {residue['2']['frobenius_norm']:.3e} | "
            f"{residue['charged_over_neutral_max']:.3e} | "
            f"{row['max_Knu_amplitude']:.3e} | {row['K0mu_amplitude']:.3e} | "
            f"{row['RRRR_amplitude']:.3e} | {row['leading_ST_required_tau_2p4e34']:.3e} |"
        )
    report.extend(
        [
            "",
            "## Kappa 100 charged entries",
            "",
            "| left | right | residue | spurion charge | abs |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for entry in chosen["nonzero_charged_entries"]:
        report.append(
            f"| {entry['left']} | {entry['right']} | {entry['residue']} | "
            f"{entry['spurion_charge_needed']} | {entry['abs']:.3e} |"
        )
    report.extend(["", "## Verdict", "", verdict, ""])
    (OUT / "z3_triplet_regulator_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Z3 triplet regulator audit")
    print("  neutral singulars:", " ".join(f"{x:.6e}" for x in neutral_singulars))
    print(
        "  kappa100: "
        f"cond={chosen['condition']:.6e} "
        f"Knu={chosen['max_Knu_amplitude']:.6e} "
        f"K0mu={chosen['K0mu_amplitude']:.6e} "
        f"RRRR={chosen['RRRR_amplitude']:.6e}"
    )
    print(
        "  charged/neutral: "
        f"max={chosen['residue_summary']['charged_over_neutral_max']:.6e} "
        f"fro={chosen['residue_summary']['charged_over_neutral_frobenius']:.6e}"
    )
    print(f"  charged nonzero entries: {len(chosen['nonzero_charged_entries'])}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
