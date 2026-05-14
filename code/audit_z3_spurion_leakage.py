#!/usr/bin/env python3
"""Audit Z3 spurion stabilization and leakage for the triplet regulator.

No web lookup is used.  The Z3 regulator needs charged spurions S1,S2 to lift
the singular neutral support.  A generic Z3 spurion sector would also allow
other charged entries.  This script enumerates the low-degree spurion monomials
and computes rigorous entrywise leakage bounds for the monitored d=5 proton
operators.

If dangerous charged entries are bounded by |delta W_AB| <= eta, then for each
linear proton row C = a . W we use

    |C| <= |C_0| + eta sum_{dangerous j} |a_j|.

This is a deterministic worst-phase bound, not a random scan.
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
import scan_triplet_mixing_nullspace as ns  # noqa: E402


INPUT = ROOT / "output" / "z3_triplet_regulator" / "z3_triplet_regulator_summary.json"
OUT = ROOT / "output" / "z3_spurion_leakage"

MONOMIAL_MAX_DEGREE = 6
DANGEROUS_ENTRY_TOL = 1.0e-10
ETA_GRID = [0.0, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2]
AMPLITUDE_CLASS_TARGET = 2.0e-4


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def coeffs_from_w(w: np.ndarray) -> np.ndarray:
    return np.array([w[i, j] for i in range(w.shape[0]) for j in range(w.shape[1])], dtype=complex)


def load_payload() -> dict[str, Any]:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    chosen = data["chosen_kappa100"]
    matrix_row = next(row for row in data["rows"] if row["label"] == chosen["label"])
    w = cmat(matrix_row["matrix"])
    charges_left = tuple(int(x) for x in data["neutral_support"]["charges_left"])
    charges_right = tuple(int(x) for x in data["neutral_support"]["charges_right"])
    names = data["neutral_support"]["basis_names"]
    return {"data": data, "chosen": chosen, "w": w, "charges_left": charges_left, "charges_right": charges_right, "names": names}


def residue(charges_left: tuple[int, ...], charges_right: tuple[int, ...], i: int, j: int) -> int:
    return int((charges_left[i] + charges_right[j]) % 3)


def enumerate_monomials(max_degree: int = MONOMIAL_MAX_DEGREE) -> list[dict[str, Any]]:
    rows = []
    for a in range(max_degree + 1):
        for b in range(max_degree + 1 - a):
            degree = a + b
            charge = (a + 2 * b) % 3
            supported_residue = (-charge) % 3
            rows.append(
                {
                    "monomial": f"S1^{a} S2^{b}",
                    "power_S1": a,
                    "power_S2": b,
                    "degree": degree,
                    "Z3_charge": charge,
                    "supports_entry_residue": supported_residue,
                }
            )
    return sorted(rows, key=lambda row: (row["degree"], row["power_S1"], row["power_S2"]))


def entry_classification(names: list[str], w: np.ndarray, charges_left: tuple[int, ...], charges_right: tuple[int, ...]) -> list[dict[str, Any]]:
    desired = []
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            r = residue(charges_left, charges_right, i, j)
            abs_w = float(abs(w[i, j]))
            desired.append(
                {
                    "left": left,
                    "right": right,
                    "i": i,
                    "j": j,
                    "residue": r,
                    "spurion_charge_needed": int((-r) % 3),
                    "abs_regulator_entry": abs_w,
                    "neutral": bool(r == 0),
                    "desired_charged": bool(r != 0 and abs_w > DANGEROUS_ENTRY_TOL),
                    "dangerous_charged_if_generated": bool(r != 0 and abs_w <= DANGEROUS_ENTRY_TOL),
                }
            )
    return desired


def row_bound_for_channel(
    coeff_base: np.ndarray,
    pair_full: list[tuple[int, int]],
    tensors: dict[tuple[int, int], np.ndarray],
    entries: list[tuple[int, int, int, int]],
    dangerous_cols: list[int],
    eta: float,
) -> dict[str, Any]:
    best = {"bound": -1.0, "base_abs": 0.0, "l1_dangerous": 0.0, "entry": entries[0], "base_value": 0.0j}
    for entry in entries:
        row = np.array([tensors[pair][entry] for pair in pair_full], dtype=complex)
        base = complex(np.dot(row, coeff_base))
        l1 = float(np.sum(np.abs(row[dangerous_cols]))) if dangerous_cols else 0.0
        bound = float(abs(base) + eta * l1)
        if bound > best["bound"]:
            best = {
                "bound": bound,
                "base_abs": float(abs(base)),
                "l1_dangerous": l1,
                "entry": entry,
                "base_value": base,
            }
    return best


def leakage_bounds(w: np.ndarray, dangerous_cols: list[int], eta: float) -> dict[str, Any]:
    _, pair_full, channel_entries, constants = lift.build_channel_tensors()
    coeff_base = coeffs_from_w(w)
    per_channel = {}
    for channel, (tensors, entries) in channel_entries.items():
        item = row_bound_for_channel(coeff_base, pair_full, tensors, entries, dangerous_cols, eta)
        life = ns.lifetime_from_amplitude(item["bound"], constants)
        per_channel[channel] = {
            "bound": item["bound"],
            "base_abs_at_selected_entry": item["base_abs"],
            "l1_dangerous_at_selected_entry": item["l1_dangerous"],
            "selected_index": "".join(str(x) for x in item["entry"]),
            "S_T_required_tau_2p4e34": life["S_T_required_tau_2p4e34"],
        }
    max_knu = max(per_channel["LLLL_upupdown_Knu"]["bound"], per_channel["LLLL_downdownup_Knu"]["bound"])
    k0mu = per_channel["LLLL_upupdown_K0mu"]["bound"]
    rrrr = per_channel["RRRR_uusd_anycharged"]["bound"]
    leading = max(per_channel.items(), key=lambda item: item[1]["bound"])
    return {
        "eta": float(eta),
        "max_Knu_bound": float(max_knu),
        "K0mu_bound": float(k0mu),
        "RRRR_bound": float(rrrr),
        "worst_monitored_bound": float(max(max_knu, k0mu, rrrr)),
        "leading_channel": leading[0],
        "leading_bound": leading[1]["bound"],
        "leading_ST_required_tau_2p4e34": leading[1]["S_T_required_tau_2p4e34"],
        "per_channel": per_channel,
    }


def eta_for_target(w: np.ndarray, dangerous_cols: list[int], target: float) -> float:
    _, pair_full, channel_entries, _ = lift.build_channel_tensors()
    coeff_base = coeffs_from_w(w)
    eta_limits = []
    for tensors, entries in channel_entries.values():
        for entry in entries:
            row = np.array([tensors[pair][entry] for pair in pair_full], dtype=complex)
            base = float(abs(np.dot(row, coeff_base)))
            l1 = float(np.sum(np.abs(row[dangerous_cols]))) if dangerous_cols else 0.0
            if l1 > 0.0:
                eta_limits.append(max(0.0, (target - base) / l1))
            elif base > target:
                return 0.0
    return float(min(eta_limits)) if eta_limits else math.inf


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = load_payload()
    names = payload["names"]
    w = payload["w"]
    charges_left = payload["charges_left"]
    charges_right = payload["charges_right"]
    entries = entry_classification(names, w, charges_left, charges_right)
    dangerous = [row for row in entries if row["dangerous_charged_if_generated"]]
    desired_charged = [row for row in entries if row["desired_charged"]]
    pair_full = [(i, j) for i in range(len(names)) for j in range(len(names))]
    dangerous_cols = [n for n, (i, j) in enumerate(pair_full) if any(row["i"] == i and row["j"] == j for row in dangerous)]

    monomials = enumerate_monomials()
    bounds = [leakage_bounds(w, dangerous_cols, eta) for eta in ETA_GRID]
    eta_target = eta_for_target(w, dangerous_cols, AMPLITUDE_CLASS_TARGET)

    with (OUT / "z3_spurion_monomials.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(monomials[0].keys()))
        writer.writeheader()
        writer.writerows(monomials)

    with (OUT / "z3_triplet_entry_classification.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "left",
            "right",
            "residue",
            "spurion_charge_needed",
            "abs_regulator_entry",
            "neutral",
            "desired_charged",
            "dangerous_charged_if_generated",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in entries:
            writer.writerow({key: row[key] for key in fieldnames})

    csv_rows = []
    for row in bounds:
        csv_rows.append(
            {
                "eta": row["eta"],
                "max_Knu_bound": row["max_Knu_bound"],
                "K0mu_bound": row["K0mu_bound"],
                "RRRR_bound": row["RRRR_bound"],
                "worst_monitored_bound": row["worst_monitored_bound"],
                "leading_channel": row["leading_channel"],
                "leading_ST_required_tau_2p4e34": row["leading_ST_required_tau_2p4e34"],
            }
        )
    with (OUT / "z3_spurion_leakage_bounds.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    stabilization = {
        "driving_superpotential": (
            "W_drive = X1(S1^3-v1^3) + X2(S2^3-v2^3) + X12(S1 S2-v12^2), "
            "with R(Xi)=2 and <Xi>=0."
        ),
        "F_flat_solution": (
            "F_X1=F_X2=F_X12=0 fixes S1^3=v1^3, S2^3=v2^3, and S1 S2=v12^2. "
            "At X1=X2=X12=0, F_S1=F_S2=0.  Consistency requires v12^6=v1^3 v2^3 up to the chosen Z3 branch."
        ),
        "remaining_alignment_need": (
            "Z3 and singlet F-flatness alone do not select which residue-1 or residue-2 triplet entries are generated. "
            "A matrix-alignment driver or additional shaping symmetry is still needed to keep the six dangerous charged entries small."
        ),
    }
    verdict = (
        "The spurion magnitudes can be fixed by a standard F-flat singlet driver, but Z3 alone is not enough: "
        "at degree one, S2 can populate all residue-1 entries and S1 can populate all residue-2 entries.  "
        "The regulator wants only four charged entries; six additional charged entries are dangerous if generated.  "
        f"The deterministic leakage bound keeps all monitored amplitudes below {AMPLITUDE_CLASS_TARGET:.1e} only for "
        f"entrywise dangerous leakage eta <= {eta_target:.3e}.  This is only a mild suppression relative to the "
        "desired 1e-2 lift if one asks solely for the proxy proton-amplitude class, but a much stronger alignment "
        "is needed if the paper wants structural support purity at the 1e-4 level.  Therefore the next "
        "model-building requirement is a matrix-alignment driver for the unused charged entries."
    )
    output = {
        "note": "No web lookup used. Z3 spurion monomial and leakage audit.",
        "charges": {
            "left": {name: int(charges_left[i]) for i, name in enumerate(names)},
            "right": {name: int(charges_right[i]) for i, name in enumerate(names)},
            "S1": 1,
            "S2": 2,
        },
        "monomial_max_degree": MONOMIAL_MAX_DEGREE,
        "monomials": monomials,
        "entry_classification": entries,
        "desired_charged_entries": desired_charged,
        "dangerous_charged_entries": dangerous,
        "leakage_bounds": bounds,
        "amplitude_class_target": AMPLITUDE_CLASS_TARGET,
        "eta_max_for_target": eta_target,
        "stabilization": stabilization,
        "verdict": verdict,
    }
    (OUT / "z3_spurion_leakage_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Z3 spurion leakage audit",
        "",
        "No web lookup was used.",
        "",
        "Driving-sector candidate:",
        "",
        "```text",
        stabilization["driving_superpotential"],
        "```",
        "",
        "F-flatness: " + stabilization["F_flat_solution"],
        "",
        "## Dangerous charged entries",
        "",
        "| left | right | residue | spurion charge |",
        "|---|---|---:|---:|",
    ]
    for row in dangerous:
        report.append(f"| {row['left']} | {row['right']} | {row['residue']} | {row['spurion_charge_needed']} |")
    report.extend(
        [
            "",
            "## Leakage bounds",
            "",
            "| eta | max Knu | K0mu | RRRR | worst | leading | ST max |",
            "|---:|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in bounds:
        report.append(
            f"| {row['eta']:.1e} | {row['max_Knu_bound']:.3e} | {row['K0mu_bound']:.3e} | "
            f"{row['RRRR_bound']:.3e} | {row['worst_monitored_bound']:.3e} | "
            f"`{row['leading_channel']}` | {row['leading_ST_required_tau_2p4e34']:.3e} |"
        )
    report.extend(
        [
            "",
            f"Target eta for worst monitored amplitude <= {AMPLITUDE_CLASS_TARGET:.1e}: `{eta_target:.3e}`.",
            "",
            "## Verdict",
            "",
            verdict,
            "",
        ]
    )
    (OUT / "z3_spurion_leakage_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Z3 spurion leakage audit")
    print(f"  desired charged entries: {len(desired_charged)}")
    print(f"  dangerous charged entries: {len(dangerous)}")
    print(f"  eta max for worst amplitude <= {AMPLITUDE_CLASS_TARGET:.1e}: {eta_target:.6e}")
    for row in bounds:
        if row["eta"] in [0.0, 1.0e-4, 1.0e-3, 1.0e-2]:
            print(
                f"  eta={row['eta']:.1e}: "
                f"Knu={row['max_Knu_bound']:.6e} "
                f"K0mu={row['K0mu_bound']:.6e} "
                f"RRRR={row['RRRR_bound']:.6e}"
            )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
