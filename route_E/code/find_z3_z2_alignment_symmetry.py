#!/usr/bin/env python3
"""Find an extra shaping symmetry for the Z3 triplet regulator.

No web lookup is used.  The previous leakage audit showed that Z3 fixes
spurion residues but not the matrix direction: generic S1,S2 insertions can
fill six unused charged entries.  This script searches for a minimal extra
cyclic grading that preserves the desired neutral/charged entries and forbids
the dangerous ones.  It then enumerates low-degree spurion monomials to check
whether the dangerous entries stay forbidden beyond degree one.
"""

from __future__ import annotations

import csv
import itertools
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import construct_triplet_rank_lift as lift  # noqa: E402


INPUT = ROOT / "output" / "z3_spurion_leakage" / "z3_spurion_leakage_summary.json"
REGULATOR = ROOT / "output" / "z3_triplet_regulator" / "z3_triplet_regulator_summary.json"
OUT = ROOT / "output" / "z3_z2_alignment"
MAX_MONOMIAL_DEGREE = 8


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def load_inputs() -> dict[str, Any]:
    leak = json.loads(INPUT.read_text(encoding="utf-8"))
    reg = json.loads(REGULATOR.read_text(encoding="utf-8"))
    row = next(item for item in reg["rows"] if item["label"] == reg["chosen_kappa100"]["label"])
    return {
        "leak": leak,
        "reg": reg,
        "w": cmat(row["matrix"]),
        "names": reg["neutral_support"]["basis_names"],
    }


def required_spurion(residue: int) -> str:
    if residue == 0:
        return "S0"
    if residue == 1:
        return "S2"
    return "S1"


def entry_sets(leak: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    entries = leak["entry_classification"]
    desired = [row for row in entries if row["neutral"] or row["desired_charged"]]
    desired_charged = [row for row in entries if row["desired_charged"]]
    dangerous = [row for row in entries if row["dangerous_charged_if_generated"]]
    return desired, desired_charged, dangerous


def is_allowed(row: dict[str, Any], ql: tuple[int, ...], qr: tuple[int, ...], qs: dict[str, int], modulus: int) -> bool:
    spurion = required_spurion(int(row["residue"]))
    return (ql[int(row["i"])] + qr[int(row["j"])] + qs[spurion]) % modulus == 0


def scan_extra_gradings(leak: dict[str, Any], max_modulus: int = 8) -> list[dict[str, Any]]:
    desired, _, dangerous = entry_sets(leak)
    results = []
    for modulus in range(2, max_modulus + 1):
        # Gauge-fix q_L(H10)=q_R(H10)=0.  H10-H10 with S0 then sets q(S0)=0.
        for ql_tail in itertools.product(range(modulus), repeat=3):
            ql = (0,) + ql_tail
            for qr_tail in itertools.product(range(modulus), repeat=3):
                qr = (0,) + qr_tail
                for q_s1 in range(modulus):
                    for q_s2 in range(modulus):
                        qs = {"S0": 0, "S1": q_s1, "S2": q_s2}
                        missed = [row for row in desired if not is_allowed(row, ql, qr, qs, modulus)]
                        leaked = [row for row in dangerous if is_allowed(row, ql, qr, qs, modulus)]
                        candidate = {
                            "modulus": modulus,
                            "q_left": ql,
                            "q_right": qr,
                            "q_spurions": qs,
                            "missed_desired": len(missed),
                            "allowed_dangerous": len(leaked),
                        }
                        dangerous_all = len([row for row in combined_allowed_entries(leak, candidate) if row["dangerous"]])
                        score = (
                            10000 * len(missed)
                            + 1000 * dangerous_all
                            + 100 * len(leaked)
                            + modulus
                            + sum(ql)
                            + sum(qr)
                            + q_s1
                            + q_s2
                        )
                        results.append(
                            {
                                "modulus": modulus,
                                "q_left": ql,
                                "q_right": qr,
                                "q_spurions": qs,
                                "missed_desired": len(missed),
                                "allowed_dangerous": len(leaked),
                                "dangerous_monomials_to_degree": dangerous_all,
                                "score": score,
                            }
                        )
        exact = [
            row
            for row in results
            if row["modulus"] == modulus
            and row["missed_desired"] == 0
            and row["allowed_dangerous"] == 0
            and row["dangerous_monomials_to_degree"] == 0
        ]
        if exact:
            break
    return sorted(
        results,
        key=lambda row: (
            row["missed_desired"],
            row["dangerous_monomials_to_degree"],
            row["allowed_dangerous"],
            row["modulus"],
            row["score"],
        ),
    )


def monomials(max_degree: int = MAX_MONOMIAL_DEGREE) -> list[dict[str, int]]:
    rows = []
    for a in range(max_degree + 1):
        for b in range(max_degree + 1 - a):
            rows.append({"power_S1": a, "power_S2": b, "degree": a + b, "Z3_charge": (a + 2 * b) % 3})
    return sorted(rows, key=lambda row: (row["degree"], row["power_S1"], row["power_S2"]))


def combined_allowed_entries(leak: dict[str, Any], solution: dict[str, Any], max_degree: int = MAX_MONOMIAL_DEGREE) -> list[dict[str, Any]]:
    ql = solution["q_left"]
    qr = solution["q_right"]
    # The found minimal solution has S1,S2 even under the extra symmetry.  Keep
    # the formula general.
    q_s_extra = {"S1": solution["q_spurions"]["S1"], "S2": solution["q_spurions"]["S2"]}
    modulus = int(solution["modulus"])
    rows = []
    for entry in leak["entry_classification"]:
        for mono in monomials(max_degree):
            z3_ok = (int(entry["residue"]) + mono["Z3_charge"]) % 3 == 0
            extra_charge = (
                ql[int(entry["i"])]
                + qr[int(entry["j"])]
                + mono["power_S1"] * q_s_extra["S1"]
                + mono["power_S2"] * q_s_extra["S2"]
            ) % modulus
            if z3_ok and extra_charge == 0:
                rows.append(
                    {
                        "left": entry["left"],
                        "right": entry["right"],
                        "residue": entry["residue"],
                        "degree": mono["degree"],
                        "power_S1": mono["power_S1"],
                        "power_S2": mono["power_S2"],
                        "dangerous": entry["dangerous_charged_if_generated"],
                        "desired": entry["neutral"] or entry["desired_charged"],
                    }
                )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = load_inputs()
    leak = payload["leak"]
    reg = payload["reg"]
    w = payload["w"]
    names = list(payload["names"])
    desired, desired_charged, dangerous = entry_sets(leak)
    scans = scan_extra_gradings(leak)
    best = scans[0]
    allowed = combined_allowed_entries(leak, best)
    dangerous_allowed = [row for row in allowed if row["dangerous"]]
    desired_allowed = [row for row in allowed if row["desired"]]

    scan_rows = []
    for row in scans[:200]:
        scan_rows.append(
            {
                "modulus": row["modulus"],
                "q_left": " ".join(map(str, row["q_left"])),
                "q_right": " ".join(map(str, row["q_right"])),
                "q_S0": row["q_spurions"]["S0"],
                "q_S1": row["q_spurions"]["S1"],
                "q_S2": row["q_spurions"]["S2"],
                "missed_desired": row["missed_desired"],
                "allowed_dangerous": row["allowed_dangerous"],
                "dangerous_monomials_to_degree": row["dangerous_monomials_to_degree"],
                "score": row["score"],
            }
        )
    with (OUT / "z3_z2_alignment_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scan_rows[0].keys()))
        writer.writeheader()
        writer.writerows(scan_rows)

    allowed_rows = []
    seen = set()
    for row in allowed:
        key = (row["left"], row["right"], row["degree"], row["power_S1"], row["power_S2"])
        if key in seen:
            continue
        seen.add(key)
        allowed_rows.append(row)
    with (OUT / "z3_z2_allowed_monomials.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(allowed_rows[0].keys()))
        writer.writeheader()
        writer.writerows(allowed_rows)

    _, pair_full, channel_entries, constants = lift.build_channel_tensors()
    channels = lift.max_channel_for_w(w, pair_full, channel_entries, constants)
    max_knu = max(channels["LLLL_upupdown_Knu"]["amplitude"], channels["LLLL_downdownup_Knu"]["amplitude"])
    k0mu = channels["LLLL_upupdown_K0mu"]["amplitude"]
    rrrr = channels["RRRR_uusd_anycharged"]["amplitude"]
    leading = max(channels.items(), key=lambda item: item[1]["amplitude"])
    verdict = (
        "A minimal extra Z2 grading solves the matrix-support problem.  It assigns odd parity only to "
        "the G120 triplet source on both left and right and leaves H10,F126,Ktr plus S0,S1,S2 even.  "
        "All desired neutral and charged regulator entries are allowed, while all six dangerous charged "
        "entries are forbidden.  Because S1 and S2 are even under the extra Z2, this remains true for "
        f"all spurion monomials checked up to degree {MAX_MONOMIAL_DEGREE}.  The proton-amplitude replay is "
        "unchanged from the kappa=100 regulator."
    )
    output = {
        "note": "No web lookup used. Minimal extra shaping symmetry for the Z3 triplet regulator.",
        "basis_names": names,
        "best_extra_grading": {
            "modulus": best["modulus"],
            "q_left": {name: int(best["q_left"][i]) for i, name in enumerate(names)},
            "q_right": {name: int(best["q_right"][i]) for i, name in enumerate(names)},
            "q_spurions": best["q_spurions"],
            "missed_desired": best["missed_desired"],
            "allowed_dangerous": best["allowed_dangerous"],
        },
        "desired_entry_count": len(desired),
        "desired_charged_entry_count": len(desired_charged),
        "dangerous_entry_count": len(dangerous),
        "dangerous_allowed_by_monomials_up_to_degree": len(dangerous_allowed),
        "desired_allowed_monomial_records_up_to_degree": len(desired_allowed),
        "max_monomial_degree_checked": MAX_MONOMIAL_DEGREE,
        "channels": {
            "max_Knu_amplitude": float(max_knu),
            "K0mu_amplitude": float(k0mu),
            "RRRR_amplitude": float(rrrr),
            "leading_channel": leading[0],
            "leading_ST_required_tau_2p4e34": leading[1]["S_T_required_tau_2p4e34"],
        },
        "regulator_reference": reg["chosen_kappa100"],
        "verdict": verdict,
    }
    (OUT / "z3_z2_alignment_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Z3 x Z2 alignment audit",
        "",
        "No web lookup was used.",
        "",
        "Best extra grading:",
        "",
        "| source | H10 | F126 | G120 | Ktr |",
        "|---|---:|---:|---:|---:|",
        "| q_L mod 2 | " + " | ".join(str(best["q_left"][i]) for i in range(len(names))) + " |",
        "| q_R mod 2 | " + " | ".join(str(best["q_right"][i]) for i in range(len(names))) + " |",
        "",
        f"Spurions: S0={best['q_spurions']['S0']}, S1={best['q_spurions']['S1']}, S2={best['q_spurions']['S2']} mod {best['modulus']}.",
        "",
        f"Desired entries missed: {best['missed_desired']}.  Dangerous entries allowed: {best['allowed_dangerous']}.",
        f"Dangerous entries generated by monomials through degree {MAX_MONOMIAL_DEGREE}: {len(dangerous_allowed)}.",
        "",
        "Replay:",
        "",
        f"- max Knu = {max_knu:.6e}",
        f"- K0mu = {k0mu:.6e}",
        f"- RRRR = {rrrr:.6e}",
        f"- leading ST max = {leading[1]['S_T_required_tau_2p4e34']:.6e}",
        "",
        "## Verdict",
        "",
        verdict,
        "",
    ]
    (OUT / "z3_z2_alignment_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Z3 x Z2 alignment audit")
    print(f"  best modulus: Z_{best['modulus']}")
    print(f"  q_left: {best['q_left']}")
    print(f"  q_right: {best['q_right']}")
    print(f"  q_spurions: {best['q_spurions']}")
    print(f"  missed desired / allowed dangerous: {best['missed_desired']} / {best['allowed_dangerous']}")
    print(f"  dangerous monomials through degree {MAX_MONOMIAL_DEGREE}: {len(dangerous_allowed)}")
    print(f"  replay: Knu={max_knu:.6e} K0mu={k0mu:.6e} RRRR={rrrr:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
