#!/usr/bin/env python3
"""Orbit-stabilizer audit for the F_54 order parameter.

No web lookup is used.  The audit checks whether the traceless-symmetric
SO(10) tensor

    S0 = diag(-2,-2,-2,-2,-2,-2, 3,3,3,3)

can consistently be treated as a fixed constrained order parameter for the
F_54 Clebsch operator used in the threshold model.

The checks are deliberately finite-dimensional:

1. Stabilizer in so(10): [A,S0]=0 gives so(6)+so(4).
2. Orbit dimension: dim SO(10)-dim(SO(6)xSO(4)) = 24.
3. Two-form Clebsch spectrum: A_ij gets eigenvalue s_i+s_j.
4. Projector check: P_X=-(9/25)(F-2)(F+4/3) isolates (6,2,2).
5. Optional dynamical realization: a single 54_H cubic superpotential has
   an F-flat and D-flat vacuum on this orbit.
6. Landau audit: compare a fixed background, one propagating 54, and the
   previously documented large 54/210 alignment sector.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "f54_orbit_stabilizer"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"

S0 = [-2.0] * 6 + [3.0] * 4
N = len(S0)

T_INDEX = {
    "45": 8.0,
    "54": 12.0,
    "210": 56.0,
}
GAUGE_CONTRIBUTION = 24.0


def trace_power(power: int) -> float:
    return float(sum(value**power for value in S0))


def minimal_polynomial_residuals() -> list[float]:
    # (S+2 I)(S-3 I)=0, equivalently S^2-S-6 I=0.
    return [value * value - value - 6.0 for value in S0]


def fterm_residuals(m: float = 1.0, lam: float = 1.0) -> dict[str, Any]:
    """Check the projected 54_H cubic F-term.

    For W = m/2 Tr S^2 + lam/3 Tr S^3, the traceless-symmetric F-term is

        F = m S + lam (S^2 - Tr(S^2)/10 I).

    Since S0^2 - 6 I = S0, the vacuum S=v S0 is F-flat for v=-m/lam.
    """
    v = -m / lam
    mean_s2 = sum((v * value) ** 2 for value in S0) / N
    residual = [
        m * v * value + lam * ((v * value) ** 2 - mean_s2)
        for value in S0
    ]
    return {
        "m": m,
        "lambda": lam,
        "v": v,
        "mean_TrS2_over_10": mean_s2,
        "residuals": residual,
        "max_abs_residual": max(abs(value) for value in residual),
        "passes": max(abs(value) for value in residual) < 1.0e-12,
    }


def dflat_residual() -> dict[str, Any]:
    """For the real diagonal representative, [S^\dagger,S]=0 exactly."""
    return {
        "commutator_norm_proxy": 0.0,
        "reason": "S0 is real diagonal, so [S0^dagger,S0]=0 and D^a=Tr([S0^dagger,S0] T^a)=0.",
        "passes": True,
    }


def generator_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for i in range(N):
        for j in range(i + 1, N):
            si = S0[i]
            sj = S0[j]
            if i < 6 and j < 6:
                sector = "(15,1,1)"
                block = "color-color"
            elif i >= 6 and j >= 6:
                sector = "(1,3,1)+(1,1,3)"
                block = "weak-weak"
            else:
                sector = "(6,2,2)"
                block = "mixed"
            commutator_weight = sj - si
            clebsch_raw = si + sj
            rows.append(
                {
                    "i": i + 1,
                    "j": j + 1,
                    "block": block,
                    "sector": sector,
                    "s_i": si,
                    "s_j": sj,
                    "commutator_weight": commutator_weight,
                    "stabilizer_generator": abs(commutator_weight) < 1.0e-12,
                    "orbit_generator": abs(commutator_weight) >= 1.0e-12,
                    "clebsch_raw": clebsch_raw,
                    "F54_normalized": clebsch_raw / 3.0,
                }
            )
    return rows


def summarize_generators(rows: list[dict[str, Any]]) -> dict[str, Any]:
    stabilizer_dim = sum(1 for row in rows if row["stabilizer_generator"])
    orbit_dim = sum(1 for row in rows if row["orbit_generator"])
    block_counts = Counter(row["block"] for row in rows)
    raw_counts = Counter(row["clebsch_raw"] for row in rows)
    normalized_counts = Counter(row["F54_normalized"] for row in rows)
    return {
        "so10_dimension": len(rows),
        "stabilizer_dimension": stabilizer_dim,
        "expected_stabilizer_dimension": 15 + 6,
        "orbit_dimension": orbit_dim,
        "expected_orbit_dimension": 45 - 21,
        "block_counts": dict(block_counts),
        "raw_clebsch_counts": {str(key): value for key, value in sorted(raw_counts.items())},
        "normalized_F54_counts": {str(key): value for key, value in sorted(normalized_counts.items())},
        "weak_two_forms_note": "The six weak-weak two-forms split as Lambda^2_+ plus Lambda^2_- = 3+3, but F54 assigns both the same value 2.",
        "passes": stabilizer_dim == 21 and orbit_dim == 24 and raw_counts == Counter({-4.0: 15, 1.0: 24, 6.0: 6}),
    }


def projector_rows() -> list[dict[str, Any]]:
    fragments = [
        ("Sigma_L", "(1,3,1)", 2.0),
        ("Sigma_R", "(1,1,3)", 2.0),
        ("Sigma8", "(15,1,1)", -4.0 / 3.0),
        ("X_622", "(6,2,2)", 1.0 / 3.0),
    ]
    out = []
    for name, sector, f54 in fragments:
        px = -(9.0 / 25.0) * (f54 - 2.0) * (f54 + 4.0 / 3.0)
        out.append(
            {
                "fragment": name,
                "sector": sector,
                "F54": f54,
                "P_X": px,
                "target": 1.0 if name == "X_622" else 0.0,
                "abs_error": abs(px - (1.0 if name == "X_622" else 0.0)),
            }
        )
    return out


def load_vacuum_rows() -> list[dict[str, Any]]:
    payload = json.loads(VACUUM.read_text(encoding="utf-8"))
    return payload["replay_rows"]


def sum_dynkin(counts: dict[str, int]) -> float:
    return float(sum(T_INDEX[rep] * count for rep, count in counts.items()))


def landau_ratio(alpha_inv: float, b10: float) -> float:
    if b10 <= 0.0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def landau_scenarios() -> list[dict[str, Any]]:
    scenarios = [
        {
            "name": "fixed_F54_background_plus_three_45",
            "counts": {"45": 3},
            "interpretation": "F54 is a constrained background/order parameter; only the 45 mediator triplet propagates above M_G.",
        },
        {
            "name": "one_propagating_54_plus_three_45",
            "counts": {"45": 3, "54": 1},
            "interpretation": "A single dynamical 54_H fixes the orbit through W54, plus the 45 mediator triplet.",
        },
        {
            "name": "large_54_210_alignment_plus_three_45",
            "counts": {"45": 3, "54": 6, "210": 2},
            "interpretation": "Previously documented large source/driver alignment tower treated as propagating.",
        },
    ]
    out = []
    for scenario in scenarios:
        sum_t = sum_dynkin(scenario["counts"])
        b10 = sum_t - GAUGE_CONTRIBUTION
        rows = []
        for vac in load_vacuum_rows():
            r = float(vac["R"])
            alpha_inv = float(vac["alphaG_inv"])
            alpha_at_r = alpha_inv - b10 * math.log(r) / (2.0 * math.pi)
            rows.append(
                {
                    "R": r,
                    "alphaG_inv_MG": alpha_inv,
                    "alpha_inv_at_RMG": alpha_at_r,
                    "landau_ratio": landau_ratio(alpha_inv, b10),
                    "passes": alpha_at_r > 0.0,
                }
            )
        out.append(
            {
                "name": scenario["name"],
                "counts": scenario["counts"],
                "sum_T": sum_t,
                "b10": b10,
                "interpretation": scenario["interpretation"],
                "rows": rows,
                "passes_R50_R200": all(row["passes"] for row in rows if row["R"] in (50.0, 200.0)),
            }
        )
    return out


def build_payload() -> dict[str, Any]:
    rows = generator_rows()
    gen_summary = summarize_generators(rows)
    proj_rows = projector_rows()
    fterm = fterm_residuals()
    dflat = dflat_residual()
    landau = landau_scenarios()
    max_projector_error = max(row["abs_error"] for row in proj_rows)
    return {
        "note": "No web lookup used. F54 orbit-stabilizer and constrained-order-parameter audit.",
        "S0": S0,
        "trace_checks": {
            "TrS": trace_power(1),
            "TrS2": trace_power(2),
            "TrS3": trace_power(3),
            "minimal_polynomial": "S0^2 - S0 - 6 I = 0",
            "minimal_polynomial_max_abs_residual": max(abs(value) for value in minimal_polynomial_residuals()),
            "passes": abs(trace_power(1)) < 1.0e-12
            and max(abs(value) for value in minimal_polynomial_residuals()) < 1.0e-12,
        },
        "orbit_stabilizer": gen_summary,
        "fterm_single_54": fterm,
        "dflatness": dflat,
        "projector_rows": proj_rows,
        "landau_scenarios": landau,
        "verdict": {
            "can_be_fixed_constrained_order_parameter": bool(
                gen_summary["passes"]
                and fterm["passes"]
                and dflat["passes"]
                and max_projector_error < 1.0e-12
            ),
            "max_projector_error": max_projector_error,
            "precise_claim": (
                "F54 is consistent as a fixed SO(10)/(SO(6)xSO(4)) constrained order parameter, "
                "and even as a single 54_H cubic F-flat/D-flat order parameter. This validates the "
                "low-index projector use of F54 but does not derive the order parameter uniquely from PSLT."
            ),
            "remaining_gap": (
                "A first-principles origin for selecting this orbit, or a fully specified low-index "
                "source constraint replacing the spurion language, remains an EFT assumption."
            ),
        },
    }


def write_generator_csv(rows: list[dict[str, Any]]) -> None:
    fields = [
        "i",
        "j",
        "block",
        "sector",
        "s_i",
        "s_j",
        "commutator_weight",
        "stabilizer_generator",
        "orbit_generator",
        "clebsch_raw",
        "F54_normalized",
    ]
    with (OUT / "f54_two_form_spectrum.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_landau_csv(scenarios: list[dict[str, Any]]) -> None:
    fields = [
        "scenario",
        "sum_T",
        "b10",
        "R",
        "alphaG_inv_MG",
        "alpha_inv_at_RMG",
        "landau_ratio",
        "passes",
        "interpretation",
    ]
    with (OUT / "f54_order_parameter_landau.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for scenario in scenarios:
            for row in scenario["rows"]:
                writer.writerow(
                    {
                        "scenario": scenario["name"],
                        "sum_T": scenario["sum_T"],
                        "b10": scenario["b10"],
                        "R": row["R"],
                        "alphaG_inv_MG": row["alphaG_inv_MG"],
                        "alpha_inv_at_RMG": row["alpha_inv_at_RMG"],
                        "landau_ratio": row["landau_ratio"],
                        "passes": row["passes"],
                        "interpretation": scenario["interpretation"],
                    }
                )


def write_report(payload: dict[str, Any]) -> None:
    orbit = payload["orbit_stabilizer"]
    fterm = payload["fterm_single_54"]
    lines: list[str] = []
    lines.append("# F54 orbit-stabilizer audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Orbit data")
    lines.append("")
    lines.append("```text")
    lines.append("S0 = diag(-2^6, 3^4)")
    lines.append(f"Tr S0 = {payload['trace_checks']['TrS']:.0f}")
    lines.append(f"Tr S0^2 = {payload['trace_checks']['TrS2']:.0f}")
    lines.append(f"Tr S0^3 = {payload['trace_checks']['TrS3']:.0f}")
    lines.append("S0^2 - S0 - 6 I = 0")
    lines.append("```")
    lines.append("")
    lines.append(
        f"The stabilizer count is {orbit['stabilizer_dimension']} = 15+6, "
        f"so the orbit dimension is {orbit['orbit_dimension']} = 45-21."
    )
    lines.append("")
    lines.append("## Two-form Clebsch spectrum")
    lines.append("")
    lines.append("| block | multiplicity | raw s_i+s_j | normalized F54 |")
    lines.append("|---|---:|---:|---:|")
    lines.append("| color-color | 15 | -4 | -4/3 |")
    lines.append("| weak-weak | 6 | 6 | 2 |")
    lines.append("| mixed | 24 | 1 | 1/3 |")
    lines.append("")
    lines.append("The weak six two-forms split into self-dual and anti-self-dual 3+3,")
    lines.append("but F54 assigns both the same Clebsch value, so F54 alone cannot split L/R.")
    lines.append("")
    lines.append("## Single-54 F/D flatness")
    lines.append("")
    lines.append("For W54 = m/2 Tr S^2 + lambda/3 Tr S^3,")
    lines.append("the projected F-term is F = m S + lambda(S^2 - Tr(S^2) I/10).")
    lines.append(
        f"With m=lambda=1 and v={fterm['v']:.0f}, max |F| = {fterm['max_abs_residual']:.3e}."
    )
    lines.append("D-flatness holds because the representative is real diagonal.")
    lines.append("")
    lines.append("## Projector")
    lines.append("")
    lines.append("P_X = -(9/25)(F54-2)(F54+4/3) gives:")
    lines.append("")
    lines.append("| fragment | F54 | P_X | target |")
    lines.append("|---|---:|---:|---:|")
    for row in payload["projector_rows"]:
        lines.append(
            f"| {row['fragment']} | {row['F54']:.12g} | {row['P_X']:.12g} | {row['target']:.0f} |"
        )
    lines.append("")
    lines.append("## Landau interpretation")
    lines.append("")
    lines.append("| scenario | sum T | b10 | R=50/200 pass |")
    lines.append("|---|---:|---:|---|")
    for scenario in payload["landau_scenarios"]:
        lines.append(
            f"| {scenario['name']} | {scenario['sum_T']:.0f} | {scenario['b10']:.0f} | "
            f"{'yes' if scenario['passes_R50_R200'] else 'no'} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"]["precise_claim"])
    lines.append(payload["verdict"]["remaining_gap"])
    lines.append("")
    (OUT / "f54_orbit_stabilizer_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = generator_rows()
    payload = build_payload()
    write_generator_csv(rows)
    write_landau_csv(payload["landau_scenarios"])
    (OUT / "f54_orbit_stabilizer_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    verdict = payload["verdict"]
    print("F54 orbit-stabilizer audit complete")
    print(f"  can_be_fixed_constrained_order_parameter={verdict['can_be_fixed_constrained_order_parameter']}")
    print(f"  stabilizer_dim={payload['orbit_stabilizer']['stabilizer_dimension']}")
    print(f"  orbit_dim={payload['orbit_stabilizer']['orbit_dimension']}")
    print(f"  single_54_Fterm_max={payload['fterm_single_54']['max_abs_residual']:.3e}")
    print(f"  max_projector_error={verdict['max_projector_error']:.3e}")
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
