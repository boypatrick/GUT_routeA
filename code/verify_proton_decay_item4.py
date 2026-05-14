#!/usr/bin/env python3
"""Item 4: proton-decay operator ledger and lifetime estimates.

The script separates unavoidable dimension-6 gauge exchange from model-dependent
dimension-5 colored-triplet operators.  It intentionally uses approximate local
benchmarks rather than web lookups.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, asdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
YUKAWA_JSON = ROOT / "output" / "yukawa_o2" / "cp1_o2_yukawa_scan.json"
OUT = ROOT / "output" / "proton_decay"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def f(num: int, den: int = 1) -> Fraction:
    return Fraction(num, den)


def frac_str(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


@dataclass(frozen=True)
class Field:
    name: str
    y: Fraction
    b: Fraction
    l: Fraction
    dimension: Fraction


FIELDS = {
    "Q": Field("Q", f(1, 6), f(1, 3), f(0), f(3, 2)),
    "L": Field("L", f(-1, 2), f(0), f(1), f(3, 2)),
    "u^c": Field("u^c", f(-2, 3), f(-1, 3), f(0), f(3, 2)),
    "d^c": Field("d^c", f(1, 3), f(-1, 3), f(0), f(3, 2)),
    "e^c": Field("e^c", f(1), f(0), f(-1), f(3, 2)),
}


@dataclass(frozen=True)
class Operator:
    name: str
    fields: tuple[str, ...]
    kind: str


OPERATORS = [
    Operator("O_QQQL", ("Q", "Q", "Q", "L"), "dimension-6 gauge / dimension-5 SUSY superpotential core"),
    Operator("O_UUDE", ("u^c", "u^c", "d^c", "e^c"), "dimension-6 gauge / dimension-5 SUSY superpotential core"),
]


@dataclass(frozen=True)
class HadronicConstants:
    m_p_GeV: float = 0.938272
    f_pi_GeV: float = 0.130
    beta_H_GeV3: float = 0.012
    D: float = 0.80
    F: float = 0.47
    A_R: float = 2.4

    @property
    def width_prefactor_GeV5(self) -> float:
        chiral = (1.0 + self.D + self.F) ** 2
        return (
            self.m_p_GeV
            / (64.0 * math.pi * self.f_pi_GeV**2)
            * self.beta_H_GeV3**2
            * chiral
            * self.A_R**2
        )


def load_yukawa_singulars(sector: str) -> list[float]:
    payload = json.loads(YUKAWA_JSON.read_text(encoding="utf-8"))
    return payload["results"][sector]["normalized_singular_values"]


def operator_quantum_numbers(op: Operator) -> dict[str, object]:
    y = sum((FIELDS[name].y for name in op.fields), f(0))
    b = sum((FIELDS[name].b for name in op.fields), f(0))
    l = sum((FIELDS[name].l for name in op.fields), f(0))
    dim = sum((FIELDS[name].dimension for name in op.fields), f(0))
    return {
        "name": op.name,
        "fields": " ".join(op.fields),
        "kind": op.kind,
        "hypercharge": frac_str(y),
        "B": frac_str(b),
        "L": frac_str(l),
        "B_minus_L": frac_str(b - l),
        "fermion_operator_dimension": frac_str(dim),
        "hypercharge_zero": y == 0,
        "B_minus_L_zero": b - l == 0,
        "four_fermion_dimension_6": dim == 6,
    }


def lifetime_years_from_c6(c6_GeV_minus2: float, constants: HadronicConstants) -> float:
    width = constants.width_prefactor_GeV5 * c6_GeV_minus2**2
    if width <= 0.0:
        return math.inf
    return HBAR_GEV_S / width / SECONDS_PER_YEAR


def c6_gauge(alpha_G_inv: float, M_X_GeV: float, flavor_factor: float = 1.0) -> float:
    alpha = 1.0 / alpha_G_inv
    g2 = 4.0 * math.pi * alpha
    return flavor_factor * g2 / M_X_GeV**2


def gauge_lifetime_table(constants: HadronicConstants) -> list[dict[str, float]]:
    rows = []
    for alpha_inv in [24.0, 30.0, 40.0]:
        for mx in [1.0e15, 3.0e15, 1.0e16, 2.0e16, 5.0e16]:
            c6 = c6_gauge(alpha_inv, mx)
            rows.append(
                {
                    "alpha_G_inv": alpha_inv,
                    "M_X_GeV": mx,
                    "C6_GeV_minus2": c6,
                    "tau_years": lifetime_years_from_c6(c6, constants),
                }
            )
    return rows


def mx_required_for_tau(alpha_G_inv: float, tau_years: float, constants: HadronicConstants, flavor_factor: float = 1.0) -> float:
    alpha = 1.0 / alpha_G_inv
    g2 = 4.0 * math.pi * alpha
    width_target = HBAR_GEV_S / (tau_years * SECONDS_PER_YEAR)
    c6_max = math.sqrt(width_target / constants.width_prefactor_GeV5)
    return math.sqrt(flavor_factor * g2 / c6_max)


def dimension5_scenarios(constants: HadronicConstants) -> list[dict[str, float]]:
    up = load_yukawa_singulars("up")
    down = load_yukawa_singulars("down")

    y_top = 0.60
    y_bottom = 0.024
    y_u1 = y_top * up[2]
    y_c = y_top * up[1]
    y_t = y_top * up[0]
    y_d1 = y_bottom * down[2]
    y_s = y_bottom * down[1]
    y_b = y_bottom * down[0]

    alpha2 = 1.0 / 25.0
    loop = alpha2 / (4.0 * math.pi)
    M_T = 1.0e16
    m_wino = 1.0e3
    tau_benchmark = 1.0e34

    scenarios = [
        ("rank_staircase_first_gen", y_u1 * y_d1, 1.0, 1.0e5),
        ("second_gen_leakage", y_c * y_s, 1.0, 1.0e5),
        ("third_gen_unsuppressed", y_t * y_b, 1.0, 1.0e5),
        ("third_gen_with_triplet_geometric_filter_1e-4", y_t * y_b, 1.0e-4, 1.0e5),
        ("third_gen_with_triplet_geometric_filter_1e-5", y_t * y_b, 1.0e-5, 1.0e5),
        ("third_gen_with_10TeV_sfermions_filter_1e-4", y_t * y_b, 1.0e-4, 1.0e4),
    ]

    rows = []
    for name, yprod, triplet_filter, m_sfermion in scenarios:
        c5 = triplet_filter * yprod / M_T
        dressing = loop * m_wino / m_sfermion**2
        c6 = c5 * dressing
        tau = lifetime_years_from_c6(c6, constants)
        c6_max = math.sqrt(
            HBAR_GEV_S / (tau_benchmark * SECONDS_PER_YEAR) / constants.width_prefactor_GeV5
        )
        filter_required = c6_max / max((yprod / M_T) * dressing, 1e-300)
        rows.append(
            {
                "scenario": name,
                "yprod": yprod,
                "triplet_filter": triplet_filter,
                "M_T_GeV": M_T,
                "m_wino_GeV": m_wino,
                "m_sfermion_GeV": m_sfermion,
                "loop_alpha2_over_4pi": loop,
                "C5_GeV_minus1": c5,
                "dressing_GeV_minus1": dressing,
                "C6_dressed_GeV_minus2": c6,
                "tau_years": tau,
                "triplet_filter_required_for_tau_1e34": filter_required,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, object]) -> None:
    constants = payload["hadronic_constants"]
    lines: list[str] = []
    lines.append("# Item 4 proton-decay operator verification")
    lines.append("")
    lines.append("This report separates dimension-6 gauge exchange from model-dependent")
    lines.append("dimension-5 colored-triplet operators.  No web lookup was used; the")
    lines.append("lifetime target `1e34 yr` is a replaceable benchmark.")
    lines.append("")
    lines.append("## Operator ledger")
    lines.append("")
    lines.append("| operator | fields | Y | B | L | B-L | dim | checks |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|")
    for row in payload["operator_checks"]:
        checks = []
        if row["hypercharge_zero"]:
            checks.append("Y=0")
        if row["B_minus_L_zero"]:
            checks.append("B-L=0")
        if row["four_fermion_dimension_6"]:
            checks.append("dim6")
        lines.append(
            f"| `{row['name']}` | `{row['fields']}` | `{row['hypercharge']}` | `{row['B']}` | "
            f"`{row['L']}` | `{row['B_minus_L']}` | `{row['fermion_operator_dimension']}` | "
            f"{', '.join(checks)} |"
        )
    lines.append("")
    lines.append("## Dimension-6 gauge exchange")
    lines.append("")
    lines.append("Integrating out a heavy vector `X` with coupling `g_G X_mu J_X^mu` gives")
    lines.append("")
    lines.append("```text")
    lines.append("L_eff = - g_G^2 J_X^dagger J_X / M_X^2 + O(M_X^-4).")
    lines.append("```")
    lines.append("")
    lines.append("The estimate uses")
    lines.append("")
    lines.append("```text")
    lines.append("Gamma ~= K |C6|^2,  C6 = g_G^2/M_X^2")
    lines.append(f"K = {constants['width_prefactor_GeV5']:.6e} GeV^5")
    lines.append("```")
    lines.append("")
    lines.append("Representative lifetimes:")
    lines.append("")
    lines.append("| alpha_G^-1 | M_X [GeV] | tau [yr] |")
    lines.append("|---:|---:|---:|")
    for row in payload["gauge_lifetime_table"]:
        if row["alpha_G_inv"] in (24.0, 40.0) and row["M_X_GeV"] in (1e15, 3e15, 1e16, 2e16):
            lines.append(f"| {row['alpha_G_inv']:.0f} | {row['M_X_GeV']:.3e} | {row['tau_years']:.3e} |")
    lines.append("")
    lines.append("Mass required for the benchmark `tau > 1e34 yr`:")
    lines.append("")
    for row in payload["mx_required_for_tau_1e34"]:
        lines.append(f"- `alpha_G^-1={row['alpha_G_inv']:.0f}`: `M_X > {row['M_X_required_GeV']:.3e} GeV`")
    lines.append("")
    lines.append("## Dimension-5 colored-triplet operators")
    lines.append("")
    lines.append("In a supersymmetric or holomorphic completion, triplet exchange gives")
    lines.append("")
    lines.append("```text")
    lines.append("W5 = (C_L/M_T) QQQL + (C_R/M_T) u^c e^c u^c d^c.")
    lines.append("C6_dressed ~= (alpha_2/4pi)(m_wino/m_sfermion^2)(S_T y_a y_b/M_T).")
    lines.append("```")
    lines.append("")
    lines.append("Here `S_T` is the proposed CP1/O(2) triplet geometric filter.")
    lines.append("")
    lines.append("| scenario | yprod | S_T | m_sfermion [GeV] | tau [yr] | S_T required for 1e34 yr |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in payload["dimension5_scenarios"]:
        lines.append(
            f"| `{row['scenario']}` | {row['yprod']:.3e} | {row['triplet_filter']:.1e} | "
            f"{row['m_sfermion_GeV']:.1e} | {row['tau_years']:.3e} | "
            f"{row['triplet_filter_required_for_tau_1e34']:.3e} |"
        )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- Dimension-6 decay is controlled mostly by `M_X`; it cannot be removed by Yukawa texture.")
    lines.append("- Dimension-5 decay is model-dependent and can be killed by absence of low triplets, heavy sfermions, or the CP1/O(2) triplet filter.")
    lines.append("- The dangerous operator ledger preserves `B-L`; the item-3 `B-L=2` Majorana trace-lift does not by itself induce proton decay.")
    lines.append("- Item 5 must lock threshold corrections to these proton constraints, because lowering `M_X` to fix unification can immediately revive dimension-6 decay.")
    (OUT / "proton_decay_item4_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    constants = HadronicConstants()
    operator_checks = [operator_quantum_numbers(op) for op in OPERATORS]
    gauge_rows = gauge_lifetime_table(constants)
    d5_rows = dimension5_scenarios(constants)
    tau_target = 1.0e34
    mx_rows = [
        {
            "alpha_G_inv": alpha_inv,
            "tau_target_years": tau_target,
            "M_X_required_GeV": mx_required_for_tau(alpha_inv, tau_target, constants),
        }
        for alpha_inv in [24.0, 30.0, 40.0]
    ]

    write_csv(OUT / "dimension6_gauge_lifetimes.csv", gauge_rows)
    write_csv(OUT / "dimension5_triplet_lifetimes.csv", d5_rows)

    checks = {
        "operator_hypercharge_zero": all(row["hypercharge_zero"] for row in operator_checks),
        "operator_B_minus_L_zero": all(row["B_minus_L_zero"] for row in operator_checks),
        "operator_dimension_6": all(row["four_fermion_dimension_6"] for row in operator_checks),
        "gauge_MX_2e16_alpha24_tau_gt_1e34": next(
            row["tau_years"]
            for row in gauge_rows
            if row["alpha_G_inv"] == 24.0 and row["M_X_GeV"] == 2.0e16
        )
        > tau_target,
        "gauge_MX_3e15_alpha24_tau_lt_1e34": next(
            row["tau_years"]
            for row in gauge_rows
            if row["alpha_G_inv"] == 24.0 and row["M_X_GeV"] == 3.0e15
        )
        < tau_target,
        "dimension5_filter_1e_minus5_can_restore_tau": next(
            row["tau_years"]
            for row in d5_rows
            if row["scenario"] == "third_gen_with_triplet_geometric_filter_1e-5"
        )
        > tau_target,
    }

    payload = {
        "input": {
            "yukawa_json": str(YUKAWA_JSON),
            "note": "Approximate local constants; no web lookup used.",
        },
        "hadronic_constants": {
            **asdict(constants),
            "width_prefactor_GeV5": constants.width_prefactor_GeV5,
        },
        "operator_checks": operator_checks,
        "gauge_lifetime_table": gauge_rows,
        "mx_required_for_tau_1e34": mx_rows,
        "dimension5_scenarios": d5_rows,
        "checks": checks,
        "all_checks_ok": all(checks.values()),
    }
    (OUT / "proton_decay_verification.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)

    print("Item 4 proton-decay verification")
    print(f"  width prefactor K: {constants.width_prefactor_GeV5:.6e} GeV^5")
    print("  operator checks:")
    for key, value in checks.items():
        print(f"    {key}: {value}")
    print("  M_X required for tau > 1e34 yr:")
    for row in mx_rows:
        print(f"    alpha_G^-1={row['alpha_G_inv']:.0f}: {row['M_X_required_GeV']:.6e} GeV")
    print("  selected d=5 scenarios:")
    for row in d5_rows:
        print(
            "    {scenario}: tau={tau:.3e} yr, S_required={sreq:.3e}".format(
                scenario=row["scenario"],
                tau=row["tau_years"],
                sreq=row["triplet_filter_required_for_tau_1e34"],
            )
        )
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")

    if not payload["all_checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
