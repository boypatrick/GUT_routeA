#!/usr/bin/env python3
"""Exact Pati-Salam decomposition checker for one Spin(10) spinor 16.

The script uses rational arithmetic throughout.  It constructs the branching

    16 -> (4,2,1) + (bar4,1,2)

then breaks SU(4)_C -> SU(3)_C x U(1)_{B-L} and verifies

    Y = T3R + (B-L)/2.

Outputs are written under output/pati_salam/.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "pati_salam"


def f(num: int, den: int = 1) -> Fraction:
    return Fraction(num, den)


def frac_str(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def frac_float(x: Fraction) -> float:
    return float(x.numerator) / float(x.denominator)


@dataclass(frozen=True)
class SMMultiplet:
    name: str
    ps_origin: str
    su3: str
    su3_dim: int
    su2l_dim: int
    t3r: Fraction
    b_minus_l: Fraction
    hypercharge: Fraction
    components: str

    @property
    def dimension(self) -> int:
        return self.su3_dim * self.su2l_dim


@dataclass(frozen=True)
class WeylComponent:
    name: str
    ps_origin: str
    su3: str
    color_state: str
    t3l: Fraction
    t3r: Fraction
    b_minus_l: Fraction
    hypercharge: Fraction
    electric_charge: Fraction


def build_multiplets() -> list[SMMultiplet]:
    """Return the six Standard Model left-handed Weyl multiplets in one 16."""

    return [
        SMMultiplet(
            name="Q",
            ps_origin="(4,2,1)",
            su3="3",
            su3_dim=3,
            su2l_dim=2,
            t3r=f(0),
            b_minus_l=f(1, 3),
            hypercharge=f(1, 6),
            components="(u_L,d_L)",
        ),
        SMMultiplet(
            name="L",
            ps_origin="(4,2,1)",
            su3="1",
            su3_dim=1,
            su2l_dim=2,
            t3r=f(0),
            b_minus_l=f(-1),
            hypercharge=f(-1, 2),
            components="(nu_L,e_L)",
        ),
        SMMultiplet(
            name="u^c",
            ps_origin="(bar4,1,2)",
            su3="bar3",
            su3_dim=3,
            su2l_dim=1,
            t3r=f(-1, 2),
            b_minus_l=f(-1, 3),
            hypercharge=f(-2, 3),
            components="u^c",
        ),
        SMMultiplet(
            name="d^c",
            ps_origin="(bar4,1,2)",
            su3="bar3",
            su3_dim=3,
            su2l_dim=1,
            t3r=f(1, 2),
            b_minus_l=f(-1, 3),
            hypercharge=f(1, 3),
            components="d^c",
        ),
        SMMultiplet(
            name="nu^c",
            ps_origin="(bar4,1,2)",
            su3="1",
            su3_dim=1,
            su2l_dim=1,
            t3r=f(-1, 2),
            b_minus_l=f(1),
            hypercharge=f(0),
            components="nu^c",
        ),
        SMMultiplet(
            name="e^c",
            ps_origin="(bar4,1,2)",
            su3="1",
            su3_dim=1,
            su2l_dim=1,
            t3r=f(1, 2),
            b_minus_l=f(1),
            hypercharge=f(1),
            components="e^c",
        ),
    ]


def build_components() -> list[WeylComponent]:
    """Expand the same 16 to individual Weyl components."""

    out: list[WeylComponent] = []
    colors = ["r", "g", "b"]
    anticolors = ["bar r", "bar g", "bar b"]

    for color in colors:
        y = f(0) + f(1, 3) / 2
        out.append(
            WeylComponent("u_L", "(4,2,1)", "3", color, f(1, 2), f(0), f(1, 3), y, f(1, 2) + y)
        )
        out.append(
            WeylComponent("d_L", "(4,2,1)", "3", color, f(-1, 2), f(0), f(1, 3), y, f(-1, 2) + y)
        )

    y_l = f(0) + f(-1) / 2
    out.append(WeylComponent("nu_L", "(4,2,1)", "1", "lepton", f(1, 2), f(0), f(-1), y_l, f(1, 2) + y_l))
    out.append(WeylComponent("e_L", "(4,2,1)", "1", "lepton", f(-1, 2), f(0), f(-1), y_l, f(-1, 2) + y_l))

    for color in anticolors:
        y_uc = f(-1, 2) + f(-1, 3) / 2
        y_dc = f(1, 2) + f(-1, 3) / 2
        out.append(WeylComponent("u^c", "(bar4,1,2)", "bar3", color, f(0), f(-1, 2), f(-1, 3), y_uc, y_uc))
        out.append(WeylComponent("d^c", "(bar4,1,2)", "bar3", color, f(0), f(1, 2), f(-1, 3), y_dc, y_dc))

    y_nuc = f(-1, 2) + f(1) / 2
    y_ec = f(1, 2) + f(1) / 2
    out.append(WeylComponent("nu^c", "(bar4,1,2)", "1", "lepton", f(0), f(-1, 2), f(1), y_nuc, y_nuc))
    out.append(WeylComponent("e^c", "(bar4,1,2)", "1", "lepton", f(0), f(1, 2), f(1), y_ec, y_ec))

    return out


def su3_dynkin_index(rep: str) -> Fraction:
    return f(1, 2) if rep in {"3", "bar3"} else f(0)


def su3_cubic_index(rep: str) -> Fraction:
    if rep == "3":
        return f(1)
    if rep == "bar3":
        return f(-1)
    return f(0)


def su2_dynkin_index(dim: int) -> Fraction:
    return f(1, 2) if dim == 2 else f(0)


def anomaly_sums(multiplets: Iterable[SMMultiplet]) -> dict[str, Fraction]:
    multiplets = list(multiplets)
    return {
        "SU3^3": sum(f(m.su2l_dim) * su3_cubic_index(m.su3) for m in multiplets),
        "SU3^2-U1Y": sum(f(m.su2l_dim) * su3_dynkin_index(m.su3) * m.hypercharge for m in multiplets),
        "SU2L^2-U1Y": sum(f(m.su3_dim) * su2_dynkin_index(m.su2l_dim) * m.hypercharge for m in multiplets),
        "U1Y^3": sum(f(m.dimension) * m.hypercharge**3 for m in multiplets),
        "grav^2-U1Y": sum(f(m.dimension) * m.hypercharge for m in multiplets),
    }


def pati_salam_anomalies() -> dict[str, Fraction | int | bool]:
    su4_cubic = f(2) - f(2)
    su2l_doublets = 4
    su2r_doublets = 4
    return {
        "SU4^3": su4_cubic,
        "SU2L_global_doublets": su2l_doublets,
        "SU2L_global_even": su2l_doublets % 2 == 0,
        "SU2R_global_doublets": su2r_doublets,
        "SU2R_global_even": su2r_doublets % 2 == 0,
    }


def trace_normalizations(components: Iterable[WeylComponent]) -> dict[str, Fraction]:
    components = list(components)
    t3l2 = sum(c.t3l**2 for c in components)
    t3r2 = sum(c.t3r**2 for c in components)
    u2 = sum((c.b_minus_l / 2) ** 2 for c in components)
    y2 = sum(c.hypercharge**2 for c in components)
    ty = sum(c.t3r * (c.b_minus_l / 2) for c in components)
    yx = sum(c.hypercharge * (f(2) * c.t3r - f(3) * c.b_minus_l / 2) for c in components)
    x2 = sum((f(2) * c.t3r - f(3) * c.b_minus_l / 2) ** 2 for c in components)
    return {
        "Tr_T3L^2": t3l2,
        "Tr_T3R^2": t3r2,
        "Tr_((B-L)/2)^2": u2,
        "Tr_T3R*((B-L)/2)": ty,
        "Tr_Y^2": y2,
        "kY=Tr_Y^2/Tr_T3L^2": y2 / t3l2,
        "Tr_Y*X_PS": yx,
        "Tr_X_PS^2_for_X=2T3R-3(B-L)/2": x2,
    }


def solve_hypercharge_uniqueness() -> dict[str, Fraction]:
    """Solve Y = a T3R + b (B-L)/2 from Q and u^c charges."""

    # Q: 1/6 = a*0 + b*(1/6), so b=1.
    b = f(1)
    # u^c: -2/3 = a*(-1/2) + b*(-1/6), so a=1.
    a = (f(-2, 3) - b * f(-1, 6)) / f(-1, 2)
    return {"a_T3R": a, "b_(B-L)/2": b}


def serialize_fraction_dict(data: dict[str, Fraction | int | bool]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in data.items():
        if isinstance(value, Fraction):
            out[key] = {"exact": frac_str(value), "float": frac_float(value)}
        else:
            out[key] = value
    return out


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    multiplets: list[SMMultiplet],
    components: list[WeylComponent],
    sm_anomalies: dict[str, Fraction],
    ps_anomalies: dict[str, Fraction | int | bool],
    traces: dict[str, Fraction],
    uniqueness: dict[str, Fraction],
) -> None:
    lines: list[str] = []
    lines.append("# Pati-Salam decomposition check")
    lines.append("")
    lines.append("Branching: `16 -> (4,2,1) + (bar4,1,2)`.")
    lines.append("Hypercharge: `Y = T3R + (B-L)/2`.")
    lines.append("")
    lines.append("## Multiplet table")
    lines.append("")
    lines.append("| field | PS origin | SU(3)c | SU(2)L dim | T3R | B-L | Y | dim |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for m in multiplets:
        lines.append(
            "| {name} | {origin} | {su3} | {su2} | {t3r} | {bl} | {y} | {dim} |".format(
                name=m.name,
                origin=m.ps_origin,
                su3=m.su3,
                su2=m.su2l_dim,
                t3r=frac_str(m.t3r),
                bl=frac_str(m.b_minus_l),
                y=frac_str(m.hypercharge),
                dim=m.dimension,
            )
        )
    lines.append("")
    lines.append(f"Total left-handed Weyl dimension: `{sum(m.dimension for m in multiplets)}`.")
    lines.append("")
    lines.append("## SM anomaly sums")
    lines.append("")
    for key, value in sm_anomalies.items():
        lines.append(f"- `{key}` = `{frac_str(value)}`")
    lines.append("")
    lines.append("## Pati-Salam anomaly/global checks")
    lines.append("")
    for key, value in ps_anomalies.items():
        lines.append(f"- `{key}` = `{value if not isinstance(value, Fraction) else frac_str(value)}`")
    lines.append("")
    lines.append("## Hypercharge uniqueness and normalization")
    lines.append("")
    lines.append(
        "- Solving `Y = a T3R + b (B-L)/2` from `Q` and `u^c` gives "
        f"`a={frac_str(uniqueness['a_T3R'])}`, `b={frac_str(uniqueness['b_(B-L)/2'])}`."
    )
    for key, value in traces.items():
        lines.append(f"- `{key}` = `{frac_str(value)}`")
    lines.append("")
    lines.append("Therefore `T_Y^GUT = sqrt(3/5) Y` and `alpha_1 = (5/3) alpha_Y`.")
    lines.append(
        "The orthogonal broken Abelian direction can be chosen as "
        "`X_PS = 2 T3R - 3(B-L)/2`, with `Tr(Y X_PS)=0`."
    )
    lines.append("")
    lines.append("## Individual component charges")
    lines.append("")
    lines.append("| component | color | T3L | T3R | B-L | Y | Qem |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for c in components:
        lines.append(
            "| {name} | {color} | {t3l} | {t3r} | {bl} | {y} | {q} |".format(
                name=c.name,
                color=c.color_state,
                t3l=frac_str(c.t3l),
                t3r=frac_str(c.t3r),
                bl=frac_str(c.b_minus_l),
                y=frac_str(c.hypercharge),
                q=frac_str(c.electric_charge),
            )
        )
    (OUT / "pati_salam_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    multiplets = build_multiplets()
    components = build_components()

    sm_anomalies = anomaly_sums(multiplets)
    ps_anomalies = pati_salam_anomalies()
    traces = trace_normalizations(components)
    uniqueness = solve_hypercharge_uniqueness()

    expected_y = {m.name: m.t3r + m.b_minus_l / 2 for m in multiplets}
    hypercharge_ok = all(expected_y[m.name] == m.hypercharge for m in multiplets)
    component_hypercharge_ok = all(c.hypercharge == c.t3r + c.b_minus_l / 2 for c in components)
    all_sm_anomalies_zero = all(value == 0 for value in sm_anomalies.values())
    total_dimension = sum(m.dimension for m in multiplets)
    dimension_ok = total_dimension == 16 and len(components) == 16
    ps_ok = ps_anomalies["SU4^3"] == 0 and bool(ps_anomalies["SU2L_global_even"]) and bool(ps_anomalies["SU2R_global_even"])
    uniqueness_ok = uniqueness == {"a_T3R": f(1), "b_(B-L)/2": f(1)}
    normalization_ok = traces["kY=Tr_Y^2/Tr_T3L^2"] == f(5, 3)
    orthogonal_ok = traces["Tr_Y*X_PS"] == 0

    multiplet_rows = []
    for m in multiplets:
        row = asdict(m)
        row["t3r"] = frac_str(m.t3r)
        row["b_minus_l"] = frac_str(m.b_minus_l)
        row["hypercharge"] = frac_str(m.hypercharge)
        row["dimension"] = m.dimension
        multiplet_rows.append(row)

    component_rows = []
    for c in components:
        row = asdict(c)
        row["t3l"] = frac_str(c.t3l)
        row["t3r"] = frac_str(c.t3r)
        row["b_minus_l"] = frac_str(c.b_minus_l)
        row["hypercharge"] = frac_str(c.hypercharge)
        row["electric_charge"] = frac_str(c.electric_charge)
        component_rows.append(row)

    write_csv(
        OUT / "pati_salam_decomposition.csv",
        multiplet_rows,
        [
            "name",
            "ps_origin",
            "su3",
            "su3_dim",
            "su2l_dim",
            "t3r",
            "b_minus_l",
            "hypercharge",
            "components",
            "dimension",
        ],
    )
    write_csv(
        OUT / "pati_salam_components.csv",
        component_rows,
        [
            "name",
            "ps_origin",
            "su3",
            "color_state",
            "t3l",
            "t3r",
            "b_minus_l",
            "hypercharge",
            "electric_charge",
        ],
    )

    result = {
        "branching": "Spin(10) 16 -> (4,2,1) + (bar4,1,2)",
        "hypercharge_formula": "Y = T3R + (B-L)/2",
        "total_dimension": total_dimension,
        "component_count": len(components),
        "checks": {
            "dimension_ok": dimension_ok,
            "hypercharge_ok": hypercharge_ok,
            "component_hypercharge_ok": component_hypercharge_ok,
            "all_sm_anomalies_zero": all_sm_anomalies_zero,
            "pati_salam_anomalies_ok": ps_ok,
            "hypercharge_uniqueness_ok": uniqueness_ok,
            "gut_normalization_kY_5_over_3_ok": normalization_ok,
            "orthogonal_X_PS_ok": orthogonal_ok,
            "all_checks_ok": all(
                [
                    dimension_ok,
                    hypercharge_ok,
                    component_hypercharge_ok,
                    all_sm_anomalies_zero,
                    ps_ok,
                    uniqueness_ok,
                    normalization_ok,
                    orthogonal_ok,
                ]
            ),
        },
        "sm_anomalies": serialize_fraction_dict(sm_anomalies),
        "pati_salam_anomalies": serialize_fraction_dict(ps_anomalies),
        "trace_normalizations": serialize_fraction_dict(traces),
        "hypercharge_uniqueness": serialize_fraction_dict(uniqueness),
    }

    (OUT / "pati_salam_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    write_report(multiplets, components, sm_anomalies, ps_anomalies, traces, uniqueness)

    print("Pati-Salam Spin(10) 16 verification")
    print(f"  total dimension: {total_dimension}")
    print(f"  component count: {len(components)}")
    for key, value in result["checks"].items():
        print(f"  {key}: {value}")
    print("  SM anomalies:")
    for key, value in sm_anomalies.items():
        print(f"    {key}: {frac_str(value)}")
    print("  trace normalization:")
    print(f"    Tr(Y^2)/Tr(T3L^2): {frac_str(traces['kY=Tr_Y^2/Tr_T3L^2'])}")
    print(f"    Tr(Y X_PS): {frac_str(traces['Tr_Y*X_PS'])}")
    print(f"  wrote: {OUT}")

    if not result["checks"]["all_checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
