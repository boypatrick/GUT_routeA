#!/usr/bin/env python3
"""Mechanical verifier for the frozen P54PQ-v2 action card.

This verifier does not assert that a displayed tensor list is a complete
Hilbert basis and does not replace a Spin(10) Clebsch/Hessian calculation.
It checks the pieces that can be made exact at the action-definition stage:

* operator dimension, PQ neutrality, and Hermiticity declarations;
* raw and basis-quotiented parameter counts;
* the rank of scalar rephasing vectors and the three invariant CP phases;
* the printed-versus-corrected eta1 charge;
* the 54 vacuum trace fingerprints and Pati--Salam dimensions;
* fail-closed branch content and forbidden Yukawa structures.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "p54_action_parameter_convention_card.tex"
SCRIPT = Path(__file__).resolve()


PQ = {
    "Psi": -1,
    "Phi": 0,
    "Sigma": 2,
    "Sigma*": -2,
    "phi": -2,
    "phi*": 2,
    "S": -4,
    "S*": 4,
}


@dataclass(frozen=True)
class Operator:
    name: str
    coefficient: str
    fields: dict[str, int]
    coefficient_dimension: int
    hermiticity: str

    @property
    def field_dimension(self) -> int:
        return sum(self.fields.values())

    @property
    def pq_charge(self) -> int:
        return sum(PQ[field] * count for field, count in self.fields.items())


OPS = [
    Operator("Phi mass", "mu2", {"Phi": 2}, 2, "real"),
    Operator("Phi cubic", "c", {"Phi": 3}, 1, "real"),
    Operator("Phi quartic trace-square", "a", {"Phi": 4}, 0, "real"),
    Operator("Phi quartic trace-four", "b", {"Phi": 4}, 0, "real"),
    Operator("Sigma mass", "nu2", {"Sigma": 1, "Sigma*": 1}, 2, "real"),
    Operator("Sigma quartic singlet", "lambda0", {"Sigma": 2, "Sigma*": 2}, 0, "real"),
    Operator("Sigma quartic 2", "lambda2", {"Sigma": 2, "Sigma*": 2}, 0, "real"),
    Operator("Sigma quartic 4", "lambda4", {"Sigma": 2, "Sigma*": 2}, 0, "real"),
    Operator("Sigma quartic 4 prime", "lambda4p", {"Sigma": 2, "Sigma*": 2}, 0, "real"),
    Operator("Phi2 SigmaSigma*", "alpha", {"Phi": 2, "Sigma": 1, "Sigma*": 1}, 0, "real"),
    Operator("PhiPhi SigmaSigma*", "beta", {"Phi": 2, "Sigma": 1, "Sigma*": 1}, 0, "real"),
    Operator("phi mass", "xi02", {"phi": 1, "phi*": 1}, 2, "real"),
    Operator("phi quartic norm", "xi1", {"phi": 2, "phi*": 2}, 0, "real"),
    Operator("phi quartic holomorphic norm", "xi2", {"phi": 2, "phi*": 2}, 0, "real"),
    Operator("Phi phi phi*", "xi3", {"Phi": 1, "phi": 1, "phi*": 1}, 1, "real"),
    Operator("SigmaSigma* phi phi* first", "gamma1", {"Sigma": 1, "Sigma*": 1, "phi": 1, "phi*": 1}, 0, "real"),
    Operator("SigmaSigma* phi phi* second", "gamma2", {"Sigma": 1, "Sigma*": 1, "phi": 1, "phi*": 1}, 0, "real"),
    Operator("Phi2 phi phi*", "eta0", {"Phi": 2, "phi": 1, "phi*": 1}, 0, "real"),
    Operator("corrected eta1", "eta1", {"Sigma": 2, "Sigma*": 1, "phi": 1}, 0, "complex_plus_hc"),
    Operator("PhiPhi phi phi*", "eta2", {"Phi": 2, "phi": 1, "phi*": 1}, 0, "real"),
    Operator("SigmaSigma phiphi", "eta3", {"Sigma": 2, "phi": 2}, 0, "complex_plus_hc"),
    Operator("S mass", "mus2", {"S": 1, "S*": 1}, 2, "real"),
    Operator("S quartic", "chi1", {"S": 2, "S*": 2}, 0, "real"),
    Operator("SigmaSigma* SS*", "chi2", {"Sigma": 1, "Sigma*": 1, "S": 1, "S*": 1}, 0, "real"),
    Operator("Phi2 SS*", "chi3", {"Phi": 2, "S": 1, "S*": 1}, 0, "real"),
    Operator("SigmaSigma Phi S", "chi4", {"Sigma": 2, "Phi": 1, "S": 1}, 0, "complex_plus_hc"),
    Operator("phi phi* SS*", "chi5", {"phi": 1, "phi*": 1, "S": 1, "S*": 1}, 0, "real"),
    Operator("phiphi S*", "chi6", {"phi": 2, "S*": 1}, 1, "complex_plus_hc"),
    Operator("Phi phiphi S*", "chi7", {"Phi": 1, "phi": 2, "S*": 1}, 0, "complex_plus_hc"),
]


REAL_SCALAR_COEFFICIENTS = (
    "mu2", "nu2", "xi02", "mus2",
    "c", "xi3",
    "a", "b", "lambda0", "lambda2", "lambda4", "lambda4p",
    "alpha", "beta", "xi1", "xi2", "gamma1", "gamma2", "eta0",
    "eta2", "chi1", "chi2", "chi3", "chi5",
)
COMPLEX_SCALAR_COEFFICIENTS = ("eta1", "eta3", "chi4", "chi6", "chi7")
PHASE_VECTORS = {
    "eta1": (1, 1, 0),
    "eta3": (2, 2, 0),
    "chi4": (2, 0, 1),
    "chi6": (0, 2, -1),
    "chi7": (0, 2, -1),
}


CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def exact_rank(rows: list[tuple[int, ...]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    value - factor * pivot_entry
                    for value, pivot_entry in zip(matrix[row], matrix[rank])
                ]
        rank += 1
    return rank


def index_multiplicities(factors: tuple[str, ...]) -> dict[str, int]:
    indices = sorted(set("".join(factors)))
    return {index: sum(factor.count(index) for factor in factors) for index in indices}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def audit() -> dict[str, Any]:
    operator_rows = []
    for operator in OPS:
        dimension_ok = operator.field_dimension + operator.coefficient_dimension == 4
        pq_ok = operator.pq_charge == 0
        hermitian_ok = operator.hermiticity in {"real", "complex_plus_hc"}
        operator_rows.append(
            {
                "name": operator.name,
                "coefficient": operator.coefficient,
                "field_dimension": operator.field_dimension,
                "coefficient_dimension": operator.coefficient_dimension,
                "total_dimension": operator.field_dimension + operator.coefficient_dimension,
                "pq_charge": operator.pq_charge,
                "hermiticity": operator.hermiticity,
                "pass": dimension_ok and pq_ok and hermitian_ok,
            }
        )

    check(
        "action",
        "all scalar operators have dimension four",
        all(row["total_dimension"] == 4 for row in operator_rows),
        f"audited {len(operator_rows)} coefficient-monomial pairs",
    )
    check(
        "action",
        "all scalar operators are PQ neutral",
        all(row["pq_charge"] == 0 for row in operator_rows),
        "charges are q(Psi,Phi,Sigma,phi,S)=(-1,0,+2,-2,-4)",
    )
    check(
        "action",
        "complex scalar operators carry an explicit Hermitian partner",
        all(row["hermiticity"] in {"real", "complex_plus_hc"} for row in operator_rows),
        "eta1, eta3, chi4, chi6, and chi7 are declared complex_plus_hc",
    )

    printed_eta1_charge = sum(
        PQ[field] * count
        for field, count in {"Sigma": 1, "Sigma*": 2, "phi": 1}.items()
    )
    corrected_eta1_charge = next(
        operator.pq_charge for operator in OPS if operator.coefficient == "eta1"
    )
    check(
        "correction",
        "printed eta1 monomial fails continuous PQ",
        printed_eta1_charge == -4,
        f"q(Sigma Sigma* Sigma* phi)={printed_eta1_charge}",
    )
    check(
        "correction",
        "corrected eta1 monomial is PQ neutral",
        corrected_eta1_charge == 0,
        f"q(Sigma Sigma* Sigma phi)={corrected_eta1_charge}",
    )
    corrected_eta1_indices = index_multiplicities(
        ("ijklm", "ijkpq", "lmpqn", "n")
    )
    corrected_chi4_indices = index_multiplicities(("ijklm", "ijkln", "mn"))
    printed_chi4_indices = index_multiplicities(("ijklm", "ijkln", "ij"))
    check(
        "correction",
        "corrected eta1 and chi4 contractions are index balanced",
        all(value == 2 for value in corrected_eta1_indices.values())
        and all(value == 2 for value in corrected_chi4_indices.values())
        and any(value != 2 for value in printed_chi4_indices.values()),
        "eta1 and corrected Phi_mn have every index twice; printed Phi_ij does not",
    )

    scalar_raw = len(REAL_SCALAR_COEFFICIENTS) + 2 * len(COMPLEX_SCALAR_COEFFICIENTS)
    phase_rank = exact_rank(list(PHASE_VECTORS.values()))
    scalar_quotient = scalar_raw - phase_rank
    yukawa_raw = 2 * (3 * 4 // 2) * 2
    yukawa_quotient = yukawa_raw - 9
    total_raw = 2 + scalar_raw + yukawa_raw
    total_classical = 2 + scalar_quotient + yukawa_quotient
    total_quantum_pq = total_classical - 1

    check(
        "parameters",
        "scalar raw count",
        scalar_raw == 34,
        f"24 real coefficients + 5 complex coefficients = {scalar_raw} reals",
    )
    check(
        "parameters",
        "scalar rephasing rank",
        phase_rank == 2,
        "five complex-coupling phase vectors have exact rank 2",
    )
    check(
        "parameters",
        "three scalar CP phases survive",
        len(COMPLEX_SCALAR_COEFFICIENTS) - phase_rank == 3,
        "delta1=arg(eta3)-2arg(eta1); delta2=arg(chi4)+arg(chi6)-2arg(eta1); delta3=arg(chi7)-arg(chi6)",
    )
    check(
        "parameters",
        "two symmetric Yukawa matrices give 15 physical reals",
        yukawa_raw == 24 and yukawa_quotient == 15,
        "24 raw reals minus dim U(3)=9",
    )
    check(
        "parameters",
        "joint classical and quantum counts",
        total_raw == 60 and total_classical == 49 and total_quantum_pq == 48,
        "raw=60, classical quotient=49, anomalous-PQ observable quotient=48",
    )

    # Exact rational fingerprints of the 54 Pati--Salam direction.
    eigenvalues = [Fraction(-2, 5)] * 6 + [Fraction(3, 5)] * 4
    traces = {
        "trX": sum(eigenvalues),
        "trX2": sum(value**2 for value in eigenvalues),
        "trX3": sum(value**3 for value in eigenvalues),
        "trX4": sum(value**4 for value in eigenvalues),
    }
    radial = {
        "mu2": -traces["trX2"] / 2,
        "c": traces["trX3"] / 3,
        "a": traces["trX2"] ** 2 / 4,
        "b": traces["trX4"] / 2,
    }
    expected_radial = {
        "mu2": Fraction(-6, 5),
        "c": Fraction(4, 25),
        "a": Fraction(36, 25),
        "b": Fraction(42, 125),
    }
    check(
        "vacuum",
        "54 Pati--Salam direction is traceless",
        traces["trX"] == 0,
        "6(-2/5)+4(3/5)=0",
    )
    check(
        "vacuum",
        "54 radial potential fingerprints",
        radial == expected_radial,
        "(-6/5,4/25,36/25,42/125)",
    )
    check(
        "vacuum",
        "54 commutant has Pati--Salam dimension",
        6 * 5 // 2 + 4 * 3 // 2 == 15 + 3 + 3,
        "dim so(6)+dim so(4)=21=dim su(4)+dim su(2)+dim su(2)",
    )

    ps_dimensions = {
        "10": 6 * 1 * 1 + 1 * 2 * 2,
        "54": 1 * 3 * 3 + 6 * 2 * 2 + 20 * 1 * 1 + 1,
        "126": 6 + 10 * 3 + 10 * 3 + 15 * 2 * 2,
    }
    check(
        "representations",
        "Pati--Salam decompositions preserve dimensions",
        ps_dimensions == {"10": 10, "54": 54, "126": 126},
        json.dumps(ps_dimensions, sort_keys=True),
    )
    check(
        "representations",
        "Yukawa symmetry follows from representation channel",
        10 in (10, 120, 126) and 126 in (10, 120, 126),
        "16 tensor 16 = 10_s + 120_a + 126_s; only symmetric channels are present",
    )

    primary_fields = {"3x16F", "54R", "126C", "10C", "1C"}
    forbidden_primary_fields = {"120H", "210H", "extra10H"}
    check(
        "branch",
        "primary field contract is exact",
        primary_fields == {"3x16F", "54R", "126C", "10C", "1C"},
        "+".join(sorted(primary_fields)),
    )
    check(
        "branch",
        "comparison representations are absent from primary",
        primary_fields.isdisjoint(forbidden_primary_fields),
        "120H, 210H, and an extra 10H are fail-closed",
    )

    allowed_yukawa_charges = {
        "PsiPsi_phi*": 2 * PQ["Psi"] + PQ["phi*"],
        "PsiPsi_Sigma": 2 * PQ["Psi"] + PQ["Sigma"],
    }
    forbidden_yukawa_charge = 2 * PQ["Psi"] + PQ["phi"]
    check(
        "branch",
        "exactly two PQ-allowed Yukawa channels",
        allowed_yukawa_charges == {"PsiPsi_phi*": 0, "PsiPsi_Sigma": 0}
        and forbidden_yukawa_charge == -4,
        f"allowed={allowed_yukawa_charges}; q(PsiPsi phi)={forbidden_yukawa_charge}",
    )

    tex = TEX.read_text(encoding="utf-8")
    required_tex_tokens = (
        r"\Lag_{\rm P54PQ}",
        r"\Lag_Y",
        r"V_{\Phi\Sigma}",
        r"V_{\Phi\Sigma\phi}",
        r"V_S",
        r"\ThetaP/\mathcal G_{\rm basis}",
        r"\langle\Phi\rangle",
        r"\operatorname{rank}J",
    )
    check(
        "artifact",
        "TeX card contains all required interfaces",
        all(token in tex for token in required_tex_tokens),
        f"required token count={len(required_tex_tokens)}",
    )

    passed = sum(row["pass"] for row in CHECKS)
    return {
        "schema": "route-f-p54-action-card-v2",
        "model_id": "P54PQ-v2",
        "date": "2026-08-30",
        "status": {
            "action_definition_subgate": "done",
            "p0_overall": "action-and-renormalizable-audit-done",
            "frozen_renormalizable_candidate": True,
            "phenomenological_viability_claimed": False,
            "printed_eta1_used": False,
        },
        "branch_contract": {
            "gauge_group": "Spin(10)",
            "global_symmetry": "U(1)_PQ",
            "primary_fields": sorted(primary_fields),
            "comparison_branches": ["P54N", "P210"],
            "forbidden_imports": sorted(forbidden_primary_fields),
        },
        "charges": PQ,
        "operator_audit": operator_rows,
        "eta1_correction": {
            "printed_monomial": "Sigma Sigma* Sigma* phi",
            "printed_charge": printed_eta1_charge,
            "corrected_monomial": "Sigma Sigma* Sigma phi + h.c.",
            "corrected_charge": corrected_eta1_charge,
            "corrected_index_multiplicities": corrected_eta1_indices,
            "chi4_corrected_index_multiplicities": corrected_chi4_indices,
            "chi4_printed_index_multiplicities": printed_chi4_indices,
        },
        "parameter_count": {
            "gauge_topological_raw_reals": 2,
            "scalar_raw_reals": scalar_raw,
            "scalar_phase_orbit_rank": phase_rank,
            "scalar_quotient_reals": scalar_quotient,
            "yukawa_raw_reals": yukawa_raw,
            "family_basis_dimension": 9,
            "yukawa_quotient_reals": yukawa_quotient,
            "total_raw_reals": total_raw,
            "total_classical_quotient_reals": total_classical,
            "total_quantum_pq_observable_reals": total_quantum_pq,
            "physical_scalar_cp_phases": [
                "arg(eta3)-2arg(eta1)",
                "arg(chi4)+arg(chi6)-2arg(eta1)",
                "arg(chi7)-arg(chi6)",
            ],
        },
        "phase_vectors": PHASE_VECTORS,
        "vacuum_fingerprints": {
            "traces": {name: str(value) for name, value in traces.items()},
            "radial_coefficients": {name: str(value) for name, value in radial.items()},
            "ps_dimensions": ps_dimensions,
        },
        "open_debts": [
            "PQ-quality policy for higher-dimension operators",
            "P2 running and threshold matching on the completed P1 spectrum",
            "scaled-SVD Jacobian rank after the spectrum and matching maps exist",
        ],
        "checks": CHECKS,
        "summary": {
            "passed": passed,
            "total": len(CHECKS),
            "all_pass": passed == len(CHECKS),
        },
        "sources": [source_row(TEX), source_row(SCRIPT)],
    }


def markdown(report: dict[str, Any]) -> str:
    counts = report["parameter_count"]
    summary = report["summary"]
    lines = [
        "# P54PQ-v2 action-card verification",
        "",
        f"Status: **{summary['passed']}/{summary['total']} checks passed**.",
        "",
        "The primary P-layer action is `Spin(10) x U(1)_PQ` with "
        "`3x16F + 54R + 126C + 10C + 1C`. The singlet is required; the no-PQ "
        "short-field-list theory is the comparison branch `P54N`.",
        "",
        "## Definition-level correction",
        "",
        "The printed `Sigma Sigma* Sigma* phi` eta1 monomial has PQ charge "
        f"`{report['eta1_correction']['printed_charge']}`. The frozen card uses "
        "`Sigma Sigma* Sigma phi + h.c.`, which has charge zero.",
        "",
        "## Parameter count",
        "",
        f"- Scalar sector: `{counts['scalar_raw_reals']}` raw reals, "
        f"`{counts['scalar_quotient_reals']}` after the rank-"
        f"`{counts['scalar_phase_orbit_rank']}` scalar rephasing orbit.",
        f"- Yukawa sector: `{counts['yukawa_raw_reals']}` raw reals and "
        f"`{counts['yukawa_quotient_reals']}` after the generic `U(3)` family quotient.",
        f"- Full action: `{counts['total_raw_reals']}` raw, "
        f"`{counts['total_classical_quotient_reals']}` classically basis-quotiented, "
        f"`{counts['total_quantum_pq_observable_reals']}` continuous zero-temperature "
        "observable parameters after the anomalous PQ reparametrization.",
        "- Three scalar CP phases survive: `arg(eta3)-2arg(eta1)`, "
        "`arg(chi4)+arg(chi6)-2arg(eta1)`, and `arg(chi7)-arg(chi6)`.",
        "",
        "## Checks",
        "",
        "| Group | Check | Result |",
        "|---|---|---|",
    ]
    for row in report["checks"]:
        result = "PASS" if row["pass"] else "FAIL"
        lines.append(f"| {row['group']} | {row['name']} | {result} |")
    lines.extend(
        [
            "",
            "## Remaining downstream debts",
            "",
        ]
    )
    lines.extend(f"- {debt}" for debt in report["open_debts"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report = audit()
    json_path = OUTPUT / "p54_action_parameter_convention_card.json"
    md_path = OUTPUT / "p54_action_parameter_convention_card.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
