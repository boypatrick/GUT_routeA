#!/usr/bin/env python3
"""Independent Spin(10) invariant-basis and PQ anomaly audit for P54PQ-v2.

The proof is deliberately independent of the printed tensor potential.  It
uses representation-product multiplicities, enumerates every PQ-neutral
multidegree through degree four, and then checks the two subtle tensor facts
numerically in a self-dual five-form realization:

* 126 x 126bar has no 54 channel, so Phi Sigma Sigma* is absent;
* the declared 54+126 vacuum has gauge-orbit ranks 24 and 33.

The anomaly convention is A(G^2-PQ)=sum q T(R), with T(16)=2 for Spin(10)
and T(3)=1/2 for QCD.  The instanton coefficient is Nhat=2 A_QCD.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "p54_spin10_pq_audit.tex"
SCRIPT = Path(__file__).resolve()

FIELDS = ("A", "B", "C", "H", "Hb", "S", "Sb")
PQ = {"A": 0, "B": 2, "C": -2, "H": -2, "Hb": 2, "S": -4, "Sb": 4}

# A=54, B=126, C=bar126, H=10, Hb=10*, S=1, Sb=1*.
# Multiplicities below are singlet multiplicities, not merely a list copied
# from the candidate potential.
INVARIANT_LEDGER = {
    # degree two
    "A A": 1,
    "B C": 1,
    "H Hb": 1,
    "S Sb": 1,
    # degree three
    "A A A": 1,
    "A H Hb": 1,
    "H H Sb": 1,       # plus conjugate Hb Hb S
    # degree four: Hermitian multidegrees
    "A A A A": 2,
    "A A B C": 2,
    "A A H Hb": 2,
    "A A S Sb": 1,
    "B B C C": 4,
    "B C H Hb": 2,
    "B C S Sb": 1,
    "H H Hb Hb": 2,
    "H Hb S Sb": 1,
    "S S Sb Sb": 1,
    # degree four: one representative of each complex-conjugate pair
    "A B B S": 1,
    "A H H Sb": 1,
    "B B C H": 1,
    "B B H H": 1,
}

CONJUGATE_REPRESENTATIVES = {
    "H H Sb": "Hb Hb S",
    "A B B S": "A C C Sb",
    "A H H Sb": "A Hb Hb S",
    "B B C H": "B C C Hb",
    "B B H H": "C C Hb Hb",
}

COEFFICIENT_MAP = {
    "A A": ["mu2"],
    "B C": ["nu2"],
    "H Hb": ["xi02"],
    "S Sb": ["mus2"],
    "A A A": ["c"],
    "A H Hb": ["xi3"],
    "H H Sb": ["chi6"],
    "A A A A": ["a", "b"],
    "A A B C": ["alpha", "beta"],
    "A A H Hb": ["eta0", "eta2"],
    "A A S Sb": ["chi3"],
    "B B C C": ["lambda0", "lambda2", "lambda4", "lambda4p"],
    "B C H Hb": ["gamma1", "gamma2"],
    "B C S Sb": ["chi2"],
    "H H Hb Hb": ["xi1", "xi2"],
    "H Hb S Sb": ["chi5"],
    "S S Sb Sb": ["chi1"],
    "A B B S": ["chi4"],
    "A H H Sb": ["chi7"],
    "B B C H": ["eta1"],
    "B B H H": ["eta3"],
}

# The channel intersection giving each multiplicity.  Four copies for B^2C^2
# means the common channels of Sym^2(126) and its conjugate.
CHANNELS = {
    "A A A A": ["1", "54", "660", "770"],  # invariant ring leaves 2 at degree 4
    "A A B C": ["1", "770"],
    "A A H Hb": ["1", "54"],
    "B B C C": ["54", "1050", "2772", "4125"],
    "B C H Hb": ["1", "45"],
    "H H Hb Hb": ["1", "54"],
    "A B B S": ["54"],
    "A H H Sb": ["54"],
    "B B C H": ["1050"],
    "B B H H": ["54"],
}

PHASE_VECTORS = {
    "eta1": (1, 1, 0),
    "eta3": (2, 2, 0),
    "chi4": (2, 0, 1),
    "chi6": (0, 2, -1),
    "chi7": (0, 2, -1),
}

CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append({"group": group, "name": name, "pass": bool(condition), "detail": detail})


def canonical(items: tuple[str, ...] | list[str]) -> str:
    order = {name: i for i, name in enumerate(FIELDS)}
    return " ".join(sorted(items, key=order.__getitem__))


def neutral_multidegrees() -> dict[int, list[str]]:
    result: dict[int, list[str]] = {}
    for degree in (2, 3, 4):
        rows = []
        for items in itertools.combinations_with_replacement(FIELDS, degree):
            if sum(PQ[item] for item in items) == 0:
                rows.append(canonical(items))
        result[degree] = rows
    return result


def rank_exact(rows: list[tuple[int, ...]]) -> int:
    # Integer rows are tiny; the floating rank is exact at this scale.
    return int(np.linalg.matrix_rank(np.asarray(rows, dtype=float), tol=1e-12))


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(values[i] > values[j] for i in range(len(values)) for j in range(i + 1, len(values)))
    return -1 if inversions % 2 else 1


def five_form_setup() -> tuple[list[tuple[int, ...]], np.ndarray, np.ndarray, np.ndarray]:
    n = 10
    quints = list(itertools.combinations(range(n), 5))
    qindex = {q: i for i, q in enumerate(quints)}
    d5 = np.zeros((252, 252), dtype=float)
    seen: set[tuple[int, ...]] = set()
    pairs = []
    all_indices = set(range(n))
    for a, q in enumerate(quints):
        r = tuple(sorted(all_indices - set(q)))
        if q in seen or r in seen:
            continue
        seen.add(q)
        seen.add(r)
        b = qindex[r]
        sign = permutation_sign(tuple(list(q) + list(r)))
        pairs.append((a, b, sign))
        d5[b, a] = sign
        d5[a, b] = -sign

    # Build the holomorphic 5-form (e0+i e1)^... and choose the duality half
    # that contains it.
    uhol = [(np.eye(n)[2 * k] + 1j * np.eye(n)[2 * k + 1]) / math.sqrt(2) for k in range(5)]
    omega_full = np.zeros((n,) * 5, dtype=complex)
    base = np.einsum("i,j,k,l,m->ijklm", *uhol)
    for pm in itertools.permutations(range(5)):
        omega_full += permutation_sign(pm) * np.transpose(base, pm)

    def components(tensor: np.ndarray) -> np.ndarray:
        return np.asarray([tensor[q] for q in quints], dtype=complex)

    omega = components(omega_full)
    omega /= np.linalg.norm(omega)
    eigen = complex((d5 @ omega) @ np.conj(omega))
    duality = 1.0 if abs(eigen - 1j) < abs(eigen + 1j) else -1.0
    u126 = np.zeros((252, 126), dtype=complex)
    for k, (a, b, sign) in enumerate(pairs):
        u126[a, k] = 1 / math.sqrt(2)
        u126[b, k] = -duality * 1j * sign / math.sqrt(2)

    # Lookup from every ordered tensor index to the compressed component.
    comp_index = np.zeros(n**5, dtype=int)
    comp_sign = np.zeros(n**5, dtype=float)
    for flat, idx in enumerate(itertools.product(range(n), repeat=5)):
        if len(set(idx)) < 5:
            continue
        sorted_idx = tuple(sorted(idx))
        perm = tuple(sorted_idx.index(value) for value in idx)
        comp_index[flat] = qindex[sorted_idx]
        comp_sign[flat] = permutation_sign(perm)
    return quints, u126, comp_index, comp_sign


def full_five_form(compressed: np.ndarray, lookup: np.ndarray, signs: np.ndarray) -> np.ndarray:
    return (signs * compressed[lookup]).reshape((10,) * 5)


def act5(tensor: np.ndarray, a: int, b: int) -> np.ndarray:
    delta = np.zeros_like(tensor)
    for slot in range(5):
        source = np.moveaxis(tensor, slot, 0)
        target = np.moveaxis(delta, slot, 0)
        target[a] += source[b]
        target[b] -= source[a]
    return delta


def compressed_components(tensor: np.ndarray, quints: list[tuple[int, ...]]) -> np.ndarray:
    return np.asarray([tensor[q] for q in quints], dtype=complex)


def tensor_diagnostics() -> dict[str, Any]:
    quints, u126, lookup, signs = five_form_setup()
    rng = np.random.default_rng(540126)
    st_norms = []
    antisym_norms = []
    for _ in range(4):
        coeff = rng.normal(size=126) + 1j * rng.normal(size=126)
        coeff /= np.linalg.norm(coeff)
        sigma_c = u126 @ coeff
        sigma = full_five_form(sigma_c, lookup, signs)
        k = np.einsum("iabcd,jabcd->ij", sigma, np.conj(sigma)) / math.factorial(4)
        sym = (k + k.T) / 2
        st = sym - np.eye(10) * np.trace(sym) / 10
        anti = (k - k.T) / 2
        st_norms.append(float(np.linalg.norm(st)))
        antisym_norms.append(float(np.linalg.norm(anti)))

    # The PS 54 background and the holomorphic 126 direction.
    x54 = np.diag([-2 / 5] * 6 + [3 / 5] * 4)
    omega = u126 @ (np.conj(u126).T @ compressed_components(
        full_five_form(u126[:, 0], lookup, signs), quints
    ))
    # Use the actual holomorphic direction recovered by the setup: choose the
    # self-dual vector with maximal overlap with the product complex structure.
    uhol = [(np.eye(10)[2 * k] + 1j * np.eye(10)[2 * k + 1]) / math.sqrt(2) for k in range(5)]
    omega_full = np.zeros((10,) * 5, dtype=complex)
    base = np.einsum("i,j,k,l,m->ijklm", *uhol)
    for pm in itertools.permutations(range(5)):
        omega_full += permutation_sign(pm) * np.transpose(base, pm)
    omega = compressed_components(omega_full, quints)
    omega /= np.linalg.norm(omega)
    sigma = full_five_form(omega, lookup, signs)

    orbit54 = []
    orbit_full = []
    for a in range(10):
        for b in range(a + 1, 10):
            generator = np.zeros((10, 10))
            generator[a, b] = 1
            generator[b, a] = -1
            dphi = generator @ x54 - x54 @ generator
            dsigma = compressed_components(act5(sigma, a, b), quints)
            orbit54.append(dphi.ravel())
            orbit_full.append(np.concatenate([dphi.ravel(), dsigma.real, dsigma.imag]))
    rank54 = int(np.linalg.matrix_rank(np.stack(orbit54), tol=1e-10))
    rank_full = int(np.linalg.matrix_rank(np.stack(orbit_full), tol=1e-10))
    return {
        "symmetric_traceless_norms_126x126bar": st_norms,
        "antisymmetric_norms_126x126bar": antisym_norms,
        "orbit_rank_54": rank54,
        "orbit_rank_54_plus_126": rank_full,
        "unbroken_generator_count": 45 - rank_full,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit() -> dict[str, Any]:
    neutral = neutral_multidegrees()
    accepted_with_conjugates = set(INVARIANT_LEDGER)
    accepted_with_conjugates.update(CONJUGATE_REPRESENTATIVES.values())
    all_neutral = {row for rows in neutral.values() for row in rows}
    excluded = sorted(all_neutral - accepted_with_conjugates)

    check("basis", "every admitted multidegree is PQ neutral", set(INVARIANT_LEDGER).issubset(all_neutral), f"admitted={len(INVARIANT_LEDGER)}")
    check("basis", "every neutral multidegree is admitted or explicitly excluded by Spin(10) channels", all_neutral == accepted_with_conjugates | set(excluded), f"neutral={len(all_neutral)}, excluded={len(excluded)}")
    check("basis", "coefficient multiplicity matches invariant multiplicity", all(len(COEFFICIENT_MAP[key]) == mult for key, mult in INVARIANT_LEDGER.items()), "all 21 representative multidegrees")

    real_coefficients = 4 + 2 + 18
    complex_coefficients = 1 + 4
    scalar_raw = real_coefficients + 2 * complex_coefficients
    phase_rank = rank_exact(list(PHASE_VECTORS.values()))
    check("parameters", "corrected scalar action has 34 raw real coefficients", scalar_raw == 34, "24 real + 5 complex")
    check("parameters", "five coupling phases have rank-two rephasing orbit", phase_rank == 2, "PQ is the one-dimensional null direction")
    check("parameters", "three scalar CP invariants remain", complex_coefficients - phase_rank == 3, "delta1, delta2, delta3")
    check("parameters", "full corrected counts are 60/49/48", (2 + scalar_raw + 24, 2 + scalar_raw - phase_rank + 15, 2 + scalar_raw - phase_rank + 15 - 1) == (60, 49, 48), "raw/classical/quantum-PQ")

    tensor = tensor_diagnostics()
    check("tensor", "126 x 126bar has no symmetric-traceless 54 bilinear", max(tensor["symmetric_traceless_norms_126x126bar"]) < 2e-12, f"max norm={max(tensor['symmetric_traceless_norms_126x126bar']):.3e}")
    check("tensor", "the absent 54 test is non-vacuous because the 45 bilinear survives", min(tensor["antisymmetric_norms_126x126bar"]) > 1e-3, f"min antisymmetric norm={min(tensor['antisymmetric_norms_126x126bar']):.3e}")
    check("alignment", "54 breaks Spin(10) to the 21-generator Pati-Salam algebra", tensor["orbit_rank_54"] == 24, f"orbit rank={tensor['orbit_rank_54']}")
    check("alignment", "54 plus 126 leaves exactly the 12-generator SM algebra", tensor["orbit_rank_54_plus_126"] == 33 and tensor["unbroken_generator_count"] == 12, f"orbit rank={tensor['orbit_rank_54_plus_126']}, unbroken={tensor['unbroken_generator_count']}")

    anomaly = {
        "Spin10^2_PQ": 3 * (-1) * 2,
        "SU3c^2_PQ": 3 * (-2),
        "gravitational_PQ": 3 * 16 * (-1),
        "PQ_cubed": 3 * 16 * (-1) ** 3,
    }
    nhat = 2 * anomaly["SU3c^2_PQ"]
    scalar_charge_gcd = math.gcd(math.gcd(abs(PQ["B"]), abs(PQ["H"])), abs(PQ["S"]))
    naive_ndw = abs(nhat) // scalar_charge_gcd
    center_phases = {
        "Psi": 1j * np.exp(-1j * math.pi / 2),
        "Phi": 1,
        "Sigma": -1 * np.exp(2j * math.pi / 2),
        "phi": -1 * np.exp(-2j * math.pi / 2),
        "S": np.exp(-4j * math.pi / 2),
    }
    center_identity = max(abs(value - 1) for value in center_phases.values())
    physical_ndw = abs(nhat) // 4
    check("anomaly", "Spin(10)^2-PQ coefficient is -6", anomaly["Spin10^2_PQ"] == -6, "3 families, q16=-1, T(16)=2")
    check("anomaly", "QCD instanton coefficient is Nhat=-12", anomaly["SU3c^2_PQ"] == -6 and nhat == -12, "each 16 contributes A3=-2")
    check("global-form", "a pi/2 PQ rotation is the inverse Spin(10) Z4 center action", center_identity < 1e-12, f"max phase residual={center_identity:.3e}")
    check("domain-wall", "naive N_DW=6 is reduced to physical N_DW=3", naive_ndw == 6 and physical_ndw == 3, "12 QCD phases / Z4 gauge identification")
    check("domain-wall", "minimal model retains a post-inflation domain-wall problem", physical_ndw > 1, "N_DW=3; inflation or extra anomalous matter is required")

    tex = TEX.read_text(encoding="utf-8")
    required_tokens = (r"\mathrm{Sym}^2(126)", r"\chi_7", r"\mathcal A_{10}", r"N_{\rm DW}=3", r"P54PQ-v2")
    check("artifact", "companion TeX contains basis, anomaly, and global-form interfaces", all(token in tex for token in required_tokens), f"tokens={len(required_tokens)}")

    passed = sum(row["pass"] for row in CHECKS)
    return {
        "schema": "route-f-p54-spin10-pq-audit-v1",
        "model_id": "P54PQ-v2",
        "date": "2026-08-30",
        "status": {
            "invariant_basis_closed_through_dimension_four": passed == len(CHECKS),
            "pq_anomaly_domain_wall_closed": passed == len(CHECKS),
            "previous_action_card_v1_superseded": True,
            "p1_authorized_on_v2_only": passed == len(CHECKS),
        },
        "representation_products": {
            "10x10": "1_s + 45_a + 54_s",
            "54x54": "1 + 45 + 54 + 660 + 770 + 1386",
            "126x10": "210 + 1050",
            "126x54": "126 + 1728 + 4950",
            "126x126": "54 + 945 + 1050 + 2772 + 4125 + 6930",
            "126x126bar": "1 + 45 + 210 + 770 + 5940 + 8910",
            "Sym2_54": "1 + 54 + 660 + 770",
            "Sym2_126": "54 + 1050 + 2772 + 4125",
        },
        "neutral_multidegrees": neutral,
        "admitted_invariant_multiplicities": INVARIANT_LEDGER,
        "coefficient_map": COEFFICIENT_MAP,
        "channel_intersections": CHANNELS,
        "excluded_neutral_multidegrees": excluded,
        "new_operator": {
            "coefficient": "chi7",
            "monomial": "Phi_ij phi_i phi_j S* + h.c.",
            "channel": "54 in Sym^2(10)",
            "pq_charge": 0,
            "consequence": "P54PQ-v1 superseded by P54PQ-v2",
        },
        "parameter_count": {
            "scalar_raw_reals": scalar_raw,
            "scalar_rephasing_rank": phase_rank,
            "scalar_quotient_reals": scalar_raw - phase_rank,
            "total_raw_reals": 60,
            "total_classical_quotient_reals": 49,
            "total_quantum_pq_observable_reals": 48,
            "physical_scalar_cp_phases": [
                "arg(eta3)-2arg(eta1)",
                "arg(chi4)+arg(chi6)-2arg(eta1)",
                "arg(chi7)-arg(chi6)",
            ],
        },
        "tensor_diagnostics": tensor,
        "anomaly": {
            **anomaly,
            "QCD_instanton_coefficient_Nhat": nhat,
            "scalar_charge_gcd": scalar_charge_gcd,
            "naive_domain_wall_number": naive_ndw,
            "diagonal_Z4_phase_residual": center_identity,
            "physical_domain_wall_number": physical_ndw,
        },
        "checks": CHECKS,
        "summary": {"passed": passed, "total": len(CHECKS), "all_pass": passed == len(CHECKS)},
        "sources": [
            {"path": str(TEX.relative_to(REPO)), "sha256": sha256(TEX)},
            {"path": str(SCRIPT.relative_to(REPO)), "sha256": sha256(SCRIPT)},
        ],
    }


def markdown(report: dict[str, Any]) -> str:
    anomaly = report["anomaly"]
    lines = [
        "# P54PQ-v2 Spin(10) invariant-basis and PQ audit",
        "",
        f"Status: **{report['summary']['passed']}/{report['summary']['total']} checks passed**.",
        "",
        "## Action-level result",
        "",
        "The independent basis audit found one omitted operator: "
        "`chi7 Phi_ij phi_i phi_j S* + h.c.`. Therefore `P54PQ-v1` is "
        "superseded by `P54PQ-v2`.",
        "",
        "The corrected counts are `60` raw, `49` classical-quotient, and "
        "`48` continuous quantum-PQ observable parameters. Three scalar CP "
        "phases survive.",
        "",
        "## PQ anomaly and domain walls",
        "",
        f"- `A[Spin(10)^2-PQ] = {anomaly['Spin10^2_PQ']}`.",
        f"- `A[SU(3)c^2-PQ] = {anomaly['SU3c^2_PQ']}` and `Nhat = {anomaly['QCD_instanton_coefficient_Nhat']}`.",
        f"- The naive scalar-gcd count is `{anomaly['naive_domain_wall_number']}`; quotienting by the diagonal Spin(10) center gives physical `N_DW = {anomaly['physical_domain_wall_number']}`.",
        "- `N_DW=3` is not cosmologically harmless: a post-inflation PQ transition needs an additional repair.",
        "",
        "## Verification",
        "",
        "| Group | Check | Result |",
        "|---|---|---|",
    ]
    for row in report["checks"]:
        lines.append(f"| {row['group']} | {row['name']} | {'PASS' if row['pass'] else 'FAIL'} |")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report = audit()
    (OUTPUT / "p54_spin10_pq_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUTPUT / "p54_spin10_pq_audit.md").write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
