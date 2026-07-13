#!/usr/bin/env python3
"""Reduced constrained-source audit for the 54/210 projector sector.

The Schur-locking audit moved the main obstruction to the large 54/210
source/driver tower: complete multiplets are one-loop matching silent, but
they create a short SO(10) Landau-pole distance if treated as elementary
propagating chiral fields at M_G.

This script formalizes the reduced alternative:

* F54 is a constrained spectral order parameter on the SO(10)/(SO(6)xSO(4))
  orbit of S0=diag(-2^6,3^4).
* Omega210 is a constrained decomposable weak-volume four-form on the same
  Grassmannian orbit.
* Only radial singlets and gauge-orbit coordinates are retained; all normal
  non-singlet components are constrained/composite rather than elementary
  propagating fields.

The audit verifies orbit dimensions, Clebsch/projector data, and the SO(10)
UV beta-function cost of the reduced branch.
"""

from __future__ import annotations

from itertools import combinations
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "constrained_54_210_source_sector"

ALPHA_G_INV_REF = 41.36308203884868
BASELINE_SUM_T = 70.0
STERILE_10G_T = 1.0


def so10_generators(n: int = 10) -> List[np.ndarray]:
    gens: List[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            mat = np.zeros((n, n), dtype=float)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            gens.append(mat)
    return gens


def symmetric_traceless_basis(n: int = 10) -> List[np.ndarray]:
    basis: List[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            mat = np.zeros((n, n), dtype=float)
            mat[i, j] = 1.0 / math.sqrt(2.0)
            mat[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(mat)
    for k in range(n - 1):
        mat = np.zeros((n, n), dtype=float)
        coeff = 1.0 / math.sqrt((k + 1) * (k + 2))
        for i in range(k + 1):
            mat[i, i] = coeff
        mat[k + 1, k + 1] = -(k + 1) * coeff
        basis.append(mat)
    assert len(basis) == 54
    return basis


def symmetric_basis(n: int = 10) -> List[np.ndarray]:
    basis: List[np.ndarray] = []
    for i in range(n):
        mat = np.zeros((n, n), dtype=float)
        mat[i, i] = 1.0
        basis.append(mat)
    for i in range(n):
        for j in range(i + 1, n):
            mat = np.zeros((n, n), dtype=float)
            mat[i, j] = 1.0 / math.sqrt(2.0)
            mat[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(mat)
    assert len(basis) == 55
    return basis


def mat_inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a.T @ b))


def rank_from_svd(mat: np.ndarray, tol: float = 1.0e-10) -> int:
    return int(np.sum(np.linalg.svd(mat, compute_uv=False) > tol))


def audit_54() -> Dict[str, object]:
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    residual = s0 @ s0 - s0 - 6.0 * np.eye(10)

    domain = symmetric_traceless_basis()
    codomain = symmetric_basis()
    jac = np.zeros((len(codomain), len(domain)))
    for a, e in enumerate(domain):
        image = s0 @ e + e @ s0 - e
        for b, c in enumerate(codomain):
            jac[b, a] = mat_inner(c, image)
    fixed_rank = rank_from_svd(jac)
    fixed_nullity = 54 - fixed_rank

    tangents = []
    for t in so10_generators():
        tangents.append(t @ s0 - s0 @ t)
    tangent_matrix = np.array([[mat_inner(b, tan) for tan in tangents] for b in domain])
    orbit_rank = rank_from_svd(tangent_matrix)

    return {
        "representation_dimension": 54,
        "fixed_orbit_dimension": fixed_nullity,
        "orbit_rank_from_generators": orbit_rank,
        "radial_singlets_retained": 1,
        "reduced_local_dimension": orbit_rank + 1,
        "normal_modes_removed_if_radial_retained": 54 - (orbit_rank + 1),
        "spectral_constraint_rank_fixed_radius": fixed_rank,
        "spectral_constraint_residual_norm": float(np.linalg.norm(residual)),
        "stabilizer_dimension": 45 - orbit_rank,
        "expected_stabilizer": "SO(6)xSO(4), dimension 15+6=21",
        "passes": bool(orbit_rank == 24 and fixed_nullity == 24 and np.linalg.norm(residual) < 1.0e-12),
    }


FourIndex = Tuple[int, int, int, int]


def four_basis(n: int = 10) -> List[FourIndex]:
    return list(combinations(range(n), 4))


def permutation_sign(seq: Iterable[int]) -> int:
    values = list(seq)
    inv = 0
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            inv += int(values[i] > values[j])
    return -1 if inv % 2 else 1


def omega_component(index: Tuple[int, int, int, int], support: FourIndex) -> float:
    if len(set(index)) < 4:
        return 0.0
    if set(index) != set(support):
        return 0.0
    return float(permutation_sign(index))


def generator_action_on_fourform(gen: np.ndarray, support: FourIndex, basis: List[FourIndex]) -> np.ndarray:
    vec = np.zeros(len(basis), dtype=float)
    for out_pos, out_idx in enumerate(basis):
        total = 0.0
        out_list = list(out_idx)
        for slot in range(4):
            i_slot = out_list[slot]
            for p in range(10):
                coeff = gen[i_slot, p]
                if coeff == 0.0:
                    continue
                in_idx = tuple(out_list[:slot] + [p] + out_list[slot + 1 :])
                total += coeff * omega_component(in_idx, support)
        vec[out_pos] = total
    return vec


def twoform_basis(n: int = 10) -> List[Tuple[int, int]]:
    return list(combinations(range(n), 2))


def d210_matrix(support: FourIndex = (6, 7, 8, 9)) -> np.ndarray:
    basis = twoform_basis()
    pos = {pair: i for i, pair in enumerate(basis)}
    mat = np.zeros((len(basis), len(basis)), dtype=float)
    for a, pair_a in enumerate(basis):
        for b, pair_b in enumerate(basis):
            i, j = pair_a
            k, l = pair_b
            mat[a, b] = 0.5 * omega_component((i, j, k, l), support)
    # The two-form basis is ordered only by i<j.  The 1/2 in D=(1/2)Phi_ijkl A_kl
    # is compensated by summing k<l once; multiply by 2 to represent the operator
    # on independent two-form coordinates.
    return 2.0 * mat


def audit_210() -> Dict[str, object]:
    support = (6, 7, 8, 9)
    basis4 = four_basis()
    omega = np.zeros(len(basis4), dtype=float)
    omega[basis4.index(support)] = 1.0
    tangent_cols = [generator_action_on_fourform(gen, support, basis4) for gen in so10_generators()]
    tangent_matrix = np.column_stack(tangent_cols)
    orbit_rank = rank_from_svd(tangent_matrix)

    dmat = d210_matrix(support)
    eig = np.linalg.eigvalsh(0.5 * (dmat + dmat.T))
    counts = {
        "-1": int(np.sum(np.isclose(eig, -1.0, atol=1.0e-10))),
        "0": int(np.sum(np.isclose(eig, 0.0, atol=1.0e-10))),
        "+1": int(np.sum(np.isclose(eig, 1.0, atol=1.0e-10))),
    }
    return {
        "representation_dimension": len(basis4),
        "decomposable_support": [index + 1 for index in support],
        "omega_norm": float(np.linalg.norm(omega)),
        "orbit_rank_from_generators": orbit_rank,
        "radial_singlets_retained": 1,
        "reduced_local_dimension": orbit_rank + 1,
        "normal_modes_removed_if_radial_retained": len(basis4) - (orbit_rank + 1),
        "stabilizer_dimension": 45 - orbit_rank,
        "expected_stabilizer": "SO(6)xSO(4), dimension 15+6=21",
        "D210_eigenvalue_counts_on_twoforms": counts,
        "D210_eigenvalues_minmax": [float(np.min(eig)), float(np.max(eig))],
        "passes": bool(orbit_rank == 24 and counts == {"-1": 3, "0": 39, "+1": 3}),
    }


def projector_checks() -> Dict[str, object]:
    f_values = {
        "Sigma_L": 2.0,
        "Sigma_R": 2.0,
        "Sigma8": -4.0 / 3.0,
        "X622": 1.0 / 3.0,
    }
    d_values = {
        "Sigma_L": 1.0,
        "Sigma_R": -1.0,
        "Sigma8": 0.0,
        "X622": 0.0,
    }
    px = {
        key: -(9.0 / 25.0) * (f - 2.0) * (f + 4.0 / 3.0)
        for key, f in f_values.items()
    }
    pr = {
        key: 0.5 * (d * d - d)
        for key, d in d_values.items()
    }
    targets_x = {"Sigma_L": 0.0, "Sigma_R": 0.0, "Sigma8": 0.0, "X622": 1.0}
    targets_r = {"Sigma_L": 0.0, "Sigma_R": 1.0, "Sigma8": 0.0, "X622": 0.0}
    max_err = max(
        [abs(px[key] - targets_x[key]) for key in px]
        + [abs(pr[key] - targets_r[key]) for key in pr]
    )
    return {
        "F54_values": f_values,
        "D210_values": d_values,
        "PX_values": px,
        "PR_values": pr,
        "max_projector_error": max_err,
        "passes": bool(max_err < 1.0e-12),
    }


def landau_ratio(alpha_inv: float, sum_t: float) -> float:
    b10 = sum_t - 24.0
    if b10 <= 0:
        return math.inf
    return math.exp(2.0 * math.pi * alpha_inv / b10)


def uv_rows() -> List[Dict[str, object]]:
    scenarios = [
        ("baseline", BASELINE_SUM_T, "minimal matter/Yukawa/Higgs branch"),
        ("baseline_plus_10G", BASELINE_SUM_T + STERILE_10G_T, "minimal Schur completion"),
        ("constrained_54_210_sources_plus_10G", BASELINE_SUM_T + STERILE_10G_T, "54/210 constrained; only singlets/orbit coordinates retained"),
        ("elementary_54_210_plus_10G", BASELINE_SUM_T + STERILE_10G_T + 12.0 + 56.0, "one elementary 54 and one elementary 210 added"),
        ("documented_drive_tower_plus_10G", 185.0, "three 54 link/driver pairs plus one 210 pair plus 10G"),
    ]
    rows: List[Dict[str, object]] = []
    for name, sum_t, note in scenarios:
        b10 = sum_t - 24.0
        lp = landau_ratio(ALPHA_G_INV_REF, sum_t)
        rows.append(
            {
                "scenario": name,
                "sum_T": sum_t,
                "b10": b10,
                "landau_ratio_Lambda_over_MG": lp,
                "reaches_R50": bool(lp > 50.0),
                "reaches_R200": bool(lp > 200.0),
                "interpretation": note,
            }
        )
    return rows


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_report(payload: Dict[str, object]) -> str:
    a54 = payload["audit_54"]
    a210 = payload["audit_210"]
    proj = payload["projector_checks"]
    lines = [
        "# Constrained 54/210 source-sector audit",
        "",
        "No web lookup was used.",
        "",
        "## 54 spectral orbit",
        "",
        f"- fixed-orbit dimension: `{a54['fixed_orbit_dimension']}`",
        f"- orbit rank from SO(10) generators: `{a54['orbit_rank_from_generators']}`",
        f"- radial singlets retained: `{a54['radial_singlets_retained']}`",
        f"- normal modes removed if radial retained: `{a54['normal_modes_removed_if_radial_retained']}`",
        f"- spectral residual norm: `{a54['spectral_constraint_residual_norm']:.3e}`",
        "",
        "## 210 decomposable four-form orbit",
        "",
        f"- representation dimension: `{a210['representation_dimension']}`",
        f"- orbit rank from SO(10) generators: `{a210['orbit_rank_from_generators']}`",
        f"- radial singlets retained: `{a210['radial_singlets_retained']}`",
        f"- normal modes removed if radial retained: `{a210['normal_modes_removed_if_radial_retained']}`",
        f"- D210 two-form eigenvalue counts: `{a210['D210_eigenvalue_counts_on_twoforms']}`",
        "",
        "## Projector preservation",
        "",
        f"`max projector error = {proj['max_projector_error']:.3e}`.",
        "",
        "## UV rows",
        "",
        "| scenario | sum T | b10 | LP/MG | reaches R=50 | reaches R=200 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in payload["uv_rows"]:
        lines.append(
            f"| `{row['scenario']}` | {row['sum_T']:.0f} | {row['b10']:.0f} | "
            f"{row['landau_ratio_Lambda_over_MG']:.3g} | {row['reaches_R50']} | {row['reaches_R200']} |"
        )
    lines.extend(["", "## Verdict", "", str(payload["verdict"])])
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a54 = audit_54()
    a210 = audit_210()
    proj = projector_checks()
    uv = uv_rows()
    payload: Dict[str, object] = {
        "note": "No web lookup used. Reduced constrained 54/210 source-sector audit.",
        "audit_54": a54,
        "audit_210": a210,
        "projector_checks": proj,
        "uv_rows": uv,
        "verdict": (
            "The reduced constrained-source branch preserves the verified F54/D210 "
            "projector data while removing 29 normal 54 modes and 185 normal 210 "
            "modes from the elementary spectrum if radial singlets are retained. "
            "Because the remaining radial modes are SM singlets and the orbit modes "
            "are gauge orientations, this branch has the same SO(10) UV cost as the "
            "minimal Schur completion, not the large 54/210 driver tower.  It is a "
            "conditional constrained/composite source sector, not an elementary "
            "renormalizable 54_H+210_H tower."
        ),
        "passes": bool(a54["passes"] and a210["passes"] and proj["passes"]),
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(make_report(payload), encoding="utf-8")
    write_csv(OUT / "uv_rows.csv", uv)
    write_csv(
        OUT / "projector_values.csv",
        [
            {
                "fragment": key,
                "F54": proj["F54_values"][key],
                "D210": proj["D210_values"][key],
                "PX": proj["PX_values"][key],
                "PR": proj["PR_values"][key],
            }
            for key in proj["F54_values"]
        ],
    )
    print("Constrained 54/210 source-sector audit complete")
    print(f"passes={payload['passes']}")
    print(f"F54_normal_removed={a54['normal_modes_removed_if_radial_retained']}")
    print(f"D210_normal_removed={a210['normal_modes_removed_if_radial_retained']}")
    print(f"projector_error={proj['max_projector_error']:.3e}")
    print(f"outputs={OUT}")


if __name__ == "__main__":
    main()
