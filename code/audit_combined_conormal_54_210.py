#!/usr/bin/env python3
"""Combined conormal action for aligned constrained 54/210 sources.

This script tests whether the constrained/composite source sector can be
written as local F-term constraints without reintroducing relative-orientation
moduli.  The proposed conormal equations are

    C54(S) = S^2 - S - 6 I = 0,
    A(S,Omega) = (K_S - 12 I) Omega = 0,
    N(Omega) = <Omega,Omega> - 1 = 0,

where K_S is the induced action of the symmetric tensor S on four-forms.  At
S0=diag(-2^6,3^4), the eigenvalue 12 subspace of K_S on Lambda^4 is the single
weak-volume line e7^e8^e9^e10.  Thus A aligns Omega with the weak four-plane
selected by S.

The rank test is the central check:

    variables: delta S in 54, delta Omega in Lambda^4 R^10 (210 total)
    constraints: C54 (55 rows), A (210 rows), N (1 row)

If the Jacobian has nullity 24, the only zero modes are the shared
SO(10)/(SO(6)xSO(4)) orbit coordinates.  Normal-bundle multipliers then pair
all normal directions; unprojected multipliers leave extra massless multiplier
zero modes and are rejected as propagating fields.
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
OUT = ROOT / "output" / "combined_conormal_54_210"


FourIndex = Tuple[int, int, int, int]


def mat_inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a.T @ b))


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


def so10_generators(n: int = 10) -> List[np.ndarray]:
    gens: List[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            mat = np.zeros((n, n), dtype=float)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            gens.append(mat)
    assert len(gens) == 45
    return gens


def four_basis(n: int = 10) -> List[FourIndex]:
    return list(combinations(range(n), 4))


def permutation_sign(seq: Iterable[int]) -> int:
    vals = list(seq)
    inv = 0
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            inv += int(vals[i] > vals[j])
    return -1 if inv % 2 else 1


def omega_component(index: Tuple[int, int, int, int], support: FourIndex) -> float:
    if len(set(index)) < 4:
        return 0.0
    if set(index) != set(support):
        return 0.0
    return float(permutation_sign(index))


def action_on_fourform(mat: np.ndarray, omega: np.ndarray, basis4: List[FourIndex]) -> np.ndarray:
    """Induced gl(10) action on independent four-form coordinates."""

    pos = {idx: i for i, idx in enumerate(basis4)}

    def comp(idx: Tuple[int, int, int, int]) -> float:
        if len(set(idx)) < 4:
            return 0.0
        ordered = tuple(sorted(idx))
        return float(permutation_sign(idx) * omega[pos[ordered]])

    out = np.zeros(len(basis4), dtype=float)
    for out_pos, out_idx in enumerate(basis4):
        total = 0.0
        out_list = list(out_idx)
        for slot in range(4):
            i_slot = out_list[slot]
            for p in range(10):
                coeff = mat[i_slot, p]
                if coeff == 0.0:
                    continue
                in_idx = tuple(out_list[:slot] + [p] + out_list[slot + 1 :])
                total += coeff * comp(in_idx)
        out[out_pos] = total
    return out


def k_matrix_on_fourforms(mat: np.ndarray, basis4: List[FourIndex]) -> np.ndarray:
    columns = []
    for j in range(len(basis4)):
        omega = np.zeros(len(basis4), dtype=float)
        omega[j] = 1.0
        columns.append(action_on_fourform(mat, omega, basis4))
    return np.column_stack(columns)


def linearized_c54(s0: np.ndarray, x: np.ndarray) -> np.ndarray:
    return s0 @ x + x @ s0 - x


def build_jacobian() -> Dict[str, object]:
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    support = (6, 7, 8, 9)
    s_domain = symmetric_traceless_basis()
    c_codomain = symmetric_basis()
    basis4 = four_basis()
    omega0 = np.zeros(len(basis4), dtype=float)
    omega0[basis4.index(support)] = 1.0
    k0 = k_matrix_on_fourforms(s0, basis4)
    a0 = k0 - 12.0 * np.eye(len(basis4))

    rows = 55 + 210 + 1
    cols = 54 + 210
    jac = np.zeros((rows, cols), dtype=float)

    # C54 rows.
    for a, e in enumerate(s_domain):
        image = linearized_c54(s0, e)
        for b, c in enumerate(c_codomain):
            jac[b, a] = mat_inner(c, image)

    # Alignment rows A(S,Omega)=(K_S-12)Omega.
    start = 55
    for a, e in enumerate(s_domain):
        jac[start : start + 210, a] = action_on_fourform(e, omega0, basis4)
    jac[start : start + 210, 54:] = a0

    # Normalization row.
    jac[-1, 54:] = 2.0 * omega0

    singular = np.linalg.svd(jac, compute_uv=False)
    rank = int(np.sum(singular > 1.0e-10))
    nullity = cols - rank

    # Independent S and Omega orbits would give 48 orientation modes.  The
    # aligned conormal constraints should leave only the diagonal/shared 24.
    s_tangents = []
    o_tangents = []
    for gen in so10_generators():
        s_tangents.append(gen @ s0 - s0 @ gen)
        o_tangents.append(action_on_fourform(gen, omega0, basis4))
    tangent_matrix = np.zeros((cols, len(s_tangents)))
    for j, (st, ot) in enumerate(zip(s_tangents, o_tangents)):
        for i, b in enumerate(s_domain):
            tangent_matrix[i, j] = mat_inner(b, st)
        tangent_matrix[54:, j] = ot
    tangent_rank = int(np.sum(np.linalg.svd(tangent_matrix, compute_uv=False) > 1.0e-10))

    full_multiplier_dim = rows
    normal_multiplier_dim = rank
    full_hessian_zeros = cols + full_multiplier_dim - 2 * rank
    normal_hessian_zeros = cols + normal_multiplier_dim - 2 * rank

    return {
        "variables_dimension": cols,
        "constraint_rows_unprojected": rows,
        "jacobian_rank": rank,
        "jacobian_nullity": nullity,
        "shared_orbit_tangent_rank": tangent_rank,
        "independent_orbit_rank_before_alignment": 48,
        "relative_orientation_modes_removed": 24,
        "full_multiplier_dim": full_multiplier_dim,
        "normal_bundle_multiplier_dim": normal_multiplier_dim,
        "full_multiplier_hessian_zero_modes": full_hessian_zeros,
        "normal_bundle_hessian_zero_modes": normal_hessian_zeros,
        "radial_extension_zero_modes": normal_hessian_zeros + 2,
        "nonzero_singular_min": float(min(x for x in singular if x > 1.0e-10)),
        "nonzero_singular_max": float(max(singular)),
        "passes": bool(rank == 240 and nullity == 24 and tangent_rank == 24 and normal_hessian_zeros == 24),
    }


def d210_counts() -> Dict[str, int]:
    basis2 = list(combinations(range(10), 2))
    support = (6, 7, 8, 9)
    mat = np.zeros((len(basis2), len(basis2)), dtype=float)
    for a, (i, j) in enumerate(basis2):
        for b, (k, l) in enumerate(basis2):
            mat[a, b] = omega_component((i, j, k, l), support)
    eig = np.linalg.eigvalsh(0.5 * (mat + mat.T))
    return {
        "-1": int(np.sum(np.isclose(eig, -1.0, atol=1.0e-10))),
        "0": int(np.sum(np.isclose(eig, 0.0, atol=1.0e-10))),
        "+1": int(np.sum(np.isclose(eig, 1.0, atol=1.0e-10))),
    }


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_report(payload: Dict[str, object]) -> str:
    j = payload["combined_jacobian"]
    lines = [
        "# Combined conormal 54/210 audit",
        "",
        "No web lookup was used.",
        "",
        "## F-term system",
        "",
        "```text",
        "C54(S) = S^2 - S - 6 I = 0",
        "A(S,Omega) = (K_S - 12 I) Omega = 0",
        "N(Omega) = <Omega,Omega> - 1 = 0",
        "```",
        "",
        "## Rank audit",
        "",
        f"- variables dimension: `{j['variables_dimension']}`",
        f"- unprojected constraint rows: `{j['constraint_rows_unprojected']}`",
        f"- Jacobian rank/nullity: `{j['jacobian_rank']}/{j['jacobian_nullity']}`",
        f"- shared orbit tangent rank: `{j['shared_orbit_tangent_rank']}`",
        f"- relative orientation modes removed: `{j['relative_orientation_modes_removed']}`",
        "",
        "## Multiplier Hessian",
        "",
        f"- full multiplier Hessian zero modes: `{j['full_multiplier_hessian_zero_modes']}`",
        f"- normal-bundle multiplier Hessian zero modes: `{j['normal_bundle_hessian_zero_modes']}`",
        f"- with two radial singlets retained: `{j['radial_extension_zero_modes']}`",
        "",
        "## 210 Clebsch check",
        "",
        f"`D210` eigenvalue counts on two-forms: `{payload['D210_eigen_counts']}`.",
        "",
        "## Verdict",
        "",
        str(payload["verdict"]),
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jac = build_jacobian()
    counts = d210_counts()
    payload: Dict[str, object] = {
        "note": "No web lookup used. Combined conormal action for aligned constrained 54/210 sources.",
        "combined_jacobian": jac,
        "D210_eigen_counts": counts,
        "passes": bool(jac["passes"] and counts == {"-1": 3, "0": 39, "+1": 3}),
        "verdict": (
            "The aligned conormal F-term system has nullity 24, exactly the "
            "shared SO(10)/(SO(6)xSO(4)) orbit.  A normal-bundle multiplier "
            "pairs all 240 normal directions and leaves only those orbit modes; "
            "unprojected multipliers would leave 50 Hessian zero modes, including "
            "26 extra multiplier/radial modes.  With two radial source singlets "
            "retained deliberately, the physical zero-mode count is 26.  Therefore "
            "the constrained branch is viable only if the multipliers are "
            "normal-bundle auxiliary/composite fields, not propagating unprojected "
            "54/210 towers."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(make_report(payload), encoding="utf-8")
    write_csv(OUT / "rank_summary.csv", [jac])
    print("Combined conormal 54/210 audit complete")
    print(f"passes={payload['passes']}")
    print(f"rank={jac['jacobian_rank']}")
    print(f"nullity={jac['jacobian_nullity']}")
    print(f"normal_hessian_zeros={jac['normal_bundle_hessian_zero_modes']}")
    print(f"radial_extension_zero_modes={jac['radial_extension_zero_modes']}")
    print(f"D210_counts={counts}")
    print(f"outputs={OUT}")


if __name__ == "__main__":
    main()
