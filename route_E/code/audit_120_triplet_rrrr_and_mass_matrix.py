#!/usr/bin/env python3
"""Audit the 120-like antisymmetric triplet null against RRRR and mass lifts.

No web lookup is used.  The previous component audit found an exact monitored
LLLL Knu null when the QQT coupling is a pure antisymmetric 120-like direction,

    p = (0,0,2,0),  q = (1,1,1,1).

That is not yet a proton-decay proof.  In a Spin(10)/Pati-Salam triplet sector
the same 10-10 triplet source also feeds a right-handed UE coupling, and the
same 10-5bar source feeds UD.  This script therefore checks:

  * locked-component RRRR leakage for the antisymmetric QQT/UE direction;
  * whether q can be moved to an RRRR null while keeping the LLLL null;
  * whether a joint rank-one p_I q_J null exists for LLLL+RRRR;
  * whether a finite 4x4 triplet mass-matrix lift of the rank-one source
    reintroduces leakage through orthogonal complete-source directions.

The calculation is a component-level obstruction/target audit, not a full
SO(10) tensor contraction.  It gives the next action-level constraint that an
explicit 120_H triplet mass matrix must satisfy.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_triplet_clebsch_nullspace_from_components as comp  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402
import fit_two_kernel_flavor_then_d5 as two_kernel  # noqa: E402


PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "triplet_120_rrrr_mass_matrix"

RNG_SEED = 202605082236
NULL_TOL = 1.0e-9


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def pack_complex(values: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(values), np.imag(values)]).astype(float)


def unpack_complex(values: np.ndarray) -> np.ndarray:
    n = len(values) // 2
    return values[:n] + 1j * values[n:]


def normalize_to(vec: np.ndarray, target_norm: float) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm <= 1.0e-30:
        return vec
    return vec * (target_norm / norm)


def proton_constants() -> dict[str, float]:
    return read_json(PROTON)["hadronic_constants"]


def target_filter() -> float:
    return float(read_json(KNU_TARGET)["S_T_target_inferred"])


def rrrr_indices() -> list[tuple[int, int, int, int]]:
    # Tensor order is (u_a, u_b, d_c, charged_d) after the helper convention.
    rows: list[tuple[int, int, int, int]] = []
    for ell in range(3):
        for u_pair in [(0, 0), (0, 1), (1, 0)]:
            rows.append((u_pair[0], u_pair[1], 1, ell))
    return rows


def build_locked_maps(
    b10: list[np.ndarray],
    b5: list[np.ndarray],
    ul: dict[str, np.ndarray],
    ur: dict[str, np.ndarray],
    u_nu: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, int, int]], list[tuple[int, int, int, int]]]:
    l_indices = comp.entry_indices()
    r_indices = rrrr_indices()
    lmap = np.zeros((len(l_indices), len(b10), len(b5)), dtype=complex)
    rmap = np.zeros((len(r_indices), len(b10), len(b5)), dtype=complex)
    for i, y10 in enumerate(b10):
        for j, y5 in enumerate(b5):
            # LLLL sees the symmetric QQ projection.
            ltensor = d5.tensor_llll(d5.sym(y10), y5, ul["up"], ul["up"], ul["down"], u_nu)
            # RRRR sees the full 10-10 source; antisymmetric 120-like pieces are not
            # projected away by the same symmetrizer in this proxy.
            rtensor = d5.tensor_rrrr(y10, y5, ur["up"], ur["charged_lepton"], ur["up"], ur["down"])
            for n, idx in enumerate(l_indices):
                lmap[n, i, j] = ltensor[idx]
            for n, idx in enumerate(r_indices):
                rmap[n, i, j] = rtensor[idx]
    return lmap, rmap, l_indices, r_indices


def bilinear_entries(mats: np.ndarray, p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.einsum("i,nij,j->n", p, mats, q, optimize=True)


def max_abs_entries(values: np.ndarray, indices: list[tuple[int, int, int, int]]) -> tuple[float, tuple[int, int, int, int], complex]:
    pick = int(np.argmax(np.abs(values)))
    return float(abs(values[pick])), indices[pick], complex(values[pick])


def physical_channel_amplitudes(
    p: np.ndarray,
    q: np.ndarray,
    b10: list[np.ndarray],
    b5: list[np.ndarray],
    ul: dict[str, np.ndarray],
    ur: dict[str, np.ndarray],
    u_nu: np.ndarray,
) -> dict[str, Any]:
    y10_raw = sum(p[i] * b10[i] for i in range(len(b10)))
    y5_raw = sum(q[j] * b5[j] for j in range(len(b5)))
    y10 = d5.scale_to_largest_singular(y10_raw, 0.60)
    y5 = d5.scale_to_largest_singular(y5_raw, 0.024)
    l_tensor = d5.tensor_llll(d5.sym(y10), y5, ul["up"], ul["up"], ul["down"], u_nu)
    r_tensor = d5.tensor_rrrr(y10, y5, ur["up"], ur["charged_lepton"], ur["up"], ur["down"])
    l_amp, l_idx = d5.max_perm_entry(l_tensor, (0, 0, 1), None)
    r_amp, r_idx = d5.max_rrrr_uusd(r_tensor, None)
    return {
        "LLLL_Knu": {"amplitude": float(l_amp), "index": "".join(str(x) for x in l_idx)},
        "RRRR_uusd": {"amplitude": float(r_amp), "index": "".join(str(x) for x in r_idx)},
    }


def add_lifetime(row: dict[str, Any], constants: dict[str, float], st: float) -> dict[str, Any]:
    out = dict(row)
    for key in ["LLLL_Knu", "RRRR_uusd"]:
        amp = out[key]["amplitude"]
        life = d5.lifetime_from_amplitude(amp, constants, triplet_filter=st)
        out[key]["tau_years_at_ST_target"] = life["tau_years"]
        out[key]["margin_1e35_at_ST_target"] = math.inf if math.isinf(life["tau_years"]) else life["tau_years"] / 1.0e35
    out["worst_margin_1e35_at_ST_target"] = min(out["LLLL_Knu"]["margin_1e35_at_ST_target"], out["RRRR_uusd"]["margin_1e35_at_ST_target"])
    return out


def svd_nullity(mat: np.ndarray, tol: float = NULL_TOL) -> tuple[np.ndarray, np.ndarray, int, int]:
    _u, s, vh = np.linalg.svd(mat, full_matrices=True)
    rank = int(np.sum(s > tol * max(float(s[0]), 1.0)))
    return s, vh, mat.shape[1] - rank, rank


def q_least_for_fixed_p(mats: np.ndarray, p: np.ndarray, q0: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, float]:
    lq = np.einsum("i,nij->nj", p, mats, optimize=True)
    s, vh, nullity, _rank = svd_nullity(lq)
    q = normalize_to(vh[-1].conjugate(), float(np.linalg.norm(q0)))
    residual = float(np.linalg.norm(lq @ q))
    return q, s, nullity, residual


def solve_joint_rank_one(
    lmap: np.ndarray,
    rmap: np.ndarray,
    p0: np.ndarray,
    q0: np.ndarray,
    p_seed: np.ndarray,
    q_seed: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    rng = np.random.default_rng(RNG_SEED)
    pnorm = float(np.linalg.norm(p0))
    qnorm = float(np.linalg.norm(q0))

    def objective(x: np.ndarray) -> np.ndarray:
        p = unpack_complex(x[:8])
        q = unpack_complex(x[8:])
        vals = np.concatenate([bilinear_entries(lmap, p, q), bilinear_entries(rmap, p, q)])
        pieces = list(np.real(vals)) + list(np.imag(vals))
        pieces.append(float(np.linalg.norm(p) / pnorm - 1.0))
        pieces.append(float(np.linalg.norm(q) / qnorm - 1.0))
        return np.array(pieces, dtype=float)

    starts = [
        np.concatenate([pack_complex(p0), pack_complex(q0)]),
        np.concatenate([pack_complex(p_seed), pack_complex(q_seed)]),
    ]
    for _ in range(24):
        p = normalize_to(rng.normal(size=4) + 1j * rng.normal(size=4), pnorm)
        q = normalize_to(rng.normal(size=4) + 1j * rng.normal(size=4), qnorm)
        starts.append(np.concatenate([pack_complex(p), pack_complex(q)]))

    best_x = starts[0]
    best_cost = math.inf
    for start in starts:
        res = least_squares(
            objective,
            start,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=1800,
        )
        cost = float(np.linalg.norm(objective(res.x)))
        if cost < best_cost:
            best_x = res.x
            best_cost = cost
    return unpack_complex(best_x[:8]), unpack_complex(best_x[8:]), best_cost


def complete_unitary(first: np.ndarray) -> np.ndarray:
    e = normalize_to(first.astype(complex), 1.0)
    # null_space returns an orthonormal basis for vectors orthogonal to e^dagger.
    rest = null_space(e.conjugate().reshape(1, -1))
    return np.column_stack([e, rest])


def linear_w_entries(mats: np.ndarray, w: np.ndarray) -> np.ndarray:
    return np.einsum("ij,nij->n", w, mats, optimize=True)


def finite_lift_rows(
    p: np.ndarray,
    q: np.ndarray,
    p0: np.ndarray,
    q0: np.ndarray,
    lmap: np.ndarray,
    rmap: np.ndarray,
    l_indices: list[tuple[int, int, int, int]],
    r_indices: list[tuple[int, int, int, int]],
) -> list[dict[str, Any]]:
    u = complete_unitary(p)
    v = complete_unitary(q)
    natural = np.outer(p0, q0.conjugate())
    l_nat = float(np.max(np.abs(linear_w_entries(lmap, natural))))
    r_nat = float(np.max(np.abs(linear_w_entries(rmap, natural))))
    rows: list[dict[str, Any]] = []
    for kappa in [30.0, 100.0, 300.0, 1000.0, 1.0e4]:
        eps = 1.0 / kappa
        w = u @ np.diag([1.0, eps, eps, eps]) @ v.conjugate().T
        l_vals = linear_w_entries(lmap, w)
        r_vals = linear_w_entries(rmap, w)
        l_amp, l_idx, _l_val = max_abs_entries(l_vals, l_indices)
        r_amp, r_idx, _r_val = max_abs_entries(r_vals, r_indices)
        rows.append(
            {
                "label": f"finite_mass_lift_kappa_{int(kappa)}",
                "kappa": float(kappa),
                "epsilon": float(eps),
                "linear_LLLL_ratio_to_natural": float(l_amp / max(l_nat, 1.0e-30)),
                "linear_RRRR_ratio_to_natural": float(r_amp / max(r_nat, 1.0e-30)),
                "linear_LLLL_amplitude": float(l_amp),
                "linear_RRRR_amplitude": float(r_amp),
                "LLLL_index": "".join(str(x) for x in l_idx),
                "RRRR_index": "".join(str(x) for x in r_idx),
                "mass_condition_number": float(kappa),
                "interpretation": "Finite lift of rank-one W by equal orthogonal singular values.",
            }
        )
    return rows


def compact_row(
    label: str,
    p: np.ndarray,
    q: np.ndarray,
    p0: np.ndarray,
    q0: np.ndarray,
    b10: list[np.ndarray],
    b5: list[np.ndarray],
    ul: dict[str, np.ndarray],
    ur: dict[str, np.ndarray],
    u_nu: np.ndarray,
    constants: dict[str, float],
    st: float,
    residual: float | None = None,
) -> dict[str, Any]:
    row = {
        "label": label,
        **physical_channel_amplitudes(p, q, b10, b5, ul, ur, u_nu),
        "p_deformation_ratio": float(np.linalg.norm(p - p0) / max(np.linalg.norm(p0), 1.0e-30)),
        "q_deformation_ratio": float(np.linalg.norm(q - q0) / max(np.linalg.norm(q0), 1.0e-30)),
        "residual": None if residual is None else float(residual),
        "p": [cjson(z) for z in p],
        "q": [cjson(z) for z in q],
    }
    return add_lifetime(row, constants, st)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    y, data = comp.load_component_card()
    b10, b5, names = comp.component_bases(data)
    ul, ur, u_nu = two_kernel.rotations_for(y)
    constants = proton_constants()
    st = target_filter()

    lmap, rmap, l_indices, r_indices = build_locked_maps(b10, b5, ul, ur, u_nu)
    joint_w = np.vstack([lmap.reshape((len(l_indices), -1)), rmap.reshape((len(r_indices), -1))])
    joint_s, _joint_vh, joint_nullity, joint_rank = svd_nullity(joint_w)

    p0 = np.ones(4, dtype=complex)
    q0 = np.ones(4, dtype=complex)
    p_anti = np.array([0.0, 0.0, 2.0, 0.0], dtype=complex)
    q_rrrr, rrrr_q_s, rrrr_q_nullity, rrrr_q_residual = q_least_for_fixed_p(rmap, p_anti, q0)
    p_joint, q_joint, joint_rank_one_cost = solve_joint_rank_one(lmap, rmap, p0, q0, p_anti, q_rrrr)
    joint_residual = float(
        np.linalg.norm(np.concatenate([bilinear_entries(lmap, p_joint, q_joint), bilinear_entries(rmap, p_joint, q_joint)]))
    )

    rows = [
        compact_row("natural_locked", p0, q0, p0, q0, b10, b5, ul, ur, u_nu, constants, st),
        compact_row("antisymmetric_120_QQT_locked_natural_QLT", p_anti, q0, p0, q0, b10, b5, ul, ur, u_nu, constants, st),
        compact_row(
            "antisymmetric_120_QQT_with_best_RRRR_q",
            p_anti,
            q_rrrr,
            p0,
            q0,
            b10,
            b5,
            ul,
            ur,
            u_nu,
            constants,
            st,
            residual=rrrr_q_residual,
        ),
        compact_row(
            "joint_rank_one_LLLL_RRRR_search",
            p_joint,
            q_joint,
            p0,
            q0,
            b10,
            b5,
            ul,
            ur,
            u_nu,
            constants,
            st,
            residual=joint_residual,
        ),
    ]
    finite_rows = finite_lift_rows(p_anti, q_rrrr, p0, q0, lmap, rmap, l_indices, r_indices)

    with (OUT / "triplet_120_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "label",
            "LLLL_amp",
            "LLLL_margin_1e35",
            "RRRR_amp",
            "RRRR_margin_1e35",
            "worst_margin_1e35",
            "p_deformation_ratio",
            "q_deformation_ratio",
            "residual",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "label": row["label"],
                    "LLLL_amp": row["LLLL_Knu"]["amplitude"],
                    "LLLL_margin_1e35": row["LLLL_Knu"]["margin_1e35_at_ST_target"],
                    "RRRR_amp": row["RRRR_uusd"]["amplitude"],
                    "RRRR_margin_1e35": row["RRRR_uusd"]["margin_1e35_at_ST_target"],
                    "worst_margin_1e35": row["worst_margin_1e35_at_ST_target"],
                    "p_deformation_ratio": row["p_deformation_ratio"],
                    "q_deformation_ratio": row["q_deformation_ratio"],
                    "residual": row["residual"],
                }
            )

    exact_joint_rank_one = joint_residual < 1.0e-8
    fixed_anti_rrrr_null = rrrr_q_residual < 1.0e-8
    verdict_text = (
        "The antisymmetric 120-like QQT direction exactly removes the monitored LLLL Knu tensor. "
        "Under locked Spin(10)-component bookkeeping the same 10-10 source feeds RRRR, but the "
        "fixed-antisymmetric-p RRRR map has a two-dimensional q-null.  The explicit joint rank-one "
        "solution is therefore p=(0,0,2,0) with q almost entirely in the fourth component, i.e. a "
        "crossed 120_A/120_B triplet projector: 10-10 uses the G_A antisymmetric source while "
        "10-5bar uses the orthogonal G_B-like source.  This is a genuine component-level "
        "LLLL+RRRR null.  It is still not an action-level proof, because a finite mass-matrix lift "
        "leaks through orthogonal directions at O(1/kappa); the next completion must derive this "
        "crossed 120 projector from the Higgs/triplet superpotential and feed the finite-lift leakage "
        "back into threshold and dressed proton scans."
    )
    output = {
        "note": "No web lookup used. 120-like triplet RRRR and finite mass-matrix lift audit.",
        "basis_names": names,
        "operator_definitions": {
            "locked_10_10_source": "same p_I B_10^I feeds QQT and UET",
            "locked_10_5bar_source": "same q_J B_5^J feeds QLT and UDT",
            "LLLL": "C_L=sum_IJ p_I q_J M_L^{IJ}, with sym(B_10^I) on QQ",
            "RRRR": "C_R=sum_IJ p_I q_J M_R^{IJ}, no sym projection on UE in this proxy",
            "finite_lift": "W_kappa=U diag(1,1/kappa,1/kappa,1/kappa) V^dagger",
        },
        "linear_algebra": {
            "joint_W_rank": joint_rank,
            "joint_W_nullity": joint_nullity,
            "joint_W_singular_values": [float(x) for x in joint_s],
            "fixed_antisymmetric_p_RRRR_q_nullity": rrrr_q_nullity,
            "fixed_antisymmetric_p_RRRR_q_singular_values": [float(x) for x in rrrr_q_s],
            "fixed_antisymmetric_p_RRRR_q_residual": float(rrrr_q_residual),
            "joint_rank_one_cost": float(joint_rank_one_cost),
            "joint_rank_one_residual": float(joint_residual),
        },
        "rows": rows,
        "finite_lift_rows": finite_rows,
        "verdict": {
            "fixed_anti_RRRR_exact_q_null": bool(fixed_anti_rrrr_null),
            "joint_rank_one_exact_LLLL_RRRR_null": bool(exact_joint_rank_one),
            "best_physical_row": min(rows, key=lambda row: row["worst_margin_1e35_at_ST_target"]),
            "most_safe_physical_row": max(rows, key=lambda row: row["worst_margin_1e35_at_ST_target"]),
            "interpretation": verdict_text,
        },
    }
    (OUT / "summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# 120-like triplet RRRR and mass-matrix audit",
        "",
        "No web lookup was used.",
        "",
        "Locked component assumption:",
        "",
        "```text",
        "Y_QQ^T and Y_UE^T use the same p_I B_10^I.",
        "Y_QL^T and Y_UD^T use the same q_J B_5^J.",
        "```",
        "",
        "## Physical p,q rows",
        "",
        "| row | LLLL amp | LLLL margin | RRRR amp | RRRR margin | worst margin | p def | q def | residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        residual = 0.0 if row["residual"] is None else row["residual"]
        report.append(
            f"| `{row['label']}` | {row['LLLL_Knu']['amplitude']:.3e} | "
            f"{row['LLLL_Knu']['margin_1e35_at_ST_target']:.3e} | "
            f"{row['RRRR_uusd']['amplitude']:.3e} | "
            f"{row['RRRR_uusd']['margin_1e35_at_ST_target']:.3e} | "
            f"{row['worst_margin_1e35_at_ST_target']:.3e} | "
            f"{row['p_deformation_ratio']:.3e} | {row['q_deformation_ratio']:.3e} | {residual:.3e} |"
        )
    report.extend(
        [
            "",
            "## Finite mass-matrix lift leakage",
            "",
            "| kappa | epsilon | LLLL ratio | RRRR ratio | LLLL index | RRRR index |",
            "|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in finite_rows:
        report.append(
            f"| {row['kappa']:.0f} | {row['epsilon']:.3e} | "
            f"{row['linear_LLLL_ratio_to_natural']:.3e} | "
            f"{row['linear_RRRR_ratio_to_natural']:.3e} | "
            f"`{row['LLLL_index']}` | `{row['RRRR_index']}` |"
        )
    report.extend(
        [
            "",
            "## Verdict",
            "",
            verdict_text,
            "",
            f"Fixed-antisymmetric-p RRRR q-nullity: `{rrrr_q_nullity}`.",
            f"Fixed-antisymmetric-p RRRR residual: `{rrrr_q_residual:.6e}`.",
            f"Joint rank-one LLLL+RRRR residual: `{joint_residual:.6e}`.",
            "",
        ]
    )
    (OUT / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print("120-like triplet RRRR/mass-matrix audit")
    print(f"  fixed anti-p RRRR q nullity: {rrrr_q_nullity}")
    print(f"  fixed anti-p RRRR residual: {rrrr_q_residual:.6e}")
    print(f"  joint rank-one residual: {joint_residual:.6e}")
    print(f"  exact joint rank-one null: {exact_joint_rank_one}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
