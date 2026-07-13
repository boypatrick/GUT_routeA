#!/usr/bin/env python3
"""Component-level triplet Clebsch nullspace audit.

No web lookup is used.  The previous two-kernel audit showed that deriving
colored-triplet tensors from the same CP1 Veronese/contact flavor operators
nearly fixes CKM and masses but does not provide enough Knu suppression for the
1e35 yr future-stress target.

This script asks the next sharper question.  Keep the successful doublet
Yukawa matrices fixed.  In the colored-triplet sector, expand

    Y_QQ^T = sum_I p_I B_Q^I,
    Y_QL^T = sum_J q_J B_L^J,

where the component basis is the current operator basis

    B_Q = (H, r_u F, a_u G_A, b_u G_B),
    B_L = r_d (H, F, a_d G_A, b_d G_B).

For the Knu operators with external flavors (u,u,s,nu_l), build the bilinear
map

    C_a = sum_{I,J} p_I M_a^{IJ} q_J,

for all three quark permutations and all three neutrino flavors.  Then audit:

  * fixed-QQ linear nullspace in q;
  * fixed-QL linear nullspace in p;
  * rank-one bilinear null, p_I q_J, with both norms held fixed.

This is still not the final Spin(10) contraction calculation.  It is the
component-level nullspace target that a real 10/bar126/120 triplet mass matrix
or symmetry must reproduce.
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
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402
import fit_two_kernel_flavor_then_d5 as two_kernel  # noqa: E402


INPUT = ROOT / "output" / "two_kernel_flavor_then_d5" / "summary.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "triplet_clebsch_nullspace"

LABEL = "d5_both_F_minus_0"
RNG_SEED = 202605082154
NULL_TOL = 1.0e-8


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cval(raw: dict[str, float]) -> complex:
    return raw["re"] + 1j * raw["im"]


def cvec(raw: list[dict[str, float]]) -> np.ndarray:
    return np.array([cval(z) for z in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


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


def load_component_card() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    payload = read_json(INPUT)
    card = payload["matrix_cards"][LABEL]
    basis = two_kernel.load_veronese_basis()
    k = two_kernel.analytic_k()
    md = card["model_data"]
    h_geo = fit.symmetric_from_coeff(cvec(md["H_coefficients"]), basis)
    f_geo = fit.symmetric_from_coeff(cvec(md["F_coefficients"]), basis)
    h = h_geo + cval(md["epsilon_H"]) * k
    fmat = f_geo + cval(md["epsilon_F"]) * k
    ga = fit.antisym_from_coeff(cvec(md["G_A_coefficients"]))
    gb = fit.antisym_from_coeff(cvec(md["G_B_coefficients"]))
    mix = {name: cval(raw) for name, raw in md["mixing_coefficients"].items()}
    y = {name: fit.cmat(raw) for name, raw in card["Yukawa_fit"].items()}
    data = {"H": h, "F": fmat, "G_A": ga, "G_B": gb, "mix": mix, "source_label": LABEL}
    return y, data


def rotations_for(y: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    ul: dict[str, np.ndarray] = {}
    for name, mat in y.items():
        left, _right, _s = two_kernel.left_right_rotation(mat)
        ul[name] = left
    return ul, two_kernel.pmns_target(ul["charged_lepton"])


def component_bases(data: dict[str, Any]) -> tuple[list[np.ndarray], list[np.ndarray], list[str]]:
    h = data["H"]
    fmat = data["F"]
    ga = data["G_A"]
    gb = data["G_B"]
    mix = data["mix"]
    names = ["H", "F", "G_A", "G_B"]
    bq = [h, mix["r_u"] * fmat, mix["a_u"] * ga, mix["b_u"] * gb]
    bl = [
        mix["r_d"] * h,
        mix["r_d"] * fmat,
        mix["r_d"] * mix["a_d"] * ga,
        mix["r_d"] * mix["b_d"] * gb,
    ]
    return bq, bl, names


def entry_indices() -> list[tuple[int, int, int, int]]:
    return [(a, b, c, ell) for (a, b, c) in sorted(set(permutations((0, 0, 1)))) for ell in range(3)]


def build_entry_matrices(
    bq: list[np.ndarray],
    bl: list[np.ndarray],
    ul: dict[str, np.ndarray],
    u_nu: np.ndarray,
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    indices = entry_indices()
    mats = np.zeros((len(indices), len(bq), len(bl)), dtype=complex)
    for i, qq in enumerate(bq):
        for j, ql in enumerate(bl):
            tensor = d5.tensor_llll(
                d5.sym(qq),
                ql,
                ul["up"],
                ul["up"],
                ul["down"],
                u_nu,
            )
            for n, idx in enumerate(indices):
                mats[n, i, j] = tensor[idx]
    return mats, indices


def physical_amplitude(
    p: np.ndarray,
    q: np.ndarray,
    bq: list[np.ndarray],
    bl: list[np.ndarray],
    ul: dict[str, np.ndarray],
    u_nu: np.ndarray,
) -> tuple[float, tuple[int, int, int, int]]:
    qq_raw = sum(p[i] * bq[i] for i in range(len(bq)))
    ql_raw = sum(q[j] * bl[j] for j in range(len(bl)))
    qq = d5.scale_to_largest_singular(qq_raw, 0.60)
    ql = d5.scale_to_largest_singular(ql_raw, 0.024)
    tensor = d5.tensor_llll(d5.sym(qq), ql, ul["up"], ul["up"], ul["down"], u_nu)
    return d5.max_perm_entry(tensor, (0, 0, 1), None)


def svd_nullity(mat: np.ndarray, tol: float = 1.0e-10) -> tuple[np.ndarray, np.ndarray, int]:
    _u, s, vh = np.linalg.svd(mat, full_matrices=True)
    rank = int(np.sum(s > tol * max(float(s[0]), 1.0)))
    nullity = mat.shape[1] - rank
    return s, vh, nullity


def row_for(
    label: str,
    p: np.ndarray,
    q: np.ndarray,
    p0: np.ndarray,
    q0: np.ndarray,
    bq: list[np.ndarray],
    bl: list[np.ndarray],
    ul: dict[str, np.ndarray],
    u_nu: np.ndarray,
    baseline_amp: float,
    future_margin0: float,
    residual: float | None = None,
) -> dict[str, Any]:
    amp, idx = physical_amplitude(p, q, bq, bl, ul, u_nu)
    ratio = amp / max(baseline_amp, 1.0e-30)
    return {
        "label": label,
        "amplitude": float(amp),
        "selected_index": "".join(str(i) for i in idx),
        "amplitude_ratio_to_natural": float(ratio),
        "future_margin_1e35": float(future_margin0 / max(ratio * ratio, 1.0e-30)),
        "p_deformation_ratio": float(np.linalg.norm(p - p0) / max(np.linalg.norm(p0), 1.0e-30)),
        "q_deformation_ratio": float(np.linalg.norm(q - q0) / max(np.linalg.norm(q0), 1.0e-30)),
        "p_norm": float(np.linalg.norm(p)),
        "q_norm": float(np.linalg.norm(q)),
        "null_residual": None if residual is None else float(residual),
        "p": [cjson(z) for z in p],
        "q": [cjson(z) for z in q],
    }


def bilinear_entries(mats: np.ndarray, p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.einsum("i,nij,j->n", p, mats, q, optimize=True)


def solve_rank_one_null(mats: np.ndarray, p0: np.ndarray, q0: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    rng = np.random.default_rng(RNG_SEED)
    pnorm = float(np.linalg.norm(p0))
    qnorm = float(np.linalg.norm(q0))

    def objective(x: np.ndarray) -> np.ndarray:
        p = unpack_complex(x[:8])
        q = unpack_complex(x[8:])
        vals = bilinear_entries(mats, p, q)
        pieces = list(np.real(vals)) + list(np.imag(vals))
        pieces.append(float(np.linalg.norm(p) / pnorm - 1.0))
        pieces.append(float(np.linalg.norm(q) / qnorm - 1.0))
        return np.array(pieces, dtype=float)

    starts = [np.concatenate([pack_complex(p0), pack_complex(q0)])]
    # Seed with fixed-linear smallest-singular directions.
    lq = np.einsum("i,nij->nj", p0, mats, optimize=True)
    lp = np.einsum("nij,j->ni", mats, q0, optimize=True)
    _s, vhq, _null = svd_nullity(lq)
    _s2, vhp, _null2 = svd_nullity(lp)
    starts.append(np.concatenate([pack_complex(p0), pack_complex(normalize_to(vhq[-1].conjugate(), qnorm))]))
    starts.append(np.concatenate([pack_complex(normalize_to(vhp[-1].conjugate(), pnorm)), pack_complex(q0)]))
    for _ in range(21):
        p = normalize_to(rng.normal(size=4) + 1j * rng.normal(size=4), pnorm)
        q = normalize_to(rng.normal(size=4) + 1j * rng.normal(size=4), qnorm)
        starts.append(np.concatenate([pack_complex(p), pack_complex(q)]))

    best_x = starts[0]
    best_cost = math.inf
    for start in starts:
        res = least_squares(
            objective,
            start,
            xtol=2.0e-12,
            ftol=2.0e-12,
            gtol=2.0e-12,
            max_nfev=1600,
        )
        cost = float(np.linalg.norm(objective(res.x)))
        if cost < best_cost:
            best_cost = cost
            best_x = res.x
    p_best = unpack_complex(best_x[:8])
    q_best = unpack_complex(best_x[8:])
    null_residual = float(np.linalg.norm(bilinear_entries(mats, p_best, q_best)))
    return p_best, q_best, null_residual, best_cost


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    y, data = load_component_card()
    ul, u_nu = rotations_for(y)
    bq, bl, names = component_bases(data)
    mats, indices = build_entry_matrices(bq, bl, ul, u_nu)
    future_margin0 = float(read_json(KNU_TARGET)["verdict"]["future_margin"])

    p0 = np.ones(4, dtype=complex)
    q0 = np.ones(4, dtype=complex)
    baseline_amp, baseline_idx = physical_amplitude(p0, q0, bq, bl, ul, u_nu)

    lq = np.einsum("i,nij->nj", p0, mats, optimize=True)
    lp = np.einsum("nij,j->ni", mats, q0, optimize=True)
    sq, vhq, q_nullity = svd_nullity(lq)
    sp, vhp, p_nullity = svd_nullity(lp)
    q_smin = normalize_to(vhq[-1].conjugate(), float(np.linalg.norm(q0)))
    p_smin = normalize_to(vhp[-1].conjugate(), float(np.linalg.norm(p0)))

    cmat = mats.reshape(len(indices), 16)
    sc, vhc, z_nullity = svd_nullity(cmat)
    p_rank1, q_rank1, rank1_residual, rank1_cost = solve_rank_one_null(mats, p0, q0)

    rows = [
        row_for("natural_doublet_like", p0, q0, p0, q0, bq, bl, ul, u_nu, baseline_amp, future_margin0),
        row_for("fixed_QQ_smallest_QL_singular", p0, q_smin, p0, q0, bq, bl, ul, u_nu, baseline_amp, future_margin0),
        row_for("fixed_QL_smallest_QQ_singular", p_smin, q0, p0, q0, bq, bl, ul, u_nu, baseline_amp, future_margin0),
        row_for(
            "rank_one_bilinear_null_search",
            p_rank1,
            q_rank1,
            p0,
            q0,
            bq,
            bl,
            ul,
            u_nu,
            baseline_amp,
            future_margin0,
            rank1_residual,
        ),
    ]

    with (OUT / "component_nullspace_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [key for key in rows[0].keys() if key not in {"p", "q"}]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})

    best = max(rows, key=lambda row: row["future_margin_1e35"])
    rank_one_passes = bool(rank1_residual < NULL_TOL and rows[-1]["future_margin_1e35"] >= 1.0)
    finite_deformation = bool(rows[-1]["p_deformation_ratio"] < 2.0 and rows[-1]["q_deformation_ratio"] < 2.0)
    antisym_qq_only = bool(
        abs(p_rank1[0]) < 1.0e-10
        and abs(p_rank1[1]) < 1.0e-10
        and abs(p_rank1[3]) < 1.0e-10
        and abs(p_rank1[2]) > 1.0
        and np.linalg.norm(q_rank1 - q0) < 1.0e-10
    )

    payload: dict[str, Any] = {
        "note": "No web lookup used. Component-level triplet Clebsch nullspace audit.",
        "source": {
            "two_kernel_summary": str(INPUT),
            "source_label": LABEL,
            "component_basis_names": names,
        },
        "math": {
            "operator_basis_QQ": "B_Q=(H,r_u F,a_u G_A,b_u G_B)",
            "operator_basis_QL": "B_L=r_d(H,F,a_d G_A,b_d G_B)",
            "external_entries": ["".join(str(i) for i in idx) for idx in indices],
            "bilinear_map": "C_a=sum_IJ p_I M_a^{IJ} q_J",
        },
        "linear_algebra": {
            "fixed_QQ_Lq_singular_values": [float(v) for v in sq],
            "fixed_QQ_Lq_nullity": int(q_nullity),
            "fixed_QL_Lp_singular_values": [float(v) for v in sp],
            "fixed_QL_Lp_nullity": int(p_nullity),
            "full_linear_map_C_rank": int(np.linalg.matrix_rank(cmat, tol=1.0e-10)),
            "full_linear_map_C_nullity_in_p_tensor_q": int(z_nullity),
            "full_linear_map_C_singular_values": [float(v) for v in sc],
            "rank_one_search_residual": rank1_residual,
            "rank_one_search_cost": rank1_cost,
        },
        "rows": rows,
        "verdict": {
            "fixed_QQ_all_Knu_null_exists": bool(q_nullity > 0),
            "fixed_QL_all_Knu_null_exists": bool(p_nullity > 0),
            "rank_one_bilinear_null_found": bool(rank1_residual < NULL_TOL),
            "rank_one_future_safe": bool(rows[-1]["future_margin_1e35"] >= 1.0),
            "rank_one_finite_deformation": finite_deformation,
            "rank_one_null_is_antisymmetric_QQ_only": antisym_qq_only,
            "component_level_candidate_exists": bool(rank_one_passes and finite_deformation),
            "best_row": {key: best[key] for key in best if key not in {"p", "q"}},
            "interpretation": (
                "Fixed-QQ and fixed-QL linear nullities test whether one side of "
                "the triplet coupling can be varied alone.  The rank-one search "
                "tests whether the full p_I q_J component map has a nontrivial "
                "all-Knu null compatible with finite triplet coefficient norms. "
                "Here the null is the antisymmetric QQ-only direction: the LLLL "
                "QQ contraction sees sym(Y_QQ), so a pure 120-like antisymmetric "
                "QQ triplet source vanishes.  This is a component-level target for "
                "a Higgs triplet mass-matrix derivation, not yet a full d=5 proof; "
                "RRRR channels and complete 120 component matching remain to be "
                "audited."
            ),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Triplet Clebsch component nullspace audit",
        "",
        "No web lookup was used.",
        "",
        "For fixed doublet flavor data, define",
        "",
        "```text",
        "Y_QQ^T = sum_I p_I B_Q^I,   B_Q=(H,r_u F,a_u G_A,b_u G_B)",
        "Y_QL^T = sum_J q_J B_L^J,   B_L=r_d(H,F,a_d G_A,b_d G_B)",
        "C_a = sum_IJ p_I M_a^{IJ} q_J",
        "```",
        "",
        "The entries `a` are the three `(u,u,s)` permutations times three neutrino flavors.",
        "",
        "## Linear algebra",
        "",
        f"- fixed QQ nullity in QL coefficients: `{q_nullity}`",
        f"- fixed QL nullity in QQ coefficients: `{p_nullity}`",
        f"- linear nullity in unrestricted `p tensor q`: `{z_nullity}`",
        f"- rank-one null residual: `{rank1_residual:.6e}`",
        f"- rank-one null is antisymmetric QQ only: `{antisym_qq_only}`",
        "",
        "| row | amp ratio | future margin | p deformation | q deformation | residual |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        residual = row["null_residual"]
        report.append(
            f"| `{row['label']}` | {row['amplitude_ratio_to_natural']:.6e} | "
            f"{row['future_margin_1e35']:.6e} | {row['p_deformation_ratio']:.6e} | "
            f"{row['q_deformation_ratio']:.6e} | "
            f"{0.0 if residual is None else residual:.6e} |"
        )
    report += [
        "",
        "## Verdict",
        "",
        payload["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Triplet Clebsch component nullspace audit")
    print(f"  fixed QQ nullity: {q_nullity}")
    print(f"  fixed QL nullity: {p_nullity}")
    print(f"  tensor nullity: {z_nullity}")
    print(f"  rank-one residual: {rank1_residual:.6e}")
    print(f"  rank-one future margin: {rows[-1]['future_margin_1e35']:.6e}")
    print(f"  candidate exists: {payload['verdict']['component_level_candidate_exists']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
