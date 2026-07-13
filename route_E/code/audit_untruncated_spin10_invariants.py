#!/usr/bin/env python3
"""Untruncated renormalizable Spin(10) invariant audit.

No web lookup is used.  This script classifies the renormalizable invariants
available to the 54_H, 210_H, and three adjoint 45 fields used by the local
mediator card.  It also performs explicit tensor sanity checks:

* allowed contractions are nonzero on random tensors;
* contractions claimed to vanish are numerically zero;
* aligned 54/210 Clebsch operators preserve the broad
  (15,1,1) + (6,2,2) + (1,3,1) + (1,1,3) block split;
* the distinct-adjoint structure-constant invariant fABC is allowed and is
  therefore a real coefficient that must be forbidden or scanned.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "untruncated_spin10_invariants"
CARD = ROOT / "output" / "spin10_vacuum_alignment" / "benchmark_R200.json"


def pair_basis(n: int = 10) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


PAIRS = pair_basis(10)
PAIR_INDEX = {pair: idx for idx, pair in enumerate(PAIRS)}
FOUR_BASIS = list(itertools.combinations(range(10), 4))
FOUR_INDEX = {basis: idx for idx, basis in enumerate(FOUR_BASIS)}


def permutation_sign(values: tuple[int, ...]) -> int:
    if len(set(values)) != len(values):
        return 0
    inv = 0
    for i, vi in enumerate(values):
        for vj in values[i + 1 :]:
            if vi > vj:
                inv += 1
    return -1 if inv % 2 else 1


def antisym_basis_matrix(pair: tuple[int, int], n: int = 10) -> np.ndarray:
    i, j = pair
    out = np.zeros((n, n), dtype=float)
    out[i, j] = 1.0
    out[j, i] = -1.0
    return out


BASIS_45 = [antisym_basis_matrix(pair) for pair in PAIRS]


def coeffs_from_antisym(a: np.ndarray) -> np.ndarray:
    return np.array([a[i, j] for i, j in PAIRS], dtype=float)


def random_antisym(rng: np.random.Generator, n: int = 10) -> np.ndarray:
    raw = rng.normal(size=(n, n))
    return raw - raw.T


def random_sym_traceless(rng: np.random.Generator, n: int = 10) -> np.ndarray:
    raw = rng.normal(size=(n, n))
    s = 0.5 * (raw + raw.T)
    return s - np.eye(n) * float(np.trace(s)) / n


def random_fourform(rng: np.random.Generator) -> np.ndarray:
    return rng.normal(size=len(FOUR_BASIS))


def four_component(vec: np.ndarray, inds: tuple[int, int, int, int]) -> float:
    if len(set(inds)) < 4:
        return 0.0
    ordered = tuple(sorted(inds))
    return float(permutation_sign(inds) * vec[FOUR_INDEX[ordered]])


def weak_volume_fourform() -> np.ndarray:
    vec = np.zeros(len(FOUR_BASIS), dtype=float)
    vec[FOUR_INDEX[(6, 7, 8, 9)]] = 1.0
    return vec


def canonical_color_twoform_6d() -> np.ndarray:
    a = np.zeros((6, 6), dtype=float)
    for i, j in [(0, 1), (2, 3), (4, 5)]:
        a[i, j] = 1.0
        a[j, i] = -1.0
    return a


def color_fourform_from_twoform(a6: np.ndarray) -> np.ndarray:
    vec = np.zeros(len(FOUR_BASIS), dtype=float)
    for basis, idx in FOUR_INDEX.items():
        if any(v >= 6 for v in basis):
            continue
        remaining = tuple(v for v in range(6) if v not in basis)
        if len(remaining) != 2:
            continue
        r, s = remaining
        vec[idx] = permutation_sign(basis + (r, s)) * a6[r, s]
    return vec


def d_matrix_from_fourform(phi: np.ndarray) -> np.ndarray:
    d = np.zeros((len(PAIRS), len(PAIRS)), dtype=float)
    for p, (i, j) in enumerate(PAIRS):
        for q, (k, l) in enumerate(PAIRS):
            d[p, q] = four_component(phi, (i, j, k, l))
    return d


def f54_operator() -> np.ndarray:
    # Normalize diag(-2^6, 3^4) by 1/3, giving
    # F54(color,color)=-4/3, F54(color,weak)=1/3, F54(weak,weak)=2.
    weights = [-2.0 / 3.0] * 6 + [1.0] * 4
    diag = [weights[i] + weights[j] for i, j in PAIRS]
    return np.diag(diag)


def color_adj_vev_10d() -> np.ndarray:
    a = np.zeros((10, 10), dtype=float)
    for i, j in [(0, 1), (2, 3), (4, 5)]:
        a[i, j] = 1.0
        a[j, i] = -1.0
    return a


def adjoint_action_operator(a0: np.ndarray) -> np.ndarray:
    op = np.zeros((len(PAIRS), len(PAIRS)), dtype=float)
    for col, basis in enumerate(BASIS_45):
        comm = a0 @ basis - basis @ a0
        op[:, col] = coeffs_from_antisym(comm)
    return op


def block_projectors() -> dict[str, np.ndarray]:
    p_color = np.zeros((45, 45), dtype=float)
    p_x = np.zeros((45, 45), dtype=float)
    p_weak = np.zeros((45, 45), dtype=float)
    for idx, (i, j) in enumerate(PAIRS):
        if i < 6 and j < 6:
            p_color[idx, idx] = 1.0
        elif i < 6 <= j:
            p_x[idx, idx] = 1.0
        else:
            p_weak[idx, idx] = 1.0
    d_weak = d_matrix_from_fourform(weak_volume_fourform())
    p_l = 0.5 * (p_weak + d_weak)
    p_r = 0.5 * (p_weak - d_weak)
    return {
        "P_color_15_1_1": p_color,
        "P_X_6_2_2": p_x,
        "P_L_1_3_1": p_l,
        "P_R_1_1_3": p_r,
    }


def off_block_norms(op: np.ndarray, projectors: dict[str, np.ndarray]) -> dict[str, float]:
    off = np.zeros_like(op)
    names = list(projectors)
    for a, name_a in enumerate(names):
        for b, name_b in enumerate(names):
            if a == b:
                continue
            off += projectors[name_a] @ op @ projectors[name_b]
    return {
        "max_abs": float(np.max(np.abs(off))),
        "frobenius": float(np.linalg.norm(off)),
    }


def trace_a3(a: np.ndarray) -> float:
    return float(np.trace(a @ a @ a))


def invariant_s3(s: np.ndarray) -> float:
    return float(np.trace(s @ s @ s))


def invariant_s2_phi(s1: np.ndarray, s2: np.ndarray, phi: np.ndarray) -> float:
    total = 0.0
    for i in range(10):
        for j in range(10):
            for k in range(10):
                for l in range(10):
                    total += s1[i, k] * s2[j, l] * four_component(phi, (i, j, k, l))
    return float(total)


def invariant_s_phi_phi(s: np.ndarray, phi: np.ndarray) -> float:
    q = np.zeros((10, 10), dtype=float)
    triples = list(itertools.combinations(range(10), 3))
    for i in range(10):
        for j in range(10):
            q[i, j] = sum(
                four_component(phi, (i,) + abc) * four_component(phi, (j,) + abc)
                for abc in triples
            )
    q -= np.eye(10) * float(np.trace(q)) / 10.0
    return float(np.sum(s * q))


def pform_action(t: np.ndarray, vec: np.ndarray) -> np.ndarray:
    out = np.zeros_like(vec)
    for basis, idx in FOUR_INDEX.items():
        value = vec[idx]
        if value == 0.0:
            continue
        raw = list(basis)
        for pos, old in enumerate(raw):
            for new in range(10):
                coeff = t[new, old]
                if coeff == 0.0:
                    continue
                replaced = raw.copy()
                replaced[pos] = new
                if len(set(replaced)) < 4:
                    continue
                sorted_tuple = tuple(sorted(replaced))
                out[FOUR_INDEX[sorted_tuple]] += (
                    coeff * permutation_sign(tuple(replaced)) * value
                )
    return out


def invariant_a_phi_phi(a: np.ndarray, phi: np.ndarray) -> float:
    return float(np.dot(phi, pform_action(a, phi)))


def invariant_a_s_s(a: np.ndarray, s: np.ndarray) -> float:
    return float(np.trace(s.T @ (a @ s - s @ a)))


def invariant_a_s_phi_metric_candidate(a: np.ndarray, s: np.ndarray, phi: np.ndarray) -> float:
    # The only metric-style attempt pairs the two-form and symmetric tensor
    # indices into the four-form.  It vanishes because the symmetric slots are
    # contracted against an antisymmetric four-form.
    total = 0.0
    for i in range(10):
        for j in range(10):
            for k in range(10):
                for l in range(10):
                    total += a[i, j] * s[k, l] * four_component(phi, (i, j, k, l))
    return float(total)


def invariant_phi3(phi: np.ndarray) -> float:
    d = d_matrix_from_fourform(phi)
    return float(np.trace(d @ d @ d))


def invariant_sxy(s: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    return float(np.trace(x @ s @ y))


def invariant_phixy(phi: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    d = d_matrix_from_fourform(phi)
    return float(coeffs_from_antisym(x) @ d @ coeffs_from_antisym(y))


def invariant_fxyz(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    return float(np.trace(x @ (y @ z - z @ y)))


def tensor_sanity_checks(seed: int = 20260506) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    a = random_antisym(rng)
    b = random_antisym(rng)
    c = random_antisym(rng)
    s = random_sym_traceless(rng)
    s2 = random_sym_traceless(rng)
    phi = random_fourform(rng)
    return {
        "seed": seed,
        "allowed_nonzero_examples": {
            "Tr_S3": invariant_s3(s),
            "S_Phi_Phi": invariant_s_phi_phi(s, phi),
            "Phi3_metric": invariant_phi3(phi),
            "S_X_Y": invariant_sxy(s, a, b),
            "Phi_X_Y": invariant_phixy(phi, a, b),
            "f_X_Y_Z": invariant_fxyz(a, b, c),
        },
        "vanishing_examples": {
            "Tr_A3_identical": trace_a3(a),
            "f_A_A_B": invariant_fxyz(a, a, b),
            "A_S_S_moment": invariant_a_s_s(a, s),
            "A_Phi_Phi_moment": invariant_a_phi_phi(a, phi),
            "S_S_Phi_metric": invariant_s2_phi(s, s2, phi),
            "A_S_Phi_metric_attempt": invariant_a_s_phi_metric_candidate(a, s, phi),
        },
    }


def operator_audit() -> dict[str, Any]:
    projectors = block_projectors()
    f54 = f54_operator()
    d_weak = d_matrix_from_fourform(weak_volume_fourform())
    d_color = d_matrix_from_fourform(color_fourform_from_twoform(canonical_color_twoform_6d()))
    ad_a0 = adjoint_action_operator(color_adj_vev_10d())
    ops = {
        "F54_from_aligned_54": f54,
        "D210_weak_volume": d_weak,
        "D210_color_fourform": d_color,
        "ad_A0_from_fABC": ad_a0,
    }
    rows: dict[str, Any] = {}
    for name, op in ops.items():
        svals = np.linalg.svd(op, compute_uv=False)
        rows[name] = {
            "broad_PS_offblock": off_block_norms(op, projectors),
            "spectral_norm": float(svals[0]),
            "frobenius_norm": float(np.linalg.norm(op)),
            "rank_numeric": int(np.sum(svals > 1.0e-10)),
            "max_abs_plane_basis_entry": float(np.max(np.abs(op))),
        }
    return rows


def invariant_table() -> list[dict[str, str]]:
    pairs = "i,j in {A,B,C}; symmetric pair count = 6"
    return [
        {
            "sector": "quadratic",
            "term": "54_H^2, 210_H^2",
            "contraction": "Tr S^2, Phi_{ijkl} Phi_{ijkl}",
            "status": "allowed",
            "off_block": "not a mediator off-block source",
            "action": "keep in full Higgs vacuum fit",
        },
        {
            "sector": "quadratic",
            "term": "45_i 45_j",
            "contraction": "<X_i,X_j>",
            "status": "allowed",
            "off_block": "block diagonal for aligned PS basis",
            "action": f"scan full 6-parameter mass matrix ({pairs}) or impose texture symmetry",
        },
        {
            "sector": "pure Higgs cubic",
            "term": "54_H^3",
            "contraction": "Tr S^3",
            "status": "allowed",
            "off_block": "controls 54 vacuum alignment",
            "action": "include in untruncated F-flatness equations",
        },
        {
            "sector": "pure Higgs cubic",
            "term": "54_H^2 210_H",
            "contraction": "S_{ik} S_{jl} Phi_{ijkl}",
            "status": "forbidden/zero",
            "off_block": "none",
            "action": "no scan coefficient; contraction vanishes by symmetric vs antisymmetric indices",
        },
        {
            "sector": "pure Higgs cubic",
            "term": "54_H 210_H^2",
            "contraction": "S_{ij} Phi_{iabc} Phi_{jabc}",
            "status": "allowed",
            "off_block": "can move Higgs alignment",
            "action": "include in full Higgs vacuum fit",
        },
        {
            "sector": "pure Higgs cubic",
            "term": "210_H^3",
            "contraction": "Tr_{Lambda^2} D_Phi^3",
            "status": "allowed",
            "off_block": "contains the verified PS color cubic branch",
            "action": "keep; match coefficient to Goldstone-locking sector",
        },
        {
            "sector": "one adjoint plus Higgs",
            "term": "45_i 54_H^2",
            "contraction": "<S,[X_i,S]>",
            "status": "zero for one 54_H",
            "off_block": "none",
            "action": "no scan coefficient unless multiple distinct 54 copies are added",
        },
        {
            "sector": "one adjoint plus Higgs",
            "term": "45_i 210_H^2",
            "contraction": "<Phi, X_i.Phi>",
            "status": "zero for one 210_H",
            "off_block": "none",
            "action": "no scan coefficient unless multiple distinct 210 copies are added",
        },
        {
            "sector": "one adjoint plus Higgs",
            "term": "45_i 54_H 210_H",
            "contraction": "no nonzero metric scalar",
            "status": "forbidden/zero",
            "off_block": "none",
            "action": "no scan coefficient",
        },
        {
            "sector": "two adjoints plus Higgs",
            "term": "54_H 45_i 45_j",
            "contraction": "Tr(X_i S X_j)",
            "status": "allowed",
            "off_block": "aligned S gives F54 and preserves broad PS blocks",
            "action": f"scan 6 coefficients ({pairs}); current card uses a tuned subset",
        },
        {
            "sector": "two adjoints plus Higgs",
            "term": "210_H 45_i 45_j",
            "contraction": "Phi_{ijkl} X_i^{ij} X_j^{kl}",
            "status": "allowed",
            "off_block": "aligned weak/color 210 preserves broad PS blocks but changes color/weak masses",
            "action": f"scan 6 coefficients ({pairs}); tie weak and color effects by one Spin(10) coupling matrix",
        },
        {
            "sector": "three adjoints",
            "term": "45_i 45_j 45_k",
            "contraction": "Tr X_i [X_j,X_k]",
            "status": "allowed only for three distinct adjoints",
            "off_block": "broad PS blocks preserved by A0, but component mixing is O(1)",
            "action": "for A,B,C this is a real rho_ABC coefficient; forbid by symmetry or include in Hessian scan",
        },
        {
            "sector": "three adjoints",
            "term": "45_i^3 or 45_i^2 45_j",
            "contraction": "Tr X^3 or Tr X[X,Y]",
            "status": "zero",
            "off_block": "none",
            "action": "no scan coefficient for repeated adjoints",
        },
    ]


def load_reference_threshold() -> dict[str, float]:
    if not CARD.exists():
        return {}
    card = json.loads(CARD.read_text(encoding="utf-8"))
    best = card["rge_proton_replay"]["best"]
    return {
        "R": float(card["R"]),
        "finite_mediator_projected_l2": float(best["mediator_projected_l2"]),
        "M_GeV": float(best["MG_GeV"]),
    }


def write_csv(table: list[dict[str, str]]) -> None:
    with (OUT / "untruncated_spin10_invariant_table.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sector", "term", "contraction", "status", "off_block", "action"],
        )
        writer.writeheader()
        writer.writerows(table)


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Untruncated Spin(10) invariant audit")
    lines.append("")
    lines.append("No web lookup was used.  This audit lists the renormalizable")
    lines.append("Spin(10)-invariant contractions involving `54_H`, `210_H`, and the")
    lines.append("three adjoint mediator fields `A,B,C` of the local projector card.")
    lines.append("")
    lines.append("## Invariant census")
    lines.append("")
    lines.append("| sector | term | status | off-block verdict | required action |")
    lines.append("|---|---|---|---|---|")
    for row in payload["invariant_table"]:
        lines.append(
            f"| {row['sector']} | `{row['term']}` | {row['status']} | "
            f"{row['off_block']} | {row['action']} |"
        )
    lines.append("")
    lines.append("## Tensor sanity checks")
    lines.append("")
    checks = payload["tensor_sanity_checks"]
    lines.append("Allowed contractions are nonzero on a generic random draw:")
    lines.append("")
    lines.append("```text")
    for key, value in checks["allowed_nonzero_examples"].items():
        lines.append(f"{key} = {value:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("Contractions classified as zero vanish numerically:")
    lines.append("")
    lines.append("```text")
    for key, value in checks["vanishing_examples"].items():
        lines.append(f"{key} = {value:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("## Aligned-operator off-block checks")
    lines.append("")
    lines.append("Projectors are the broad PS projectors onto `(15,1,1)`, `(6,2,2)`,")
    lines.append("`(1,3,1)`, and `(1,1,3)`.  For an operator `O`, the displayed off-block")
    lines.append("norm is `sum_{p != q} P_p O P_q`.")
    lines.append("")
    lines.append("| operator | broad off max | broad off Frobenius | spectral norm | rank |")
    lines.append("|---|---:|---:|---:|---:|")
    for name, row in payload["operator_audit"].items():
        off = row["broad_PS_offblock"]
        lines.append(
            f"| `{name}` | {off['max_abs']:.3e} | {off['frobenius']:.3e} | "
            f"{row['spectral_norm']:.3e} | {row['rank_numeric']} |"
        )
    lines.append("")
    rho = payload["rho_ABC_risk"]
    lines.append("The important nontrivial result is the last row.  The distinct-adjoint")
    lines.append("structure-constant invariant `rho Tr A[B,C]` is allowed.  It does not mix")
    lines.append("the broad PS projectors for the aligned Cartan color vev, but it has")
    lines.append(f"operator norm `{rho['ad_A0_spectral_norm']:.3e}`.  A unit coefficient is")
    lines.append("therefore an order-one deformation of the mediator Hessian, not a harmless")
    lines.append("roundoff effect.")
    if rho.get("rough_rho_less_than_for_mediator_l2"):
        lines.append("")
        lines.append("Using the finite-mediator threshold residual only as a scale marker,")
        lines.append("")
        lines.append("```text")
        lines.append(
            "rho * ||ad_A0|| < ||P Delta_med|| would require "
            f"rho < {rho['rough_rho_less_than_for_mediator_l2']:.3e}."
        )
        lines.append("```")
        lines.append("")
        lines.append("This is not a physical exclusion bound; it is a warning that `rho_ABC`")
        lines.append("must be symmetry-forbidden or included in the next spectrum/RGE scan.")
    lines.append("")
    lines.append("## Next model-building decision")
    lines.append("")
    lines.append("There are two honest routes:")
    lines.append("")
    lines.append("1. **Full scan route:** add the allowed coefficient matrices")
    lines.append("   `m_ij`, `sigma_ij`, `tau_ij` and `rho_ABC` to the component Hessian,")
    lines.append("   then rerun the RGE/proton scan.")
    lines.append("2. **Symmetry route:** impose a mediator grading or R-symmetry that forbids")
    lines.append("   `rho_ABC` and the unwanted entries of `m_ij,sigma_ij,tau_ij`, while")
    lines.append("   leaving the projector-generating `AA`, `BC`, `AB`, and `AC` terms.")
    lines.append("")
    lines.append("The present local card is therefore structurally consistent only after one")
    lines.append("of these two choices is made explicitly.")
    (OUT / "untruncated_spin10_invariant_audit_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    table = invariant_table()
    write_csv(table)
    op = operator_audit()
    ref = load_reference_threshold()
    ad_norm = op["ad_A0_from_fABC"]["spectral_norm"]
    rough = None
    if ref:
        rough = ref["finite_mediator_projected_l2"] / ad_norm
    payload = {
        "note": "No web lookup used. Untruncated 54/210/45 invariant audit.",
        "field_content": {
            "S": "54_H, symmetric traceless two-tensor",
            "Phi": "210_H, antisymmetric four-form",
            "adjoints": ["A=45_Sigma", "B=45_med", "C=45_med"],
        },
        "invariant_table": table,
        "tensor_sanity_checks": tensor_sanity_checks(),
        "operator_audit": op,
        "reference_threshold": ref,
        "rho_ABC_risk": {
            "term": "rho_ABC Tr A[B,C]",
            "allowed_by_spin10": True,
            "must_scan_or_forbid": True,
            "ad_A0_spectral_norm": ad_norm,
            "rough_rho_less_than_for_mediator_l2": rough,
        },
        "verdict": {
            "forbidden_by_tensor_symmetry": [
                "54_H^2 210_H",
                "45_i 54_H^2 for one 54_H",
                "45_i 210_H^2 for one 210_H",
                "45_i 54_H 210_H",
                "45_i^3 and 45_i^2 45_j",
            ],
            "allowed_and_must_enter_next_scan_unless_symmetry_forbids": [
                "general 45_i45_j mass matrix",
                "54_H 45_i45_j matrix",
                "210_H 45_i45_j matrix",
                "rho_ABC Tr A[B,C]",
            ],
            "allowed_higgs_vacuum_terms": [
                "54_H^3",
                "54_H 210_H^2",
                "210_H^3",
            ],
        },
    }
    (OUT / "untruncated_spin10_invariant_audit_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)

    print("Untruncated Spin(10) invariant audit")
    print(f"  invariant rows: {len(table)}")
    print(f"  F54 broad off-block max: {op['F54_from_aligned_54']['broad_PS_offblock']['max_abs']:.3e}")
    print(f"  D210 weak broad off-block max: {op['D210_weak_volume']['broad_PS_offblock']['max_abs']:.3e}")
    print(f"  D210 color broad off-block max: {op['D210_color_fourform']['broad_PS_offblock']['max_abs']:.3e}")
    print(f"  fABC ad_A0 broad off-block max: {op['ad_A0_from_fABC']['broad_PS_offblock']['max_abs']:.3e}")
    print(f"  fABC ad_A0 spectral norm: {ad_norm:.3e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
