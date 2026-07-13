#!/usr/bin/env python3
"""Build mass-basis dimension-five Wilson-tensor proxies.

No web lookup is used.  This script upgrades the scalar d=5 leakage proxy to a
four-index mass-basis tensor ledger.  It is still not the final SUSY-GUT
proton-decay calculation: triplet mixing, full wino/higgsino dressing, and
chiral reduction are represented by one common filter S_T.  The purpose is to
make the flavor part reproducible and to identify which triplet-Clebsch
hypotheses are immediately dangerous.
"""

from __future__ import annotations

import csv
import json
import math
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "dimension5_wilson_tensors"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def sym(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat + mat.T)


def anti(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat - mat.T)


def scale_to_largest_singular(mat: np.ndarray, target: float) -> np.ndarray:
    singulars = np.linalg.svd(mat, compute_uv=False)
    return mat * (target / max(float(singulars[0]), 1.0e-30))


def load_inputs() -> tuple[dict[str, Any], dict[str, float]]:
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    proton = json.loads(PROTON.read_text(encoding="utf-8"))
    return payload, proton["hadronic_constants"]


def lifetime_from_amplitude(
    amplitude: float,
    constants: dict[str, float],
    triplet_filter: float = 1.0,
    m_triplet_gev: float = 1.0e16,
    m_wino_gev: float = 1.0e3,
    m_sfermion_gev: float = 1.0e5,
    alpha2_inv: float = 25.0,
) -> dict[str, float]:
    loop = (1.0 / alpha2_inv) / (4.0 * math.pi)
    c5 = triplet_filter * amplitude / m_triplet_gev
    dressing = loop * m_wino_gev / (m_sfermion_gev * m_sfermion_gev)
    c6 = c5 * dressing
    width = constants["width_prefactor_GeV5"] * c6 * c6
    tau = math.inf if width <= 0.0 else HBAR_GEV_S / width / SECONDS_PER_YEAR
    return {
        "amplitude": float(amplitude),
        "triplet_filter": float(triplet_filter),
        "M_T_GeV": float(m_triplet_gev),
        "m_wino_GeV": float(m_wino_gev),
        "m_sfermion_GeV": float(m_sfermion_gev),
        "loop_alpha2_over_4pi": float(loop),
        "C5_GeV_minus1": float(c5),
        "dressing_GeV_minus1": float(dressing),
        "C6_dressed_GeV_minus2": float(c6),
        "tau_years": float(tau),
    }


def filter_required(row_at_st1: dict[str, float], tau_bound_years: float) -> float:
    if row_at_st1["amplitude"] <= 0.0:
        return math.inf
    return math.sqrt(row_at_st1["tau_years"] / tau_bound_years)


def tensor_llll(
    y_qq: np.ndarray,
    y_ql: np.ndarray,
    uq1: np.ndarray,
    uq2: np.ndarray,
    uq3: np.ndarray,
    ul: np.ndarray,
) -> np.ndarray:
    return np.einsum("ij,kl,ia,jb,kc,ld->abcd", y_qq, y_ql, uq1, uq2, uq3, ul, optimize=True)


def tensor_rrrr(
    y_ue: np.ndarray,
    y_ud: np.ndarray,
    uu1: np.ndarray,
    ue: np.ndarray,
    uu2: np.ndarray,
    ud: np.ndarray,
) -> np.ndarray:
    return np.einsum("ij,kl,ia,jd,kb,lc->abcd", y_ue, y_ud, uu1, ue, uu2, ud, optimize=True)


def max_perm_entry(tensor: np.ndarray, q_flavors: tuple[int, int, int], lepton: int | None) -> tuple[float, tuple[int, int, int, int]]:
    best_amp = -1.0
    best_idx = (0, 0, 0, 0)
    for qidx in sorted(set(permutations(q_flavors))):
        leptons = range(3) if lepton is None else [lepton]
        for ell in leptons:
            idx = (qidx[0], qidx[1], qidx[2], ell)
            amp = float(abs(tensor[idx]))
            if amp > best_amp:
                best_amp = amp
                best_idx = idx
    return best_amp, best_idx


def max_rrrr_uusd(tensor: np.ndarray, charged: int | None) -> tuple[float, tuple[int, int, int, int]]:
    # Tensor order is u^c_a e^c_d u^c_b d^c_c, stored as (a,b,c,d).
    best_amp = -1.0
    best_idx = (0, 0, 0, 0)
    charged_set = range(3) if charged is None else [charged]
    for ell in charged_set:
        for u_pair in [(0, 0), (0, 1), (1, 0)]:
            idx = (u_pair[0], u_pair[1], 1, ell)
            amp = float(abs(tensor[idx]))
            if amp > best_amp:
                best_amp = amp
                best_idx = idx
    return best_amp, best_idx


def row_with_lifetime(
    hypothesis: str,
    operator: str,
    channel: str,
    amplitude: float,
    index: tuple[int, int, int, int],
    constants: dict[str, float],
    note: str,
) -> dict[str, Any]:
    base = lifetime_from_amplitude(amplitude, constants)
    return {
        "hypothesis": hypothesis,
        "operator": operator,
        "channel_proxy": channel,
        "selected_index": "".join(str(i) for i in index),
        "amplitude": base["amplitude"],
        "tau_years_ST_1": base["tau_years"],
        "S_T_required_tau_1e34": filter_required(base, 1.0e34),
        "S_T_required_tau_2p4e34": filter_required(base, 2.4e34),
        "S_T_required_tau_1e35": filter_required(base, 1.0e35),
        "C6_dressed_GeV_minus2_ST_1": base["C6_dressed_GeV_minus2"],
        "note": note,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload, constants = load_inputs()

    y_raw = {name: cmat(raw) for name, raw in payload["Yukawa_matrices"].items()}
    rot = payload["biunitary_rotations"]
    ul = {name: cmat(rot[name]["left_rotation"]) for name in rot}
    ur = {name: cmat(rot[name]["right_rotation"]) for name in rot}
    u_nu = cmat(payload["PMNS_convention"]["neutrino_left_rotation_target"])

    # Local item-4 normalization.  The charged-lepton value is only used in
    # RRRR proxy hypotheses and is kept explicit in the output.
    y_phys = {
        "up": scale_to_largest_singular(y_raw["up"], 0.60),
        "down": scale_to_largest_singular(y_raw["down"], 0.024),
        "charged_lepton": scale_to_largest_singular(y_raw["charged_lepton"], 0.010),
        "neutrino_dirac": scale_to_largest_singular(y_raw["neutrino_dirac"], 0.35),
    }

    k_tr = np.array(
        [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        dtype=complex,
    ) / math.sqrt(3.0)
    k_tr = scale_to_largest_singular(k_tr, 0.60)

    hypotheses = [
        {
            "name": "minimal_down_aligned",
            "Y_QQ": sym(y_phys["up"]),
            "Y_QL": y_phys["down"],
            "Y_UE": sym(y_phys["up"]),
            "Y_UD": y_phys["down"],
            "note": "SU(5)-like proxy with Y_QQ from symmetric up texture and Y_QL from down texture.",
        },
        {
            "name": "lepton_transposed",
            "Y_QQ": sym(y_phys["up"]),
            "Y_QL": y_phys["charged_lepton"].T,
            "Y_UE": sym(y_phys["up"]),
            "Y_UD": y_phys["charged_lepton"].T,
            "note": "Triplet QL coupling aligned with Y_e^T instead of Y_d.",
        },
        {
            "name": "geometric_average_10bar5",
            "Y_QQ": sym(y_phys["up"]),
            "Y_QL": 0.5 * (y_phys["down"] + y_phys["charged_lepton"].T),
            "Y_UE": sym(y_phys["up"]),
            "Y_UD": 0.5 * (y_phys["down"] + y_phys["charged_lepton"].T),
            "note": "A symmetric interpolation between the fitted down and charged-lepton 10*5bar textures.",
        },
        {
            "name": "antisymmetric_up_piece",
            "Y_QQ": anti(y_phys["up"]),
            "Y_QL": y_phys["down"],
            "Y_UE": anti(y_phys["up"]),
            "Y_UD": y_phys["down"],
            "note": "Stress test for a 120-like antisymmetric triplet piece.",
        },
        {
            "name": "transvectant_contact_QQ",
            "Y_QQ": k_tr,
            "Y_QL": y_phys["down"],
            "Y_UE": k_tr,
            "Y_UD": y_phys["down"],
            "note": "Novel contact-filter hypothesis: the QQ triplet source is the spin-zero transvectant.",
        },
    ]

    rows: list[dict[str, Any]] = []
    tensor_summaries: dict[str, Any] = {}
    for hyp in hypotheses:
        name = hyp["name"]
        # Two assignments expose the otherwise hidden SU(2)-doublet basis choice.
        cl_upupdown_nu = tensor_llll(hyp["Y_QQ"], hyp["Y_QL"], ul["up"], ul["up"], ul["down"], u_nu)
        cl_downdownup_nu = tensor_llll(hyp["Y_QQ"], hyp["Y_QL"], ul["down"], ul["down"], ul["up"], u_nu)
        cl_upupdown_e = tensor_llll(hyp["Y_QQ"], hyp["Y_QL"], ul["up"], ul["up"], ul["down"], ul["charged_lepton"])
        cr = tensor_rrrr(hyp["Y_UE"], hyp["Y_UD"], ur["up"], ur["charged_lepton"], ur["up"], ur["down"])

        for label, tensor, lepton_note in [
            ("LLLL_upupdown_Knu", cl_upupdown_nu, "neutrino mass-basis lepton"),
            ("LLLL_downdownup_Knu", cl_downdownup_nu, "neutrino mass-basis lepton"),
        ]:
            amp, idx = max_perm_entry(tensor, (0, 0, 1), None)
            rows.append(row_with_lifetime(name, "LLLL", label, amp, idx, constants, hyp["note"] + " " + lepton_note))

        amp_mu, idx_mu = max_perm_entry(cl_upupdown_e, (0, 0, 1), 1)
        rows.append(row_with_lifetime(name, "LLLL", "LLLL_upupdown_K0mu", amp_mu, idx_mu, constants, hyp["note"] + " charged-lepton proxy"))

        amp_r, idx_r = max_rrrr_uusd(cr, None)
        rows.append(row_with_lifetime(name, "RRRR", "RRRR_uusd_anycharged", amp_r, idx_r, constants, hyp["note"] + " right-handed proxy"))

        tensor_summaries[name] = {
            "norms": {
                "Y_QQ_F": float(np.linalg.norm(hyp["Y_QQ"])),
                "Y_QL_F": float(np.linalg.norm(hyp["Y_QL"])),
                "C_L_upupdown_nu_F": float(np.linalg.norm(cl_upupdown_nu)),
                "C_L_downdownup_nu_F": float(np.linalg.norm(cl_downdownup_nu)),
                "C_R_F": float(np.linalg.norm(cr)),
            },
            "max_abs": {
                "C_L_upupdown_nu": float(np.max(np.abs(cl_upupdown_nu))),
                "C_L_downdownup_nu": float(np.max(np.abs(cl_downdownup_nu))),
                "C_R": float(np.max(np.abs(cr))),
            },
        }

    rows_sorted = sorted(rows, key=lambda row: row["amplitude"], reverse=True)
    most_dangerous = rows_sorted[0]
    viable_at_st1e_minus5 = [
        row for row in rows
        if 1.0e-5 <= row["S_T_required_tau_2p4e34"]
    ]

    with (OUT / "dimension5_wilson_proxy.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows_sorted)

    payload_out = {
        "note": "No web lookup used. Four-index mass-basis d=5 proxy tensors; not a full chiral proton-decay calculation.",
        "input": str(INPUT),
        "normalization": {
            "y_t": 0.60,
            "y_b": 0.024,
            "y_tau_proxy": 0.010,
            "y_nuD_proxy": 0.35,
            "M_T_GeV": 1.0e16,
            "m_wino_GeV": 1.0e3,
            "m_sfermion_GeV": 1.0e5,
            "alpha2_inv": 25.0,
        },
        "operator_definitions": {
            "LLLL": "C_L^{abcd}=Y_QQ^{ij}Y_QL^{kl}U_1^{ia}U_2^{jb}U_3^{kc}U_L^{ld}",
            "RRRR": "C_R^{abcd}=Y_UE^{ij}Y_UD^{kl}U_u^{ia}U_e^{jd}U_u^{kb}U_d^{lc}",
            "index_order_LLLL": "(q_a,q_b,q_c,lepton_d)",
            "index_order_RRRR": "(u_a,u_b,d_c,e_d), built from u^c e^c u^c d^c",
        },
        "tensor_summaries": tensor_summaries,
        "rows": rows_sorted,
        "verdict": {
            "most_dangerous": most_dangerous,
            "count_safe_for_ST_1e_minus5_against_2p4e34": len(viable_at_st1e_minus5),
            "total_rows": len(rows),
            "interpretation": (
                "The scalar CKM proxy was optimistic for some explicit tensor assignments. "
                "A fitted triplet sector must therefore specify its Clebsch orientation; "
                "the transvectant-contact QQ hypothesis is testable because it changes the "
                "dangerous tensor entries rather than merely multiplying all channels by S_T."
            ),
        },
    }
    (OUT / "dimension5_wilson_tensors.json").write_text(json.dumps(payload_out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Dimension-five mass-basis Wilson tensor audit",
        "",
        "No web lookup was used.",
        "",
        "The entries below use the common item-4 dressing proxy with `M_T=1e16 GeV`,",
        "`m_wino=1e3 GeV`, `m_sfermion=1e5 GeV`, and `alpha2^{-1}=25`.",
        "",
        "| hypothesis | operator | channel | index | amp | S_T max 1e34 | S_T max 2.4e34 |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in rows_sorted:
        report.append(
            f"| `{row['hypothesis']}` | `{row['operator']}` | `{row['channel_proxy']}` | "
            f"`{row['selected_index']}` | {row['amplitude']:.3e} | "
            f"{row['S_T_required_tau_1e34']:.3e} | {row['S_T_required_tau_2p4e34']:.3e} |"
        )
    report.extend(
        [
            "",
            "## Verdict",
            "",
            payload_out["verdict"]["interpretation"],
            "",
            f"Most dangerous row: `{most_dangerous['hypothesis']}` / "
            f"`{most_dangerous['channel_proxy']}` with amplitude "
            f"`{most_dangerous['amplitude']:.6e}` and "
            f"`S_T^max(2.4e34 yr)={most_dangerous['S_T_required_tau_2p4e34']:.6e}`.",
            "",
        ]
    )
    (OUT / "dimension5_wilson_tensors_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Dimension-five Wilson tensor audit")
    print(f"  rows: {len(rows_sorted)}")
    print(f"  most dangerous: {most_dangerous['hypothesis']} {most_dangerous['channel_proxy']}")
    print(f"  amplitude: {most_dangerous['amplitude']:.6e}")
    print(f"  S_T max for 2.4e34 yr: {most_dangerous['S_T_required_tau_2p4e34']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
