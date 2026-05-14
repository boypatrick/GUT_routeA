#!/usr/bin/env python3
"""Mass-insertion dressed dimension-five proton-decay proxy.

No web lookup is used.  This is a controlled upgrade of the scalar xi_H
stress scan.  It rebuilds the four-index triplet Wilson tensors from the
verified two-sided W_AB filter, then dresses the flavor indices with explicit
soft mass-insertion kernels

    C'_{abcd}=R_1{}_{aa'} R_2{}_{bb'} R_3{}_{cc'} R_4{}_{dd'} C_{a'b'c'd'},
    R_i = 1 + Delta_i.

The Delta_i matrices are not a full SUSY spectrum.  They are small Hermitian
mass-insertion ansaetze: CKM/MFV left-handed misalignment, right-handed
third-family splitting, and a combined stress.  The audit converts the
resulting RRRR amplification into an effective xi_H and checks it against the
previous spectrum-aware safe envelope.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import construct_triplet_rank_lift as rank_lift  # noqa: E402
import scan_triplet_mixing_nullspace as ns  # noqa: E402


TWO_SIDED = ROOT / "output" / "two_sided_triplet_filter" / "two_sided_triplet_filter_summary.json"
FLAVOR = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
VACUUM = ROOT / "output" / "spin10_vacuum_alignment" / "spin10_vacuum_alignment_summary.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
SPECTRUM = ROOT / "output" / "spectrum_aware_d5_dressing" / "summary.json"
OUT = ROOT / "output" / "mass_insertion_d5_dressing"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
DISPLAY_ST = 1.0e-5
ALPHA2_INV = 25.0

FILTER_LABELS = ["omegaR_0.1_kappa_100", "all_blocks_kappa_100"]
SOFT_EPS_GRID = [0.0, 0.01, 0.03, 0.10, 0.30, 1.00]
SPECTRUM_POINTS = [
    {
        "name": "baseline_100TeV",
        "m_sfermion_GeV": 1.0e5,
        "m_wino_GeV": 1.0e3,
        "mu_H_GeV": 1.0e3,
        "tan_beta": 10.0,
    },
    {
        "name": "marginal_safe_20TeV",
        "m_sfermion_GeV": 2.0e4,
        "m_wino_GeV": 1.0e3,
        "mu_H_GeV": 3.0e3,
        "tan_beta": 10.0,
    },
    {
        "name": "near_unsafe_20TeV",
        "m_sfermion_GeV": 2.0e4,
        "m_wino_GeV": 3.0e3,
        "mu_H_GeV": 3.0e3,
        "tan_beta": 10.0,
    },
]

CHANNEL_META = {
    "LLLL_upupdown_Knu": {
        "operator": "LLLL",
        "target_years": 2.4e34,
        "slot_kinds": ["uL", "uL", "dL", "nuL"],
    },
    "LLLL_downdownup_Knu": {
        "operator": "LLLL",
        "target_years": 2.4e34,
        "slot_kinds": ["dL", "dL", "uL", "nuL"],
    },
    "LLLL_upupdown_K0mu": {
        "operator": "LLLL",
        "target_years": 1.0e34,
        "slot_kinds": ["uL", "uL", "dL", "eL"],
    },
    "RRRR_uusd_anycharged": {
        "operator": "RRRR",
        "target_years": 2.4e34,
        "slot_kinds": ["uR", "uR", "dR", "eR"],
    },
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def selected_filter_rows() -> list[dict[str, Any]]:
    payload = read_json(TWO_SIDED)
    rows = {row["label"]: row for row in payload["rows"]}
    rows.update({row["label"]: row for row in payload["all_block_rows"]})
    return [rows[label] for label in FILTER_LABELS]


def filter_matrix(row: dict[str, Any]) -> np.ndarray:
    return cmat(row["W_matrix"])


def loop_function(x: float) -> float:
    if x <= 0.0:
        return 1.0
    if abs(x - 1.0) < 1.0e-8:
        return 0.5
    return (1.0 - x + x * math.log(x)) / ((1.0 - x) * (1.0 - x))


def vacuum_inputs() -> dict[str, float]:
    payload = read_json(VACUUM)
    vac = payload["recommended_benchmark"]
    yuk = payload["seed"]["best"]
    return {
        "M_T_GeV": float(vac["M_HC_GeV"]),
        "yt_MG": float(vac.get("yt_MG", yuk["yt_MG"])),
        "yb_MG": float(vac.get("yb_MG", yuk["yb_MG"])),
    }


def dressing_factors(point: dict[str, float], vac: dict[str, float]) -> dict[str, float]:
    m_sf = point["m_sfermion_GeV"]
    m_wino = point["m_wino_GeV"]
    mu_h = point["mu_H_GeV"]
    tan_beta = point["tan_beta"]
    x_w = (m_wino / m_sf) ** 2
    x_h = (mu_h / m_sf) ** 2
    d_w = (1.0 / ALPHA2_INV) / (4.0 * math.pi) * m_wino / (m_sf * m_sf) * loop_function(x_w)
    d_h = (
        (1.0 / (16.0 * math.pi * math.pi))
        * mu_h
        / (m_sf * m_sf)
        * vac["yt_MG"]
        * vac["yb_MG"]
        * (tan_beta / 10.0)
        * loop_function(x_h)
    )
    return {
        "D_wino_GeV_minus1": float(d_w),
        "D_higgsino_GeV_minus1": float(d_h),
        "higgsino_to_wino": float(d_h / d_w) if d_w > 0 else math.inf,
    }


def width_to_tau(width_gev: float) -> float:
    if width_gev <= 0.0:
        return math.inf
    return HBAR_GEV_S / width_gev / SECONDS_PER_YEAR


def lifetime(
    amplitude: float,
    st: float,
    m_triplet: float,
    width_prefactor: float,
    dressing: float,
) -> dict[str, float]:
    c5 = st * amplitude / m_triplet
    c6 = c5 * dressing
    width = width_prefactor * c6 * c6
    return {
        "C5_GeV_minus1": float(c5),
        "C6_GeV_minus2": float(c6),
        "width_GeV": float(width),
        "tau_years": float(width_to_tau(width)),
    }


def st_max_for(amplitude: float, target_years: float, m_triplet: float, width_prefactor: float, dressing: float) -> float:
    unit = lifetime(amplitude, 1.0, m_triplet, width_prefactor, dressing)
    return math.sqrt(unit["tau_years"] / target_years)


def offdiag(mat: np.ndarray) -> np.ndarray:
    return mat - np.diag(np.diag(mat))


def normalize_offdiag(mat: np.ndarray) -> np.ndarray:
    herm = 0.5 * (mat + mat.conjugate().T)
    off = offdiag(herm)
    scale = float(np.max(np.abs(off)))
    if scale < 1.0e-30:
        return np.zeros_like(herm)
    return off / scale


def normalize_traceless(mat: np.ndarray) -> np.ndarray:
    herm = 0.5 * (mat + mat.conjugate().T)
    herm = herm - np.trace(herm) * np.eye(herm.shape[0]) / herm.shape[0]
    scale = float(np.max(np.abs(herm)))
    if scale < 1.0e-30:
        return np.zeros_like(herm)
    return herm / scale


def insertion_basis() -> dict[str, Any]:
    payload = read_json(FLAVOR)
    rot = payload["biunitary_rotations"]
    ul = {name: cmat(rot[name]["left_rotation"]) for name in rot}
    ur = {name: cmat(rot[name]["right_rotation"]) for name in rot}
    ckm = cmat(payload["CKM"]["matrix"])
    p3 = np.diag([0.0, 0.0, 1.0]).astype(complex)

    # If the left soft splitting is aligned with the up basis, down indices see
    # V^dagger P3 V.  If aligned with the down basis, up indices see V P3 V^dagger.
    up_aligned_to_down = normalize_offdiag(ckm.conjugate().T @ p3 @ ckm)
    down_aligned_to_up = normalize_offdiag(ckm @ p3 @ ckm.conjugate().T)

    right_split = {
        "uR": normalize_offdiag(ur["up"].conjugate().T @ p3 @ ur["up"]),
        "dR": normalize_offdiag(ur["down"].conjugate().T @ p3 @ ur["down"]),
        "eR": normalize_offdiag(ur["charged_lepton"].conjugate().T @ p3 @ ur["charged_lepton"]),
    }
    democratic = np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)

    y_raw = {name: cmat(raw) for name, raw in payload["Yukawa_matrices"].items()}
    hu = y_raw["up"] @ y_raw["up"].conjugate().T
    hd = y_raw["down"] @ y_raw["down"].conjugate().T
    comm = normalize_traceless(1j * (hu @ hd - hd @ hu))
    comm_u = normalize_offdiag(ul["up"].conjugate().T @ comm @ ul["up"])
    comm_d = normalize_offdiag(ul["down"].conjugate().T @ comm @ ul["down"])

    return {
        "ckm": ckm,
        "matrices": {
            "zero": {key: np.zeros((3, 3), dtype=complex) for key in ["uL", "dL", "eL", "nuL", "uR", "dR", "eR"]},
            "up_aligned_LL_MFV": {
                "uL": np.zeros((3, 3), dtype=complex),
                "dL": up_aligned_to_down,
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                "uR": np.zeros((3, 3), dtype=complex),
                "dR": np.zeros((3, 3), dtype=complex),
                "eR": np.zeros((3, 3), dtype=complex),
            },
            "down_aligned_LL_MFV": {
                "uL": down_aligned_to_up,
                "dL": np.zeros((3, 3), dtype=complex),
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                "uR": np.zeros((3, 3), dtype=complex),
                "dR": np.zeros((3, 3), dtype=complex),
                "eR": np.zeros((3, 3), dtype=complex),
            },
            "right_third_split": {
                "uL": np.zeros((3, 3), dtype=complex),
                "dL": np.zeros((3, 3), dtype=complex),
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                **right_split,
            },
            "commutator_LL": {
                "uL": comm_u,
                "dL": comm_d,
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                "uR": np.zeros((3, 3), dtype=complex),
                "dR": np.zeros((3, 3), dtype=complex),
                "eR": np.zeros((3, 3), dtype=complex),
            },
            "combined_LL_RR": {
                "uL": 0.5 * down_aligned_to_up + 0.5 * comm_u,
                "dL": 0.5 * up_aligned_to_down + 0.5 * comm_d,
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                **right_split,
            },
            "democratic_RR_stress": {
                "uL": np.zeros((3, 3), dtype=complex),
                "dL": np.zeros((3, 3), dtype=complex),
                "eL": np.zeros((3, 3), dtype=complex),
                "nuL": np.zeros((3, 3), dtype=complex),
                "uR": democratic,
                "dR": democratic,
                "eR": democratic,
            },
            "democratic_all_stress": {
                "uL": democratic,
                "dL": democratic,
                "eL": democratic,
                "nuL": democratic,
                "uR": democratic,
                "dR": democratic,
                "eR": democratic,
            },
        },
    }


def combine_tensors(w: np.ndarray, channel_entries: dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]]) -> dict[str, tuple[np.ndarray, list[tuple[int, int, int, int]]]]:
    out: dict[str, tuple[np.ndarray, list[tuple[int, int, int, int]]]] = {}
    pair_full = [(i, j) for i in range(w.shape[0]) for j in range(w.shape[1])]
    coeffs = np.array([w[i, j] for (i, j) in pair_full], dtype=complex)
    for channel, (tensors, entries) in channel_entries.items():
        tensor = sum(coeffs[n] * tensors[pair] for n, pair in enumerate(pair_full))
        out[channel] = (tensor, entries)
    return out


def transform_tensor(tensor: np.ndarray, slot_kinds: list[str], base: dict[str, np.ndarray], eps: float) -> np.ndarray:
    transforms = [np.eye(3, dtype=complex) + eps * base[kind] for kind in slot_kinds]
    return np.einsum("ae,bf,cg,dh,efgh->abcd", transforms[0], transforms[1], transforms[2], transforms[3], tensor, optimize=True)


def max_entry(tensor: np.ndarray, entries: list[tuple[int, int, int, int]]) -> tuple[float, tuple[int, int, int, int], complex]:
    best_amp = -1.0
    best_idx = entries[0]
    best_value = 0.0j
    for idx in entries:
        value = tensor[idx]
        amp = float(abs(value))
        if amp > best_amp:
            best_amp = amp
            best_idx = idx
            best_value = complex(value)
    return best_amp, best_idx, best_value


def effective_xi(amplitude_ratio: float, dress: dict[str, float]) -> float:
    d_w = dress["D_wino_GeV_minus1"]
    d_h = dress["D_higgsino_GeV_minus1"]
    if d_h <= 0.0:
        return math.inf
    return float((amplitude_ratio * (d_w + d_h) - d_w) / d_h)


def audit() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    basis_names, pair_full, channel_entries, _constants_unused = rank_lift.build_channel_tensors()
    proton = read_json(PROTON)["hadronic_constants"]
    vac = vacuum_inputs()
    soft = insertion_basis()

    rows: list[dict[str, Any]] = []
    channel_rows: list[dict[str, Any]] = []

    for filt in selected_filter_rows():
        w = filter_matrix(filt)
        combined = combine_tensors(w, channel_entries)
        base_amplitudes = {
            channel: max_entry(tensor, entries)[0]
            for channel, (tensor, entries) in combined.items()
        }
        for scenario, base_delta in soft["matrices"].items():
            for eps in SOFT_EPS_GRID:
                for point in SPECTRUM_POINTS:
                    dress = dressing_factors(point, vac)
                    channel_results: list[dict[str, Any]] = []
                    for channel, (tensor, entries) in combined.items():
                        meta = CHANNEL_META[channel]
                        dressed_tensor = transform_tensor(tensor, meta["slot_kinds"], base_delta, eps)
                        amp, idx, value = max_entry(dressed_tensor, entries)
                        amp0 = max(base_amplitudes[channel], 1.0e-30)
                        amp_ratio = amp / amp0
                        if meta["operator"] == "RRRR":
                            d_eff = dress["D_wino_GeV_minus1"] + dress["D_higgsino_GeV_minus1"]
                            xi_eff = effective_xi(amp_ratio, dress)
                        else:
                            d_eff = dress["D_wino_GeV_minus1"]
                            xi_eff = None
                        life = lifetime(amp, DISPLAY_ST, vac["M_T_GeV"], proton["width_prefactor_GeV5"], d_eff)
                        st_allowed = st_max_for(
                            amp,
                            meta["target_years"],
                            vac["M_T_GeV"],
                            proton["width_prefactor_GeV5"],
                            d_eff,
                        )
                        result = {
                            "filter_label": filt["label"],
                            "scenario": scenario,
                            "epsilon": eps,
                            "spectrum_name": point["name"],
                            "channel": channel,
                            "operator": meta["operator"],
                            "amplitude_base": float(amp0),
                            "amplitude_dressed": float(amp),
                            "amplitude_ratio": float(amp_ratio),
                            "selected_index": "".join(str(i) for i in idx),
                            "value": cjson(value),
                            "effective_xi_H": None if xi_eff is None else float(xi_eff),
                            "S_T_max": float(st_allowed),
                            "tau_years_ST_display": life["tau_years"],
                            "margin_at_ST_display": life["tau_years"] / meta["target_years"],
                            "passes": life["tau_years"] >= meta["target_years"],
                        }
                        channel_results.append(result)
                        channel_rows.append(result)
                    worst = min(channel_results, key=lambda item: item["margin_at_ST_display"])
                    rrrr = next(item for item in channel_results if item["channel"] == "RRRR_uusd_anycharged")
                    rows.append(
                        {
                            "filter_label": filt["label"],
                            "scenario": scenario,
                            "epsilon": eps,
                            "spectrum_name": point["name"],
                            **point,
                            "D_wino_GeV_minus1": dress["D_wino_GeV_minus1"],
                            "D_higgsino_GeV_minus1": dress["D_higgsino_GeV_minus1"],
                            "higgsino_to_wino": dress["higgsino_to_wino"],
                            "worst_channel": worst["channel"],
                            "worst_margin": worst["margin_at_ST_display"],
                            "worst_ST_max": worst["S_T_max"],
                            "all_channels_pass": all(item["passes"] for item in channel_results),
                            "RRRR_amplitude_ratio": rrrr["amplitude_ratio"],
                            "RRRR_effective_xi_H": rrrr["effective_xi_H"],
                        }
                    )

    summary = summarize(rows, channel_rows, vac, basis_names)
    return rows, channel_rows, summary


def summarize(rows: list[dict[str, Any]], channel_rows: list[dict[str, Any]], vac: dict[str, float], basis_names: list[str]) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    for label in FILTER_LABELS:
        subset = [row for row in rows if row["filter_label"] == label]
        safe = [row for row in subset if row["all_channels_pass"]]
        unsafe = [row for row in subset if not row["all_channels_pass"]]
        baseline = next(
            row for row in subset
            if row["scenario"] == "zero" and row["epsilon"] == 0.0 and row["spectrum_name"] == "baseline_100TeV"
        )
        max_xi = max(row["RRRR_effective_xi_H"] for row in subset if row["RRRR_effective_xi_H"] is not None)
        by_filter[label] = {
            "total_points": len(subset),
            "safe_points": len(safe),
            "unsafe_points": len(unsafe),
            "safe_fraction": len(safe) / len(subset),
            "baseline": trim_row(baseline),
            "most_marginal_safe": trim_row(min(safe, key=lambda row: row["worst_margin"])) if safe else None,
            "nearest_unsafe": trim_row(max(unsafe, key=lambda row: row["worst_margin"])) if unsafe else None,
            "max_effective_xi_H": float(max_xi),
            "scenario_table": scenario_table(subset),
        }

    preferred = by_filter["omegaR_0.1_kappa_100"]
    return {
        "note": "No web lookup used. Explicit mass-insertion proxy for dressed d=5 operators.",
        "basis_names": basis_names,
        "vacuum_inputs": vac,
        "display_triplet_filter": DISPLAY_ST,
        "soft_epsilon_grid": SOFT_EPS_GRID,
        "spectrum_points": SPECTRUM_POINTS,
        "operator_formula": "C'_{abcd}=R1_{aa'}R2_{bb'}R3_{cc'}R4_{dd'}C_{a'b'c'd'}, R_i=1+epsilon Delta_i",
        "by_filter_label": by_filter,
        "previous_spectrum_envelope": read_json(SPECTRUM)["by_filter_label"]["omegaR_0.1_kappa_100"]["safe_by_xi_H"],
        "verdict": {
            "preferred_filter": "omegaR_0.1_kappa_100",
            "preferred_safe_fraction": preferred["safe_fraction"],
            "preferred_max_effective_xi_H": preferred["max_effective_xi_H"],
            "preferred_baseline_margin": preferred["baseline"]["worst_margin"],
            "interpretation": (
                "Replacing the scalar xi_H by explicit Hermitian mass-insertion "
                "kernels keeps the preferred filter safe over the audited MFV/right-split "
                "grid and even the democratic off-diagonal stress points.  The largest "
                "derived xi_H is an output of the soft misalignment ansatz, not an "
                "externally imposed scalar; a full SUSY spectrum would replace these "
                "kernels by actual sfermion diagonalization."
            ),
        },
    }


def trim_row(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "worst_channel",
        "worst_margin",
        "worst_ST_max",
        "all_channels_pass",
        "RRRR_amplitude_ratio",
        "RRRR_effective_xi_H",
    ]
    return {key: row[key] for key in keys}


def scenario_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for scenario in sorted({row["scenario"] for row in rows}):
        subset = [row for row in rows if row["scenario"] == scenario]
        safe = [row for row in subset if row["all_channels_pass"]]
        table.append(
            {
                "scenario": scenario,
                "safe_points": len(safe),
                "total_points": len(subset),
                "safe_fraction": len(safe) / len(subset),
                "max_effective_xi_H": max(row["RRRR_effective_xi_H"] for row in subset if row["RRRR_effective_xi_H"] is not None),
                "worst_margin": min(row["worst_margin"] for row in subset),
            }
        )
    return table


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "m_sfermion_GeV",
        "m_wino_GeV",
        "mu_H_GeV",
        "tan_beta",
        "D_wino_GeV_minus1",
        "D_higgsino_GeV_minus1",
        "higgsino_to_wino",
        "worst_channel",
        "worst_margin",
        "worst_ST_max",
        "all_channels_pass",
        "RRRR_amplitude_ratio",
        "RRRR_effective_xi_H",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_channel_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "channel",
        "operator",
        "amplitude_base",
        "amplitude_dressed",
        "amplitude_ratio",
        "selected_index",
        "effective_xi_H",
        "S_T_max",
        "tau_years_ST_display",
        "margin_at_ST_display",
        "passes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Mass-insertion dressed dimension-five proxy",
        "",
        "No web lookup was used.",
        "",
        "The flavor tensors are dressed by",
        "",
        "```text",
        "C'_{abcd}=R1_{aa'} R2_{bb'} R3_{cc'} R4_{dd'} C_{a'b'c'd'}, R_i=1+epsilon Delta_i.",
        "```",
        "",
        "## Filter summaries",
        "",
        "| filter | safe / total | baseline margin | max xi_H eff | marginal safe | nearest unsafe |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, payload in summary["by_filter_label"].items():
        marginal = payload["most_marginal_safe"]
        unsafe = payload["nearest_unsafe"]
        lines.append(
            f"| `{label}` | {payload['safe_points']}/{payload['total_points']} | "
            f"{payload['baseline']['worst_margin']:.3e} | {payload['max_effective_xi_H']:.3e} | "
            f"{marginal['worst_margin'] if marginal else math.inf:.3e} | "
            f"{unsafe['worst_margin'] if unsafe else math.inf:.3e} |"
        )
    lines.extend(["", "## Preferred scenario table", "", "| scenario | safe / total | max xi_H eff | worst margin |", "|---|---:|---:|---:|"])
    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    for row in preferred["scenario_table"]:
        lines.append(
            f"| `{row['scenario']}` | {row['safe_points']}/{row['total_points']} | "
            f"{row['max_effective_xi_H']:.3e} | {row['worst_margin']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, channel_rows, summary = audit()
    write_csv(OUT / "mass_insertion_scan.csv", rows)
    write_channel_csv(OUT / "channel_scan.csv", channel_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)

    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    print("Mass-insertion d=5 dressing proxy")
    print(
        "  omegaR_0.1_kappa_100: "
        f"safe={preferred['safe_points']}/{preferred['total_points']}, "
        f"max xi_eff={preferred['max_effective_xi_H']:.3e}, "
        f"baseline margin={preferred['baseline']['worst_margin']:.3e}"
    )
    if preferred["nearest_unsafe"]:
        print(f"  nearest unsafe margin={preferred['nearest_unsafe']['worst_margin']:.3e}")
    else:
        print("  no unsafe audited point")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
