#!/usr/bin/env python3
"""Generate a publication-card channel table for dimension-five proton decay.

No web lookup is used.  This audit is intentionally conservative about what is
and is not closed.  It rebuilds flavor rotations and channel-level mass-basis
Wilson entries from the single publication closure card, including the embedded
source-consistent H,F,G_A,G_B triplet-source tensors.  The table is still marked
PROXY_ONLY because the paper card does not yet contain final SUSY dressing
matrices, triplet mass eigenstate mixings, or channel-specific chiral reduction
operators.  Its role is to prevent the class-level crossed-120 d=5 margin from
being mistaken for a final proton-decay theorem.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import fit_two_kernel_flavor_then_d5 as two  # noqa: E402
import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
OUT = ROOT / "output" / "publication_channel_d5_tables"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path) if path.exists() else None,
            }
        )
    return rows


def cnum(raw: dict[str, float]) -> complex:
    return complex(float(raw["re"]), float(raw["im"]))


def cvec(raw: list[dict[str, float]]) -> np.ndarray:
    return np.array([cnum(item) for item in raw], dtype=complex)


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cnum(cell) for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def finite_or_none(x: float) -> float | None:
    return float(x) if math.isfinite(x) else None


def rebuild_source_data(card: dict[str, Any]) -> dict[str, Any]:
    source = card["triplet_tensor_source"]["model_data"]
    basis = two.load_veronese_basis()
    k = two.analytic_k()
    h_geo = fit.symmetric_from_coeff(cvec(source["H_coefficients"]), basis)
    f_geo = fit.symmetric_from_coeff(cvec(source["F_coefficients"]), basis)
    h = h_geo + cnum(source["epsilon_H"]) * k
    fmat = f_geo + cnum(source["epsilon_F"]) * k
    ga = fit.antisym_from_coeff(cvec(source["G_A_coefficients"]))
    gb = fit.antisym_from_coeff(cvec(source["G_B_coefficients"]))
    # The publication row is a doublet-only CKM repair with the crossed triplet
    # projector frozen.  Therefore triplet-source component bases must use the
    # source-consistent triplet mix stored with H,F,G_A,G_B, not the repaired
    # doublet mix from selected_row.
    mix = {key: cnum(val) for key, val in source["mixing_coefficients"].items()}
    return {
        "H_geo": h_geo,
        "F_geo": f_geo,
        "H": h,
        "F": fmat,
        "G_A": ga,
        "G_B": gb,
        "mix": mix,
        "doublet_mix": {key: cnum(val) for key, val in card["selected_row"]["mixing_coefficients"].items()},
    }


def rotations_for_card(card: dict[str, Any]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    y = {name: cmat(raw) for name, raw in card["selected_row"]["Yukawa_fit"].items()}
    return two.rotations_for(y)


def profile_by_name(name: str) -> two.TripletProfile:
    for profile in two.PROFILES:
        if profile.name == name:
            return profile
    raise KeyError(name)


def max_rrrr_channel(tensor: np.ndarray, d_flavor: int, charged: int | None) -> tuple[float, tuple[int, int, int, int]]:
    """Tensor order is (u_a,u_b,d_c,e_d), matching d5.tensor_rrrr."""
    best_amp = -1.0
    best_idx = (0, 0, d_flavor, 0)
    charged_set = range(3) if charged is None else [charged]
    for ell in charged_set:
        for u_pair in [(0, 0), (0, 1), (1, 0)]:
            idx = (u_pair[0], u_pair[1], d_flavor, ell)
            amp = float(abs(tensor[idx]))
            if amp > best_amp:
                best_amp = amp
                best_idx = idx
    return best_amp, best_idx


def constants_for(channel_class: str, chiral: dict[str, Any]) -> dict[str, float]:
    key = "e_pi" if channel_class in {"e_pi", "mu_pi"} else channel_class
    return {"width_prefactor_GeV5": float(chiral["central_width_prefactors"][key])}


def lifetime_from_amp(amplitude: float, channel_class: str, chiral: dict[str, Any]) -> dict[str, float | None]:
    base = d5.lifetime_from_amplitude(amplitude, constants_for(channel_class, chiral))
    width_pref = constants_for(channel_class, chiral)["width_prefactor_GeV5"]
    return {
        "central_width_prefactor_GeV5": width_pref,
        "tau_years_ST_1": finite_or_none(base["tau_years"]),
        "C6_dressed_GeV_minus2_ST_1": finite_or_none(base["C6_dressed_GeV_minus2"]),
    }


def filter_required_from_tau(tau_years: float | None, bound_years: float) -> float | None:
    if tau_years is None:
        return None
    return math.sqrt(tau_years / bound_years)


def channel_row(
    hypothesis: str,
    operator: str,
    channel: str,
    channel_class: str,
    amplitude: float,
    index: tuple[int, int, int, int],
    class_margin: float,
    leakage_ratio: float,
    chiral: dict[str, Any],
    note: str,
) -> dict[str, Any]:
    life = lifetime_from_amp(amplitude, channel_class, chiral)
    tau = life["tau_years_ST_1"]
    return {
        "hypothesis": hypothesis,
        "operator": operator,
        "channel": channel,
        "channel_class": channel_class,
        "selected_index": "".join(str(i) for i in index),
        "raw_proxy_amplitude": float(amplitude),
        "central_width_prefactor_GeV5": life["central_width_prefactor_GeV5"],
        "tau_years_ST_1_central": tau,
        "S_T_required_tau_2p4e34_central": filter_required_from_tau(tau, 2.4e34),
        "S_T_required_tau_1e35_central": filter_required_from_tau(tau, 1.0e35),
        "crossed_finite_lift_leakage_ratio": float(leakage_ratio),
        "class_level_future_margin_1e35": float(class_margin),
        "class_level_passes_1e35": bool(class_margin >= 1.0),
        "status": "PROXY_ONLY_CLASS_LEVEL_PASS" if class_margin >= 1.0 else "PROXY_ONLY_CLASS_LEVEL_FAIL",
        "note": note,
    }


def tensor_rows(
    hypothesis: str,
    y_qq: np.ndarray,
    y_ql: np.ndarray,
    y_ue: np.ndarray,
    y_ud: np.ndarray,
    ul: dict[str, np.ndarray],
    ur: dict[str, np.ndarray],
    u_nu: np.ndarray,
    crossed_row: dict[str, float],
    chiral: dict[str, Any],
    note: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cl_nu = d5.tensor_llll(d5.sym(y_qq), y_ql, ul["up"], ul["up"], ul["down"], u_nu)
    cl_ch = d5.tensor_llll(d5.sym(y_qq), y_ql, ul["up"], ul["up"], ul["down"], ul["charged_lepton"])
    cr = d5.tensor_rrrr(y_ue, y_ud, ur["up"], ur["charged_lepton"], ur["up"], ur["down"])

    llll_margin = float(crossed_row["LLLL_margin_1e35_replayed"])
    rrrr_margin = float(crossed_row["RRRR_margin_1e35_replayed"])
    llll_leak = float(crossed_row["LLLL_leakage_ratio"])
    rrrr_leak = float(crossed_row["RRRR_leakage_ratio"])

    amp, idx = d5.max_perm_entry(cl_nu, (0, 0, 1), None)
    rows.append(channel_row(hypothesis, "LLLL", "p_to_Kplus_nubar", "Knu", amp, idx, llll_margin, llll_leak, chiral, note))
    amp, idx = d5.max_perm_entry(cl_ch, (0, 0, 1), 1)
    rows.append(channel_row(hypothesis, "LLLL", "p_to_K0_muplus", "K0mu", amp, idx, llll_margin, llll_leak, chiral, note))
    amp, idx = d5.max_perm_entry(cl_ch, (0, 0, 0), 0)
    rows.append(channel_row(hypothesis, "LLLL", "p_to_eplus_pi0", "e_pi", amp, idx, llll_margin, llll_leak, chiral, note))
    amp, idx = d5.max_perm_entry(cl_ch, (0, 0, 0), 1)
    rows.append(channel_row(hypothesis, "LLLL", "p_to_muplus_pi0", "mu_pi", amp, idx, llll_margin, llll_leak, chiral, note))

    amp, idx = max_rrrr_channel(cr, d_flavor=1, charged=None)
    rows.append(channel_row(hypothesis, "RRRR", "p_to_K_like_any_charged", "K0mu", amp, idx, rrrr_margin, rrrr_leak, chiral, note))
    amp, idx = max_rrrr_channel(cr, d_flavor=0, charged=0)
    rows.append(channel_row(hypothesis, "RRRR", "p_to_eplus_pi0", "e_pi", amp, idx, rrrr_margin, rrrr_leak, chiral, note))
    amp, idx = max_rrrr_channel(cr, d_flavor=0, charged=1)
    rows.append(channel_row(hypothesis, "RRRR", "p_to_muplus_pi0", "mu_pi", amp, idx, rrrr_margin, rrrr_leak, chiral, note))
    return rows


def build() -> dict[str, Any]:
    card = read_json(CARD)
    chiral = read_json(CHIRAL)
    data = rebuild_source_data(card)
    ul, ur, u_nu = rotations_for_card(card)
    crossed_row = next(
        row
        for row in read_json(ROOT / "output" / "crossed_120_triplet_projector" / "summary.json")["finite_leakage_replay"]
        if float(row["kappa"]) == float(card["crossed_projector"]["kappa"])
    )

    profile = profile_by_name("F_quarter_phase")
    y_qq_profile, y_ql_profile = two.triplet_matrices(data, profile)
    y_qq_crossed = d5.scale_to_largest_singular(data["G_A"], 0.60)
    y_ql_crossed = d5.scale_to_largest_singular(data["G_B"], 0.024)

    rows: list[dict[str, Any]] = []
    rows.extend(
        tensor_rows(
            "finite_Clebsch_profile_F_quarter_phase",
            y_qq_profile,
            y_ql_profile,
            y_qq_profile,
            y_ql_profile,
            ul,
            ur,
            u_nu,
            crossed_row,
            chiral,
            "Triplet matrices are the finite Clebsch profile rebuilt from the publication card.",
        )
    )
    rows.extend(
        tensor_rows(
            "crossed_pure_120_source_proxy",
            y_qq_crossed,
            y_ql_crossed,
            y_qq_crossed,
            y_ql_crossed,
            ul,
            ur,
            u_nu,
            crossed_row,
            chiral,
            "Pure crossed 120_A/120_B source proxy; LLLL uses sym(Y_QQ), so exact antisymmetric QQ nulls are visible.",
        )
    )

    rows_sorted = sorted(rows, key=lambda row: (row["operator"], row["channel"], row["hypothesis"]))
    exact_complete = all(
        key in card
        for key in [
            "final_triplet_mass_eigenstates",
            "final_susy_dressing_matrices",
            "final_chiral_reduction_operators",
            "C5L_tensor",
            "C5R_tensor",
        ]
    )
    worst_margin = min(float(row["class_level_future_margin_1e35"]) for row in rows_sorted)
    min_tau_values = [row["tau_years_ST_1_central"] for row in rows_sorted if row["tau_years_ST_1_central"] is not None]
    payload = {
        "note": "No web lookup used. Channel-level d=5 table scaffold from the publication closure card.",
        "input_manifest": manifest(
            [
                CARD,
                CHIRAL,
                ROOT / "output" / "crossed_120_triplet_projector" / "summary.json",
            ]
        ),
        "operator_definitions": {
            "LLLL": "C_L^{abcd}=Y_QQ^{ij}Y_QL^{kl}U_{q1}^{ia}U_{q2}^{jb}U_{q3}^{kc}U_l^{ld}",
            "RRRR": "C_R^{abcd}=Y_UE^{ij}Y_UD^{kl}U_u^{ia}U_e^{jd}U_u^{kb}U_d^{lc}",
            "index_order_LLLL": "(q_a,q_b,q_c,lepton_d)",
            "index_order_RRRR": "(u_a,u_b,d_c,e_d)",
        },
        "source_norms": {
            "G_A_antisym_residual": float(np.linalg.norm(data["G_A"] + data["G_A"].T)),
            "G_B_antisym_residual": float(np.linalg.norm(data["G_B"] + data["G_B"].T)),
            "profile_Y_QQ_largest_singular": float(np.linalg.svd(y_qq_profile, compute_uv=False)[0]),
            "profile_Y_QL_largest_singular": float(np.linalg.svd(y_ql_profile, compute_uv=False)[0]),
            "crossed_Y_QQ_largest_singular": float(np.linalg.svd(y_qq_crossed, compute_uv=False)[0]),
            "crossed_Y_QL_largest_singular": float(np.linalg.svd(y_ql_crossed, compute_uv=False)[0]),
        },
        "rows": rows_sorted,
        "verdict": {
            "proxy_channel_rows": len(rows_sorted),
            "exact_channel_wilson_complete": exact_complete,
            "all_class_level_margins_pass_1e35": all(row["class_level_passes_1e35"] for row in rows_sorted),
            "worst_class_level_margin_1e35": worst_margin,
            "min_tau_years_ST_1_central_proxy": min(min_tau_values) if min_tau_values else None,
            "missing_for_publication_exactness": [
                "triplet mass-eigenstate mixing matrices",
                "short-distance and long-distance dressing matrices by channel",
                "final C5L and C5R tensor arrays after the crossed inverse block",
                "channel-specific chiral/lattice reduction operators with cited input uncertainties",
            ],
            "interpretation": (
                "The card now generates channel rows for K+nu, K0mu, e+pi0, and mu+pi0, "
                "and the crossed-120 class-level kappa=30 margin is enormous for both "
                "LLLL and RRRR rows.  However this is still a proxy table: exact channel "
                "Wilson tensors and final SUSY/chiral dressing data are not yet encoded, "
                "so the OPEN publication-level d=5 item remains open."
            ),
        },
    }
    return payload


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "hypothesis",
        "operator",
        "channel",
        "channel_class",
        "selected_index",
        "raw_proxy_amplitude",
        "central_width_prefactor_GeV5",
        "tau_years_ST_1_central",
        "S_T_required_tau_2p4e34_central",
        "S_T_required_tau_1e35_central",
        "crossed_finite_lift_leakage_ratio",
        "class_level_future_margin_1e35",
        "class_level_passes_1e35",
        "status",
        "note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_manifest_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = ["path", "exists", "size_bytes", "sha256"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = [
        "# Publication-card channel d=5 table scaffold",
        "",
        "No web lookup was used.",
        "",
        "The table is generated from the single publication closure card.  It is",
        "still marked `PROXY_ONLY`: the final channel Wilson tensors and full",
        "SUSY/chiral dressing operators are not yet encoded in the card.",
        "",
        "## Class-level margins",
        "",
        f"Worst crossed-120 class-level margin at `1e35 yr`: `{payload['verdict']['worst_class_level_margin_1e35']:.6e}`.",
        f"Exact channel Wilson complete: `{payload['verdict']['exact_channel_wilson_complete']}`.",
        "",
        "## Channel rows",
        "",
        "| hypothesis | op | channel | index | amp | tau(ST=1) | class margin |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in payload["rows"]:
        tau = row["tau_years_ST_1_central"]
        tau_s = "inf" if tau is None else f"{tau:.3e}"
        lines.append(
            f"| `{row['hypothesis']}` | `{row['operator']}` | `{row['channel']}` | "
            f"`{row['selected_index']}` | {row['raw_proxy_amplitude']:.3e} | "
            f"{tau_s} | {row['class_level_future_margin_1e35']:.3e} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        payload["verdict"]["interpretation"],
        "",
        "Missing before a publication-final d=5 theorem:",
        "",
    ]
    for item in payload["verdict"]["missing_for_publication_exactness"]:
        lines.append(f"- {item}")
    lines.append("")
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build()
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "channel_table.csv", payload["rows"])
    write_manifest_csv(OUT / "input_manifest.csv", payload["input_manifest"])
    write_report(payload)
    verdict = payload["verdict"]
    print("Publication-card channel d=5 table scaffold")
    print(f"  rows: {verdict['proxy_channel_rows']}")
    print(f"  exact channel Wilson complete: {verdict['exact_channel_wilson_complete']}")
    print(f"  all class-level margins pass 1e35: {verdict['all_class_level_margins_pass_1e35']}")
    print(f"  worst class-level margin: {verdict['worst_class_level_margin_1e35']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
