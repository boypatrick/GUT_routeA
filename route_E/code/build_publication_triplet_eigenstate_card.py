#!/usr/bin/env python3
"""Build a source-basis triplet eigenstate card from the publication card.

No web lookup is used.

This is the next step after the publication channel scaffold.  It turns the
crossed 120_A/120_B class-level inverse block into an explicit finite
source-basis matrix W, computes its singular/eigenstate data, and contracts it
into four-index C5L/C5R proxy tensors:

    C5L = sum_{I,J} W_{IJ} C5L[I,J],
    C5R = sum_{I,J} W_{IJ} C5R[I,J].

The output is still not a publication-final proton-decay theorem.  It does not
include full short/long-distance SUSY dressing matrices or final chiral
operator reductions.  Its purpose is narrower and falsifiable: check whether
the crossed inverse block can be represented as an actual tensor-valued source
map instead of only as a class-level leakage ratio.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import null_space


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_publication_channel_d5_tables as pubd5  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
TRIPLET_120 = ROOT / "output" / "triplet_120_rrrr_mass_matrix" / "summary.json"
CROSSED = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
OUT = ROOT / "output" / "publication_triplet_eigenstate_card"
REFERENCE_TRIPLET_FILTER = 1.0e-5


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict[str, Any]]:
    return [
        {
            "path": str(path.relative_to(ROOT)),
            "exists": path.exists(),
            "size_bytes": path.stat().st_size if path.exists() else None,
            "sha256": sha256(path) if path.exists() else None,
        }
        for path in paths
    ]


def cnum(raw: dict[str, float]) -> complex:
    return complex(float(raw["re"]), float(raw["im"]))


def cvec(raw: list[dict[str, float]]) -> np.ndarray:
    return np.array([cnum(item) for item in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(mat: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in mat]


def tensor4_json(tensor: np.ndarray) -> list[list[list[list[dict[str, float]]]]]:
    return [
        [
            [
                [cjson(tensor[a, b, c, d]) for d in range(tensor.shape[3])]
                for c in range(tensor.shape[2])
            ]
            for b in range(tensor.shape[1])
        ]
        for a in range(tensor.shape[0])
    ]


def finite_or_none(x: float) -> float | None:
    return float(x) if math.isfinite(x) else None


def normalize(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm <= 1.0e-30:
        raise ValueError("cannot normalize zero vector")
    return vec / norm


def complete_unitary(first: np.ndarray) -> np.ndarray:
    e = normalize(first.astype(complex))
    rest = null_space(e.conjugate().reshape(1, -1))
    return np.column_stack([e, rest])


def component_bases(data: dict[str, Any]) -> tuple[list[np.ndarray], list[np.ndarray], list[str]]:
    h = data["H"]
    fmat = data["F"]
    ga = data["G_A"]
    gb = data["G_B"]
    mix = data["mix"]
    names = ["H", "F", "G_A", "G_B"]
    b10 = [h, mix["r_u"] * fmat, mix["a_u"] * ga, mix["b_u"] * gb]
    b5 = [
        mix["r_d"] * h,
        mix["r_d"] * fmat,
        mix["r_d"] * mix["a_d"] * ga,
        mix["r_d"] * mix["b_d"] * gb,
    ]
    return b10, b5, names


def inverse_blocks(p: np.ndarray, q: np.ndarray, kappa: float) -> dict[str, np.ndarray]:
    u = complete_unitary(p)
    v = complete_unitary(q)
    eps = 1.0 / kappa
    return {
        "U_left": u,
        "V_right": v,
        "rank_one": u @ np.diag([1.0, 0.0, 0.0, 0.0]) @ v.conjugate().T,
        "finite": u @ np.diag([1.0, eps, eps, eps]) @ v.conjugate().T,
    }


def contract_tensors(
    w: np.ndarray,
    b10: list[np.ndarray],
    b5: list[np.ndarray],
    ul: dict[str, np.ndarray],
    ur: dict[str, np.ndarray],
    u_nu: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    c_l = np.zeros((3, 3, 3, 3), dtype=complex)
    c_r = np.zeros((3, 3, 3, 3), dtype=complex)
    for i, y10 in enumerate(b10):
        for j, y5 in enumerate(b5):
            c_l += w[i, j] * d5.tensor_llll(
                d5.sym(y10),
                y5,
                ul["up"],
                ul["up"],
                ul["down"],
                u_nu,
            )
            c_r += w[i, j] * d5.tensor_rrrr(
                y10,
                y5,
                ur["up"],
                ur["charged_lepton"],
                ur["up"],
                ur["down"],
            )
    return c_l, c_r


def lifetime_from_amp(amplitude: float, channel_class: str, chiral: dict[str, Any], m_triplet: float) -> dict[str, float | None]:
    constants = pubd5.constants_for(channel_class, chiral)
    base = d5.lifetime_from_amplitude(amplitude, constants, m_triplet_gev=m_triplet)
    return {
        "central_width_prefactor_GeV5": constants["width_prefactor_GeV5"],
        "tau_years_unfiltered": finite_or_none(base["tau_years"]),
        "C6_dressed_GeV_minus2": finite_or_none(base["C6_dressed_GeV_minus2"]),
    }


def filter_required(tau_years: float | None, bound: float) -> float | None:
    if tau_years is None:
        return None
    return math.sqrt(tau_years / bound)


def row(
    block: str,
    operator: str,
    channel: str,
    channel_class: str,
    amplitude: float,
    idx: tuple[int, int, int, int],
    chiral: dict[str, Any],
    m_triplet: float,
) -> dict[str, Any]:
    life = lifetime_from_amp(amplitude, channel_class, chiral, m_triplet)
    tau = life["tau_years_unfiltered"]
    return {
        "inverse_block": block,
        "operator": operator,
        "channel": channel,
        "channel_class": channel_class,
        "selected_index": "".join(str(i) for i in idx),
        "amplitude": float(amplitude),
        "M_triplet_reference_GeV": float(m_triplet),
        "central_width_prefactor_GeV5": life["central_width_prefactor_GeV5"],
        "tau_years_unfiltered": tau,
        "S_T_required_tau_2p4e34": filter_required(tau, 2.4e34),
        "S_T_required_tau_1e35": filter_required(tau, 1.0e35),
    }


def channel_rows(block: str, c_l: np.ndarray, c_r: np.ndarray, chiral: dict[str, Any], m_triplet: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    amp, idx = d5.max_perm_entry(c_l, (0, 0, 1), None)
    rows.append(row(block, "LLLL", "p_to_Kplus_nubar", "Knu", amp, idx, chiral, m_triplet))
    amp, idx = d5.max_perm_entry(c_l, (0, 0, 1), 1)
    rows.append(row(block, "LLLL", "p_to_K0_muplus", "K0mu", amp, idx, chiral, m_triplet))
    amp, idx = d5.max_perm_entry(c_l, (0, 0, 0), 0)
    rows.append(row(block, "LLLL", "p_to_eplus_pi0", "e_pi", amp, idx, chiral, m_triplet))
    amp, idx = d5.max_perm_entry(c_l, (0, 0, 0), 1)
    rows.append(row(block, "LLLL", "p_to_muplus_pi0", "mu_pi", amp, idx, chiral, m_triplet))
    amp, idx = pubd5.max_rrrr_channel(c_r, d_flavor=1, charged=None)
    rows.append(row(block, "RRRR", "p_to_K_like_any_charged", "K0mu", amp, idx, chiral, m_triplet))
    amp, idx = pubd5.max_rrrr_channel(c_r, d_flavor=0, charged=0)
    rows.append(row(block, "RRRR", "p_to_eplus_pi0", "e_pi", amp, idx, chiral, m_triplet))
    amp, idx = pubd5.max_rrrr_channel(c_r, d_flavor=0, charged=1)
    rows.append(row(block, "RRRR", "p_to_muplus_pi0", "mu_pi", amp, idx, chiral, m_triplet))
    return rows


def svd_card(w: np.ndarray, m_lock: float) -> dict[str, Any]:
    u, s, vh = np.linalg.svd(w, full_matrices=True)
    return {
        "singular_values_inverse_block": [float(x) for x in s],
        "mass_singular_values_GeV": [None if x <= 1.0e-30 else float(m_lock / x) for x in s],
        "left_singular_vectors": matrix_json(u),
        "right_singular_vectors_dagger": matrix_json(vh),
        "condition_number_nonzero": float(max(s) / min(x for x in s if x > 1.0e-30)),
    }


def build() -> tuple[dict[str, Any], dict[str, Any]]:
    card = read_json(CARD)
    trip = read_json(TRIPLET_120)
    crossed = read_json(CROSSED)
    chiral = read_json(CHIRAL)
    data = pubd5.rebuild_source_data(card)
    b10, b5, names = component_bases(data)
    ul, ur, u_nu = pubd5.rotations_for_card(card)

    target = trip["verdict"]["most_safe_physical_row"]
    p = cvec(target["p"])
    q = cvec(target["q"])
    kappa = float(card["crossed_projector"]["kappa"])
    m_lock = float(crossed["threshold_locked_dilation"]["M_lock_GeV"])
    blocks = inverse_blocks(p, q, kappa)

    tensors: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    rows: list[dict[str, Any]] = []
    for block_name in ["rank_one", "finite"]:
        c_l, c_r = contract_tensors(blocks[block_name], b10, b5, ul, ur, u_nu)
        tensors[block_name] = (c_l, c_r)
        rows.extend(channel_rows(block_name, c_l, c_r, chiral, m_lock))

    finite_taus = [r["tau_years_unfiltered"] for r in rows if r["inverse_block"] == "finite" and r["tau_years_unfiltered"] is not None]
    rank_taus = [r["tau_years_unfiltered"] for r in rows if r["inverse_block"] == "rank_one" and r["tau_years_unfiltered"] is not None]
    finite_min_tau = min(finite_taus) if finite_taus else None
    rank_min_tau = min(rank_taus) if rank_taus else None

    eigen_card = {
        "basis_names": names,
        "kappa": kappa,
        "M_lock_GeV": m_lock,
        "p_direction": [cjson(z) for z in p],
        "q_direction": [cjson(z) for z in q],
        "U_left_completion": matrix_json(blocks["U_left"]),
        "V_right_completion": matrix_json(blocks["V_right"]),
        "rank_one_inverse_block": matrix_json(blocks["rank_one"]),
        "finite_inverse_block": matrix_json(blocks["finite"]),
        "rank_one_svd": svd_card(blocks["rank_one"], m_lock),
        "finite_svd": svd_card(blocks["finite"], m_lock),
        "C5L_rank_one": tensor4_json(tensors["rank_one"][0]),
        "C5R_rank_one": tensor4_json(tensors["rank_one"][1]),
        "C5L_finite": tensor4_json(tensors["finite"][0]),
        "C5R_finite": tensor4_json(tensors["finite"][1]),
    }
    summary = {
        "note": "No web lookup used. Source-basis triplet eigenstate card with explicit C5L/C5R proxy tensors.",
        "input_manifest": manifest([CARD, TRIPLET_120, CROSSED, CHIRAL]),
        "operator_definitions": {
            "inverse_block": "C5 tensors use W_{IJ}/M_lock, with W exported here and M_lock used as the triplet reference mass.",
            "C5L": "C_L^{abcd}=sum_{I,J} W_{IJ} B10_I^{ij} B5_J^{kl} U_u^{ia} U_u^{jb} U_d^{kc} U_l^{ld}, with sym(B10_I) on QQ.",
            "C5R": "C_R^{abcd}=sum_{I,J} W_{IJ} B10_I^{ij} B5_J^{kl} U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}.",
        },
        "source": {
            "triplet_direction_label": target["label"],
            "triplet_direction_residual": target["residual"],
            "G_A_antisym_residual": float(np.linalg.norm(data["G_A"] + data["G_A"].T)),
            "G_B_antisym_residual": float(np.linalg.norm(data["G_B"] + data["G_B"].T)),
        },
        "rows": rows,
        "verdict": {
            "source_basis_C5_tensors_exported": True,
            "publication_level_d5_complete": False,
            "rank_one_min_tau_years_unfiltered": rank_min_tau,
            "finite_kappa_min_tau_years_unfiltered": finite_min_tau,
            "finite_kappa_min_margin_1e35_unfiltered": None if finite_min_tau is None else finite_min_tau / 1.0e35,
            "reference_triplet_filter": REFERENCE_TRIPLET_FILTER,
            "finite_kappa_min_tau_years_at_reference_filter": (
                None if finite_min_tau is None else finite_min_tau / (REFERENCE_TRIPLET_FILTER * REFERENCE_TRIPLET_FILTER)
            ),
            "finite_kappa_min_margin_1e35_at_reference_filter": (
                None if finite_min_tau is None else finite_min_tau / (REFERENCE_TRIPLET_FILTER * REFERENCE_TRIPLET_FILTER) / 1.0e35
            ),
            "finite_block_condition_number": eigen_card["finite_svd"]["condition_number_nonzero"],
            "rank_one_has_infinite_orthogonal_masses": True,
            "interpretation": (
                "The crossed triplet inverse block is now represented as an explicit finite source-basis "
                "matrix and contracted into C5L/C5R proxy tensors.  This advances the previous class-level "
                "leakage table to a tensor-valued artifact.  It is not a final proton-decay theorem because "
                "the exported tensors still lack full SUSY dressing matrices and final chiral reductions."
            ),
        },
    }
    return summary, eigen_card


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "inverse_block",
        "operator",
        "channel",
        "channel_class",
        "selected_index",
        "amplitude",
        "M_triplet_reference_GeV",
        "central_width_prefactor_GeV5",
        "tau_years_unfiltered",
        "S_T_required_tau_2p4e34",
        "S_T_required_tau_1e35",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "exists", "size_bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Publication triplet eigenstate card",
        "",
        "No web lookup was used.",
        "",
        "This card diagonalizes the source-basis crossed triplet inverse block and",
        "contracts it into explicit C5L/C5R proxy tensors.  It is still not the",
        "final dressed proton-decay theorem.",
        "",
        "## Eigenstate summary",
        "",
        f"Triplet direction: `{summary['source']['triplet_direction_label']}`.",
        f"Finite block condition number: `{summary['verdict']['finite_block_condition_number']:.6e}`.",
        f"Finite-block minimum unfiltered lifetime: `{summary['verdict']['finite_kappa_min_tau_years_unfiltered']:.6e}` yr.",
        f"Finite-block minimum margin at 1e35 yr: `{summary['verdict']['finite_kappa_min_margin_1e35_unfiltered']:.6e}`.",
        f"Reference-filter minimum margin at 1e35 yr: `{summary['verdict']['finite_kappa_min_margin_1e35_at_reference_filter']:.6e}` for `S_T=1e-5`.",
        "",
        "## Channel rows",
        "",
        "| block | op | channel | index | amp | tau unfiltered | S_T max 1e35 |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in summary["rows"]:
        tau = row["tau_years_unfiltered"]
        tau_s = "inf" if tau is None else f"{tau:.3e}"
        st = row["S_T_required_tau_1e35"]
        st_s = "inf" if st is None else f"{st:.3e}"
        lines.append(
            f"| `{row['inverse_block']}` | `{row['operator']}` | `{row['channel']}` | "
            f"`{row['selected_index']}` | {row['amplitude']:.3e} | {tau_s} | {st_s} |"
        )
    lines += [
        "",
        "## Verdict",
        "",
        summary["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary, eigen_card = build()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "triplet_eigenstate_card.json").write_text(json.dumps(eigen_card, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "channel_table.csv", summary["rows"])
    write_manifest(OUT / "input_manifest.csv", summary["input_manifest"])
    write_report(summary)
    verdict = summary["verdict"]
    print("Publication triplet eigenstate card")
    print(f"  C5 tensors exported: {verdict['source_basis_C5_tensors_exported']}")
    print(f"  publication d5 complete: {verdict['publication_level_d5_complete']}")
    print(f"  finite block condition: {verdict['finite_block_condition_number']:.6e}")
    print(f"  finite min tau: {verdict['finite_kappa_min_tau_years_unfiltered']:.6e} yr")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
