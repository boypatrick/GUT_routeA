#!/usr/bin/env python3
"""Dressed d=5 replay sourced directly from the publication eigenstate C5 card.

No web lookup is used.

This audit is the bridge that the roadmap currently asks for:

    triplet_eigenstate_card.json
      -> explicit C5L/C5R tensors
      -> chargino/neutralino + sfermion-eigenstate dressing kernels
      -> channel-specific lifetime rows.

It is still labelled DRESSED_PROXY, not a final proton-decay theorem, because
the chiral reductions and soft spectrum are the local audited proxies already
used elsewhere in the project.  The point is narrower: no scalar leakage
shortcut is allowed here.  Every channel row must start from the exported
eigenstate-card tensors.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_eigenstate_d5_dressing as eig  # noqa: E402
import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import audit_mssm_mixing_d5_dressing as mssm  # noqa: E402
import audit_publication_channel_d5_tables as pubd5  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


CARD = ROOT / "output" / "publication_triplet_eigenstate_card" / "triplet_eigenstate_card.json"
CARD_SUMMARY = ROOT / "output" / "publication_triplet_eigenstate_card" / "summary.json"
CHIRAL = ROOT / "output" / "chiral_lattice_d5" / "summary.json"
OUT = ROOT / "output" / "publication_dressed_c5_from_eigenstate_card"

DISPLAY_ST = 1.0e-5
PAIRINGS = list(combinations(range(4), 2))

WIDTH_CASES = {
    "central": {
        "Knu": 0.00039531226931756774,
        "K0mu": 0.00006965259084548448,
        "e_pi": 0.0011318398546159686,
        "mu_pi": 0.0011318398546159686,
    },
    "max_width": {
        "Knu": 0.0019122754090593987,
        "K0mu": 0.00040366230812257124,
        "e_pi": 0.005713123603904858,
        "mu_pi": 0.005713123603904858,
    },
    "min_width": {
        "Knu": 0.00007673089919519867,
        "K0mu": 0.00001346200606010078,
        "e_pi": 0.0002088057718555193,
        "mu_pi": 0.0002088057718555193,
    },
}

CHANNELS = [
    {
        "operator": "LLLL",
        "channel": "p_to_Kplus_nubar",
        "channel_class": "Knu",
        "q_flavors": (0, 0, 1),
        "lepton": None,
        "slot_kinds": ["uL", "uL", "dL", "nuL"],
        "present_bound_years": 2.4e34,
    },
    {
        "operator": "LLLL",
        "channel": "p_to_K0_muplus",
        "channel_class": "K0mu",
        "q_flavors": (0, 0, 1),
        "lepton": 1,
        "slot_kinds": ["uL", "uL", "dL", "eL"],
        "present_bound_years": 1.0e34,
    },
    {
        "operator": "LLLL",
        "channel": "p_to_eplus_pi0",
        "channel_class": "e_pi",
        "q_flavors": (0, 0, 0),
        "lepton": 0,
        "slot_kinds": ["uL", "uL", "dL", "eL"],
        "present_bound_years": 2.4e34,
    },
    {
        "operator": "LLLL",
        "channel": "p_to_muplus_pi0",
        "channel_class": "mu_pi",
        "q_flavors": (0, 0, 0),
        "lepton": 1,
        "slot_kinds": ["uL", "uL", "dL", "eL"],
        "present_bound_years": 2.4e34,
    },
    {
        "operator": "RRRR",
        "channel": "p_to_K_like_any_charged",
        "channel_class": "K0mu",
        "d_flavor": 1,
        "charged": None,
        "slot_kinds": ["uR", "uR", "dR", "eR"],
        "present_bound_years": 2.4e34,
    },
    {
        "operator": "RRRR",
        "channel": "p_to_eplus_pi0",
        "channel_class": "e_pi",
        "d_flavor": 0,
        "charged": 0,
        "slot_kinds": ["uR", "uR", "dR", "eR"],
        "present_bound_years": 2.4e34,
    },
    {
        "operator": "RRRR",
        "channel": "p_to_muplus_pi0",
        "channel_class": "mu_pi",
        "d_flavor": 0,
        "charged": 1,
        "slot_kinds": ["uR", "uR", "dR", "eR"],
        "present_bound_years": 2.4e34,
    },
]


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


def tensor4(raw: list[list[list[list[dict[str, float]]]]]) -> np.ndarray:
    return np.array(
        [
            [
                [[cnum(raw[a][b][c][d]) for d in range(3)] for c in range(3)]
                for b in range(3)
            ]
            for a in range(3)
        ],
        dtype=complex,
    )


def width_to_tau(width_gev: float) -> float:
    return mi.width_to_tau(width_gev)


def lifetime(
    amplitude_with_dressing: float,
    st: float,
    m_triplet: float,
    width_prefactor: float,
) -> dict[str, float]:
    c6 = st * amplitude_with_dressing / m_triplet
    width = width_prefactor * c6 * c6
    return {
        "C6_GeV_minus2": float(c6),
        "width_GeV": float(width),
        "tau_years": float(width_to_tau(width)),
    }


def st_max_for(
    amplitude_with_dressing: float,
    target_years: float,
    m_triplet: float,
    width_prefactor: float,
) -> float:
    unit = lifetime(amplitude_with_dressing, 1.0, m_triplet, width_prefactor)
    return math.sqrt(unit["tau_years"] / target_years)


def select_entry(tensor: np.ndarray, channel: dict[str, Any]) -> tuple[float, tuple[int, int, int, int], complex]:
    if channel["operator"] == "LLLL":
        amp, idx = d5.max_perm_entry(tensor, channel["q_flavors"], channel["lepton"])
        return amp, idx, tensor[idx]
    amp, idx = pubd5.max_rrrr_channel(tensor, channel["d_flavor"], channel["charged"])
    return amp, idx, tensor[idx]


def channel_tensor(card: dict[str, Any], block: str, operator: str) -> np.ndarray:
    key = f"C5{'L' if operator == 'LLLL' else 'R'}_{block}"
    return tensor4(card[key])


def dressed_channel_rows_for_point(
    card: dict[str, Any],
    block: str,
    spec: dict[str, Any],
    point: dict[str, float],
    vac: dict[str, float],
    ewkino: dict[str, Any],
    scenario: str,
    epsilon: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    m_triplet = float(card["M_lock_GeV"])
    for channel in CHANNELS:
        base = channel_tensor(card, block, channel["operator"])
        pair_rows: list[dict[str, Any]] = []
        for pair in PAIRINGS:
            kind_i = channel["slot_kinds"][pair[0]]
            kind_j = channel["slot_kinds"][pair[1]]
            kernel, parts = mssm.pair_kernel(
                spec["by_kind"][kind_i],
                spec["by_kind"][kind_j],
                point,
                vac,
                ewkino,
                kind_i,
                kind_j,
                channel["operator"],
            )
            dressed_tensor = eig.apply_pair_kernel(base, kernel, pair)
            amp, idx, value = select_entry(dressed_tensor, channel)
            for width_case, prefactors in WIDTH_CASES.items():
                width_prefactor = prefactors[channel["channel_class"]]
                life_display = lifetime(amp, DISPLAY_ST, m_triplet, width_prefactor)
                st_present = st_max_for(amp, channel["present_bound_years"], m_triplet, width_prefactor)
                st_1e35 = st_max_for(amp, 1.0e35, m_triplet, width_prefactor)
                pair_rows.append(
                    {
                        "status": "DRESSED_PROXY",
                        "inverse_block": block,
                        "scenario": scenario,
                        "epsilon": float(epsilon),
                        "spectrum_name": point["name"],
                        "operator": channel["operator"],
                        "channel": channel["channel"],
                        "channel_class": channel["channel_class"],
                        "normalization_case": width_case,
                        "pair": f"{pair[0]}{pair[1]}",
                        "selected_index": "".join(str(i) for i in idx),
                        "amplitude_with_dressing": float(amp),
                        "selected_value_re": float(np.real(value)),
                        "selected_value_im": float(np.imag(value)),
                        "M_triplet_GeV": m_triplet,
                        "S_T_display": DISPLAY_ST,
                        "width_prefactor_GeV5": width_prefactor,
                        "C6_GeV_minus2_at_ST_display": life_display["C6_GeV_minus2"],
                        "tau_years_at_ST_display": life_display["tau_years"],
                        "present_bound_years": channel["present_bound_years"],
                        "margin_present_at_ST_display": life_display["tau_years"] / channel["present_bound_years"],
                        "margin_1e35_at_ST_display": life_display["tau_years"] / 1.0e35,
                        "S_T_max_present": st_present,
                        "S_T_max_1e35": st_1e35,
                        "passes_present_at_ST_display": life_display["tau_years"] >= channel["present_bound_years"],
                        "passes_1e35_at_ST_display": life_display["tau_years"] >= 1.0e35,
                        "avg_chargino_part": parts["chargino"],
                        "avg_neutralino_part": parts["neutralino"],
                        "avg_higgsino_part": parts["higgsino"],
                        "avg_total_part": parts["total"],
                        "note": "All rows are dressed from eigenstate-card C5 tensors; no scalar leakage shortcut is used.",
                    }
                )
        rows.extend(pair_rows)
    return rows


def trim(row: dict[str, Any] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    keys = [
        "status",
        "inverse_block",
        "scenario",
        "epsilon",
        "spectrum_name",
        "operator",
        "channel",
        "normalization_case",
        "pair",
        "selected_index",
        "amplitude_with_dressing",
        "S_T_max_1e35",
        "tau_years_at_ST_display",
        "margin_1e35_at_ST_display",
        "passes_1e35_at_ST_display",
    ]
    return {key: row[key] for key in keys}


def summarize(rows: list[dict[str, Any]], point_rows: list[dict[str, Any]], card_summary: dict[str, Any]) -> dict[str, Any]:
    by_block_case: dict[str, Any] = {}
    for block in ["rank_one", "finite"]:
        block_rows = [row for row in rows if row["inverse_block"] == block]
        by_case: dict[str, Any] = {}
        for case in WIDTH_CASES:
            case_rows = [row for row in block_rows if row["normalization_case"] == case]
            unsafe_present = [row for row in case_rows if not row["passes_present_at_ST_display"]]
            unsafe_1e35 = [row for row in case_rows if not row["passes_1e35_at_ST_display"]]
            by_case[case] = {
                "channel_rows": len(case_rows),
                "unsafe_present_rows": len(unsafe_present),
                "unsafe_1e35_rows": len(unsafe_1e35),
                "global_worst_present": trim(min(case_rows, key=lambda row: row["margin_present_at_ST_display"])) if case_rows else None,
                "global_worst_1e35": trim(min(case_rows, key=lambda row: row["margin_1e35_at_ST_display"])) if case_rows else None,
                "global_S_T_max_present": min((row["S_T_max_present"] for row in case_rows), default=math.inf),
                "global_S_T_max_1e35": min((row["S_T_max_1e35"] for row in case_rows), default=math.inf),
            }
        by_block_case[block] = by_case

    baseline_rows = [
        row for row in point_rows
        if row["inverse_block"] == "finite"
        and row["scenario"] == "zero"
        and row["epsilon"] == 0.0
        and row["spectrum_name"] == "baseline_100TeV"
    ]
    preferred = by_block_case["finite"]["central"]
    max_width = by_block_case["finite"]["max_width"]
    complete = (
        preferred["unsafe_present_rows"] == 0
        and preferred["unsafe_1e35_rows"] == 0
        and max_width["unsafe_present_rows"] == 0
        and max_width["unsafe_1e35_rows"] == 0
    )
    return {
        "note": "No web lookup used. Dressed d=5 proxy sourced only from publication triplet eigenstate C5 tensors.",
        "status": "DRESSED_PROXY",
        "input_manifest": manifest([CARD, CARD_SUMMARY, CHIRAL]),
        "display_triplet_filter": DISPLAY_ST,
        "source_card_verdict": card_summary["verdict"],
        "audited_blocks": ["rank_one", "finite"],
        "normalization_cases": list(WIDTH_CASES),
        "soft_scenarios": sorted({row["scenario"] for row in point_rows}),
        "spectrum_points": mi.SPECTRUM_POINTS,
        "by_block_case": by_block_case,
        "baseline_finite_point_rows": baseline_rows,
        "verdict": {
            "publication_level_d5_complete": False,
            "dressed_proxy_complete_from_eigenstate_card": True,
            "all_finite_central_rows_pass_1e35_at_ST_1e_minus_5": preferred["unsafe_1e35_rows"] == 0,
            "all_finite_max_width_rows_pass_1e35_at_ST_1e_minus_5": max_width["unsafe_1e35_rows"] == 0,
            "can_promote_to_PASS_CONDITIONAL": complete,
            "interpretation": (
                "The audit removes the scalar leakage shortcut: every channel row is generated by applying "
                "chargino/neutralino and sfermion-eigenstate dressing kernels to the exported C5L/C5R tensors. "
                "The result remains a DRESSED_PROXY because the soft spectrum and chiral reductions are still "
                "the project's local proxy layer rather than a final external phenomenology package."
            ),
        },
    }


def point_summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        key = (
            row["inverse_block"],
            row["scenario"],
            row["epsilon"],
            row["spectrum_name"],
            row["normalization_case"],
        )
        grouped.setdefault(key, []).append(row)
    out: list[dict[str, Any]] = []
    for key, items in sorted(grouped.items()):
        worst_present = min(items, key=lambda row: row["margin_present_at_ST_display"])
        worst_1e35 = min(items, key=lambda row: row["margin_1e35_at_ST_display"])
        out.append(
            {
                "inverse_block": key[0],
                "scenario": key[1],
                "epsilon": float(key[2]),
                "spectrum_name": key[3],
                "normalization_case": key[4],
                "channel_rows": len(items),
                "all_present_pass": all(row["passes_present_at_ST_display"] for row in items),
                "all_1e35_pass": all(row["passes_1e35_at_ST_display"] for row in items),
                "worst_present": trim(worst_present),
                "worst_1e35": trim(worst_1e35),
            }
        )
    return out


def audit() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    card = read_json(CARD)
    card_summary = read_json(CARD_SUMMARY)
    vac = mi.vacuum_inputs()
    soft = mi.insertion_basis()
    rows: list[dict[str, Any]] = []

    for block in ["rank_one", "finite"]:
        for scenario, base_delta in soft["matrices"].items():
            for eps in eig.SOFT_EPS_GRID:
                spec = eig.spectrum_for_scenario(base_delta, eps)
                if not spec["positive_definite"]:
                    continue
                for point in mi.SPECTRUM_POINTS:
                    ewkino = mssm.electroweakino_spectrum(point)
                    rows.extend(
                        dressed_channel_rows_for_point(
                            card,
                            block,
                            spec,
                            point,
                            vac,
                            ewkino,
                            scenario,
                            eps,
                        )
                    )
    points = point_summaries(rows)
    summary = summarize(rows, points, card_summary)
    return rows, points, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Publication dressed C5 replay from eigenstate card",
        "",
        "No web lookup was used.",
        "",
        "This replay applies the existing chargino/neutralino and positive",
        "sfermion-eigenstate dressing kernels directly to the exported",
        "`C5L`/`C5R` tensors from the publication triplet eigenstate card.",
        "It is marked `DRESSED_PROXY`, not a final proton-decay theorem.",
        "",
        "## Global gates",
        "",
    ]
    for block in ["rank_one", "finite"]:
        lines.append(f"### `{block}` inverse block")
        lines.append("")
        lines.append("| width case | rows | unsafe present | unsafe 1e35 | global S_T max 1e35 | worst 1e35 row |")
        lines.append("|---|---:|---:|---:|---:|---|")
        for case, item in summary["by_block_case"][block].items():
            worst = item["global_worst_1e35"]
            worst_s = (
                "none" if worst is None else
                f"{worst['channel']} {worst['operator']} {worst['scenario']} eps={worst['epsilon']} "
                f"{worst['spectrum_name']} pair={worst['pair']} margin={worst['margin_1e35_at_ST_display']:.3e}"
            )
            lines.append(
                f"| `{case}` | {item['channel_rows']} | {item['unsafe_present_rows']} | "
                f"{item['unsafe_1e35_rows']} | {item['global_S_T_max_1e35']:.6e} | {worst_s} |"
            )
        lines.append("")
    lines += [
        "## Verdict",
        "",
        summary["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, points, summary = audit()
    write_csv(OUT / "dressed_channel_rows.csv", rows)
    write_csv(OUT / "dressed_point_rows.csv", points)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    print("Publication dressed C5 replay from eigenstate card")
    print(f"  status: {summary['status']}")
    print(f"  finite central unsafe 1e35 rows: {summary['by_block_case']['finite']['central']['unsafe_1e35_rows']}")
    print(f"  finite max-width unsafe 1e35 rows: {summary['by_block_case']['finite']['max_width']['unsafe_1e35_rows']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
