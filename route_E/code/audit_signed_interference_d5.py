#!/usr/bin/env python3
"""Signed-interference dimension-five proton-decay proxy.

No web lookup is used.  This audit keeps the complex Wilson tensors and scans
relative CP phases between the triplet filter and the chargino/neutralino/
higgsino dressing components.  It is the next controlled step after the
coherent-absolute MSSM-mixing proxy:

* the soft spectra are diagonalized as before;
* chargino and neutralino mass matrices are resolved as before;
* the component Wilson tensors are kept complex;
* a small phase grid is used to ask whether the preferred branch needs
  accidental destructive interference.

The output is a robustness ledger, not a final proton-decay prediction.
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

import audit_eigenstate_d5_dressing as eig  # noqa: E402
import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import audit_mssm_mixing_d5_dressing as mix  # noqa: E402
import construct_triplet_rank_lift as rank_lift  # noqa: E402


OUT = ROOT / "output" / "signed_interference_d5"

DISPLAY_ST = 1.0e-5
PHASE_GRID = [0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi]
TRIPLET_PHASE_GRID = [0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0]
BASIS_PHASE_CHARGES = np.array([0.0, 0.0, 1.0, 2.0])


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def phase_deform_w(w: np.ndarray, phi: float) -> np.ndarray:
    phase = np.exp(1j * phi * BASIS_PHASE_CHARGES)
    return np.diag(phase) @ w @ np.diag(phase.conjugate())


def component_pair_kernels(
    spec_i: dict[str, Any],
    spec_j: dict[str, Any],
    point: dict[str, float],
    vac: dict[str, float],
    ewkino: dict[str, Any],
    kind_i: str,
    kind_j: str,
    operator: str,
) -> dict[str, np.ndarray]:
    ui = spec_i["eigvecs"]
    uj = spec_j["eigvecs"]
    ei = spec_i["eigvals"]
    ej = spec_j["eigvals"]
    kernels = {
        "chargino": np.zeros((3, 3, 3, 3), dtype=complex),
        "neutralino": np.zeros((3, 3, 3, 3), dtype=complex),
        "higgsino": np.zeros((3, 3, 3, 3), dtype=complex),
    }
    for r in range(3):
        for s in range(3):
            parts = mix.pair_dressing(
                point,
                vac,
                ewkino,
                float(ei[r]),
                float(ej[s]),
                kind_i,
                kind_j,
                operator,
            )
            projector = np.einsum(
                "a,A,b,B->abAB",
                ui[:, r],
                ui[:, r].conjugate(),
                uj[:, s],
                uj[:, s].conjugate(),
                optimize=True,
            )
            for key in kernels:
                kernels[key] += parts[key] * projector
    return kernels


def apply_components(
    tensor: np.ndarray,
    kernels: dict[str, np.ndarray],
    pair: tuple[int, int],
) -> dict[str, np.ndarray]:
    return {
        key: eig.apply_pair_kernel(tensor, kernel, pair)
        for key, kernel in kernels.items()
    }


def component_signs(operator: str) -> dict[str, float]:
    if operator == "LLLL":
        return {"chargino": 1.0, "neutralino": -1.0, "higgsino": 0.0}
    return {"chargino": 0.0, "neutralino": 1.0, "higgsino": -1.0}


def channel_value(tensor: np.ndarray, entries: list[tuple[int, int, int, int]]) -> tuple[float, str, str, complex]:
    values = [complex(tensor[idx]) for idx in entries]
    max_pos = max(range(len(entries)), key=lambda i: abs(values[i]))
    max_entry_amp = abs(values[max_pos])
    normalized_sum = sum(values) / math.sqrt(max(len(values), 1))
    if abs(normalized_sum) > max_entry_amp:
        return float(abs(normalized_sum)), "normalized_chiral_sum", "sum", normalized_sum
    return float(max_entry_amp), "entry_max", "".join(str(x) for x in entries[max_pos]), values[max_pos]


def expi(phi: float) -> complex:
    return complex(math.cos(phi), math.sin(phi))


def phase_scan(
    component_tensors: dict[str, np.ndarray],
    entries: list[tuple[int, int, int, int]],
    operator: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    signs = component_signs(operator)
    min_row: dict[str, Any] | None = None
    max_row: dict[str, Any] | None = None
    for phi_n in PHASE_GRID:
        for phi_h in PHASE_GRID:
            tensor = (
                signs["chargino"] * component_tensors["chargino"]
                + signs["neutralino"] * expi(phi_n) * component_tensors["neutralino"]
                + signs["higgsino"] * expi(phi_h) * component_tensors["higgsino"]
            )
            amp, reduction, selected, value = channel_value(tensor, entries)
            row = {
                "amplitude": amp,
                "phi_neutralino": phi_n,
                "phi_higgsino": phi_h,
                "reduction": reduction,
                "selected": selected,
                "value": mi.cjson(value),
            }
            if min_row is None or amp < min_row["amplitude"]:
                min_row = row
            if max_row is None or amp > max_row["amplitude"]:
                max_row = row
    if min_row is None or max_row is None:
        raise RuntimeError("empty phase scan")
    return min_row, max_row


def custom_channel_specs(
    combined: dict[str, tuple[np.ndarray, list[tuple[int, int, int, int]]]]
) -> dict[str, dict[str, Any]]:
    specs: dict[str, dict[str, Any]] = {}
    for channel, (tensor, entries) in combined.items():
        specs[channel] = {
            "tensor": tensor,
            "entries": entries,
            **mi.CHANNEL_META[channel],
        }
    # Minimal e pi proxies using existing charged-lepton and RRRR tensors.
    specs["LLLL_uude_epi"] = {
        "tensor": combined["LLLL_upupdown_K0mu"][0],
        "entries": [(0, 0, 0, 0)],
        "operator": "LLLL",
        "target_years": 2.4e34,
        "slot_kinds": ["uL", "uL", "dL", "eL"],
    }
    specs["RRRR_uude_epi"] = {
        "tensor": combined["RRRR_uusd_anycharged"][0],
        "entries": [(0, 0, 0, 0)],
        "operator": "RRRR",
        "target_years": 2.4e34,
        "slot_kinds": ["uR", "uR", "dR", "eR"],
    }
    return specs


def c6_lifetime(amplitude: float, st: float, m_triplet: float, width_prefactor: float) -> dict[str, float]:
    return eig.c6_lifetime(amplitude, st, m_triplet, width_prefactor)


def st_max_for(amplitude: float, target_years: float, m_triplet: float, width_prefactor: float) -> float:
    return eig.st_max_for(amplitude, target_years, m_triplet, width_prefactor)


def audit() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    _basis_names, _pair_full, channel_entries, _constants_unused = rank_lift.build_channel_tensors()
    proton = read_json(mi.PROTON)["hadronic_constants"]
    vac = mi.vacuum_inputs()
    soft = mi.insertion_basis()
    mssm_summary = read_json(mix.OUT / "summary.json")

    rows: list[dict[str, Any]] = []
    channel_rows: list[dict[str, Any]] = []

    for filt in mi.selected_filter_rows():
        w0 = mi.filter_matrix(filt)
        for scenario, base_delta in soft["matrices"].items():
            for eps in eig.SOFT_EPS_GRID:
                spec = eig.spectrum_for_scenario(base_delta, eps)
                for point in mi.SPECTRUM_POINTS:
                    ewkino = mix.electroweakino_spectrum(point)
                    row_channel_results: list[dict[str, Any]] = []
                    for phi_t in TRIPLET_PHASE_GRID:
                        w = phase_deform_w(w0, phi_t)
                        combined = mi.combine_tensors(w, channel_entries)
                        specs = custom_channel_specs(combined)
                        for channel, meta in specs.items():
                            if not spec["positive_definite"]:
                                continue
                            pair_results: list[dict[str, Any]] = []
                            for pair in eig.PAIRINGS:
                                kind_i = meta["slot_kinds"][pair[0]]
                                kind_j = meta["slot_kinds"][pair[1]]
                                kernels = component_pair_kernels(
                                    spec["by_kind"][kind_i],
                                    spec["by_kind"][kind_j],
                                    point,
                                    vac,
                                    ewkino,
                                    kind_i,
                                    kind_j,
                                    meta["operator"],
                                )
                                tensors = apply_components(meta["tensor"], kernels, pair)
                                min_phase, max_phase = phase_scan(tensors, meta["entries"], meta["operator"])
                                life = c6_lifetime(
                                    max_phase["amplitude"],
                                    DISPLAY_ST,
                                    vac["M_T_GeV"],
                                    proton["width_prefactor_GeV5"],
                                )
                                st_allowed = st_max_for(
                                    max_phase["amplitude"],
                                    meta["target_years"],
                                    vac["M_T_GeV"],
                                    proton["width_prefactor_GeV5"],
                                )
                                pair_results.append(
                                    {
                                        "filter_label": filt["label"],
                                        "scenario": scenario,
                                        "epsilon": eps,
                                        "spectrum_name": point["name"],
                                        "triplet_phase": phi_t,
                                        "channel": channel,
                                        "operator": meta["operator"],
                                        "pair": f"{pair[0]}{pair[1]}",
                                        "min_phase_amplitude": min_phase["amplitude"],
                                        "max_phase_amplitude": max_phase["amplitude"],
                                        "phase_spread": (
                                            max_phase["amplitude"] / max(min_phase["amplitude"], 1.0e-300)
                                        ),
                                        "phi_neutralino_at_max": max_phase["phi_neutralino"],
                                        "phi_higgsino_at_max": max_phase["phi_higgsino"],
                                        "reduction_at_max": max_phase["reduction"],
                                        "selected_at_max": max_phase["selected"],
                                        "S_T_max": float(st_allowed),
                                        "tau_years_ST_display": life["tau_years"],
                                        "margin_at_ST_display": life["tau_years"] / meta["target_years"],
                                        "passes": life["tau_years"] >= meta["target_years"],
                                    }
                                )
                            worst_pair = min(pair_results, key=lambda item: item["margin_at_ST_display"])
                            row_channel_results.append(worst_pair)
                            channel_rows.append(worst_pair)
                    if not spec["positive_definite"]:
                        rows.append(
                            {
                                "filter_label": filt["label"],
                                "scenario": scenario,
                                "epsilon": eps,
                                "spectrum_name": point["name"],
                                **point,
                                "positive_definite": False,
                                "min_soft_eigenvalue": spec["min_eigenvalue"],
                                "worst_channel": "",
                                "worst_margin": 0.0,
                                "worst_ST_max": 0.0,
                                "max_phase_spread": math.inf,
                                "all_channels_pass": False,
                            }
                        )
                        continue
                    worst = min(row_channel_results, key=lambda item: item["margin_at_ST_display"])
                    rows.append(
                        {
                            "filter_label": filt["label"],
                            "scenario": scenario,
                            "epsilon": eps,
                            "spectrum_name": point["name"],
                            **point,
                            "positive_definite": True,
                            "min_soft_eigenvalue": spec["min_eigenvalue"],
                            "worst_channel": worst["channel"],
                            "worst_triplet_phase": worst["triplet_phase"],
                            "worst_pair": worst["pair"],
                            "worst_margin": worst["margin_at_ST_display"],
                            "worst_ST_max": worst["S_T_max"],
                            "max_phase_spread": max(item["phase_spread"] for item in row_channel_results),
                            "all_channels_pass": all(item["passes"] for item in row_channel_results),
                        }
                    )

    summary = summarize(rows, channel_rows, mssm_summary)
    return rows, channel_rows, summary


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
        "positive_definite",
        "min_soft_eigenvalue",
        "worst_channel",
        "worst_triplet_phase",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_phase_spread",
        "all_channels_pass",
    ]
    return {key: row[key] for key in keys}


def scenario_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for scenario in sorted({row["scenario"] for row in rows}):
        subset = [row for row in rows if row["scenario"] == scenario]
        positive = [row for row in subset if row["positive_definite"]]
        safe = [row for row in positive if row["all_channels_pass"]]
        table.append(
            {
                "scenario": scenario,
                "safe_points": len(safe),
                "positive_points": len(positive),
                "total_points": len(subset),
                "safe_fraction_among_positive": len(safe) / len(positive) if positive else 0.0,
                "worst_margin": min((row["worst_margin"] for row in positive), default=0.0),
                "max_phase_spread": max((row["max_phase_spread"] for row in positive), default=math.nan),
            }
        )
    return table


def summarize(
    rows: list[dict[str, Any]],
    channel_rows: list[dict[str, Any]],
    mssm_summary: dict[str, Any],
) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    aligned_names = {
        "zero",
        "up_aligned_LL_MFV",
        "down_aligned_LL_MFV",
        "right_third_split",
        "commutator_LL",
        "combined_LL_RR",
    }
    for label in mi.FILTER_LABELS:
        subset = [row for row in rows if row["filter_label"] == label]
        positive = [row for row in subset if row["positive_definite"]]
        safe = [row for row in positive if row["all_channels_pass"]]
        unsafe = [row for row in positive if not row["all_channels_pass"]]
        aligned = [row for row in positive if row["scenario"] in aligned_names]
        aligned_safe = [row for row in aligned if row["all_channels_pass"]]
        by_filter[label] = {
            "total_points": len(subset),
            "positive_points": len(positive),
            "safe_points": len(safe),
            "unsafe_points": len(unsafe),
            "safe_fraction_among_positive": len(safe) / len(positive) if positive else 0.0,
            "aligned_positive_points": len(aligned),
            "aligned_safe_points": len(aligned_safe),
            "aligned_safe_fraction": len(aligned_safe) / len(aligned) if aligned else 0.0,
            "most_marginal_safe": trim_row(min(safe, key=lambda row: row["worst_margin"])) if safe else None,
            "nearest_unsafe": trim_row(max(unsafe, key=lambda row: row["worst_margin"])) if unsafe else None,
            "max_phase_spread": max((row["max_phase_spread"] for row in positive), default=math.nan),
            "scenario_table": scenario_table(subset),
        }

    preferred = by_filter["omegaR_0.1_kappa_100"]
    verdict = (
        "The preferred filter remains safe under the signed phase-grid replay, "
        "so the branch does not require accidental destructive interference "
        "within the audited CP-phase grid.  A final calculation still needs "
        "the real chiral and lattice matrix-element normalization."
    )
    if preferred["nearest_unsafe"] is not None:
        verdict = (
            "The signed phase-grid replay finds unsafe positive points in the "
            "preferred filter.  Proton safety would then require a phase or "
            "soft-alignment condition beyond the current triplet filter."
        )
    return {
        "note": "No web lookup used. Signed complex d=5 interference proxy with CP-phase scan.",
        "display_triplet_filter": DISPLAY_ST,
        "phase_grid": PHASE_GRID,
        "triplet_phase_grid": TRIPLET_PHASE_GRID,
        "basis_phase_charges": BASIS_PHASE_CHARGES.tolist(),
        "channels_include_minimal_epi_proxy": True,
        "reference_mssm_mixing": {
            "preferred_safe": mssm_summary["by_filter_label"]["omegaR_0.1_kappa_100"]["safe_points"],
            "preferred_positive": mssm_summary["by_filter_label"]["omegaR_0.1_kappa_100"]["positive_points"],
        },
        "by_filter_label": by_filter,
        "verdict": {
            "preferred_filter": "omegaR_0.1_kappa_100",
            "preferred_aligned_safe_fraction": preferred["aligned_safe_fraction"],
            "preferred_safe_fraction_among_positive": preferred["safe_fraction_among_positive"],
            "preferred_has_unsafe_positive_point": preferred["nearest_unsafe"] is not None,
            "interpretation": verdict,
        },
    }


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
        "positive_definite",
        "min_soft_eigenvalue",
        "worst_channel",
        "worst_triplet_phase",
        "worst_pair",
        "worst_margin",
        "worst_ST_max",
        "max_phase_spread",
        "all_channels_pass",
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
        "triplet_phase",
        "channel",
        "operator",
        "pair",
        "min_phase_amplitude",
        "max_phase_amplitude",
        "phase_spread",
        "phi_neutralino_at_max",
        "phi_higgsino_at_max",
        "reduction_at_max",
        "selected_at_max",
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
        "# Signed-interference d=5 proxy",
        "",
        "No web lookup was used.",
        "",
        "Complex chargino, neutralino and higgsino component tensors are combined",
        "over a small CP phase grid.  The reported margin is the worst phase point.",
        "",
        "## Filter summaries",
        "",
        "| filter | aligned safe | positive safe | max phase spread | marginal safe | nearest unsafe |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, payload in summary["by_filter_label"].items():
        marginal = payload["most_marginal_safe"]
        unsafe = payload["nearest_unsafe"]
        lines.append(
            f"| `{label}` | {payload['aligned_safe_points']}/{payload['aligned_positive_points']} | "
            f"{payload['safe_points']}/{payload['positive_points']} | {payload['max_phase_spread']:.3e} | "
            f"{marginal['worst_margin'] if marginal else math.inf:.3e} | "
            f"{unsafe['worst_margin'] if unsafe else math.inf:.3e} |"
        )
    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    lines.extend(
        [
            "",
            "## Preferred scenario table",
            "",
            "| scenario | safe/positive | worst margin | max phase spread |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in preferred["scenario_table"]:
        lines.append(
            f"| `{row['scenario']}` | {row['safe_points']}/{row['positive_points']} | "
            f"{row['worst_margin']:.3e} | {row['max_phase_spread']:.3e} |"
        )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, channel_rows, summary = audit()
    write_csv(OUT / "signed_interference_scan.csv", rows)
    write_channel_csv(OUT / "channel_scan.csv", channel_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    preferred = summary["by_filter_label"]["omegaR_0.1_kappa_100"]
    print("Signed-interference d=5 proxy")
    print(
        "  omegaR_0.1_kappa_100: "
        f"aligned_safe={preferred['aligned_safe_points']}/{preferred['aligned_positive_points']}, "
        f"positive_safe={preferred['safe_points']}/{preferred['positive_points']}, "
        f"max_phase_spread={preferred['max_phase_spread']:.3e}"
    )
    if preferred["nearest_unsafe"]:
        print(
            "  nearest unsafe: "
            f"margin={preferred['nearest_unsafe']['worst_margin']:.3e}, "
            f"scenario={preferred['nearest_unsafe']['scenario']}, "
            f"eps={preferred['nearest_unsafe']['epsilon']}"
        )
    else:
        print("  no unsafe positive-definite audited point")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
