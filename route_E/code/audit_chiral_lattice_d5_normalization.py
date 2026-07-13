#!/usr/bin/env python3
"""Chiral/lattice normalization replay for signed d=5 proton proxy.

No web lookup is used.  The previous d=5 audits used a single local hadronic
prefactor.  This script factors that prefactor into

    K = m_p/(64 pi f_pi^2) beta_H^2 A_R^2 C_chiral^2 Phi_2,

and replays the signed-interference worst rows over a conservative local
uncertainty box for beta_H, A_R, D and F.  The goal is not to replace a
lattice calculation; it is to quantify whether the displayed S_T=1e-5 filter
survives once channel-specific chiral factors are not collapsed into one
number.
"""

from __future__ import annotations

import csv
import json
import math
from itertools import product
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SIGNED = ROOT / "output" / "signed_interference_d5" / "channel_scan.csv"
SIGNED_SUMMARY = ROOT / "output" / "signed_interference_d5" / "summary.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "chiral_lattice_d5"

DISPLAY_ST = 1.0e-5
HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0

MASS = {
    "p": 0.938272,
    "pi": 0.134977,
    "K": 0.493677,
    "mu": 0.105658,
}

PARAM_BOX = {
    "beta_H_GeV3": [0.008, 0.012, 0.018],
    "A_R": [1.8, 2.4, 3.2],
    "D": [0.60, 0.80, 1.00],
    "F": [0.35, 0.47, 0.55],
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def kallen(a: float, b: float, c: float) -> float:
    return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c)


def phase_space(m_meson: float, m_lepton: float = 0.0) -> float:
    mp = MASS["p"]
    lam = max(kallen(mp * mp, m_meson * m_meson, m_lepton * m_lepton), 0.0)
    # Normalize to m_p^4 so massless pion channels are close to
    # (1-m_pi^2/m_p^2)^2.
    return math.sqrt(lam) / (mp * mp) * max(1.0 - (m_meson + m_lepton) ** 2 / (mp * mp), 0.0)


def chiral_factor(channel_class: str, d: float, f: float) -> float:
    if channel_class == "e_pi":
        return abs(1.0 + d + f)
    if channel_class == "Knu":
        direct = 1.0 + (d + 3.0 * f) / 3.0
        indirect = 2.0 * d / 3.0
        return math.sqrt(direct * direct + indirect * indirect)
    if channel_class == "K0mu":
        return math.sqrt((1.0 - d + f) ** 2 + (2.0 * d / 3.0) ** 2)
    raise ValueError(f"unknown channel class {channel_class}")


def classify_channel(channel: str) -> str:
    if "epi" in channel:
        return "e_pi"
    if "K0mu" in channel:
        return "K0mu"
    if "Knu" in channel:
        return "Knu"
    if "RRRR_uusd_anycharged" in channel:
        # Conservative interpretation of the generic charged RRRR monitor.
        return "Knu"
    raise ValueError(f"cannot classify channel {channel}")


def channel_target(channel: str) -> float:
    if "K0mu" in channel:
        return 1.0e34
    return 2.4e34


def channel_phase(channel_class: str) -> float:
    if channel_class == "e_pi":
        return phase_space(MASS["pi"], 0.0)
    if channel_class == "Knu":
        return phase_space(MASS["K"], 0.0)
    if channel_class == "K0mu":
        return phase_space(MASS["K"], MASS["mu"])
    raise ValueError(channel_class)


def width_prefactor(params: dict[str, float], channel_class: str) -> float:
    mp = MASS["p"]
    f_pi = 0.130
    c = chiral_factor(channel_class, params["D"], params["F"])
    return (
        mp
        / (64.0 * math.pi * f_pi * f_pi)
        * params["beta_H_GeV3"] ** 2
        * params["A_R"] ** 2
        * c * c
        * channel_phase(channel_class)
    )


def param_grid() -> list[dict[str, float]]:
    keys = list(PARAM_BOX)
    return [dict(zip(keys, values)) for values in product(*(PARAM_BOX[key] for key in keys))]


def read_channel_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with SIGNED.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, Any] = dict(raw)
            for key in [
                "epsilon",
                "triplet_phase",
                "min_phase_amplitude",
                "max_phase_amplitude",
                "phase_spread",
                "phi_neutralino_at_max",
                "phi_higgsino_at_max",
                "S_T_max",
                "tau_years_ST_display",
                "margin_at_ST_display",
            ]:
                row[key] = float(row[key])
            row["passes"] = row["passes"] == "True"
            rows.append(row)
    return rows


def replay_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    proton = read_json(PROTON)["hadronic_constants"]
    old_k = float(proton["width_prefactor_GeV5"])
    signed_summary = read_json(SIGNED_SUMMARY)
    rows = read_channel_rows()
    params = param_grid()
    central = {
        "beta_H_GeV3": float(proton["beta_H_GeV3"]),
        "A_R": float(proton["A_R"]),
        "D": float(proton["D"]),
        "F": float(proton["F"]),
    }

    replay: list[dict[str, Any]] = []
    for row in rows:
        cls = classify_channel(row["channel"])
        k_values = [width_prefactor(p, cls) for p in params]
        k_min = min(k_values)
        k_max = max(k_values)
        k_cent = width_prefactor(central, cls)
        target = channel_target(row["channel"])
        tau_old = row["tau_years_ST_display"]
        for label, k in [("central", k_cent), ("min_width", k_min), ("max_width", k_max)]:
            tau = tau_old * old_k / k
            margin = tau / target
            st_max = row["S_T_max"] * math.sqrt(old_k / k)
            replay.append(
                {
                    **row,
                    "channel_class": cls,
                    "target_years": target,
                    "normalization_case": label,
                    "old_width_prefactor": old_k,
                    "new_width_prefactor": k,
                    "prefactor_ratio_new_over_old": k / old_k,
                    "S_T_max_replayed": st_max,
                    "tau_years_replayed": tau,
                    "margin_replayed": margin,
                    "passes_replayed": margin >= 1.0,
                }
            )

    summary = summarize(replay, signed_summary, central, old_k)
    return replay, summary


def summarize(replay: list[dict[str, Any]], signed_summary: dict[str, Any], central: dict[str, float], old_k: float) -> dict[str, Any]:
    by_filter: dict[str, Any] = {}
    for filt in sorted({row["filter_label"] for row in replay}):
        subset = [row for row in replay if row["filter_label"] == filt]
        by_case: dict[str, Any] = {}
        for case in ["central", "min_width", "max_width"]:
            case_rows = [row for row in subset if row["normalization_case"] == case]
            safe = [row for row in case_rows if row["passes_replayed"]]
            unsafe = [row for row in case_rows if not row["passes_replayed"]]
            by_case[case] = {
                "channel_rows": len(case_rows),
                "safe_channel_rows": len(safe),
                "unsafe_channel_rows": len(unsafe),
                "worst_safe": trim(min(safe, key=lambda row: row["margin_replayed"])) if safe else None,
                "nearest_unsafe": trim(max(unsafe, key=lambda row: row["margin_replayed"])) if unsafe else None,
                "global_S_T_max": min(row["S_T_max_replayed"] for row in case_rows),
                "global_worst_margin": min(row["margin_replayed"] for row in case_rows),
            }
        by_filter[filt] = by_case

    pref = by_filter["omegaR_0.1_kappa_100"]
    verdict = (
        "The central channel-specific chiral replay preserves the preferred "
        "S_T=1e-5 branch, but the conservative max-width envelope does not.  "
        "The filter must be tightened to the max-width global S_T allowance "
        "before the d=5 claim can be called robust."
    )
    if pref["max_width"]["global_worst_margin"] >= 1.0:
        verdict = (
            "The preferred S_T=1e-5 branch survives even the conservative "
            "channel-specific chiral width envelope."
        )
    return {
        "note": "No web lookup used. Chiral/lattice normalization replay for signed d=5 proxy.",
        "old_single_width_prefactor": old_k,
        "central_params": central,
        "parameter_box": PARAM_BOX,
        "phase_space": {
            "e_pi": channel_phase("e_pi"),
            "Knu": channel_phase("Knu"),
            "K0mu": channel_phase("K0mu"),
        },
        "central_width_prefactors": {
            cls: width_prefactor(central, cls)
            for cls in ["e_pi", "Knu", "K0mu"]
        },
        "signed_reference": signed_summary["by_filter_label"]["omegaR_0.1_kappa_100"]["most_marginal_safe"],
        "by_filter_label": by_filter,
        "verdict": {
            "preferred_filter": "omegaR_0.1_kappa_100",
            "central_global_S_T_max": pref["central"]["global_S_T_max"],
            "max_width_global_S_T_max": pref["max_width"]["global_S_T_max"],
            "central_worst_margin": pref["central"]["global_worst_margin"],
            "max_width_worst_margin": pref["max_width"]["global_worst_margin"],
            "S_T_1e_minus_5_robust_under_max_width": pref["max_width"]["global_worst_margin"] >= 1.0,
            "interpretation": verdict,
        },
    }


def trim(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "triplet_phase",
        "channel",
        "channel_class",
        "operator",
        "pair",
        "normalization_case",
        "S_T_max_replayed",
        "margin_replayed",
        "prefactor_ratio_new_over_old",
    ]
    return {key: row[key] for key in keys}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "filter_label",
        "scenario",
        "epsilon",
        "spectrum_name",
        "triplet_phase",
        "channel",
        "channel_class",
        "operator",
        "pair",
        "normalization_case",
        "old_width_prefactor",
        "new_width_prefactor",
        "prefactor_ratio_new_over_old",
        "S_T_max",
        "S_T_max_replayed",
        "tau_years_ST_display",
        "tau_years_replayed",
        "margin_at_ST_display",
        "margin_replayed",
        "passes_replayed",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_report(summary: dict[str, Any]) -> None:
    lines = [
        "# Chiral/lattice normalization replay",
        "",
        "No web lookup was used.",
        "",
        "The width prefactor is split as",
        "",
        "```text",
        "K = m_p/(64 pi f_pi^2) beta_H^2 A_R^2 C_chiral^2 Phi_2",
        "```",
        "",
        "## Channel prefactors",
        "",
        "| channel class | central K | K/K_old | phase space |",
        "|---|---:|---:|---:|",
    ]
    old = summary["old_single_width_prefactor"]
    for cls, k in summary["central_width_prefactors"].items():
        lines.append(f"| `{cls}` | {k:.6e} | {k / old:.3e} | {summary['phase_space'][cls]:.3e} |")
    lines.extend(
        [
            "",
            "## Filter replay",
            "",
            "| filter | case | global S_T max | worst margin | unsafe rows |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for filt, payload in summary["by_filter_label"].items():
        for case in ["central", "min_width", "max_width"]:
            row = payload[case]
            lines.append(
                f"| `{filt}` | `{case}` | {row['global_S_T_max']:.6e} | "
                f"{row['global_worst_margin']:.6e} | {row['unsafe_channel_rows']} |"
            )
    lines.extend(["", "## Verdict", "", summary["verdict"]["interpretation"], ""])
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = replay_rows()
    write_csv(OUT / "chiral_replay.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(summary)
    verdict = summary["verdict"]
    print("Chiral/lattice d=5 normalization replay")
    print(
        "  omegaR_0.1_kappa_100: "
        f"central_STmax={verdict['central_global_S_T_max']:.3e}, "
        f"max_width_STmax={verdict['max_width_global_S_T_max']:.3e}, "
        f"max_width_margin={verdict['max_width_worst_margin']:.3e}"
    )
    print(f"  robust S_T=1e-5: {verdict['S_T_1e_minus_5_robust_under_max_width']}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
