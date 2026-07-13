#!/usr/bin/env python3
"""Pareto/homotopy scan for the CP1/O(-4) Clebsch flavor branch.

No web lookup is used.  The previous baseline showed that the restricted
10_H + overline{126}_H + 120_H ansatz can make CKM small, but at the price of a
large mass-hierarchy score.  This script starts from that CKM-heavy point and
adiabatically increases the mass-hierarchy weight, recording the Pareto curve.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402


BASELINE = ROOT / "output" / "flavor_clebsch_scan" / "scan_clebsch_flavor_fit.json"
OUT = ROOT / "output" / "flavor_clebsch_homotopy"


STAGES = [
    {"name": "ckm_lock", "ckm": 4.0, "jarlskog": 0.25, "mass": 0.6, "reg": 0.003},
    {"name": "mass_1", "ckm": 4.0, "jarlskog": 0.25, "mass": 1.0, "reg": 0.003},
    {"name": "mass_2", "ckm": 4.2, "jarlskog": 0.25, "mass": 2.0, "reg": 0.003},
    {"name": "mass_4", "ckm": 4.5, "jarlskog": 0.25, "mass": 4.0, "reg": 0.004},
    {"name": "mass_8", "ckm": 5.0, "jarlskog": 0.25, "mass": 8.0, "reg": 0.005},
]


def load_basis() -> np.ndarray:
    card = json.loads(fit.CARD.read_text(encoding="utf-8"))
    return np.array(
        [[[cell["re"] + 1j * cell["im"] for cell in row] for row in plane] for plane in card["product_clebsch_C_ij_m"]],
        dtype=complex,
    )


def run_stage(
    stage: dict[str, float | str],
    seed: np.ndarray,
    basis: np.ndarray,
    rng: np.random.Generator,
    bound: float,
) -> dict[str, Any]:
    weights = {
        "ckm": float(stage["ckm"]),
        "jarlskog": float(stage["jarlskog"]),
        "mass": float(stage["mass"]),
        "reg": float(stage["reg"]) / bound,
    }
    candidates = [np.clip(seed, -0.95 * bound, 0.95 * bound)]
    for scale in [0.015, 0.040]:
        candidates.append(np.clip(seed + rng.normal(0.0, scale, size=seed.shape), -0.95 * bound, 0.95 * bound))

    best = None
    for candidate in candidates:
        res = least_squares(
            fit.objective,
            candidate,
            bounds=(-bound * np.ones_like(seed), bound * np.ones_like(seed)),
            args=(basis, weights),
            xtol=3.0e-10,
            ftol=3.0e-10,
            gtol=3.0e-10,
            max_nfev=850,
        )
        item = fit.evaluate(res.x, basis, str(stage["name"]), bound)
        item["stage_weights"] = weights
        item["optimizer"] = {
            "success": bool(res.success),
            "nfev": int(res.nfev),
            "weighted_cost": float(np.sum(fit.objective(res.x, basis, weights) ** 2)),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(res.x), bound, rtol=0.0, atol=1.0e-5))),
        }
        if best is None or item["optimizer"]["weighted_cost"] < best["optimizer"]["weighted_cost"]:
            best = item
    assert best is not None
    return best


def flat_row(item: dict[str, Any]) -> dict[str, Any]:
    obs = item["CKM_observables"]
    scores = item["scores"]
    replay = item["seesaw_replay"]
    return {
        "stage": item["profile"],
        "mass_weight": item["stage_weights"]["mass"],
        "ckm_weight": item["stage_weights"]["ckm"],
        "Vus": obs["Vus"],
        "Vcb": obs["Vcb"],
        "Vub": obs["Vub"],
        "J_abs": obs["J"],
        "ckm_score": scores["ckm_magnitude_log_score"],
        "mass_score": scores["mass_log_score"],
        "j_score": scores["jarlskog_log_score"],
        "theta_norm": replay["theta_norm"],
        "seesaw_residual": replay["seesaw_matrix_residual"],
        "M1_GeV": replay["heavy_neutrino_masses_GeV"][0],
        "M2_GeV": replay["heavy_neutrino_masses_GeV"][1],
        "M3_GeV": replay["heavy_neutrino_masses_GeV"][2],
        "weighted_cost": item["optimizer"]["weighted_cost"],
    }


def write_csv(items: list[dict[str, Any]]) -> None:
    rows = [flat_row(item) for item in items]
    with (OUT / "scan_clebsch_flavor_homotopy.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Clebsch flavor homotopy scan")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("| stage | mass weight | Vus | Vcb | Vub | J | CKM score | mass score | theta |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for item in payload["stages"]:
        row = flat_row(item)
        lines.append(
            f"| {row['stage']} | {row['mass_weight']:.2f} | {row['Vus']:.4e} | "
            f"{row['Vcb']:.4e} | {row['Vub']:.4e} | {row['J_abs']:.4e} | "
            f"{row['ckm_score']:.3e} | {row['mass_score']:.3e} | {row['theta_norm']:.3e} |"
        )
    best = payload["best_pareto"]
    obs = best["CKM_observables"]
    scores = best["scores"]
    replay = best["seesaw_replay"]
    lines.append("")
    lines.append("## Best Pareto point")
    lines.append("")
    lines.append(f"Stage: `{best['profile']}`.")
    lines.append(
        f"CKM: `|Vus|={obs['Vus']:.6e}, |Vcb|={obs['Vcb']:.6e}, "
        f"|Vub|={obs['Vub']:.6e}, |J|={obs['J']:.6e}`."
    )
    lines.append(
        f"Scores: `CKM={scores['ckm_magnitude_log_score']:.6e}`, "
        f"`mass={scores['mass_log_score']:.6e}`."
    )
    heavy = replay["heavy_neutrino_masses_GeV"]
    lines.append(
        f"Seesaw: `theta={replay['theta_norm']:.6e}`, "
        f"`residual={replay['seesaw_matrix_residual']:.6e}`, "
        f"`M=({heavy[0]:.6e}, {heavy[1]:.6e}, {heavy[2]:.6e}) GeV`."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "scan_clebsch_flavor_homotopy_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    basis = load_basis()
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    seed = np.array(baseline["best"]["x"], dtype=float)
    rng = np.random.default_rng(2026050720)
    bound = 2.8

    stages = []
    current = seed
    print("Clebsch flavor homotopy scan", flush=True)
    for stage in STAGES:
        item = run_stage(stage, current, basis, rng, bound)
        stages.append(item)
        current = np.array(item["x"], dtype=float)
        obs = item["CKM_observables"]
        scores = item["scores"]
        print(
            f"  {stage['name']:8s} mass_w={stage['mass']:.1f} "
            f"Vus={obs['Vus']:.4e} Vcb={obs['Vcb']:.4e} Vub={obs['Vub']:.4e} "
            f"ckm={scores['ckm_magnitude_log_score']:.3e} "
            f"mass={scores['mass_log_score']:.3e}",
            flush=True,
        )

    viable_candidates = [
        item
        for item in stages
        if item["scores"]["ckm_magnitude_log_score"] < 5.0e-2
        and item["scores"]["mass_log_score"] < 2.0e-1
    ]
    if viable_candidates:
        best = min(viable_candidates, key=lambda item: item["scores"]["mass_log_score"])
    else:
        # Rank the Pareto curve by a conservative scalarized criterion that keeps
        # CKM meaningful while still rewarding mass improvement.
        best = min(stages, key=lambda item: item["scores"]["mass_log_score"] + 4.0 * item["scores"]["ckm_magnitude_log_score"])

    best_obs = best["CKM_observables"]
    best_scores = best["scores"]
    verdict = (
        "The homotopy scan quantifies the Clebsch branch as a Pareto tension.  "
        f"The best scalarized point has |Vus|={best_obs['Vus']:.6e}, "
        f"|Vcb|={best_obs['Vcb']:.6e}, |Vub|={best_obs['Vub']:.6e}, "
        f"CKM score={best_scores['ckm_magnitude_log_score']:.6e}, and "
        f"mass score={best_scores['mass_log_score']:.6e}.  If no viability check "
        "passes, the next representation-controlled move is to add the second "
        "independent 120_H doublet Clebsch direction or relax the CP1/O(-4) "
        "restriction only for one controlled symmetric tensor."
    )
    payload = {
        "note": "No web lookup used. Homotopy scan from the CKM-heavy Clebsch seed.",
        "input_baseline": str(BASELINE),
        "stages": stages,
        "best_pareto": best,
        "checks": {
            "scan_completed": True,
            "any_ckm_score_lt_0p05": any(item["scores"]["ckm_magnitude_log_score"] < 5.0e-2 for item in stages),
            "any_mass_score_lt_0p20": any(item["scores"]["mass_log_score"] < 2.0e-1 for item in stages),
            "any_joint_viable": bool(viable_candidates),
            "best_seesaw_residual_lt_1e_minus_10": bool(best["seesaw_replay"]["seesaw_matrix_residual"] < 1.0e-10),
            "best_theta_lt_1e_minus_8": bool(best["seesaw_replay"]["theta_norm"] < 1.0e-8),
        },
        "verdict": verdict,
    }
    payload["phenomenology_viable"] = payload["checks"]["any_joint_viable"]
    (OUT / "scan_clebsch_flavor_homotopy.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(stages)
    write_report(payload)
    print(f"  best pareto: {best['profile']}", flush=True)
    print(f"  phenomenology viable: {payload['phenomenology_viable']}", flush=True)
    print(f"  wrote: {OUT}", flush=True)


if __name__ == "__main__":
    main()
