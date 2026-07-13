#!/usr/bin/env python3
"""Local CKM repair with the crossed-120 triplet projector frozen.

No web lookup is used.  The previous Clebsch-null bridge showed that a
crossed 120_A/120_B triplet projector conditionally closes the dressed d=5
future-stress gap for the source-consistent flavor row
``d5_both_F_minus_0``.  The remaining numerical obstruction is purely in the
doublet flavor sector: the CKM magnitude score is 1.453203e-3, above the
strict 1e-3 target.

This script performs a deliberately local audit.  It keeps the CP1 tensors
H,F,G_A,G_B and the crossed triplet projector fixed, then refits only the
doublet-Higgs mixing coefficients

    r_u,r_d,a_u,b_u,a_d,b_d,a_e,b_e,a_nu,b_nu.

Thus a passing row would mean that the proton-side closure survives a
doublet-only CKM repair; a failing row sharpens the obstruction without
retuning the triplet Wilson tensor.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402
import fit_two_kernel_flavor_then_d5 as two  # noqa: E402


SOURCE = ROOT / "output" / "two_kernel_flavor_then_d5" / "summary.json"
CROSSED = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"
OUT = ROOT / "output" / "source_consistent_ckm_crossed120"

SOURCE_LABEL = "d5_both_F_minus_0"
STRICT_CKM = 1.0e-3
LOOSE_MASS = 2.0e-1
SEESAW_RESIDUAL_MAX = 1.0e-10
RNG_SEED = 202605090911

MIX_KEYS = ["r_u", "r_d", "a_u", "b_u", "a_d", "b_d", "a_e", "b_e", "a_nu", "b_nu"]
MINIMAL_KEYS = ["r_u", "r_d", "a_u", "a_d", "a_e", "a_nu"]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def c(z: dict[str, float]) -> complex:
    return complex(float(z["re"]), float(z["im"]))


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def cvec(raw: list[dict[str, float]]) -> np.ndarray:
    return np.array([c(item) for item in raw], dtype=complex)


def pack_complex(v: np.ndarray) -> np.ndarray:
    return np.concatenate([np.real(v), np.imag(v)]).astype(float)


def unpack_complex(x: np.ndarray) -> np.ndarray:
    half = len(x) // 2
    return x[:half] + 1j * x[half:]


def load_crossed_margin(kappa: float = 30.0) -> tuple[float, dict[str, Any]]:
    payload = read_json(CROSSED)
    row = next(item for item in payload["finite_leakage_replay"] if float(item["kappa"]) == float(kappa))
    return float(row["worst_margin_1e35_replayed"]), payload


def source_card() -> dict[str, Any]:
    payload = read_json(SOURCE)
    return payload["matrix_cards"][SOURCE_LABEL]


@dataclass(frozen=True)
class FixedTensors:
    h: np.ndarray
    fmat: np.ndarray
    ga: np.ndarray
    gb: np.ndarray
    seed_mix: dict[str, complex]


def load_fixed_tensors() -> FixedTensors:
    card = source_card()
    model = card["model_data"]
    basis = two.load_veronese_basis()
    k = two.analytic_k()
    h = fit.symmetric_from_coeff(cvec(model["H_coefficients"]), basis) + c(model["epsilon_H"]) * k
    fmat = fit.symmetric_from_coeff(cvec(model["F_coefficients"]), basis) + c(model["epsilon_F"]) * k
    ga = fit.antisym_from_coeff(cvec(model["G_A_coefficients"]))
    gb = fit.antisym_from_coeff(cvec(model["G_B_coefficients"]))
    seed_mix = {key: c(model["mixing_coefficients"][key]) for key in MIX_KEYS}
    return FixedTensors(h=h, fmat=fmat, ga=ga, gb=gb, seed_mix=seed_mix)


def y_from_mix(tensors: FixedTensors, mix: dict[str, complex]) -> dict[str, np.ndarray]:
    h, fmat, ga, gb = tensors.h, tensors.fmat, tensors.ga, tensors.gb
    return {
        "up": h + mix["r_u"] * fmat + mix["a_u"] * ga + mix["b_u"] * gb,
        "down": mix["r_d"] * (h + fmat + mix["a_d"] * ga + mix["b_d"] * gb),
        "charged_lepton": mix["r_d"] * (h - 3.0 * fmat + mix["a_e"] * ga + mix["b_e"] * gb),
        "neutrino_dirac": h - 3.0 * mix["r_u"] * fmat + mix["a_nu"] * ga + mix["b_nu"] * gb,
    }


def score_flavor(y: dict[str, np.ndarray]) -> dict[str, Any]:
    v, ckm, _mass_quark = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = fit.log10_ratio(ckm["J"], fit.TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    return {
        "CKM_abs": fit.matrix_abs_json(v),
        "CKM_observables": ckm,
        "mass_ratios": {name: {"small": vals[0], "mid": vals[1]} for name, vals in masses.items()},
        "scores": {
            "ckm_magnitude_log_score": float(ckm_score),
            "jarlskog_log_score": float(j_score),
            "mass_log_score": float(mass_score),
            "total_unweighted": float(ckm_score + j_score + mass_score),
        },
    }


def mix_to_vector(seed: dict[str, complex], keys: list[str]) -> np.ndarray:
    return pack_complex(np.array([seed[key] for key in keys], dtype=complex))


def vector_to_mix(tensors: FixedTensors, keys: list[str], x: np.ndarray) -> dict[str, complex]:
    mix = dict(tensors.seed_mix)
    vals = unpack_complex(x)
    if len(keys) != len(vals):
        raise ValueError(f"key/value length mismatch: {len(keys)} != {len(vals)}")
    for key, val in zip(keys, vals):
        mix[key] = complex(val)
    return mix


def residuals(
    x: np.ndarray,
    tensors: FixedTensors,
    keys: list[str],
    x0: np.ndarray,
    radius: float,
    weights: dict[str, float],
) -> np.ndarray:
    mix = vector_to_mix(tensors, keys, x)
    y = y_from_mix(tensors, mix)
    flv = score_flavor(y)
    ckm = flv["CKM_observables"]
    masses = flv["mass_ratios"]
    pieces = [
        weights["ckm"] * fit.log10_ratio(ckm["Vus"], fit.TARGETS["Vus"]),
        weights["ckm"] * fit.log10_ratio(ckm["Vcb"], fit.TARGETS["Vcb"]),
        weights["ckm"] * fit.log10_ratio(ckm["Vub"], fit.TARGETS["Vub"]),
        weights["jarlskog"] * fit.log10_ratio(ckm["J"], fit.TARGETS["J"]),
    ]
    for sector, vals in masses.items():
        pieces.append(weights["mass"] * fit.log10_ratio(vals["small"], fit.TARGETS[f"{sector}_small"]))
        pieces.append(weights["mass"] * fit.log10_ratio(vals["mid"], fit.TARGETS[f"{sector}_mid"]))
    pieces.append(weights["reg"] * float(np.linalg.norm((x - x0) / max(radius, 1.0e-30))))
    return np.array(pieces, dtype=float)


def evaluate(
    *,
    label: str,
    mode: str,
    tensors: FixedTensors,
    keys: list[str],
    x: np.ndarray,
    x0: np.ndarray,
    radius: float,
    d5_margin: float,
    optimizer: dict[str, Any] | None = None,
) -> dict[str, Any]:
    mix = vector_to_mix(tensors, keys, x)
    y = y_from_mix(tensors, mix)
    flv = score_flavor(y)
    replay = fit.seesaw_replay(y["neutrino_dirac"], y["charged_lepton"])
    delta = x - x0
    scores = flv["scores"]
    row = {
        "label": label,
        "mode": mode,
        "varied_keys": keys,
        "radius": radius,
        "ckm_score": scores["ckm_magnitude_log_score"],
        "jarlskog_score": scores["jarlskog_log_score"],
        "mass_score": scores["mass_log_score"],
        "future_margin_1e35_crossed120_kappa30": d5_margin,
        "passes_strict_ckm": scores["ckm_magnitude_log_score"] < STRICT_CKM,
        "passes_loose_mass": scores["mass_log_score"] <= LOOSE_MASS,
        "passes_seesaw": replay["seesaw_matrix_residual"] < SEESAW_RESIDUAL_MAX,
        "passes_future_d5": d5_margin >= 1.0,
        "passes_conditional_strict_closure": bool(
            scores["ckm_magnitude_log_score"] < STRICT_CKM
            and scores["mass_log_score"] <= LOOSE_MASS
            and replay["seesaw_matrix_residual"] < SEESAW_RESIDUAL_MAX
            and d5_margin >= 1.0
        ),
        "CKM_observables": flv["CKM_observables"],
        "mass_ratios": flv["mass_ratios"],
        "seesaw_replay": replay,
        "doublet_deformation": {
            "delta_l2": float(np.linalg.norm(delta)),
            "delta_linf": float(np.max(np.abs(delta))) if len(delta) else 0.0,
            "relative_to_seed_var_l2": float(np.linalg.norm(delta) / max(np.linalg.norm(x0), 1.0e-30)),
            "active_bound_fraction": float(np.mean(np.isclose(np.abs(delta), radius, rtol=0.0, atol=1.0e-5))),
        },
        "mixing_coefficients": {key: cjson(val) for key, val in mix.items()},
        "Yukawa_fit": {name: fit.matrix_json(mat) for name, mat in y.items()},
        "optimizer": optimizer or {},
    }
    return row


def run_mode(
    *,
    mode: str,
    tensors: FixedTensors,
    keys: list[str],
    radius: float,
    weights: dict[str, float],
    d5_margin: float,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    x0 = mix_to_vector(tensors.seed_mix, keys)
    lower = x0 - radius
    upper = x0 + radius
    starts = [x0]
    # Small deterministic local nudges are enough to avoid a bad finite
    # difference basin without turning this into a broad random scan.
    for scale in [0.15, -0.15, 0.35, -0.35]:
        starts.append(np.clip(x0 + scale * radius * rng.normal(size=x0.shape), lower, upper))

    rows: list[dict[str, Any]] = []
    for idx, start in enumerate(starts):
        res = least_squares(
            residuals,
            start,
            bounds=(lower, upper),
            args=(tensors, keys, x0, radius, weights),
            xtol=2.0e-11,
            ftol=2.0e-11,
            gtol=2.0e-11,
            max_nfev=650,
        )
        rows.append(
            evaluate(
                label=f"{mode}_{idx}",
                mode=mode,
                tensors=tensors,
                keys=keys,
                x=res.x,
                x0=x0,
                radius=radius,
                d5_margin=d5_margin,
                optimizer={
                    "success": bool(res.success),
                    "nfev": int(res.nfev),
                    "cost": float(2.0 * res.cost),
                    "message": str(res.message),
                },
            )
        )
    return rows


def flat(row: dict[str, Any]) -> dict[str, Any]:
    ckm = row["CKM_observables"]
    deform = row["doublet_deformation"]
    return {
        "label": row["label"],
        "mode": row["mode"],
        "ckm_score": row["ckm_score"],
        "jarlskog_score": row["jarlskog_score"],
        "mass_score": row["mass_score"],
        "future_margin": row["future_margin_1e35_crossed120_kappa30"],
        "Vus": ckm["Vus"],
        "Vcb": ckm["Vcb"],
        "Vub": ckm["Vub"],
        "J_abs": ckm["J"],
        "seesaw_residual": row["seesaw_replay"]["seesaw_matrix_residual"],
        "delta_l2": deform["delta_l2"],
        "delta_linf": deform["delta_linf"],
        "relative_delta": deform["relative_to_seed_var_l2"],
        "active_bound_fraction": deform["active_bound_fraction"],
        "passes_conditional_strict_closure": row["passes_conditional_strict_closure"],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    tensors = load_fixed_tensors()
    d5_margin, crossed_payload = load_crossed_margin(30.0)

    source_x0 = mix_to_vector(tensors.seed_mix, MIX_KEYS)
    source = evaluate(
        label="source_seed",
        mode="seed",
        tensors=tensors,
        keys=MIX_KEYS,
        x=source_x0,
        x0=source_x0,
        radius=0.0,
        d5_margin=d5_margin,
    )

    configs = [
        {
            "mode": "minimal_6_radius_0p08",
            "keys": MINIMAL_KEYS,
            "radius": 0.08,
            "weights": {"ckm": 16.0, "jarlskog": 0.18, "mass": 34.0, "reg": 0.018},
        },
        {
            "mode": "minimal_6_radius_0p18",
            "keys": MINIMAL_KEYS,
            "radius": 0.18,
            "weights": {"ckm": 18.0, "jarlskog": 0.16, "mass": 32.0, "reg": 0.015},
        },
        {
            "mode": "all10_radius_0p08",
            "keys": MIX_KEYS,
            "radius": 0.08,
            "weights": {"ckm": 16.0, "jarlskog": 0.16, "mass": 34.0, "reg": 0.020},
        },
        {
            "mode": "all10_radius_0p18",
            "keys": MIX_KEYS,
            "radius": 0.18,
            "weights": {"ckm": 18.0, "jarlskog": 0.14, "mass": 30.0, "reg": 0.016},
        },
        {
            "mode": "all10_radius_0p30_ckmheavy",
            "keys": MIX_KEYS,
            "radius": 0.30,
            "weights": {"ckm": 26.0, "jarlskog": 0.10, "mass": 24.0, "reg": 0.018},
        },
    ]

    rows = [source]
    for cfg in configs:
        rows.extend(
            run_mode(
                mode=cfg["mode"],
                tensors=tensors,
                keys=cfg["keys"],
                radius=float(cfg["radius"]),
                weights=cfg["weights"],
                d5_margin=d5_margin,
                rng=rng,
            )
        )

    flat_rows = [flat(row) for row in rows]
    with (OUT / "scan_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0].keys()))
        writer.writeheader()
        writer.writerows(flat_rows)

    best_ckm = min(rows, key=lambda row: row["ckm_score"])
    valid_mass = [row for row in rows if row["mass_score"] <= LOOSE_MASS]
    best_valid_mass_ckm = min(valid_mass, key=lambda row: row["ckm_score"]) if valid_mass else None
    closures = [row for row in rows if row["passes_conditional_strict_closure"]]
    best_closure = min(closures, key=lambda row: row["ckm_score"]) if closures else None

    payload = {
        "note": "No web lookup used. Source-consistent CKM repair with crossed-120 triplet projector frozen.",
        "mathematical_setup": {
            "fixed_tensors": "H,F,G_A,G_B are reconstructed from output/two_kernel_flavor_then_d5 matrix card d5_both_F_minus_0.",
            "varied_doublet_coefficients": "r_u,r_d,a_u,b_u,a_d,b_d,a_e,b_e,a_nu,b_nu, or the minimal 6-key subset.",
            "frozen_triplet_projector": "crossed 120_A/120_B finite lift kappa=30 from output/crossed_120_triplet_projector.",
            "closure_test": "CKM<1e-3, mass<=0.2, seesaw residual<1e-10, crossed d5 future margin>=1.",
        },
        "config": {
            "rng_seed": RNG_SEED,
            "strict_ckm": STRICT_CKM,
            "loose_mass": LOOSE_MASS,
            "seesaw_residual_max": SEESAW_RESIDUAL_MAX,
            "configs": configs,
        },
        "selection_rule_caveat": {
            "field_only_unbroken_spin10_projector_possible": crossed_payload["verdict"]["field_only_unbroken_spin10_projector_possible"],
            "ps_eft_or_constrained_projector_possible": crossed_payload["verdict"]["ps_eft_or_constrained_projector_possible"],
            "interpretation": crossed_payload["verdict"]["interpretation"],
        },
        "source_seed": source,
        "rows": rows,
        "verdict": {
            "row_count": len(rows),
            "strict_conditional_closure_count": len(closures),
            "best_ckm": flat(best_ckm),
            "best_mass_valid_ckm": flat(best_valid_mass_ckm) if best_valid_mass_ckm else None,
            "best_closure": flat(best_closure) if best_closure else None,
            "d5_gap_closed_by_frozen_crossed_projector": d5_margin >= 1.0,
            "source_seed_ckm_factor_to_strict": source["ckm_score"] / STRICT_CKM,
            "best_mass_valid_ckm_factor_to_strict": (
                best_valid_mass_ckm["ckm_score"] / STRICT_CKM if best_valid_mass_ckm else math.inf
            ),
            "local_conditional_closure_complete": bool(closures),
            "publication_complete": False,
            "interpretation": (
                "This is a local doublet-only audit.  It does not repair the "
                "field-only Spin(10) no-go for the crossed triplet projector. "
                "It asks whether strict CKM can be reached without retuning "
                "the d=5 projector or leaving the source-consistent tensor card."
            ),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Source-consistent CKM repair with crossed-120 frozen",
        "",
        "No web lookup was used.",
        "",
        "The triplet side is frozen to the crossed 120_A/120_B finite-lift projector at kappa=30.",
        "Only doublet-Higgs mixing coefficients are varied.",
        "",
        "| role | mode | CKM score | mass score | d5 margin | seesaw residual | rel. delta | closure |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    selected = [("seed", source), ("best_ckm", best_ckm)]
    if best_valid_mass_ckm is not None:
        selected.append(("best_mass_valid_ckm", best_valid_mass_ckm))
    if best_closure is not None:
        selected.append(("best_closure", best_closure))
    seen: set[str] = set()
    for role, row in selected:
        if row["label"] in seen:
            continue
        seen.add(row["label"])
        report.append(
            f"| `{role}` | `{row['mode']}` | {row['ckm_score']:.6e} | "
            f"{row['mass_score']:.6e} | {row['future_margin_1e35_crossed120_kappa30']:.6e} | "
            f"{row['seesaw_replay']['seesaw_matrix_residual']:.3e} | "
            f"{row['doublet_deformation']['relative_to_seed_var_l2']:.3e} | "
            f"`{row['passes_conditional_strict_closure']}` |"
        )
    report += [
        "",
        "## Verdict",
        "",
        f"Strict conditional closures: `{len(closures)}`.",
        f"Frozen crossed-projector future d=5 margin: `{d5_margin:.6e}`.",
        "",
        payload["verdict"]["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Source-consistent CKM repair with crossed-120 frozen")
    print(f"  rows: {len(rows)}")
    print(f"  strict closures: {len(closures)}")
    print(f"  best CKM: {best_ckm['ckm_score']:.6e}")
    if best_valid_mass_ckm is not None:
        print(f"  best mass-valid CKM: {best_valid_mass_ckm['ckm_score']:.6e}")
    print(f"  crossed d5 margin: {d5_margin:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
