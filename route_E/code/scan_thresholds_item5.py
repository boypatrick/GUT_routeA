#!/usr/bin/env python3
"""Item 5: one-loop threshold-correction scan with proton lower bound.

No web lookup is used.  Input low-energy constants are the same local
benchmarks used in earlier scans.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROTON_JSON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "thresholds"

MZ_GEV = 91.1876
ALPHA_EM_INV_MZ = 127.955
SIN2_THETA_W_MZ = 0.23122
ALPHA_S_MZ = 0.1184
HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


@dataclass(frozen=True)
class BetaModel:
    name: str
    b: tuple[float, float, float]
    note: str


MODELS = [
    BetaModel("SM", (41.0 / 10.0, -19.0 / 6.0, -7.0), "one Higgs doublet, non-SUSY one-loop beta"),
    BetaModel("MSSM", (33.0 / 5.0, 1.0, -3.0), "MSSM one-loop beta above the superpartner threshold"),
]


def alpha_inverse_mz() -> np.ndarray:
    alpha_em = 1.0 / ALPHA_EM_INV_MZ
    cos2 = 1.0 - SIN2_THETA_W_MZ
    alpha1 = (5.0 / 3.0) * alpha_em / cos2
    alpha2 = alpha_em / SIN2_THETA_W_MZ
    alpha3 = ALPHA_S_MZ
    return np.array([1.0 / alpha1, 1.0 / alpha2, 1.0 / alpha3], dtype=float)


def load_proton_constants() -> dict[str, float]:
    payload = json.loads(PROTON_JSON.read_text(encoding="utf-8"))
    return {
        "K_GeV5": float(payload["hadronic_constants"]["width_prefactor_GeV5"]),
        "tau_target_years": 1.0e34,
    }


def mx_required(alpha_g_inv: float, constants: dict[str, float]) -> float:
    width_target = HBAR_GEV_S / (constants["tau_target_years"] * SECONDS_PER_YEAR)
    c6_max = math.sqrt(width_target / constants["K_GeV5"])
    g2 = 4.0 * math.pi / alpha_g_inv
    return math.sqrt(g2 / c6_max)


def required_threshold(alpha_inv_mz: np.ndarray, b: np.ndarray, mg: float) -> tuple[float, np.ndarray, np.ndarray]:
    log = math.log(mg / MZ_GEV)
    unification_estimates = alpha_inv_mz - b * log / (2.0 * math.pi)
    alpha_g_inv = float(np.mean(unification_estimates))
    delta = unification_estimates - alpha_g_inv
    diff = np.array([delta[0] - delta[1], delta[1] - delta[2]], dtype=float)
    return alpha_g_inv, delta, diff


def pairwise_crossings(alpha_inv_mz: np.ndarray, b: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {}
    pairs = [("12", 0, 1), ("23", 1, 2), ("13", 0, 2)]
    for label, i, j in pairs:
        log = 2.0 * math.pi * (alpha_inv_mz[i] - alpha_inv_mz[j]) / (b[i] - b[j])
        out[f"M{label}_GeV"] = MZ_GEV * math.exp(log)
        out[f"alpha{label}_inv"] = float(alpha_inv_mz[i] - b[i] * log / (2.0 * math.pi))
    return out


def scalar_threshold_basis_fit(delta: np.ndarray, mg: float) -> dict[str, object]:
    # Complex scalar fragments in the SM normalization:
    # H_C: (3,1,-1/3), Sigma_3: real (1,3,0), Sigma_8: real (8,1,0).
    # Delta_i = sum_r b_i^r log(M_G/M_r)/(2pi).
    basis = np.array(
        [
            [1.0 / 15.0, 0.0, 0.0],
            [0.0, 1.0 / 3.0, 0.0],
            [1.0 / 6.0, 0.0, 1.0 / 2.0],
        ],
        dtype=float,
    )
    logs, *_ = np.linalg.lstsq(basis / (2.0 * math.pi), delta, rcond=None)
    masses = mg * np.exp(-logs)
    reconstructed = basis @ logs / (2.0 * math.pi)
    residual = float(np.linalg.norm(reconstructed - delta))
    return {
        "basis": ["H_C_scalar_(3,1,-1/3)", "Sigma3_real_(1,3,0)", "Sigma8_real_(8,1,0)"],
        "logs_ln_MG_over_Mr": [float(x) for x in logs],
        "masses_GeV": [float(x) for x in masses],
        "max_abs_log": float(np.max(np.abs(logs))),
        "relative_residual": residual / max(float(np.linalg.norm(delta)), 1e-30),
    }


def scan_model(model: BetaModel, constants: dict[str, float]) -> dict[str, object]:
    alpha_inv = alpha_inverse_mz()
    b = np.array(model.b, dtype=float)
    rows = []
    best = None
    best_proton_safe = None
    best_split_safe = None

    for log10_mg in np.linspace(14.0, 18.0, 1601):
        mg = 10.0**log10_mg
        alpha_g_inv, delta, diff = required_threshold(alpha_inv, b, mg)
        if alpha_g_inv <= 0:
            continue
        mreq = mx_required(alpha_g_inv, constants)
        rho_min = max(1.0, mreq / mg)
        ln_rho_min = math.log(rho_min)
        l2 = float(np.linalg.norm(delta))
        span = float(np.max(delta) - np.min(delta))
        row = {
            "model": model.name,
            "log10_MG": log10_mg,
            "MG_GeV": mg,
            "alphaG_inv": alpha_g_inv,
            "Delta1": float(delta[0]),
            "Delta2": float(delta[1]),
            "Delta3": float(delta[2]),
            "Delta12": float(diff[0]),
            "Delta23": float(diff[1]),
            "Delta_l2": l2,
            "Delta_span": span,
            "M_X_required_GeV": mreq,
            "rho_X_min": rho_min,
            "ln_rho_X_min": ln_rho_min,
            "proton_safe_if_MX_equals_MG": mg >= mreq,
            "proton_safe_within_factor10_split": ln_rho_min <= math.log(10.0),
        }
        rows.append(row)

        objective = l2 + 0.50 * max(0.0, ln_rho_min - math.log(10.0)) ** 2
        if best is None or objective < best[0]:
            best = (objective, row)
        if row["proton_safe_if_MX_equals_MG"] and (best_proton_safe is None or l2 < best_proton_safe[0]):
            best_proton_safe = (l2, row)
        if row["proton_safe_within_factor10_split"] and (best_split_safe is None or l2 < best_split_safe[0]):
            best_split_safe = (l2, row)

    crossings = pairwise_crossings(alpha_inv, b)
    selected = {
        "global_best_with_split_penalty": best[1],
        "best_MX_equals_MG_proton_safe": best_proton_safe[1] if best_proton_safe else None,
        "best_with_MX_within_factor10": best_split_safe[1] if best_split_safe else None,
    }
    for key, row in list(selected.items()):
        if row is None:
            continue
        delta = np.array([row["Delta1"], row["Delta2"], row["Delta3"]], dtype=float)
        row["scalar_threshold_basis_fit"] = scalar_threshold_basis_fit(delta, row["MG_GeV"])

    return {
        "model": model.name,
        "note": model.note,
        "beta": list(model.b),
        "pairwise_crossings": crossings,
        "selected": selected,
        "rows": rows,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    flattened = []
    for row in rows:
        r = {key: value for key, value in row.items() if key != "scalar_threshold_basis_fit"}
        flattened.append(r)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flattened[0].keys()))
        writer.writeheader()
        writer.writerows(flattened)


def write_report(payload: dict[str, object]) -> None:
    lines: list[str] = []
    lines.append("# Item 5 threshold-correction scan")
    lines.append("")
    lines.append("No web lookup was used. Low-energy inputs are local benchmark constants:")
    lines.append("")
    lines.append("```text")
    lines.append(f"alpha_em^-1(MZ) = {ALPHA_EM_INV_MZ}")
    lines.append(f"sin^2 theta_W(MZ) = {SIN2_THETA_W_MZ}")
    lines.append(f"alpha_s(MZ) = {ALPHA_S_MZ}")
    lines.append("alpha_1 = (5/3) alpha_Y")
    lines.append("```")
    lines.append("")
    lines.append("Matching convention:")
    lines.append("")
    lines.append("```text")
    lines.append("alpha_i^-1(MZ) = alpha_G^-1 + b_i/(2pi) log(M_G/MZ) + Delta_i")
    lines.append("sum_i Delta_i = 0 fixes the unobservable common threshold shift.")
    lines.append("```")
    lines.append("")
    lines.append("The proton bound from item 4 is imposed as")
    lines.append("")
    lines.append("```text")
    lines.append("M_X >= M_X^min(alpha_G),  tau(p -> e+ pi0) > 1e34 yr.")
    lines.append("```")
    lines.append("")
    for model in payload["models"]:
        lines.append(f"## {model['model']} baseline")
        lines.append("")
        lines.append(f"Beta coefficients: `{model['beta']}`.")
        c = model["pairwise_crossings"]
        lines.append("")
        lines.append("Pairwise one-loop crossings without thresholds:")
        lines.append("")
        lines.append("```text")
        lines.append(f"M12 = {c['M12_GeV']:.6e} GeV, alpha12^-1 = {c['alpha12_inv']:.6f}")
        lines.append(f"M23 = {c['M23_GeV']:.6e} GeV, alpha23^-1 = {c['alpha23_inv']:.6f}")
        lines.append(f"M13 = {c['M13_GeV']:.6e} GeV, alpha13^-1 = {c['alpha13_inv']:.6f}")
        lines.append("```")
        lines.append("")
        for label, title in [
            ("global_best_with_split_penalty", "Best allowing M_X above M_G"),
            ("best_MX_equals_MG_proton_safe", "Best with M_X = M_G and proton safe"),
            ("best_with_MX_within_factor10", "Best with M_X/M_G <= 10"),
        ]:
            row = model["selected"][label]
            if row is None:
                lines.append(f"{title}: none.")
                continue
            lines.append(f"### {title}")
            lines.append("")
            lines.append("```text")
            lines.append(f"M_G = {row['MG_GeV']:.6e} GeV")
            lines.append(f"alpha_G^-1 = {row['alphaG_inv']:.6f}")
            lines.append(
                f"Delta = ({row['Delta1']:+.6f}, {row['Delta2']:+.6f}, {row['Delta3']:+.6f})"
            )
            lines.append(f"||Delta||_2 = {row['Delta_l2']:.6f}, span = {row['Delta_span']:.6f}")
            lines.append(f"M_X^min = {row['M_X_required_GeV']:.6e} GeV")
            lines.append(f"rho_X^min = M_X^min/M_G clipped at 1 = {row['rho_X_min']:.6f}")
            lines.append(f"ln rho_X^min = {row['ln_rho_X_min']:.6f}")
            lines.append("```")
            fit = row["scalar_threshold_basis_fit"]
            lines.append("")
            lines.append("Illustrative scalar-fragment threshold fit:")
            lines.append("")
            lines.append("```text")
            for name, logv, mass in zip(fit["basis"], fit["logs_ln_MG_over_Mr"], fit["masses_GeV"]):
                lines.append(f"{name}: ln(M_G/M_r)={logv:+.4f}, M_r={mass:.6e} GeV")
            lines.append(f"max |log| = {fit['max_abs_log']:.4f}, residual = {fit['relative_residual']:.3e}")
            lines.append("```")
            lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- Threshold corrections live in a two-dimensional difference plane; complete GUT multiplets move only the common alpha_G direction.")
    lines.append("- The proton-safe cone removes threshold solutions that lower the broken gauge boson below the item-4 bound.")
    lines.append("- The MSSM-like beta coefficients need only very small threshold corrections near 2e16 GeV and are proton safe.")
    lines.append("- The SM beta coefficients need multi-unit threshold differences even after the proton bound; this is possible only with sizable split multiplets or an Ed gauge-kinetic contribution.")
    lines.append("- The scalar-fragment fit is illustrative, not a final spectrum: enormous logs mean the chosen fragment basis is an unnatural explanation for that point.")
    (OUT / "threshold_item5_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def json_clean(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): json_clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_clean(v) for v in value]
    if isinstance(value, tuple):
        return [json_clean(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    constants = load_proton_constants()
    models = [scan_model(model, constants) for model in MODELS]
    all_rows = []
    for model in models:
        all_rows.extend(model["rows"])
    write_csv(OUT / "threshold_scan.csv", all_rows)
    summary = {
        "input": {
            "proton_json": str(PROTON_JSON),
            "MZ_GeV": MZ_GEV,
            "alpha_em_inv_MZ": ALPHA_EM_INV_MZ,
            "sin2_theta_W_MZ": SIN2_THETA_W_MZ,
            "alpha_s_MZ": ALPHA_S_MZ,
            "alpha_inverse_MZ_GUT_normalized": [float(x) for x in alpha_inverse_mz()],
        },
        "proton_bound": constants,
        "models": [
            {key: value for key, value in model.items() if key != "rows"}
            for model in models
        ],
    }
    summary = json_clean(summary)
    (OUT / "threshold_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary)

    print("Item 5 threshold scan")
    print(f"  alpha_i^-1(MZ): {', '.join(f'{x:.6f}' for x in alpha_inverse_mz())}")
    for model in summary["models"]:
        print(f"  {model['model']}:")
        for label, row in model["selected"].items():
            if row is None:
                print(f"    {label}: none")
                continue
            print(
                "    {label}: MG={mg:.3e} alphaG^-1={ag:.3f} "
                "Delta=({d1:+.3f},{d2:+.3f},{d3:+.3f}) "
                "rhoX={rho:.3f} l2={l2:.3f}".format(
                    label=label,
                    mg=row["MG_GeV"],
                    ag=row["alphaG_inv"],
                    d1=row["Delta1"],
                    d2=row["Delta2"],
                    d3=row["Delta3"],
                    rho=row["rho_X_min"],
                    l2=row["Delta_l2"],
                )
            )
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
