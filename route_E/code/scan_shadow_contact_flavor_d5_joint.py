#!/usr/bin/env python3
"""Joint shadow/contact flavor and Knu d=5 scan.

No web lookup is used.  This script tests the next nontrivial idea suggested
by the roadmap: the same spin-zero contact/shadow direction that repairs the
Veronese flavor obstruction might also reduce the dangerous Knu dimension-five
amplitude.

Starting from the current transvectant flavor point, scan

    Y_a -> Y_a + eps_a K,

where

    K = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]]

is the CP1 O(2) second transvectant.  For each deformed point, compute

  * CKM and mass-ratio scores,
  * a seesaw replay,
  * the mass-basis Knu Wilson amplitude in the same proxy convention as
    scan_dimension5_wilson_tensors.py,
  * an inferred 1e35 yr Knu future margin by rescaling the calibrated local
    target-map margin.

This is not the final proton-decay calculation.  It is a Pareto-direction
audit: does a common contact/shadow deformation move CKM and Knu in the right
direction at all?
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

import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


TRANSVECTANT = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "shadow_contact_flavor_d5_joint"

RNG_SEED = 202605082123
RANDOM_PER_MODE = 2500
BUDGETS = [0.0, 0.015, 0.03, 0.06, 0.10, 0.16]
MODES = ["down_only", "up_down", "down_lepton_tied", "all_independent"]
DISPLAY_ST = 7.5e-6


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def analytic_k() -> np.ndarray:
    return np.array(
        [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        dtype=complex,
    ) / math.sqrt(3.0)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def load_seed() -> dict[str, np.ndarray]:
    payload = read_json(TRANSVECTANT)
    return {
        name: fit.cmat(raw)
        for name, raw in payload["Yukawa_matrices"].items()
    }


def random_complex_vector(rng: np.random.Generator, n: int, budget: float) -> np.ndarray:
    if budget == 0.0:
        return np.zeros(n, dtype=complex)
    raw = rng.normal(size=2 * n)
    norm = float(np.linalg.norm(raw))
    if norm == 0.0:
        return np.zeros(n, dtype=complex)
    # Uniform radius in the real 2n-ball.
    radius = budget * float(rng.random() ** (1.0 / (2 * n)))
    vals = radius * raw / norm
    return vals[:n] + 1j * vals[n:]


def coefficients_for(mode: str, raw: np.ndarray) -> dict[str, complex]:
    zero = 0.0 + 0.0j
    if mode == "down_only":
        return {
            "up": zero,
            "down": raw[0],
            "charged_lepton": zero,
            "neutrino_dirac": zero,
        }
    if mode == "up_down":
        return {
            "up": raw[0],
            "down": raw[1],
            "charged_lepton": zero,
            "neutrino_dirac": zero,
        }
    if mode == "down_lepton_tied":
        return {
            "up": zero,
            "down": raw[0],
            "charged_lepton": raw[0],
            "neutrino_dirac": raw[1],
        }
    if mode == "all_independent":
        return {
            "up": raw[0],
            "down": raw[1],
            "charged_lepton": raw[2],
            "neutrino_dirac": raw[3],
        }
    raise ValueError(mode)


def n_coeff(mode: str) -> int:
    return {
        "down_only": 1,
        "up_down": 2,
        "down_lepton_tied": 2,
        "all_independent": 4,
    }[mode]


def deform(seed: dict[str, np.ndarray], coeffs: dict[str, complex], k: np.ndarray) -> dict[str, np.ndarray]:
    return {name: mat + coeffs[name] * k for name, mat in seed.items()}


def scores(y: dict[str, np.ndarray]) -> dict[str, float]:
    _v, ckm, _mass_quark = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(
        fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2
        for key in ["Vus", "Vcb", "Vub"]
    )
    j_score = fit.log10_ratio(ckm["J"], fit.TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    return {
        "Vus": ckm["Vus"],
        "Vcb": ckm["Vcb"],
        "Vub": ckm["Vub"],
        "J": ckm["J"],
        "ckm_score": float(ckm_score),
        "jarlskog_score": float(j_score),
        "mass_score": float(mass_score),
    }


def left_right_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, s, vh = np.linalg.svd(y, full_matrices=True)
    order = np.argsort(s)
    s = s[order]
    left = u[:, order]
    right = vh.conjugate().T[:, order]
    diag = left.conjugate().T @ y @ right
    phases = np.array(
        [
            1.0 if abs(diag[i, i]) < 1.0e-30 else diag[i, i] / abs(diag[i, i])
            for i in range(3)
        ],
        dtype=complex,
    )
    right = right @ np.diag(phases.conjugate())
    return left, right, s


def knu_proxy(y: dict[str, np.ndarray]) -> dict[str, Any]:
    """Return the dangerous Knu proxy amplitude after item-4 normalization."""
    y_phys = {
        "up": d5.scale_to_largest_singular(y["up"], 0.60),
        "down": d5.scale_to_largest_singular(y["down"], 0.024),
        "charged_lepton": d5.scale_to_largest_singular(y["charged_lepton"], 0.010),
        "neutrino_dirac": d5.scale_to_largest_singular(y["neutrino_dirac"], 0.35),
    }
    ul: dict[str, np.ndarray] = {}
    ur: dict[str, np.ndarray] = {}
    singulars: dict[str, list[float]] = {}
    for name, mat in y.items():
        left, right, s = left_right_rotation(mat)
        ul[name] = left
        ur[name] = right
        singulars[name] = [float(v) for v in s]

    # Use the same local PMNS target convention as the transvectant export: the
    # neutrino mass-basis left rotation is U_e U_PMNS.
    import verify_seesaw_item3 as seesaw  # local import avoids a global dependency cycle

    benchmark = seesaw.LightNuBenchmark(
        m1_eV=1.0e-3,
        dm21_eV2=7.42e-5,
        dm31_eV2=2.517e-3,
        sin2_theta12=0.304,
        sin2_theta13=0.0222,
        sin2_theta23=0.573,
        delta_cp_rad=1.20 * math.pi,
        alpha21_rad=0.35 * math.pi,
        alpha31_rad=1.10 * math.pi,
    )
    u_pmns = seesaw.standard_pmns(
        benchmark.sin2_theta12,
        benchmark.sin2_theta13,
        benchmark.sin2_theta23,
        benchmark.delta_cp_rad,
        benchmark.alpha21_rad,
        benchmark.alpha31_rad,
    )
    u_nu = ul["charged_lepton"] @ u_pmns

    y_qq = d5.sym(y_phys["up"])
    y_ql = y_phys["down"]
    cl_upupdown_nu = d5.tensor_llll(
        y_qq,
        y_ql,
        ul["up"],
        ul["up"],
        ul["down"],
        u_nu,
    )
    amp, idx = d5.max_perm_entry(cl_upupdown_nu, (0, 0, 1), None)
    return {
        "amplitude": float(amp),
        "selected_index": "".join(str(i) for i in idx),
        "singular_values": singulars,
    }


def row_from(
    *,
    label: str,
    mode: str,
    budget: float,
    coeffs: dict[str, complex],
    y: dict[str, np.ndarray],
    baseline: dict[str, Any],
    target_future_margin: float,
) -> dict[str, Any]:
    s = scores(y)
    kproxy = knu_proxy(y)
    amp = float(kproxy["amplitude"])
    amp_ratio = amp / baseline["knu_amplitude"]
    inferred_future_margin = target_future_margin / max(amp_ratio * amp_ratio, 1.0e-30)
    coeff_norm = math.sqrt(sum(abs(v) ** 2 for v in coeffs.values()))
    return {
        "label": label,
        "mode": mode,
        "budget": budget,
        "coeff_norm": coeff_norm,
        "eps_up": coeffs["up"],
        "eps_down": coeffs["down"],
        "eps_charged_lepton": coeffs["charged_lepton"],
        "eps_neutrino_dirac": coeffs["neutrino_dirac"],
        **s,
        "knu_amplitude": amp,
        "knu_selected_index": kproxy["selected_index"],
        "knu_amplitude_ratio_to_baseline": amp_ratio,
        "inferred_knu_future_margin_1e35": inferred_future_margin,
        "inferred_suppression_needed": max(1.0 / math.sqrt(inferred_future_margin), 1.0),
        "ckm_improvement_factor": baseline["ckm_score"] / max(s["ckm_score"], 1.0e-30),
        "knu_amplitude_reduction_factor": baseline["knu_amplitude"] / max(amp, 1.0e-30),
        "mass_score_ratio": s["mass_score"] / max(baseline["mass_score"], 1.0e-30),
    }


def csv_row(row: dict[str, Any]) -> dict[str, Any]:
    out = dict(row)
    for key in ["eps_up", "eps_down", "eps_charged_lepton", "eps_neutrino_dirac"]:
        z = out.pop(key)
        out[key + "_re"] = float(np.real(z))
        out[key + "_im"] = float(np.imag(z))
        out[key + "_abs"] = float(abs(z))
    return out


def matrix_payload(y: dict[str, np.ndarray]) -> dict[str, Any]:
    return {name: fit.matrix_json(mat) for name, mat in y.items()}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    seed = load_seed()
    k = analytic_k()
    target = read_json(KNU_TARGET)
    target_future_margin = float(target["verdict"]["future_margin"])

    baseline_y = seed
    baseline_scores = scores(baseline_y)
    baseline_knu = knu_proxy(baseline_y)
    baseline = {
        **baseline_scores,
        "knu_amplitude": float(baseline_knu["amplitude"]),
    }
    rows: list[dict[str, Any]] = [
        row_from(
            label="baseline_transvectant",
            mode="baseline",
            budget=0.0,
            coeffs={name: 0.0 + 0.0j for name in seed},
            y=baseline_y,
            baseline=baseline,
            target_future_margin=target_future_margin,
        )
    ]
    y_cache: dict[str, dict[str, np.ndarray]] = {"baseline_transvectant": baseline_y}

    for mode in MODES:
        n = n_coeff(mode)
        for budget in BUDGETS:
            samples = 1 if budget == 0.0 else RANDOM_PER_MODE
            for sample in range(samples):
                raw = random_complex_vector(rng, n, budget)
                coeffs = coefficients_for(mode, raw)
                y = deform(seed, coeffs, k)
                label = f"{mode}_b{budget:g}_{sample}"
                row = row_from(
                    label=label,
                    mode=mode,
                    budget=budget,
                    coeffs=coeffs,
                    y=y,
                    baseline=baseline,
                    target_future_margin=target_future_margin,
                )
                rows.append(row)
                # Keep only potentially useful matrix payloads to avoid huge JSON.
                if (
                    row["ckm_score"] < baseline["ckm_score"]
                    or row["knu_amplitude"] < baseline["knu_amplitude"]
                    or row["inferred_knu_future_margin_1e35"] > 1.0
                ):
                    y_cache[label] = y

    # Pareto fronts under two physically distinct partial orders.
    best_ckm = min(rows, key=lambda r: r["ckm_score"])
    best_knu = min(rows, key=lambda r: r["knu_amplitude"])
    best_balanced = min(
        rows,
        key=lambda r: (
            r["ckm_score"] / max(baseline["ckm_score"], 1.0e-30)
            + 0.35 * r["mass_score"] / max(baseline["mass_score"], 1.0e-30)
            + max(0.0, 1.0 - r["inferred_knu_future_margin_1e35"])
        ),
    )
    joint_improvers = [
        r
        for r in rows
        if r["ckm_score"] < baseline["ckm_score"]
        and r["knu_amplitude"] < baseline["knu_amplitude"]
        and r["mass_score"] < 1.25 * baseline["mass_score"]
    ]
    future_safe = [
        r for r in rows if r["inferred_knu_future_margin_1e35"] >= 1.0
    ]
    future_safe_with_ckm_gain = [
        r
        for r in future_safe
        if r["ckm_score"] < baseline["ckm_score"]
        and r["mass_score"] < 1.50 * baseline["mass_score"]
    ]
    best_joint = (
        min(joint_improvers, key=lambda r: r["ckm_score"] + r["mass_score"])
        if joint_improvers
        else None
    )

    keep_labels = {
        "baseline_transvectant",
        best_ckm["label"],
        best_knu["label"],
        best_balanced["label"],
    }
    if best_joint is not None:
        keep_labels.add(best_joint["label"])
    if future_safe_with_ckm_gain:
        keep_labels.add(
            min(future_safe_with_ckm_gain, key=lambda r: r["ckm_score"])["label"]
        )

    csv_rows = [csv_row(r) for r in sorted(rows, key=lambda r: r["ckm_score"])[:250]]
    fields = list(csv_rows[0].keys())
    with (OUT / "shadow_contact_joint_scan_top_ckm.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(csv_rows)

    summary_rows = []
    for name, row in [
        ("baseline", rows[0]),
        ("best_ckm", best_ckm),
        ("best_knu", best_knu),
        ("best_balanced", best_balanced),
        ("best_joint", best_joint),
    ]:
        if row is not None:
            summary_rows.append({"role": name, **csv_row(row)})
    with (OUT / "shadow_contact_joint_summary_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    matrix_cards = {
        label: matrix_payload(y_cache[label])
        for label in sorted(keep_labels)
        if label in y_cache
    }
    verdict = {
        "samples": len(rows),
        "baseline_ckm_score": baseline["ckm_score"],
        "baseline_mass_score": baseline["mass_score"],
        "baseline_knu_amplitude": baseline["knu_amplitude"],
        "baseline_future_margin": target_future_margin,
        "joint_improver_count": len(joint_improvers),
        "future_safe_count": len(future_safe),
        "future_safe_with_ckm_gain_count": len(future_safe_with_ckm_gain),
        "best_ckm": csv_row(best_ckm),
        "best_knu": csv_row(best_knu),
        "best_balanced": csv_row(best_balanced),
        "best_joint": csv_row(best_joint) if best_joint else None,
        "interpretation": (
            "The contact/shadow direction has a real joint gradient if "
            "joint_improver_count > 0: at least one small deformation improves "
            "both the CKM score and the Knu proxy while keeping the mass score "
            "within 25 percent of the transvectant baseline.  A point with "
            "future_safe_with_ckm_gain_count > 0 would close the local Knu "
            "future-stress proxy at the same time; otherwise the scan only "
            "identifies a useful Pareto direction, not a full solution."
        ),
    }

    output = {
        "note": "No web lookup used. Joint shadow/contact flavor and Knu d=5 Pareto scan.",
        "deformation": {
            "formula": "Y_a -> Y_a + eps_a K",
            "K": fit.matrix_json(k),
            "modes": MODES,
            "budgets": BUDGETS,
            "random_per_mode": RANDOM_PER_MODE,
            "rng_seed": RNG_SEED,
        },
        "target_reference": {
            "knu_target_map": str(KNU_TARGET),
            "baseline_future_margin": target_future_margin,
            "display_S_T": DISPLAY_ST,
        },
        "matrix_cards": matrix_cards,
        "verdict": verdict,
    }
    (OUT / "summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Shadow/contact joint flavor and d=5 scan",
        "",
        "No web lookup was used.",
        "",
        "The audited deformation is",
        "",
        "```text",
        "Y_a -> Y_a + eps_a K,",
        "K = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]].",
        "```",
        "",
        "| role | mode | budget | CKM score | mass score | Knu amp ratio | inferred future margin |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        report.append(
            f"| `{row['role']}` | `{row['mode']}` | {row['budget']:.3g} | "
            f"{row['ckm_score']:.6e} | {row['mass_score']:.6e} | "
            f"{row['knu_amplitude_ratio_to_baseline']:.6e} | "
            f"{row['inferred_knu_future_margin_1e35']:.6e} |"
        )
    report += [
        "",
        "## Verdict",
        "",
        f"Joint improvers: `{len(joint_improvers)}`.",
        f"Future-safe points: `{len(future_safe)}`.",
        f"Future-safe with CKM gain: `{len(future_safe_with_ckm_gain)}`.",
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Shadow/contact joint flavor+d5 scan")
    print(f"  samples: {len(rows)}")
    print(f"  baseline CKM score: {baseline['ckm_score']:.6e}")
    print(f"  best CKM score: {best_ckm['ckm_score']:.6e}")
    print(f"  best Knu amp ratio: {best_knu['knu_amplitude_ratio_to_baseline']:.6e}")
    print(f"  joint improvers: {len(joint_improvers)}")
    print(f"  future-safe with CKM gain: {len(future_safe_with_ckm_gain)}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
