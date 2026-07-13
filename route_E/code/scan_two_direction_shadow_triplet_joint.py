#!/usr/bin/env python3
"""Two-direction shadow/contact plus triplet-alignment scan.

No web lookup is used.  This is the next audit after
scan_shadow_contact_flavor_d5_joint.py.  The previous one-line transvectant
deformation found a real CKM--Knu Pareto direction, but no point was both
CKM-improved and future-safe for the 1e35 yr Knu stress test.

Here we separate two logically different deformations:

  1. a doublet/flavor shadow direction

         Y_a -> Y_a + eps_a K,

     where K is the CP1 O(2) spin-zero second transvectant;

  2. a triplet-only colored-Higgs direction

         Y_QL^T -> Y_QL^T + delta_T,

     used only in the dimension-five Knu Wilson tensor.  This tests whether a
     locked-link/triplet alignment can protect proton decay while leaving the
     successful doublet flavor fit approximately unchanged.

Two triplet directions are audited:

  * contact_K: the same K direction, scaled to the down-type triplet size;
  * linear_cancel: the minimum-norm direction that cancels the currently
    dominant Knu linear functional.  This is a diagnostic/tuned fallback, not
    by itself a symmetry derivation.

The output states whether any sampled point satisfies all local targets:
CKM score < 1e-3, mass score within 25 percent of the baseline, and inferred
future Knu margin > 1.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_dimension5_wilson_tensors as d5  # noqa: E402


TRANSVECTANT = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
KNU_TARGET = ROOT / "output" / "knu_target_map" / "summary.json"
OUT = ROOT / "output" / "two_direction_shadow_triplet_joint"

RNG_SEED = 202605082131
FLAVOR_RANDOM_PER_MODE = 1800
FLAVOR_BUDGETS = [0.0, 0.015, 0.03, 0.06, 0.10, 0.16]
FLAVOR_MODES = ["down_only", "up_down", "down_lepton_tied", "all_independent"]
KEEP_TOP = 260
CONTACT_PHASES = 12
CONTACT_STRENGTHS = [0.0, 0.20, 0.40, 0.65, 0.90, 1.20, 1.60]
LINEAR_CANCEL_FRACTIONS = [0.0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50]
CKM_TARGET_SCORE = 1.0e-3
MASS_TOLERANCE = 1.25


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
    return {name: fit.cmat(raw) for name, raw in payload["Yukawa_matrices"].items()}


def random_complex_vector(rng: np.random.Generator, n: int, budget: float) -> np.ndarray:
    if budget == 0.0:
        return np.zeros(n, dtype=complex)
    raw = rng.normal(size=2 * n)
    norm = float(np.linalg.norm(raw))
    if norm == 0.0:
        return np.zeros(n, dtype=complex)
    radius = budget * float(rng.random() ** (1.0 / (2 * n)))
    vals = radius * raw / norm
    return vals[:n] + 1j * vals[n:]


def n_coeff(mode: str) -> int:
    return {
        "down_only": 1,
        "up_down": 2,
        "down_lepton_tied": 2,
        "all_independent": 4,
    }[mode]


def coefficients_for(mode: str, raw: np.ndarray) -> dict[str, complex]:
    zero = 0.0 + 0.0j
    if mode == "down_only":
        return {"up": zero, "down": raw[0], "charged_lepton": zero, "neutrino_dirac": zero}
    if mode == "up_down":
        return {"up": raw[0], "down": raw[1], "charged_lepton": zero, "neutrino_dirac": zero}
    if mode == "down_lepton_tied":
        return {"up": zero, "down": raw[0], "charged_lepton": raw[0], "neutrino_dirac": raw[1]}
    if mode == "all_independent":
        return {"up": raw[0], "down": raw[1], "charged_lepton": raw[2], "neutrino_dirac": raw[3]}
    raise ValueError(mode)


def deform(seed: dict[str, np.ndarray], coeffs: dict[str, complex], k: np.ndarray) -> dict[str, np.ndarray]:
    return {name: mat + coeffs[name] * k for name, mat in seed.items()}


def coeff_norm(coeffs: dict[str, complex]) -> float:
    return math.sqrt(sum(abs(v) ** 2 for v in coeffs.values()))


def left_right_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, s, vh = np.linalg.svd(y, full_matrices=True)
    order = np.argsort(s)
    s = s[order]
    left = u[:, order]
    right = vh.conjugate().T[:, order]
    diag = left.conjugate().T @ y @ right
    phases = np.array(
        [1.0 if abs(diag[i, i]) < 1.0e-30 else diag[i, i] / abs(diag[i, i]) for i in range(3)],
        dtype=complex,
    )
    right = right @ np.diag(phases.conjugate())
    return left, right, s


def scores(y: dict[str, np.ndarray]) -> dict[str, float]:
    _v, ckm, _mass_quark = fit.ckm_matrix(y["up"], y["down"])
    masses = fit.all_mass_ratios(y)
    ckm_score = sum(fit.log10_ratio(ckm[key], fit.TARGETS[key]) ** 2 for key in ["Vus", "Vcb", "Vub"])
    j_score = fit.log10_ratio(ckm["J"], fit.TARGETS["J"]) ** 2
    mass_score = 0.0
    for sector, (small, mid) in masses.items():
        mass_score += fit.log10_ratio(small, fit.TARGETS[f"{sector}_small"]) ** 2
        mass_score += fit.log10_ratio(mid, fit.TARGETS[f"{sector}_mid"]) ** 2
    return {
        "Vus": float(ckm["Vus"]),
        "Vcb": float(ckm["Vcb"]),
        "Vub": float(ckm["Vub"]),
        "J": float(ckm["J"]),
        "ckm_score": float(ckm_score),
        "jarlskog_score": float(j_score),
        "mass_score": float(mass_score),
    }


def physical_yukawas(y: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {
        "up": d5.scale_to_largest_singular(y["up"], 0.60),
        "down": d5.scale_to_largest_singular(y["down"], 0.024),
        "charged_lepton": d5.scale_to_largest_singular(y["charged_lepton"], 0.010),
        "neutrino_dirac": d5.scale_to_largest_singular(y["neutrino_dirac"], 0.35),
    }


def pmns_target(u_e: np.ndarray) -> np.ndarray:
    import verify_seesaw_item3 as seesaw

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
    return u_e @ u_pmns


@dataclass
class FlavorCard:
    label: str
    mode: str
    budget: float
    coeffs: dict[str, complex]
    y: dict[str, np.ndarray]
    y_phys: dict[str, np.ndarray]
    ul: dict[str, np.ndarray]
    ur: dict[str, np.ndarray]
    u_nu: np.ndarray
    score: dict[str, float]


def build_card(label: str, mode: str, budget: float, coeffs: dict[str, complex], y: dict[str, np.ndarray]) -> FlavorCard:
    ul: dict[str, np.ndarray] = {}
    ur: dict[str, np.ndarray] = {}
    for name, mat in y.items():
        left, right, _s = left_right_rotation(mat)
        ul[name] = left
        ur[name] = right
    return FlavorCard(
        label=label,
        mode=mode,
        budget=budget,
        coeffs=coeffs,
        y=y,
        y_phys=physical_yukawas(y),
        ul=ul,
        ur=ur,
        u_nu=pmns_target(ul["charged_lepton"]),
        score=scores(y),
    )


def knu_amplitude(card: FlavorCard, y_ql_triplet: np.ndarray) -> tuple[float, tuple[int, int, int, int], np.ndarray]:
    y_qq = d5.sym(card.y_phys["up"])
    tensor = d5.tensor_llll(
        y_qq,
        y_ql_triplet,
        card.ul["up"],
        card.ul["up"],
        card.ul["down"],
        card.u_nu,
    )
    amp, idx = d5.max_perm_entry(tensor, (0, 0, 1), None)
    return float(amp), idx, tensor


def linear_functional_for_idx(card: FlavorCard, idx: tuple[int, int, int, int]) -> tuple[complex, np.ndarray]:
    """Return A and L where A=sum_kl Y_QL[kl] L[kl] for a selected Knu index."""
    y_qq = d5.sym(card.y_phys["up"])
    a, b, c, ell = idx
    pref = np.einsum("ij,i,j->", y_qq, card.ul["up"][:, a], card.ul["up"][:, b], optimize=True)
    linear = pref * np.outer(card.ul["down"][:, c], card.u_nu[:, ell])
    y_ql = card.y_phys["down"]
    amp_complex = np.einsum("kl,kl->", y_ql, linear, optimize=True)
    return complex(amp_complex), linear


def contact_rows_for_card(
    card: FlavorCard,
    baseline_amp: float,
    target_future_margin: float,
    k_triplet: np.ndarray,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    y_base = card.y_phys["down"]
    phase_grid = [np.exp(2j * math.pi * p / CONTACT_PHASES) for p in range(CONTACT_PHASES)]
    for strength in CONTACT_STRENGTHS:
        phases = [1.0 + 0.0j] if strength == 0.0 else phase_grid
        for phase in phases:
            delta = strength * phase * k_triplet
            rows.append(row_from_triplet(card, y_base + delta, "contact_K", strength, phase, delta, baseline_amp, target_future_margin))
    return rows


def linear_cancel_rows_for_card(card: FlavorCard, baseline_amp: float, target_future_margin: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    y_base = card.y_phys["down"]
    _amp, idx, _tensor = knu_amplitude(card, y_base)
    amp_complex, linear = linear_functional_for_idx(card, idx)
    norm = float(np.linalg.norm(linear))
    if norm <= 1.0e-30:
        return []
    direction = np.conjugate(linear) / norm
    zeta_cancel = -amp_complex / norm
    for frac in LINEAR_CANCEL_FRACTIONS:
        delta = frac * zeta_cancel * direction
        rows.append(row_from_triplet(card, y_base + delta, "linear_cancel", frac, zeta_cancel, delta, baseline_amp, target_future_margin))
    return rows


def row_from_triplet(
    card: FlavorCard,
    y_ql_triplet: np.ndarray,
    triplet_mode: str,
    triplet_strength: float,
    triplet_phase_or_zeta: complex,
    delta: np.ndarray,
    baseline_amp: float,
    target_future_margin: float,
) -> dict[str, Any]:
    amp, idx, _tensor = knu_amplitude(card, y_ql_triplet)
    amp_ratio = amp / max(baseline_amp, 1.0e-30)
    future_margin = target_future_margin / max(amp_ratio * amp_ratio, 1.0e-30)
    row = {
        "label": card.label,
        "flavor_mode": card.mode,
        "flavor_budget": card.budget,
        "flavor_coeff_norm": coeff_norm(card.coeffs),
        "triplet_mode": triplet_mode,
        "triplet_strength": float(triplet_strength),
        "triplet_phase_re": float(np.real(triplet_phase_or_zeta)),
        "triplet_phase_im": float(np.imag(triplet_phase_or_zeta)),
        "triplet_delta_norm": float(np.linalg.norm(delta)),
        "triplet_delta_over_down_norm": float(np.linalg.norm(delta) / max(np.linalg.norm(card.y_phys["down"]), 1.0e-30)),
        "knu_amplitude": amp,
        "knu_selected_index": "".join(str(i) for i in idx),
        "knu_amplitude_ratio_to_baseline": amp_ratio,
        "inferred_knu_future_margin_1e35": future_margin,
        "inferred_suppression_needed": max(1.0 / math.sqrt(future_margin), 1.0),
        **card.score,
    }
    row["passes_local_closure"] = bool(
        row["ckm_score"] < CKM_TARGET_SCORE
        and row["mass_score"] <= MASS_TOLERANCE * BASELINE_MASS_SCORE
        and row["inferred_knu_future_margin_1e35"] >= 1.0
    )
    row["ckm_improvement"] = float(BASELINE_CKM_SCORE / max(row["ckm_score"], 1.0e-30))
    row["mass_score_ratio"] = float(row["mass_score"] / max(BASELINE_MASS_SCORE, 1.0e-30))
    row["knu_reduction"] = float(baseline_amp / max(amp, 1.0e-30))
    return row


def csv_safe(row: dict[str, Any]) -> dict[str, Any]:
    return {k: (float(v) if isinstance(v, np.floating) else v) for k, v in row.items()}


def matrix_payload(card: FlavorCard) -> dict[str, Any]:
    return {
        "Yukawa_matrices": {name: fit.matrix_json(mat) for name, mat in card.y.items()},
        "coefficients": {name: cjson(z) for name, z in card.coeffs.items()},
    }


# These are set in main and used only to keep row construction compact.
BASELINE_CKM_SCORE = math.nan
BASELINE_MASS_SCORE = math.nan


def main() -> None:
    global BASELINE_CKM_SCORE, BASELINE_MASS_SCORE

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    seed = load_seed()
    k = analytic_k()
    target_future_margin = float(read_json(KNU_TARGET)["verdict"]["future_margin"])

    zero_coeffs = {name: 0.0 + 0.0j for name in seed}
    baseline_card = build_card("baseline_transvectant", "baseline", 0.0, zero_coeffs, seed)
    baseline_amp, baseline_idx, _tensor = knu_amplitude(baseline_card, baseline_card.y_phys["down"])
    BASELINE_CKM_SCORE = baseline_card.score["ckm_score"]
    BASELINE_MASS_SCORE = baseline_card.score["mass_score"]

    flavor_cards = [baseline_card]
    for mode in FLAVOR_MODES:
        n = n_coeff(mode)
        for budget in FLAVOR_BUDGETS:
            samples = 1 if budget == 0.0 else FLAVOR_RANDOM_PER_MODE
            for sample in range(samples):
                raw = random_complex_vector(rng, n, budget)
                coeffs = coefficients_for(mode, raw)
                y = deform(seed, coeffs, k)
                flavor_cards.append(build_card(f"{mode}_b{budget:g}_{sample}", mode, budget, coeffs, y))

    def balanced_key(card: FlavorCard) -> float:
        return (
            card.score["ckm_score"] / max(BASELINE_CKM_SCORE, 1.0e-30)
            + 0.35 * card.score["mass_score"] / max(BASELINE_MASS_SCORE, 1.0e-30)
        )

    keep: dict[str, FlavorCard] = {baseline_card.label: baseline_card}
    for card in sorted(flavor_cards, key=lambda c: c.score["ckm_score"])[:KEEP_TOP]:
        keep[card.label] = card
    for card in sorted(flavor_cards, key=balanced_key)[:KEEP_TOP]:
        keep[card.label] = card
    for card in flavor_cards:
        if card.score["ckm_score"] < BASELINE_CKM_SCORE and card.score["mass_score"] <= 1.25 * BASELINE_MASS_SCORE:
            keep[card.label] = card
    kept_cards = list(keep.values())

    k_triplet = d5.scale_to_largest_singular(k, 0.024)
    triplet_rows: list[dict[str, Any]] = []
    for card in kept_cards:
        triplet_rows.extend(contact_rows_for_card(card, baseline_amp, target_future_margin, k_triplet))
        triplet_rows.extend(linear_cancel_rows_for_card(card, baseline_amp, target_future_margin))

    best_ckm = min(triplet_rows, key=lambda r: r["ckm_score"])
    best_knu = min(triplet_rows, key=lambda r: r["knu_amplitude"])
    best_balanced = min(
        triplet_rows,
        key=lambda r: (
            r["ckm_score"] / max(CKM_TARGET_SCORE, 1.0e-30)
            + 0.25 * r["mass_score"] / max(BASELINE_MASS_SCORE, 1.0e-30)
            + max(0.0, 1.0 - r["inferred_knu_future_margin_1e35"])
            + 0.10 * r["triplet_delta_over_down_norm"]
        ),
    )
    closure_rows = [r for r in triplet_rows if r["passes_local_closure"]]
    natural_closure_rows = [r for r in closure_rows if r["triplet_delta_over_down_norm"] <= 1.0]
    contact_closure_rows = [r for r in closure_rows if r["triplet_mode"] == "contact_K"]
    linear_closure_rows = [r for r in closure_rows if r["triplet_mode"] == "linear_cancel"]
    ckm_future_rows = [
        r
        for r in triplet_rows
        if r["ckm_score"] < CKM_TARGET_SCORE and r["inferred_knu_future_margin_1e35"] >= 1.0
    ]
    best_closure = min(closure_rows, key=lambda r: r["triplet_delta_over_down_norm"]) if closure_rows else None

    summary_rows = [
        {"role": "baseline", **triplet_rows[0]},
        {"role": "best_ckm", **best_ckm},
        {"role": "best_knu", **best_knu},
        {"role": "best_balanced", **best_balanced},
    ]
    if best_closure is not None:
        summary_rows.append({"role": "best_closure_min_delta", **best_closure})

    top_rows = sorted(
        triplet_rows,
        key=lambda r: (
            0 if r["passes_local_closure"] else 1,
            r["ckm_score"] / max(CKM_TARGET_SCORE, 1.0e-30)
            + max(0.0, 1.0 - r["inferred_knu_future_margin_1e35"])
            + 0.10 * r["mass_score"] / max(BASELINE_MASS_SCORE, 1.0e-30),
        ),
    )[:500]
    with (OUT / "two_direction_scan_top.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_safe(top_rows[0]).keys()))
        writer.writeheader()
        writer.writerows(csv_safe(r) for r in top_rows)

    with (OUT / "two_direction_summary_rows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_safe(summary_rows[0]).keys()))
        writer.writeheader()
        writer.writerows(csv_safe(r) for r in summary_rows)

    keep_labels = {r["label"] for r in summary_rows}
    if best_closure is not None:
        keep_labels.add(best_closure["label"])
    card_by_label = {card.label: card for card in kept_cards}
    matrix_cards = {label: matrix_payload(card_by_label[label]) for label in sorted(keep_labels) if label in card_by_label}

    verdict = {
        "flavor_samples": len(flavor_cards),
        "kept_flavor_cards": len(kept_cards),
        "triplet_rows": len(triplet_rows),
        "baseline_ckm_score": BASELINE_CKM_SCORE,
        "baseline_mass_score": BASELINE_MASS_SCORE,
        "baseline_knu_amplitude": baseline_amp,
        "baseline_knu_index": "".join(str(i) for i in baseline_idx),
        "baseline_future_margin": target_future_margin,
        "local_closure_count": len(closure_rows),
        "natural_local_closure_count": len(natural_closure_rows),
        "contact_closure_count": len(contact_closure_rows),
        "linear_cancel_closure_count": len(linear_closure_rows),
        "ckm_future_rows_count": len(ckm_future_rows),
        "best_ckm": csv_safe(best_ckm),
        "best_knu": csv_safe(best_knu),
        "best_balanced": csv_safe(best_balanced),
        "best_closure_min_delta": csv_safe(best_closure) if best_closure else None,
        "interpretation": (
            "A nonzero local_closure_count means the two-direction completion "
            "can satisfy the local CKM, mass, and future-Knu proxy targets in "
            "this reduced audit.  contact_closure_count is the conservative "
            "symmetry-like result.  linear_cancel_closure_count is a diagnostic: "
            "it proves that a triplet-only functional direction exists, but it "
            "must still be derived from a symmetry or mediator Hessian before it "
            "can close the theory."
        ),
    }

    payload = {
        "note": "No web lookup used. Two-direction shadow/contact flavor and triplet-only d=5 scan.",
        "mathematical_setup": {
            "flavor_deformation": "Y_a -> Y_a + eps_a K",
            "triplet_deformation": "Y_QL^T -> Y_QL^T + delta_T",
            "K": fit.matrix_json(k),
            "linear_cancel_derivation": (
                "For a fixed dangerous Knu index A=sum_kl Y_QL[kl] L[kl]. "
                "The minimum-Frobenius-norm direction changing only this "
                "linear functional is D=conj(L)/||L||, with zeta=-A/||L||."
            ),
        },
        "scan_config": {
            "rng_seed": RNG_SEED,
            "flavor_random_per_mode": FLAVOR_RANDOM_PER_MODE,
            "flavor_budgets": FLAVOR_BUDGETS,
            "flavor_modes": FLAVOR_MODES,
            "keep_top": KEEP_TOP,
            "contact_phases": CONTACT_PHASES,
            "contact_strengths": CONTACT_STRENGTHS,
            "linear_cancel_fractions": LINEAR_CANCEL_FRACTIONS,
            "ckm_target_score": CKM_TARGET_SCORE,
            "mass_tolerance": MASS_TOLERANCE,
        },
        "matrix_cards": matrix_cards,
        "verdict": verdict,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Two-direction shadow/triplet joint scan",
        "",
        "No web lookup was used.",
        "",
        "The scan separates the doublet flavor shadow direction from a triplet-only colored-Higgs direction:",
        "",
        "```text",
        "Y_a -> Y_a + eps_a K,",
        "Y_QL^T -> Y_QL^T + delta_T.",
        "```",
        "",
        "The diagnostic linear-cancel direction is derived from the local linear functional",
        "`A=sum_kl Y_QL[kl] L[kl]`; the minimum-norm cancellation direction is",
        "`D=conj(L)/||L||`, `zeta=-A/||L||`.",
        "",
        "| role | flavor mode | triplet mode | CKM score | mass score | Knu amp ratio | future margin | delta/down | closure |",
        "|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        report.append(
            f"| `{row['role']}` | `{row['flavor_mode']}` | `{row['triplet_mode']}` | "
            f"{row['ckm_score']:.6e} | {row['mass_score']:.6e} | "
            f"{row['knu_amplitude_ratio_to_baseline']:.6e} | "
            f"{row['inferred_knu_future_margin_1e35']:.6e} | "
            f"{row['triplet_delta_over_down_norm']:.6e} | `{row['passes_local_closure']}` |"
        )
    report += [
        "",
        "## Verdict",
        "",
        f"Flavor samples: `{len(flavor_cards)}`.",
        f"Kept flavor cards: `{len(kept_cards)}`.",
        f"Triplet rows: `{len(triplet_rows)}`.",
        f"Local closures: `{len(closure_rows)}`.",
        f"Natural closures with `delta/down <= 1`: `{len(natural_closure_rows)}`.",
        f"Contact-K closures: `{len(contact_closure_rows)}`.",
        f"Linear-cancel closures: `{len(linear_closure_rows)}`.",
        "",
        verdict["interpretation"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(report), encoding="utf-8")

    print("Two-direction shadow/triplet joint scan")
    print(f"  flavor samples: {len(flavor_cards)}")
    print(f"  kept flavor cards: {len(kept_cards)}")
    print(f"  triplet rows: {len(triplet_rows)}")
    print(f"  local closures: {len(closure_rows)}")
    print(f"  natural local closures: {len(natural_closure_rows)}")
    print(f"  contact closures: {len(contact_closure_rows)}")
    print(f"  linear-cancel closures: {len(linear_closure_rows)}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
