#!/usr/bin/env python3
"""Scan how predictive the Majorana contact lift is.

No web lookup is used.  The previous audit showed that the source-consistent
inverse-seesaw Majorana matrix decomposes as

    M_R / M_* = M_V + zeta K,

where M_V is the five-complex-dimensional CP1/O(2) Veronese component and K is
the normalized spin-zero contact direction.  This script tests whether the
contact coefficient behaves like a harmless small correction or a rigid source
parameter by scanning

    M_R(s) = M_* (M_V + s zeta K)

and a phase ring

    M_R(delta) = M_* (M_V + exp(i delta) zeta K).

For each point it computes the type-I seesaw light-neutrino matrix, Takagi
mixing, PMNS angles in the same charged-lepton basis, light mass splittings,
and conditioning.  The scan is deliberately one-dimensional/two-dimensional in
the contact coefficient; it is not a new unconstrained fit.
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


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import audit_source_majorana_texture_rank as rank_audit  # noqa: E402
import verify_seesaw_item3 as seesaw  # noqa: E402


OUT = ROOT / "output" / "majorana_contact_sensitivity"
CLOSURE_CARD = ROOT / "output" / "publication_closure_card" / "publication_closure_card.json"
RANK_SUMMARY = ROOT / "output" / "source_majorana_texture_rank" / "summary.json"


TARGET_ANGLES = {
    "sin2_theta12": 0.304,
    "sin2_theta13": 0.0222,
    "sin2_theta23": 0.573,
}
TARGET_DM = {
    "dm21_eV2": 7.42e-5,
    "dm31_eV2": 2.517e-3,
}
LOOSE_ANGLE_TOL = 3.0e-2
LOOSE_DM_REL_TOL = 3.0e-1
TIGHT_ANGLE_TOL = 1.0e-3
TIGHT_DM_REL_TOL = 1.0e-2


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cnum(raw: dict[str, float]) -> complex:
    return complex(raw["re"], raw["im"])


def mixing_angles(u: np.ndarray) -> dict[str, float]:
    uabs2 = np.abs(u) ** 2
    s13 = float(uabs2[0, 2])
    denom = max(1.0 - s13, 1.0e-30)
    return {
        "sin2_theta12": float(uabs2[0, 1] / denom),
        "sin2_theta13": s13,
        "sin2_theta23": float(uabs2[1, 2] / denom),
    }


def setup() -> dict[str, Any]:
    card = read_json(CLOSURE_CARD)
    rank_summary = read_json(RANK_SUMMARY)
    y = {name: rank_audit.cmat(raw) for name, raw in card["selected_row"]["Yukawa_fit"].items()}
    u_e, _ = seesaw.left_rotation(y["charged_lepton"])
    d = rank_summary["decomposition"]
    return {
        "card": card,
        "rank_summary": rank_summary,
        "u_e": u_e,
        "mD_eV": y["neutrino_dirac"] * 100.0e9,
        "scale_GeV": float(d["scale_largest_singular_GeV"]),
        "veronese": cmat(d["veronese_projection"]),
        "contact_basis": cmat(d["contact_basis"]),
        "zeta": cnum(d["contact_coefficient"]),
    }


def eval_coeff(data: dict[str, Any], coeff: complex) -> dict[str, Any]:
    mr_gev = data["scale_GeV"] * (data["veronese"] + coeff * data["contact_basis"])
    mr_ev = mr_gev * 1.0e9
    try:
        heavy = np.linalg.svd(mr_gev, compute_uv=False)
        cond = float(heavy[0] / max(heavy[-1], 1.0e-300))
        inv_mr = np.linalg.inv(mr_ev)
        m_light = -(data["mD_eV"] @ inv_mr @ data["mD_eV"].T)
        masses, u_nu = seesaw.takagi_by_h(m_light)
        u_pmns = data["u_e"].conjugate().T @ u_nu
        angles = mixing_angles(u_pmns)
        dm = {
            "dm21_eV2": float(masses[1] ** 2 - masses[0] ** 2),
            "dm31_eV2": float(masses[2] ** 2 - masses[0] ** 2),
        }
        angle_diffs = {key: angles[key] - TARGET_ANGLES[key] for key in TARGET_ANGLES}
        dm_rel = {
            key: (dm[key] - TARGET_DM[key]) / TARGET_DM[key]
            for key in TARGET_DM
        }
        angle_max = max(abs(x) for x in angle_diffs.values())
        dm_rel_max = max(abs(x) for x in dm_rel.values())
        # A dimensionless local score used only for ordering scan points.
        score = math.sqrt(
            sum((angle_diffs[key] / max(TARGET_ANGLES[key], 1.0e-30)) ** 2 for key in TARGET_ANGLES)
            + sum(dm_rel[key] ** 2 for key in TARGET_DM)
        )
        theta_norm = float(np.linalg.norm(data["mD_eV"] @ inv_mr, ord=2))
        valid = True
        error = ""
    except np.linalg.LinAlgError as exc:
        heavy = np.array([math.nan, math.nan, math.nan])
        cond = math.inf
        masses = np.array([math.nan, math.nan, math.nan])
        angles = {key: math.nan for key in TARGET_ANGLES}
        dm = {key: math.nan for key in TARGET_DM}
        angle_diffs = {key: math.nan for key in TARGET_ANGLES}
        dm_rel = {key: math.nan for key in TARGET_DM}
        angle_max = math.inf
        dm_rel_max = math.inf
        score = math.inf
        theta_norm = math.inf
        valid = False
        error = str(exc)
    return {
        "valid": valid,
        "error": error,
        "coeff_re": float(np.real(coeff)),
        "coeff_im": float(np.imag(coeff)),
        "coeff_abs": float(abs(coeff)),
        "heavy_min_GeV": float(np.min(heavy)),
        "heavy_mid_GeV": float(np.sort(heavy)[1]),
        "heavy_max_GeV": float(np.max(heavy)),
        "heavy_condition": cond,
        "theta_norm": theta_norm,
        "m1_eV": float(masses[0]),
        "m2_eV": float(masses[1]),
        "m3_eV": float(masses[2]),
        "dm21_eV2": dm["dm21_eV2"],
        "dm31_eV2": dm["dm31_eV2"],
        "dm21_rel_diff": dm_rel["dm21_eV2"],
        "dm31_rel_diff": dm_rel["dm31_eV2"],
        "sin2_theta12": angles["sin2_theta12"],
        "sin2_theta13": angles["sin2_theta13"],
        "sin2_theta23": angles["sin2_theta23"],
        "d_sin2_theta12": angle_diffs["sin2_theta12"],
        "d_sin2_theta13": angle_diffs["sin2_theta13"],
        "d_sin2_theta23": angle_diffs["sin2_theta23"],
        "angle_max_absdiff": float(angle_max),
        "dm_rel_max_absdiff": float(dm_rel_max),
        "score": float(score),
        "passes_loose": bool(angle_max <= LOOSE_ANGLE_TOL and dm_rel_max <= LOOSE_DM_REL_TOL),
        "passes_tight": bool(angle_max <= TIGHT_ANGLE_TOL and dm_rel_max <= TIGHT_DM_REL_TOL),
    }


def line_row(data: dict[str, Any], s: float) -> dict[str, Any]:
    row = eval_coeff(data, s * data["zeta"])
    row["mode"] = "real_scale"
    row["s"] = float(s)
    row["delta_phase_rad"] = 0.0
    return row


def phase_row(data: dict[str, Any], delta: float) -> dict[str, Any]:
    row = eval_coeff(data, data["zeta"] * np.exp(1j * delta))
    row["mode"] = "phase_ring"
    row["s"] = 1.0
    row["delta_phase_rad"] = float(delta)
    return row


def finite_derivatives(data: dict[str, Any]) -> dict[str, Any]:
    h = 1.0e-4
    plus = line_row(data, 1.0 + h)
    minus = line_row(data, 1.0 - h)
    derivs = {}
    for key in [
        "sin2_theta12",
        "sin2_theta13",
        "sin2_theta23",
        "dm21_eV2",
        "dm31_eV2",
        "m1_eV",
        "m2_eV",
        "m3_eV",
    ]:
        derivs[f"d_{key}_ds_at_1"] = float((plus[key] - minus[key]) / (2.0 * h))
    return derivs


def interval_around(rows: list[dict[str, Any]], center_key: str, predicate_key: str) -> dict[str, Any]:
    valid = [row for row in rows if bool(row[predicate_key])]
    if not valid:
        return {"exists": False, "min": None, "max": None, "width": 0.0}
    # Keep the connected component containing the target point closest to zero
    # phase or s=1, because distant disconnected passes are not a source lock.
    target = 1.0 if center_key == "s" else 0.0
    values = sorted(float(row[center_key]) for row in valid)
    step_candidates = [abs(values[i + 1] - values[i]) for i in range(len(values) - 1) if values[i + 1] > values[i]]
    step = min(step_candidates) if step_candidates else 0.0
    tol = max(1.5 * step, 1.0e-12)
    closest = min(values, key=lambda x: abs(x - target))
    component = [closest]
    changed = True
    while changed:
        changed = False
        lo, hi = min(component), max(component)
        for value in values:
            if value < lo and lo - value <= tol:
                component.append(value)
                changed = True
            if value > hi and value - hi <= tol:
                component.append(value)
                changed = True
    return {
        "exists": True,
        "min": float(min(component)),
        "max": float(max(component)),
        "width": float(max(component) - min(component)),
        "grid_step": float(step),
    }


def adaptive_interval(
    make_row: Any,
    center: float,
    lower_bound: float,
    upper_bound: float,
    predicate_key: str,
) -> dict[str, Any]:
    center_row = make_row(center)
    if not bool(center_row[predicate_key]):
        return {"exists": False, "min": None, "max": None, "width": 0.0, "method": "adaptive"}

    def upper_edge() -> float:
        step = 1.0e-9
        passed = center
        candidate = min(center + step, upper_bound)
        while candidate < upper_bound and bool(make_row(candidate)[predicate_key]):
            passed = candidate
            step *= 2.0
            candidate = min(center + step, upper_bound)
        if bool(make_row(candidate)[predicate_key]):
            return candidate
        failed = candidate
        for _ in range(80):
            mid = 0.5 * (passed + failed)
            if bool(make_row(mid)[predicate_key]):
                passed = mid
            else:
                failed = mid
        return passed

    def lower_edge() -> float:
        step = 1.0e-9
        passed = center
        candidate = max(center - step, lower_bound)
        while candidate > lower_bound and bool(make_row(candidate)[predicate_key]):
            passed = candidate
            step *= 2.0
            candidate = max(center - step, lower_bound)
        if bool(make_row(candidate)[predicate_key]):
            return candidate
        failed = candidate
        for _ in range(80):
            mid = 0.5 * (passed + failed)
            if bool(make_row(mid)[predicate_key]):
                passed = mid
            else:
                failed = mid
        return passed

    lo = lower_edge()
    hi = upper_edge()
    return {
        "exists": True,
        "min": float(lo),
        "max": float(hi),
        "width": float(hi - lo),
        "half_width": float(0.5 * (hi - lo)),
        "method": "adaptive",
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def report(summary: dict[str, Any]) -> str:
    baseline = summary["representative_rows"]["veronese_only_s0"]
    target = summary["representative_rows"]["target_s1"]
    best = summary["representative_rows"]["best_real_scale"]
    phase_best = summary["representative_rows"]["best_phase_ring"]
    derivs = summary["finite_derivatives"]
    return "\n".join(
        [
            "# Majorana Contact Sensitivity Scan",
            "",
            "No web lookup was used.  This scan keeps the source-consistent Veronese projection fixed and varies only the contact coefficient.",
            "",
            "## Real contact scale",
            "",
            f"- Veronese-only s=0 angle max residual: `{baseline['angle_max_absdiff']:.6e}`",
            f"- Veronese-only s=0 splitting max relative residual: `{baseline['dm_rel_max_absdiff']:.6e}`",
            f"- target s=1 angle max residual: `{target['angle_max_absdiff']:.6e}`",
            f"- target s=1 splitting max relative residual: `{target['dm_rel_max_absdiff']:.6e}`",
            f"- best real-scale s: `{best['s']:.6f}` with score `{best['score']:.6e}`",
            f"- loose interval around s=1: `{summary['real_scale_loose_interval']}`",
            f"- tight interval around s=1: `{summary['real_scale_tight_interval']}`",
            "",
            "## Contact phase ring",
            "",
            f"- best phase delta: `{phase_best['delta_phase_rad']:.6e}` rad with score `{phase_best['score']:.6e}`",
            f"- loose phase interval around zero: `{summary['phase_loose_interval']}`",
            f"- tight phase interval around zero: `{summary['phase_tight_interval']}`",
            "",
            "## Local derivatives at s=1",
            "",
            f"- d sin2(theta12)/ds: `{derivs['d_sin2_theta12_ds_at_1']:.6e}`",
            f"- d sin2(theta13)/ds: `{derivs['d_sin2_theta13_ds_at_1']:.6e}`",
            f"- d sin2(theta23)/ds: `{derivs['d_sin2_theta23_ds_at_1']:.6e}`",
            f"- d dm21/ds: `{derivs['d_dm21_eV2_ds_at_1']:.6e}` eV^2",
            f"- d dm31/ds: `{derivs['d_dm31_eV2_ds_at_1']:.6e}` eV^2",
            "",
            "## Verdict",
            "",
            summary["verdict"]["interpretation"],
            "",
        ]
    )


def build() -> dict[str, Any]:
    data = setup()
    real_scales = [i / 1000.0 for i in range(0, 2001)]
    phases = [-math.pi + (2.0 * math.pi) * i / 1440.0 for i in range(0, 1441)]
    real_rows = [line_row(data, s) for s in real_scales]
    phase_rows = [phase_row(data, delta) for delta in phases]
    best_real = min((row for row in real_rows if row["valid"]), key=lambda row: row["score"])
    best_phase = min((row for row in phase_rows if row["valid"]), key=lambda row: row["score"])
    baseline = real_rows[0]
    target = min(real_rows, key=lambda row: abs(row["s"] - 1.0))
    derivs = finite_derivatives(data)
    grid_loose_real = interval_around(real_rows, "s", "passes_loose")
    grid_tight_real = interval_around(real_rows, "s", "passes_tight")
    grid_loose_phase = interval_around(phase_rows, "delta_phase_rad", "passes_loose")
    grid_tight_phase = interval_around(phase_rows, "delta_phase_rad", "passes_tight")
    loose_real = adaptive_interval(lambda x: line_row(data, x), 1.0, 0.0, 2.0, "passes_loose")
    tight_real = adaptive_interval(lambda x: line_row(data, x), 1.0, 0.0, 2.0, "passes_tight")
    loose_phase = adaptive_interval(lambda x: phase_row(data, x), 0.0, -math.pi, math.pi, "passes_loose")
    tight_phase = adaptive_interval(lambda x: phase_row(data, x), 0.0, -math.pi, math.pi, "passes_tight")
    contact_fraction = read_json(RANK_SUMMARY)["decomposition"]["contact_fraction"]
    phase_tolerance_rad = loose_phase["width"] / 2.0 if loose_phase["exists"] else 0.0
    scale_tolerance = loose_real["width"] / 2.0 if loose_real["exists"] else 0.0
    verdict = {
        "veronese_only_fails_loose_pmns_dm": not baseline["passes_loose"],
        "contact_scale_must_be_locked": scale_tolerance < 0.1,
        "contact_phase_must_be_locked": phase_tolerance_rad < 0.1,
        "predictive_source_action_still_needed": True,
        "interpretation": (
            "The contact coefficient is not a decorative perturbation: removing it "
            "spoils the PMNS/mass-splitting target, while the target point is recovered "
            "only in a narrow neighborhood of the reconstructed complex contact "
            "coefficient.  A source action must therefore predict both the magnitude "
            "and phase of the contact lift, or the PMNS result remains an inverse-seesaw "
            "compatibility statement rather than a predictive theorem."
        ),
    }
    return {
        "note": "No web lookup used. Majorana contact sensitivity scan.",
        "input_manifest": [
            {
                "label": "publication_closure_card",
                "path": str(CLOSURE_CARD.relative_to(ROOT)),
                "sha256": sha256(CLOSURE_CARD),
            },
            {
                "label": "source_majorana_texture_rank",
                "path": str(RANK_SUMMARY.relative_to(ROOT)),
                "sha256": sha256(RANK_SUMMARY),
            },
        ],
        "scan_definition": {
            "real_scale": "M_R(s)=M_*(M_V+s zeta K), s in [0,2]",
            "phase_ring": "M_R(delta)=M_*(M_V+exp(i delta) zeta K), delta in [-pi,pi]",
            "loose_angle_tolerance": LOOSE_ANGLE_TOL,
            "loose_dm_relative_tolerance": LOOSE_DM_REL_TOL,
            "tight_angle_tolerance": TIGHT_ANGLE_TOL,
            "tight_dm_relative_tolerance": TIGHT_DM_REL_TOL,
            "contact_fraction": contact_fraction,
        },
        "representative_rows": {
            "veronese_only_s0": baseline,
            "target_s1": target,
            "best_real_scale": best_real,
            "best_phase_ring": best_phase,
        },
        "real_scale_loose_interval": loose_real,
        "real_scale_tight_interval": tight_real,
        "phase_loose_interval": loose_phase,
        "phase_tight_interval": tight_phase,
        "grid_resolution_intervals": {
            "real_scale_loose_interval": grid_loose_real,
            "real_scale_tight_interval": grid_tight_real,
            "phase_loose_interval": grid_loose_phase,
            "phase_tight_interval": grid_tight_phase,
        },
        "finite_derivatives": derivs,
        "verdict": verdict,
        "real_scale_rows": real_rows,
        "phase_rows": phase_rows,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    write_csv(OUT / "real_scale_scan.csv", summary["real_scale_rows"])
    write_csv(OUT / "phase_ring_scan.csv", summary["phase_rows"])
    light_summary = {
        key: value
        for key, value in summary.items()
        if key not in {"real_scale_rows", "phase_rows"}
    }
    (OUT / "summary.json").write_text(json.dumps(light_summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(light_summary), encoding="utf-8")
    print("Majorana contact sensitivity scan")
    print(f"  s=0 angle max: {summary['representative_rows']['veronese_only_s0']['angle_max_absdiff']:.6e}")
    print(f"  s=0 dm rel max: {summary['representative_rows']['veronese_only_s0']['dm_rel_max_absdiff']:.6e}")
    print(f"  best s: {summary['representative_rows']['best_real_scale']['s']:.6f}")
    print(f"  loose s interval: {summary['real_scale_loose_interval']}")
    print(f"  loose phase interval: {summary['phase_loose_interval']}")


if __name__ == "__main__":
    main()
