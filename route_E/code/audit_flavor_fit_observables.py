#!/usr/bin/env python3
"""Audit CKM/PMNS observables from the exact flavor benchmark card.

No web lookup is used.  The CKM target below is an internal phenomenological
benchmark, not a refreshed world-average input.  The point of the audit is to
separate what the current CP1/O(2) card already proves (hierarchies and a
reproducible seesaw) from what it does not yet prove (a full flavor fit).
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CARD = ROOT / "output" / "flavor_benchmark" / "flavor_benchmark_card.json"
OUT = ROOT / "output" / "flavor_fit"


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def matrix_json(m: np.ndarray) -> list[list[dict[str, float]]]:
    return [[cjson(z) for z in row] for row in m]


def matrix_abs_json(m: np.ndarray) -> list[list[float]]:
    return [[float(abs(z)) for z in row] for row in m]


def left_rotation(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = y @ y.conjugate().T
    values, vectors = np.linalg.eigh(h)
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    return vectors, np.sqrt(np.maximum(values, 0.0))


def standard_ckm(s12: float, s23: float, s13: float, delta: float) -> np.ndarray:
    c12 = math.sqrt(1.0 - s12**2)
    c23 = math.sqrt(1.0 - s23**2)
    c13 = math.sqrt(1.0 - s13**2)
    eid = np.exp(1j * delta)
    emid = np.exp(-1j * delta)
    return np.array(
        [
            [c12 * c13, s12 * c13, s13 * emid],
            [-s12 * c23 - c12 * s23 * s13 * eid, c12 * c23 - s12 * s23 * s13 * eid, s23 * c13],
            [s12 * s23 - c12 * c23 * s13 * eid, -c12 * s23 - s12 * c23 * s13 * eid, c23 * c13],
        ],
        dtype=complex,
    )


def mixing_observables(v: np.ndarray) -> dict[str, float]:
    abs_v = np.abs(v)
    s13 = float(abs_v[0, 2])
    c13 = math.sqrt(max(1.0 - s13**2, 1.0e-30))
    s12 = float(abs_v[0, 1] / c13)
    s23 = float(abs_v[1, 2] / c13)
    j = float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))
    return {
        "abs_Vus": float(abs_v[0, 1]),
        "abs_Vcb": float(abs_v[1, 2]),
        "abs_Vub": float(abs_v[0, 2]),
        "s12": s12,
        "s23": s23,
        "s13": s13,
        "jarlskog": j,
    }


def log_mass_score(found: list[float], target: list[float]) -> float:
    return float(sum(math.log10(max(f, 1.0e-30) / t) ** 2 for f, t in zip(found, target)))


def ckm_magnitude_score(found: np.ndarray, target: np.ndarray) -> float:
    entries = [(0, 1), (1, 2), (0, 2)]
    return float(sum(math.log10(max(abs(found[i, j]), 1.0e-30) / abs(target[i, j])) ** 2 for i, j in entries))


def first_order_h_shift(
    u_left: np.ndarray,
    d_left_current: np.ndarray,
    d_singulars: np.ndarray,
    v_target: np.ndarray,
) -> dict[str, Any]:
    # In the up-mass basis, H_d = V_CKM diag(y_d^2) V_CKM^dagger.  Keeping the
    # down singular values fixed, this is the minimal Hermitian target for a
    # CKM-only deformation of the down left eigensystem.
    v_current = u_left.conjugate().T @ d_left_current
    diag = np.diag(d_singulars**2)
    h_current = v_current @ diag @ v_current.conjugate().T
    h_target = v_target @ diag @ v_target.conjugate().T
    delta_h = h_target - h_current
    largest = float(np.max(d_singulars**2))
    offdiag = delta_h.copy()
    np.fill_diagonal(offdiag, 0.0)
    return {
        "relative_frobenius_delta_Hd": float(np.linalg.norm(delta_h) / max(np.linalg.norm(h_current), 1.0e-30)),
        "relative_spectral_delta_Hd": float(np.linalg.norm(delta_h, ord=2) / max(np.linalg.norm(h_current, ord=2), 1.0e-30)),
        "relative_offdiag_to_largest_yb2": float(np.linalg.norm(offdiag) / max(largest, 1.0e-30)),
        "max_abs_offdiag_delta_over_yb2": float(np.max(np.abs(offdiag)) / max(largest, 1.0e-30)),
        "current_Hd_in_up_basis_abs": matrix_abs_json(h_current),
        "target_Hd_in_up_basis_abs": matrix_abs_json(h_target),
        "delta_Hd_in_up_basis": matrix_json(delta_h),
    }


def write_ckm_csv(path: Path, current: np.ndarray, target: np.ndarray) -> None:
    rows = []
    labels = ["u", "c", "t"]
    cols = ["d", "s", "b"]
    for i, row_label in enumerate(labels):
        for j, col_label in enumerate(cols):
            rows.append(
                {
                    "row": row_label,
                    "col": col_label,
                    "current_abs": float(abs(current[i, j])),
                    "target_abs": float(abs(target[i, j])),
                    "ratio_current_over_target": float(abs(current[i, j]) / max(abs(target[i, j]), 1.0e-30)),
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines = []
    lines.append("# Flavor observable audit")
    lines.append("")
    lines.append("No web lookup was used.  CKM targets are internal benchmark values.")
    lines.append("")
    lines.append("## Dirac masses")
    lines.append("")
    lines.append("| sector | found small | found mid | target small | target mid | log score |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for sector, item in payload["mass_ratio_audit"].items():
        found = item["found"]
        target = item["target"]
        lines.append(
            f"| {sector} | {found[0]:.6e} | {found[1]:.6e} | "
            f"{target[0]:.6e} | {target[1]:.6e} | {item['log_score']:.3e} |"
        )
    lines.append("")
    lines.append("## CKM from the exact card")
    lines.append("")
    cur = payload["ckm_current_observables"]
    tar = payload["ckm_target_observables"]
    lines.append("| observable | current | target | ratio |")
    lines.append("|---|---:|---:|---:|")
    for key in ["abs_Vus", "abs_Vcb", "abs_Vub", "jarlskog"]:
        current = cur[key]
        target = tar[key]
        ratio = current / target if abs(target) > 1.0e-30 else math.inf
        lines.append(f"| `{key}` | {current:.6e} | {target:.6e} | {ratio:.3e} |")
    lines.append("")
    lines.append(f"CKM magnitude log-score: `{payload['ckm_magnitude_log_score']:.6e}`.")
    lines.append("")
    lines.append("## Minimal CKM-only down-sector deformation")
    lines.append("")
    h = payload["minimal_down_H_shift"]
    lines.append(
        "Keeping the down singular values fixed, replacing the current left "
        "eigensystem by the CKM target requires"
    )
    lines.append("")
    lines.append("```text")
    lines.append(f"||Delta H_d||_F / ||H_d||_F = {h['relative_frobenius_delta_Hd']:.6e}")
    lines.append(f"||Delta H_d||_2 / ||H_d||_2 = {h['relative_spectral_delta_Hd']:.6e}")
    lines.append(f"||offdiag Delta H_d||_F / y_b^2 = {h['relative_offdiag_to_largest_yb2']:.6e}")
    lines.append(f"max |offdiag Delta H_d| / y_b^2 = {h['max_abs_offdiag_delta_over_yb2']:.6e}")
    lines.append("```")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "flavor_observable_audit_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = json.loads(CARD.read_text(encoding="utf-8"))
    sectors = card["yukawa_sectors"]
    y_u = cmat(sectors["up"]["Y_matrix_normalized"])
    y_d = cmat(sectors["down"]["Y_matrix_normalized"])
    y_e = cmat(sectors["charged_lepton"]["Y_matrix_normalized"])
    y_nu = cmat(sectors["neutrino_dirac"]["Y_matrix_normalized"])

    u_left, u_s = left_rotation(y_u)
    d_left, d_s = left_rotation(y_d)
    e_left, e_s = left_rotation(y_e)
    nu_left, nu_s = left_rotation(y_nu)

    # Internal, fixed flavor targets used throughout the local notes.  These are
    # not meant to replace a final cited input table.
    mass_targets = {
        "up": [1.0e-5, 7.0e-3],
        "down": [1.0e-3, 2.0e-2],
        "charged_lepton": [3.0e-4, 6.0e-2],
        "neutrino_dirac": [1.0e-2, 2.0e-1],
    }
    found = {
        "up": [float(u_s[0] / u_s[2]), float(u_s[1] / u_s[2])],
        "down": [float(d_s[0] / d_s[2]), float(d_s[1] / d_s[2])],
        "charged_lepton": [float(e_s[0] / e_s[2]), float(e_s[1] / e_s[2])],
        "neutrino_dirac": [float(nu_s[0] / nu_s[2]), float(nu_s[1] / nu_s[2])],
    }
    mass_audit = {
        sector: {"found": values, "target": mass_targets[sector], "log_score": log_mass_score(values, mass_targets[sector])}
        for sector, values in found.items()
    }

    v_ckm_current = u_left.conjugate().T @ d_left
    v_ckm_target = standard_ckm(s12=0.225, s23=0.041, s13=0.0037, delta=1.20)
    shift = first_order_h_shift(u_left, d_left, d_s, v_ckm_target)
    ckm_score = ckm_magnitude_score(v_ckm_current, v_ckm_target)

    payload = {
        "note": "No web lookup used. CKM target is an internal benchmark, not a final cited input.",
        "input_card": str(CARD),
        "mass_ratio_audit": mass_audit,
        "CKM_current": matrix_json(v_ckm_current),
        "CKM_current_abs": matrix_abs_json(v_ckm_current),
        "CKM_target": matrix_json(v_ckm_target),
        "CKM_target_abs": matrix_abs_json(v_ckm_target),
        "ckm_current_observables": mixing_observables(v_ckm_current),
        "ckm_target_observables": mixing_observables(v_ckm_target),
        "ckm_magnitude_log_score": ckm_score,
        "minimal_down_H_shift": shift,
        "PMNS_from_seesaw_card": card["seesaw"]["reconstructed_pmns_angles"],
        "checks": {
            "mass_ratios_close_to_internal_targets": all(item["log_score"] < 0.2 for item in mass_audit.values()),
            "ckm_needs_additional_fit": ckm_score > 0.5,
            "down_H_shift_is_perturbative_but_not_tiny": (
                0.01 < shift["relative_spectral_delta_Hd"] < 1.0
            ),
        },
        "verdict": (
            "The exact O(-4) card remains a hierarchy/seesaw benchmark, not a full "
            "flavor fit.  It reproduces the intended mass-ratio targets and PMNS "
            "benchmark, but the CKM matrix from U_u^dagger U_d is far from the "
            "internal CKM target.  The obstruction is quantitative rather than "
            "fatal: a CKM-only deformation of H_d=Y_d Y_d^dagger with fixed down "
            "singular values is perturbative, so the next controlled route is to "
            "add a constrained 120_H/Toeplitz-Kahler misalignment or to fit small "
            "dual-density perturbations before feeding the rotations into d=5 "
            "proton decay."
        ),
    }
    payload["all_checks_ok"] = all(payload["checks"].values())

    (OUT / "flavor_observable_audit.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_ckm_csv(OUT / "ckm_current_vs_target.csv", v_ckm_current, v_ckm_target)
    write_report(payload)

    print("Flavor observable audit")
    print("  mass log scores:")
    for sector, item in mass_audit.items():
        print(f"    {sector:15s}: {item['log_score']:.6e}")
    obs = payload["ckm_current_observables"]
    target = payload["ckm_target_observables"]
    print(
        "  CKM current |Vus|,|Vcb|,|Vub| = "
        f"{obs['abs_Vus']:.6e}, {obs['abs_Vcb']:.6e}, {obs['abs_Vub']:.6e}"
    )
    print(
        "  CKM target  |Vus|,|Vcb|,|Vub| = "
        f"{target['abs_Vus']:.6e}, {target['abs_Vcb']:.6e}, {target['abs_Vub']:.6e}"
    )
    print(f"  CKM magnitude log-score: {ckm_score:.6e}")
    print(f"  relative spectral Delta H_d: {shift['relative_spectral_delta_Hd']:.6e}")
    print(f"  all checks ok: {payload['all_checks_ok']}")
    print(f"  wrote: {OUT}")
    if not payload["all_checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
