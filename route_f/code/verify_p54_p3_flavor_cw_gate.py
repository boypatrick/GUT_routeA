#!/usr/bin/env python3
"""Diagnostic input helpers and compatibility entry point for the P54 audit.

The 2026-09-02 tree conjugation/exclusion and real-slice CW implementation
has been superseded. The authoritative calculation is verify_p54_theory_audit.
SM-only transport below is not the Pati-Salam matching likelihood.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OVERLAP_JSON = ROUTE_F / "output" / "p54_p3_light_doublet_overlaps.json"
P2_JSON = ROUTE_F / "output" / "p54_p2_running_thresholds_cosmology.json"
P2_CW_JSON = ROUTE_F / "output" / "p54_p2_bosonic_cw.json"
OUTPUT_JSON = ROUTE_F / "output" / "p54_p3_flavor_cw_gate.json"
OUTPUT_MD = ROUTE_F / "output" / "p54_p3_flavor_cw_gate.md"
SCRIPT = Path(__file__).resolve()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_json(value: complex) -> dict[str, float]:
    return {"re": float(value.real), "im": float(value.imag)}


def matrix_json(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    return [[complex_json(complex(value)) for value in row] for row in matrix]


def singular_values(matrix: np.ndarray) -> list[float]:
    return np.sort(np.linalg.svd(matrix, compute_uv=False)).astype(float).tolist()


def standard_ckm(s12: float, s23: float, s13: float, delta: float) -> np.ndarray:
    c12, c23, c13 = (math.sqrt(1.0 - value**2) for value in (s12, s23, s13))
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


def diagnostic_sm_run(mu_high: float) -> dict[str, Any]:
    """One-loop SM transport for diagnostic target scales only.

    Inputs: Mummidi--Patel Table II and PDG 2024 CKM parametrisation.
    PS Yukawa running and matching are not included; no full exclusion follows.
    """

    mu0 = 173.1
    vev = 174.0
    masses_gev = {
        "u": 1.21e-3,
        "c": 0.61,
        "t": 163.35,
        "d": 2.58e-3,
        "s": 52.74e-3,
        "b": 2.72,
        "e": 0.499e-3,
        "mu": 0.104,
        "tau": 1.759,
    }
    ckm = standard_ckm(0.22501, 0.04183, 0.003732, 1.147)
    yu0 = np.diag([masses_gev[key] / vev for key in ("u", "c", "t")]).astype(complex)
    yd0 = ckm @ np.diag([masses_gev[key] / vev for key in ("d", "s", "b")])
    ye0 = np.diag([masses_gev[key] / vev for key in ("e", "mu", "tau")]).astype(complex)

    def pack(gauge: np.ndarray, matrices: list[np.ndarray]) -> np.ndarray:
        data: list[float] = gauge.astype(float).tolist()
        for matrix in matrices:
            data.extend(matrix.real.ravel().tolist())
            data.extend(matrix.imag.ravel().tolist())
        return np.asarray(data, dtype=float)

    def unpack(data: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
        gauge = data[:3]
        matrices = []
        cursor = 3
        for _ in range(3):
            real = data[cursor : cursor + 9].reshape(3, 3)
            imag = data[cursor + 9 : cursor + 18].reshape(3, 3)
            matrices.append(real + 1j * imag)
            cursor += 18
        return gauge, matrices

    initial = pack(np.array([0.4632, 0.6540, 1.1630]), [yu0, yd0, ye0])

    def beta(_log_mu: float, data: np.ndarray) -> np.ndarray:
        (g1, g2, g3), (yu, yd, ye) = unpack(data)
        trace = float(
            3.0 * np.trace(yu.conj().T @ yu).real
            + 3.0 * np.trace(yd.conj().T @ yd).real
            + np.trace(ye.conj().T @ ye).real
        )
        ident = np.eye(3)
        bu = (
            1.5 * (yu @ yu.conj().T - yd @ yd.conj().T) @ yu
            + (trace - 17.0 / 20.0 * g1**2 - 9.0 / 4.0 * g2**2 - 8.0 * g3**2) * ident @ yu
        )
        bd = (
            1.5 * (yd @ yd.conj().T - yu @ yu.conj().T) @ yd
            + (trace - 1.0 / 4.0 * g1**2 - 9.0 / 4.0 * g2**2 - 8.0 * g3**2) * ident @ yd
        )
        be = 1.5 * (ye @ ye.conj().T) @ ye + (
            trace - 9.0 / 4.0 * g1**2 - 9.0 / 4.0 * g2**2
        ) * ident @ ye
        loop = 1.0 / (16.0 * math.pi**2)
        bg = loop * np.array([41.0 / 10.0 * g1**3, -19.0 / 6.0 * g2**3, -7.0 * g3**3])
        return pack(bg, [loop * bu, loop * bd, loop * be])

    solution = solve_ivp(
        beta,
        (0.0, math.log(mu_high / mu0)),
        initial,
        rtol=2.0e-10,
        atol=2.0e-12,
        method="DOP853",
    )
    if not solution.success:
        raise RuntimeError(f"Diagnostic RGE integration failed: {solution.message}")
    gauge, (yu, yd, ye) = unpack(solution.y[:, -1])
    u_left, _, _ = np.linalg.svd(yu)
    d_left, _, _ = np.linalg.svd(yd)
    # SVD orders heavy to light; publish the conventional (u,c,t)/(d,s,b) order.
    ckm_high = u_left[:, ::-1].conj().T @ d_left[:, ::-1]
    return {
        "status": "diagnostic_only_missing_Pati_Salam_Yukawa_RGE_and_threshold_matching",
        "mu_low_GeV": mu0,
        "mu_high_GeV": mu_high,
        "gauge_at_high": gauge.tolist(),
        "yukawa_singular_values_at_high": {
            "up": singular_values(yu),
            "down": singular_values(yd),
            "charged_lepton": singular_values(ye),
        },
        "CKM_abs_at_high": np.abs(ckm_high).tolist(),
        "source_common_scale_inputs": "Mummidi--Patel JHEP 12 (2021) 042, Table II",
        "source_ckm": "PDG 2024 CKM review, Eq. (12.28)",
    }


def phase_uniform_bound(r_abs: float, s_abs: float, yb: float, ytau: float) -> dict[str, float]:
    """Return ||Yu|| <= A ||Yd|| + B ||Ye|| for arbitrary phases."""

    coefficient_d = r_abs * (3.0 + s_abs) / 4.0
    coefficient_e = r_abs * (1.0 + s_abs) / 4.0
    upper = coefficient_d * yb + coefficient_e * ytau
    return {"A_down": coefficient_d, "B_charged_lepton": coefficient_e, "yt_upper": upper}


def gaussian_halfspace_chi2_lower(
    yt: float,
    yb: float,
    ytau: float,
    coefficient_d: float,
    coefficient_e: float,
    fractional_sigma: float = 0.10,
) -> float:
    """Mahalanobis distance from the target to t-A b-B tau <= 0."""

    violation = yt - coefficient_d * yb - coefficient_e * ytau
    if violation <= 0.0:
        return 0.0
    denominator = (
        (fractional_sigma * yt) ** 2
        + (coefficient_d * fractional_sigma * yb) ** 2
        + (coefficient_e * fractional_sigma * ytau) ** 2
    )
    return violation**2 / denominator


def appendix_b_matrices() -> tuple[np.ndarray, np.ndarray]:
    """Published effective H' and F' matrices, not a P54 fit result."""

    h_eff = np.diag(np.array([0.00023, -0.04811, -5.79504]) * 1.0e-3).astype(complex)
    f_eff = np.array(
        [
            [-0.0088 + 0.0178j, 0.0475 - 0.0889j, 0.4635 + 0.6797j],
            [0.0475 - 0.0889j, 1.1279 + 0.5108j, -1.2218 - 2.5921j],
            [0.4635 + 0.6797j, -1.2218 - 2.5921j, 5.4683 - 5.9856j],
        ],
        dtype=complex,
    ) * 1.0e-4
    return h_eff, f_eff


def run():
    """Return the corrected audit; retained for caller compatibility."""
    from verify_p54_theory_audit import run as audited_run
    return audited_run()


def main() -> None:
    from verify_p54_theory_audit import main as audited_main
    audited_main()


if __name__ == "__main__":
    main()
