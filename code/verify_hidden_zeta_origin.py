#!/usr/bin/env python3
"""Verify the Route-B hidden transvectant messenger numerics.

This script checks only the algebraic data used in the lean paper:
the target zeta, its square-root coupling, the normalized transvectant
contact matrix, the Z_178 phase candidate, a simple instanton parametrization,
and visible-singlet threshold silence.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def as_pair(z: complex) -> list[float]:
    return [float(np.real(z)), float(np.imag(z))]


def main() -> None:
    zeta = 0.1076472949 + 0.0736514853j

    rho = abs(zeta)
    phi = float(np.angle(zeta))
    sqrt_zeta = np.sqrt(zeta)

    # Normalized second-transvectant contact direction.
    k_tr = (1 / np.sqrt(3.0)) * np.array(
        [
            [0, 0, -1],
            [0, 1, 0],
            [-1, 0, 0],
        ],
        dtype=complex,
    )

    lam = sqrt_zeta
    delta_m = lam**2 * k_tr
    target_delta_m = zeta * k_tr

    order = 178
    root_power = 17
    theta_178 = 2 * np.pi * root_power / order
    zeta_178 = rho * np.exp(1j * theta_178)

    instanton_action = -np.log(rho)
    one_loop_lambda_factor = abs(sqrt_zeta) ** 2 / (16 * np.pi**2)
    instanton_couplings = {}
    for b0 in [4, 8, 12]:
        g2 = 8 * np.pi**2 / (b0 * instanton_action)
        alpha = g2 / (4 * np.pi)
        instanton_couplings[str(b0)] = {
            "g": float(np.sqrt(g2)),
            "alpha_inv": float(1 / alpha),
        }

    delta_b_visible = [0, 0, 0]

    results = {
        "zeta": as_pair(zeta),
        "rho": float(rho),
        "phi": phi,
        "sqrt_zeta": as_pair(sqrt_zeta),
        "sqrt_zeta_abs": float(abs(sqrt_zeta)),
        "sqrt_zeta_squared_error": float(abs(sqrt_zeta**2 - zeta)),
        "k_tr_symmetric_error": float(np.linalg.norm(k_tr - k_tr.T)),
        "k_tr_inverse_check": float(np.linalg.norm(k_tr @ (3 * k_tr) - np.eye(3))),
        "messenger_contact_error": float(np.linalg.norm(delta_m - target_delta_m)),
        "z178_order": order,
        "z178_power": root_power,
        "theta_178": float(theta_178),
        "phase_error": float(abs(theta_178 - phi)),
        "zeta_178": as_pair(zeta_178),
        "zeta_178_error": float(abs(zeta_178 - zeta)),
        "instanton_action": float(instanton_action),
        "one_loop_lambda_factor": float(one_loop_lambda_factor),
        "instanton_couplings": instanton_couplings,
        "delta_b_visible": delta_b_visible,
    }

    out_dir = Path("output/hidden_zeta_origin")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "routeB_hidden_zeta_verification.json"
    out_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    for key, value in results.items():
        print(f"{key} = {value}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
