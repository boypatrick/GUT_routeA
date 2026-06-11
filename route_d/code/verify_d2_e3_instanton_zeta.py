#!/usr/bin/env python3
"""Route-D D2-E3 numerical checks for the instanton zeta audit.

This script checks only convention-independent arithmetic attached to the
benchmark coefficient.  It does not verify the existence of an E3 divisor, its
charged zero-mode spectrum, or a complete Green-Schwarz selection rule in a
global compactification.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def main() -> None:
    zeta = complex(0.1076472949, 0.0736514853)
    rho = abs(zeta)
    phase = math.atan2(zeta.imag, zeta.real)
    s_eff_a1 = -math.log(rho)
    two_instanton_scale = rho**2

    delta_s_window = 3.208798e-5
    delta_phase_window = 5.424119e-5
    two_inst_over_delta_s = two_instanton_scale / delta_s_window

    # Convention: zeta = A exp(2 pi i T), T = axion + i volume.
    # If |A|=1, Im(T)=S_eff/(2 pi).
    im_t_if_prefactor_one = s_eff_a1 / (2.0 * math.pi)

    # Minimal charge bookkeeping for a Majorana-only instanton insertion:
    # if q(N)=+1, the perturbatively forbidden NN operator has charge +2.
    # The instanton factor must transform with charge -2.  This is not a
    # zero-mode proof; it is only the algebraic charge target.
    q_n = 1
    q_nn = 2 * q_n
    q_instanton_factor = -q_nn
    charge_neutral_majorana = q_nn + q_instanton_factor == 0

    report = {
        "zeta_real": zeta.real,
        "zeta_imag": zeta.imag,
        "abs_zeta": rho,
        "arg_zeta_rad": phase,
        "effective_action_if_prefactor_one": s_eff_a1,
        "im_T_if_zeta_equals_exp_2pi_i_T_and_prefactor_one": im_t_if_prefactor_one,
        "two_instanton_scale_abs_zeta_squared": two_instanton_scale,
        "delta_s_window": delta_s_window,
        "delta_phase_window_rad": delta_phase_window,
        "two_instanton_scale_over_delta_s_window": two_inst_over_delta_s,
        "q_N_toy": q_n,
        "q_NN_toy": q_nn,
        "q_instanton_factor_required": q_instanton_factor,
        "charge_neutral_majorana_toy_check_pass": charge_neutral_majorana,
        "zero_mode_spectrum_checked": False,
        "global_e3_divisor_constructed": False,
        "ten_digit_zeta_prediction_claimed": False,
    }

    out_dir = Path(__file__).resolve().parents[1] / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "d2_e3_instanton_zeta.json"
    md_path = out_dir / "d2_e3_instanton_zeta.md"

    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    md_path.write_text(
        "# D2-E3 E3-Instanton Zeta Verification\n\n"
        f"`|zeta| = {rho}`.\n"
        f"`arg(zeta) = {phase}` rad.\n"
        f"`S_eff = -log(|zeta|) = {s_eff_a1}` if the prefactor has unit modulus.\n"
        f"If `zeta = exp(2 pi i T)`, then `Im(T) = {im_t_if_prefactor_one}` "
        "for unit prefactor.\n"
        f"`|zeta|^2 = {two_instanton_scale}`.\n"
        f"`|zeta|^2 / Delta_s = {two_inst_over_delta_s}` for "
        f"`Delta_s = {delta_s_window}`.\n\n"
        "Toy GS charge check: if `q(N)=+1`, then `q(NN)=+2`; an instanton "
        "factor of charge `-2` makes the Majorana insertion neutral.  This is "
        "only charge bookkeeping, not a global zero-mode proof.\n\n"
        "Boundary: this verification does not construct an E3 divisor, does "
        "not count charged zero modes, and does not claim a ten-digit "
        "`zeta` prediction.\n"
    )

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    print(
        "D2-E3 checks: "
        f"|zeta|={rho:.16g}, S_eff={s_eff_a1:.16g}, "
        f"|zeta|^2/Delta_s={two_inst_over_delta_s:.6g}"
    )


if __name__ == "__main__":
    main()
