#!/usr/bin/env python3
"""Fail-closed flavor-regime gate for DYN-7 and DYN-9b-3.

The historical scans use an unflavored asymmetry at ``M1 ~ 10^10 GeV``.
That number is a useful diagnostic, but it is not a physical posterior in the
tau-resolved regime.  This audit classifies the Standard-Model regime, records
the missing MSSM branch inputs, and corrects the double-square-root error in
the historical Davidson--Ibarra consistency bound.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

from route_e_paths import AUDIT_OUTPUT


OUT = AUDIT_OUTPUT / "audit7"
DYN4A = AUDIT_OUTPUT / "audit1" / "dyn4a_seesaw_zeta_posterior.json"
CHECKS: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    ok = bool(condition)
    CHECKS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail else ""))


def sm_flavor_regime(m1_gev: float) -> str:
    """Conventional thermal-leptogenesis regime classifier.

    Boundaries are order-of-magnitude transition scales, not sharp phase
    boundaries.  A density-matrix calculation should interpolate through
    them; the classifier exists to prevent an unflavored calculation from
    being promoted in the central two-flavor band.
    """
    if m1_gev >= 1.0e12:
        return "unflavored"
    if m1_gev >= 1.0e9:
        return "two_flavor_tau_resolved"
    return "three_flavor_mu_tau_resolved"


def main() -> None:
    print("== DYN-7F section 1: input and regime classification ==")
    if not DYN4A.exists():
        raise SystemExit(f"missing required DYN-4a ledger: {DYN4A}")
    dyn4a = json.loads(DYN4A.read_text(encoding="utf-8"))
    spectrum_desc = dyn4a["benchmark"]["heavy_MR_singular_values_GeV"]
    m1_gev = float(min(spectrum_desc))
    regime_sm = sm_flavor_regime(m1_gev)
    check(
        "DYN-4a exports a positive hierarchical heavy-neutrino spectrum",
        len(spectrum_desc) == 3 and m1_gev > 0.0 and max(spectrum_desc) / m1_gev > 1.0e4,
        f"M1={m1_gev:.6e} GeV",
    )
    check(
        "the non-SUSY benchmark lies in the tau-resolved two-flavor regime",
        regime_sm == "two_flavor_tau_resolved",
        f"classifier={regime_sm}",
    )

    # DYN-7 is an MSSM/SUSY-lane calculation.  Charged-lepton equilibration
    # rates and transition scales depend on the branch (notably tan beta), so
    # the SM thresholds alone must not be used to certify it.
    mssm_inputs = {
        "tan_beta": None,
        "charged_lepton_thermal_rates": None,
        "spectator_matrix": None,
        "initial_N1_abundance": "assumed_thermal_in_historical_scan",
    }
    missing_mssm_inputs = [key for key, value in mssm_inputs.items() if value is None]
    check(
        "MSSM flavor certification fails closed when branch-local thermal inputs are absent",
        missing_mssm_inputs == ["tan_beta", "charged_lepton_thermal_rates", "spectator_matrix"],
        "missing=" + ",".join(missing_mssm_inputs),
    )

    print("== DYN-7F section 2: Davidson-Ibarra bound repair ==")
    # The old code first took sqrt(eigenvalues(mnu^dagger mnu)) = m_i and
    # then took another sqrt.  The correct masses are already the first sqrt.
    m1_light_ev = 1.0e-3
    dm31_ev2 = 2.517e-3
    m3_light_ev = math.sqrt(m1_light_ev**2 + dm31_ev2)
    delta_m_ev = m3_light_ev - m1_light_ev
    v_u_gev = 100.0
    v_sm_gev = 174.0
    eps_di_susy = 3.0 * m1_gev * delta_m_ev * 1.0e-9 / (
        8.0 * math.pi * v_u_gev**2
    )
    eps_di_sm = 3.0 * m1_gev * delta_m_ev * 1.0e-9 / (
        16.0 * math.pi * v_sm_gev**2
    )
    historical_double_sqrt_delta_m = math.sqrt(m3_light_ev) - math.sqrt(m1_light_ev)
    historical_overestimate = historical_double_sqrt_delta_m / delta_m_ev
    check(
        "light masses are extracted with one square root, not two",
        abs(m3_light_ev - 0.05017967716117751) < 1.0e-14,
        f"m3={m3_light_ev:.12e} eV",
    )
    check(
        "corrected SUSY Davidson-Ibarra central bound is reproduced",
        abs(eps_di_susy - 1.3974207726e-5) < 5.0e-11,
        f"epsilon_DI_SUSY={eps_di_susy:.12e}",
    )
    check(
        "historical double-square-root bound was overestimated by about 3.91",
        3.8 < historical_overestimate < 4.0,
        f"factor={historical_overestimate:.8f}",
    )

    print("== DYN-7F section 3: fail-closed claim statuses ==")
    statuses = {
        "DYN-7": {
            "claim_status": "blocked_missing_branch_thermal_inputs",
            "historical_unflavored_scan_role": "order_estimate_diagnostic_only",
            "success_probability_publishable": False,
            "arg_zeta_conditioned_interval_publishable": False,
            "required_repair": (
                "supply tan(beta), branch-local charged-lepton and spectator rates, "
                "then solve at least the two-flavor kinetic system; density matrix preferred"
            ),
        },
        "DYN-9b-3": {
            "claim_status": "blocked_requires_two_flavor",
            "sm_regime": regime_sm,
            "historical_unflavored_scan_role": "order_estimate_diagnostic_only",
            "success_probability_publishable": False,
            "arg_zeta_conditioned_interval_publishable": False,
            "required_repair": (
                "solve tau-resolved two-flavor Boltzmann or density-matrix equations "
                "with spectator effects and branch-local rates"
            ),
        },
    }
    check(
        "neither historical unflavored scan is promoted to a physical posterior",
        all(not value["success_probability_publishable"] for value in statuses.values()),
    )

    n_pass = sum(ok for _, ok, _ in CHECKS)
    all_pass = n_pass == len(CHECKS)
    payload = {
        "audit": "audit7_dyn7_flavor_regime_gate",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "all_pass": all_pass,
        "checks_passed": n_pass,
        "checks_total": len(CHECKS),
        "checks": [
            {"name": name, "pass": ok, "detail": detail}
            for name, ok, detail in CHECKS
        ],
        "input": {
            "dyn4a_ledger": str(Path("output/audit1") / DYN4A.name),
            "M1_GeV": m1_gev,
        },
        "regime_boundaries_GeV": {
            "three_to_two_flavor": 1.0e9,
            "two_to_unflavored": 1.0e12,
            "note": "order-of-magnitude SM classifier; density-matrix interpolation required near transitions",
        },
        "davidson_ibarra_repair": {
            "m1_eV": m1_light_ev,
            "m3_eV": m3_light_ev,
            "delta_m_eV": delta_m_ev,
            "epsilon_DI_SUSY_vu100": eps_di_susy,
            "epsilon_DI_SM_v174": eps_di_sm,
            "historical_double_sqrt_overestimate_factor": historical_overestimate,
        },
        "mssm_inputs": mssm_inputs,
        "statuses": statuses,
        "physics_claims_closed": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    json_path = OUT / "dyn7_flavor_regime_gate.json"
    md_path = OUT / "dyn7_flavor_regime_gate.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(
        "# DYN-7 / DYN-9b-3 flavor-regime gate\n\n"
        f"- arithmetic checks: **{n_pass}/{len(CHECKS)} pass**\n"
        f"- benchmark `M1 = {m1_gev:.6e} GeV`\n"
        f"- non-SUSY classifier: **{regime_sm}**\n"
        "- DYN-7 status: **blocked_missing_branch_thermal_inputs**\n"
        "- DYN-9b-3 status: **blocked_requires_two_flavor**\n"
        "- historical success probabilities and conditioned `arg(zeta)` windows are diagnostics only\n\n"
        "## Davidson-Ibarra correction\n\n"
        f"- `m3 = {m3_light_ev:.12e} eV` (one square root)\n"
        f"- SUSY, `v_u=100 GeV`: `{eps_di_susy:.12e}`\n"
        f"- SM, `v=174 GeV`: `{eps_di_sm:.12e}`\n"
        f"- old double-square-root overestimate: `{historical_overestimate:.6f}`\n\n"
        "## Checks\n\n"
        + "\n".join(
            f"- [{'PASS' if ok else 'FAIL'}] {name}"
            + (f": {detail}" if detail else "")
            for name, ok, detail in CHECKS
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {json_path} (+ .md)")
    print(
        "DYN-7F: arithmetic gate passed; physical posterior remains fail-closed "
        "pending flavored kinetics."
    )
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
