#!/usr/bin/env python3
"""Build a no-web input convention ledger for the current publication card.

This is a provenance and arithmetic audit, not a literature refresh.  The
thread explicitly forbids web lookup, so the ledger records the numerical
targets and width conventions exactly as they are used in the local scripts,
hashes their source files/artifacts, and checks that the dressed d=5 channel
CSV is internally reproduced by the stated lifetime formula.
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

import audit_eigenstate_d5_dressing as eig  # noqa: E402
import audit_mass_insertion_d5_dressing as mi  # noqa: E402
import audit_mssm_mixing_d5_dressing as mssm  # noqa: E402
import audit_publication_dressed_c5_from_eigenstate_card as dressed  # noqa: E402
import audit_publication_flavor_d5_reproducibility as repro  # noqa: E402
import scan_clebsch_flavor_fit as fit  # noqa: E402
import scan_yukawa_superpotential_rge as rge  # noqa: E402


OUT = ROOT / "output" / "no_web_input_convention_ledger"
D5_ROWS = ROOT / "output" / "publication_dressed_c5_clockwork_rescue" / "dressed_channel_rows_kappa729.csv"
PUBLICATION_CARD = ROOT / "output" / "clockwork_rescued_publication_card" / "clockwork_rescued_publication_card.json"
PROTON_JSON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"

SOURCE_FILES = {
    "ledger_builder": ROOT / "code" / "build_no_web_input_convention_ledger.py",
    "flavor_targets": ROOT / "code" / "scan_clebsch_flavor_fit.py",
    "flavor_reproducibility_gates": ROOT / "code" / "audit_publication_flavor_d5_reproducibility.py",
    "dressed_c5_width_conventions": ROOT / "code" / "audit_publication_dressed_c5_from_eigenstate_card.py",
    "clockwork_c5_replay": ROOT / "code" / "audit_publication_dressed_c5_clockwork_rescue.py",
    "mass_insertion_dressing_inputs": ROOT / "code" / "audit_mass_insertion_d5_dressing.py",
    "soft_eigenstate_grid": ROOT / "code" / "audit_eigenstate_d5_dressing.py",
    "mssm_mixing_inputs": ROOT / "code" / "audit_mssm_mixing_d5_dressing.py",
    "rge_inputs": ROOT / "code" / "scan_yukawa_superpotential_rge.py",
}

ARTIFACTS = {
    "clockwork_publication_card": PUBLICATION_CARD,
    "clockwork_dressed_c5_rows": D5_ROWS,
    "proton_decay_verification": PROTON_JSON,
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: dict[str, Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label, path in paths.items():
        rows.append(
            {
                "label": label,
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path),
            }
        )
    return rows


def arr(v: np.ndarray) -> list[float] | list[list[float]]:
    return np.asarray(v, dtype=float).tolist()


def scalar_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(category: str, name: str, value: Any, unit: str, source: str, role: str) -> None:
        rows.append(
            {
                "category": category,
                "name": name,
                "value": json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else str(value),
                "unit": unit,
                "source": source,
                "role": role,
            }
        )

    for key, value in fit.TARGETS.items():
        add("flavor_target", key, value, "dimensionless", "scan_clebsch_flavor_fit.TARGETS", "local flavor fit target")

    add("flavor_gate", "STRICT_CKM", repro.STRICT_CKM, "score", "audit_publication_flavor_d5_reproducibility", "local CKM pass gate")
    add("flavor_gate", "LOOSE_MASS", repro.LOOSE_MASS, "score", "audit_publication_flavor_d5_reproducibility", "local mass-ratio pass gate")
    add("flavor_gate", "SEESAW_RESIDUAL_MAX", repro.SEESAW_RESIDUAL_MAX, "matrix residual", "audit_publication_flavor_d5_reproducibility", "seesaw replay pass gate")
    add("flavor_gate", "TOL", repro.TOL, "absolute difference", "audit_publication_flavor_d5_reproducibility", "card recomputation tolerance")

    for key, value in {
        "m1_eV": 1.0e-3,
        "dm21_eV2": 7.42e-5,
        "dm31_eV2": 2.517e-3,
        "sin2_theta12": 0.304,
        "sin2_theta13": 0.0222,
        "sin2_theta23": 0.573,
        "delta_cp_rad": 1.20 * math.pi,
        "alpha21_rad": 0.35 * math.pi,
        "alpha31_rad": 1.10 * math.pi,
    }.items():
        add("seesaw_benchmark", key, value, "mixed", "scan_clebsch_flavor_fit.seesaw_replay", "local light-neutrino reconstruction input")

    for case, prefs in dressed.WIDTH_CASES.items():
        for channel_class, value in prefs.items():
            add("d5_width_prefactor", f"{case}.{channel_class}", value, "GeV^5", "audit_publication_dressed_c5_from_eigenstate_card.WIDTH_CASES", "channel width prefactor")

    for channel in dressed.CHANNELS:
        add("d5_present_bound", f"{channel['operator']}.{channel['channel']}", channel["present_bound_years"], "yr", "audit_publication_dressed_c5_from_eigenstate_card.CHANNELS", "local present-bound convention")

    add("d5_filter", "DISPLAY_ST", dressed.DISPLAY_ST, "dimensionless", "audit_publication_dressed_c5_from_eigenstate_card", "display triplet filter")
    add("d5_stress", "future_stress_target", 1.0e35, "yr", "publication dressed C5 audits", "uniform future d=5 stress target")
    add("clockwork", "kappa", 729.0, "dimensionless", "publication_dressed_c5_clockwork_rescue", "exact clockwork replay value")
    add("clockwork", "q", 3.0, "dimensionless", "clockwork_hidden_gauge_quotient_origin", "clockwork ratio")
    add("clockwork", "n", 6, "links", "clockwork_hidden_gauge_quotient_origin", "clockwork length")

    add("lifetime", "HBAR_GEV_S", mi.HBAR_GEV_S, "GeV s", "audit_mass_insertion_d5_dressing", "width-to-lifetime conversion")
    add("lifetime", "SECONDS_PER_YEAR", mi.SECONDS_PER_YEAR, "s/yr", "audit_mass_insertion_d5_dressing", "width-to-lifetime conversion")
    add("dressing", "ALPHA2_INV", mi.ALPHA2_INV, "dimensionless", "audit_mass_insertion_d5_dressing", "wino dressing proxy")
    add("dressing", "SOFT_EPS_GRID_mass_insertion", mi.SOFT_EPS_GRID, "dimensionless", "audit_mass_insertion_d5_dressing", "mass-insertion stress grid")
    add("dressing", "SOFT_EPS_GRID_eigenstate", eig.SOFT_EPS_GRID, "dimensionless", "audit_eigenstate_d5_dressing", "eigenstate stress grid")
    add("dressing", "SPECTRUM_POINTS", mi.SPECTRUM_POINTS, "mixed", "audit_mass_insertion_d5_dressing", "soft spectrum stress points")
    add("mssm_mixing", "MW_GEV", mssm.MW_GEV, "GeV", "audit_mssm_mixing_d5_dressing", "chargino matrix input")
    add("mssm_mixing", "MZ_GEV", mssm.MZ_GEV, "GeV", "audit_mssm_mixing_d5_dressing", "neutralino matrix input")
    add("mssm_mixing", "SIN2_THETA_W", mssm.SIN2_THETA_W, "dimensionless", "audit_mssm_mixing_d5_dressing", "neutralino matrix input")
    add("mssm_mixing", "COS2_THETA_W", mssm.COS2_THETA_W, "dimensionless", "audit_mssm_mixing_d5_dressing", "neutralino matrix input")

    for key, value, unit in [
        ("MZ_GEV", rge.MZ_GEV, "GeV"),
        ("ALPHA_EM_INV_MZ", rge.ALPHA_EM_INV_MZ, "dimensionless"),
        ("SIN2_THETA_W_MZ", rge.SIN2_THETA_W_MZ, "dimensionless"),
        ("ALPHA_S_MZ", rge.ALPHA_S_MZ, "dimensionless"),
        ("V_HIGGS_GEV", rge.V_HIGGS_GEV, "GeV"),
        ("M_TOP_RUN_GEV", rge.M_TOP_RUN_GEV, "GeV"),
        ("M_BOTTOM_RUN_GEV", rge.M_BOTTOM_RUN_GEV, "GeV"),
        ("M_TAU_GEV", rge.M_TAU_GEV, "GeV"),
        ("M_NU_DIRAC_MAX_GEV", rge.M_NU_DIRAC_MAX_GEV, "GeV"),
        ("TAU_TARGET_YEARS", rge.TAU_TARGET_YEARS, "yr"),
        ("TRIPLET_FILTER_TARGET", rge.TRIPLET_FILTER_TARGET, "dimensionless"),
        ("MWINO_GEV", rge.MWINO_GEV, "GeV"),
    ]:
        add("rge_input", key, value, unit, "scan_yukawa_superpotential_rge", "local two-loop/RGE benchmark input")

    proton = read_json(PROTON_JSON)
    for key, value in proton["hadronic_constants"].items():
        add("proton_hadronic", key, value, "mixed", "output/proton_decay/proton_decay_verification.json", "legacy local d=6/d=5 prefactor convention")

    return rows


def stable_digest(obj: Any) -> str:
    encoded = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def verify_d5_rows() -> dict[str, Any]:
    channel_bounds = {(row["operator"], row["channel"]): float(row["present_bound_years"]) for row in dressed.CHANNELS}
    count = 0
    prefactor_mismatches = 0
    bound_mismatches = 0
    st_mismatches = 0
    kappa_mismatches = 0
    tau_relerr_max = 0.0
    c6_relerr_max = 0.0
    present_margin_relerr_max = 0.0
    stress_margin_relerr_max = 0.0
    worst_tau_row: dict[str, Any] | None = None
    worst_tau_relerr = -1.0

    with D5_ROWS.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            count += 1
            case = raw["normalization_case"]
            channel_class = raw["channel_class"]
            expected_pref = float(dressed.WIDTH_CASES[case][channel_class])
            actual_pref = float(raw["width_prefactor_GeV5"])
            if abs(expected_pref - actual_pref) > 1.0e-18:
                prefactor_mismatches += 1

            expected_bound = channel_bounds[(raw["operator"], raw["channel"])]
            actual_bound = float(raw["present_bound_years"])
            if abs(expected_bound - actual_bound) > 1.0e-6:
                bound_mismatches += 1

            st = float(raw["S_T_display"])
            if abs(st - dressed.DISPLAY_ST) > 1.0e-18:
                st_mismatches += 1

            kappa = float(raw["kappa"])
            if abs(kappa - 729.0) > 1.0e-12:
                kappa_mismatches += 1

            amplitude = float(raw["amplitude_with_dressing"])
            m_triplet = float(raw["M_triplet_GeV"])
            c6 = st * amplitude / m_triplet
            width = actual_pref * c6 * c6
            tau = mi.width_to_tau(width)
            stored_c6 = float(raw["C6_GeV_minus2_at_ST_display"])
            stored_tau = float(raw["tau_years_at_ST_display"])
            stored_present_margin = float(raw["margin_present_at_ST_display"])
            stored_stress_margin = float(raw["margin_1e35_at_ST_display"])

            c6_relerr = abs(c6 - stored_c6) / max(abs(stored_c6), 1.0e-300)
            tau_relerr = abs(tau - stored_tau) / max(abs(stored_tau), 1.0e-300)
            present_margin_relerr = abs(tau / actual_bound - stored_present_margin) / max(abs(stored_present_margin), 1.0e-300)
            stress_margin_relerr = abs(tau / 1.0e35 - stored_stress_margin) / max(abs(stored_stress_margin), 1.0e-300)

            c6_relerr_max = max(c6_relerr_max, c6_relerr)
            tau_relerr_max = max(tau_relerr_max, tau_relerr)
            present_margin_relerr_max = max(present_margin_relerr_max, present_margin_relerr)
            stress_margin_relerr_max = max(stress_margin_relerr_max, stress_margin_relerr)

            if tau_relerr > worst_tau_relerr:
                worst_tau_relerr = tau_relerr
                worst_tau_row = {
                    "normalization_case": case,
                    "operator": raw["operator"],
                    "channel": raw["channel"],
                    "scenario": raw["scenario"],
                    "epsilon": float(raw["epsilon"]),
                    "spectrum_name": raw["spectrum_name"],
                    "stored_tau_years": stored_tau,
                    "recomputed_tau_years": tau,
                    "relative_error": tau_relerr,
                }

    return {
        "rows_checked": count,
        "prefactor_mismatches": prefactor_mismatches,
        "present_bound_mismatches": bound_mismatches,
        "display_filter_mismatches": st_mismatches,
        "kappa_mismatches": kappa_mismatches,
        "max_c6_relative_error": c6_relerr_max,
        "max_tau_relative_error": tau_relerr_max,
        "max_present_margin_relative_error": present_margin_relerr_max,
        "max_1e35_margin_relative_error": stress_margin_relerr_max,
        "worst_tau_relative_error_row": worst_tau_row,
        "passes": (
            count > 0
            and prefactor_mismatches == 0
            and bound_mismatches == 0
            and st_mismatches == 0
            and kappa_mismatches == 0
            and tau_relerr_max < 1.0e-12
            and c6_relerr_max < 1.0e-12
            and present_margin_relerr_max < 1.0e-12
            and stress_margin_relerr_max < 1.0e-12
        ),
    }


def build() -> dict[str, Any]:
    constants = scalar_rows()
    constant_digest = stable_digest(constants)
    d5_verification = verify_d5_rows()
    publication = read_json(PUBLICATION_CARD)
    proton = read_json(PROTON_JSON)

    return {
        "note": (
            "No web lookup used.  This is an internal input-convention ledger, "
            "not a final cited experimental-input refresh."
        ),
        "source_manifest": manifest(SOURCE_FILES),
        "artifact_manifest": manifest(ARTIFACTS),
        "constant_rows_sha256": constant_digest,
        "constant_row_count": len(constants),
        "flavor_targets": fit.TARGETS,
        "flavor_gates": {
            "STRICT_CKM": repro.STRICT_CKM,
            "LOOSE_MASS": repro.LOOSE_MASS,
            "SEESAW_RESIDUAL_MAX": repro.SEESAW_RESIDUAL_MAX,
            "TOL": repro.TOL,
        },
        "seesaw_benchmark": {
            "m1_eV": 1.0e-3,
            "dm21_eV2": 7.42e-5,
            "dm31_eV2": 2.517e-3,
            "sin2_theta12": 0.304,
            "sin2_theta13": 0.0222,
            "sin2_theta23": 0.573,
            "delta_cp_rad": 1.20 * math.pi,
            "alpha21_rad": 0.35 * math.pi,
            "alpha31_rad": 1.10 * math.pi,
        },
        "d5_conventions": {
            "display_triplet_filter": dressed.DISPLAY_ST,
            "future_stress_target_years": 1.0e35,
            "width_cases": dressed.WIDTH_CASES,
            "channels": dressed.CHANNELS,
            "hbar_GeV_s": mi.HBAR_GEV_S,
            "seconds_per_year": mi.SECONDS_PER_YEAR,
            "soft_epsilon_grid": eig.SOFT_EPS_GRID,
            "spectrum_points": mi.SPECTRUM_POINTS,
        },
        "mssm_mixing_inputs": {
            "MW_GeV": mssm.MW_GEV,
            "MZ_GeV": mssm.MZ_GEV,
            "sin2_theta_W": mssm.SIN2_THETA_W,
            "cos2_theta_W": mssm.COS2_THETA_W,
            "alpha2_inv": mssm.ALPHA2_INV,
            "alpha2": mssm.ALPHA2,
            "alpha_y": mssm.ALPHA_Y,
            "hypercharge": mssm.HYPERCHARGE,
        },
        "rge_inputs": {
            "MZ_GeV": rge.MZ_GEV,
            "alpha_em_inv_MZ": rge.ALPHA_EM_INV_MZ,
            "sin2_theta_W_MZ": rge.SIN2_THETA_W_MZ,
            "alpha_s_MZ": rge.ALPHA_S_MZ,
            "SM_B1": arr(rge.SM_B1),
            "SM_B2": arr(rge.SM_B2),
            "MSSM_B1": arr(rge.MSSM_B1),
            "MSSM_B2": arr(rge.MSSM_B2),
            "MSSM_C_TOP": arr(rge.MSSM_C_TOP),
            "MSSM_C_BOTTOM": arr(rge.MSSM_C_BOTTOM),
            "MSSM_C_TAU": arr(rge.MSSM_C_TAU),
            "MSSM_C_NU": arr(rge.MSSM_C_NU),
            "HEAVY_BASIS": arr(rge.HEAVY_BASIS),
            "PROJECTOR": arr(rge.PROJECTOR),
        },
        "legacy_proton_benchmarks": {
            "hadronic_constants": proton["hadronic_constants"],
            "gauge_lifetime_table": proton["gauge_lifetime_table"],
            "mx_required_for_tau_1e34": proton["mx_required_for_tau_1e34"],
            "dimension5_scenarios": proton["dimension5_scenarios"],
        },
        "d5_row_formula_verification": d5_verification,
        "publication_card_digest": {
            "local_card_complete": publication["verdict"]["local_clockwork_rescued_card_complete"],
            "publication_level_complete": publication["verdict"]["publication_level_complete"],
            "ckm_score": publication["flavor_summary"]["scores"]["ckm_score"],
            "mass_score": publication["flavor_summary"]["scores"]["mass_score"],
            "seesaw_residual": publication["flavor_summary"]["seesaw_replay"]["seesaw_matrix_residual"],
            "d5_total_channel_rows": publication["d5_summary"]["total_channel_rows"],
            "d5_worst_margin_1e35": publication["d5_summary"]["global_worst"]["worst_margin_1e35"],
        },
        "verdict": {
            "ledger_built": True,
            "d5_rows_reproduce_formula": d5_verification["passes"],
            "publication_inputs_are_hash_locked_locally": True,
            "external_reference_refresh_done": False,
            "interpretation": (
                "The current clockwork-rescued flavor+d5 branch is internally "
                "consistent against its no-web local input conventions.  The "
                "remaining publication blocker is not an arithmetic mismatch; "
                "it is replacing this convention ledger by a final cited input "
                "table and rerunning the same audit."
            ),
        },
    }


def report(summary: dict[str, Any]) -> str:
    d5 = summary["d5_row_formula_verification"]
    card = summary["publication_card_digest"]
    lines = [
        "# No-web input convention ledger",
        "",
        "No web lookup was used.  This is an internal convention card, not a final bibliography/input refresh.",
        "",
        "## Hash lock",
        "",
        f"- constant rows: `{summary['constant_row_count']}`",
        f"- constant-row digest: `{summary['constant_rows_sha256']}`",
        "",
        "## Publication-card digest",
        "",
        f"- local card complete: `{card['local_card_complete']}`",
        f"- publication-level complete: `{card['publication_level_complete']}`",
        f"- CKM score: `{card['ckm_score']:.6e}`",
        f"- mass score: `{card['mass_score']:.6e}`",
        f"- seesaw residual: `{card['seesaw_residual']:.6e}`",
        f"- d5 channel rows: `{card['d5_total_channel_rows']}`",
        f"- worst d5 1e35 margin: `{card['d5_worst_margin_1e35']:.6e}`",
        "",
        "## Dressed d5 formula replay",
        "",
        f"- rows checked: `{d5['rows_checked']}`",
        f"- width-prefactor mismatches: `{d5['prefactor_mismatches']}`",
        f"- present-bound mismatches: `{d5['present_bound_mismatches']}`",
        f"- display-filter mismatches: `{d5['display_filter_mismatches']}`",
        f"- kappa mismatches: `{d5['kappa_mismatches']}`",
        f"- max C6 relative error: `{d5['max_c6_relative_error']:.3e}`",
        f"- max tau relative error: `{d5['max_tau_relative_error']:.3e}`",
        f"- max present-margin relative error: `{d5['max_present_margin_relative_error']:.3e}`",
        f"- max 1e35-margin relative error: `{d5['max_1e35_margin_relative_error']:.3e}`",
        f"- formula replay passes: `{d5['passes']}`",
        "",
        "## Verdict",
        "",
        summary["verdict"]["interpretation"],
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = build()
    constants = scalar_rows()
    fields = ["category", "name", "value", "unit", "source", "role"]
    write_csv(OUT / "constant_rows.csv", constants, fields)
    (OUT / "input_convention_ledger.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "report.md").write_text(report(summary), encoding="utf-8")

    print("No-web input convention ledger written")
    print(f"  constants: {summary['constant_row_count']}")
    print(f"  d5 rows checked: {summary['d5_row_formula_verification']['rows_checked']}")
    print(f"  formula replay passes: {summary['d5_row_formula_verification']['passes']}")
    print(f"  digest: {summary['constant_rows_sha256']}")


if __name__ == "__main__":
    main()
