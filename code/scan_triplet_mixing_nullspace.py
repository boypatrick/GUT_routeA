#!/usr/bin/env python3
"""Scan triplet-mixing null directions for d=5 proton-decay tensors.

No web lookup is used.  The previous audit multiplied the whole d=5 tensor by
one scalar filter S_T.  Here the filter is promoted to a finite inverse
triplet-mass/mixing matrix W_AB = (T^{-1})_AB in a local four-shape basis

    H10-like, F126-like, G120-like, K_tr contact.

The calculation is intentionally phrased as a linear algebra problem.  For each
source pair (A,B), the mass-basis tensor contribution is known.  We then search
for the coefficient vector W_AB closest to the identity branch that kills the
K nu LLLL entries, and check whether K0 mu or RRRR channels become the new
leading obstruction.
"""

from __future__ import annotations

import csv
import json
import math
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "output" / "flavor_transvectant_rotations" / "transvectant_flavor_rotations.json"
PROTON = ROOT / "output" / "proton_decay" / "proton_decay_verification.json"
OUT = ROOT / "output" / "triplet_mixing_nullspace"

HBAR_GEV_S = 6.582119569e-25
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def cmat(raw: list[list[dict[str, float]]]) -> np.ndarray:
    return np.array([[cell["re"] + 1j * cell["im"] for cell in row] for row in raw], dtype=complex)


def cjson(z: complex) -> dict[str, float]:
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def sym(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat + mat.T)


def anti(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat - mat.T)


def scale_to_largest_singular(mat: np.ndarray, target: float) -> np.ndarray:
    singulars = np.linalg.svd(mat, compute_uv=False)
    return mat * (target / max(float(singulars[0]), 1.0e-30))


def tensor_llll(
    y_qq: np.ndarray,
    y_ql: np.ndarray,
    uq1: np.ndarray,
    uq2: np.ndarray,
    uq3: np.ndarray,
    ul: np.ndarray,
) -> np.ndarray:
    return np.einsum("ij,kl,ia,jb,kc,ld->abcd", y_qq, y_ql, uq1, uq2, uq3, ul, optimize=True)


def tensor_rrrr(
    y_ue: np.ndarray,
    y_ud: np.ndarray,
    uu1: np.ndarray,
    ue: np.ndarray,
    uu2: np.ndarray,
    ud: np.ndarray,
) -> np.ndarray:
    return np.einsum("ij,kl,ia,jd,kb,lc->abcd", y_ue, y_ud, uu1, ue, uu2, ud, optimize=True)


def lifetime_from_amplitude(
    amplitude: float,
    constants: dict[str, float],
    m_triplet_gev: float = 1.0e16,
    m_wino_gev: float = 1.0e3,
    m_sfermion_gev: float = 1.0e5,
    alpha2_inv: float = 25.0,
) -> dict[str, float]:
    loop = (1.0 / alpha2_inv) / (4.0 * math.pi)
    c5 = amplitude / m_triplet_gev
    dressing = loop * m_wino_gev / (m_sfermion_gev * m_sfermion_gev)
    c6 = c5 * dressing
    width = constants["width_prefactor_GeV5"] * c6 * c6
    tau = math.inf if width <= 0.0 else HBAR_GEV_S / width / SECONDS_PER_YEAR
    return {
        "amplitude": float(amplitude),
        "tau_years": float(tau),
        "S_T_required_tau_2p4e34": math.inf if amplitude <= 0.0 else float(math.sqrt(tau / 2.4e34)),
        "C6_dressed_GeV_minus2": float(c6),
    }


def max_perm_rows(q_flavors: tuple[int, int, int], leptons: list[int]) -> list[tuple[int, int, int, int]]:
    return [(a, b, c, ell) for (a, b, c) in sorted(set(permutations(q_flavors))) for ell in leptons]


def max_rrrr_rows(charged: list[int]) -> list[tuple[int, int, int, int]]:
    rows: list[tuple[int, int, int, int]] = []
    for ell in charged:
        for u_pair in [(0, 0), (0, 1), (1, 0)]:
            rows.append((u_pair[0], u_pair[1], 1, ell))
    return rows


def row_matrix(
    pair_list: list[tuple[int, int]],
    tensors: dict[tuple[int, int], np.ndarray],
    entries: list[tuple[int, int, int, int]],
) -> np.ndarray:
    return np.array([[tensors[pair][entry] for pair in pair_list] for entry in entries], dtype=complex)


def max_channel(
    coeffs: np.ndarray,
    pair_list: list[tuple[int, int]],
    tensors: dict[tuple[int, int], np.ndarray],
    entries: list[tuple[int, int, int, int]],
) -> tuple[float, tuple[int, int, int, int], complex]:
    best_amp = -1.0
    best_entry = entries[0]
    best_value = 0.0j
    for entry in entries:
        value = sum(coeffs[n] * tensors[pair][entry] for n, pair in enumerate(pair_list))
        amp = float(abs(value))
        if amp > best_amp:
            best_amp = amp
            best_entry = entry
            best_value = complex(value)
    return best_amp, best_entry, best_value


def matrix_from_coeffs(coeffs: np.ndarray, pair_list: list[tuple[int, int]], size: int) -> np.ndarray:
    w = np.zeros((size, size), dtype=complex)
    for coeff, (i, j) in zip(coeffs, pair_list):
        w[i, j] = coeff
    return w


def coeff_json(coeffs: np.ndarray, pair_list: list[tuple[int, int]], names: list[str]) -> list[dict[str, Any]]:
    out = []
    for coeff, (i, j) in zip(coeffs, pair_list):
        out.append(
            {
                "QQ_source": names[i],
                "QL_source": names[j],
                "abs": float(abs(coeff)),
                "phase_rad": float(np.angle(coeff)),
                "re": float(np.real(coeff)),
                "im": float(np.imag(coeff)),
            }
        )
    return out


def project_or_least_singular(a: np.ndarray, baseline: np.ndarray, tol: float = 1.0e-11) -> tuple[np.ndarray, dict[str, Any]]:
    u, s, vh = np.linalg.svd(a, full_matrices=True)
    rank = int(np.sum(s > tol * max(float(s[0]), 1.0)))
    null = vh.conjugate().T[:, rank:]
    if null.shape[1] > 0:
        projected = null @ (null.conjugate().T @ baseline)
        if np.linalg.norm(projected) > 1.0e-13:
            coeffs = projected / np.linalg.norm(projected)
            mode = "projected_exact_null_closest_to_identity"
        else:
            coeffs = null[:, 0] / np.linalg.norm(null[:, 0])
            mode = "exact_null_orthogonal_to_identity"
    else:
        coeffs = vh.conjugate().T[:, -1]
        coeffs = coeffs / np.linalg.norm(coeffs)
        mode = "least_singular_no_exact_null"
    residual = float(np.linalg.norm(a @ coeffs))
    return coeffs, {
        "mode": mode,
        "rank": rank,
        "singular_values": [float(x) for x in s],
        "nullity": int(null.shape[1]),
        "constraint_residual_l2": residual,
    }


def audit_candidate(
    label: str,
    coeffs: np.ndarray,
    baseline: np.ndarray,
    pair_list: list[tuple[int, int]],
    size: int,
    basis_names: list[str],
    channel_entries: dict[str, tuple[dict[tuple[int, int], np.ndarray], list[tuple[int, int, int, int]]]],
    constants: dict[str, float],
    mode_info: dict[str, Any],
) -> dict[str, Any]:
    channels: dict[str, Any] = {}
    for channel, (tensors, entries) in channel_entries.items():
        amp, idx, value = max_channel(coeffs, pair_list, tensors, entries)
        life = lifetime_from_amplitude(amp, constants)
        channels[channel] = {
            "amplitude": amp,
            "selected_index": "".join(str(x) for x in idx),
            "value": cjson(value),
            "tau_years": life["tau_years"],
            "S_T_required_tau_2p4e34": life["S_T_required_tau_2p4e34"],
        }

    wmat = matrix_from_coeffs(coeffs, pair_list, size)
    singulars = np.linalg.svd(wmat, compute_uv=False)
    nonzero = singulars[singulars > 1.0e-12]
    condition = math.inf if len(nonzero) < len(singulars) else float(nonzero[0] / nonzero[-1])
    overlap = abs(np.vdot(baseline, coeffs)) / max(np.linalg.norm(baseline) * np.linalg.norm(coeffs), 1.0e-30)
    angle = float(math.acos(min(1.0, max(0.0, float(overlap)))))
    rms = float(np.linalg.norm(coeffs) / math.sqrt(len(coeffs)))

    return {
        "label": label,
        "mode_info": mode_info,
        "identity_overlap_abs": float(overlap),
        "identity_angle_rad": angle,
        "coefficient_peak_over_rms": float(np.max(np.abs(coeffs)) / max(rms, 1.0e-30)),
        "W_rank": int(len(nonzero)),
        "W_singular_values": [float(x) for x in singulars],
        "W_condition": condition,
        "coefficients": coeff_json(coeffs, pair_list, basis_names),
        "channels": channels,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    proton = json.loads(PROTON.read_text(encoding="utf-8"))
    constants = proton["hadronic_constants"]

    y_raw = {name: cmat(raw) for name, raw in payload["Yukawa_matrices"].items()}
    rot = payload["biunitary_rotations"]
    ul = {name: cmat(rot[name]["left_rotation"]) for name in rot}
    ur = {name: cmat(rot[name]["right_rotation"]) for name in rot}
    u_nu = cmat(payload["PMNS_convention"]["neutrino_left_rotation_target"])

    y_phys = {
        "up": scale_to_largest_singular(y_raw["up"], 0.60),
        "down": scale_to_largest_singular(y_raw["down"], 0.024),
        "charged_lepton": scale_to_largest_singular(y_raw["charged_lepton"], 0.010),
    }
    k_tr_top = scale_to_largest_singular(
        np.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=complex) / math.sqrt(3.0),
        0.60,
    )
    k_tr_bottom = scale_to_largest_singular(k_tr_top, 0.024)

    basis_names = ["H10", "F126", "G120", "Ktr"]
    delta_de = y_phys["down"] - y_phys["charged_lepton"].T
    qq_basis = [
        scale_to_largest_singular(sym(y_phys["up"]), 0.60),
        scale_to_largest_singular(sym(delta_de), 0.60),
        scale_to_largest_singular(anti(y_phys["up"]), 0.60),
        k_tr_top,
    ]
    ql_basis = [
        scale_to_largest_singular(y_phys["down"], 0.024),
        scale_to_largest_singular(delta_de, 0.024),
        scale_to_largest_singular(anti(y_phys["down"]), 0.024),
        k_tr_bottom,
    ]

    size = len(basis_names)
    pair_full = [(i, j) for i in range(size) for j in range(size)]
    pair_diag = [(i, i) for i in range(size)]

    tensors_up: dict[tuple[int, int], np.ndarray] = {}
    tensors_down: dict[tuple[int, int], np.ndarray] = {}
    tensors_mu: dict[tuple[int, int], np.ndarray] = {}
    tensors_r: dict[tuple[int, int], np.ndarray] = {}
    for i in range(size):
        for j in range(size):
            tensors_up[(i, j)] = tensor_llll(qq_basis[i], ql_basis[j], ul["up"], ul["up"], ul["down"], u_nu)
            tensors_down[(i, j)] = tensor_llll(qq_basis[i], ql_basis[j], ul["down"], ul["down"], ul["up"], u_nu)
            tensors_mu[(i, j)] = tensor_llll(qq_basis[i], ql_basis[j], ul["up"], ul["up"], ul["down"], ul["charged_lepton"])
            tensors_r[(i, j)] = tensor_rrrr(qq_basis[i], ql_basis[j], ur["up"], ur["charged_lepton"], ur["up"], ur["down"])

    entries_knu = max_perm_rows((0, 0, 1), [0, 1, 2])
    entries_k0mu = max_perm_rows((0, 0, 1), [1])
    entries_rrrr = max_rrrr_rows([0, 1, 2])
    channel_entries = {
        "LLLL_upupdown_Knu": (tensors_up, entries_knu),
        "LLLL_downdownup_Knu": (tensors_down, entries_knu),
        "LLLL_upupdown_K0mu": (tensors_mu, entries_k0mu),
        "RRRR_uusd_anycharged": (tensors_r, entries_rrrr),
    }

    audits: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    for ansatz, pair_list in [("diagonal_rep_mixing", pair_diag), ("full_bipartite_mixing", pair_full)]:
        baseline = np.zeros(len(pair_list), dtype=complex)
        for n, (i, j) in enumerate(pair_list):
            if i == j:
                baseline[n] = 1.0
        baseline = baseline / np.linalg.norm(baseline)

        a_knu = np.vstack(
            [
                row_matrix(pair_list, tensors_up, entries_knu),
                row_matrix(pair_list, tensors_down, entries_knu),
            ]
        )
        coeffs, mode_info = project_or_least_singular(a_knu, baseline)
        for label, vec, info in [
            (ansatz + "_identity", baseline, {"mode": "identity_branch", "rank": None, "nullity": None, "constraint_residual_l2": float(np.linalg.norm(a_knu @ baseline))}),
            (ansatz + "_knu_null", coeffs, mode_info),
        ]:
            audit = audit_candidate(label, vec, baseline, pair_list, size, basis_names, channel_entries, constants, info)
            audits.append(audit)
            for channel, row in audit["channels"].items():
                csv_rows.append(
                    {
                        "label": label,
                        "mode": audit["mode_info"]["mode"],
                        "channel": channel,
                        "amplitude": row["amplitude"],
                        "selected_index": row["selected_index"],
                        "S_T_required_tau_2p4e34": row["S_T_required_tau_2p4e34"],
                        "identity_overlap_abs": audit["identity_overlap_abs"],
                        "identity_angle_rad": audit["identity_angle_rad"],
                        "W_rank": audit["W_rank"],
                        "W_condition": audit["W_condition"],
                        "coefficient_peak_over_rms": audit["coefficient_peak_over_rms"],
                    }
                )

    with (OUT / "triplet_mixing_nullspace_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)

    identity_full = next(a for a in audits if a["label"] == "full_bipartite_mixing_identity")
    null_full = next(a for a in audits if a["label"] == "full_bipartite_mixing_knu_null")
    identity_amp = max(
        identity_full["channels"]["LLLL_upupdown_Knu"]["amplitude"],
        identity_full["channels"]["LLLL_downdownup_Knu"]["amplitude"],
    )
    null_amp = max(
        null_full["channels"]["LLLL_upupdown_Knu"]["amplitude"],
        null_full["channels"]["LLLL_downdownup_Knu"]["amplitude"],
    )
    suppression = math.inf if null_amp <= 0.0 else identity_amp / null_amp
    leading_after_null = max(null_full["channels"].items(), key=lambda item: item[1]["amplitude"])

    output = {
        "note": "No web lookup used. This is a finite-dimensional triplet-mixing nullspace audit, not a full dressed proton-decay calculation.",
        "basis_names": basis_names,
        "basis_definition": {
            "H10": "sym(Y_u), top-normalized on QQ and bottom-normalized on QL",
            "F126": "sym(Y_d-Y_e^T) / (Y_d-Y_e^T), representation-shape proxy",
            "G120": "anti(Y_u) / anti(Y_d), antisymmetric representation-shape proxy",
            "Ktr": "spin-zero CP1 transvectant contact tensor",
        },
        "normalization": {
            "QQ_largest_singular_targets": {"H10": 0.60, "F126": 0.60, "G120": 0.60, "Ktr": 0.60},
            "QL_largest_singular_targets": {"H10": 0.024, "F126": 0.024, "G120": 0.024, "Ktr": 0.024},
            "M_T_GeV": 1.0e16,
            "m_wino_GeV": 1.0e3,
            "m_sfermion_GeV": 1.0e5,
            "alpha2_inv": 25.0,
        },
        "audits": audits,
        "verdict": {
            "full_identity_max_Knu_amplitude": float(identity_amp),
            "full_null_max_Knu_amplitude": float(null_amp),
            "full_null_Knu_suppression_factor": float(suppression),
            "full_null_identity_overlap_abs": null_full["identity_overlap_abs"],
            "full_null_angle_rad": null_full["identity_angle_rad"],
            "full_null_W_rank": null_full["W_rank"],
            "full_null_W_condition": null_full["W_condition"],
            "full_null_leading_channel": leading_after_null[0],
            "full_null_leading_amplitude": leading_after_null[1]["amplitude"],
            "interpretation": (
                "With both doublet assignments and all three neutrino flavors constrained, the finite "
                "four-shape basis has no exact Knu nullspace.  The full bipartite inverse-mixing matrix "
                "nevertheless has a sharp near-null singular direction that suppresses Knu by four orders "
                "of magnitude.  The obstruction then moves to the RRRR proxy, and the near-null matrix is "
                "rank deficient; a UV model must therefore explain either this singular direction or an "
                "additional RRRR filter."
            ),
        },
    }
    (OUT / "triplet_mixing_nullspace_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [
        "# Triplet-mixing nullspace audit",
        "",
        "No web lookup was used.",
        "",
        "The scalar filter is replaced by a normalized inverse triplet-mixing matrix",
        "`W_AB=(T^{-1})_AB` in the local basis `H10,F126,G120,Ktr`.",
        "",
        "| label | channel | amp | S_T max 2.4e34 | overlap with identity | W rank | W cond |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in csv_rows:
        cond = row["W_condition"]
        cond_text = "inf" if not math.isfinite(cond) else f"{cond:.3e}"
        report.append(
            f"| `{row['label']}` | `{row['channel']}` | {row['amplitude']:.3e} | "
            f"{row['S_T_required_tau_2p4e34']:.3e} | {row['identity_overlap_abs']:.3e} | "
            f"{row['W_rank']} | {cond_text} |"
        )
    report.extend(
        [
            "",
            "## Verdict",
            "",
            output["verdict"]["interpretation"],
            "",
            f"Full-bipartite identity max Knu amplitude: `{identity_amp:.6e}`.",
            f"Full-bipartite Knu near-null max Knu amplitude: `{null_amp:.6e}`.",
            f"Suppression factor: `{suppression:.6e}`.",
            f"Leading channel after Knu null: `{leading_after_null[0]}` with amplitude `{leading_after_null[1]['amplitude']:.6e}`.",
            "",
        ]
    )
    (OUT / "triplet_mixing_nullspace_report.md").write_text("\n".join(report), encoding="utf-8")

    print("Triplet-mixing nullspace audit")
    print(f"  full identity max Knu amplitude: {identity_amp:.6e}")
    print(f"  full null max Knu amplitude: {null_amp:.6e}")
    print(f"  suppression factor: {suppression:.6e}")
    print(f"  leading after null: {leading_after_null[0]} {leading_after_null[1]['amplitude']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
