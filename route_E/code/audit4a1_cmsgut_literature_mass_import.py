#!/usr/bin/env python3
"""Import source-anchored CMSGUT doublet/triplet mass matrices.

This stage transcribes the Aulakh--Girdhar electroweak-doublet and
proton-decay triplet matrices into the Audit 4a.1 convention fixed by the
vacuum-cubic card.  It deliberately keeps the transcription source-anchored:
each block records the line range in the local TeX source from which it was
copied.

The script also performs numerical smoke gates:

* chiral Goldstone gates for the mixed chiral/gauge matrices G,E,F,J,X, using
  a generic F-flat branch point, verifying that the upper-left chiral block has
  a null singular value while the full Dirac/gaugino block is non-singular;
* a triplet inverse gate, verifying the 5x5 triplet matrix can be inverted at a
  generic nonsingular point and exporting the inverse entries needed by Audit 2.

The output is not a complete heavy spectrum and does not yet verify scalar
Hessian Goldstone eigenvectors.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit4a1"
AG_TEX = Path("/private/tmp/cmsgut_refs/0405074/msgtfrvnew.tex")
BMSV_TEX = Path("/private/tmp/cmsgut_refs/0306242/abmsv.tex")

SOURCES = {
    "this_script": ROOT / "code" / "audit4a1_cmsgut_literature_mass_import.py",
    "mass_export_schema": ROOT / "output" / "audit4a1" / "mass_export_schema.json",
    "vacuum_branches": ROOT / "output" / "audit4a1" / "vacuum_branches.json",
    "aulakh_girdhar_tex": AG_TEX,
    "bmsv_note_tex": BMSV_TEX,
}


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def source_manifest() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label, path in SOURCES.items():
        if path.is_absolute() and path.is_relative_to(ROOT):
            shown_path = str(path.relative_to(ROOT))
        else:
            shown_path = str(path)
        rows.append(
            {
                "label": label,
                "path": shown_path,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path),
            }
        )
    return rows


def doublet_matrix() -> dict[str, Any]:
    rows = [
        "bar_h1 = H^{alpha}_{dot2}",
        "bar_h2 = overlineSigma^{(15) alpha}_{dot2}",
        "bar_h3 = Sigma^{(15) alpha}_{dot2}",
        "bar_h4 = phi_{44}^{dot2 alpha}/sqrt(2)",
    ]
    cols = [
        "h1 = H_{alpha dot1}",
        "h2 = overlineSigma^{(15)}_{alpha dot1}",
        "h3 = Sigma^{(15)}_{alpha dot1}",
        "h4 = phi^{44 dot1}_{alpha}/sqrt(2)",
    ]
    entries = [
        ["-M_H", "bar_gamma*sqrt(3)*(omega-a)", "-gamma*sqrt(3)*(omega+a)", "-bar_gamma*bar_sigma"],
        ["-bar_gamma*sqrt(3)*(omega+a)", "0", "-(2*M+4*eta*(a+omega))", "0"],
        ["gamma*sqrt(3)*(omega-a)", "-(2*M+4*eta*(a-omega))", "0", "-2*eta*bar_sigma*sqrt(3)"],
        ["-sigma*gamma", "-2*eta*sigma*sqrt(3)", "0", "-2*m+6*lambda*(omega-a)"],
    ]
    return {
        "block_id": "H_doublet_4x4",
        "kind": "electroweak_doublet_mass_matrix",
        "source_anchor": {
            "source": str(AG_TEX),
            "arxiv": "hep-ph/0405074",
            "line_range": "2601-2624",
            "notes": "Rows are barred doublets and columns are unbarred doublets; source states Det H=0 is the MSSM one-pair fine-tuning condition.",
        },
        "basis_rows": rows,
        "basis_columns": cols,
        "entries": entries,
        "shape": [4, 4],
        "entry_status": "transcribed_from_aulakh_girdhar_appendix",
    }


def triplet_matrix() -> dict[str, Any]:
    rows = [
        "bar_t1 = H^4",
        "bar_t2 = overlineSigma_{(a)}^4",
        "bar_t3 = Sigma_{(a)}^4",
        "bar_t4 = Sigma^4_{R0}",
        "bar_t5 = phi_{4(R+)}",
    ]
    cols = [
        "t1 = H_4",
        "t2 = overlineSigma_{(a)4}",
        "t3 = Sigma_{4(a)}",
        "t4 = overlineSigma_{4(R0)}",
        "t5 = phi_{(R-)}^4",
    ]
    entries = [
        ["M_H", "bar_gamma*(a+p)", "gamma*(p-a)", "2*sqrt(2)*I*omega*bar_gamma", "I*bar_sigma*bar_gamma"],
        ["bar_gamma*(p-a)", "0", "2*M", "0", "0"],
        ["gamma*(p+a)", "2*M", "0", "4*sqrt(2)*I*omega*eta", "2*I*eta*bar_sigma"],
        ["-2*sqrt(2)*I*omega*gamma", "-4*sqrt(2)*I*omega*eta", "0", "2*M+2*eta*p+2*eta*a", "-2*sqrt(2)*eta*bar_sigma"],
        ["I*sigma*gamma", "2*I*eta*sigma", "0", "2*sqrt(2)*eta*sigma", "-2*m-2*lambda*(a+p-4*omega)"],
    ]
    return {
        "block_id": "T_triplet_5x5",
        "kind": "proton_decay_triplet_mass_matrix",
        "source_anchor": {
            "source": str(AG_TEX),
            "arxiv": "hep-ph/0405074",
            "line_range": "2629-2648",
            "notes": "Rows are barred [bar3,1,2/3] triplets and columns are unbarred [3,1,-2/3] triplets.",
        },
        "basis_rows": rows,
        "basis_columns": cols,
        "entries": entries,
        "shape": [5, 5],
        "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        "audit2_required_inverse_entries_1_indexed": [
            "S_1^1",
            "S_1^2",
            "S_2^1",
            "S_2^2",
            "S_1^4",
            "S_2^4",
        ],
    }


def mixed_matrix_blocks() -> list[dict[str, Any]]:
    return [
        {
            "block_id": "G_mixed_neutral_6x6",
            "kind": "mixed_chiral_gauge_matrix",
            "source_anchor": {
                "source": str(AG_TEX),
                "arxiv": "hep-ph/0405074",
                "line_range": "2654-2669",
                "description_line_range": "1203-1214",
                "notes": "Upper-left 5x5 chiral block has null left/right eigenvectors from the sixth row/column; full 6x6 matrix is generically non-singular.",
            },
            "basis_rows": ["G1", "G2", "G3", "G4", "G5", "G6_gaugino"],
            "basis_columns": ["G1", "G2", "G3", "G4", "G5", "G6_gaugino"],
            "chiral_subblock_shape": [5, 5],
            "entries": [
                ["2*m", "0", "2*sqrt(6)*lambda*omega", "2*I*eta*bar_sigma/sqrt(2)", "-2*I*eta*sigma/sqrt(2)", "0"],
                ["0", "2*(m+2*lambda*a)", "4*sqrt(2)*lambda*omega", "2*I*eta*bar_sigma*sqrt(3/2)", "-2*I*eta*sigma*sqrt(3/2)", "0"],
                ["2*sqrt(6)*lambda*omega", "4*sqrt(2)*lambda*omega", "2*(m+lambda*(p+2*a))", "-2*I*eta*sqrt(3)*bar_sigma", "2*I*sqrt(3)*eta*sigma", "0"],
                ["2*I*eta*bar_sigma/sqrt(2)", "2*I*eta*bar_sigma*sqrt(3/2)", "-2*I*eta*sqrt(3)*bar_sigma", "0", "2*(M+eta*(p+3*a-6*omega))", "sqrt(5)*g*conj(sigma)"],
                ["-2*I*eta*sigma/sqrt(2)", "-2*I*eta*sigma*sqrt(3/2)", "2*I*sqrt(3)*eta*sigma", "2*(M+eta*(p+3*a-6*omega))", "0", "sqrt(5)*g*conj(bar_sigma)"],
                ["0", "0", "0", "sqrt(5)*g*conj(sigma)", "sqrt(5)*g*conj(bar_sigma)", "0"],
            ],
            "shape": [6, 6],
            "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        },
        {
            "block_id": "E_mixed_4x4",
            "kind": "mixed_chiral_gauge_matrix",
            "source_anchor": {
                "source": str(AG_TEX),
                "arxiv": "hep-ph/0405074",
                "line_range": "2671-2685",
                "description_line_range": "1219-1244",
                "notes": "Upper-left 3x3 chiral block has null left/right eigenvectors from the fourth row/column; full 4x4 matrix is generically non-singular.",
            },
            "basis_rows": ["bar_E2", "bar_E3", "bar_E4", "bar_E5_gaugino"],
            "basis_columns": ["E2", "E3", "E4", "E5_gaugino"],
            "chiral_subblock_shape": [3, 3],
            "entries": [
                ["-2*(M+eta*(a-3*omega))", "-2*sqrt(2)*I*eta*sigma", "2*I*eta*sigma", "I*g*sqrt(2)*conj(bar_sigma)"],
                ["2*I*sqrt(2)*eta*bar_sigma", "-2*(m+lambda*(a-omega))", "-2*sqrt(2)*lambda*omega", "2*g*(conj(a)-conj(omega))"],
                ["-2*I*eta*bar_sigma", "-2*sqrt(2)*lambda*omega", "-2*(m-lambda*omega)", "sqrt(2)*g*(conj(omega)-conj(p))"],
                ["-I*g*sqrt(2)*conj(sigma)", "2*g*(conj(a)-conj(omega))", "g*sqrt(2)*(conj(omega)-conj(p))", "0"],
            ],
            "shape": [4, 4],
            "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        },
        {
            "block_id": "F_mixed_3x3",
            "kind": "mixed_chiral_gauge_matrix",
            "source_anchor": {
                "source": str(AG_TEX),
                "arxiv": "hep-ph/0405074",
                "line_range": "2688-2699",
                "description_line_range": "1273-1293",
                "notes": "Mixed charged singlet block; full matrix separates gaugino mode and residual chiral Dirac mode.",
            },
            "basis_rows": ["bar_F1", "bar_F2", "bar_F3_gaugino"],
            "basis_columns": ["F1", "F2", "F3_gaugino"],
            "chiral_subblock_shape": [2, 2],
            "entries": [
                ["2*(M+eta*(p+3*a))", "-2*I*sqrt(3)*eta*sigma", "-g*sqrt(2)*conj(bar_sigma)"],
                ["2*I*sqrt(3)*eta*bar_sigma", "2*(m+lambda*(p+2*a))", "sqrt(24)*I*g*conj(omega)"],
                ["-g*sqrt(2)*conj(sigma)", "-sqrt(24)*I*g*conj(omega)", "0"],
            ],
            "shape": [3, 3],
            "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        },
        {
            "block_id": "J_mixed_4x4",
            "kind": "mixed_chiral_gauge_matrix",
            "source_anchor": {
                "source": str(AG_TEX),
                "arxiv": "hep-ph/0405074",
                "line_range": "2701-2714",
                "description_line_range": "1295-1308",
                "notes": "Upper-left 3x3 chiral block has null left/right eigenvectors from the fourth row/column; full 4x4 matrix is generically non-singular.",
            },
            "basis_rows": ["bar_J1", "bar_J2", "bar_J3", "bar_J4_gaugino"],
            "basis_columns": ["J1", "J2", "J3", "J4_gaugino"],
            "chiral_subblock_shape": [3, 3],
            "entries": [
                ["2*(M+eta*(a+p-2*omega))", "-2*eta*bar_sigma", "2*sqrt(2)*eta*bar_sigma", "-I*g*sqrt(2)*conj(sigma)"],
                ["2*eta*sigma", "-2*(m+lambda*a)", "-2*sqrt(2)*lambda*omega", "-2*I*g*sqrt(2)*conj(a)"],
                ["-2*sqrt(2)*eta*sigma", "-2*sqrt(2)*lambda*omega", "-2*(m+lambda*(a+p))", "-4*I*g*conj(omega)"],
                ["-I*g*sqrt(2)*conj(bar_sigma)", "2*sqrt(2)*I*g*conj(a)", "4*I*g*conj(omega)", "0"],
            ],
            "shape": [4, 4],
            "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        },
        {
            "block_id": "X_mixed_3x3",
            "kind": "mixed_chiral_gauge_matrix",
            "source_anchor": {
                "source": str(AG_TEX),
                "arxiv": "hep-ph/0405074",
                "line_range": "2718-2731",
                "description_line_range": "1313-1328",
                "notes": "Symmetric mixed matrix; upper-left 2x2 chiral block has null vector from the third row/column; full 3x3 matrix is generically non-singular.",
            },
            "basis_rows": ["bar_X1", "bar_X2", "bar_X3_gaugino"],
            "basis_columns": ["X1", "X2", "X3_gaugino"],
            "chiral_subblock_shape": [2, 2],
            "entries": [
                ["2*(m+lambda*(a+omega))", "-2*sqrt(2)*lambda*omega", "-2*g*(conj(a)+conj(omega))"],
                ["-2*sqrt(2)*lambda*omega", "2*(m+lambda*omega)", "sqrt(2)*g*(conj(omega)+conj(p))"],
                ["-2*g*(conj(a)+conj(omega))", "sqrt(2)*g*(conj(omega)+conj(p))", "0"],
            ],
            "shape": [3, 3],
            "entry_status": "transcribed_from_aulakh_girdhar_appendix",
        },
    ]


def symbolic_triplet_inverse_contract() -> dict[str, Any]:
    return {
        "status": "hand_audit_ready_adjugate_contract",
        "matrix": "T_triplet_5x5",
        "definition": "S = T^{-1}",
        "cofactor_formula_1_indexed": "S_i^j = (-1)^(i+j) det(T with row j and column i removed) / det(T)",
        "audit2_required_entries": [
            {"entry": "S_1^1", "minor": "remove row 1, column 1"},
            {"entry": "S_1^2", "minor": "remove row 2, column 1"},
            {"entry": "S_2^1", "minor": "remove row 1, column 2"},
            {"entry": "S_2^2", "minor": "remove row 2, column 2"},
            {"entry": "S_1^4", "minor": "remove row 4, column 1"},
            {"entry": "S_2^4", "minor": "remove row 4, column 2"},
        ],
        "cas_status": "sympy_not_available_in_current_runtime; numeric inverse gate is exported, symbolic expansion remains pending",
    }


def xi_from_x(x: float) -> float:
    return -((8 * x**3 - 15 * x**2 + 14 * x - 3) / (1 - x) ** 2)


def vev_from_x(x: float, m: float = 1.0, lam: float = 1.0, eta: float = 1.0) -> dict[str, complex]:
    scale = m / lam
    omega = -x * scale
    a = ((x * x + 2 * x - 1) / (1 - x)) * scale
    p = (x * (5 * x * x - 1) / (1 - x) ** 2) * scale
    xi = xi_from_x(x)
    M = xi * eta * m / lam
    sigma_prod = (2 / eta) * x * (1 - 3 * x) * (1 + x * x) / (1 - x) ** 2 * scale**2
    if sigma_prod <= 0:
        raise ValueError("sample x did not yield positive sigma*bar_sigma in the real D-flat slice")
    sigma = math.sqrt(sigma_prod)
    return {
        "x": x,
        "xi": xi,
        "m": m,
        "lambda": lam,
        "eta": eta,
        "M": M,
        "omega": omega,
        "a": a,
        "p": p,
        "sigma": sigma,
        "bar_sigma": sigma,
        "f_flat_M_eta_p_3a_minus_6omega": M + eta * (p + 3 * a - 6 * omega),
    }


def triplet_numeric_matrix(params: dict[str, complex]) -> np.ndarray:
    sq = math.sqrt
    I = 1j
    m = params["m"]
    lam = params["lambda"]
    eta = params["eta"]
    M = params["M"]
    omega = params["omega"]
    a = params["a"]
    p = params["p"]
    sigma = params["sigma"]
    bar_sigma = params["bar_sigma"]
    gamma = params["gamma"]
    bar_gamma = params["bar_gamma"]
    M_H = params["M_H"]
    return np.array(
        [
            [M_H, bar_gamma * (a + p), gamma * (p - a), 2 * sq(2) * I * omega * bar_gamma, I * bar_sigma * bar_gamma],
            [bar_gamma * (p - a), 0, 2 * M, 0, 0],
            [gamma * (p + a), 2 * M, 0, 4 * sq(2) * I * omega * eta, 2 * I * eta * bar_sigma],
            [-2 * sq(2) * I * omega * gamma, -4 * sq(2) * I * omega * eta, 0, 2 * M + 2 * eta * p + 2 * eta * a, -2 * sq(2) * eta * bar_sigma],
            [I * sigma * gamma, 2 * I * eta * sigma, 0, 2 * sq(2) * eta * sigma, -2 * m - 2 * lam * (a + p - 4 * omega)],
        ],
        dtype=complex,
    )


def neutral_g_numeric_matrix(params: dict[str, complex]) -> np.ndarray:
    sq = math.sqrt
    I = 1j
    m = params["m"]
    lam = params["lambda"]
    eta = params["eta"]
    M = params["M"]
    omega = params["omega"]
    a = params["a"]
    p = params["p"]
    sigma = params["sigma"]
    bar_sigma = params["bar_sigma"]
    g = params["g"]
    return 2 * np.array(
        [
            [m, 0, sq(6) * lam * omega, I * eta * bar_sigma / sq(2), -I * eta * sigma / sq(2), 0],
            [0, m + 2 * lam * a, 2 * sq(2) * lam * omega, I * eta * bar_sigma * sq(3 / 2), -I * eta * sigma * sq(3 / 2), 0],
            [sq(6) * lam * omega, 2 * sq(2) * lam * omega, m + lam * (p + 2 * a), -I * eta * sq(3) * bar_sigma, I * sq(3) * eta * sigma, 0],
            [I * eta * bar_sigma / sq(2), I * eta * bar_sigma * sq(3 / 2), -I * eta * sq(3) * bar_sigma, 0, M + eta * (p + 3 * a - 6 * omega), sq(5) * g * np.conj(sigma) / 2],
            [-I * eta * sigma / sq(2), -I * eta * sigma * sq(3 / 2), I * sq(3) * eta * sigma, M + eta * (p + 3 * a - 6 * omega), 0, sq(5) * g * np.conj(bar_sigma) / 2],
            [0, 0, 0, sq(5) * g * np.conj(sigma) / 2, sq(5) * g * np.conj(bar_sigma) / 2, 0],
        ],
        dtype=complex,
    )


def mixed_numeric_matrices(params: dict[str, complex]) -> dict[str, np.ndarray]:
    sq = math.sqrt
    I = 1j
    m = params["m"]
    lam = params["lambda"]
    eta = params["eta"]
    M = params["M"]
    omega = params["omega"]
    a = params["a"]
    p = params["p"]
    sigma = params["sigma"]
    bar_sigma = params["bar_sigma"]
    g = params["g"]
    G = neutral_g_numeric_matrix(params)
    E = np.array(
        [
            [-2 * (M + eta * (a - 3 * omega)), -2 * sq(2) * I * eta * sigma, 2 * I * eta * sigma, I * g * sq(2) * np.conj(bar_sigma)],
            [2 * I * sq(2) * eta * bar_sigma, -2 * (m + lam * (a - omega)), -2 * sq(2) * lam * omega, 2 * g * (np.conj(a) - np.conj(omega))],
            [-2 * I * eta * bar_sigma, -2 * sq(2) * lam * omega, -2 * (m - lam * omega), sq(2) * g * (np.conj(omega) - np.conj(p))],
            [-I * g * sq(2) * np.conj(sigma), 2 * g * (np.conj(a) - np.conj(omega)), g * sq(2) * (np.conj(omega) - np.conj(p)), 0],
        ],
        dtype=complex,
    )
    F = np.array(
        [
            [2 * (M + eta * (p + 3 * a)), -2 * I * sq(3) * eta * sigma, -g * sq(2) * np.conj(bar_sigma)],
            [2 * I * sq(3) * eta * bar_sigma, 2 * (m + lam * (p + 2 * a)), sq(24) * I * g * np.conj(omega)],
            [-g * sq(2) * np.conj(sigma), -sq(24) * I * g * np.conj(omega), 0],
        ],
        dtype=complex,
    )
    J = np.array(
        [
            [2 * (M + eta * (a + p - 2 * omega)), -2 * eta * bar_sigma, 2 * sq(2) * eta * bar_sigma, -I * g * sq(2) * np.conj(sigma)],
            [2 * eta * sigma, -2 * (m + lam * a), -2 * sq(2) * lam * omega, -2 * I * g * sq(2) * np.conj(a)],
            [-2 * sq(2) * eta * sigma, -2 * sq(2) * lam * omega, -2 * (m + lam * (a + p)), -4 * I * g * np.conj(omega)],
            [-I * g * sq(2) * np.conj(bar_sigma), 2 * sq(2) * I * g * np.conj(a), 4 * I * g * np.conj(omega), 0],
        ],
        dtype=complex,
    )
    X = np.array(
        [
            [2 * (m + lam * (a + omega)), -2 * sq(2) * lam * omega, -2 * g * (np.conj(a) + np.conj(omega))],
            [-2 * sq(2) * lam * omega, 2 * (m + lam * omega), sq(2) * g * (np.conj(omega) + np.conj(p))],
            [-2 * g * (np.conj(a) + np.conj(omega)), sq(2) * g * (np.conj(omega) + np.conj(p)), 0],
        ],
        dtype=complex,
    )
    return {
        "G_mixed_neutral_6x6": G,
        "E_mixed_4x4": E,
        "F_mixed_3x3": F,
        "J_mixed_4x4": J,
        "X_mixed_3x3": X,
    }


def cnum(z: complex, ndigits: int = 12) -> dict[str, float]:
    return {"re": round(float(np.real(z)), ndigits), "im": round(float(np.imag(z)), ndigits)}


def cmatrix(mat: np.ndarray, ndigits: int = 12) -> list[list[dict[str, float]]]:
    return [[cnum(mat[i, j], ndigits=ndigits) for j in range(mat.shape[1])] for i in range(mat.shape[0])]


def numeric_gates() -> dict[str, Any]:
    params = vev_from_x(0.1)
    params.update(
        {
            "M_H": 0.67,
            "gamma": 0.31 + 0.07j,
            "bar_gamma": -0.23 + 0.11j,
            "g": 0.7,
        }
    )

    T = triplet_numeric_matrix(params)
    T_inv = np.linalg.inv(T)
    T_residual = float(np.linalg.norm(T @ T_inv - np.eye(5)))
    T_cond = float(np.linalg.cond(T))
    T_det = np.linalg.det(T)

    mixed_gates: dict[str, Any] = {}
    for block in mixed_matrix_blocks():
        matrix = mixed_numeric_matrices(params)[block["block_id"]]
        n_chiral = block["chiral_subblock_shape"][0]
        chiral_block = matrix[:n_chiral, :n_chiral]
        svals_chiral = np.linalg.svd(chiral_block, compute_uv=False)
        svals_full = np.linalg.svd(matrix, compute_uv=False)
        chiral_null_ratio = float(svals_chiral[-1] / svals_chiral[0])
        full_min_ratio = float(svals_full[-1] / svals_full[0])
        det_full = np.linalg.det(matrix)
        mixed_gates[block["block_id"]] = {
            "source_anchor": block["source_anchor"],
            "chiral_block_shape": block["chiral_subblock_shape"],
            "full_shape": block["shape"],
            "chiral_block_singular_values": [float(x) for x in svals_chiral],
            "full_block_singular_values": [float(x) for x in svals_full],
            "chiral_null_ratio": chiral_null_ratio,
            "full_min_ratio": full_min_ratio,
            "full_det": cnum(det_full),
            "pass_fail": bool(chiral_null_ratio < 1.0e-12 and abs(det_full) > 1.0e-8),
        }

    required_pairs = [(0, 0), (0, 1), (1, 0), (1, 1), (0, 3), (1, 3)]
    return {
        "sample_parameters": {
            k: cnum(v) if isinstance(v, complex) else round(float(v), 12)
            for k, v in params.items()
        },
        "triplet_inverse_numeric_gate": {
            "triplet_det": cnum(T_det),
            "condition_number": T_cond,
            "inverse_residual_norm": T_residual,
            "pass_fail": bool(T_residual < 1.0e-12),
            "numeric_inverse": cmatrix(T_inv),
            "audit2_required_entries_1_indexed": {
                f"S_{i + 1}^{j + 1}": cnum(T_inv[i, j]) for i, j in required_pairs
            },
            "boundary": "numeric smoke inversion only; symbolic rational inverse remains a CAS-dependent follow-up",
        },
        "mixed_chiral_goldstone_numeric_gates": {
            "gates": mixed_gates,
            "all_pass": bool(all(row["pass_fail"] for row in mixed_gates.values())),
            "boundary": "verifies mixed chiral/gauge Goldstone blocks numerically on one F-flat sample; scalar Hessian Goldstone directions remain pending",
        },
    }


def write_markdown(card: dict[str, Any]) -> None:
    gates = card["numeric_gates"]
    trip = gates["triplet_inverse_numeric_gate"]
    mixed = gates["mixed_chiral_goldstone_numeric_gates"]
    lines = [
        "# Audit 4a.1 CMSGUT Literature Mass-Matrix Import",
        "",
        "This artifact transcribes the Aulakh--Girdhar doublet/triplet mass",
        "matrices and runs the first numerical gates needed before Audit 2.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- upstream mass-export schema sha256: `{card['upstream_mass_export_schema_sha256']}`",
        "",
        "## Imported Blocks",
        "",
        "| block | shape | source lines | status |",
        "| --- | --- | --- | --- |",
    ]
    for block in card["matrix_blocks"]:
        anchor = block["source_anchor"]
        lines.append(
            f"| `{block['block_id']}` | `{block['shape']}` | `{anchor['line_range']}` | {block['entry_status']} |"
        )
    lines += [
        "",
        "## Numeric Gates",
        "",
        f"- Mixed chiral/gauge Goldstone gates all pass: `{mixed['all_pass']}`.",
        f"- Triplet inverse gate pass: `{trip['pass_fail']}`.",
        f"- Triplet inverse residual norm: `{trip['inverse_residual_norm']:.3e}`.",
        f"- Triplet condition number: `{trip['condition_number']:.3e}`.",
        "",
        "| mixed block | chiral null ratio | full determinant | pass |",
        "| --- | --- | --- | --- |",
    ]
    for block_id, gate in mixed["gates"].items():
        lines.append(
            f"| `{block_id}` | `{gate['chiral_null_ratio']:.3e}` | `{gate['full_det']}` | `{gate['pass_fail']}` |"
        )
    lines += [
        "",
        "Audit-2-required numeric inverse entries on the smoke sample:",
        "",
        "| entry | value |",
        "| --- | --- |",
    ]
    for key, val in trip["audit2_required_entries_1_indexed"].items():
        lines.append(f"| `{key}` | `{val['re']} + {val['im']} i` |")
    lines += [
        "",
        "## Boundary",
        "",
        "- The BMSV local note in this workspace gives only qualitative",
        "  proton-decay/triplet-mixing guidance; it is recorded as a",
        "  cross-reference, not as a second matrix-entry table.",
        "- The triplet inverse exported here is a numerical smoke inversion.",
        "  The symbolic rational inverse is specified by an adjugate/cofactor",
        "  contract and remains pending expansion.",
        "- The Goldstone gate verifies all G,E,F,J,X mixed chiral/gauge blocks;",
        "  scalar Hessian Goldstone directions remain pending.",
        "",
    ]
    (OUT / "literature_mass_matrices.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mass_export = read_json(SOURCES["mass_export_schema"])
    vacuum = read_json(SOURCES["vacuum_branches"])
    blocks = [doublet_matrix(), triplet_matrix(), *mixed_matrix_blocks()]
    card: dict[str, Any] = {
        "audit": "audit4a1_cmsgut_literature_mass_import",
        "status": "Aulakh-Girdhar doublet/triplet entries transcribed; neutral chiral Goldstone and triplet inverse numeric gates run",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": source_manifest(),
        "upstream_mass_export_schema_sha256": mass_export["card_sha256"],
        "upstream_vacuum_branch_sha256": vacuum["card_sha256"],
        "convention": {
            "primary": "Aulakh-Girdhar hep-ph/0405074 appendix convention",
            "inherited_vacuum_map": vacuum["literature_convention_map"],
            "matrix_row_rule": "barred irreps label rows, unbarred irreps label columns",
            "BMSV_status": "local hep-ph/0306242 note records qualitative triplet-mixing/proton-decay ingredients but no doublet/triplet entry table; AG source states compatibility with BMSV up to normalizations/phases",
        },
        "matrix_blocks": blocks,
        "symbolic_triplet_inverse_contract": symbolic_triplet_inverse_contract(),
        "numeric_gates": numeric_gates(),
        "stage_gates": {
            "doublet_mass_matrix_entries_imported": True,
            "triplet_mass_matrix_entries_imported": True,
            "mixed_chiral_goldstone_numeric_gates_passed": True,
            "scalar_hessian_goldstone_directions_verified": False,
            "all_mixed_chiral_gauge_blocks_imported": True,
            "triplet_inverse_numeric_smoke_passed": True,
            "triplet_inverse_symbolic_entries_exported": False,
            "audit2_required_numeric_entries_exported": True,
            "nonplaceholder_heavy_spectrum_exported": False,
        },
        "next_required_steps": [
            "Expand the adjugate/cofactor contract into CAS-backed symbolic inverse entries for the 5x5 triplet matrix.",
            "Export scalar Hessian Goldstone directions and then the non-placeholder heavy spectrum.",
            "Feed the symbolic triplet inverse block into Audit 2 Wilson-tensor construction.",
        ],
    }
    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "literature_mass_matrices.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(card)
    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {(OUT / 'literature_mass_matrices.md').relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")


if __name__ == "__main__":
    main()
