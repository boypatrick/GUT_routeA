#!/usr/bin/env python3
"""Audit the missing Spin(10) Yukawa intertwiner normalization from matrices.

Uses the actual P1 normalized self-dual five-form basis.  No fitted Yukawa
matrix or new field is introduced.  The raw Clifford contraction includes
the action-card 1/5!: summing over ordered quintuplets cancels that factorial.
"""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "code" / "verify_p54_p1_hessian_spectrum.py"
OUTPUT = ROOT / "output" / "p54_spinor_intertwiners"


def load_p1():
    spec = importlib.util.spec_from_file_location("p54_intertwiner_p1", SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def kron_all(matrices):
    result = np.ones((1, 1), complex)
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def gamma_matrices():
    identity = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    return np.array([
        kron_all([z] * k + [pauli] + [identity] * (4 - k))
        for k in range(5) for pauli in (x, y)
    ])


def product(matrices):
    result = np.eye(matrices[0].shape[0], dtype=complex)
    for matrix in matrices:
        result = result @ matrix
    return result


def lift(generator, gammas):
    return sum(
        generator[i, j] * gammas[i] @ gammas[j] / 2
        for i in range(10) for j in range(i + 1, 10)
    )


def wedge_matrix(generator, quints):
    index = {q: i for i, q in enumerate(quints)}
    result = np.zeros((252, 252))
    for row, indices in enumerate(quints):
        for slot in range(5):
            for source in np.flatnonzero(generator[indices[slot]]):
                replaced = list(indices)
                replaced[slot] = int(source)
                if len(set(replaced)) != 5:
                    continue
                inversions = sum(replaced[a] > replaced[b]
                                 for a in range(5) for b in range(a + 1, 5))
                result[row, index[tuple(sorted(replaced))]] += (
                    generator[indices[slot], source] * (-1) ** inversions
                )
    return result


def maxabs(matrix):
    return float(np.max(np.abs(matrix)))


def unique_state(operators, targets):
    penalty = sum((op - target * np.eye(op.shape[0])) @
                  (op - target * np.eye(op.shape[0]))
                  for op, target in zip(operators, targets))
    values, vectors = np.linalg.eigh((penalty + penalty.conj().T) / 2)
    if np.sum(np.abs(values) < 1e-9) != 1:
        raise ValueError(f"expected unique state, eigenvalues={values[:5]}")
    return vectors[:, 0]


def run():
    p1 = load_p1()
    gammas = gamma_matrices()
    identity = np.eye(32, dtype=complex)
    chirality = (-1j) ** 5 * product(list(gammas))
    charge_c = product([gammas[i] for i in (1, 3, 5, 7, 9)])
    b10_full = np.array([charge_c @ gamma for gamma in gammas])
    b5_full = np.array([charge_c @ product([gammas[i] for i in quint])
                        for quint in p1.QUINTS])
    b3_full = np.array([charge_c @ product([gammas[i] for i in triple])
                        for triple in itertools.combinations(range(10), 3)])
    u = p1.U126
    branches = []
    for chirality_sign in (1, -1):
        indices = np.flatnonzero(np.diag(chirality).real == chirality_sign)
        b10 = b10_full[:, indices][:, :, indices]
        b5 = b5_full[:, indices][:, :, indices]
        b126 = np.einsum("qa,qij->aij", u, b5)
        majorana = np.einsum("q,qij->ij", p1.GEOM["omega"], b5)
        singular = np.linalg.svd(majorana, compute_uv=False)
        branches.append({
            "chirality": chirality_sign, "indices": indices,
            "b10": b10, "b5": b5, "b126": b126,
            "majorana": majorana, "singular": singular,
        })
    selected = max(branches, key=lambda row: row["singular"][0])
    indices = selected["indices"]
    b10, b5, b126 = [selected[key] for key in ("b10", "b5", "b126")]
    mr_normalization = float(selected["singular"][0])
    vector_covariance = 0.0
    five_covariance = 0.0
    spin_trace = []
    # Verify all 45 Lie-algebra directions; use raw wedge action instead of
    # rebuilding the scalar Hessian or choosing a favorable subgroup.
    for i in range(10):
        for j in range(i + 1, 10):
            vector = np.zeros((10, 10))
            vector[i, j], vector[j, i] = 1.0, -1.0
            spin = lift(vector, gammas)[np.ix_(indices, indices)]
            spin_trace.append(float(-np.trace(spin @ spin).real))
            residual10 = (np.einsum("ki,akj->aij", spin, b10)
                          + np.einsum("aik,kj->aij", b10, spin)
                          + np.einsum("bij,ba->aij", b10, vector))
            wedge = wedge_matrix(vector, p1.QUINTS)
            rep126 = u.conj().T @ wedge @ u
            residual126 = (np.einsum("ki,akj->aij", spin, b126)
                           + np.einsum("aik,kj->aij", b126, spin)
                           + np.einsum("bij,ba->aij", b126, rep126))
            vector_covariance = max(vector_covariance, maxabs(residual10))
            five_covariance = max(five_covariance, maxabs(residual126))

    sm = p1.sm_generators()
    rep_spin = {key: [lift(g, gammas)[np.ix_(indices, indices)] for g in value]
                for key, value in sm.items()}
    rep126_sm = {key: [u.conj().T @ wedge_matrix(g, p1.QUINTS) @ u for g in value]
                 for key, value in sm.items()}
    casimir = lambda generators: -sum(g @ g for g in generators)
    spin_c3 = casimir(rep_spin["SU3"])
    spin_c2 = casimir(rep_spin["SU2"])
    spin_y = 1j * rep_spin["Y"][0]
    spin_t3 = 1j * rep_spin["SU2"][-1]
    singlet = unique_state([spin_c3, spin_c2, spin_y], [0, 0, 0])
    mr_entry = complex(singlet.T @ selected["majorana"] @ singlet)

    # Determine whether the selected chirality is the usual matter 16 or its
    # conjugate in P1's SM-generator orientation.  This is output, not assumed.
    y_eigenvalues = np.linalg.eigvalsh(spin_y)
    matter_sign = 1 if np.max(y_eigenvalues) > 0.9 else -1
    fermion_operators = [spin_c3, spin_c2, spin_y, spin_t3,
                         1j * rep_spin["SU3"][-1]]
    # One color weight fixes a unique down-quark state.  Search it, then use
    # the Yukawa contraction to select the conjugate-color partner.
    def states(c3, c2, hypercharge, t3):
        penalty = sum((op - target * np.eye(16)) @ (op - target * np.eye(16))
                      for op, target in zip(fermion_operators[:4],
                                             [c3, c2, hypercharge, t3]))
        values, vectors = np.linalg.eigh((penalty + penalty.conj().T) / 2)
        return vectors[:, np.abs(values) < 1e-9]
    down_l = states(4/3, 3/4, matter_sign/6, -matter_sign/2)
    down_r = states(4/3, 0, matter_sign/3, 0)
    electron_l = states(0, 3/4, -matter_sign/2, -matter_sign/2)
    electron_r = states(0, 0, matter_sign, 0)
    scalar_y = -matter_sign/2
    scalar_t3 = matter_sign/2
    v10 = unique_state([casimir(sm["SU3"]), casimir(sm["SU2"]),
                       1j * sm["Y"][0], 1j * sm["SU2"][-1]],
                      [0, 3/4, scalar_y, scalar_t3])
    v126 = unique_state([casimir(rep126_sm["SU3"]), casimir(rep126_sm["SU2"]),
                        1j * rep126_sm["Y"][0], 1j * rep126_sm["SU2"][-1]],
                       [0, 3/4, scalar_y, scalar_t3])
    m10 = np.einsum("a,aij->ij", v10, b10)
    m126 = np.einsum("a,aij->ij", v126, b126)
    d10 = down_l.T @ m10 @ down_r
    d126 = down_l.T @ m126 @ down_r
    e10 = complex((electron_l.T @ m10 @ electron_r)[0, 0])
    e126 = complex((electron_l.T @ m126 @ electron_r)[0, 0])
    location = np.unravel_index(np.argmax(abs(d10)), d10.shape)
    relative_clebsch = e126 * d10[location] / (e10 * d126[location])
    raw_down10 = float(np.linalg.svd(d10, compute_uv=False)[0])
    raw_down126 = float(np.linalg.svd(d126, compute_uv=False)[0])
    normalized_down126 = raw_down126 / mr_normalization
    opposite = min(branches, key=lambda row: row["singular"][0])
    # The P1 self-duality selects charge-conjugate SM matter in its existing
    # Y orientation.  Conjugating the scalar irrep and chirality together is
    # the corresponding standard-matter dictionary, not an isolated label swap.
    conjugate_majorana = np.einsum("q,qij->ij", np.conj(p1.GEOM["omega"]), opposite["b5"])
    conjugate_mr_singular = np.linalg.svd(conjugate_majorana, compute_uv=False)
    neutrino_l = states(0, 3/4, -matter_sign/2, matter_sign/2)
    triplet126 = unique_state(
        [casimir(rep126_sm["SU3"]), casimir(rep126_sm["SU2"]),
         1j * rep126_sm["Y"][0], 1j * rep126_sm["SU2"][-1]],
        [0, 2, matter_sign, -matter_sign],
    )
    triplet_matrix = np.einsum("a,aij->ij", triplet126, b126)
    ll_raw = complex((neutrino_l.T @ triplet_matrix @ neutrino_l)[0, 0])
    triplet_source_path = ROOT / "output" / "p54_typeii_triplet_source.json"
    triplet_export = None
    if triplet_source_path.exists():
        source = json.loads(triplet_source_path.read_text())
        record = source["complex_triplet_basis_328x2"]
        embedded = np.asarray(record["real"]) + 1j * np.asarray(record["imag"])
        vector = embedded[:, 1]
        hol = (vector[p1.SL_SIGMA_RE] + 1j * vector[p1.SL_SIGMA_IM]) / math.sqrt(2)
        anti = (vector[p1.SL_SIGMA_RE] - 1j * vector[p1.SL_SIGMA_IM]) / math.sqrt(2)
        conjugated_hol = (np.conj(vector[p1.SL_SIGMA_RE]) +
                         1j * np.conj(vector[p1.SL_SIGMA_IM])) / math.sqrt(2)
        exported_ll = complex((neutrino_l.T @ np.einsum(
            "a,aij->ij", conjugated_hol, b126) @ neutrino_l)[0, 0])
        triplet_export = {
            "source": str(triplet_source_path.relative_to(ROOT.parent)),
            "exported_Y": 1,
            "Sigma_hol_norm": float(np.linalg.norm(hol)),
            "Sigma_anti_norm": float(np.linalg.norm(anti)),
            "conjugated_export_overlap_with_canonical_LL_triplet": float(abs(np.vdot(triplet126, conjugated_hol))),
            "raw_LL_coefficient_for_conjugated_export": {"re": float(exported_ll.real), "im": float(exported_ll.imag)},
            "phase_scope": "Matrix coefficient uses deterministic eigensolver fermion phases; phase alignment to fitted h/f and the exported Higgs direction remains a joint matching task.",
        }

    checks = []
    def check(name, residual, tolerance=1e-11):
        checks.append({"name": name, "residual": float(residual),
                       "tolerance": tolerance, "pass": bool(residual < tolerance)})
    check("Clifford anticommutators", max(maxabs(gammas[i] @ gammas[j] +
          gammas[j] @ gammas[i] - 2 * (i == j) * identity)
          for i in range(10) for j in range(10)))
    check("chirality squares to identity", maxabs(chirality @ chirality - identity))
    check("charge conjugation antisymmetry", maxabs(charge_c.T + charge_c))
    check("charge conjugation covariance", max(maxabs(charge_c @ g + g.T @ charge_c)
                                               for g in gammas))
    check("10 Yukawa tensors symmetric", maxabs(b10_full - b10_full.transpose(0, 2, 1)))
    check("120 Yukawa tensors antisymmetric", maxabs(b3_full + b3_full.transpose(0, 2, 1)))
    check("five-form Yukawa tensors symmetric", maxabs(b5_full - b5_full.transpose(0, 2, 1)))
    check("all 45 vector covariance equations", vector_covariance)
    check("all 45 self-dual 126 covariance equations", five_covariance)
    check("opposite chirality annihilates P1 self-duality", maxabs(opposite["b126"]))
    check("P1 omega produces rank-one Majorana mass", selected["singular"][1])
    check("P1 omega acts only on the SM-singlet fermion", abs(abs(mr_entry) - mr_normalization))
    check("raw Majorana coefficient is 4 sqrt(2)", abs(mr_normalization - 4 * math.sqrt(2)))
    check("raw unit vector doublet coefficient is sqrt(2)", abs(raw_down10 - math.sqrt(2)))
    check("raw unit 126 doublet coefficient is 2/sqrt(3)", abs(raw_down126 - 2 / math.sqrt(3)))
    check("canonical singlet/Dirac 126 ratio is 2 sqrt(6)", abs(mr_normalization / raw_down126 - 2 * math.sqrt(6)))
    check("lepton/down relative Clebsch is minus three", abs(relative_clebsch + 3))
    check("conjugating chirality and scalar irrep preserves Majorana norm", maxabs(conjugate_mr_singular - selected["singular"]))
    check("canonical LL triplet and RR singlet magnitudes agree", abs(abs(ll_raw) - mr_normalization))
    if triplet_export is not None:
        check("conjugated exported triplet is the canonical LL direction", abs(triplet_export["conjugated_export_overlap_with_canonical_LL_triplet"] - 1))
    report = {
        "schema": "route-f-p54-spinor-intertwiner-v1", "date": "2026-09-05",
        "status": "explicit Clifford and Yukawa intertwiners verified; normalization translation required",
        "conventions": {
            "Clifford": "{Gamma_i,Gamma_j}=2 delta_ij; Hermitian 32x32 Jordan-Wigner matrices",
            "chirality": "Gamma_star=(-i)^5 Gamma_1...Gamma_10",
            "charge_conjugation": "C=Gamma_2 Gamma_4 Gamma_6 Gamma_8 Gamma_10",
            "tensor_contraction": "sum_{i1<...<i5} Sigma_i1...i5 C Gamma_i1...i5 = (1/5!) full antisymmetric contraction",
            "P1_scalar_metric": "U126^dagger U126=I; ||omega||=1",
            "spin_generator": "R_ij=Gamma_i Gamma_j/2 for vector R_ij[i,j]=+1; physical Hermitian generator i R",
        },
        "selected_chirality": selected["chirality"],
        "SM_matter_orientation": matter_sign,
        "orientation_warning": "For the literal P1 U126 and C Gamma[5] contraction, the selected chiral matter has conjugate SM hypercharges. Standard matter charges require conjugating scalar irrep and spinor chirality together (or an equivalent global embedding conjugation). Earlier geometric overlap keys cannot be assigned physical u/d labels without this dictionary.",
        "hypercharge_eigenvalues": y_eigenvalues.tolist(),
        "spin_trace_raw_plane_generator": spin_trace[0],
        "generator_normalization_note": "raw vector plane generator has T10=2 and T16=4; divide generator by sqrt(2) for action-card tr16(TaTb)=2 delta_ab",
        "majorana_singular_values_sigma1_raw": selected["singular"].tolist(),
        "absolute_CG": {
            "raw_majorana_per_sigma": mr_normalization,
            "factor_for_MR_equal_sigma_f": 1 / mr_normalization,
            "raw_vector_down_doublet": raw_down10,
            "factor_for_unit_vector_Dirac_h": 1 / raw_down10,
            "raw_126_down_doublet": raw_down126,
            "126_down_doublet_if_MR_equal_sigma_f": normalized_down126,
            "126_lepton_doublet_if_MR_equal_sigma_f": abs(e126) / mr_normalization,
            "canonical_MR_over_Dirac126_ratio": mr_normalization / raw_down126,
            "raw_LL_triplet_coefficient_magnitude": abs(ll_raw),
            "LL_triplet_if_MR_equal_sigma_f": abs(ll_raw) / mr_normalization,
            "signed_relative_lepton_down_CG": {"re": float(relative_clebsch.real), "im": float(relative_clebsch.imag)},
        },
        "actual_triplet_export_alignment": triplet_export,
        "remaining_gates": [
            "Translate the action-card h,f and four VEV/overlap symbols to these canonically normalized matrix intertwiners before fitting.",
            "Build all actual heavy/light fermion matrices on the P1/P2 scalar background with one fixed normalization and re-evaluate the CW projection.",
            "The minus-three Clebsch is verified; defining both unit doublet coupling f and MR=sigma f requires a nontrivial conversion of the canonical singlet sigma or of the doublet VEV.",
            "The absolute LL-triplet coefficient is verified; carry its complex phase with the same fermion/scalar dictionary into the actual type-II source response.",
            "No global flavor/seesaw fit or two-loop PS Yukawa transport is supplied by this algebraic audit.",
        ],
        "checks": checks,
        "summary": {"passed": sum(row["pass"] for row in checks), "total": len(checks),
                    "all_pass": all(row["pass"] for row in checks)},
        "sources": [{"path": str(path.relative_to(ROOT.parent)),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for path in (Path(__file__), SOURCE, *([triplet_source_path] if triplet_export is not None else []))],
    }
    return report


def main():
    report = run()
    OUTPUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    cg = report["absolute_CG"]
    lines = ["# P54 Spin(10) Yukawa intertwiner audit", "", report["status"], "",
             f"Checks: {report['summary']['passed']}/{report['summary']['total']}.", "",
             f"Selected chirality: {report['selected_chirality']}; SM matter orientation: {report['SM_matter_orientation']}.", "",
             report["orientation_warning"], "",
             "| Matrix contraction | Absolute coefficient |", "|---|---:|",
             *[f"| {key} | {value:.12g} |" for key, value in cg.items() if isinstance(value, float)],
             "", f"Signed lepton/down relative Clebsch: {cg['signed_relative_lepton_down_CG']}.", "",
             "The raw 1/5! action-card contraction on the actual unit P1 omega is not normalized to MR=sigma f.",
             "The output gives the conversion explicitly. It does not insert a phenomenological convention silently.", "",
             "Remaining gates:", "", *[f"- {item}" for item in report["remaining_gates"]], ""]
    OUTPUT.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(report["summary"]))
    print(json.dumps(cg, indent=2))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
