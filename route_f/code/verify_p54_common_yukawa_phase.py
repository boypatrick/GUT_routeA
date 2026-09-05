#!/usr/bin/env python3
"""Match the actual P54 scalar bases to one complex Spin(10) fermion basis.

The global conjugation is explicit: physical scalar coordinates are K x_old,
the physical self-dual basis is conjugate(U126), and the physical spinor has
standard SM hypercharges.  All Yukawa coefficients below use the same states,
including MR and the induced LL triplet response.  This is algebraic matching,
not a flavor fit or a finite-threshold calculation.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

RF = Path(__file__).resolve().parents[1]
INTERTWINER = RF / "code/verify_p54_spinor_intertwiners.py"
BOSON = RF / "output/p54_full_doublet_cw.json"
TRIPLET = RF / "output/p54_typeii_triplet_source.json"
OUTPUT = RF / "output/p54_common_yukawa_phase"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(obj)
    return obj


def cjson(value):
    array = np.asarray(value)
    return {"real": array.real.tolist(), "imag": array.imag.tolist()}


def from_cjson(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def build_geometry():
    """Return raw, unnormalized intertwiners for standard-hypercharge 16.

    Stable consumer API: b10[a,i,j], b126[A,i,j], U126_physical[q,A],
    spin_sm/group lists, p1, lift, wedge_matrix, conjugation_diagonal.
    Scalar coefficients z126 multiply U126_physical with no conjugation.
    The action is h Psi Psi phi_physical* + f Psi Psi Sigma_physical.
    """
    audit = module("p54_common_intertwiner", INTERTWINER)
    p1 = audit.load_p1()
    gammas = audit.gamma_matrices()
    chirality = (-1j) ** 5 * audit.product(list(gammas))
    indices = np.flatnonzero(np.diag(chirality).real > 0)
    charge_c = audit.product([gammas[i] for i in (1, 3, 5, 7, 9)])
    b10 = np.array([(charge_c @ gamma)[np.ix_(indices, indices)] for gamma in gammas])
    b5 = np.array([(charge_c @ audit.product([gammas[i] for i in q]))[
        np.ix_(indices, indices)] for q in p1.QUINTS])
    u = np.conj(p1.U126)
    b126 = np.einsum("qa,qij->aij", u, b5)
    sm = p1.sm_generators()
    spin_sm = {key: [audit.lift(g, gammas)[np.ix_(indices, indices)] for g in generators]
               for key, generators in sm.items()}
    scalar126_sm = {key: [u.conj().T @ audit.wedge_matrix(g, p1.QUINTS) @ u
                         for g in generators] for key, generators in sm.items()}
    k = np.ones(p1.N_REAL)
    k[p1.SL_SIGMA_IM] = -1
    k[p1.SL_H_IM] = -1
    k[327] = -1
    return {"p1": p1, "gammas": gammas, "C": charge_c, "indices": indices,
            "b10": b10, "b5": b5, "b126": b126, "U126_physical": u,
            "spin_sm": spin_sm, "scalar126_sm": scalar126_sm, "scalar10_sm": sm,
            "conjugation_diagonal": k, "lift": audit.lift,
            "wedge_matrix": audit.wedge_matrix}


def real_scalar_yukawa_tensors(geometry=None, *, old_coordinates=False):
    """Return (Zh,Zf), each shape (328,16,16), in the raw action convention.

    M_fermion(x)=h_raw tensor sum_A Zh[A] x[A] +
                 f_raw tensor sum_A Zf[A] x[A], with family tensor products
    supplied by the caller. old_coordinates=True pulls the tensors back by K
    for direct use on cached historical scalar Hessian coordinates.
    """
    g = build_geometry() if geometry is None else geometry
    p1 = g["p1"]
    zh = np.zeros((p1.N_REAL, 16, 16), complex)
    zf = np.zeros_like(zh)
    zh[p1.SL_H_RE] = g["b10"] / math.sqrt(2)
    zh[p1.SL_H_IM] = -1j * g["b10"] / math.sqrt(2)
    zf[p1.SL_SIGMA_RE] = g["b126"] / math.sqrt(2)
    zf[p1.SL_SIGMA_IM] = 1j * g["b126"] / math.sqrt(2)
    if old_coordinates:
        zh *= g["conjugation_diagonal"][:, None, None]
        zf *= g["conjugation_diagonal"][:, None, None]
    return zh, zf


def phase_positive(vector):
    return vector * np.exp(-1j * np.angle(vector[np.argmax(abs(vector))]))


def run():
    g = build_geometry()
    p1, sm = g["p1"], g["spin_sm"]
    boson = json.loads(BOSON.read_text())
    triplet = json.loads(TRIPLET.read_text())
    record = boson["basis"]
    candidates = [value for value in record.values() if isinstance(value, dict)
                  and "real" in value and np.asarray(value["real"]).shape == (328, 4)]
    if len(candidates) != 1:
        raise ValueError("Expected exactly one exported scalar doublet embedding")
    old_basis = from_cjson(candidates[0])
    k = g["conjugation_diagonal"]
    basis = k[:, None] * old_basis
    light = from_cjson(boson["retuned_bosonic_eigenpair"]["light_coefficients"])
    tree_light = from_cjson(record["tree_light_coefficients"])
    old_triplet = from_cjson(triplet["complex_triplet_basis_328x2"])
    triplet_basis = k[:, None] * old_triplet
    induced = from_cjson(triplet["complex_induced_amplitude_over_v2_times_omega"])

    def scalar_matrices(direction):
        phi_star = (direction[p1.SL_H_RE] - 1j * direction[p1.SL_H_IM]) / math.sqrt(2)
        sigma = (direction[p1.SL_SIGMA_RE] + 1j * direction[p1.SL_SIGMA_IM]) / math.sqrt(2)
        return (np.einsum("a,aij->ij", phi_star, g["b10"]),
                np.einsum("a,aij->ij", sigma, g["b126"]))
    up_h, up_f, down_h, down_f = [], [], [], []
    for column in basis.T:
        mh, mf = scalar_matrices(column)
        up_h.append(mh)
        up_f.append(mf)
        mh, mf = scalar_matrices(np.conj(column))
        down_h.append(mh)
        down_f.append(mf)
    up_h, up_f, down_h, down_f = map(np.asarray, (up_h, up_f, down_h, down_f))

    casimir = lambda matrices: -sum(matrix @ matrix for matrix in matrices)
    c3, c2 = casimir(sm["SU3"]), casimir(sm["SU2"])
    y, t3 = 1j * sm["Y"][0], 1j * sm["SU2"][-1]
    color3, color8 = 1j * sm["SU3"][-2], 1j * sm["SU3"][-1]
    # These commuting Cartans are diagonal in the Jordan-Wigner occupation basis.
    diagonality = max(np.linalg.norm(op - np.diag(np.diag(op)))
                      for op in (c3, c2, y, t3, color3, color8))
    if diagonality > 1e-10:
        raise ValueError("Fock-basis state selector requires diagonal commuting operators")
    eye = np.eye(16, dtype=complex)
    def states(targets):
        mask = np.ones(16, dtype=bool)
        for operator, target in zip((c3, c2, y, t3), targets):
            mask &= abs(np.diag(operator).real - target) < 1e-10
        indices = np.flatnonzero(mask)
        # Sort quarks by fixed color Cartans, not eigensolver ordering.
        indices = sorted(indices, key=lambda index: (color3[index, index].real,
                                                      color8[index, index].real))
        return eye[:, indices]
    states_by_name = {
        "uL": states((4/3, 3/4, 1/6, 1/2)),
        "dL": states((4/3, 3/4, 1/6, -1/2)),
        "uc": states((4/3, 0, -2/3, 0)),
        "dc": states((4/3, 0, 1/3, 0)),
        "nuL": states((0, 3/4, -1/2, 1/2)),
        "eL": states((0, 3/4, -1/2, -1/2)),
        "nuc": states((0, 0, 0, 0)),
        "ec": states((0, 0, 1, 0)),
    }
    # Fix otherwise unphysical component phases once, by making each 10_H
    # reference Yukawa coefficient positive. This also fixes the relative
    # phases carried into f, MR and ML; no later entry-by-entry sign choices.
    pairing = {}
    for name, left, right, matrix in (
        ("u", "uL", "uc", up_h[0]), ("d", "dL", "dc", down_h[1]),
        ("nu", "nuL", "nuc", up_h[0]), ("e", "eL", "ec", down_h[1]),
    ):
        vl, vr = states_by_name[left], states_by_name[right]
        preliminary = vl.T @ matrix @ vr
        permutation = np.argmax(abs(preliminary), axis=1)
        if len(set(permutation.tolist())) != len(permutation):
            raise ValueError("Color partners are not a bijection")
        vr = vr[:, permutation]
        phases = np.exp(-1j * np.angle(np.diag(vl.T @ matrix @ vr)))
        vr = vr * phases[None, :]
        states_by_name[right] = vr
        pairing[name] = {"left": left, "right": right,
                         "right_permutation": permutation.tolist(),
                         "right_rephasing": cjson(phases)}

    coefficients = {}
    color_error = 0.0
    for name, hmatrices, fmatrices in (
        ("u", up_h, up_f), ("d", down_h, down_f),
        ("e", down_h, down_f), ("nu", up_h, up_f),
    ):
        pair = pairing[name]
        vl, vr = states_by_name[pair["left"]], states_by_name[pair["right"]]
        values = {}
        for label, matrices in (("h_raw", hmatrices), ("f_raw", fmatrices)):
            blocks = np.asarray([vl.T @ matrix @ vr for matrix in matrices])
            scalars = np.trace(blocks, axis1=1, axis2=2) / vl.shape[1]
            color_error = max(color_error, float(np.linalg.norm(
                blocks - scalars[:, None, None] * np.eye(vl.shape[1]))))
            values[label] = scalars
        coefficients[name] = values

    # MR uses the actual canonical P1 unit vacuum, globally conjugated.
    mr = np.einsum("a,aij->ij", np.conj(p1.GEOM["omega126"]), g["b126"])
    n = states_by_name["nuc"]
    l = states_by_name["nuL"]
    kappa_r = complex((n.T @ mr @ n)[0, 0])
    triplet_ll = []
    triplet_ll_conjugate = []
    for vector in triplet_basis.T:
        _, matrix = scalar_matrices(vector)
        _, conjugate_matrix = scalar_matrices(np.conj(vector))
        triplet_ll.append(complex((l.T @ matrix @ l)[0, 0]))
        triplet_ll_conjugate.append(complex((l.T @ conjugate_matrix @ l)[0, 0]))
    triplet_ll = np.asarray(triplet_ll)
    triplet_ll_conjugate = np.asarray(triplet_ll_conjugate)
    kappa_h = complex(coefficients["d"]["h_raw"][1])
    kappa_d = complex(coefficients["d"]["f_raw"][2])
    projected = {}
    for label, vector in (("tree", tree_light), ("bosonic", light)):
        values = {}
        for species, rows in coefficients.items():
            weights = vector if species in ("u", "nu") else np.conj(vector)
            values[species] = {coupling: complex(row @ weights)
                               for coupling, row in rows.items()}
        a = values["d"]["h_raw"] / kappa_h
        b = values["u"]["h_raw"] / kappa_h
        d = values["d"]["f_raw"] / kappa_d
        e = values["u"]["f_raw"] / kappa_d
        projected[label] = {
            "raw_species_coefficients": {species: {key: cjson(value) for key, value in rows.items()}
                                         for species, rows in values.items()},
            "Dirac_unit_abde": {key: cjson(value) for key, value in (("a", a), ("b", b), ("d", d), ("e", e))},
            "complex_r": cjson(b / a), "complex_s": cjson(a * e / (b * d)),
            "sum_rule_residual": float(max(abs(values["nu"]["h_raw"] - values["u"]["h_raw"]),
                                           abs(values["e"]["h_raw"] - values["d"]["h_raw"]),
                                           abs(values["nu"]["f_raw"] + 3 * values["u"]["f_raw"]),
                                           abs(values["e"]["f_raw"] + 3 * values["d"]["f_raw"]))),
        }
    typeii_raw = complex(triplet_ll @ induced + triplet_ll_conjugate @ np.conj(induced))
    sigma = float(triplet["vacuum"]["sigma"])
    source_real = np.asarray(triplet["canonical_induced_z_over_v2_times_omega"])
    real_basis = np.column_stack([math.sqrt(2) * triplet_basis.real,
                                  -math.sqrt(2) * triplet_basis.imag])
    actual_response = real_basis @ source_real
    _, actual_ml = scalar_matrices(actual_response)
    direct_ml = complex((l.T @ actual_ml @ l)[0, 0])

    covariance10, covariance126 = 0.0, 0.0
    for i in range(10):
        for j in range(i+1, 10):
            scalar10 = np.zeros((10, 10))
            scalar10[i, j], scalar10[j, i] = 1.0, -1.0
            spin = g["lift"](scalar10, g["gammas"])[np.ix_(g["indices"], g["indices"])]
            u = g["U126_physical"]
            scalar126 = u.conj().T @ g["wedge_matrix"](scalar10, p1.QUINTS) @ u
            for tensors, scalar, label in ((g["b10"], scalar10, "10"),
                                           (g["b126"], scalar126, "126")):
                residual = (np.einsum("ki,akj->aij", spin, tensors)
                            + np.einsum("aik,kj->aij", tensors, spin)
                            + np.einsum("bij,ba->aij", tensors, scalar))
                error = float(np.max(abs(residual)))
                if label == "10":
                    covariance10 = max(covariance10, error)
                else:
                    covariance126 = max(covariance126, error)
    physical_ry = k[:, None] * p1.representation_matrix(p1.sm_generators()["Y"][0]) * k[None, :]
    physical_r3 = k[:, None] * p1.representation_matrix(p1.sm_generators()["SU2"][-1]) * k[None, :]
    standard_y = sorted([1/6]*6 + [-2/3]*3 + [1/3]*3 + [-1/2]*2 + [1, 0])
    test_phases = np.exp(1j * np.array([.17, -.29, .43, -.61]))
    rephased_light = light / test_phases
    phase_error = 0.0
    for species, rows in coefficients.items():
        if species in ("u", "nu"):
            left_factor, right_factor = test_phases, rephased_light
            original_weights = light
        else:
            left_factor, right_factor = np.conj(test_phases), np.conj(rephased_light)
            original_weights = np.conj(light)
        for row in rows.values():
            phase_error = max(phase_error, abs((row * left_factor) @ right_factor - row @ original_weights))
    triplet_phase = np.exp(1j * np.array([.37, -.51]))
    rephased_typeii = (triplet_ll * triplet_phase) @ (induced / triplet_phase)
    rephased_typeii += (triplet_ll_conjugate * np.conj(triplet_phase)) @ np.conj(induced / triplet_phase)
    zh, zf = real_scalar_yukawa_tensors(g, old_coordinates=True)
    old_response = k * actual_response
    tensor_ml = complex((l.T @ np.einsum("a,aij->ij", old_response, zf) @ l)[0, 0])
    spin_frame = np.column_stack(list(states_by_name.values()))
    # Non-diagonal copy changes must be tested on the actual embedded scalar
    # vectors: B -> B U, c -> U^dagger c. In the down sector the scalar enters
    # conjugated; using c instead of c* would fail this test by order one.
    rng = np.random.default_rng(540926)
    copy_unitary, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)))
    transformed_basis = basis @ copy_unitary
    transformed_light = copy_unitary.conj().T @ light
    copy_unitary_error = 0.0
    wrong_down_weight_error = 0.0
    for species, rows in coefficients.items():
        pair = pairing[species]
        vl, vr = states_by_name[pair["left"]], states_by_name[pair["right"]]
        is_up = species in ("u", "nu")
        for coupling_index, (coupling, row) in enumerate(rows.items()):
            transformed_row = []
            for column in transformed_basis.T:
                matrices = scalar_matrices(column if is_up else np.conj(column))
                block = vl.T @ matrices[coupling_index] @ vr
                transformed_row.append(np.trace(block) / vl.shape[1])
            transformed_row = np.asarray(transformed_row)
            old_value = row @ (light if is_up else np.conj(light))
            new_value = transformed_row @ (transformed_light if is_up else np.conj(transformed_light))
            copy_unitary_error = max(copy_unitary_error, abs(new_value-old_value))
            if not is_up:
                wrong_down_weight_error = max(wrong_down_weight_error,
                                              abs(transformed_row @ transformed_light-old_value))

    # Synthetic matrices here test tensor covariance only; they are not fit
    # points and cannot produce a flavor likelihood or a claimed spectrum.
    htest = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    htest = (htest+htest.T)/10
    ftest = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    ftest = np.eye(3) + (ftest+ftest.T)/10
    ynu_h = coefficients["nu"]["h_raw"] @ light
    ynu_f = coefficients["nu"]["f_raw"] @ light
    def seesaw(ht, ft):
        ynu = ynu_h*ht+ynu_f*ft
        return typeii_raw*ft - ynu @ np.linalg.solve(ft, ynu.T)/(2*sigma*kappa_r)
    mass_test = seesaw(htest, ftest)
    family_unitary, _ = np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
    family_phase = np.diag(np.exp(1j*np.array([.23, -.71, 1.19])))
    family_errors = []
    for rotation in (family_unitary, family_phase):
        moved = seesaw(rotation.T @ htest @ rotation, rotation.T @ ftest @ rotation)
        expected = rotation.T @ mass_test @ rotation
        family_errors.append(float(np.linalg.norm(moved-expected)/max(np.linalg.norm(expected), 1e-15)))
    alpha_l, alpha_r = .41, -.63
    ynu_test = ynu_h*htest+ynu_f*ftest
    ynu_rephased = np.exp(1j*(alpha_l+alpha_r))*ynu_test
    mr_coefficient_rephased = np.exp(2j*alpha_r)*kappa_r
    ml_coefficient_rephased = np.exp(2j*alpha_l)*typeii_raw
    mass_rephased = ml_coefficient_rephased*ftest
    mass_rephased -= ynu_rephased @ np.linalg.solve(ftest, ynu_rephased.T)/(2*sigma*mr_coefficient_rephased)
    neutrino_rephase_error = float(np.linalg.norm(mass_rephased-np.exp(2j*alpha_l)*mass_test)
                                  / max(np.linalg.norm(mass_test), 1e-15))
    checks = []
    def check(name, error, tolerance=1e-10):
        checks.append({"name": name, "residual": float(error), "tolerance": tolerance,
                       "pass": bool(error < tolerance)})
    check("standard matter 16 hypercharge spectrum", np.max(abs(np.linalg.eigvalsh(y) - standard_y)))
    check("all spinor component wavefunctions have unit norm", max(np.linalg.norm(v.conj().T @ v - np.eye(v.shape[1])) for v in states_by_name.values()))
    check("all 16 fixed spinor states form one orthonormal basis", np.linalg.norm(spin_frame.conj().T @ spin_frame - np.eye(16)))
    check("commuting Fock-basis quantum numbers are diagonal", diagonality)
    check("global scalar conjugation is an orthogonal involution", np.max(abs(k*k - 1)))
    check("four complex scalar-copy wavefunctions are orthonormal", np.linalg.norm(basis.conj().T @ basis - np.eye(4)))
    check("global conjugation preserves scalar Y=+1/2", np.linalg.norm(1j * physical_ry @ basis - basis / 2))
    check("scalar neutral copies have T3=-1/2", np.linalg.norm(1j * physical_r3 @ basis + basis / 2))
    check("all 45 Spin(10) covariance equations of physical 10", covariance10)
    check("all 45 Spin(10) covariance equations of conjugated 126", covariance126)
    check("all color copies yield one common complex Yukawa coefficient", color_error)
    check("raw h reference coefficient sqrt(2)", abs(kappa_h - math.sqrt(2)))
    check("raw f Dirac reference magnitude 2/sqrt(3)", abs(abs(kappa_d) - 2/math.sqrt(3)))
    check("raw MR reference magnitude 4 sqrt(2)", abs(abs(kappa_r) - 4*math.sqrt(2)))
    check("MR/Dirac magnitude 2 sqrt(6)", abs(abs(kappa_r/kappa_d) - 2*math.sqrt(6)))
    check("all four species share the same phase-resolved Clebsch sum rules", max(value["sum_rule_residual"] for value in projected.values()))
    check("54 triplet has no renormalizable fermion coupling", abs(triplet_ll[0]) + abs(triplet_ll_conjugate[0]))
    check("126 LL triplet and MR share absolute normalization", abs(abs(triplet_ll[1]) - abs(kappa_r)))
    check("LL triplet couples only to physical Y=+1 amplitude", np.linalg.norm(triplet_ll_conjugate))
    check("direct real-coordinate response equals complex type-II contraction", abs(typeii_raw - direct_ml))
    check("bosonic-improved light vector is normalized", abs(np.vdot(light, light).real - 1))
    check("arbitrary scalar-copy rephasings preserve all species projections", phase_error)
    check("triplet basis rephasing preserves the complex type-II response", abs(rephased_typeii-typeii_raw))
    check("public real-scalar tensor interface reproduces type-II response", abs(tensor_ml-typeii_raw))
    check("random complex U4 copy rotation preserves actual Yukawa projections", copy_unitary_error)
    check("complete type I plus II covaries under random family U3", family_errors[0])
    check("complete type I plus II covaries under family rephasings", family_errors[1])
    check("independent nuL and nuc rephasings preserve seesaw covariance", neutrino_rephase_error)

    report = {
        "schema": "p54-common-complex-yukawa-phase-v1", "date": "2026-09-05",
        "status": "common complex scalar/spinor matching interface; no global fit performed",
        "conventions": {
            "global_conjugation": "physical x=K old x, flipping Im(Sigma), Im(phi), Im(S); U126physical=conj(U126old); standard spinor Gamma_star=+1; scalar and spinor dictionaries conjugated together",
            "scalar_action_pullback": "Vphysical(xphysical)=Vold(K xphysical); physical scalar Hessian=K Hold K, so the existing spectrum and complex copy eigenvector are transformed consistently",
            "couplings": "h_raw,f_raw denote coefficients in the globally conjugated action; they are complex conjugates of old-action coefficients if old couplings had been fixed",
            "Yukawa_action": "-L_Y=(1/2) h_raw Psi^T C Gamma_i Psi phi_physical_i* +(1/(2*5!)) f_raw Psi^T C Gamma_[5] Psi Sigma_physical_[5]+h.c.",
            "Higgs_definition": "x_physical=2 Re(B_physical c H0); H0 has Y=+1/2,T3=-1/2 and <H0>=v/sqrt(2)",
            "projection": "u,nu use c; d,e use conjugate(c), because H0 versus H0*",
            "spinor_phase": "Fock component vectors initially have largest entry real positive; right-handed-conjugate states then rephased once to make their reference 10 Yukawa entry positive; all f,MR,ML use those same states",
            "scalar_basis": "actual exported P1/P2 complex basis phases are preserved under real K, with no independent absolute-value replacement",
            "Dirac_unit": "h_D=kappa_h h_raw; f_D=kappa_d f_raw, with complex kappa_d recorded below",
            "Majorana_unit": "f_M=kappa_R f_raw; M_R=sigma_canonical*omega*f_M; f_M=(kappa_R/kappa_d) f_D",
            "phase_interpretation": "The displayed i in f_M/f_D is a component-basis convention, not an extra physical CP parameter. Independent nuL/nuc rephasings are transported through Ynu,MR,ML and verified on the full type I+II sum.",
            "mass_relations": "Yd=a h_D+d f_D; Ye=a h_D-3d f_D; Yu=b h_D+e f_D; Ynu=b h_D-3e f_D",
        },
        "scalar_copy_order_old_geometric": record["order"],
        "scalar_copy_basis_physical_328x4": cjson(basis),
        "bosonic_light_coefficients": cjson(light),
        "fermion_states_16": {key: cjson(value) for key, value in states_by_name.items()},
        "pairing_and_once_only_rephasings": pairing,
        "raw_copy_coefficients": {species: {key: cjson(value) for key, value in rows.items()}
                                  for species, rows in coefficients.items()},
        "normalization": {"kappa_h": cjson(kappa_h), "kappa_d": cjson(kappa_d),
                          "kappa_R": cjson(kappa_r), "f_M_over_f_D": cjson(kappa_r/kappa_d)},
        "projected_light_maps": projected,
        "typeII_common_phase": {
            "triplet_order": ["54", "126"], "canonical_LL_raw_coefficients": cjson(triplet_ll),
            "canonical_LL_conjugate_coefficients": cjson(triplet_ll_conjugate),
            "actual_complex_induced_amplitude_over_v2_times_omega": cjson(induced),
            "ML_over_v2_over_omega_per_f_raw": cjson(typeii_raw),
            "ML_over_v2_over_omega_per_f_D": cjson(typeii_raw/kappa_d),
            "ML_over_v2_over_omega_per_f_M": cjson(typeii_raw/kappa_r),
            "typeI_coefficient_in_v2_over_omega_units_f_D": cjson(-1/(2*sigma*(kappa_r/kappa_d))),
            "seesaw_formula": "Mnu=(v^2/omega)[CII_D f_D + CI_D Ynu f_D^{-1} Ynu^T], assuming invertible f_D; coefficients are complex and include the common MR/LL phase convention",
            "order": "tree triplet K and cubic source evaluated along fixed-VEV bosonic-CW light direction; this is mixed-order, not complete one-loop Weinberg matching",
        },
        "covariance_test_diagnostics": {
            "seed": 540926,
            "wrong_down_c_instead_of_conjugate_c_error": float(wrong_down_weight_error),
            "synthetic_family_matrices_scope": "used only for algebraic covariance tests, never for fit or matching-nuisance promotion",
        },
        "scope": {"complex_phase_matching_computed": True, "new_field_content": False,
                  "global_flavor_fit_performed": False, "finite_PS_Yukawa_thresholds": False,
                  "all_group_covariance": "all 45 Spin(10) generators checked again after the global scalar/spinor conjugation",
                  "consumer_API": "build_geometry() returns raw intertwiners; real_scalar_yukawa_tensors(g,old_coordinates=True) supplies (328,16,16) tensors directly on cached old scalar coordinates"},
        "checks": checks,
        "summary": {"passed": sum(row["pass"] for row in checks), "total": len(checks),
                    "all_pass": all(row["pass"] for row in checks)},
        "sources": [{"path": str(path.relative_to(RF.parent)),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for path in (Path(__file__), INTERTWINER, RF/"code/verify_p54_p1_hessian_spectrum.py", BOSON, TRIPLET)],
    }
    return report


def main():
    report = run()
    OUTPUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Common complex P54 Yukawa phase matching", "", report["status"], "",
             f"Checks: {report['summary']['passed']}/{report['summary']['total']}.", "",
             "The physical scalar basis and spinor chirality are conjugated together. Standard matter hypercharges are verified.", "",
             "The JSON exports the actual scalar basis, fixed spinor states, complex copy coefficients, complex light maps, and the common-phase type-II response.", "",
             "Normalization:", "", *[f"- {key}: {value}" for key, value in report["normalization"].items()], "",
             "Bosonic light map in Dirac-unit convention:", "",
             *[f"- {key}: {value}" for key, value in report["projected_light_maps"]["bosonic"]["Dirac_unit_abde"].items()], "",
             "The seesaw coefficients use these same phases and canonical P1 sigma. Finite PS matching and the global fit remain separate calculations.", ""]
    OUTPUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(report["summary"]))
    print(json.dumps(report["normalization"]))
    print(json.dumps(report["projected_light_maps"]["bosonic"], indent=2))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
