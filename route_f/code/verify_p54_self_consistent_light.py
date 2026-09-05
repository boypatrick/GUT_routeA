#!/usr/bin/env python3
"""Clebsch-constrained local scalar/flavor feedback on two synthetic families.

The bosonic background and curvature are the actual existing P54 artifacts.
The family matrices are deterministic examples, NOT fitted parameters. The
fermionic fixed-VEV correction uses the common complex Clifford dictionary.
For each new light direction the tree triplet source is recomputed as its
quadratic contraction with same-action cubic jets. No old CII is frozen in.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm

RF = Path(__file__).resolve().parents[1]
P1 = RF / "code/verify_p54_p1_hessian_spectrum.py"
PHASE = RF / "output/p54_common_yukawa_phase.json"
BOSON = RF / "output/p54_full_doublet_cw.json"
TRIPLET = RF / "output/p54_typeii_triplet_source.json"
FERMION = RF / "code/verify_p54_fermion_tadpole.py"
CACHE = RF.parent / "tmp/p54_full_doublet_cw"
OUT = RF / "output/p54_self_consistent_light"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(obj)
    return obj


def cj(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def decode(x):
    return np.asarray(x["real"]) + 1j*np.asarray(x["imag"])


def herm(x):
    return (x+x.conj().T)/2


def align(vector, reference):
    return vector*np.exp(-1j*np.angle(np.vdot(reference, vector)))


def cubic_triplet_tensor(p1, boson, phase, triplet):
    """Build T[A,p,q]=V'''[triplet_A,doublet_p,doublet_q] from cached jets.

    P54's potential is quartic, hence centered differences of its Hessian
    yield the exact cubic tensor algebraically (up to floating-point errors).
    Hypercharge rotates real into imaginary Higgs components at the invariant
    background, supplying the other four directional jets exactly.
    """
    parameters = boson["tree_parameters"]
    vevs = boson["vacuum"]
    x0 = p1.vacuum_vector(vevs["omega"], vevs["sigma"], vevs["vs"])
    physical_basis = decode(phase["scalar_copy_basis_physical_328x4"])
    k = np.ones(p1.N_REAL)
    k[p1.SL_SIGMA_IM] = k[p1.SL_H_IM] = -1
    k[327] = -1
    basis = k[:, None]*physical_basis
    qr, qi = math.sqrt(2)*basis.real, -math.sqrt(2)*basis.imag
    doublet_real = np.column_stack([qr, qi])
    triplet_complex = decode(triplet["complex_triplet_basis_328x2"])
    triplet_real = np.column_stack([math.sqrt(2)*triplet_complex.real,
                                    -math.sqrt(2)*triplet_complex.imag])
    keybase = hashlib.sha256(P1.read_bytes()+json.dumps(parameters, sort_keys=True).encode()).digest()
    counts = {"cache_hits": 0, "evaluated": 0}
    hessian_function = None
    def hessian(x):
        nonlocal hessian_function
        key = hashlib.sha256(keybase+np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = CACHE/(key+".npz")
        if path.exists():
            with np.load(path) as data:
                result = data["h"]
            counts["cache_hits"] += 1
            return result
        if hessian_function is None:
            hessian_function = p1.jax.jit(p1.hessian(p1.potential_factory(parameters)))
        start = time.time()
        result = np.asarray(hessian_function(p1.anp.asarray(x)), dtype=float)
        result = (result+result.T)/2
        CACHE.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, h=result)
        counts["evaluated"] += 1
        print(f"New scalar Hessian jet: {time.time()-start:.2f}s", flush=True)
        return result
    step = .02
    jets = [(hessian(x0+step*q)-hessian(x0-step*q))/(2*step) for q in qr.T]
    ry = p1.representation_matrix(p1.sm_generators()["Y"][0])
    rotation = expm(-math.pi*ry)
    rotation_error = float(np.linalg.norm(rotation@qr-qi))
    jets.extend([rotation@jet@rotation.T for jet in jets[:4]])
    tensor = np.einsum("iA,pij,jq->Apq", triplet_real, np.asarray(jets), doublet_real)
    symmetry_error = float(np.linalg.norm(tensor-tensor.transpose(0, 2, 1)))
    tensor = (tensor+tensor.transpose(0, 2, 1))/2
    c0 = decode(phase["bosonic_light_coefficients"])
    real_c = np.r_[c0.real, c0.imag]
    source = np.einsum("Apq,p,q->A", tensor, real_c, real_c)
    old_source = np.asarray(triplet["canonical_cubic_source_over_omega"])
    source_error = float(np.linalg.norm(source-old_source)/max(np.linalg.norm(old_source), 1e-15))
    return tensor, {"cache": counts, "step": step,
                    "real_imaginary_rotation_error": rotation_error,
                    "cubic_tensor_permutation_error": symmetry_error,
                    "existing_actual_source_reconstruction_relative_error": source_error}


def run():
    p1 = module("p54_feedback_p1", P1)
    fermion = module("p54_feedback_fermion", FERMION)
    audit = fermion.helpers()
    phase, boson, triplet = [json.loads(path.read_text()) for path in (PHASE, BOSON, TRIPLET)]
    tensor, jet_audit = cubic_triplet_tensor(p1, boson, phase, triplet)
    m0 = decode(boson["tree_doublet_matrix_over_omega2"])
    correction_b = decode(boson["bosonic_CW"]["total_matrix_over_omega2"])
    correction_b += decode(boson["tadpoles"]["doublet_CT_over_omega2"])
    cb = decode(phase["bosonic_light_coefficients"])
    p126 = np.diag([0., 0., 1., 1.]).astype(complex)
    p10 = np.eye(4)-p126
    sigma = float(boson["vacuum"]["sigma"])
    mu = float(boson["scheme"]["mu_over_omega"])
    kappa_r = complex(decode(phase["normalization"]["kappa_R"]))
    kappa_d = complex(decode(phase["normalization"]["kappa_d"]))
    coefficients = {species: {key: decode(value) for key, value in rows.items()}
                    for species, rows in phase["raw_copy_coefficients"].items()}
    kt = np.asarray(triplet["real_triplet_K_over_omega2"])
    ll = decode(phase["typeII_common_phase"]["canonical_LL_raw_coefficients"])
    llc = decode(phase["typeII_common_phase"]["canonical_LL_conjugate_coefficients"])

    def source_and_typeii(c):
        cr = np.r_[c.real, c.imag]
        source = np.einsum("Apq,p,q->A", tensor, cr, cr)
        response = -.5*np.linalg.solve(kt, source)
        amplitude = (response[:2]+1j*response[2:])/math.sqrt(2)
        cii = complex(ll@amplitude+llc@np.conj(amplitude))
        return source, response, amplitude, cii

    def solve(hraw, fraw, *, tree=m0, bosonic=correction_b, projector=p10,
              p_126=p126, reference=cb, coeff=coefficients["nu"], k_r=kappa_r):
        ys = [a*hraw+b*fraw for a, b in zip(coeff["h_raw"], coeff["f_raw"])]
        fixed = fermion.fixed_vev_fermion(ys, k_r*fraw, sigma, mu, p_126, audit)
        tuned = fermion.schur_retune(tree, bosonic+fixed["fixed"], projector)
        if not tuned["success"]:
            raise RuntimeError("Benchmark lacks a stable Schur-retuned light eigenpair")
        c = align(decode(tuned["light_coefficients"]), reference)
        return {"ys": ys, "fixed": fixed, "retuning": tuned, "c": c,
                "matrix": decode(tuned["matrix"])}

    def exact_matrix_derivative(hraw, fraw, result):
        """d/dtheta at zero for f_raw(theta)=(1+theta) f_raw.

        X'=2X, W=X(log(X/mu^2)-1), W'=2X log(X/mu^2).
        No finite difference is used in this derivative of the CW matrix.
        """
        mr = result["fixed"]["mr"]
        eigenvalues, u = np.linalg.eigh(mr.conj().T@mr)
        logs = np.log(eigenvalues/mu**2)
        w = (u*(eigenvalues*(logs-1)))@u.conj().T
        dw = (u*(2*eigenvalues*logs))@u.conj().T
        ys = result["ys"]
        dys = [b*fraw for b in coefficients["nu"]["f_raw"]]
        dpi = np.array([[-np.trace((da.conj().T@b+a.conj().T@db)@w+
                                   a.conj().T@b@dw)/(8*math.pi**2)
                         for b, db in zip(ys, dys)] for a, da in zip(ys, dys)])
        dnu2 = -float(np.sum(eigenvalues**2*(4*logs-2)))/(8*math.pi**2*sigma**2)
        return herm(dpi-dnu2*p126/2)

    rng = np.random.default_rng(540930)
    def unitary(n):
        return np.linalg.qr(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))[0]
    rows = []
    checks = []
    def check(name, error, tolerance=1e-10):
        checks.append({"name": name, "residual": float(error), "tolerance": tolerance,
                       "pass": bool(error < tolerance)})
    check("common phase prerequisite passed", 0 if phase["summary"]["all_pass"] else 1)
    check("cached same-action jets reproduce actual previous source", jet_audit["existing_actual_source_reconstruction_relative_error"], 1e-9)
    check("cubic tensor is symmetric in its Higgs indices", jet_audit["cubic_tensor_permutation_error"], 1e-9)
    check("hypercharge transports all imaginary directional jets", jet_audit["real_imaginary_rotation_error"])
    old_cii = source_and_typeii(cb)[-1]
    for case, singular_values in enumerate(([.025, .04, .065], [.035, .055, .085])):
        u = unitary(3)
        fraw = u.conj()@np.diag(singular_values)@u.conj().T
        hraw = .12*(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
        hraw = (hraw+hraw.T)/2
        result = solve(hraw, fraw)
        c, matrix = result["c"], result["matrix"]
        source, response, amplitude, cii = source_and_typeii(c)
        ynu = sum(ci*yi for ci, yi in zip(c, result["ys"]))
        species_yukawas = {}
        species_projection = {}
        for species, row in coefficients.items():
            weights = c if species in ("u", "nu") else c.conj()
            ah, af = complex(row["h_raw"]@weights), complex(row["f_raw"]@weights)
            species_yukawas[species] = ah*hraw+af*fraw
            species_projection[species] = {"h_raw": cj(ah), "f_raw": cj(af)}
        ci = -1/(2*sigma*kappa_r)
        ml = cii*fraw
        mi = ci*ynu@np.linalg.solve(fraw, ynu.T)
        mnu = ml+mi
        derivative = exact_matrix_derivative(hraw, fraw, result)
        dxi = float(np.vdot(c, derivative@c).real/np.vdot(c, p10@c).real)
        eigenvalues, vectors = np.linalg.eigh(matrix)
        pseudo = (vectors[:, 1:]*(1/eigenvalues[1:]))@vectors[:, 1:].conj().T
        dc = -pseudo@(derivative-dxi*p10)@c
        difference_checks = []
        for step in (1e-3, 3e-4):
            plus, minus = solve(hraw, (1+step)*fraw), solve(hraw, (1-step)*fraw)
            dc_fd = (align(plus["c"], c)-align(minus["c"], c))/(2*step)
            dxi_fd = (plus["retuning"]["delta_xi02"]-minus["retuning"]["delta_xi02"])/(2*step)
            dm_fd = (plus["fixed"]["fixed"]-minus["fixed"]["fixed"])/(2*step)
            difference_checks.append({
                "step": step,
                "implicit_vector_derivative_relative_error": float(np.linalg.norm(dc_fd-dc)/max(np.linalg.norm(dc), 1e-14)),
                "retuning_derivative_relative_error": abs(dxi_fd-dxi)/max(abs(dxi), 1e-14),
                "CW_matrix_derivative_relative_error": float(np.linalg.norm(dm_fd-derivative)/max(np.linalg.norm(derivative), 1e-14)),
            })
        # The new source is evaluated quadratically in c. Its implicit
        # derivative comes from the same tensor and the horizontal dc.
        cr, dcr = np.r_[c.real, c.imag], np.r_[dc.real, dc.imag]
        dj = 2*np.einsum("Apq,p,q->A", tensor, dcr, cr)
        test_step = 1e-3
        cplus = align(solve(hraw, (1+test_step)*fraw)["c"], c)
        cminus = align(solve(hraw, (1-test_step)*fraw)["c"], c)
        dj_fd = (source_and_typeii(cplus)[0]-source_and_typeii(cminus)[0])/(2*test_step)
        source_derivative_error = float(np.linalg.norm(dj_fd-dj)/max(np.linalg.norm(dj), 1e-14))

        family = unitary(3)
        moved = solve(family.T@hraw@family, family.T@fraw@family)
        moved_c = align(moved["c"], c)
        moved_source, _, _, moved_cii = source_and_typeii(moved_c)
        moved_ynu = sum(z*y for z, y in zip(moved_c, moved["ys"]))
        moved_f = family.T@fraw@family
        moved_mass = moved_cii*moved_f+ci*moved_ynu@np.linalg.solve(moved_f, moved_ynu.T)
        family_mass_error = float(np.linalg.norm(moved_mass-family.T@mnu@family)/max(np.linalg.norm(mnu), 1e-14))
        cp_coeff = {key: value.conj() for key, value in coefficients["nu"].items()}
        cp = solve(hraw.conj(), fraw.conj(), tree=m0.conj(), bosonic=correction_b.conj(),
                   projector=p10.conj(), p_126=p126.conj(), reference=c.conj(),
                   coeff=cp_coeff, k_r=kappa_r.conjugate())
        cp_vector_error = float(np.linalg.norm(cp["c"]-c.conj()))
        # This is covariance under conjugation of every complex input, not
        # an assertion that the original parameter point is CP invariant.
        cp_doublet_sign = np.r_[np.ones(4), -np.ones(4)]
        cp_triplet_sign = np.r_[np.ones(2), -np.ones(2)]
        cp_tensor = (tensor * cp_triplet_sign[:, None, None]
                     * cp_doublet_sign[None, :, None] * cp_doublet_sign[None, None, :])
        cp_kt = cp_triplet_sign[:, None]*kt*cp_triplet_sign[None, :]
        cp_cr = np.r_[cp["c"].real, cp["c"].imag]
        cp_source = np.einsum("Apq,p,q->A", cp_tensor, cp_cr, cp_cr)
        cp_response = -.5*np.linalg.solve(cp_kt, cp_source)
        cp_amplitude = (cp_response[:2]+1j*cp_response[2:])/math.sqrt(2)
        cp_cii = complex(ll.conj()@cp_amplitude+llc.conj()@cp_amplitude.conj())
        cp_ynu = sum(z*y for z, y in zip(cp["c"], cp["ys"]))
        cp_mass = cp_cii*fraw.conj()+ci.conjugate()*cp_ynu@np.linalg.solve(fraw.conj(), cp_ynu.T)
        cp_mass_error = float(np.linalg.norm(cp_mass-mnu.conj())/max(np.linalg.norm(mnu), 1e-14))

        prefix = f"case {case}: "
        check(prefix+"stable unique Schur-retuned light mode", 0 if eigenvalues[1]>1e-3 else 1)
        check(prefix+"full eigenpair equation", np.linalg.norm(matrix@c))
        check(prefix+"exact Schur complement", result["retuning"]["Schur_residual"])
        check(prefix+"horizontal derivative c-dagger dc=0", abs(np.vdot(c, dc)))
        check(prefix+"differentiated eigenpair equation", np.linalg.norm(matrix@dc+(derivative-dxi*p10)@c))
        check(prefix+"analytic CW derivative agrees with both differences", max(v["CW_matrix_derivative_relative_error"] for v in difference_checks), 1e-5)
        check(prefix+"implicit c derivative agrees with both differences", max(v["implicit_vector_derivative_relative_error"] for v in difference_checks), 1e-4)
        check(prefix+"implicit retuning derivative agrees with both differences", max(v["retuning_derivative_relative_error"] for v in difference_checks), 1e-4)
        check(prefix+"new quadratic triplet source solves its field equation", np.linalg.norm(kt@response+source/2))
        check(prefix+"new triplet source derivative follows implicit c", source_derivative_error, 1e-4)
        check(prefix+"family covariance of corrected light direction", np.linalg.norm(moved_c-c))
        check(prefix+"family covariance of I plus II mass", family_mass_error)
        check(prefix+"full-input complex conjugation of light direction", cp_vector_error)
        check(prefix+"CP-transformed cubic tensor recomputes conjugate source", np.linalg.norm(cp_source-cp_triplet_sign*source))
        check(prefix+"CP-transformed triplet response recomputes conjugate CII", abs(cp_cii-cii.conjugate()))
        check(prefix+"full-input complex conjugation of seesaw expression", cp_mass_error)
        rows.append({
            "case": case, "scope": "synthetic family matrices constrained by actual P54 Clebsches; not a fit",
            "h_raw": cj(hraw), "f_raw": cj(fraw), "f_raw_Takagi_values": list(singular_values),
            "h_raw_operator_norm": float(np.linalg.norm(hraw, 2)), "f_raw_operator_norm": float(np.linalg.norm(fraw, 2)),
            "actual_neutrino_copy_Yukawas": [cj(value) for value in result["ys"]],
            "MR_over_omega": cj(result["fixed"]["mr"]),
            "fermionic_CW_over_omega2": cj(result["fixed"]["pi"]),
            "fermionic_tadpole_CT_over_omega2": cj(result["fixed"]["ct"]),
            "fermionic_fixed_VEV_correction_over_omega2": cj(result["fixed"]["fixed"]),
            "retuning": result["retuning"], "common_phase_light_coefficients": cj(c),
            "rotation_from_bosonic_light_radians": math.acos(min(1., abs(np.vdot(cb, c)))),
            "implicit_parameter": "f_raw(theta)=(1+theta) f_raw, h_raw held fixed; derivative at theta=0",
            "d_fixed_CW_dtheta_over_omega2": cj(derivative), "d_delta_xi02_dtheta": dxi,
            "dc_dtheta_horizontal": cj(dc), "finite_difference_checks": difference_checks,
            "tree_triplet_source_at_new_c_over_omega": source.tolist(),
            "tree_triplet_real_response_over_v2_times_omega": response.tolist(),
            "tree_triplet_complex_response_over_v2_times_omega": cj(amplitude),
            "CII_raw_at_new_c": cj(cii), "CII_D_at_new_c": cj(cii/kappa_d),
            "CII_raw_at_old_bosonic_c": cj(old_cii), "CII_change_absolute": abs(cii-old_cii),
            "Ynu_at_new_c": cj(ynu),
            "all_four_Yukawa_matrices_at_new_c": {name: cj(value) for name, value in species_yukawas.items()},
            "all_four_complex_scalar_projections_at_new_c": species_projection,
            "ML_over_v2_over_omega": cj(ml), "MI_over_v2_over_omega": cj(mi),
            "Mnu_over_v2_over_omega": cj(mnu),
            "Mnu_Takagi_values_in_v2_over_omega_units": np.linalg.svd(mnu, compute_uv=False).tolist(),
            "source_derivative_relative_error": source_derivative_error,
        })
    return {
        "schema": "p54-clebsch-local-scalar-flavor-feedback-v1", "date": "2026-09-05",
        "scope": "two deterministic perturbative synthetic families with actual complex Spin(10) Clebsches and actual bosonic curvature; local feedback only, no finite thresholds or global fit",
        "fixed_background": {"sigma_over_omega": sigma, "mu_over_omega": mu,
                             "scheme": "fixed mu and fixed VEV fermionic radial tadpole subtraction"},
        "formulae": {
            "CW_derivative": "X'=2X; W=X(log(X/mu^2)-1); W'=2X log(X/mu^2); differentiate each actual Ya and the same radial counterterm",
            "retuning_derivative": "xi'=(c^dagger Delta_F' c)/(c^dagger P10 c)",
            "implicit_light_derivative": "c'=-D^+ (Delta_F'-xi' P10)c, with c^dagger c'=0 and D^+ on the positive heavy subspace",
            "triplet_source": "J_A(c)=T_Apq c_real[p]c_real[q]; J_A'=2T_Apq c_real'[p]c_real[q]",
            "triplet_response": "response=-K_tree^(-1)J(c)/2; complex amplitude=(response_R+i response_I)/sqrt(2)",
        },
        "cubic_jet_audit": jet_audit, "cubic_tensor_4x8x8": tensor.tolist(),
        "cases": rows,
        "remaining_physical_gates": [
            "The synthetic families are not fitted to charged masses, CKM, PMNS or seesaw observations.",
            "The complete finite PS matching, full multi-parent Yukawa flow and sequential Weinberg thresholds remain necessary.",
            "The updated triplet source is tree-level at each new c; triplet K is the baseline tree block. Bosonic and fermionic loop triplet masses/cubic vertices are not claimed.",
            "All background VEVs are held by the declared tadpole prescription; this is not a new globally minimized quantum vacuum.",
            "Complex-conjugation tests transform all complex inputs; they do not establish CP symmetry of a benchmark.",
        ],
        "physical_fit_performed": False,
        "checks": checks, "summary": {"passed": sum(row["pass"] for row in checks),
                                        "total": len(checks), "all_pass": all(row["pass"] for row in checks)},
        "sources": [{"path": str(path.relative_to(RF.parent)),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for path in (Path(__file__), P1, FERMION, RF/"code/verify_p54_theory_audit.py", PHASE, BOSON, TRIPLET)],
    }


def main():
    result = run()
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# P54 local scalar/flavor feedback", "", result["scope"], "",
             f"Checks: {result['summary']['passed']}/{result['summary']['total']}.", "",
             "Both examples use actual four-copy neutrino Clebsches, MR normalization and fixed-VEV fermionic corrections. The Schur-retuned light vector is recomputed for each family.", "",
             "The type-II tree source is recomputed from same-action cubic jets at the new light direction; it is not frozen to its previous value.", ""]
    for row in result["cases"]:
        lines.extend([f"Case {row['case']}: heavy gap {row['retuning']['nearest_heavy_gap']:.10g}; total delta_xi02 {row['retuning']['delta_xi02']:.10g}; rotation {row['rotation_from_bosonic_light_radians']:.6g} rad; |delta CII_raw| {row['CII_change_absolute']:.6g}.", ""])
    lines.extend(["Physical gates:", "", *[f"- {text}" for text in result["remaining_physical_gates"]], ""])
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"]))
    print(json.dumps(result["cubic_jet_audit"]))
    for row in result["cases"]:
        print(json.dumps({"case": row["case"], "gap": row["retuning"]["nearest_heavy_gap"],
                          "dc_dtheta_norm": float(np.linalg.norm(decode(row["dc_dtheta_horizontal"]))),
                          "FD": row["finite_difference_checks"], "CII_change": row["CII_change_absolute"]}))
    for row in result["checks"]:
        if not row["pass"]:
            print("FAIL", row)
    if not result["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
