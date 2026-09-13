#!/usr/bin/env python3
"""Upper-PS heavy-vector finite Yukawa/kinetic SUBSET in DR, Landau gauge.

Actual P54 generators and cached stationary PS Hessian. Retained scalar
masses are expanded to leading power in the hard region; eliminated scalar
masses are kept exactly. This is not full upper/lower matching or a fit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.integrate import quad

import verify_p54_upper_yukawa_thresholds as yuk

RF = Path(__file__).resolve().parents[1]
LOOP = 16 * math.pi**2
OUT = RF / "output/p54_upper_vector_matching"


def b0(x, y, mu):
    """Finite B0(0;x,y), without 1/(16 pi^2), via stable quadrature."""
    if x <= 0 or y < 0 or mu <= 0:
        raise ValueError("hard x>0, internal y>=0 and mu>0 required")
    return -quad(lambda t: math.log(((1-t)*x+t*y)/mu**2), 0, 1,
                 epsabs=3e-12, epsrel=3e-12)[0]


def build_vector_geometry(geometry):
    ps = geometry["ps"]
    g = ps["common_geometry"]
    p1, ub, ids = g["p1"], ps["fermion_basis"], g["indices"]
    rs, ts = [], []
    for a in range(6):
        for b in range(6, 10):
            gen = np.zeros((10, 10))
            gen[a, b], gen[b, a] = 1/math.sqrt(2), -1/math.sqrt(2)
            rs.append(p1.representation_matrix(gen))
            spin = g["lift"](gen, g["gammas"])[np.ix_(ids, ids)]
            ts.append(1j * ub.conj().T @ spin @ ub)
    rs, ts = np.asarray(rs), np.asarray(ts)
    old = json.loads((RF / "output/p54_p2_two_site_matching.json").read_text())
    coupling = old["gauge_coupling_iteration"][-1]["g_output"]
    orbit = np.einsum("aij,j->ai", rs, geometry["x"])
    mass = coupling**2 * orbit @ orbit.T
    values, u = np.linalg.eigh(mass)
    # All masses are degenerate at this PS saddle. Keep the original
    # orthonormal generators; arbitrary diagonalizer phases cannot enter.
    if yuk.error(mass, np.eye(24)*np.mean(values)) > 1e-11:
        raise ValueError("This implementation requires the actual degenerate upper coset")
    active = geometry["active"]
    ra = np.einsum("Aij,ja->Aia", rs, active, optimize=True)
    cf = sum(t@t for t in ts)
    cs = np.einsum("Aia,Aib->ab", ra, ra, optimize=True)
    return {"R": rs, "T": ts, "RA": ra, "CF": cf, "CS": cs,
            "g": coupling, "mass2": float(np.mean(values)), "mass_matrix": mass}


def kinetic(geometry, vector, mu):
    mv = vector["mass2"]
    heavy, ra = geometry["heavy"], vector["RA"]
    w0 = 3*b0(mv, 0., mu)-.5
    wh = np.array([3*b0(mv, float(m), mu)-.5 for m in geometry["masses2"]])
    # Complete state sum: complement is the retained/Goldstone/PQ block.
    # Its mass expansion is explicit, not a claim that all these masses vanish.
    projected = np.einsum("ij,Aia->Aja", heavy, ra, optimize=True)
    weighted = w0*vector["CS"] + np.einsum(
        "j,Aja,Ajb->ab", wh-w0, projected, projected, optimize=True)
    kphi = -vector["g"]**2/LOOP * weighted
    kpsi = -1.5*vector["g"]**2/LOOP * vector["CF"]
    return kpsi, (kphi+kphi.T)/2


def threshold(tree, geometry, vector, mu):
    nf = tree.shape[-1]//16
    kp, ks = kinetic(geometry, vector, mu)
    kp = np.kron(kp, np.eye(nf))
    # The -2 and -1/2 constants arise from d-dimensional numerator factors
    # times the UV pole; no 4D-numerator shortcut is made here.
    weight = 3*b0(vector["mass2"], 0., mu)-2
    vertex = np.zeros_like(tree)
    for t in vector["T"]:
        tf = np.kron(t, np.eye(nf))
        vertex -= vector["g"]**2/LOOP * weight * (tf.T@tree@tf)
    fermion_legs = -.5*(kp.T@tree+tree@kp)
    scalar_legs = -.5*np.einsum("ab,bij->aij", ks, tree, optimize=True)
    return {"Kpsi": kp, "Kphi": ks, "vertex": vertex,
            "fermion_legs": fermion_legs, "scalar_legs": scalar_legs,
            "delta": vertex+fermion_legs+scalar_legs}


def require_complete_matching(payload):
    """No partial diagram package, or mismatched scheme, may enable a fit."""
    # One shared consumer guard, not a competing list of fit requirements.
    from verify_p54_finite_yukawa_interface import require_physical_inputs
    require_physical_inputs(payload)


def run(cache_dir):
    geometry = yuk.build_upper_geometry(cache_dir)
    vector = build_vector_geometry(geometry)
    mu = geometry["mu"]
    mirror = yuk.build_mirror_basis(geometry)
    checks = []

    def check(name, a, b=0., tol=2e-10):
        e = yuk.error(np.asarray(a), np.asarray(b))
        checks.append({"name": name, "residual": e, "tolerance": tol, "pass": e < tol})

    def rejects(name, value):
        try:
            require_complete_matching(value)
        except ValueError:
            check(name, 0.)
        else:
            check(name, 1.)

    check("upper_coset_24_degenerate_masses", vector["mass_matrix"], np.eye(24)*vector["mass2"])
    check("coset_spinor_Casimir_3", vector["CF"], 3*np.eye(16))
    check("all_upper_scalar_masses_positive", float(geometry["masses2"].min() <= 0))
    check("active_heavy_tree_mass_block_zero", geometry["active"].T@geometry["hessian"]@geometry["heavy"])
    check("Landau_Goldstone_Yukawa_h_zero", np.einsum("ab,aij->bij", geometry["goldstone"], geometry["all_h"]))
    check("Landau_Goldstone_Yukawa_f_zero", np.einsum("ab,aij->bij", geometry["goldstone"], geometry["all_f"]))
    # A single active scalar has no VV or ghost mass vertex at sigma=0:
    # the background is in 54+S, the active parents are in 10+126.
    rx = np.einsum("Aij,j->Ai", vector["R"], geometry["x"])
    check("active_linear_vector_mass_vertex_zero",
          np.einsum("Ai,Bia->ABa", rx, vector["RA"]))
    check("active_vector_current_no_eaten_scalar",
          np.einsum("ij,Aia->Aja", geometry["goldstone"], vector["RA"]))
    expected = np.repeat([3., 7., 6., 6.], [8, 120, 60, 60])
    check("actual_parent_coset_Casimirs", vector["CS"], np.diag(expected))

    for x, y in ((.37, 0.), (.37, .37), (.37, 2.4), (2.4, .37)):
        exact = 1-math.log(x/mu**2) if y == 0 else (
            -math.log(x/mu**2) if x == y else
            1-(x*math.log(x/mu**2)-y*math.log(y/mu**2))/(x-y))
        check(f"B0_quadrature_{x}_{y}", b0(x,y,mu), exact)
    # Laurent finite parts reconstructed at nonzero epsilon, independently
    # of the group sums. Coefficients correspond to d=4-2 epsilon.
    e, b = 1e-5, .731
    d = 4-2*e
    check("d_dim_vertex_finite_minus2", (d-1)*(1/e+b)-3/e, 3*b-2, tol=2e-5)
    check("d_dim_scalar_finite_minus_half", 4*(d-1)/d*(1/e+b)-3/e, 3*b-.5, tol=1e-5)
    check("d_dim_fermion_finite_minus_three_halves", -(5-d-4/d)*(1/e+b), -1.5, tol=1e-5)
    # Independent neutral Dirac mass insertion: Sigma_m+Sigma_slash.
    # The longitudinal contribution cancels in the matching mass.
    for xi in (0., .3, 1., 2.):
        log = .417
        long = 0. if xi == 0 else xi*(1-log-math.log(xi))
        sigma_mass = 1-3*log+long
        sigma_slash = 1.5-long
        check(f"Abelian_neutral_mass_Rxi_cancellation_{xi}", sigma_mass+sigma_slash, 2.5-3*log)

    kpsi, kphi = kinetic(geometry, vector, mu)
    check("Kpsi_Hermitian", kpsi, kpsi.conj().T)
    check("Kphi_real_symmetric", kphi, kphi.T)
    check("positive_fermion_kinetic_metric", float(np.linalg.eigvalsh(np.eye(16)+kpsi).min() <= 0))
    check("positive_scalar_kinetic_metric", float(np.linalg.eigvalsh(np.eye(248)+kphi).min() <= 0))
    ps = geometry["ps"]
    parents = yuk.module("vector_parent_generators", RF / "code/verify_p54_p2_two_site_matching.py")
    commutator = 0.
    for group in parents.ps_generators(ps["common_geometry"]["p1"]).values():
        for gen in group:
            r = ps["common_geometry"]["p1"].representation_matrix(gen)
            ar = geometry["active"].T@r@geometry["active"]
            commutator = max(commutator, yuk.error(kphi@ar, ar@kphi))
    check("PS_covariance_of_scalar_kinetic", commutator)

    coefficients, norm_rows = {}, {}
    for spurion in ("h", "f"):
        h, f = np.array([[float(spurion == "h")]], complex), np.array([[float(spurion == "f")]], complex)
        tree = geometry["flow"].assemble_real_yukawas(geometry["flow"].spin10_boundary(h,f), ps)
        actual_tree = np.einsum("ab,aij->bij", geometry["active"], geometry["all_"+spurion])
        check(f"{spurion}_tree_intertwiner_dictionary", tree, actual_tree)
        # Ward-derived Casimir identity, independent of loop kernels.
        x = sum(t.T@tree@t for t in vector["T"])
        ward = np.einsum("ab,bij->aij",vector["CS"],tree)-vector["CF"].T@tree-tree@vector["CF"]
        check(f"{spurion}_coset_Yukawa_Ward", 2*x, ward)
        result = threshold(tree,geometry,vector,mu)
        projection = yuk.project_declared(result["delta"],geometry)
        mir = yuk.project_mirror(projection["residual"],geometry,mirror)
        reconstruction = projection["reconstruction"]+yuk.assemble_mirror(mir,geometry,mirror)
        check(f"{spurion}_six_invariant_closure",result["delta"],reconstruction)
        coefficients[spurion] = {"direct":{k:yuk.cjson(v) for k,v in projection["couplings"].items()},
                                 "mirror":{k:yuk.cjson(v) for k,v in mir.items()}}
        norm_rows[spurion] = {k:float(np.linalg.norm(result[k])) for k in ("vertex","fermion_legs","scalar_legs","delta")}
        step = 1e-3
        dp = threshold(tree,geometry,vector,mu*math.exp(step))["delta"]
        dm = threshold(tree,geometry,vector,mu*math.exp(-step))["delta"]
        # beta_PS-beta_Spin10=+3 g^2 {C_coset,Y}/(16 pi^2).
        target = 3*vector["g"]**2/LOOP*(vector["CF"].T@tree+tree@vector["CF"])
        check(f"{spurion}_full_gauge_matching_scale_derivative",(dp-dm)/(2*step),target)
        no_scalar = result["vertex"]+result["fermion_legs"]
        check(f"{spurion}_omitting_scalar_legs_is_detectable",float(np.linalg.norm(result["delta"]-no_scalar)<1e-7))

    scheme = {"action":"P54PQ-v2","subtraction":"MSbar-DR",
              "gauge":"background-field Landau","tadpoles":"fixed-VEV"}
    rejects("fit_rejects_this_partial_package", scheme)
    rejects("fit_rejects_scheme_switch", {**scheme,"subtraction":"DRbar"})
    old_flags = {key:True for key in ("upper_finite_gauge","lower_finite_gauge",
        "upper_finite_Yukawa","lower_finite_Yukawa","all_active_PS_flow",
        "lower_EFT_flow","sequential_Weinberg_matching","same_action_scalar_feedback")}
    rejects("old_boolean_flags_do_not_bypass_KJC", {**scheme,**old_flags})
    active_m = np.linalg.eigvalsh(geometry["active"].T@geometry["hessian"]@geometry["active"])
    max_ratio = float(np.max(abs(active_m))/vector["mass2"])
    sources = [Path(__file__),Path(yuk.__file__),RF/"code/verify_p54_p1_hessian_spectrum.py",
               RF/"code/verify_p54_common_yukawa_phase.py",RF/"output/p54_ps_finite_thresholds.json",
               RF/"output/p54_p2_two_site_matching.json",RF/"code/verify_p54_finite_yukawa_interface.py"]
    report = {"schema":"p54-upper-vector-subset-v1", "date":"2026-09-11",
              "status":"upper vector Yukawa/kinetic hard subset; no full matching or fit",
              "scheme":scheme,"vacuum":geometry["upper"]["vacuum"],
              "cache_path":geometry["cache_path"],
              "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
              "g":vector["g"],"mu_over_omega_reference":mu,"vector_mass2":vector["mass2"],
              "coset_spinor_Casimir":3,"coset_parent_Casimirs":{"H":3,"F":7,"L":6,"R":6},
              "Kpsi":yuk.cjson(kpsi),"Kphi":yuk.cjson(kphi),
              "scalar_kinetic_eigenvalue_range":[float(x) for x in np.linalg.eigvalsh(kphi)[[0,-1]]],
              "fermion_kinetic_eigenvalue_range":[float(x) for x in np.linalg.eigvalsh(kpsi)[[0,-1]]],
              "finite_linear_spurion_coefficients":coefficients,"tensor_norms":norm_rows,
              "retained_mass_expansion":{"order":"leading retained mass / vector mass in hard region",
                 "largest_abs_m2_over_MV2":max_ratio,"error_bound_established":False,
                 "all_retained_masses_zero_claim":False},
              "full_matching_complete":False,"physical_fit_enabled":False,
              "missing":["upper scalar-cubic kinetic, potential/tadpole and Wilson matching",
                 "upper full source/contact hard-minus-EFT matching and operator mixing",
                 "controlled retained-mass corrections for the stated accuracy",
                 "lower covariant complete finite matching","finite sequential seesaw matching",
                 "fermionic light eigenstate feedback with matched families"],
              "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),
                 "all_pass":all(c["pass"] for c in checks)}}
    return report, geometry, vector


def markdown(report):
    lines = ["# P54 upper heavy-vector finite matching subset", "",
             "The existing upper PS saddle and generators are used. **This is not full (K,J,C) matching and not a physical fit.**", "",
             "DR, MSbar, background-field Landau; fixed-VEV prescription. One-loop vertices and kinetic terms only.", "",
             f"Checks: {report['summary']['passed']}/{report['summary']['total']}. "
             f"M_V^2 = {report['vector_mass2']:.12g}; mu = {report['mu_over_omega_reference']:.12g}; g = {report['g']:.12g}.", "",
             f"Kpsi eigenvalue range: {report['fermion_kinetic_eigenvalue_range']}. "
             f"Kphi eigenvalue range: {report['scalar_kinetic_eigenvalue_range']}.", "",
             "Eliminated scalar masses are kept exactly in the vector-scalar bubble; retained scalar masses are expanded at leading hard order. "
             f"The largest |m_retained^2|/M_V^2 is {report['retained_mass_expansion']['largest_abs_m2_over_MV2']:.6g}; no precision error bound is claimed.", "",
             "The -2, -1/2 and -3/2 finite terms in the vertex, scalar and fermion formulas are derived before taking d=4. "
             "The two spurion tensors satisfy the coset Ward identity and the complete gauge scale derivative; a scale check alone cannot fix finite constants.", "",
             "## Finite linear coefficients in raw family units", "",
             "| Spurion | Parent | Direct | Mirror |", "|---|---|---:|---:|"]
    for s, c in report["finite_linear_spurion_coefficients"].items():
        for key, value in c["direct"].items():
            v = complex(value["real"][0][0],value["imag"][0][0])
            m = c["mirror"].get(key,{"real":[[0]],"imag":[[0]]})
            mv = complex(m["real"][0][0],m["imag"][0][0])
            lines.append(f"| {s}_raw | {key} | {v:.10g} | {mv:.10g} |")
    lines += ["", "These are action-derived group contractions, not fitted observables or additional UV family parameters.", "",
              "## Remaining requirements", ""]+["- "+x for x in report["missing"]]
    lines += ["", "## Regressions", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in report["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir",type=Path)
    args = parser.parse_args()
    report, _, _ = run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(report))
    print(json.dumps({"summary":report["summary"],"mV2":report["vector_mass2"],
                      "Kphi":report["scalar_kinetic_eigenvalue_range"],
                      "failed":[c for c in report["checks"] if not c["pass"]]},indent=2))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)
