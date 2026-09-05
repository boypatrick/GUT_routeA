#!/usr/bin/env python3
"""Controlled finite gauge-threshold diagnostics for the four-parent PS EFT.

The upper PS-symmetric stationary saddle uses the same scalar parameters as
the historical broken vacuum. Tachyons in retained active fields are not
mistaken for a failure of the integrated heavy block. A separate one-step
SO(10)->SM threshold is a regression, not a change of physical branch.
The lower staged finite matching is not replaced by parent-weighted logs.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
CACHE = ROOT / "tmp/p54_full_doublet_cw"
OUT = RF / "output/p54_ps_finite_thresholds.json"
ACTIVE = ["phi10:(1,2,2)", "Sigma126:(15,2,2)",
          "Sigma126:(10-pair,3,1)", "Sigma126:(10-pair,1,3)"]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def sym(a):
    return (a + a.T) / 2


def image_basis(a, tol=1e-10):
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    return u[:, s > tol]


def range_basis(p):
    d, u = np.linalg.eigh(sym(p))
    return u[:, d > .5]


def log_matrix(h, scale):
    d, u = np.linalg.eigh(sym(h))
    if d.min() <= 0:
        raise ValueError("Cannot use a nonpositive heavy mass in a finite logarithm")
    return (u * np.log(np.sqrt(d) / scale)) @ u.T


def indices(ops, p):
    return np.array([np.trace(op @ p) for op in ops], dtype=float)


def heavy_valley_api(h, active_basis, heavy_basis, upper_goldstone_projector):
    """Quadratic tree EFT data in the upper unitary-gauge slice.

    Inputs are canonical full-field tensors. The returned matrices are an
    API for a future lower matcher, not one-loop Wilson coefficients. The
    active coordinates retain their original parent normalization; therefore
    G0 is generally not the identity. Full momentum dependence is needed
    beyond the derivative expansion K(p)=S+p^2 G+O(p^4/C).
    """
    embedding = (np.eye(h.shape[0]) - upper_goldstone_projector) @ active_basis
    heavy = sym(heavy_basis.T @ h @ heavy_basis)
    source = heavy_basis.T @ h @ embedding
    response = -np.linalg.solve(heavy, source)
    g0 = sym(embedding.T @ embedding)
    schur = sym(embedding.T @ h @ embedding + source.T @ response)
    metric = sym(g0 + response.T @ response)
    tangent = embedding + heavy_basis @ response
    return {"embedding":embedding, "heavy_basis":heavy_basis, "C":heavy,
            "linear_heavy_source":source, "heavy_response_jacobian":response,
            "G0":g0, "S":schur, "G":metric, "valley_tangent":tangent}


def finite_lambda(h, heavy_scalar_basis, mv2, heavy_vector_basis, scalar_ops, vector_ops, scale):
    scalar_log = heavy_scalar_basis @ log_matrix(heavy_scalar_basis.T @ h @ heavy_scalar_basis, scale) @ heavy_scalar_basis.T
    vector_log = heavy_vector_basis @ log_matrix(heavy_vector_basis.T @ mv2 @ heavy_vector_basis, scale) @ heavy_vector_basis.T
    ps = heavy_scalar_basis @ heavy_scalar_basis.T
    pv = heavy_vector_basis @ heavy_vector_basis.T
    constant = indices(vector_ops, pv)
    scalar = np.array([np.trace(op @ scalar_log) for op in scalar_ops])
    vector_log_term = -21 * np.array([np.trace(op @ vector_log) for op in vector_ops])
    return {"scalar_log": scalar.tolist(), "vector_log": vector_log_term.tolist(),
            "vector_finite_constant": constant.tolist(),
            "lambda_total": (scalar + vector_log_term + constant).tolist(),
            "scalar_real_indices": indices(scalar_ops, ps).tolist(),
            "vector_real_indices": constant.tolist(),
            "dlambda_dlogmu": (21 * constant - indices(scalar_ops, ps)).tolist()}


def run():
    start = time.time()
    p1path = RF / "code/verify_p54_p1_hessian_spectrum.py"
    fullpath = RF / "output/p54_full_doublet_cw.json"
    censuspath = RF / "output/p54_ps_active_census.json"
    helperpath = RF / "code/verify_p54_p2_two_site_matching.py"
    p1 = module("p54_finite_p1", p1path)
    helper = module("p54_finite_helpers", helperpath)
    full = json.loads(fullpath.read_text())
    census = json.loads(censuspath.read_text())
    pars = full["tree_parameters"]
    v = full["vacuum"]
    gu = full["scheme"]["mu_over_omega"]
    x0 = p1.vacuum_vector(v["omega"], v["sigma"], v["vs"])
    potential = p1.potential_factory(pars)
    hf = p1.jax.jit(p1.hessian(potential))
    gf = p1.jax.jit(p1.grad(potential))
    keybase = hashlib.sha256(p1path.read_bytes() + json.dumps(pars, sort_keys=True).encode()).digest()
    CACHE.mkdir(parents=True, exist_ok=True)
    counts = {"cache_hits": 0, "evaluated": 0}

    def hessian(x, name):
        tick = time.time()
        key = hashlib.sha256(keybase + np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = CACHE / (key + ".npz")
        if path.exists():
            with np.load(path) as f:
                h = f["h"]
            counts["cache_hits"] += 1
        else:
            h = sym(np.asarray(hf(p1.anp.asarray(x)), float))
            np.savez_compressed(path, h=h)
            counts["evaluated"] += 1
        print(f"{name}: {time.time() - tick:.2f}s", flush=True)
        return h

    h0 = hessian(x0, "historical broken Hessian")
    psgen = helper.ps_generators(p1)
    casimir = helper.casimir_operators(p1, psgen)
    parents = helper.parent_projectors(p1, casimir)
    pa = sum(row["projector"] for row in parents if row["label"] in ACTIVE)
    ba = range_basis(pa)
    psops = helper.ps_index_operators(casimir)
    smops = helper.sm_index_operators(p1)
    gauge = [helper.so_generator(a, b) for a in range(10) for b in range(a + 1, 10)]
    high_columns = [i for i, (a, b) in enumerate(( (a, b) for a in range(10) for b in range(a + 1, 10))) if a < 6 <= b]
    pu_adj = np.zeros((45, 45)); pu_adj[high_columns, high_columns] = 1
    psadj = helper.adjoint_index_operators(psgen, gauge, {"SU4": 15, "SU2L": 3, "SU2R": 3})
    smgen = p1.sm_generators()
    smadj = helper.adjoint_index_operators(smgen, gauge, {"SU3": 8, "SU2": 3, "Y": 1}, hypercharge=True)

    # A PS background must solve its heavy radial equations with the same
    # action parameters. Re-deriving quadratic masses would change the model.
    def radial_eq(logr):
        w, s = np.exp(logr)
        grad = p1.radial_gradient(np.array([w, 0., s]), pars)
        return grad[[0, 2]]
    solved = root(radial_eq, np.log([v["omega"], v["vs"]]), tol=1e-11)
    wp, sp = np.exp(solved.x)
    if np.linalg.norm(solved.fun) > 1e-9:
        raise ValueError("No nearby same-action PS radial stationary point was recovered")
    xp = p1.vacuum_vector(float(wp), 0., float(sp))
    hp = hessian(xp, "same-action PS stationary Hessian")
    gradp = np.asarray(gf(p1.anp.asarray(xp)), float)
    orbitp = p1.gauge_orbit(xp)
    bgp = image_basis(orbitp)
    pgp = bgp @ bgp.T
    pq = p1.pq_direction(xp)
    pq -= pgp @ pq
    pq /= np.linalg.norm(pq)
    ppq = np.outer(pq, pq)
    ph = sym(np.eye(p1.N_REAL) - pa - pgp - ppq)
    bh = range_basis(ph)
    heavy = sym(bh.T @ hp @ bh)
    active = sym(ba.T @ hp @ ba)
    mixed = ba.T @ hp @ bh
    active_values = np.linalg.eigvalsh(active)
    heavy_values = np.linalg.eigvalsh(heavy)
    heavy_positive = bool(heavy_values.min() > 1e-8)
    mvps = gu**2 * orbitp.T @ orbitp
    bvp = image_basis(mvps)
    mu_upper = gu * wp
    upper_finite = finite_lambda(hp, bh, mvps, bvp, psops, psadj, mu_upper) if heavy_positive else None
    upper_scale_rows = []
    if upper_finite is not None:
        for factor in (.5, 1., 2.):
            row = finite_lambda(hp, bh, mvps, bvp, psops, psadj, factor*mu_upper)
            upper_scale_rows.append({"mu_over_reference":factor, "lambda_4_L_R":row["lambda_total"]})
    schur = sym(active - mixed @ np.linalg.solve(heavy, mixed.T)) if heavy_positive else None
    parent_spectrum = []
    for row in parents:
        br = range_basis(row["projector"])
        dr = np.linalg.eigvalsh(sym(br.T @ hp @ br))
        parent_spectrum.append({"label": row["label"], "active": row["label"] in ACTIVE,
                                "min_m2": float(dr.min()), "max_m2": float(dr.max()),
                                "negative_count": int(np.sum(dr < -1e-8)), "dimension": len(dr)})

    # At the simultaneous broken vacuum the U/I orbit projectors are
    # constructed from the very same generator sites as the vector spectrum.
    orbit0 = p1.gauge_orbit(x0)
    bgu = image_basis(orbit0 @ pu_adj)
    pgu = bgu @ bgu.T
    bgi = image_basis(orbit0 @ (np.eye(45) - pu_adj))
    pgi = bgi @ bgi.T
    pg = pgu + pgi
    bg = image_basis(orbit0)
    # The 55 upper heavy scalar directions stay SM-covariant and orthogonal
    # to the actual gauge orbit. The active parent coordinates must first be
    # pulled into the upper unitary-gauge slice, which gives a nontrivial G0.
    valley = heavy_valley_api(h0, ba, bh, pgu)
    valley_values = eigh(valley["S"], valley["G"], eigvals_only=True)
    pq0 = p1.pq_direction(x0)
    pq0 -= pg @ pq0
    pq0 /= np.linalg.norm(pq0)
    image = valley["embedding"]
    image_projector = image @ np.linalg.solve(valley["G0"], image.T)
    valley_partition_error = np.linalg.norm(image_projector + bh@bh.T + pgu + np.outer(pq0,pq0) - np.eye(p1.N_REAL))
    valley_heavy_min = float(np.linalg.eigvalsh(valley["C"]).min())
    valley_metric_min = float(np.linalg.eigvalsh(valley["G"]).min())
    lower_orbit_coordinates = ba.T @ bgi
    lower_ward_error = float(np.linalg.norm(valley["S"]@lower_orbit_coordinates))
    sm_reps = p1.sm_representation_matrices()
    heavy_sm_comm = max(np.linalg.norm((bh@bh.T)@t-t@(bh@bh.T)) for group in sm_reps.values() for t in group)
    valley_metric_comm = max(np.linalg.norm(valley["G"]@(ba.T@t@ba)-(ba.T@t@ba)@valley["G"]) for group in sm_reps.values() for t in group)
    valley_mass_comm = max(np.linalg.norm(valley["S"]@(ba.T@t@ba)-(ba.T@t@ba)@valley["S"]) for group in sm_reps.values() for t in group)
    d0, u0 = np.linalg.eigh(h0)
    bphys = u0[:, d0 > 1e-7]
    mv0 = gu**2 * orbit0.T @ orbit0
    bv0 = image_basis(mv0)
    mu_full = gu * v["omega"]
    one_step = finite_lambda(h0, bphys, mv0, bv0, smops, smadj, mu_full)
    scale_rows = []
    for factor in (.5, 1., 2.):
        r = finite_lambda(h0, bphys, mv0, bv0, smops, smadj, factor * mu_full)
        scale_rows.append({"mu_over_reference": factor, "lambda_SM_3_2_1": r["lambda_total"]})
    fd = (np.asarray(scale_rows[2]["lambda_SM_3_2_1"]) - np.asarray(scale_rows[0]["lambda_SM_3_2_1"])) / math.log(4)

    repair = census["minimal_stage_consistent_repair"]
    get = lambda key: np.array([a["value"] for a in repair[key]])
    vi, vu = get("heavy_vector_real_indices_MI_SM"), get("heavy_vector_real_indices_MU_PS")
    si, su = get("physical_heavy_scalar_real_indices_MI_SM"), get("physical_heavy_scalar_real_indices_MU_PS")
    a_ps = np.array([r["value"] for r in census["reconstructed_a_PS"]])
    a_sm = np.array([-7., -19/6, 41/10])
    a10 = -34/3
    project = lambda a: np.array([a[0], a[1], .4*a[0] + .6*a[2]])
    slope_i, slope_u = 21*vi-si, 21*vu-su
    required_i, required_u = -6*(project(a_ps)-a_sm), -6*(a10-a_ps)
    log_checks = []
    for stage, observed, required, labels in (("I", slope_i, required_i, ["3","2","1"]), ("U", slope_u, required_u, ["4","L","R"])):
        for label, got, want in zip(labels, observed, required):
            log_checks.append({"stage": stage, "group": label, "derived_slope": float(got),
                               "beta_required_slope": float(want), "residual": float(got-want)})
    actual_active_comm = float(np.linalg.norm(pa @ h0 - h0 @ pa))
    upper_active_comm = float(np.linalg.norm(pa @ hp - hp @ pa))
    goldstone_comm = float(np.linalg.norm(pa @ pgu - pgu @ pa))
    upper_goldstone_comm = float(np.linalg.norm(pa @ pgp - pgp @ pa))
    # An explicit equal-charge scalar example shows the issue is more than
    # basis dependence: parent-weighted log is not the EFT Schur determinant.
    toy = np.array([[2., .5], [.5, 5.]])
    toy_log = 2*log_matrix(toy, 1.)
    toy_s0 = 2.-.5**2/5.
    toy_p2 = .7
    toy_sp = toy_p2+2.-.5**2/(toy_p2+5.)
    toy_det_error = abs(np.linalg.slogdet(toy_p2*np.eye(2)+toy)[1] - math.log(toy_p2+5.) - math.log(toy_sp))
    toy_weighted_error = float(toy_log[0,0]-math.log(toy_s0))
    checks = [
        ("four-parent active projector is orthogonal and has rank 248", np.linalg.norm(pa@pa-pa)<1e-10 and ba.shape[1]==248),
        ("same-action upper PS background is fully stationary", np.linalg.norm(gradp)<1e-8),
        ("upper PS orbit has 24 gauge Goldstones", bgp.shape[1]==24),
        ("upper PS orbit obeys stationary Ward identity", np.linalg.norm(hp@bgp)<1e-9),
        ("upper retained plus integrated plus gauge/PQ planes partition fields", np.linalg.norm(ph@ph-ph)<1e-9 and bh.shape[1]+ba.shape[1]+bgp.shape[1]+1==328),
        ("PS-preserving Hessian commutes with active projector", upper_active_comm<1e-9),
        ("upper heavy block positive despite retained tachyons", heavy_positive),
        ("upper Goldstone projector commutes with active projector", upper_goldstone_comm<1e-9),
        ("broken-vacuum site orbit ranks are 24 and 9", bgu.shape[1]==24 and bgi.shape[1]==9),
        ("site orbit projectors are orthogonal and sum to full orbit", np.linalg.norm(pgu@pgi)<1e-9 and np.linalg.norm(pg-bg@bg.T)<1e-9),
        ("actual 55-dimensional heavy valley is positive and gauge-orthogonal", valley_heavy_min>1e-8 and np.linalg.norm(bh.T@bg)<1e-9),
        ("upper-gauge-slice active plus heavy and symmetry planes are complete", valley_partition_error<1e-9),
        ("heavy valley pullback metric is positive", valley_metric_min>1e-8),
        ("heavy valley kinetic and mass tensors are SM covariant", max(heavy_sm_comm,valley_metric_comm,valley_mass_comm)<1e-8),
        ("heavy valley preserves nine lower gauge plus four Higgs zeros", int(np.sum(abs(valley_values)<1e-7))==13 and int(np.sum(valley_values< -1e-7))==0 and lower_ward_error<1e-9),
        ("momentum Schur determinant identity holds while weighted-log shortcut fails", toy_det_error<1e-14 and abs(toy_weighted_error)>1e-4),
        ("full broken scalar physical hard rank is 290", bphys.shape[1]==290),
        ("one-step unbroken generators commute with scalar mass", max(np.linalg.norm(op@h0-h0@op) for op in smops)<1e-9),
        ("one-step SM indices are physical scalars 79,77,377/5", np.linalg.norm(np.asarray(one_step["scalar_real_indices"])-[79,77,377/5])<1e-8),
        ("one-step finite vector constants are 5,6,8", np.linalg.norm(np.asarray(one_step["vector_finite_constant"])-[5,6,8])<1e-9),
        ("one-step finite threshold scale derivative agrees with beta difference", np.linalg.norm(fd+6*(a10-a_sm))<1e-8),
        ("all six staged logarithmic identities close", max(abs(row["residual"]) for row in log_checks)<1e-12),
        ("staged slopes compose to one-step slope", np.linalg.norm(project(slope_u)+slope_i-fd)<1e-8),
    ]
    if upper_finite is not None:
        checks.append(("actual upper finite scalar/vector indices realize repaired upper slopes", np.linalg.norm(np.asarray(upper_finite["dlambda_dlogmu"])-required_u)<1e-8))
        upper_fd = (np.asarray(upper_scale_rows[2]["lambda_4_L_R"])-np.asarray(upper_scale_rows[0]["lambda_4_L_R"]))/math.log(4)
        checks.append(("upper finite threshold scale variation reproduces all three slopes", np.linalg.norm(upper_fd-required_u)<1e-8))
    checks = [{"name": name, "pass": bool(ok)} for name, ok in checks]
    sources = [p1path, helperpath, fullpath, censuspath, Path(__file__).resolve()]
    return {
        "schema": "p54-ps-finite-thresholds-v1", "date": "2026-09-05",
        "status": {"upper_PS_finite_matching_computed": upper_finite is not None,
                   "lower_PS_to_SM_finite_matching_computed": False, "six_stage_log_identities_closed": True,
                   "one_step_SM_finite_regression_computed": True, "physical_branch_changed": False,
                   "scales_refitted": False, "whole_model_excluded": False},
        "scheme": {
            "matching": "alpha_low^-1=alpha_high^-1-lambda/(12pi), MSbar one-loop zero-momentum gauge matching",
            "finite_formula": "lambda_i=Tr_V I_i -21 Tr_V[I_i log(MV/mu)] +Tr_physicalS[I_i log(MS/mu)]",
            "Goldstone": "one real Goldstone per massive vector removed from scalar trace; vector coefficient -21 includes vector/Goldstone/ghost package",
            "quadratic_mass_scheme": "same declared tree masses; finite fixed-VEV invariant CTs inherited from full-doublet ledger, not rederived at sigma=0",
            "upper_stationary_saddle_scheme": "tree heavy-radial stationary saddle at sigma=0 in the same frozen renormalized parameters; not the one-loop fixed-VEV retuned action's stationary point",
            "background_shift_order": "a one-loop shift of this upper saddle changes one-loop gauge mass logarithms only at two-loop order; the finite tree displacement from the nonzero-sigma vacuum is retained exactly",
            "nu_CT": full["tadpoles"]["finite_mass_counterterms_over_omega2"]["delta_nu2"],
            "nu_CT_normalization": "deltaH=-delta_nu2 P126/2; inserting this one-loop CT into one-loop gauge mass logs would add selected two-loop terms and is not done",
            "neutral_fermions": "heavy RH neutrinos are unbroken-SM singlets; no gauge threshold at this order",
        },
        "active_parent_labels": ACTIVE,
        "upper_PS": {"vacuum": {"omega":float(wp), "sigma":0., "vs":float(sp)},
            "same_action_parameters": pars, "gradient_norm":float(np.linalg.norm(gradp)),
            "gauge_rank":bgp.shape[1], "active_real_dimension":ba.shape[1], "physical_integrated_dimension":bh.shape[1],
            "active_mass_eigenvalue_min":float(active_values.min()), "active_mass_eigenvalue_max":float(active_values.max()),
            "active_tachyon_count":int(np.sum(active_values< -1e-8)),
            "heavy_mass_eigenvalue_min":float(heavy_values.min()), "heavy_mass_eigenvalue_max":float(heavy_values.max()),
            "heavy_mass_eigenvalues":heavy_values.tolist(),
            "largest_retained_to_smallest_integrated_abs_m2_ratio":float(max(abs(active_values))/heavy_values.min()),
            "heavy_tachyon_count":int(np.sum(heavy_values< -1e-8)),
            "active_heavy_mass_mixing_norm":float(np.linalg.norm(mixed)),
            "Schur_mass_difference_from_active_norm":float(np.linalg.norm(schur-active)) if schur is not None else None,
            "parent_spectrum":parent_spectrum, "finite_threshold":upper_finite, "finite_scale_variation":upper_scale_rows,
            "mu_over_omega_reference":float(mu_upper)},
        "simultaneous_background_obstruction": {
            "active_mass_commutator_norm":actual_active_comm, "upper_active_mass_commutator_norm":upper_active_comm,
            "active_upperGoldstone_commutator_norm":goldstone_comm,
            "active_upperGoldstone_overlap_trace":float(np.trace(pa@pgu)),
            "site_Goldstone_ranks":[bgu.shape[1],bgi.shape[1]],
            "site_Goldstone_index_MI_SM":indices(smops,pgi).tolist(),
            "site_Goldstone_index_MU_SM":indices(smops,pgu).tolist(),
            "required_lower_EFT": "derive the PS-covariant functional with all heavy fields eliminated and transport it to the broken background; its scalar Schur operator is momentum-dependent",
            "Schur_operator": "in orthonormal block coordinates Keff(p)=p^2+A-B(p^2+C)^-1 Bdagger, so Keff(0)=A-BC^-1 Bdagger and Zeff=I+BC^-2 Bdagger",
            "why_weighted_log_is_insufficient": "Tr(Pactive I log H) is not log det of the active EFT; replacing it by log of zero-momentum Schur mass also omits induced kinetic/gauge vertices",
            "upper_background_limit": "upper matching is the dimension-four PS coefficient at sigma=0; finite sigma generates higher-dimensional gauge operators and must not be silently equated to the broken weighted-log partition",
            "equal_charge_2x2_counterexample": {"mass_squared_matrix":toy.tolist(), "parent_weighted_log":float(toy_log[0,0]),
                "log_zero_momentum_Schur_mass":math.log(toy_s0), "difference":toy_weighted_error,
                "induced_kinetic_coefficient":1.+.5**2/5.**2,"full_momentum_determinant_residual":toy_det_error},
        },
        "actual_vacuum_tree_heavy_valley": {
            "scope":"same-action quadratic tree heavy-valley and derivative-expansion diagnostic, not finite lower one-loop matching",
            "coordinate_formula":"B=(I-PGupper)Bactive; W=BH^T H B; Keff(p^2)=p^2 G0+B^T H B-W^T(p^2 I+C)^-1 W; G0=B^T B",
            "heavy_real_dimension":bh.shape[1], "active_real_dimension":ba.shape[1],
            "heavy_min_m2":valley_heavy_min,"heavy_max_m2":float(np.linalg.eigvalsh(valley["C"]).max()),
            "source_linear_norm":float(np.linalg.norm(valley["linear_heavy_source"])),
            "heavy_response_Jacobian_norm":float(np.linalg.norm(valley["heavy_response_jacobian"],2)),
            "slice_metric_G0_eigenvalue_min":float(np.linalg.eigvalsh(valley["G0"]).min()),
            "pullback_metric_eigenvalue_min":valley_metric_min,
            "pullback_metric_eigenvalue_max":float(np.linalg.eigvalsh(valley["G"]).max()),
            "Schur_generalized_eigenvalue_groups":p1.group_eigenvalues(valley_values,abs_tol=1e-7),
            "Schur_zero_count":int(np.sum(abs(valley_values)<1e-7)),"Schur_negative_count":int(np.sum(valley_values< -1e-7)),
            "lower_Goldstone_Ward_residual":lower_ward_error,"slice_partition_residual":float(valley_partition_error),
            "SM_covariance_residuals":{"heavy_projector":float(heavy_sm_comm),"Schur_mass":float(valley_mass_comm),"pullback_metric":float(valley_metric_comm)},
            "public_source_API":"heavy_valley_api(H0, Bactive, Bheavy_upper, PGupper_actual) returns embedding,C,linear_heavy_source,heavy_response_jacobian,G0,S,G,valley_tangent",
            "basis_construction":"Bactive=range(sum of four exact PS tensor parent projectors); Bheavy_upper=range(I-Pactive-PGupper_PS-PQupper); PGupper_actual=image(O(actual) Padjoint_SO10/PS)",
            "missing_for_lower_finite_matching":"upper-matched scalar interactions and gauge higher-dimensional operators, full induced covariant kinetic vertices, lower vector/Goldstone/ghost functional, and consistent matching-scale transport; generalized(S,G) masses are low-momentum approximations, not full poles for threshold logs",
        },
        "one_step_SM_regression": {"scope":"same broken tree spectrum; complete one-step matching only, not a replacement of the selected two-stage branch",
            "mu_over_omega":float(mu_full), "finite_threshold":one_step, "scale_variation":scale_rows,
            "measured_scale_derivative":fd.tolist(), "required_scale_derivative":(-6*(a10-a_sm)).tolist()},
        "six_stage_log_identities":log_checks,
        "checks":checks, "summary":{"passed":sum(x["pass"] for x in checks),"total":len(checks),"all_pass":all(x["pass"] for x in checks)},
        "cache":counts,"runtime_seconds":time.time()-start,
        "sources":[{"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
        "primary_threshold_references": [{"title":"Threshold effects in SO(10) models with one intermediate breaking scale", "url":"https://doi.org/10.1140/epjc/s10052-020-8308-9"},
            {"title":"Grand Unification of Effective Gauge Theories", "url":"https://doi.org/10.1016/0550-3213(81)90498-3"}],
    }


def markdown(r):
    up=r["upper_PS"]; full=r["one_step_SM_regression"]
    lines=["# PS finite gauge-threshold audit", "",
        "The four-parent active EFT is unchanged. The upper PS-symmetric coefficient and a full one-step SM regression are computed; the finite lower staged matching is still open.", "",
        f"Checks: **{r['summary']['passed']}/{r['summary']['total']}**.", "",
        "## Upper same-action PS background", "",
        f"- VEVs: `{up['vacuum']}`; full gradient norm `{up['gradient_norm']:.3e}`.",
        f"- Retained active fields: `{up['active_real_dimension']}` real, `{up['active_tachyon_count']}` tachyonic modes.",
        f"- Integrated physical fields: `{up['physical_integrated_dimension']}` real; smallest mass squared `{up['heavy_mass_eigenvalue_min']:.12g}`.",
        "- Active tachyons describe the lower symmetry breaking and are not included in upper hard logarithms.",
        "- This is a tree stationary saddle of the frozen parameters. It is not claimed to be stationary for the one-loop fixed-VEV retuned action.",
        f"- Largest retained / smallest integrated squared-mass ratio: `{up['largest_retained_to_smallest_integrated_abs_m2_ratio']:.6g}`. There is no uniform hierarchy between every retained and integrated scalar."]
    if up["finite_threshold"]:
        lines += [f"- Finite upper lambda in `(4,L,R)`: `{up['finite_threshold']['lambda_total']}`.",
                  f"- Upper scale slopes: `{up['finite_threshold']['dlambda_dlogmu']}`."]
    ob=r["simultaneous_background_obstruction"]
    lines += ["", "## Lower matching remains a functional calculation", "",
        f"At simultaneous nonzero VEVs, `||[Pactive,H]||={ob['active_mass_commutator_norm']:.6g}` and `||[Pactive,P_G,U]||={ob['active_upperGoldstone_commutator_norm']:.6g}`.",
        "The same site-resolved gauge orbit is used for vector and Goldstone projectors. Its upper Goldstone image has active-parent support, so a naive complementary scalar projector is not a staged EFT.",
        "The exact scalar Schur operator contains momentum dependence and induced gauge vertices. A weighted mass logarithm, or only the zero-momentum Schur mass, does not compute that finite matching.", "",
        "## Actual-vacuum tree heavy valley", ""]
    valley=r["actual_vacuum_tree_heavy_valley"]
    lines += [f"- The same 55 heavy scalar directions have positive `C`, with minimum mass squared `{valley['heavy_min_m2']:.12g}`.",
        "- The upper gauge slice uses `B=(I-PGupper)Bactive`, with source `W=BH^T H B` and response `D=-C^-1 W`.",
        "- Its tree Schur mass and pullback metric are `S=B^T H B-W^T C^-1 W` and `G=B^T B+D^T D`.",
        f"- Metric eigenvalues range from `{valley['pullback_metric_eigenvalue_min']:.12g}` to `{valley['pullback_metric_eigenvalue_max']:.12g}`; `||D||2={valley['heavy_response_Jacobian_norm']:.12g}`.",
        f"- The generalized pair `(S,G)` has `{valley['Schur_zero_count']}` zeros (nine lower Goldstones and four Higgs coordinates) and no negative eigenvalues. Its lower-orbit Ward residual is `{valley['lower_Goldstone_Ward_residual']:.3e}`.",
        "- `heavy_valley_api(...)` exports all quadratic source, metric and embedding matrices. These derivative-expansion masses do not replace the finite lower functional or its gauge vertices.", "",
        "## One-step unbroken-SM regression", "",
        f"- Complete finite lambda `(3,2,1)`: `{full['finite_threshold']['lambda_total']}`.",
        f"- Finite vector constants: `{full['finite_threshold']['vector_finite_constant']}`.",
        f"- Measured scale slopes: `{full['measured_scale_derivative']}`.",
        "- All six staged logarithmic identities close and sum to this one-step slope. No physical two-stage scale refit is claimed.", "",
        "The fixed-VEV mass counterterms retain the original invariant normalization, including `deltaH_nu=-delta_nu2 P126/2`. One-loop counterterm insertion into these one-loop threshold mass logs is beyond the stated order.", "",
        "Finite-threshold convention: [Hall](https://doi.org/10.1016/0550-3213(81)90498-3), [SO(10) threshold calculation](https://doi.org/10.1140/epjc/s10052-020-8308-9).", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    report=run()
    OUT.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    OUT.with_suffix(".md").write_text(markdown(report))
    print(f"PS finite thresholds: {report['summary']['passed']}/{report['summary']['total']}")
    print("upper:",report["upper_PS"]["finite_threshold"])
    print("one-step:",report["one_step_SM_regression"]["finite_threshold"])
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)
