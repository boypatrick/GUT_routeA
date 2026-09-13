#!/usr/bin/env python3
"""Same-action, nonuniform P54 scalar candidate evaluation.

This module does not select a point, fit data, or mutate a benchmark.  Its
caller supplies exact ambient Hessian jets assembled from the declared
invariant basis.  Every candidate gets its OWN stationary mass inputs,
tree light-doublet tuning, spectrum, Coleman--Weinberg weights and scalar
kinetic matrix.  Frozen-propagator attribution is not reused as a physical
loop prediction at a changed point.

The result is a radial hard-potential/kinetic engineering audit, NOT a
full gauge/ghost/fermion momentum calculation or a physical pole test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

import verify_p54_common_renormalization as common


LOOP = 16 * np.pi**2


def bubble_slope_matrix(mass_squared, soft_count=38):
    r"""Return d/ds int_0^1 log[(1-x)a+xb+x(1-x)s] dx at s=0.

    All entries containing at least one hard state are retained.  The
    soft--soft derivative is nonlocal and is excluded, not approximated
    by zero as a physical self-energy.  For a,b>0 the exact weight is

        [(a+b)/2 - ab log(b/a)/(b-a)] / (b-a)^2.

    Equal masses give 1/(6a), and (0,b) gives 1/(2b).  The convergent
    even series in t=(b-a)/(a+b) avoids cancellation near degeneracy.
    No mass-degeneracy grouping approximation is needed.
    """
    masses = np.asarray(mass_squared, dtype=float).copy()
    if masses.ndim != 1 or not np.all(np.isfinite(masses)):
        raise ValueError("Expected a finite one-dimensional mass-squared array")
    if not 0 <= soft_count < len(masses):
        raise ValueError("soft_count must leave a nonempty hard spectrum")
    if np.any(masses[soft_count:] <= 0):
        raise ValueError("Every retained hard squared mass must be positive")
    # The caller independently verifies the exact soft cluster.  Roundoff
    # in symmetry/critical eigenvalues must not create regulator masses.
    masses[:soft_count] = 0.0
    aa, bb = np.broadcast_arrays(masses[:, None], masses[None, :])
    weights = np.zeros_like(aa)
    mixed = (aa == 0) ^ (bb == 0)
    weights[mixed] = 1.0 / (2.0 * (aa + bb)[mixed])
    hard = (aa > 0) & (bb > 0)
    a, b = aa[hard], bb[hard]
    mean = (a + b) / 2.0
    t = (b - a) / (a + b)
    near = np.abs(t) < 0.05
    value = np.empty_like(a)
    # O(t^12) with |t|<.05 has relative error below double precision.
    value[near] = sum(
        t[near] ** (2 * n) / ((2 * n + 1) * (2 * n + 3))
        for n in range(6)
    ) / (2.0 * mean[near])
    far = ~near
    d = b[far] - a[far]
    value[far] = (
        mean[far] - a[far] * b[far] * (np.log(b[far]) - np.log(a[far])) / d
    ) / d**2
    weights[hard] = value
    return (weights + weights.T) / 2.0


def verify_bubble_weights(checks):
    """Persist the independent 80-pair analytic-versus-quadrature check.

    Scaling the integration denominator by max(a,b) keeps quadrature's
    absolute tolerance meaningful over the sixteen-decade mass ratio.
    The loosest relative errors occur at extreme ratios, where adaptive
    quadrature resolves a narrow endpoint less accurately; this does not
    replace the exact equal-mass and mixed hard--soft limit checks.
    """
    masses = np.array([0., 1e-8, .001, 1., 1.000000000001,
                       1.099, 1.101, 4., 1e8])
    weights = bubble_slope_matrix(masses, soft_count=1)
    rows = []
    for i, a in enumerate(masses):
        for j, b in enumerate(masses):
            if i + j == 0:
                continue
            scale = max(a, b)
            integral, quadrature_error = quad(
                lambda x: x * (1 - x) / ((1 - x) * (a / scale) + x * (b / scale)),
                0., 1., epsabs=1e-13, epsrel=1e-12,
            )
            numerical = integral / scale
            error = float(abs(weights[i, j] - numerical) / abs(numerical))
            tolerance = 2e-9
            check = {
                "name": f"analytic_scalar_bubble_slope_pair_{i}_{j}",
                "residual": error, "tolerance": tolerance,
                "pass": bool(np.isfinite(error) and error < tolerance),
            }
            checks.append(check)
            rows.append({
                "squared_masses": [float(a), float(b)],
                "analytic_slope": float(weights[i, j]),
                "quadrature_slope": float(numerical),
                "quadrature_estimated_absolute_error": float(quadrature_error / scale),
                "relative_error": error,
            })
    symmetry_error = float(np.linalg.norm(weights - weights.T))
    checks.append({"name": "analytic_scalar_bubble_slope_symmetric",
                   "residual": symmetry_error, "tolerance": 1e-14,
                   "pass": symmetry_error < 1e-14})
    negative_part = float(max(0., -weights.min()))
    checks.append({"name": "analytic_scalar_bubble_slope_nonnegative",
                   "residual": negative_part, "tolerance": 1e-14,
                   "pass": negative_part < 1e-14})
    equal_error = float(np.max(abs(weights.diagonal()[1:] * (6 * masses[1:]) - 1)))
    checks.append({"name": "analytic_scalar_bubble_slope_equal_mass_limit",
                   "residual": equal_error, "tolerance": 1e-14,
                   "pass": equal_error < 1e-14})
    mixed_error = float(np.max(abs(weights[0, 1:] * (2 * masses[1:]) - 1)))
    checks.append({"name": "analytic_scalar_bubble_slope_massless_limit",
                   "residual": mixed_error, "tolerance": 1e-14,
                   "pass": mixed_error < 1e-14})
    return {
        "quadrature_pair_count": len(rows),
        "maximum_relative_error": max(row["relative_error"] for row in rows),
        "pair_checks": rows,
        "soft_soft_derivative_excluded_as_nonlocal": True,
    }


def evaluate_candidate(
    p1, cw, search, pars, r, h, t, second, mu, rad, checks, label,
    *, engineering_margin=0.3, gauge_coupling=None,
):
    r"""Reanchor, retune and fully reweight one independently chosen point.

    Parameters ``pars`` and h=V_,AB must describe the SAME untuned action
    at r=(w,sigma,vs).  t_r=partial_r h and second_rs=partial_r partial_s h
    are fixed-Lagrangian-parameter derivatives, not derivatives along a
    stationary-mass curve.  Their shapes are (3,N,N) and (3,3,N,N).

    The exact affine dependent-mass update is
      dp = -A(r)^(-1) grad_r V,  A=-diag(12w/5,sigma,vs),
      h -> h + D(dp).
    The existing -xi02 P10 Hessian is then removed and the first physical
    crossing of h-xi02 P10 is found.  Both updates have zero radial jets,
    so t and second remain unchanged.  One-loop tadpoles use a SEPARATE
    finite common mass CT, which is not folded into the tree propagators.

    The historical benchmark has g_U=mu/omega.  ``gauge_coupling`` lets a
    caller distinguish them; omission explicitly retains that convention.
    ``checks`` is the caller's list of numerical identity checks.
    """
    pars = dict(pars)
    r = np.asarray(r, dtype=float)
    h = np.array(h, dtype=float, copy=True)
    t = np.asarray(t, dtype=float)
    second = np.asarray(second, dtype=float)
    rad = np.asarray(rad, dtype=float)
    n = p1.N_REAL
    if r.shape != (3,) or np.any(r <= 0):
        raise ValueError("This evaluator requires three strictly positive radii")
    if h.shape != (n, n) or t.shape != (3, n, n) or second.shape != (3, 3, n, n):
        raise ValueError("Incorrect exact Hessian-jet shapes")
    if rad.shape != (n, 3):
        raise ValueError("Incorrect radial embedding shape")
    if not all(np.all(np.isfinite(a)) for a in (r, h, t, second, rad)):
        raise ValueError("Nonfinite candidate input")
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError("The renormalization scale must be positive")
    gauge = float(mu if gauge_coupling is None else gauge_coupling)
    if not np.isfinite(gauge) or gauge <= 0:
        raise ValueError("The gauge coupling must be positive")
    if not 0 < engineering_margin < 1:
        raise ValueError("An explicit engineering margin in (0,1) is required")
    check_start = len(checks)

    def check(name, actual, expected=0.0, tolerance=2e-8):
        actual, expected = np.asarray(actual), np.asarray(expected)
        err = float(np.linalg.norm(actual - expected) /
                    max(1.0, np.linalg.norm(actual), np.linalg.norm(expected)))
        passed = bool(np.isfinite(err) and err < tolerance)
        checks.append({"name": str(label) + "_" + name,
                       "residual": err, "tolerance": tolerance, "pass": passed})
        return passed

    metric = rad.T @ rad
    x = rad @ r
    initial_parameters = dict(pars)
    result = {
        "label": str(label), "radii": r.tolist(),
        "input_tree_parameters": initial_parameters,
        "mu_over_omega": float(mu), "gauge_coupling": gauge,
        "gauge_equals_mu_by_default": gauge_coupling is None,
        "engineering_margin": float(engineering_margin),
        "radial_metric": metric.tolist(),
        "accepted_for_next_radial_audit": False,
        "is_new_physical_benchmark": False,
        "full_background_stability_decided": False,
        "full_gauge_ghost_fermion_momentum_kernel": False,
        "one_loop_light_doublet_retuned": False,
        "physical_fit_enabled": False,
        "default_matching_inputs_mutated": False,
        "missing": [
            "bounded-from-below quartic-action audit and global-vacuum comparison; radial Hessian jets cannot constrain pure-H quartics xi1/xi2",
            "complete gauge/ghost/fermion momentum terms and gauge/Nielsen consistency",
            "nonradial one-loop stability, momentum-dependent poles and loop light-doublet tuning",
            "new upper PS saddle, running and threshold verification at this point",
            "complete Wilson/box and lower-end finite gauge/ghost/Yukawa matching",
            "complete finite seesaw including C_HN and mass-by-mass running",
        ],
    }

    def reject(reason):
        result["failure_reason"] = reason
        result["numerical_identity_checks_pass"] = all(c["pass"] for c in checks[check_start:])
        result["tree_parameters"] = dict(pars)
        return result

    # On the pure complex-10 ray, q=|H^T H|^2/(H^dagger H)^2 spans [0,1].
    # V4/(H^dagger H)^2=xi1+xi2*q, so BOTH endpoint inequalities are
    # necessary.  At the H=0 radial background these quartics contribute
    # identically zero to h,t,second: a radial pass cannot detect this
    # exact unbounded direction.  This is a candidate gate, not a failed
    # algebraic identity.  Flat quartic endpoints need lower-degree tests.
    xi1, xi2 = float(pars["xi1"]), float(pars["xi2"])
    pure_h_bfb = bool(xi1 >= 0. and xi1 + xi2 >= 0.)
    result["pure_H_quartic_necessary_BFB"] = {
        "xi1": xi1, "xi1_plus_xi2": xi1 + xi2,
        "minimum_normalized_quartic": min(xi1, xi1 + xi2),
        "conditions": ["xi1 >= 0", "xi1 + xi2 >= 0"],
        "pass": pure_h_bfb,
        "sufficient_for_full_action_boundedness": False,
        "equality_requires_lower_degree_flat_direction_audit": bool(xi1 == 0. or xi1 + xi2 == 0.),
        "full_mixed_field_BFB_completed": False,
    }
    if not pure_h_bfb:
        return reject("A necessary pure-H quartic BFB condition fails, despite radial-Hessian invisibility")

    check("ambient_Hessian_symmetric", h, h.T)
    check("radial_first_jets_symmetric", t, t.swapaxes(-1, -2))
    check("radial_second_jets_symmetric", second, second.swapaxes(-1, -2))
    check("radial_second_jets_Schwarz", second, second.swapaxes(0, 1))
    check("canonical_radial_metric", metric, np.diag([12 / 5, 2., 1.]))
    if not all(c["pass"] for c in checks[check_start:]):
        return reject("Input exact-action jet identities failed")

    dmass = -np.linalg.solve(common.radial_mass_map(r), p1.radial_gradient(r, pars))
    doperator = common.mass_operator(dmass)
    check("mass_update_normalization", rad.T @ doperator @ x,
          common.radial_mass_map(r) @ dmass)
    for key, value in zip(common.MASS_KEYS, dmass):
        pars[key] += float(value)
    h += doperator
    result["stationarity_mass_reanchor"] = dmass.tolist()

    p10 = np.zeros((n, n))
    p10[p1.SL_H_RE, p1.SL_H_RE] = np.eye(p1.SL_H_RE.stop - p1.SL_H_RE.start)
    p10[p1.SL_H_IM, p1.SL_H_IM] = np.eye(p1.SL_H_IM.stop - p1.SL_H_IM.start)
    h += float(pars["xi02"]) * p10
    pars["xi02"] = 0.0
    stationary = check("reanchored_tree_stationary", p1.radial_gradient(r, pars))
    symmetry, symmetry_basis, complement = search.symmetry_complement(p1, x)
    sym_rank = symmetry_basis.shape[1]
    result["tree_symmetry_rank"] = int(sym_rank)
    ward = check("tree_Goldstone_Ward_identity", h @ symmetry_basis)
    if not stationary or not ward or sym_rank != 34:
        return reject("Stationarity or the expected rank-34 symmetry orbit failed")
    tuning = search.tune_first_instability(h, p10, complement)
    result["tree_tuning"] = tuning
    if not tuning.get("tunable"):
        return reject("No admissible first tree light-mode crossing")
    pars["xi02"] = float(tuning["xi02"])
    h -= pars["xi02"] * p10
    lam, u = np.linalg.eigh((h + h.T) / 2)
    tolerance = float(tuning["zero_tolerance"])
    zero_count = int(np.sum(np.abs(lam) < tolerance))
    negative_count = int(np.sum(lam < -tolerance))
    result.update({"tree_scalar_eigenvalues": lam.tolist(),
                   "tree_zero_count": zero_count,
                   "tree_negative_count": negative_count,
                   "tree_parameters": dict(pars)})
    if zero_count != 38 or negative_count != 0 or lam[38] <= tolerance:
        return reject("Tree spectrum is not exactly 34 symmetry plus four critical modes with a positive hard gap")
    classification = p1.classify_extra_zero_modes(lam, u, tolerance, symmetry)
    result["extra_zero_mode_classification"] = classification
    first_is_doublet = classification["rank"] == 4
    for group, expected in (("SU3", 0.0), ("SU2", 0.75), ("Y", 0.25)):
        values = np.asarray(classification["casimir_eigenvalues"][group])
        first_is_doublet &= values.shape == (4,) and bool(np.max(abs(values - expected)) < 3e-6)
        first_is_doublet &= classification["invariance_leakage"][group] < 3e-6
    result["tree_first_crossing_is_one_complex_doublet"] = bool(first_is_doublet)
    if not first_is_doublet:
        return reject("The first crossing is not one invariant SM complex doublet")
    result["tree_min_hard_mass2"] = float(lam[38])

    gens = cw.generators()
    orbit = cw.field_orbit(p1, x, gens)
    orad = [cw.field_orbit(p1, q, gens) for q in rad.T]
    mv = gauge**2 * orbit.T @ orbit
    vl, vu = np.linalg.eigh(mv)
    vector_tolerance = 1e-10 * max(1., float(vl[-1]))
    result["tree_vector_eigenvalues"] = vl.tolist()
    if np.sum(vl > vector_tolerance) != 33 or np.min(vl) < -vector_tolerance:
        return reject("The vector spectrum does not have the expected 33 hard and 12 massless modes")
    vt = [gauge**2 * (oi.T @ orbit + orbit.T @ oi) for oi in orad]
    mu2 = float(mu)**2
    te = np.array([u.T @ v @ u for v in t])
    ts = np.array([
        np.dot(lam[38:] * (np.log(lam[38:] / mu2) - 1), np.diag(v)[38:]) / (2 * LOOP)
        for v in te
    ])
    tv = np.array([
        3 * np.dot(vl[-33:] * (np.log(vl[-33:] / mu2) - 1 / 3),
                   np.diag(vu[:, -33:].T @ v @ vu[:, -33:])) / (2 * LOOP)
        for v in vt
    ])
    dp_scalar = -np.linalg.solve(common.radial_mass_map(r), ts)
    dp_vector = -np.linalg.solve(common.radial_mass_map(r), tv)
    dp = dp_scalar + dp_vector
    ct_scalar = rad.T @ common.mass_operator(dp_scalar) @ rad
    ct_vector = rad.T @ common.mass_operator(dp_vector) @ rad
    ct = common.mass_operator(dp)
    check("common_fixed_VEV_condition", ts + tv + rad.T @ ct @ x)
    hs, hv = np.zeros((3, 3)), np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            hs[i, j] = cw.hard_trace_hessian(
                lam, u, t[i], t[j], second[i, j], 290, mu2, 1.5, 1.)
            vv = gauge**2 * (orad[i].T @ orad[j] + orad[j].T @ orad[i])
            hv[i, j] = cw.hard_trace_hessian(
                vl, vu, vt[i], vt[j], vv, 33, mu2, 5 / 6, 3.)
    tree = rad.T @ h @ rad
    hard = tree + hs + hv + ct_scalar + ct_vector
    weights = bubble_slope_matrix(lam, 38)
    dz = np.einsum("ab,iab,jab->ij", weights, te, te, optimize=True) / (2 * LOOP)
    check("scalar_hard_Hessian_symmetric", hs, hs.T)
    check("vector_hard_Hessian_symmetric", hv, hv.T)
    check("scalar_kinetic_symmetric", dz, dz.T)
    tree_masses = eigh(tree, metric, eigvals_only=True)
    masses = eigh(hard, metric, eigvals_only=True)
    wave, wave_vectors = eigh(dz, metric)
    check("positive_scalar_bubble_weights", min(0., float(weights.min())))
    check("scalar_kinetic_Gram_PSD", min(0., float(wave.min())))
    result.update({
        "radial_tree_Hessian": tree.tolist(),
        "radial_scalar_CW_Hessian_before_CT": hs.tolist(),
        "radial_vector_CW_Hessian_before_CT": hv.tolist(),
        "radial_scalar_mass_CT": ct_scalar.tolist(),
        "radial_vector_mass_CT": ct_vector.tolist(),
        "radial_scalar_fixed_VEV_insertion": (hs + ct_scalar).tolist(),
        "radial_vector_fixed_VEV_insertion": (hv + ct_vector).tolist(),
        "radial_hard_Hessian": hard.tolist(),
        "radial_hard_scalar_kinetic": dz.tolist(),
        "radial_total_metric_scalar_hard_only": (metric + dz).tolist(),
        "tree_radial_eigenvalues": tree_masses.tolist(),
        "radial_hard_eigenvalues": masses.tolist(),
        "hard_scalar_kinetic_eigenvalues": wave.tolist(),
        "leading_kinetic_radial_vector": wave_vectors[:, -1].tolist(),
        "leading_kinetic_canonical_radial_fractions":
            (np.diag(metric) * wave_vectors[:, -1]**2).tolist(),
        "common_finite_mass_CT": dp.tolist(),
        "scalar_finite_mass_CT": dp_scalar.tolist(),
        "vector_finite_mass_CT": dp_vector.tolist(),
        "radial_scalar_tadpoles": ts.tolist(),
        "radial_vector_tadpoles": tv.tolist(),
        "minimum_hard_radial_mass2": float(masses.min()),
        "rho_Z": float(np.max(abs(wave))),
        "tree_radial_positive": bool(tree_masses.min() > 0),
        "kinetic_scope": "hard-hard and hard-soft scalar pairs; nonlocal soft-soft and all gauge/ghost/fermion momentum terms excluded",
    })
    if tree_masses.min() <= 0:
        return reject("Relative-insertion norm is undefined against a nonpositive tree radial Hessian")
    relative = eigh(hard - tree, tree, eigvals_only=True)
    rho_h = float(np.max(abs(relative)))
    rho_z = float(np.max(abs(wave)))
    identities_pass = all(c["pass"] for c in checks[check_start:])
    result.update({
        "relative_radial_insertions": relative.tolist(),
        "rho_H": rho_h,
        "radial_hard_positive": bool(masses.min() > 0),
        "mass_insertion_below_margin": rho_h < engineering_margin,
        "kinetic_insertion_below_margin": rho_z < engineering_margin,
        "numerical_identity_checks_pass": identities_pass,
        "accepted_for_next_radial_audit": bool(
            identities_pass and masses.min() > 0
            and rho_h < engineering_margin and rho_z < engineering_margin),
        "engineering_margin_is_convergence_theorem": False,
    })
    if not result["accepted_for_next_radial_audit"]:
        result["failure_reason"] = "One or more independent radial positivity/control gates failed"
    return result


def run(cache_dir):
    import verify_p54_invariant_decomposition as decomposition
    rf=Path(__file__).resolve().parents[1]
    oldpath=rf/'output/p54_full_doublet_cw.json'
    ledgerpath=rf/'output/p54_invariant_decomposition.json'
    old=json.loads(oldpath.read_text());ledger=json.loads(ledgerpath.read_text())
    if not ledger['summary']['all_pass']:
        raise ValueError('Invariant reconstruction gate is not closed')
    p1=common.yuk.module('nonuniform_p1',rf/'code/verify_p54_p1_hessian_spectrum.py')
    cw=common.yuk.module('nonuniform_cw',rf/'code/verify_p54_p2_bosonic_cw.py')
    search=common.yuk.module('nonuniform_search',rf/'code/search_p54_hierarchical_p1.py')
    p=old['tree_parameters'];r0=np.asarray(ledger['reference_radii']);mu=ledger['mu']
    rad=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    inv=decomposition.InvariantJets(p1,p,r0,cache_dir)
    checks=[];bubble_checks=verify_bubble_weights(checks)
    def check(name,a,b=0.,tol=3e-8):
        err=common.yuk.error(np.asarray(a),np.asarray(b))
        checks.append(dict(name=name,residual=err,tolerance=tol,pass_=bool(err<tol)))
    def passed(c):return c.get('pass',c.get('pass_',False))
    # These finite, coarse cards were selected AFTER attribution identified
    # the three transverse Sigma quartics, but BEFORE their candidate loops.
    # No experimental masses or flavor objectives enter this design.
    designs=[
      ('T30',.30,1.,1.,r0[1]),
      ('T10',.10,1.,1.,r0[1]),
      ('T03',.03,1.,1.,r0[1]),
      ('T03_L2',.03,2.,1.,r0[1]),
      ('T03_X',.03,1.,.03,r0[1]),
      ('T03_X_S30',.03,1.,.03,.30),
    ]
    rows=[]
    for label,tau,radial_factor,chi4_factor,sigma in designs:
        pars=dict(p)
        for key in ('lambda2','lambda4','lambda4p'):pars[key]*=tau
        pars['lambda0']*=radial_factor;pars['chi4']*=chi4_factor
        r=np.array([1.,sigma,.25]);h,t,second=inv.assemble(pars,r)
        row=evaluate_candidate(p1,cw,search,pars,r,h,t,second,mu,rad,checks,label,
                               gauge_coupling=mu)
        row['design']=dict(transverse_quartic_multiplier=tau,lambda0_multiplier=radial_factor,
                           chi4_multiplier=chi4_factor,sigma_over_omega=sigma)
        rows.append(row)
        print(label,{k:row.get(k) for k in ('failure_reason','minimum_hard_radial_mass2','rho_H','rho_Z','accepted_for_next_radial_audit')},flush=True)
    # Retained points must independently reproduce the unmodified action
    # using their FINAL reanchored and retuned parameter card. Also test one
    # rejected but evaluable point, so this regression is not vacuous.
    retained=[i for i,row in enumerate(rows) if row['accepted_for_next_radial_audit']]
    targets=retained.copy()
    if not targets:
        candidates=[i for i,row in enumerate(rows) if 'rho_H' in row]
        if candidates:targets=[min(candidates,key=lambda i:rows[i]['rho_H'])]
    for i in targets:
        row=rows[i];r=np.array(row['radii']);pars=row['tree_parameters']
        h,t,u=inv.assemble(pars,r);store=common.scalar.ActionJets(p1,pars,cache_dir);x=rad@r
        start=len(checks)
        direct=store.hessian(x,row['label']+'_independent_retuned_action')
        check(row['label']+'_independent_same_action_H',h,direct)
        for j in range(3):
            tj,uj=store.jets(x,direct,rad[:,j],row['label']+f'_radial_{j}')
            check(row['label']+f'_independent_same_action_T_{j}',t[j],tj)
            check(row['label']+f'_independent_same_action_U_{j}{j}',u[j,j],uj)
        direction=np.array([.37,-.23,.41]);tj,uj=store.jets(x,direct,rad@direction,row['label']+'_mixed_direction')
        check(row['label']+'_independent_same_action_mixed_T',np.einsum('r,rij->ij',direction,t),tj)
        check(row['label']+'_independent_same_action_mixed_U',np.einsum('r,s,rsij->ij',direction,direction,u),uj)
        row['independent_final_action_validation']=dict(passed=all(passed(c) for c in checks[start:]),
            cache_hits=store.hits,new_Hessians=store.evaluated)
        row['accepted_for_next_radial_audit'] &= row['independent_final_action_validation']['passed']
        p1.jax.clear_caches()
    # Strong negative controls for the radial-blind pure-H quartic.
    h,t,u=inv.assemble(p,r0)
    for name,xi1,xi2 in (('negative_null_H',-.1,.2),('negative_real_H',.1,-.2)):
        bad=dict(p,xi1=xi1,xi2=xi2)
        bh,bt,bu=inv.assemble(bad,r0)
        check(name+'_radial_H_unchanged_despite_unbounded_action',bh,h)
        check(name+'_radial_T_unchanged_despite_unbounded_action',bt,t)
        badrow=evaluate_candidate(p1,cw,search,bad,r0,bh,bt,bu,mu,rad,checks,name,gauge_coupling=mu)
        check(name+'_necessary_BFB_gate_rejects',float(badrow['accepted_for_next_radial_audit']))
    sources=[Path(__file__),Path(decomposition.__file__),Path(p1.__file__),Path(cw.__file__),
             Path(search.__file__),Path(common.__file__),Path(common.scalar.__file__),
             Path(common.yuk.__file__),oldpath,ledgerpath]
    return dict(schema='p54-nonuniform-scalar-reselection-v1',date='2026-09-14',
        scope='six attribution-guided engineering cards; no experimental mass fit',
        candidate_design_policy='only existing transverse Sigma quartics, lambda0, chi4 and one radius; dependent masses/light-doublet retuned',
        engineering_margin=.3,candidates=rows,bubble_weight_regression=bubble_checks,
        accepted_candidate_labels=[row['label'] for row in rows if row['accepted_for_next_radial_audit']],
        selected_physical_benchmark=None,default_parameters_changed=False,physical_fit_enabled=False,
        complete_mixed_BFB=False,complete_gauged_pole_kernel=False,
        source_sha256={str(p.relative_to(rf)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        checks=checks,summary=dict(passed=sum(passed(c) for c in checks),total=len(checks),all_pass=all(passed(c) for c in checks)))


def markdown(report):
    lines=['# Attribution-guided nonuniform P54 scalar cards','',
      'No default benchmark is changed. Passing a radial engineering margin is NOT a physical vacuum certificate.',
      f"Checks: {report['summary']['passed']}/{report['summary']['total']}.",'',
      '| Card | tau | lambda0 multiplier | chi4 multiplier | sigma | min hard radial mass2 | rho_H | rho_Z | Radial margin |',
      '|---|---:|---:|---:|---:|---:|---:|---:|:---:|']
    for row in report['candidates']:
        d=row['design'];prefix=f"| {row['label']} | {d['transverse_quartic_multiplier']} | {d['lambda0_multiplier']} | {d['chi4_multiplier']} | {d['sigma_over_omega']} |"
        if 'rho_H' in row:
            lines.append(prefix+f" {row['minimum_hard_radial_mass2']:.9g} | {row['rho_H']:.9g} | {row['rho_Z']:.9g} | {row['accepted_for_next_radial_audit']} |")
        else:lines.append(prefix+' tree gate rejected | - | - | False |')
    lines+=['','Full parameter cards, 328-state spectra, split radial CW/CT matrices, kinetic matrices and independent exact-action checks are retained in JSON.',
       'Remaining: full mixed BFB, competing stationary branches, full nonradial/momentum stability, loop light-doublet retuning, new running/thresholds, Wilson/box, lower matching and finite seesaw.']
    for row in report['candidates']:
        if 'rho_H' not in row:lines.append(f"- {row['label']}: {row.get('failure_reason')}; {row.get('tree_tuning')}")
    return '\n'.join(lines)+'\n'


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cache-dir',type=Path,required=True)
    args=ap.parse_args();report=run(args.cache_dir)
    output=Path(__file__).resolve().parents[1]/'output/p54_nonuniform_reselection'
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    output.with_suffix('.md').write_text(markdown(report))
    print(markdown(report))
    if not report['summary']['all_pass']:
        print('failed',[c for c in report['checks'] if not c.get('pass',c.get('pass_',False))])
        raise SystemExit(1)
