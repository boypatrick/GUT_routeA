#!/usr/bin/env python3
"""Local, gauge-covariant tree PS EFT with the upper scalar/vector fields removed.

This is a constructive two-derivative EFT, not a lower one-loop gauge
threshold.  The retained chart consists of the four declared PS scalar
parents (248 real coordinates); the gauge-neutral axion is held fixed.
The 24 SO(10)/PS directions are a coset, NOT a 24-dimensional subgroup.
Their algebraic vector response is therefore a PS-covariant one-form, not
a principal connection for an invented upper gauge subgroup.

Public functions below return the full computable local tensors rather than
parent-weighted mass logarithms.  Tensor inputs use the P1 canonical real
328-field convention.  All scalar jets are from the frozen tree action;
the selected Higgs ray is the exported bosonic-improved direction only.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
CACHE = ROOT / "tmp/p54_full_doublet_cw"
OUT = RF / "output/p54_lower_covariant_eft.json"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def sym(a):
    return (a + a.T) / 2


def from_cjson(x):
    return np.asarray(x["real"]) + 1j * np.asarray(x["imag"])


def invsqrt(a):
    d, u = np.linalg.eigh(sym(a))
    if d.min() <= 0:
        raise ValueError("Positive metric required")
    return (u / np.sqrt(d)) @ u.T


def scalar_valley_api(h, active_basis, heavy_basis):
    """Jets of BH^T grad V(x0+BA y+BH z(y))=0, in a PS-invariant slice.

    Unlike an already horizontal embedding, BA is the raw PS tensor
    embedding.  The upper-vector algebraic equation supplies the quotient
    metric subsequently. C>0 proves existence of the local analytic valley
    by the implicit-function theorem, not global validity of the chart.
    """
    ba, bh = active_basis, heavy_basis
    c = sym(bh.T @ h @ bh)
    if np.linalg.eigvalsh(c).min() <= 0:
        raise ValueError("The integrated scalar block must be positive")
    w = bh.T @ h @ ba
    d = -np.linalg.solve(c, w)
    t = ba + bh @ d
    return {"BA": ba, "BH": bh, "C": c, "W": w, "D": d, "T": t,
            "A": sym(ba.T @ h @ ba),
            "S": sym(t.T @ h @ t), "raw_metric": sym(t.T @ t),
            "scalar_p4_kernel": -sym(w.T @ np.linalg.solve(c, np.linalg.solve(c, np.linalg.solve(c, w))))}


def scalar_nonlocal_kernel_api(p_euclidean_squared, valley, vector):
    """Exact quadratic scalar Schur resolvent, before a derivative truncation.

    This is the scalar block with PS gauge fields set to zero, after the
    upper longitudinal-vector response is eliminated. The upper transverse
    vector determinant and the lower gauge/ghost blocks are NOT included.
    Continuation to a trial Minkowski mass uses p_E^2=-m^2 and is singular
    at eigenvalues of C; do not replace its poles by O(p^2) mass logs.
    """
    ba, c, w = valley["BA"], valley["C"], valley["W"]
    g0 = ba.T @ vector["horizontal_projector"] @ ba
    return sym(valley["A"] + p_euclidean_squared*g0
               - w.T @ np.linalg.solve(c+p_euclidean_squared*np.eye(len(c)), w))


def scalar_valley_directional_jet_api(valley, hessian_direction):
    """X_,qa = -BH C^-1 BH^T V'''[X_,q,X_,a].

    hessian_direction is the 328x328 directional derivative of the FULL
    tree Hessian along X_,q. It is not a selected principal mass block.
    """
    bh, c, t = valley["BH"], valley["C"], valley["T"]
    z_qa = -np.linalg.solve(c, bh.T @ hessian_direction @ t)
    return {"z_qa": z_qa, "X_qa": bh @ z_qa}


def upper_vector_api(x, tangent, upper_representations):
    """Eliminate 24 vector fields at O(D^2), covariantly under PS.

    With R=-iT and D_PS X+g O A_U, O_A=R_A X, I=O^T O and K=O^T T,
    A_U*=-g^-1 I^-1 K D y. The returned one_form is I^-1 K,
    before the explicit -1/g. M_U^2=g^2 I; the kinetic correction is
    -1/2 j^T M_U^-2 j, j=g K D y. Positive I makes this algebraic
    elimination well-defined locally. SO(10)/PS is not a subgroup.
    """
    o = np.column_stack([r @ x for r in upper_representations])
    inertia = sym(o.T @ o)
    if np.linalg.eigvalsh(inertia).min() <= 0:
        raise ValueError("Upper-vector inertia is not invertible")
    k = o.T @ tangent
    connection = np.linalg.solve(inertia, k)
    horizontal = sym(np.eye(len(x)) - o @ np.linalg.solve(inertia, o.T))
    raw = sym(tangent.T @ tangent)
    subtraction = sym(k.T @ connection)
    return {"O": o, "I": inertia, "K": k, "one_form": connection,
            "horizontal_projector": horizontal, "raw_metric": raw,
            "current_current_metric": subtraction, "metric": sym(raw - subtraction)}


def upper_vector_directional_jet_api(vector, x_q, tangent, tangent_q, upper_representations):
    """Exact field derivative of the upper-vector response and metric."""
    o, inertia, k, a = (vector[key] for key in ("O", "I", "K", "one_form"))
    oq = np.column_stack([r @ x_q for r in upper_representations])
    iq = sym(oq.T @ o + o.T @ oq)
    kq = oq.T @ tangent + o.T @ tangent_q
    aq = np.linalg.solve(inertia, kq - iq @ a)
    rawq = sym(tangent_q.T @ tangent + tangent.T @ tangent_q)
    metricq = sym(rawq - kq.T @ a - k.T @ aq)
    return {"O_q": oq, "I_q": iq, "K_q": kq, "one_form_q": aq,
            "metric_q": metricq, "raw_metric_q": rawq}


def spectrum(a, tolerance=1e-9):
    d = np.linalg.eigvalsh(sym(a))
    return {"min": float(d.min()), "max": float(d.max()),
            "rank": int(np.sum(abs(d) > tolerance)), "trace": float(d.sum()),
            "nonzero_eigenvalues": d[abs(d) > tolerance].tolist()}


def run(include_tensors=False):
    """Run diagnostics; optionally expose the same full tensors to a matcher.

    ``run(include_tensors=True)`` returns ``{"report": ..., "tensors": ...}``.
    The tensor payload deliberately stays in memory (NumPy arrays and the
    exact potential callable) and is not reduced to JSON spectral weights.
    """
    started = time.time()
    p1path = RF / "code/verify_p54_p1_hessian_spectrum.py"
    upperpath = RF / "code/verify_p54_ps_finite_thresholds.py"
    helperspath = RF / "code/verify_p54_p2_two_site_matching.py"
    fullpath = RF / "output/p54_full_doublet_cw.json"
    finitepath = RF / "output/p54_ps_finite_thresholds.json"
    tripletpath = RF / "output/p54_typeii_triplet_source.json"
    p1 = module("p54_lower_p1", p1path)
    upper = module("p54_lower_upper", upperpath)
    helper = module("p54_lower_helpers", helperspath)
    full, finite, triplet = [json.loads(p.read_text()) for p in (fullpath, finitepath, tripletpath)]
    pars, v = full["tree_parameters"], full["vacuum"]
    x0 = p1.vacuum_vector(v["omega"], v["sigma"], v["vs"])
    vu = finite["upper_PS"]["vacuum"]
    xu = p1.vacuum_vector(vu["omega"], vu["sigma"], vu["vs"])
    gu = full["scheme"]["mu_over_omega"]
    potential = p1.potential_factory(pars)
    hf = p1.jax.jit(p1.hessian(potential))
    keybase = hashlib.sha256(p1path.read_bytes() + json.dumps(pars, sort_keys=True).encode()).digest()
    counts = {"cache_hits": 0, "evaluated": 0}

    def hessian(x, label):
        tick = time.time()
        key = hashlib.sha256(keybase + np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = CACHE / (key + ".npz")
        if path.exists():
            with np.load(path) as cached:
                result = cached["h"]
            counts["cache_hits"] += 1
        else:
            result = sym(np.asarray(hf(p1.anp.asarray(x)), float))
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, h=result)
            counts["evaluated"] += 1
        print(f"{label}: {time.time()-tick:.2f}s", flush=True)
        return result

    h0 = hessian(x0, "actual tree background")
    psgen = helper.ps_generators(p1)
    ps = [r for rs in psgen.values() for r in rs]
    casimir = helper.casimir_operators(p1, psgen)
    parents = helper.parent_projectors(p1, casimir)
    pa = sum(row["projector"] for row in parents if row["label"] in upper.ACTIVE)
    ba = upper.range_basis(pa)
    gold_upper = upper.image_basis(p1.gauge_orbit(xu))
    pgu_upper = gold_upper @ gold_upper.T
    pq = p1.pq_direction(xu)
    pq -= pgu_upper @ pq
    pq /= np.linalg.norm(pq)
    bh = upper.range_basis(np.eye(328) - pa - pgu_upper - np.outer(pq, pq))
    print(f"PS chart: {ba.shape[1]} active, {bh.shape[1]} heavy", flush=True)
    ug = [helper.so_generator(a, b) for a in range(6) for b in range(6, 10)]
    ru = [p1.representation_matrix(t) for t in ug]
    rp = [p1.representation_matrix(t) for t in ps]
    valley = scalar_valley_api(h0, ba, bh)
    c, t, s = (valley[key] for key in ("C", "T", "S"))
    vector = upper_vector_api(x0, t, ru)
    metric, o = vector["metric"], vector["O"]
    pgu = np.eye(328) - vector["horizontal_projector"]
    old = upper.heavy_valley_api(h0, ba, bh, pgu)

    # Lie brackets identify a symmetric coset, not an upper subgroup.
    projection = lambda m, basis: sum(np.sum(b*m) * b for b in basis)
    lie_pu_error, lie_uu_ps_error, uu_ps_norm = 0., 0., 0.
    for a in ps:
        for b in ug:
            bracket = a @ b - b @ a
            lie_pu_error = max(lie_pu_error, np.linalg.norm(bracket-projection(bracket, ug)))
    for i, a in enumerate(ug):
        for b in ug[i+1:]:
            bracket = a @ b - b @ a
            lie_uu_ps_error = max(lie_uu_ps_error, np.linalg.norm(bracket-projection(bracket, ps)))
            uu_ps_norm = max(uu_ps_norm, np.linalg.norm(projection(bracket, ps)))

    # Manifest covariance of the raw chart and the actual two-derivative metric.
    ph = bh @ bh.T
    ps_chart_error = max(np.linalg.norm(pa @ r-r @ pa) for r in rp)
    heavy_ps_error = max(np.linalg.norm(ph @ r-r @ ph) for r in rp)
    smreps = p1.sm_representation_matrices()
    sm = [r for rows in smreps.values() for r in rows]
    sm_tangent_error = max(np.linalg.norm(r @ t-t @ (ba.T @ r @ ba)) for r in sm)
    sm_metric_error = max(np.linalg.norm(metric @ (ba.T @ r @ ba)-(ba.T @ r @ ba) @ metric) for r in sm)
    sm_mass_error = max(np.linalg.norm(s @ (ba.T @ r @ ba)-(ba.T @ r @ ba) @ s) for r in sm)
    # A finite broken-PS rotation must leave the tensor construction covariant.
    ops = np.column_stack([r @ x0 for r in rp])
    broken = int(np.argmax(np.linalg.norm(ops, axis=0)))
    rotation = expm(.173 * rp[broken])
    rotated = upper_vector_api(rotation @ x0, rotation @ t, ru)
    finite_covariance_error = np.linalg.norm(rotated["metric"]-metric)
    lower_coordinates = ba.T @ ops
    lower_tangent_error = np.linalg.norm(t @ lower_coordinates-ops)
    lower_mass_full = gu**2 * ops.T @ ops
    lower_mass_eft = gu**2 * lower_coordinates.T @ metric @ lower_coordinates
    lower_mass_error = np.linalg.norm(lower_mass_full-lower_mass_eft)
    lower_goldstones = upper.image_basis(lower_coordinates)
    lower_orth = lower_goldstones @ invsqrt(lower_goldstones.T @ metric @ lower_goldstones)
    lower_goldstone_norm_error = np.linalg.norm(lower_orth.T @ metric @ lower_orth-np.eye(lower_orth.shape[1]))
    lower_goldstone_ward = np.linalg.norm(s @ lower_orth)

    # Two existing finite-difference steps yield exact cubic/quartic jets for
    # the quartic tree polynomial, with a numerical cancellation check.
    candidates = [row for row in full["basis"].values() if isinstance(row, dict)
                  and "real" in row and np.asarray(row["real"]).shape == (328, 4)]
    if len(candidates) != 1:
        raise ValueError("Unique exported complex doublet embedding required")
    bd = from_cjson(candidates[0])
    coeff = from_cjson(full["retuned_bosonic_eigenpair"]["light_coefficients"])
    q = math.sqrt(2) * (bd @ coeff).real
    aq = ba.T @ q
    dhs, v4s = [], []
    for step in (.02, .01):
        hp, hm = hessian(x0+step*q, f"q +{step}"), hessian(x0-step*q, f"q -{step}")
        dhs.append((hp-hm)/(2*step))
        v4s.append(float(q @ ((hp+hm-2*h0)/step**2) @ q))
    dhq = dhs[0]
    scalarjet = scalar_valley_directional_jet_api(valley, dhq)
    tq = scalarjet["X_qa"]
    vectorjet = upper_vector_directional_jet_api(vector, q, t, tq, ru)
    j = dhq @ q
    jh = bh.T @ j
    zqq = -np.linalg.solve(c, jh)
    xqq = bh @ zqq
    heavy_exchange = float(jh @ np.linalg.solve(c, jh))
    lambda_direct = v4s[0] / 6
    lambda_relaxed = lambda_direct-heavy_exchange/2
    step_error = np.linalg.norm(dhs[0]-dhs[1])/max(np.linalg.norm(dhs[0]), 1e-14)
    quartic_step_error = abs(v4s[0]-v4s[1])/max(abs(v4s[0]), 1e-14)
    # This identity probes the actual V''' and the nonlinear valley map.
    nonlinear_ward = []
    potential_ward = []
    for r, orbit in zip(rp, ops.T):
        ta = ba.T @ r @ ba
        x_qk = -bh @ np.linalg.solve(c, bh.T @ dhq @ orbit)
        nonlinear_ward.append(np.linalg.norm(r @ q-t @ ta @ aq-x_qk))
        potential_ward.append(np.linalg.norm(t.T @ (dhq @ orbit+h0 @ (r @ q)-r @ (h0 @ q))))

    # Independent finite difference of the algebraic vector map checks its jet.
    eps = 1e-5
    plus = upper_vector_api(x0+eps*q, t+eps*tq, ru)
    minus = upper_vector_api(x0-eps*q, t-eps*tq, ru)
    metricjet_fd_error = np.linalg.norm((plus["metric"]-minus["metric"])/(2*eps)-vectorjet["metric_q"])
    connectionjet_fd_error = np.linalg.norm((plus["one_form"]-minus["one_form"])/(2*eps)-vectorjet["one_form_q"])
    # Along the color-singlet Higgs/heavy-relaxation ray, no upper-color
    # vector current exists to this order. This is a selection rule, not fit.
    radial_current = max(np.linalg.norm(np.column_stack([r @ xx for r in ru]).T @ yy)
                         for xx, yy in ((x0,q), (x0,xqq), (q,xqq)))
    radial_kinetic_h2 = float(zqq @ zqq)

    # Actual trilinear in the retained 126 triplet after ALL 55 upper heavy
    # fields have been removed.  Canonical normalization is not omitted.
    bt = from_cjson(triplet["complex_triplet_basis_328x2"])
    tr126 = np.column_stack((math.sqrt(2)*bt[:,1].real, -math.sqrt(2)*bt[:,1].imag))
    ad = ba.T @ tr126
    td = t @ ad
    gd = sym(ad.T @ metric @ ad)
    sd = sym(ad.T @ s @ ad)
    jd = td.T @ j
    normd = invsqrt(gd)
    yresponse = -.5*np.linalg.solve(sd, jd)
    canonical_mass = sym(normd.T @ sd @ normd)
    canonical_source = normd.T @ jd
    reconstructed = td @ yresponse + .5*xqq
    trfull = np.column_stack((math.sqrt(2)*bt.real, -math.sqrt(2)*bt.imag))
    triplet_response = trfull.T @ reconstructed
    reference_response = np.asarray(triplet["canonical_induced_z_over_v2_times_omega"])
    triplet_response_error = np.linalg.norm(triplet_response-reference_response)
    # Direct H^4 and the upper scalar p^4 kernel are distinct Wilson data.
    p4 = valley["scalar_p4_kernel"]
    p4_values = np.linalg.eigvalsh(p4)
    metric_values = np.linalg.eigvalsh(metric)
    eigenvalues, eigenvectors = eigh(s, metric)
    # A global smallest-gap comparison is overconservative when selection
    # rules kill W. Resolve the actual support of every scalar source.
    dc, uc = np.linalg.eigh(c)
    source_weights = uc.T @ valley["W"] @ eigenvectors
    coupling_groups = {}
    coupled_count, beyond_radius = 0, 0
    near_mode = None
    for jmode, mass in enumerate(eigenvalues):
        source = source_weights[:, jmode]
        use = abs(source) > max(1e-10, 1e-8*np.linalg.norm(source))
        if not np.any(use):
            continue
        coupled_count += 1
        ratio = float(mass/dc[use].min())
        beyond_radius += int(ratio >= 1)
        if ratio >= 1 and near_mode is None:
            near_mode = jmode
        key = (round(float(mass), 9), round(float(dc[use].min()), 9))
        if key not in coupling_groups:
            coupling_groups[key] = {"retained_m2":float(mass), "multiplicity":0,
                "smallest_actually_coupled_C_eigenvalue":float(dc[use].min()),
                "mass_to_coupled_C_ratio":ratio,"source_norm":float(np.linalg.norm(source)),
                "source_support_C_eigenvalues":sorted(set(round(float(v), 12) for v in dc[use]))}
        coupling_groups[key]["multiplicity"] += 1
    resolvent_zero_error = np.linalg.norm(scalar_nonlocal_kernel_api(0., valley, vector)-s)
    p4_remainders = []
    for momentum in (1e-3, 5e-4):
        remainder = scalar_nonlocal_kernel_api(momentum,valley,vector)-s-momentum*metric-momentum**2*p4
        p4_remainders.append(float(np.linalg.norm(remainder)))
    # Resolve the actual near-crossing by a closed canonical 2x2 block,
    # rather than treating Taylor nonconvergence as a physical obstruction.
    aa = eigenvectors[:, near_mode]
    active_pair = vector["horizontal_projector"] @ ba @ aa
    active_pair /= np.linalg.norm(active_pair)
    heavy_pair = bh @ valley["W"] @ aa
    heavy_pair /= np.linalg.norm(heavy_pair)
    pair_basis = np.column_stack((active_pair, heavy_pair))
    pair_mass = sym(pair_basis.T @ h0 @ pair_basis)
    pair_closure_error = np.linalg.norm(h0@pair_basis-pair_basis@pair_mass)
    pair_poles = np.linalg.eigvalsh(pair_mass)
    selected_connection = np.unravel_index(np.argmax(abs(vectorjet["one_form_q"])), vectorjet["one_form_q"].shape)
    tensor_residuals = {
        "active_PS_projector": float(ps_chart_error), "heavy_PS_projector": float(heavy_ps_error),
        "SM_valley_tangent": float(sm_tangent_error), "SM_metric": float(sm_metric_error),
        "SM_mass": float(sm_mass_error), "finite_broken_PS_metric_covariance": float(finite_covariance_error),
        "lower_orbit_tangent": float(lower_tangent_error), "lower_vector_mass": float(lower_mass_error),
        "lower_Goldstone_metric_normalization": float(lower_goldstone_norm_error),
        "lower_Goldstone_mass_Ward": float(lower_goldstone_ward),
        "nonlinear_valley_q_orbit_Ward": float(max(nonlinear_ward)),
        "full_V3_q_orbit_Ward": float(max(potential_ward)),
        "metric_directional_jet_finite_difference": float(metricjet_fd_error),
        "connection_directional_jet_finite_difference": float(connectionjet_fd_error),
        "cubic_step_relative": float(step_error), "quartic_step_relative": float(quartic_step_error),
        "triplet_source_response_reconstruction": float(triplet_response_error)}
    checks = [
        ("raw PS-invariant chart has 248 active and 55 integrated real scalars", ba.shape[1]==248 and bh.shape[1]==55 and max(ps_chart_error,heavy_ps_error)<1e-10),
        ("positive scalar block gives a local analytic heavy valley", np.linalg.eigvalsh(c).min()>1e-3),
        ("upper 24-vector inertia is positive", np.linalg.eigvalsh(vector["I"]).min()>1e-3),
        ("SO10/PS is a symmetric coset, not a subgroup", max(lie_pu_error,lie_uu_ps_error)<1e-12 and uu_ps_norm>.1),
        ("scalar-heavy directions are orthogonal to actual upper gauge orbit", np.linalg.norm(bh.T @ o)<1e-11),
        ("raw scalar valley followed by vector elimination reproduces prior horizontal G and S", max(np.linalg.norm(metric-old["G"]),np.linalg.norm(s-old["S"]))<1e-10),
        ("two-derivative metric is positive", metric_values.min()>.9),
        ("actual unbroken SM tangent, mass and metric Ward identities", max(sm_tangent_error,sm_mass_error,sm_metric_error)<1e-10),
        ("finite broken-PS rotation preserves vector-eliminated metric", finite_covariance_error<1e-10),
        ("lower PS/SM orbit has rank nine and correct tangent", lower_goldstones.shape[1]==9 and lower_tangent_error<1e-10),
        ("lower gauge mass normalization agrees before and after elimination", lower_mass_error<1e-10),
        ("all nine lower Goldstones are metric-normalized zero modes", max(lower_goldstone_norm_error,lower_goldstone_ward)<1e-10),
        ("retained spectrum has thirteen zero and no negative tree modes", np.sum(abs(eigenvalues)<1e-8)==13 and np.sum(eigenvalues < -1e-8)==0),
        ("exported Higgs lies in the raw valley tangent with no linear heavy response", np.linalg.norm(t@aq-q)<1e-10 and abs(q@q-1)<1e-10),
        ("full cubic and quartic jets agree at two cached steps", step_error<1e-8 and quartic_step_error<1e-8),
        ("nonlinear heavy-valley equivariance passes q-orbit V3 test", max(nonlinear_ward)<1e-9),
        ("full scalar cubic obeys differentiated PS potential Ward identity", max(potential_ward)<1e-9),
        ("covariant response one-form and metric derivatives pass independent finite difference", max(metricjet_fd_error,connectionjet_fd_error)<1e-8),
        ("Higgs-ray upper vector current vanishes by color selection", radial_current<1e-10),
        ("retained canonical triplet block remains positive", np.linalg.eigvalsh(canonical_mass).min()>0),
        ("upper-relaxed triplet cubic reproduces full 54/126 mixed response", triplet_response_error<1e-9),
        ("upper scalar exchange produces negative semidefinite p4 quadratic kernel", p4_values.max()<1e-9 and p4_values.min() < -1e-5),
        ("exact scalar Schur resolvent reproduces S, G and p4 with cubic remainder", resolvent_zero_error<1e-10 and .10<p4_remainders[1]/p4_remainders[0]<.15),
        ("near-crossing is an invariant canonical 2x2 block with positive exact poles", pair_closure_error<1e-10 and pair_poles.min()>0),
    ]
    checkrows = [{"name":name,"pass":bool(passed)} for name,passed in checks]
    sourcepaths = [p1path, upperpath, helperspath, fullpath, finitepath, tripletpath, Path(__file__).resolve()]
    report = {
        "schema":"p54-lower-covariant-tree-eft-v1", "date":"2026-09-06",
        "status":{"local_tree_two_derivative_EFT_constructed":True,
                  "lower_one_loop_finite_gauge_matching_completed":False,
                  "upper_one_loop_Wilson_vertices_completed":False,
                  "global_valley_continuation_proved":False,
                  "fermionic_matching_completed":False},
        "scope":{"chart":"Four declared PS scalar parents, upper Phi54 unitary gauge; gauge-neutral axion fixed as spectator",
                 "action":"Same frozen tree parameters and actual stationary vacuum as P1/P2; no scalar CT inserted in these tree jets",
                 "loop_scheme":"One-loop fixed-VEV retuning needs CT and loop shifts of Veff and metric at O(hbar); these are not included. Insertion of such shifts into an already one-loop gauge threshold is O(hbar^2), but tree EFT coefficients themselves shift at O(hbar).",
                 "Higgs_ray":"Exported bosonic-improved complex light eigenvector; evaluating frozen tree tensors along it is a mixed-order diagnostic, not a loop-complete Higgs Wilson coefficient",
                 "local_validity":"Analytic implicit-function germ at actual vacuum while C and upper inertia remain positive; no global upper-to-lower path or all-momentum separation is inferred",
                 "derivative_counting":"Scalar two-covariant-derivative action plus leading PS Yang-Mills. Induced F_U and [A_U,A_U]_PS terms are covariant p4, not computed finite matching."},
        "vacuum":v,"gauge_coupling":gu,"active_parent_labels":upper.ACTIVE,
        "functional_definition":{
            "valley":"X(y)=x0+BA y+BH z(y); BH^T grad V(X)=0; z(0)=0",
            "valley_first_jet":"C=BH^T H BH, W=BH^T H BA, D=-C^-1 W, T=BA+BH D",
            "valley_second_jet":"X_,ab=-BH C^-1 BH^T V'''[T_a,T_b]",
            "PS_transformation":"k_p(y)=BA^T R_p(x0+BA y); z transforms in BH; k is affine, not purely linear at broken vacuum",
            "gauge_convention":"R=-iT, D=partial-ig A T=partial+g A R; scalar kinetic is |D_PS X+g O A_U|^2/2",
            "upper_vector":"O_A=R_A X, I=O^T O, j=g O^T D_PS X, M_U^2=g^2 I; A_U*=-g^-1 I^-1 O^T D_PS X",
            "scalar_metric":"g_ab=T_a^T (1-O I^-1 O^T) T_b; L2=1/2 g_ab D y^a D y^b-V(X)",
            "current_current":"Delta L2=-1/2 j^T M_U^-2 j",
            "coset_warning":"[PS,U] subset U and [U,U] subset PS is nonzero: I^-1 O^T T is a PS-covariant heavy-vector one-form, not a principal connection for an upper subgroup",
            "effective_cubic":"Veff,abc=V'''[T_a,T_b,T_c] at the stationary valley, since BH^T H T=0",
            "Higgs_quartic":"lambda_upper=V''''[q,q,q,q]/6-(BH^T V'''[q,q])^T C^-1(BH^T V'''[q,q])/2",
            "Higgs_kinetic":"Along X(h)=x0+q h+BH z_qq h^2/2+..., g_hh=1+||z_qq||^2 h^2+O(h^3), for this zero-upper-current Higgs ray",
            "scalar_p4":"K_eff(p)=S+p^2 G-p^4 W^T C^-3 W+... before additional vector-kinetic p4 terms",
            "vector_p4":"Reinsert A_U* into F_U=D_PS A_U* and F_PS=F_PS^(0)+g[A_U*,A_U*]_PS; their squares/cross terms are required beyond the declared scalar O(D^2) truncation"},
        "local_tensors":{
            "C_eigenvalues":np.linalg.eigvalsh(c).tolist(),"D_operator_norm":float(np.linalg.norm(valley["D"],2)),
            "upper_inertia_eigenvalues":np.linalg.eigvalsh(vector["I"]).tolist(),
            "upper_vector_m2_eigenvalues":(gu**2*np.linalg.eigvalsh(vector["I"])).tolist(),
            "raw_metric":spectrum(valley["raw_metric"]),
            "current_current_metric":spectrum(vector["current_current_metric"]),
            "quotient_metric_min":float(metric_values.min()),"quotient_metric_max":float(metric_values.max()),
            "prior_horizontal_metric_error":float(np.linalg.norm(metric-old["G"])),
            "prior_horizontal_mass_error":float(np.linalg.norm(s-old["S"])),
            "upper_one_form_operator_norm":float(np.linalg.norm(vector["one_form"],2)),
            "upper_one_form_q_operator_norm":float(np.linalg.norm(vectorjet["one_form_q"],2)),
            "metric_q_operator_norm":float(np.linalg.norm(vectorjet["metric_q"],2)),
            "metric_q_Frobenius_norm":float(np.linalg.norm(vectorjet["metric_q"])),
            "selected_one_form_q_entry":{"coset_generator_pair":[int(selected_connection[0]//4),int(6+selected_connection[0]%4)],
                 "active_chart_column":int(selected_connection[1]),"coefficient":float(vectorjet["one_form_q"][selected_connection]),
                 "warning":"Raw basis component depends on deterministic tensor eigenspace basis; operator norms are basis-independent."},
            "scalar_p4_kernel_min":float(p4_values.min()),"scalar_p4_kernel_max":float(p4_values.max())},
        "derivative_expansion_gate":{
            "retained_generalized_m2_max":float(eigenvalues.max()),
            "retained_to_integrated_min_m2_ratio":float(eigenvalues.max()/np.linalg.eigvalsh(c).min()),
            "retained_generalized_modes_above_smallest_integrated_m2":int(np.sum(eigenvalues>np.linalg.eigvalsh(c).min())),
            "linearly_coupled_retained_modes":coupled_count,
            "linearly_decoupled_retained_modes":len(eigenvalues)-coupled_count,
            "coupled_modes_at_or_beyond_Taylor_radius":beyond_radius,
            "coupling_resolved_groups":list(coupling_groups.values()),
            "exact_resolvent_zero_error":float(resolvent_zero_error),
            "p4_subtracted_remainder_at_pE2_0p001_0p0005":p4_remainders,
            "near_crossing_exact_block":{
                "canonical_mass_matrix":pair_mass.tolist(),"full_Hessian_closure_error":float(pair_closure_error),
                "exact_tree_poles_m2":pair_poles.tolist(),"retained_O_p2_m2":float(eigenvalues[near_mode]),
                "relative_retained_O_p2_error":float(abs(pair_poles[1]-eigenvalues[near_mode])/pair_poles[1]),
                "interpretation":"The source-selected two-dimensional block is invariant under the full 328-field Hessian. Exact diagonalization resolves the local near-crossing; no one-loop threshold is inferred."},
            "interpretation":"The global smallest-gap comparison is only sufficient and is overconservative for W-decoupled sectors. Actual source support leaves 25 linearly coupled retained modes; 12 lie just beyond their coupled C Taylor radius and require the exact resolvent, despite small mixing. This does not invalidate the local two-derivative EFT or prove a model no-go. Generalized masses must not be blindly promoted to finite-threshold logarithms."},
        "Lie_algebra":{"PS_U_closure_error":float(lie_pu_error),"U_U_PS_closure_error":float(lie_uu_ps_error),
                       "nonzero_U_U_PS_bracket_max_norm":float(uu_ps_norm)},
        "lower_gauge":{"m2_eigenvalues":np.linalg.eigvalsh(sym(lower_mass_eft)).tolist(),
                       "massive_rank":lower_goldstones.shape[1],"metric_normalized_Goldstone_count":lower_orth.shape[1]},
        "Higgs_local_vertices":{"coordinate_convention":"h is the specified UV-projection ray coordinate with unit metric at the origin, not a globally/geodesically canonical coordinate; all 248 active parents are still retained, so this is not SM lambda",
                    "tree_ray_m2":float(q @ h0 @ q),"tree_direct_V4":v4s[0],
                    "tree_direct_lambda":lambda_direct,"upper_heavy_exchange_J_Cinv_J":heavy_exchange,
                    "upper_relaxed_lambda":lambda_relaxed,"heavy_valley_curvature_norm":float(np.linalg.norm(zqq)),
                    "g_hh_h2_coefficient":radial_kinetic_h2,"L_h2_dh2_coefficient":radial_kinetic_h2/2,
                    "upper_current_selection_rule_residual":float(radial_current)},
        "retained_126_neutral_triplet_vertices":{
                    "real_basis_order":["sqrt2 Re B126","-sqrt2 Im B126"],
                    "metric":gd.tolist(),"Schur_m2":sd.tolist(),"cubic_J_hh":jd.tolist(),
                    "canonical_m2":canonical_mass.tolist(),"canonical_cubic_J_hh":canonical_source.tolist(),
                    "retained_coordinate_over_h2":yresponse.tolist(),
                    "reconstructed_full_triplet_over_h2":triplet_response.tolist(),
                    "no_fermionic_Clebsch_or_neutrino_matrix_claim":True},
        "tensor_residuals":tensor_residuals,
        "reusable_API":{"module":str(Path(__file__).resolve().relative_to(ROOT)),
                    "functions":["run(include_tensors=True)","scalar_valley_api","scalar_nonlocal_kernel_api","scalar_valley_directional_jet_api","upper_vector_api","upper_vector_directional_jet_api"],
                    "input_basis":"BA is raw sum-of-four-PS-parent image; BH is upper physical heavy image exported by the prior finite-threshold construction; RU are the actual 24 canonical SO10/PS scalar representation matrices",
                    "available_full_local_tensors":"C,W,D,T,S,G; upper inertia I and current K; covariant response one-form; X_,qa, one-form derivative and metric derivative; scalar p4 quadratic kernel; exact scalar quadratic Schur resolvent at supplied p_E^2",
                    "lower_matcher_still_needs":["upper one-loop Wilson action including higher-dimensional scalar/vector/fermion operators in a common scheme",
                        "full lower covariant gauge-scalar-ghost quadratic fluctuation operator using this metric and connection",
                        "consistent hard-minus-EFT subtraction, ghosts and site-resolved Goldstone measure",
                        "regulator and momentum/derivative-order declaration beyond the displayed local jets"]},
        "references":[{"title":"Henning, Lu, Murayama, One-loop Matching and Running with Covariant Derivative Expansion", "url":"https://arxiv.org/abs/1604.01019", "use":"Gauge-covariant heavy-field functional matching; present calculation is tree-level only"},
                      {"title":"Henning, Lu, Murayama, How to use the Standard Model effective field theory", "url":"https://arxiv.org/abs/1412.1837", "use":"Covariant derivative expansion and operator-level EFT matching"}],
        "checks":checkrows,"summary":{"passed":sum(r["pass"] for r in checkrows),"total":len(checkrows)},
        "cache":counts,"runtime_seconds":time.time()-started,
        "sources":[{"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sourcepaths]}
    if include_tensors:
        return {"report":report, "tensors":{
            "x0":x0, "H0":h0, "potential":potential, "p1":p1, "parameters":pars,
            "BA":ba, "BH":bh, "PS_representations":rp, "upper_representations":ru,
            "active_vacuum_coordinates":ba.T@x0, "heavy_vacuum_coordinates":bh.T@x0,
            "valley":valley, "upper_vector":vector, "Higgs_direction":q,
            "Hessian_derivative_q":dhq, "valley_directional_jet_q":scalarjet,
            "upper_vector_directional_jet_q":vectorjet,
            "lower_orbit_coordinates":lower_coordinates, "lower_Goldstone_basis":lower_orth,
            "retained_triplet_coordinates":ad, "canonical_triplet_coordinates":ad@normd}}
    return report


def markdown(r):
    a, h, d = r["local_tensors"], r["Higgs_local_vertices"], r["retained_126_neutral_triplet_vertices"]
    lines = ["# P54 lower covariant tree EFT", "",
        "Constructed a local PS-covariant tree scalar EFT through two covariant derivatives, including 55 scalar-valley responses and algebraic elimination of the 24 upper vectors. This is **not** completed lower one-loop finite gauge matching.", "",
        "## Exact local construction", "",
        r"The raw 248-dimensional PS tensor chart is $X(y)=x_0+B_Ay+B_Hz(y)$ with $B_H^T\nabla V(X)=0$. Its positive heavy block defines an analytic germ by the implicit-function theorem. The axion is held fixed as a gauge-neutral spectator.", "",
        r"$$C=B_H^THB_H,\quad D=-C^{-1}B_H^THB_A,\quad T=B_A+B_HD,\quad X_{,ab}=-B_HC^{-1}B_H^TV^{(3)}[T_a,T_b].$$", "",
        r"With $R=-iT$, the card's $D=\partial-igAT$ is $D=\partial+gAR$. For $O_A=R_AX$, $I=O^TO$ and $j=gO^TD_{PS}X$,", "",
        r"$$A_U^*=-g^{-1}I^{-1}O^TD_{PS}X,\qquad \Delta\mathcal L_2=-\frac12j^T(g^2I)^{-1}j,$$", "",
        r"$$\mathcal L_2=\frac12(Dy)^T\underbrace{T^T(1-OI^{-1}O^T)T}_{g(y)}(Dy)-V(X(y)).$$", "",
        r"The 24 directions are the symmetric coset $SO(10)/PS$, not a subgroup: $[PS,U]\subset U$, $[U,U]\subset PS\ne0$. Thus $I^{-1}O^TT$ is a PS-covariant response one-form, not a principal connection for a fictitious 24-dimensional group.", "",
        "## Actual numerical tensors", "",
        f"- Heavy scalar gap: `{min(a['C_eigenvalues']):.12g}`; quotient metric range: `[{a['quotient_metric_min']:.12g}, {a['quotient_metric_max']:.12g}]`.",
        f"- Upper-vector current-current metric subtraction has rank `{a['current_current_metric']['rank']}`, trace `{a['current_current_metric']['trace']:.12g}`.",
        f"- Higgs-direction response derivative norm: `{a['upper_one_form_q_operator_norm']:.12g}`; metric-derivative norm: `{a['metric_q_operator_norm']:.12g}`.",
        f"- Previous horizontal Schur metric is reproduced with error `{a['prior_horizontal_metric_error']:.3e}`. All nine lower gauge masses and Goldstone normalizations agree.", "",
        "## Beyond the Hessian", "",
        r"For the exported corrected neutral direction $q$, $J_H=B_H^TV^{(3)}[q,q]$,", "",
        r"$$\lambda_{\rm upper}=\frac16V^{(4)}[q^4]-\frac12J_H^TC^{-1}J_H,\qquad z_{,hh}=-C^{-1}J_H.$$",
        "",
        f"- Direct tree quartic `{h['tree_direct_lambda']:.12g}` becomes `{h['upper_relaxed_lambda']:.12g}` after all 55 upper scalars relax.",
        f"- The induced radial derivative operator is `{h['L_h2_dh2_coefficient']:.12g} h^2 (D h)^2` in omega units. The upper-vector current on this color-singlet ray vanishes by selection rules.",
        f"- Canonically normalized retained 126 triplet masses: `{np.linalg.eigvalsh(np.asarray(d['canonical_m2'])).tolist()}`.",
        f"- Its actual effective hh source: `{d['canonical_cubic_J_hh']}`; reconstruction agrees with the full mixed 54/126 response to `{r['tensor_residuals']['triplet_source_response_reconstruction']:.3e}`.",
        "- These are frozen tree tensors evaluated on a bosonic-improved ray, not loop-complete Higgs/type-II Wilson coefficients. The ray coordinate has unit metric only at the origin; all 248 active parents remain. Thus this is neither an SM quartic nor an all-orders geodesic-coordinate coupling. No fermion Clebsch is inferred.", "",
        "## Scope and remaining matching", "",
        r"Scalar elimination already produces $-p^4W^TC^{-3}W$. Substitution of $A_U^*$ into $F_U=D_{PS}A_U^*$ and $F_{PS}=F_{PS}^{(0)}+g[A_U^*,A_U^*]_{PS}$ generates further covariant order-p4 operators. They cannot be silently discarded inside a claimed full one-loop finite threshold.", "",
        "The retained PS transformation is affine about the broken vacuum. The verifier checks projector covariance, SM Ward identities, a finite broken-PS rotation, nonlinear V3 q-orbit identities, all nine lower Goldstone normalizations, and two independent jet steps. Positivity proves only a local valley germ, not global continuation or uniformly small retained momenta.", "",
        f"The largest retained O(p2) generalized mass squared is `{r['derivative_expansion_gate']['retained_to_integrated_min_m2_ratio']:.6g}` times the smallest integrated scalar mass squared; `{r['derivative_expansion_gate']['retained_generalized_modes_above_smallest_integrated_m2']}` retained real modes lie above that smallest gap. Thus these generalized masses must not be promoted to uniformly accurate poles or finite-threshold logarithms.", "",
        f"That global gap test is only sufficient, not a new no-go: actual source support leaves `{r['derivative_expansion_gate']['linearly_coupled_retained_modes']}` linearly coupled modes and `{r['derivative_expansion_gate']['linearly_decoupled_retained_modes']}` decoupled modes. The `{r['derivative_expansion_gate']['coupled_modes_at_or_beyond_Taylor_radius']}` modes at m2 = 0.268823846583 couple weakly to C = 0.261600225 (ratio 1.0276132); for this near-crossing the Taylor series fails while the exact quadratic Schur resolvent remains available. The public `scalar_nonlocal_kernel_api` computes that resolvent rather than guessing a mass correction.", "",
        f"The actual source-selected canonical 2x2 block is invariant under the full Hessian to `{r['derivative_expansion_gate']['near_crossing_exact_block']['full_Hessian_closure_error']:.3e}` and has exact tree poles `{r['derivative_expansion_gate']['near_crossing_exact_block']['exact_tree_poles_m2']}`. Its retained O(p2) mass error is only `{100*r['derivative_expansion_gate']['near_crossing_exact_block']['relative_retained_O_p2_error']:.6g}%`: a resolvent repair, not a physical model obstruction.", "",
        "The action is frozen at tree level. A fixed-VEV one-loop CT affects tree EFT coefficients at one-loop order and must accompany a future loop EFT; inserting that shift inside an already one-loop gauge threshold is a higher-order operation. No upper saddle displacement or physical branch change is made here.", "",
        "The public tensor functions return full local valley, current, metric, directional-vertex tensors and the exact quadratic scalar resolvent. Calling `run(include_tensors=True)` also exports the exact common bases, all generator representations, the potential callable and the same tested full tensor payload. The next matcher must construct the lower covariant gauge-scalar-ghost fluctuation operator and perform a hard-minus-EFT subtraction with the upper Wilson action in the same scheme.", "",
        f"Validation: **{r['summary']['passed']}/{r['summary']['total']}**; Hessian cache: `{r['cache']}`.", "",
        "References: [covariant functional matching](https://arxiv.org/abs/1604.01019), [covariant derivative expansion](https://arxiv.org/abs/1412.1837).", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    report = run()
    OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    OUT.with_suffix(".md").write_text(markdown(report))
    print(json.dumps(report["summary"], indent=2))
    print(json.dumps(report["Higgs_local_vertices"], indent=2))
    print(json.dumps(report["retained_126_neutral_triplet_vertices"], indent=2))
    for check in report["checks"]:
        if not check["pass"]:
            print("FAILED:",check["name"])
    if report["summary"]["passed"] != report["summary"]["total"]:
        raise SystemExit(1)
