#!/usr/bin/env python3
"""Frozen-action scalar poles seen by actual P54 interaction vertices.

Not a new action, flavor fit, loop pole calculation or extra spatial dimension.
The public Gaussian reduction transports the kernel, sources AND contact
response. Cached Hessians are content-addressed by action, parameters and x.
Run with the pinned requirements_p54_matching.txt environment. --cache-dir
may point to a read-only prior checkout; this consumer never writes its cache.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

RF = Path(__file__).resolve().parents[1]
ROOT = RF.parent
OUT = RF / "output/p54_spectral_channels"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def decode(x):
    return np.asarray(x["real"]) + 1j*np.asarray(x["imag"])


def encode(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def sym(x):
    return (x+x.T)/2


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-30))


def image(a, tol=1e-9):
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    return u[:, s > tol]


def range_basis(p):
    d, u = np.linalg.eigh(sym(p))
    return u[:, d > .5]


def clusters(values, tol=1e-9):
    groups = []
    for i, value in enumerate(values):
        if not groups or abs(value-values[groups[-1][0]]) > tol:
            groups.append([i])
        else:
            groups[-1].append(i)
    return groups


def eliminate(kernel, sources, contact, keep):
    """Exact Gaussian response map (K,J,C), analytic in complex momentum.

    Response is C + J.T solve(K,J). The transpose is NOT a dagger: the
    original real Hermitian-current sources are analytically continued.
    No inverse of a singular eliminated block is silently regularized.
    """
    keep = np.asarray(keep, dtype=int)
    if len(np.unique(keep)) != len(keep):
        raise ValueError("Duplicate retained coordinates")
    drop = np.setdiff1d(np.arange(len(kernel)), keep)
    a = kernel[np.ix_(keep, keep)]
    b = kernel[np.ix_(keep, drop)]
    c = kernel[np.ix_(drop, drop)]
    jr, jh = sources[keep], sources[drop]
    solved = np.linalg.solve(c, np.column_stack((b.T, jh)))
    nr = len(keep)
    return (sym(a-b@solved[:, :nr]), jr-b@solved[:, nr:],
            sym(contact+jh.T@solved[:, nr:]))


def response(kernel, sources, contact):
    return sym(contact+sources.T@np.linalg.solve(kernel, sources))


def run(cache_dir):
    inputs = [RF/"code/verify_p54_common_yukawa_phase.py",
              RF/"code/verify_p54_spinor_intertwiners.py",
              RF/"code/verify_p54_p1_hessian_spectrum.py",
              RF/"code/verify_p54_p2_two_site_matching.py",
              RF/"code/verify_p54_ps_finite_thresholds.py",
              RF/"output/p54_full_doublet_cw.json",
              RF/"output/p54_common_yukawa_phase.json",
              RF/"output/p54_ps_finite_thresholds.json",
              RF/"output/p54_scalar_chn.json",
              RF/"output/p54_p2_two_site_matching.json"]
    common = module("spectral_common", inputs[0])
    helper = module("spectral_helper", inputs[3])
    upper = module("spectral_upper", inputs[4])
    full, phase, finite, chn, p2 = [json.loads(p.read_text()) for p in inputs[5:]]
    geom = common.build_geometry()
    p1 = geom["p1"]
    pars, vev = full["tree_parameters"], full["vacuum"]
    x0 = p1.vacuum_vector(vev["omega"], vev["sigma"], vev["vs"])
    keybase = hashlib.sha256(inputs[2].read_bytes()+json.dumps(pars, sort_keys=True).encode()).digest()
    cache_ledger = []

    def hessian(x):
        key = hashlib.sha256(keybase+np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = cache_dir/(key+".npz")
        if not path.exists():
            raise FileNotFoundError("Missing same-action Hessian; run the existing P54 cache producer: "+str(path))
        with np.load(path) as data:
            result = np.asarray(data["h"])
        cache_ledger.append({"key": key, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        return sym(result)

    h0 = hessian(x0)
    checks = []

    def check(name, residual, tolerance=2e-9):
        checks.append({"name": name, "residual": float(residual), "tolerance": tolerance,
                       "passed": bool(np.isfinite(residual) and residual < tolerance)})

    # Gauge directions are not particles. The physical PQ mode is retained.
    gauge = image(p1.gauge_orbit(x0))
    phys = range_basis(np.eye(328)-gauge@gauge.T)
    physical_projector = phys@phys.T
    kphys = sym(phys.T@h0@phys)
    masses, u = np.linalg.eigh(kphys)
    modes = phys@u
    check("33 gauge directions, 295 physical real scalar directions", abs(gauge.shape[1]-33)+abs(phys.shape[1]-295))
    check("tree gauge Ward identity", np.linalg.norm(h0@gauge))
    check("physical slice is invariant under the same Hessian", np.linalg.norm(h0@phys-phys@kphys))
    check("no negative physical tree eigenvalues", max(0., -masses.min()))
    check("four tree Higgs zeros plus one physical PQ zero retained", abs(np.sum(abs(masses)<1e-8)-5))

    sm = p1.sm_representation_matrices()
    sm_all = [r for rs in sm.values() for r in rs]
    casimir = -sum(r@r for r in sm_all)
    ce, cu = np.linalg.eigh(sym(casimir))
    bs = cu[:, abs(ce)<1e-9]
    psing = bs@bs.T
    check("five real SM singlets before gauge removal", abs(bs.shape[1]-5))
    check("unbroken SM commutes with physical scalar Hessian",
          max(np.linalg.norm(h0@r-r@h0) for r in sm_all))

    # The exact TREE light vector, not the improved mixed-order ray.
    bd = geom["conjugation_diagonal"][:, None]*decode(phase["scalar_copy_basis_physical_328x4"])
    cr = decode(full["basis"]["tree_light_coefficients"])
    qr, qi = math.sqrt(2)*bd.real, -math.sqrt(2)*bd.imag
    real_c = np.r_[cr.real, cr.imag]
    doublet = np.column_stack((qr, qi))
    q = doublet@real_c
    jets_r = [(hessian(x0+.02*a)-hessian(x0-.02*a))/.04 for a in qr.T]
    phase_transport = expm(-np.pi*sm["Y"][0])
    jets = np.asarray(jets_r+[phase_transport@a@phase_transport.T for a in jets_r])
    dh = np.einsum("a,aij->ij", real_c, jets)
    jqq = dh@q
    jrho = psing@jqq
    check("common complex tree light direction solves Hq=0", np.linalg.norm(h0@q), 2e-11)
    check("hypercharge transport reconstructs imaginary-copy jets", np.linalg.norm(phase_transport@qr-qi), 2e-11)
    # Average the four real Higgs components using exact SM covariance of V3.
    transports = [np.eye(328)]+[expm(np.pi*r) for r in sm["SU2"]]
    quartet = np.column_stack([r@q for r in transports])
    check("four transported Higgs components are canonical", np.linalg.norm(quartet.T@quartet-np.eye(4)))
    averaged = sum(r@jqq for r in transports)/4
    check("HdaggerH source equals full doublet trace of cubic tensor", np.linalg.norm(averaged-jrho))
    check("HdaggerH source is an SM singlet", max(np.linalg.norm(r@jrho) for r in sm_all))
    check("HdaggerH source has no massless scalar component", np.linalg.norm(modes[:, abs(masses)<1e-8].T@jrho))

    names = ["HdaggerH"]
    columns = [jrho]
    channels = [{"name": names[0], "definition": "J_A=(1/4) sum_alpha V3[e_A,q_alpha,q_alpha]", "unit": "omega", "family_factor": "none"}]

    # True Yukawa vertices of declared Spin(10) spinor states. h and f are
    # kept symbolic: the unit-spurion coefficients are not a family fit.
    zh, zf = common.real_scalar_yukawa_tensors(geom, old_coordinates=True)
    for group in sm:
        for index, (rs, rf) in enumerate(zip(sm[group], geom["spin_sm"][group])):
            for label, tensor in (("h", zh), ("f", zf)):
                ward = (np.einsum("ik,akj->aij", rf.T, tensor)
                        +np.einsum("aik,kj->aij", tensor, rf)
                        +np.einsum("ba,bij->aij", rs, tensor))
                check(f"full {label} Yukawa {group}[{index}] intertwiner Ward identity", np.linalg.norm(ward))
    states = {name: decode(value) for name, value in phase["fermion_states_16"].items()}
    pairs = [("uLuc", "uL", "uc"), ("dLdc", "dL", "dc"),
             ("eLec", "eL", "ec"), ("nuLN", "nuL", "nuc"),
             ("NN", "nuc", "nuc"), ("nuLnuL", "nuL", "nuL"),
             ("ucdc", "uc", "dc")]
    raw_yukawa = {}
    selection = []
    for pair, left, right in pairs:
        l, r = states[left][:, 0], states[right][:, 0]
        yf = 1j*geom["spin_sm"]["Y"][0]
        charge = float(np.real(l.conj()@yf@l+r.conj()@yf@r))
        for label, tensor in (("h", zh), ("f", zf)):
            source = np.einsum("i,aij,j->a", l, tensor, r)
            raw_yukawa[pair+"_"+label] = source
            ward = np.linalg.norm(1j*sm["Y"][0]@source-charge*source)
            check(pair+" "+label+" hypercharge selection rule", ward)
            selection.append({"pair": pair, "spurion": label, "fermion_pair_Y": charge,
                              "raw_vertex_norm": float(np.linalg.norm(source)),
                              "hypercharge_Ward_residual": float(ward)})
            for component, value in (("Re", source.real), ("Im", source.imag)):
                name = pair+"_"+label+"_"+component
                names.append(name)
                columns.append(value)
                channels.append({"name": name, "definition": "Hermitian-current component of l^T Z_"+label+"[A] r", "unit": "dimensionless", "family_factor": label+"_raw; unfitted"})
    check("10_H has no NN Yukawa vertex", np.linalg.norm(raw_yukawa["NN_h"]), 2e-12)
    check("actual NN scalar coupling reproduces previous common-phase tensor",
          np.linalg.norm(raw_yukawa["NN_f"]-decode(chn["raw_NN_scalar_couplings_328"])), 2e-12)

    # Scalar-vector-vector vertices from |D X|^2/2. Average over each full
    # degenerate gauge-mass eigenspace; no preferred vector component.
    generators = [helper.so_generator(a, b) for a in range(10) for b in range(a+1, 10)]
    reps = np.asarray([p1.representation_matrix(r) for r in generators])
    orbit = np.column_stack([r@x0 for r in reps])
    check("canonical generator orbit reproduces P1 orbit", np.linalg.norm(orbit-p1.gauge_orbit(x0)))
    g = float(p2["gauge_coupling_iteration"][-1]["g_output"])
    check("historical gauge coupling agrees with declared muU=gU omega card", abs(g-float(full["scheme"]["mu_over_omega"])), 2e-12)
    mv = sym(g*g*orbit.T@orbit)
    dv, uv = np.linalg.eigh(mv)
    vector_rows = []
    for group in clusters(dv):
        mass = float(np.mean(dv[group]))
        vb = uv[:, group]
        rp = np.einsum("ac,aij->cij", vb, reps)
        vertex = 2*g*g*sum(r.T@(r@x0) for r in rp)/len(group)
        check("vector group SM-singlet vertex m2="+str(round(mass, 6)), max(np.linalg.norm(r@vertex) for r in sm_all))
        if mass < 1e-9:
            check("unbroken massless vectors have no tree single-scalar mass vertex", np.linalg.norm(vertex))
            continue
        name = "VV_"+str(len(vector_rows))
        # Independent derivative of trace(P_alpha M_V^2(x))/d_alpha.
        direction = jrho/np.linalg.norm(jrho)
        step = 1e-4
        def mean_mass(x):
            ox = np.column_stack([r@x for r in reps])@vb
            return g*g*np.sum(ox*ox)/len(group)
        fd = (mean_mass(x0+step*direction)-mean_mass(x0-step*direction))/(2*step)
        check(name+" kinetic-term derivative normalization", abs(fd-vertex@direction), 2e-10)
        names.append(name)
        columns.append(vertex)
        channels.append({"name": name, "definition": "(1/d) trace(P_vector d M_V^2/d x_A); vertex multiplies metric polarization contraction", "unit": "omega", "family_factor": "none"})
        vector_rows.append({"channel": name, "vector_mass_squared_over_omega2": mass,
                            "real_multiplicity": len(group), "scalar_vertex_norm_over_omega": float(np.linalg.norm(vertex))})

    raw_j = np.column_stack(columns)
    jphys = phys.T@raw_j
    jfull = phys@jphys
    amplitudes = modes.T@raw_j
    residue_rows = []
    for group in clusters(masses):
        a = amplitudes[group]
        residue = sym(a.T@a)
        value = float(np.mean(masses[group]))
        # Do not hide the precise raw eigenvalue residual in the zero cluster.
        value_display = 0. if abs(value)<1e-9 else value
        check("residue positivity at m2="+str(round(value_display, 6)), max(0., -np.linalg.eigvalsh(residue).min()))
        residue_rows.append({"mass_squared_over_omega2": value_display,
                             "raw_mean_eigenvalue": value, "multiplicity": len(group),
                             "coupling_squared_diagonal": np.diag(residue).tolist(),
                             "residue_matrix": residue.tolist()})
    check("zeroth spectral sum rule includes massless modes", rel(sum(np.asarray(r["residue_matrix"]) for r in residue_rows), jphys.T@jphys))
    check("first spectral moment equals J^T H J", rel(sum(r["raw_mean_eigenvalue"]*np.asarray(r["residue_matrix"]) for r in residue_rows), jphys.T@kphys@jphys))

    # Physical axion vs eaten Goldstone: both are zero modes of the raw H,
    # but only the former survives the physical scalar slice.
    axion = physical_projector@p1.pq_direction(x0)
    axion /= np.linalg.norm(axion)
    gnn = raw_yukawa["NN_f"]
    axion_nn = complex(gnn@axion)
    check("physical PQ vector is a retained tree zero mode", np.linalg.norm(h0@axion))
    # The first attempted positive NN-PQ coupling hypothesis FAILED. In this
    # vacuum Sigma's phase is eaten; the remaining PQ tangent is the S phase.
    # Verify its representation support before recording the zero selection.
    singlet_field = np.zeros((328, 328))
    singlet_field[p1.SL_S, p1.SL_S] = np.eye(2)
    check("physical PQ tangent lies entirely in gauge-singlet S", np.linalg.norm(axion-singlet_field@axion))
    check("S has no renormalizable Yukawa tensor", np.linalg.norm(zh[p1.SL_S])+np.linalg.norm(zf[p1.SL_S]))
    check("NN physical PQ coupling vanishes by representation support", abs(axion_nn))
    check("HdaggerH has no physical PQ pole residue", abs(jrho@axion))
    check("all declared single-scalar probe vertices miss the physical PQ mode", np.linalg.norm(axion@raw_j))
    massive = masses>1e-8
    chn_value = .5*(gnn@modes[:, massive])@((modes[:, massive].T@jrho)/masses[massive])
    prior_chn = complex(decode(chn["rays"][0]["CHN_raw_times_omega"]))
    check("zero-momentum massive spectral moment reproduces actual tree CHN", abs(chn_value-prior_chn), 2e-11)

    # Reconstruct the declared upper 55-scalar split; keep PQ explicitly.
    parents = helper.parent_projectors(p1, helper.casimir_operators(p1, helper.ps_generators(p1)))
    pa = sum(row["projector"] for row in parents if row["label"] in upper.ACTIVE)
    vu = finite["upper_PS"]["vacuum"]
    xu = p1.vacuum_vector(vu["omega"], vu["sigma"], vu["vs"])
    gu = image(p1.gauge_orbit(xu))
    pqu = (np.eye(328)-gu@gu.T)@p1.pq_direction(xu)
    pqu /= np.linalg.norm(pqu)
    bh = range_basis(np.eye(328)-pa-gu@gu.T-np.outer(pqu, pqu))
    check("55 upper heavy scalars lie in actual physical slice", abs(bh.shape[1]-55)+np.linalg.norm(gauge.T@bh))
    # Final coordinate chart: actual phi10 bidoublet parent + physical PQ.
    # This is a factorization chart, NOT a new low-energy EFT field census.
    p10 = next(row["projector"] for row in parents if row["label"]=="phi10:(1,2,2)")
    bl = image(np.column_stack((range_basis(p10), axion)))
    bm = range_basis(physical_projector-bh@bh.T-bl@bl.T)
    basis = np.column_stack((bl, bm, bh))
    nl, nm, nh = bl.shape[1], bm.shape[1], bh.shape[1]
    check("nested physical factorization is complete and orthonormal", np.linalg.norm(basis.T@basis-np.eye(295))+np.linalg.norm(basis@basis.T-physical_projector))
    ktree = sym(basis.T@h0@basis)
    source = basis.T@raw_j
    reduction_rows = []
    points = [.001, .03, .2, complex(-.02, .003), complex(-.265, .001)]
    for t in points:
        k = ktree+t*np.eye(295)
        contact = np.zeros((len(names), len(names)))
        original = response(k, source, contact)
        first = eliminate(k, source, contact, np.arange(nl+nm))
        staged = eliminate(*first, np.arange(nl))
        direct = eliminate(k, source, contact, np.arange(nl))
        # Reverse order: remove the middle block first, then the upper one.
        alternate_first = eliminate(k, source, contact, np.r_[np.arange(nl), np.arange(nl+nm, 295)])
        alternate = eliminate(*alternate_first, np.arange(nl))
        exact = response(*staged)
        missing_contact = staged[1].T@np.linalg.solve(staged[0], staged[1])
        check("HH final linear source vanishes by SM/PQ selection at pE2="+str(t), np.linalg.norm(staged[1][:, 0]))
        hh_missing_error = abs(original[0, 0]-missing_contact[0, 0])/abs(original[0, 0])
        check("negative control loses full HH response at pE2="+str(t), abs(hh_missing_error-1.))
        spectral = amplitudes.T@((1/(masses+t))[:, None]*amplitudes)
        errors = {"full_vs_staged_response": rel(original, exact),
                  "direct_vs_staged_kernel": rel(direct[0], staged[0]),
                  "direct_vs_staged_sources": rel(direct[1], staged[1]),
                  "direct_vs_staged_contact": rel(direct[2], staged[2]),
                  "opposite_order_response": rel(original, response(*alternate)),
                  "spectral_vs_full_response": rel(original, spectral)}
        for name, error in errors.items():
            check(name+" at pE2="+str(t), error, 5e-8)
        # A passive non-orthonormal coordinate change must transform both
        # metric/kernel and vertices; it is not physical object rotation.
        scale = np.linspace(.8, 1.2, 295)
        transformed = response(k*scale[:, None]*scale[None, :], source*scale[:, None], contact)
        check("noncanonical coordinate covariance at pE2="+str(t), rel(original, transformed), 5e-8)
        reduction_rows.append({"pE2_over_omega2": encode(t), "errors": errors,
                               "relative_error_if_induced_contact_omitted": rel(original, missing_contact),
                               "HH_relative_error_if_induced_contact_omitted": float(hh_missing_error),
                               "final_kernel_condition_number": float(np.linalg.cond(staged[0]))})
    check("negative control detects lost response without induced contact",
          0 if max(r["relative_error_if_induced_contact_omitted"] for r in reduction_rows)>.01 else 1)

    # A change inside a degenerate eigenspace cannot change its residue.
    zero_group = next(group for group in clusters(masses) if abs(masses[group[0]])<1e-9)
    a = amplitudes[zero_group]
    swap = np.eye(len(zero_group))[::-1]
    check("degenerate pole residue independent of eigenvector labels", rel(a.T@a, (swap@a).T@(swap@a)))
    # Source rank is diagnostic only; the physical vertices were fixed first.
    singular = np.linalg.svd(jphys, compute_uv=False)
    source_rank = int(np.sum(singular>1e-9))
    observable_basis = image(np.column_stack([modes[:, group]@(modes[:, group].T@raw_j)
                                              for group in clusters(masses)]))
    check("channel-generated spectral subspace is H-invariant",
          np.linalg.norm((np.eye(328)-observable_basis@observable_basis.T)@h0@observable_basis), 2e-8)
    # Minimal observable realization reconstructed from the actual sources.
    ko = sym(observable_basis.T@h0@observable_basis)
    jo = observable_basis.T@raw_j
    for t in points:
        check("minimal interaction-generated realization at pE2="+str(t),
              rel(response(kphys+t*np.eye(295), jphys, np.zeros((len(names), len(names)))),
                  response(ko+t*np.eye(len(ko)), jo, np.zeros((len(names), len(names))))), 5e-8)

    report = {"schema": "p54-physical-scalar-spectral-channels-v1", "date": "2026-09-10",
              "all_checks_pass": all(r["passed"] for r in checks),
              "checks_passed": sum(r["passed"] for r in checks), "checks_total": len(checks), "checks": checks,
              "scope": {"action": "P54PQ-v2 unchanged; frozen tree polynomial and stationary background",
                        "kinematic_object": "tree scalar-exchange form factors; no full cross sections or loop pole masses",
                        "family_inputs": "h_raw and f_raw remain symbolic, correlated UV matrices; unit-spurion geometric vertices only",
                        "spacetime": "3+1 unchanged; layer labels denote retained field-coordinate descriptions",
                        "gauge": "33 eaten directions excluded from physical scalar response; PQ axion retained",
                        "units": "omega=1; J_HH and J_VV in omega, Yukawa vertices dimensionless; residue entries carry corresponding product units",
                        "source_basis": "scalar action and Yukawa/gauge vertices fixed before any source-rank decomposition; quark pairs use first matched colour in common-phase dictionary, not a colour-averaged bilinear",
                        "IR": "axion is massless only in this perturbative tree action; no QCD potential included",
                        "not_completed": ["full finite gauge/ghost/scalar matching", "fermion-inclusive fitted light state", "physical flavor fit", "momentum-dependent loop poles and widths", "extra-dimensional completion"]},
              "vacuum": vev, "gauge_coupling": g,
              "channels": channels, "channel_order": names,
              "raw_interaction_sources_328": raw_j.tolist(),
              "physical_projected_sources_328": jfull.tolist(),
              "channel_source_rank": source_rank, "channel_generated_spectral_subspace_dimension": observable_basis.shape[1],
              "physical_scalar_dimension": len(masses), "physical_mass_squared_over_omega2": masses.tolist(),
              "spectral_clusters": residue_rows, "yukawa_selection_rules": selection,
              "vector_pair_channels": vector_rows,
              "physical_axion": {"NN_raw_coupling": encode(axion_nn), "NN_raw_residue": abs(axion_nn)**2,
                                  "HdaggerH_coupling": float(jrho@axion),
                                  "all_declared_vertices_norm": float(np.linalg.norm(axion@raw_j)),
                                  "interpretation": "Physical PQ tangent is the S phase; all declared linear probe residues vanish. It remains a physical zero mode, not a gauge artifact or absent state.",
                                  "rejected_initial_hypothesis": "NN would see a nonzero physical PQ scalar pole. The first run returned zero; support in S, absent from both Yukawa tensors, explains the failure without parameter changes."},
              "CHN_spectral_moment": {"omega_CHN_over_f_raw": encode(chn_value),
                                       "previous_same_action_coefficient": encode(prior_chn)},
              "nested_elimination": {"dimensions_L_M_U": [nl, nm, nh],
                                     "chart": "phi10 bidoublet plus physical PQ | other retained physical parents | upper 55",
                                     "low_energy_EFT_census_claimed": False, "points": reduction_rows},
              "cache_recomputed": False, "cache_reads": cache_ledger,
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in inputs+[Path(__file__).resolve()]},
              "references": ["https://arxiv.org/abs/2105.02058", "https://arxiv.org/abs/1604.01019"]}
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    write_markdown(report)
    write_tex_data(report)
    print(json.dumps({"passed": report["checks_passed"], "total": len(checks),
                      "source_rank": source_rank, "observable_dimension": observable_basis.shape[1],
                      "dimensions": [nl, nm, nh], "axion": report["physical_axion"],
                      "failed": [r for r in checks if not r["passed"]]}, indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit(1)
    return report


def write_markdown(r):
    names = r["channel_order"]
    hh = names.index("HdaggerH")
    nr, ni = names.index("NN_f_Re"), names.index("NN_f_Im")
    lines = ["# P54 common scalar spectrum and interaction-defined channels", "",
             "Date: 2026-09-10. Action and background unchanged. This is a tree scalar-exchange audit, not a fitted particle spectrum or loop-complete scattering prediction.", "",
             f"Validation: **{r['checks_passed']}/{r['checks_total']}**. No cached Hessian was recomputed.", "",
             "## Physical vertices, not arbitrary observation projectors", "",
             "The HdaggerH vertex is the singlet trace of the cubic tensor on the exact tree light doublet. Seven named fermion bilinears use actual Clifford/Yukawa tensors; h_raw and f_raw are factored out, not fitted. Massive vector-pair vertices are derivatives of the actual kinetic-term mass matrix, averaged over complete degenerate vector eigenspaces.", "",
             f"There are {len(names)} real current slots (including selection-rule zeros), source rank {r['channel_source_rank']}, and a {r['channel_generated_spectral_subspace_dimension']}-dimensional channel-generated invariant subspace inside all 295 physical real scalar coordinates. The 33 gauge directions are excluded; all five physical tree zero modes remain.", "",
             r"For physical scalar H and interaction vertices J, $$F_E(t)=J^T(H+tI)^{-1}J=\sum_\lambda R_\lambda/(t+\lambda),\quad R_\lambda=J^T\Pi_\lambda J\succeq0.$$", "",
             "Residues are aggregated over exact degenerate eigenspaces, not assigned to arbitrary eigenvectors. The JSON contains every pole cluster and its full cross-channel residue matrix. Diagonal entries are squared vertex strengths, not branching fractions.", "",
             "## Singlet-channel poles and squared couplings", "",
             "HH weights are in omega^2; NN weights have the family spurion removed and are dimensionless. Other channels can see additional poles listed in JSON.", "",
             "| m^2 / omega^2 | HH weight | NN weight per unit f_raw |", "|---:|---:|---:|"]
    for row in r["spectral_clusters"]:
        diag = row["coupling_squared_diagonal"]
        if diag[hh]+diag[nr]+diag[ni]>1e-12:
            lines.append(f"| {row['mass_squared_over_omega2']:.9g} | {diag[hh]:.9g} | {diag[nr]+diag[ni]:.9g} |")
    a = r["physical_axion"]
    lines += ["", "The physical PQ mode lies in the S phase, which has no direct Yukawa tensor. Its residues vanish in ALL declared linear channels, including NN and HdaggerH. An initial nonzero-NN-PQ hypothesis was rejected by the actual vertices, not repaired by fitting. This is a physical mode invisible to this current set, not a deleted state, a gauge artifact, a globally decoupled axion, or a claim of zero QCD axion mass. The massive inverse spectral moment independently reproduces the previous nonzero CHN coefficient, including its common complex phase.", "",
              f"The symmetry-enforced zero NN axion residue has numerical residual `{a['NN_raw_residue']:.3g}`; this is not a predicted tiny nonzero coupling. The massive moment gives `omega CHN / f_raw = {complex(decode(r['CHN_spectral_moment']['omega_CHN_over_f_raw'])).imag:.12g} i`, with real part consistent with zero.", "",
              "## Exact nested elimination must also transport vertices and contact terms", "",
              r"For retained r and eliminated h, $$K'=K_{rr}-K_{rh}K_{hh}^{-1}K_{hr},\quad J'=J_r-K_{rh}K_{hh}^{-1}J_h,\quad C'=C+J_h^TK_{hh}^{-1}J_h.$$", "",
              r"The full response is $$F=C+J^TK^{-1}J=C'+J'^TK'^{-1}J'.$$", "",
              f"The actual nested dimensions are `{r['nested_elimination']['dimensions_L_M_U']}`. The final chart is a factorization device, not an asserted SM-only EFT. Direct, staged, and opposite-order elimination agree at spacelike and complex near-pole momenta, and under noncanonical coordinate rescaling.", "",
              "| p_E^2 / omega^2 | full vs staged relative residual | error if induced contact omitted |", "|---:|---:|---:|"]
    for row in r["nested_elimination"]["points"]:
        lines.append(f"| {complex(decode(row['pE2_over_omega2'])):.5g} | {row['errors']['full_vs_staged_response']:.3g} | {row['relative_error_if_induced_contact_omitted']:.3g} |")
    lines += ["", "The induced contact is nonlocal before a derivative expansion and can carry poles. Omitting it loses 100% of the HH-to-HH scalar-exchange response at every tested point. The matrix-norm errors above use the declared unit-spurion current normalization, not fitted cross sections. Keeping K alone is not an observable-preserving projection. Feshbach-Schur transitivity is existing mathematics, not evidence of a fourth spatial direction.", "",
              "## Scope and next gates", "",
              "- Preserve the action and symbolic UV family matrices. The selected currents do not exhaust all physical channels; darkness is always relative to the declared current set.",
              "- Use these source/contact identities in the common finite matching workflow. Full vector/Goldstone/ghost/scalar terms, running, widths, and a physical flavor fit remain open.",
              "- Extra spatial dimensions are comparison-only. A fixed interval spectrum would obey m_n^2=M_5^2+(n*pi/L)^2 and therefore (m_2^2-m_0^2)/(m_1^2-m_0^2)=4. No extra-dimensional parameters or particle assignments are fitted here.", "",
              "## Reproduction", "",
              "Use the pinned `code/requirements_p54_matching.txt` environment and the existing full-doublet/kinetic Hessian cache. This consumer is read-only on the cache. The exact keys and file hashes are recorded in JSON; an absent or incompatible cache fails explicitly.", "",
              "```sh", "python3 route_f/code/verify_p54_spectral_channels.py --cache-dir /path/to/tmp/p54_full_doublet_cw", "```", "",
              "Derivation: [TeX](../tex/p54_spectral_channels_nested_elimination.tex), [PDF](pdf/p54_spectral_channels_nested_elimination.pdf).", "",
              "References: [Feshbach-Schur map](https://arxiv.org/abs/2105.02058), [nonlocal covariant EFT matching](https://arxiv.org/abs/1604.01019).", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))


def write_tex_data(r):
    """Generated numerical companion; all derivations live in the main TeX."""
    names = r["channel_order"]
    nr, ni = names.index("NN_f_Re"), names.index("NN_f_Im")
    lines = ["% Generated by verify_p54_spectral_channels.py; do not hand edit.",
             rf"\newcommand{{\ChecksPassed}}{{{r['checks_passed']}}}",
             rf"\newcommand{{\ChecksTotal}}{{{r['checks_total']}}}"]
    rows = []
    def tex_number(value):
        text = f"{value:.8g}"
        if "e" in text:
            base, exponent = text.split("e")
            return rf"\ensuremath{{{base}\times10^{{{int(exponent)}}}}}"
        return text
    for row in r["spectral_clusters"]:
        d = row["coupling_squared_diagonal"]
        if d[0]+d[nr]+d[ni]>1e-12:
            rows.append(f"{row['mass_squared_over_omega2']:.8f} & {tex_number(d[0])} & {tex_number(d[nr]+d[ni])}"+r" \\")
    lines += [r"\newcommand{\SingletPoleRows}{", *rows, "}"]
    OUT.with_name(OUT.name+"_tables.tex").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT/"tmp/p54_full_doublet_cw")
    args = parser.parse_args()
    run(args.cache_dir)
