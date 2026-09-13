#!/usr/bin/env python3
"""Source-preserving one-loop jets on the actual upper-PS P54 background.

Tests the computed vector vertex/kinetic contribution, including induced C.
The other hard diagrams are explicitly UNKNOWN; their absent arrays here
mean a labelled additive sub-contribution, never a full physical matching.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import verify_p54_spectral_channels as spectral
import verify_p54_upper_vector_matching as vector_matching
import verify_p54_upper_yukawa_thresholds as yuk

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_kjc_finite_matching"


def reduce_jet(tree, loop, keep):
    """Frechet derivative of Gaussian (K,J,C) elimination; no resummation.

    K and C are complex symmetric analytic kernels, NOT Hermitian at
    complex momentum. A singular eliminated block is deliberately rejected.
    """
    k, j, c = tree
    dk, dj, dc = loop
    keep = np.asarray(keep, dtype=int)
    if len(np.unique(keep)) != len(keep) or np.any(keep < 0) or np.any(keep >= len(k)):
        raise ValueError("invalid retained coordinates")
    if any(np.shape(x) != np.shape(y) for x,y in zip(tree,loop)):
        raise ValueError("tree/loop shape mismatch")
    for a in (k,c,dk,dc):
        if yuk.error(a,a.T) > 1e-10:
            raise ValueError("analytic real-field kernel requires transpose symmetry")
    drop = np.setdiff1d(np.arange(len(k)),keep)
    b, d = k[np.ix_(keep,drop)], k[np.ix_(drop,drop)]
    db, dd = dk[np.ix_(keep,drop)], dk[np.ix_(drop,drop)]
    jh, djh = j[drop],dj[drop]
    # q = D^-1 B^T, w = D^-1 J_h; transposes use symmetry of D.
    q, w = np.linalg.solve(d,b.T), np.linalg.solve(d,jh)
    dkr = dk[np.ix_(keep,keep)]-db@q-q.T@db.T+q.T@dd@q
    djr = dj[keep]-db@w-q.T@djh+q.T@dd@w
    dcr = dc+djh.T@w+w.T@djh-w.T@dd@w
    return spectral.eliminate(k,j,c,keep), (dkr,djr,dcr)


def response_jet(tree, loop):
    k,j,c = tree
    dk,dj,dc = loop
    w = np.linalg.solve(k,j)
    return spectral.response(k,j,c), dc+dj.T@w+w.T@dj-w.T@dk@w


def reparameterize_jet(tree,loop,s):
    """phi_old=(I+ell*s) phi_new. Sources transform with the kernel."""
    k,j,c = tree
    dk,dj,dc = loop
    return tree,(dk+s.T@k+k@s,dj+s.T@j,dc)


def merge_contributions(contributions):
    """Add disjoint computed diagram jets only at identical matching context.

    This operation does not certify that the diagram list is complete.
    Background, scale, momentum, basis and truncation are part of context.
    """
    if not contributions:
        raise ValueError("No computed contributions")
    required = ("action","subtraction","gauge","tadpoles","background",
                "mu","momentum","basis","loop_order","mass_expansion")
    context = contributions[0]["context"]
    if any(key not in context for key in required):
        raise ValueError("Missing matching context")
    seen, total = set(), None
    for piece in contributions:
        if piece["context"] != context:
            raise ValueError("Cannot combine different matching prescriptions/backgrounds")
        diagrams = piece["diagrams"]
        if not diagrams or len(set(diagrams)) != len(diagrams) or seen.intersection(diagrams):
            raise ValueError("Duplicated or undeclared diagrams")
        seen.update(diagrams)
        if len(piece["loop"]) != 3 or any(x is None for x in piece["loop"]):
            raise ValueError("A computed contribution must specify all K,J,C components")
        if total is None:
            total = tuple(np.array(x,copy=True) for x in piece["loop"])
        else:
            if any(x.shape != y.shape for x,y in zip(total,piece["loop"])):
                raise ValueError("Mismatched source/field dimensions")
            total = tuple(x+y for x,y in zip(total,piece["loop"]))
    return total


def source_columns(tensors,geometry):
    """P-SPEC1 pairs plus two different-colour diquark probes in PS chart.

    Family spurions stay factored out. Re/Im are Hermitian current slots.
    """
    phase = json.loads((RF/"output/p54_common_yukawa_phase.json").read_text())
    ub = geometry["ps"]["fermion_basis"]
    states = {key:ub.conj().T@(np.asarray(v["real"])+1j*np.asarray(v["imag"]))
              for key,v in phase["fermion_states_16"].items()}
    pairs = [("uLuc","uL","uc"),("dLdc","dL","dc"),("eLec","eL","ec"),
             ("nuLN","nuL","nuc"),("NN","nuc","nuc"),
             ("nuLnuL","nuL","nuL"),("ucdc","uc","dc")]
    indexed = [(name,left,0,right,0) for name,left,right in pairs]
    indexed += [("ucdc_crosscolour","uc",0,"dc",1),
                ("uLdL_crosscolour","uL",0,"dL",1)]
    cols,names = [],[]
    for name,left,il,right,ir in indexed:
        for spurion,tensor in tensors.items():
            z = np.einsum("i,aij,j->a",states[left][:,il],tensor,states[right][:,ir])
            for component,value in (("Re",z.real),("Im",z.imag)):
                cols.append(value)
                names.append(name+"_"+spurion+"_"+component)
    return np.column_stack(cols), names


def run(cache_dir):
    geometry = yuk.build_upper_geometry(cache_dir)
    vector = vector_matching.build_vector_geometry(geometry)
    active,heavy,pq = geometry["active"],geometry["heavy"],geometry["pq"]
    basis = np.column_stack((active,pq,heavy)) # 248+1+55 = 304
    h = basis.T@geometry["hessian"]@basis
    kp,ks = vector_matching.kinetic(geometry,vector,geometry["mu"])
    dkin = np.zeros((304,304)); dkin[:248,:248] = ks
    ys, dys = {}, {}
    for s in ("h","f"):
        raw = geometry["all_"+s]
        tree_y = np.einsum("ab,aij->bij",active,raw)
        finite = vector_matching.threshold(tree_y,geometry,vector,geometry["mu"])
        # Fermion legs apply also to heavy-source endpoints. Scalar legs
        # are NOT inserted here: K already carries their normalization.
        dy_raw = -.5*(kp.T@raw+raw@kp)+np.einsum("ab,bij->aij",active,finite["vertex"])
        ys[s] = np.einsum("ab,aij->bij",basis,raw)
        dys[s] = np.einsum("ab,aij->bij",basis,dy_raw)
    j,names = source_columns(ys,geometry)
    dj,_ = source_columns(dys,geometry)
    c = np.zeros((len(names),len(names)))
    # dc=0 ONLY for this specified external-leg/active-vertex/active-kinetic
    # additive contribution before elimination. Unknown 1PI boxes etc are
    # not inferred to vanish. Induced delta C below is computed, not omitted.
    keep_final = np.r_[np.arange(128),248] # H,F,PQ, full PS parents
    keep_upper = np.arange(249) # integrate only the actual upper 55
    keep_reverse = np.r_[keep_final,np.arange(249,304)]
    checks, rows = [], []

    def check(name,a,b=0.,tol=3e-10):
        err = yuk.error(np.asarray(a),np.asarray(b))
        checks.append({"name":name,"residual":err,"tolerance":tol,"pass":err < tol})

    check("304_physical_upper_chart_orthonormal",basis.T@basis,np.eye(304))
    check("chart_excludes_24_eaten_upper_directions",geometry["goldstone"].T@basis)
    check("PQ_source_zero",j[248])
    check("PQ_loop_source_zero_in_computed_subset",dj[248])
    check("different_colour_probes_detect_upper_scalar_exchange",float(np.linalg.norm(j[249:])<1e-8))
    for z in (.013, .27, -.11+.07j, .03+.12j):
        tree = (h+z*np.eye(304),j,c)
        context = {"action":"P54PQ-v2","subtraction":"MSbar-DR",
            "gauge":"background-field Landau","tadpoles":"fixed-VEV",
            "background":geometry["upper"]["vacuum"],"mu":geometry["mu"],
            "momentum":[float(np.real(z)),float(np.imag(z))],
            "basis": "upper-HFLR-PQ-heavy/"+hashlib.sha256(j.tobytes()).hexdigest(),
            "loop_order":1,"mass_expansion":"leading retained m2/MV2; exact hard scalar masses"}
        pieces = [{"context":context,"diagrams":["active_scalar_vector_kinetic"],
                   "loop":(z*dkin,np.zeros_like(j),np.zeros_like(c))},
                  {"context":context,"diagrams":["active_vector_Y_vertex","vector_fermion_external_legs"],
                   "loop":(np.zeros_like(h),dj,np.zeros_like(c))}]
        loop = merge_contributions(pieces)
        check(f"z={z}_context_checked_sum",loop[1],dj)
        full = response_jet(tree,loop)
        direct = reduce_jet(tree,loop,keep_final)
        upper = reduce_jet(tree,loop,keep_upper)
        check(f"z={z}_upper_induced_delta_C_nonzero",float(np.linalg.norm(upper[1][2])<1e-7))
        staged = reduce_jet(*upper,keep_final)
        reverse0 = reduce_jet(tree,loop,keep_reverse)
        reverse = reduce_jet(*reverse0,np.arange(len(keep_final)))
        for order in (0,1):
            for idx,label in enumerate(("K","J","C")):
                check(f"z={z}_order{order}_{label}_staged",direct[order][idx],staged[order][idx])
                check(f"z={z}_order{order}_{label}_reverse",direct[order][idx],reverse[order][idx])
            check(f"z={z}_order{order}_response",full[order],response_jet(*direct)[order])
        # Derivative of the independently implemented exact nonlinear map.
        eps = 2e-4
        plus = spectral.eliminate(*(a+eps*b for a,b in zip(tree,loop)),keep_final)
        minus = spectral.eliminate(*(a-eps*b for a,b in zip(tree,loop)),keep_final)
        for index,label in enumerate(("K","J","C")):
            check(f"z={z}_{label}_finite_difference",(plus[index]-minus[index])/(2*eps),direct[1][index],tol=2e-8)
        canonical = reparameterize_jet(tree,loop,-.5*dkin)
        check(f"z={z}_canonical_response",response_jet(*canonical)[1],full[1])
        check(f"z={z}_canonical_mass_insertion",canonical[1][0],-.5*(dkin@h+h@dkin))
        # Incorrectly normalizing J twice while retaining its K correction.
        double_counted = (loop[0],loop[1]-.5*dkin@j,loop[2])
        double_error = yuk.error(response_jet(tree,double_counted)[1],full[1])
        check(f"z={z}_double_scalar_leg_detected",float(double_error<1e-5))
        no_c = (direct[1][0],direct[1][1],np.zeros_like(c))
        omission = yuk.error(response_jet(direct[0],no_c)[1],full[1])
        check(f"z={z}_dropping_induced_loop_C_detected",float(omission<1e-5))
        rows.append({"z":{"real":float(np.real(z)),"imag":float(np.imag(z))},
                     "upper_delta_C_norm":float(np.linalg.norm(upper[1][2])),
                     "final_delta_C_norm":float(np.linalg.norm(direct[1][2])),
                     "drop_delta_C_response_error":omission,"double_count_scalar_leg_error":double_error,
                     "full_loop_response_norm":float(np.linalg.norm(full[1]))})
    for name,bad in (
        ("different_background",{**context,"background":{**context["background"],"sigma":.1265}}),
        ("different_scale",{**context,"mu":2*context["mu"]}),
        ("different_scheme",{**context,"subtraction":"DRbar"}),
        ("different_current_basis",{**context,"basis":"other"})):
        try:
            merge_contributions([pieces[0],{**pieces[1],"context":bad}])
        except ValueError:
            check("reject_merge_"+name,0.)
        else:
            check("reject_merge_"+name,1.)
    try:
        merge_contributions([pieces[0],pieces[0]])
    except ValueError:
        check("reject_double_counted_diagram",0.)
    else:
        check("reject_double_counted_diagram",1.)
    # Additional dense complex algebra case exercises delta D and delta B:
    # these need not be nonzero in the sparse actual upper vector subpiece.
    rng = np.random.default_rng(20260910)
    r = rng.normal(size=(12,12))
    k = r@r.T+np.eye(12)*(1+.4j)
    jt = rng.normal(size=(12,5))
    dr = rng.normal(size=(12,12))+1j*rng.normal(size=(12,12))
    dk = .02*(dr+dr.T)
    dj_test = .02*(rng.normal(size=(12,5))+1j*rng.normal(size=(12,5)))
    dc = np.diag(np.arange(5))*.003
    test_tree, test_loop = (k,jt,np.eye(5)*.1),(dk,dj_test,dc)
    reduced = reduce_jet(test_tree,test_loop,np.arange(4))
    eps = 1e-4
    plus = spectral.eliminate(*(a+eps*b for a,b in zip(test_tree,test_loop)),np.arange(4))
    minus = spectral.eliminate(*(a-eps*b for a,b in zip(test_tree,test_loop)),np.arange(4))
    for i,label in enumerate(("K","J","C")):
        check("synthetic_dense_deltaB_deltaD_"+label,(plus[i]-minus[i])/(2*eps),reduced[1][i],tol=1e-8)
    check("synthetic_dense_response",response_jet(test_tree,test_loop)[1],response_jet(*reduced)[1])
    try:
        reduce_jet((np.zeros((2,2)),np.ones((2,1)),np.zeros((1,1))),
                   (np.zeros((2,2)),np.zeros((2,1)),np.zeros((1,1))),[0])
    except np.linalg.LinAlgError:
        check("singular_elimination_rejected",0.)
    else:
        check("singular_elimination_rejected",1.)
    prescription = {"action":"P54PQ-v2","subtraction":"MSbar-DR",
                    "gauge":"background-field Landau","tadpoles":"fixed-VEV"}
    try:
        vector_matching.require_complete_matching(prescription)
    except ValueError as exc:
        blocker = str(exc)
        check("partial_KJC_cannot_enable_fit",0.)
    else:
        blocker = "ERROR"
        check("partial_KJC_cannot_enable_fit",1.)
    sources = [Path(__file__), Path(vector_matching.__file__),Path(spectral.__file__),
               RF/"output/p54_common_yukawa_phase.json", RF/"output/p54_ps_finite_thresholds.json"]
    report = {"schema":"p54-kjc-one-loop-subset-transport-v1","date":"2026-09-11",
              "scheme":prescription,"vacuum":geometry["upper"]["vacuum"],
              "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
              "scope":"exact one-loop jet algebra; actual computed vector subset, not all finite diagrams",
              "physical_chart_dimension":304,"upper_retained_dimension":249,
              "final_chart_dimension":129,"channel_order":names,
              "source_definition":"seven P-SPEC1 pairs plus uc(colour1)dc(colour2) and uL(colour1)dL(colour2); h/f factored out",
              "J0":j.tolist(),"delta_J_computed_subset":dj.tolist(),
              "rows":rows,"finite_matching_complete":False,"physical_fit_enabled":False,
              "fit_rejection":blocker,"checks":checks,
              "summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),
                         "all_pass":all(c["pass"] for c in checks)}}
    return report


def markdown(r):
    lines = ["# P54 one-loop (K,J,C) transport checkpoint", "",
             "**Complete algebra, partial diagram content. No physical fit.**", "",
             f"{r['summary']['passed']}/{r['summary']['total']} checks pass. The actual stationary upper PS chart has 304 physical real scalars; "
             "248 active parents plus PQ remain after eliminating 55 upper scalars. The final 129-dimensional chart retains whole H/F parents and PQ. "
             "Further elimination is an off-shell algebra regression, not a claim that L/R can physically be decoupled at the upper scale.", "",
             "The current slots retain seven P-SPEC1 pairs and add two action-defined different-colour diquark pairs, with the unknown two UV families factored out. "
             "The earlier same-colour pairs do not see the upper triplets: antisymmetric colour contractions vanish. The new probes detect their finite induced contact response. "
             "One-loop contributions are the computed active vector vertex, universal vector fermion external legs (including heavy-source endpoints), "
             "and active vector-scalar kinetic term. Missing self-energies, heavy-source 1PI vertices, boxes, operator mixing and lower matching remain unknown.", "",
             "The initial delta C is zero **for this selected additive set of graphs only**; the induced delta C is nonzero. "
             "This is not a zero boundary condition for the complete finite contact matching.", "",
             "| z=p_E^2 | Induced upper delta C norm | Drop final delta C response error | Double-count scalar leg error |",
             "|---|---:|---:|---:|"]
    for row in r["rows"]:
        z = complex(row["z"]["real"],row["z"]["imag"])
        lines.append(f"| {z} | {row['upper_delta_C_norm']:.7g} | {row['drop_delta_C_response_error']:.7g} | {row['double_count_scalar_leg_error']:.7g} |")
    lines += ["", "Errors use ||A-B||/max(1,||A||,||B||); they are not physical cross-section errors.", "",
              "The Frechet derivative agrees with direct finite differences and both elimination orders. "
              "Canonicalizing the field metric also transforms the mass matrix and sources; normalizing J while retaining the unnormalized K double counts external scalar legs.", "",
              "## Fit gate", "", r["fit_rejection"], "",
              "## Regression ledger", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir",type=Path)
    args = parser.parse_args()
    r = run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"rows":r["rows"],
                      "failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]:
        raise SystemExit(1)
