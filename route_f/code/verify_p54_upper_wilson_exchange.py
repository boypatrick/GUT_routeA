#!/usr/bin/env python3
"""Factorizable one-loop upper scalar-exchange four-fermion Wilson response.

Actual 24 Yukawa-coupled real upper six modes, complete PS-invariant mass
subblock, hard bosonic self-energy and common-CT valley response; actual
pure-Yukawa/vector vertex and fermion-leg contractions. Nonfactorizable
boxes and the fermionic fixed-VEV counterterms are NOT supplied.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

import verify_p54_upper_yukawa_thresholds as yuk
import verify_p54_upper_scalar_matching as scalar
import verify_p54_upper_vector_matching as vec
import verify_p54_kjc_finite_matching as kjc

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_upper_wilson_exchange"
LOOP=16*np.pi**2


def six_space(geometry):
    norms=np.linalg.norm(geometry["heavy_h"],axis=(1,2))+np.linalg.norm(geometry["heavy_f"],axis=(1,2))
    indices=np.flatnonzero(norms>1e-9)
    if len(indices)!=24:raise ValueError("Expected four real six-dimensional PS copies")
    mass=geometry["masses2"][indices]
    groups=[]
    for i,m in enumerate(mass):
        if not groups or abs(m-mass[groups[-1][0]])>1e-9:groups.append([i])
        else:groups[-1].append(i)
    if [len(g) for g in groups]!=[6,6,6,6]:raise ValueError("Unexpected upper six-copy degeneracy")
    bh=geometry["heavy"][:,indices]
    p1=geometry["ps"]["common_geometry"]["p1"]
    helper=yuk.module("wilson_PS_generators",RF/"code/verify_p54_p2_two_site_matching.py")
    gens=[bh.T@p1.representation_matrix(g)@bh for group in helper.ps_generators(p1).values() for g in group]
    mats,names=[],[]
    max_intertwiner_error=0.
    for i,group in enumerate(groups):
        p=np.zeros((24,24));p[np.ix_(group,group)]=np.eye(6)/math.sqrt(6)
        mats.append(p);names.append(f"copy{i}_norm")
        for j in range(i):
            other=groups[j]
            # Row-major vectorization: Ri X - X Rj = 0.
            equations=[]
            for r in gens:
                ri,rj=r[np.ix_(group,group)],r[np.ix_(other,other)]
                equations.append(np.kron(ri,np.eye(6))-np.kron(np.eye(6),rj.T))
            ns=null_space(np.vstack(equations),rcond=1e-9)
            if ns.shape[1]!=1:raise ValueError("Six copies must have one real intertwiner")
            inter=ns[:,0].reshape(6,6)
            inter*=np.sign(inter.ravel()[np.argmax(abs(inter))])
            q=np.zeros((24,24));q[np.ix_(group,other)]=inter/math.sqrt(2)
            q[np.ix_(other,group)]=inter.T/math.sqrt(2)
            max_intertwiner_error=max(max_intertwiner_error,max(yuk.error(r@q,q@r) for r in gens))
            mats.append(q);names.append(f"copy{i}_copy{j}_mixed")
    return indices,bh,np.array(mats),names,max_intertwiner_error


def run(cache_dir):
    g=yuk.build_upper_geometry(cache_dir)
    p1=g["ps"]["common_geometry"]["p1"]
    indices,bh,mats,names,ward=six_space(g)
    full=json.loads((RF/"output/p54_full_doublet_cw.json").read_text())
    upper=json.loads((RF/"output/p54_upper_scalar_matching.json").read_text())
    mu=full["scheme"]["mu_over_omega"]
    if upper["common_CT_upper_tadpole"]["mu"]!=mu:raise ValueError("Tadpole anchor mismatch")
    store=scalar.ActionJets(p1,g["upper"]["same_action_parameters"],cache_dir)
    vectors=vec.build_vector_geometry(g)
    basis,m,slope,bubble,a0=scalar.loop_weights(g,mu)
    directions,measurement=scalar.directions_from_basis(mats)
    checks=[]
    def check(name,a,b=0.,tol=3e-8):
        err=yuk.error(np.asarray(a),np.asarray(b))
        checks.append({"name":name,"residual":err,"tolerance":tol,"pass":err<tol})
    check("24_Yukawa_coupled_upper_real_modes",len(indices),24)
    check("ten_real_symmetric_six_copy_invariants",np.einsum("aij,bij->ab",mats,mats),np.eye(10))
    check("all_six_copy_intertwiners_PS_covariant",ward)
    x,h0=g["x"],g["hessian"]
    def mass_contraction(t,u):
        te=basis.T@t@basis;ue=basis.T@u@basis
        return (np.dot(a0,np.diag(ue))-np.sum(te*te*bubble))/(2*LOOP)
    values=[]
    for i,q in enumerate(directions):
        t,u=store.jets(x,h0,bh@q,"six_mass_probe_"+str(i))
        values.append(mass_contraction(t,u))
    coeff=np.linalg.solve(measurement,values)
    ms=np.einsum("a,aij->ij",coeff,mats)
    rng=np.random.default_rng(20260912)
    q=rng.normal(size=24);q/=np.linalg.norm(q)
    t,u=store.jets(x,h0,bh@q,"six_independent_probe")
    check("independent_full_six_mass_reconstruction",mass_contraction(t,u),q@ms@q,tol=5e-7)
    ra=np.einsum("Aij,ja->Aia",vectors["R"],bh)
    rx=np.einsum("Aij,j->Ai",vectors["R"],x)
    check("six_linear_VV_mass_vertices_vanish",np.einsum("Ai,Bia->ABa",rx,ra))
    cs=np.einsum("Aia,Aib->ab",ra,ra)
    mv=3*vectors["g"]**2*vectors["mass2"]*(math.log(vectors["mass2"]/mu**2)-1/3)*cs/LOOP
    cts=full["tadpoles"]["finite_mass_counterterms_over_omega2"]
    ct=np.zeros((328,328));ct[:54,:54]=-cts["delta_mu2"]*np.eye(54)
    ct[54:306,54:306]=-cts["delta_nu2"]*.5*np.eye(252)
    ct[326:,326:]=-cts["delta_mus2"]*np.eye(2)
    shift=upper["common_CT_upper_tadpole"]["hard_plus_CT_linear_displacement"]
    radial=[p1.vacuum_vector(1,0,0),p1.vacuum_vector(0,0,1)]
    ms_shift=np.zeros((24,24))
    for i,r in enumerate(radial):
        t,u=store.jets(x,h0,r,"six_radial_reuse_"+str(i))
        ms_shift+=shift[i]*(bh.T@t@bh)
    dm=ms+mv+bh.T@ct@bh+ms_shift
    d=np.diag(g["masses2"][indices])
    check("same_action_tree_six_mass_block",bh.T@h0@bh,d)
    cases=[]
    kp_v,_=vec.kinetic(g,vectors,mu)
    weight=3*vec.b0(vectors["mass2"],0.,mu)-2
    for i,(h,f) in enumerate(((.20+.04j,.025-.01j),(.13-.07j,.041+.013j))):
        # Deliberately labelled one-family benchmarks. No fit is performed.
        y=g["heavy_h"][indices]*h+g["heavy_f"][indices]*f
        all_heavy=g["heavy_h"]*h+g["heavy_f"]*f
        py=yuk.threshold_from_tensors(y,all_heavy,g["masses2"],mu)
        gv=-vectors["g"]**2/LOOP*weight*sum(t.T@y@t for t in vectors["T"])
        dy=py["delta"]+gv-.5*(kp_v.T@y+y@kp_v)
        j,order=kjc.source_columns({"combined":y},g)
        dj,_=kjc.source_columns({"combined":dy},g)
        invj=np.linalg.solve(d,j)
        c0=j.T@invj
        c1=dj.T@invj+invj.T@dj-invj.T@dm@invj
        check(f"case{i}_tree_C_positive_semidefinite",max(0.,-np.linalg.eigvalsh(c0).min()))
        check(f"case{i}_finite_C_transpose_symmetry",c1,c1.T)
        eps=1e-4
        def exact(e):return (j+e*dj).T@np.linalg.solve(d+e*dm,j+e*dj)
        check(f"case{i}_exact_resolvent_finite_difference",(exact(eps)-exact(-eps))/(2*eps),c1,tol=1e-6)
        # Arbitrary scalar reparameterization cancels between vertices and mass.
        raw=rng.normal(size=(24,24));s=.01*(raw+raw.T)
        transformed_d=dm+s.T@d+d@s
        transformed_j=dj+s.T@j
        cc=transformed_j.T@invj+invj.T@transformed_j-invj.T@transformed_d@invj
        check(f"case{i}_heavy_kinetic_redefinition_cancels_in_C",cc,c1)
        check(f"case{i}_self_energy_insertion_is_not_zero",float(np.linalg.norm(invj.T@dm@invj)<1e-7))
        cases.append({"synthetic_h_raw":yuk.cjson(h),"synthetic_f_raw":yuk.cjson(f),
            "source_order":order,"tree_C":c0.tolist(),"delta_C_factorizable":c1.tolist(),
            "tree_C_norm":float(np.linalg.norm(c0)),"delta_C_norm":float(np.linalg.norm(c1)),
            "delta_C_over_tree_norm":float(np.linalg.norm(c1)/np.linalg.norm(c0)),
            "full_finite_Wilson":False})
    invroot=np.diag(1/np.sqrt(np.diag(d)))
    mass_control,mass_modes=np.linalg.eigh(invroot@dm@invroot)
    radius=float(max(abs(mass_control)))
    dominant=mass_modes[:,np.argmax(abs(mass_control))]
    mode_contributions={name:float(dominant@invroot@value@invroot@dominant)
        for name,value in {"scalar":ms,"vector":mv,"common_bosonic_CT":bh.T@ct@bh,
                           "hard_plus_CT_valley_response":ms_shift}.items()}
    paths=[Path(__file__),Path(scalar.__file__),Path(yuk.__file__),Path(vec.__file__),Path(kjc.__file__),
           Path(p1.__file__),RF/"output/p54_upper_scalar_matching.json",RF/"output/p54_full_doublet_cw.json"]
    return {"schema":"p54-upper-factorizable-Wilson-v1","date":"2026-09-11",
        "scope":"scalar-exchange four-fermion Wilson factorizable one-loop subset; no boxes or fitted families",
        "vacuum":g["upper"]["vacuum"],"mu":mu,"heavy_real_dimension":24,
        "invariant_names":names,"scalar_mass_coefficients":coeff.tolist(),
        "tree_mass":d.tolist(),"delta_mass_scalar":ms.tolist(),"delta_mass_vector":mv.tolist(),
        "delta_mass_CT_and_valley":(bh.T@ct@bh+ms_shift).tolist(),"delta_mass_total_subset":dm.tolist(),
        "dominant_relative_mass_mode_contributions":mode_contributions,
        "relative_mass_insertion_eigenvalues":mass_control.tolist(),"cases":cases,
        "propagator_expansion_control":{
            "scope":"computed bosonic/common-CT mass subset only; not full-model exclusion",
            "spectral_radius":radius,
            "Neumann_series_converges_at_unit_loop_parameter":radius<1.,
            "formal_loop_parameter_convergence_radius":1/radius if radius else None,
            "resummed_subset_positive_relative_mass":bool(min(1+mass_control)>0),
            "small_source_projection_does_not_certify_full_mass_control":True},
        "complete_upper_Wilson_matching":False,"physical_fit_enabled":False,
        "missing":["nonfactorizable scalar/vector/fermion boxes and direct vector exchange",
            "complete Lorentz/gauge/family operator basis and mixing","fermionic fixed-VEV counterterms and final light retuning",
            "upper scalar three/four-point potential Wilson vertices","controlled retained-mass corrections"],
        "cache":{"hits":store.hits,"new_vertex_evaluations":store.evaluated,"keys":store.keys},
        "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# Upper factorizable scalar-exchange Wilson matching", "",
        "**One-loop scalar-exchange subset, not a complete four-fermion matching.**", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "All 24 Yukawa-coupled real upper six modes are retained, including the ten real PS-invariant mass mixings between their four copies. Actual scalar/vector mass diagrams, the broken-anchored bosonic counterterms and the hard heavy-valley response enter the same propagator. Pure-Yukawa/vector proper vertices and external fermion legs enter J. Heavy scalar external-leg normalization is not inserted twice.", "",
        "C0=J^T D^-1 J; delta C=delta J^T D^-1 J+J^T D^-1 delta J-J^T D^-1 delta D D^-1 J.", "",
        "| Synthetic one-family case | ||C0|| | ||delta C|| | Ratio |", "|---|---:|---:|---:|"]
    lines += [f"| {i} | {c['tree_C_norm']:.8g} | {c['delta_C_norm']:.8g} | {c['delta_C_over_tree_norm']:.8g} |" for i,c in enumerate(r["cases"])]
    control=r["propagator_expansion_control"]
    lines += ["", "These are diagnostic source kernels, not proton lifetimes or fitted physical amplitudes. The ratio alone is not a full perturbativity assessment because important diagrams and parameter feedback remain absent.", "",
        "## Non-small mass insertion: do not promote to a fit", "",
        f"The spectral radius of D^(-1/2) deltaD D^(-1/2) is {control['spectral_radius']:.8g}. Its Neumann expansion does not converge at unit loop parameter for this computed subset, even though the resummed subset mass is positive. Smaller corrections in selected source projections cannot certify the whole propagator. Missing diagrams, common parameter feedback and the declared truncation must be resolved before judging the complete model. The algebra checks below do not assert perturbative control.", "",
        "Contributions on the dominant normalized mass direction (not separate eigenvalues):", "",
        "```json",json.dumps(r["dominant_relative_mass_mode_contributions"],indent=2),"```", "",
        "## Missing contributions", ""]+["- "+s for s in r["missing"]]
    lines += ["", "## Checks", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--cache-dir",type=Path,required=True)
    args=ap.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"ratios":[x["delta_C_over_tree_norm"] for x in r["cases"]],
        "mass_insertion_range":[min(r["relative_mass_insertion_eigenvalues"]),max(r["relative_mass_insertion_eigenvalues"])],
        "failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
