#!/usr/bin/env python3
"""Bounded action-preserving scalar reselection after the P-IR1 failure.

Uniformly weaken existing scalar couplings and relax the VEV hierarchy.
Dependent masses are re-solved, and the tree light-doublet condition is
reimposed. No measured masses or flavor data enter. Candidates are NOT
certified physical backgrounds and do not replace matching inputs.
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

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_scalar_reselection"
LOOP=16*np.pi**2


def run(cache_dir):
    oldpath=RF/"output/p54_full_doublet_cw.json"
    irpath=RF/"output/p54_goldstone_ir.json"
    old=json.loads(oldpath.read_text());ir=json.loads(irpath.read_text())
    if ir["hard_scalar_kinetic"]["perturbation_smaller_than_tree_metric"]:
        raise ValueError("Reselection is conditional on the diagnosed control failure")
    p1=common.yuk.module("reselect_p1",RF/"code/verify_p54_p1_hessian_spectrum.py")
    cw=common.yuk.module("reselect_cw",RF/"code/verify_p54_p2_bosonic_cw.py")
    search=common.yuk.module("reselect_search",RF/"code/search_p54_hierarchical_p1.py")
    p=old["tree_parameters"];r0=np.array([old["vacuum"][k] for k in ("omega","sigma","vs")])
    mu=old["scheme"]["mu_over_omega"];mu2=mu*mu
    store=common.scalar.ActionJets(p1,p,cache_dir)
    rad=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)]);metric=rad.T@rad
    x0=rad@r0;h0=store.hessian(x0,"old_base")
    # Only sigma changes, so its cached quadratic jet determines the
    # complete ambient Hessian at EVERY tested new radial point exactly.
    tr=[];ur=[]
    for i in range(3):
        ti,ui=store.jets(x0,h0,rad[:,i],f"reselect_radial_{i}")
        tr.append(ti);ur.append(ui)
    second={}
    for i in range(3):
        for j in range(i,3):
            if i==j:second[i,j]=ur[i]
            else:
                step=.02
                hp=store.hessian(x0+step*(rad[:,i]+rad[:,j]),f"reselect_pair_{i}{j}")
                second[i,j]=(hp-h0-step*(tr[i]+tr[j])-.5*step**2*(ur[i]+ur[j]))/step**2
            second[j,i]=second[i,j]
    p10=np.zeros_like(h0);p10[306:326,306:326]=np.eye(20)
    gens=cw.generators();orad=[cw.field_orbit(p1,q,gens) for q in rad.T]
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        a,b=np.asarray(a),np.asarray(b)
        err=float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append({"name":name,"residual":err,"tolerance":tol,"pass":err<tol})
    rows=[]
    # Four declared coarse engineering probes, not a minimization against
    # observations. No fine scan or extra free Lagrangian coefficient.
    for z,sigma in ((.03,.3),(.03,.5),(.1,.3),(.1,.5)):
        r=np.array([1.,sigma,.25]);x=rad@r;ds=sigma-r0[1]
        pars={key:z*val for key,val in p.items()}
        # Same normalization as the actual action; CT map also describes
        # an exact change in its linearly occurring dependent mass inputs.
        dmass=-np.linalg.solve(common.radial_mass_map(r),p1.radial_gradient(r,pars))
        for key,value in zip(common.MASS_KEYS,dmass):pars[key]+=value
        h=z*(h0+ds*tr[1]+.5*ds*ds*ur[1])+common.mass_operator(dmass)+z*p["xi02"]*p10
        pars["xi02"]=0.
        symmetry,sg,comp=search.symmetry_complement(p1,x)
        tuning=search.tune_first_instability(h,p10,comp)
        tag=f"z{z}_sigma{sigma}"
        check(tag+"_stationary",p1.radial_gradient(r,pars))
        check(tag+"_exact_polynomial_tree_Goldstones",h@sg)
        if not tuning.get("tunable"):
            rows.append({"z":z,"sigma":sigma,"tree_tuning":tuning,"accepted_for_next_radial_audit":False})
            continue
        pars["xi02"]=tuning["xi02"];h-=pars["xi02"]*p10
        lam,u=np.linalg.eigh(h)
        check(tag+"_38_zero_modes",lam[:38])
        check(tag+"_positive_tree_complement",float(lam[38]<=0))
        classification=p1.classify_extra_zero_modes(lam,u,tuning["zero_tolerance"],symmetry)
        check(tag+"_one_complex_doublet_rank",classification["rank"],4)
        for group,value in (("SU3",0.),("SU2",.75),("Y",.25)):
            check(tag+"_doublet_Casimir_"+group,classification["casimir_eigenvalues"][group],np.full(4,value))
        t=[z*(tr[i]+ds*second[i,1]) for i in range(3)]
        orbit=cw.field_orbit(p1,x,gens);mv=mu2*orbit.T@orbit
        vl,vu=np.linalg.eigh(mv)
        vt=[mu2*(oi.T@orbit+orbit.T@oi) for oi in orad]
        ts=np.array([np.dot(lam[38:]*(np.log(lam[38:]/mu2)-1),np.diag(u[:,38:].T@v@u[:,38:]))/(2*LOOP) for v in t])
        tv=np.array([3*np.dot(vl[-33:]*(np.log(vl[-33:]/mu2)-1/3),np.diag(vu[:,-33:].T@v@vu[:,-33:]))/(2*LOOP) for v in vt])
        dp=-np.linalg.solve(common.radial_mass_map(r),ts+tv)
        ct=common.mass_operator(dp)
        hs=np.zeros((3,3));hv=hs.copy()
        for i in range(3):
            for j in range(3):
                hs[i,j]=cw.hard_trace_hessian(lam,u,t[i],t[j],z*second[i,j],290,mu2,1.5,1.)
                vv=mu2*(orad[i].T@orad[j]+orad[j].T@orad[i])
                hv[i,j]=cw.hard_trace_hessian(vl,vu,vt[i],vt[j],vv,33,mu2,5/6,3.)
        tree=rad.T@h@rad;hard=tree+hs+hv+rad.T@ct@rad
        check(tag+"_common_fixed_VEV_condition",ts+tv+rad.T@ct@x)
        # Group-degenerate state sums for the hard scalar kinetic matrix.
        groups=[np.arange(38)];m=[0.];start=38
        while start<328:
            stop=start+1
            while stop<328 and abs(lam[stop]-lam[start])<1e-8:stop+=1
            groups.append(np.arange(start,stop));m.append(float(lam[start:stop].mean()));start=stop
        te=[u.T@v@u for v in t];dz=np.zeros((3,3))
        for k,a in enumerate(m):
            for l,b in enumerate(m):
                if k+l==0:continue
                slope=quad(lambda xx:xx*(1-xx)/((1-xx)*a+xx*b),0,1,epsabs=1e-12)[0]
                for i in range(3):
                    for j in range(3):
                        dz[i,j]+=slope*np.sum(te[i][np.ix_(groups[k],groups[l])]*te[j][np.ix_(groups[k],groups[l])])/(2*LOOP)
        masses=eigh(hard,metric,eigvals_only=True)
        rel=eigh(hard-tree,tree,eigvals_only=True)
        wave=eigh(dz,metric,eigvals_only=True)
        # Two independent CONTROL tests, neither implies a whole-vacuum or
        # pole certificate. A strict 0.3 engineering margin is explicit.
        accepted=bool(masses.min()>0 and np.max(abs(rel))<.3 and wave.max()<.3)
        rows.append({"z":z,"sigma":sigma,"vs":.25,"tree_parameters":pars,"tree_tuning":tuning,
            "extra_zero_mode_classification":classification,"tree_min_hard_mass2":float(lam[38]),
            "tree_radial_eigenvalues":eigh(tree,metric,eigvals_only=True).tolist(),
            "radial_hard_Hessian":hard.tolist(),"radial_hard_eigenvalues":masses.tolist(),
            "relative_radial_insertions":rel.tolist(),"hard_scalar_kinetic_eigenvalues":wave.tolist(),
            "common_finite_mass_CT":dp.tolist(),"accepted_for_next_radial_audit":accepted,
            "is_new_physical_benchmark":False})
    accepted=[i for i,row in enumerate(rows) if row["accepted_for_next_radial_audit"]]
    sources=[Path(__file__),Path(p1.__file__),Path(cw.__file__),Path(search.__file__),Path(common.__file__),Path(common.scalar.__file__),oldpath,irpath]
    return {"schema":"p54-bounded-scalar-reselection-v1","date":"2026-09-13",
        "scope":"four declared weak-coupling/less-hierarchical engineering probes; no measured input fit or default mutation",
        "coupling_scale_rule":"scale every existing scalar coefficient by z; re-solve mu2,nu2,mus2 at new radii; retune xi02 at first physical doublet crossing",
        "gauge_and_mu_held_at_reference":mu,"engineering_margin":.3,"candidates":rows,
        "accepted_candidate_indices":accepted,"selected_physical_benchmark":None,
        "missing":["full nonradial one-loop stability and complete gauge/fermion momentum consistency",
            "one-loop light-doublet retuning and actual new upper PS saddle/matching",
            "two-loop gauge running and full thresholds at any reselected point"],
        "old_matching_inputs_changed":False,"physical_fit_enabled":False,
        "cache":{"hits":store.hits,"new_Hessian_evaluations":store.evaluated},
        "source_sha256":{str(path.relative_to(RF)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# Bounded P54 scalar parameter reselection","",
        "**Engineering probes only; none is automatically a viable model or a replacement matching input.**","",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.","",
        "| Scalar multiplier | sigma/omega | Lowest hard radial eigenvalue | Max relative radial insertion | Max hard scalar kinetic correction | Pass 0.3 radial-control margin |",
        "|---:|---:|---:|---:|---:|:---:|"]
    for c in r["candidates"]:
        if "radial_hard_eigenvalues" in c:
            lines.append(f"| {c['z']} | {c['sigma']} | {min(c['radial_hard_eigenvalues']):.8g} | {max(abs(v) for v in c['relative_radial_insertions']):.8g} | {max(c['hard_scalar_kinetic_eigenvalues']):.8g} | {c['accepted_for_next_radial_audit']} |")
        else:lines.append(f"| {c['z']} | {c['sigma']} | tree tuning failed | - | - | False |")
    lines += ["","The exact quartic-action polynomial jets avoid new Hessian evaluations. All accepted points would still require full nonradial, gauge/momentum, light-doublet and threshold verification. No experimental mass was used. Parameter cards, tree classifications and all finite tadpole shifts are retained in JSON.","","## Missing",""]
    lines += ["- "+s for s in r["missing"]]
    return "\n".join(lines)+"\n"


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--cache-dir",type=Path,required=True)
    args=ap.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(markdown(r));print("failed",[c for c in r["checks"] if not c["pass"]])
    if not r["summary"]["all_pass"]:raise SystemExit(1)
