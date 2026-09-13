#!/usr/bin/env python3
"""Linear-CHN finite vertex and operator closure at sterile thresholds.

Direct CHN--Yukawa bubble, associated terminal C5 subset and Higgs tadpole;
exact field-dependent Schur generation of the dimension-six LNH operator.
Partial-threshold RG residual is identified with its mass-dependent
down-mixing. This is NOT a complete nuSMEFT anomalous dimension or matcher.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import verify_p54_sequential_seesaw as seq
import verify_p54_finite_seesaw_matching as finite

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_chn_threshold_closure"
LOOP=16*np.pi**2


def proper_vertex(y,masses,c,mu,internal=None):
    """L contains -Y LNH + C NN HdaggerH; Takagi mass coordinates.

    CHN vertex factor 2 and chirality-flip propagator M give
    deltaY=2 Y diag[M(1+log(mu^2/M^2))] C / (16 pi^2).
    External scalar/fermion wave contributions with one CHN vanish at
    this order qH=0; the nonzero Higgs MASS tadpole is separate.
    """
    m=np.asarray(masses)
    if np.any(m<=0) or mu<=0:raise ValueError("Positive Takagi masses and scale required")
    ids=np.arange(len(m)) if internal is None else np.asarray(internal,dtype=int)
    weight=m[ids]*(1+np.log(mu**2/m[ids]**2))
    return 2*(y[:,ids]*weight)@c[ids,:]/LOOP


def terminal_linear_chn(state,mu):
    """Only linear-CHN finite subset with all N removed and leading qH=0."""
    if state["qH"]!=0:raise ValueError("Declared leading light-mass truncation requires qH=0")
    m,u=seq.takagi(state["MR"]);y=state["Ynu"]@u;c=u.T@state["CHN"]@u
    dy=proper_vertex(y,m,c,mu)
    dc=-dy@(y/m).T-(y/m)@dy.T
    dq=4*np.sum(m**3*np.diag(c).real*(np.log(m**2/mu**2)-1))/LOOP
    return {"delta_C5":seq.sym(dc),"delta_qH":float(dq),"delta_Y_in_Takagi_chart":dy}


def partial_data(state,h,u,mu):
    """Removed-internal-line hard subset plus tree LNH Wilson coefficient.

    At the matching point M_hr=0 but C_hr need not vanish. The remaining
    internal-N graphs have an EFT realization through D6 and must not be
    silently added to a dimension-five-only matching coefficient.
    """
    m=np.diag(u.T@state["MR"]@u).real
    y=state["Ynu"]@u;c=u.T@state["CHN"]@u
    h=np.asarray(h,dtype=int);r=np.array([i for i in range(len(m)) if i not in h],dtype=int)
    dy=proper_vertex(y,m,c,mu,internal=h)
    yh,yr=y[:,h],y[:,r]
    dc=-dy[:,h]@(yh/m[h]).T-(yh/m[h])@dy[:,h].T
    d6=2*(yh/m[h])@c[np.ix_(h,r)]
    # Contribution of the generated O_LNH to beta_C5, in the declared
    # no-half CHN convention. Other d6 operators and mixings are not given.
    downmix=2*((d6*m[r])@yr.T+yr@(d6*m[r]).T)/LOOP
    delta_effective=dc-dy[:,r]@(yr/m[r]).T-(yr/m[r])@dy[:,r].T
    return {"m":m,"y":y,"c":c,"h":h,"r":r,"D6":d6,"delta_C5_hard_HH":seq.sym(dc),
            "delta_Yr_hard":dy[:,r],"delta_effective_C5_subset":seq.sym(delta_effective),
            "beta_C5_D6_mass_downmix":seq.sym(downmix)}


def chn_beta_effective(state):
    off=seq.state_copy(state);off["CHN"]*=0
    # Isolate CHN feedback, while retaining the actual other couplings.
    return finite.effective_beta(state)-finite.effective_beta(off)


def run():
    local=json.loads(seq.INPUT.read_text());chn=json.loads(seq.CHN_INPUT.read_text())
    base=seq.load_local_light(local["cases"][0],seq.decode(chn["rays"][2]["CHN_matrix_times_omega"]))
    base["qH"]=0.
    rng=np.random.default_rng(20260912)
    checks=[]
    def check(name,a,b=0.,tol=2e-9):
        a,b=np.asarray(a),np.asarray(b)
        residual=float(np.linalg.norm(a-b)/max(1,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append({"name":name,"residual":residual,"tolerance":tol,"pass":residual<tol})
    mu=.043;m,u=seq.takagi(base["MR"]);y=base["Ynu"]@u
    c=u.T@base["CHN"]@u
    for i,mi in enumerate(m):
        integral=-quad(lambda t:np.log(t*mi**2/mu**2),0,1,epsabs=1e-12)[0]
        check(f"massive_massless_vertex_integral_{i}",integral,1+np.log(mu**2/mi**2))
    result=terminal_linear_chn(base,mu)
    step=1e-4
    scale_derivative=(terminal_linear_chn(base,mu*np.exp(step))["delta_C5"]-
                      terminal_linear_chn(base,mu*np.exp(-step))["delta_C5"])/(2*step)
    out,_,_=seq.decouple_block(base,list(range(3)))
    check("all_removed_CHN_C5_scale_derivative_beta_difference",scale_derivative,
          chn_beta_effective(out)-chn_beta_effective(base))
    rotations=[seq.unitary(rng,3) for _ in range(6)]
    rotated=seq.family_transform(base,*rotations)
    transformed=terminal_linear_chn(rotated,mu)
    check("terminal_CHN_C5_family_covariance",transformed["delta_C5"],
          rotations[3].T@result["delta_C5"]@rotations[3])
    check("terminal_CHN_Higgs_tadpole_family_invariance",transformed["delta_qH"],result["delta_qH"])
    # Independent finite fermion determinant derivative: M(rho)=M-2Crho.
    def potential(rho):
        mm=base["MR"]-2*rho*base["CHN"]
        x=np.linalg.eigvalsh(mm.conj().T@mm)
        return -np.sum(x*x*(np.log(x/mu**2)-1.5))/(2*LOOP)
    dr=1e-6
    check("CHN_finite_Higgs_mass_direct_determinant",(potential(dr)-potential(-dr))/(2*dr),result["delta_qH"],tol=1e-8)
    dscale=(terminal_linear_chn(base,mu*np.exp(step))["delta_qH"]-
            terminal_linear_chn(base,mu*np.exp(-step))["delta_qH"])/(2*step)
    off=seq.state_copy(base);off["CHN"]*=0
    check("CHN_finite_Higgs_mass_scale_derivative",dscale,-(seq.beta(base)["qH"]-seq.beta(off)["qH"]))
    deg=seq.state_copy(base);deg["MR"]=.035*np.eye(3,dtype=complex)
    orth=np.linalg.qr(rng.normal(size=(3,3)))[0]
    dy=proper_vertex(deg["Ynu"],np.full(3,.035),deg["CHN"],mu)
    rdy=proper_vertex(deg["Ynu"]@orth,np.full(3,.035),orth.T@deg["CHN"]@orth,mu)
    check("degenerate_CHN_vertex_O3_covariance",rdy,dy@orth)
    # Nonaligned stress case tests mixed blocks that align only at the UV
    # boundary of the actual scalar-tree dictionary.
    stress=seq.state_copy(base)
    noise=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    stress["CHN"]+=.03*seq.sym(noise)
    tests=[]
    for h in ([2],[1,2],list(range(3))):
        data=partial_data(stress,h,u,mu)
        after,_,_=seq.decouple_block(stress,h,u)
        deriv=(partial_data(stress,h,u,mu*np.exp(step))["delta_effective_C5_subset"]-
               partial_data(stress,h,u,mu*np.exp(-step))["delta_effective_C5_subset"])/(2*step)
        naive=chn_beta_effective(after)-chn_beta_effective(stress)
        corrected=naive+data["beta_C5_D6_mass_downmix"]
        check(f"partial_CHN_matching_RG_closes_with_D6_{len(h)}",deriv,corrected)
        rr=data["r"];hh=data["h"];yy=data["y"];cc=data["c"];mm=np.diag(data["m"])
        def schur_y(rho):
            shifted=mm-2*rho*cc
            return yy[:,rr]-yy[:,hh]@np.linalg.solve(shifted[np.ix_(hh,hh)],shifted[np.ix_(hh,rr)])
        if len(rr):
            check(f"tree_LNH_equals_exact_field_Schur_derivative_{len(h)}",(schur_y(dr)-schur_y(-dr))/(2*dr),data["D6"],tol=1e-8)
            check(f"dimension5_only_mixed_block_residual_detected_{len(h)}",float(np.linalg.norm(deriv-naive)<1e-9))
        tests.append({"removed":len(h),"naive_RG_residual_norm":float(np.linalg.norm(deriv-naive)),
            "D6":seq.cjson(data["D6"]),"mass_downmix":seq.cjson(data["beta_C5_D6_mass_downmix"])})
    # Inspect the ACTUAL old trajectories; do not promote or silently
    # modify them by applying an incomplete new matcher.
    trajectories=[]
    for i,case in enumerate(local["cases"]):
        state=seq.load_local_light(case,seq.decode(chn["rays"][2+i]["CHN_matrix_times_omega"]))
        rows=[]
        def inspect_and_tree_match(s,h,uu,scale):
            data=partial_data(s,h,uu,scale)
            rows.append({"mu":scale,"removed":len(h),"retained":len(data["r"]),
                "D6":seq.cjson(data["D6"]),"D6_norm":float(np.linalg.norm(data["D6"])),
                "CHN_hr_norm":float(np.linalg.norm(data["c"][np.ix_(data["h"],data["r"])])),
                "mass_downmix_norm":float(np.linalg.norm(data["beta_C5_D6_mass_downmix"])),
                "finite_Yr_hard_norm":float(np.linalg.norm(data["delta_Yr_hard"])),
                "new_terms_applied_to_trajectory":False})
            return seq.decouple_block(s,h,uu)
        seq.evolve(state,.1,1e-4,threshold_matcher=inspect_and_tree_match)
        check(f"actual_trajectory_three_events_{i}",len(rows),3)
        check(f"actual_running_mixed_CHN_generates_D6_{i}",float(max(x["D6_norm"] for x in rows)==0))
        trajectories.append({"synthetic_case_with_actual_group_dictionary":i,"events":rows})
    paths=[Path(__file__),Path(seq.__file__),Path(finite.__file__),seq.INPUT,seq.CHN_INPUT]
    return {"schema":"p54-linear-CHN-finite-and-operator-closure-v1","date":"2026-09-12",
        "scope":"linear CHN, leading qH=0 finite vertex/C5/Higgs-mass subset; induced D6 required in sequential mass-dependent RG",
        "terminal_linear_CHN":{"delta_C5":seq.cjson(result["delta_C5"]),"delta_qH":result["delta_qH"],
                               "complete_finite_seesaw":False},
        "partial_block_closure_tests":tests,"actual_trajectory_operator_audit":trajectories,
        "old_trajectories_modified":False,"full_D6_anomalous_dimensions":False,
        "complete_finite_CHN_matching":False,"physical_fit_enabled":False,
        "missing":["full dimension-six operator matching/running and its mass-suppressed feedback",
            "nonzero-qH terms and finite correlated MR,lambda and field maps",
            "remaining pre-existing-C5/mixed-sterile/box diagrams and lower scalar/gauge matching"],
        "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# CHN finite vertex and sequential operator closure", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "**This is not a complete finite seesaw matcher. Existing diagnostic trajectories are not modified.**", "",
        "The no-half CHN interaction gives M(rho)=M-2 C rho. A direct Yukawa/CHN bubble has finite vertex 2 Y diag[M(1+log(mu^2/M^2))] C/(16 pi^2). Its all-removed linear-CHN C5 and finite Higgs-mass subsets pass independent RG, integral, determinant and covariance checks at leading qH=0.", "",
        "At a partially removed Takagi block the exact tree Schur map generates Y_eff(rho)=Y_r+D6 rho, D6=2 Y_h M_h^-1 C_hr. The operator -(D6 LHN_r)(HdaggerH)+h.c. has dimension six, but retained Majorana masses allow it to mix down into C5. The explicit mixed-block beta residual equals this missing down-mixing. No hierarchy suppression is silently assumed.", "",
        "| Actual diagnostic case | Threshold | D6 norm | Mixed CHN norm |", "|---|---:|---:|---:|"]
    for case in r["actual_trajectory_operator_audit"]:
        for i,row in enumerate(case["events"]):
            lines.append(f"| {case['synthetic_case_with_actual_group_dictionary']} | {i} | {row['D6_norm']:.8g} | {row['CHN_hr_norm']:.8g} |")
    lines += ["", "The scalar-tree CHN can align with MR at the input, but actual running generates mixed Takagi entries. These numbers use the pre-existing synthetic families, not fitted data.", "",
        "An additional controlled expansion in retained/removed Majorana masses can consistently place the D6 feedback beyond leading dimension-five accuracy. Keeping retained-mass effects without that expansion instead needs the generated operator or an explicit remainder bound; no fitted-trajectory bound is claimed here.", "",
        "## Missing", ""]+["- "+s for s in r["missing"]]
    lines += ["", "## Checks", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__=="__main__":
    r=run();OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"failed":[c for c in r["checks"] if not c["pass"]],
        "actual_D6_norms":[[e["D6_norm"] for e in c["events"]] for c in r["actual_trajectory_operator_audit"]]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
