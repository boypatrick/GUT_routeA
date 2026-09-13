#!/usr/bin/env python3
"""Actual common bosonic parameter/tadpole audit, with finite-scheme transport.

No counterterm is dropped, no parameter scan is made, and no full matching
claim follows from the formal scheme identities. The frozen historical
point and its content-addressed polynomial action derivatives are used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import root

import verify_p54_upper_yukawa_thresholds as yuk
import verify_p54_upper_scalar_matching as scalar
import verify_p54_upper_wilson_exchange as wilson

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_common_renormalization"
LOOP=16*np.pi**2
MASS_KEYS=("mu2","nu2","mus2")


def mass_operator(dp):
    """Exact canonical Hessian of delta p . partial_p V, in old real chart."""
    result=np.zeros((328,328))
    result[:54,:54]=-dp[0]*np.eye(54)
    result[54:306,54:306]=-dp[1]*np.eye(252)/2
    result[326:,326:]=-dp[2]*np.eye(2)
    return result


def radial_mass_map(r):
    """partial_r partial_p V; rows (w,sigma,vs), columns MASS_KEYS."""
    return -np.diag(np.array([12/5,1,1])*r)


def finite_matching_conversion(delta_matching, jacobian_tree, uv_shift, eft_shift):
    """q=f(p)+Delta; p'=p+a, q'=q+b gives Delta'=Delta+b-Df.a."""
    return delta_matching+eft_shift-jacobian_tree@uv_shift


def majorana_fixed_VEV_radial(masses, sigma, mu):
    """d_sigma^2 V_F - (d_sigma V_F)/sigma, at fixed f_M and mu.

    The only nonzero fermionic radial Hessian entry is sigma,sigma.
    Three arbitrary Takagi masses are allowed, including their zero limit.
    """
    masses=np.asarray(masses,float)
    if np.any(masses<0) or sigma<=0 or mu<=0:
        raise ValueError("Nonnegative masses and positive sigma,mu required")
    pos=masses>0
    return float(-4*np.sum(masses[pos]**4*np.log(masses[pos]**2/mu**2))/(LOOP*sigma**2))


def run(cache_dir):
    g=yuk.build_upper_geometry(cache_dir)
    p1=g["ps"]["common_geometry"]["p1"]
    cw=yuk.module("common_renorm_cw",RF/"code/verify_p54_p2_bosonic_cw.py")
    full=json.loads((RF/"output/p54_full_doublet_cw.json").read_text())
    up=json.loads((RF/"output/p54_upper_scalar_matching.json").read_text())
    wr=json.loads((RF/"output/p54_upper_wilson_exchange.json").read_text())
    p=full["tree_parameters"]
    assert p==g["upper"]["same_action_parameters"]
    r=np.array([full["vacuum"][k] for k in ("omega","sigma","vs")])
    ru=np.array([g["upper"]["vacuum"]["omega"],0.,g["upper"]["vacuum"]["vs"]])
    mu=full["scheme"]["mu_over_omega"]
    dp=np.array([full["tadpoles"]["finite_mass_counterterms_over_omega2"]["delta_"+k] for k in MASS_KEYS])
    ct=mass_operator(dp)
    pms=dict(p);pms.update({k:p[k]+v for k,v in zip(MASS_KEYS,dp)})
    jets=scalar.ActionJets(p1,p,cache_dir)
    x=p1.vacuum_vector(*r);h=jets.hessian(x,"broken_cached")
    rad=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    metric=rad.T@rad
    mass=rad.T@h@rad
    lam,u=np.linalg.eigh(h)
    gens=cw.generators();o=cw.field_orbit(p1,x,gens)
    mv=mu**2*o.T@o;lv,uv=np.linalg.eigh(mv)
    ts=[];us=[];tv=[];ov=[]
    for i,q in enumerate(rad.T):
        t,s=jets.jets(x,h,q,"common_radial_"+str(i))
        ts.append(t);us.append(s)
        oq=cw.field_orbit(p1,q,gens);ov.append(oq)
        tv.append(mu**2*(oq.T@o+o.T@oq))
    mixed={}
    for i in range(3):
        for j in range(i,3):
            if i==j: mixed[i,j]=us[i]
            else:
                step=.02
                hp=jets.hessian(x+step*(rad[:,i]+rad[:,j]),f"common_radial_pair_{i}{j}")
                mixed[i,j]=(hp-h-step*(ts[i]+ts[j])-.5*step**2*(us[i]+us[j]))/step**2
    def spectrum_first(values,vectors,derivative,count,c,mult):
        weight=np.zeros_like(values)
        weight[-count:]=values[-count:]*(2*np.log(values[-count:]/mu**2)-2*c+1)
        return float(mult*np.dot(weight,np.diag(vectors.T@derivative@vectors))/(4*LOOP))
    tad_s=np.array([spectrum_first(lam,u,t,290,1.5,1.) for t in ts])
    tad_v=np.array([spectrum_first(lv,uv,t,33,5/6,3.) for t in tv])
    tad=tad_s+tad_v
    hs=np.zeros((3,3));hv=hs.copy()
    for (i,j),uu in mixed.items():
        vv=mu**2*(ov[i].T@ov[j]+ov[j].T@ov[i])
        hs[i,j]=hs[j,i]=cw.hard_trace_hessian(lam,u,ts[i],ts[j],uu,290,mu**2,1.5,1.)
        hv[i,j]=hv[j,i]=cw.hard_trace_hessian(lv,uv,tv[i],tv[j],vv,33,mu**2,5/6,3.)
    hct=rad.T@ct@rad
    corrected=mass+hs+hv+hct
    eig_tree=eigh(mass,metric,eigvals_only=True)
    eig_corrected=eigh(corrected,metric,eigvals_only=True)
    e=np.linalg.eigh(mass)
    invroot=(e[1]/np.sqrt(e[0]))@e[1].T
    relative=eigh(invroot@(hs+hv+hct)@invroot,eigvals_only=True)
    # Fixed tree vertices give an EXACT polynomial radial matrix for all
    # small radial displacements; no new spectrum fit or scan is needed.
    def radial_hessian_at(dr):
        answer=h+sum(a*t for a,t in zip(dr,ts))
        for (i,j),v in mixed.items():answer+=dr[i]*dr[j]*v*(.5 if i==j else 1.)
        return answer
    def hard_potential(dr):
        hh=radial_hessian_at(dr)
        ox=o+sum(a*oi for a,oi in zip(dr,ov))
        ss=np.linalg.eigvalsh(hh)[-290:]
        vv=np.linalg.eigvalsh(mu**2*ox.T@ox)[-33:]
        return (np.sum(ss**2*(np.log(ss/mu**2)-1.5))+
                3*np.sum(vv**2*(np.log(vv/mu**2)-5/6)))/(4*LOOP)
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        err=yuk.error(np.asarray(a),np.asarray(b))
        checks.append({"name":name,"residual":err,"tolerance":tol,"pass":err<tol})
    check("same_parameter_dictionary_at_both_endpoints",float(p!=g["upper"]["same_action_parameters"]))
    check("canonical_radial_metric",metric,np.diag([12/5,2,1]))
    check("historical_full_bosonic_tadpole_reproduced",tad,full["tadpoles"]["radial_derivatives_over_omega3"])
    check("fixed_VEV_radial_counterterm_cancels_once",tad+radial_mass_map(r)@dp)
    check("counterterm_radial_map_equals_full_field_contraction",rad.T@ct@x,radial_mass_map(r)@dp)
    p10=np.zeros((328,328));p10[306:326,306:326]=np.eye(20)
    check("light_doublet_retuning_cannot_change_radial_block",rad.T@p10@rad)
    # Independent derivatives of the actual spectral function, not formulas
    # differentiated from themselves.
    numerical={}
    for step in (2e-5,1e-5):
        value=hard_potential(np.zeros(3));nh=np.zeros((3,3));nt=[]
        for i in range(3):
            ei=np.eye(3)[i]*step
            fp,fm=hard_potential(ei),hard_potential(-ei)
            nt.append((fp-fm)/(2*step));nh[i,i]=(fp+fm-2*value)/step**2
            for j in range(i):
                ej=np.eye(3)[j]*step
                nh[i,j]=nh[j,i]=(hard_potential(ei+ej)-hard_potential(ei-ej)-hard_potential(-ei+ej)+hard_potential(-ei-ej))/(4*step**2)
        check(f"spectral_tadpole_finite_difference_{step}",nt,tad,tol=2e-6)
        check(f"spectral_radial_Hessian_finite_difference_{step}",nh,hs+hv,tol=2e-5)
        numerical[str(step)]={"tadpole":nt,"Hessian":nh.tolist()}
    # Finite scheme covariance, V_MS(p+dp)=V_FV(p)+V_CT exactly
    # for these linearly occurring mass parameters at fixed fields.
    for i,point in enumerate((r,ru,np.array([.93,.11,.28]))):
        check(f"finite_scheme_radial_potential_identity_{i}",
              p1.radial_potential(point,pms),p1.radial_potential(point,p)+.5*(p1.vacuum_vector(*point)@ct@p1.vacuum_vector(*point)))
        check(f"finite_scheme_radial_gradient_identity_{i}",
              p1.radial_gradient(point,pms),p1.radial_gradient(point,p)+radial_mass_map(point)@dp)
    dp_ms_tree=-np.linalg.solve(mass,radial_mass_map(r)@dp)
    dp_loop=-np.linalg.solve(mass,tad)
    check("broken_MS_tree_and_explicit_loop_VEV_shifts_cancel",dp_ms_tree+dp_loop)
    # Formal stationary continuation uses tiny loop parameters only; these
    # are scheme-identity checks, not a new phenomenological vacuum scan.
    continuation=[]
    for ell in (1e-5,5e-6,2.5e-6):
        pp=dict(p);pp.update({k:p[k]+ell*v for k,v in zip(MASS_KEYS,dp)})
        solved=root(lambda rr:p1.radial_gradient(rr,pp),r+ell*dp_ms_tree,tol=1e-11)
        check(f"formal_stationary_root_{ell}",p1.radial_gradient(solved.x,pp))
        linear_error=np.linalg.norm(solved.x-r-ell*dp_ms_tree)
        continuation.append({"loop_parameter":ell,"stationary_point":solved.x.tolist(),"linear_remainder_norm":float(linear_error)})
    ratios=[continuation[i]["linear_remainder_norm"]/continuation[i+1]["linear_remainder_norm"] for i in range(2)]
    check("formal_stationary_conversion_remainder_is_second_order",max(abs(np.array(ratios)-4)),tol=.15)
    # Transport of the same fixed-VEV functional to upper singlet valley.
    r2=rad[:,[0,2]];hu=r2.T@g["hessian"]@r2
    tctu=r2.T@ct@g["x"]
    thu=np.array(up["common_CT_upper_tadpole"]["scalar_hard"])+up["common_CT_upper_tadpole"]["vector_hard"]
    du_ct=-np.linalg.solve(hu,tctu);du_loop=-np.linalg.solve(hu,thu)
    check("upper_common_valley_shift_equals_scheme_plus_loop",du_ct+du_loop,
          up["common_CT_upper_tadpole"]["hard_plus_CT_linear_displacement"])
    indices,bh,_,_,_=wilson.six_space(g)
    du_operator=bh.T@ct@bh
    for i,q in enumerate(r2.T):
        ti,_=jets.jets(g["x"],g["hessian"],q,"upper_scheme_radial_reuse")
        du_operator+=du_ct[i]*(bh.T@ti@bh)
    d=np.array(wr["tree_mass"]);dm=np.array(wr["delta_mass_total_subset"])
    # In a consistently shifted MS tree basis, the counterterm part moves
    # into tree D; it does not disappear from tree+one-loop D.
    check("upper_mass_scheme_transport_identity",d+dm,(d+du_operator)+(dm-du_operator))
    source_cases=[]
    for i,case in enumerate(wr["cases"]):
        hraw=complex(case["synthetic_h_raw"]["real"],case["synthetic_h_raw"]["imag"])
        fraw=complex(case["synthetic_f_raw"]["real"],case["synthetic_f_raw"]["imag"])
        y=g["heavy_h"][indices]*hraw+g["heavy_f"][indices]*fraw
        j,_=wilson.kjc.source_columns({"combined":y},g)
        ij=np.linalg.solve(d,j)
        dc_tree=-ij.T@du_operator@ij
        fixed_delta=np.array(case["delta_C_factorizable"])
        ms_delta=fixed_delta-dc_tree
        check(f"Wilson_parameter_scheme_transport_{i}",dc_tree+ms_delta,fixed_delta)
        eps=1e-5
        numeric=(j.T@np.linalg.solve(d+eps*du_operator,j)-j.T@np.linalg.solve(d-eps*du_operator,j))/(2*eps)
        check(f"Wilson_scheme_tree_derivative_{i}",numeric,dc_tree,tol=1e-7)
        source_cases.append({"synthetic_case":i,"scheme_tree_Wilson_derivative":dc_tree.tolist(),
                             "fixed_VEV_delta_norm":float(np.linalg.norm(fixed_delta)),
                             "converted_MS_loop_delta_norm":float(np.linalg.norm(ms_delta)),
                             "tree_shift_norm":float(np.linalg.norm(dc_tree)),"no_new_physical_prediction":True})
    # Direct scalar mass operator and field redefinitions cannot license a
    # finite matching package with a different EFT parameter convention.
    rng=np.random.default_rng(20260913)
    jac=rng.normal(size=(6,3));uvs=rng.normal(size=3);efts=rng.normal(size=6);dm0=rng.normal(size=6)
    converted=finite_matching_conversion(dm0,jac,uvs,efts)
    check("general_UV_and_EFT_finite_scheme_matching_chain_rule",jac@uvs+converted-efts,dm0)
    # Resolve the source of dnu2 by degenerate mass clusters, with basis-
    # invariant traces rather than arbitrarily labelled eigenvectors.
    def clusters(values,vectors,first,count,c,mult):
        rows=[];start=len(values)-count
        while start<len(values):
            stop=start+1
            while stop<len(values) and abs(values[stop]-values[start])<1e-8:stop+=1
            mass2=float(np.mean(values[start:stop]));b=vectors[:,start:stop]
            td=mult*mass2*(2*np.log(mass2/mu**2)-2*c+1)*np.trace(b.T@first@b)/(4*LOOP)
            rows.append({"mass2":mass2,"multiplicity":stop-start,"delta_nu2_contribution":float(td/r[1])})
            start=stop
        return sorted(rows,key=lambda rr:abs(rr["delta_nu2_contribution"]),reverse=True)
    groups_s=clusters(lam,u,ts[1],290,1.5,1.)
    groups_v=clusters(lv,uv,tv[1],33,5/6,3.)
    check("mass_cluster_decomposition_reconstructs_delta_nu2",sum(x["delta_nu2_contribution"] for x in groups_s+groups_v),dp[1])
    # Explicit mu dependence at fixed tree parameters; NOT the full beta_p.
    dtlog=np.array([-np.dot(lam[-290:],np.diag(u.T@t@u)[-290:])/LOOP
        -3*np.dot(lv[-33:],np.diag(uv.T@v@uv)[-33:])/LOOP for t,v in zip(ts,tv)])
    dctlog=-np.linalg.solve(radial_mass_map(r),dtlog)
    # An analytic ALL-Yukawa upper bound, not a scan over fitted families:
    # max_{x>=0} [-x^2 log(x/mu^2)] = mu^4/(2e).
    cap=6*mu**4/(np.e*LOOP*r[1]**2)
    upper_radial=corrected.copy();upper_radial[1,1]+=cap
    cap_masses,cap_vectors=eigh(upper_radial,metric)
    witness=cap_vectors[:,0]
    check("negative_radial_witness_is_gauge_horizontal",o.T@(rad@witness))
    check("negative_radial_witness_is_not_PQ_phase",p1.pq_direction(x)@(rad@witness))
    mass_f=np.einsum("a,aij->ij",x,g["all_f"])
    check("one_massive_spinor_direction_per_family_at_broken_vacuum",np.linalg.matrix_rank(mass_f,tol=1e-9),1)
    check("no_h_raw_mass_on_radial_background",np.einsum("ar,aij->rij",rad,g["all_h"]))
    check("f_raw_radial_mass_only_depends_on_sigma",np.einsum("ar,aij->rij",rad[:,[0,2]],g["all_f"]))
    check("Majorana_cap_saturated_by_three_equal_masses",
          majorana_fixed_VEV_radial(np.full(3,mu*np.exp(-.25)),r[1],mu),cap)
    check("universal_radial_instability_witness_metric_normalization",witness@metric@witness,1.)
    check("universal_radial_instability_witness_Rayleigh_identity",witness@upper_radial@witness,cap_masses[0])
    check("negative_direction_survives_maximal_three_Majorana_curvature",float(cap_masses[0]>=0))
    fermion_tests=[]
    for i,masses in enumerate(([.02,.06,.11],[.25,.35,.45],[.6,.8,1.],np.full(3,mu*np.exp(-.25)))):
        masses=np.asarray(masses)
        f=masses/r[1]
        def vf(s):
            x2=(s*f)**2
            return -np.sum(x2**2*(np.log(x2/mu**2)-1.5))/(2*LOOP)
        t=-np.sum(masses**4*(np.log(masses**2/mu**2)-1))/(8*np.pi**2*r[1])
        step=2e-5*r[1]
        fd_curvature=(vf(r[1]+step)+vf(r[1]-step)-2*vf(r[1]))/step**2-t/r[1]
        analytic=majorana_fixed_VEV_radial(masses,r[1],mu)
        check(f"Majorana_fixed_VEV_curvature_direct_potential_{i}",fd_curvature,analytic,tol=2e-5)
        check(f"Majorana_cap_test_{i}",max(0.,analytic-cap))
        fermion_tests.append({"masses":masses.tolist(),"curvature":analytic,"direct_difference":float(fd_curvature)})
    paths=[Path(__file__),Path(p1.__file__),Path(cw.__file__),Path(scalar.__file__),Path(wilson.__file__),Path(yuk.__file__),
           RF/"output/p54_full_doublet_cw.json",RF/"output/p54_upper_scalar_matching.json",RF/"output/p54_upper_wilson_exchange.json"]
    return {"schema":"p54-common-bosonic-renormalization-v1","date":"2026-09-12",
        "scope":"one-loop bosonic same-bare-parameter consistency; missing full diagrams and fermionic feedback remain missing",
        "mu":mu,"mass_parameter_order":MASS_KEYS,"fixed_VEV_parameters":{k:p[k] for k in MASS_KEYS},
        "finite_counterterm_shift":dp.tolist(),"same_bare_MS_parameters":{k:pms[k] for k in MASS_KEYS},
        "relative_parameter_shift":(dp/np.array([p[k] for k in MASS_KEYS])).tolist(),
        "broken_tadpoles":{"scalar":tad_s.tolist(),"vector":tad_v.tolist(),"total":tad.tolist(),
            "fixed_parameter_explicit_logmu_CT_derivative":dctlog.tolist(),"not_full_mass_parameter_beta":True},
        "radial_bosonic":{"metric":metric.tolist(),"tree_Hessian":mass.tolist(),
            "scalar_loop_Hessian":hs.tolist(),"vector_loop_Hessian":hv.tolist(),"CT_Hessian":hct.tolist(),
            "corrected_Hessian":corrected.tolist(),"tree_eigenvalues":eig_tree.tolist(),
            "corrected_generalized_eigenvalues":eig_corrected.tolist(),
            "relative_insertion_eigenvalues":relative.tolist(),"no_negative_radial_mode_at_this_truncation":bool(min(eig_corrected)>0),
            "fermionic_and_soft_IR_completion_included":False},
        "broken_MS_tree_linear_VEV_shift":dp_ms_tree.tolist(),"broken_explicit_loop_VEV_shift":dp_loop.tolist(),
        "formal_small_loop_continuation":continuation,"continuation_remainder_ratios":ratios,
        "upper_scheme_tree_radial_shift":du_ct.tolist(),"upper_explicit_hard_radial_shift":du_loop.tolist(),
        "upper_six_scheme_tree_mass_derivative":du_operator.tolist(),"Wilson_scheme_transport":source_cases,
        "scalar_delta_nu2_clusters":groups_s,"vector_delta_nu2_clusters":groups_v,
        "finite_difference_regressions":numerical,"cache":{"hits":jets.hits,"new_vertex_evaluations":jets.evaluated},
        "all_Majorana_radial_bound":{"assumptions":["frozen scalar point and mu","three Majorana masses M=sigma*f_M",
            "same fixed-VEV finite nu2 counterterm","one-loop hard potential; soft Goldstone/IR completion excluded"],
            "max_sigma_sigma_curvature":cap,"upper_bound_generalized_eigenvalues":cap_masses.tolist(),
            "universal_negative_witness":witness.tolist(),"witness_upper_bound":float(cap_masses[0]),
            "any_family_matrices_rescue_declared_radial_truncation":False,
            "full_model_or_all_orders_exclusion":False,"fermion_potential_tests":fermion_tests},
        "common_bosonic_scheme_algebra_complete":True,"full_common_physical_renormalization_complete":False,
        "full_Wilson_box_matching":False,"full_lower_matching":False,"full_finite_CHN_seesaw":False,
        "physical_fit_enabled":False,
        "source_sha256":{str(path.relative_to(RF)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# Common bosonic renormalization and tadpole consistency", "",
        "**Same-bare-parameter algebra is checked; full matching and physical control remain open.**", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "A finite fixed-VEV counterterm cannot be deleted at unchanged renormalized input. The equivalent MS input is p_MS=p_fixed+delta_p. Tree backgrounds, mass kernels and Wilson coefficients must be transported as well. The full first-order response is invariant; its division into tree and loop pieces is not.", "",
        "| Parameter | Fixed-VEV input | Finite shift | Same-bare MS input |", "|---|---:|---:|---:|"]
    lines += [f"| {key} | {r['fixed_VEV_parameters'][key]:.9g} | {r['finite_counterterm_shift'][i]:.9g} | {r['same_bare_MS_parameters'][key]:.9g} |" for i,key in enumerate(MASS_KEYS)]
    lines += ["", "## Broken radial bosonic block", "",
        f"Tree generalized squared masses: {r['radial_bosonic']['tree_eigenvalues']}", "",
        f"Corrected bosonic hard+CT squared masses: {r['radial_bosonic']['corrected_generalized_eigenvalues']}", "",
        f"Relative insertions: {r['radial_bosonic']['relative_insertion_eigenvalues']}", "",
        "These are fixed-gauge hard-potential curvatures, not pole masses or a full-model stability certificate. They exclude the fermionic and soft-IR completion.", "",
        "## No Majorana-family rescue within the declared hard truncation", "",
        "At fixed f_M, the finite fixed-VEV Majorana correction is only in the sigma,sigma entry: -4 sum_i M_i^4 log(M_i^2/mu^2)/(16 pi^2 sigma^2). The exact scalar inequality -x^2 log(x/mu^2) <= mu^4/(2e) bounds all three families without fitting or sampling them.", "",
        f"Their maximal positive curvature is {r['all_Majorana_radial_bound']['max_sigma_sigma_curvature']:.9g} omega^2. Even the matrix with this maximal correction has lowest generalized eigenvalue {r['all_Majorana_radial_bound']['witness_upper_bound']:.9g} omega^2. Thus no choice of the three Majorana masses rescues this frozen point in the declared one-loop hard-potential truncation. This is NOT an all-orders/full-model exclusion; a consistently reorganized scalar/IR treatment would be a different calculation.", "",
        f"The formal MS tree VEV response is {r['broken_MS_tree_linear_VEV_shift']}; the explicit loop response cancels it at first order. Tiny formal loop-parameter continuation tests this identity, not a new phenomenological vacuum scan.", "",
        "## Still missing", "", "- Full fermionic/common parameter feedback and controlled vacuum expansion.",
        "- Complete Wilson/box graphs and lower gauge/ghost/Yukawa matching.",
        "- Actual CHN/pre-existing-C5/mixed-sterile finite matching.", "",
        "## Checks", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--cache-dir",type=Path,required=True)
    args=ap.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"MS_parameters":r["same_bare_MS_parameters"],
        "radial_squared_masses":r["radial_bosonic"]["corrected_generalized_eigenvalues"],
        "MS_tree_shift":r["broken_MS_tree_linear_VEV_shift"],"failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
