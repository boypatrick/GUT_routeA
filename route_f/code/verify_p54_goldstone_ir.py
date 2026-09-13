#!/usr/bin/env python3
"""Same-action Goldstone Ward, soft-IR and scalar momentum audit.

Actual polynomial radial vertices; no parameter changes or new fit. Scalar
Euclidean momentum dependence is complete at one loop on this radial slice.
The vector contribution remains the old zero-momentum hard potential;
therefore the combined diagnostic is NOT a complete pole or stability test.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

import verify_p54_common_renormalization as common

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_goldstone_ir"
LOOP=16*np.pi**2


@lru_cache(maxsize=8)
def quadrature(n):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)/2,w/2


def log_bubble(a,b,s,mu2,n=128):
    """Euclidean MS finite log integral, -B0; arguments are squared masses.

    Analytic zero/one-massless cases remove endpoint logarithms. Gaussian
    integration is used only for two strictly positive masses.
    """
    if min(a,b,s)<0 or mu2<=0:raise ValueError("Nonnegative masses/momentum required")
    if a==0 and b==0:
        if s==0:raise ValueError("Massless zero-momentum bubble is IR divergent")
        return float(np.log(s/mu2)-2)
    if a==0 or b==0:
        m=max(a,b)
        if s==0:return float(np.log(m/mu2)-1)
        z=s/m
        # int log[t*(m+(1-t)s)] = log(m)-2+(1+1/z)log(1+z)
        return float(np.log(m/mu2)-2+(1+1/z)*np.log1p(z))
    if s==0:
        if abs(a-b)<1e-10*max(a,b):return float(np.log((a+b)/(2*mu2)))
        return float((a*np.log(a/mu2)-b*np.log(b/mu2))/(a-b)-1)
    x,w=quadrature(n)
    return float(w@np.log(((1-x)*a+x*b+x*(1-x)*s)/mu2))


def run(cache_dir):
    full_path=RF/"output/p54_full_doublet_cw.json"
    ren_path=RF/"output/p54_common_renormalization.json"
    full=json.loads(full_path.read_text());ren=json.loads(ren_path.read_text())
    p1=common.yuk.module("goldstone_p1",RF/"code/verify_p54_p1_hessian_spectrum.py")
    cw=common.yuk.module("goldstone_cw",RF/"code/verify_p54_p2_bosonic_cw.py")
    p=full["tree_parameters"];r=np.array([full["vacuum"][k] for k in ("omega","sigma","vs")])
    mu2=ren["mu"]**2
    store=common.scalar.ActionJets(p1,p,cache_dir)
    x=p1.vacuum_vector(*r);h=store.hessian(x,"same_broken_background")
    radial=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    metric=radial.T@radial
    lam,u=np.linalg.eigh(h);soft=u[:,:38];heavy=u[:,38:]
    t=[];uu={}
    for i,q in enumerate(radial.T):
        first,second=store.jets(x,h,q,"ir_radial_"+str(i))
        t.append(first);uu[i,i]=second
    for i in range(3):
        for j in range(i+1,3):
            step=.02
            hp=store.hessian(x+step*(radial[:,i]+radial[:,j]),f"ir_pair_{i}{j}")
            uu[i,j]=uu[j,i]=(hp-h-step*(t[i]+t[j])-.5*step**2*(uu[i,i]+uu[j,j]))/step**2
    t=np.array(t)
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        a,b=np.asarray(a),np.asarray(b)
        error=float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append({"name":name,"residual":error,"tolerance":tol,"pass":error<tol})
    check("frozen_tree_stationarity",p1.radial_gradient(r,p))
    check("38_tree_soft_modes",lam[:38])
    check("positive_hard_gap",float(lam[38]<=0))
    gens=cw.generators();orbit=cw.field_orbit(p1,x,gens)
    ug,sv,_=np.linalg.svd(orbit,full_matrices=False);gauge=ug[:,sv>1e-9]
    pq=p1.pq_direction(x);pq-=gauge@(gauge.T@pq);pq=pq/np.linalg.norm(pq)
    gold=np.column_stack([gauge,pq])
    rem=soft-gold@(gold.T@soft)
    uh,sh,_=np.linalg.svd(rem,full_matrices=False);higgs=uh[:,sh>1e-9]
    check("33_gauge_1_PQ_4_light_dimensions",[gauge.shape[1],np.linalg.matrix_rank(pq[:,None]),higgs.shape[1]],[33,1,4])
    check("actual_soft_projector_decomposition",soft@soft.T,gold@gold.T+higgs@higgs.T)
    check("goldstone_tree_zero_modes",h@gold)
    # Ward identity for the invariant hard functional: H1 T x = T grad V1.
    tad=np.array(ren["broken_tadpoles"]["total"])
    grad=radial@np.linalg.solve(metric,tad)
    b=np.column_stack([orbit,p1.pq_direction(x)])
    tg=np.column_stack([cw.field_orbit(p1,grad,gens),p1.pq_direction(grad)])
    ub,sb,vb=np.linalg.svd(b,full_matrices=False);keep=sb>1e-9
    transform=vb[keep,:].T/sb[keep]
    pred=transform.T@b.T@tg@transform
    ct=common.mass_operator(np.array(ren["finite_counterterm_shift"]))
    gold_ct=ub[:,keep].T@ct@ub[:,keep]
    check("hard_Ward_mass_matrix_symmetric",pred,pred.T)
    check("common_CT_cancels_all_34_hard_Ward_masses",pred+gold_ct)
    for i in range(3):
        th=np.column_stack([cw.field_orbit(p1,h@radial[:,i],gens),p1.pq_direction(h@radial[:,i])])
        check(f"actual_cubic_Goldstone_Ward_{i}",soft.T@t[i]@b,soft.T@th)
    _,_,vg=np.linalg.svd(orbit,full_matrices=False);unbroken=vg[33:,:].T
    for i in range(3):
        check(f"unbroken_vectors_stay_massless_along_radial_{i}",cw.field_orbit(p1,radial[:,i],gens)@unbroken)
    # Actual Gram coefficient of the IR logarithm, including all soft modes.
    blocks={"gauge":gauge,"PQ":pq[:,None],"Goldstone":gold,"Higgs":higgs,"all_soft":soft}
    grams={};a_soft=None
    for key,q in blocks.items():
        a=np.array([q.T@v@q for v in t])
        grams[key]=np.einsum("rij,sji->rs",a,a)/(2*LOOP)
        check(f"{key}_IR_Gram_positive",max(0.,-np.linalg.eigvalsh(grams[key]).min()))
        if key=="all_soft":a_soft=a
    # Gauge and physical PQ can mix in a chosen orthogonal basis; retain it.
    grams["gauge_PQ_mixed"]=grams["Goldstone"]-grams["gauge"]-grams["PQ"]
    check("Goldstone_Higgs_radial_mixing_vanishes",np.array([gold.T@v@higgs for v in t]))
    check("all_soft_IR_Gram_equals_Goldstone_plus_Higgs",grams["all_soft"],grams["Goldstone"]+grams["Higgs"])
    e_soft=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            ti=soft.T@t[i]@heavy;tj=soft.T@t[j]@heavy
            e_soft[i,j]=(np.trace(soft.T@uu[i,j]@soft)-2*np.sum(ti*tj/lam[38:]))/(2*LOOP)
    trace_a=np.trace(a_soft,axis1=1,axis2=2)/(2*LOOP)
    def regulated_soft(eps):
        logarithm=np.log(eps/mu2)
        return grams["all_soft"]*logarithm+eps*(logarithm-1)*e_soft
    regulators=[]
    q=np.array(ren["all_Majorana_radial_bound"]["universal_negative_witness"])
    hb=np.array(ren["radial_bosonic"]["corrected_Hessian"])
    for eps in (1e-4,1e-6,1e-8,1e-10):
        soft_h=regulated_soft(eps);soft_t=eps*(np.log(eps/mu2)-1)*trace_a
        regulators.append({"regulator_mass2":eps,"soft_tadpole":soft_t.tolist(),
            "soft_Hessian":soft_h.tolist(),"witness_soft_curvature":float(q@soft_h@q),
            "static_diagnostic_eigenvalues":eigh(hb+soft_h,metric,eigvals_only=True).tolist()})
    check("soft_tadpole_vanishes_with_regulator",float(np.linalg.norm(regulators[-1]["soft_tadpole"])>=1e-7))
    check("IR_log_coefficient_nonzero_on_old_negative_witness",float(q@grams["all_soft"]@q<=1e-10))
    for eps in (1e-6,1e-8):
        step=.01
        slope=(regulated_soft(eps*np.exp(step))-regulated_soft(eps*np.exp(-step)))/(2*step)
        expected=grams["all_soft"]+eps*np.log(eps/mu2)*e_soft
        check(f"regulated_Hessian_log_slope_{eps}",slope,expected,tol=1e-7)
    # Direct spectral differentiation of the local soft eigenvalue branches.
    def soft_potential(dr,eps):
        hh=h+np.einsum("i,ijk->jk",dr,t)
        for i in range(3):
            for j in range(3):hh=hh+.5*dr[i]*dr[j]*uu[i,j]
        masses=np.linalg.eigvalsh(hh)[:38]+eps
        if np.any(masses<=0):raise ValueError("Finite-difference point left positive regulator neighborhood")
        return np.sum(masses**2*(np.log(masses/mu2)-1.5))/(4*LOOP)
    eps=1e-3;step=2e-6;fd=np.zeros((3,3));f0=soft_potential(np.zeros(3),eps)
    for i in range(3):
        ei=np.eye(3)[i]*step
        fd[i,i]=(soft_potential(ei,eps)+soft_potential(-ei,eps)-2*f0)/step**2
        for j in range(i):
            ej=np.eye(3)[j]*step
            fd[i,j]=fd[j,i]=(soft_potential(ei+ej,eps)-soft_potential(ei-ej,eps)-soft_potential(-ei+ej,eps)+soft_potential(-ei-ej,eps))/(4*step**2)
    check("soft_regulated_Hessian_independent_spectral_difference",fd,regulated_soft(eps),tol=2e-5)
    # All scalar one-loop bubbles at Euclidean momentum, grouped by actual
    # degenerate tree eigenvalues. The zero block includes the four Higgs modes.
    groups=[np.arange(38)];masses=[0.];start=38
    while start<328:
        stop=start+1
        while stop<328 and abs(lam[stop]-lam[start])<1e-8:stop+=1
        groups.append(np.arange(start,stop));masses.append(float(np.mean(lam[start:stop])));start=stop
    te=np.array([u.T@v@u for v in t]);ng=len(groups)
    weight=np.empty((3,3,ng,ng));seagull=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            for k,gg in enumerate(groups):
                for l,hh in enumerate(groups):
                    weight[i,j,k,l]=np.sum(te[i][np.ix_(gg,hh)]*te[j][np.ix_(gg,hh)])/(2*LOOP)
            seagull[i,j]=np.dot(lam[38:]*(np.log(lam[38:]/mu2)-1),np.diag(heavy.T@uu[i,j]@heavy))/(2*LOOP)
    logs0=np.array([[log_bubble(a,b,0,mu2) if i+j else 0. for j,b in enumerate(masses)] for i,a in enumerate(masses)])
    scalar0=seagull+np.einsum("rsij,ij->rs",weight,logs0)
    check("all_scalar_zero_momentum_hard_part_reproduces_old_Frechet",scalar0,ren["radial_bosonic"]["scalar_loop_Hessian"])
    check("grouped_soft_bubble_equals_actual_IR_Gram",weight[:,:,0,0],grams["all_soft"])
    for k,(a,b,s) in enumerate(((0.,0.,.03),(0.,.02,.04),(.02,.02,.07),(.013,.19,.11))):
        direct=quad(lambda z:np.log(((1-z)*a+z*b+z*(1-z)*s)/mu2),0,1,epsabs=1e-12)[0]
        check(f"bubble_independent_integral_{k}",log_bubble(a,b,s,mu2),direct,tol=1e-10)
    # At s>0 the regulator log cancels before taking eps -> 0.
    s_test=.02
    for eps in (1e-4,1e-7,1e-10):
        l=log_bubble(eps,eps,s_test,mu2,n=256)
        assembled=regulated_soft(eps)+grams["all_soft"]*(l-np.log(eps/mu2))
        check(f"momentum_IR_cancellation_identity_{eps}",assembled,
              grams["all_soft"]*l+eps*(np.log(eps/mu2)-1)*e_soft)
    def momentum(s,n=128):
        logs=np.array([[log_bubble(a,b,s,mu2,n) for b in masses] for a in masses])
        diff=logs-logs0;diff[0,0]=0.
        hard_motion=np.einsum("rsij,ij->rs",weight,diff)
        soft_motion=grams["all_soft"]*(np.log(s/mu2)-2)
        return hard_motion,soft_motion,hb+s*metric+hard_motion+soft_motion
    momenta=[]
    # These are kinematic probes, not scalar-parameter or vacuum scans.
    for s in (.001,.01,.05,.1,.5,1.):
        dh,ss,k=momentum(s)
        check(f"scalar_hard_momentum_difference_positive_Gram_{s}",max(0.,-np.linalg.eigvalsh(dh).min()))
        momenta.append({"pE2":s,"hard_scalar_momentum_difference":dh.tolist(),"soft_scalar_bubble":ss.tolist(),
            "radial_Euclidean_subset":k.tolist(),"generalized_eigenvalues":eigh(k,metric,eigvals_only=True).tolist(),
            "fixed_witness":float(q@k@q),"is_full_gauged_pole_kernel":False})
    dh,ss,k=momentum(.05,64);dh2,_,k2=momentum(.05,256)
    check("scalar_momentum_quadrature_64_vs_256",k,k2,tol=1e-9)
    # Finite hard scalar wave-function coefficient; soft-soft remains
    # nonlocal and cannot be represented by this derivative expansion.
    slopes=np.zeros((ng,ng))
    for i,a in enumerate(masses):
        for j,bm in enumerate(masses):
            if i+j:
                slopes[i,j]=quad(lambda z:z*(1-z)/((1-z)*a+z*bm),0,1,epsabs=1e-12)[0]
    dz=np.einsum("rsij,ij->rs",weight,slopes)
    dz_eig=eigh(dz,metric,eigvals_only=True)
    for small_s in (1e-5,5e-6):
        dd,_,_=momentum(small_s)
        check(f"hard_scalar_kinetic_matches_momentum_derivative_{small_s}",dd/small_s,dz,tol=2e-4)
    check("hard_scalar_kinetic_positive",max(0.,-dz_eig.min()))
    # A uniform rescaling is an exact family of EXISTING scalar parameters,
    # not a new operator, new free family matrix or fitted mass formula.
    # V0 -> z V0 preserves the tree stationary radii and all tree inertia.
    # Only inspect this repair path if the old hard scalar derivative
    # expansion already fails its small-insertion control criterion.
    hs=np.array(ren["radial_bosonic"]["scalar_loop_Hessian"])
    hv=np.array(ren["radial_bosonic"]["vector_loop_Hessian"])
    tree=np.array(ren["radial_bosonic"]["tree_Hessian"])
    def radial_ct(tadpole):
        shift=-np.linalg.solve(common.radial_mass_map(r),tadpole)
        return radial.T@common.mass_operator(shift)@radial
    scalar_fixed=hs+radial_ct(np.array(ren["broken_tadpoles"]["scalar"]))
    vector_fixed=hv+radial_ct(np.array(ren["broken_tadpoles"]["vector"]))
    log_tad=np.array([np.dot(lam[38:],np.diag(heavy.T@v@heavy)) for v in t])/(2*LOOP)
    log_hess=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            log_hess[i,j]=(np.dot(lam[38:],np.diag(heavy.T@uu[i,j]@heavy))+
                np.sum(te[i]*te[j])-np.sum(te[i,:38,:38]*te[j,:38,:38]))/(2*LOOP)
    scalar_log=log_hess+radial_ct(log_tad)
    check("scalar_vector_CT_split_reconstructs_old_radial",tree+scalar_fixed+vector_fixed,hb)
    # For 0<z<=1, the nonnegative L_sigma,sigma times log(z) can be
    # dropped in an upper bound. The remaining concave quadratic has a
    # negative global maximum, excluding this entire bosonic repair ray.
    bound_uniform=vector_fixed[1,1]+tree[1,1]**2/(-4*scalar_fixed[1,1])
    check("uniform_ray_bound_sign_assumptions",float(scalar_fixed[1,1]>=0 or scalar_log[1,1]<0))
    check("entire_bosonic_weakening_ray_has_negative_sigma_witness",float(bound_uniform>=0))
    def scaled_hard(z):return z*tree+z*z*(scalar_fixed+np.log(z)*scalar_log)+vector_fixed
    scale_rows=[]
    if dz_eig.max()>1:
        for z in (.01,.03,.1,.3):
            hz=scaled_hard(z)
            scale_rows.append({"scalar_parameter_multiplier":z,
                "radial_hard_eigenvalues":eigh(hz,metric,eigvals_only=True).tolist(),
                "relative_radial_insertions":eigh(hz-z*tree,z*tree,eigvals_only=True).tolist(),
                "hard_scalar_kinetic_eigenvalues":(z*dz_eig).tolist(),
                "full_background_candidate_certified":False})
            pz={key:z*val for key,val in p.items()}
            check(f"uniform_scalar_rescaling_stationarity_{z}",p1.radial_gradient(r,pz))
            for rr in (r,np.array([.93,.11,.28])):
                check(f"uniform_scalar_rescaling_exact_action_radial_{z}_{rr[0]}",p1.radial_potential(rr,pz),z*p1.radial_potential(rr,p))
            # Direct use of rescaled matrices in the old Frechet functional,
            # independent of the derived z^2 log(z) polynomial above.
            hz_s=np.zeros((3,3));tz=[]
            for i in range(3):
                lz=z*lam[38:]
                tz.append(np.dot(lz*(np.log(lz/mu2)-1),z*np.diag(heavy.T@t[i]@heavy))/(2*LOOP))
                for j in range(3):
                    hz_s[i,j]=cw.hard_trace_hessian(z*lam,u,z*t[i],z*t[j],z*uu[i,j],290,mu2,1.5,1.)
            check(f"uniform_scalar_rescaling_Frechet_identity_{z}",hz_s+radial_ct(np.array(tz)),z*z*(scalar_fixed+np.log(z)*scalar_log))
    # Since q^T B q >=0, a low-momentum soft scalar bubble cannot be
    # a positive rescue for s <= exp(2)*mu^2 in this MS convention.
    check("soft_bubble_negative_semidefinite_below_exp2_mu2",max(0.,np.linalg.eigvalsh(ss).max()))
    # Perturbative control is tested separately from a pole-mass declaration.
    p_shift=np.array(ren["finite_counterterm_shift"])
    h_shift=common.mass_operator(p_shift)
    ratios=[]
    for gg,m2 in zip(groups[1:],masses[1:]):
        ev=np.linalg.eigvalsh(u[:,gg].T@h_shift@u[:,gg])/m2
        ratios.append({"tree_mass2":m2,"multiplicity":len(gg),"CT_only_relative_eigenvalues":ev.tolist()})
    sources=[Path(__file__),Path(common.__file__),Path(p1.__file__),Path(cw.__file__),Path(common.scalar.__file__),full_path,ren_path]
    return {"schema":"p54-goldstone-ir-radial-v1","date":"2026-09-13",
        "scope":"one-loop same-action Goldstone Ward and scalar momentum slice; full gauge/fermion momentum and two-loop reorganization absent",
        "mu2":mu2,"radial_metric":metric.tolist(),"soft_dimensions":{"gauge":33,"physical_PQ":1,"light_Higgs":4},
        "soft_gap":{"max_soft_abs_mass2":float(np.max(abs(lam[:38]))),"min_hard_mass2":float(lam[38])},
        "Goldstone_Ward":{"hard_projected_mass":pred.tolist(),"common_CT_projected_mass":gold_ct.tolist(),
            "net_norm":float(np.linalg.norm(pred+gold_ct)),"bare_hard_mass_eigenvalues":np.linalg.eigvalsh(pred).tolist(),
            "extra_positive_Goldstone_regulator_is_not_a_scheme_preserving_repair":True},
        "IR_Gram_matrices":{key:value.tolist() for key,value in grams.items()},"soft_reduced_second_trace":e_soft.tolist(),
        "negative_witness_IR_coefficients":{key:float(q@value@q) for key,value in grams.items()},
        "regulator_audit":regulators,"scalar_momentum_probes":momenta,"scalar_tree_mass_clusters":masses,
        "hard_scalar_kinetic":{"matrix":dz.tolist(),"generalized_eigenvalues":dz_eig.tolist(),
            "perturbation_smaller_than_tree_metric":bool(dz_eig.max()<1),"soft_nonlocal_part_excluded":True},
        "uniform_scalar_parameter_repair":{"triggered_by_large_hard_kinetic":bool(dz_eig.max()>1),
            "formula":"H_hard(z)=z H_tree+z^2(S+log(z) L_S)+V; all existing scalar coefficients scaled, gauge held fixed",
            "scalar_fixed_matrix":scalar_fixed.tolist(),"vector_fixed_matrix":vector_fixed.tolist(),
            "scalar_log_matrix":scalar_log.tolist(),"diagnostics":scale_rows,
            "all_0_to_1_sigma_curvature_upper_bound":float(bound_uniform),
            "uniform_weakening_bosonic_ray_excluded":bool(bound_uniform<0),
            "does_not_exclude_added_fermion_corrections_at_new_points":True,
            "updates_default_inputs":False,"changes_gauge_coupling":False},
        "CT_only_mass_ratio_diagnostics":ratios,
        "soft_tadpole_limit_zero":True,"zero_momentum_full_Hessian_is_IR_finite":False,
        "scalar_Goldstone_IR_alone_repairs_hard_benchmark":False,
        "full_momentum_gauge_ghost_fermion_kernel_complete":False,"two_loop_Goldstone_resummation_complete":False,
        "full_background_stability_decided":False,"new_scalar_point_selected":False,"physical_fit_enabled":False,
        "missing":["same-prescription vector/Goldstone/ghost mixed and fermionic momentum maps with Nielsen consistency",
            "controlled scalar parameter feedback or an explicitly reorganized perturbative counting",
            "complete Wilson/box, lower gauge/ghost/Yukawa matching and finite CHN seesaw"],
        "cache":{"hits":store.hits,"new_Hessian_evaluations":store.evaluated},
        "source_sha256":{str(path.relative_to(RF)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# P54 same-action Goldstone / IR background audit","",
        "**The one-loop scalar IR and momentum subgate is computed; full gauged background stability is NOT decided.**","",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.","",
        "All 34 Goldstone hard mass shifts are canceled by the existing common fixed-VEV counterterms. Adding a positive Goldstone mass by hand is not a same-scheme repair. Actual soft tadpoles vanish as the regulator tends to zero, but the radial potential Hessian has a nonzero positive-Gram coefficient multiplying log(epsilon/mu^2). It is not a finite physical mass matrix.","",
        "The full scalar one-loop momentum dependence on the radial slice is computed from the actual cubic/quartic action tensors. External Euclidean momentum replaces the divergent log by log(pE^2/mu^2)-2 for the massless bubble. The hard scalar momentum difference is positive semidefinite, but this is not the full gauged pole kernel.","",
        "| pE^2 / omega^2 | Lowest generalized eigenvalue of diagnostic kernel |","|---:|---:|"]
    for row in r["scalar_momentum_probes"]:lines.append(f"| {row['pE2']:.4g} | {row['generalized_eigenvalues'][0]:.9g} |")
    lines += ["","The diagnostic includes the old zero-momentum vector hard term, not its dynamical gauge/ghost completion. Its negative values cannot be promoted to a physical tachyon; positive values at high Euclidean momentum cannot establish a stable vacuum either. No default scalar parameters or old matching trajectories were changed.","",
        "## Background control and a ruled-out repair ray","",
        f"The hard scalar kinetic correction has generalized eigenvalues {r['hard_scalar_kinetic']['generalized_eigenvalues']}. The largest exceeds one in the declared MS convention. A finite canonical field redefinition can absorb this metric, but does not establish a controlled loop remainder or change the inertia of a fixed curvature matrix.","",
        "For the exact existing-parameter ray V0 -> z V0 at fixed radii/gauge/scale, H_hard(z)=z H_tree+z^2(S+log(z) L_S)+V and deltaZ_h(z)=z deltaZ_h(1). The sigma coordinate gives an analytic concave-quadratic upper bound for every 0<z<=1:","",
        f"`H_sigma,sigma(z) <= {r['uniform_scalar_parameter_repair']['all_0_to_1_sigma_curvature_upper_bound']:.9g} omega^2 < 0`.","",
        "This excludes that entire BOSONIC repair ray, not arbitrary scalar parameters or fermionic corrections at new points. The separate scalar_reselection report tests a few less-hierarchical engineering points; none overwrites the old inputs.","","## Missing",""]
    lines += ["- "+s for s in r["missing"]]
    lines += ["","## Checks","","| Check | Residual | Pass |","|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--cache-dir",type=Path,required=True)
    args=ap.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"failed":[c for c in r["checks"] if not c["pass"]],
        "IR_coefficients":r["negative_witness_IR_coefficients"],
        "momentum_eigenvalues":[[p["pE2"],p["generalized_eigenvalues"][0]] for p in r["scalar_momentum_probes"]],
        "hard_kinetic_eigenvalues":r["hard_scalar_kinetic"]["generalized_eigenvalues"],
        "uniform_scalar_repair":r["uniform_scalar_parameter_repair"]["diagnostics"],
        "cache":r["cache"]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
