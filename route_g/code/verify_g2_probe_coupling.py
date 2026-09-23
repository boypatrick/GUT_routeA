#!/usr/bin/env python3
"""Candidate circle-uniform Higgs-portal normalization tests, not a fit.

All numerical cards are synthetic and dimensionless in one arbitrary mass
unit. The candidate portal is NOT activated in the existing G2 kernels.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"output"/"g2_probe_coupling"
CARD=dict(R=1.0,alpha=0.25,bare_M_squared=0.25,portal_lambda=0.02,
          v=2.0,mh=4.0,mN=0.8,fN=0.3,eft_cutoff=8.0,
          synthetic_tower_width_allowance=1e-4)


def mass_squared(n,physical_M_squared,alpha=CARD["alpha"],R=CARD["R"]):
    return physical_M_squared+(np.asarray(n)+alpha)**2/R**2


def physical_base(bare_mass_squared,portal_lambda,v):
    return bare_mass_squared+portal_lambda*v*v/2


def width_one(mass,portal_lambda=CARD["portal_lambda"],v=CARD["v"],mh=CARD["mh"]):
    if mh<=2*mass:return 0.0
    beta=np.sqrt((1-2*mass/mh)*(1+2*mass/mh))
    return float(portal_lambda**2*v*v*beta/(16*np.pi*mh))


def width_tower(physical_M_squared,alpha=CARD["alpha"],R=CARD["R"],
                portal_lambda=CARD["portal_lambda"],v=CARD["v"],mh=CARD["mh"],
                eft_cutoff=CARD["eft_cutoff"]):
    # The loose interval follows m_n >= |n+alpha|/R. Exact thresholds filter it.
    low,high=int(np.ceil(-R*mh/2-alpha)),int(np.floor(R*mh/2-alpha))
    rows=[]
    for n in range(low,high+1):
        mass=float(np.sqrt(mass_squared(n,physical_M_squared,alpha,R)))
        if 2*mass>=mh:continue
        beta=float(np.sqrt((1-2*mass/mh)*(1+2*mass/mh)))
        rows.append(dict(n=n,mass=mass,beta=beta,width=width_one(mass,portal_lambda,v,mh),
                         below_eft_cutoff=bool(mass<eft_cutoff)))
    return dict(rigorous_candidate_interval=[low,high],open_modes=rows,
                beta_sum=sum(row["beta"] for row in rows),
                total_width=sum(row["width"] for row in rows),
                all_open_modes_below_eft_cutoff=all(row["below_eft_cutoff"] for row in rows),
                sufficient_eft_energy_condition=bool(eft_cutoff>mh),
                caveat="All kinematically open signed modes are counted; this prediction is qualified only when the external Higgs energy and relevant states lie within a independently controlled EFT")


def kinematics_lab(mS,mN,speed):
    if not 0<=speed<1:raise ValueError("0<=beam speed<1 required")
    root=np.sqrt(1-speed*speed)
    # Rationalized E-m avoids cancellation in the nonrelativistic limit.
    energy_minus_mass=mS*speed*speed/(root*(1+root))
    gap=2*mN*energy_minus_mass
    s=(mS+mN)**2+gap
    kallen=gap*(gap+4*mS*mN)
    p2=kallen/(4*s)
    return dict(s=float(s),kallen=float(kallen),com_momentum_squared=float(p2),
                lab_momentum=float(mS*speed/root),lab_energy=float(mS+energy_minus_mass))


def amplitude_squared(t,portal_lambda,mN,fN,mh):
    return portal_lambda**2*fN**2*mN**2*(4*mN*mN-t)/(mh*mh-t)**2


def nr_cross_section(mS,portal_lambda=CARD["portal_lambda"],mN=CARD["mN"],
                     fN=CARD["fN"],mh=CARD["mh"]):
    return float(portal_lambda**2*fN**2*mN**4/(4*np.pi*mh**4*(mS+mN)**2))


def log_remainder(x):
    """Stable log(1+x)-x/(1+x), including its O(x^2) low-speed limit."""
    if abs(x)<1e-4:
        return sum((-1)**k*(k-1)*x**k/k for k in range(2,14))
    return float(np.log1p(x)-x/(1+x))


def total_cross_section(mS,speed,portal_lambda=CARD["portal_lambda"],mN=CARD["mN"],
                        fN=CARD["fN"],mh=CARD["mh"]):
    kin=kinematics_lab(mS,mN,speed)
    if speed==0:
        value=nr_cross_section(mS,portal_lambda,mN,fN,mh)
        return dict(**kin,speed=speed,total_cross_section=value,
                    direct_t_integral=value,analytic_t_integral=value,nr_limit=value)
    p2,s=kin["com_momentum_squared"],kin["s"]
    # dt=4p^2 du cancels Kallen(s)=4 s p^2, making slow-beam integration stable.
    value=quad(lambda u:amplitude_squared(-4*p2*u,portal_lambda,mN,fN,mh)/(16*np.pi*s),
               0,1,epsabs=1e-25,epsrel=1e-12)[0]
    direct=quad(lambda t:amplitude_squared(t,portal_lambda,mN,fN,mh)/(16*np.pi*kin["kallen"]),
                -4*p2,0,epsabs=1e-25,epsrel=1e-12)[0]
    length=4*p2;a=mh*mh;b=4*mN*mN;x=length/a
    integral=b*length/(a*(a+length))+log_remainder(x)
    analytic=portal_lambda**2*fN**2*mN*mN*integral/(16*np.pi*kin["kallen"])
    return dict(**kin,speed=speed,total_cross_section=float(value),direct_t_integral=float(direct),
                analytic_t_integral=float(analytic),nr_limit=nr_cross_section(mS,portal_lambda,mN,fN,mh))


def gamma_matrices():
    identity=np.eye(2);zero=np.zeros((2,2))
    paulis=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
    return [np.block([[identity,zero],[zero,-identity]])]+[np.block([[zero,s],[-s,zero]]) for s in paulis]


def width_compatibility(beta_sum,mS,width_allowance,card=CARD):
    if beta_sum<=0:
        return dict(beta_sum=beta_sum,bound_available=False,
                    reason="No kinematically open tower mode: no constraint follows from this Higgs partial-width channel")
    coupling_squared_max=16*np.pi*card["mh"]*width_allowance/(card["v"]**2*beta_sum)
    sigma_max=4*card["fN"]**2*card["mN"]**4*width_allowance/(card["v"]**2*card["mh"]**3
                  *(mS+card["mN"])**2*beta_sum)
    return dict(beta_sum=beta_sum,bound_available=True,synthetic_width_allowance=width_allowance,
                portal_lambda_squared_max=float(coupling_squared_max),nr_sigma_max=float(sigma_max),
                assumptions="Fixed physical spectrum; the allowance bounds this entire tower partial width. An invisible-width interpretation additionally requires all counted final states to be experimentally invisible")


def run():
    unchanged_paths=[ROOT/"code"/name for name in ("verify_g2_local_conversion.py","verify_g2_resolution.py","verify_g2_driven.py")]
    before={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in unchanged_paths}
    checks=[]

    def check(name,passed,measured=None,tolerance=None):
        row=dict(name=name,passed=bool(passed))
        if measured is not None:row["measured"]=float(measured)
        if tolerance is not None:row["tolerance"]=float(tolerance)
        checks.append(row)

    labels=np.arange(-3,4);R=CARD["R"];theta0=.7
    overlap=np.empty((len(labels),len(labels)),dtype=complex)
    for i,m in enumerate(labels):
        for j,n in enumerate(labels):
            value=quad(lambda theta:np.cos((n-m)*theta)/(2*np.pi),0,2*np.pi,epsabs=1e-13)[0]
            imaginary=quad(lambda theta:np.sin((n-m)*theta)/(2*np.pi),0,2*np.pi,epsabs=1e-13)[0]
            overlap[i,j]=complex(value,imaginary)
    normalized_local=np.exp(1j*(labels[None,:]-labels[:,None])*theta0)
    check("uniform_circle_overlap:identity",np.max(abs(overlap-np.eye(len(labels))))<1e-13)
    singular=np.linalg.svd(normalized_local,compute_uv=False)
    check("normalized_brane_overlap:rank_one",np.sum(singular>1e-12)==1)
    check("normalized_brane_overlap:nonzero_off_diagonal",abs(normalized_local[0,1])>.99)
    check("uniform_probe:no_mode_conversion",np.max(abs(overlap-np.diag(np.diag(overlap))))<1e-13)
    shifted=physical_base(CARD["bare_M_squared"],CARD["portal_lambda"],CARD["v"])
    original=mass_squared(labels,CARD["bare_M_squared"])
    new=mass_squared(labels,shifted)
    check("common_mass_shift",np.max(abs(new-original-CARD["portal_lambda"]*CARD["v"]**2/2))<1e-14)
    check("shift_preserves_second_difference",np.max(abs(np.diff(new,2)-2/R**2))<1e-14)
    tower=width_tower(shifted)
    open_labels=[row["n"] for row in tower["open_modes"]]
    exhaustive=[n for n in range(-100,101) if 2*np.sqrt(mass_squared(n,shifted))<CARD["mh"]]
    check("width:complete_signed_mode_enumeration",open_labels==exhaustive)
    check("width:sum_once_per_complex_field",abs(tower["total_width"]-CARD["portal_lambda"]**2*CARD["v"]**2*tower["beta_sum"]/(16*np.pi*CARD["mh"]))<1e-18)
    for row in tower["open_modes"]:
        energy=CARD["mh"]/2
        momentum=np.sqrt((energy-row["mass"])*(energy+row["mass"]))
        delta_jacobian=2*momentum/energy
        radial=momentum**2/(energy**2*delta_jacobian)
        phase_space=quad(lambda cosine:2*np.pi*radial/(16*np.pi**2),-1,1)[0]
        independent_width=(CARD["portal_lambda"]*CARD["v"])**2*phase_space/(2*CARD["mh"])
        row["independent_two_body_phase_space"]=float(phase_space)
        row["independent_partial_width"]=float(independent_width)
        check(f"width:n={row['n']}:independent_distinct_pair_phase_space",
              abs(independent_width/row["width"]-1)<1e-13)
    light_neutral=width_tower(shifted,alpha=0)
    neutral_zero=next(row["width"] for row in light_neutral["open_modes"] if row["n"]==0)
    neutral_plus=next(row["width"] for row in light_neutral["open_modes"] if row["n"]==1)
    check("width:neutral_signed_labels_are_distinct_channels",
          [row["n"] for row in light_neutral["open_modes"]]==[-1,0,1]
          and abs(light_neutral["total_width"]-neutral_zero-2*neutral_plus)<1e-18)
    check("width:synthetic_EFT_energy_condition",tower["sufficient_eft_energy_condition"] and tower["all_open_modes_below_eft_cutoff"])
    check("width:threshold_closed",width_one(CARD["mh"]/2)==0)
    check("width:zero_coupling",width_tower(shifted,portal_lambda=0)["total_width"]==0)
    doubled=width_tower(shifted,portal_lambda=2*CARD["portal_lambda"])
    check("width:lambda_squared_scaling_fixed_physical_masses",abs(doubled["total_width"]/tower["total_width"]-4)<1e-13)
    check("width:low_cutoff_not_qualified",not width_tower(shifted,eft_cutoff=.8)["sufficient_eft_energy_condition"])
    mS=float(np.sqrt(mass_squared(0,shifted)))
    trace_rows=[];gammas=gamma_matrices()
    slash=lambda p:p[0]*gammas[0]-sum(p[i]*gammas[i] for i in (1,2,3))
    speed_rows=[]
    for speed in (.3,.1,.03,.01,.003,.001):
        row=total_cross_section(mS,speed)
        row["relative_departure_from_nr"]=abs(row["total_cross_section"]/row["nr_limit"]-1)
        speed_rows.append(row)
        check(f"sigma:v={speed}:independent_t_integral",abs(row["direct_t_integral"]/row["total_cross_section"]-1)<1e-12)
        check(f"sigma:v={speed}:analytic_integral",abs(row["analytic_t_integral"]/row["total_cross_section"]-1)<1e-11)
        momentum=np.sqrt(row["com_momentum_squared"])
        check(f"kinematics:v={speed}:independent_COM_momentum",
              abs(momentum-CARD["mN"]*row["lab_momentum"]/np.sqrt(row["s"]))<1e-14)
        en=np.hypot(CARD["mN"],momentum)
        incoming=np.array([en,0,0,-momentum])
        for cosine in (-1.,-.3,.7,1.):
            outgoing=np.array([en,-momentum*np.sqrt(1-cosine*cosine),0,-momentum*cosine])
            t=-2*momentum**2*(1-cosine)
            trace=np.trace((slash(outgoing)+CARD["mN"]*np.eye(4))@(slash(incoming)+CARD["mN"]*np.eye(4)))/2
            expected=4*CARD["mN"]**2-t
            check(f"Dirac_trace:v={speed}:cos={cosine}",abs(trace-expected)<1e-13)
            trace_rows.append(dict(speed=speed,cosine=cosine,t=t,spin_averaged_trace=float(trace.real),expected=expected))
    check("sigma:NR_convergence",speed_rows[-1]["relative_departure_from_nr"]<1e-6,
          speed_rows[-1]["relative_departure_from_nr"],1e-6)
    check("sigma:monotonic_low_speed_closure",all(a["relative_departure_from_nr"]>b["relative_departure_from_nr"] for a,b in zip(speed_rows,speed_rows[1:])))
    zero=total_cross_section(mS,.1,portal_lambda=0)
    check("sigma:zero_coupling",zero["total_cross_section"]==0)
    base_sigma=total_cross_section(mS,.1)["total_cross_section"]
    check("sigma:lambda_squared_scaling",abs(total_cross_section(mS,.1,portal_lambda=2*CARD["portal_lambda"])["total_cross_section"]/base_sigma-4)<1e-13)
    mediator_scan=[dict(mh=mh,sigma=total_cross_section(mS,.1,mh=mh)["total_cross_section"]) for mh in (2.,4.,8.)]
    check("sigma:heavier_mediator_suppression",all(a["sigma"]>b["sigma"] for a,b in zip(mediator_scan,mediator_scan[1:])))
    neutral_base=physical_base(25.,.03,CARD["v"])
    degeneracy=[]
    for j in (1,2,3):
        plus=float(np.sqrt(mass_squared(j,neutral_base,alpha=0)))
        minus=float(np.sqrt(mass_squared(-j,neutral_base,alpha=0)))
        sp=total_cross_section(plus,.1,portal_lambda=.03)["total_cross_section"]
        sm=total_cross_section(minus,.1,portal_lambda=.03)["total_cross_section"]
        check(f"uniform_neutral_portal:j={j}:sign_blind",plus==minus and sp==sm)
        degeneracy.append(dict(abs_j=j,mass_plus=plus,mass_minus=minus,sigma_plus=sp,sigma_minus=sm))
    heavy_tower=width_tower(neutral_base,alpha=0,portal_lambda=.03)
    no_bound=width_compatibility(heavy_tower["beta_sum"],degeneracy[0]["mass_plus"],CARD["synthetic_tower_width_allowance"])
    check("closed_tower:no_width_constraint",heavy_tower["beta_sum"]==0 and not no_bound["bound_available"])
    bound=width_compatibility(tower["beta_sum"],mS,CARD["synthetic_tower_width_allowance"])
    allowed_lambda=np.sqrt(bound["portal_lambda_squared_max"])
    saturated=width_tower(shifted,portal_lambda=allowed_lambda)["total_width"]
    saturated_sigma=nr_cross_section(mS,portal_lambda=allowed_lambda)
    check("compatibility:width_saturation",abs(saturated/CARD["synthetic_tower_width_allowance"]-1)<1e-13)
    check("compatibility:eliminated_lambda_sigma",abs(saturated_sigma/bound["nr_sigma_max"]-1)<1e-13)
    ratio=4*CARD["fN"]**2*CARD["mN"]**4/(CARD["v"]**2*CARD["mh"]**3*(mS+CARD["mN"])**2*tower["beta_sum"])
    check("compatibility:sigma_over_tower_width",abs(nr_cross_section(mS)/tower["total_width"]/ratio-1)<1e-13)
    attenuation=[]
    for number_density,path_length in ((0.,1.),(.4,3.),(1e3,1e3)):
        depth=number_density*path_length*base_sigma
        probability=float(-np.expm1(-depth))
        attenuation.append(dict(synthetic_number_density=number_density,synthetic_path_length=path_length,
                                optical_depth=depth,at_least_one_interaction_probability=probability))
        check(f"attenuation:n={number_density}:d={path_length}",0<=probability<=1 and probability<=depth)
    check("attenuation:dilute_linear_limit",abs(attenuation[1]["at_least_one_interaction_probability"]/attenuation[1]["optical_depth"]-1)<1e-7)
    matching_rows=[]
    for scale in (.1,1.,10.):
        frozen_physical_m2=25*scale**2
        maximum_lambda=2*frozen_physical_m2/CARD["v"]**2
        for factor in (0.,.5,1.,1.1):
            coupling=factor*maximum_lambda
            bare_m2=frozen_physical_m2-coupling*CARD["v"]**2/2
            tolerance=1e-12*frozen_physical_m2
            gate=bare_m2>=-tolerance
            expected_sign=1 if factor<1 else 0 if factor==1 else -1
            sign_correct=abs(bare_m2)<tolerance if expected_sign==0 else expected_sign*bare_m2>0
            check(f"X_only_fixed_mass:Lambda={scale}:factor={factor}:mass_matching",
                  abs(physical_base(bare_m2,coupling,CARD["v"])-frozen_physical_m2)<tolerance)
            check(f"X_only_fixed_mass:Lambda={scale}:factor={factor}:bare_sign",sign_correct)
            check(f"X_only_fixed_mass:Lambda={scale}:factor={factor}:necessary_gate",gate==(factor<=1))
            matching_rows.append(dict(synthetic_scale_Lambda=scale,frozen_physical_MD_squared=frozen_physical_m2,
                                      v=CARD["v"],lambda_X_max=maximum_lambda,factor_of_maximum=factor,
                                      lambda_X=coupling,bare_MD_squared=bare_m2,
                                      necessary_no_runaway_gate=bool(gate),
                                      H_zero_X_constant_potential_differences=[bare_m2*amplitude**2 for amplitude in (1.,10.,100.)]))
    after={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in unchanged_paths}
    check("existing_G2_kernels_unchanged",before==after)
    result=dict(status="Synthetic candidate portal normalization and compatibility tests; not activated in G2 or empirically qualified",
        synthetic_card=CARD,physical_mass_base_squared=shifted,probe_scalar_mass=mS,
        action=dict(candidate="-HdaggerH[lambdaPhi sum_n|phi_n|^2+lambdaX sum_l|x_l|^2]",
                    higgs_expansion="HdaggerH=(v+h)^2/2; common Delta M_S^2=lambdaS v^2/2; h Sdagger S vertex=-i lambdaS v",
                    nucleon_matching="L_hNN=-(fN mN/v) h Nbar N, treated as an explicit postulated low-energy matching coefficient",
                    activation="Candidate only. Existing G2 and powered kernels remain unchanged.",
                    both_portals="If both lambdaPhi and lambdaX are nonzero, tree h-mediated Phi-X diagonal elastic scattering is added and interferes with existing q=0 elastic scattering. Old total branch fractions cannot be retained.",
                    x_only_option="lambdaPhi=0 removes that additional tree h-exchange Phi-X amplitude, but does not prove absence of radiative corrections or a complete apparatus"),
        overlap=dict(labels=labels.tolist(),uniform_real=overlap.real.tolist(),uniform_imag=overlap.imag.tolist(),
                     normalized_brane_real=normalized_local.real.tolist(),normalized_brane_imag=normalized_local.imag.tolist(),
                     brane_singular_values=singular.tolist(),
                     normalization="Uniform integral int R dtheta u_m* u_n=delta_mn; normalized brane matrix=(2 pi R)u_m*(theta0)u_n(theta0). The rescaling removes dimensions for rank comparison, not a physical equality of couplings."),
        tower_width=tower,neutral_width_counting_regression=light_neutral,
        neutral_sign_degeneracy=degeneracy,closed_heavy_tower=heavy_tower,
        nucleon_scattering=dict(amplitude_squared="lambdaS^2 fN^2 mN^2 (4mN^2-t)/(mh^2-t)^2",
                    differential="d sigma/dt=|M|^2/[16 pi Kallen(s,mS^2,mN^2)], -4pCOM^2<=t<=0",
                    nr_limit="lambdaS^2 fN^2 mN^4/[4 pi mh^4(mS+mN)^2]",nr_value=nr_cross_section(mS),
                    speed_rows=speed_rows,dirac_trace_checks=trace_rows,mediator_scan=mediator_scan,
                    matching_scope="Pointlike effective scalar nucleon coupling with fixed synthetic fN; form-factor variation, nuclear coherence, target composition and empirical matching are absent"),
        coupling_elimination=dict(open_tower=bound,closed_tower=no_bound,saturated_width=saturated,saturated_nr_sigma=saturated_sigma,
                    formula="lambdaS^2_max=16 pi mh Gamma_allow/(v^2 Bsum); sigma_SN_max=4 fN^2 mN^4 Gamma_allow/[v^2 mh^3(mS+mN)^2 Bsum]",
                    mass_convention="These lambda^2 scalings and elimination hold a physical spectrum fixed. If bare masses are fixed instead, lambda changes masses, beta_sum and open thresholds, so the constraint is implicit; no tuning/naturalness claim is made."),
        per_pass=attenuation,
        x_only_fixed_mass_matching=dict(formula="MD_bare^2=25 Lambda^2-lambdaX v^2/2; lambdaX<=50 Lambda^2/v^2",
                rows=matching_rows,
                assumptions="No X self-quartic or other stabilizing X interactions; X is circle-neutral and admits a constant j=0 mode. Set H=Phi=0.",
                proof="Along constant X, gradients and all portal terms vanish. V(X)-V(0)=MD_bare^2 |X|^2; negative MD_bare^2 is unbounded below. Equality is a flat direction of this restricted configuration, not strict isolation.",
                scope="A necessary tree-level consistency gate for matching a fixed physical mass, not a sufficient vacuum/perturbativity/EFT certificate. Large synthetic lambda values are not recommendations.",
                charged_phi_qualification="Do not transfer this argument to Phi: nontrivial holonomy can enforce a positive minimum covariant-gradient contribution, so a negative Phi bare mass coefficient alone need not imply the same constant-mode runaway."),
        limits=[
            "All cards are synthetic in one arbitrary mass unit: no GeV map, measured Higgs/nucleon inputs, target performance, or empirical exclusion is supplied.",
            "Each open signed Fourier label represents one complex scalar field; its particle-antiparticle final state has no extra factor two. Distinct plus/minus labels must each be counted when open.",
            "Finite tower enumeration is exact kinematics for the declared spectrum, not a UV completion. The relevant energies and states must lie below an independently controlled EFT cutoff.",
            "Uniform coupling is KK diagonal and cannot directly distinguish neutral +j from -j. Brane-local rank-one off-diagonal couplings are a different, unadopted physical assumption.",
            "A partial-width allowance is not automatically an invisible-width bound: invisibility and applicable experimental selection must be established separately.",
            "Nonnegative positive-lambda mass shifts preserve tower curvature; maintaining the original physical G2 masses requires an explicit matching prescription, not silently adding shifts.",
            "With no X self-quartic, X-only matching at fixed MD^2=25 Lambda^2 requires lambdaX<=50 Lambda^2/v^2 to avoid the H=Phi=0 constant-X runaway. This is necessary, not a full stability or weak-coupling guarantee.",
            "Nonzero portals for both species modify the elastic amplitude and hence total rates and conditional fractions; this test does not activate them in the old kernels.",
            "P=1-exp(-nT d sigma) is only an independent dilute-target first-interaction model. It supplies neither momentum resolution, recoil-sign inference, exposure, timing, nor a real detector efficiency.",
            "Dark-matter abundance, cosmology, nuclear many-body matching, physical momentum calibration, loops and UV completion are not inferred."],
        checks=checks,summary=dict(checks=len(checks),passed=sum(row["passed"] for row in checks),
                                  failed=[row["name"] for row in checks if not row["passed"]]),
        source_sha256={**after,Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    write_report(result)
    print(json.dumps(result["summary"],indent=2))
    print(json.dumps(dict(open_n=open_labels,beta_sum=tower["beta_sum"],width=tower["total_width"],
                         nr_sigma=nr_cross_section(mS),last_NR_error=speed_rows[-1]["relative_departure_from_nr"],
                         bound=bound),indent=2))
    if result["summary"]["failed"]:raise SystemExit(1)


def write_report(result):
    tower=result["tower_width"];scattering=result["nucleon_scattering"];bound=result["coupling_elimination"]["open_tower"]
    lines=["# Route G2 — candidate circle-uniform probe coupling tests","",
        "This is a synthetic EFT normalization test, not a physical-scale calibration, particle fit, or detector qualification. No existing G2 scattering kernel is changed. Every number uses one arbitrary engineering mass unit; none is a measured Higgs, nucleon or target input.","",
        "## Candidate interaction and selection rule","",
        r"$$\mathcal L_{\rm portal}=-H^\dagger H\left[\lambda_\Phi\sum_n|\phi_n|^2+\lambda_X\sum_l|x_l|^2\right],\qquad H^\dagger H=\frac{(v+h)^2}{2}.$$","",
        r"$$\int_0^{2\pi}R\,d\theta\,u_m^*u_n=\delta_{mn},\quad \Delta M_S^2=\lambda_Sv^2/2,\quad hS^*S:\ -i\lambda_Sv.$$","",
        "The normalized uniform-circle overlap is the identity. A different, brane-local probe has matrix u_m*(theta0)u_n(theta0), rank one with off-diagonal entries; this alternative is not adopted. Its matrix is rescaled by 2 pi R only for a dimensionless rank comparison.","",
        "A common mass shift preserves the signed-label second difference 2/R^2. A neutral tower remains exactly degenerate under j -> -j, and this uniform scalar probe has identical masses and cross sections for both signs. It cannot serve as a direct sign tag.","",
        "If both species' portals are nonzero, an extra tree h-mediated diagonal Phi-X elastic amplitude interferes with the original q=0 elastic amplitude. Old total rates and branch fractions must then be recomputed. The lambdaPhi=0, X-only option avoids this additional tree amplitude, but is not a radiative-stability proof.","",
        "## Complex-scalar Higgs widths","",
        r"$$\Gamma(h\to S_nS_n^*)=\frac{\lambda_S^2v^2}{16\pi m_h}\sqrt{1-\frac{4m_n^2}{m_h^2}},\qquad 2m_n<m_h.$$","",
        "Each signed label is one complex field. Its distinguishable particle-antiparticle pair is already included: do not multiply by another two. Different +n and -n fields are separate channels when kinematically open. An independent radial two-body phase-space integration verifies the width normalization, and a neutral light tower checks Gamma_total=Gamma_0+2 Gamma_1. The finite signed-label bound is followed by the exact threshold; no infinite degeneracy factor is inserted.","",
        "Synthetic card: R=1, alpha=0.25, bare M^2=0.25, lambda=0.02, v=2, mh=4, mN=0.8, fN=0.3, stated EFT cutoff=8. After the common shift, physical M^2=0.29.","",
        "| signed n | physical scalar mass | beta | partial width |","|---:|---:|---:|---:|"]
    for row in tower["open_modes"]:
        lines.append(f"| {row['n']} | {row['mass']:.9g} | {row['beta']:.9g} | {row['width']:.9g} |")
    lines += ["",f"Bsum={tower['beta_sum']:.9g}; total width={tower['total_width']:.9g}. The synthetic neutral heavy tower has no open modes at mh=4, hence no bound from this partial-width channel. The listed cutoff is an assumed scope condition, not an independently derived UV cutoff.","",
        "## Explicit low-energy nucleon matching","",
        r"$$\mathcal L_{hNN}=-\frac{f_Nm_N}{v}h\bar NN,\qquad\overline{|\mathcal M|^2}=\frac{\lambda_S^2f_N^2m_N^2(4m_N^2-t)}{(m_h^2-t)^2}.$$","",
        r"$$\frac{d\sigma}{dt}=\frac{\overline{|\mathcal M|^2}}{16\pi\mathcal K(s,m_S^2,m_N^2)},\quad -4p_{\rm COM}^2\leq t\leq0,\quad\sigma_{\rm NR}=\frac{\lambda_S^2f_N^2m_N^4}{4\pi m_h^4(m_S+m_N)^2}.$$","",
        "The initial nucleon spin is averaged and the final spin summed. Explicit Dirac-matrix traces verify 4mN^2-t. Direct t integration, a rescaled integration, and an analytic antiderivative agree. This is a fixed synthetic pointlike fN matching ansatz; momentum-dependent form factors and nuclear composition are absent.","",
        "| scalar speed in nucleon-rest frame | integrated sigma | relative departure from NR |","|---:|---:|---:|"]
    for row in scattering["speed_rows"]:
        lines.append(f"| {row['speed']:g} | {row['total_cross_section']:.10g} | {row['relative_departure_from_nr']:.6g} |")
    lines += ["",f"NR limit={scattering['nr_value']:.10g} in inverse-square engineering mass units. Tests also cover zero coupling, lambda-squared scaling at fixed physical masses, increasing-mediator suppression, and exact neutral +/-j degeneracy.","",
        "## Width versus interaction compatibility","",
        r"For $B=\sum_{\rm open}\sqrt{1-4m_n^2/m_h^2}>0$ and a hypothetical upper allowance $\Gamma_{\rm allow}$ on this entire tower partial width:","",
        r"$$\lambda_S^2\leq\frac{16\pi m_h\Gamma_{\rm allow}}{v^2B},\qquad \sigma_{SN}^{\rm NR}\leq\frac{4f_N^2m_N^4\Gamma_{\rm allow}}{v^2m_h^3(m_S+m_N)^2B}.$$","",
        f"Using only the synthetic allowance Gamma_allow=0.0001 gives lambda^2_max={bound['portal_lambda_squared_max']:.9g} and sigma_max={bound['nr_sigma_max']:.9g}. Substitution saturates both algebraic bounds. If B=0, this decay channel supplies no such limit. These are not experimental bounds.","",
        "The elimination and lambda-squared tests keep the physical spectrum fixed. Holding bare masses fixed instead changes the common portal mass shift, beta sum and thresholds, making the constraint implicit. Interpreting a partial-width allowance as an invisible-width constraint additionally requires the final states really to be invisible under the relevant selection.","",
        "## A necessary X-only fixed-mass matching gate","",
        r"$$M_{D,\rm bare}^2=25\Lambda^2-\frac{\lambda_Xv^2}{2},\qquad \lambda_X\leq\frac{50\Lambda^2}{v^2}.$$","",
        "This condition follows only for the stated action with no X self-quartic or other stabilizing X interaction. At H=Phi=0, the neutral j=0 X mode may be constant, so its gradients and portal terms vanish: V(X)-V(0)=MD_bare^2 |X|^2. A negative coefficient runs to minus infinity as |X| grows. At equality this restricted direction is flat, not strictly isolated.","",
        "Synthetic checks use Lambda=0.1,1,10, v=2 and lambdaX at 0,0.5,1,1.1 times the bound. All preserve the frozen physical MD^2=25 Lambda^2; the last factor gives a negative bare coefficient and fails this necessary gate. Large algebraic coupling values are not weak-coupling or EFT recommendations.","",
        "This constant-X proof must not be generalized to the charged Phi field: nontrivial circle holonomy can supply an unavoidable covariant-gradient term, so a negative Phi bare mass coefficient alone need not imply the same runaway. The X bound itself does not certify the full vacuum, loop stability, perturbativity, or detector performance.","",
        "## A per-pass probability is not a resolution model","",
        r"$$P(\text{at least one interaction})=1-e^{-n_Td\sigma}.$$","",
        "This expression is tested only as an independent dilute-target first-interaction law with synthetic number density and path length. It does not determine momentum resolution, sign-classification error, trigger acceptance, luminosity, or a real detector efficiency. A larger coupling is therefore not itself a readout specification.","",
        f"**{result['summary']['passed']}/{result['summary']['checks']} checks pass.** Existing G2, G2-S and driven verifier hashes are unchanged. No GeV map, physical hardware number, cosmology, dark-matter abundance, experimental exclusion, or UV completion has been added.","",
        "Reproduce: python3 route_g/code/verify_g2_probe_coupling.py. JSON retains overlap matrices, all widths, spin traces, velocity checks, hypothetical compatibility bounds and source hashes."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":run()
