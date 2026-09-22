#!/usr/bin/env python3
"""Bounded physics-first LO scale and gauge-proton screen, not a P54 fit.

Regenerates the four-parent PS beta coefficients from actual Casimirs.
Nine declared whole-parent mass cards are leading-log hypotheses, not
action-derived spectra or phenomenologically accepted benchmarks.
"""
from __future__ import annotations
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

import verify_p54_ps_active_census as census

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_scales_lo_screen"
PARENTS = RF / "output/p54_p2_two_site_matching.json"
MZ = 91.1876
# A frozen, explicitly approximate screening card, not a 2026 correlated fit.
INPUT = dict(alpha_em_inverse=127.955, sin2_theta_w=.23122, alpha_s=.1180)
A_SM = np.array([-7., -19/6, 41/10])
P = np.array([[1.,0.,0.],[0.,1.,0.],[2/5,0.,3/5]])
ACTIVE = ("phi10:(1,2,2)", "Sigma126:(15,2,2)",
          "Sigma126:(10-pair,3,1)", "Sigma126:(10-pair,1,3)")
HI20 = "Phi54:(20prime,1,1)"
HI33 = "Phi54:(1,3,3)"
MID15 = "Sigma126:(15,2,2)"
YEAR_S = 365.25*24*3600
HBAR_GEV_S = 6.582119569e-25
MP = .93827208816
MPI = .1349768
TAU_LIMIT = 2.4e34
# Yoo et al. 2111.01608 Table 8, continuum, MSbar(2 GeV), Q^2=0.
W_PLUS = -.159
W_PLUS_ERRORS = np.array([.015,.020,.025])
W_PI = abs(W_PLUS)/np.sqrt(2)
W_ERR = np.linalg.norm(W_PLUS_ERRORS)/np.sqrt(2)


def parent_beta(row):
    """Real canonical trace convention: delta b=T_real/6."""
    return row["dimension"]*np.array([row["C4"],row["C2L"],row["C2R"]])/np.array([15,3,3])/6


def alpha_sm0():
    a,s = INPUT["alpha_em_inverse"],INPUT["sin2_theta_w"]
    return np.array([1/INPUT["alpha_s"],a*s,3*a*(1-s)/5])


def segment_rg(a0, b, length, gamma):
    """Gauge-only dlog C downward; exact one-loop segment integral."""
    a1 = a0-b*length/(2*np.pi)
    if min(a0.min(),a1.min()) <= 0:
        raise ValueError("Gauge Landau pole in a declared interval")
    result = np.zeros(3)
    for k in range(3):
        result[k] = gamma[k]/b[k]*np.log(a0[k]/a1[k]) if abs(b[k])>1e-13 else gamma[k]*length/(2*np.pi*a0[k])
    return a1,float(np.sum(result))


def low_qcd_factor():
    inv = 1/INPUT["alpha_s"]
    logfactor = 0.
    for hi,lo,nf in ((MZ,4.18,5),(4.18,2.,4)):
        b0 = 11-2*nf/3
        new = inv-b0*np.log(hi/lo)/(2*np.pi)
        logfactor += 2/b0*np.log(inv/new)
        inv = new
    return float(np.exp(logfactor)),float(1/inv)


def solve_card(card_id, factors, rows, aps):
    # Lower and upper logarithms are tied to actual whole-parent masses.
    shift_ps = sum((parent_beta(rows[label])*np.log(v[1]) for label,v in factors.items()),start=np.zeros(3))
    delta = -P@shift_ps/(2*np.pi)
    design = np.column_stack([np.ones(3),A_SM/(2*np.pi),P@aps/(2*np.pi)])
    alpha_u,t,u = np.linalg.solve(design,alpha_sm0()-delta)
    mi,mu = MZ*np.exp(t),MZ*np.exp(t+u)
    inv_sm_i = alpha_sm0()-A_SM*t/(2*np.pi)
    inv_ps_i = np.linalg.solve(P,inv_sm_i)
    g4 = np.sqrt(4*np.pi/inv_ps_i[0]);gu = np.sqrt(4*np.pi/alpha_u)
    sigma,omega = mi/g4,mu/gu
    rho = sigma/omega
    mx,mxp = mu,mu*np.sqrt(1+rho*rho)

    # Piecewise PS running, with no group-by-group adjustable threshold.
    masses = {k:(mi if site=="MI" else mu)*fac for k,(site,fac) in factors.items()}
    cuts = sorted(set([mi,mu]+[m for m in masses.values() if mi<m<mu]))
    ai = inv_ps_i.copy();log_ps=0.;segments=[]
    for left,right in zip(cuts[:-1],cuts[1:]):
        middle=np.sqrt(left*right);beta=aps.copy()
        for label,(site,fac) in factors.items():
            if site=="MI" and middle<masses[label]: beta-=parent_beta(rows[label])
            if site=="MU" and middle>masses[label]: beta+=parent_beta(rows[label])
        nxt,lg = segment_rg(ai,beta,np.log(right/left),np.array([15/4,9/4,9/4]))
        segments.append(dict(low_GeV=float(left),high_GeV=float(right),b_PS=beta.tolist(),alpha_inverse_low=ai.tolist(),alpha_inverse_high=nxt.tolist()))
        ai=nxt;log_ps+=lg
    above = sum((parent_beta(rows[k])*np.log(m/mu) for k,m in masses.items() if m>mu),start=np.zeros(3))
    residual = ai+above/(2*np.pi)-alpha_u
    _,log_sm_l=segment_rg(alpha_sm0(),A_SM,t,np.array([2,9/4,23/20]))
    _,log_sm_r=segment_rg(alpha_sm0(),A_SM,t,np.array([2,9/4,11/20]))
    alow,alpha2=low_qcd_factor()
    al,ar=alow*np.exp(log_ps+log_sm_l),alow*np.exp(log_ps+log_sm_r)

    # Define scalar-three-quark Wilson convention explicitly. F includes
    # Fierz/Clebsch and flavor contractions, not just an arbitrary unit number.
    k2=gu*gu/(2*mx*mx)
    phase=MP/(32*np.pi)*(1-MPI**2/MP**2)**2
    width_unit=phase*W_PI**2*k2**2
    tau_unit=HBAR_GEV_S/(width_unit*YEAR_S)
    fmax=np.sqrt(tau_unit/TAU_LIMIT)
    out=dict(id=card_id,declared_parent_mass_ratios={k:dict(site=v[0],ratio=v[1]) for k,v in factors.items()},
             parent_masses_GeV=masses,log10_MI=float(np.log10(mi)),log10_MU=float(np.log10(mu)),
             log10_MX=float(np.log10(mx)),MI_GeV=float(mi),MU_GeV=float(mu),MX_GeV=float(mx),MXprime_GeV=float(mxp),
             alphaU_inverse=float(alpha_u),gU=float(gu),alpha_PS_inverse_at_MI=inv_ps_i.tolist(),g4I=float(g4),
             sigma_GeV=float(sigma),omega_GeV=float(omega),sigma_over_omega=float(rho),
             delta_alpha_inverse_from_declared_logs=delta.tolist(),PS_segments=segments,
             PS_upper_matching_residual=residual.tolist(),gauge_only_A_L=float(al),gauge_only_A_R=float(ar),
             low_QCD_A=float(alow),alpha_s_2GeV_LO=alpha2,MX2_over_MXprime2=float(mx*mx/(mxp*mxp)),
             gauge_operator_k2_GeVminus2=float(k2),tau_e_pi_unit_low_energy_F_years=float(tau_unit),
             allowed_effective_flavor_norm=float(fmax),allowed_UV_F_L_if_F_R_zero=float(fmax/al),
             allowed_UV_F_R_if_F_L_zero=float(fmax/ar),
             allowed_effective_flavor_norm_hadronic_1sigma=[float(fmax*W_PI/(W_PI+W_ERR)),float(fmax*W_PI/(W_PI-W_ERR))],
             hierarchy_ok=bool(t>0 and u>0),perturbative_gauge_ok=bool(alpha_u>1 and min(inv_ps_i)>1),
             action_realization_tested=False,flavor_rotations_fitted=False,scalar_mediated_amplitude_included=False,
             physical_candidate=False,
             conditions=["A_L^2 |F_L|^2 + A_R^2 |F_R|^2 < allowed_effective_flavor_norm^2 in the stated gauge-only approximation",
                         "F_L,R must be calculated from the same flavor solution and both vector masses; not freely chosen nuisance counterterms",
                         "Whole-parent mass card must be realized by a stable same-action scalar point before promotion",
                         "Other proton channels and scalar-mediated amplitudes remain untested"])
    return out


def actual_vector_audit(rho):
    """Actual canonical orbit from Phi and normalized decomposable Sigma.

    For Omega=u1 wedge ... wedge u5 and P=U Udag, the induced wedge
    inner product is Tr[P Ta^dag(1-P)Tb]+Tr[P Ta^dag]Tr[P Tb].
    This avoids JAX and a 328-real Hessian entirely, while preserving the
    original action normalization and five complex planes.
    """
    basis=[]
    for a in range(10):
        for b in range(a+1,10):
            q=np.zeros((10,10));q[a,b]=1/np.sqrt(2);q[b,a]=-1/np.sqrt(2);basis.append(q)
    phi=np.diag([-2/5]*6+[3/5]*4)
    u=np.column_stack([(np.eye(10)[2*k]+1j*np.eye(10)[2*k+1])/np.sqrt(2) for k in range(5)])
    proj=u@u.conj().T;comp=np.eye(10)-proj
    dphi=[t@phi-phi@t for t in basis]
    mass=np.array([[np.trace(dp.T@dq)+2*rho*rho*np.real(np.trace(proj@ta.T@comp@tb)+np.trace(proj@ta.T)*np.trace(proj@tb))
                    for dq,tb in zip(dphi,basis)] for dp,ta in zip(dphi,basis)])
    yg10=np.zeros((10,10))
    for k,y in enumerate([-1/3]*3+[1/2]*2):yg10[2*k,2*k+1]=y;yg10[2*k+1,2*k]=-y
    yg=np.array([[np.trace(a.T@(yg10@b-b@yg10)) for b in basis] for a in basis]);y2=-yg@yg
    eigen,v=np.linalg.eigh(mass);q=np.diag(v.T@y2@v)
    selections={"X_Y5over6":np.abs(q-25/36)<1e-8,"Xprime_Y1over6":np.abs(q-1/36)<1e-8}
    return dict(rho=rho,masses_over_gUomega={k:np.sqrt(eigen[sel]).tolist() for k,sel in selections.items()},
                multiplicities={k:int(sum(sel)) for k,sel in selections.items()},
                full_mass2_eigenvalues_over_gU2omega2=eigen.tolist(),
                intermediate_charged_vector_mass2_over_g2sigma2=float(np.mean(eigen[12:20])/rho**2))


def run():
    data=json.loads(PARENTS.read_text());rows={r["label"]:r for r in data["parent_census"]}
    active=[]
    for label in ACTIVE:
        row=rows[label];c=[census.rational(row[k]) for k in ("C4","C2L","C2R")]
        dim=census.Q(row["dimension"],2)
        active.append(dict(C=c,T_complex=[dim*z/d for z,d in zip(c,(15,3,3))]))
    af,_=census.beta_coefficients(active);aps=np.array([float(z) for z in af])
    baseline=solve_card("four_parent_degenerate_LO",{},rows,aps)
    cards=[]
    for k20,k15 in itertools.product((1/3,1.,3.),(1.,2.,3.)):
        factors={HI20:("MU",k20),HI33:("MU",1/k20),MID15:("MI",k15)}
        cards.append(solve_card(f"P20_{k20:g}_P33_{1/k20:g}_B15_{k15:g}",factors,rows,aps))
    vector=actual_vector_audit(baseline["sigma_over_omega"])
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        a,b=np.asarray(a),np.asarray(b);err=float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append(dict(name=name,residual=err,tolerance=tol,passed=err<tol))
    check("correct_four_parent_PS_beta",aps,[2/3,26/3,26/3])
    check("real_Phi20_parent_beta",parent_beta(rows[HI20]),[4/3,0,0])
    check("real_Phi33_parent_beta",parent_beta(rows[HI33]),[0,1,1])
    check("complex_Sigma15_parent_beta",parent_beta(rows[MID15]),[16/3,5,5])
    check("baseline_independent_LO_reference",[baseline["MI_GeV"],baseline["MU_GeV"]],[5.082267472797e13,1.899353581072e15])
    invariant=np.array([-2.,-3.,5.])
    t_parity=2*np.pi*(invariant@alpha_sm0())/(invariant@A_SM)
    check("independent_D_parity_MI_formula",np.log(baseline["MI_GeV"]/MZ),t_parity)
    check("D_parity_SM_beta_combination_is_44",invariant@A_SM,44.)
    check("actual_X_multiplicity",vector["multiplicities"]["X_Y5over6"],12)
    check("actual_Xprime_multiplicity",vector["multiplicities"]["Xprime_Y1over6"],12)
    check("actual_X_mass_normalization",vector["masses_over_gUomega"]["X_Y5over6"],np.ones(12))
    check("actual_Xprime_mass_normalization",vector["masses_over_gUomega"]["Xprime_Y1over6"],np.full(12,np.sqrt(1+baseline["sigma_over_omega"]**2)))
    check("actual_intermediate_charged_vector_normalization",vector["intermediate_charged_vector_mass2_over_g2sigma2"],1.)
    for card in cards:
        check(card["id"]+"_piecewise_threshold_reconstruction",card["PS_upper_matching_residual"])
        check(card["id"]+"_D_parity_fixes_MI",card["MI_GeV"],baseline["MI_GeV"])
        check(card["id"]+"_sigma_uses_g4_at_MI",card["g4I"]*card["sigma_GeV"],card["MI_GeV"])
        check(card["id"]+"_vector_mass_uses_unification_coupling",card["gU"]*card["omega_GeV"],card["MX_GeV"])
    phase=MP/(32*np.pi)*(1-MPI**2/MP**2)**2
    cmax=np.sqrt(HBAR_GEV_S/(YEAR_S*TAU_LIMIT*phase*W_PI**2))
    check("coefficient_limit_reproduces_unit_amplitude_lifetime",cmax,baseline["gauge_operator_k2_GeVminus2"]*baseline["allowed_effective_flavor_norm"])
    # A scalar ODE quadrature provides an independent sign/normalization
    # regression for the exponent, including both signs and zero beta.
    from scipy.integrate import quad
    for beta in (-7.,0.,26/3):
        a0=np.array([38.,42.,45.]);bb=np.array([beta,beta,beta]);gamma=np.array([2.,9/4,15/4]);length=2.3
        _,analytic=segment_rg(a0,bb,length,gamma)
        direct=quad(lambda x:np.sum(gamma/(2*np.pi*(a0-bb*x/(2*np.pi)))),0.,length,epsabs=1e-12)[0]
        check("RG_exponent_independent_integral_beta_"+str(beta),analytic,direct)
    paths=[Path(__file__),Path(census.__file__),PARENTS,RF/"code/verify_p54_p1_hessian_spectrum.py",RF/"code/verify_p54_p2_two_site_matching.py"]
    return dict(schema="p54-scales-lo-physics-screen-v1",date="2026-09-22",input_card=INPUT,MZ_GeV=MZ,
                ordering=dict(SM=["3","2","1_GUT"],PS=["4","L","R"]),a_SM=A_SM.tolist(),a_PS=aps.tolist(),a_PS_exact=[str(z) for z in af],
                active_parent_labels=list(ACTIVE),baseline=baseline,scenarios=cards,actual_vector_audit=vector,
                normalization=dict(MI="charged PS leptoquark mass = g4(MI)*sigma",MU="X(Y=5/6) mass = gU*omega",
                    Xprime="gU*sqrt(omega^2+sigma^2)",proton_coefficient="C_L,R(2GeV)=[gU^2/(2 MX^2)] A_L,R F_L,R; F includes Fierz/Clebsch/flavor and Xprime ratio"),
                proton_input=dict(mp_GeV=MP,mpi0_GeV=MPI,W_pi0_LR_GeV2=W_PI,W_pi0_quadrature_error_GeV2=W_ERR,
                    W_plus_LR_GeV2=W_PLUS,W_plus_errors=W_PLUS_ERRORS.tolist(),hadronic_scheme="MSbar 2 GeV",tau_epi_limit_years=TAU_LIMIT,confidence_level=.90,
                    scalar_operator_coefficient_norm_bound_GeVminus2=float(cmax)),
                precision_contract=dict(one_loop_gauge_running=True,tree_gauge_matching=True,whole_parent_leading_logs=True,
                    finite_gauge_matching_constants=False,two_loop_running=False,finite_yukawa_matching=False,
                    dimension_six_RG="one-loop gauge-only diagonal factors, PS/SM plus piecewise low-QCD to 2 GeV; no Yukawa anomalous mixing",
                    hadronic_errors="quadrature sensitivity only, no joint experimental/lattice confidence region",
                    no_universal_proton_exclusion=True,no_accepted_action_benchmark=True),
                falsification_conditions=["Reject a declared scale card if ln(MI/MZ)<=0, ln(MU/MI)<=0 or a gauge Landau pole occurs below MU",
                    "Reject a specified gauge-only flavor realization if A_L^2|F_L|^2+A_R^2|F_R|^2 exceeds the saved limit squared",
                    "A full-model exclusion needs every allowed flavor/threshold realization to fail, including other-channel and scalar-amplitude consistency",
                    "A scale solution alone never certifies scalar vacuum, light-doublet or flavor viability"],
                references=[dict(topic="gauge operator normalization and gauge-only anomalous exponents",url="https://arxiv.org/html/1507.06712v2",location="Sections 5 and 10; old beta table not imported"),
                    dict(topic="physical-pion lattice form factor",url="https://arxiv.org/html/2111.01608v1",location="Table 8 and isospin relation; continuum MSbar(2 GeV)"),
                    dict(topic="p to e+ pi0 experimental limit",url="https://arxiv.org/abs/2010.16098",location="published 90 percent CL partial lifetime, not a newly measured 2026 value")],
                checks=checks,passed=sum(x["passed"] for x in checks),total=len(checks),
                source_hashes={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def markdown(r):
    b=r["baseline"]
    lines=["# P54 bounded LO scale / gauge-proton screen", "", "**Conditional scale solutions, not accepted physical benchmarks.** The historical frozen scalar point and its old P2 scales are not inputs.","",
        f"Checks: **{r['passed']}/{r['total']}**. Correct four-parent PS beta: `(2/3,26/3,26/3)`, not historical `(1,26/3,26/3)`.","",
        "## Declared approximation", "", "One-loop gauge running, tree gauge matching, one light SM Higgs doublet, three 16F families and the four saved active PS parents. The gauge-singlet axion does not affect these gauge betas. The central input card is a deliberately approximate inherited screening card, not a current correlated electroweak fit. Finite gauge constants, two-loop running and scalar-action realization are not included.","",
        "For `t=log(MI/MZ)`, `u=log(MU/MI)`, solve the three linear equations", "",
        "`alpha_SM^-1(MZ) = alphaU^-1 1 + b_SM t/(2 pi) + P b_PS u/(2 pi) - P sum_A delta_b_A log(kappa_A)/(2 pi)`,", "",
        "with `P(4,L,R)=(4,L,2*4/5+3*R/5)` and `delta_b_A=T_real,A/6`. Each logarithm belongs to an entire declared PS parent; no independent per-gauge compensator is introduced. The nine cards vary whole `Phi54:(20prime,1,1)` and `Phi54:(1,3,3)` reciprocally by `(1/3,1,3)`, and whole `Sigma126:(15,2,2)` by `(1,2,3)` at MI. Other stage masses stay at their reference scale; no Goldstone parent is moved. These are provisional mass hypotheses, not proofs that all nine spectra come from a stable scalar action.","",
        "All cards preserve left-right parity. Therefore `5 alpha1^-1-3 alpha2^-1-2 alpha3^-1` fixes MI independently of these logs. A narrow/vertical MI locus is a structural LO result, not a lack of scan coverage.","",
        "## Scale results", "",f"Baseline: `MI={b['MI_GeV']:.6g} GeV`, `MX=MU={b['MU_GeV']:.6g} GeV`, `alphaU^-1={b['alphaU_inverse']:.6f}`; `g4(MI)={b['g4I']:.6f}`, `sigma={b['sigma_GeV']:.6g} GeV`. The actual 45-generator orbit verifies `MX=gU*omega`, `MXprime=gU*sqrt(omega^2+sigma^2)` and charged PS-vector mass `g4(MI)*sigma`. The latter uses the intermediate coupling, not gU.","",
        "| Card | log10 MI | log10 MX | alphaU inverse | effective flavor norm limit |", "|---|---:|---:|---:|---:|"]
    lines += [f"| {c['id']} | {c['log10_MI']:.5f} | {c['log10_MX']:.5f} | {c['alphaU_inverse']:.4f} | {c['allowed_effective_flavor_norm']:.4f} |" for c in r["scenarios"]]
    lines += ["", "## What the proton number means", "",
        "Use `C_L,R(2 GeV)=gU^2/(2 MX^2) A_L,R F_L,R` for scalar three-quark operators. F includes the Fierz/Clebsch factors, fitted flavor rotations and the coherent second-vector term weighted by `MX^2/MXprime^2`; F=1 is merely a coefficient normalization, not a model prediction. The one-loop gauge RG factors are integrated with the same PS parent thresholds and SM running, followed by QCD to 2 GeV. [Gauge normalization and exponents](https://arxiv.org/html/1507.06712v2).", "",
        "The direct continuum lattice value is `W_pi+=-0.159(15)(20)(25) GeV^2`; isospin gives `|W_pi0|=|W_pi+|/sqrt(2)` in MSbar at 2 GeV. No extra chiral-Lagrangian `(1+D+F)/f_pi` multiplier is applied. [Lattice Table 8](https://arxiv.org/html/2111.01608v1).", "",
        "`Gamma_e_pi = mp/(32 pi) (1-mpi^2/mp^2)^2 |W_pi0|^2 [|C_L|^2+|C_R|^2]`.", "",
        f"Using the published `tau/Br > {TAU_LIMIT:.2g} yr` limit gives `sqrt(|C_L|^2+|C_R|^2)<{r['proton_input']['scalar_operator_coefficient_norm_bound_GeVminus2']:.6g} GeV^-2` at the central hadronic input. [Super-K original measurement](https://arxiv.org/abs/2010.16098).", "",
        "Equivalently every card must obey `sqrt(A_L^2 |F_L|^2 + A_R^2 |F_R|^2) < saved limit`. The table is this necessary condition in the stated gauge-only approximation, not a green/excluded model map. A pure-channel diagnostic and one-standard-error hadronic sensitivity are exported in JSON, but no joint confidence region is implied.", "",
        "## Explicit decisions", "", "- A nonordered scale solution or gauge Landau pole fails the declared scale card.", "- A specified flavor realization that violates the saved inequality fails this gauge-only channel test; F cannot be reset to zero to declare success.", "- Whole-model exclusion needs the correlated allowed flavor domain and all relevant amplitudes/channels. Scalar exchange and interference are currently uncomputed, so no unconditional exclusion is made.", "- All nine cards remain `physical_candidate=false`: vacuum realization, scalar-mediated proton decay, other channels and the common constrained flavor fit have not passed.", "", "Next use the saved sigma and correlated vector masses in the LO flavor screen, or test whether the action can realize a promising parent-mass pattern. Do not restore the obsolete P2 beta table or old MU threshold result.",""]
    return "\n".join(lines)


if __name__=="__main__":
    r=run();OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n");OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps(dict(passed=r["passed"],total=r["total"],baseline={k:r["baseline"][k] for k in ("MI_GeV","MX_GeV","alphaU_inverse","sigma_GeV","allowed_effective_flavor_norm")})))
    if r["passed"]!=r["total"]:raise SystemExit(1)
