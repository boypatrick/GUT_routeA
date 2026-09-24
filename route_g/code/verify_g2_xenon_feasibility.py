#!/usr/bin/env python3
"""Conditional boosted-X/Xe feasibility with the public XENON1T NR response.

No halo assumption, observed-event likelihood, production flux, or exclusion
claim is supplied. Only the X-only matching-scale Higgs-portal candidate is used.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import spherical_jn

import verify_g2_resolution as source
import verify_g2_driven as driven

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"data"/"xenon1t_s2only"
OUT = ROOT/"output"/"g2_xenon_feasibility"
CARD = dict(mh_GeV=125.08, v_GeV=246.0, Gamma_SM_GeV=.00410,
            invisible_branching_allowance=.107, fN=.308, fN_error=.018,
            mN_GeV=.939, atomic_mass_GeV=.93149410242,
            molar_mass_g_mol=131.293, Avogadro=6.02214076e23,
            GeV_minus2_to_cm2=3.893793721e-28, hbarc_GeV_fm=.1973269804,
            exposure_kg_day=356770.0, seconds_per_day=86400.0,
            recoil_window_keV=[.7, 50.0], requested_S2_window_PE=[150.0, 3000.0],
            lambda_X_scan_min=1e-5, lambda_X_scan_cap=1.0, expected_event_goal=3.0,
            source_sigma_i=.1, lambda_Phi_at_matching=0.0)
ISOTOPES = np.array([124,126,128,129,130,131,132,134,136])
ABUNDANCE = np.array([.00095,.00089,.01910,.26401,.04071,.21232,.26909,.10436,.08857])
ANCHORS = [1.0,2.0,5.0,10.0,12.0,12.508,15.0,30.0]
RAW_HASHES = {"s2_response_nr.csv":"76098a792eaf576ac835fcab7c412b2cc7763f2808c93db93998e256684dc202",
              "s2_binning_info.csv":"6ee893b8dff0a92848369026ec34ddfbeeac583ed6b9cdb5a9b2e1331cdbd341"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def helm_squared(q_GeV, A):
    q = np.asarray(q_GeV)/CARD["hbarc_GeV_fm"]
    c, a, skin = 1.23*A**(1/3)-.60, .52, .9
    radius = np.sqrt(c*c+7*np.pi*np.pi*a*a/3-5*skin*skin)
    z = q*radius
    ratio = np.ones_like(z, dtype=float)
    np.divide(3*spherical_jn(1,z), z, out=ratio, where=abs(z)>1e-8)
    ratio = np.where(abs(z)<=1e-8, 1-z*z/10, ratio)
    return ratio*ratio*np.exp(-(q*skin)**2)


def endpoint(p, mX, mA):
    energy = np.hypot(mX,p)
    s = mX*mX+mA*mA+2*mA*energy
    return 2*mA*p*p/s


def nuclear_H(ER_GeV, A, use_helm=True, fN=None):
    """d sigma/dER = lambda_X^2 / p_X^2 * H(ER) within the endpoint."""
    fN = CARD["fN"] if fN is None else fN
    mA = A*CARD["atomic_mass_GeV"]
    ff = helm_squared(np.sqrt(2*mA*np.asarray(ER_GeV)), A) if use_helm else 1.0
    return fN*fN*CARD["mN_GeV"]**2*A*A*mA*ff/(8*np.pi*(CARD["mh_GeV"]**2+2*mA*np.asarray(ER_GeV))**2)


def dsigma_dER(ER, p, mX, A, coupling=1.0, use_helm=True):
    mA = A*CARD["atomic_mass_GeV"]
    return np.where((np.asarray(ER)>=0)&(np.asarray(ER)<=endpoint(p,mX,mA)),
                    coupling*coupling*nuclear_H(ER,A,use_helm)/p**2,0.0)


def response_integrals(energy_rows, order=24):
    nodes, weights = leggauss(order)
    out = np.zeros((len(ISOTOPES),len(energy_rows)))
    for i,row in enumerate(energy_rows):
        lo = max(float(row["energy_bin_start_kev"]),CARD["recoil_window_keV"][0])*1e-6
        hi = min(float(row["energy_bin_end_kev"]),CARD["recoil_window_keV"][1])*1e-6
        if hi<=lo:
            continue
        energies = (lo+hi)/2+(hi-lo)/2*nodes
        for k,A in enumerate(ISOTOPES):
            out[k,i] = float(np.dot(weights,nuclear_H(energies,int(A)))*(hi-lo)/2)
    return out


def tower_width_unit(scale):
    """Width/lambda_X^2 for all signed complex X modes, with fixed physical MD."""
    if scale>=CARD["mh_GeV"]/10:
        return dict(open_signed_modes=[], beta_sum=0.0, unit_width_GeV=0.0)
    radius = np.sqrt(max((CARD["mh_GeV"]/(2*scale))**2-25,0))
    labels = range(-int(np.ceil(radius)),int(np.ceil(radius))+1)
    rows = []
    for j in labels:
        mass = scale*np.sqrt(25+j*j)
        if 2*mass<CARD["mh_GeV"]:
            beta = np.sqrt((1-2*mass/CARD["mh_GeV"])*(1+2*mass/CARD["mh_GeV"]))
            rows.append(dict(j=j,mass_GeV=mass,beta=float(beta)))
    beta_sum = sum(r["beta"] for r in rows)
    return dict(open_signed_modes=rows,beta_sum=beta_sum,
                unit_width_GeV=CARD["v_GeV"]**2*beta_sum/(16*np.pi*CARD["mh_GeV"]))


def run():
    checks = []

    def check(name,passed,value=None,tolerance=None):
        row = dict(name=name,passed=bool(passed))
        if value is not None:row["measured"]=float(value)
        if tolerance is not None:row["tolerance"]=float(tolerance)
        checks.append(row)

    dependencies = {"resolution_source":Path(source.__file__),"driven_source":Path(driven.__file__),
                    "g2_source":Path(source.g2.__file__),
                    "resolution_artifact":ROOT/"output"/"g2_resolution.json"}
    dependency_hashes = {key:sha(path) for key,path in dependencies.items()}
    for name,expected in RAW_HASHES.items():
        check("data_hash:"+name,sha(DATA/name)==expected)
    old = json.loads(dependencies["resolution_artifact"].read_text())
    check("source_hash:resolution",old["verification"]["source_sha256"]==dependency_hashes["resolution_source"])
    check("source_hash:g2",old["verification"]["g2_source_sha256"]==dependency_hashes["g2_source"])
    energy_rows = np.genfromtxt(DATA/"s2_response_nr.csv",delimiter=",",names=True)
    bins = np.genfromtxt(DATA/"s2_binning_info.csv",delimiter=",",names=True)
    response = np.column_stack([energy_rows[f"s2_bin_{i:03d}"] for i in range(len(bins))])
    selected = (bins["start_pe"]>=150)&(bins["end_pe"]<=3000)
    columns = np.flatnonzero(selected)
    check("response:nonnegative",np.min(response)>=0)
    check("response:probability_not_efficiency_twice",np.max(response.sum(axis=1))<=1+1e-12)
    check("response:one_to_one_bin_keys",np.array_equal(bins["s2_bin_number"],np.arange(len(bins))))
    check("response:contiguous_bins",np.max(abs(bins["end_pe"][:-1]-bins["start_pe"][1:]))<1e-8)
    check("response:energy_edges_contiguous",np.max(abs(energy_rows["energy_bin_end_kev"][:-1]-energy_rows["energy_bin_start_kev"][1:]))<1e-12)
    check("response:ROI_supported",energy_rows["energy_bin_start_kev"].min()<=.7 and energy_rows["energy_bin_end_kev"].max()>=50)
    check("response:actual_S2_lower_edge",abs(bins["start_pe"][columns[0]]-150.027)<.001)
    check("isotopes:number_fractions_sum_one",abs(ABUNDANCE.sum()-1)<1e-14)
    for A in ISOTOPES:
        check(f"Helm:A={A}:zero_momentum",abs(float(helm_squared(0,int(A)))-1)<1e-14)
        vals=helm_squared(np.linspace(0,.25,501),int(A))
        check(f"Helm:A={A}:bounded",np.min(vals)>=0 and np.max(vals)<=1+1e-13)
    integrals = response_integrals(energy_rows)
    coarse_integrals = response_integrals(energy_rows,order=12)
    nuclear_bins = ABUNDANCE@integrals@response[:,selected]
    coarse_nuclear_bins = ABUNDANCE@coarse_integrals@response[:,selected]
    G = float(nuclear_bins.sum())
    check("response:quadrature_convergence",np.max(abs(nuclear_bins-coarse_nuclear_bins))/G<1e-12)
    # Independent integration of the row-summed acceptance, including keV->GeV.
    direct = 0.0
    acceptance = response[:,selected].sum(axis=1)
    for row,eps in zip(energy_rows,acceptance):
        lo=max(float(row["energy_bin_start_kev"]),.7)*1e-6
        hi=min(float(row["energy_bin_end_kev"]),50)*1e-6
        if hi<=lo:continue
        direct += eps*quad(lambda ER:sum(float(frac)*float(nuclear_H(ER,int(A))) for A,frac in zip(ISOTOPES,ABUNDANCE)),lo,hi,epsabs=1e-24,epsrel=1e-10)[0]
    check("response:independent_acceptance_integral",abs(direct/G-1)<1e-10,abs(direct/G-1),1e-10)
    # Constant spectrum integration checks clipping rather than center-based cuts.
    widths=np.maximum(0,np.minimum(energy_rows["energy_bin_end_kev"],50)-np.maximum(energy_rows["energy_bin_start_kev"],.7))
    check("response:clipped_energy_width",abs(widths.sum()-49.3)<1e-12)
    check("response:lower_edge_partial_bin",np.any((energy_rows["energy_bin_start_kev"]<.7)&(energy_rows["energy_bin_end_kev"]>.7)))
    source_rows,_,_,tail = source.components(.1,True)
    coarse_rows,_,_,_=source.components(.1,True,order=24)
    Z=sum(float(r["rates"].sum()) for r in source_rows)
    def moments(rows):
        z=sum(float(r["rates"].sum()) for r in rows)
        return np.array([sum(float(np.dot(r["rates"],1/r["means"]**2)) for r in rows if r["sign"]==sign)/z for sign in [1,-1]])
    B=moments(source_rows)
    Bcoarse=moments(coarse_rows)
    check("source:inverse_momentum_moment_convergence",np.max(abs(B-Bcoarse))<1e-13)
    oldpoint=next(r for r in old["scan"] if r["sigma_i"]==.1 and r["sigma_d"]==.05 and r["driven"])
    check("source:rate_recovery",abs(Z/oldpoint["total_pair_rate"]-1)<1e-13)
    source_priors=np.array([sum(float(r["rates"].sum()) for r in source_rows if r["sign"]==sign)/Z for sign in [1,-1]])
    detected_priors=B/B.sum()
    profile=nuclear_bins/G
    joint=detected_priors[:,None]*profile[None,:]
    s2_error=float(np.minimum(joint[0],joint[1]).sum())
    check("S2:factorization_normalized",abs(joint.sum()-1)<1e-13)
    check("S2:same_conditional_shape",np.max(abs(joint[0]/detected_priors[0]-joint[1]/detected_priors[1]))<1e-14)
    check("S2:no_sign_information",abs(s2_error-min(detected_priors))<1e-14)
    minimum_rows=[driven.driven_channel(1,l,0,q,0) for l in [0,-2] for q in [-1,0,1]]
    pmin=min(r["outgoing_momentum"] for r in minimum_rows)
    exposure_factor=CARD["seconds_per_day"]*CARD["exposure_kg_day"]*(1000*CARD["Avogadro"]/CARD["molar_mass_g_mol"])
    unit_sigma_cm2=float(B.sum()*G*CARD["GeV_minus2_to_cm2"])
    unit_count_coefficient=exposure_factor*unit_sigma_cm2
    Gamma_max=CARD["invisible_branching_allowance"]/(1-CARD["invisible_branching_allowance"])*CARD["Gamma_SM_GeV"]
    scale_grid=np.unique(np.r_[np.geomspace(1,30,121),ANCHORS])
    scans=[]
    for scale in scale_grid:
        tower=tower_width_unit(float(scale))
        width_bound=np.sqrt(Gamma_max/tower["unit_width_GeV"]) if tower["unit_width_GeV"]>0 else None
        stability=50*scale*scale/CARD["v_GeV"]**2
        limits={"scan_cap":1.0,"bare_mass_stability":float(stability)}
        if width_bound is not None:limits["conditional_invisible_width"]=float(width_bound)
        active=min(limits,key=limits.get)
        maximum=limits[active]
        coefficient=unit_count_coefficient/scale**2
        minimum_endpoint=min(endpoint(pmin*scale,np.sqrt(26)*scale,float(A)*CARD["atomic_mass_GeV"]) for A in ISOTOPES)*1e6
        check(f"Lambda={scale:.10g}:all_endpoints_above_ROI",minimum_endpoint>50)
        check(f"Lambda={scale:.10g}:matched_bare_mass_nonnegative",25*scale**2-maximum*CARD["v_GeV"]**2/2>=-1e-10)
        check(f"Lambda={scale:.10g}:width_bound",maximum**2*tower["unit_width_GeV"]<=Gamma_max*(1+1e-12))
        # Compare scale with the exact threshold expression to avoid promoting
        # floating-point multiplication noise into a spurious zero-velocity mode.
        check(f"Lambda={scale:.10g}:signed_tower_complete",len(tower["open_signed_modes"])==sum(scale<CARD["mh_GeV"]/(2*np.sqrt(25+j*j)) for j in range(-1000,1001)))
        scans.append(dict(Lambda_GeV=float(scale),mX_absj1_GeV=float(np.sqrt(26)*scale),
                          minimum_recoil_endpoint_keV=float(minimum_endpoint),
                          width_unit_GeV=tower["unit_width_GeV"],width_beta_sum=tower["beta_sum"],
                          width_open_signed_mode_count=len(tower["open_signed_modes"]),
                          lambda_stability_max=float(stability),lambda_width_max=width_bound,
                          lambda_allowed_conditional_max=maximum,active_constraint=active,
                          accepted_sigma_cm2_per_lambda_squared=unit_sigma_cm2/scale**2,
                          counts_per_incident_flux_per_lambda_squared=coefficient,
                          minimum_flux_for_three_cm2_s=3/(coefficient*maximum**2),
                          required_flux_for_three_at_lambda1_cm2_s=3/coefficient,
                          pump_energy_quantum_GeV=.2*float(scale)))
    threshold=CARD["mh_GeV"]/10
    check("Higgs:tower_closed_at_threshold",tower_width_unit(threshold)["beta_sum"]==0)
    check("Higgs:tower_open_below_threshold",tower_width_unit(threshold*(1-1e-8))["beta_sum"]>0)
    # Independent finite-velocity phase-space normalization and NR limit.
    kin_examples=[]
    for A in [1,132]:
        mA=A*CARD["atomic_mass_GeV"]
        mX=np.sqrt(26)
        for speed in [1e-4,.2]:
            p=mX*speed/np.sqrt(1-speed**2)
            Emax=endpoint(p,mX,mA)
            analytic=(CARD["fN"]**2*CARD["mN_GeV"]**2*A*A*mA/(8*np.pi*p*p)
                      *Emax/(CARD["mh_GeV"]**2*(CARD["mh_GeV"]**2+2*mA*Emax)))
            integrated=quad(lambda u:float(dsigma_dER(Emax*u,p,mX,A,use_helm=False))*Emax,0,1,epsabs=1e-25,epsrel=1e-11)[0]
            reduced=mX*mA/(mX+mA)
            nr=CARD["fN"]**2*CARD["mN_GeV"]**2*A*A*reduced**2/(4*np.pi*CARD["mh_GeV"]**4*mX*mX)
            check(f"pointlike:A={A}:v={speed}:analytic_integral",abs(integrated/analytic-1)<1e-11)
            if speed==1e-4:check(f"NR:A={A}:limit",abs(integrated/nr-1)<1e-7)
            ER=Emax*.4
            s=mX*mX+mA*mA+2*mA*np.hypot(mX,p)
            kallen=4*mA*mA*p*p
            M2=4*CARD["fN"]**2*CARD["mN_GeV"]**2*A*A*mA*mA/(CARD["mh_GeV"]**2+2*mA*ER)**2
            full=2*mA*M2/(16*np.pi*kallen)
            check(f"finite_velocity:A={A}:v={speed}:normalization",abs(float(dsigma_dER(ER,p,mX,A,use_helm=False))/full-1)<1e-13)
            kin_examples.append(dict(A=A,speed=speed,exact_pointlike_sigma_GeV_minus2=integrated,nr_sigma_GeV_minus2=nr,exact_over_NR=integrated/nr))
    for coupling in [.001,.1,1.0]:
        ratio=float(dsigma_dER(1e-5,1.1,np.sqrt(26),132,coupling)/dsigma_dER(1e-5,1.1,np.sqrt(26),132))
        check(f"coupling:{coupling}:square_scaling",abs(ratio-coupling*coupling)<1e-13)
    for A in [124,132,136]:
        for recoil in [.7e-6,12e-6,50e-6]:
            left=float(dsigma_dER(recoil,.809144464424409,np.sqrt(26),A))*.809144464424409**2
            right=float(dsigma_dER(recoil,1.6,np.sqrt(26),A))*1.6**2
            check(f"S2:direct_kernel_factorization:A={A}:ER={recoil}",abs(left/right-1)<1e-13)
    fN_nuisance=dict(fN=CARD["fN"],delta_fN=CARD["fN_error"],
                     yield_multipliers=[((CARD["fN"]+sign*CARD["fN_error"])/CARD["fN"])**2 for sign in [-1,1]],
                     flux_requirement_multipliers=[(CARD["fN"]/(CARD["fN"]+sign*CARD["fN_error"]))**2 for sign in [-1,1]],
                     scope="Only the supplied nucleon scalar-form-factor uncertainty; not total nuclear, detector, EFT or flux uncertainty.")
    check("dependencies:unchanged",all(sha(path)==dependency_hashes[key] for key,path in dependencies.items()))
    check("data:unchanged",all(sha(DATA/name)==expected for name,expected in RAW_HASHES.items()))
    failed=[r["name"] for r in checks if not r["passed"]]
    anchors=[next(r for r in scans if r["Lambda_GeV"]==v) for v in ANCHORS]
    result=dict(status="Conditional X-only candidate feasibility with a real public NR response, not a realized beam, observed signal, exclusion or discovery sensitivity.",
                card=CARD,source_drive=source.DRIVE,isotopes=dict(A=ISOTOPES.tolist(),number_fractions=ABUNDANCE.tolist(),mass_approximation="mA=A*0.93149410242 GeV; isotope/electron/binding mass corrections neglected."),
                input_response=dict(directory="route_g/data/xenon1t_s2only",sha256=RAW_HASHES,
                                    interpretation="NR matrix already contains all full-volume averaged analysis selections. No second fiducial/efficiency factor.",
                                    energy_interpolation="Published response row is piecewise constant over its supplied energy-bin edges; integrate theory within clipped [.7,50] keV bin.",
                                    actual_selected_S2_PE=[float(bins["start_pe"][columns[0]]),float(bins["end_pe"][columns[-1]])],
                                    selected_bin_numbers=columns.tolist(),energy_row_count=len(energy_rows),full_S2_bin_count=len(bins)),
                normalization=dict(exposure_target_seconds=exposure_factor,unit_lambda1_Lambda1_sigma_cm2=unit_sigma_cm2,
                                   counts_per_flux_at_lambda1_Lambda1=unit_count_coefficient,
                                   flux_definition="Stationary uniform fluence rate of the selected m=0,|j|=1 X population at the detector in particles/cm^2/s; no halo density or source production flux inferred.",
                                   event_equation="N=Phi_pair *86400*356770*(1000*N_A/131.293)*<sigma_acc_cm2>; <sigma_acc>=lambda_X^2/Lambda_GeV^2 * unit_sigma.",
                                   flux_goal="Phi_3=3/[counts_per_flux_per_lambda_squared*lambda_X^2]; three expected events, not a confidence limit or discovery criterion."),
                factorization=dict(pmin_hat=pmin,minimum_endpoint_over_scan_keV=min(r["minimum_recoil_endpoint_keV"] for r in scans),
                                   source_event_priors=source_priors.tolist(),inverse_squared_momentum_joint_moments=B.tolist(),
                                   post_detection_sign_priors=detected_priors.tolist(),S2_only_bayes_error=s2_error,
                                   S2_only_mutual_information_bits=0.0,
                                   theorem="For the entire accepted recoil window, all exact endpoints are above50keV. Thus d sigma/dER=lambda_X^2*H_A(ER)/p_X^2. Folding an ER-only response yields P(S2bin|j=+)=P(S2bin|j=-), so S2 carries no sign information; only rate-weighted priors change.",
                                   continuous_scale_coverage="At fixed p_hat,m_hat and nucleus, endpoint(Lambda)=2mA*p_hat^2*Lambda^2/(mA^2+2mA*Ehat*Lambda+m_hat^2*Lambda^2), whose derivative is4mA*p_hat^2*Lambda*(mA^2+mA*Ehat*Lambda)/denominator^2>0. Checking Lambda=1 for every isotope bounds the whole1..30GeV interval, not merely scan points.",
                                   no_promotion="The old ideal momentum-classifier error2.20% is not an S2 recoil-energy performance prediction.",
                                   beam_tail_probability=tail,beam_tail_qualification="Formal-tree Gaussian quadrature z<=10 without renormalization. Endpoint bound uses exact p_i>=0 minimum, independent of that truncation; unknown UV validity is not certified."),
                S2_spectrum=dict(bin_numbers=columns.tolist(),start_PE=bins["start_pe"][selected].tolist(),end_PE=bins["end_pe"][selected].tolist(),
                                 center_PE=bins["log_center_pe"][selected].tolist(),common_conditional_probability=profile.tolist(),
                                 conditional_probability_by_true_sign=[profile.tolist(),profile.tolist()],
                                 joint_sign_bin_probability=joint.tolist(),nuclear_integral_per_bin=nuclear_bins.tolist()),
                Higgs_constraint=dict(Gamma_tower_max_GeV=Gamma_max,last_open_threshold_Lambda_GeV=threshold,
                                      formula="Gamma_X=lambda_X^2*v^2/(16*pi*mh)*sum_signed_j sqrt(1-4*Lambda^2*(25+j^2)/mh^2), open modes only.",
                                      scope="ATLAS invisible branching allowance imported conditionally: SM production/visible width and invisible/escaping tower final states must apply; EFT must support external mh and all contributing modes. Above threshold no constraint follows from this width channel."),
                fN_nuisance=fN_nuisance,kinematic_checks=kin_examples,scan=scans,anchors=anchors,
                assumptions=["Physical scale1..30GeV and X-only portal are explicit candidate choices, not values derived by Route G.",
                             "Physical MD=5Lambda is retained by bare mass matching. With no X self-quartic, lambda_X<=50Lambda^2/v^2 avoids negative bare mass; the scan cap1 is a chosen range, not a theorem of all-orders perturbativity.",
                             "lambda_Phi=0 holds at a stated matching scale, not to all orders. Loop, gauge/UV completion and portal backreaction are not certified.",
                             "Source COM is identified with target lab, the output ensemble is incoherent and selected externally, flux is stationary/spatially uniform at the detector, and attenuation/multiple scattering are neglected.",
                             "Source compact labels are ideal selections. No actual mode-selecting preparation or event-level j measurement is provided.",
                             "The nuclear scalar-current/Helm approximation keeps finite incident velocity and exact elastic endpoint. It is not a full nuclear EFT treatment; isotope spin-dependent effects and binding uncertainties are not included.",
                             "The published NR matrix is used for nuclear recoil energy only. No actual source timing, directional information, backgrounds, counts or statistical likelihood are fitted.",
                             "Events outside.7..50keV and selected complete S2 bins are not modeled as accepted here. This is a deliberately specified response window, not all possible detector information."],
                checks=checks,summary=dict(checks=len(checks),passed=len(checks)-len(failed),failed=failed),
                verification=dict(source_sha256=sha(__file__),dependencies=dependency_hashes,data_sha256=RAW_HASHES))
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    lines=["# G2 xenon feasibility: counts per incident flux and a failed S2 sign readout", "",
           f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**. The matrix is a published experimental response; the source and portal remain a hypothetical candidate. No actual beam flux or signal confidence is claimed.", "",
           "The chosen physical card is Lambda=1..30GeV, mX=√26Lambda, mh=125.08GeV, v=246GeV, fN=.308±.018. The portal couples X only at matching. Physical MD=5Lambda is held fixed by bare matching; no new X stabilizing quartic is introduced.", "",
           "## The physical obstruction", "",
           f"The source minimum momentum is {pmin:.6f}Lambda. Every exact recoil endpoint in the scan exceeds {result['factorization']['minimum_endpoint_over_scan_keV']:.1f}keV, above the entire0.7–50keV window. Consequently dσ/dER=lambda_X² H_A(ER)/pX² throughout the window. The two signs have **identical normalized S2 spectra** after the real response is folded.", "",
           f"Source sign priors: +{100*source_priors[0]:.3f}%, −{100*source_priors[1]:.3f}%. Accepted-event priors: +{100*detected_priors[0]:.3f}%, −{100*detected_priors[1]:.3f}%. The best S2-only error is **{100*s2_error:.3f}%**, simply guessing the accepted majority sign; no S2-bin classifier improves it in this model/window. The previous2.20% ideal momentum result does not describe a xenon energy-deposit measurement.", "",
           "## Response and rate normalization", "",
           f"Use complete S2 bins from {result['input_response']['actual_selected_S2_PE'][0]:.3f} to {result['input_response']['actual_selected_S2_PE'][1]:.0f}PE. Nuclear-recoil response rows are piecewise constant over supplied energy-bin edges, clipped to0.7–50keV; the theoretical kernel is integrated within each bin. All full-volume averaged analysis selections are already in the matrix. Do not multiply another fiducial efficiency.", "",
           "N = Phi_pair ×86400×356770kg-days×(1000 N_A/131.293)×<sigma_acc_cm²>. Phi_pair is the stationary uniform selected-particle flux at the detector in cm⁻²s⁻¹, not source production or a halo prediction. No additional velocity factor is used. Three expected events require Phi_3=3/(N/Phi_pair); this is not a discovery or exclusion threshold.", "",
           "## Conditional coupling budget and required flux", "",
           "lambda_max=min(1,50Lambda²/v²,lambda_Higgs when a tower mode is open). The entire signed complex-X tower is counted in the width. The ATLAS B_inv=.107 allowance implies Gamma_X<=B_inv/(1-B_inv)×.00410GeV only under its invisible-channel and SM-width premises. At Lambda>=12.508GeV this particular Higgs-decay constraint disappears; stability and the chosen scan cap remain.", "",
           "| Lambda [GeV] | mX [GeV] | Conditional lambda_max | Active bound | Minimum flux for3 expected events [cm⁻²s⁻¹] |",
           "|---:|---:|---:|---|---:|"]
    for row in anchors:
        lines.append(f"| {row['Lambda_GeV']:g} | {row['mX_absj1_GeV']:.3f} | {row['lambda_allowed_conditional_max']:.5g} | {row['active_constraint']} | {row['minimum_flux_for_three_cm2_s']:.5g} |")
    lines += ["",f"At lambda_X=1,Lambda=1GeV the accepted effective cross section is {unit_sigma_cm2:.6g}cm²; scale it by lambda_X²/Lambda_GeV². This normalization point itself need not satisfy the physical coupling constraints.", "",
              f"The fN-only yield factors are {fN_nuisance['yield_multipliers'][0]:.4f} and {fN_nuisance['yield_multipliers'][1]:.4f}; reciprocal flux factors are {fN_nuisance['flux_requirement_multipliers'][0]:.4f} and {fN_nuisance['flux_requirement_multipliers'][1]:.4f}. These are not total uncertainties.", "",
              "## What is and is not established", "",
              "The calculation checks coherent scalar-current normalization, the low-speed limit, exact finite-velocity endpoints, isotope-number weighting, Helm F(0)=1, response clipping, independent energy integration, source regression and full signed-tower thresholds. Natural-xenon isotope masses use A×.93149410242GeV; binding/electron corrections and nuclear theory uncertainties are not promoted into a precision model.", "",
              "Useful event counts still require an independently attainable detector-incident flux, a controlled UV/EFT and an applicable invisible-width interpretation. Public detector calibration does not provide that source. An energy-only ROI below all endpoints is not a momentum spectrometer: recovering sign information needs an independent observable or physically accessible endpoint information, not further fitting of a Gaussian error width.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"]))
    print(json.dumps(dict(priors=detected_priors.tolist(),S2error=s2_error,unit_sigma_cm2=unit_sigma_cm2,anchors=anchors)))
    if failed:raise SystemExit("Failed checks: "+", ".join(failed))


if __name__=="__main__":
    run()
