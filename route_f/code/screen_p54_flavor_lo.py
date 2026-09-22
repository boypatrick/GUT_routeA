#!/usr/bin/env python3
"""Bounded physics-first three-family, two-matrix leading-order screen.

The common-scale tree relations are a deliberately permissive flavor
hypothesis, NOT a matched P54 UV prediction.  Charged targets use one-loop
SM transport only below MI.  Neutrino angle/ratio running, PS Yukawa logs,
finite thresholds and type II are omitted.  No scalar quartic is touched.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
import verify_p54_p3_flavor_cw_gate as inputs
import verify_p54_sequential_seesaw as seesaw

RF=Path(__file__).resolve().parents[1]
OUT=RF/'output/p54_flavor_lo_screen'
MI=5.082267472797e13
G4I=.570524280319144
SIGMA=MI/G4I
V174=174.
TRIU=np.triu_indices(3)
HSCALE=np.array([1e-6,1e-4,.005])
FSCALE=np.array([3e-6,1e-5,1e-4,1e-4,3e-4,1e-3])
ANGLES=np.array([.308,.470,.02215]) # theta12,theta23,theta13
NU_BANDS=np.array([[.275,.345],[.435,.585],[.02030,.02388]])
DM21=7.49e-5
DM31=2.513e-3
DM21_BAND=np.array([6.92e-5,8.05e-5])
DM31_BAND=np.array([2.451e-3,2.578e-3])
SOURCES=[
 {'url':'https://arxiv.org/html/2109.04050v2','use':'Table II MSbar charged masses and gauge inputs at 173.1 GeV; Appendix B provides numerical starts ONLY, originally a different 2HDM fit'},
 {'url':'https://pdg.lbl.gov/2024/reviews/rpp2024-rev-ckm-matrix.pdf','use':'frozen PDG 2024 standard CKM central inputs already encoded in the legacy transport helper'},
 {'url':'https://www.nu-fit.org/sites/default/files/v60.tbl-parameters.pdf','use':'NuFIT 6.0 IC24 with SK, normal ordering: three mixing angles and two mass-squared differences; frozen 2024 input, not claimed to be the latest fit'},
 {'url':'https://arxiv.org/abs/2410.05380','use':'primary NuFIT 6.0 methodology and conditional ordering choice'},
]

def cjson(a):
    a=np.asarray(a);return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def decode(a):return np.asarray(a['real'])+1j*np.asarray(a['imag'])

def unpack(x):
    h=np.diag(x[:3]*HSCALE).astype(complex)
    f=np.zeros((3,3),complex);f[TRIU]=(x[3:9]+1j*x[9:15])*FSCALE
    f=f+f.T-np.diag(np.diag(f))
    return h,f,np.exp(x[15]),x[16]+1j*x[17]

def pack(h,f,r,s):
    return np.r_[h.diagonal().real/HSCALE,f[TRIU].real/FSCALE,f[TRIU].imag/FSCALE,np.log(r),s.real,s.imag]

def species(h,f,r,s):
    return {'up':r*(h+s*f),'down':h+f,'charged_lepton':h-3*f,'neutrino_dirac':r*(h-3*s*f)}

def diagonalize(y):
    u,m,_=np.linalg.svd(y)
    return m[::-1],u[:,::-1]

def mixing(v):
    mod=abs(v);s13=mod[0,2];c13=np.sqrt(max(1-s13*s13,1e-30))
    s12=np.clip(mod[0,1]/c13,0,1);s23=np.clip(mod[1,2]/c13,0,1)
    c12=np.sqrt(1-s12*s12);c23=np.sqrt(1-s23*s23)
    j=float(np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()))
    denominator=max(s12*s23*s13*c12*c23*c13*c13,1e-30)
    sin_delta=np.clip(j/denominator,-1,1)
    cos_delta=np.clip((s12*s12*s23*s23+c12*c12*c23*c23*s13*s13-mod[2,0]**2)/
                      max(2*s12*s23*c12*c23*s13,1e-30),-1,1)
    return np.array([s12,s23,s13]),float(np.arctan2(sin_delta,cos_delta)),j

def observables(x):
    h,f,r,s=unpack(x);ys=species(h,f,r,s)
    masses={};left={}
    for key in ('up','down','charged_lepton'):
        masses[key],left[key]=diagonalize(ys[key])
    # Consistent with the repository's two-left-Weyl common phase dictionary.
    ckm=left['up'].T@left['down'].conj()
    ckm_angles,ckm_delta,ckm_j=mixing(ckm)
    shape=ys['neutrino_dirac']@np.linalg.solve(f,ys['neutrino_dirac'].T)
    q,nu_left=diagonalize(shape)
    pmns=left['charged_lepton'].T@nu_left.conj()
    pmns_angles,pmns_delta,pmns_j=mixing(pmns)
    ratio=(q[1]**2-q[0]**2)/max(q[2]**2-q[0]**2,1e-100)
    return dict(h=h,f=f,r=r,s=s,ys=ys,masses=masses,ckm=ckm,ckm_angles=ckm_angles,
                ckm_delta=ckm_delta,ckm_j=ckm_j,shape=shape,q=q,pmns=pmns,
                pmns_sin2=pmns_angles**2,pmns_delta=pmns_delta,pmns_j=pmns_j,ratio=ratio)

def targets(mu):
    target=inputs.diagnostic_sm_run(mu)
    mod=np.asarray(target['CKM_abs_at_high'])
    # Initial PDG J is positive; reconstruct its physical branch from moduli.
    s13=mod[0,2];s12=mod[0,1]/np.sqrt(1-s13*s13);s23=mod[1,2]/np.sqrt(1-s13*s13)
    c12=np.sqrt(1-s12*s12);c23=np.sqrt(1-s23*s23)
    cd=(s12*s12*s23*s23+c12*c12*c23*c23*s13*s13-mod[2,0]**2)/(2*s12*s23*c12*c23*s13)
    target['CKM_s12_s23_s13']=np.array([s12,s23,s13]).tolist()
    target['CKM_delta_rad']=float(np.arccos(np.clip(cd,-1,1)))
    return target

def residual(x,target):
    try:
        obs=observables(x)
        values=np.r_[obs['masses']['up'],obs['masses']['down'],obs['masses']['charged_lepton']]
        desired=np.concatenate([target['yukawa_singular_values_at_high'][key] for key in ('up','down','charged_lepton')])
        # These are declared LO discrepancy widths, NOT measured 1-sigma
        # errors and not free additive threshold nuisance parameters.
        widths=np.array([.3,.1,.1,.3,.3,.1,.1,.1,.1])
        mass_res=np.log(np.maximum(values,1e-30)/desired)/widths
        angles=np.log(np.maximum(obs['ckm_angles'],1e-30)/target['CKM_s12_s23_s13'])/.10
        delta=np.angle(np.exp(1j*(obs['ckm_delta']-target['CKM_delta_rad'])))/.12
        nures=(obs['pmns_sin2']-ANGLES)/((NU_BANDS[:,1]-NU_BANDS[:,0])/6)
        ratio=np.log(max(obs['ratio'],1e-30)/(DM21/DM31))/.05
        return np.r_[mass_res,angles,delta,nures,ratio]
    except (np.linalg.LinAlgError,ValueError,FloatingPointError):
        return np.full(17,1e6)

def normalization_family(obs,kappa,sigma,hcap=1.,fcap=1.,comparison_scale=MI):
    h,f,r,s=(obs[key] for key in ('h','f','r','s'))
    hn=float(np.linalg.norm(h,2));fn=float(np.linalg.norm(f,2))
    a_weight=1+r*r;d_weight=1+abs(r*s)**2
    dmin=fn/fcap
    room=1-a_weight*(hn/hcap)**2
    dmax=np.sqrt(max(room,0)/d_weight)
    q=obs['q'];factor=np.sqrt(DM31/(q[2]**2-q[0]**2))
    sigma_per_d=V174**2*1e9/(abs(kappa)*factor)
    interval=np.array([dmin,dmax])*sigma_per_d
    required_d=sigma/sigma_per_d
    feasible=bool(room>0 and 0<dmin<=dmax)
    baseline=bool(feasible and dmin<=required_d<=dmax)
    chosen_d=required_d if baseline else (np.sqrt(dmin*dmax) if feasible else .5/np.sqrt(d_weight))
    a=np.sqrt(max(1-d_weight*chosen_d**2,0)/a_weight)
    b=r*a;d=chosen_d;e=r*s*d
    hd=h/max(a,1e-100);fd=f/max(d,1e-100)
    selected_sigma=sigma_per_d*d
    mr=kappa*selected_sigma*fd
    mnu=-V174**2*obs['ys']['neutrino_dirac']@np.linalg.solve(mr,obs['ys']['neutrino_dirac'].T)
    nu_eV=np.linalg.svd(mnu,compute_uv=False)[::-1]*1e9
    actual_mr=kappa*sigma*fd
    actual_mr_values=np.linalg.svd(actual_mr,compute_uv=False)[::-1]
    ynu_norm=float(np.linalg.norm(obs['ys']['neutrino_dirac'],2))
    log_estimate=ynu_norm**2/(16*np.pi**2)*np.maximum(np.log(comparison_scale/actual_mr_values),0.)
    return {
        'norm_cap_convention':'spectral norms of Dirac-unit h_D,f_D; engineering filter, not complete Spin(10) perturbative control',
        'h_D_norm_cap':hcap,'f_D_norm_cap':fcap,
        'exact_norm_overlap_existence_lhs':float(a_weight*(hn/hcap)**2+d_weight*(fn/fcap)**2),
        'normalized_overlap_family_exists_under_caps':feasible,
        'd_interval_under_caps':[float(dmin),float(dmax)] if feasible else None,
        'theta_interval_rad':[float(np.arcsin(np.clip(dmin*np.sqrt(d_weight),0,1))),float(np.arcsin(np.clip(dmax*np.sqrt(d_weight),0,1)))] if feasible else None,
        'sigma_interval_GeV_central_dm31':interval.tolist() if feasible else None,
        'sigma_interval_GeV_dm31_3sigma_envelope':
            [float(interval[0]*np.sqrt(DM31/DM31_BAND[1])),float(interval[1]*np.sqrt(DM31/DM31_BAND[0]))] if feasible else None,
        'sigma_over_d_GeV':float(sigma_per_d),
        'sigma_supremum_from_overlap_normalization_without_any_norm_caps_GeV':float(sigma_per_d/np.sqrt(d_weight)),
        'baseline_sigma_GeV':sigma,'baseline_required_d':float(required_d),
        'baseline_sigma_inside_this_witness_family':baseline,
        'baseline_failure_does_not_exclude_other_flavor_textures':True,
        'selected_sigma_GeV':float(selected_sigma),
        'selected_overlaps_abde':{key:cjson(value) for key,value in zip(('a','b','d','e'),(a,b,d,e))},
        'selected_h_D':cjson(hd),'selected_f_D':cjson(fd),
        'selected_h_D_norm':float(np.linalg.norm(hd,2)),'selected_f_D_norm':float(np.linalg.norm(fd,2)),
        'selected_raw_h_norm':float(np.linalg.norm(hd,2)/np.sqrt(2)),
        'selected_raw_f_norm':float(np.linalg.norm(fd,2)*np.sqrt(3)/2),
        'selected_Majorana_f_M_norm':float(abs(kappa)*np.linalg.norm(fd,2)),
        'selected_MR_GeV':cjson(mr),'selected_MR_Takagi_masses_GeV':np.linalg.svd(mr,compute_uv=False)[::-1].tolist(),
        'selected_Mnu_GeV':cjson(mnu),'selected_neutrino_masses_eV':nu_eV.tolist(),
        'selected_neutrino_sum_eV':float(nu_eV.sum()),
        'baseline_at_same_selected_overlap':{
            'MR_GeV':cjson(actual_mr),'MR_Takagi_masses_GeV':actual_mr_values.tolist(),
            'MN_over_comparison_scale':(actual_mr_values/comparison_scale).tolist(),
            'any_MR_above_comparison_scale':bool(np.any(actual_mr_values>comparison_scale)),
            'neutrino_masses_eV':(nu_eV*selected_sigma/sigma).tolist(),
            'Ynu_spectral_norm':ynu_norm,
            'Ynu_squared_log_MI_over_MN_over_16pi2_diagnostic':log_estimate.tolist(),
            'diagnostic_does_not_replace_matrix_RGE_or_thresholds':True,
            'sterile_threshold_EFT_assignment_completed':False,
            'same_overlap_chosen_for_representation_not_a_baseline_fit':True,
        },
        'baseline_max_dm31_eV2_over_this_capped_overlap_family':float(DM31*(interval[1]/sigma)**2) if feasible else None,
        'overlap_action_realization_demonstrated':False,
        'theta_is_a_permissive_LO_profile_variable_not_a_new_action_parameter':True,
    }

def fixed_sigma_residual(x,target,kappa,sigma):
    """One bounded refinement; actual acceptance still uses exact gates.

    Atmospheric central value fixes d. Remaining inequalities are
    (1+r^2)||H'||^2+(1+|rs|^2)d^2<=1 and ||F'||/d<=1.
    A 1% interior target avoids treating a soft-penalty boundary as a
    feasibility certificate. No observable discrepancy width is changed.
    """
    try:
        o=observables(x);h,f,r,s=(o[k] for k in ('h','f','r','s'))
        q=o['q'];factor=np.sqrt(DM31/(q[2]**2-q[0]**2))
        required_d=sigma*abs(kappa)*factor/(V174**2*1e9)
        lhs=(1+r*r)*np.linalg.norm(h,2)**2+(1+abs(r*s)**2)*required_d**2
        f_norm=np.linalg.norm(f,2)/required_d
        penalty=np.maximum(np.array([lhs,f_norm])-.99,0.)/.01
        return np.r_[residual(x,target),penalty]
    except (np.linalg.LinAlgError,ValueError,FloatingPointError):
        return np.full(19,1e6)


def export_fixed_candidate(fit,target,kappa,sigma,mu):
    o=observables(fit.x)
    masses=np.concatenate([o['masses'][key] for key in ('up','down','charged_lepton')])
    desired=np.concatenate([target['yukawa_singular_values_at_high'][key] for key in ('up','down','charged_lepton')])
    widths=np.array([.3,.1,.1,.3,.3,.1,.1,.1,.1]);mass_relative=masses/desired-1
    ckm_relative=o['ckm_angles']/target['CKM_s12_s23_s13']-1
    delta_error=abs(np.angle(np.exp(1j*(o['ckm_delta']-target['CKM_delta_rad']))))
    ratio_band=[DM21_BAND[0]/DM31_BAND[1],DM21_BAND[1]/DM31_BAND[0]]
    gates={'charged_masses_within_declared_LO_widths':bool(np.all(abs(mass_relative)<=widths)),
       'three_CKM_sines_within_10_percent':bool(np.all(abs(ckm_relative)<=.1)),
       'CKM_delta_within_0p12_rad':bool(delta_error<=.12),
       'three_PMNS_sin2_inside_frozen_NuFIT_3sigma_NO':bool(np.all(o['pmns_sin2']>=NU_BANDS[:,0]) and np.all(o['pmns_sin2']<=NU_BANDS[:,1])),
       'delta_m21_over_delta_m31_inside_3sigma_rectangle':bool(ratio_band[0]<=o['ratio']<=ratio_band[1])}
    family=normalization_family(o,kappa,sigma,comparison_scale=mu)
    names=['u','c','t','d','s','b','e','mu','tau'];worst=int(np.argmax(abs(mass_relative)/widths))
    return dict(label='one_fixed_sigma_refinement_from_best_shape',x=fit.x.tolist(),
        objective=float(np.sum(fit.fun[:17]**2)),penalized_objective=float(np.sum(fit.fun**2)),
        scaled_residuals=fit.fun[:17].tolist(),constraint_penalty_residuals=fit.fun[17:].tolist(),
        nfev=fit.nfev,njev=fit.njev,optimizer_success=fit.success,message=fit.message,
        fixed_sigma_GeV=sigma,observable_widths_unchanged=True,extra_start_or_retry_count=0,
        constraint_policy='one-percent interior target and 0.01 penalty scale for exact normalized-overlap/norm inequalities; final acceptance checked separately',
        Hprime=cjson(o['h']),Fprime=cjson(o['f']),r=float(o['r']),s=cjson(o['s']),
        Yukawa_matrices={k:cjson(v) for k,v in o['ys'].items()},
        charged_singular_values={k:v.tolist() for k,v in o['masses'].items()},
        charged_fractional_residuals=mass_relative.tolist(),CKM_moduli=abs(o['ckm']).tolist(),
        CKM_s12_s23_s13=o['ckm_angles'].tolist(),CKM_delta_rad=o['ckm_delta'],CKM_J=o['ckm_j'],
        CKM_fractional_residuals=ckm_relative.tolist(),CKM_delta_error_rad=float(delta_error),
        PMNS_moduli=abs(o['pmns']).tolist(),PMNS_sin2_theta12_theta23_theta13=o['pmns_sin2'].tolist(),
        PMNS_delta_rad_unfitted=o['pmns_delta'],PMNS_J=o['pmns_j'],mass_squared_ratio=float(o['ratio']),
        shape_screen_gates=gates,normalized_overlap_seesaw_family=family,
        joint_LO_flavor_and_baseline_sigma_screen_pass=bool(all(gates.values()) and family['baseline_sigma_inside_this_witness_family']),
        worst_charged_discrepancy=dict(species=names[worst],fractional_residual=float(mass_relative[worst]),
             declared_fractional_width=float(widths[worst]),width_units=float(abs(mass_relative[worst])/widths[worst])))


def run(mu=MI,sigma=SIGMA,max_nfev=1000):
    start=time.time();target=targets(mu)
    phasepath=RF/'output/p54_common_yukawa_phase.json'
    phases=json.loads(phasepath.read_text());kappa=complex(decode(phases['normalization']['f_M_over_f_D']))
    h0,f0=inputs.appendix_b_matrices()
    scale=target['yukawa_singular_values_at_high']['down'][-1]/.00530
    h0*=scale;f0*=scale
    r0=77.4189*target['yukawa_singular_values_at_high']['up'][-1]/(.437*scale)
    seeds=[('published_numbers_CP_conjugated_then_refitted',h0,f0.conj(),r0,.3140+.0282j),
           ('published_numbers_then_refitted',h0,f0,r0,.3140-.0282j),
           ('fixed_phase_perturbation_then_refitted',h0,f0.conj()*np.exp(.25j),r0,.30+.06j)]
    lower=np.r_[np.full(15,-100.),np.log(2.),-3.,-3.]
    upper=np.r_[np.full(15,100.),np.log(500.),3.,3.]
    fits=[]
    for label,h,f,r,s in seeds:
        initial=pack(h,f,r,s)
        fit=least_squares(residual,initial,args=(target,),bounds=(lower,upper),
                          max_nfev=max_nfev,xtol=1e-10,ftol=1e-10,gtol=1e-9,
                          x_scale='jac',diff_step=1e-5)
        fits.append({'label':label,'x':fit.x.tolist(),'initial_x':initial.tolist(),
                     'sum_squared_scaled_residuals':float(np.sum(fit.fun**2)),
                     'max_abs_scaled_residual':float(max(abs(fit.fun))),
                     'nfev':fit.nfev,'njev':fit.njev,'optimizer_success':fit.success,
                     'message':fit.message,'scaled_residuals':fit.fun.tolist()})
        print(label,'objective',fits[-1]['sum_squared_scaled_residuals'],'nfev',fit.nfev,flush=True)
    best=min(fits,key=lambda row:row['sum_squared_scaled_residuals']);obs=observables(np.array(best['x']))
    fixed_fit=least_squares(fixed_sigma_residual,np.asarray(best['x']),args=(target,kappa,sigma),
        bounds=(lower,upper),max_nfev=max_nfev,xtol=1e-10,ftol=1e-10,gtol=1e-9,
        x_scale='jac',diff_step=1e-5)
    fixed=export_fixed_candidate(fixed_fit,target,kappa,sigma,mu)
    print('one fixed-sigma refinement',fixed['objective'],fixed['penalized_objective'],
          'nfev',fixed['nfev'],'joint pass',fixed['joint_LO_flavor_and_baseline_sigma_screen_pass'],flush=True)
    masses=np.concatenate([obs['masses'][key] for key in ('up','down','charged_lepton')])
    desired=np.concatenate([target['yukawa_singular_values_at_high'][key] for key in ('up','down','charged_lepton')])
    widths=np.array([.3,.1,.1,.3,.3,.1,.1,.1,.1]);mass_relative=masses/desired-1
    ckm_relative=obs['ckm_angles']/target['CKM_s12_s23_s13']-1
    delta_error=abs(np.angle(np.exp(1j*(obs['ckm_delta']-target['CKM_delta_rad']))))
    nuok=bool(np.all(obs['pmns_sin2']>=NU_BANDS[:,0]) and np.all(obs['pmns_sin2']<=NU_BANDS[:,1]))
    ratio_band=[DM21_BAND[0]/DM31_BAND[1],DM21_BAND[1]/DM31_BAND[0]]
    gates={'charged_masses_within_declared_LO_widths':bool(np.all(abs(mass_relative)<=widths)),
           'three_CKM_sines_within_10_percent':bool(np.all(abs(ckm_relative)<=.1)),
           'CKM_delta_within_0p12_rad':bool(delta_error<=.12),
           'three_PMNS_sin2_inside_frozen_NuFIT_3sigma_NO':nuok,
           'delta_m21_over_delta_m31_inside_3sigma_rectangle':bool(ratio_band[0]<=obs['ratio']<=ratio_band[1])}
    family=normalization_family(obs,kappa,sigma,comparison_scale=mu)
    strict_family=normalization_family(obs,kappa,sigma,fcap=np.sqrt(4*np.pi)/(2*np.sqrt(6)),comparison_scale=mu)
    loose_family=normalization_family(obs,kappa,sigma,hcap=np.sqrt(4*np.pi),fcap=np.sqrt(4*np.pi),comparison_scale=mu)
    checks=[]
    def check(name,a,b,tol=2e-9):
        a,b=np.asarray(a),np.asarray(b);err=float(np.linalg.norm(a-b)/max(1e-30,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append(dict(name=name,residual=err,tolerance=tol,pass_=bool(err<tol)))
    check('common_Majorana_Dirac_magnitude',abs(kappa),2*np.sqrt(6))
    c={k:complex(decode(v)) for k,v in family['selected_overlaps_abde'].items()}
    hd=decode(family['selected_h_D']);fd=decode(family['selected_f_D'])
    check('four_overlaps_unit_normalized',sum(abs(v)**2 for v in c.values()),1.)
    for key,value in dict(up=c['b']*hd+c['e']*fd,down=c['a']*hd+c['d']*fd,
                          charged_lepton=c['a']*hd-3*c['d']*fd,neutrino_dirac=c['b']*hd-3*c['e']*fd).items():
        check('same_h_f_reconstruct_'+key,value,obs['ys'][key])
    check('MR_proportional_to_same_fD',decode(family['selected_MR_GeV']),kappa*family['selected_sigma_GeV']*fd)
    nm=np.asarray(family['selected_neutrino_masses_eV'])
    check('atmospheric_scale_fixed_by_single_sigma',nm[2]**2-nm[0]**2,DM31)
    _,un=seesaw.takagi(decode(family['selected_Mnu_GeV']))
    _,ue=diagonalize(obs['ys']['charged_lepton'])
    check('PMNS_moduli_SVD_vs_independent_Takagi',abs(ue.T@un),abs(obs['pmns']))
    # An algebraic family rotation check uses the fitted matrices, not toy data.
    rotation=np.linalg.qr(np.array([[1+.2j,.3,-.2j],[.1j,1-.3j,.2],[.4,.1j,1+.1j]]))[0]
    moved={key:rotation.T@value@rotation for key,value in obs['ys'].items()}
    mf=rotation.T@obs['f']@rotation
    ms=moved['neutrino_dirac']@np.linalg.solve(mf,moved['neutrino_dirac'].T)
    check('seesaw_family_covariance',ms,rotation.T@obs['shape']@rotation)
    _,mu_left=diagonalize(moved['up']);_,md_left=diagonalize(moved['down'])
    check('CKM_moduli_family_covariance',abs(mu_left.T@md_left.conj()),abs(obs['ckm']))
    ff=fixed['normalized_overlap_seesaw_family']
    fc={k:complex(decode(v)) for k,v in ff['selected_overlaps_abde'].items()}
    check('fixed_sigma_candidate_four_overlaps_normalized',sum(abs(v)**2 for v in fc.values()),1.)
    check('fixed_sigma_candidate_same_f_MR',decode(ff['selected_MR_GeV']),
          kappa*ff['selected_sigma_GeV']*decode(ff['selected_f_D']))
    fhd=decode(ff['selected_h_D']);ffd=decode(ff['selected_f_D'])
    for key,value in dict(up=fc['b']*fhd+fc['e']*ffd,down=fc['a']*fhd+fc['d']*ffd,
          charged_lepton=fc['a']*fhd-3*fc['d']*ffd,neutrino_dirac=fc['b']*fhd-3*fc['e']*ffd).items():
        check('fixed_sigma_candidate_same_h_f_'+key,value,decode(fixed['Yukawa_matrices'][key]))
    fnm=np.asarray(ff['selected_neutrino_masses_eV'])
    check('fixed_sigma_candidate_single_atmospheric_scale',fnm[2]**2-fnm[0]**2,DM31)
    status='LO_shape_witness_found' if all(gates.values()) else 'no_LO_shape_witness_in_this_bounded_search'
    local_sources=[Path(__file__),Path(inputs.__file__),Path(seesaw.__file__),phasepath]
    return dict(schema='p54-physics-first-flavor-lo-screen-v1',date='2026-09-22',status=status,
       scales=dict(comparison_mu_GeV=mu,baseline_MI_GeV=MI,baseline_g4I=G4I,baseline_sigma_GeV=sigma,
                   sigma_equals_MI_over_g4=True,v174_GeV=V174),
       hypotheses=['tree two-matrix sum rules imposed at MI as a permissive LO ansatz',
                   'common complex h_D,f_D and unit-normalized (a,b,d,e); MR=kappa*sigma*f_D',
                   'pure type I, normal ordering; type II set to zero as an approximation, not proved by scalar decoupling',
                   'charged targets transported with SM one-loop RG below MI; no free finite-threshold shifts',
                   'neutrino mixing and mass-ratio running neglected; no statistical global likelihood'],
       target_record=target,neutrino_targets=dict(sin2_theta12_theta23_theta13=ANGLES.tolist(),bands_3sigma=NU_BANDS.tolist(),
                  dm21_eV2=DM21,dm31_eV2=DM31,ratio_3sigma_rectangle=ratio_band),
       discrepancy_policy=dict(charged_mass_order=['u','c','t','d','s','b','e','mu','tau'],
           charged_fractional_widths=widths.tolist(),CKM_fractional_width=.10,CKM_delta_rad_width=.12,
           widths_are_not_experimental_errors=True,threshold_nuisance_parameters_added=False),
       search=dict(algorithm='bounded scipy least_squares, three predeclared literature-based shape starts plus one fixed-sigma refinement',
           max_nfev_per_start=max_nfev,parameter_lower=lower.tolist(),parameter_upper=upper.tolist(),
           no_global_minimum_or_exclusion_claim=True,starts=fits),
       best=dict(label=best['label'],objective=best['sum_squared_scaled_residuals'],
           Hprime=cjson(obs['h']),Fprime=cjson(obs['f']),r=float(obs['r']),s=cjson(obs['s']),
           Yukawa_matrices={k:cjson(v) for k,v in obs['ys'].items()},
           charged_singular_values={k:v.tolist() for k,v in obs['masses'].items()},
           charged_fractional_residuals=mass_relative.tolist(),CKM_moduli=abs(obs['ckm']).tolist(),
           CKM_s12_s23_s13=obs['ckm_angles'].tolist(),CKM_delta_rad=obs['ckm_delta'],CKM_J=obs['ckm_j'],
           CKM_fractional_residuals=ckm_relative.tolist(),CKM_delta_error_rad=float(delta_error),
           PMNS_moduli=abs(obs['pmns']).tolist(),PMNS_sin2_theta12_theta23_theta13=obs['pmns_sin2'].tolist(),
           PMNS_delta_rad_unfitted=obs['pmns_delta'],PMNS_J=obs['pmns_j'],mass_squared_ratio=float(obs['ratio']),
           shape_screen_gates=gates,normalized_overlap_seesaw_family=family,
           norm_cap_sensitivity_without_refitting={
               'Majorana_component_below_sqrt4pi':strict_family,
               'Dirac_units_below_sqrt4pi_not_full_perturbativity':loose_family},
           charged_tolerance_policy_sensitivity_without_refitting={
               'largest_light_quark_fractional_deviation':float(np.max(abs(mass_relative[[0,3,4]]))),
               'all_charged_pass_if_light_quark_width_35percent':bool(np.all(abs(mass_relative)<=np.array([.35,.1,.1,.35,.35,.1,.1,.1,.1]))),
               'official_declared_screen_not_changed':True}),
       fixed_sigma_refinement=fixed,
       joint_LO_flavor_and_scale_witness_found=fixed['joint_LO_flavor_and_baseline_sigma_screen_pass'],
       joint_flavor_seesaw_proton_matched_UV_witness_found=False,
       untested_observables=['leptonic Dirac CP phase and Majorana phases','absolute lightest-neutrino mass bounds and cosmology likelihood',
           'neutrinoless double beta decay likelihood','baryon asymmetry/leptogenesis','proton decay and gauge unification in this subtask',
           'precision electroweak/Higgs observables'],
       missing=['MU-to-MI Pati-Salam Yukawa running','neutrino-Yukawa feedback and sequential sterile thresholds',
           'finite upper/lower gauge/ghost/Yukawa matching, Wilson/box and C_HN seesaw matching',
           'actual scalar action realization of overlaps and suppression of type II','full perturbative control of collective Yukawa tensors'],
       scalar_benchmark_reused_for_physical_fit=False,full_physical_fit_completed=False,
       ultraviolet_model_viability_established=False,default_parameters_changed=False,
       necessary_falsification_policy='Only reject a specified LO witness/profile family under the stated inputs and caps; bounded optimizer failure is not a model no-go.',
       sources=SOURCES,source_sha256={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in local_sources},
       checks=checks,summary=dict(passed=sum(c['pass_'] for c in checks),total=len(checks),all_pass=all(c['pass_'] for c in checks)),
       runtime_seconds=time.time()-start)

def markdown(report):
    b=report['best'];f=b['normalized_overlap_seesaw_family'];g=b['shape_screen_gates']
    fixed=report['fixed_sigma_refinement'];ff=fixed['normalized_overlap_seesaw_family']
    lines=['# Physics-first two-matrix LO flavor screen','',report['status'],'',
      '**Conditional common-scale flavor screen, not a matched P54 physical fit.**','',
      f"Comparison scale: {report['scales']['comparison_mu_GeV']:.8g} GeV; sigma=MI/g4={report['scales']['baseline_sigma_GeV']:.8g} GeV.",
      'Nine charged masses and four CKM parameters are checked after one-loop SM target transport. Three PMNS angles and the solar/atmospheric ratio use frozen NuFIT 6.0 normal-ordering inputs without neutrino RG transport.',
      'The charged 10%/30% widths and CKM 10%/0.12-rad widths are declared LO discrepancy filters, not experimental uncertainties or adjustable threshold shifts.','',
      f"Bounded search objective: {b['objective']:.8g}; numerical checks {report['summary']['passed']}/{report['summary']['total']}.",'']
    lines += ['- '+key+': '+str(value) for key,value in g.items()]
    lines += ['',f"Normalized-overlap sigma interval (central atmospheric scale, ||h_D||,||f_D||<=1): {f['sigma_interval_GeV_central_dm31']} GeV.",
      f"Baseline sigma belongs to this witness family: {f['baseline_sigma_inside_this_witness_family']}.",
      f"Shape-minimum down Yukawa residual: {100*b['charged_fractional_residuals'][3]:.5g}%, compared with a predeclared 30% width. The formal gate is not relaxed after seeing this near miss.",
      f"Even removing norm caps leaves the normalized-overlap supremum sigma={f['sigma_supremum_from_overlap_normalization_without_any_norm_caps_GeV']:.8g} GeV for THIS shape texture.",
      'A failure excludes only this texture/overlap profile under these approximations, not all two-matrix models.','',
      '## One bounded fixed-sigma refinement','',
      f"Exactly one refinement from the best shape start used at most {report['search']['max_nfev_per_start']} evaluations. Flavor objective={fixed['objective']:.7g}; joint LO screen pass={fixed['joint_LO_flavor_and_baseline_sigma_screen_pass']}.",
      f"Its sigma interval is {ff['sigma_interval_GeV_central_dm31']} GeV; the baseline is inside, but the charged/CKM/neutrino gates must also pass.",
      f"Largest charged discrepancy: {fixed['worst_charged_discrepancy']['species']} {100*fixed['worst_charged_discrepancy']['fractional_residual']:+.5g}% versus the fixed {100*fixed['worst_charged_discrepancy']['declared_fractional_width']:.3g}% width.",
      f"At actual sigma, MN/MI={ff['baseline_at_same_selected_overlap']['MN_over_comparison_scale']}. Sterile threshold/EFT assignment remains uncomputed, including a state above MI if indicated.",
      'No extra starts, wider tolerances, or scalar retuning are used to turn the result green. Bounded failure is not a global no-go.','',
      '## Exact algebraic interface','',
      'The matrices are Hprime=a h_D and Fprime=d f_D. The common relations are Yd=Hprime+Fprime, Ye=Hprime-3Fprime, Yu=r(Hprime+s Fprime), and Ynu=r(Hprime-3s Fprime).',
      'Choose a=cos(theta)/sqrt(1+r^2), d=sin(theta)/sqrt(1+|r s|^2), b=r a, e=r s d. This gives unit norm exactly and MR=(i 2 sqrt(6)) sigma Fprime/d in the stored common phase convention.',
      'Writing H=||Hprime||_2, F=||Fprime||_2, the cap-compatible interval is F/fcap <= d <= sqrt([1-(1+r^2)H^2/hcap^2]/[1+|r s|^2]). It exists iff (1+r^2)H^2/hcap^2+(1+|r s|^2)F^2/fcap^2 <= 1.',
      'For Q=Ynu Fprime^(-1) Ynu^T with singular values q_i, matching the atmospheric central value gives sigma/d = 174^2 * 10^9 * sqrt(q_3^2-q_1^2)/(2 sqrt(6) sqrt(Delta m31^2)). This is the only neutrino scale adjustment.',
      'Caps on Dirac-unit matrices are screening conventions, not a perturbativity theorem. Raw h/f and Majorana f_M norms, a stricter Majorana-component cap, and a looser Dirac-unit cap sensitivity are all retained in JSON without refitting.',
      'The four normalized overlaps have NOT been obtained from a scalar stationary point or mass eigenvector.','',
      '## Still untested','']
    lines += ['- '+s for s in report['missing']+report['untested_observables']]
    lines += ['','## Primary sources','']+['- ['+s['use']+']('+s['url']+')' for s in SOURCES]
    return '\n'.join(lines)+'\n'

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mu',type=float,default=MI)
    parser.add_argument('--sigma',type=float,default=SIGMA);parser.add_argument('--max-nfev',type=int,default=1000)
    args=parser.parse_args();report=run(args.mu,args.sigma,args.max_nfev)
    OUT.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    OUT.with_suffix('.md').write_text(markdown(report));print(markdown(report))
    if not report['summary']['all_pass']:raise SystemExit(1)
