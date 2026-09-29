#!/usr/bin/env python3
"""G-R5 observation model and data-eligibility checks. No experimental fit.

Fixtures test algebra only; they are not surrogate clock observations.
The full release audit is separately reproducible by audit_gr5_clock_data.py.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/gr5_external_reference'


def run():
    cfg=json.loads((ROOT/'data/BBR_CLOCK_INPUTS.json').read_text())
    card=json.loads((ROOT/'data/EXTERNAL_REFERENCE_CARD.json').read_text())
    audit=json.loads((ROOT/'output/gr5_clock_data_full_audit.json').read_text())
    k=cfg['constants']
    h,kb,c,eps=[k[s] for s in ['h_J_s','kB_J_K','c_m_s','epsilon0_F_m']]
    arad=8*np.pi**5*kb**4/(15*h**3*c**3)
    u0=arad*300**4
    nu=np.array([v['frequency_Hz'] for v in cfg['transitions']])
    alpha=np.array([v['delta_alpha_SI'] for v in cfg['transitions']])
    sigma=np.array([v['published_sigma_alpha_SI'] for v in cfg['transitions']])
    s=-alpha*u0/(2*h*eps*nu)
    sig_s=sigma*u0/(2*h*eps*nu)
    # Units rescaled for matrix-rank diagnostics, not extra physical input.
    v=s/1e-16
    checks=[]
    def check(name, ok):
        checks.append(dict(name=name,passed=bool(ok)))
    def close(name,a,b,tol=1e-12):
        check(name,np.allclose(a,b,rtol=tol,atol=tol))
    rho=s[0]/s[1]
    weights=np.array([1,-rho])/(1-rho)
    close('composite:retains_common_factor',weights.sum(),1)
    close('composite:cancels_leading_target_BBR',weights@v,0)
    check('composite:requires_negative_weight',weights[0]<0<weights[1])
    d=np.array([-1.,1.])
    close('differential:common_factor_null',d@np.ones(2),0)
    M=np.column_stack([np.ones(2),v])
    check('design:external_two_channels_rank_two',np.linalg.matrix_rank(M)==2)
    check('design:internal_ratio_has_no_common_column',np.linalg.matrix_rank((d@M)[None,:])==1)
    N=np.column_stack([M,np.ones(2)])
    check('design:unknown_common_bias_not_identifiable',np.linalg.matrix_rank(N)==2)
    close('design:common_bias_exact_null_vector',N@np.array([1.,0.,-1.]),np.zeros(2))
    # Ideal fixtures are intentionally far above rounding error in scaled units.
    for chi,x in [(0.,.14),(.3,0.),(-.2,.09)]:
        y=M@np.array([chi,x])
        close(f'fixture:{chi}:{x}:recover_two_parameters',np.linalg.solve(M,y),[chi,x])
        close(f'fixture:{chi}:{x}:common_composite',weights@y,chi)
        close(f'fixture:{chi}:{x}:differential',d@y,(v[1]-v[0])*x)
    # Calibration of radiation is an independent input, not extracted from y
    # and then treated as an independent validation of the same response.
    radiometry_design=np.vstack([M,[0,1]])
    check('design:independent_radiometry_adds_overconstraint',
          radiometry_design.shape[0]-np.linalg.matrix_rank(radiometry_design)==1)
    x=(310/300)**4-1
    source_delta=s*x
    log_delta=np.log1p(s*(310/300)**4)-np.log1p(s)
    close('static:linear_vs_log_at_retained_order',log_delta/1e-16,source_delta/1e-16)
    # Exact log observation model (within the declared multiplicative model).
    target_b=np.array([-.07,-.01])
    ref_b=-.03
    chiS,chiR,g,link=.04,.015,.002,.003
    logq=chiS-chiR+np.log1p(target_b)-np.log1p(ref_b)+g+link
    residual=logq-np.log1p(target_b)+np.log1p(ref_b)-g-link
    close('observation:reference_response_retained',residual,np.full(2,chiS-chiR))
    close('observation:ratio_cancels_reference',d@logq,np.log1p(target_b[1])-np.log1p(target_b[0]))
    close('observation:globally_common_factor_cancels',
          (chiS+.2)-(chiR+.2),chiS-chiR)
    # A constructed triangle supplies no third independent physical datum.
    triangle=np.vstack([np.eye(2),d])
    check('observation:derived_triangle_rank_two',np.linalg.matrix_rank(triangle)==2)
    close('observation:derived_triangle_closure',np.array([1,-1,1])@triangle,np.zeros(2))
    # ABBA cancels offset and linear drift, not quadratic/heater-locked bias.
    t=np.array([-3.,-1.,1.,3.])
    w=np.array([-.5,.5,.5,-.5])
    state=np.array([0.,1.,1.,0.])
    close('ABBA:offset_cancels',w@np.ones(4),0)
    close('ABBA:linear_drift_cancels',w@t,0)
    close('ABBA:unit_state_contrast',w@state,1)
    close('ABBA:quadratic_drift_survives',w@t**2,-8)
    # Same reference error cancels only with identical effective windows.
    W2=np.array([.5,.5,0.])
    W3=np.array([0.,.5,.5])
    drift=np.array([0.,1.,2.])
    close('sampling:matched_reference_cancels',(W2-W2)@drift,0)
    close('sampling:interleaved_reference_leakage',(W3-W2)@drift,1)
    # Shared-reference covariance must not be reduced as independent noise.
    common_sigma=3.
    C=np.diag([4.,1.])+common_sigma**2*np.ones((2,2))
    close('covariance:difference_removes_common',d@C@d,5)
    close('covariance:composite_keeps_common',
          weights@C@weights,4*weights[0]**2+weights[1]**2+common_sigma**2)
    inv=np.linalg.inv(C)
    wgls=inv@np.ones(2)/(np.ones(2)@inv@np.ones(2))
    close('covariance:GLS_normalization',wgls.sum(),1)
    check('covariance:GLS_common_floor',1/(np.ones(2)@inv@np.ones(2))>=common_sigma**2)
    # Illustrative control budgets, not an achieved experimental uncertainty.
    ref_leak=s[1]*((300.1/300)**4-1)
    gravity_1cm=9.80665*.01/c**2
    coefficient_errors=abs(weights)*sig_s*x
    coefficient_bounds=[abs(coefficient_errors[0]-coefficient_errors[1]),sum(coefficient_errors)]
    check('budget:reference_thermal_leak_is_not_zero',ref_leak<0)
    check('budget:centimeter_gravity_is_relevant',gravity_1cm>1e-18)
    check('budget:old_polarizabilities_not_sub_1e18_common_calibration',coefficient_bounds[0]>1e-18)
    # Re-read packaged primary bytes, not reconstructed observations.
    archive=(ROOT/'data/clock_comparison_audit/ptb_lpi_2021.zip').read_bytes()
    check('data:PTB_archive_hash',hashlib.sha256(archive).hexdigest()==audit['ptb_2021']['archive_sha256'])
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        a=np.loadtxt(io.StringIO(z.read('Data/E3_E2_frequency_ratio_measurement.txt').decode()))
    check('data:PTB_eleven_published_points',a.shape==(11,5))
    check('data:PTB_absolute_not_fractional_offset',not audit['ptb_2021']['ratio_is_fractional'])
    for name, p in audit['rocit_2025']['files'].items():
        local=ROOT/'data/clock_comparison_audit/rocit_2025'/Path(name).name
        if local.exists():
            check('data:bundled_ROCIT_day_hash',hashlib.sha256(local.read_bytes()).hexdigest()==p['sha256'])
    check('data:full_audit_34_files',audit['rocit_2025']['selected_files']==34)
    check('data:full_audit_has_2324201_rows',audit['rocit_2025']['rows']==2324201)
    check('data:thermal_test_remains_open',not audit['empirical_differential_test_ready'])
    check('data:no_common_time_measurement',not audit['common_time_effect_measured'])
    check('card:no_fitted_time_law',card['no_new_time_law'])
    failed=[r['name'] for r in checks if not r['passed']]
    result=dict(stage='G-R5',status='reference-model-bounded-done; empirical-thermal-test-open',
                measurement_claim=False, fitted_parameters=False,
                composite_weights_E2_E3=weights.tolist(),
                static_target_300_to_310_fractional_shifts=source_delta.tolist(),
                reference_300_to_300p1_fractional_shift=ref_leak,
                gravity_1cm_fractional_shift=gravity_1cm,
                composite_coefficient_sigma_bounds_unknown_correlation=coefficient_bounds,
                dataset_audit='gr5_clock_data_full_audit.json',
                checks=checks,summary=dict(checks=len(checks),passed=len(checks)-len(failed),failed=failed))
    OUT.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    report=f'''# G-R5: data acquisition and external-reference observation model

Status: {result['status']}. {len(checks)-len(failed)}/{len(checks)} checks pass.
No new thermal measurement, no fitted universal-time parameter.

## Acquired evidence and empirical gate

The PTB 2021 archive contains 11 published ratio points. The full ROCIT
2025 release contains 34 selected same-package E2/E3 files with 2,324,201
rows, all marked valid, zero nonfinite values and zero duplicate timestamps.
The selected fields do not include paired bath-state labels, independent
radiometry or a per-run BBR correction ledger. Data integrity is not
eligibility for a thermal-response test. No such fit is run.
Source attribution, version/hash and original-vs-processed distinctions:
data/clock_comparison_audit/PROVENANCE.md.

## Selected reference and physical model

Use a second independent Yb+ E3 ion clock outside the modulated target
enclosure, connected by a comb and stabilized optical paths. Its own
radiation response, drift, common-factor response and path/gravity terms
remain in the model. The measurable candidate is Delta log(F_S/F_R),
not an absolute F_S. Two target-to-reference channels add common-mode
information; their internally derived ratio is not a third independent datum.

At leading static order a composite with E2/E3 weights
({weights[0]:.6f}, {weights[1]:.6f}) cancels the target's T^4 response
and retains common target-reference response with coefficient one.
It still retains reference drift and common path/gravity biases.
An unconstrained such bias is exactly degenerate with the desired signal.

For the illustrative 300 -> 310 K target change, uncertainties of the
pinned polarizabilities alone give a composite coefficient standard-error
range [{coefficient_bounds[0]:.4g}, {coefficient_bounds[1]:.4g}] over unknown
correlation. This is not an achieved total uncertainty or a confidence
interval. These old inputs do not support a sub-1e-18 common-mode claim.
A 0.1 K reference-bath change contributes {ref_leak:.4g}; a 1 cm relative
height change gives about {gravity_1cm:.4g} in the weak-field model.

## Concrete next step

Obtain the missing clock/thermal/correction records described in
data/EXTERNAL_REFERENCE_CARD.json. Use matched effective time weights
and low-high-high-low blocks to suppress linear drift; nonlinear and
heater-correlated bias do not cancel automatically. Keep radiometry
independent of the frequency law under test. No lab contact has been sent.
Without those records, the empirical test stays open; do not replace it
with another fit, dataset-wide Gaussian scan, or a universal-time claim.

All numerical fixtures test algebra, not fabricated observations.
Full derivation: tex/route_g_external_reference.tex.
'''
    OUT.with_suffix('.md').write_text(report)
    print(json.dumps(result['summary']))
    print('composite weights:',weights,'coefficient-error bounds:',coefficient_bounds)
    if failed:
        raise SystemExit(1)


if __name__=='__main__':
    run()
