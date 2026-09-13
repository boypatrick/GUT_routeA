#!/usr/bin/env python3
"""Exact scalar determinant subtraction on the actual broken PS valley.

Computes the finite MSbar scalar determinant difference missed by replacing
the nonlocal PS Schur operator with its two-derivative truncation. This is
a constant-background scalar result, NOT the complete lower gauge/ghost
matching or its full source derivatives.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, eigh, expm

import verify_p54_upper_yukawa_thresholds as yuk
import verify_p54_upper_vector_matching as vec
import verify_p54_lower_covariant_eft as lower

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_lower_scalar_finite_subtraction"


def finite_trace(masses2,mu):
    values=np.asarray(masses2)
    if values.min() < -1e-8:
        raise ValueError("Finite stable scalar determinant requires nonnegative mass squares")
    hard=values[values>1e-8]
    return float(np.sum(hard**2*(np.log(hard/mu**2)-1.5))/(64*np.pi**2))


def calculate(h,x,ba,bh,reps,mu):
    valley=lower.scalar_valley_api(h,ba,bh)
    vector=lower.upper_vector_api(x,valley["T"],reps)
    raw=np.column_stack((ba,bh))
    full_metric=raw.T@vector["horizontal_projector"]@raw
    full_mass=raw.T@h@raw
    a,c,w=valley["A"],valley["C"],valley["W"]
    g0=ba.T@vector["horizontal_projector"]@ba
    g=vector["metric"]
    mf=eigh(full_mass,full_metric,eigvals_only=True)
    mc=np.linalg.eigvalsh(c)
    ml=eigh(valley["S"],g,eigvals_only=True)
    vf,vc,vl=[finite_trace(v,mu) for v in (mf,mc,ml)]
    return {"valley":valley,"vector":vector,"full_metric":full_metric,"full_mass":full_mass,
            "G0":g0,"m_full":mf,"m_C":mc,"m_local":ml,
            "V_full":vf,"V_C":vc,"V_PS_nonlocal":vf-vc,"V_PS_local":vl,
            "finite_local_truncation_difference":vf-vc-vl}


def run(cache_dir):
    geometry=yuk.build_upper_geometry(cache_dir)
    p1=geometry["ps"]["common_geometry"]["p1"]
    full=json.loads((RF/"output/p54_full_doublet_cw.json").read_text())
    pars,v=full["tree_parameters"],full["vacuum"]
    x=p1.vacuum_vector(v["omega"],v["sigma"],v["vs"])
    base=hashlib.sha256(Path(p1.__file__).read_bytes()+json.dumps(pars,sort_keys=True).encode()).digest()
    key=hashlib.sha256(base+np.asarray(x,dtype="<f8").tobytes()).hexdigest()
    path=Path(cache_dir)/(key+".npz")
    with np.load(path) as data:h=data["h"]
    ba,bh=geometry["active"],geometry["heavy"]
    reps=vec.build_vector_geometry(geometry)["R"]
    mu=full["scheme"]["mu_over_omega"]
    result=calculate(h,x,ba,bh,reps,mu)
    valley=result["valley"]
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        e=yuk.error(np.asarray(a),np.asarray(b))
        checks.append({"name":name,"residual":e,"tolerance":tol,"pass":e<tol})
    check("upper_heavy_basis_horizontal_at_broken_background",bh.T@result["vector"]["O"])
    check("full_kinetic_block_structure",result["full_metric"],block_diag(result["G0"],np.eye(55)))
    original=np.linalg.eigvalsh(h)
    check("quotient_full_290_positive_scalar_poles",result["m_full"][-290:],original[-290:])
    check("full_303_chart_has_13_zero_modes",np.sum(abs(result["m_full"])<1e-8),13)
    check("local_248_chart_has_13_zero_modes",np.sum(abs(result["m_local"])<1e-8),13)
    check("local_metric_includes_heavy_response",result["vector"]["metric"],
          result["G0"]+valley["W"].T@np.linalg.solve(valley["C"],np.linalg.solve(valley["C"],valley["W"])))
    rows=[]
    for z in (1e-4,.003,.03,.3,3.,30.):
        k=result["full_mass"]+z*result["full_metric"]
        d=valley["C"]+z*np.eye(55)
        s=lower.scalar_nonlocal_kernel_api(z,valley,result["vector"])
        sl=valley["S"]+z*result["vector"]["metric"]
        logs=[np.linalg.slogdet(a) for a in (k,d,s,sl)]
        check(f"p2={z}_positive_Euclidean_kernels",[a[0] for a in logs],np.ones(4))
        defect=logs[0][1]-logs[1][1]-logs[2][1]
        check(f"p2={z}_exact_nonlocal_logdet_subtraction",defect,tol=5e-8)
        rows.append({"p_E2":z,"exact_logdet_defect":float(defect),
                     "nonlocal_minus_local_logdet":float(logs[2][1]-logs[3][1])})
    # Independent normal-mode evaluation uses det(K)=det(G)*prod(z+m_i^2).
    z=.077
    spectral=np.linalg.slogdet(result["full_metric"])[1]+np.log(z+result["m_full"]).sum()
    check("full_generalized_spectrum_logdet",spectral,
          np.linalg.slogdet(result["full_mass"]+z*result["full_metric"])[1])
    # Scale derivative of the finite subtraction is fixed by the UV polynomial.
    step=1e-4
    def difference(scale):
        return finite_trace(result["m_full"],scale)-finite_trace(result["m_C"],scale)-finite_trace(result["m_local"],scale)
    derivative=(difference(mu*np.exp(step))-difference(mu*np.exp(-step)))/(2*step)
    delta_m4=np.sum(result["m_full"]**2)-np.sum(result["m_C"]**2)-np.sum(result["m_local"]**2)
    check("finite_subtraction_scale_derivative",derivative,-delta_m4/(32*np.pi**2))
    check("local_determinant_is_not_exact",float(abs(result["finite_local_truncation_difference"])<1e-9))
    # Broken PS rotations transform the background, not arbitrary eigenvectors.
    parent=yuk.module("lower_finite_PS",RF/"code/verify_p54_p2_two_site_matching.py")
    rg=[p1.representation_matrix(g) for group in parent.ps_generators(p1).values() for g in group]
    broken=max(rg,key=lambda r:np.linalg.norm(r@x))
    rotation=expm(.173*broken)
    rotated=calculate(rotation@h@rotation.T,rotation@x,ba,bh,reps,mu)
    check("finite_scalar_subtraction_broken_PS_covariance",
          rotated["finite_local_truncation_difference"],result["finite_local_truncation_difference"])
    # Wilson contact C remains needed even at the quadratic determinant stage.
    paths=[Path(__file__),Path(lower.__file__),Path(yuk.__file__),Path(vec.__file__),Path(p1.__file__),RF/"output/p54_full_doublet_cw.json"]
    return {"schema":"p54-lower-finite-scalar-subtraction-v1","date":"2026-09-11",
        "scheme":"same broken background, MSbar scalar c=3/2, exact quadratic PS Schur kernel; not full gauge/ghost functional",
        "vacuum":v,"mu":mu,"cache_key":key,
        "dimensions":{"full_quotient":303,"upper_heavy":55,"PS_local":248,"soft_zeros":13},
        "finite_values":{key:result[key] for key in ("V_full","V_C","V_PS_nonlocal","V_PS_local","finite_local_truncation_difference")},
        "UV_delta_trace_M4":float(delta_m4),"momentum_rows":rows,
        "complete_lower_gauge_covariant_matching":False,"source_derivatives_completed":False,"physical_fit_enabled":False,
        "missing":["background/source derivatives of the finite functional","transverse gauge and ghost/Goldstone package",
                   "upper loop Wilson insertion and lower EFT subtraction","finite kinetic/Yukawa/box operator matching"],
        "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# Lower finite scalar determinant subtraction", "",
        "**Actual finite constant-background scalar result; full lower covariant matching remains open.**", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "The 303-dimensional upper-vector quotient has 290 positive scalar poles and 13 retained zeros. Its exact quadratic Schur determinant factorizes into 55 upper scalar poles and the nonlocal 248-dimensional PS operator. The local two-derivative PS operator is not interchangeable with it inside a loop trace.", "",
        "```json",json.dumps(r["finite_values"],indent=2),"```", "",
        "Values are in omega^4 reference units. The nonlocal-minus-local difference is a finite matching contribution in the declared scalar determinant, not a particle mass prediction. A constant vacuum term alone does not provide its source derivatives or the required Wilson coefficients.", "",
        "| p_E^2 | Exact logdet defect | Nonlocal minus local logdet |", "|---|---:|---:|"]
    lines += [f"| {v['p_E2']} | {v['exact_logdet_defect']:.3e} | {v['nonlocal_minus_local_logdet']:.8g} |" for v in r["momentum_rows"]]
    lines += ["", "## Missing components", ""]+["- "+s for s in r["missing"]]
    lines += ["", "## Checks", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    p=argparse.ArgumentParser();p.add_argument("--cache-dir",type=Path,required=True)
    args=p.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"finite_values":r["finite_values"],
                      "failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
