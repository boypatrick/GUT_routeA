#!/usr/bin/env python3
"""Actual upper PS scalar kinetic/mass hard kernels and common-CT tadpoles.

The action is unchanged. New directional Hessian evaluations are polynomial
vertex evaluations, NOT vacuum/lattice scans. Leading retained-mass hard
expansion, exact heavy masses. Full upper Wilson matching is not asserted.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm

import verify_p54_upper_yukawa_thresholds as yuk
import verify_p54_scalar_kinetic_matching as scalar
import verify_p54_upper_vector_matching as vector

RF = Path(__file__).resolve().parents[1]
ROOT = RF.parent
OUT = RF/"output/p54_upper_scalar_matching"
LOOP = 16*np.pi**2


class ActionJets:
    def __init__(self, p1, parameters, read_cache):
        self.p1, self.pars = p1, parameters
        self.read_cache = Path(read_cache)
        self.write_cache = ROOT/"tmp/p54_upper_scalar_jets"
        self.keybase = hashlib.sha256(Path(p1.__file__).read_bytes()+json.dumps(parameters,sort_keys=True).encode()).digest()
        self.hf = None
        self.keys, self.hits, self.evaluated = [],0,0

    def hessian(self,x,label):
        key = hashlib.sha256(self.keybase+np.asarray(x,dtype="<f8").tobytes()).hexdigest()
        self.keys.append(key)
        for directory in (self.read_cache,self.write_cache):
            path = directory/(key+".npz")
            if path.exists():
                self.hits += 1
                with np.load(path) as data:
                    return data["h"]
        tick = time.time()
        if self.hf is None:
            self.hf = self.p1.jax.jit(self.p1.hessian(self.p1.potential_factory(self.pars)))
        result = np.asarray(self.hf(self.p1.anp.asarray(x)),float)
        self.write_cache.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(self.write_cache/(key+".npz"),h=result)
        self.evaluated += 1
        print(f"{label}: new polynomial Hessian jet, {time.time()-tick:.1f}s",flush=True)
        return result

    def jets(self,x,h0,q,label,step=.02):
        hp,hm = self.hessian(x+step*q,label+"+"),self.hessian(x-step*q,label+"-")
        return (hp-hm)/(2*step),(hp+hm-2*h0)/step**2


def invariant_basis(geometry):
    """Eight real PS-invariant symmetric bilinears: H(3),F(3),L(1),R(1).

    Real-type complex bidoublets have two equivalent real copies. Mirror
    tensors construct their actual real structure; no fitted projector.
    """
    mirror = yuk.build_mirror_basis(geometry)
    flow, ps = geometry["flow"],geometry["ps"]
    mats,names = [],[]
    offset = 0
    for key in flow.KEYS:
        size = 2*len(ps["tensors"][key])
        p = np.zeros((248,248)); p[offset:offset+size,offset:offset+size]=np.eye(size)
        mats.append(p/np.linalg.norm(p)); names.append(key+"_norm")
        if key in mirror:
            couplings = {k:np.zeros((1,1),complex) for k in flow.KEYS}
            couplings[key][:] = 1
            tree = flow.assemble_real_yukawas(couplings,ps)[offset:offset+size]
            flat = np.column_stack((tree.real.reshape(size,-1),tree.imag.reshape(size,-1)))
            for phase,label in ((1.,"mirror_re"),(1j,"mirror_im")):
                mc = {k:np.array([[phase if k == key else 0]],complex) for k in mirror}
                yy = yuk.assemble_mirror(mc,geometry,mirror)[offset:offset+size]
                target = np.column_stack((yy.real.reshape(size,-1),yy.imag.reshape(size,-1)))
                s = target@np.linalg.pinv(flat)
                if yuk.error(s,s.T)>1e-10 or yuk.error(s@flat,target)>1e-10:
                    raise ValueError("Mirror field real structure failed")
                full = np.zeros((248,248)); full[offset:offset+size,offset:offset+size]=s
                for old in mats:
                    full -= np.sum(full*old)*old
                full /= np.linalg.norm(full)
                mats.append(full); names.append(key+"_"+label)
        offset += size
    return np.array(mats),names


def directions_from_basis(mats):
    """Greedy exact tensor-basis probes; includes complex-phase directions."""
    candidates=[]
    for mat in mats:
        values,u=np.linalg.eigh(mat)
        candidates.extend([u[:,0],u[:,-1]])
    selected,rows=[],[]
    while len(rows)<len(mats):
        best=None
        for i,q in enumerate(candidates):
            row=np.einsum("i,aij,j->a",q,mats,q)
            residual=row.copy()
            if rows:
                old=np.asarray(rows)
                residual-=row@np.linalg.pinv(old)@old
            score=np.linalg.norm(residual)
            if best is None or score>best[0]:
                best=(score,i,row)
        if best[0]<1e-10:
            raise ValueError("Insufficient invariant probes")
        selected.append(candidates.pop(best[1]));rows.append(best[2])
    return np.array(selected),np.array(rows)


def loop_weights(geometry,mu):
    heavy=geometry["heavy"]
    # Keep all upper Goldstone/PQ directions in the Landau scalar sum.
    d,u=np.linalg.eigh(np.eye(328)-heavy@heavy.T)
    soft=u[:,d>.5]
    basis=np.column_stack((heavy,soft))
    m=np.r_[geometry["masses2"],np.zeros(len(soft.T))]
    hard=np.arange(328)<55
    include=hard[:,None]|hard[None,:]
    xx,yy=np.broadcast_arrays(m[:,None],m[None,:])
    slopes=np.zeros((328,328));b0=np.zeros_like(slopes)
    slopes[include]=scalar.bubble_slope(xx[include],yy[include])
    for i in range(328):
        for j in range(i,328):
            if include[i,j]:
                x,y=max(m[i],m[j]),min(m[i],m[j])
                if y == 0:
                    b=1-math.log(x/mu**2)
                elif abs(x-y)<1e-8*max(x,y):
                    b=-math.log((x+y)/(2*mu**2))
                else:
                    b=1-(x*math.log(x/mu**2)-y*math.log(y/mu**2))/(x-y)
                b0[i,j]=b0[j,i]=b
    a0=np.zeros(328)
    a0[hard]=m[hard]*(np.log(m[hard]/mu**2)-1)
    return basis,m,slopes,b0,a0


def run(read_cache):
    geometry=yuk.build_upper_geometry(read_cache)
    p1=geometry["ps"]["common_geometry"]["p1"]
    pars=geometry["upper"]["same_action_parameters"]
    store=ActionJets(p1,pars,read_cache)
    mats,names=invariant_basis(geometry)
    directions,measurement=directions_from_basis(mats)
    vgeom=vector.build_vector_geometry(geometry)
    full=json.loads((RF/"output/p54_full_doublet_cw.json").read_text())
    # Evaluate all terms at the fixed-VEV anchor scale, not an almost-equal
    # scale with silently frozen finite counterterms.
    mu=float(full["scheme"]["mu_over_omega"])
    basis,m,slopes,bubbles,a0=loop_weights(geometry,mu)
    x,h0,active=geometry["x"],geometry["hessian"],geometry["active"]
    checks=[]
    def check(name,a,b=0.,tol=2e-8):
        e=yuk.error(np.asarray(a),np.asarray(b))
        checks.append({"name":name,"residual":e,"tolerance":tol,"pass":e<tol})
    check("eight_PS_invariant_real_symmetric_bilinears",np.einsum("aij,bij->ab",mats,mats),np.eye(8))
    check("full_heavy_soft_state_sum",basis@basis.T,np.eye(328))
    check("probe_matrix_invertible",float(np.linalg.matrix_rank(measurement)!=8))

    def contractions(t,u):
        te=basis.T@t@basis
        ue=basis.T@u@basis
        kinetic=np.sum(te*te*slopes)/(2*LOOP)
        seagull=np.dot(a0,np.diag(ue))/(2*LOOP)
        bubble=-np.sum(te*te*bubbles)/(2*LOOP)
        return np.array([kinetic,seagull,bubble])

    values=[]
    for i,q in enumerate(directions):
        t,u=store.jets(x,h0,active@q,"upper_active_"+str(i))
        values.append(contractions(t,u))
    coefficients=np.linalg.solve(measurement,np.array(values))
    ks,ms_seagull,ms_bubble=[np.einsum("a,aij->ij",coefficients[:,i],mats) for i in range(3)]
    mass=ms_seagull+ms_bubble
    check("scalar_kinetic_positive_Gram",max(0.,-np.linalg.eigvalsh(ks).min()))
    check("matched_scalar_kinetic_metric_positive",float(np.linalg.eigvalsh(np.eye(248)+ks).min()<=0))
    # Independent mixed-parent and complex-phase probes, not reconstruction inputs.
    rng=np.random.default_rng(20260911)
    for i in range(2):
        q=rng.normal(size=248); q/=np.linalg.norm(q)
        t,u=store.jets(x,h0,active@q,"independent_probe_"+str(i))
        got=contractions(t,u)
        check("independent_kinetic_mass_invariant_reconstruction_"+str(i),got,
              [q@ks@q,q@ms_seagull@q,q@ms_bubble@q],tol=5e-7)
    # Exact quartic polynomial: central Hessian differences are exact jets,
    # not an extrapolated fit to a tabulated mass formula.
    q=directions[0]
    t,u=store.jets(x,h0,active@q,"halfstep_probe",step=.01)
    check("quartic_action_two_step_jet_identity",contractions(t,u),values[0],tol=5e-7)
    parents=yuk.module("upper_scalar_PS_generators",RF/"code/verify_p54_p2_two_site_matching.py")
    ward=0.
    for group in parents.ps_generators(p1).values():
        for gen in group:
            r=active.T@p1.representation_matrix(gen)@active
            ward=max(ward,yuk.error(r@ks,ks@r),yuk.error(r@mass,mass@r))
    check("full_PS_Ward_for_scalar_kinetic_and_mass",ward)

    # Common UV counterterms are anchored at the historical broken vacuum.
    # They cannot be fixed independently at the upper stationary point.
    ctrow=full["tadpoles"]["finite_mass_counterterms_over_omega2"]
    ct=np.zeros((328,328))
    ct[:54,:54]=-ctrow["delta_mu2"]*np.eye(54)
    ct[54:306,54:306]=-ctrow["delta_nu2"]*.5*np.eye(252)
    ct[326:,326:]=-ctrow["delta_mus2"]*np.eye(2)
    radial=np.column_stack((p1.vacuum_vector(1,0,0),p1.vacuum_vector(0,0,1)))
    radial_jets=[]; tad_s=[];tad_v=[]
    rs=vgeom["R"]
    rx=np.einsum("Aij,j->Ai",rs,x)
    mv=vgeom["mass2"]
    # Use one common scale for these upper scalar/vector hard contributions.
    for i,q in enumerate(radial.T):
        t,u=store.jets(x,h0,q,"upper_radial_"+str(i))
        radial_jets.append(t)
        tad_s.append(float(np.dot(a0,np.diag(basis.T@t@basis))/(2*LOOP)))
        rq=np.einsum("Aij,j->Ai",rs,q)
        m1=vgeom["g"]**2*(rq@rx.T+rx@rq.T)
        tad_v.append(float(3*mv*(math.log(mv/mu**2)-1/3)*np.trace(m1)/(2*LOOP)))
    tct=radial.T@ct@x
    hard_t=np.array(tad_s)+np.array(tad_v)
    # Broken-anchored CT is quoted at mu0. Its one-loop scale dependence
    # must be included rather than silently reading it at muU.
    mu0=float(full["scheme"]["mu_over_omega"])
    if abs(mu/mu0-1)>1e-14:
        # Anchor the tadpole diagnostic at exactly the CT reference scale.
        _,_,_,_,anchor_a0=loop_weights(geometry,mu0)
        for i,t in enumerate(radial_jets):
            tad_s[i]=float(np.dot(anchor_a0,np.diag(basis.T@t@basis))/(2*LOOP))
            q=radial[:,i];rq=np.einsum("Aij,j->Ai",rs,q)
            m1=vgeom["g"]**2*(rq@rx.T+rx@rq.T)
            tad_v[i]=float(3*mv*(math.log(mv/mu0**2)-1/3)*np.trace(m1)/(2*LOOP))
        hard_t=np.array(tad_s)+np.array(tad_v)
    residual=hard_t+tct
    hr=radial.T@h0@radial
    displacement=-np.linalg.solve(hr,residual)
    check("upper_hard_tadpole_linear_shift_equation",hr@displacement+residual)
    # This is the heavy contribution plus common counterterm, NOT the
    # full upper quantum saddle: retained EFT loops are not included.
    delta_mass_shift=sum(a*(active.T@t@active) for a,t in zip(displacement,radial_jets))
    # Active VV-scalar linear vertices vanish at this PS background. The
    # vector CW curvature is therefore its quadratic mass seagull alone.
    vector_mass=3*vgeom["g"]**2*mv*(math.log(mv/mu**2)-1/3)*vgeom["CS"]/LOOP
    proper_mass=mass+vector_mass+active.T@ct@active+delta_mass_shift
    _,kv=vector.kinetic(geometry,vgeom,mu)
    ksum=ks+kv
    tree_mass=active.T@h0@active
    canonical_mass=proper_mass-.5*(ksum@tree_mass+tree_mass@ksum)
    check("combined_bosonic_kinetic_metric_positive",float(np.linalg.eigvalsh(np.eye(248)+ksum).min()<=0))
    check("canonical_mass_is_real_symmetric",canonical_mass,canonical_mass.T)
    phases=yuk.build_mirror_basis(geometry)
    scalar_y={}
    for s in ("h","f"):
        y0=np.einsum("ab,aij->bij",active,geometry["all_"+s])
        dy=-.5*np.einsum("ab,bij->aij",ks,y0)
        proj=yuk.project_declared(dy,geometry)
        mir=yuk.project_mirror(proj["residual"],geometry,phases)
        check(s+"_scalar_leg_six_invariant_closure",dy,
              proj["reconstruction"]+yuk.assemble_mirror(mir,geometry,phases))
        scalar_y[s]={"direct":{k:yuk.cjson(v) for k,v in proj["couplings"].items()},
                     "mirror":{k:yuk.cjson(v) for k,v in mir.items()}}
    source_paths=[Path(__file__),Path(yuk.__file__),Path(scalar.__file__),Path(vector.__file__),
                  Path(p1.__file__),RF/"output/p54_full_doublet_cw.json",RF/"output/p54_ps_finite_thresholds.json"]
    report={"schema":"p54-upper-scalar-hard-subset-v1","date":"2026-09-11",
        "scope":"actual upper scalar-cubic kinetic and hard seagull/bubble mass; common-CT hard tadpole diagnostic",
        "scheme":{"action":"P54PQ-v2","gauge":"background-field Landau","subtraction":"MSbar-DR",
                  "tadpole_anchor":"historical broken fixed-VEV; never independently retuned at upper PS"},
        "vacuum":geometry["upper"]["vacuum"],"mass_kinetic_matching_mu":mu,
        "invariant_names":names,"invariant_coefficients_kinetic_seagull_bubble":coefficients.tolist(),
        "Kscalar":ks.tolist(),"hard_scalar_mass":mass.tolist(),
        "Kscalar_eigenvalue_range":np.linalg.eigvalsh(ks)[[0,-1]].tolist(),
        "hard_scalar_mass_eigenvalue_range":np.linalg.eigvalsh(mass)[[0,-1]].tolist(),
        "common_anchor_two_point":{"mu":mu,"vector_CW_curvature":vector_mass.tolist(),
            "scalar_plus_vector_kinetic":ksum.tolist(),
            "hard_plus_CT_plus_valley_mass":proper_mass.tolist(),
            "canonical_mass_increment":canonical_mass.tolist(),
            "additional_light_doublet_retuning_included":False,
            "fermionic_fixed_VEV_counterterms_included":False,
            "full_physical_stationarity_claim":False},
        "scalar_leg_Yukawa_coefficients":scalar_y,
        "common_CT_upper_tadpole":{"mu":mu0,"scalar_hard":tad_s,"vector_hard":tad_v,
            "counterterm":tct.tolist(),"residual":residual.tolist(),
            "radial_hessian":hr.tolist(),"hard_plus_CT_linear_displacement":displacement.tolist(),
            "induced_active_mass_shift":delta_mass_shift.tolist(),
            "full_upper_quantum_stationary_point":False,
            "reason":"retained EFT loops and fermionic fixed-VEV counterterms remain absent"},
        "cache":{"hits":store.hits,"new_vertex_evaluations":store.evaluated,"keys":store.keys},
        "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "full_upper_Wilson_matching":False,"full_lower_matching":False,"physical_fit_enabled":False,
        "missing":["source-dependent three/four-point upper Wilson loops","retained-mass expansion control",
                   "upper EFT soft-loop tadpoles in common global prescription","lower gauge-covariant finite matching",
                   "finite matching with nonzero CHN and partially retained sterile fields"],
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}
    return report


def markdown(r):
    lines=["# Upper scalar hard matching and common tadpole anchor", "",
        "**Actual scalar kinetic/mass subset, not full upper Wilson or lower finite matching.**", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "The eight real PS-invariant scalar bilinears are derived from the existing intertwiners. Eight polynomial-action probes reconstruct the full 248x248 metric and mass correction; independent complex/mixed-parent directions test the reconstruction.", "",
        f"Scalar-cubic kinetic eigenvalue range: {r['Kscalar_eigenvalue_range']}. "
        f"Hard scalar mass eigenvalue range: {r['hard_scalar_mass_eigenvalue_range']}.", "",
        "Heavy-heavy and mixed-heavy-soft propagators are included at leading retained-mass hard order. Pure-soft graphs are not included in this hard coefficient. The gauge-fixed upper Goldstone directions are kept in the scalar state sum, not mistaken for physical retained PS scalars.", "",
        "## Common fixed-VEV prescription", "",
        "The pre-existing invariant mass counterterms are anchored at the broken vacuum and cannot be independently adjusted at the upper PS saddle. The following diagnostic uses the exact counterterm reference scale, not the nearby upper matching scale:", "",
        "```json",json.dumps({k:v for k,v in r['common_CT_upper_tadpole'].items() if k!='induced_active_mass_shift'},indent=2),"```", "",
        "This displacement is only the hard-plus-counterterm contribution. Retained EFT-loop tadpoles are still necessary for a full upper quantum stationary background. It must not be substituted as a completed background into a fit.", "",
        "## Missing items", ""]+["- "+s for s in r["missing"]]
    lines += ["", "## Checks", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--cache-dir",type=Path,required=True)
    args=ap.parse_args();r=run(args.cache_dir)
    OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"K_range":r["Kscalar_eigenvalue_range"],
         "tadpole_residual":r["common_CT_upper_tadpole"]["residual"],
         "displacement":r["common_CT_upper_tadpole"]["hard_plus_CT_linear_displacement"],
         "failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]:raise SystemExit(1)
