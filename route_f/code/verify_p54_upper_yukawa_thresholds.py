#!/usr/bin/env python3
"""Actual upper-PS pure-Yukawa finite thresholds, not full finite matching.

The scalar saddle/masses and Spin(10) intertwiners are the existing P54 action.
Only the physical upper heavy scalar block is integrated. All fermions are
massless at sigma=0. Family matrices in the executable test are synthetic.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
OUT = RF / "output/p54_upper_yukawa_thresholds"
UPPER = RF / "output/p54_ps_finite_thresholds.json"
FLOW = RF / "code/verify_p54_ps_yukawa_flow.py"
COMMON = RF / "code/verify_p54_common_yukawa_phase.py"
P1 = RF / "code/verify_p54_p1_hessian_spectrum.py"
LOOP = 16 * math.pi**2


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def cjson(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def error(x, y):
    return float(np.linalg.norm(x-y)/max(1., np.linalg.norm(x), np.linalg.norm(y)))


def image(x, tol=1e-10):
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, s > tol]


def build_upper_geometry(cache_dir=None):
    flow = module("upper_ps_flow", FLOW)
    common = module("upper_ps_common", COMMON)
    geometry = flow.build_ps_geometry()
    g, p1 = geometry["common_geometry"], geometry["common_geometry"]["p1"]
    upper = json.loads(UPPER.read_text())["upper_PS"]
    pars, v = upper["same_action_parameters"], upper["vacuum"]
    x = p1.vacuum_vector(v["omega"], 0., v["vs"])
    keybase = hashlib.sha256(P1.read_bytes()+json.dumps(pars, sort_keys=True).encode()).digest()
    key = hashlib.sha256(keybase+np.asarray(x, dtype="<f8").tobytes()).hexdigest()
    cache_root = ROOT / "tmp/p54_full_doublet_cw" if cache_dir is None else Path(cache_dir)
    cache = cache_root / (key+".npz")
    if not cache.exists():
        raise FileNotFoundError("Run the same-action PS finite-threshold prerequisite: "+str(cache))
    with np.load(cache) as data:
        hessian = data["h"]
    active = geometry["scalar_real_embedding_old"]
    goldstone = image(p1.gauge_orbit(x))
    pq = p1.pq_direction(x)
    pq -= goldstone @ (goldstone.T @ pq)
    pq /= np.linalg.norm(pq)
    ph = np.eye(328)-active@active.T-goldstone@goldstone.T-np.outer(pq,pq)
    d, u = np.linalg.eigh((ph+ph.T)/2)
    heavy_basis = u[:, d > .5]
    mh = heavy_basis.T @ hessian @ heavy_basis
    masses, rotation = np.linalg.eigh((mh+mh.T)/2)
    heavy = heavy_basis @ rotation
    zh, zf = common.real_scalar_yukawa_tensors(g, old_coordinates=True)
    ub = geometry["fermion_basis"]
    zh = np.einsum("pi,apq,qj->aij", ub, zh, ub)
    zf = np.einsum("pi,apq,qj->aij", ub, zf, ub)
    yh = np.einsum("ab,aij->bij", heavy, zh)
    yf = np.einsum("ab,aij->bij", heavy, zf)
    return {"flow": flow, "common": common, "ps": geometry, "upper": upper,
            "x": x, "hessian": hessian, "active": active, "heavy": heavy,
            "goldstone": goldstone, "pq": pq, "masses2": masses,
            "heavy_h": yh, "heavy_f": yf, "all_h": zh, "all_f": zf,
            "cache_path": str(cache.relative_to(ROOT)) if cache.is_relative_to(ROOT) else str(cache),
            "mu": upper["mu_over_omega_reference"]}


def family_tensor(tensor, family):
    n = family.shape[0]
    return np.einsum("aij,pq->aipjq", tensor, family).reshape(len(tensor),16*n,16*n)


def kernels(masses2, mu):
    """One heavy scalar, one massless fermion; squared masses in common units."""
    m = np.asarray(masses2)
    if mu <= 0 or np.any(m <= 0):
        raise ValueError("Positive heavy masses and scale required")
    log = np.log(m/mu**2)
    return (1-log)/LOOP, (.5*log-.25)/LOOP


def threshold_from_tensors(active_y, heavy_y, masses2, mu):
    """Finite hard 1PI vertex and fermion kinetic diagrams, in full Weyl units.

    L=-Y^a psi psi s_a/2+h.c.; Y^a symmetric; Zpsi=I+Kpsi.
    Scalar pure-Yukawa hard kinetic term is zero here: its loop has only
    massless, retained fermions and cancels against the identical EFT graph.
    """
    if np.ndim(masses2) == 2:
        masses2, rotation = np.linalg.eigh(masses2)
        heavy_y = np.einsum("ab,aij->bij",rotation,heavy_y)
    f, h = kernels(masses2, mu)
    kinetic = -np.einsum("b,bki,bkj->ij", h, heavy_y.conj(), heavy_y)
    vertex = np.zeros_like(active_y)
    adjoint = active_y.conj().transpose(0,2,1)
    for weight, y in zip(f, heavy_y):
        if np.linalg.norm(y) > 1e-14:
            vertex -= weight*(y@adjoint@y)
    legs = -(kinetic.T@active_y+active_y@kinetic)/2
    return {"Kpsi": kinetic, "vertex": vertex, "legs": legs,
            "delta": vertex+legs, "matched": active_y+vertex+legs}


def calculate(h_raw, f_raw, geometry, mu=None):
    for name, family in (("h_raw",h_raw),("f_raw",f_raw)):
        if family.shape[0] != family.shape[1] or np.linalg.norm(family-family.T)>1e-10:
            raise ValueError(name+" must be complex symmetric")
    flow = geometry["flow"]
    c = flow.spin10_boundary(h_raw, f_raw)
    active_y = flow.assemble_real_yukawas(c, geometry["ps"])
    heavy_y = family_tensor(geometry["heavy_h"],h_raw)+family_tensor(geometry["heavy_f"],f_raw)
    result = threshold_from_tensors(active_y,heavy_y,geometry["masses2"],geometry["mu"] if mu is None else mu)
    result.update({"tree":active_y,"heavy_y":heavy_y,"tree_couplings":c})
    return result


def project_declared(y, geometry):
    """Orthogonal projection onto the four declared holomorphic PS tensors.

    Return residual rather than discard any finite mirror Yukawa operators.
    """
    flow, ps = geometry["flow"], geometry["ps"]
    nf = y.shape[-1]//16
    result, offsets, offset = {}, {}, 0
    for key in flow.KEYS:
        k = ps["tensors"][key]
        part = y[offset:offset+2*len(k)]
        offsets[key] = (offset,offset+2*len(k))
        offset += 2*len(k)
        w = (part[0::2]-ps["complex_signs"][key]*1j*part[1::2])/math.sqrt(2)
        w = w.reshape(len(k),16,nf,16,nf)
        if key in ("H","F"):
            kg, wg = k[:,:8,8:],w[:,:8,:,8:,:]
        elif key == "L":
            kg,wg = k[:,:8,:8],w[:,:8,:,:8,:]
        else:
            kg,wg = k[:,8:,8:],w[:,8:,:,8:,:]
        result[key] = np.einsum("aij,aipjq->pq",kg.conj(),wg)/np.vdot(kg,kg).real
    reconstruction = flow.assemble_real_yukawas(result,ps)
    residual = y-reconstruction
    return {"couplings":result,"reconstruction":reconstruction,"residual":residual,
            "residual_norms":{key:float(np.linalg.norm(residual[a:b])) for key,(a,b) in offsets.items()},
            "offsets":offsets}


def build_mirror_basis(geometry):
    """Construct the second PS invariant from ACTUAL six-mode contractions.

    No guessed SU(5) or synthetic group tensor is used. A six-fold real
    phi mass eigenspace maps K_a^dagger to the unique conjugate-bidoublet
    intertwiner. Fix its norm to the original parent norm and its largest
    entry to positive real. The phase fixing is a recorded basis choice.
    """
    nh = np.linalg.norm(geometry["heavy_h"],axis=(1,2))
    first = np.flatnonzero(nh>1.)[0]
    mask = abs(geometry["masses2"]-geometry["masses2"][first])<1e-9
    generators = geometry["heavy_h"][mask]
    out = {}
    for key in ("H","F"):
        original = geometry["ps"]["tensors"][key]
        raw = sum(a@original.conj().transpose(0,2,1)@a for a in generators)
        if np.linalg.norm(raw)<1e-10:
            raise ValueError("No actual mirror intertwiner found for "+key)
        phase = np.exp(-1j*np.angle(raw.ravel()[np.argmax(abs(raw))]))
        out[key] = raw*phase*np.linalg.norm(original)/np.linalg.norm(raw)
    return out


def assemble_mirror(couplings, geometry, mirror):
    nf = next(iter(couplings.values())).shape[0]
    out = np.zeros((248,16*nf,16*nf),complex)
    offset = 0
    for key in geometry["flow"].KEYS:
        number = len(geometry["ps"]["tensors"][key])
        if key in mirror:
            k = mirror[key].copy()
            k[:,8:,:8] = 0
            family = couplings[key]
            w = (np.einsum("aij,pq->aipjq",k,family)+
                 np.einsum("aji,pq->aipjq",k,family.T)).reshape(number,16*nf,16*nf)
            out[offset:offset+2*number:2] = w/math.sqrt(2)
            out[offset+1:offset+2*number:2] = -geometry["ps"]["complex_signs"][key]*1j*w/math.sqrt(2)
        offset += 2*number
    return out


def project_mirror(y, geometry, mirror):
    nf,offset,result = y.shape[-1]//16,0,{}
    for key in geometry["flow"].KEYS:
        number = len(geometry["ps"]["tensors"][key])
        if key in mirror:
            part = y[offset:offset+2*number]
            w = (part[0::2]+geometry["ps"]["complex_signs"][key]*1j*part[1::2])/math.sqrt(2)
            w = w.reshape(number,16,nf,16,nf)[:,:8,:,8:,:]
            k = mirror[key][:,:8,8:]
            result[key] = np.einsum("aij,aipjq->pq",k.conj(),w)/np.vdot(k,k).real
        offset += 2*number
    return result


def mass_groups(geometry):
    values = geometry["masses2"]
    output = []
    begin = 0
    while begin<len(values):
        end = begin+1
        while end<len(values) and abs(values[end]-values[begin])<1e-9:
            end += 1
        hs = float(np.linalg.norm(geometry["heavy_h"][begin:end])**2)
        fs = float(np.linalg.norm(geometry["heavy_f"][begin:end])**2)
        output.append({"m2":float(np.mean(values[begin:end])),"multiplicity":end-begin,
                       "raw_h_tensor_norm2":hs,"raw_f_tensor_norm2":fs,
                       "Yukawa_active":hs+fs>1e-10})
        begin = end
    return output


def polynomial_coefficients(geometry, mirror, mu):
    """Exact finite cubic-spurion coefficients, obtained by group contraction.

    Their family reconstruction is tested on independent generic complex
    symmetric matrices; these are not coefficients fitted to observables.
    """
    empty = {key:np.zeros((1,1),complex) for key in geometry["flow"].KEYS}
    coeffs = {}
    fweights,hweights = kernels(geometry["masses2"],mu)
    for raw_name in ("h","f"):
        y = geometry["heavy_"+raw_name]
        k = -np.einsum("b,bki,bkj->ij",hweights,y.conj(),y)
        coeffs["K_"+raw_name] = complex(np.trace(k)/16)
        coeffs["K_"+raw_name+"_group_residual"] = error(k,np.eye(16)*np.trace(k)/16)
        for key in ("H","F"):
            c = {**empty,key:np.ones((1,1),complex)}
            a = geometry["flow"].assemble_real_yukawas(c,geometry["ps"])
            vertex = threshold_from_tensors(a,y,geometry["masses2"],mu)["vertex"]
            coeff = project_mirror(vertex,geometry,mirror)[key][0,0]
            # External parent coefficients H,F equal sqrt2*h_raw,4*f_raw.
            coeffs["mirror_"+key+"_from_"+raw_name] = complex(coeff*geometry["ps"]["scales_from_raw"][key])
    return coeffs


def full_real_weyl_beta(real_y, gauge_casimirs=None, gauge=None):
    """Full LWX real-Weyl one-loop coefficient, without pair cancellation.

    A matrix contraction evaluates V_a=sum_b Y_b Y_a^dagger Y_b. The
    non-conjugated Gram matrix is essential. This is valid after finite
    mirror matching, unlike the original holomorphic-pair shortcut.
    Returns 16*pi^2 beta, in canonical real-scalar coordinates.
    """
    count,n,_ = real_y.shape
    flat = real_y.reshape(count,n*n)
    s = np.einsum("aki,akj->ij",real_y.conj(),real_y)
    result = (s.T@real_y+real_y@s)/2
    # T_ij,kl=sum_b Y_b,ik Y_b,lj; no complex conjugation in this tensor.
    tensor = (flat.T@flat).reshape(n,n,n,n).transpose(0,3,1,2).reshape(n*n,n*n)
    adjoint = real_y.conj().transpose(0,2,1).reshape(count,n*n)
    result += 2*(adjoint@tensor.T).reshape(count,n,n)
    gram = (flat.conj()@flat.T).real
    result += np.einsum("ab,bij->aij",gram,real_y)
    if gauge is not None:
        for casimir,coupling in zip(gauge_casimirs,gauge):
            result -= 3*coupling**2*(casimir.T@real_y+real_y@casimir)
    return result


def six_parent_flow_audit(geometry, mirror, actual_matched):
    """Check complete PS-covariant dim4 Yukawa closure with existing scalars.

    General H,F,mirrorH,mirrorF and symmetric L,R are allowed. They are
    tested off the parity locus and away from the tiny finite benchmark.
    No new scalar, independent UV flavor parameter, or fit is introduced.
    """
    flow,ps = geometry["flow"],geometry["ps"]
    rng = np.random.default_rng(540966)
    general = lambda scale:scale*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    couplings = {key:general(.04) for key in flow.KEYS}
    for key in ("L","R"):
        couplings[key] = (couplings[key]+couplings[key].T)/2
    mirrors = {key:general(.025) for key in mirror}
    gauge = np.array([.57,.54,.61])
    casimirs = [np.kron(ps["spin_casimirs"][key],np.eye(3)) for key in ("SU4","SU2L","SU2R")]
    original = flow.assemble_real_yukawas(couplings,ps)
    y = original+assemble_mirror(mirrors,geometry,mirror)
    beta = full_real_weyl_beta(y,casimirs,gauge)
    declared = project_declared(beta,geometry)
    mirrored = project_mirror(beta,geometry,mirror)
    reconstructed = declared["reconstruction"]+assemble_mirror(mirrored,geometry,mirror)
    checks = []
    def check(name,residual,tolerance=2e-10):
        checks.append({"name":name,"residual":float(residual),"tolerance":tolerance,
                       "pass":bool(residual<tolerance)})
    check("full real-Weyl beta agrees with original fast path when mirrors vanish",error(
        full_real_weyl_beta(original,casimirs,gauge),flow.generic_beta(original,casimirs,gauge,paired_holomorphic=True)))
    # Direct vertex sums for three different actual parent external states
    # independently validate the fast matrix-index contraction above.
    s = np.einsum("aki,akj->ij",y.conj(),y)
    flat = y.reshape(len(y),-1)
    scalar_gram = (flat.conj()@flat.T).real
    for a in (0,19,181):
        direct = (s.T@y[a]+y[a]@s)/2+2*sum(b@y[a].conj().T@b for b in y)
        direct += np.einsum("b,bij->ij",scalar_gram[a],y)
        for c,g in zip(casimirs,gauge):
            direct -= 3*g*g*(c.T@y[a]+y[a]@c)
        check("full mirror beta matches literal Weyl diagram sum at scalar "+str(a),error(beta[a],direct))
    check("six Yukawa invariant families close under generic complex off-parity one-loop flow",error(beta,reconstructed))
    check("original four-matrix projection is not a complete beta after mirror matching",0 if
          np.linalg.norm(declared["residual"])>1e-3 else 1)
    try:
        flow.generic_beta(y,casimirs,gauge,paired_holomorphic=True)
        rejected = False
    except ValueError:
        rejected = True
    check("old complex-pair fast path fails closed on actual mirror tensor structure",0 if rejected else 1)
    # The actual finite-matched boundary is also checked, not only a large
    # generic benchmark chosen to make missing terms visible.
    actual_beta = full_real_weyl_beta(actual_matched,casimirs,gauge)
    ad = project_declared(actual_beta,geometry)
    am = project_mirror(actual_beta,geometry,mirror)
    check("actual finite-matched boundary beta is reconstructed by all six invariants",error(
        actual_beta,ad["reconstruction"]+assemble_mirror(am,geometry,mirror)))
    um = flow.unitary(rng,3);un = flow.unitary(rng,3)
    changed = flow.family_transform(couplings,um,un)
    changed_m = {key:um.T@value@un for key,value in mirrors.items()}
    yr = flow.assemble_real_yukawas(changed,ps)+assemble_mirror(changed_m,geometry,mirror)
    ef = np.zeros((48,48),complex)
    ef[:24,:24] = np.kron(np.eye(8),um)
    ef[24:,24:] = np.kron(np.eye(8),un)
    check("full mirror flow has independent U3L times U3R covariance",error(
        full_real_weyl_beta(yr,casimirs,gauge),ef.T@beta@ef))
    extra = {"H":mirrors["H"],"F":mirrors["F"],"L":np.zeros((3,3),complex),"R":np.zeros((3,3),complex)}
    measured_y4 = flow.tensor_gauge_y4(y,casimirs,[15,3,3])
    predicted_y4 = flow.gauge_y4(couplings)+flow.gauge_y4(extra)
    check("full mirror gauge Y4 includes both conjugate bidoublet norms",error(measured_y4,predicted_y4))
    return {"scope":"complete one-loop PS-invariant dimension-four Yukawa tensor-space closure test on existing248real scalars; not finite gauge/scalar matching or global fit",
            "complex_parameter_count_for_three_families":48,
            "representation_proof":"LR=(1+15,2,2), each retained H and F complex parent allows its two conjugate couplings; LL symmetric gauge=(10,3,1)+(6,1,1), only the bar10 triplet is retained; RR conjugate. Therefore H,F,Hmirror,Fmirror are general family3x3, L,R symmetric:4*9+2*6=48 complex entries. No further PS-invariant renormalizable Yukawa tensors exist for these matter/scalar fields.",
            "original_four_parent_beta_complete_after_nonzero_mirror_matching":False,
            "full_real_Weyl_beta_required":True,
            "six_invariant_closure_residual":error(beta,reconstructed),
            "four_invariant_projection_residual_norm":float(np.linalg.norm(declared["residual"])),
            "actual_matched_four_invariant_projection_residual_norm":float(np.linalg.norm(ad["residual"])),
            "actual_matched_six_invariant_closure_residual":error(actual_beta,ad["reconstruction"]+assemble_mirror(am,geometry,mirror)),
            "synthetic_generic_couplings":{key:cjson(value) for key,value in couplings.items()},
            "synthetic_generic_mirror_couplings":{key:cjson(value) for key,value in mirrors.items()},
            "beta_original_parent":{key:cjson(value) for key,value in declared["couplings"].items()},
            "beta_mirror_parent":{key:cjson(value) for key,value in mirrored.items()},
            "gauge_Y4":measured_y4.tolist(),"checks":checks}


def run():
    geometry = build_upper_geometry()
    rng = np.random.default_rng(540951)
    sym = lambda x:(x+x.T)/2
    h = sym(.10*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))))
    f = sym(.018*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))))
    result = calculate(h,f,geometry)
    projection = project_declared(result["delta"],geometry)
    mirror = build_mirror_basis(geometry)
    mirror_c = project_mirror(result["delta"],geometry,mirror)
    mirror_y = assemble_mirror(mirror_c,geometry,mirror)
    coeffs = polynomial_coefficients(geometry,mirror,geometry["mu"])
    kfamily = coeffs["K_h"]*h.conj().T@h+coeffs["K_f"]*f.conj().T@f
    predicted_mirror = {
        "H":coeffs["mirror_H_from_h"]*(h@h.conj().T@h)+coeffs["mirror_H_from_f"]*(f@h.conj().T@f),
        "F":coeffs["mirror_F_from_h"]*(h@f.conj().T@h)+coeffs["mirror_F_from_f"]*(f@f.conj().T@f)}
    checks = []
    def check(name,residual,tolerance=2e-10):
        checks.append({"name":name,"residual":float(residual),"tolerance":tolerance,
                       "pass":bool(residual<tolerance)})
    check("same-action upper finite prerequisite",0 if json.loads(UPPER.read_text())["summary"]["all_pass"] else 1)
    check("heavy55 spectrum reproduces stored upper spectrum",error(geometry["masses2"],np.asarray(geometry["upper"]["heavy_mass_eigenvalues"])))
    check("actual canonical active and heavy planes are orthogonal",np.linalg.norm(geometry["active"].T@geometry["heavy"]))
    check("all 55 integrated scalar eigenvalues are positive",0 if len(geometry["masses2"])==55 and min(geometry["masses2"])>0 else 1)
    check("only 24 real six-parent scalar modes have Yukawa vertices",0 if np.sum(
        np.linalg.norm(geometry["heavy_h"],axis=(1,2))+np.linalg.norm(geometry["heavy_f"],axis=(1,2))>1e-8)==24 else 1)
    mass_tensor = family_tensor(geometry["all_h"],h)+family_tensor(geometry["all_f"],f)
    check("all 48 Weyl fermions are massless at the actual PS saddle",np.linalg.norm(
        np.einsum("a,aij->ij",geometry["x"],mass_tensor)))
    check("heavy Goldstone and PQ planes have zero Yukawa couplings",max(
        np.linalg.norm(np.einsum("ab,aij->bij",geometry["goldstone"],mass_tensor)),
        np.linalg.norm(np.einsum("a,aij->ij",geometry["pq"],mass_tensor))))
    check("actual fermion kinetic correction is Hermitian",error(result["Kpsi"],result["Kpsi"].conj().T))
    check("heavy six contractions give a common Spin10-family kinetic matrix",error(result["Kpsi"],np.kron(np.eye(16),kfamily)))
    check("two noncommuting complex family spurions reconstruct finite mirror coefficients",max(
        error(mirror_c[key],predicted_mirror[key]) for key in mirror))
    check("four original plus two actual mirror invariants reconstruct full finite result",error(
        result["delta"],projection["reconstruction"]+mirror_y))
    check("finite mirror operators are genuinely outside the original four-matrix space",0 if min(
        projection["residual_norms"][key] for key in mirror)>1e-6 else 1)
    check("original holomorphic F=sqrt2 R survives this pure-Yukawa subset",np.linalg.norm(
        projection["couplings"]["F"]-math.sqrt(2)*projection["couplings"]["R"]))
    check("original holomorphic L=R survives this pure-Yukawa subset",error(
        projection["couplings"]["L"],projection["couplings"]["R"]))
    check("both generated mirror family matrices are symmetric",max(error(a,a.T) for a in mirror_c.values()))
    # An independent full-tensor UV-minus-EFT beta contraction fixes the
    # overall kinetic coefficient/sign, without using finite kernels.
    flow = geometry["flow"]
    # The historical coordinates store Re/Im blocks separately; reorder
    # actual UV tensors into complete adjacent pairs for the checked fast API.
    p1 = geometry["ps"]["common_geometry"]["p1"]
    order = list(range(54))
    for re,im in ((p1.SL_SIGMA_RE,p1.SL_SIGMA_IM),(p1.SL_H_RE,p1.SL_H_IM)):
        order += [index for pair in zip(range(re.start,re.stop),range(im.start,im.stop)) for index in pair]
    order += [326,327]
    buv_ordered = flow.generic_beta(mass_tensor[order],paired_holomorphic=True)
    buv = np.empty_like(buv_ordered)
    buv[order] = buv_ordered
    blo = flow.generic_beta(result["tree"],paired_holomorphic=True)
    heavy_beta = np.einsum("as,aij->sij",geometry["active"],buv)-blo
    eps = 2e-4
    plus = calculate(h,f,geometry,geometry["mu"]*math.exp(eps))
    minus = calculate(h,f,geometry,geometry["mu"]*math.exp(-eps))
    slope = (plus["delta"]-minus["delta"])/(2*eps)
    check("finite matching-scale slope equals minus actual UV-minus-EFT beta",error(slope,-heavy_beta/LOOP),1e-9)
    check("finite mirror contributions are matching-scale independent",error(plus["vertex"],minus["vertex"]))
    # Each scalar by itself has a nonzero vertex slope; this explicitly
    # checks the factor two before complete complex pairs cancel it.
    bindex = int(np.flatnonzero(np.linalg.norm(result["heavy_y"],axis=(1,2))>1e-8)[0])
    one = result["heavy_y"][bindex:bindex+1]
    one_m = geometry["masses2"][bindex:bindex+1]
    vp = threshold_from_tensors(result["tree"],one,one_m,geometry["mu"]*math.exp(eps))["vertex"]
    vm = threshold_from_tensors(result["tree"],one,one_m,geometry["mu"]*math.exp(-eps))["vertex"]
    check("single real heavy mode has the Weyl vertex factor two and sign",error((vp-vm)/(2*eps),
        -2*(one[0]@result["tree"].conj().transpose(0,2,1)@one[0])/LOOP))
    degenerate = threshold_from_tensors(result["tree"],result["heavy_y"],
        np.full(55,.31),geometry["mu"])
    check("degenerate complete heavy complex pairs cancel finite vertices",np.linalg.norm(degenerate["vertex"]))
    integral_error = 0.
    for mass in (.19,.2,.37,.42):
        got = kernels(np.array([mass]),geometry["mu"])
        fnum = -quad(lambda t:math.log((1-t)*mass/geometry["mu"]**2),0,1)[0]/LOOP
        hnum = quad(lambda t:(1-t)*math.log((1-t)*mass/geometry["mu"]**2),0,1)[0]/LOOP
        integral_error = max(integral_error,abs(got[0][0]-fnum),abs(got[1][0]-hnum))
    check("independent zero-momentum Feynman-parameter integration",integral_error)
    # Covariance transports the entire complex scalar/spinor convention.
    u,_ = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    family_rotated = calculate(u.T@h@u,u.T@f@u,geometry)
    uf = np.kron(np.eye(16),u)
    check("arbitrary complex family U3 covariance of finite vertices and legs",error(
        family_rotated["delta"],uf.T@result["delta"]@uf))
    check("arbitrary complex family U3 covariance of kinetic matching",error(
        family_rotated["Kpsi"],uf.conj().T@result["Kpsi"]@uf))
    cp = threshold_from_tensors(result["tree"].conj(),result["heavy_y"].conj(),
        geometry["masses2"],geometry["mu"])
    check("CP covariance with conjugation of couplings AND group tensors",error(cp["delta"],result["delta"].conj()))
    oh,_ = np.linalg.qr(rng.normal(size=(55,55)))
    rh = np.einsum("ab,aij->bij",oh,result["heavy_y"])
    rm = oh.T@np.diag(geometry["masses2"])@oh
    rotated = threshold_from_tensors(result["tree"],rh,rm,geometry["mu"])
    check("nondegenerate heavy real-basis covariance with full mass rotation",error(rotated["delta"],result["delta"]))
    # Complex scalar phase rotations act as real O(2) rotations. Retain the
    # exact scalar-coordinate transformation even for the generated mirrors.
    angles = rng.normal(size=124)
    oa = np.zeros((248,248))
    for j,a in enumerate(angles):
        oa[2*j:2*j+2,2*j:2*j+2] = [[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]]
    ya = np.einsum("ab,bij->aij",oa,result["tree"])
    rotated = threshold_from_tensors(ya,result["heavy_y"],geometry["masses2"],geometry["mu"])
    check("all 124 complex scalar phase choices transport finite mirror vertices",error(
        rotated["delta"],np.einsum("ab,bij->aij",oa,result["delta"])))
    # A global PQ transformation moves the saddle along its axion orbit.
    # In historical coordinates q(Sigma,h,S)=(2,-2,-4); after common global
    # conjugation the selected standard-16 fermion has q=+1 for this generator.
    # Thus every raw old-coordinate Yukawa direction has effective charge -2.
    angle = .173
    oq = np.eye(328)
    for re,im,charge in ((p1.SL_SIGMA_RE,p1.SL_SIGMA_IM,2),
                         (p1.SL_H_RE,p1.SL_H_IM,-2),(slice(326,327),slice(327,328),-4)):
        for a,b in zip(range(re.start,re.stop),range(im.start,im.stop)):
            oq[np.ix_([a,b],[a,b])] = [[math.cos(charge*angle),-math.sin(charge*angle)],
                                      [math.sin(charge*angle),math.cos(charge*angle)]]
    moved_y = np.einsum("ab,aij->bij",oq@geometry["heavy"],mass_tensor)
    pq_result = threshold_from_tensors(result["tree"],moved_y,geometry["masses2"],geometry["mu"])
    check("actual PQ orbit of heavy basis transports all Yukawa phases",error(moved_y,
        np.exp(-2j*angle)*result["heavy_y"]))
    check("generated mirrors carry the required PQ spurion phase",error(pq_result["vertex"],
        np.exp(-4j*angle)*result["vertex"]))
    check("PQ axion motion leaves fermion kinetic matching invariant",error(pq_result["Kpsi"],result["Kpsi"]))
    # A second independent family point tests the polynomial API, not a fit.
    h2 = sym(.07*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))))
    f2 = sym(.026*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))))
    second = calculate(h2,f2,geometry)
    second_mc = project_mirror(second["delta"],geometry,mirror)
    second_pred = {"H":coeffs["mirror_H_from_h"]*(h2@h2.conj().T@h2)+coeffs["mirror_H_from_f"]*(f2@h2.conj().T@f2),
                   "F":coeffs["mirror_F_from_h"]*(h2@f2.conj().T@h2)+coeffs["mirror_F_from_f"]*(f2@f2.conj().T@f2)}
    check("independent second complex family point verifies universal cubic coefficients",max(
        error(second_mc[key],second_pred[key]) for key in mirror))
    check("PS intertwiner prerequisite has no hidden normalization gate",max(
        geometry["ps"]["casimir_residual"],max(geometry["ps"]["block_leakage"].values())))
    # Ward covariance of the ACTUAL finite correction, with scalar generators
    # pulled back to the historical real coordinates used by the Hessian.
    p1 = geometry["ps"]["common_geometry"]["p1"]
    helper = module("upper_yukawa_ps_generators",RF/"code/verify_p54_p2_two_site_matching.py")
    gg = helper.ps_generators(p1)
    ward = 0.
    for matrices in gg.values():
        for generator in matrices:
            spin = geometry["ps"]["common_geometry"]["lift"](generator,geometry["ps"]["common_geometry"]["gammas"])
            ix = geometry["ps"]["common_geometry"]["indices"]
            ub = geometry["ps"]["fermion_basis"]
            spin = np.kron(ub.conj().T@spin[np.ix_(ix,ix)]@ub,np.eye(3))
            scalar = geometry["active"].T@p1.representation_matrix(generator)@geometry["active"]
            residual = spin.T@result["delta"]+result["delta"]@spin+np.einsum("ba,bij->aij",scalar,result["delta"])
            ward = max(ward,float(np.linalg.norm(residual)))
    check("all 21 PS generators preserve full finite Yukawa tensors",ward)
    # The exact group factors, computed from Clifford contractions above,
    # multiply finite logarithms of the two split six-parent masses. These
    # checks are not numerical fits of the physical family benchmark.
    yukawa_groups = [row for row in mass_groups(geometry) if row["Yukawa_active"]]
    mh = [row["m2"] for row in yukawa_groups if row["raw_h_tensor_norm2"]>1]
    mf = [row["m2"] for row in yukawa_groups if row["raw_f_tensor_norm2"]>1]
    dh,df = math.log(mh[1]/mh[0])/LOOP,math.log(mf[1]/mf[0])/LOOP
    magnitudes = {"mirror_H_from_h":3*math.sqrt(2)*dh,
                  "mirror_H_from_f":6*math.sqrt(2)*df,
                  "mirror_F_from_h":4*dh,"mirror_F_from_f":8*df}
    check("actual Clifford group factors equal 3sqrt2,6sqrt2,4,8 times split logs",max(
        abs(abs(coeffs[key])-value) for key,value in magnitudes.items()))
    for name,mass,factor in (("h",mh,3),("f",mf,6)):
        expected = -factor*sum(kernels(np.array(mass),geometry["mu"])[1])
        check("actual "+name+" kinetic group factor and finite -1/4 constant",abs(coeffs["K_"+name]-expected))
    flow_closure = six_parent_flow_audit(geometry,mirror,result["matched"])
    checks.extend(flow_closure["checks"])
    report = {
        "schema":"p54-upper-pure-yukawa-threshold-v1","date":"2026-09-06",
        "scope":"Actual same-action PS saddle and heavy55 spectrum; finite pure-Yukawa fermion kinetic and 1PI vertex graphs; synthetic family benchmarks, not full upper/lower matching or fit",
        "vacuum":geometry["upper"]["vacuum"],"mu_over_reference_omega":geometry["mu"],
        "mass_groups":mass_groups(geometry),
        "formulae":{"Kpsi":"-sum_b Yb^dagger Yb h(Mb^2,0)",
            "delta_vertex":"-sum_b Yb Ya^dagger Yb f(Mb^2,0)",
            "delta_canonical":"delta_vertex-(Kpsi^T Ya+Ya Kpsi)/2",
            "f":"[1-log(M^2/mu^2)]/(16pi^2)","h":"[log(M^2/mu^2)/2-1/4]/(16pi^2)",
            "mirror_H":"A_Hh h_raw h_raw^dagger h_raw + A_Hf f_raw h_raw^dagger f_raw",
            "mirror_F":"A_Fh h_raw f_raw^dagger h_raw + A_Ff f_raw f_raw^dagger f_raw",
            "Kfamily":"k_h h_raw^dagger h_raw+k_f f_raw^dagger f_raw",
            "original_holomorphic_delta":"-1/2(Kfamily^T C+C Kfamily), C=H,F,L,R on the tree Spin10 boundary",
            "scale_identity":"d(delta_Y)/dlog(mu)=-(beta_UV-beta_PS)_pureYukawa",
            "mirror_definition":"opposite real/imaginary scalar sign to original parent; group tensors are actual six-mode contractions, normalized to original K and largest entry real-positive"},
        "finite_group_coefficients":{key:cjson(value) if isinstance(value,complex) else value for key,value in coeffs.items()},
        "split_log_parameters":{"Delta_h":dh,"Delta_f":df,
            "definition":"Delta_x=log(m_x,high^2/m_x,low^2)/(16pi^2)",
            "actual_group_factor_magnitudes":{"H_from_h":"3sqrt2","H_from_f":"6sqrt2","F_from_h":"4","F_from_f":"8"}},
        "mirror_group_tensors":{key:cjson(value) for key,value in mirror.items()},
        "synthetic_benchmark":{"h_raw":cjson(h),"f_raw":cjson(f),"Kfamily":cjson(kfamily),
            "Kpsi_norm":float(np.linalg.norm(result["Kpsi"])),"vertex_norm":float(np.linalg.norm(result["vertex"])),
            "delta_original_parent":{key:cjson(value) for key,value in projection["couplings"].items()},
            "generated_mirror_parent":{key:cjson(value) for key,value in mirror_c.items()},
            "outside_original_four_parent_norms":projection["residual_norms"],
            "delta_total_norm":float(np.linalg.norm(result["delta"])),
            "kinetic_metric_min":float(np.linalg.eigvalsh(np.eye(48)+result["Kpsi"]).min()),
            "matching_scale_slope_norm":float(np.linalg.norm(slope))},
        "second_synthetic_benchmark":{"h_raw":cjson(h2),"f_raw":cjson(f2),
            "generated_mirror_parent":{key:cjson(value) for key,value in second_mc.items()}},
        "closure_result":{"original_holomorphic_F_sqrt2R_preserved_in_this_subset":True,
            "original_four_Yukawa_tensor_space_closed_under_finite_matching":False,
            "generated_mirror_H_and_F_nonzero":True,
            "new_fundamental_fields_or_independent_UV_flavor_parameters_added":False},
        "diagram_ledger":[
            {"graph":"heavy physical scalar plus massless fermion self-energy","status":"computed with all55 eigenstates;24 nonzero Yukawa states"},
            {"graph":"two massless fermions plus one heavy scalar triangle","status":"computed with split real masses and complete relative phases"},
            {"graph":"pure-Yukawa light scalar kinetic graph","status":"zero hard matching: both fermions are retained and massless; UV minus EFT cancels identically"},
            {"graph":"scalar-cubic plus two-Yukawa triangle","status":"zero at dimension-four and zero momentum with massless internal fermions: requires chirality-flip mass or external momentum; not a missing pure-Yukawa cubic term"},
            {"graph":"heavy-vector fermion kinetic and Yukawa vertex graphs","status":"not computed here; gauge-fixing package and BFG-Landau conventions required"},
            {"graph":"heavy-vector/scalar, Goldstone and ghost scalar kinetic package","status":"not computed here"},
            {"graph":"scalar-cubic contribution to upper scalar kinetic metric","status":"not computed here; separate from pure-Yukawa matching"},
            {"graph":"upper heavy-radial/tadpole shift and full scalar matching","status":"not computed here; tree stationary background used, one-loop mass insertion omitted consistently"},
            {"graph":"lower staged PS-to-SM functional and sequential Weinberg finite matching","status":"not computed here"}],
        "PS_flow_implication":"The previous four-parent one-loop flow remains a valid zeroth-order holomorphic flow, but its strict finite-matched boundary requires mirror operators. Since generated mirrors are one-loop, inserting them in one-loop running produces order-two-loop effects; a claimed full finite-matched coupled system must extend the tensor basis rather than discard them.",
        "PQ_and_CP":"Mirror operators have calculable PQ spurion/axion dependence. Under the historical q(S)=-4 orbit theta, mirror coefficients multiply exp(-4i theta); they are not explicit PQ breaking. The displayed i in the chosen mirror invariant basis is conventional, not a new independent CP parameter. Family couplings and all group tensors must be conjugated together under CP.",
        "six_invariant_PS_flow_closure":flow_closure,
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)},
        "cache_path":geometry["cache_path"],
        "sources":[{"path":str(path.relative_to(ROOT)),"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
                   for path in (Path(__file__),UPPER,FLOW,COMMON,P1)],
        "primary_references":[
            {"url":"https://arxiv.org/html/2310.16563v2","use":"Weyl canonical matching Eq22, finite scalar kernels Eq23-24, pure scalar diagrams Eq8-9; no SU5 group coefficients imported"},
            {"url":"https://arxiv.org/abs/hep-ph/0211440","use":"generic Weyl one-loop Yukawa beta and independent full UV-minus-EFT slope check"},
            {"url":"https://doi.org/10.1140/epjc/s10052-019-6570-5","use":"zero-momentum UV-minus-EFT hard matching and general diagram classification"}]
    }
    return report


def markdown(report):
    c = report["finite_group_coefficients"]
    scalar = lambda key: complex(c[key]["real"],c[key]["imag"])
    lines = ["# P54 upper-PS finite pure-Yukawa threshold", "",
        f"{report['summary']['passed']}/{report['summary']['total']} executable checks pass. "
        "The calculation uses the actual same-action stationary PS saddle, all 55 physical upper heavy scalar eigenstates, and the common Spin(10) Clifford intertwiners. "
        "Family benchmarks are synthetic. This is not full upper/lower finite matching or a flavor fit.", "",
        "## Main result: the four-matrix finite boundary is not closed", "",
        "The original holomorphic H,F,L,R corrections preserve F=sqrt(2)R and L=R in this pure-Yukawa subset. "
        "However, the actual split heavy real scalars generate nonzero conjugate-bidoublet Yukawa operators for H and F. "
        "Discarding the component orthogonal to the original four tensors would lose genuine finite matching terms. "
        "No fundamental field or independent UV family parameter has been added.", "",
        "## Diagram derivation and normalization", "",
        r"Use $\mathcal L_Y=-\frac12Y^a_{IJ}\psi_I\psi_Js_a+\mathrm{h.c.}$ with canonical real scalars, $Y^{aT}=Y^a$, "
        r"and $Z_\psi=I+K_\psi$. All fermions are massless at $\sigma=0$. A physical heavy scalar eigenstate $b$ has "
        r"$Y_b=\sum_A(E_H)_{Ab}Y_A$ and squared mass $M_b^2>0$.", "",
        r"The massless fermion/heavy-scalar self-energy has numerator $\slashed{k}$. Combining "
        r"$k^2[(k-p)^2-M_b^2]$ with parameter $t$ and shifting $k\mapsto\ell+(1-t)p$ gives the "
        r"finite $\slashed p$ coefficient $-h(M_b^2,0)Y_b^\dagger Y_b$, where", "",
        r"$$h(M^2,0)=\frac1{16\pi^2}\int_0^1(1-t)\log\frac{(1-t)M^2}{\mu^2}\,dt"
        r"=\frac{\frac12\log(M^2/\mu^2)-\frac14}{16\pi^2}.$$", "",
        r"In the zero-momentum triangle, two massless Weyl numerators contract to $k^2$; the remaining integral is "
        r"$[k^2(k^2-M_b^2)]^{-1}$. The oriented Yukawa chain is $Y_bY_a^\dagger Y_b$. "
        r"Its finite 1PI correction relative to the declared $-Y_a/2$ vertex is $-f(M_b^2,0)Y_bY_a^\dagger Y_b$, with", "",
        r"$$f(M^2,0)=-\frac1{16\pi^2}\int_0^1\log\frac{(1-t)M^2}{\mu^2}\,dt"
        r"=\frac{1-\log(M^2/\mu^2)}{16\pi^2}.$$", "",
        r"The $1/2$ in the symmetric Weyl action is canceled by differentiating its two fermions; no extra Majorana factor is inserted. "
        r"Combining the two external fermion legs gives", "",
        r"$$K_\psi=-\sum_{b\in H}Y_b^\dagger Y_b h_b,\qquad "
        r"\delta Y_a^{\rm 1PI}=-\sum_{b\in H}Y_bY_a^\dagger Y_b f_b,$$", "",
        r"$$\Delta Y_a=\delta Y_a^{\rm 1PI}-\frac12(K_\psi^TY_a+Y_aK_\psi).$$", "",
        "The canonical field operation and scalar kernels agree with "
        "[Patel and Shukla, Eqs. (22)–(24)](https://arxiv.org/html/2310.16563v2). "
        "No SU(5) group coefficient is imported. The actual full UV and PS Weyl tensors independently check the matching-scale slope against "
        "[Luo, Wang and Xiao](https://arxiv.org/abs/hep-ph/0211440). "
        "The subtraction follows zero-momentum hard amplitude matching as in "
        "[Gabelmann, Mühlleitner and Staub](https://doi.org/10.1140/epjc/s10052-019-6570-5).", "",
        r"Indeed $\partial_{\log\mu}f=2/(16\pi^2)$ and $\partial_{\log\mu}h=-1/(16\pi^2)$. "
        r"Consequently", "",
        r"$$\partial_{\log\mu}\Delta Y_a=-\frac1{16\pi^2}\left["
        r"\frac12(S_H^TY_a+Y_aS_H)+2\sum_{b\in H}Y_bY_a^\dagger Y_b\right],\quad "
        r"S_H=\sum_{b\in H}Y_b^\dagger Y_b.$$", "",
        "The bracket is precisely the heavy part of the pure-Yukawa beta function here. "
        "The scalar-wavefunction fermion loop contains only retained massless fermions and cancels between UV and EFT. "
        "Trace mixing between a heavy six scalar and an active parent vanishes by PS symmetry. "
        "The code checks the complete UV-minus-EFT contraction, not only derivatives of assumed kernels.", "",
        "## Actual upper spectrum", "",
        f"The matching scale is mu/reference-omega = {report['mu_over_reference_omega']:.12g}. "
        f"The same-action saddle is {report['vacuum']}.", "",
        "| Squared mass / reference omega² | Real multiplicity | Raw h norm² | Raw f norm² |",
        "|---:|---:|---:|---:|"]
    for row in report["mass_groups"]:
        lines.append(f"| {row['m2']:.12g} | {row['multiplicity']} | {row['raw_h_tensor_norm2']:.6g} | {row['raw_f_tensor_norm2']:.6g} |")
    lines += ["", "Only 24 of the 55 real heavy modes couple to fermions: the split real components of the 10 and 126 PS six parents. "
        "Gauge Goldstones and the retained PQ phase have zero Yukawa couplings at this saddle. "
        "The retained active tachyonic directions are not inserted into upper heavy logarithms.", "",
        "## Finite mirror coefficients and exact split-log structure", "",
        r"Set $h=h_{\rm raw}$, $f=f_{\rm raw}$ and $H=\sqrt2h$, $F=4f$, $L=R=2\sqrt2f$ at tree level. "
        r"The computed common fermion metric is $K_\psi=I_{16}\otimes k$,", "",
        r"$$k=k_h h^\dagger h+k_f f^\dagger f,\qquad "
        r"k_h=-3[h(m_{h,-}^2,0)+h(m_{h,+}^2,0)],\quad "
        r"k_f=-6[h(m_{f,-}^2,0)+h(m_{f,+}^2,0)].$$", "",
        f"At the actual scale: k_h={scalar('K_h').real:.14g}, k_f={scalar('K_f').real:.14g}.", "",
        "The mirror invariant basis is constructed from the actual six-mode Clifford contraction, normalized to the original parent tensor norm, "
        "and its largest entry made positive real. Its full complex tensors are exported in JSON. In that recorded basis:", "",
        r"$$\widetilde H=A_{Hh}hh^\dagger h+A_{Hf}fh^\dagger f,\qquad "
        r"\widetilde F=A_{Fh}hf^\dagger h+A_{Ff}ff^\dagger f.$$", "",
        "| Coefficient | Actual complex value | Magnitude from actual group contraction |",
        "|---|---:|---|",
        f"| A_Hh | {scalar('mirror_H_from_h'):.12g} | 3 sqrt(2) Delta_h |",
        f"| A_Hf | {scalar('mirror_H_from_f'):.12g} | 6 sqrt(2) Delta_f |",
        f"| A_Fh | {scalar('mirror_F_from_h'):.12g} | 4 Delta_h |",
        f"| A_Ff | {scalar('mirror_F_from_f'):.12g} | 8 Delta_f |", "",
        r"Here $\Delta_x=\log(m_{x,+}^2/m_{x,-}^2)/(16\pi^2)$. "
        "All four coefficients are derived from the group tensors, not fitted to family data. "
        "They are matching-scale independent and vanish when the two real masses of each complete complex six parent coincide. "
        "Both random complex symmetric family benchmarks verify the universal cubic polynomial map.", "",
        "The displayed i is a tensor-basis phase, not a new CP parameter. The PQ/axion phase is also transported: "
        "under the historical q(S)=-4 orbit by theta, the mirror coefficients acquire exp(-4i theta). "
        "Thus these operators are not explicit PQ breaking; their background spurion must be retained if the axion is restored as a dynamical field.", "",
        "## What changes and what does not", "",
        "The holomorphic four-matrix projection alone still preserves the protected F=sqrt(2)R relation in this subset. "
        "The genuinely new result is that this projection is not the entire finite boundary. "
        "The new mirror matrices are calculable nonlinear functions of the original two UV family matrices; they are not free nuisance matrices. "
        "They can change the tree two-spurion sum-rule geometry, but their observed magnitudes do not establish a successful flavor fit. "
        "Their effect inside a one-loop beta function begins at two-loop order; a resummed system that includes them must extend its tensor basis.", "",
        "## Complete dimension-four PS tensor-space closure", "",
        "Once the nonzero mirror boundary is included, the original beta_closed(H,F,L,R) is not the complete flow. "
        "Its paired-holomorphic vertex cancellation is no longer valid. The new full_real_weyl_beta implementation retains the complete real-Weyl contraction:", "",
        r"$$16\pi^2\beta_{Y_a}=\frac12(S^TY_a+Y_aS)+2\sum_bY_bY_a^\dagger Y_b"
        r"+\sum_bY_b\operatorname{Re}\operatorname{Tr}(Y_a^\dagger Y_b)"
        r"-3\sum_i g_i^2(C_i^TY_a+Y_aC_i),\qquad S=\sum_bY_b^\dagger Y_b.$$", "",
        "The existing 248 real scalars are unchanged. The code adds the two already-derived mirror intertwiners and tests completely general "
        "complex H,F,mirrorH,mirrorF, with symmetric L,R, at an unequal-gL/gR and unequal-L/R point. "
        "It compares the full beta tensor with its projection onto all six invariant families. "
        "Three literal vertex sums independently check the optimized tensor-index contraction; the zero-mirror limit reproduces the previous verifier.", "",
        "The absence of further dimension-four Yukawa tensors also follows directly from representation theory. With "
        r"$\psi_L=(4,2,1)$ and $\psi_R=(\overline4,1,2)$, their mixed product is "
        r"$(1\oplus15,2,2)$. Both H and F are real-type PS representations carried by complex scalar fields, "
        "so each admits two conjugate couplings. The symmetric left-left gauge product is "
        r"$(10,3,1)\oplus(6,1,1)$; only its conjugate triplet parent is retained. The right-right statement is conjugate. "
        "Thus the complete family tensor space has four general matrices and two symmetric matrices:", "",
        r"$$\dim_\mathbb C\mathcal Y_{\rm PS}=4n_f^2+2\frac{n_f(n_f+1)}2=48\quad(n_f=3).$$", "",
        "This counts the allowed low-energy tensor space, not independent fundamental UV parameters. "
        "PS symmetry and renormalizability force the dimension-four Yukawa counterterms to remain in this space. "
        "Higher-dimensional operators and their feedback are outside this statement.", "",
        f"Generic six-invariant closure residual: {report['six_invariant_PS_flow_closure']['six_invariant_closure_residual']:.3e}. "
        f"At the actual finite-matched benchmark the residual is {report['six_invariant_PS_flow_closure']['actual_matched_six_invariant_closure_residual']:.3e}; "
        f"discarding the two mirror directions leaves beta-tensor norm {report['six_invariant_PS_flow_closure']['actual_matched_four_invariant_projection_residual_norm']:.9g}.", "",
        "The full gauge Yukawa trace is checked too: the mirror H/F norm contributions must be included in Y4. "
        "The old fast path now rejects, rather than silently mis-evolves, this non-holomorphic-pair input. "
        "The new function supplies the full beta tensor and a verified six-invariant projection; it does not claim a completed global trajectory/fit.", "",
        "## Missing-diagram ledger", "",
        "| Graph | Status |", "|---|---|"]
    for row in report["diagram_ledger"]:
        lines.append(f"| {row['graph']} | {row['status']} |")
    lines += ["", "The pure-Yukawa subset is independent of gauge fixing, but this does not provide the missing "
        "background-field Landau vector/Goldstone/ghost package. No unknown finite contribution is silently set to zero. "
        "One-loop mass-counterterm insertions in these one-loop diagrams would be selected two-loop terms and are not included. "
        "The separate broken-background four-doublet scalar-cubic kinetic calculation is a one-step hard-scalar subset; "
        "it is neither this upper-only Yukawa calculation nor complete lower staged matching.", "",
        "## Executable validation", "", "| Check | Residual | Pass |", "|---|---:|:---:|"]
    for check in report["checks"]:
        lines.append(f"| {check['name']} | {check['residual']:.3e} | {'yes' if check['pass'] else 'NO'} |")
    lines += ["", "API: build_upper_geometry(); calculate(h_raw,f_raw,geometry,mu); "
        "threshold_from_tensors accepts either diagonal heavy squared masses or a full real symmetric heavy mass matrix; "
        "build_mirror_basis/project_mirror/assemble_mirror preserve the generated invariant directions. "
        "The source hashes and original heavy-Hessian cache key are recorded in JSON.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    report = run()
    OUT.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(report))
    print(json.dumps({"summary":report["summary"],"coefficients":report["finite_group_coefficients"],
                      "failed":[c for c in report["checks"] if not c["pass"]]},indent=2))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)
