#!/usr/bin/env python3
"""Full four-parent P54 Pati--Salam Yukawa flow, in canonical tensor units.

This is an MS-bar one-loop matrix beta-function implementation, not a flavor
fit and not a finite threshold calculation. Actual Spin(10) Clifford tensors
are projected onto all four active PS scalar parents. Both Majorana parents
are retained. Bidoublet matrices are general complex off the parity locus.
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
from scipy.integrate import solve_ivp

RF = Path(__file__).resolve().parents[1]
OUTPUT = RF / "output/p54_ps_yukawa_flow"
COMMON = RF / "code/verify_p54_common_yukawa_phase.py"
PARENTS = RF / "code/verify_p54_p2_two_site_matching.py"
KEYS = ("H", "F", "L", "R")
LOOP = 16 * np.pi**2
GAUGE_A = np.array([2/3, 26/3, 26/3])
GAUGE_B = np.array([[3551/6, 249/2, 249/2],
                    [1245/2, 779/3, 48], [1245/2, 48, 779/3]])


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(obj)
    return obj


def norm_error(a, b):
    return float(np.linalg.norm(a-b) / max(1.0, np.linalg.norm(a), np.linalg.norm(b)))


def cjson(value):
    a = np.asarray(value)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def build_ps_geometry():
    """Return canonical PS intertwiners and their actual scalar embeddings.

    tensors[key] has shape (complex_parent_dimension,16,16), with fermion
    order (eight L, eight R). normalizations refer to raw Spin(10) family
    matrices: H=sqrt(2) h_raw, F=4 f_raw, L=R=2sqrt(2) f_raw.
    scalar_real_embedding_old maps interleaved parent real coordinates into
    the old 328-coordinate bosonic basis; fermion_basis maps ordered PS
    fermions to the common, standard-hypercharge 16-component spinor.
    """
    common = module("p54_flow_common", COMMON)
    g = common.build_geometry()
    p1 = g["p1"]
    parents = module("p54_flow_parents", PARENTS)
    gens = parents.ps_generators(p1)
    spin = {key: [g["lift"](a, g["gammas"])[np.ix_(g["indices"], g["indices"])]
                  for a in value] for key, value in gens.items()}
    spin_c = {key: -sum(a@a for a in value) for key, value in spin.items()}
    eig, vec = np.linalg.eigh(spin_c["SU2L"])
    ul, ur = vec[:, abs(eig-.75) < 1e-10], vec[:, abs(eig) < 1e-10]
    assert ul.shape == ur.shape == (16, 8)
    fermion_basis = np.column_stack((ul, ur))
    u = g["U126_physical"]
    reps = {key: [u.conj().T @ g["wedge_matrix"](a, p1.QUINTS) @ u for a in value]
            for key, value in gens.items()}
    casimir = {key: -sum(a@a for a in value) for key, value in reps.items()}
    combined = casimir["SU4"] + math.sqrt(2)*casimir["SU2L"] + np.pi*casimir["SU2R"]
    ev, vv = np.linalg.eigh(combined)
    targets = {"F": (4, .75, .75, 60), "L": (4.5, 2, 0, 30), "R": (4.5, 0, 2, 30)}
    scalar_basis = {"H": np.eye(10, dtype=complex)[:, 6:10]}
    casimir_residual = 0.0
    for key, (c4, cl, cr, dimension) in targets.items():
        value = c4 + math.sqrt(2)*cl + np.pi*cr
        basis = vv[:, abs(ev-value)<1e-9]
        assert basis.shape == (126, dimension)
        scalar_basis[key] = basis
        for group, c in zip(("SU4", "SU2L", "SU2R"), (c4, cl, cr)):
            casimir_residual = max(casimir_residual, norm_error(casimir[group]@basis, c*basis))
    scales = {"H": math.sqrt(2), "F": 4., "L": 2*math.sqrt(2), "R": 2*math.sqrt(2)}
    tensors, embeddings, signs = {}, {}, {}
    block_leakage = {}
    for key in KEYS:
        sb = scalar_basis[key]
        raw = g["b10"] if key == "H" else g["b126"]
        # H couples to phi*, whereas the three 126 parents couple to Sigma.
        local = np.einsum("as,aij->sij", sb.conj() if key == "H" else sb, raw)
        local = np.einsum("pi,apq,qj->aij", fermion_basis, local, fermion_basis) / scales[key]
        mask = np.zeros((16, 16))
        if key in ("H", "F"):
            mask[:8, 8:] = 1
            mask[8:, :8] = 1
        elif key == "L":
            mask[:8, :8] = 1
        else:
            mask[8:, 8:] = 1
        block_leakage[key] = float(np.linalg.norm(local*(1-mask)))
        tensors[key] = local*mask
        emb = np.zeros((328, 2*sb.shape[1]))
        re = p1.SL_H_RE if key == "H" else p1.SL_SIGMA_RE
        im = p1.SL_H_IM if key == "H" else p1.SL_SIGMA_IM
        emb[re, 0::2], emb[re, 1::2] = sb.real, -sb.imag
        emb[im, 0::2], emb[im, 1::2] = sb.imag, sb.real
        embeddings[key] = g["conjugation_diagonal"][:, None]*emb
        signs[key] = -1 if key == "H" else 1
    embedding = np.column_stack([embeddings[key] for key in KEYS])
    return {"tensors": tensors, "scalar_complex_bases": scalar_basis,
            "scalar_real_embedding_old": embedding, "scalar_embeddings_old": embeddings,
            "fermion_basis": fermion_basis, "complex_signs": signs,
            "spin_casimirs": {key: fermion_basis.conj().T@value@fermion_basis
                              for key, value in spin_c.items()},
            "scales_from_raw": scales, "casimir_residual": casimir_residual,
            "block_leakage": block_leakage, "common_geometry": g}


def complex_yukawas(couplings, geometry):
    """Parent-ordered complex scalar Yukawa tensors, gauge-major family-minor."""
    nf = couplings["H"].shape[0]
    values = {}
    for key in KEYS:
        k = geometry["tensors"][key]
        family = couplings[key]
        if key in ("H", "F"):
            klr = k.copy()
            klr[ :, 8:, :8] = 0
            y = np.einsum("aij,pq->aipjq", klr, family)
            y += np.einsum("aji,pq->aipjq", klr, family.T)
        else:
            if np.linalg.norm(family-family.T)>1e-9:
                raise ValueError(f"{key} must be complex symmetric")
            y = np.einsum("aij,pq->aipjq", k, family)
        values[key] = y.reshape(len(k), 16*nf, 16*nf)
    return values


def assemble_real_yukawas(couplings, geometry):
    """248 real-scalar tensors, ordered H,F,L,R and real,imag per column.

    Scalar coordinates are exactly scalar_real_embedding_old. Fermion basis
    is geometry['fermion_basis'] tensor family identity. The 4 family inputs
    are PS normalized H,F,L,R, NOT all the same raw Spin(10) matrix.
    """
    complex_values = complex_yukawas(couplings, geometry)
    output = []
    for key in KEYS:
        value = complex_values[key]
        pair = np.empty((2*len(value), *value.shape[1:]), complex)
        pair[0::2] = value/math.sqrt(2)
        pair[1::2] = geometry["complex_signs"][key]*1j*value/math.sqrt(2)
        output.append(pair)
    return np.concatenate(output)


def generic_beta(real_y, gauge_casimirs=None, gauge=None, paired_holomorphic=False):
    """Luo--Wang--Xiao Eq.(33), two-component Weyl convention kappa=1/2.

    Return 16 pi^2 beta. For symmetric full Weyl Yukawa matrices the left
    wavefunction is S^T, not S, with S=sum Y^dagger Y. The optional paired
    flag is valid only after Y_im=+/-i Y_re has been checked. Such complete
    complex pairs cancel the vertex term exactly, reducing cubic cost.
    """
    if paired_holomorphic:
        if len(real_y)%2:
            raise ValueError("Holomorphic fast path requires complete adjacent real/imaginary scalar pairs.")
        even,odd = real_y[0::2],real_y[1::2]
        axes = tuple(range(1,even.ndim))
        mismatch = np.minimum(np.linalg.norm(odd-1j*even,axis=axes),
                              np.linalg.norm(odd+1j*even,axis=axes))
        scale = np.maximum(np.linalg.norm(even,axis=axes),np.linalg.norm(odd,axis=axes))
        if np.any(mismatch>2e-12*np.maximum(scale,np.finfo(float).tiny)):
            raise ValueError("Holomorphic fast path requires Y_im=+/-i Y_re pairwise; "
                             "use the full vertex contraction after unpaired or mixed real-mode decoupling.")
    s = np.einsum("aki,akj->ij", real_y.conj(), real_y)
    answer = (s.T@real_y + real_y@s)/2
    if not paired_holomorphic:
        for a in range(len(real_y)):
            answer[a] += 2*sum(y@real_y[a].conj().T@y for y in real_y)
    flat = real_y.reshape(len(real_y), -1)
    gram = (flat.conj()@flat.T).real
    answer += np.einsum("ab,bij->aij", gram, real_y)
    if gauge is not None:
        for c, coupling in zip(gauge_casimirs, gauge):
            answer -= 3*coupling**2*(c.T@real_y+real_y@c)
    return answer


def anomalous_matrices(c):
    h, f, l, r = (c[key] for key in KEYS)
    return (h@h.conj().T + 15/4*f@f.conj().T + 15/2*l@l.conj().T,
            h.conj().T@h + 15/4*f.conj().T@f + 15/2*r.conj().T@r)


def beta_closed(c, gauge):
    """Full four-parent one-loop coefficient, i.e. 16 pi^2 dY/dlog(mu)."""
    h, f, l, r = (c[key] for key in KEYS)
    al, ar = anomalous_matrices(c)
    g4, gl, gr = gauge
    gd = 9/4*(gl**2+gr**2+5*g4**2)
    gleft = 9/4*(2*gl**2+5*g4**2)
    gright = 9/4*(2*gr**2+5*g4**2)
    trace = lambda y: float(np.vdot(y, y).real)
    return {"H": al@h+h@ar+(4*trace(h)-gd)*h,
            "F": al@f+f@ar+(trace(f)-gd)*f,
            "L": al@l+l@al.T+(2*trace(l)-gleft)*l,
            "R": ar.T@r+r@ar+(2*trace(r)-gright)*r}


def gauge_y4(c):
    """Weyl Yukawa contribution to two-loop gauge running, order (4,L,R)."""
    al, ar = anomalous_matrices(c)
    return np.array([2*np.trace(al+ar).real, 4*np.trace(al).real, 4*np.trace(ar).real])


def beta_gauge(c, gauge):
    g = np.asarray(gauge)
    return g**3*(GAUGE_A/LOOP + (GAUGE_B@(g*g)-gauge_y4(c))/LOOP**2)


def spin10_boundary(h_raw, f_raw):
    """Tree Spin(10)-normalized boundary; finite matching must be added later."""
    return {"H": math.sqrt(2)*h_raw, "F": 4*f_raw,
            "L": 2*math.sqrt(2)*f_raw, "R": 2*math.sqrt(2)*f_raw}


def integrate(c, gauge, interval=(0., math.log(8))):
    n = c["H"].shape[0]
    m = n*n
    def pack(values, gs):
        z = np.concatenate([values[key].ravel() for key in KEYS])
        return np.r_[z.real, z.imag, gs]
    def unpack(y):
        z = y[:4*m]+1j*y[4*m:8*m]
        return {key: z[k*m:(k+1)*m].reshape(n,n) for k,key in enumerate(KEYS)}, y[-3:]
    def rhs(t,y):
        matrices, gs = unpack(y)
        return pack({k:v/LOOP for k,v in beta_closed(matrices,gs).items()}, beta_gauge(matrices,gs))
    solution = solve_ivp(rhs, interval, pack(c,gauge), method="DOP853", rtol=2e-10, atol=2e-12)
    if not solution.success:
        raise RuntimeError(solution.message)
    out, gs = unpack(solution.y[:,-1])
    return out, gs, {"nfev": solution.nfev, "log_scale_interval": list(interval), "success": True}


def tensor_gauge_y4(y, casimirs, dimensions):
    s = np.einsum("aik,ajk->ij", y, y.conj())
    return np.array([np.trace(c@s).real/d for c,d in zip(casimirs,dimensions)])


def unitary(rng, n):
    z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,r = np.linalg.qr(z)
    return q*np.exp(-1j*np.angle(np.diag(r)))[None,:]


def family_transform(c, ul, ur):
    return {"H": ul.T@c["H"]@ur, "F": ul.T@c["F"]@ur,
            "L": ul.T@c["L"]@ul, "R": ur.T@c["R"]@ur}


def parity_transform(c):
    return {"H": c["H"].T, "F": c["F"].T, "L": c["R"], "R": c["L"]}


def literature_no_delta_l(c, gauge):
    """Literal regression of arXiv:1612.07973 (76)--(78); NOT the P54 EFT."""
    h,f,r = c["H"],c["F"],c["R"]
    left = h@h.conj().T+15/4*f@f.conj().T
    right = h.conj().T@h+15/4*(f.conj().T@f+2*r.conj()@r)
    g4,gl,gr = gauge
    return {"H": left@h+h@right+4*np.trace(h@h.conj().T)*h
                  -9/4*(gl*gl+gr*gr+5*g4*g4)*h,
            "F": left@f+f@right+np.trace(f@f.conj().T)*f
                  -9/4*(gl*gl+gr*gr+5*g4*g4)*f,
            "R": (h.T@h.conj()+15/4*(f.T@f.conj()+2*r@r.conj()))@r
                  +r@right+2*np.trace(r@r.conj())*r
                  -9/4*(2*gr*gr+5*g4*g4)*r}


def sm_top_regression():
    """Independent SM top-only tensor, including the full vertex contraction."""
    coupling = .63
    k = np.zeros((2,9,9),complex)
    for weak in range(2):
        for color in range(3):
            q,u = 2*color+weak,6+color
            k[weak,q,u] = k[weak,u,q] = coupling
    y = np.empty((4,9,9),complex)
    y[0::2],y[1::2] = k/math.sqrt(2),1j*k/math.sqrt(2)
    c3 = np.eye(9)*4/3
    c2 = np.diag([.75]*6+[0]*3)
    c1 = np.diag([3/5*(1/6)**2]*6+[3/5*(2/3)**2]*3)
    gs = np.array([1.01,.63,.47])  # g3,g2,GUT-normalized g1
    beta = generic_beta(y,[c3,c2,c1],gs)
    coefficient = 9/2*coupling**2-8*gs[0]**2-9/4*gs[1]**2-17/20*gs[2]**2
    expected_y4 = coupling**2*np.array([2.,1.5,1.7])
    return {"beta_relative_residual": norm_error(beta,coefficient*y),
            "gauge_Y4_relative_residual": norm_error(tensor_gauge_y4(y,[c3,c2,c1],[8,3,1]),expected_y4),
            "expected_cubic_coefficient": 4.5,
            "gauge_coefficient_order_g3_g2_g1": [8.,2.25,.85]}


def run():
    geometry = build_ps_geometry()
    rng = np.random.default_rng(543126)
    general = lambda: .045*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    c = {key: general() for key in KEYS}
    for key in ("L", "R"):
        c[key] = (c[key]+c[key].T)/2
    gauge = np.array([.57,.51,.62])
    y = assemble_real_yukawas(c,geometry)
    nf = 3
    gauge_c = [np.kron(geometry["spin_casimirs"][key],np.eye(nf))
               for key in ("SU4","SU2L","SU2R")]
    generic = generic_beta(y,gauge_c,gauge,paired_holomorphic=True)
    expected = assemble_real_yukawas(beta_closed(c,gauge),geometry)
    # Sector contraction diagnostics expose any normalization error.
    contractions = {}
    for key in KEYS:
        t = geometry["tensors"][key]
        s = np.einsum("aki,akj->ij", t.conj(),t)
        contractions[key] = {"left_S": float(np.trace(s[:8,:8]).real/8),
                             "right_S": float(np.trace(s[8:,8:]).real/8),
                             "scalar_norm": float(np.vdot(t[0],t[0]).real)}
    checks = []
    def check(name, residual, tolerance=2e-11):
        checks.append({"name": name, "residual": float(residual), "tolerance": tolerance,
                       "passed": bool(residual<tolerance)})
    check("actual PS Casimir eigenspaces",geometry["casimir_residual"])
    check("four parents have only the allowed chiral fermion blocks",max(geometry["block_leakage"].values()))
    sm_scalar = geometry["common_geometry"]["scalar126_sm"]
    scalar_c3 = -sum(a@a for a in sm_scalar["SU3"])
    scalar_y = 1j*sm_scalar["Y"][0]
    colorless_y = {}
    for key in ("L","R"):
        sb = geometry["scalar_complex_bases"][key]
        color_values,color_vectors = np.linalg.eigh(sb.conj().T@scalar_c3@sb)
        singlets = sb@color_vectors[:,abs(color_values)<1e-10]
        colorless_y[key] = np.linalg.eigvalsh(singlets.conj().T@scalar_y@singlets).tolist()
    check("physical conjugated triplet labels from actual colorless hypercharges",max(
        norm_error(np.asarray(colorless_y["L"]),np.array([1.,1.,1.])),
        norm_error(np.asarray(colorless_y["R"]),np.array([-2.,-1.,0.]))))
    embedding = geometry["scalar_real_embedding_old"]
    check("248 canonical real scalar directions are orthonormal",norm_error(embedding.T@embedding,np.eye(248)))
    check("full Weyl Yukawa matrices are complex symmetric",norm_error(y,y.transpose(0,2,1)))
    check("all-active generic tensor beta equals closed four-matrix system",norm_error(generic,expected))
    target_contractions = {"H":(2,2,8),"F":(7.5,7.5,2),"L":(15,0,4),"R":(0,15,4)}
    check("actual raw normalization and all parent contraction constants",max(
        max(abs(contractions[key][field]-target) for field,target in zip(
            ("left_S","right_S","scalar_norm"),target_contractions[key])) for key in KEYS))
    a = y[17].conj().T
    vertex = y[0::2]@a@y[0::2]+y[1::2]@a@y[1::2]
    check("complete complex scalar pairs cancel the vertex contraction",float(np.linalg.norm(vertex)))
    invalid_pair_cases = {}
    broken_pair = y[:2].copy()
    broken_pair[1] += .17*broken_pair[0]
    for name,invalid_y in (("odd_real_scalar_count",y[:3]),("nonholomorphic_pair",broken_pair)):
        try:
            generic_beta(invalid_y,paired_holomorphic=True)
        except ValueError:
            invalid_pair_cases[name] = True
        else:
            invalid_pair_cases[name] = False
    check("holomorphic fast-path guard rejects incomplete and nonholomorphic pairs",0 if
          all(invalid_pair_cases.values()) else 1)
    ul,ur = unitary(rng,3),unitary(rng,3)
    changed = family_transform(c,ul,ur)
    transformed_beta = family_transform(beta_closed(c,gauge),ul,ur)
    check("independent U(3)L times U(3)R family covariance",max(norm_error(
        beta_closed(changed,gauge)[key],transformed_beta[key]) for key in KEYS))
    ef = np.zeros((48,48),complex)
    ef[:24,:24] = np.kron(np.eye(8),ul)
    ef[24:,24:] = np.kron(np.eye(8),ur)
    check("actual intertwiners implement family covariance",norm_error(
        assemble_real_yukawas(changed,geometry),ef.T@y@ef))
    pg = gauge[[0,2,1]]
    pb = parity_transform(beta_closed(c,gauge))
    check("left-right exchange covariance including unequal gauge couplings",max(norm_error(
        beta_closed(parity_transform(c),pg)[key],pb[key]) for key in KEYS))
    check("Majorana beta matrices remain complex symmetric",max(norm_error(
        beta_closed(c,gauge)[key],beta_closed(c,gauge)[key].T) for key in ("L","R")))
    angles = rng.normal(size=124)
    rot = np.zeros((248,248))
    for j,angle in enumerate(angles):
        rot[2*j:2*j+2,2*j:2*j+2] = [[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
    ry = np.einsum("ab,bij->aij",rot,y)
    check("canonical scalar phase-basis covariance of generic beta",norm_error(
        generic_beta(ry,gauge_c,gauge,paired_holomorphic=True),np.einsum("ab,bij->aij",rot,generic)))
    scalar_beta = generic_beta(y,None,None,paired_holomorphic=True)
    zero_gauge_expected = assemble_real_yukawas(beta_closed(c,np.zeros(3)),geometry)
    check("Yukawa-only tensor contraction independently separates gauge factors",norm_error(scalar_beta,zero_gauge_expected))
    no_l = {**c,"L":np.zeros((3,3),complex)}
    reduced = literature_no_delta_l(no_l,gauge)
    check("DeltaL-absent limit agrees with corrected literature equations",max(norm_error(
        beta_closed(no_l,gauge)[key],reduced[key]) for key in ("H","F","R")))
    y4_actual = tensor_gauge_y4(y,gauge_c,[15,3,3])
    check("two-loop gauge Yukawa trace agrees with closed Y4 formula",norm_error(y4_actual,gauge_y4(c)))
    check("two-loop gauge Yukawa trace is family-basis invariant",norm_error(gauge_y4(changed),gauge_y4(c)))
    sm = sm_top_regression()
    check("SM top-only 9/2 cubic and 8,9/4,17/20 gauge coefficients",sm["beta_relative_residual"])
    check("SM top-only two-loop gauge Yukawa coefficients",sm["gauge_Y4_relative_residual"])
    final,gfinal,ode = integrate(c,gauge)
    final_cov,gcov,_ = integrate(changed,gauge)
    transformed_final = family_transform(final,ul,ur)
    check("nonzero numerical ODE respects independent family basis covariance",max(
        max(norm_error(final_cov[key],transformed_final[key]) for key in KEYS),norm_error(gfinal,gcov)))
    check("numerical ODE preserves Majorana symmetry",max(norm_error(final[key],final[key].T) for key in ("L","R")))
    ode["initial_couplings"] = {key:cjson(c[key]) for key in KEYS}
    ode["final_couplings"] = {key:cjson(final[key]) for key in KEYS}
    ode["initial_gauge"] = gauge.tolist()
    ode["final_gauge"] = gfinal.tolist()
    ode["change_norms"] = {key:float(np.linalg.norm(final[key]-c[key])) for key in KEYS}
    check("short ODE actually evolves all four nonzero matrices",0 if min(ode["change_norms"].values())>1e-5 else 1)
    hraw,fraw = general(),.25*general()
    hraw,fraw = (hraw+hraw.T)/2,(fraw+fraw.T)/2
    boundary = spin10_boundary(hraw,fraw)
    common = module("p54_flow_common_check",COMMON)
    zh,zf = common.real_scalar_yukawa_tensors(geometry["common_geometry"],old_coordinates=True)
    zh = np.einsum("as,aij->sij",embedding,zh)
    zf = np.einsum("as,aij->sij",embedding,zf)
    ub = geometry["fermion_basis"]
    zh = np.einsum("pi,apq,qj->aij",ub,zh,ub)
    zf = np.einsum("pi,apq,qj->aij",ub,zf,ub)
    raw = (np.einsum("aij,pq->aipjq",zh,hraw)+np.einsum("aij,pq->aipjq",zf,fraw)).reshape(248,48,48)
    check("Spin10 boundary exactly reconstructs actual raw common-phase action",norm_error(
        assemble_real_yukawas(boundary,geometry),raw))
    gp = np.array([.58,.56,.56])
    parity_final,parity_g,parity_ode = integrate(boundary,gp,(0.,-math.log(10)))
    parity_ode.update({"initial_raw_h":cjson(hraw),"initial_raw_f":cjson(fraw),
                       "initial_couplings":{key:cjson(boundary[key]) for key in KEYS},
                       "final_couplings":{key:cjson(parity_final[key]) for key in KEYS},
                       "initial_gauge":gp.tolist(),"final_gauge":parity_g.tolist()})
    check("actual Spin10 boundary flow preserves symmetric bidoublets",max(norm_error(
        parity_final[key],parity_final[key].T) for key in ("H","F")))
    check("actual Spin10 boundary flow preserves L=R and gL=gR",max(
        norm_error(parity_final["L"],parity_final["R"]),abs(parity_g[1]-parity_g[2])))
    parity_ode["protected_126_ratio_residual_F_minus_sqrt2_R"] = float(np.linalg.norm(
        parity_final["F"]-math.sqrt(2)*parity_final["R"]))
    check("parity locus protects F=sqrt2 R at one loop",parity_ode["protected_126_ratio_residual_F_minus_sqrt2_R"])
    asym_case = {**c,"H":(c["H"]+c["H"].T)/2,"F":(c["F"]+c["F"].T)/2}
    asym_beta = beta_closed(asym_case,gauge)
    antisymmetric_source = {key:float(np.linalg.norm(asym_beta[key]-asym_beta[key].T)) for key in ("H","F")}
    check("unequal L/R produces a genuine antisymmetric bidoublet beta source",0 if
          min(antisymmetric_source.values())>1e-6 else 1)
    broken_gauge = np.array([.58,.53,.61])
    broken_final,_,broken_ode = integrate(boundary,broken_gauge,(0.,-math.log(10)))
    broken_ode["F_minus_sqrt2_R_norm"] = float(np.linalg.norm(broken_final["F"]-math.sqrt(2)*broken_final["R"]))
    broken_ode["L_minus_R_norm"] = float(np.linalg.norm(broken_final["L"]-broken_final["R"]))
    check("unequal gL/gR negative control breaks protected 126 ratios",0 if min(
        broken_ode["F_minus_sqrt2_R_norm"],broken_ode["L_minus_R_norm"])>1e-5 else 1)
    result = {"status":"derived_and_numerically_verified_not_a_flavor_fit",
              "all_checks_pass":all(x["passed"] for x in checks),"checks":checks,
              "checks_passed":sum(x["passed"] for x in checks),"checks_total":len(checks),
              "scheme":"MS-bar, one-loop Yukawa; two-loop gauge including Yukawa term; t=log(mu)",
              "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),COMMON,PARENTS)},
              "parent_dimensions_complex":{"H":4,"F":60,"L":30,"R":30},
              "physical_scalar_labels":{"H":"(1,2,2)","F":"(15,2,2)",
                                       "L":"(bar10,3,1)","R":"(10,1,3)"},
              "triplet_colorless_hypercharges":colorless_y,
              "holomorphic_fast_path_negative_controls":invalid_pair_cases,
              "normalization_from_raw":geometry["scales_from_raw"],
              "parent_contraction_constants":contractions,
              "gauge_a_order_4_L_R":GAUGE_A.tolist(),"gauge_b_without_yukawa":GAUGE_B.tolist(),
              "synthetic_point_Y4":y4_actual.tolist(),"sm_top_regression":sm,
              "unequal_LR_antisymmetric_bidoublet_beta":antisymmetric_source,
              "general_LR_ODE":ode,"Spin10_boundary_parity_ODE":parity_ode,
              "broken_parity_negative_control_ODE":broken_ode,
              "API":{"geometry":"build_ps_geometry()",
                     "tensors":"assemble_real_yukawas({H,F,L,R},geometry) -> (248,48,48)",
                     "scalar_coordinates":"geometry['scalar_real_embedding_old'] -> (328,248)",
                     "fermion_coordinates":"geometry['fermion_basis'] -> (16,16), ordered L then R",
                     "beta":"beta_closed(couplings,gauge) returns 16*pi^2*beta_Y",
                     "gauge_beta":"beta_gauge returns full dg/dlog(mu), order 4,L,R"},
              "limitations":["Synthetic family entries are not fitted observables.",
                             "Tree Spin10 boundary is checked, not finite Yukawa threshold matching.",
                             "EFT valid only while all four complete PS parents are active.",
                             "No two-loop Yukawa beta, physical fit, pole mass or determinant promotion."]}
    OUTPUT.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    write_report(result)
    print(json.dumps({"all_checks_pass":result["all_checks_pass"],"checks_passed":result["checks_passed"],
                      "checks_total":len(checks),"tensor_residual":norm_error(generic,expected),
                      "failed":[x for x in checks if not x["passed"]]},indent=2))
    if not result["all_checks_pass"]:
        raise SystemExit(1)
    return result


def write_report(result):
    lines = [
        "# P54: complete four-parent Pati–Salam Yukawa evolution", "",
        f"Status: {result['checks_passed']}/{result['checks_total']} executable checks pass. "
        "This closes the all-active one-loop Yukawa contraction/ODE task, not a physical fit or finite matching.", "",
        "## Convention and derivation", "",
        r"Use canonical real scalars and two-component Weyl fields: "
        r"$\mathcal L_Y=-\frac12Y^a_{IJ}\psi_I\psi_J\varphi_a+\mathrm{h.c.}$, $Y^{aT}=Y^a$. "
        r"For $S=\sum_bY^{b\dagger}Y^b$ the one-loop coefficient is", "",
        r"$$\mathcal B^a=\tfrac12(S^TY^a+Y^aS)+2\sum_bY^bY^{a\dagger}Y^b"
        r"+\sum_bY^b\operatorname{Re}\operatorname{Tr}(Y^{a\dagger}Y^b)"
        r"-3\sum_i g_i^2(C_i^TY^a+Y^aC_i),\qquad16\pi^2\dot Y^a=\mathcal B^a.$$", "",
        "The tensor formula and Weyl factor are from "
        "[Luo, Wang and Xiao, Eqs. (19), (30)–(34)](https://arxiv.org/pdf/hep-ph/0211440). "
        "The transposed left wavefunction follows from symmetric Weyl tensors; replacing it by S fails complex family covariance.", "",
        r"Actual Clifford tensors are projected by PS Casimirs. Retain "
        r"$H:(1,2,2)$, $F:(15,2,2)$, $L:(\overline{10},3,1)$, $R:(10,1,3)$ "
        r"with complex dimensions $4,60,30,30$. Each complex scalar gives "
        r"$Y_R=K/\sqrt2$, $Y_I=s iK/\sqrt2$, with $s=-1$ for $\phi^*$ and $s=+1$ for $\Sigma$. "
        r"Thus $Y_RX Y_R+Y_IXY_I=0$ for any $X$: the full vertex sum cancels. "
        "This cancellation requires complete complex pairs. The executable fast path checks even scalar count "
        "and the pairwise relation with a scale-relative tolerance, rejecting incomplete or nonholomorphic pairs. "
        "After arbitrary real-mode mass decoupling the full vertex contraction must be retained.", "",
        r"These are physical labels after the common global conjugation, with "
        r"$4=(3,B-L=1/3)+(1,B-L=-1)$. The actual color-singlet scalar hypercharges are "
        r"$(1,1,1)$ for L and $(-2,-1,0)$ for R. Thus the historical literal-P1 "
        r"$10/\overline{10}$ names must be swapped; their Casimirs, indices and beta coefficients are unchanged.", "",
        "The following group contractions are computed from the actual normalized K tensors, with family matrices factored out. "
        "The first two columns are eigenvalues of sum K†K on eight L/R gauge states; the last is the norm of any complex scalar component.", "",
        "| Parent | left S | right S | Tr K†K per complex scalar |",
        "|---|---:|---:|---:|",
        "| H | 2 | 2 | 8 |", "| F | 15/2 | 15/2 | 2 |",
        "| L | 15 | 0 | 4 |", "| R | 0 | 15 | 4 |", "",
        r"This fixes the PS units $H=\sqrt2h_{\rm raw}$, $F=4f_{\rm raw}$, "
        r"$L=R=2\sqrt2f_{\rm raw}$ at the tree Spin(10) boundary. "
        r"In the common component-phase convention $f_D=2f_{\rm raw}/\sqrt3$ and $f_M=+i4\sqrt2f_{\rm raw}$: "
        r"$F=2\sqrt3f_D$, $R=\sqrt6f_D$, $f_M=+2iR$. "
        "These are normalization-aware boundary maps, not four new UV family matrices.", "",
        "## Complete one-loop matrix system", "",
        r"For general complex H,F and symmetric L,R define", "",
        r"$$A_L=HH^\dagger+\frac{15}{4}FF^\dagger+\frac{15}{2}LL^\dagger,\qquad "
        r"A_R=H^\dagger H+\frac{15}{4}F^\dagger F+\frac{15}{2}R^\dagger R,$$", "",
        r"$$G_D=\frac94(g_L^2+g_R^2+5g_4^2),\quad G_L=\frac94(2g_L^2+5g_4^2),"
        r"\quad G_R=\frac94(2g_R^2+5g_4^2).$$", "",
        r"$$\begin{aligned}16\pi^2\dot H&=A_LH+HA_R+[4\operatorname{tr}(H^\dagger H)-G_D]H,\\"
        r"16\pi^2\dot F&=A_LF+FA_R+[\operatorname{tr}(F^\dagger F)-G_D]F,\\"
        r"16\pi^2\dot L&=A_LL+LA_L^T+[2\operatorname{tr}(L^\dagger L)-G_L]L,\\"
        r"16\pi^2\dot R&=A_R^TR+RA_R+[2\operatorname{tr}(R^\dagger R)-G_R]R.\end{aligned}$$", "",
        r"Realification of the scalar norm column yields trace coefficients $(4,1,2,2)$; the fermion columns yield "
        r"$S_{LL}=2I_8\otimes A_L^T$, $S_{RR}=2I_8\otimes A_R$. "
        r"Gauge terms follow from $C_4(4)=15/8$, $C_2(2)=3/4$ and the $-3$ anticommutator. "
        "The code compares all entries of 248 real-scalar 48×48 beta matrices with the compact system, at a generic complex unequal-L/R point.", "",
        "Deleting L reproduces the corrected "
        "[Meloni, Ohlsson and Riad, Eqs. (76)–(78)](https://arxiv.org/pdf/1612.07973). "
        "This is a regression only; P54 retains DeltaL. An independent SM-top tensor test recovers 9/2 and the "
        "GUT-normalized gauge coefficients 17/20, 9/4, 8.", "",
        "## Symmetry and a protected Spin(10) subflow", "",
        r"The equations are covariant under $H,F\mapsto U_L^T(H,F)U_R$, "
        r"$L\mapsto U_L^TLU_L$, $R\mapsto U_R^TRU_R$, and under "
        r"$H,F\mapsto H^T,F^T$, $L\leftrightarrow R$, $g_L\leftrightarrow g_R$. "
        r"So symmetric $H,F$ with $L=R,g_L=g_R$ form an invariant parity locus.", "",
        r"There is a stronger result at one loop. Impose additionally $F=\sqrt2R$. Then "
        r"$A_R^T=A_L$, $G_D=G_L=G_R$, and $\operatorname{tr}F^\dagger F=2\operatorname{tr}R^\dagger R$. "
        r"Substitution gives", "",
        r"$$\mathcal B_F-\sqrt2\mathcal B_R=0,\quad \mathcal B_L-\mathcal B_R=0,\quad"
        r"\mathcal B_H^T=\mathcal B_H,\quad\mathcal B_F^T=\mathcal B_F.$$",
        "",
        r"The gauge equations also preserve $g_L=g_R$, because $a_L=a_R$, $b_{L4}=b_{R4}$, "
        r"$b_{LL}=b_{RR}$, $b_{LR}=b_{RL}$, and $Y_{4,L}=Y_{4,R}$. "
        "Uniqueness of the smooth perturbative ODE then proves invariance of this entire submanifold. "
        "The actual tree Spin(10) tensor boundary lies on it; an explicit downward numerical flow preserves it. "
        "This is conditional one-loop protection, not a theorem about finite threshold matching or two-loop Yukawa flow.", "",
        "An unequal-gL/gR negative-control flow breaks the proportionality. With unequal Majorana parents, even initially symmetric H,F "
        "acquire antisymmetric beta sources. Therefore the general implementation retains four independently evolved matrices; "
        "it does not enforce a possibly false low-scale boundary relation.", "",
        "## Upper boundary versus lower matching", "",
        r"At the upper tree boundary use $H=\sqrt2h_{\rm raw}$, $F=4f_{\rm raw}$, $L=R=2\sqrt2f_{\rm raw}$. "
        r"At the lower scale use the evolved $H$, $f_D=F/(2\sqrt3)$ for Dirac projections; in the fixed common component phase "
        r"$M_R=+2i\,\sigma_{\rm dimful}R$, $M_L=-2i\,\Delta_{L,\rm dimful}L$. "
        r"If the scalar card writes a dimensionless $\sigma$ times $\omega$, then $\sigma_{\rm dimful}=\sigma\omega$. "
        r"The light-doublet coefficients $(a,b,d,e)$ multiply H and $F/(2\sqrt3)$; type I and II use R and L. "
        "The ratio fM/fD is preserved only on the protected locus just proved. Finite-matched off-locus data must be propagated by the full system.", "",
        "## Two-loop gauge Yukawa term", "",
        r"The actual Weyl trace $Y_{4,i}=\operatorname{Tr}[C_i\sum_aY^aY^{a\dagger}]/d(G_i)$ gives", "",
        r"$$Y_{4,4}=4\|H\|_F^2+15\|F\|_F^2+15\|L\|_F^2+15\|R\|_F^2,\quad "
        r"Y_{4,L}=4\|H\|_F^2+15\|F\|_F^2+30\|L\|_F^2,\quad "
        r"Y_{4,R}=4\|H\|_F^2+15\|F\|_F^2+30\|R\|_F^2.$$", "",
        r"$$\dot g_i=\frac{a_ig_i^3}{16\pi^2}+\frac{g_i^3}{(16\pi^2)^2}"
        r"\left[\sum_jb_{ij}g_j^2-Y_{4,i}\right],\qquad a=(2/3,26/3,26/3).$$", "",
        "The b matrix is the corrected four-parent census, with b44=3551/6, not the historical extra-six table. "
        "The Yukawa gauge trace is checked against the full tensor and independently against SM-top.", "",
        "## Executable interface", "",
        "build_ps_geometry() exports the scalar embeddings, actual group tensors and fermion basis. "
        "assemble_real_yukawas(couplings, geometry) returns 248×48×48 tensors in real-scalar order H,F,L,R, "
        "with real/imaginary coordinates interleaved. The 328×248 old-coordinate embedding and 16×16 fermion basis "
        "are included. beta_closed returns 16π² times the one-loop Yukawa beta. beta_gauge returns dg/dlog(mu), "
        "including gauge two-loop and Yukawa traces. integrate performs a coupled ODE without a fit.", "",
        "| Check | Residual | Pass |", "|---|---:|:---:|"]
    for item in result["checks"]:
        lines.append(f"| {item['name']} | {item['residual']:.3e} | {'yes' if item['passed'] else 'NO'} |")
    lines += ["", "## Remaining physical gates", "",
              "A beta function does not provide finite Yukawa threshold matching at the actual upper/lower spectra. "
              "The common action reconstruction verifies the tree boundary only. Sequential Majorana/type-II matching, "
              "the actual loop-corrected light-doublet projection and the global flavor/seesaw fit remain separate. "
              "All family entries in these ODE tests are synthetic; no physical fit, Higgs pole, determinant or portal is promoted.", ""]
    OUTPUT.with_suffix(".md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
