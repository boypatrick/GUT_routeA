#!/usr/bin/env python3
"""Moving-threshold SM+N seesaw evolution using actual local-light inputs.

One-loop MS-bar running including the 2024 dimension-five feedback correction,
the actual scalar-tree C_HN, and TREE sterile threshold matching. C_BN=0 is
retained as an RG-invariant boundary; its finite loop matching is not asserted.
Not a finite one-loop matching calculation, dimension-six truncation, or fit.
Dimensionful quantities are M/omega, omega*C5, omega*C_HN and qH=m_Z^2/omega^2.
V=+m_Z^2 HdaggerH+lambda(HdaggerH)^2. The omega reference is a
declared benchmark reference scale, not a solved P54 unification scale.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

RF = Path(__file__).resolve().parents[1]
INPUT = RF/"output/p54_self_consistent_light.json"
CHN_INPUT = RF/"output/p54_scalar_chn.json"
OUTPUT = RF/"output/p54_sequential_seesaw"
LOOP = 16*np.pi**2
MATRIX_KEYS = ("Yu","Yd","Ye","Ynu","MR","C5","C5II","CHN")
GAUGE_B = np.array([41/10,-19/6,-7.])


def cjson(value):
    a = np.asarray(value)
    return {"real":a.real.tolist(),"imag":a.imag.tolist()}


def decode(value):
    return np.asarray(value["real"])+1j*np.asarray(value["imag"])


def sym(a):
    return (a+a.T)/2


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))


def state_copy(state):
    return {key:np.array(value,copy=True) if isinstance(value,np.ndarray) else value
            for key,value in state.items()}


def takagi(matrix, *, allow_zero=False):
    """Return positive ascending masses, U with U.T M U=diag(m).

    The real symmetric 2n construction avoids phase square roots and remains
    regular in exactly degenerate positive-mass subspaces. For light-neutrino
    observables allow_zero admits a unitary kernel basis; sterile threshold
    callers reject zeros because massless states have no decoupling threshold.
    """
    n = len(matrix)
    if n == 0:
        return np.empty(0),np.empty((0,0),complex)
    scale = np.linalg.norm(matrix)
    if scale == 0:
        if allow_zero:
            return np.zeros(n),np.eye(n,dtype=complex)
        raise ValueError("Massless states do not define positive decoupling thresholds.")
    if relative(matrix,matrix.T)>2e-10:
        raise ValueError("Takagi input must be complex symmetric.")
    m = matrix/scale
    real = np.block([[m.real,-m.imag],[-m.imag,-m.real]])
    values,vectors = np.linalg.eigh(real)
    values,vectors = values[-n:],vectors[:,-n:]
    positive = values>2e-13
    if not all(positive) and not allow_zero:
        raise ValueError("Rank-deficient Takagi block: keep genuinely massless sterile states active.")
    u = vectors[:n,positive]+1j*vectors[n:,positive]
    masses = values[positive]*scale
    if not all(positive):
        zero_count = n-sum(positive)
        _,_,vh = np.linalg.svd(m)
        kernel = vh.conj().T[:,-zero_count:]
        kernel -= u@(u.conj().T@kernel)
        kernel = np.linalg.qr(kernel)[0]
        u = np.column_stack((kernel,u))
        masses = np.r_[np.zeros(zero_count),masses]
    if relative(u.conj().T@u,np.eye(n))>3e-10 or relative(u.T@matrix@u,np.diag(masses))>3e-10:
        raise ValueError("Takagi congruence failed.")
    return masses,u


def beta(state, *, include_dim5_feedback=True):
    """One-loop d/dlog(mu), gauge order (GUT-normalized g1,g2,g3).

    All Dirac Yukawas have doublet-family rows and singlet-family columns.
    V(H)=+qH*omega^2 HdaggerH + lambda (HdaggerH)^2. C5 and C5II have mass
    convention Mnu=v^2/2*(C5-Ynu MR^-1 Ynu.T), before omega rescaling.
    Default includes Zhang 2024 Eq.(8), with C5_here=-C5_Zhang.
    False removes the three Weinberg feedback corrections as a negative control,
    while retaining the C_HN and Higgs-mass terms in both systems.
    """
    u,d,e,n = (state[key] for key in ("Yu","Yd","Ye","Ynu"))
    hu,hd,he,hn = (a@a.conj().T for a in (u,d,e,n))
    tr = lambda a: float(np.trace(a).real)
    trace = tr(3*hu+3*hd+he+hn)
    quartic_trace = tr(3*hu@hu+3*hd@hd+he@he+hn@hn)
    g1,g2,g3 = state["g"]
    lam = float(state["lambda"])
    eye = np.eye(3)
    out = {
        "Yu": (1.5*(hu-hd)+(trace-17/20*g1*g1-9/4*g2*g2-8*g3*g3)*eye)@u,
        "Yd": (1.5*(hd-hu)+(trace-1/4*g1*g1-9/4*g2*g2-8*g3*g3)*eye)@d,
        "Ye": (1.5*(he-hn)+(trace-9/4*g1*g1-9/4*g2*g2)*eye)@e,
        "Ynu":(1.5*(hn-he)+(trace-9/20*g1*g1-9/4*g2*g2)*eye)@n,
    }
    sn = n.conj().T@n
    m,chn,qh = state["MR"],state["CHN"],state["qH"]
    alpha_h = 12*lam+2*trace-.9*g1*g1-4.5*g2*g2
    out["MR"] = sn.T@m+m@sn-8*qh*chn
    out["CHN"] = alpha_h*chn+4*(chn@sn+sn.T@chn)
    out["qH"] = (alpha_h*qh-4*np.trace(m.conj().T@m@sn).real
                 +8*np.trace(m.conj().T@m@m.conj().T@chn).real)
    out["Ynu"] -= 4*n@m.conj().T@chn
    p = -1.5*he+(.5+3*include_dim5_feedback)*hn
    alpha = 4*lam-3*g2*g2+2*trace
    for key in ("C5","C5II"):
        out[key] = p@state[key]+state[key]@p.T+alpha*state[key]
    out["g"] = GAUGE_B*state["g"]**3
    out["lambda"] = (24*lam**2-3*lam*(3*g2*g2+3/5*g1*g1)
                     +3/4*g2**4+3/8*(3/5*g1*g1+g2*g2)**2
                     +4*lam*trace-2*quartic_trace
                     +16*np.trace(chn.conj().T@m@sn).real)
    if include_dim5_feedback:
        out["Ynu"] += 3*state["C5"]@n.conj()@state["MR"]
        out["lambda"] -= 4*np.trace(state["C5"]@n.conj()@state["MR"]@n.conj().T).real
    return {key:value/LOOP for key,value in out.items()}


def effective_c5(state):
    n = state["Ynu"]
    return state["C5"]-n@np.linalg.solve(state["MR"],n.T) if n.shape[1] else state["C5"].copy()


def decouple_block(state, heavy_indices, sterile_basis=None):
    """Exact algebraic block-tree Schur matching; dimension-six terms omitted.

    U is unitary and acts as N_old=U N_new, hence Ynew=Yold U,
    Mnew=U.T Mold U. A non-diagonal retained/heavy mass block is allowed.
    """
    out = state_copy(state)
    count = state["Ynu"].shape[1]
    u = np.eye(count,dtype=complex) if sterile_basis is None else sterile_basis
    if relative(u.conj().T@u,np.eye(count))>2e-10:
        raise ValueError("Threshold change of sterile basis must be unitary.")
    h = np.asarray(heavy_indices,dtype=int)
    if len(set(h.tolist()))!=len(h) or np.any(h<0) or np.any(h>=count) or not len(h):
        raise ValueError("Heavy block must contain distinct valid indices.")
    r = np.asarray([j for j in range(count) if j not in h],dtype=int)
    y = state["Ynu"]@u
    m = u.T@state["MR"]@u
    a,b,d = m[np.ix_(r,r)],m[np.ix_(r,h)],m[np.ix_(h,h)]
    yh,yr = y[:,h],y[:,r]
    out["C5"] = sym(state["C5"]-yh@np.linalg.solve(d,yh.T))
    out["Ynu"] = yr-yh@np.linalg.solve(d,b.T)
    out["MR"] = sym(a-b@np.linalg.solve(d,b.T))
    # M(rho)=M-2*C_HN*rho: differentiate its Schur complement at rho=0.
    # This is T.T C_HN T, not just its retained block unless b=0.
    pullback = np.zeros((count,len(r)),complex)
    pullback[r,:] = np.eye(len(r))
    pullback[h,:] = -np.linalg.solve(d,b.T)
    rotated_chn = u.T@state["CHN"]@u
    out["CHN"] = sym(pullback.T@rotated_chn@pullback)
    record = {"removed":len(h),"retained":len(r),
              "tree_C5_increment":cjson(out["C5"]-state["C5"]),
              "tree_effective_C5_continuity":relative(effective_c5(out),effective_c5(state)),
              "offdiagonal_mass_block_norm":float(np.linalg.norm(b)),
              "retained_yukawa_schur_shift_norm":float(np.linalg.norm(out["Ynu"]-yr)),
              "CHN_before_times_omega":cjson(state["CHN"]),
              "CHN_after_times_omega":cjson(out["CHN"]),
              "CHN_pullback_correction_norm":float(np.linalg.norm(out["CHN"]-rotated_chn[np.ix_(r,r)]))}
    return out,record,r


def packer(n):
    shapes = {"Yu":(3,3),"Yd":(3,3),"Ye":(3,3),"Ynu":(3,n),"MR":(n,n),
              "C5":(3,3),"C5II":(3,3),"CHN":(n,n)}
    lengths = [np.prod(shapes[key],dtype=int) for key in MATRIX_KEYS]
    size = sum(lengths)
    def pack(state):
        z = np.concatenate([state[key].ravel() for key in MATRIX_KEYS])
        return np.r_[z.real,z.imag,state["g"],state["lambda"],state["qH"]]
    def unpack(values):
        z = values[:size]+1j*values[size:2*size]
        result = {}
        start = 0
        for key,length in zip(MATRIX_KEYS,lengths):
            result[key] = z[start:start+length].reshape(shapes[key])
            start += length
        result["g"],result["lambda"],result["qH"] = values[-5:-2],float(values[-2]),float(values[-1])
        return result
    return pack,unpack


def evolve(state, mu_start, mu_stop, *, cluster_log_tolerance=1e-7,
           rtol=3e-10, atol=2e-12, max_step=.30, include_dim5_feedback=True,
           threshold_matcher=None, threshold_scale_factor=1.):
    """Run downward in omega units, solving mu=max singular(MR(mu)) live.

    Ordering is recomputed from the running matrix at every event. Eigenvector
    derivatives and labels fixed at the initial mass spectrum are never used.
    Degenerate/near-degenerate threshold clusters are removed in one block.
    """
    if not 0<mu_stop<mu_start or threshold_scale_factor<=0:
        raise ValueError("Need 0 < mu_stop < mu_start.")
    state = state_copy(state)
    n0 = state["Ynu"].shape[1]
    initial_masses = np.linalg.svd(state["MR"],compute_uv=False) if n0 else np.empty(0)
    if n0 and mu_start<=threshold_scale_factor*initial_masses[0]*(1+1e-8):
        raise ValueError("Initial scale must be above all retained running Majorana masses.")
    t,tstop = math.log(mu_start),math.log(mu_stop)
    stages,events = [],[]
    sterile_embedding = np.eye(n0,dtype=complex)
    while t>tstop+1e-13:
        n = state["Ynu"].shape[1]
        pack,unpack = packer(n)
        begin_state = state_copy(state)
        def rhs(time,values):
            return pack(beta(unpack(values),include_dim5_feedback=include_dim5_feedback))
        def threshold(time,values):
            matrix = unpack(values)["MR"]
            highest = np.linalg.svd(matrix,compute_uv=False)[0]
            return time-math.log(threshold_scale_factor*highest)
        threshold.terminal = True
        threshold.direction = -1
        solution = solve_ivp(rhs,(t,tstop),pack(state),method="DOP853",rtol=rtol,atol=atol,
                             events=threshold if n else None,max_step=max_step)
        if not solution.success:
            raise RuntimeError(solution.message)
        state = unpack(solution.y[:,-1])
        new_t = float(solution.t[-1])
        stages.append({"active_N":n,"mu_start_over_omega":math.exp(t),"mu_end_over_omega":math.exp(new_t),
                       "nfev":solution.nfev,
                       "MR_start_over_omega":cjson(begin_state["MR"]),
                       "MR_end_over_omega":cjson(state["MR"]),
                       "MR_end_offdiagonal_norm":float(np.linalg.norm(
                           state["MR"]-np.diag(np.diag(state["MR"])))) if n else 0.,
                       "running_mass_start":np.linalg.svd(begin_state["MR"],compute_uv=False).tolist() if n else [],
                       "running_mass_end":np.linalg.svd(state["MR"],compute_uv=False).tolist() if n else [],
                       "C5_change_norm":float(np.linalg.norm(state["C5"]-begin_state["C5"])),
                       "CHN_start_times_omega":cjson(begin_state["CHN"]),
                       "CHN_end_times_omega":cjson(state["CHN"]),
                       "CHN_change_norm":float(np.linalg.norm(state["CHN"]-begin_state["CHN"])),
                       "qH_start":float(begin_state["qH"]),"qH_end":float(state["qH"]),
                       "all_Yukawa_change_norms":{key:float(np.linalg.norm(state[key]-begin_state[key]))
                                                 for key in ("Yu","Yd","Ye","Ynu")}})
        t = new_t
        if not n or not len(solution.t_events[0]):
            break
        masses,u = takagi(state["MR"])
        maximum = masses[-1]
        heavy = np.flatnonzero(abs(np.log(masses/maximum))<=cluster_log_tolerance)
        heavy_embedding = sterile_embedding@u[:,heavy]
        old_y = state["Ynu"].copy()
        if threshold_matcher is None:
            state,record,retained = decouple_block(state,heavy,u)
        else:
            state,record,retained = threshold_matcher(state,heavy,u,math.exp(t))
        record.update({"mu_over_omega":math.exp(t),"running_masses_at_event":masses.tolist(),
                       "root_log_residual":t-math.log(threshold_scale_factor*maximum),"cluster_log_tolerance":cluster_log_tolerance,
                       "cluster_spread_log":float(np.ptp(np.log(masses[heavy]))),
                       "removed_projector_in_initial_basis":cjson(heavy_embedding@heavy_embedding.conj().T),
                       "removed_yukawa_matrix":cjson(old_y@u[:,heavy])})
        events.append(record)
        sterile_embedding = sterile_embedding@u[:,retained]
    return state,{"stages":stages,"events":events,"initial_running_masses":initial_masses.tolist(),
                  "initial_mu_over_omega":mu_start,"final_mu_over_omega":math.exp(t),
                  "remaining_N":state["Ynu"].shape[1],
                  "matching_order":"tree" if threshold_matcher is None else "explicit callback; see per-event completeness ledger",
                  "threshold_scale_factor":threshold_scale_factor,
                  "running_order":"one-loop through dimension five, nonzero C_HN, C_BN=0",
                  "Zhang_2024_dimension5_feedback_included":include_dim5_feedback}


def family_transform(state, uq, uu, ud, ul, ue, un):
    out = state_copy(state)
    for key,left,right in (("Yu",uq,uu),("Yd",uq,ud),("Ye",ul,ue),("Ynu",ul,un)):
        out[key] = left.T@state[key]@right
    for key in ("MR","CHN"):
        out[key] = un.T@state[key]@un
    for key in ("C5","C5II"):
        out[key] = ul.T@state[key]@ul
    return out


def unitary(rng,n):
    q,r = np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    return q*np.exp(-1j*np.angle(np.diag(r)))[None,:]


def terminal_observables(state, omega_GeV, v_GeV=246.22):
    """MS-bar mass proxies at fixed reference v, not pole observables."""
    c = sym(effective_c5(state))
    masses,unu = takagi(c,allow_zero=True)
    hu = state["Yu"]@state["Yu"].conj().T
    hd = state["Yd"]@state["Yd"].conj().T
    he = state["Ye"]@state["Ye"].conj().T
    eu,vu = np.linalg.eigh(hu)
    ed,vd = np.linalg.eigh(hd)
    ee,ve = np.linalg.eigh(he)
    pmns = ve.T@unu
    ckm = vu.T@vd.conj()
    mod = abs(pmns)
    s13 = mod[0,2]
    s12 = mod[0,1]/math.sqrt(max(1-s13*s13,1e-30))
    s23 = mod[1,2]/math.sqrt(max(1-s13*s13,1e-30))
    angles = np.degrees(np.arcsin(np.clip([s12,s13,s23],0,1)))
    neutrino_eV = masses*v_GeV**2/(2*omega_GeV)*1e9
    return {"fixed_v_GeV":v_GeV,"benchmark_omega_GeV":omega_GeV,
            "physical_endpoint":False,
            "neutrino_masses_eV":neutrino_eV.tolist(),
            "delta_m21_sq_eV2":float(neutrino_eV[1]**2-neutrino_eV[0]**2),
            "delta_m31_sq_eV2":float(neutrino_eV[2]**2-neutrino_eV[0]**2),
            "PMNS_moduli":mod.tolist(),"angles_theta12_theta13_theta23_deg":angles.tolist(),
            "PMNS_J":float(np.imag(pmns[0,0]*pmns[1,1]*pmns[0,1].conj()*pmns[1,0].conj())),
            "CKM_moduli":abs(ckm).tolist(),
            "CKM_J":float(np.imag(ckm[0,0]*ckm[1,1]*ckm[0,1].conj()*ckm[1,0].conj())),
            "charged_running_mass_proxies_GeV":{"u":(np.sqrt(np.maximum(eu,0))*v_GeV/math.sqrt(2)).tolist(),
                                               "d":(np.sqrt(np.maximum(ed,0))*v_GeV/math.sqrt(2)).tolist(),
                                               "e":(np.sqrt(np.maximum(ee,0))*v_GeV/math.sqrt(2)).tolist()},
            "terminal_gauge_g1_g2_g3":state["g"].tolist(),"terminal_lambda":state["lambda"],
            "terminal_qH_mZ2_over_omega2":state["qH"],
            "terminal_Higgs_mZ2_GeV2_diagnostic":state["qH"]*omega_GeV**2,
            "C5_total_times_omega":cjson(state["C5"]),"C5_typeII_times_omega":cjson(state["C5II"]),
            "Mnu_GeV":cjson(c*v_GeV**2/(2*omega_GeV))}


def load_local_light(case, chn_matrix=None):
    """Use exported actual local-light matrices; never invent four UV inputs."""
    row = case["all_four_Yukawa_matrices_at_new_c"]
    out = {key:decode(row[species]) for key,species in (("Yu","u"),("Yd","d"),("Ye","e"),("Ynu","nu"))}
    out["MR"] = decode(case["MR_over_omega"])
    out["C5"] = 2*decode(case["ML_over_v2_over_omega"])
    out["C5II"] = out["C5"].copy()
    out["CHN"] = np.zeros_like(out["MR"]) if chn_matrix is None else np.array(chn_matrix,copy=True)
    out["qH"] = 0.
    out["g"],out["lambda"] = np.array([.53,.55,.60]),.12
    return out


def serialize_state(state):
    return {key:cjson(state[key]) for key in MATRIX_KEYS}|{
        "g":state["g"].tolist(),"lambda":float(state["lambda"]),"qH":float(state["qH"])}


def run():
    source = json.loads(INPUT.read_text())
    chn_source = json.loads(CHN_INPUT.read_text())
    if not source["summary"]["all_pass"]:
        raise ValueError("Actual local-light prerequisite has failing checks.")
    if not chn_source["all_checks_pass"]:
        raise ValueError("Actual scalar C_HN prerequisite has failing checks.")
    checks = []
    def check(name,residual,tolerance=2e-8):
        checks.append({"name":name,"residual":float(residual),"tolerance":tolerance,
                       "passed":bool(residual<tolerance)})
    rng = np.random.default_rng(541031)
    cases = []
    omega = 1e14
    stop = 246.22/omega
    for index,case in enumerate(source["cases"]):
        chn_row = chn_source["rays"][2+index]
        state = load_local_light(case,decode(chn_row["CHN_matrix_times_omega"]))
        check(f"case {index}: actual local I+II normalization reconstruction",relative(
            effective_c5(state)/2,decode(case["Mnu_over_v2_over_omega"])),2e-12)
        final,history = evolve(state,.1,stop)
        higgs_ratio = float(abs(final["qH"])/stop**2)
        applicability = {
            "terminal_mu_over_omega":float(stop),
            "abs_qH_over_terminal_mu_squared":higgs_ratio,
            "terminal_light_Higgs_EFT_condition":bool(higgs_ratio < 1.),
            "physical_endpoint":False,
            "reason":"Unmatched quadratic boundary; mass-independent MS-bar diagnostic continuation, not a physical light-Higgs EFT endpoint."}
        check(f"case {index}: actual scalar-tree CHN is nonzero and runs",0 if
              np.linalg.norm(state["CHN"])>1e-4 and history["stages"][0]["CHN_change_norm"]>1e-7 else 1)
        check(f"case {index}: Higgs qH generated from diagnostic zero boundary",0 if
              state["qH"]==0 and abs(final["qH"])>1e-9
              and not applicability["terminal_light_Higgs_EFT_condition"] else 1)
        rotations = [unitary(rng,3) for _ in range(6)]
        transformed = family_transform(state,*rotations)
        rotated_final,rotated_history = evolve(transformed,.1,stop)
        expected = family_transform(final,*rotations[:5],np.empty((0,0),complex))
        check(f"case {index}: full running plus moving thresholds family covariance",
              max([relative(rotated_final[key],expected[key]) for key in ("Yu","Yd","Ye","C5","C5II")]
                  +[abs(rotated_final["qH"]-expected["qH"]),abs(rotated_final["lambda"]-expected["lambda"])]))
        check(f"case {index}: event scales are sterile-basis independent",relative(
            np.array([row["mu_over_omega"] for row in history["events"]]),
            np.array([row["mu_over_omega"] for row in rotated_history["events"]])))
        check(f"case {index}: all three moving thresholds and below-last running",0 if
              history["remaining_N"]==0 and len(history["events"])==3 and history["stages"][-1]["active_N"]==0 else 1)
        check(f"case {index}: live threshold roots",max(abs(row["root_log_residual"]) for row in history["events"]),2e-10)
        check(f"case {index}: exact tree continuity at every threshold",
              max(row["tree_effective_C5_continuity"] for row in history["events"]),2e-11)
        initial_m = np.asarray(history["initial_running_masses"])
        event_scales = np.array([row["mu_over_omega"] for row in history["events"]])
        drift = (event_scales-initial_m)/initial_m
        check(f"case {index}: thresholds are not frozen initial masses",0 if max(abs(drift))>1e-4 else 1)
        observables = terminal_observables(final,omega)
        rotated_observables = terminal_observables(rotated_final,omega)
        check(f"case {index}: terminal PMNS and neutrino mass covariance",max(
            relative(np.asarray(observables[key]),np.asarray(rotated_observables[key]))
            for key in ("PMNS_moduli","neutrino_masses_eV")))
        check(f"case {index}: active and below-last stages evolve nontrivially",0 if all(
            row["all_Yukawa_change_norms"]["Yu"]>1e-5 for row in history["stages"]) else 1)
        legacy,legacy_history = evolve(state,.1,stop,include_dim5_feedback=False)
        legacy_observables = terminal_observables(legacy,omega)
        legacy_difference = {"scope":"omit Weinberg feedback only; CHN/qH still active; negative control, not prediction",
                             "C5_relative_difference":relative(final["C5"],legacy["C5"]),
                             "lambda_difference":float(final["lambda"]-legacy["lambda"]),
                             "threshold_relative_difference":(
                                 np.array([x["mu_over_omega"] for x in history["events"]])/
                                 np.array([x["mu_over_omega"] for x in legacy_history["events"]])-1).tolist(),
                             "neutrino_masses_eV":legacy_observables["neutrino_masses_eV"]}
        without_chn = state_copy(state)
        without_chn["CHN"] *= 0
        no_chn_final,no_chn_history = evolve(without_chn,.1,stop)
        no_chn_comparison = {"C5_relative_difference":relative(final["C5"],no_chn_final["C5"]),
                             "lambda_difference":float(final["lambda"]-no_chn_final["lambda"]),
                             "qH_difference":float(final["qH"]-no_chn_final["qH"]),
                             "event_scales_without_CHN":[r["mu_over_omega"] for r in no_chn_history["events"]]}
        check(f"case {index}: actual nonzero CHN changes the full trajectory",0 if
              no_chn_comparison["C5_relative_difference"]>1e-6 else 1)
        cases.append({"input_case":index,"input_state":serialize_state(state),"history":history,
                      "CHN_source":chn_row["name"],
                      "CHN_matching_order_notice":"actual corrected external light ray, frozen tree scalar vertices/Hessian; mixed-order projection",
                      "relative_threshold_drift_from_initial_descending_masses":drift.tolist(),
                      "terminal_state":serialize_state(final),"terminal_observables":observables,
                      "endpoint_applicability":applicability,
                      "without_scalar_CHN_comparison":no_chn_comparison,
                      "legacy_without_dim5_feedback_comparison":legacy_difference})
    # General non-diagonal tree Schur elimination agrees with one-shot.
    state = load_local_light(source["cases"][0])
    # Generic complex symmetric CHN tests both covariance and non-diagonal pullback.
    state["CHN"] = decode(chn_source["rays"][2]["CHN_matrix_times_omega"])
    one_shot = effective_c5(state)
    first,record,_ = decouple_block(state,[1])
    second,record2,_ = decouple_block(first,[0])
    third,record3,_ = decouple_block(second,[0])
    check("non-diagonal exact block-tree sequential equals one-shot Schur",relative(third["C5"],one_shot),2e-12)
    check("non-diagonal Schur test exercises retained-Yukawa shift",0 if
          record["retained_yukawa_schur_shift_norm"]>1e-3 else 1)
    # Differentiate the full field-dependent Schur complement independently.
    def schur_at_rho(rho):
        m = state["MR"]-2*rho*state["CHN"]
        r,h = [0,2],[1]
        return m[np.ix_(r,r)]-m[np.ix_(r,h)]@np.linalg.solve(m[np.ix_(h,h)],m[np.ix_(h,r)])
    step = 1e-4
    fd_chn = -(schur_at_rho(step)-schur_at_rho(-step))/(4*step)
    check("CHN block pullback equals field-dependent Schur derivative",relative(first["CHN"],fd_chn),2e-9)
    check("non-diagonal CHN matching is not naive corner truncation",0 if
          record["CHN_pullback_correction_norm"]>1e-5 else 1)
    # Two successive partial eliminations leave one sterile state.
    combined,_,_ = decouple_block(state,[0,1])
    check("CHN sequential pullback equals one-shot retained Schur derivative",
          relative(second["CHN"],combined["CHN"]),2e-12)
    # Degenerate Takagi subspace: no preferred eigenvectors or drop order.
    deg = state_copy(state)
    deg.update({"Yu":.12*np.eye(3,dtype=complex),"Yd":.09*np.eye(3,dtype=complex),
                "Ye":.07*np.eye(3,dtype=complex),"Ynu":.30*np.eye(3,dtype=complex),
                "MR":.035*np.eye(3,dtype=complex),"C5":.003j*np.eye(3,dtype=complex),
                "C5II":.003j*np.eye(3,dtype=complex),"CHN":.004j*np.eye(3,dtype=complex)})
    deg_final,deg_history = evolve(deg,.1,stop)
    rotations = [unitary(rng,3) for _ in range(6)]
    deg_rot = family_transform(deg,*rotations)
    dr_final,dr_history = evolve(deg_rot,.1,stop)
    expected = family_transform(deg_final,*rotations[:5],np.empty((0,0),complex))
    check("exact degenerate Takagi cluster removed as one block",0 if
          len(deg_history["events"])==1 and deg_history["events"][0]["removed"]==3 else 1)
    check("degenerate threshold and final C5 are fully family covariant",max(
        relative(dr_final["C5"],expected["C5"]),
        abs(deg_history["events"][0]["mu_over_omega"]-dr_history["events"][0]["mu_over_omega"])))
    masses,u = takagi(deg_rot["MR"])
    q,_ = np.linalg.qr(rng.normal(size=(3,3)))
    b1,_,_ = decouple_block(deg_rot,[0,1,2],u)
    b2,_,_ = decouple_block(deg_rot,[0,1,2],u@q)
    check("degenerate Takagi O(3) basis choice leaves block matching invariant",relative(b1["C5"],b2["C5"]),2e-12)
    # Level crossing: initially larger M_2 runs down fastest; N_1 leaves first.
    crossing = state_copy(state)
    crossing["Yu"] = .12*np.eye(3,dtype=complex)
    crossing["Yd"] = .09*np.eye(3,dtype=complex)
    crossing["Ye"] = .07*np.eye(3,dtype=complex)
    crossing["C5"] = .001*np.eye(3,dtype=complex)
    crossing["C5II"] = .001*np.eye(3,dtype=complex)
    crossing["MR"] = np.diag([.049,.05]).astype(complex)
    crossing["CHN"] = np.zeros((2,2),complex)
    crossing["Ynu"] = np.array([[.03,0],[0,1.2],[0,0]],complex)
    crossing_final,crossing_history = evolve(crossing,1.,.005)
    projector = decode(crossing_history["events"][0]["removed_projector_in_initial_basis"])
    check("live event ordering detects a running sterile-mass level crossing",
          relative(projector,np.diag([1.,0.])),2e-9)
    # Differential covariance before any matching is separately checked.
    rotations = [unitary(rng,3) for _ in range(6)]
    changed = family_transform(state,*rotations)
    expected_beta = family_transform(beta(state),*rotations)
    check("one-loop full matrix beta family covariance",max(
        [relative(beta(changed)[key],expected_beta[key]) for key in MATRIX_KEYS]
        +[abs(beta(changed)["qH"]-expected_beta["qH"]),
          abs(beta(changed)["lambda"]-expected_beta["lambda"])]),2e-12)
    # Directly isolate all CHN feedback signs, with a nonzero signed Higgs mass.
    probe = state_copy(state)
    probe["qH"] = -.003
    empty_chn = state_copy(probe)
    empty_chn["CHN"] *= 0
    bp,b0 = beta(probe),beta(empty_chn)
    m,c,y = probe["MR"],probe["CHN"],probe["Ynu"]
    sn = y.conj().T@y
    check("CHN feedback uses Zhang positive-potential-mass signs",max(
        relative(bp["MR"]-b0["MR"],-8*probe["qH"]*c/LOOP),
        relative(bp["Ynu"]-b0["Ynu"],-4*y@m.conj().T@c/LOOP),
        abs((bp["lambda"]-b0["lambda"])-16*np.trace(c.conj().T@m@sn).real/LOOP),
        abs((bp["qH"]-b0["qH"])-8*np.trace(m.conj().T@m@m.conj().T@c).real/LOOP)),2e-11)
    check("MR and CHN one-loop derivatives preserve complex symmetry",max(
        relative(bp[key],bp[key].T) for key in ("MR","CHN")),2e-12)
    check("zero CHN remains RG-invariant but is not the actual scalar boundary",
          np.linalg.norm(b0["CHN"]),2e-12)
    # Independent one-family CW-divergence polynomial, not an RGE-to-itself test.
    cw_y,cw_m,cw_c = .37,.63,.04
    h_grid = np.linspace(-.8,.8,25)
    def cw_polynomial(c):
        values = []
        for hv in h_grid:
            mf = np.array([[0.,cw_y*hv/np.sqrt(2)],
                           [cw_y*hv/np.sqrt(2),cw_m-c*hv*hv]])
            x = mf.T@mf
            values.append(np.trace(x@x))
        return np.polynomial.polynomial.polyfit(h_grid,values,8)
    baseline = cw_polynomial(0.)
    odd = (cw_polynomial(cw_c)-cw_polynomial(-cw_c))/2
    check("independent CW polynomial fixes CHN Higgs mass/quartic feedback signs",max(
        abs(baseline[2]-2*cw_m**2*cw_y**2),
        abs(baseline[4]-.5*cw_y**4),
        abs(-2*odd[2]-8*cw_m**3*cw_c),
        abs(-4*odd[4]-16*cw_y**2*cw_m*cw_c)),2e-11)
    cw_diagnostic = {"y":cw_y,"M":cw_m,"CHN":cw_c,
                     "trace_MdaggerM_squared_baseline_h2_h4":baseline[[2,4]].tolist(),
                     "odd_in_CHN_h2_h4":odd[[2,4]].tolist(),
                     "inferred_beta_mZ2_CHN":float(-2*odd[2]),
                     "inferred_beta_lambda_CHN":float(-4*odd[4])}
    # Canonical-subsector composite identities are deliberately tested at CHN=0.
    state["CHN"] *= 0
    # Composite seesaw identity separates C5 and active-N alpha terms.
    b = beta(state)
    y,m = state["Ynu"],state["MR"]
    q = y@np.linalg.solve(m,y.T)
    dq = (b["Ynu"]@np.linalg.solve(m,y.T)+y@np.linalg.solve(m,b["Ynu"].T)
          -y@np.linalg.solve(m,b["MR"]@np.linalg.solve(m,y.T)))
    h = {key:state[key]@state[key].conj().T for key in ("Yu","Yd","Ye","Ynu")}
    tr = np.trace(3*h["Yu"]+3*h["Yd"]+h["Ye"]+h["Ynu"]).real
    p = -1.5*h["Ye"]+.5*h["Ynu"]
    g1,g2,_ = state["g"]
    expected_dq = (p@q+q@p.T+(2*tr-.9*g1*g1-4.5*g2*g2)*q
                   +3*(h["Ynu"]@state["C5"]+state["C5"]@h["Ynu"].T))/LOOP
    check("active seesaw composite beta includes corrected dimension-five feedback",relative(dq,expected_dq),2e-12)
    total = effective_c5(state)
    alpha_c = 4*state["lambda"]-3*g2*g2+2*tr
    alpha_gap = 4*state["lambda"]+.9*g1*g1+1.5*g2*g2
    expected_total = (p@total+total@p.T+alpha_c*total+alpha_gap*q)/LOOP
    check("I-plus-II coefficient has the required active-seesaw source term",
          relative(b["C5"]-dq,expected_total),2e-12)
    missing_source = float(np.linalg.norm(alpha_gap*q/LOOP))
    check("one Weinberg coefficient alone is not closed above all thresholds",
          0 if missing_source>1e-5 else 1)
    old = beta(state,include_dim5_feedback=False)
    correction_y = 3*state["C5"]@state["Ynu"].conj()@state["MR"]/LOOP
    correction_c = 3*(h["Ynu"]@state["C5"]+state["C5"]@h["Ynu"].T)/LOOP
    correction_lambda = -4*np.trace(state["C5"]@state["Ynu"].conj()@
                                    state["MR"]@state["Ynu"].conj().T).real/LOOP
    check("three Zhang 2024 corrections have the declared C5 sign",max(
        relative(b["Ynu"]-old["Ynu"],correction_y),
        relative(b["C5"]-old["C5"],correction_c),
        abs(b["lambda"]-old["lambda"]-correction_lambda)),2e-11)
    old_dq = (old["Ynu"]@np.linalg.solve(m,y.T)+y@np.linalg.solve(m,old["Ynu"].T)
              -y@np.linalg.solve(m,old["MR"]@np.linalg.solve(m,y.T)))
    check("new dimension-five terms cancel instantaneously in beta(C5-Q)",
          relative(b["C5"]-dq,old["C5"]-old_dq),2e-12)
    # Below the last threshold, charged-lepton-free C5 scales multiplicatively.
    zero = state_copy(state)
    zero.update({"Yu":np.zeros((3,3),complex),"Yd":np.zeros((3,3),complex),
                 "Ye":np.zeros((3,3),complex),"Ynu":np.zeros((3,0),complex),
                 "MR":np.empty((0,0),complex),"CHN":np.empty((0,0),complex)})
    alpha = (4*zero["lambda"]-3*zero["g"][1]**2)/LOOP
    check("canonical Higgs quartic gives +4lambda in below-threshold C5 beta",
          relative(beta(zero)["C5"],alpha*zero["C5"]),2e-12)
    minimal = state_copy(state)
    minimal["Ynu"] = minimal["Ynu"][:,:2]
    minimal["MR"] = np.diag([.025,.055]).astype(complex)
    minimal["CHN"] = np.zeros((2,2),complex)
    minimal["C5"] = np.zeros((3,3),complex)
    minimal["C5II"] = np.zeros((3,3),complex)
    minimal_final,minimal_history = evolve(minimal,.1,stop)
    singular = np.linalg.svd(minimal_final["C5"],compute_uv=False)
    check("minimal two-sterile model preserves a massless neutrino at one loop",
          singular[-1]/singular[0],2e-10)
    minimal_observables = terminal_observables(minimal_final,omega)
    check("terminal Takagi interface admits the physical rank-two light spectrum",
          minimal_observables["neutrino_masses_eV"][0],1e-12)
    result = {"schema":"p54-live-sequential-seesaw-v2-CHN","date":"2026-09-06",
              "all_checks_pass":all(row["passed"] for row in checks),
              "checks_passed":sum(row["passed"] for row in checks),"checks_total":len(checks),"checks":checks,
              "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),INPUT,CHN_INPUT)},
              "scope":"one-loop SM+N+C5+CHN+qH including all 2024 dim5 feedback at C_BN=0; exact block TREE matching; not a fit",
              "benchmark_scale_notice":"omega=1e14 GeV and initial g/lambda are illustrative boundary data, not a P54 scale solution",
              "cases":cases,"degenerate_cluster_history":deg_history,"level_crossing_history":crossing_history,
              "active_seesaw_alpha_C5_minus_alpha_Q":alpha_gap,
              "naive_single_C5_missing_source_norm":missing_source,
              "minimal_two_sterile_rank_test":{"smallest_over_largest_C5_singular_value":float(singular[-1]/singular[0]),
                                               "history":minimal_history,"terminal_observables":minimal_observables},
              "tree_Schur_records":[record,record2,record3],
              "independent_one_family_CW_feedback_sign_check":cw_diagnostic,
              "finite_one_loop_matching_implemented":False,
              "limitations":["Upper/lower finite Yukawa and scalar thresholds remain open.",
                             "Actual scalar-tree CHN is nonzero; corrected external rays use frozen tree vertices/Hessian (mixed order).",
                             "CBN=0 is one-loop RG invariant and has no scalar-tree source; finite dipole matching remains open.",
                             "qH=m_Zhang^2/omega^2 starts at zero as a diagnostic, not an electroweak-matched Higgs mass.",
                             "The PQ axion is retained as a spectator; its loop EFT is not included in this nuSMEFT subsystem.",
                             "No dimension-six nuSMEFT basis, double-dim5 insertions, or its feedback is included.",
                             "Tree Schur matching uses leading 1/M seesaw EFT, not exact finite-v pole diagonalization.",
                             "No two-loop running or one-loop finite type-I C5 matching.",
                             "Fixed-v MS-bar proxies are not physical pole observables.",
                             "Synthetic family inputs are constrained by actual local-light Clebsches, but are not fitted."]}
    OUTPUT.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    write_report(result)
    print(json.dumps({"all_checks_pass":result["all_checks_pass"],"checks_passed":result["checks_passed"],
                      "checks_total":len(checks),"failed":[row for row in checks if not row["passed"]]},indent=2))
    if not result["all_checks_pass"]:
        raise SystemExit(1)
    return result


def write_report(result):
    lines = [
        "# P54: moving-threshold sequential seesaw", "",
        "Date: 2026-09-06.", "",
        f"Status: {result['checks_passed']}/{result['checks_total']} checks pass. "
        "An operational SM+N matrix flow now solves the running mass thresholds, performs exact block-tree matching, "
        "and continues below the last sterile neutrino. The actual nonzero scalar-tree CHN now enters the complete "
        "one-loop dimension-five running subsystem with CBN=0 and a running Higgs mass parameter, plus tree matching, "
        "not a completed finite one-loop EFT calculation or a physical flavor fit.", "",
        "## 1. Input contract and conventions", "",
        "Both runs read the existing local-light JSON. Yu,Yd,Ye,Ynu are the four matrices evaluated on its actual "
        "self-consistently corrected scalar light vector. MR/omega is its common-phase Majorana matrix. "
        "The two synthetic h_raw,f_raw matrices and actual Clifford projections determine all these inputs; "
        "the code does not invent four free UV family matrices. CHN is read from p54_scalar_chn.json, evaluated "
        "on the same corrected external light ray with the frozen tree scalar Hessian and vertices. "
        "It is correlated with the same f_raw, not a new family matrix. This external-ray improvement is "
        "explicitly mixed order, not a completed loop-level scalar matching calculation.", "",
        r"All Dirac matrices have doublet-family rows and singlet-family columns. With "
        r"$H^0=v/\sqrt2$, define", "",
        r"$$V(H)=m_Z^2H^\dagger H+\lambda(H^\dagger H)^2,\qquad q_H=m_Z^2/\omega^2,\qquad "
        r"M_\nu=\frac{v^2}{2}\left(C_5-Y_\nu M_R^{-1}Y_\nu^T\right).$$", "",
        r"The program stores $\widehat M_R=M_R/\omega$, $\widehat C_5=\omega C_5$, "
        r"$\widehat C_{HN}=\omega C_{HN}$, and $q_H$, so "
        r"$M_\nu=v^2[\widehat C_5-Y_\nu\widehat M_R^{-1}Y_\nu^T]/(2\omega)$. "
        r"The initial type-II coefficient is $\widehat C_5^{II}=2C_{II,\rm raw}f_{\rm raw}$. "
        "This factor of two is tested against the pre-existing local I+II mass matrix. "
        r"The no-half operator convention is $\mathcal L\supset-[N^TM_RN/2+\mathrm{h.c.}]"
        r"+[C_{HN}NN(H^\dagger H)+\mathrm{h.c.}]$, hence $M_R(\rho)=M_R-2C_{HN}\rho$. "
        r"Here $m_Z^2$ is Zhang's signed Lagrangian parameter, not the Z-boson mass. "
        r"If another card writes $V=-m^2H^\dagger H+\lambda(H^\dagger H)^2$, then $q_H=-m^2/\omega^2$.", "",
        "The illustrative boundary g=(0.53,0.55,0.60), lambda=0.12 and mu_start/omega=0.1 are declared diagnostic "
        "inputs, not solved P54 matching data. The additional qH=0 boundary is a diagnostic; it is not "
        "an electroweak-matched Higgs mass and the running generally regenerates nonzero qH. "
        "In particular, no scalar-card lambda_eff with another normalization "
        "is silently reused. Dimensionless thresholds are the primary results. The reference omega=10^14 GeV and "
        "v=246.22 GeV only convert them into explicitly synthetic terminal mass proxies.", "",
        "## 2. One-loop equations with every surviving sterile state", "",
        r"Write $H_x=Y_xY_x^\dagger$, $S_N=Y_\nu^\dagger Y_\nu$, "
        r"$T=\operatorname{tr}(3H_u+3H_d+H_e+H_\nu)$, and "
        r"$T_4=\operatorname{tr}(3H_u^2+3H_d^2+H_e^2+H_\nu^2)$. "
        r"Let $\mathcal B=16\pi^2d/d\log\mu$. In GUT hypercharge normalization:", "",
        r"$$\begin{aligned}"
        r"\mathcal B_{Y_u}&=\left[\frac32(H_u-H_d)+T-\frac{17}{20}g_1^2-\frac94g_2^2-8g_3^2\right]Y_u,\\"
        r"\mathcal B_{Y_d}&=\left[\frac32(H_d-H_u)+T-\frac14g_1^2-\frac94g_2^2-8g_3^2\right]Y_d,\\"
        r"\mathcal B_{Y_e}&=\left[\frac32(H_e-H_\nu)+T-\frac94g_1^2-\frac94g_2^2\right]Y_e,\\"
        r"\mathcal B_{Y_\nu}&=\left[\frac32(H_\nu-H_e)+T-\frac9{20}g_1^2-\frac94g_2^2\right]Y_\nu"
        r"+3C_5Y_\nu^*M_R-4Y_\nu M_R^\dagger C_{HN},\\"
        r"\mathcal B_{M_R}&=S_N^TM_R+M_RS_N-8m_Z^2C_{HN}.\end{aligned}$$", "",
        r"With $P_C=-3H_e/2+7H_\nu/2$ and $\alpha_C=4\lambda-3g_2^2+2T$:", "",
        r"$$\mathcal B_{C_5}=P_CC_5+C_5P_C^T+\alpha_CC_5,\qquad"
        r"\mathcal B_{C_5^{II}}=P_CC_5^{II}+C_5^{II}P_C^T+\alpha_CC_5^{II}.$$", "",
        r"$$\mathcal B_{g_i}=b_i g_i^3,\quad b=(41/10,-19/6,-7),$$", "",
        r"$$\mathcal B_\lambda=24\lambda^2-3\lambda(3g_2^2+\tfrac35g_1^2)"
        r"+\tfrac34g_2^4+\tfrac38(g_2^2+\tfrac35g_1^2)^2+4\lambda T-2T_4"
        r"-4\operatorname{Re}\operatorname{tr}(C_5Y_\nu^*M_RY_\nu^\dagger)"
        r"+16\operatorname{Re}\operatorname{tr}(C_{HN}^\dagger M_RS_N).$$", "",
        r"Define $A_H=12\lambda+2T-\frac9{10}g_1^2-\frac92g_2^2$. The remaining nonzero coefficients obey", "",
        r"$$\mathcal B_{C_{HN}}=A_HC_{HN}+4(C_{HN}S_N+S_N^TC_{HN}),$$", "",
        r"$$\mathcal B_{m_Z^2}=A_Hm_Z^2-4\operatorname{tr}(M_R^\dagger M_RS_N)"
        r"+8\operatorname{Re}\operatorname{tr}(M_R^\dagger M_RM_R^\dagger C_{HN}).$$", "",
        r"All these formulas are dimensionful; the code divides by powers of the fixed reference $\omega$ "
        r"to evolve the hatted variables. It does not identify $\omega$ with the running scale. "
        r"The beta function of $C_{BN}$ is homogeneous, so its zero boundary remains zero at this order.", "",
        "The CHN feedback signs have an independent CW-divergence test: for one real family take "
        r"$\mathscr M(h)=\left(\begin{smallmatrix}0&yh/\sqrt2\\yh/\sqrt2&M-Ch^2\end{smallmatrix}\right)$. "
        r"A polynomial regression of $\operatorname{tr}[(\mathscr M^\dagger\mathscr M)^2]$ yields "
        r"$[h^2]=2M^2y^2-4M^3C$ and $[h^4]=y^4/2-4My^2C+O(C^2)$. "
        "Calibrating against the renormalizable mass/quartic terms gives the displayed +8 M^3 C and "
        "+16 M y^2 C corrections. Taking the odd part in C removes the double-insertion contribution.", "",
        "The renormalizable baseline uses the transposed-row convention of "
        "[Antusch et al., Appendix D.1](https://arxiv.org/pdf/hep-ph/0501272). "
        "Their quartic is lambda_A/4, so lambda_A=4 lambda here. The +4 lambda coefficient is independently "
        "confirmed in the canonical convention by "
        "[Wang, Zhang and Zhou, Eq. (2.15)](https://arxiv.org/pdf/2302.08140). "
        "No sterile gauge contribution is added because the retained neutrinos are SM singlets.", "",
        "Crucially, the default equations include the three corrections in "
        "[Di Zhang (2024), Eqs. (5)–(6), (8)](https://arxiv.org/html/2405.18017v2): "
        "Weinberg feedback into Ynu and lambda and the replacement 1/2 by 7/2 in its own anomalous dimension. "
        "Our C5 is minus that paper's coefficient, explaining the feedback signs. "
        "An explicitly named negative control removes only these three Weinberg corrections while retaining "
        "the CHN and qH equations. The CHN signs do not flip with the Weinberg convention. "
        "The actual tree scalar exchange falsifies the old CHN=0 boundary assumption. CBN has no source "
        "from this nonderivative scalar-tree graph, but its finite loop matching remains uncomputed.", "",
        "The state space changes from 3 to 2 to 1 to 0 sterile columns. The full complex symmetric MR matrix "
        "runs within every interval. C5II runs continuously, while the total C5 receives type-I threshold increments. "
        "C5II is a tagged initial-source component along the common nonlinear trajectory, not a separately fitted "
        "or scheme-independent decomposition of loop effects.", "",
        "## 3. Why one frozen Weinberg matrix is insufficient", "",
        r"First restrict to $C_{HN}=0$ to isolate the Weinberg issue, and set $Q=Y_\nu M_R^{-1}Y_\nu^T$. "
        r"The product rule, including "
        r"$dM_R^{-1}=-M_R^{-1}(dM_R)M_R^{-1}$, gives", "",
        r"With $P=-3H_e/2+H_\nu/2$,", "",
        r"$$\mathcal B_Q=PQ+QP^T+\alpha_QQ+3(H_\nu C_5+C_5H_\nu^T),\qquad"
        r"\alpha_Q=2T-\frac9{10}g_1^2-\frac92g_2^2.$$", "",
        r"Thus for $C_{\rm eff}=C_5-Q$:", "",
        r"$$\mathcal B_{C_{\rm eff}}=PC_{\rm eff}+C_{\rm eff}P^T+\alpha_CC_{\rm eff}"
        r"+\underbrace{\left(4\lambda+\frac9{10}g_1^2+\frac32g_2^2\right)Q}_{\text{active-seesaw source}}.$$",
        "",
        "This is a structural, not numerical, distinction: above thresholds the type-I and type-II pieces have "
        "different flavor-universal running. Evolving only their initial sum with the below-threshold Weinberg "
        "equation loses the displayed source. The new 2024 contributions to beta C5 and beta Q cancel "
        "instantaneously in their difference, so this last equation retains its form. "
        "Nevertheless Ynu and lambda follow changed trajectories; the JSON quantifies the resulting difference "
        "from the system with only these corrections omitted. Both identities are checked directly in the "
        "canonical CHN=0 subsector; neither is claimed without its extra terms at nonzero CHN.", "",
        r"For actual nonzero $C=C_{HN}$, the additional contribution to $\mathcal B_Q$ is", "",
        r"$$\Delta_{HN}\mathcal B_Q=-4Y_\nu M_R^\dagger C M_R^{-1}Y_\nu^T"
        r"-4Y_\nu M_R^{-1}C M_R^*Y_\nu^T"
        r"+8m_Z^2Y_\nu M_R^{-1}C M_R^{-1}Y_\nu^T.$$",
        "Thus beta(C5-Q) acquires minus this expression. Keeping only a single total Weinberg coefficient "
        "would also lose this independent scalar-exchange effect. The implementation evolves every matrix "
        "and obtains the composite by matrix inversion, avoiding such a closure assumption.", "",
        "## 4. Derivation of the sign and block matching", "",
        r"At a chosen sterile split, rotate by $N_{\rm old}=U N_{\rm new}$, so "
        r"$Y'=YU$ and $M'=U^TMU$. Write", "",
        r"$$M'=\begin{pmatrix}A&B\\B^T&D\end{pmatrix},\qquad Y'=(Y_r,Y_h).$$", "",
        r"The algebraic heavy equation of motion is $N_h=-D^{-1}(B^TN_r+Y_h^T\ell H)$. "
        "Substitution, or the block-inverse identity, gives", "",
        r"$$\boxed{\ C_5^-=C_5^+-Y_hD^{-1}Y_h^T,\quad "
        r"Y_r^-=Y_r-Y_hD^{-1}B^T,\quad M_r^-=A-BD^{-1}B^T\ }.$$", "",
        r"Consequently $C_5^--Y_r^-(M_r^-)^{-1}(Y_r^-)^T=C_5^+-YM^{-1}Y^T$ "
        "exactly at tree level. The minus sign follows from the declared mass convention, not a convention-free "
        "choice of Weinberg operator sign. Iterated Schur complements equal the one-shot Schur complement; "
        "a non-diagonal numerical test also exercises the retained-Yukawa correction.", "",
        r"For $M(\rho)=M-2C_{HN}\rho$, differentiate the same Schur complement. Writing "
        r"$T_N=\binom{I}{-D^{-1}B^T}$ in retained/heavy order gives", "",
        r"$$C_{HN}^-=T_N^T C_{HN}'T_N"
        r"=C_{rr}-C_{rh}D^{-1}B^T-BD^{-1}C_{hr}"
        r"+BD^{-1}C_{hh}D^{-1}B^T.$$",
        "This formula is checked against an independent finite difference of the field-dependent Schur "
        "complement and against two successive eliminations. Merely truncating the retained corner is "
        "wrong for a general off-diagonal mass split. A CHN insertion in the induced Weinberg term instead "
        "belongs to dimension seven and is not included in this dimension-five truncation.", "",
        "At a Takagi threshold the cross block B vanishes up to roundoff. The general block formula is nevertheless "
        "used, so roundoff cannot silently spoil the tree identity. This is exact matrix algebra at leading seesaw "
        "order; it does not assert exact light eigenvalues of the full electroweak mass matrix or include dimension-six kinetic effects.", "",
        "## 5. Moving events and degeneracies", "",
        r"The solver detects $f(t)=t-\log\sigma_{\max}[\widehat M_R(t)]=0$, where $t=\log(\mu/\omega)$. "
        "The singular values are recomputed inside the event function. Neither the initial eigenvalues nor an "
        "initial ordering determine the threshold sequence.", "",
        r"For symmetric $M=X+iY$, the real symmetric matrix "
        r"$\begin{pmatrix}X&-Y\\-Y&-X\end{pmatrix}$ has positive eigenvectors $(a,b)$ yielding "
        r"$u=a+ib$ with $Mu=m u^*$. This gives $U^TMU=\operatorname{diag}(m)$ without square-root phase choices. "
        "A complete near-degenerate positive-mass cluster is removed together, using a declared relative log-mass "
        "tolerance of 10^-7. Exact degeneracy is insensitive to the remaining real-orthogonal Takagi freedom. "
        "A separate test uses arbitrary complex sterile basis transformations. Rank-deficient Takagi blocks "
        "are rejected for sterile threshold events rather than assigning a threshold to a genuinely massless state. "
        "The terminal light-neutrino Takagi interface separately allows a kernel. A minimal two-sterile test "
        "preserves its massless light mode through the whole corrected flow.", "",
        "The level-crossing negative control begins with M1=0.049 < M2=0.050 but a much larger second-column "
        "Yukawa. The second mass runs downward faster, and the state initially called N1 actually leaves first. "
        "The removed projector, rather than an ambiguous eigenvector phase, is saved and tested.", "",
        "## 6. Numerical checkpoints", "",
        "| Local-light case | Initial descending masses / omega | Moving event scales / omega |",
        "|---|---|---|"]
    for row in result["cases"]:
        initial = ", ".join(f"{x:.8f}" for x in row["history"]["initial_running_masses"])
        events = ", ".join(f"{x['mu_over_omega']:.8f}" for x in row["history"]["events"])
        lines.append(f"| {row['input_case']} | {initial} | {events} |")
    lines += ["", "Every event is followed by continued running, including a final zero-sterile interval. "
              "The JSON records all initial/final matrices, moving mass spectra, matching continuity, terminal "
              "PMNS/CKM moduli, Jarlskog invariants and fixed-v running mass proxies. They are synthetic "
              "diagnostics, not agreement with experiment.", "",
              "| Check | Residual | Pass |", "|---|---:|:---:|"]
    for row in result["checks"]:
        lines.append(f"| {row['name']} | {row['residual']:.3e} | {'yes' if row['passed'] else 'NO'} |")
    lines += ["", "## 7. Finite matching is still a separate gate", "",
              "[Zhang and Zhou, Eqs. (87)–(88)](https://arxiv.org/pdf/2107.12133) provide finite "
              "one-loop matching when the full type-I sector is integrated out. This includes hard vertex "
              "terms and field-normalization contributions; changing C5 alone is not the entire matching operation. "
              "That all-sterile formula is not silently applied to a partially active sterile EFT with a pre-existing "
              "type-II coefficient. The current solver deliberately uses tree matching and one-loop running. "
              "Finite logarithms/constant terms, all lower scalar thresholds, and the relevant higher-dimensional "
              "operator closure remain open.", "",
              "The running subsystem now includes all one-loop dimension-five feedback in the dipole-free "
              "nuSMEFT, with actual nonzero scalar-tree CHN. The retained PQ axion is a spectator here; its "
              "interactions and loop thresholds are not incorporated, so this is not the complete P54 EFT. "
              "Finite scalar/Yukawa/dipole matching and a matched Higgs mass boundary remain open. "
              "Dimension-six operators, double dimension-five insertions and exact finite-v pole matching are "
              "outside the implementation. The modest mass ratios in the benchmarks do not establish that "
              "omitted terms are negligible. A finite and operator-complete matching audit is required before "
              "promoting these diagnostic terminal quantities to physical predictions.", "",
              "The terminal applicability flags are explicitly false: the generated positive qH is much larger "
              "than the terminal scale squared. A mass-independent MS-bar ODE can be continued as a mathematical "
              "diagnostic, but retaining a light Higgs down to that endpoint requires quadratic finite matching "
              "and retuning. Fixed-v mass proxies are not a physically valid endpoint in these runs.", "",
              "Next: connect the properly matched lower scalar EFT and independently evolved six-invariant "
              "PS inputs, including the calculable conjugate-bidoublet matrices, "
              "then add finite sequential type-I/type-II matching with a consistent operator basis. A full flavor "
              "fit, Higgs pole mass, determinant or portal is not promoted by this checkpoint.", ""]
    OUTPUT.with_suffix(".md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
