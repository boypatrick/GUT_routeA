#!/usr/bin/env python3
"""Finite C5 matching with explicit applicability and live sterile thresholds.

Canonical type-I full C5 formula: all singlets removed, CHN=C5_UV=0,
leading light mass/M. Actual P54 partially retained, nonzero-CHN trajectories
receive ONLY the universal gauge/quartic subset, never a complete-fit flag.
"""
from __future__ import annotations
import json
import hashlib
import math
from pathlib import Path

import numpy as np
import verify_p54_sequential_seesaw as seq

RF=Path(__file__).resolve().parents[1]
OUT=RF/"output/p54_finite_seesaw_matching"
LOOP=16*np.pi**2


def finite_c5_parts(y,masses,mu,lam,g):
    """Renormalized d<=5 canonical type-I coefficients, row-Weyl convention.

    C5_here = -C5_ZhangZhou; g[0] is GUT normalized, gY^2=3*g1^2/5.
    Separate pieces so a subset cannot silently stand in for the full map.
    """
    masses=np.asarray(masses,float)
    if y.shape[1]!=len(masses) or np.any(masses<=0) or mu<=0:
        raise ValueError("Positive removed masses, common scale and matching columns required")
    logs=np.log(mu**2/masses**2)
    q=(y/masses)@y.T
    zh=float(np.dot(np.sum(abs(y)**2,axis=0),.5+logs))/LOOP
    # This row-family tensor is K_L^T in a literal two-left-Weyl convention.
    klrow=(y*((3+2*logs)/4))@y.conj().T/LOOP
    weights=2*lam*(1+logs)+(.6*g[0]**2+g[1]**2)*(1+3*logs)/4
    universal=-(y*(weights/masses))@y.T/LOOP
    tree=-q
    legs=-zh*tree-.5*(klrow@tree+tree@klrow.T)
    return {"tree":seq.sym(tree),"gauge_quartic":seq.sym(universal),
            "ZH":zh,"KL_row":klrow,"canonical_legs":seq.sym(legs),
            "canonical_typeI_C5":seq.sym(tree+universal+legs)}


def canonical_typeI_terminal(state,mu):
    """Restricted full C5 result; does not claim all SMEFT coefficients."""
    if np.linalg.norm(state["CHN"])>1e-14 or np.linalg.norm(state["C5"])>1e-14:
        raise ValueError("Canonical type-I terminal formula requires CHN=C5_UV=0; P54 is more general")
    if abs(state["qH"])>1e-14:
        raise ValueError("This implementation requires the declared leading qH/M^2=0 truncation")
    masses,u=seq.takagi(state["MR"])
    y=state["Ynu"]@u
    return finite_c5_parts(y,masses,mu,state["lambda"],state["g"])


def universal_threshold(state,heavy,u,mu):
    """One finite gauge/quartic C5 subset at an actual running Takagi block.

    Keeps the old nonzero CHN in running/tree pullback; does not set it to
    zero to borrow the canonical-model formula. CHN finite loops are missing.
    """
    out,record,retained=seq.decouple_block(state,heavy,u)
    mass=u.T@state["MR"]@u
    if seq.relative(mass,np.diag(np.diag(mass)))>1e-10:
        raise ValueError("Finite spectral kernel requires a Takagi mass basis")
    y=state["Ynu"]@u[:,heavy]
    masses=np.real(np.diag(mass)[heavy])
    part=finite_c5_parts(y,masses,mu,state["lambda"],state["g"])
    out["C5"]=seq.sym(out["C5"]+part["gauge_quartic"])
    record.update({"finite_C5_gauge_quartic_increment":seq.cjson(part["gauge_quartic"]),
        "finite_increment_norm":float(np.linalg.norm(part["gauge_quartic"])),
        "finite_matching_complete":False,"CHN_retained_not_zeroed":float(np.linalg.norm(state["CHN"])),
        "missing":["finite CHN/pre-existing-C5 operator insertions","mixed removed/retained-N graphs",
                   "correlated finite Yukawa/MR/Higgs matching and field normalization"]})
    return out,record,retained


def effective_beta(state):
    b=seq.beta(state)
    y,m=state["Ynu"],state["MR"]
    if not y.shape[1]: return b["C5"]
    inverse_y=np.linalg.solve(m,y.T)
    return seq.sym(b["C5"]-b["Ynu"]@inverse_y-y@np.linalg.solve(m,b["Ynu"].T)
        +y@np.linalg.solve(m,b["MR"]@inverse_y))


def run():
    local=json.loads(seq.INPUT.read_text())
    chn=json.loads(seq.CHN_INPUT.read_text())
    state=seq.load_local_light(local["cases"][0],seq.decode(chn["rays"][2]["CHN_matrix_times_omega"]))
    rng=np.random.default_rng(20260911)
    checks=[]
    def check(name,a,b=0.,tol=2e-9):
        a,b=np.asarray(a),np.asarray(b)
        e=float(np.linalg.norm(a-b)/max(1,np.linalg.norm(a),np.linalg.norm(b)))
        checks.append({"name":name,"residual":e,"tolerance":tol,"pass":e<tol})
    try:canonical_typeI_terminal(state,.05)
    except ValueError:check("actual_CHN_and_typeII_point_rejected_by_canonical_formula",0.)
    else:check("actual_CHN_and_typeII_point_rejected_by_canonical_formula",1.)

    # Explicit canonical model LIMIT for validation, not a replacement P54 point.
    canonical=seq.state_copy(state)
    canonical["CHN"]*=0;canonical["C5"]*=0;canonical["C5II"]*=0
    mu=.043
    parts=canonical_typeI_terminal(canonical,mu)
    tree_out,_,_=seq.decouple_block(canonical,list(range(3)))
    step=1e-4
    derivative=(canonical_typeI_terminal(canonical,mu*math.exp(step))["canonical_typeI_C5"]
               -canonical_typeI_terminal(canonical,mu*math.exp(-step))["canonical_typeI_C5"])/(2*step)
    check("full_canonical_C5_matching_scale_derivative",derivative,effective_beta(tree_out)-effective_beta(canonical))
    check("canonical_finite_C5_complex_symmetric",parts["canonical_typeI_C5"],parts["canonical_typeI_C5"].T)
    rotations=[seq.unitary(rng,3) for _ in range(6)]
    rotated=seq.family_transform(canonical,*rotations)
    rpart=canonical_typeI_terminal(rotated,mu)
    check("canonical_finite_C5_full_family_covariance",rpart["canonical_typeI_C5"],
          rotations[3].T@parts["canonical_typeI_C5"]@rotations[3])
    check("omitted_finite_legs_detected",float(np.linalg.norm(parts["canonical_legs"])<1e-7))
    masses,u=seq.takagi(canonical["MR"])
    y=canonical["Ynu"]@u
    # At mu=M, the independent wave-function calculation gives 1/2 and 3/4.
    equal=np.full(3,.035)
    eq=finite_c5_parts(y,equal,.035,canonical["lambda"],canonical["g"])
    check("degenerate_Higgs_wave_constant",eq["ZH"],np.linalg.norm(y)**2/(2*LOOP))
    check("degenerate_lepton_wave_constant",eq["KL_row"],3*y@y.conj().T/(4*LOOP))
    orth=np.linalg.qr(rng.normal(size=(3,3)))[0]
    eqrot=finite_c5_parts(y@orth,equal,.035,canonical["lambda"],canonical["g"])
    check("degenerate_Majorana_block_O3_invariance",eqrot["canonical_typeI_C5"],eq["canonical_typeI_C5"])
    # Universal gauge/quartic part is additive by removed spectral block.
    total=sum(finite_c5_parts(y[:,i:i+1],masses[i:i+1],mu,canonical["lambda"],canonical["g"])["gauge_quartic"] for i in range(3))
    check("nondegenerate_universal_block_additivity",total,parts["gauge_quartic"])
    for heavy in ([2],[1,2],[0,1,2]):
        before=seq.state_copy(state)
        after,_,_=seq.decouple_block(before,heavy,u)
        dplus=universal_threshold(before,heavy,u,mu*math.exp(step))[1]
        dminus=universal_threshold(before,heavy,u,mu*math.exp(-step))[1]
        d=(seq.decode(dplus["finite_C5_gauge_quartic_increment"])-seq.decode(dminus["finite_C5_gauge_quartic_increment"]))/ (2*step)
        def gauge_lambda_beta(s):
            off=seq.state_copy(s);off["g"]=np.zeros(3);off["lambda"]=0.
            return effective_beta(s)-effective_beta(off)
        check("partial_block_gauge_quartic_RG_identity_"+str(len(heavy)),d,
              gauge_lambda_beta(after)-gauge_lambda_beta(before))
    cases=[]
    for i,case in enumerate(local["cases"]):
        actual=seq.load_local_light(case,seq.decode(chn["rays"][2+i]["CHN_matrix_times_omega"]))
        # Stop below all sterile events but do not invent a physical EW endpoint.
        stop=1e-4
        base,base_history=seq.evolve(actual,.1,stop)
        finite,history=seq.evolve(actual,.1,stop,threshold_matcher=universal_threshold)
        check(f"case{i}_three_live_finite_events",len(history["events"]),3)
        check(f"case{i}_all_event_roots",max(abs(e["root_log_residual"]) for e in history["events"]))
        check(f"case{i}_CHN_kept_nonzero_at_all_thresholds",float(any(e["CHN_retained_not_zeroed"]<=0 for e in history["events"])))
        check(f"case{i}_no_event_claims_complete_matching",float(any(e["finite_matching_complete"] for e in history["events"])))
        check(f"case{i}_finite_result_differs_from_tree",float(np.linalg.norm(finite["C5"]-base["C5"])<1e-8))
        transformed=seq.family_transform(actual,*rotations)
        rf,rh=seq.evolve(transformed,.1,stop,threshold_matcher=universal_threshold)
        check(f"case{i}_finite_event_trajectory_covariance",rf["C5"],rotations[3].T@finite["C5"]@rotations[3],tol=1e-7)
        scale_rows=[]
        for factor in (.8,1.25):
            varied,vh=seq.evolve(actual,.1,stop,threshold_matcher=universal_threshold,threshold_scale_factor=factor)
            scale_rows.append({"factor":factor,"C5_difference_norm":float(np.linalg.norm(varied["C5"]-finite["C5"])),
                "events":[e["mu_over_omega"] for e in vh["events"]]})
        cases.append({"input_case":i,"input_scope":"actual P54 group/light dictionary on pre-existing synthetic families; no physical fit",
            "history":history,"terminal_C5":seq.cjson(finite["C5"]),
            "tree_vs_finite_relative_C5":seq.relative(finite["C5"],base["C5"]),
            "matching_scale_variation":scale_rows,"physical_endpoint":False})
    paths=[Path(__file__),Path(seq.__file__),seq.INPUT,seq.CHN_INPUT]
    return {"schema":"p54-finite-seesaw-subsectors-v1","date":"2026-09-11",
        "canonical_terminal_limit":{"C5_finite_complete_at_declared_truncation":True,
            "assumptions":["all N removed","CHN=0","C5_UV=0","leading qH/M^2=0"],
            "result":seq.cjson(parts["canonical_typeI_C5"]),"not_actual_P54_boundary":True},
        "P54_sequential_scope":"universal gauge/quartic C5 finite subset only; CHN/tree-II and retained N never silently discarded",
        "cases":cases,"full_P54_finite_seesaw_matching":False,"physical_fit_enabled":False,
        "references":[{"url":"https://arxiv.org/html/2107.12133v3","equations":"49, 51, 89, 91; C5_here=-C5_paper; gY2=3/5 g1_GUT2"},
                      {"url":"https://arxiv.org/html/2201.00840v2","equations":"42, 45; independent degenerate wave-function constants"}],
        "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "checks":checks,"summary":{"passed":sum(c["pass"] for c in checks),"total":len(checks),"all_pass":all(c["pass"] for c in checks)}}


def markdown(r):
    lines=["# Finite seesaw matching: canonical limit and actual sequential subset", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}.", "",
        "**Full P54 finite seesaw matching remains incomplete.** The canonical all-N-decoupled C5 formula is implemented only for CHN=C5_UV=0 and leading qH/M^2=0. Actual P54 input has nonzero CHN and type-II C5, so that restricted formula explicitly rejects it.", "",
        "The universal gauge/quartic finite C5 term is now applied at each moving Takagi mass block while retaining the existing nonzero CHN in the trajectory. Finite CHN/pre-existing-C5 insertions, mixed removed/retained-N graphs and correlated finite parameter matching are still absent.", "",
        "| Diagnostic case | Live finite thresholds | Relative C5 shift from tree thresholds |", "|---|---:|---:|"]
    lines += [f"| {c['input_case']} | {len(c['history']['events'])} | {c['tree_vs_finite_relative_C5']:.8g} |" for c in r["cases"]]
    lines += ["", "These shifts use the old synthetic families in the actual P54 group/light dictionary. They are not fitted neutrino predictions. The matching-scale variations are diagnostics of an incomplete map, not calibrated theory errors.", "",
        "The canonical full-C5 scale derivative and the partial-block universal scale derivatives reproduce the corresponding beta differences. Degenerate blocks and generic family transformations are checked.", "",
        "Sources: [Zhang-Zhou, Eqs. 49/51/91](https://arxiv.org/html/2107.12133v3); [Ohlsson-Pernow, Eqs. 42/45](https://arxiv.org/html/2201.00840v2).", "",
        "| Check | Residual | Pass |", "|---|---:|:---:|"]
    lines += [f"| {c['name']} | {c['residual']:.3e} | {c['pass']} |" for c in r["checks"]]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    r=run();OUT.with_suffix(".json").write_text(json.dumps(r,indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(r))
    print(json.dumps({"summary":r["summary"],"finite_C5_shifts":[c["tree_vs_finite_relative_C5"] for c in r["cases"]],
                      "failed":[c for c in r["checks"] if not c["pass"]]},indent=2))
    if not r["summary"]["all_pass"]: raise SystemExit(1)
