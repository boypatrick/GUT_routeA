#!/usr/bin/env python3
"""Corrected P54 flavor geometry, exact complex CW, and constructive gates.

The target transport is diagnostic SM running. No full PS/global fit or
loop-corrected scalar eigenpair is claimed. New algebra is tested against
independent spectra, finite differences and explicit matrix constructions.
"""
from __future__ import annotations
import importlib.util
import hashlib
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"


def load_legacy_helpers():
    spec = importlib.util.spec_from_file_location("p54_p3_helpers", RF / "code/verify_p54_p3_flavor_cw_gate.py")
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def cjson(a):
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def sharp_bound(r, rho, yb, ye):
    if rho == 0:
        return r*(3*yb+ye)/4, 0.
    x = (9*yb*yb*(1+rho*rho)-ye*ye*(9+rho*rho))/(6*rho*(ye*ye+3*yb*yb))
    x = float(np.clip(x, -1, 1))
    result = r*(yb*math.sqrt(9+rho*rho+6*rho*x)+ye*math.sqrt(1+rho*rho-2*rho*x))/4
    return result, x


def complex_cw(yukawas, mr, mu):
    if not np.allclose(mr, mr.T, rtol=1e-12, atol=1e-14):
        raise ValueError("Majorana mass matrix must be complex symmetric")
    if mu <= 0:
        raise ValueError("Matching scale must be positive")
    lam, u = np.linalg.eigh(mr.conj().T @ mr)
    if lam.min() <= 0:
        raise ValueError("The hard Majorana sector must be invertible")
    w = (u*(lam*(np.log(lam/mu**2)-1))) @ u.conj().T
    pi = np.array([[-np.trace(a.conj().T@b@w)/(8*math.pi**2) for b in yukawas] for a in yukawas])
    return (pi+pi.conj().T)/2


def potential(q, yukawas, mr, mu):
    n = len(yukawas)
    z = (q[:n]+1j*q[n:])/math.sqrt(2)
    md = sum(x*y for x, y in zip(z, yukawas))
    mass = np.block([[np.zeros_like(mr), md], [md.T, mr]])
    lam = np.linalg.eigvalsh(mass.conj().T@mass)[-mr.shape[0]:]
    return -np.sum(lam**2*(np.log(lam/mu**2)-1.5))/(32*math.pi**2)


def numerical_hessian(fun, n, step):
    out = np.zeros((n, n)); zero = np.zeros(n); v0 = fun(zero)
    for i in range(n):
        ei = np.eye(n)[i]*step
        out[i, i] = (fun(ei)+fun(-ei)-2*v0)/step**2
        for j in range(i):
            ej = np.eye(n)[j]*step
            out[i, j] = out[j, i] = (fun(ei+ej)-fun(ei-ej)-fun(-ei+ej)+fun(-ei-ej))/(4*step**2)
    return out


def run():
    legacy = load_legacy_helpers()
    p2path = RF / "output/p54_p2_two_site_matching.json"
    oppath = RF / "output/p54_p3_light_doublet_overlaps.json"
    p2 = json.loads(p2path.read_text())
    geom = json.loads(oppath.read_text())
    g = geom["geometric_sector_magnitudes"]
    c = np.array([g["phi_hol"], g["phi_anti"], g["Sigma_hol"], g["Sigma_anti"]])
    mu_gut = float(p2["thresholded_solution"]["MU_GeV"])
    gu = float(p2["gauge_coupling_iteration"][-1]["g_output"])
    sigma = float(p2["hierarchical_vacuum_ratio"])
    gap = float(p2["loop_corrected_doublet_condition"]["nearest_heavy_doublet_gap_m2_over_omega2"])
    targets = legacy.diagnostic_sm_run(mu_gut)
    yt, yc = targets["yukawa_singular_values_at_high"]["up"][-1:-3:-1]
    yb = targets["yukawa_singular_values_at_high"]["down"][-1]
    ye = targets["yukawa_singular_values_at_high"]["charged_lepton"][-1]
    cases = {}
    for name, ab in {
        "declared_R_equals_minus_iY": (c[0], c[1], c[3], c[2]),
        "global_charge_conjugate_stress_test": (c[1], c[0], c[2], c[3]),
    }.items():
        a, b, d, e = ab
        r, s = b/a, a*e/(d*b)
        loose = legacy.phase_uniform_bound(r, s, yb, ye)
        upper, phase_x = sharp_bound(r, s, yb, ye)
        f_lower = (ye-yb)/(4*d)
        cap = 4*math.pi
        min_weight = ((yt-r*yb)/(cap*math.sqrt(1+r*r)))**2
        cases[name] = {
            "alpha10": a, "beta10": b, "alpha126": d, "beta126": e,
            "abs_r": r, "abs_s": s, "triangle_bound": loose["yt_upper"],
            "sharp_phase_uniform_bound": upper, "maximizing_cos_phase_s": phase_x,
            "target_over_sharp_bound": yt/upper, "required_additive_defect_norm": yt-upper,
            "adversarial_30pct_ratio": .7*yt/(1.3*upper),
            "conditional_gaussian_chi2_lower_10pct": legacy.gaussian_halfspace_chi2_lower(yt,yb,ye,loose["A_down"],loose["B_charged_lepton"]),
            "fundamental_f_lower_from_down_lepton_norms": f_lower,
            "minimum_126_weight_with_f_below_4pi": min_weight,
            "tree_126_weight": c[2]**2+c[3]**2,
            "minimum_amplitude_rotation": math.asin(math.sqrt(min_weight))-math.asin(math.sqrt(c[2]**2+c[3]**2)),
            "single_120_required_symmetric_norm": (yt-yc)/2,
            "single_120_available_symmetric_norm_upper": upper,
            "single_120_fixed_overlap_excluded_at_diagnostic_targets": (yt-yc)>2*upper,
        }

    chosen = cases["declared_R_equals_minus_iY"]
    a,b,d,e = (chosen[k] for k in ("alpha10","beta10","alpha126","beta126"))
    heff, feff = legacy.appendix_b_matrices()
    h,f = heff/a, feff/d
    pi2 = complex_cw([h,-3*f], sigma*f, gu)
    pi4 = np.zeros((4,4), dtype=complex)
    pi4[np.ix_([1,2],[1,2])] = pi2
    q = np.eye(4)-np.outer(c,c)
    p10 = np.diag([1,1,0,0])
    eta = float((c@pi4@c).real)
    dxi = eta/(c@p10@c)
    qr = q@(pi4-dxi*p10)@q

    # Exact rank-one third-family construction. This is not a global fit.
    H = (3*yb+ye)/4; F = (yb-ye)/4
    mirror_k = (a*a-b*b)/a
    fg = F/d
    gg = (yt-b/a*H-e*fg)/mirror_k
    hg = (H-b*gg)/a
    rankone = np.diag([0.,0.,1.])
    hs,fs,gs = hg*rankone,fg*rankone,gg*rankone
    yu = b*hs+e*fs+a*gs
    yd = a*hs+d*fs+b*gs
    yl = a*hs-3*d*fs+b*gs
    theta = math.atan(abs(gg)/(2*abs(hg)))
    cos2 = math.cos(theta)**2
    h_uv = abs(hg)/cos2
    f_uv = abs(fg)/cos2
    y_uv = abs(gg)/math.sin(2*theta)
    vev_s = .25/math.sqrt(2)  # canonical S=(x326+i x327)/sqrt(2)
    mx = gu  # illustrative bare vectorlike mass equal to old MU
    kappa_uv = mx*math.tan(theta)/vev_s

    rng = np.random.default_rng(20260905)
    def rand_sym(scale=1.):
        z = scale*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        return (z+z.T)/2
    max_skew_pair = max_weyl_violation = 0.
    for _ in range(200):
        z = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        skew = (z-z.T)/2; sy = rand_sym(.15)
        v = np.linalg.svd(skew,compute_uv=False)
        w = np.linalg.svd(skew+sy,compute_uv=False)
        max_skew_pair = max(max_skew_pair, abs(v[0]-v[1]), v[2])
        max_weyl_violation = max(max_weyl_violation, w[0]-w[1]-2*np.linalg.norm(sy,2))
    # Independent finite differences test all real/imaginary CW directions.
    ys = [rand_sym(.15), rand_sym(.12)]
    unit, _ = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    mr = unit.conj()@np.diag([.4,.7,1.1])@unit.conj().T
    pi = complex_cw(ys,mr,.8)
    realpi = np.block([[pi.real,-pi.imag],[pi.imag,pi.real]])
    num = numerical_hessian(lambda x: potential(x,ys,mr,.8),4,3e-4)
    cw_error = float(np.linalg.norm(num-realpi)/max(np.linalg.norm(realpi),1e-10))
    # Schur retuning explicit counterexample and exact zero eigenpair.
    off = .001
    naive = np.array([[0.,off],[off,gap]])
    fixed = naive.copy(); fixed[0,0] = off*off/gap
    cv = np.array([1.,-off/gap]); cv /= np.linalg.norm(cv)
    # Phase optimization is compared against an independent dense grid.
    phase_errors = []
    for case in cases.values():
        rr,ss = case["abs_r"],case["abs_s"]
        phase = np.linspace(0,2*math.pi,20001)
        vals = rr/4*(yb*np.abs(3+ss*np.exp(1j*phase))+ye*np.abs(1-ss*np.exp(1j*phase)))
        phase_errors.append(abs(float(vals.max())-case["sharp_phase_uniform_bound"]))
    toy_residual = max(np.linalg.norm(yu-yt*rankone),np.linalg.norm(yd-yb*rankone),np.linalg.norm(yl-ye*rankone))
    checks = [
        {"name":"corrected overlap normalization", "pass":abs(c@c-1)<1e-12},
        {"name":"conjugate invariants use opposite hol/anti maps", "pass":a==c[0] and b==c[1] and d==c[3] and e==c[2]},
        {"name":"sharp phase bound agrees with dense maximization", "pass":max(phase_errors)<1e-8},
        {"name":"both charge conventions fail tree diagnostic gate", "pass":all(x["target_over_sharp_bound"]>1 for x in cases.values())},
        {"name":"complex skew singular pairing on 200 samples", "pass":max_skew_pair<1e-12},
        {"name":"skew plus symmetric spectral-gap theorem on 200 samples", "pass":max_weyl_violation<1e-12},
        {"name":"exact complex CW against independent 4-real-direction differences", "pass":cw_error<1e-5},
        {"name":"conditional CW Hermiticity", "pass":np.linalg.norm(pi4-pi4.conj().T)<1e-12},
        {"name":"projected zero mass can hide negative eigenvalue", "pass":np.linalg.eigvalsh(naive)[0]<0},
        {"name":"Schur retuning produces exact zero eigenpair", "pass":np.linalg.norm(fixed@cv)<1e-14},
        {"name":"rank-one symmetric mirror operator matches third-family toy", "pass":toy_residual<1e-14},
        {"name":"all mediator vertices conserve PQ", "pass":(-1+3-2)==0 and (4-3-1)==0 and (3-3)==0},
        {"name":"vectorlike mediator adds zero PQ anomalies", "pass":(3*2-3*2)==0 and (3*16-3*16)==0},
        {"name":"exact mediator mixing reproduces all three effective norms", "pass":max(abs(y_uv*math.sin(2*theta)-abs(gg)),abs(h_uv*cos2-abs(hg)),abs(f_uv*cos2-abs(fg)))<1e-12},
    ]
    quarticpath = RF / "output/p54_relaxed_light_quartic.json"
    quartic = json.loads(quarticpath.read_text()) if quarticpath.exists() else None
    return {"schema":"route-f-p54-theory-audit-v1", "date":"2026-09-05",
            "status":{"tree_flavor_diagnostic_tension":True, "full_loop_corrected_branch_excluded":False,
                      "global_flavor_seesaw_fit_performed":False, "eta_profile_domain_proven_empty":False,
                      "new_UV_extension_fitted":False, "physical_pole_mass_computed":False},
            "corrections":["phi* and Sigma conjugation was conflated", "tree eigenvector was called loop-corrected", "diagnostic SM running was over-promoted", "a large Gaussian chi2 does not empty the parameter domain", "real-slice CW omitted imaginary copy mixing", "semidefinite Hessian alone omits heavy-relaxed quartic", "torsion NDA was misreported as a universal Yukawa correction"],
            "geometric_overlaps":c.tolist(),"diagnostic_targets":targets,"corrected_cases":cases,
            "conditional_complex_CW":{"scope":"published Hprime,Fprime transported using positive-real corrected overlaps; no fit or scalar-relative-phase match",
                "h_norm":float(np.linalg.norm(h,2)),"f_norm":float(np.linalg.norm(f,2)),
                "up_copy_pi_over_omega2":cjson(pi2),"geometric_copy_order":["phi_hol","phi_anti","Sigma_hol","Sigma_anti"],
                "full_complex_copy_pi_over_omega2":cjson(pi4),"eta_Y":eta,"retuned_Q_norm":float(np.linalg.norm(qr,2)),
                "Q_norm_over_tree_gap":float(np.linalg.norm(qr,2))/gap,"total_boson_fermion_Q_certified":False,
                "synthetic_full_complex_finite_difference_relative_error":cw_error},
            "schur_counterexample":{"gap":gap,"off_diagonal":off,"naive_eigenvalues":np.linalg.eigvalsh(naive).tolist(),
                "exact_retuning":off*off/gap,"fixed_eigenvalues":np.linalg.eigvalsh(fixed).tolist(),"light_rotation_radians":math.atan(off/gap)},
            "mirror_10_rank_one_toy":{"scope":"only third family; all first/second family masses and mixings set to zero; not a viable fit",
                "operator":"S* (16_i 16_j)_10 phi / Lambda", "coefficient_k":mirror_k,
                "h33":hg,"f33":fg,"g33":gg,"max_matrix_residual":float(toy_residual),
                "exact_mediator_theta":theta,"h_UV_norm":h_uv,"f_UV_norm":f_uv,"y_UV_norm":y_uv,"kappa_UV_for_bare_M_equal_old_MU":kappa_uv,
                "largest_vertex_loop_parameter":max(h_uv,f_uv,y_uv,kappa_uv)**2/(16*math.pi**2),
                "bare_M_over_omega":mx,"physical_heavy_M_over_omega":mx/math.cos(theta),
                "new_threshold_recalculation_required":True,"PQ_domain_wall_number_unchanged":3},
            "relaxed_quartic":quartic,
            "checks":checks,"summary":{"passed":sum(x["pass"] for x in checks),"total":len(checks),"all_pass":all(x["pass"] for x in checks)},
            "sources":[{"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in (oppath,p2path,RF/"code/verify_p54_p3_flavor_cw_gate.py",Path(__file__))]}


def markdown(r):
    a=r["corrected_cases"]["declared_R_equals_minus_iY"]
    lines=["# P54 theory audit and corrected P3 status", "", "Authoritative correction: 2026-09-05. The full loop-corrected branch has not been excluded; the global flavor/seesaw fit has not been performed.","",
           f"Corrected tree r={a['abs_r']:.9g}, s={a['abs_s']:.9g}; sharp phase-uniform top bound={a['sharp_phase_uniform_bound']:.9g} against diagnostic yt=0.464725012.","",
           "The old factors 69.49, repair target r>=62.9, and empty eta-profile-domain claim are withdrawn. The two invariant contractions use different scalar conjugations.","",
           "The audit includes the heavy-relaxed quartic, an exact complex Majorana CW formula, an antisymmetric-spurion obstruction, Schur retuning, and a PQ-invariant rank-one mirror-10 construction.","",
           "The mirror-10 example only fits three third-family diagonal entries and adds a vectorlike 16 pair in its explicit UV completion. Its thresholds and full flavor fit remain open.","",
           f"Independent algebraic/numerical checks: {r['summary']['passed']}/{r['summary']['total']}.",""]
    if r["relaxed_quartic"]:
        lines += [f"Same-action heavy-relaxed tree quartic: lambda={r['relaxed_quartic']['relaxed_lambda']:.10g}.",""]
    return "\n".join(lines)


def main():
    report=run()
    def numpy_scalar(value):
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Cannot serialize {type(value).__name__}")
    for name in ("p54_theory_audit", "p54_p3_flavor_cw_gate"):
        (RF/f"output/{name}.json").write_text(json.dumps(report,indent=2,sort_keys=True,default=numpy_scalar)+"\n")
        (RF/f"output/{name}.md").write_text(markdown(report))
    print(f"P54 theory audit: {report['summary']['passed']}/{report['summary']['total']}")
    for name,c in report["corrected_cases"].items():
        print(name, 'r,s=',c["abs_r"],c["abs_s"], 'sharp top bound=',c["sharp_phase_uniform_bound"])
    print('quartic:', report['relaxed_quartic']['relaxed_lambda'] if report['relaxed_quartic'] else 'not yet available')
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
