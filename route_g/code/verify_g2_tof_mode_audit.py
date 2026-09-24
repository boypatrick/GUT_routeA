#!/usr/bin/env python3
"""G2-T: separate mass metrology, source-conditioned signs and event budgets.

Algebra/conditioning tests, not a detector simulation or a Gaussian scan.
The at-rest Higgs examples and error boxes are requirements, not calibration.
No action, source strength, or earlier output is changed by this verifier.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import verify_g2_local_conversion as g2
import verify_g2_xenon_feasibility as xe

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_tof_mode_audit"
C_M_NS = .299792458
CARD = dict(
    Lambda_GeV=5., baseline_m=10., energy_GeV=xe.CARD["mh_GeV"]/2,
    signed_examples=[0, 1, -1, 5], target_A=132, recoil_keV=10.,
    bounded_relative_p_error=.001, bounded_path_error_m=.001,
    bounded_time_error_ns=.020,
    example_scope="Higgs at rest for illustration only, not an LHC boost distribution, flux or calibrated instrument.",
    action_changed=False, Gaussian_scan=False, new_source_or_sensor_claimed=False,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mass_squared(p, time_ns, length_m):
    """GeV^2 for p in GeV; retain negative estimates instead of clipping."""
    return p*p*((C_M_NS*time_ns/length_m)**2-1)


def tof_ns(mass, p, length_m):
    return length_m*np.hypot(mass, p)/(C_M_NS*p)


def mass_box(p, time_ns, length_m, dp_fraction, dt_ns, dl_m):
    """Exact extrema in a rectangular error box entirely inside T>L/c."""
    p_lo, p_hi = p*(1-dp_fraction), p*(1+dp_fraction)
    t_lo, t_hi = time_ns-dt_ns, time_ns+dt_ns
    l_lo, l_hi = length_m-dl_m, length_m+dl_m
    if min(p_lo, t_lo, l_lo) <= 0 or C_M_NS*t_lo <= l_hi:
        raise ValueError("Box touches the unphysical/no-real-mass domain")
    return (np.sqrt(mass_squared(p_lo, t_lo, l_hi)),
            np.sqrt(mass_squared(p_hi, t_hi, l_lo)))


def recoil_inverse(beta, recoil_GeV, theta, target_mass):
    """Free elastic target at rest; not the S2 energy-only likelihood."""
    q = np.sqrt(recoil_GeV*(recoil_GeV+2*target_mass))
    denominator = beta*q*np.cos(theta)-recoil_GeV
    if not 0 < beta < 1 or denominator <= 0:
        raise ValueError("No positive-energy massive solution")
    energy = target_mass*recoil_GeV/denominator
    if energy <= recoil_GeV:
        raise ValueError("No positive outgoing energy")
    return energy, beta*energy, energy*np.sqrt(1-beta*beta)


def run():
    checks = []

    def check(name, ok, measured=None):
        row = dict(name=name, passed=bool(ok))
        if measured is not None:
            row["measured"] = float(measured)
        checks.append(row)

    paths = dict(
        local_action=Path(g2.__file__),
        xenon_kernel=Path(xe.__file__),
        higgs_source_code=ROOT/"code"/"verify_g2_higgs_source.py",
        higgs_source_result=ROOT/"output"/"g2_higgs_source.json",
    )
    hashes = {k: sha(p) for k, p in paths.items()}
    prior = json.loads(paths["higgs_source_result"].read_text())
    check("prior:no_failed_checks", not prior["summary"]["failed"])
    check("prior:code_hash", prior["verification"]["source_sha256"]
          == hashes["higgs_source_code"])
    check("prior:xenon_dependency_hash",
          prior["verification"]["dependencies"]["xenon_source"] == hashes["xenon_kernel"])
    bound = prior["global_bound"]["event_ceiling"]
    E, L = CARD["energy_GeV"], CARD["baseline_m"]
    examples = []
    for j in CARD["signed_examples"]:
        mass = CARD["Lambda_GeV"]*np.sqrt(25+j*j)
        p = np.sqrt((E-mass)*(E+mass))
        beta, gamma2 = p/E, E*E/(mass*mass)
        T = tof_ns(mass, p, L)
        m2 = mass_squared(p, T, L)
        check(f"j={j}:mass_shell_reconstruction", abs(m2/mass**2-1) < 2e-14)
        check(f"j={j}:proper_time", abs((T/np.sqrt(gamma2))/(L*mass/(C_M_NS*p))-1) < 2e-14)
        # Independent check of the logarithmic Jacobian, not a fit.
        x = np.array([p, T, L])
        jacobian = np.array([1., gamma2, -gamma2])
        step = 1e-6
        numerical = []
        for k in range(3):
            shift = np.zeros(3)
            shift[k] = step
            high = .5*np.log(mass_squared(*(x*np.exp(shift))))
            low = .5*np.log(mass_squared(*(x*np.exp(-shift))))
            numerical.append((high-low)/(2*step))
        check(f"j={j}:log_jacobian", np.max(abs(np.array(numerical)-jacobian)) < 2e-8)
        # Eight non-Gaussian, correlated small-error points realize S S^T.
        S = 1e-6*np.array([[2., 0., 0.], [.3, 1., 0.], [-.2, .4, .7]])
        errors = np.array([S@np.array(s) for s in itertools.product([-1., 1.], repeat=3)])
        log_m = np.array([.5*np.log(mass_squared(*(x*np.exp(e)))) for e in errors])
        propagated = jacobian@(S@S.T)@jacobian
        check(f"j={j}:correlated_nongaussian_variance",
              abs(np.var(log_m)/propagated-1) < 1e-7)
        box = mass_box(p, T, L, CARD["bounded_relative_p_error"],
                       CARD["bounded_time_error_ns"], CARD["bounded_path_error_m"])
        corners = []
        for sp, st, sl in itertools.product([-1, 1], repeat=3):
            corners.append(np.sqrt(mass_squared(
                p*(1+sp*CARD["bounded_relative_p_error"]),
                T+st*CARD["bounded_time_error_ns"],
                L+sl*CARD["bounded_path_error_m"])))
        check(f"j={j}:exact_box", np.allclose([min(corners), max(corners)], box, atol=1e-12))
        examples.append(dict(j=j, mass_GeV=float(mass), p_GeV=float(p),
                             beta=float(beta), gamma_squared=float(gamma2),
                             time_ns=float(T), proper_time_ns=float(T/np.sqrt(gamma2)),
                             reconstructed_mass_interval_GeV=[float(v) for v in box]))
    a, b, minus = examples[:3]
    check("signed_pair:identical_kinematics",
          all(b[k] == minus[k] for k in ["mass_GeV", "p_GeV", "beta", "time_ns"]))
    check("mass_box:adjacent_modes_disjoint",
          a["reconstructed_mass_interval_GeV"][1] < b["reconstructed_mass_interval_GeV"][0])
    check("mass_estimator:negative_not_clipped", mass_squared(10., .99*L/C_M_NS, L) < 0)
    for scale in [.25, 2., 7.]:
        check(f"tof_only:scale_degeneracy:{scale}",
              abs(tof_ns(scale*b["mass_GeV"], scale*b["p_GeV"], L)/b["time_ns"]-1) < 1e-14)

    def separation(dt):
        low_b = mass_box(b["p_GeV"], b["time_ns"], L, .001, dt, .001)[0]
        high_a = mass_box(a["p_GeV"], a["time_ns"], L, .001, dt, .001)[1]
        return low_b-high_a

    touching_dt = brentq(separation, 0, .1, xtol=1e-14)
    check("box:touching_threshold", abs(separation(touching_dt)) < 1e-11)
    check("box:20ps_is_below_touching_threshold", CARD["bounded_time_error_ns"] < touching_dt)

    # A discrete correlated source, NOT an actual collider distribution.
    source_p = np.array([20., 40., 60., 80.])
    source_t0 = np.array([0., 8., -4., 0.])
    source_weights = np.array([.1, .2, .3, .4])
    efficiency = np.array([.8, .6, .7, .5])
    arrival = source_t0+tof_ns(b["mass_GeV"], source_p, L)
    gate = (arrival < 40.).astype(float)
    before = source_weights*efficiency
    after = before*gate
    mean_before = float(before@source_p/before.sum())
    mean_after = float(after@source_p/after.sum())
    check("gate:source_normalized", abs(source_weights.sum()-1) < 1e-15)
    check("gate:acceptance_reduces_count", 0 < after.sum() < before.sum() <= 1)
    check("gate:observed_momentum_changes", abs(mean_after-mean_before) > 1)
    for weight in [0., .2, .75, 1.]:
        check(f"gate:old_event_ceiling:{weight}", weight*bound <= bound)
    time_gate = dict(scope="Synthetic four-atom joint source; no Gaussian or physical flux claim.",
                     p_GeV=source_p.tolist(), emission_ns=source_t0.tolist(),
                     production_weights=source_weights.tolist(), efficiency=efficiency.tolist(),
                     arrival_ns=arrival.tolist(), gate=gate.tolist(),
                     fraction_before=float(before.sum()), fraction_after=float(after.sum()),
                     conditional_mean_p_before_GeV=mean_before,
                     conditional_mean_p_after_GeV=mean_after)

    # Exact elastic inversion and its physical boundaries.
    target = CARD["target_A"]*xe.CARD["atomic_mass_GeV"]
    recoil_cases = []
    for mass, p in [(5., 2.), (b["mass_GeV"], b["p_GeV"]), (40., 100.)]:
        energy = np.hypot(mass, p)
        beta = p/energy
        endpoint = xe.endpoint(p, mass, target)
        for fraction in [1e-6, .1, .8]:
            er = endpoint*fraction
            q = np.sqrt(er*(er+2*target))
            costheta = (energy+target)*er/(p*q)
            theta = np.arccos(costheta)
            inferred = recoil_inverse(beta, er, theta, target)
            tag = f"recoil:m={mass:.4g}:fraction={fraction:g}"
            check(tag+":physical", 0 < costheta <= 1 and energy-er > 0)
            check(tag+":inverse", np.allclose(inferred, [energy, p, mass], rtol=2e-11))
            pout2 = p*p+q*q-2*p*q*costheta
            check(tag+":outgoing_shell", abs(((energy-er)**2-pout2)/mass**2-1) < 2e-13)
            recoil_cases.append(dict(mass_GeV=float(mass), p_GeV=float(p),
                                     recoil_GeV=float(er), theta_rad=float(theta),
                                     scope="Algebra validation; not a nuclear response extrapolation."))
    er = CARD["recoil_keV"]*1e-6
    q = np.sqrt(er*(er+2*target))
    theta = np.arccos((E+target)*er/(b["p_GeV"]*q))
    D = b["beta"]*q*np.cos(theta)-er
    angular_jacobian = b["beta"]*q*np.sin(theta)/D
    angular_budget = .01/angular_jacobian
    step = 1e-8
    numerical_angle_jac = (np.log(recoil_inverse(b["beta"], er, theta+step, target)[2])
                          - np.log(recoil_inverse(b["beta"], er, theta-step, target)[2]))/(2*step)
    check("directional:angular_jacobian", abs(numerical_angle_jac/angular_jacobian-1) < 1e-7)
    try:
        recoil_inverse(b["beta"], er, np.pi/2, target)
    except ValueError:
        check("directional:reject_negative_denominator", True)
    else:
        check("directional:reject_negative_denominator", False)
    directional = dict(target_mass_GeV=float(target), recoil_keV=CARD["recoil_keV"],
                       theta_degrees=float(np.degrees(theta)), denominator_GeV=float(D),
                       abs_dlogm_dtheta_per_radian=float(angular_jacobian),
                       angular_budget_for_one_percent_mass_rad=float(angular_budget),
                       angular_budget_for_one_percent_mass_arcsec=float(np.degrees(angular_budget)*3600),
                       scope="Local angular-error budget alone with beta and ER exact; not an attainable resolution or joint error forecast.",
                       other_missing_inputs=["known target isotope per event or mixture likelihood",
                                             "source timestamp", "flight length and incoming direction",
                                             "recoil direction", "joint energy/timing/direction response"])

    # Revalidate existing G2-R information, not a newly invented sign sensor.
    channels, recovered = [], []
    n, pi = g2.CARD["incoming_phi_n"], g2.CARD["incoming_com_momentum"]
    for l in range(g2.PACKET["min_l"], g2.PACKET["max_l"]+1):
        channels += g2.enumerate_channels(n, l, pi)["open_channels"]
    for row in channels:
        d, j = row["m"]-n, row["j"]
        if d == 0:
            continue
        ef = row["final_phi_energy"]+row["final_detector_energy"]
        en = np.hypot(g2.phi_mass(n), pi)
        l2 = g2.CARD["R"]**2*((ef-en)**2-g2.CARD["MD"]**2-pi*pi)
        inferred_j = (l2-d*d-abs(j)**2)/(2*d)
        check(f"source_sign:l={row['l']}:m={row['m']}", abs(inferred_j-j) < 2e-13)
        recovered.append(dict(l=row["l"], m=row["m"], j=j, reconstructed_j=float(inferred_j)))
    check("source_sign:old_channel_count", len(channels) == 25)
    # Non-degenerate Phi labels require the actual alpha, not an assumed sensor.
    for u, v in [(-2, -1), (-1, 0), (-1, 1), (0, 1), (1, 2)]:
        check(f"phi_labels:quarter_holonomy:{u}:{v}", g2.phi_mass(u) != g2.phi_mass(v))
    half_card = dict(g2.CARD, alpha=.5)
    check("phi_labels:half_holonomy_counterexample", g2.phi_mass(0, half_card) == g2.phi_mass(-1, half_card))
    interference = []
    for name, even, odd in [("pure_odd_is_blind", 0j, 1+0j),
                            ("quadrature_is_blind", 1+0j, 1j),
                            ("aligned_interferes", 1+0j, .2+0j)]:
        plus, minus_rate = abs(even+odd)**2, abs(even-odd)**2
        contrast = plus-minus_rate
        check("interference:"+name, abs(contrast-4*(even.conjugate()*odd).real) < 1e-14)
        if name != "aligned_interferes":
            check("interference:blind:"+name, plus == minus_rate)
        interference.append(dict(case=name, plus_rate=plus, minus_rate=minus_rate,
                                 difference=contrast, synthetic_amplitudes=True))
    check("sign:no_change_to_higgs_theorem", prior["sign_result"]["optimal_error"] == .5)
    check("dependencies:unchanged", all(sha(p) == hashes[k] for k, p in paths.items()))
    failed = [r["name"] for r in checks if not r["passed"]]
    result = dict(
        status="Bounded mathematical measurement audit complete; no feasible instrument, new interaction, or revived Higgs-xenon chain.",
        card=CARD, mass_examples=examples,
        adjacent_mass_tof_separation_ns=b["time_ns"]-a["time_ns"],
        exact_box_touching_time_halfwidth_ns=float(touching_dt),
        time_gate_example=time_gate, directional_reconstruction=directional,
        recoil_algebra_checks=recoil_cases, existing_source_sign_recovery=recovered,
        interference_examples=interference, retained_event_ceiling=bound,
        summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
        checks=checks, verification=dict(source_sha256=sha(__file__), dependencies=hashes),
    )
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# G2-T: mass measurement and signed-mode information are separate", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**.", "",
             "No action, coupling, physical source, response matrix or old result is changed. No Gaussian scan or actual hardware performance is claimed.", "",
             "## Different masses: TOF plus a genuinely independent constraint", "",
             "m²=p²[(cT/L)²-1]. TOF alone fixes m/p, not m. Neutral X has no ordinary magnetic-curvature momentum readout. Two time hits still only determine velocity.", "",
             "The logarithmic Jacobian for (p,T,L) is (1,gamma²,-gamma²); the full covariance must be used. Negative noisy m² estimates are not silently clipped.", "",
             "At-rest Higgs illustration only: Lambda=5 GeV, E=62.54 GeV, L=10 m. This is not a real LHC source distribution.", "",
             "| j | mass [GeV] | p [GeV] | TOF [ns] | Exact reconstructed-mass range [GeV] |",
             "|---:|---:|---:|---:|---:|"]
    for r in examples:
        lo, hi = r["reconstructed_mass_interval_GeV"]
        lines.append(f"| {r['j']} | {r['mass_GeV']:.6f} | {r['p_GeV']:.6f} | {r['time_ns']:.6f} | [{lo:.6f}, {hi:.6f}] |")
    lines += ["", "Intervals assume bounded |delta p|/p<=0.1%, |delta T|<=20 ps and |delta L|<=1 mm. These are uncalibrated requirements, not Gaussian sigmas, confidence intervals or event-by-event guarantees for the collider.",
              f"The j=0 and abs(j)=1 nominal TOFs differ by {result['adjacent_mass_tof_separation_ns']*1000:.6f} ps. Their exact error boxes first touch at a timing halfwidth of {touching_dt*1000:.6f} ps with the other two bounds fixed.", "",
              "## Same action, alternative constraint: directional elastic recoil", "",
              "For a known target mass M at rest, q=sqrt(ER(ER+2M)), D=beta q cos(theta)-ER, E=M ER/D, p=beta E and m=E sqrt(1-beta²). Require 0<beta<1, D>0 and positive outgoing energy. The reconstructed p is correlated with beta and ER, not an additional independent datum.",
              f"At A=132, ER=10 keV and the abs(j)=1 example, theta={directional['theta_degrees']:.8f} degrees. The local angular budget alone for 1% mass error is {directional['angular_budget_for_one_percent_mass_arcsec']:.6f} arcsec, assuming exact beta and ER. This severe conditioning is a diagnostic, not detector performance.",
              "Natural xenon also has an unobserved isotope mixture. The existing S2-only response has neither recoil-vector information nor a tagged source-to-hit X flight time. No directional or timing calibration is imported.", "",
              "## Same mass, opposite sign", "",
              "Existing undriven Phi-X source correlations can label j=n+l-m if the needed incoming/output labels are measured. With unobserved l, G2-R's full COM kinematic inversion works only for d=m-n nonzero and known incoming conditions. This is revalidated, not a constructed source or independent conservation-law test. The driven branch needs the work exchange too.",
              "A new sign-odd amplitude alone is insufficient: |Me+j Mo|²-|Me-j Mo|²=4j Re(Me* Mo). There must be coherent even/odd contributions in the same observed channel, a calibrated signed external response, or an informative preparation. No such new interaction is activated.",
              "Direct tree Higgs production plus the diagonal portal retains equal signed distributions and 50% optimal sign error.", "",
              "## Signal budget and next gate", "",
              f"The previous arbitrary-momentum, efficiency-one ceiling stays {bound:.8g}. Any added timing/momentum/direction gate <=1 can only lower it. No timing improvement revives the rejected on-shell-Higgs/single-pass-xenon chain.",
              "Next require an actual production tag and independent momentum/energy mechanism before instrument forecasts; for sign choose source correlation or specify a new sign-sensitive interaction for separate authorization/audit. Production observables avoiding rescattering remain preferable for the rejected source. Full response, background and coherent finite-pulse likelihoods are not supplied."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")
    print(json.dumps(result["summary"]))
    print(f"Adjacent TOF gap={result['adjacent_mass_tof_separation_ns']*1000:.6f} ps; box touching halfwidth={touching_dt*1000:.6f} ps")
    print(f"Directional angular-only 1% budget={directional['angular_budget_for_one_percent_mass_arcsec']:.6f} arcsec; retained count bound={bound:.8g}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
