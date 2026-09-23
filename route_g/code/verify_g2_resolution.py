#!/usr/bin/env python3
"""G2-S: rate-weighted recoil-sign inference with explicit resolution cards.

These are conditional engineering forecasts, not detector calibration or an
absolute collision probability. The incoming distribution is an incoherent
ensemble; the drive is the long-time weak-scattering Floquet rate limit.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import ndtr

import verify_g2_local_conversion as g2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_resolution"
SIGMA_I = [0.0, 0.02, 0.1, 0.2, 0.4]
SIGMA_D = [0.02, 0.05, 0.1, 0.2, 0.4]
DRIVE = dict(epsilon=0.5, Omega=0.2)
Z_CUT = 10.0
RESOLVED_ERROR = 0.05
SQRT2PI = np.sqrt(2 * np.pi)


def incoming_quadrature(sigma, order=48):
    """Truncated-normal quadrature, without renormalizing the discarded tail."""
    p0 = g2.CARD["incoming_com_momentum"]
    if sigma == 0:
        return np.array([p0]), np.array([1.0]), 0.0
    lower = -p0 / sigma
    norm = ndtr(p0 / sigma)
    # Piecewise GL integrates the Gaussian accurately even for a narrow beam.
    edges = np.linspace(lower, Z_CUT, int(np.ceil(Z_CUT - lower)) + 1)
    nodes, weights = leggauss(order)
    zs, ws = [], []
    for left, right in zip(edges[:-1], edges[1:]):
        z = (left + right) / 2 + (right - left) / 2 * nodes
        zs.extend(z)
        ws.extend(weights * (right - left) / 2 * np.exp(-z*z/2) / (SQRT2PI * norm))
    return p0 + sigma * np.array(zs), np.array(ws), float(ndtr(-Z_CUT) / norm)


def kernel(p, l, q=0, driven=False):
    n, m = g2.CARD["incoming_phi_n"], 0
    j = n + l - m
    en, el = np.hypot(g2.phi_mass(n), p), np.hypot(g2.detector_mass(l), p)
    energy = en + el + (q * DRIVE["Omega"] if driven else 0)
    ma, mb = g2.phi_mass(m), g2.detector_mass(j)
    rad = (energy-ma-mb)*(energy+ma+mb)*(energy-ma+mb)*(energy+ma-mb)
    pf = np.sqrt(np.maximum(rad, 0)) / (2*energy)
    pf = np.where(energy > ma+mb, pf, 0)
    coupling = g2.CARD["g"] * (DRIVE["epsilon"]/2 if q else 1)
    rate = coupling**2 * pf / (16*np.pi*en*el*energy)
    return pf, rate


def components(sigma_i, driven, order=48):
    p, iw, tail = incoming_quadrature(sigma_i, order)
    labels = np.arange(g2.PACKET["min_l"], g2.PACKET["max_l"]+1)
    w = np.exp(-labels**2/(2*g2.PACKET["width"]**2))
    w /= w.sum()
    rows = []
    for l, sign in [(0, 1), (-2, -1)]:
        for q in ([-1, 0, 1] if driven else [0]):
            pf, rate = kernel(p, l, q, driven)
            rows.append(dict(l=l, sign=sign, q=q, means=pf,
                             rates=iw*w[labels.tolist().index(l)]*rate))
    return rows, p, iw, tail


def density(x, means, weights, sigma):
    values = np.asarray(x)
    return np.exp(-0.5*((values[..., None]-means)/sigma)**2) @ weights / (SQRT2PI*sigma)


def interval_mass(left, right, means, weights, sigma):
    a, b = (left-means)/sigma, (right-means)/sigma
    # Use survival probabilities in the positive tail: no 1-Phi cancellation.
    probability = np.where(a >= 0, ndtr(-a)-ndtr(-b), ndtr(b)-ndtr(a))
    return float(np.dot(weights, probability))


def classify(rows, sigma_d, grid_divisor=4):
    means = [np.concatenate([r["means"] for r in rows if r["sign"] == s]) for s in [1, -1]]
    weights = [np.concatenate([r["rates"] for r in rows if r["sign"] == s]) for s in [1, -1]]
    total = sum(float(w.sum()) for w in weights)
    weights = [w/total for w in weights]
    priors = [float(w.sum()) for w in weights]
    low = min(float(m.min()) for m in means)-Z_CUT*sigma_d
    high = max(float(m.max()) for m in means)+Z_CUT*sigma_d
    grid = np.linspace(low, high, int(np.ceil((high-low)*grid_divisor/sigma_d))+1)
    # Chunking limits temporary matrices in the broadest incoming card.
    delta = np.concatenate([density(xs, means[0], weights[0], sigma_d)
                            - density(xs, means[1], weights[1], sigma_d)
                            for xs in np.array_split(grid, max(1, len(grid)//128))])

    def difference(x):
        return float(density(x, means[0], weights[0], sigma_d)
                     - density(x, means[1], weights[1], sigma_d))

    roots = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], delta[:-1], delta[1:]):
        if fa*fb < 0:
            roots.append(float(brentq(difference, a, b, xtol=1e-13)))
    # Exact component CDF integration between numerically resolved boundaries.
    # Outside [low,high] the omitted response probability is <=2 Phi(-10).
    wrong = [0.0, 0.0]
    for left, right in zip([low]+roots, roots+[high]):
        choose = 0 if difference((left+right)/2) >= 0 else 1
        other = 1-choose
        wrong[other] += interval_mass(left, right, means[other], weights[other], sigma_d)
    error = sum(wrong)
    moments = []
    for ms, ws, prior in zip(means, weights, priors):
        mean = float(np.dot(ws, ms)/prior)
        var = float(np.dot(ws, (ms-mean)**2)/prior)
        moments.append(dict(mean_momentum=mean, intrinsic_variance=var,
                            measured_standard_deviation=float(np.sqrt(var+sigma_d**2))))
    return dict(bayes_error=error, majority_only_error=min(priors), priors=priors,
                false_minus_given_plus=wrong[0]/priors[0],
                false_plus_given_minus=wrong[1]/priors[1],
                bayes_boundaries=roots, total_pair_rate=total, sign_moments=moments,
                resolved_at_declared_5pct=error <= RESOLVED_ERROR,
                response_tail_error_bound=float(2*ndtr(-Z_CUT)))


def evaluate(sigma_i, sigma_d, driven, order=48, grid_divisor=4):
    rows, p, iw, tail = components(sigma_i, driven, order)
    result = classify(rows, sigma_d, grid_divisor)
    if driven:
        tagged = []
        for q in [-1, 0, 1]:
            qr = classify([r for r in rows if r["q"] == q], sigma_d, grid_divisor)
            tagged.append(dict(q=q, rate=qr["total_pair_rate"], error=qr["bayes_error"]))
        result["ideal_q_tagged_error"] = sum(r["rate"]*r["error"] for r in tagged)/result["total_pair_rate"]
        result["q_tagged_components"] = tagged
    # pf <= Eout/2 implies Kq <= |gq|²/(32*pi*mn*mu_l).
    labels = np.arange(g2.PACKET["min_l"], g2.PACKET["max_l"]+1)
    w = np.exp(-labels**2/(2*g2.PACKET["width"]**2)); w /= w.sum()
    bound = 0.0
    for l in [0, -2]:
        bound += w[labels.tolist().index(l)]*g2.CARD["g"]**2/(32*np.pi*g2.phi_mass(1)*g2.detector_mass(l))
    bound *= tail*(1+DRIVE["epsilon"]**2/2 if driven else 1)
    result.update(sigma_i=sigma_i, sigma_d=sigma_d, driven=driven,
                  incident_tail_probability=tail, omitted_pair_rate_bound=float(bound),
                  normalized_incident_tail_error_bound=float(2*bound/result["total_pair_rate"]),
                  actual_incident_mean=float(np.dot(iw, p)/iw.sum()),
                  actual_incident_sd=float(np.sqrt(np.dot(iw, (p-np.dot(iw,p)/iw.sum())**2)/iw.sum())))
    return result


def run():
    checks = []

    def check(name, passed, measured=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if measured is not None: row["measured"] = float(measured)
        if tolerance is not None: row["tolerance"] = float(tolerance)
        checks.append(row)

    joint_path = ROOT/"output"/"g2_recoil_joint.json"
    joint = json.loads(joint_path.read_text())
    g2_path = Path(g2.__file__)
    check("G2-R:source_hash", joint["verification"]["g2_source_sha256"] == hashlib.sha256(g2_path.read_bytes()).hexdigest())
    check("G2-R:artifact_hash", joint["verification"]["g2_artifact_sha256"] == hashlib.sha256((ROOT/"output"/"g2_local_conversion.json").read_bytes()).hexdigest())
    for si in SIGMA_I:
        ps, iw, tail = incoming_quadrature(si)
        check(f"incoming:quadrature_normalization:{si}", abs(iw.sum()+tail-1) < 1e-14)
        check(f"incoming:nonnegative_support_weights:{si}", np.all(ps >= 0) and np.all(iw > 0))
        if si:
            a = -g2.CARD["incoming_com_momentum"]/si
            ratio = np.exp(-a*a/2)/(SQRT2PI*ndtr(-a))
            exact_mean = g2.CARD["incoming_com_momentum"]+si*ratio
            exact_var = si*si*(1+a*ratio-ratio*ratio)
            mean = np.dot(iw, ps)
            variance = np.dot(iw, (ps-mean)**2)
            check(f"incoming:analytic_truncated_mean:{si}", abs(mean-exact_mean) < 1e-13)
            check(f"incoming:analytic_truncated_variance:{si}", abs(variance-exact_var) < 1e-13)
    rows, _, _, _ = components(0, False)
    original = [r for r in joint["events"] if r["m"] == 0 and abs(r["j"]) == 1]
    for row in rows:
        old = next(r for r in original if r["j"] == row["sign"])
        check(f"sharp:momentum:j={row['sign']}", abs(row["means"][0]-old["outgoing_momentum"]) < 2e-14)
        check(f"sharp:weighted_rate:j={row['sign']}", abs(row["rates"][0]/old["weighted_rate"]-1) < 2e-13)
    n, m, p0 = 1, 0, g2.CARD["incoming_com_momentum"]
    slopes = []
    for l in [0, -2]:
        pf = kernel(np.array([p0]), l)[0][0]
        en, el = np.hypot(g2.phi_mass(n), p0), np.hypot(g2.detector_mass(l), p0)
        ef, er = np.hypot(g2.phi_mass(m), pf), np.hypot(g2.detector_mass(n+l-m), pf)
        slope = p0*(1/en+1/el)/(pf*(1/ef+1/er))
        h = 1e-5
        finite = (kernel(np.array([p0+h]), l)[0][0]-kernel(np.array([p0-h]), l)[0][0])/(2*h)
        check(f"dp_f_dp:l={l}", abs(slope-finite) < 2e-9, abs(slope-finite), 2e-9)
        slopes.append(dict(l=l, j=n+l-m, derivative=float(slope), finite_difference=float(finite)))

    scans = []
    for driven in [False, True]:
        for si in SIGMA_I:
            previous = -1.0
            for sd in SIGMA_D:
                row = evaluate(si, sd, driven)
                coarse = evaluate(si, sd, driven, order=24, grid_divisor=8)
                delta = abs(row["bayes_error"]-coarse["bayes_error"])
                tolerance = 3e-9*max(row["bayes_error"], 1e-14)
                check(f"convergence:drive={driven}:si={si}:sd={sd}", delta < tolerance, delta, tolerance)
                row["quadrature_order_doubling_error"] = delta
                check(f"positive_and_bounded:{driven}:{si}:{sd}", 0 <= row["bayes_error"] <= row["majority_only_error"]+2e-12)
                check(f"normalized:{driven}:{si}:{sd}", abs(sum(row["priors"])-1) < 1e-13)
                check(f"noise_data_processing:{driven}:{si}:{sd}", row["bayes_error"] >= previous-1e-12)
                previous = row["bayes_error"]
                if driven:
                    check(f"q_tag_data_processing:{si}:{sd}", row["ideal_q_tagged_error"] <= row["bayes_error"]+1e-12)
                if si == 0 and not driven:
                    a, b = row["priors"]
                    ma, mb = [s["mean_momentum"] for s in row["sign_moments"]]
                    boundary = (ma+mb)/2+sd**2/(mb-ma)*np.log(a/b)
                    exact = a*ndtr((ma-boundary)/sd)+b*ndtr((boundary-mb)/sd)
                    row["analytic_gaussian_error"] = float(exact)
                    relative = abs(row["bayes_error"]-exact)/exact
                    check(f"sharp:analytic_gaussian:sd={sd}", relative < 1e-6, relative, 1e-6)
                scans.append(row)

    for idx, slope in enumerate(slopes):
        narrow = next(r for r in scans if not r["driven"] and r["sigma_i"] == .02 and r["sigma_d"] == .02)
        actual = narrow["sign_moments"][idx]["intrinsic_variance"]
        linear = (slope["derivative"]*.02)**2
        relative = abs(actual/linear-1)
        check(f"linearized_width:j={slope['j']}", relative < .03, relative, .03)
        slope.update(sigma_i=.02, exact_rate_weighted_variance=actual, linear_variance=linear,
                     relative_linearization_error=relative)

    profiles = []
    for driven in [False, True]:
        cr, _, _, _ = components(.1, driven)
        z = sum(float(r["rates"].sum()) for r in cr)
        x = np.linspace(.4, 2.1, 341)
        profiles.append(dict(sigma_i=.1, sigma_d=.1, driven=driven, x=x.tolist(),
                             plus_density=(sum(density(x, r["means"], r["rates"], .1) for r in cr if r["sign"] == 1)/z).tolist(),
                             minus_density=(sum(density(x, r["means"], r["rates"], .1) for r in cr if r["sign"] == -1)/z).tolist(),
                             note="Plot-ready samples only; all reported integrals use CDF boundary integration, not this finite profile window."))

    # Floquet q=0 must retain the original local kernel exactly.
    ptest = np.array([0.0, .2, .7, 2.0])
    for l in [0, -2]:
        undriven, driven = kernel(ptest, l), kernel(ptest, l, 0, True)
        check(f"Floquet:q0_baseline:l={l}", all(np.array_equal(a, b) for a, b in zip(undriven, driven)))
        for q in [-1, 0, 1]:
            pf, rate = kernel(ptest, l, q, True)
            energy = np.hypot(g2.phi_mass(0), pf)+np.hypot(g2.detector_mass(1+l), pf)
            work = energy-np.hypot(g2.phi_mass(1), ptest)-np.hypot(g2.detector_mass(l), ptest)
            check(f"Floquet:work:l={l}:q={q}", np.max(abs(work-q*DRIVE["Omega"])) < 3e-14)
            check(f"Floquet:nonnegative:l={l}:q={q}", np.all(rate >= 0))

    failed = [r["name"] for r in checks if not r["passed"]]
    result = dict(status="Conditional idealized finite-resolution feasibility forecast; no absolute collision probability or empirical calibration.",
                  card=g2.CARD, compact_packet=g2.PACKET, drive=DRIVE,
                  readout_card=dict(selected_mode=0, selected_abs_recoil=1, signed_labels=[1,-1], incoming_labels=[0,-2],
                                    sigma_i=SIGMA_I, sigma_d=SIGMA_D, error_threshold=RESOLVED_ERROR,
                                    incoming="Truncated normal on p>=0 centered before truncation at p0=.2. sigma_i is underlying normal width, not truncated SD; zero width means delta.",
                                    measured="Gaussian reconstructed momentum estimator x in R, not a nonnegative exact norm; perfect m,|j| labeling and equal acceptance.",
                                    geometry="Common fixed COM frame and direction-independent response; no angular/boost spread or coherent 4D wave packet.",
                                    drive="Homogeneous periodic coupling g(t)=g[1+epsilon cos(Omega t)]; tree, weak scattering, long-time resolved Floquet sidebands. q is unobserved unless explicitly ideal-tagged.",
                                    truncation="Incident standardized z<=10 retained without renormalization; explicit rate bound included. Bayes grid spans every component mean +/-10 sigma_d; omitted response mass <=2 Phi(-10).",
                                    EFT_limit="The incident-tail rate bound assumes formal continuation of the tree kernel. No EFT validity scale, loop error, or UV cross-section bound has been specified; small numerical tails do not certify physical Gaussian support to infinite momentum."),
                  equations=dict(component="A_s(x)=w_l sum_q integral_0^infty dp f_i(p) K_lmq(p) R(x|pf_lmq(p)); l=0 for s=+, -2 for s=-.",
                                 rate="K_lmq=|g_q|^2 pf/[16*pi*E_n*E_l*(E_n+E_l+q Omega)], open channel only.",
                                 bayes="e_B=integral dx min(A_+,A_-)/integral dx(A_++A_-)",
                                 tagged="e_B,tag=sum_q integral dx min(A_+,q,A_-,q)/Z <= e_B,unobserved.",
                                 slope="dp_f/dp = p(1/E_n+1/E_l)/[pf(1/E_f+1/E_r)]",
                                 work="E_final-E_initial=q Omega; compact momentum remains n+l=m+j."),
                  slopes=slopes, scan=scans, profiles=profiles, checks=checks,
                  summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
                  verification=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                    g2_source_sha256=hashlib.sha256(g2_path.read_bytes()).hexdigest(),
                                    g2_artifact_sha256=hashlib.sha256((ROOT/"output"/"g2_local_conversion.json").read_bytes()).hexdigest(),
                                    g2r_artifact_sha256=hashlib.sha256(joint_path.read_bytes()).hexdigest()))
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Route G2-S: resolution, incoming spread and powered-sideband comparison", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**.", "",
             "This is a conditional recoil-sign classification forecast for selected scattered events with m=0 and |j|=1, not the probability of a collision. Units are the unchanged G2 card. No measured detector performance is assumed.", "",
             "The incident ensemble has p>=0 and density exp[-(p-.2)^2/(2 sigma_i^2)]/[sqrt(2 pi) sigma_i Phi(.2/sigma_i)]. sigma_i is the pre-truncation width. The response is a Gaussian real-valued reconstructed momentum estimator with width sigma_d. Both compact labels are otherwise perfectly resolved, efficiency is identical, and all components share a COM frame and dilute external overlap.", "",
             "Rates are integrated before normalization: A_s(x)=w_l sum_q integral f_i K_lmq R dp; the optimal error is integral min(A_+,A_-) dx / Z. The 5% boundary is an explicitly chosen engineering criterion, not a theorem or confidence level.", "",
             "## Optimal unobserved-sign errors", "",
             "Each cell is an error percentage; **bold** means <=5%. Incident spread changes the event-selected priors as well as the line shape.", ""]
    for powered in [False, True]:
        lines += ["### "+("Powered: epsilon=.5, Omega=.2; q unobserved" if powered else "Stationary G2 baseline"), "",
                  "| sigma_i \\ sigma_d | "+" | ".join(map(str, SIGMA_D))+" |",
                  "|---|"+"---:|"*len(SIGMA_D)]
        for si in SIGMA_I:
            values = [r for r in scans if r["driven"] == powered and r["sigma_i"] == si]
            def format_error(r):
                value = f"{100*r['bayes_error']:.6g}%"
                return "**"+value+"**" if r["resolved_at_declared_5pct"] else value
            lines.append("| "+str(si)+" | "+" | ".join(format_error(r) for r in values)+" |")
        lines.append("")
    lines += ["## Selected comparison and event priors", "",
              "| sigma_i | sigma_d | stationary error | powered error | powered ideal q-tag error | stationary majority-only error |",
              "|---:|---:|---:|---:|---:|---:|"]
    for si, sd in [(0,.1),(.1,.1),(.2,.1),(.4,.1),(.1,.2)]:
        a = next(r for r in scans if not r["driven"] and r["sigma_i"] == si and r["sigma_d"] == sd)
        b = next(r for r in scans if r["driven"] and r["sigma_i"] == si and r["sigma_d"] == sd)
        vals = [a["bayes_error"], b["bayes_error"], b["ideal_q_tagged_error"], a["majority_only_error"]]
        lines.append(f"| {si} | {sd} | "+" | ".join(f"{100*v:.6g}%" for v in vals)+" |")
    lines += ["", "The ideal q-tag is a mathematical information benchmark, not a built pump detector. A real tag must include its noise and work readout. The drive supplies/removes q Omega energy, does not alter compact selection n+l=m+j, and produces additional momentum lines rather than a coordinate change of mass.", "",
              "## Numerical controls and limitations", "",
              "- Piecewise Gauss-Legendre orders 24 and 48 and Bayes-boundary grids sigma_d/8 and sigma_d/4 are compared for every scan point. Gaussian component CDF differences are evaluated with survival tails to avoid cancellation; no error floor is imposed.",
              "- Incident z>10 is omitted without reweighting; Kq<=|gq|^2/(32 pi m_n mu_l) bounds the omitted rate. The response-window loss is <=2 Phi(-10). Bounds and convergence differences are stored per point.",
              "- Grid/order agreement is a numerical convergence test, not certified exhaustive isolation of all mixture roots or a rigorous global error bar. The analytic tail bounds do not certify that separate boundary search.",
              "- Tests cover exact sharp-beam G2-R momenta, analytic unequal-prior Gaussian Bayes errors, kinematic derivatives, narrow-spread widths, positivity, normalized priors, noise data processing, and the advantage of ideal sideband information.",
              "- JSON records both false-sign rates, priors, decision boundaries, selected momentum moments, all scans, and representative profile samples. The finite sampled profiles are not used to compute errors.",
              "- Unknown boosts, angular spread, mass-label confusion, unequal acceptance, correlated recoil measurement, finite pump duration and coherent packet interference are outside this card. The powered rate limit presumes resolved sidebands; no finite-time absolute transition probability is claimed.",
              "- The small Gaussian-tail bounds certify numerical integration of the declared tree kernel, not the physical UV validity of that kernel. The EFT cutoff and loop errors are unspecified; an unbounded physical Gaussian beam is not validated by this calculation. A physical beam card must be restricted to a justified EFT domain or supplemented by UV control.",
              "", "Next: replace this engineering response card with specified or calibrated joint momentum/energy response. Specify four-dimensional incoming packets, timing, geometry and pump envelope only if an absolute collision probability is required.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"]))
    for si, sd in [(0,.1),(.1,.1),(.2,.1),(.4,.1)]:
        selected = [r for r in scans if r["sigma_i"] == si and r["sigma_d"] == sd]
        print(json.dumps([dict(sigma_i=si, sigma_d=sd, driven=r["driven"], error=r["bayes_error"], tagged=r.get("ideal_q_tagged_error")) for r in selected]))
    if failed:
        raise SystemExit("Failed checks: "+", ".join(failed))


if __name__ == "__main__":
    run()
