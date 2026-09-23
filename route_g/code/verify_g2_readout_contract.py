#!/usr/bin/env python3
"""Frozen momentum-only engineering readout contract on unchanged G2-S rates.

No detector calibration, physical fit, absolute collision probability, or new
interaction is supplied. Oracle retuning is reported separately from deployment.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr

import verify_g2_resolution as resolution

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_readout_contract"
CARD = dict(sigma_i=.1, nominal_sigma_d=.05, driven=True,
            scale_errors=[-.01, 0.0, .01], offset_errors=[-.02, 0.0, .02],
            true_sigmas=[.04, .05, .06], chosen_error_ceiling=.05)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def transform(rows, scale_error=0.0, offset_error=0.0):
    return [dict(row, means=(1+scale_error)*row["means"]+offset_error) for row in rows]


def class_arrays(rows):
    means = [np.concatenate([r["means"] for r in rows if r["sign"] == s]) for s in [1, -1]]
    rates = [np.concatenate([r["rates"] for r in rows if r["sign"] == s]) for s in [1, -1]]
    z = sum(float(r.sum()) for r in rates)
    return means, [r/z for r in rates]


def nominal_rule(rows, sigma, grid_divisor=8):
    fitted = resolution.classify(rows, sigma, grid_divisor)
    boundaries = fitted["bayes_boundaries"]
    means, weights = class_arrays(rows)
    edges = [-np.inf]+boundaries+[np.inf]
    decisions = []
    for left, right in zip(edges[:-1], edges[1:]):
        if np.isneginf(left):
            x = right-10*sigma
        elif np.isposinf(right):
            x = left+10*sigma
        else:
            x = (left+right)/2
        values = [float(resolution.density(x, m, w, sigma)) for m, w in zip(means, weights)]
        decisions.append(1 if values[0] >= values[1] else -1)
    return dict(boundaries=boundaries, interval_predictions=decisions,
                interval_convention="[-infinity,t1), [t1,t2), ..., [tk,infinity); ties have zero Gaussian probability",
                origin="Numerically resolved nominal Bayes boundaries; frozen for every calibration stress case.")


def assess_rule(rows, sigma, rule):
    means, weights = class_arrays(rows)
    matrix = np.zeros((2, 2))
    edges = [-np.inf]+rule["boundaries"]+[np.inf]
    for left, right, prediction in zip(edges[:-1], edges[1:], rule["interval_predictions"]):
        predicted = 0 if prediction == 1 else 1
        for truth in [0, 1]:
            matrix[truth, predicted] += resolution.interval_mass(left, right, means[truth], weights[truth], sigma)
    priors, predicted = matrix.sum(axis=1), matrix.sum(axis=0)
    recall = np.diag(matrix)/priors
    purity = np.diag(matrix)/predicted
    return dict(confusion_joint=matrix.tolist(), confusion_rows_truth_columns_prediction=[1, -1],
                priors=priors.tolist(), predicted_fractions=predicted.tolist(),
                conditional_confusion=(matrix/priors[:, None]).tolist(),
                recall=recall.tolist(), purity=purity.tolist(),
                overall_error=float(matrix[0, 1]+matrix[1, 0]),
                false_minus_given_plus=float(matrix[0, 1]/priors[0]),
                false_plus_given_minus=float(matrix[1, 0]/priors[1]))


def continuous_box_bound(rows, threshold, box=None, epsabs=1e-12):
    """Analytic componentwise envelope; numerical integration is not interval-certified."""
    means, weights = class_arrays(rows)
    box = box or dict(a=[min(CARD["scale_errors"]), max(CARD["scale_errors"])],
                      b=[min(CARD["offset_errors"]), max(CARD["offset_errors"])],
                      sigma=[min(CARD["true_sigmas"]), max(CARD["true_sigmas"])])
    amin, amax = box["a"]
    bmin, bmax = box["b"]
    smin, smax = box["sigma"]
    # pf >= 0: largest mean maximizes plus->minus, smallest maximizes minus->plus.
    delta_plus = (1+amax)*means[0]+bmax-threshold
    delta_minus = threshold-((1+amin)*means[1]+bmin)
    tails = [np.maximum(ndtr(delta/smin), ndtr(delta/smax)) for delta in [delta_plus, delta_minus]]
    sampled_bounds = [float(np.dot(w, tail)) for w, tail in zip(weights, tails)]
    # max over sigma has a kink at delta=0. Split there rather than mistaking
    # smooth-beam quadrature convergence for convergence of this new envelope.
    bounds = [0.0, 0.0]
    estimated_error = 0.0
    split_points = []
    packet = resolution.g2.PACKET
    labels = np.arange(packet["min_l"], packet["max_l"]+1)
    compact_weights = np.exp(-labels**2/(2*packet["width"]**2))
    compact_weights /= compact_weights.sum()
    ztotal = sum(float(row["rates"].sum()) for row in rows)
    p0, si = resolution.g2.CARD["incoming_com_momentum"], CARD["sigma_i"]
    lower, upper = -p0/si, resolution.Z_CUT
    for row in rows:
        index = 0 if row["sign"] == 1 else 1
        target_pf = (threshold-bmax)/(1+amax) if index == 0 else (threshold-bmin)/(1+amin)

        def target(z):
            return float(resolution.kernel(np.array([p0+si*z]), row["l"], row["q"], True)[0][0])-target_pf

        points = [brentq(target, lower, upper)] if target(lower)*target(upper) < 0 else []
        split_points.append(dict(l=row["l"], q=row["q"], kink_z=points))
        wl = compact_weights[labels.tolist().index(row["l"])]

        def integrand(z):
            pf, rate = resolution.kernel(np.array([p0+si*z]), row["l"], row["q"], True)
            delta = ((1+amax)*pf[0]+bmax-threshold if index == 0
                     else threshold-((1+amin)*pf[0]+bmin))
            tail = max(ndtr(delta/smin), ndtr(delta/smax))
            fz = np.exp(-z*z/2)/(resolution.SQRT2PI*ndtr(p0/si))
            return float(wl*fz*rate[0]*tail/ztotal)

        integral, error = quad(integrand, lower, upper, points=points, epsabs=epsabs, epsrel=1e-11)
        bounds[index] += integral
        estimated_error += error
    return dict(box=box, overall_error_upper_envelope=float(sum(bounds)),
                wrong_joint_mass_upper_envelope=bounds,
                conditional_wrong_rate_upper_envelope=[v/float(w.sum()) for v, w in zip(bounds, weights)],
                kink_aware_quadrature_error_estimate=estimated_error,
                kink_splits=split_points, unsplit_beam_quadrature_value=float(sum(sampled_bounds)),
                formula="U=sum_plus weights max_{sigma=smin,smax} Phi(((1+amax)pf+bmax-t)/sigma) + sum_minus weights max_{sigma=smin,smax} Phi((t-(1+amin)pf-bmin)/sigma)",
                interpretation="For each p,q component, monotonicity in mean and endpoint maximization in sigma gives a pointwise bound on every fixed calibration point in the continuous box. Summing separately maximized components is conservative; no single calibration must attain U.",
                numerical_status="The inequality is analytic; its displayed integral uses kink-split adaptive floating-point beam quadrature, not certified interval arithmetic or a global numerical proof.")


def run():
    checks = []

    def check(name, passed, value=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if value is not None: row["measured"] = float(value)
        if tolerance is not None: row["tolerance"] = float(tolerance)
        checks.append(row)

    parent_path = ROOT/"output"/"g2_resolution.json"
    parent = json.loads(parent_path.read_text())
    deps = dict(g2s_source=Path(resolution.__file__), g2s_artifact=parent_path,
                g2_source=ROOT/"code"/"verify_g2_local_conversion.py",
                g2_artifact=ROOT/"output"/"g2_local_conversion.json",
                g2r_artifact=ROOT/"output"/"g2_recoil_joint.json")
    before_hashes = {key: sha(path) for key, path in deps.items()}
    check("dependency:G2-S_source", before_hashes["g2s_source"] == parent["verification"]["source_sha256"])
    for key in ["g2_source", "g2_artifact", "g2r_artifact"]:
        check("dependency:"+key, before_hashes[key] == parent["verification"][key+"_sha256"])
    check("dependency:G2-S_all_checks_passed", not parent["summary"]["failed"])
    rows, _, _, _ = resolution.components(CARD["sigma_i"], True)
    coarse_rows, _, _, _ = resolution.components(CARD["sigma_i"], True, order=24)
    rule = nominal_rule(rows, CARD["nominal_sigma_d"])
    coarse_rule = nominal_rule(coarse_rows, CARD["nominal_sigma_d"], grid_divisor=16)
    check("nominal:single_threshold_verified", len(rule["boundaries"]) == 1 and rule["interval_predictions"] == [1, -1])
    if len(rule["boundaries"]) != 1 or rule["interval_predictions"] != [1, -1]:
        raise RuntimeError("This card does not admit the declared one-threshold bound")
    t = rule["boundaries"][0]
    check("nominal:boundary_grid_and_quadrature_convergence", len(coarse_rule["boundaries"]) == 1 and abs(t-coarse_rule["boundaries"][0]) < 1e-10)
    nominal = assess_rule(rows, CARD["nominal_sigma_d"], rule)
    old = next(r for r in parent["scan"] if r["sigma_i"] == CARD["sigma_i"] and r["sigma_d"] == CARD["nominal_sigma_d"] and r["driven"])
    check("nominal:parent_error_recovery", abs(nominal["overall_error"]-old["bayes_error"]) < 1e-13)
    check("nominal:parent_false_minus_recovery", abs(nominal["false_minus_given_plus"]-old["false_minus_given_plus"]) < 1e-13)
    check("nominal:parent_false_plus_recovery", abs(nominal["false_plus_given_minus"]-old["false_plus_given_minus"]) < 1e-13)
    grid = []
    for a, b, sigma in itertools.product(CARD["scale_errors"], CARD["offset_errors"], CARD["true_sigmas"]):
        changed = transform(rows, a, b)
        frozen = assess_rule(changed, sigma, rule)
        oracle_rule = nominal_rule(changed, sigma)
        oracle = assess_rule(changed, sigma, oracle_rule)
        coarse = assess_rule(transform(coarse_rows, a, b), sigma, rule)
        error = abs(frozen["overall_error"]-coarse["overall_error"])
        label = f"a={a}:b={b}:sigma={sigma}"
        check("stress:normalization:"+label, abs(np.asarray(frozen["confusion_joint"]).sum()-1) < 1e-13)
        check("stress:nonnegative_confusion:"+label, np.min(frozen["confusion_joint"]) >= 0)
        check("stress:oracle_not_worse:"+label, oracle["overall_error"] <= frozen["overall_error"]+1e-13)
        affine_oracle = resolution.classify(rows, sigma/(1+a), grid_divisor=8)["bayes_error"]
        check("stress:oracle_affine_identity:"+label, abs(oracle["overall_error"]-affine_oracle) < 1e-12)
        check("stress:quadrature:"+label, error < 1e-12, error, 1e-12)
        grid.append(dict(scale_error=a, offset_error=b, sigma_true=sigma,
                         frozen=frozen, reoptimized_oracle=oracle,
                         oracle_rule=oracle_rule, quadrature_order_doubling_error=error))
    worst = max(grid, key=lambda row: row["frozen"]["overall_error"])
    bound = continuous_box_bound(rows, t)
    coarse_bound = continuous_box_bound(coarse_rows, t, epsabs=1e-10)
    bound["quadrature_refinement_difference"] = abs(bound["overall_error_upper_envelope"]-coarse_bound["overall_error_upper_envelope"])
    bound["formal_tree_incident_tail_error_bound"] = old["normalized_incident_tail_error_bound"]
    check("box:all_sampled_points_below_component_envelope", max(r["frozen"]["overall_error"] for r in grid) <= bound["overall_error_upper_envelope"]+1e-13)
    check("box:quadrature_convergence", bound["quadrature_refinement_difference"] < 1e-10)
    subboxes = []
    for a, b, sigma in itertools.product([[-.01, 0.0], [0.0, .01]], [[-.02, 0.0], [0.0, .02]], [[.04, .05], [.05, .06]]):
        sub = continuous_box_bound(rows, t, dict(a=a, b=b, sigma=sigma))
        subboxes.append(sub)
    partitioned_bound = max(sub["overall_error_upper_envelope"] for sub in subboxes)
    check("box:partition_no_worse_than_unpartitioned", partitioned_bound <= bound["overall_error_upper_envelope"]+1e-12)
    check("box:partition_above_sampled_risks", partitioned_bound >= worst["frozen"]["overall_error"]-1e-12)
    partition = dict(predeclared_partition="One equal 2x2x2 split of shared a,b,sigma intervals; no adaptive refinement or fitted tolerances.",
                     subboxes=subboxes, overall_error_upper_envelope=partitioned_bound,
                     below_chosen_5pct_ceiling=partitioned_bound < CARD["chosen_error_ceiling"],
                     interpretation="Every continuous calibration point belongs to a subbox. Each subbox has a componentwise analytic upper envelope; max over all eight covers the original box. Numerical integrals retain the same non-interval-certified qualification.")
    # The nominal Bayes error is nondecreasing under extra independent Gaussian noise.
    def optimized_error(sigma):
        return resolution.classify(rows, sigma, grid_divisor=8)["bayes_error"]

    ceiling = CARD["chosen_error_ceiling"]
    sigma_limit = brentq(lambda sigma: optimized_error(sigma)-ceiling, .05, .1, xtol=1e-10)
    frozen_limit = brentq(lambda sigma: assess_rule(rows, sigma, rule)["overall_error"]-ceiling, .05, .1, xtol=1e-10)
    check("resolution_limit:root", abs(optimized_error(sigma_limit)-ceiling) < 1e-9)
    check("resolution_limit:frozen_not_better_than_oracle", frozen_limit <= sigma_limit+1e-9)
    check("resolution_limit:bracket", optimized_error(.05) < ceiling < optimized_error(.1))
    check("resolution_limit:local_monotonicity", optimized_error(sigma_limit-.001) < ceiling < optimized_error(sigma_limit+.001))
    check("dependency:no_parent_artifact_changed", all(sha(path) == before_hashes[key] for key, path in deps.items()))
    failed = [row["name"] for row in checks if not row["passed"]]
    result = dict(status="Frozen conditional momentum-readout engineering contract, not calibrated hardware or an experimental claim.",
                  contract=CARD, drive=resolution.DRIVE, rule=rule, nominal=nominal,
                  calibration_response="x=(1+a) pf+b+Normal(0,sigma_true^2); all rates and incident distributions stay unchanged.",
                  assumptions=["Only selected scattered events with m=0 and |j|=1; both labels and acceptance are ideal.",
                               "sigma_i=.1 is the underlying Gaussian width before p>=0 truncation; the incoherent beam has common COM and external overlap.",
                               "Sideband q is not measured; all three Floquet rates are summed before classification.",
                               "The decision rule is trained once at a=b=0,sigma=.05 and remains fixed throughout the stress card.",
                               "Oracle results retune the classifier using true calibration; they are not deployed performance.",
                               "The 5% criterion concerns overall event-weighted error. It does not require each class error <=5%; nominal minority j=-1 error is about6.23%.",
                               "Scale, offset and momentum widths are in arbitrary model units; none is a hardware calibration or confidence interval.",
                               "Beam quadrature retains standardized z<=10; its omitted-rate bound assumes formal tree-kernel continuation. EFT cutoff, loop error and unknown UV tails remain unspecified.",
                               "Decision intervals extend over all real estimator values, so Gaussian response tails are integrated by exact component CDFs. Numerical root completeness is tested by refined scans, not formally certified."],
                  stress_grid=grid, grid_worst_frozen=worst,
                  worst_reoptimized_oracle=max(grid, key=lambda row: row["reoptimized_oracle"]["overall_error"]),
                  continuous_calibration_box=bound,
                  partitioned_continuous_calibration_box=partition,
                  nominal_resolution_limit=dict(optimal_retrained_sigma_max=sigma_limit,
                                                frozen_nominal_rule_sigma_max=frozen_limit,
                                                error_ceiling=ceiling,
                                                scope="Only a=b=0, unchanged beam and drive. Retrained Bayes limit is not a frozen-rule or calibration-box guarantee."),
                  checks=checks, summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
                  verification=dict(source_sha256=sha(__file__), dependencies=before_hashes))
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    c = nominal["conditional_confusion"]
    lines = ["# G2 readout contract: frozen momentum-only rule", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**. Engineering forecast only; no calibrated hardware or absolute collision probability.", "",
             "The main card keeps the existing driven action (epsilon=.5, Omega=.2), sigma_i=.1, perfect m=0 and |j|=1 selection, and a Gaussian real-valued momentum estimator with nominal sigma_d=.05. Sideband q is unobserved. All quantities use unchanged arbitrary model units.", "",
             "## Deployable frozen decision", "",
             f"Predict j=+1 when x < **{t:.6f}**; predict j=-1 otherwise. One boundary is found in both ordinary and refined searches. The rule is trained once at a=b=0 and sigma=.05.", "",
             f"Nominal overall error: **{100*nominal['overall_error']:.4f}%**. Majority-only error: {100*min(nominal['priors']):.4f}%.", "",
             "| True sign | Predict + | Predict - | Recall | Purity of this predicted label |",
             "|---|---:|---:|---:|---:|"]
    for i, sign in enumerate([1, -1]):
        lines.append(f"| {sign:+d} | {100*c[i][0]:.4f}% | {100*c[i][1]:.4f}% | {100*nominal['recall'][i]:.4f}% | {100*nominal['purity'][i]:.4f}% |")
    lines += ["", "The two middle columns are conditional on the true sign. Purity instead conditions on the predicted sign; it uses the event-selected prior and is not the same as recall. The chosen 5% ceiling applies to the overall event-weighted error, not each class: the nominal j=-1 error is 6.23%, so this contract does not claim 95% recall for both signs.", "",
              "## Calibration stress card", "",
              "True response: x=(1+a)pf+b+Normal(0,sigma_true^2), with a in {-0.01,0,0.01}, b in {-0.02,0,0.02}, sigma_true in {0.04,0.05,0.06}. The frozen classifier is never adjusted to these true values. The oracle is retrained and reported only as a comparison.", "",
              "| a | b | sigma_true | Frozen error | Retuned oracle error |",
              "|---:|---:|---:|---:|---:|"]
    for row in grid:
        lines.append(f"| {row['scale_error']:+.2f} | {row['offset_error']:+.2f} | {row['sigma_true']:.2f} | {100*row['frozen']['overall_error']:.4f}% | {100*row['reoptimized_oracle']['overall_error']:.4f}% |")
    lines += ["", f"Worst sampled frozen error: **{100*worst['frozen']['overall_error']:.4f}%**, at a={worst['scale_error']:+.2f}, b={worst['offset_error']:+.2f}, sigma_true={worst['sigma_true']:.2f}. Its oracle error is {100*worst['reoptimized_oracle']['overall_error']:.4f}%.", "",
              "## Continuous-box bound and numerical qualification", "",
              "For the fixed threshold t, plus-to-minus error is Phi((mu-t)/sigma), maximized at mu=(1+a_max)pf+b_max. Minus-to-plus error is Phi((t-mu)/sigma), maximized at mu=(1+a_min)pf+b_min. For either signed difference d, max_{sigma in [s_min,s_max]} Phi(d/sigma) occurs at an endpoint. Sum these componentwise maxima using the unchanged normalized rate weights.", "",
              f"The resulting conservative continuous-box envelope is **{100*bound['overall_error_upper_envelope']:.4f}%**. Unlike a grid maximum, the inequality covers all fixed calibration points in the stated continuous box. Different components may maximize at different calibrations, so this envelope need not be attainable.", "",
              f"A single predeclared equal 2x2x2 partition of the shared a,b,sigma box reduces the covering envelope to **{100*partitioned_bound:.4f}%**. This maximum of eight subbox bounds differs both from the loose global envelope and from the 27-point sampled maximum; it does not retune the classifier.", "",
              f"Each bound integral splits the kinks where component mean crosses threshold. Refined adaptive integration and beam-normalization difference {bound['quadrature_refinement_difference']:.2g}; formal tree-model omitted-beam-tail error bound {bound['formal_tree_incident_tail_error_bound']:.3g}. This is a converged numerical value of an analytic envelope, not interval-certified arithmetic or a rigorous global numerical proof. No EFT cutoff or UV completion is inferred from the tiny Gaussian tail.", "",
              "## Available nominal resolution", "",
              f"At exact a=b=0, the 5% optimal-Bayes ceiling crosses near sigma_d={sigma_limit:.4f}; a rounded-down requirement is **sigma_d <= {np.floor(1000*sigma_limit)/1000:.3f}** if the classifier is retrained for that resolution. Keeping the actual frozen rule crosses near {frozen_limit:.4f}, with rounded-down requirement **sigma_d <= {np.floor(1000*frozen_limit)/1000:.3f}**. Neither number applies unchanged to unknown calibration errors or altered beam/drive parameters.", "",
              "Next: specify or calibrate a joint momentum response and mass-label confusion/acceptance. Treat finite noisy sideband tagging as a separate comparison; it is not used to rescue this momentum-only result. Four-dimensional packet geometry, pump timing and an EFT-valid beam domain are required before an absolute collision probability is meaningful.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"]))
    print(json.dumps(dict(threshold=t, nominal=nominal, grid_worst=dict(a=worst["scale_error"], b=worst["offset_error"], sigma=worst["sigma_true"], error=worst["frozen"]["overall_error"]), bound=bound, resolution_limit=result["nominal_resolution_limit"])))
    if failed:
        raise SystemExit("Failed checks: "+", ".join(failed))


if __name__ == "__main__":
    run()
