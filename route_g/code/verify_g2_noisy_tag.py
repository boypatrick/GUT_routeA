#!/usr/bin/env python3
"""G2-M comparison: a postulated noisy ternary sideband tag, not apparatus.

Reuse the unchanged rate-weighted G2-S densities. Only an independent
classical response T_eta(z|q) is added; no new scattering or wavepackets.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import verify_g2_resolution as resolution

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_noisy_tag"
TAGS = (-1, 0, 1)
ETAS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
SIGMA_I = 0.1
SIGMA_D = (0.1, 0.05)
TARGET_ERROR = 0.05
THRESHOLD_BRACKET_WIDTH = 0.001


def response_matrix(eta):
    if not 0 <= eta <= 1:
        raise ValueError("eta is a stochastic mixing strength in [0,1]")
    return (1-eta)*np.eye(3)+eta*np.ones((3,3))/3


def tag_rows(rows, eta, tag):
    matrix = response_matrix(eta)
    return [dict(row, rates=row["rates"]*matrix[TAGS.index(tag),TAGS.index(row["q"])])
            for row in rows]


def evaluate_tag(eta, sigma_d, order=48, grid_divisor=4):
    rows, _, _, tail = resolution.components(SIGMA_I, True, order)
    total = sum(float(row["rates"].sum()) for row in rows)
    tagged = []
    for tag in TAGS:
        result = resolution.classify(tag_rows(rows,eta,tag),sigma_d,grid_divisor)
        probability = result["total_pair_rate"]/total
        tagged.append(dict(z=tag,label_probability=probability,
                           conditional_sign_priors=result["priors"],
                           conditional_bayes_error=result["bayes_error"],
                           contribution_to_total_error=probability*result["bayes_error"],
                           bayes_boundaries=result["bayes_boundaries"],
                           false_minus_given_plus=result["false_minus_given_plus"],
                           false_plus_given_minus=result["false_plus_given_minus"],
                           pair_rate_for_label=result["total_pair_rate"]))
    error=sum(row["contribution_to_total_error"] for row in tagged)
    priors=[sum(row["label_probability"]*row["conditional_sign_priors"][index]
                for row in tagged) for index in (0,1)]
    return dict(eta=eta,actual_wrong_tag_probability=2*eta/3,
                sigma_i=SIGMA_I,sigma_d=sigma_d,response=response_matrix(eta).tolist(),
                bayes_error=error,passes_declared_5pct=error<=TARGET_ERROR,
                total_pair_rate=total,sign_priors=priors,label_rows=tagged,
                incident_tail_probability=tail,
                response_tail_error_bound=float(2*resolution.ndtr(-resolution.Z_CUT)))


def direct_density_integral(eta,sigma_d,order=48):
    """Independent adaptive x integration; do not reuse boundary/CDF masses."""
    rows,_,_,_=resolution.components(SIGMA_I,True,order)
    total=sum(float(row["rates"].sum()) for row in rows)
    value,error=0.0,0.0
    for tag in TAGS:
        transformed=tag_rows(rows,eta,tag)
        means=[np.concatenate([row["means"] for row in transformed if row["sign"]==sign])
               for sign in (1,-1)]
        weights=[np.concatenate([row["rates"] for row in transformed if row["sign"]==sign])/total
                 for sign in (1,-1)]
        low=min(float(mean.min()) for mean in means)-resolution.Z_CUT*sigma_d
        high=max(float(mean.max()) for mean in means)+resolution.Z_CUT*sigma_d
        boundaries=resolution.classify(transformed,sigma_d)["bayes_boundaries"]
        integral,estimate=quad(lambda x:min(float(resolution.density(x,means[0],weights[0],sigma_d)),
                                          float(resolution.density(x,means[1],weights[1],sigma_d))),
                               low,high,points=boundaries,epsabs=1e-11,epsrel=1e-10,limit=200)
        value+=integral;error+=estimate
    return dict(bayes_error=value,adaptive_quadrature_error_estimate=error,
                note="Adaptive density integration independent of CDF integration; boundary roots are shared as optional break points, not integral values")


def threshold_bracket(sigma_d):
    ideal,untagged=evaluate_tag(0,sigma_d),evaluate_tag(1,sigma_d)
    if ideal["bayes_error"]>TARGET_ERROR:
        return dict(sigma_d=sigma_d,status="No eta passes, even an ideal tag fails",exists=False)
    if untagged["bayes_error"]<=TARGET_ERROR:
        return dict(sigma_d=sigma_d,status="Every eta in [0,1] passes; no tag is required in this model",
                    exists=True,all_eta_pass=True,max_eta=1.0,
                    actual_wrong_tag_probability_at_max=2/3)
    low,high=0.0,1.0
    trace=[]
    while high-low>THRESHOLD_BRACKET_WIDTH:
        mid=(low+high)/2
        result=evaluate_tag(mid,sigma_d)
        trace.append(dict(eta=mid,bayes_error=result["bayes_error"],passes=result["passes_declared_5pct"]))
        if result["passes_declared_5pct"]:low=mid
        else:high=mid
    endpoints=[]
    for eta in (low,high):
        standard=evaluate_tag(eta,sigma_d)
        refined=evaluate_tag(eta,sigma_d,order=72,grid_divisor=8)
        coarse=evaluate_tag(eta,sigma_d,order=24,grid_divisor=8)
        endpoints.append(dict(eta=eta,bayes_error=standard["bayes_error"],
                              refined_error=refined["bayes_error"],coarse_error=coarse["bayes_error"],
                              max_refinement_change=max(abs(standard["bayes_error"]-refined["bayes_error"]),
                                                        abs(standard["bayes_error"]-coarse["bayes_error"])),
                              threshold_margin=abs(standard["bayes_error"]-TARGET_ERROR)))
    return dict(sigma_d=sigma_d,status="Numerically bracketed supremum of passing eta under the exact postulated response",
                exists=True,all_eta_pass=False,eta_bracket=[low,high],bracket_width=high-low,
                actual_wrong_tag_probability_bracket=[2*low/3,2*high/3],
                endpoints=endpoints,bisection_trace=trace,
                interpretation="A numerical threshold bracket, not a statistical confidence interval, physical calibration uncertainty, or apparatus tolerance")


def run():
    source_paths=[Path(__file__),Path(resolution.__file__),Path(resolution.g2.__file__)]
    hashes_before={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    checks=[]

    def check(name,passed,measured=None,tolerance=None):
        row=dict(name=name,passed=bool(passed))
        if measured is not None:row["measured"]=float(measured)
        if tolerance is not None:row["tolerance"]=float(tolerance)
        checks.append(row)

    composition=[]
    for eta in ETAS:
        matrix=response_matrix(eta)
        check(f"response:eta={eta}:nonnegative",np.min(matrix)>=0)
        check(f"response:eta={eta}:column_stochastic",np.max(abs(matrix.sum(axis=0)-1))<1e-14)
        check(f"response:eta={eta}:wrong_tag_probability",np.max(abs(1-np.diag(matrix)-2*eta/3))<1e-14)
    for first in ETAS[:-1]:
        for second in (value for value in ETAS if value>=first):
            lam=(second-first)/(1-first)
            error=float(np.max(abs(response_matrix(lam)@response_matrix(first)-response_matrix(second))))
            composition.append(dict(eta1=first,eta2=second,lambda_additional_noise=lam,max_matrix_error=error))
            check(f"Blackwell_degradation:{first}->{second}",error<1e-14,error,1e-14)
    cards=[]
    for sd in SIGMA_D:
        baseline=resolution.evaluate(SIGMA_I,sd,True)
        scans=[]
        previous=-1.0
        for eta in ETAS:
            result=evaluate_tag(eta,sd)
            refined=evaluate_tag(eta,sd,order=24,grid_divisor=8)
            error=abs(result["bayes_error"]-refined["bayes_error"])
            result["quadrature_grid_refinement_change"]=error
            scans.append(result)
            check(f"si=.1:sd={sd}:eta={eta}:refinement",error<1e-10,error,1e-10)
            check(f"si=.1:sd={sd}:eta={eta}:label_normalization",
                  abs(sum(row["label_probability"] for row in result["label_rows"])-1)<1e-13)
            check(f"si=.1:sd={sd}:eta={eta}:sign_prior_invariance",
                  np.max(abs(np.array(result["sign_priors"])-baseline["priors"]))<1e-13)
            check(f"si=.1:sd={sd}:eta={eta}:rate_invariance",
                  abs(result["total_pair_rate"]/baseline["total_pair_rate"]-1)<1e-13)
            check(f"si=.1:sd={sd}:eta={eta}:data_processing",
                  result["bayes_error"]>=previous-1e-12)
            check(f"si=.1:sd={sd}:eta={eta}:bounded_risk",
                  baseline["ideal_q_tagged_error"]-1e-12<=result["bayes_error"]<=baseline["bayes_error"]+1e-12)
            previous=result["bayes_error"]
        check(f"sd={sd}:eta_zero_ideal_limit",abs(scans[0]["bayes_error"]-baseline["ideal_q_tagged_error"])<1e-12)
        check(f"sd={sd}:eta_one_no_tag_limit",abs(scans[-1]["bayes_error"]-baseline["bayes_error"])<1e-12)
        check(f"sd={sd}:eta_one_label_independent",
              max(abs(row["label_probability"]-1/3) for row in scans[-1]["label_rows"])<1e-13)
        direct=direct_density_integral(.5,sd)
        midpoint=next(row for row in scans if row["eta"]==.5)
        difference=abs(direct["bayes_error"]-midpoint["bayes_error"])
        check(f"sd={sd}:independent_density_integral",difference<1e-9,difference,1e-9)
        threshold=threshold_bracket(sd)
        if threshold.get("all_eta_pass") is False:
            lower,upper=threshold["endpoints"]
            check(f"sd={sd}:threshold_bracket_pass_fail",lower["bayes_error"]<=TARGET_ERROR<upper["bayes_error"])
            check(f"sd={sd}:threshold_margins_exceed_refinement",
                  all(row["threshold_margin"]>100*row["max_refinement_change"] for row in threshold["endpoints"]))
        cards.append(dict(sigma_i=SIGMA_I,sigma_d=sd,
                          role="primary momentum-resolution plan" if sd==.05 else "broader-momentum comparison",
                          no_tag_error=baseline["bayes_error"],ideal_tag_error=baseline["ideal_q_tagged_error"],
                          sign_priors=baseline["priors"],total_pair_rate=baseline["total_pair_rate"],
                          normalized_incident_tail_error_bound=baseline["normalized_incident_tail_error_bound"],
                          scans=scans,independent_integral_at_eta_half=direct,threshold=threshold))
    check("narrower_momentum:all_tested_eta_improve",
          all(second["bayes_error"]<=first["bayes_error"]+1e-12
              for first,second in zip(cards[0]["scans"],cards[1]["scans"])))
    # Raw tag accuracy is not an information guarantee for an imbalanced q prior.
    original_rows,_,_,_=resolution.components(SIGMA_I,True)
    total=sum(float(row["rates"].sum()) for row in original_rows)
    q_priors={str(q):sum(float(row["rates"].sum()) for row in original_rows if row["q"]==q)/total
              for q in TAGS}
    constant_response=np.zeros((3,3));constant_response[1,:]=1
    constant_error=resolution.classify(original_rows,.1)["bayes_error"]
    check("constant_q_zero:response_stochastic",np.max(abs(constant_response.sum(axis=0)-1))<1e-14)
    check("constant_q_zero:no_information",abs(constant_error-cards[0]["no_tag_error"])<1e-14)
    check("constant_q_zero:high_accuracy_yet_fails",
          q_priors["0"]>.889 and constant_error>TARGET_ERROR)
    accuracy_counterexample=dict(response=constant_response.tolist(),q_priors=q_priors,
                                 prediction="Always announce z=0, independently of true q and x",
                                 raw_correct_tag_probability=q_priors["0"],
                                 bayes_error_sigma_d_point_one=constant_error,
                                 information="None: z is constant, so the error equals momentum-only",
                                 lesson="The approximately 88.2 percent threshold applies only to the full symmetric T_eta confusion matrix, not arbitrary aggregate tag accuracy")
    hashes_after={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    check("unchanged_imported_verifiers",hashes_before==hashes_after)
    result=dict(status="G2-M conditional candidate classical noisy-tag comparison; no apparatus or empirical calibration",
        card=dict(resolution.g2.CARD),drive=dict(resolution.DRIVE),cards=cards,
        misleading_raw_accuracy_counterexample=accuracy_counterexample,
        tag_model=dict(labels=list(TAGS),mixing_strengths=list(ETAS),
                       conditional_response="T_eta(z|q)=(1-eta) delta_zq+eta/3",
                       actual_wrong_tag_probability="2 eta/3, not eta",
                       joint_factorization="P(x,z|j,q,p)=N(x;p_f(j,q,p),sigma_d^2) T_eta(z|q)",
                       response_hypothesis="Extra independent tag conditional on true q and p; class-independent response; postulated rather than derived from the classical pump",
                       unnormalized_density="h_jz(x)=sum_q T_eta(z|q) h_jq(x)",
                       risk="E_eta=sum_z integral min(h_+z(x),h_-z(x)) dx / Z",
                       classifier_contract="Bayes-optimal boundaries are recomputed separately for each known calibrated T_eta and sigma_d; this is not the risk of one frozen classifier under unknown tag-response drift",
                       selection="Perfect m=0 and |j|=1 event labeling, equal acceptance, known rate priors, all as in unchanged G2-S"),
        monotonicity=dict(channel_family="T_eta=(1-eta)I+eta U; U_ij=1/3; U^2=U",
                          composition="T_eta2=T_lambda T_eta1 for 0<=eta1<=eta2<=1 and eta1<1; lambda=(eta2-eta1)/(1-eta1)",
                          proof="Extra tag noise is stochastic post-processing leaving x unchanged. min(sum_z A_z,sum_z B_z)>=sum_z min(A_z,B_z) for nonnegative weighted terms; column sums one then give E_eta2>=E_eta1.",
                          endpoints="eta=0: exact q label; eta=1: independent uniform z gives exactly the momentum-only error",
                          numerical_composition_checks=composition),
        limits=[
            "The independent tag is a candidate observation model, not a demonstrated physical sideband detector.",
            "Aggregate tag accuracy is not sufficient: always guessing the dominant q=0 has high accuracy and zero information; the full confusion matrix matters.",
            "A tag inferred from the same measured x cannot be counted as independent evidence; that would double-count information.",
            "Knowing the pump phase, frequency or a clock does not by itself reveal the event's exchanged work q Omega.",
            "The prescribed classical pump has no modeled energy-counting apparatus or quantized pump readout.",
            "Gaussian momentum noise and ternary tag factorize only by the explicit trusted-model assumption; correlated or class-dependent errors require a new calibrated joint response.",
            "These risks assume the true tag response is known and the classifier is optimized for it; frozen-classifier calibration robustness is a separate calculation.",
            "Reported errors are conditional Bayes risks for selected m=0, |j|=1 scattering, not absolute probabilities for the whole beam or all channels.",
            "No SI calibration, quantitative hardware cost, device tolerance, finite 4D collision packet, or absolute luminosity is supplied.",
            "The five-percent boundary is a declared engineering criterion; the numerical eta bracket is not a confidence interval or physical calibration uncertainty.",
            "The G2-S long-time, weak-scattering and common-overlap assumptions remain unchanged."],
        checks=checks,summary=dict(checks=len(checks),passed=sum(row["passed"] for row in checks),
                                  failed=[row["name"] for row in checks if not row["passed"]]),
        source_sha256=hashes_after)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    write_report(result)
    print(json.dumps(result["summary"],indent=2))
    for card in cards:
        print(json.dumps(dict(sigma_d=card["sigma_d"],errors=[row["bayes_error"] for row in card["scans"]],
                              threshold=card["threshold"]),indent=2))
    if result["summary"]["failed"]:raise SystemExit(1)


def write_report(result):
    lines=["# Route G2-M — noisy sideband-tag comparison","",
        "This adds a candidate classical readout model to the unchanged G2-S joint densities. It does not build a detector, modify the scattering action, or claim empirical calibration.","",
        "## A deliberately explicit additional assumption","",
        r"$$z,q\in\{-1,0,1\},\qquad T_\eta(z|q)=(1-\eta)\delta_{zq}+\eta/3,\qquad 0\leq\eta\leq1.$$","",
        "Here eta is a mixing strength: the actual wrong-tag probability is **2 eta/3**, while the correct-tag probability is 1-2 eta/3. eta=0 is an ideal tag; eta=1 produces an independent, uniformly random label.","",
        r"$$P(x,z\mid j,q,p)=\mathcal N(x;p_f(j,q,p),\sigma_d^2)T_\eta(z|q),\quad h_{jz}(x)=\sum_qT_\eta(z|q)h_{jq}(x).$$","",
        r"$$E_\eta=\frac1Z\sum_z\int_{\mathbb R}\min\{h_{+z}(x),h_{-z}(x)\}\,dx,\qquad Z=\sum_{jq}\int h_{jq}(x)dx.$$","",
        "The extra tag is independent of the momentum readout conditional on the physical state, with a response depending only on q. This factorization and its class-independence are **postulated**, not derived from the homogeneous classical pump. A label inferred from the same x cannot supply independent evidence. Knowing the drive phase, frequency or a clock does not tell us the individual event's exchanged work.","",
        "Each reported error is the optimum for a **known, calibrated T_eta**: decision boundaries are recomputed for that response. It is not the error of one frozen classifier when the real tag response drifts or is unknown. Such calibration robustness must be assessed separately; the optimum does not certify it.","",
        "The frozen preparation is sigma_i=0.1, with the same truncated incoming ensemble, finite internal packet, epsilon=0.5, Omega=0.2, perfect m=0/|j|=1 selection, equal acceptance and rate-weighted priors as G2-S. Momentum x is a Gaussian real-valued estimator, not an exact nonnegative norm. All widths remain in the existing engineering units.","",
        "## Conditional error comparison","",
        "| eta | actual wrong-tag probability | error at sigma_d=0.10 | error at sigma_d=0.05 |","|---:|---:|---:|---:|"]
    broad,narrow=result["cards"]
    for first,second in zip(broad["scans"],narrow["scans"]):
        lines.append(f"| {first['eta']:g} | {100*first['actual_wrong_tag_probability']:.2f}% | {100*first['bayes_error']:.4f}% | {100*second['bayes_error']:.4f}% |")
    lines += ["",f"The event priors are fixed at P(j=+1)={broad['sign_priors'][0]:.8g}, P(j=-1)={broad['sign_priors'][1]:.8g}. Tagging changes the posterior, not the production rate or prior. JSON records every tag-label probability and the conditional sign priors; both normalizations are checked.","",
        "## Five-percent criterion: a finite requirement, not a hardware claim",""
    ]
    threshold=broad["threshold"]
    if threshold.get("all_eta_pass") is False:
        low,high=threshold["eta_bracket"]
        wrong_low,wrong_high=threshold["actual_wrong_tag_probability_bracket"]
        lines += [f"For sigma_d=0.10 the largest admissible mixing strength is numerically bracketed by **eta in [{low:.4f}, {high:.4f}]**. This corresponds to an actual wrong-tag probability between **{100*wrong_low:.2f}% and {100*wrong_high:.2f}%** (about {100*(1-(wrong_low+wrong_high)/2):.1f}% correct tagging under this symmetric confusion matrix). Values at the lower bracket endpoint pass; values at the upper endpoint fail the declared five-percent criterion.","",
            "The bracket width is an intentional numerical stopping tolerance, not a statistical confidence interval. Endpoint quadrature/grid refinements are recorded; no experimental response uncertainty or apparatus feasibility has been inferred.",""]
    counterexample=result["misleading_raw_accuracy_counterexample"]
    lines += [f"**Raw tag accuracy is not sufficient.** Here P(q=0)={counterexample['raw_correct_tag_probability']:.8g}. A useless readout that always announces z=0 has **{100*counterexample['raw_correct_tag_probability']:.2f}% correct labels**, yet zero extra information: its recoil-sign error remains **{100*counterexample['bayes_error_sigma_d_point_one']:.4f}%**. The complete confusion matrix, including sideband sensitivity, must be established; an aggregate '88% accurate' specification would be misleading.",""]
    lines += [f"For the primary sigma_d=0.05 plan, momentum alone gives **{100*narrow['no_tag_error']:.4f}%** error. Thus every eta in [0,1] meets five percent: a sideband tag is not required by this trusted-model criterion. An ideal tag would reduce the model error to {100*narrow['ideal_tag_error']:.4f}%, but its physical realization and cost are unknown.","",
        "This is a comparison between tighter momentum resolution and a supplementary observation channel, not a quantitative hardware-cost tradeoff. The main plan can remain momentum-only; the broader-resolution alternative would need a separately established tag response.","",
        "## Why worse tagging cannot help","",
        r"Write $T_\eta=(1-\eta)I+\eta U$, with $U_{zq}=1/3$ and $U^2=U$. For $0\leq\eta_1\leq\eta_2\leq1$ and $\eta_1<1$:","",
        r"$$T_{\eta_2}=T_\lambda T_{\eta_1},\qquad\lambda=\frac{\eta_2-\eta_1}{1-\eta_1}.$$","",
        r"The higher-noise experiment is stochastic post-processing of the lower-noise tag, with x unchanged. For nonnegative terms, $\min(\sum_z a_z,\sum_z b_z)\geq\sum_z\min(a_z,b_z)$; using column sums one proves $E_{\eta_2}\geq E_{\eta_1}$. The eta=1 endpoint is already fully degraded. Consequently the passing set is an interval and the threshold bracket is meaningful.","",
        "The exact eta=0 and eta=1 limits recover ideal q tagging and the original momentum-only result. All listed noise levels satisfy this data-processing order numerically.","",
        "## Verification and limits","",
        f"**{result['summary']['passed']}/{result['summary']['checks']} checks pass.** Tests cover response stochasticity, actual mislabel rates, channel-composition identity, unchanged production priors/rates, tag probabilities, ideal/no-tag endpoints, monotonic Bayes risk, quadrature/grid refinement, and independent adaptive density integration.","",
        "This remains a conditional selected-event forecast. It is not an absolute beam conversion probability, a device design, a quantized-pump apparatus, SI calibration, or a finite 4D collision model. Correlated tag/momentum noise or a tag response depending on recoil sign requires a calibrated joint model before these numbers can be used.","",
        "Reproduce: python3 route_g/code/verify_g2_noisy_tag.py. The G2-S and G2 source files are imported unchanged; their SHA-256 values are recorded and checked before/after execution."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":
    run()
