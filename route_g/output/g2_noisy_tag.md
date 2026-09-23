# Route G2-M — noisy sideband-tag comparison

This adds a candidate classical readout model to the unchanged G2-S joint densities. It does not build a detector, modify the scattering action, or claim empirical calibration.

## A deliberately explicit additional assumption

$$z,q\in\{-1,0,1\},\qquad T_\eta(z|q)=(1-\eta)\delta_{zq}+\eta/3,\qquad 0\leq\eta\leq1.$$

Here eta is a mixing strength: the actual wrong-tag probability is **2 eta/3**, while the correct-tag probability is 1-2 eta/3. eta=0 is an ideal tag; eta=1 produces an independent, uniformly random label.

$$P(x,z\mid j,q,p)=\mathcal N(x;p_f(j,q,p),\sigma_d^2)T_\eta(z|q),\quad h_{jz}(x)=\sum_qT_\eta(z|q)h_{jq}(x).$$

$$E_\eta=\frac1Z\sum_z\int_{\mathbb R}\min\{h_{+z}(x),h_{-z}(x)\}\,dx,\qquad Z=\sum_{jq}\int h_{jq}(x)dx.$$

The extra tag is independent of the momentum readout conditional on the physical state, with a response depending only on q. This factorization and its class-independence are **postulated**, not derived from the homogeneous classical pump. A label inferred from the same x cannot supply independent evidence. Knowing the drive phase, frequency or a clock does not tell us the individual event's exchanged work.

Each reported error is the optimum for a **known, calibrated T_eta**: decision boundaries are recomputed for that response. It is not the error of one frozen classifier when the real tag response drifts or is unknown. Such calibration robustness must be assessed separately; the optimum does not certify it.

The frozen preparation is sigma_i=0.1, with the same truncated incoming ensemble, finite internal packet, epsilon=0.5, Omega=0.2, perfect m=0/|j|=1 selection, equal acceptance and rate-weighted priors as G2-S. Momentum x is a Gaussian real-valued estimator, not an exact nonnegative norm. All widths remain in the existing engineering units.

## Conditional error comparison

| eta | actual wrong-tag probability | error at sigma_d=0.10 | error at sigma_d=0.05 |
|---:|---:|---:|---:|
| 0 | 0.00% | 3.7759% | 0.0641% |
| 0.1 | 6.67% | 4.6534% | 0.7269% |
| 0.25 | 16.67% | 5.2396% | 1.3208% |
| 0.5 | 33.33% | 5.7090% | 1.9380% |
| 0.75 | 50.00% | 5.8960% | 2.1347% |
| 1 | 66.67% | 5.9449% | 2.2022% |

The event priors are fixed at P(j=+1)=0.77469551, P(j=-1)=0.22530449. Tagging changes the posterior, not the production rate or prior. JSON records every tag-label probability and the conditional sign priors; both normalizations are checked.

## Five-percent criterion: a finite requirement, not a hardware claim

For sigma_d=0.10 the largest admissible mixing strength is numerically bracketed by **eta in [0.1758, 0.1768]**. This corresponds to an actual wrong-tag probability between **11.72% and 11.78%** (about 88.2% correct tagging under this symmetric confusion matrix). Values at the lower bracket endpoint pass; values at the upper endpoint fail the declared five-percent criterion.

The bracket width is an intentional numerical stopping tolerance, not a statistical confidence interval. Endpoint quadrature/grid refinements are recorded; no experimental response uncertainty or apparatus feasibility has been inferred.

**Raw tag accuracy is not sufficient.** Here P(q=0)=0.88983936. A useless readout that always announces z=0 has **88.98% correct labels**, yet zero extra information: its recoil-sign error remains **5.9449%**. The complete confusion matrix, including sideband sensitivity, must be established; an aggregate '88% accurate' specification would be misleading.

For the primary sigma_d=0.05 plan, momentum alone gives **2.2022%** error. Thus every eta in [0,1] meets five percent: a sideband tag is not required by this trusted-model criterion. An ideal tag would reduce the model error to 0.0641%, but its physical realization and cost are unknown.

This is a comparison between tighter momentum resolution and a supplementary observation channel, not a quantitative hardware-cost tradeoff. The main plan can remain momentum-only; the broader-resolution alternative would need a separately established tag response.

## Why worse tagging cannot help

Write $T_\eta=(1-\eta)I+\eta U$, with $U_{zq}=1/3$ and $U^2=U$. For $0\leq\eta_1\leq\eta_2\leq1$ and $\eta_1<1$:

$$T_{\eta_2}=T_\lambda T_{\eta_1},\qquad\lambda=\frac{\eta_2-\eta_1}{1-\eta_1}.$$

The higher-noise experiment is stochastic post-processing of the lower-noise tag, with x unchanged. For nonnegative terms, $\min(\sum_z a_z,\sum_z b_z)\geq\sum_z\min(a_z,b_z)$; using column sums one proves $E_{\eta_2}\geq E_{\eta_1}$. The eta=1 endpoint is already fully degraded. Consequently the passing set is an interval and the threshold bracket is meaningful.

The exact eta=0 and eta=1 limits recover ideal q tagging and the original momentum-only result. All listed noise levels satisfy this data-processing order numerically.

## Verification and limits

**125/125 checks pass.** Tests cover response stochasticity, actual mislabel rates, channel-composition identity, unchanged production priors/rates, tag probabilities, ideal/no-tag endpoints, monotonic Bayes risk, quadrature/grid refinement, and independent adaptive density integration.

This remains a conditional selected-event forecast. It is not an absolute beam conversion probability, a device design, a quantized-pump apparatus, SI calibration, or a finite 4D collision model. Correlated tag/momentum noise or a tag response depending on recoil sign requires a calibrated joint model before these numbers can be used.

Reproduce: python3 route_g/code/verify_g2_noisy_tag.py. The G2-S and G2 source files are imported unchanged; their SHA-256 values are recorded and checked before/after execution.
