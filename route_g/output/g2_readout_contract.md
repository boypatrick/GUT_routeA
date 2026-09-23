# G2 readout contract: frozen momentum-only rule

Verification: **154/154 checks passed**. Engineering forecast only; no calibrated hardware or absolute collision probability.

The main card keeps the existing driven action (epsilon=.5, Omega=.2), sigma_i=.1, perfect m=0 and |j|=1 selection, and a Gaussian real-valued momentum estimator with nominal sigma_d=.05. Sideband q is unobserved. All quantities use unchanged arbitrary model units.

## Deployable frozen decision

Predict j=+1 when x < **1.258929**; predict j=-1 otherwise. One boundary is found in both ordinary and refined searches. The rule is trained once at a=b=0 and sigma=.05.

Nominal overall error: **2.2022%**. Majority-only error: 22.5304%.

| True sign | Predict + | Predict - | Recall | Purity of this predicted label |
|---|---:|---:|---:|---:|
| +1 | 98.9685% | 1.0315% | 98.9685% | 98.2028% |
| -1 | 6.2279% | 93.7721% | 93.7721% | 96.3556% |

The two middle columns are conditional on the true sign. Purity instead conditions on the predicted sign; it uses the event-selected prior and is not the same as recall. The chosen 5% ceiling applies to the overall event-weighted error, not each class: the nominal j=-1 error is 6.23%, so this contract does not claim 95% recall for both signs.

## Calibration stress card

True response: x=(1+a)pf+b+Normal(0,sigma_true^2), with a in {-0.01,0,0.01}, b in {-0.02,0,0.02}, sigma_true in {0.04,0.05,0.06}. The frozen classifier is never adjusted to these true values. The oracle is retrained and reported only as a comparison.

| a | b | sigma_true | Frozen error | Retuned oracle error |
|---:|---:|---:|---:|---:|
| -0.01 | -0.02 | 0.04 | 2.0972% | 1.7041% |
| -0.01 | -0.02 | 0.05 | 2.9692% | 2.2296% |
| -0.01 | -0.02 | 0.06 | 3.8648% | 2.7783% |
| -0.01 | +0.00 | 0.04 | 1.7108% | 1.7041% |
| -0.01 | +0.00 | 0.05 | 2.3324% | 2.2296% |
| -0.01 | +0.00 | 0.06 | 3.0718% | 2.7783% |
| -0.01 | +0.02 | 0.04 | 1.8539% | 1.7041% |
| -0.01 | +0.02 | 0.05 | 2.2542% | 2.2296% |
| -0.01 | +0.02 | 0.06 | 2.7848% | 2.7783% |
| +0.00 | -0.02 | 0.04 | 1.7517% | 1.6849% |
| +0.00 | -0.02 | 0.05 | 2.4492% | 2.2022% |
| +0.00 | -0.02 | 0.06 | 3.2454% | 2.7451% |
| +0.00 | +0.00 | 0.04 | 1.7376% | 1.6849% |
| +0.00 | +0.00 | 0.05 | 2.2022% | 2.2022% |
| +0.00 | +0.00 | 0.06 | 2.7997% | 2.7451% |
| +0.00 | +0.02 | 0.04 | 2.1110% | 1.6849% |
| +0.00 | +0.02 | 0.05 | 2.3891% | 2.2022% |
| +0.00 | +0.02 | 0.06 | 2.7831% | 2.7451% |
| +0.01 | -0.02 | 0.04 | 1.6709% | 1.6663% |
| +0.01 | -0.02 | 0.05 | 2.2021% | 2.1755% |
| +0.01 | -0.02 | 0.06 | 2.8650% | 2.7126% |
| +0.01 | +0.00 | 0.04 | 1.9235% | 1.6663% |
| +0.01 | +0.00 | 0.05 | 2.2592% | 2.1755% |
| +0.01 | +0.00 | 0.06 | 2.7153% | 2.7126% |
| +0.01 | +0.02 | 0.04 | 2.4612% | 1.6663% |
| +0.01 | +0.02 | 0.05 | 2.6338% | 2.1755% |
| +0.01 | +0.02 | 0.06 | 2.9123% | 2.7126% |

Worst sampled frozen error: **3.8648%**, at a=-0.01, b=-0.02, sigma_true=0.06. Its oracle error is 2.7783%.

## Continuous-box bound and numerical qualification

For the fixed threshold t, plus-to-minus error is Phi((mu-t)/sigma), maximized at mu=(1+a_max)pf+b_max. Minus-to-plus error is Phi((t-mu)/sigma), maximized at mu=(1+a_min)pf+b_min. For either signed difference d, max_{sigma in [s_min,s_max]} Phi(d/sigma) occurs at an endpoint. Sum these componentwise maxima using the unchanged normalized rate weights.

The resulting conservative continuous-box envelope is **5.3902%**. Unlike a grid maximum, the inequality covers all fixed calibration points in the stated continuous box. Different components may maximize at different calibrations, so this envelope need not be attainable.

A single predeclared equal 2x2x2 partition of the shared a,b,sigma box reduces the covering envelope to **4.4571%**. This maximum of eight subbox bounds differs both from the loose global envelope and from the 27-point sampled maximum; it does not retune the classifier.

Each bound integral splits the kinks where component mean crosses threshold. Refined adaptive integration and beam-normalization difference 6.9e-18; formal tree-model omitted-beam-tail error bound 4.69e-23. This is a converged numerical value of an analytic envelope, not interval-certified arithmetic or a rigorous global numerical proof. No EFT cutoff or UV completion is inferred from the tiny Gaussian tail.

## Available nominal resolution

At exact a=b=0, the 5% optimal-Bayes ceiling crosses near sigma_d=0.0907; a rounded-down requirement is **sigma_d <= 0.090** if the classifier is retrained for that resolution. Keeping the actual frozen rule crosses near 0.0891, with rounded-down requirement **sigma_d <= 0.089**. Neither number applies unchanged to unknown calibration errors or altered beam/drive parameters.

Next: specify or calibrate a joint momentum response and mass-label confusion/acceptance. Treat finite noisy sideband tagging as a separate comparison; it is not used to rescue this momentum-only result. Four-dimensional packet geometry, pump timing and an EFT-valid beam domain are required before an absolute collision probability is meaningful.
