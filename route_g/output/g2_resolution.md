# Route G2-S: resolution, incoming spread and powered-sideband comparison

Verification: **272/272 checks passed**.

This is a conditional recoil-sign classification forecast for selected scattered events with m=0 and |j|=1, not the probability of a collision. Units are the unchanged G2 card. No measured detector performance is assumed.

The incident ensemble has p>=0 and density exp[-(p-.2)^2/(2 sigma_i^2)]/[sqrt(2 pi) sigma_i Phi(.2/sigma_i)]. sigma_i is the pre-truncation width. The response is a Gaussian real-valued reconstructed momentum estimator with width sigma_d. Both compact labels are otherwise perfectly resolved, efficiency is identical, and all components share a COM frame and dilute external overlap.

Rates are integrated before normalization: A_s(x)=w_l sum_q integral f_i K_lmq R dp; the optimal error is integral min(A_+,A_-) dx / Z. The 5% boundary is an explicitly chosen engineering criterion, not a theorem or confidence level.

## Optimal unobserved-sign errors

Each cell is an error percentage; **bold** means <=5%. Incident spread changes the event-selected priors as well as the line shape.

### Stationary G2 baseline

| sigma_i \ sigma_d | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 |
|---|---:|---:|---:|---:|---:|
| 0.0 | **7.63657e-16%** | **0.027595%** | **3.55155%** | 14.6067% | 21.453% |
| 0.02 | **2.11653e-15%** | **0.0283822%** | **3.55928%** | 14.6095% | 21.4528% |
| 0.1 | **3.16498e-06%** | **0.0634668%** | **3.77384%** | 14.6873% | 21.4477% |
| 0.2 | **0.0721838%** | **0.634898%** | **4.84993%** | 15.0523% | 21.438% |
| 0.4 | **4.96822%** | 7.13216% | 11.0294% | 17.3989% | 21.5736% |

### Powered: epsilon=.5, Omega=.2; q unobserved

| sigma_i \ sigma_d | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 |
|---|---:|---:|---:|---:|---:|
| 0.0 | **1.12187%** | **2.06843%** | 5.76199% | 15.3207% | 21.5409% |
| 0.02 | **1.12271%** | **2.07423%** | 5.76838% | 15.3232% | 21.5407% |
| 0.1 | **1.15327%** | **2.20225%** | 5.94485% | 15.39% | 21.5341% |
| 0.2 | **1.47078%** | **2.74121%** | 6.83823% | 15.7082% | 21.5187% |
| 0.4 | 6.53743% | 8.6925% | 12.271% | 17.8235% | 21.6272% |

## Selected comparison and event priors

| sigma_i | sigma_d | stationary error | powered error | powered ideal q-tag error | stationary majority-only error |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.1 | 3.55155% | 5.76199% | 3.55366% | 22.5425% |
| 0.1 | 0.1 | 3.77384% | 5.94485% | 3.77589% | 22.5195% |
| 0.2 | 0.1 | 4.84993% | 6.83823% | 4.8521% | 22.44% |
| 0.4 | 0.1 | 11.0294% | 12.271% | 11.0331% | 22.2082% |
| 0.1 | 0.2 | 14.6873% | 15.39% | 14.6854% | 22.5195% |

The ideal q-tag is a mathematical information benchmark, not a built pump detector. A real tag must include its noise and work readout. The drive supplies/removes q Omega energy, does not alter compact selection n+l=m+j, and produces additional momentum lines rather than a coordinate change of mass.

## Numerical controls and limitations

- Piecewise Gauss-Legendre orders 24 and 48 and Bayes-boundary grids sigma_d/8 and sigma_d/4 are compared for every scan point. Gaussian component CDF differences are evaluated with survival tails to avoid cancellation; no error floor is imposed.
- Incident z>10 is omitted without reweighting; Kq<=|gq|^2/(32 pi m_n mu_l) bounds the omitted rate. The response-window loss is <=2 Phi(-10). Bounds and convergence differences are stored per point.
- Grid/order agreement is a numerical convergence test, not certified exhaustive isolation of all mixture roots or a rigorous global error bar. The analytic tail bounds do not certify that separate boundary search.
- Tests cover exact sharp-beam G2-R momenta, analytic unequal-prior Gaussian Bayes errors, kinematic derivatives, narrow-spread widths, positivity, normalized priors, noise data processing, and the advantage of ideal sideband information.
- JSON records both false-sign rates, priors, decision boundaries, selected momentum moments, all scans, and representative profile samples. The finite sampled profiles are not used to compute errors.
- Unknown boosts, angular spread, mass-label confusion, unequal acceptance, correlated recoil measurement, finite pump duration and coherent packet interference are outside this card. The powered rate limit presumes resolved sidebands; no finite-time absolute transition probability is claimed.
- The small Gaussian-tail bounds certify numerical integration of the declared tree kernel, not the physical UV validity of that kernel. The EFT cutoff and loop errors are unspecified; an unbounded physical Gaussian beam is not validated by this calculation. A physical beam card must be restricted to a justified EFT domain or supplemented by UV control.

Next: replace this engineering response card with specified or calibrated joint momentum/energy response. Specify four-dimensional incoming packets, timing, geometry and pump envelope only if an absolute collision probability is required.
