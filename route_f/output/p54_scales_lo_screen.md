# P54 bounded LO scale / gauge-proton screen

**Conditional scale solutions, not accepted physical benchmarks.** The historical frozen scalar point and its old P2 scales are not inputs.

Checks: **52/52**. Correct four-parent PS beta: `(2/3,26/3,26/3)`, not historical `(1,26/3,26/3)`.

## Declared approximation

One-loop gauge running, tree gauge matching, one light SM Higgs doublet, three 16F families and the four saved active PS parents. The gauge-singlet axion does not affect these gauge betas. The central input card is a deliberately approximate inherited screening card, not a current correlated electroweak fit. Finite gauge constants, two-loop running and scalar-action realization are not included.

For `t=log(MI/MZ)`, `u=log(MU/MI)`, solve the three linear equations

`alpha_SM^-1(MZ) = alphaU^-1 1 + b_SM t/(2 pi) + P b_PS u/(2 pi) - P sum_A delta_b_A log(kappa_A)/(2 pi)`,

with `P(4,L,R)=(4,L,2*4/5+3*R/5)` and `delta_b_A=T_real,A/6`. Each logarithm belongs to an entire declared PS parent; no independent per-gauge compensator is introduced. The nine cards vary whole `Phi54:(20prime,1,1)` and `Phi54:(1,3,3)` reciprocally by `(1/3,1,3)`, and whole `Sigma126:(15,2,2)` by `(1,2,3)` at MI. Other stage masses stay at their reference scale; no Goldstone parent is moved. These are provisional mass hypotheses, not proofs that all nine spectra come from a stable scalar action.

All cards preserve left-right parity. Therefore `5 alpha1^-1-3 alpha2^-1-2 alpha3^-1` fixes MI independently of these logs. A narrow/vertical MI locus is a structural LO result, not a lack of scan coverage.

## Scale results

Baseline: `MI=5.08227e+13 GeV`, `MX=MU=1.89935e+15 GeV`, `alphaU^-1=38.222412`; `g4(MI)=0.570524`, `sigma=8.90807e+13 GeV`. The actual 45-generator orbit verifies `MX=gU*omega`, `MXprime=gU*sqrt(omega^2+sigma^2)` and charged PS-vector mass `g4(MI)*sigma`. The latter uses the intermediate coupling, not gU.

| Card | log10 MI | log10 MX | alphaU inverse | effective flavor norm limit |
|---|---:|---:|---:|---:|
| P20_0.333333_P33_3_B15_1 | 13.70606 | 15.41777 | 37.9553 | 3.6252 |
| P20_0.333333_P33_3_B15_2 | 13.70606 | 15.40522 | 38.5467 | 3.4751 |
| P20_0.333333_P33_3_B15_3 | 13.70606 | 15.39789 | 38.8927 | 3.3898 |
| P20_1_P33_1_B15_1 | 13.70606 | 15.27861 | 38.2224 | 1.9234 |
| P20_1_P33_1_B15_2 | 13.70606 | 15.26606 | 38.8138 | 1.8435 |
| P20_1_P33_1_B15_3 | 13.70606 | 15.25873 | 39.1598 | 1.7981 |
| P20_3_P33_0.333333_B15_1 | 13.70606 | 15.13945 | 38.4895 | 1.0204 |
| P20_3_P33_0.333333_B15_2 | 13.70606 | 15.12690 | 39.0810 | 0.9779 |
| P20_3_P33_0.333333_B15_3 | 13.70606 | 15.11957 | 39.4269 | 0.9538 |

## What the proton number means

Use `C_L,R(2 GeV)=gU^2/(2 MX^2) A_L,R F_L,R` for scalar three-quark operators. F includes the Fierz/Clebsch factors, fitted flavor rotations and the coherent second-vector term weighted by `MX^2/MXprime^2`; F=1 is merely a coefficient normalization, not a model prediction. The one-loop gauge RG factors are integrated with the same PS parent thresholds and SM running, followed by QCD to 2 GeV. [Gauge normalization and exponents](https://arxiv.org/html/1507.06712v2).

The direct continuum lattice value is `W_pi+=-0.159(15)(20)(25) GeV^2`; isospin gives `|W_pi0|=|W_pi+|/sqrt(2)` in MSbar at 2 GeV. No extra chiral-Lagrangian `(1+D+F)/f_pi` multiplier is applied. [Lattice Table 8](https://arxiv.org/html/2111.01608v1).

`Gamma_e_pi = mp/(32 pi) (1-mpi^2/mp^2)^2 |W_pi0|^2 [|C_L|^2+|C_R|^2]`.

Using the published `tau/Br > 2.4e+34 yr` limit gives `sqrt(|C_L|^2+|C_R|^2)<8.76417e-32 GeV^-2` at the central hadronic input. [Super-K original measurement](https://arxiv.org/abs/2010.16098).

Equivalently every card must obey `sqrt(A_L^2 |F_L|^2 + A_R^2 |F_R|^2) < saved limit`. The table is this necessary condition in the stated gauge-only approximation, not a green/excluded model map. A pure-channel diagnostic and one-standard-error hadronic sensitivity are exported in JSON, but no joint confidence region is implied.

## Explicit decisions

- A nonordered scale solution or gauge Landau pole fails the declared scale card.
- A specified flavor realization that violates the saved inequality fails this gauge-only channel test; F cannot be reset to zero to declare success.
- Whole-model exclusion needs the correlated allowed flavor domain and all relevant amplitudes/channels. Scalar exchange and interference are currently uncomputed, so no unconditional exclusion is made.
- All nine cards remain `physical_candidate=false`: vacuum realization, scalar-mediated proton decay, other channels and the common constrained flavor fit have not passed.

Next use the saved sigma and correlated vector masses in the LO flavor screen, or test whether the action can realize a promising parent-mass pattern. Do not restore the obsolete P2 beta table or old MU threshold result.
