# P54: moving-threshold sequential seesaw

Date: 2026-09-06.

Status: 46/46 checks pass. An operational SM+N matrix flow now solves the running mass thresholds, performs exact block-tree matching, and continues below the last sterile neutrino. The actual nonzero scalar-tree CHN now enters the complete one-loop dimension-five running subsystem with CBN=0 and a running Higgs mass parameter, plus tree matching, not a completed finite one-loop EFT calculation or a physical flavor fit.

## 1. Input contract and conventions

Both runs read the existing local-light JSON. Yu,Yd,Ye,Ynu are the four matrices evaluated on its actual self-consistently corrected scalar light vector. MR/omega is its common-phase Majorana matrix. The two synthetic h_raw,f_raw matrices and actual Clifford projections determine all these inputs; the code does not invent four free UV family matrices. CHN is read from p54_scalar_chn.json, evaluated on the same corrected external light ray with the frozen tree scalar Hessian and vertices. It is correlated with the same f_raw, not a new family matrix. This external-ray improvement is explicitly mixed order, not a completed loop-level scalar matching calculation.

All Dirac matrices have doublet-family rows and singlet-family columns. With $H^0=v/\sqrt2$, define

$$V(H)=m_Z^2H^\dagger H+\lambda(H^\dagger H)^2,\qquad q_H=m_Z^2/\omega^2,\qquad M_\nu=\frac{v^2}{2}\left(C_5-Y_\nu M_R^{-1}Y_\nu^T\right).$$

The program stores $\widehat M_R=M_R/\omega$, $\widehat C_5=\omega C_5$, $\widehat C_{HN}=\omega C_{HN}$, and $q_H$, so $M_\nu=v^2[\widehat C_5-Y_\nu\widehat M_R^{-1}Y_\nu^T]/(2\omega)$. The initial type-II coefficient is $\widehat C_5^{II}=2C_{II,\rm raw}f_{\rm raw}$. This factor of two is tested against the pre-existing local I+II mass matrix. The no-half operator convention is $\mathcal L\supset-[N^TM_RN/2+\mathrm{h.c.}]+[C_{HN}NN(H^\dagger H)+\mathrm{h.c.}]$, hence $M_R(\rho)=M_R-2C_{HN}\rho$. Here $m_Z^2$ is Zhang's signed Lagrangian parameter, not the Z-boson mass. If another card writes $V=-m^2H^\dagger H+\lambda(H^\dagger H)^2$, then $q_H=-m^2/\omega^2$.

The illustrative boundary g=(0.53,0.55,0.60), lambda=0.12 and mu_start/omega=0.1 are declared diagnostic inputs, not solved P54 matching data. The additional qH=0 boundary is a diagnostic; it is not an electroweak-matched Higgs mass and the running generally regenerates nonzero qH. In particular, no scalar-card lambda_eff with another normalization is silently reused. Dimensionless thresholds are the primary results. The reference omega=10^14 GeV and v=246.22 GeV only convert them into explicitly synthetic terminal mass proxies.

## 2. One-loop equations with every surviving sterile state

Write $H_x=Y_xY_x^\dagger$, $S_N=Y_\nu^\dagger Y_\nu$, $T=\operatorname{tr}(3H_u+3H_d+H_e+H_\nu)$, and $T_4=\operatorname{tr}(3H_u^2+3H_d^2+H_e^2+H_\nu^2)$. Let $\mathcal B=16\pi^2d/d\log\mu$. In GUT hypercharge normalization:

$$\begin{aligned}\mathcal B_{Y_u}&=\left[\frac32(H_u-H_d)+T-\frac{17}{20}g_1^2-\frac94g_2^2-8g_3^2\right]Y_u,\\\mathcal B_{Y_d}&=\left[\frac32(H_d-H_u)+T-\frac14g_1^2-\frac94g_2^2-8g_3^2\right]Y_d,\\\mathcal B_{Y_e}&=\left[\frac32(H_e-H_\nu)+T-\frac94g_1^2-\frac94g_2^2\right]Y_e,\\\mathcal B_{Y_\nu}&=\left[\frac32(H_\nu-H_e)+T-\frac9{20}g_1^2-\frac94g_2^2\right]Y_\nu+3C_5Y_\nu^*M_R-4Y_\nu M_R^\dagger C_{HN},\\\mathcal B_{M_R}&=S_N^TM_R+M_RS_N-8m_Z^2C_{HN}.\end{aligned}$$

With $P_C=-3H_e/2+7H_\nu/2$ and $\alpha_C=4\lambda-3g_2^2+2T$:

$$\mathcal B_{C_5}=P_CC_5+C_5P_C^T+\alpha_CC_5,\qquad\mathcal B_{C_5^{II}}=P_CC_5^{II}+C_5^{II}P_C^T+\alpha_CC_5^{II}.$$

$$\mathcal B_{g_i}=b_i g_i^3,\quad b=(41/10,-19/6,-7),$$

$$\mathcal B_\lambda=24\lambda^2-3\lambda(3g_2^2+\tfrac35g_1^2)+\tfrac34g_2^4+\tfrac38(g_2^2+\tfrac35g_1^2)^2+4\lambda T-2T_4-4\operatorname{Re}\operatorname{tr}(C_5Y_\nu^*M_RY_\nu^\dagger)+16\operatorname{Re}\operatorname{tr}(C_{HN}^\dagger M_RS_N).$$

Define $A_H=12\lambda+2T-\frac9{10}g_1^2-\frac92g_2^2$. The remaining nonzero coefficients obey

$$\mathcal B_{C_{HN}}=A_HC_{HN}+4(C_{HN}S_N+S_N^TC_{HN}),$$

$$\mathcal B_{m_Z^2}=A_Hm_Z^2-4\operatorname{tr}(M_R^\dagger M_RS_N)+8\operatorname{Re}\operatorname{tr}(M_R^\dagger M_RM_R^\dagger C_{HN}).$$

All these formulas are dimensionful; the code divides by powers of the fixed reference $\omega$ to evolve the hatted variables. It does not identify $\omega$ with the running scale. The beta function of $C_{BN}$ is homogeneous, so its zero boundary remains zero at this order.

The CHN feedback signs have an independent CW-divergence test: for one real family take $\mathscr M(h)=\left(\begin{smallmatrix}0&yh/\sqrt2\\yh/\sqrt2&M-Ch^2\end{smallmatrix}\right)$. A polynomial regression of $\operatorname{tr}[(\mathscr M^\dagger\mathscr M)^2]$ yields $[h^2]=2M^2y^2-4M^3C$ and $[h^4]=y^4/2-4My^2C+O(C^2)$. Calibrating against the renormalizable mass/quartic terms gives the displayed +8 M^3 C and +16 M y^2 C corrections. Taking the odd part in C removes the double-insertion contribution.

The renormalizable baseline uses the transposed-row convention of [Antusch et al., Appendix D.1](https://arxiv.org/pdf/hep-ph/0501272). Their quartic is lambda_A/4, so lambda_A=4 lambda here. The +4 lambda coefficient is independently confirmed in the canonical convention by [Wang, Zhang and Zhou, Eq. (2.15)](https://arxiv.org/pdf/2302.08140). No sterile gauge contribution is added because the retained neutrinos are SM singlets.

Crucially, the default equations include the three corrections in [Di Zhang (2024), Eqs. (5)–(6), (8)](https://arxiv.org/html/2405.18017v2): Weinberg feedback into Ynu and lambda and the replacement 1/2 by 7/2 in its own anomalous dimension. Our C5 is minus that paper's coefficient, explaining the feedback signs. An explicitly named negative control removes only these three Weinberg corrections while retaining the CHN and qH equations. The CHN signs do not flip with the Weinberg convention. The actual tree scalar exchange falsifies the old CHN=0 boundary assumption. CBN has no source from this nonderivative scalar-tree graph, but its finite loop matching remains uncomputed.

The state space changes from 3 to 2 to 1 to 0 sterile columns. The full complex symmetric MR matrix runs within every interval. C5II runs continuously, while the total C5 receives type-I threshold increments. C5II is a tagged initial-source component along the common nonlinear trajectory, not a separately fitted or scheme-independent decomposition of loop effects.

## 3. Why one frozen Weinberg matrix is insufficient

First restrict to $C_{HN}=0$ to isolate the Weinberg issue, and set $Q=Y_\nu M_R^{-1}Y_\nu^T$. The product rule, including $dM_R^{-1}=-M_R^{-1}(dM_R)M_R^{-1}$, gives

With $P=-3H_e/2+H_\nu/2$,

$$\mathcal B_Q=PQ+QP^T+\alpha_QQ+3(H_\nu C_5+C_5H_\nu^T),\qquad\alpha_Q=2T-\frac9{10}g_1^2-\frac92g_2^2.$$

Thus for $C_{\rm eff}=C_5-Q$:

$$\mathcal B_{C_{\rm eff}}=PC_{\rm eff}+C_{\rm eff}P^T+\alpha_CC_{\rm eff}+\underbrace{\left(4\lambda+\frac9{10}g_1^2+\frac32g_2^2\right)Q}_{\text{active-seesaw source}}.$$

This is a structural, not numerical, distinction: above thresholds the type-I and type-II pieces have different flavor-universal running. Evolving only their initial sum with the below-threshold Weinberg equation loses the displayed source. The new 2024 contributions to beta C5 and beta Q cancel instantaneously in their difference, so this last equation retains its form. Nevertheless Ynu and lambda follow changed trajectories; the JSON quantifies the resulting difference from the system with only these corrections omitted. Both identities are checked directly in the canonical CHN=0 subsector; neither is claimed without its extra terms at nonzero CHN.

For actual nonzero $C=C_{HN}$, the additional contribution to $\mathcal B_Q$ is

$$\Delta_{HN}\mathcal B_Q=-4Y_\nu M_R^\dagger C M_R^{-1}Y_\nu^T-4Y_\nu M_R^{-1}C M_R^*Y_\nu^T+8m_Z^2Y_\nu M_R^{-1}C M_R^{-1}Y_\nu^T.$$
Thus beta(C5-Q) acquires minus this expression. Keeping only a single total Weinberg coefficient would also lose this independent scalar-exchange effect. The implementation evolves every matrix and obtains the composite by matrix inversion, avoiding such a closure assumption.

## 4. Derivation of the sign and block matching

At a chosen sterile split, rotate by $N_{\rm old}=U N_{\rm new}$, so $Y'=YU$ and $M'=U^TMU$. Write

$$M'=\begin{pmatrix}A&B\\B^T&D\end{pmatrix},\qquad Y'=(Y_r,Y_h).$$

The algebraic heavy equation of motion is $N_h=-D^{-1}(B^TN_r+Y_h^T\ell H)$. Substitution, or the block-inverse identity, gives

$$\boxed{\ C_5^-=C_5^+-Y_hD^{-1}Y_h^T,\quad Y_r^-=Y_r-Y_hD^{-1}B^T,\quad M_r^-=A-BD^{-1}B^T\ }.$$

Consequently $C_5^--Y_r^-(M_r^-)^{-1}(Y_r^-)^T=C_5^+-YM^{-1}Y^T$ exactly at tree level. The minus sign follows from the declared mass convention, not a convention-free choice of Weinberg operator sign. Iterated Schur complements equal the one-shot Schur complement; a non-diagonal numerical test also exercises the retained-Yukawa correction.

For $M(\rho)=M-2C_{HN}\rho$, differentiate the same Schur complement. Writing $T_N=\binom{I}{-D^{-1}B^T}$ in retained/heavy order gives

$$C_{HN}^-=T_N^T C_{HN}'T_N=C_{rr}-C_{rh}D^{-1}B^T-BD^{-1}C_{hr}+BD^{-1}C_{hh}D^{-1}B^T.$$
This formula is checked against an independent finite difference of the field-dependent Schur complement and against two successive eliminations. Merely truncating the retained corner is wrong for a general off-diagonal mass split. A CHN insertion in the induced Weinberg term instead belongs to dimension seven and is not included in this dimension-five truncation.

At a Takagi threshold the cross block B vanishes up to roundoff. The general block formula is nevertheless used, so roundoff cannot silently spoil the tree identity. This is exact matrix algebra at leading seesaw order; it does not assert exact light eigenvalues of the full electroweak mass matrix or include dimension-six kinetic effects.

## 5. Moving events and degeneracies

The solver detects $f(t)=t-\log\sigma_{\max}[\widehat M_R(t)]=0$, where $t=\log(\mu/\omega)$. The singular values are recomputed inside the event function. Neither the initial eigenvalues nor an initial ordering determine the threshold sequence.

For symmetric $M=X+iY$, the real symmetric matrix $\begin{pmatrix}X&-Y\\-Y&-X\end{pmatrix}$ has positive eigenvectors $(a,b)$ yielding $u=a+ib$ with $Mu=m u^*$. This gives $U^TMU=\operatorname{diag}(m)$ without square-root phase choices. A complete near-degenerate positive-mass cluster is removed together, using a declared relative log-mass tolerance of 10^-7. Exact degeneracy is insensitive to the remaining real-orthogonal Takagi freedom. A separate test uses arbitrary complex sterile basis transformations. Rank-deficient Takagi blocks are rejected for sterile threshold events rather than assigning a threshold to a genuinely massless state. The terminal light-neutrino Takagi interface separately allows a kernel. A minimal two-sterile test preserves its massless light mode through the whole corrected flow.

The level-crossing negative control begins with M1=0.049 < M2=0.050 but a much larger second-column Yukawa. The second mass runs downward faster, and the state initially called N1 actually leaves first. The removed projector, rather than an ambiguous eigenvector phase, is saved and tested.

## 6. Numerical checkpoints

| Local-light case | Initial descending masses / omega | Moving event scales / omega |
|---|---|---|
| 0 | 0.04651348, 0.02862368, 0.01788980 | 0.04648960, 0.02861195, 0.01786944 |
| 1 | 0.06082533, 0.03935756, 0.02504572 | 0.06080995, 0.03932527, 0.02503881 |

Every event is followed by continued running, including a final zero-sterile interval. The JSON records all initial/final matrices, moving mass spectra, matching continuity, terminal PMNS/CKM moduli, Jarlskog invariants and fixed-v running mass proxies. They are synthetic diagnostics, not agreement with experiment.

| Check | Residual | Pass |
|---|---:|:---:|
| case 0: actual local I+II normalization reconstruction | 2.118e-16 | yes |
| case 0: actual scalar-tree CHN is nonzero and runs | 0.000e+00 | yes |
| case 0: Higgs qH generated from diagnostic zero boundary | 0.000e+00 | yes |
| case 0: full running plus moving thresholds family covariance | 1.965e-15 | yes |
| case 0: event scales are sterile-basis independent | 1.812e-16 | yes |
| case 0: all three moving thresholds and below-last running | 0.000e+00 | yes |
| case 0: live threshold roots | 0.000e+00 | yes |
| case 0: exact tree continuity at every threshold | 2.718e-16 | yes |
| case 0: thresholds are not frozen initial masses | 0.000e+00 | yes |
| case 0: terminal PMNS and neutrino mass covariance | 1.583e-15 | yes |
| case 0: active and below-last stages evolve nontrivially | 0.000e+00 | yes |
| case 0: actual nonzero CHN changes the full trajectory | 0.000e+00 | yes |
| case 1: actual local I+II normalization reconstruction | 2.283e-16 | yes |
| case 1: actual scalar-tree CHN is nonzero and runs | 0.000e+00 | yes |
| case 1: Higgs qH generated from diagnostic zero boundary | 0.000e+00 | yes |
| case 1: full running plus moving thresholds family covariance | 2.944e-15 | yes |
| case 1: event scales are sterile-basis independent | 4.272e-16 | yes |
| case 1: all three moving thresholds and below-last running | 0.000e+00 | yes |
| case 1: live threshold roots | 8.882e-16 | yes |
| case 1: exact tree continuity at every threshold | 2.911e-16 | yes |
| case 1: thresholds are not frozen initial masses | 0.000e+00 | yes |
| case 1: terminal PMNS and neutrino mass covariance | 2.007e-15 | yes |
| case 1: active and below-last stages evolve nontrivially | 0.000e+00 | yes |
| case 1: actual nonzero CHN changes the full trajectory | 0.000e+00 | yes |
| non-diagonal exact block-tree sequential equals one-shot Schur | 2.593e-16 | yes |
| non-diagonal Schur test exercises retained-Yukawa shift | 0.000e+00 | yes |
| CHN block pullback equals field-dependent Schur derivative | 7.383e-13 | yes |
| non-diagonal CHN matching is not naive corner truncation | 0.000e+00 | yes |
| CHN sequential pullback equals one-shot retained Schur derivative | 1.325e-16 | yes |
| exact degenerate Takagi cluster removed as one block | 0.000e+00 | yes |
| degenerate threshold and final C5 are fully family covariant | 6.992e-16 | yes |
| degenerate Takagi O(3) basis choice leaves block matching invariant | 3.552e-16 | yes |
| live event ordering detects a running sterile-mass level crossing | 0.000e+00 | yes |
| one-loop full matrix beta family covariance | 2.507e-14 | yes |
| CHN feedback uses Zhang positive-potential-mass signs | 1.380e-14 | yes |
| MR and CHN one-loop derivatives preserve complex symmetry | 8.895e-17 | yes |
| zero CHN remains RG-invariant but is not the actual scalar boundary | 0.000e+00 | yes |
| independent CW polynomial fixes CHN Higgs mass/quartic feedback signs | 1.346e-15 | yes |
| active seesaw composite beta includes corrected dimension-five feedback | 1.781e-16 | yes |
| I-plus-II coefficient has the required active-seesaw source term | 2.835e-16 | yes |
| one Weinberg coefficient alone is not closed above all thresholds | 0.000e+00 | yes |
| three Zhang 2024 corrections have the declared C5 sign | 2.155e-13 | yes |
| new dimension-five terms cancel instantaneously in beta(C5-Q) | 1.888e-16 | yes |
| canonical Higgs quartic gives +4lambda in below-threshold C5 beta | 1.428e-16 | yes |
| minimal two-sterile model preserves a massless neutrino at one loop | 3.041e-17 | yes |
| terminal Takagi interface admits the physical rank-two light spectrum | 0.000e+00 | yes |

## 7. Finite matching is still a separate gate

[Zhang and Zhou, Eqs. (87)–(88)](https://arxiv.org/pdf/2107.12133) provide finite one-loop matching when the full type-I sector is integrated out. This includes hard vertex terms and field-normalization contributions; changing C5 alone is not the entire matching operation. That all-sterile formula is not silently applied to a partially active sterile EFT with a pre-existing type-II coefficient. The current solver deliberately uses tree matching and one-loop running. Finite logarithms/constant terms, all lower scalar thresholds, and the relevant higher-dimensional operator closure remain open.

The running subsystem now includes all one-loop dimension-five feedback in the dipole-free nuSMEFT, with actual nonzero scalar-tree CHN. The retained PQ axion is a spectator here; its interactions and loop thresholds are not incorporated, so this is not the complete P54 EFT. Finite scalar/Yukawa/dipole matching and a matched Higgs mass boundary remain open. Dimension-six operators, double dimension-five insertions and exact finite-v pole matching are outside the implementation. The modest mass ratios in the benchmarks do not establish that omitted terms are negligible. A finite and operator-complete matching audit is required before promoting these diagnostic terminal quantities to physical predictions.

The terminal applicability flags are explicitly false: the generated positive qH is much larger than the terminal scale squared. A mass-independent MS-bar ODE can be continued as a mathematical diagnostic, but retaining a light Higgs down to that endpoint requires quadratic finite matching and retuning. Fixed-v mass proxies are not a physically valid endpoint in these runs.

Next: connect the properly matched lower scalar EFT and independently evolved six-invariant PS inputs, including the calculable conjugate-bidoublet matrices, then add finite sequential type-I/type-II matching with a consistent operator basis. A full flavor fit, Higgs pole mass, determinant or portal is not promoted by this checkpoint.
