# P54 generates a sterile-Higgs dimension-five operator at tree level

Date: 2026-09-06.

18/18 checks pass. The actual frozen P54 action does not match onto CHN=0. The primary result uses its exact tree light eigenvector, with no new parameters:

$$\boxed{\ \omega C_{HN}=-0.0848692011507\,i\, f_{\rm raw}=-0.0150028969119\,f_M\ }.$$

The displayed phase belongs to the fixed common fermion convention. This is a new required Wilson coefficient, not a new arbitrary family matrix: it is proportional to the same Majorana Yukawa matrix at this tree matching point.

## Derivation

Let $s_A$ denote canonical massive real SM-singlet fluctuations, $\rho=H^\dagger H$, and write

$$\mathcal L\supset-\frac12s^TKs-s^TJ\rho-\left[\frac12N^T(M_R+G_As_A)N+\mathrm{h.c.}\right].$$

Here $K$ is the same-action scalar Hessian, $J_A=V'''[s_A,q,q]$ for a unit canonical real Higgs component $q$, and $G_A=\partial M_R/\partial s_A$. SM invariance makes the result the coefficient of $\rho=\frac12\sum_{\alpha=1}^4h_\alpha^2$, not a neutral-component-only interaction. At leading derivative order the scalar equation gives $s=-K^{-1}J\rho+\cdots$. Therefore

$$\delta M_R=-G^TK^{-1}J\rho,\qquad\mathcal L_{\rm eff}\supset C_{HN}^{ij}N_iN_jH^\dagger H+\mathrm{h.c.},\qquad\boxed{\ C_{HN}=\frac12G^TK^{-1}J\ }.$$

The factor 1/2 is required because the Majorana mass term has it while the CHN operator does not. This agrees with the operator normalization in [Di Zhang, Eqs. (2)–(3)](https://arxiv.org/html/2405.18017v2).

In program units $\widehat K=K/\omega^2$, $\widehat J=J/\omega$, and $G_A=g_Af_{\rm raw}$, so $\omega C_{HN}=\frac12g^T\widehat K^{-1}\widehat J f_{\rm raw}$. The actual common-phase nuc spinor projects the real scalar Yukawa tensors to g. Its contraction with the vacuum reproduces the exported MR, including the absolute Clifford factor and phase.

## Which modes were integrated out

The exact SM Casimir kernel has five canonical real singlets. Its Hessian has three positive radial modes and two zero modes associated with gauge/PQ directions. Only the positive block is inverted. The gauge/PQ modes are retained (the axion held as a spectator); no infrared zero eigenvalue is divided by.

| Radial mode | mass squared / omega squared | J / omega | g_raw |
|---|---:|---:|---:|
| 0 | 0.0135204700328 | -0.000485127490028 | 3.99339938521 i |
| 1 | 0.0997156901857 | -0.0144644689488 | 0.221945656786 i |
| 2 | 2.82916586478 | -0.274601275142 | -0.0591732678958 i |

The exact tree Higgs has vanishing same-action mass residual and does not source either zero mode. The coefficient is invariant under arbitrary real orthogonal changes of the three heavy coordinates.

## Independent finite-field check

The verifier also restricts the original potential to the Higgs plus three radial directions, solves its actual radial stationarity equations at finite Higgs amplitude, and reads the induced Majorana mass. With rho=h²/2, minus DeltaMR/h² tends to CHN:

| h / omega | extracted Im(omega CHN / f_raw) | relative error |
|---|---:|---:|
| 0.01 | -0.0848806763395 | 1.352e-04 |
| 0.005 | -0.0848720707143 | 3.381e-05 |
| 0.0025 | -0.0848699190465 | 8.459e-06 |

The error falls quadratically when h is halved, as expected for the next local higher-order term. The source and radial Hessian are separately checked against this restricted original potential, independently of the cached directional-Hessian contractions.

## Improved light rays are a separate order statement

| External Higgs ray | Im(omega CHN / f_raw) | Scope |
|---|---:|---|
| exact_tree_light | -0.0848692011507 | tree matching |
| bosonic_improved_projection | -0.252600617021 | frozen tree vertices on improved light ray; mixed order |
| local_case_0_mixed_order_projection | -0.255045272061 | frozen tree vertices on improved light ray; mixed order |
| local_case_1_mixed_order_projection | -0.256536567164 | frozen tree vertices on improved light ray; mixed order |

The last three entries use the actual improved light directions but frozen tree radial masses and vertices. They are useful matched-input diagnostics at the same external ray as the local Yukawas and type-II projection, not loop-complete Wilson coefficients. The roughly threefold change from the tree ray reflects Higgs alignment dependence, not a numerical precision fit.

## Consequence for sequential seesaw

The CHN=0 running subsector is mathematically invariant, but it is not the tree matching boundary of this P54 action. CHN and the Higgs quadratic coefficient must therefore be included in the sequential EFT flow. The initial CHN matrix remains correlated with the existing f_raw; it is not a new profiling nuisance.

This nonderivative scalar tree graph cannot generate the sterile hypercharge dipole CBN: it contains no external field strength or spin-tensor vertex. The declared model has no additional heavy fermion mixing that would create such a tree dipole. CBN=0 is also closed by the one-loop dimension-five RGE. Finite one-loop dipole matching, axion-loop effects, scalar finite thresholds and complete pole matching are not asserted here.

| Check | Residual | Pass |
|---|---:|:---:|
| actual common-phase Majorana mass reconstructed by scalar derivative | 2.220e-16 | yes |
| ten-Higgs tensor has no NN scalar-singlet coupling | 0.000e+00 | yes |
| SM-singlet subspace has three radial modes and two massless modes | 0.000e+00 | yes |
| positive radial inverse excludes the exact massless subspace | 4.599e-19 | yes |
| cached imaginary jets are exact hypercharge rotations | 2.545e-15 | yes |
| exact tree Higgs vector is a zero mode of the same scalar action | 8.168e-16 | yes |
| exact tree source does not source gauge/PQ massless modes | 4.079e-20 | yes |
| tree Wilson coefficient is heavy-basis invariant | 2.914e-16 | yes |
| same-action radial stationary equations vanish at the reference vacuum | 6.028e-16 | yes |
| cached radial Hessian agrees with independently restricted action | 2.273e-15 | yes |
| restricted-potential cubic source at h=0.01 | 1.142e-11 | yes |
| actual heavy radial valley is stationary at h=0.01 | 1.993e-16 | yes |
| restricted-potential cubic source at h=0.005 | 1.651e-11 | yes |
| actual heavy radial valley is stationary at h=0.005 | 1.695e-17 | yes |
| restricted-potential cubic source at h=0.0025 | 1.931e-10 | yes |
| actual heavy radial valley is stationary at h=0.0025 | 5.310e-16 | yes |
| finite-field Majorana response converges to tree CHN coefficient | 8.459e-06 | yes |
| actual P54 tree CHN is nonzero | 0.000e+00 | yes |
