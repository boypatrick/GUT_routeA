# P54 Pati-Salam Yukawa running and matching audit

Date: 2026-09-05. Status: research/implementation contract; no flavor fit or replacement P2 scale solution is claimed.

The present action is sufficient in principle to define the calculation, but the existing gauge beta table, intermediate active-field projector, and Yukawa normalization do not yet form one consistent executable EFT. The numerical scalar Hessian remains useful input. The immediate repair is to close these interfaces before using a Pati-Salam (PS) evolution in a physical likelihood.

## 1. Independently reproducible beta-table discrepancy

The local PS parent census is in `p54_p2_two_site_matching.json`. Its real dimensions and Casimirs determine the complex scalar indices without any mass fitting:

\[
T_i(R)=\frac{d_{\mathbb C}(R)C_i(R)}{d(G_i)},\qquad
d_{\mathbb C}=\tfrac12d_{\mathbb R}
\]

for the parents belonging to the complex 10 and 126. Group order is \((4,L,R)\).

| Active parent | Complex dimension | \(T_4,T_L,T_R\) | \(C_4,C_L,C_R\) |
|---|---:|---|---|
| \(\phi_{10}:(1,2,2)\) | 4 | \((0,1,1)\) | \((0,3/4,3/4)\) |
| \(\Sigma:(15,2,2)\) | 60 | \((16,15,15)\) | \((4,3/4,3/4)\) |
| \(\Delta_L:(10,3,1)\) | 30 | \((9,20,0)\) | \((9/2,2,0)\) |
| \(\Delta_R:(\overline{10},1,3)\) | 30 | \((9,0,20)\) | \((9/2,0,2)\) |

For three generations of Weyl \((4,2,1)+(\overline4,1,2)\), \(\sum_F T_i=(6,6,6)\). With complex scalar counting,

\[
a_i=-\frac{11}{3}C_i(G)+\frac23\sum_FT_i(F)+\frac13\sum_ST_i(S),
\]

\[
b_{ij}=-\frac{34}{3}C_i(G)^2\delta_{ij}
+\sum_F\left[2C_j(F)+\frac{10}{3}C_i(G)\delta_{ij}\right]T_i(F)
+\sum_S\left[4C_j(S)+\frac23C_i(G)\delta_{ij}\right]T_i(S).
\]

These are the gauge-and-matter terms; Yukawa-dependent two-loop gauge terms must be added when the Yukawa matrices are specified. The general renormalization framework is given by [Luo, Wang and Xiao](https://arxiv.org/abs/hep-ph/0211440).

Direct rational evaluation gives

\[
a_{\rm four\ parents}=\left(\frac23,\frac{26}{3},\frac{26}{3}\right),
\quad
b_{\rm four\ parents}=\begin{pmatrix}
3551/6&249/2&249/2\\
1245/2&779/3&48\\
1245/2&48&779/3
\end{pmatrix}.
\]

The historical code instead stores \(a_4=1\), \(b_{44}=1209/2\). Exactly one extra complex \((6,1,1)\) produces the difference

\[
\Delta a=(1/3,0,0),\qquad \Delta b_{44}=38/3,
\]

with all other entries zero. Gauge indices alone cannot identify whether that six belongs to the 10 or 126. Either choice brings additional Yukawa interactions and requires an explicit EFT decision.

The discrepancy is inherited rather than introduced by the local numerical solver: [Babu and Khan, Eq. (7)](https://arxiv.org/pdf/1507.06712) prints the historical table, while their Table 1 assigns both sixes to the GUT scale. The local audit reproduces the discrepancy; it does not treat agreement with a printed table as an independent physics validation.

## 2. A stronger matching-scale covariance failure

This finding does not depend on the extra-six ambiguity. The existing `verify_p54_p2_two_site_matching.py` selects only `Sigma126:(10-pair,1,3)` as its intermediate scalar projector and assigns its complement to the upper threshold. That parent is an \(SU(2)_L\) singlet. The intermediate heavy vectors are also \(SU(2)_L\) singlets. Therefore its entire lower-scale threshold obeys

\[
\lambda_{I,L}(\mu)=0,\qquad
\frac{d\lambda_{I,L}}{d\log\mu}=0.
\]

For the convention already used by P2,

\[
\alpha_{2,\rm SM}^{-1}=\alpha_{L,\rm PS}^{-1}
-\frac{\lambda_{I,L}}{12\pi},\qquad
\frac{d\alpha_i^{-1}}{d\log\mu}=-\frac{a_i}{2\pi}+O(\alpha),
\]

scale covariance requires

\[
\lambda'_{I,L}=-6(a_L^{\rm PS}-a_2^{\rm SM})
=-6\left(\frac{26}{3}+\frac{19}{6}\right)=-71.
\]

The historical interface instead leaves

\[
\frac{d}{d\log\mu}\left[
\alpha_{2,\rm SM}^{-1}-\alpha_{L,\rm PS}^{-1}
+\frac{\lambda_{I,L}}{12\pi}\right]
=\frac{71}{12\pi}=1.88333349325.
\]

An independent scalar count gives the same result: remove the two complete PS bidoublets and \(\Delta_L\), retain one SM complex Higgs doublet, and obtain

\[
\Delta a_L=\frac13+5+\frac{20}{3}-\frac16=\frac{71}{6}.
\]

Thus a root of the current unification equations is not yet a matching-scale-independent EFT result. This is a repairable interface error, not a proof that the full P54 action has no viable point.

The new verifier `route_f/code/verify_p54_ps_active_census.py` reconstructs the indices from the exported Casimirs, evaluates the beta coefficients exactly with rational arithmetic, and compares a finite-difference matching-scale variation against this analytic result. Its ledger is `route_f/output/p54_ps_active_census.json`: **14/14 audit regressions pass; the two historical physical consistency gates fail**. It leaves every historical P2 output unchanged.

### Minimal repair in the same action: all gauge directions

Choose the four active parents in the first table. Their Casimir triples, together with the original field identifier, define the PS projector unambiguously:

\[
P_A=P_{\phi;(0,3/4,3/4)}+P_{\Sigma;(4,3/4,3/4)}
+P_{\Sigma;(9/2,2,0)}+P_{\Sigma;(9/2,0,2)}.
\]

Keep these complete parents in the intermediate theory. The upper threshold integrates out the complementary PS representations. At the lower threshold integrate out the retained parents except for the actual light Higgs doublet. This is a possible EFT organization within the same action, conditional on a controlled staged matching of the actual spectrum; it is not a new scalar field content.

For logarithmic counting, write real scalar indices as \(I_S=2T_{\mathbb C}\). The active and complete indices are

\[
I_A^{\rm PS}=(68,72,72),\qquad I_{\rm full}^{\rm PS}=(84,84,84).
\]

The PS-to-SM index map is

\[
\mathcal R(v_4,v_L,v_R)=\left(v_4,v_L,\frac25v_4+\frac35v_R\right),
\qquad I_h^{\rm SM}=(0,1,3/5).
\]

The massive-vector real trace indices follow from the adjoint differences:

\[
I_{V,U}^{\rm PS}=(8-4,8-2,8-2)=(4,6,6),
\]

\[
I_{V,I}^{\rm SM}=\left(4-3,2-2,\frac25\,4+\frac35\,2\right)
=(1,0,14/5).
\]

The last line also follows by counting the PS leptoquark vectors \((3,1,2/3)+\mathrm{h.c.}\) and \(W_R^\pm\). A massive vector removes one real Goldstone from the physical scalar trace, so use the same site-resolved gauge orbit in both sectors. Then

\[
I_{S,U}^{\rm phys}=I_{\rm full}^{\rm PS}-I_A^{\rm PS}-I_{V,U}^{\rm PS}
=(12,6,6),
\]

\[
I_{S,I}^{\rm phys}=\mathcal R I_A^{\rm PS}-I_h^{\rm SM}-I_{V,I}^{\rm SM}
=(67,71,67).
\]

In the conventional one-loop gauge matching trace, a physical real scalar contributes \(+I_S\log(M_S/\mu)\), and the massive vector plus its gauge/ghost sector contributes \(-21I_V\log(M_V/\mu)\). The associated Goldstone must not be counted again as an independent physical scalar. Scale differentiation therefore gives

\[
\lambda'_U=21I_{V,U}-I_{S,U}^{\rm phys}=(72,120,120),
\]

\[
\lambda'_I=21I_{V,I}-I_{S,I}^{\rm phys}=(-46,-71,-41/5).
\]

For the full original Spin(10) action,

\[
a_{10}=-\frac{11}{3}\,8+\frac23\,6+\frac16\,84=-\frac{34}{3}.
\]

The six cancellation identities close exactly:

\[
\lambda'_U=-6\left[a_{10}(1,1,1)-a_{\rm PS}\right],\qquad
\lambda'_I=-6\left[\mathcal R a_{\rm PS}-a_{\rm SM}\right].
\]

These are verified with rational arithmetic in the new ledger. Finite constants, finite logarithmic mass matrices, and new unification scales are still to be calculated. At the simultaneous nonzero \(\omega,\sigma\) background, parent, scalar mass, and Goldstone projectors need not commute. Merely substituting \(P_A\) into the old broken-background weighted trace is not a proof of finite matching. One must derive the two-stage EFT in a common gauge/tadpole scheme, resolve the Goldstone orbits with the vector site projectors, and check that the finite implementation reproduces the six identities above.

## 3. Literature RGE that can be used as a regression, with its actual scope

Use [Meloni, Ohlsson and Riad, Appendix A, Eqs. (70)-(78)](https://arxiv.org/pdf/1612.07973), whose corrected equations supersede earlier minimal-model formulas. Their intermediate theory contains \(\Phi=(1,2,2)\), \(\Sigma=(15,2,2)\), and \(\Delta_R\), without \(\Delta_L\). Their lower EFT has two Higgs doublets. Consequently their equations are a regression target, not a direct P54 replacement.

In their matrix convention write \(H=Y_F^{(10)}\), \(F=Y_F^{(126)}\), and symmetric \(R=Y_R^{(126)}\). Define

\[
A_L=HH^\dagger+\frac{15}{4}FF^\dagger,\qquad
A_R=H^\dagger H+\frac{15}{4}F^\dagger F+\frac{15}{2}R^*R,
\]

\[
G_D=\frac94(g_L^2+g_R^2+5g_4^2),\qquad
G_R=\frac94(2g_R^2+5g_4^2).
\]

Their one-loop equations become

\[
16\pi^2\dot H=A_LH+HA_R+4\operatorname{tr}(HH^\dagger)H-G_DH,
\]

\[
16\pi^2\dot F=A_LF+FA_R+\operatorname{tr}(FF^\dagger)F-G_DF,
\]

\[
16\pi^2\dot R=A_R^TR+RA_R+2\operatorname{tr}(RR^*)R-G_RR.
\]

The corresponding gauge coefficients are \((-7/3,2,26/3)\). This makes the absence of \(\Delta_L\) explicit.

For a parity-preserving P54 interval, \(L=Y_L^{(126)}\) must also be retained. It contributes to the left fermion anomalous dimension and has its own matrix RGE. A parity completion should generate \(L\)-dependent terms and satisfy the left/right exchange identity. If a six is retained to justify the historical beta table, its left-left and right-right Yukawa tensors must also be included. **The complete P54 PS Yukawa coefficients have not yet been generated or independently checked.** Blindly importing the three equations above would omit those interactions.

The conventional upper-scale matching used in the cited PS literature is

\[
H(M_U)=\sqrt2\,h(M_U),\quad
F(M_U)=4\sqrt2\,f(M_U),\quad
R(M_U)=4f(M_U).
\]

Those factors depend on the tensor and kinetic normalization. [Djouadi, Fonseca, Ouyang and Raidal, Eq. (52)](https://arxiv.org/pdf/2212.11315) explicitly discusses this dependence. The local card declares a canonically normalized self-dual tensor, factors \(1/5!\), and the mass convention \(M_R=f\sigma\). These must be matched by an explicit Clifford/intertwiner projection before importing numerical factors from another paper. The common Clebsch ratio \(-3\) is not enough to fix the absolute normalization.

**Same-day implementation update:** `verify_p54_spinor_intertwiners.py` now supplies this absolute audit on the actual P1 `U126`, with 20/20 checks. The raw action-card coefficients give `h_D=sqrt(2)Y10`, `f_D=2Y126/sqrt(3)`, and `MR=4sqrt(2)sigma Y126=2sqrt(6)sigma f_D`. Thus the mass-normalized matrix is `f_M=2sqrt(6)f_D`, and `MR=sigma f_M`. The action card is corrected. The raw LL-triplet coefficient is also `4sqrt(2)`. This closes the absolute-normalization subgate, not the complete PS matching: the actual scalar/spinor global-conjugation dictionary and all four-copy relative phases must still be exported consistently, and finite kinetic/vertex terms remain open.

## 4. Required threshold chain and mathematical matching contract

At \(M_U\), decompose the actual two Spin(10) matrices into every active PS Yukawa tensor, with matching of scalar/fermion kinetic terms and finite one-loop vertices. At \(M_I\), match PS onto **one** light SM doublet and any right-handed neutrinos whose thresholds have not yet been crossed. Retained light-Higgs coefficients must come from the same quantum stationary point and complex mass matrix used for the spectrum.

For an amputated coupling \(Y\) and hard full-minus-EFT corrections, canonical normalization gives

\[
Y^-_a=Y^{(0)}_a+\Delta\Gamma_a^{\rm hard}
-\frac12\left[(\Delta Z_L)^T Y^{(0)}_a
+Y^{(0)}_a\Delta Z_R
+\sum_b(\Delta Z_\phi)_{ba}Y^{(0)}_b\right]+O(\hbar^2).
\]

Here \(Y^{(0)}_a\) is the coupling projected with the complete scalar eigenvector, including conjugate-field assignments. This identity follows by substituting the three inverse square roots of the kinetic matrices into the vertex. It defines the matching convention without guessing finite threshold coefficients.

Writing the tree projection as \(\mathcal P(\mu)Y^+\), the required scale check is

\[
\frac{d\Delta Y}{d\log\mu}
=\beta^-_Y-\mathcal P\beta^+_Y
-\frac{d\mathcal P}{d\log\mu}Y^+ +O(\hbar^2).
\]

The derivative of \(\mathcal P\) is relevant: a scale-dependent doublet mixture is not a fixed numerical constant. A scalar projected CW curvature \(\eta_Y\) alone cannot supply the vertex and fermion wave-function terms.

For sequential type-I thresholds, define \(M_R=U^*\operatorname{diag}(M_i)U^\dagger\), \(y_i=(Y_\nu U)_{:i}\), and choose \(m_\nu^{I}=-v^2\kappa/2\). Tree matching then gives

\[
\kappa(\mu_i^-)=\kappa(\mu_i^+)+\frac{y_i y_i^T}{M_i},
\qquad \mu_i\simeq M_i,
\]

with each removed column deleted from the remaining dynamical Yukawa matrix. The threshold order follows the actual Takagi masses, not an assumption that all three equal \(M_I\). [Antusch et al.](https://arxiv.org/abs/hep-ph/0501272) treat this sequential evolution; [Zhang and Zhou](https://arxiv.org/abs/2107.12133) provide complete one-loop type-I matching, including renormalizable couplings and the Weinberg coefficient. Their SMEFT formula applies only after additional heavy PS fields are removed or their effects are matched separately.

The local action retains type-II seesaw. Its scalar-triplet Yukawa and Higgs-triplet trilinear must be obtained from the actual potential and eigenvectors before solving the heavy-triplet equation of motion. [Li, Zhang and Zhou](https://arxiv.org/abs/2201.05082) provide the corresponding one-loop type-II matching and a machine-readable ancillary notebook. A type-I-only fit is an additional approximation requiring a demonstrated small triplet contribution.

## 5. Existing assets, missing objects, and next implementation order

Available locally: the canonical 328-real scalar potential and its derivatives; SM and PS generators on the scalar space; the PS parent Casimirs; tree stationary vectors; the complete mass-eigenvalue census; the gauge ODE integrator; diagnostic SM matrix Yukawa evolution; complex fermion mass algebra. These provide substantial reusable infrastructure.

Missing before a physical fit:

1. A single active PS census shared by beta generation, the upper/lower matching projectors, and scalar-Yukawa RGEs. Recompute the two-loop gauge coefficients from this census, rather than selecting a printed table by model name.
2. An explicit 16-spinor basis and normalized Yukawa intertwiners for all retained scalar components. Verify the opposite conjugation of the card's \(\phi^*\) and \(\Sigma\) couplings, the \(-3\) ratio, and the absolute \(h,f,R,L\) normalization.
3. The full complex, quantum-corrected light eigenvector and tadpole convention, together with a consistent definition of every renormalized parameter.
4. The complete PS matrix beta system, including \(\Delta_L\) and any chosen six, and Yukawa terms in the two-loop gauge running. The general symbolic machinery in [PyR@TE 3](https://arxiv.org/abs/2007.12700) and its [official repository](https://github.com/LSartore/pyrate) can generate this system from an explicit model file; no model file was generated or installed in this audit.
5. Finite hard vertex and kinetic matching at \(M_U\), \(M_I\), and the actual seesaw thresholds. Each interface must pass its logarithmic derivative identity before covariance propagation is interpreted physically.

The smallest useful next step is the active-census and logarithmic matching repair. It determines which Yukawa RGE is actually required. A PS-aware flavor fit becomes justified after these objects exist; replacing the missing evolution by a fixed 10% or 30% box does not prove the needed error bound.
