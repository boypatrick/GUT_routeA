# P54: complete four-parent Pati–Salam Yukawa evolution

Status: 29/29 executable checks pass. This closes the all-active one-loop Yukawa contraction/ODE task, not a physical fit or finite matching.

## Convention and derivation

Use canonical real scalars and two-component Weyl fields: $\mathcal L_Y=-\frac12Y^a_{IJ}\psi_I\psi_J\varphi_a+\mathrm{h.c.}$, $Y^{aT}=Y^a$. For $S=\sum_bY^{b\dagger}Y^b$ the one-loop coefficient is

$$\mathcal B^a=\tfrac12(S^TY^a+Y^aS)+2\sum_bY^bY^{a\dagger}Y^b+\sum_bY^b\operatorname{Re}\operatorname{Tr}(Y^{a\dagger}Y^b)-3\sum_i g_i^2(C_i^TY^a+Y^aC_i),\qquad16\pi^2\dot Y^a=\mathcal B^a.$$

The tensor formula and Weyl factor are from [Luo, Wang and Xiao, Eqs. (19), (30)–(34)](https://arxiv.org/pdf/hep-ph/0211440). The transposed left wavefunction follows from symmetric Weyl tensors; replacing it by S fails complex family covariance.

Actual Clifford tensors are projected by PS Casimirs. Retain $H:(1,2,2)$, $F:(15,2,2)$, $L:(\overline{10},3,1)$, $R:(10,1,3)$ with complex dimensions $4,60,30,30$. Each complex scalar gives $Y_R=K/\sqrt2$, $Y_I=s iK/\sqrt2$, with $s=-1$ for $\phi^*$ and $s=+1$ for $\Sigma$. Thus $Y_RX Y_R+Y_IXY_I=0$ for any $X$: the full vertex sum cancels. This cancellation requires complete complex pairs. The executable fast path checks even scalar count and the pairwise relation with a scale-relative tolerance, rejecting incomplete or nonholomorphic pairs. After arbitrary real-mode mass decoupling the full vertex contraction must be retained.

These are physical labels after the common global conjugation, with $4=(3,B-L=1/3)+(1,B-L=-1)$. The actual color-singlet scalar hypercharges are $(1,1,1)$ for L and $(-2,-1,0)$ for R. Thus the historical literal-P1 $10/\overline{10}$ names must be swapped; their Casimirs, indices and beta coefficients are unchanged.

The following group contractions are computed from the actual normalized K tensors, with family matrices factored out. The first two columns are eigenvalues of sum K†K on eight L/R gauge states; the last is the norm of any complex scalar component.

| Parent | left S | right S | Tr K†K per complex scalar |
|---|---:|---:|---:|
| H | 2 | 2 | 8 |
| F | 15/2 | 15/2 | 2 |
| L | 15 | 0 | 4 |
| R | 0 | 15 | 4 |

This fixes the PS units $H=\sqrt2h_{\rm raw}$, $F=4f_{\rm raw}$, $L=R=2\sqrt2f_{\rm raw}$ at the tree Spin(10) boundary. In the common component-phase convention $f_D=2f_{\rm raw}/\sqrt3$ and $f_M=+i4\sqrt2f_{\rm raw}$: $F=2\sqrt3f_D$, $R=\sqrt6f_D$, $f_M=+2iR$. These are normalization-aware boundary maps, not four new UV family matrices.

## Complete one-loop matrix system

For general complex H,F and symmetric L,R define

$$A_L=HH^\dagger+\frac{15}{4}FF^\dagger+\frac{15}{2}LL^\dagger,\qquad A_R=H^\dagger H+\frac{15}{4}F^\dagger F+\frac{15}{2}R^\dagger R,$$

$$G_D=\frac94(g_L^2+g_R^2+5g_4^2),\quad G_L=\frac94(2g_L^2+5g_4^2),\quad G_R=\frac94(2g_R^2+5g_4^2).$$

$$\begin{aligned}16\pi^2\dot H&=A_LH+HA_R+[4\operatorname{tr}(H^\dagger H)-G_D]H,\\16\pi^2\dot F&=A_LF+FA_R+[\operatorname{tr}(F^\dagger F)-G_D]F,\\16\pi^2\dot L&=A_LL+LA_L^T+[2\operatorname{tr}(L^\dagger L)-G_L]L,\\16\pi^2\dot R&=A_R^TR+RA_R+[2\operatorname{tr}(R^\dagger R)-G_R]R.\end{aligned}$$

Realification of the scalar norm column yields trace coefficients $(4,1,2,2)$; the fermion columns yield $S_{LL}=2I_8\otimes A_L^T$, $S_{RR}=2I_8\otimes A_R$. Gauge terms follow from $C_4(4)=15/8$, $C_2(2)=3/4$ and the $-3$ anticommutator. The code compares all entries of 248 real-scalar 48×48 beta matrices with the compact system, at a generic complex unequal-L/R point.

Deleting L reproduces the corrected [Meloni, Ohlsson and Riad, Eqs. (76)–(78)](https://arxiv.org/pdf/1612.07973). This is a regression only; P54 retains DeltaL. An independent SM-top tensor test recovers 9/2 and the GUT-normalized gauge coefficients 17/20, 9/4, 8.

## Symmetry and a protected Spin(10) subflow

The equations are covariant under $H,F\mapsto U_L^T(H,F)U_R$, $L\mapsto U_L^TLU_L$, $R\mapsto U_R^TRU_R$, and under $H,F\mapsto H^T,F^T$, $L\leftrightarrow R$, $g_L\leftrightarrow g_R$. So symmetric $H,F$ with $L=R,g_L=g_R$ form an invariant parity locus.

There is a stronger result at one loop. Impose additionally $F=\sqrt2R$. Then $A_R^T=A_L$, $G_D=G_L=G_R$, and $\operatorname{tr}F^\dagger F=2\operatorname{tr}R^\dagger R$. Substitution gives

$$\mathcal B_F-\sqrt2\mathcal B_R=0,\quad \mathcal B_L-\mathcal B_R=0,\quad\mathcal B_H^T=\mathcal B_H,\quad\mathcal B_F^T=\mathcal B_F.$$

The gauge equations also preserve $g_L=g_R$, because $a_L=a_R$, $b_{L4}=b_{R4}$, $b_{LL}=b_{RR}$, $b_{LR}=b_{RL}$, and $Y_{4,L}=Y_{4,R}$. Uniqueness of the smooth perturbative ODE then proves invariance of this entire submanifold. The actual tree Spin(10) tensor boundary lies on it; an explicit downward numerical flow preserves it. This is conditional one-loop protection, not a theorem about finite threshold matching or two-loop Yukawa flow.

An unequal-gL/gR negative-control flow breaks the proportionality. With unequal Majorana parents, even initially symmetric H,F acquire antisymmetric beta sources. Therefore the general implementation retains four independently evolved matrices; it does not enforce a possibly false low-scale boundary relation.

## Upper boundary versus lower matching

At the upper tree boundary use $H=\sqrt2h_{\rm raw}$, $F=4f_{\rm raw}$, $L=R=2\sqrt2f_{\rm raw}$. At the lower scale use the evolved $H$, $f_D=F/(2\sqrt3)$ for Dirac projections; in the fixed common component phase $M_R=+2i\,\sigma_{\rm dimful}R$, $M_L=-2i\,\Delta_{L,\rm dimful}L$. If the scalar card writes a dimensionless $\sigma$ times $\omega$, then $\sigma_{\rm dimful}=\sigma\omega$. The light-doublet coefficients $(a,b,d,e)$ multiply H and $F/(2\sqrt3)$; type I and II use R and L. The ratio fM/fD is preserved only on the protected locus just proved. Finite-matched off-locus data must be propagated by the full system.

## Two-loop gauge Yukawa term

The actual Weyl trace $Y_{4,i}=\operatorname{Tr}[C_i\sum_aY^aY^{a\dagger}]/d(G_i)$ gives

$$Y_{4,4}=4\|H\|_F^2+15\|F\|_F^2+15\|L\|_F^2+15\|R\|_F^2,\quad Y_{4,L}=4\|H\|_F^2+15\|F\|_F^2+30\|L\|_F^2,\quad Y_{4,R}=4\|H\|_F^2+15\|F\|_F^2+30\|R\|_F^2.$$

$$\dot g_i=\frac{a_ig_i^3}{16\pi^2}+\frac{g_i^3}{(16\pi^2)^2}\left[\sum_jb_{ij}g_j^2-Y_{4,i}\right],\qquad a=(2/3,26/3,26/3).$$

The b matrix is the corrected four-parent census, with b44=3551/6, not the historical extra-six table. The Yukawa gauge trace is checked against the full tensor and independently against SM-top.

## Executable interface

build_ps_geometry() exports the scalar embeddings, actual group tensors and fermion basis. assemble_real_yukawas(couplings, geometry) returns 248×48×48 tensors in real-scalar order H,F,L,R, with real/imaginary coordinates interleaved. The 328×248 old-coordinate embedding and 16×16 fermion basis are included. beta_closed returns 16π² times the one-loop Yukawa beta. beta_gauge returns dg/dlog(mu), including gauge two-loop and Yukawa traces. integrate performs a coupled ODE without a fit.

| Check | Residual | Pass |
|---|---:|:---:|
| actual PS Casimir eigenspaces | 6.280e-16 | yes |
| four parents have only the allowed chiral fermion blocks | 0.000e+00 | yes |
| physical conjugated triplet labels from actual colorless hypercharges | 9.133e-16 | yes |
| 248 canonical real scalar directions are orthonormal | 1.552e-16 | yes |
| full Weyl Yukawa matrices are complex symmetric | 0.000e+00 | yes |
| all-active generic tensor beta equals closed four-matrix system | 2.277e-16 | yes |
| actual raw normalization and all parent contraction constants | 1.776e-15 | yes |
| complete complex scalar pairs cancel the vertex contraction | 0.000e+00 | yes |
| holomorphic fast-path guard rejects incomplete and nonholomorphic pairs | 0.000e+00 | yes |
| independent U(3)L times U(3)R family covariance | 1.811e-16 | yes |
| actual intertwiners implement family covariance | 1.911e-16 | yes |
| left-right exchange covariance including unequal gauge couplings | 2.861e-17 | yes |
| Majorana beta matrices remain complex symmetric | 0.000e+00 | yes |
| canonical scalar phase-basis covariance of generic beta | 1.620e-16 | yes |
| Yukawa-only tensor contraction independently separates gauge factors | 1.716e-15 | yes |
| DeltaL-absent limit agrees with corrected literature equations | 2.733e-16 | yes |
| two-loop gauge Yukawa trace agrees with closed Y4 formula | 7.693e-16 | yes |
| two-loop gauge Yukawa trace is family-basis invariant | 2.683e-16 | yes |
| SM top-only 9/2 cubic and 8,9/4,17/20 gauge coefficients | 1.337e-16 | yes |
| SM top-only two-loop gauge Yukawa coefficients | 2.266e-16 | yes |
| nonzero numerical ODE respects independent family basis covariance | 1.537e-16 | yes |
| numerical ODE preserves Majorana symmetry | 0.000e+00 | yes |
| short ODE actually evolves all four nonzero matrices | 0.000e+00 | yes |
| Spin10 boundary exactly reconstructs actual raw common-phase action | 1.022e-16 | yes |
| actual Spin10 boundary flow preserves symmetric bidoublets | 0.000e+00 | yes |
| actual Spin10 boundary flow preserves L=R and gL=gR | 0.000e+00 | yes |
| parity locus protects F=sqrt2 R at one loop | 3.292e-17 | yes |
| unequal L/R produces a genuine antisymmetric bidoublet beta source | 0.000e+00 | yes |
| unequal gL/gR negative control breaks protected 126 ratios | 0.000e+00 | yes |

## Remaining physical gates

A beta function does not provide finite Yukawa threshold matching at the actual upper/lower spectra. The common action reconstruction verifies the tree boundary only. Sequential Majorana/type-II matching, the actual loop-corrected light-doublet projection and the global flavor/seesaw fit remain separate. All family entries in these ODE tests are synthetic; no physical fit, Higgs pole, determinant or portal is promoted.
