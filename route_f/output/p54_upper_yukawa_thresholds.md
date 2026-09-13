# P54 upper-PS finite pure-Yukawa threshold

44/44 executable checks pass. The calculation uses the actual same-action stationary PS saddle, all 55 physical upper heavy scalar eigenstates, and the common Spin(10) Clifford intertwiners. Family benchmarks are synthetic. This is not full upper/lower finite matching or a flavor fit.

## Main result: the four-matrix finite boundary is not closed

The original holomorphic H,F,L,R corrections preserve F=sqrt(2)R and L=R in this pure-Yukawa subset. However, the actual split heavy real scalars generate nonzero conjugate-bidoublet Yukawa operators for H and F. Discarding the component orthogonal to the original four tensors would lose genuine finite matching terms. No fundamental field or independent UV family parameter has been added.

## Diagram derivation and normalization

Use $\mathcal L_Y=-\frac12Y^a_{IJ}\psi_I\psi_Js_a+\mathrm{h.c.}$ with canonical real scalars, $Y^{aT}=Y^a$, and $Z_\psi=I+K_\psi$. All fermions are massless at $\sigma=0$. A physical heavy scalar eigenstate $b$ has $Y_b=\sum_A(E_H)_{Ab}Y_A$ and squared mass $M_b^2>0$.

The massless fermion/heavy-scalar self-energy has numerator $\slashed{k}$. Combining $k^2[(k-p)^2-M_b^2]$ with parameter $t$ and shifting $k\mapsto\ell+(1-t)p$ gives the finite $\slashed p$ coefficient $-h(M_b^2,0)Y_b^\dagger Y_b$, where

$$h(M^2,0)=\frac1{16\pi^2}\int_0^1(1-t)\log\frac{(1-t)M^2}{\mu^2}\,dt=\frac{\frac12\log(M^2/\mu^2)-\frac14}{16\pi^2}.$$

In the zero-momentum triangle, two massless Weyl numerators contract to $k^2$; the remaining integral is $[k^2(k^2-M_b^2)]^{-1}$. The oriented Yukawa chain is $Y_bY_a^\dagger Y_b$. Its finite 1PI correction relative to the declared $-Y_a/2$ vertex is $-f(M_b^2,0)Y_bY_a^\dagger Y_b$, with

$$f(M^2,0)=-\frac1{16\pi^2}\int_0^1\log\frac{(1-t)M^2}{\mu^2}\,dt=\frac{1-\log(M^2/\mu^2)}{16\pi^2}.$$

The $1/2$ in the symmetric Weyl action is canceled by differentiating its two fermions; no extra Majorana factor is inserted. Combining the two external fermion legs gives

$$K_\psi=-\sum_{b\in H}Y_b^\dagger Y_b h_b,\qquad \delta Y_a^{\rm 1PI}=-\sum_{b\in H}Y_bY_a^\dagger Y_b f_b,$$

$$\Delta Y_a=\delta Y_a^{\rm 1PI}-\frac12(K_\psi^TY_a+Y_aK_\psi).$$

The canonical field operation and scalar kernels agree with [Patel and Shukla, Eqs. (22)–(24)](https://arxiv.org/html/2310.16563v2). No SU(5) group coefficient is imported. The actual full UV and PS Weyl tensors independently check the matching-scale slope against [Luo, Wang and Xiao](https://arxiv.org/abs/hep-ph/0211440). The subtraction follows zero-momentum hard amplitude matching as in [Gabelmann, Mühlleitner and Staub](https://doi.org/10.1140/epjc/s10052-019-6570-5).

Indeed $\partial_{\log\mu}f=2/(16\pi^2)$ and $\partial_{\log\mu}h=-1/(16\pi^2)$. Consequently

$$\partial_{\log\mu}\Delta Y_a=-\frac1{16\pi^2}\left[\frac12(S_H^TY_a+Y_aS_H)+2\sum_{b\in H}Y_bY_a^\dagger Y_b\right],\quad S_H=\sum_{b\in H}Y_b^\dagger Y_b.$$

The bracket is precisely the heavy part of the pure-Yukawa beta function here. The scalar-wavefunction fermion loop contains only retained massless fermions and cancels between UV and EFT. Trace mixing between a heavy six scalar and an active parent vanishes by PS symmetry. The code checks the complete UV-minus-EFT contraction, not only derivatives of assumed kernels.

## Actual upper spectrum

The matching scale is mu/reference-omega = 0.563058770439. The same-action saddle is {'omega': 1.0008103059215108, 'sigma': 0.0, 'vs': 0.2542776970837077}.

| Squared mass / reference omega² | Real multiplicity | Raw h norm² | Raw f norm² |
|---:|---:|---:|---:|
| 0.102882040959 | 1 | 0 | 8.1852e-34 |
| 0.190903147974 | 6 | 48 | 0 |
| 0.200245831815 | 6 | 48 | 0 |
| 0.260340433542 | 20 | 0 | 1.20713e-28 |
| 0.372275090654 | 6 | 0 | 96 |
| 0.415462413403 | 6 | 0 | 96 |
| 0.540956581209 | 9 | 0 | 2.55677e-28 |
| 2.83317056101 | 1 | 0 | 4.46816e-31 |

Only 24 of the 55 real heavy modes couple to fermions: the split real components of the 10 and 126 PS six parents. Gauge Goldstones and the retained PQ phase have zero Yukawa couplings at this saddle. The retained active tachyonic directions are not inserted into upper heavy logarithms.

## Finite mirror coefficients and exact split-log structure

Set $h=h_{\rm raw}$, $f=f_{\rm raw}$ and $H=\sqrt2h$, $F=4f$, $L=R=2\sqrt2f$ at tree level. The computed common fermion metric is $K_\psi=I_{16}\otimes k$,

$$k=k_h h^\dagger h+k_f f^\dagger f,\qquad k_h=-3[h(m_{h,-}^2,0)+h(m_{h,+}^2,0)],\quad k_f=-6[h(m_{f,-}^2,0)+h(m_{f,+}^2,0)].$$

At the actual scale: k_h=0.018681538032514, k_f=0.010809708789426.

The mirror invariant basis is constructed from the actual six-mode Clifford contraction, normalized to the original parent tensor norm, and its largest entry made positive real. Its full complex tensors are exported in JSON. In that recorded basis:

$$\widetilde H=A_{Hh}hh^\dagger h+A_{Hf}fh^\dagger f,\qquad \widetilde F=A_{Fh}hf^\dagger h+A_{Ff}ff^\dagger f.$$

| Coefficient | Actual complex value | Magnitude from actual group contraction |
|---|---:|---|
| A_Hh | -7.86030569023e-20+0.00128368533616j | 3 sqrt(2) Delta_h |
| A_Hf | 9.01333892534e-19-0.00589775824991j | 6 sqrt(2) Delta_f |
| A_Fh | -7.41076727434e-20-0.00121027014148j | 4 Delta_h |
| A_Ff | -2.10620271079e-21+0.00556045980308j | 8 Delta_f |

Here $\Delta_x=\log(m_{x,+}^2/m_{x,-}^2)/(16\pi^2)$. All four coefficients are derived from the group tensors, not fitted to family data. They are matching-scale independent and vanish when the two real masses of each complete complex six parent coincide. Both random complex symmetric family benchmarks verify the universal cubic polynomial map.

The displayed i is a tensor-basis phase, not a new CP parameter. The PQ/axion phase is also transported: under the historical q(S)=-4 orbit by theta, the mirror coefficients acquire exp(-4i theta). Thus these operators are not explicit PQ breaking; their background spurion must be retained if the axion is restored as a dynamical field.

## What changes and what does not

The holomorphic four-matrix projection alone still preserves the protected F=sqrt(2)R relation in this subset. The genuinely new result is that this projection is not the entire finite boundary. The new mirror matrices are calculable nonlinear functions of the original two UV family matrices; they are not free nuisance matrices. They can change the tree two-spurion sum-rule geometry, but their observed magnitudes do not establish a successful flavor fit. Their effect inside a one-loop beta function begins at two-loop order; a resummed system that includes them must extend its tensor basis.

## Complete dimension-four PS tensor-space closure

Once the nonzero mirror boundary is included, the original beta_closed(H,F,L,R) is not the complete flow. Its paired-holomorphic vertex cancellation is no longer valid. The new full_real_weyl_beta implementation retains the complete real-Weyl contraction:

$$16\pi^2\beta_{Y_a}=\frac12(S^TY_a+Y_aS)+2\sum_bY_bY_a^\dagger Y_b+\sum_bY_b\operatorname{Re}\operatorname{Tr}(Y_a^\dagger Y_b)-3\sum_i g_i^2(C_i^TY_a+Y_aC_i),\qquad S=\sum_bY_b^\dagger Y_b.$$

The existing 248 real scalars are unchanged. The code adds the two already-derived mirror intertwiners and tests completely general complex H,F,mirrorH,mirrorF, with symmetric L,R, at an unequal-gL/gR and unequal-L/R point. It compares the full beta tensor with its projection onto all six invariant families. Three literal vertex sums independently check the optimized tensor-index contraction; the zero-mirror limit reproduces the previous verifier.

The absence of further dimension-four Yukawa tensors also follows directly from representation theory. With $\psi_L=(4,2,1)$ and $\psi_R=(\overline4,1,2)$, their mixed product is $(1\oplus15,2,2)$. Both H and F are real-type PS representations carried by complex scalar fields, so each admits two conjugate couplings. The symmetric left-left gauge product is $(10,3,1)\oplus(6,1,1)$; only its conjugate triplet parent is retained. The right-right statement is conjugate. Thus the complete family tensor space has four general matrices and two symmetric matrices:

$$\dim_\mathbb C\mathcal Y_{\rm PS}=4n_f^2+2\frac{n_f(n_f+1)}2=48\quad(n_f=3).$$

This counts the allowed low-energy tensor space, not independent fundamental UV parameters. PS symmetry and renormalizability force the dimension-four Yukawa counterterms to remain in this space. Higher-dimensional operators and their feedback are outside this statement.

Generic six-invariant closure residual: 4.954e-15. At the actual finite-matched benchmark the residual is 4.942e-15; discarding the two mirror directions leaves beta-tensor norm 0.000308811384.

The full gauge Yukawa trace is checked too: the mirror H/F norm contributions must be included in Y4. The old fast path now rejects, rather than silently mis-evolves, this non-holomorphic-pair input. The new function supplies the full beta tensor and a verified six-invariant projection; it does not claim a completed global trajectory/fit.

## Missing-diagram ledger

| Graph | Status |
|---|---|
| heavy physical scalar plus massless fermion self-energy | computed with all55 eigenstates;24 nonzero Yukawa states |
| two massless fermions plus one heavy scalar triangle | computed with split real masses and complete relative phases |
| pure-Yukawa light scalar kinetic graph | zero hard matching: both fermions are retained and massless; UV minus EFT cancels identically |
| scalar-cubic plus two-Yukawa triangle | zero at dimension-four and zero momentum with massless internal fermions: requires chirality-flip mass or external momentum; not a missing pure-Yukawa cubic term |
| heavy-vector fermion kinetic and Yukawa vertex graphs | not computed here; gauge-fixing package and BFG-Landau conventions required |
| heavy-vector/scalar, Goldstone and ghost scalar kinetic package | not computed here |
| scalar-cubic contribution to upper scalar kinetic metric | not computed here; separate from pure-Yukawa matching |
| upper heavy-radial/tadpole shift and full scalar matching | not computed here; tree stationary background used, one-loop mass insertion omitted consistently |
| lower staged PS-to-SM functional and sequential Weinberg finite matching | not computed here |

The pure-Yukawa subset is independent of gauge fixing, but this does not provide the missing background-field Landau vector/Goldstone/ghost package. No unknown finite contribution is silently set to zero. One-loop mass-counterterm insertions in these one-loop diagrams would be selected two-loop terms and are not included. The separate broken-background four-doublet scalar-cubic kinetic calculation is a one-step hard-scalar subset; it is neither this upper-only Yukawa calculation nor complete lower staged matching.

## Executable validation

| Check | Residual | Pass |
|---|---:|:---:|
| same-action upper finite prerequisite | 0.000e+00 | yes |
| heavy55 spectrum reproduces stored upper spectrum | 4.044e-16 | yes |
| actual canonical active and heavy planes are orthogonal | 0.000e+00 | yes |
| all 55 integrated scalar eigenvalues are positive | 0.000e+00 | yes |
| only 24 real six-parent scalar modes have Yukawa vertices | 0.000e+00 | yes |
| all 48 Weyl fermions are massless at the actual PS saddle | 0.000e+00 | yes |
| heavy Goldstone and PQ planes have zero Yukawa couplings | 0.000e+00 | yes |
| actual fermion kinetic correction is Hermitian | 3.796e-19 | yes |
| heavy six contractions give a common Spin10-family kinetic matrix | 2.000e-18 | yes |
| two noncommuting complex family spurions reconstruct finite mirror coefficients | 1.235e-19 | yes |
| four original plus two actual mirror invariants reconstruct full finite result | 2.521e-17 | yes |
| finite mirror operators are genuinely outside the original four-matrix space | 0.000e+00 | yes |
| original holomorphic F=sqrt2 R survives this pure-Yukawa subset | 2.211e-18 | yes |
| original holomorphic L=R survives this pure-Yukawa subset | 6.776e-21 | yes |
| both generated mirror family matrices are symmetric | 5.971e-20 | yes |
| finite matching-scale slope equals minus actual UV-minus-EFT beta | 8.070e-15 | yes |
| finite mirror contributions are matching-scale independent | 8.063e-19 | yes |
| single real heavy mode has the Weyl vertex factor two and sign | 6.007e-16 | yes |
| degenerate complete heavy complex pairs cancel finite vertices | 3.097e-19 | yes |
| independent zero-momentum Feynman-parameter integration | 2.602e-18 | yes |
| arbitrary complex family U3 covariance of finite vertices and legs | 3.210e-18 | yes |
| arbitrary complex family U3 covariance of kinetic matching | 4.321e-18 | yes |
| CP covariance with conjugation of couplings AND group tensors | 0.000e+00 | yes |
| nondegenerate heavy real-basis covariance with full mass rotation | 7.903e-18 | yes |
| all 124 complex scalar phase choices transport finite mirror vertices | 1.148e-18 | yes |
| actual PQ orbit of heavy basis transports all Yukawa phases | 4.813e-17 | yes |
| generated mirrors carry the required PQ spurion phase | 8.037e-19 | yes |
| PQ axion motion leaves fermion kinetic matching invariant | 7.307e-19 | yes |
| independent second complex family point verifies universal cubic coefficients | 2.581e-20 | yes |
| PS intertwiner prerequisite has no hidden normalization gate | 6.280e-16 | yes |
| all 21 PS generators preserve full finite Yukawa tensors | 6.626e-19 | yes |
| actual Clifford group factors equal 3sqrt2,6sqrt2,4,8 times split logs | 8.674e-18 | yes |
| actual h kinetic group factor and finite -1/4 constant | 3.469e-18 | yes |
| actual f kinetic group factor and finite -1/4 constant | 6.051e-22 | yes |
| full real-Weyl beta agrees with original fast path when mirrors vanish | 4.074e-19 | yes |
| full mirror beta matches literal Weyl diagram sum at scalar 0 | 7.201e-17 | yes |
| full mirror beta matches literal Weyl diagram sum at scalar 19 | 8.122e-17 | yes |
| full mirror beta matches literal Weyl diagram sum at scalar 181 | 5.558e-17 | yes |
| six Yukawa invariant families close under generic complex off-parity one-loop flow | 4.954e-15 | yes |
| original four-matrix projection is not a complete beta after mirror matching | 0.000e+00 | yes |
| old complex-pair fast path fails closed on actual mirror tensor structure | 0.000e+00 | yes |
| actual finite-matched boundary beta is reconstructed by all six invariants | 4.942e-15 | yes |
| full mirror flow has independent U3L times U3R covariance | 2.704e-16 | yes |
| full mirror gauge Y4 includes both conjugate bidoublet norms | 1.259e-15 | yes |

API: build_upper_geometry(); calculate(h_raw,f_raw,geometry,mu); threshold_from_tensors accepts either diagonal heavy squared masses or a full real symmetric heavy mass matrix; build_mirror_basis/project_mirror/assemble_mirror preserve the generated invariant directions. The source hashes and original heavy-Hessian cache key are recorded in JSON.
