# Route G G1: common circle spectrum and fixed probe channels

This is a bounded free Gaussian-model verification, not a particle fit. No measured masses or Route-F benchmark values were used.

## Declared model

For a periodic complex scalar on fixed $M_4\times S^1$, $D_\theta=\partial_\theta+i\alpha$ and $\Phi=\sum_n\phi_n e^{in\theta}/\sqrt{2\pi R}$ give $m_n^2=M_5^2+(n+\alpha)^2/R^2$.

The common source operators are $O_\pm=[\Phi(0)\pm W(0,\pi)\Phi(\pi)]/\sqrt2$, $W=e^{i\pi\alpha}$. The interaction $\sum_a\lambda J_a^\dagger O_a+\mathrm{h.c.}$ gives $g_{\pm,n}=\lambda[1\pm e^{i\pi(n+\alpha)}]/\sqrt{4\pi R}$.

Every residue $\rho_n=g_ng_n^\dagger$ is positive semidefinite. Its diagonal obeys $|g_+|^2+|g_-|^2=\lambda^2/(\pi R)$ for each chosen Fourier momentum branch. Normalized relative probe weights are $(1\pm\cos\pi(n+\alpha))/2$, not probabilities or branching fractions.

At $\alpha=0$, these same probes select even/odd mode numbers. The signed-branch second difference is $m_{n+1}^2-2m_n^2+m_{n-1}^2=2/R^2$; this is not a statement about consecutively sorted degenerate levels.

The source $J_a$ and charged, gauge-covariant operator $O_a$ transform with the same charge at theta=0, so $J_a^\dagger O_a$ is invariant. With dimensionless lambda, $[J_a]=5/2$. This is a fixed-background source response, not a dynamical gauge-theory S-matrix. The paired nonlocal probes have not been derived from a local detector theory.

Here R is compactification radius, not the cited paper's dimension label R; alpha is a modeling hypothesis, not a proved identification with Ei, Ep, or En.

## Four predeclared engineering cards

One arbitrary common mass unit is used; R is inverse mass. M5=0.5 and common source normalization lambda=1 are fixed.

| Card | R | alpha | Lowest four mass-squared values (multiplicity retained) |
|---|---:|---:|---|
| circle_trivial | 1 | 0 | 0.25, 1.25, 1.25, 4.25 |
| circle_quarter | 1 | 0.25 | 0.3125, 0.8125, 1.8125, 3.3125 |
| circle_half | 1 | 0.5 | 0.5, 0.5, 2.5, 2.5 |
| circle_larger | 2 | 0.25 | 0.265625, 0.390625, 0.640625, 1.015625 |

These are background cards, not measured particles or a dynamical transition. Negative mode number is not negative energy.

## Independent numerical checks

Link-covariant finite-difference Laplacians are built and diagonalized at 32, 64, and 128 sites. Arbitrary site phases transform both links and the same physical probe Wilson line. Passive real and complex mode-basis changes transform both the kernel and the probe map. Under alpha -> alpha+1, the mode and cutoff window are relabelled n -> n-1.

An arbitrary eigenbasis within a degenerate eigenspace is not fixed by mass alone: unitary basis changes preserve its summed residue/projector. Particular states can still be distinguished by other operators or preparation. The residue sums are explicitly tested at alpha=0 and alpha=1/2.

**Important visibility qualification:** at $\alpha=0$, the columns for $n=\pm k$ are equal, so their sine standing-wave combination is dark to both antipodal probes. The coupling map on that two-dimensional eigenspace has rank one, although its pole is visible through the other combination. At $\alpha=1/2$, the tested degenerate pair has rank two. Therefore the branch-weight sum rule is not a proof that every degenerate state is visible. The invariant diagnostic is $\dim\ker G_E=\dim E-\operatorname{rank}G_E$.

| Card | Low-spectrum error: 32 sites | 64 sites | 128 sites | Last error ratio |
|---|---:|---:|---:|---:|
| circle_trivial | 0.805741 | 0.204563 | 0.0513382 | 3.98461 |
| circle_quarter | 0.623965 | 0.158118 | 0.0396637 | 3.98647 |
| circle_half | 0.474589 | 0.120055 | 0.0301026 | 3.98821 |
| circle_larger | 0.155991 | 0.0395295 | 0.00991592 | 3.98647 |

An error ratio approaching four is second-order discretization convergence, not a particle-physics uncertainty.

## Infinite image sum and rigorous truncation bound

For $p_E^2\ge0$, $G_{ab}=\sum_n g_{a,n}g_{b,n}^*/(p_E^2+m_n^2)$. Set $\beta=R\sqrt{M_5^2+p_E^2}>0$, $a=e^{-\pi\beta}$, $t=2\pi\alpha$, $D=1+a^4-2a^2\cos t$. Independent image sums give

$$S_0=\frac\pi\beta\frac{1-a^4}{D},\qquad S_\pi=\frac\pi\beta\frac{a(1-a^2)(1+e^{it})}{D},$$
$$G=\frac{\lambda^2R}{2\pi}\begin{pmatrix}S_0+\Re S_\pi&i\Im S_\pi\\-i\Im S_\pi&S_0-\Re S_\pi\end{pmatrix}.$$

Here $S_0=\sum_n[(n+\alpha)^2+\beta^2]^{-1}$ and $S_\pi=\sum_n e^{i\pi(n+\alpha)}[(n+\alpha)^2+\beta^2]^{-1}$.

The script compares this expression with finite mode sums at cutoffs 16,32,64,128,256 and p_E^2=0,0.7,4 for every card.

Since $|g_{a,n}g_{b,n}^*|\le\lambda^2/(\pi R)$ and $p_E^2+m_n^2\ge(|n|-|\alpha|)^2/R^2$, every entry for $N>|\alpha|$ obeys

$$|G_{ab}-G^{(N)}_{ab}|\le\frac{2\lambda^2R}{\pi}\sum_{n=N+1}^\infty\frac1{(n-|\alpha|)^2}\le\frac{2\lambda^2R}{\pi(N-|\alpha|)}.$$

This is a rigorous bound for the declared free response, not a statistical interval or interacting UV regulator.

| Card | p_E^2 | Largest entry error at cutoff 256 | Rigorous bound |
|---|---:|---:|---:|
| circle_trivial | 0 | 0.00124339 | 0.0024868 |
| circle_trivial | 0.7 | 0.00124339 | 0.0024868 |
| circle_trivial | 4 | 0.00124336 | 0.0024868 |
| circle_quarter | 0 | 0.00124268 | 0.00248923 |
| circle_quarter | 0.7 | 0.00124268 | 0.00248923 |
| circle_quarter | 4 | 0.00124266 | 0.00248923 |
| circle_half | 0 | 0.00124098 | 0.00249166 |
| circle_half | 0.7 | 0.00124097 | 0.00249166 |
| circle_half | 4 | 0.00124095 | 0.00249166 |
| circle_larger | 0 | 0.00248536 | 0.00497845 |
| circle_larger | 0.7 | 0.00248532 | 0.00497845 |
| circle_larger | 4 | 0.00248515 | 0.00497845 |

## Result and boundary of the claim

**331/331 numerical/code checks passed.** This does not validate a physical theory or establish an empirical mass spectrum.

Passive coordinate/field-basis changes leave poles and complete probe response invariant. Physically inequivalent holonomy or radius cards change poles. Fixed physical probes can weight modes differently without moving a pole.

Still absent: dynamical stabilization, local detector realization, transition amplitudes, spin/chirality/three families, a Standard Model spectrum, and a derivation of Route-F Yukawa matrices. None is implicitly promoted.

Reproduce: python3 route_g/code/verify_g1_circle_spectrum.py. JSON retains complex residues, responses, every check, seed, and source hash.

References: [Tong, section 8](https://arxiv.org/abs/0908.0333); [Yamada, section 2](https://doi.org/10.1093/ptep/ptab085).
