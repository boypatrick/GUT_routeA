# P54 actual scalar-cubic kinetic matching subset

Date: 2026-09-06.

The calculation reuses the frozen tree action, 290 massive scalar propagators and all 38 soft directions. It includes mixed hard-soft loops and subtracts purely soft loops. It is a one-step broken-background scalar contribution, not the completed staged lower threshold. Massless Goldstone propagators inherit background-field Landau xi=0; this scalar subset is not gauge independent. With exactly zero soft masses, the mixed bubble's soft-region momentum expansion is scaleless in dimensional regularization. Finite soft masses or a staged EFT require a new subtraction, not this shortcut.

## Derivation

For $V\supset\tfrac12\eta_i(M^2_{ij}+T_{aij}q_a)\eta_j$, expand $\Gamma_1=\tfrac12\mathrm{Tr}\log(-D^2+M^2+T_aq_a+\cdots)$. The bubble is $-\tfrac14\mathrm{Tr}(GTqGTq)$; the inverse two-point function therefore has the following finite Euclidean momentum difference:

$$\Pi_{ab}(p_E^2)-\Pi_{ab}(0)=\frac1{32\pi^2}\sum_{ij}^{\rm hard}T_{aij}T_{bji}\int_0^1\log\!\left[1+\frac{t(1-t)p_E^2}{(1-t)m_i^2+tm_j^2}\right]dt.$$

Thus $K_{ab}=\partial_{p_E^2}\Pi_{ab}(0)$ is a positive weighted Gram matrix. The ordered pair sum contains at least one hard index. Quartic scalar tadpoles have no momentum derivative. A one-loop mass counterterm inside this already one-loop bubble would be a selected two-loop contribution and is not inserted.

$$I(x,y)=\int_0^1\frac{t(1-t)dt}{(1-t)x+ty}=\frac{x^2-y^2-2xy\log(x/y)}{2(x-y)^3},\quad I(x,x)=\frac1{6x},\quad I(x,0)=\frac1{2x}.$$

The stable near-degenerate series is $I=\frac1m\sum_{n\ge0}d^{2n}/[2(2n+1)(2n+3)]$, with $m=(x+y)/2$ and $d=(x-y)/(x+y)$. A single real-heavy field with interaction $ghS^2/2$ gives $K_h=g^2/(192\pi^2M^2)$.

## Actual tensors and covariant completion

Cubic jets are exact central differences of the quartic action's Hessian. Hypercharge transports the real jets to imaginary jets; an independent imaginary cache check verifies this. SU(2) generates the charged partners, and all 16 real external directions are assembled. The resulting full matrix obeys the four electroweak generator Ward identities.

At quadratic order in the doublets, the covariant kinetic operator is $(D_\mu H_a)^\dagger K_{ab}(D^\mu H_b)$. Its three- and four-point gauge couplings are fixed by this operator. This completion does not replace an independent calculation of gauge-loop diagrams or the Nielsen identity, and higher-field derivative operators remain separate.

```text
[[ 2.0683806247e-03 -3.3315753981e-04 -2.0783805116e-06 -9.0740536523e-05]
 [-3.3315753981e-04  2.1359111983e-03 -1.7663312370e-05  3.7882581286e-06]
 [-2.0783805116e-06 -1.7663312370e-05  7.1006598164e-02 -4.5734524453e-04]
 [-9.0740536523e-05  3.7882581286e-06 -4.5734524453e-04  8.5994596130e-02]]
```

The imaginary matrix norm is 3.45469e-18; the four eigenvalues are [0.0017672291465608018, 0.0024369597547190575, 0.07099266045263997, 0.08600863676333875].

On the actual bosonic-improved light vector, delta Z = 0.00243676715822; the scalar-leg-only normalisation factor is 0.998783838597.

Canonical matching uses $Z=1+K$, $D_c=Z^{-1/2}DZ^{-1/2}$ and $c_c=Z^{1/2}c/\sqrt{c^\dagger Zc}$. It follows that $Y_c(c_c)=Y(c)/\sqrt{c^\dagger Zc}$, with conjugation for down/e. The exact square roots test the algebra; their higher powers are not a two-loop prediction.

## Validation and scope

24/24 checks pass; 13 old Hessian cache hits, zero new Hessians.

| Check | Residual | Pass |
|---|---:|:---:|
| exactly 290 positive and 38 soft scalar directions | 0.000e+00 | yes |
| all soft masses vanish within the tree tolerance | 9.421e-16 | yes |
| sixteen canonical real components of four weak doublets | 3.123e-15 | yes |
| hypercharge transports real to imaginary cubic jets | 3.773e-16 | yes |
| independent cubic derivative steps agree | 0.000e+00 | yes |
| charged rotation fixes the SM-invariant vacuum | 5.926e-16 | yes |
| complex kinetic matrix is Hermitian | 5.431e-18 | yes |
| eight-real kinetic matrix has the full complex block structure | 1.804e-16 | yes |
| charged and neutral kinetic tensors coincide | 7.491e-16 | yes |
| no neutral-charged kinetic mixing | 5.795e-17 | yes |
| hard-heavy and mixed hard-soft pieces exhaust the result | 4.945e-16 | yes |
| scalar-cubic contribution is a positive weighted Gram matrix | 0.000e+00 | yes |
| all four electroweak Ward identities hold on the actual real doublet plane | 8.655e-16 | yes |
| independent Feynman-parameter integration including degenerate and massless limits | 1.166e-15 | yes |
| single real-heavy toy has g^2/(192 pi^2 M^2) | 0.000e+00 | yes |
| purely soft IR bubble is rejected by the hard-kernel API | 0.000e+00 | yes |
| nonzero Euclidean momentum difference tends to the computed slope | 1.606e-05 | yes |
| halving momentum squared improves the derivative approximation | 0.000e+00 | yes |
| u: scalar leg and light-vector normalisation cancel copy rotations | 4.291e-16 | yes |
| d: scalar leg and light-vector normalisation cancel copy rotations | 1.077e-16 | yes |
| e: scalar leg and light-vector normalisation cancel copy rotations | 1.248e-16 | yes |
| nu: scalar leg and light-vector normalisation cancel copy rotations | 4.437e-16 | yes |
| canonical congruence preserves the tuned zero mode | 1.049e-16 | yes |
| canonically normalised light vector has unit norm | 0.000e+00 | yes |

The nonzero-momentum Feynman-parameter integral is independently evaluated at two momenta. This closes a scalar momentum-dependence subset, not Higgs pole matching. In particular the heavy-doublet gap used in the original curvature analysis is not reinterpreted as a pole gap.
