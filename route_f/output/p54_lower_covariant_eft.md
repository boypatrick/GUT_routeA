# P54 lower covariant tree EFT

Constructed a local PS-covariant tree scalar EFT through two covariant derivatives, including 55 scalar-valley responses and algebraic elimination of the 24 upper vectors. This is **not** completed lower one-loop finite gauge matching.

## Exact local construction

The raw 248-dimensional PS tensor chart is $X(y)=x_0+B_Ay+B_Hz(y)$ with $B_H^T\nabla V(X)=0$. Its positive heavy block defines an analytic germ by the implicit-function theorem. The axion is held fixed as a gauge-neutral spectator.

$$C=B_H^THB_H,\quad D=-C^{-1}B_H^THB_A,\quad T=B_A+B_HD,\quad X_{,ab}=-B_HC^{-1}B_H^TV^{(3)}[T_a,T_b].$$

With $R=-iT$, the card's $D=\partial-igAT$ is $D=\partial+gAR$. For $O_A=R_AX$, $I=O^TO$ and $j=gO^TD_{PS}X$,

$$A_U^*=-g^{-1}I^{-1}O^TD_{PS}X,\qquad \Delta\mathcal L_2=-\frac12j^T(g^2I)^{-1}j,$$

$$\mathcal L_2=\frac12(Dy)^T\underbrace{T^T(1-OI^{-1}O^T)T}_{g(y)}(Dy)-V(X(y)).$$

The 24 directions are the symmetric coset $SO(10)/PS$, not a subgroup: $[PS,U]\subset U$, $[U,U]\subset PS\ne0$. Thus $I^{-1}O^TT$ is a PS-covariant response one-form, not a principal connection for a fictitious 24-dimensional group.

## Actual numerical tensors

- Heavy scalar gap: `0.0994502574109`; quotient metric range: `[0.984249788817, 1.0180072705]`.
- Upper-vector current-current metric subtraction has rank `12`, trace `0.189002534197`.
- Higgs-direction response derivative norm: `0.00141167025192`; metric-derivative norm: `0.00716183316691`.
- Previous horizontal Schur metric is reproduced with error `1.042e-15`. All nine lower gauge masses and Goldstone normalizations agree.

## Beyond the Hessian

For the exported corrected neutral direction $q$, $J_H=B_H^TV^{(3)}[q,q]$,

$$\lambda_{\rm upper}=\frac16V^{(4)}[q^4]-\frac12J_H^TC^{-1}J_H,\qquad z_{,hh}=-C^{-1}J_H.$$

- Direct tree quartic `0.42782363928` becomes `0.407929900088` after all 55 upper scalars relax.
- The induced radial derivative operator is `0.0247875490729 h^2 (D h)^2` in omega units. The upper-vector current on this color-singlet ray vanishes by selection rules.
- Canonically normalized retained 126 triplet masses: `[0.26883282614585813, 0.2688328261458583]`.
- Its actual effective hh source: `[0.002527431093495623, 1.7899039645361203e-16]`; reconstruction agrees with the full mixed 54/126 response to `6.960e-17`.
- These are frozen tree tensors evaluated on a bosonic-improved ray, not loop-complete Higgs/type-II Wilson coefficients. The ray coordinate has unit metric only at the origin; all 248 active parents remain. Thus this is neither an SM quartic nor an all-orders geodesic-coordinate coupling. No fermion Clebsch is inferred.

## Scope and remaining matching

Scalar elimination already produces $-p^4W^TC^{-3}W$. Substitution of $A_U^*$ into $F_U=D_{PS}A_U^*$ and $F_{PS}=F_{PS}^{(0)}+g[A_U^*,A_U^*]_{PS}$ generates further covariant order-p4 operators. They cannot be silently discarded inside a claimed full one-loop finite threshold.

The retained PS transformation is affine about the broken vacuum. The verifier checks projector covariance, SM Ward identities, a finite broken-PS rotation, nonlinear V3 q-orbit identities, all nine lower Goldstone normalizations, and two independent jet steps. Positivity proves only a local valley germ, not global continuation or uniformly small retained momenta.

The largest retained O(p2) generalized mass squared is `3.71306` times the smallest integrated scalar mass squared; `224` retained real modes lie above that smallest gap. Thus these generalized masses must not be promoted to uniformly accurate poles or finite-threshold logarithms.

That global gap test is only sufficient, not a new no-go: actual source support leaves `25` linearly coupled modes and `223` decoupled modes. The `12` modes at m2 = 0.268823846583 couple weakly to C = 0.261600225 (ratio 1.0276132); for this near-crossing the Taylor series fails while the exact quadratic Schur resolvent remains available. The public `scalar_nonlocal_kernel_api` computes that resolvent rather than guessing a mass correction.

The actual source-selected canonical 2x2 block is invariant under the full Hessian to `5.448e-12` and has exact tree poles `[0.26135949490123755, 0.269078530098762]`. Its retained O(p2) mass error is only `0.0946503%`: a resolvent repair, not a physical model obstruction.

The action is frozen at tree level. A fixed-VEV one-loop CT affects tree EFT coefficients at one-loop order and must accompany a future loop EFT; inserting that shift inside an already one-loop gauge threshold is a higher-order operation. No upper saddle displacement or physical branch change is made here.

The public tensor functions return full local valley, current, metric, directional-vertex tensors and the exact quadratic scalar resolvent. Calling `run(include_tensors=True)` also exports the exact common bases, all generator representations, the potential callable and the same tested full tensor payload. The next matcher must construct the lower covariant gauge-scalar-ghost fluctuation operator and perform a hard-minus-EFT subtraction with the upper Wilson action in the same scheme.

Validation: **24/24**; Hessian cache: `{'cache_hits': 5, 'evaluated': 0}`.

References: [covariant functional matching](https://arxiv.org/abs/1604.01019), [covariant derivative expansion](https://arxiv.org/abs/1412.1837).
