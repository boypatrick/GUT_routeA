# P54 gauge / ghost / fermion momentum audit

Date: 2026-09-13. Scope: the unchanged P54PQ-v2 action and original radial
background. This is a prescription and vertex audit plus an **evaluated
radial fermion subset**. It is not a completed gauged pole calculation.

## 1. What is already fixed, and what is not

The frozen action supplies 328 canonical real scalar coordinates, 45 real
antisymmetric gauge generators, the scalar Yukawa intertwiners and the
background `X=R r`, with

\[
r=(w,\sigma,v_s)=(1,0.1265,0.25)\omega,
\qquad R^T R=G=\operatorname{diag}(12/5,2,1).
\]

The original subtraction is dimensional regularization / MSbar, with a
background-field Landau label and common finite invariant mass
counterterms fixing the VEV. This label fixes the constant-background CW
calculation, but a momentum implementation must also write down its
gauge-fixing functional and distinguish background from quantum external
fields. Off-shell two-point functions in different such conventions must
not be spliced together without their finite field/parameter conversion.

The existing scalar result includes all scalar internal pairs, not merely
290 hard scalars. Its zero-momentum vector term does **not** include the
scalar-vector derivative bubble. The one-loop Landau diagram inventory is
independently confirmed by the three tensor structures in Braathen and
Goodsell, Eqs. (3.30)-(3.31): scalar-vector, vector seagull, and vector-vector
bubble. The first vanishes at zero momentum, so differentiating the vector
mass spectrum alone cannot recover it. Their two-loop mass treatment has
additional approximation boundaries and is not an already-complete P54
implementation. [Primary source](https://arxiv.org/html/1609.06977v2).

## 2. An explicit background-field completion of the scheme

The following equations are derived directly from the Euclidean kinetic
action; they specify a concrete completion to implement, not a claim that
the existing code has already implemented it. Choose

\[
D_\mu X=\partial_\mu X+g a_\mu^a T^a X,
\quad T^{aT}=-T^a,\quad X=\bar X+\eta,
\quad v_a=T^a\bar X,
\]
\[
F^a=\partial_\mu a_\mu^a+\xi g v_a^T\eta,
\qquad {\cal L}_{\rm gf}=\frac{F^aF^a}{2\xi},
\qquad \bar A_\mu=0.
\]

At a nonconstant scalar background, the quadratic bosonic operator is

\[
{\cal Q}_B=\begin{pmatrix}
Q_{SS}&W^T\\W&Q_{VV}
\end{pmatrix},
\]
\[
Q_{SS}=-\partial^2+V''(\bar X)+\xi g^2\sum_a v_av_a^T,
\]
\[
(Q_{VV})^{ab}_{\mu\nu}
=\left[-\partial^2\delta_{\mu\nu}
 +(1-\xi^{-1})\partial_\mu\partial_\nu\right]\delta^{ab}
 +M_{ab}^2(\bar X)\delta_{\mu\nu},
\quad M_{ab}^2=g^2v_a^Tv_b,
\]
\[
W^a_{\mu I}=-2g(T^a\partial_\mu\bar X)_I,
\qquad Q_c=-\partial^2\mathbf1+\xi M^2.
\]

The factor two in `W` is important. Before gauge fixing the mixed kinetic
term is `g a_mu [v dot partial_mu eta - (partial_mu v) dot eta]`.
Integrating the cross term of `F^2/(2 xi)` by parts cancels the first
piece and contributes another copy of the second. Therefore a constant
background diagonalizes the quadratic operator, but an external radial
field with nonzero momentum does not.

The full one-loop background functional has the unambiguous operator form

\[
\Gamma_1=\frac12\operatorname{Tr}\log Q_B
 -\operatorname{Tr}\log Q_c
 -\frac12\operatorname{Tr}\log D_F+\Gamma_{\rm CT}.
\]

The last factor `1/2` is for a Majorana/Nambu description without double
counting. For any bosonic radial source `r`, differentiation gives

\[
\Pi^B_{rs}=\frac12\operatorname{Tr}
 [Q_B^{-1}Q_{B,rs}-Q_B^{-1}Q_{B,r}Q_B^{-1}Q_{B,s}],
\]
\[
\Pi^c_{rs}=-\operatorname{Tr}
 [Q_c^{-1}Q_{c,rs}-Q_c^{-1}Q_{c,r}Q_c^{-1}Q_{c,s}].
\]

Products contain the appropriately routed external momentum, and local
counterterms are added only once. These formulas are an exact diagram
generator; they are **not** numerical results for the omitted integrals.

### Actual gauge tensors, requiring no flavor parameters

At the frozen point define `d_ar=T^a R_r`. Then

\[
B_{r,ab}=\partial_rM_{ab}^2
=g^2(d_{ar}^Tv_b+v_a^Td_{br}),
\]
\[
C_{rs,ab}=\partial_r\partial_sM_{ab}^2
=g^2(d_{ar}^Td_{bs}+d_{as}^Td_{br}),
\qquad W_{r,\mu I}^a(p)=-2igp_\mu(d_{ar})_I.
\]

The new verifier evaluates all `B_r`, `C_rs` and the derivative-vertex Gram
from the existing representation, saving them in its JSON. It verifies
`r_r B_r=2 M^2`, `r_r r_s C_rs=2 M^2`, independent finite differences, 33
broken generators, and the vanishing gauge vertex for `v_s`.

At `xi -> 0`, the vector propagator in a vector mass eigenbasis is

\[
D_{a,\mu\nu}(k)=
\frac{\delta_{\mu\nu}-k_\mu k_\nu/k^2}{k^2+m_a^2}.
\]

The remaining integrals are a vector seagull with `C_rs`, a two-vector
bubble with `B_r B_s`, and a mixed bubble from `W_r W_s`. In particular
the mixed contribution has the operator-integral structure

\[
\Pi^{SV}_{rs}(p)=-4g^2\sum_{a,k}
 [(T^aR_r)\cdot u_k][(T^aR_s)\cdot u_k]
 \int_k\frac{p_\mu D_{a,\mu\nu}(k)p_\nu}
 {(k+p)^2+m_k^2},
\]

with gauge generators rotated by the same vector mass diagonalization as
`D_a`. The formula denotes a dimensionally regularized integral followed
by the declared MS subtraction. Its zero at `p=0` does not justify
omitting it from a kinetic or pole calculation.

**Dimensional rational terms must be retained:** reducing a transverse
vector trace to the integer 3 before subtracting a UV pole discards the
finite product of its `-2 epsilon` component with that pole. This is why
the MS vector CW constant is `5/6`, rather than the scalar `3/2`, and why a
naive replacement of scalar bubble weights by three vector polarizations
is insufficient.

### Ghost qualification

There is no license to add a free ghost correction. In the explicit
functional above, the direct background-scalar ghost vertices are
`xi B_r` and `xi C_rs`; at strict Landau gauge the scalar-dependent ghost
operator is absent. At nonexceptional external momentum these direct
one-loop scalar ghost contributions have a vanishing Landau limit.
Nevertheless the finite-`xi` formulation, gauge-fixing scalar terms and
order of infrared limits must be retained for a Nielsen/BRST check.
Ghosts do not decouple from gauge-background correlators, general vertex
matching or box amplitudes. Consequently the old statement “massless
Landau ghosts have no hard CW term” is correct but is not a complete
momentum, matching or gauge-independence audit.

## 3. Exact radial fermion simplification

The already-verified P54 Yukawa tensors imply, on this radial family,

\[
M_F(r)=\sigma F_M,\qquad \partial_wM_F=\partial_{v_s}M_F=0,
\]

with only one massive spinor direction per family. Takagi diagonalizing
`F_M` leaves three nonnegative masses `M_i=sigma f_i`. The radial vertex
is simultaneously diagonal: `partial_sigma M_i=f_i`. Thus the **complete
one-loop fermion contribution on this radial slice** needs three Takagi
singular values, not a global flavor fit or a free matrix-valued
`eta_Y(p)` function. Full nonradial Higgs/PS amplitudes still require the
complex `h_raw`, `f_raw` matrices and their actual intertwiners.

Let `s=p_E^2 >= 0`, and use exactly the scalar verifier's convention

\[
L(m^2,m^2;s)=\int_0^1dx\,
 \log\frac{m^2+x(1-x)s}{\mu^2}.
\]

For a single four-component Majorana field, the quadratic expansion of
`-1/2 Tr log(slash partial+m+y delta sigma)` yields

\[
\Pi_F^{\rm raw}(s)
=\frac{y^2}{2}\int_k
 \operatorname{tr}\frac{(-i\not k+m)(-i\not k-i\not p+m)}
 {(k^2+m^2)((k+p)^2+m^2)}
\]
\[
=2y^2\int_k\frac{m^2-k\cdot(k+p)}
 {(k^2+m^2)((k+p)^2+m^2)}.
\]

With `A_E(m^2)=m^2[log(m^2/mu^2)-1]/(16 pi^2)` and the finite bubble
`I_B=-L/(16 pi^2)`, the numerator identity gives

\[
\Pi_F^{\rm raw}(s)=
-\frac{y^2}{16\pi^2}
 \left[(4m^2+s)L(m^2,m^2;s)
 +2m^2\left(\log\frac{m^2}{\mu^2}-1\right)\right].
\]

The same determinant gives

\[
V_F=-\frac{m^4}{32\pi^2}
 \left(\log\frac{m^2}{\mu^2}-\frac32\right),
\quad \frac{t_\sigma^F}{\sigma}
=-\frac{2y^2m^2}{16\pi^2}
 \left(\log\frac{m^2}{\mu^2}-1\right).
\]

The invariant fixed-VEV prescription adds the momentum-independent mass
counterterm `-t_sigma^F/sigma`, not a momentum-dependent subtraction.
Therefore

\[
\boxed{\Pi^{F,\rm FV}_{\sigma\sigma}(s)=
-\frac{1}{16\pi^2}\sum_{i=1}^3
 \frac{M_i^2}{\sigma^2}(4M_i^2+s)L(M_i^2,M_i^2;s).}
\]

All other radial entries vanish. At zero momentum this reproduces
the existing exact Majorana potential/tadpole result:

\[
\Pi^{F,\rm FV}_{\sigma\sigma}(0)=
-\frac{4}{16\pi^2\sigma^2}\sum_iM_i^4
 \log\frac{M_i^2}{\mu^2}
\leq\frac{6\mu^4}{e\,16\pi^2\sigma^2}.
\]

The inequality follows by maximizing `-x^2 log(x/mu^2)` at
`x=mu^2 exp(-1/2)` for each squared mass. It is **only** the old
zero-momentum bound and is not asserted at arbitrary momentum.

Differentiating the bubble at zero momentum gives the finite kinetic term

\[
\boxed{\Delta Z^F_{\sigma\sigma}=
-\frac{1}{16\pi^2}\sum_i\frac{M_i^2}{\sigma^2}
 \left(\log\frac{M_i^2}{\mu^2}+\frac23\right).}
\]

The canonical/generalized sigma component is half this coordinate value,
because `G_sigma,sigma=2`. A Dirac fermion with otherwise identical mass
and coupling contributes twice the Majorana result. Zero Yukawa masses
have the continuous zero limit, and no spurious `log(0)` is evaluated.

### Independent numerical verification

`code/verify_p54_radial_fermion_momentum.py` checks the closed bubble
formula against the separately reduced numerator integral

\[
\Pi_F^{\rm raw}(s)=\frac{2y^2}{16\pi^2}
 \int_0^1dx\,D(x)[1-3\log(D(x)/\mu^2)],
\quad D(x)=m^2+x(1-x)s,
\]

then adds the same counterterm. It also checks the kinetic derivative by
an independent one-sided finite difference, the old zero-momentum
function, its maximal bound, and zero Yukawa limits. Its five mass cards
are explicitly **synthetic regression inputs**, not fitted physical
families. No gauge momentum integral is silently inferred from these
tests.

## 4. Tadpole, background and physical-pole consistency

The fermionic radial expression must be combined with its own common
invariant mass shift, while the already-applied bosonic mass shift is not
applied a second time. At any new scalar benchmark all spectra, tadpoles
and their dependent counterterms must be recomputed. A fixed set of bare
parameters must be used under finite scheme changes. VEV/tadpole
renormalization can generate gauge dependences in parameter definitions;
a named tadpole prescription does not remove the need to transport those
definitions consistently. [Primary source: Dudenäs and Löschner](https://arxiv.org/abs/2010.15076).

Gauge covariance of a constant-background invariant potential gives
`H T X=T grad V`, not gauge independence of its off-shell eigenvalues.
For the full background effective action, the Nielsen identity has the
functional form

\[
\partial_\xi\Gamma+\int C_I\frac{\delta\Gamma}{\delta\bar X_I}=0.
\]

On a stationary background transported with `xi`, differentiation gives
the homogeneous two-point relation, schematically on the properly closed
physical block,

\[
\frac{d}{d\xi}\Gamma^{(2)}(p)
=\Gamma^{(2)}(p)\Lambda(p)
 +\Lambda^T(p)\Gamma^{(2)}(p).
\]

If the Nielsen insertion is regular at an isolated physical pole, this
relation preserves the location of that pole. It does not preserve every
off-shell matrix entry or wave-function coefficient. The published
all-orders analysis includes mixed propagators and explains why the
complex pole, rather than an arbitrary curvature definition, is the
appropriate gauge-independence target. Extending these identities to the
full P54 implementation remains an actual check to perform.
[Primary source: Gambino and Grassi](https://arxiv.org/abs/hep-ph/9907254).

The physical calculation must retain a closed SM/BRST sector (or its
justified Schur complement), perform the Euclidean-to-Minkowski
continuation with one sign convention, and check residues and perturbative
control. A Euclidean zero at positive `s` would correspond to a negative
Minkowski mass squared **only for that completed physical kernel**. The
existing partial scalar-plus-static-vector roots decide none of this.

## 5. Concrete next options and non-promotion gates

1. **Minimal radial completion:** use the already-evaluated `B,C,d` tensors
   to integrate the massive/massless vector and mixed loops, retaining
   dimensional rational terms. Reproduce the old zero-momentum vector
   CW derivative, then compare at small nonzero gauge parameter. Insert
   the exact three-mass fermion expression as a profile subset, with
   explicitly constrained Yukawa perturbativity. This closes a useful
   radial subproblem before attempting all 328 external fields.
2. **Gauge-invariant probe cross-check:** couple sources to the same
   invariant radial composites (e.g. `Tr Phi^2`, `Sigma† Sigma`, `S†S`)
   and compute their connected response with the required vertex/contact
   terms. Agreement of physical singularities offers an independent
   check; the composite correlator is not simply a rotated elementary
   Hessian and must not be substituted without those terms.
3. **Continue the scalar repair independently:** the nonuniform
   invariant analysis can reject or identify controlled candidates
   without calling a partial momentum kernel a physical stability proof.

No option promotes complete Wilson/box diagrams, lower gauge/ghost/Yukawa
matching, the finite seesaw including `C_HN`, or global flavor fit.
The new result removes unnecessary family-matrix complexity from one
well-defined radial fermion subset; it does not cure the large original
scalar insertion or finish the theory.
