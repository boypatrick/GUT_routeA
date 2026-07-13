# Route A--E Theory Audit and Route-F Synthesis

Date: 2026-07-13

## Executive verdict

The project has a valuable conditional mathematical core, but it is not yet a
closed physical theory.  The robust part is:

\[
\operatorname{Spin}(10):16
\quad+\quad
H^0(\mathbb P^1,\mathcal O(2))
\quad+\quad
\text{the unique }SL(2)\text{-invariant bilinear direction}.
\]

These statements organize one family, a three-dimensional carrier, and one
Majorana tensor direction.  They do **not** yet derive a unique action, the
family geometry, the coefficient \(\zeta\), a physical flavor fit, proton
lifetimes, or a UV compactification.

Route F should therefore do four things in order:

1. repair the theorem-level logical gaps and freeze a reproducible evidence
   set;
2. select one explicit action and discard incompatible mixtures of SUSY and
   non-SUSY assumptions;
3. connect the actual vacuum spectrum to two-loop running, flavor, amplitudes,
   and proton decay;
4. treat string geometry and the hidden messenger as optional completions with
   hard promotion gates.

## Route A: representation/index/contact skeleton

### What is established

The standard Pati--Salam restriction of the half-spinor is correct:

\[
16\to(4,2,1)\oplus(\bar4,1,2),\qquad
Y=T_{3R}+\frac{B-L}{2},
\]

and it yields \(Q,L,u^c,d^c,\nu^c,e^c\).  The hypercharge normalization
also checks directly:

\[
\operatorname{Tr}Y^2
=6\left(\frac16\right)^2+2\left(\frac12\right)^2
+3\left(\frac23\right)^2+3\left(\frac13\right)^2+1
=\frac{10}{3},
\]

\[
\operatorname{Tr}T_{3L}^2=2,
\qquad k_Y=\frac{\operatorname{Tr}Y^2}{\operatorname{Tr}T_{3L}^2}
=\frac53.
\]

The family-index calculation is also correct after its assumptions are made:

\[
\chi(\mathcal O(2))=h^0-h^1=3,
\qquad
h^1(\mathcal O(2))=h^0(\mathcal O(-4))=0,
\]

so \(h^0=3\).  The paper itself correctly lists \(C=\mathbb P^1\), the
degree-two divisor, the carrier choice, and absence of other chiral sectors as
assumptions (`paper/gut_framework.tex:403-466`).

Finally,

\[
\operatorname{Sym}^2(\operatorname{Sym}^2\mathbb C^2)
\simeq \operatorname{Sym}^4\mathbb C^2\oplus\mathbf1
\]

does give one invariant symmetric direction.  This fixes a direction, not its
coefficient or its physical origin (`paper/gut_framework.tex:595-648`).

### Missing pieces

P0:

- A1--A7 are labels, not one normalized action.  SUSY status, Higgs content,
  breaking potential, source/projector, doublet--triplet splitting, and scales
  must be specified simultaneously.
- If the project claims a first-principles geometric origin, the index is not
  yet a compactification theorem: there is no flux quantization, anomaly
  inflow, KK spectrum, or proof of exactly three complete chiral `16`s with no
  other charged zero modes.  This is optional, rather than P0, for a strictly
  conditional four-dimensional EFT whose carrier is declared as an axiom.
- The theorem boundary invokes companion phenomenology that the current Route-A
  manifest explicitly calls deferred (`paper/gut_framework.tex:1248-1332`).

P1:

- The phrase "standard minimal" should always be qualified by the candidate
  class and minimal-dimension assumption.  Conjugating the subgroup by a gauge
  transformation is normally gauge redundancy, not a distinct physical face;
  distinct faces require vacuum or global data.
- Any symmetric \(3\times3\) tensor lies in the \(5\oplus1\) decomposition.
  A machine-zero `Veronese + contact` residual is therefore representation
  bookkeeping, not a phenomenological prediction.

## Route B: hidden messenger

### Correct algebra

For

\[
\frac{W}{M_*}=\frac12N^TM_VN+\lambda X^TN
-\frac12X^TK_{\rm tr}^{-1}X,
\]

the \(X\) equation gives \(X=\lambda K_{\rm tr}N\), hence

\[
\frac{W_{\rm eff}}{M_*}
=\frac12N^T(M_V+\lambda^2K_{\rm tr})N.
\]

This Schur-complement statement is sound.  Setting
\(\lambda^2=\zeta\) realizes an already supplied coefficient; it does not
predict it (`paper/gut_framework.tex:725-810`).

### P0 selection-rule counterexample

The paper assigns

\[
R(X)=R(16)=1,\qquad R(H_u)=0,\qquad R(W)=2.
\]

After \(B-L\) breaking, \(X\) is an SM singlet.  Therefore

\[
W\supset y_X X_iL_jH_u,
\qquad R(XLH_u)=1+1+0=2,
\]

is both SM-gauge invariant and allowed by the stated `U(1)_R`.  The proof in
`paper/gut_framework.tex:988-1040` only forbids adding extra \(X\)'s to an
existing `16 16 H` operator; it does not forbid replacing the sterile leg by
\(X\).  Consequently the "tree-level matching silence" assumption in
`paper/gut_framework.tex:1053-1094` is not derived.  Route F must add a genuine
gauge/discrete selection rule and enumerate all operators through a declared
dimension.

### P0 canonical-normalization correction

If \(X\) has canonical Kahler potential, substituting the low-energy solution
also yields

\[
K_{\rm eff}=N^\dagger N+X^\dagger X
=N^\dagger\!\left(I+|\lambda|^2K_{\rm tr}^\dagger K_{\rm tr}\right)N.
\]

For the explicit matrix used here, \(K_{\rm tr}\) is real symmetric, so
\(K_{\rm tr}^\dagger K_{\rm tr}=K_{\rm tr}^2=I/3\).  Therefore

\[
\delta Z_{N,\text{tree}}=\frac{|\lambda|^2}{3}
=\frac{|\zeta|}{3}=0.0434773011.
\]

The loop estimate printed in the paper is

\[
\delta Z_{N,\text{loop}}\sim\frac{|\lambda|^2}{16\pi^2}
=8.2597\times10^{-4},
\]

about 52.6 times smaller.  If \(N\) and \(X\) are at comparable scales, the
proper calculation is the full \(6\times6\) Takagi/block-propagator problem

\[
\mathcal M_{NX}=M_*
\begin{pmatrix}M_V&\lambda I\\ \lambda I&-K_{\rm tr}^{-1}\end{pmatrix},
\]

not a superpotential-only elimination followed by a loop-only wavefunction
estimate.  The zero-momentum Schur complement may survive, but canonical
normalization and the Weinberg operator must be verified together.

There is a second issue with the quoted loop estimate: the displayed Route-B
superpotential is purely quadratic.  A Gaussian bilinear mass mixing by itself
does not generate a one-loop anomalous dimension proportional to
\(|\lambda|^2/(16\pi^2)\).  Such a loop requires additional propagating
trilinear or gauge interactions, which must be written before this estimate can
be used.  In the displayed model the tree-level kinetic mixing is the first
effect to audit.

For the canonical low-energy interpretation,
\(Z_N=1+|\zeta|/3=1.0434773011\).  Without a compensating redefinition of the
input parameters, canonical normalization changes the contact magnitude from
\(0.1304319033\) to \(|\zeta|/Z_N=0.1249973556\), a \(4.17\%\) shift.  The
shift is a basis-parameter change, not by itself an observable correction:
when the Dirac coupling and Majorana mass are transformed consistently, the
physical combination \(Y M^{-1}Y^T\) is invariant.  The physical gate is
therefore the complete Weinberg coefficient, not \(\zeta\) in isolation.

### Other missing pieces

- \(P_{\nu^c}\) is an input.  A \(\overline{126}_H\) source can itself carry
  an arbitrary symmetric family tensor, so Route B must show why the source
  cannot already generate the contact it is meant to explain.
- \(\lambda=\sqrt\zeta\) and \(S=-\log|\zeta|\) are reparameterizations unless
  the hidden parameters are frozen before viewing flavor data.
- The physical phase must be stated through rephasing-invariant CP observables.
- `Delta b_vis = 0` is only the singlet's direct one-loop gauge-beta
  contribution.  It does not establish silence of the `B-L` source thresholds,
  two-loop Yukawa effects, or heavy-neutrino matching.

## Route C: amplitude bootstrap

### What is useful

Route C has good charge/anomaly bookkeeping, a transparent debt ledger, and a
correct staged neutral-vector matrix.  For example,

\[
M_0^2=\frac{v_Y^2}{4}
\begin{pmatrix}g_R^2&-g_Rg_{BL}\\-g_Rg_{BL}&g_{BL}^2\end{pmatrix}
\]

has eigenvalues \(0\) and
\(v_Y^2(g_R^2+g_{BL}^2)/4\).  With
\((g_R,g_{BL},v_Y)=(0.73,0.51,2.4)\), they are numerically
\((-1.24\times10^{-16},1.14192)\), consistent with the analytic result.

### What the present passes do not prove

- P3 defines a residue matrix as \(R=gg^\dagger\) and then checks that it is
  positive semidefinite.  This is true by construction, not a residue extracted
  from a helicity amplitude.
- P9 checks counts and subset relations but not
  \([T_A,T_B]=if_{AB}{}^CT_C\), Hermiticity, or trace normalization.
- P10 assembles positive blocks to obtain 33 positive vector masses; it does
  not derive them from a chosen Higgs potential.
- P16 forms Cartesian products of reverse maps selected by charges.  Its 144
  "canonical BNV" rows are not 144 independent physical Wilson operators
  until root, color, weak, propagator, Fierz, and flavor contractions are done.
- The scalar/source and vector branches mix dimension-five SUSY and
  dimension-six non-SUSY language without first choosing an action.

For an uncancelled longitudinal-vector amplitude,

\[
a_0\sim\frac{g^2E^2}{16\pi M^2}.
\]

The tree-unitarity condition \(|\Re a_0|\le1/2\) gives

\[
\frac{E}{M}\lesssim\frac{\sqrt{8\pi}}{g}.
\]

For \(g=0.7\), the scale is only \(E/M\simeq7.16\).  Charge filtering and a
formal PSD matrix do not test the cancellations needed to avoid this growth.
Route F must extract residues from complete crossing-symmetric amplitudes and
check Ward/Goldstone identities pointwise.

## Route D: string-constrained lift

D1 correctly checks the local dimension identity

\[
78=45+1+16+\overline{16}
\]

and the local curve cohomology

\[
K_{\mathbb P^1}^{1/2}\otimes\mathcal O(3)=\mathcal O(2),
\qquad(h^0,h^1)=(3,0).
\]

D2 correctly reproduces the benchmark arithmetic

\[
|\zeta|=0.1304319033,
\quad\arg\zeta=0.6000380203,
\quad-\log|\zeta|=2.0369040.
\]

However, D1's script does not construct a resolved global geometry, quantized
\(G_4\) flux, D3 tadpole, massless hypercharge, or an exotic-free spectrum.
D2's script does not construct an instanton divisor or calculate universal and
charged zero modes, Freed--Witten flux, or a Pfaffian.  The Majorana-only
zero-mode saturation is assumed.  D1/D2 are therefore checked local
placements/interpretations, not a global F-theory/E3 completion; their
promotion to a global completion has not been performed.

## Route E: first-principles reconstruction

### Strong part

Route E improves the explicitness of the assumptions, performs a useful finite
representation scan, derives the genus ladder, identifies the Killing form
once \(\mathfrak{sl}_2\) is chosen, and proves an honest underdetermination
statement for \(\zeta\) inside a specified \(\zeta\)-blind witness family.

### P0 counterexample to the three-family selection as written

For a compact curve,

\[
h^0(\Sigma_g,T_{\Sigma_g})=
\begin{cases}3,&g=0,\\1,&g=1,\\0,&g\ge2.\end{cases}
\]

The ladder is correct.  The next inference is not.  H3 asks for a
"canonical invariant symmetric bilinear form", but the proof uses the Cartan
criterion to reject the one-dimensional abelian algebra because its Killing
form is zero (`route_E/tex/route_e_first_principles.tex:228-234,590-604` and
`route_E/tex/three_families_killing_contact.tex:238-250`).

Let \(\mathfrak g=\mathbb C e\), \([e,e]=0\), and define \(B(e,e)=1\).
Then

\[
B([x,y],z)+B(y,[x,z])=0
\]

identically and \(\det B=1\).  Thus a non-degenerate symmetric ad-invariant
form exists.  The Cartan criterion distinguishes the **Killing form**; it does
not say all invariant forms on a non-semisimple algebra are degenerate.
This counterexample invalidates that step of the printed proof; it does not by
itself make \(B=[1]\) a geometrically canonical choice, because "canonical"
has not yet been defined.

Two valid repairs are available:

1. strengthen H3 to require the bracket-generated Killing form explicitly; or
2. define "canonical" functorially as invariance under every automorphism of
   the Lie algebra with no scale/polarization datum, and prove that this
   naturality kills the abelian branch.

For the second repair, the relevant group must be stated precisely.  The full
Lie-algebra automorphism group of \(\mathbb C e\) is \(\mathbb C^*\); invariance
under every rescaling forces \(B=0\).  Requiring invariance only under the
finite automorphisms of a generic elliptic curve (often just a small subgroup)
would not by itself remove the counterexample.

Until one repair is written and independently checked, the theorem proves
\(N_{\rm fam}\le3\), not \(N_{\rm fam}=3\).

### Other theorem-domain issues

- The "two centers" result is only a degree-two zero divisor in general.  Two
  distinct zeros require a regular semisimple vector field; a nilpotent field
  can have a double zero.
- The forced \(\nu^c\) result is a minimal-dimension statement.  The theorem
  and Letter must keep the minimality assumption visible and must not imply an
  enumeration of every higher-dimensional irrep and every SM embedding.
- The boundary theorem shows that a \(\zeta\)-blind principle set valid across
  the witness family cannot select one \(\zeta\).  It does not forbid a new
  physical axiom from selecting a subfamily; the paper itself acknowledges
  this.

### P0 evidence and status integrity

The paper says the completed DYN-0--DYN-5, DYN-7--DYN-9b and Route-E
string-pricing-card `RE-SC3/4/5` lanes were executed and collected
(`route_E/tex/route_e_first_principles.tex:926-1037,1277-1285`).  The user has
now restored all 19 DYN source scripts under the non-canonical
`route_E/code_dyn/` location.  This corrects the earlier source-existence
finding.  It does not yet restore the cited `code/...` paths, any generated
`dyn*.json/md` ledger, clean-clone provenance, or the Route-E string cards
`RE-SC3/4/5` historically named as D3--D5 files under `route_d/`.

An isolated compatibility replay separates source quality from deployment:
16 scripts execute with `219/219` of their own checks green; DYN-9b-2 fails on
the missing `RE-SC3` ledger, DYN-9b-3 then fails on the missing 9b-2 ledger,
and DYN-8 fails on `RE-SC3` before reaching its remaining dependencies.  The detailed
method, numerical results, and scientific reclassification are in
`route_f/CODE_DYN_REAUDIT.md`.  The evidence audit covers a targeted 29-item
dynamics subset plus 4 core artifacts, not a complete automatic extraction of
every path in the roadmap.

The recovered calculations support a conditional CMSGUT benchmark spectrum
and a no-go for that specified SUSY slice.  They do not close the active
non-SUSY theory.  In particular:

- DYN-5 treats the quadratic mass mixing `lambda X N` as a dimensionless
  Yukawa interaction, omits the larger tree-level kinetic correction and the
  allowed `X L H_u` operator, and therefore cannot establish one-loop silence;
- DYN-4a's displayed 68% intervals are fixed-coordinate slices whose grids are
  under-resolved by factors about `4.96e3` and `1.66e3`; DYN-4b remains
  conditional on archival Dirac kernels and nuisance priors; DYN-4c's exact
  absorption is the inverse-seesaw identity, not a global Yukawa fit;
- DYN-7 and DYN-9b-3 use an unflavored leptogenesis calculation at
  `M_1 ~ 2.4e10 GeV`, where a flavored treatment is required;
- DYN-9b-1b/c/d are tree-level random scans with disclosed omitted invariants
  and stationarity work.  The `8/8800` LR hits are an existence hint in the
  sampled measure, not a global vacuum theorem;
- DYN-8's K5r check accepts rare 210-only LR minima while its emitted branch
  map still publishes the superseded claim that a `45_H` route is required.

At the same time:

- the branch map says the non-SUSY vacuum/flavor work ran;
- the status ledger still lists that same work as open
  (`route_E/tex/route_e_first_principles.tex:1224-1237`);
- the conclusion again says the deferred audits remain open;
- the existing core JSON flags `companion_phenomenology_audits_replayed=false`;
  that means the core 44-check card cannot support the later dynamics claims,
  not that it proves those dynamics were never run elsewhere;
- the text freezes experimental status at 2026-07-05 even though later sections
  incorporate 2026-07-06 corrections.

The existing Route-E core and recovered DYN source scripts are real and should
be kept.  Source recovery is not evidence closure: no DYN claim should be
cited beyond its reclassified conditional domain until the exact paths,
inputs, outputs, commands, digests, scientific gates, and clean-clone replay
are repaired.

## Unifying proposal for Route F

### Three candidate action branches

1. **F-210 continuity branch (recommended first audit):** non-SUSY
   `Spin(10)` with `210_H` for flexible Pati--Salam/left-right breaking and
   `10_H + 120_H + overline{126}_H` for flavor.  It best reuses the current
   Route-E geometry and the modern minimal-Yukawa phenomenology, but it needs a
   complete non-SUSY scalar potential and two-loop thresholds.
2. **F-54 conservative branch:** `54_H + 126_H + 10_H,C`, close to a published
   non-SUSY two-step model.  It is a strong cross-check and may be simpler, but
   its parity/intermediate-scale constraints must be tested against proton
   decay and thresholds.
3. **F-string branch:** keep D1/D2 as optional interpretations until a global
   compactification and zero-mode audit pass.  It must not carry any core
   phenomenology before then.

Do not choose the minimal `45_H + 126_H + 10_H,C` model as the default without
addressing the 2023 analysis that finds a low-GUT-scale problem and difficulty
obtaining a perturbative SM-like Higgs doublet.

Because the recommended branches are non-supersymmetric, Route B cannot be
inserted literally as a superpotential/Kahler/`U(1)_R` sector.  Route F must
either translate it into ordinary Weyl-fermion mass and kinetic matrices with
a genuine gauge or discrete symmetry, or retain the superspace mechanism only
as a SUSY comparison branch.  The Schur-complement algebra survives this
translation; the holomorphy and `R`-symmetry claims do not automatically do so.

### One identifiability map

Every quantitative claim should be organized as an observable map

\[
F:\Theta/\mathcal G_{\rm basis}\longrightarrow\mathcal O,
\]

where \(\Theta\) is the action parameter space and
\(\mathcal G_{\rm basis}\) contains field redefinitions.  The Jacobian rank

\[
r(\theta)=\operatorname{rank}\frac{\partial F}{\partial\theta}
\]

separates identifiable combinations from flat directions.  This converts the
Route-E boundary idea into a practical test: a parameter is predicted only if
it is fixed on the quotient after the fit, not merely reconstructed by choosing
hidden couplings.

### Why thresholds require an uncertainty calculation

For gauge-mediated dimension-six proton decay,

\[
C_6\sim\frac{g_U^2}{M_X^2},\qquad
\Gamma_p\propto\frac{\alpha_U^2}{M_X^4},\qquad
\tau_p\propto\frac{M_X^4}{\alpha_U^2}.
\]

A factor two shift in \(M_X\) changes \(\tau_p\) by \(16\); one decade changes
it by \(10^4\).  A one-loop central scale without correlated threshold and
lattice uncertainties cannot support a sharp alive/dead verdict.  Route F
therefore requires two-loop running, one-loop matching, the actual spectrum,
and an uncertainty budget before applying experimental bounds.

## Final recommendation

The theory is not ready to be called complete.  The mathematical skeleton can
be retained, and several no-go/boundary statements are valuable.  The next
decisive work is not another texture scan: it is (i) repair H3, (ii) repair the
Route-B operator and Kähler matching, (iii) recover the missing evidence, and
(iv) freeze one non-SUSY action whose vacuum, running, flavor, amplitudes, and
proton decay are computed in a single chain.  The concrete execution order and
kill criteria are in `ROADMAP.md`.
