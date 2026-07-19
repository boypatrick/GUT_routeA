# Route E Roadmap: Dependency Closure

Goal (2026-07-04): remove or fully explain every dependency of
`tex/route_e_first_principles.tex` on Route A/B/D material and on classical
citations; keep the proved-vs-retained boundary explicit.  Items are ordered;
each is `open` / `done` / `permanently-open (by theorem)` / `out-of-scope`.

Categories: **[remove]** = prove in-repo, dependency eliminated;
**[explain]** = retained input, justified and tagged, not eliminable by
principles of the P0-P3 kind; **[impossible]** = provably not closable;
**[scope]** = inherited Route-A obligations, untouched by Route E.

## AP-E12 endpoint density/regularity update (2026-07-20; current authority)

No Route-E dependency or portal is promoted.  AP-E12 proves that degree is
continuous in the complete-minor graph norm, so any target-valued endpoint
smooth recovery is automatically eventually fixed-degree.  It also proves
that the AP-E11 continuum density has uniform rank-one ellipticity with
constant one and an exact strong-quasiconvexity identity obtained from the
null-Lagrangian means of all first, second, and third minors.

The prioritized density proof does not close.  Malý's Cartesian theorem gives
strong minor approximation only after lowering exponents, while area-graph
approximation is an `L1` result.  A sharp audit excludes the direct radial
Malý dipole: `M_2 in L2` requires `alpha>1/2`, but non-vanishing filling cost
requires `alpha<=1/2`, including the logarithmic borderline.  This is not a
general density theorem or counterexample.

The requested relaxed-minimizer fallback also remains open.  The formal
stress is explicit, but local fixed-trace recovery is required before a
relaxed minimizer can be assigned the naive weak Euler--Lagrange system.  The
available 2024 relaxed-quasiconvex partial-regularity theorem requires `q<3`
for `(n,p)=(3,2)`; the present pointwise growth reaches `q=6` and adds an `S3`
constraint and degree relaxation.  Therefore endpoint density is neither
proved nor disproved, and weak EL, partial regularity, classicality,
continuum isolation, Hessian, determinant, and the degree-one portal all
remain false.

The `9/9` production card finds decreasing one-tetrahedron concentration and
a very tight perturbed-start basin after translation/target-`SO(3)`
quotienting.  This is only finite-grid evidence.  Ordered work is now an
annular target-valued endpoint replacement/equiintegrability lemma or a
reverse-Hölder `L^(2+delta)` estimate for the complete-minor vector, followed
by controlled target projection.  No Riemann Hessian is authorized first.

## AP-E11 compatible cochain update (2026-07-19; superseded above)

No Route-E dependency or portal is promoted.  AP-E11 replaces the refuted
one-corner stencil by a genuine tetrahedral cochain complex.  The
Alexander--Whitney cup has an exact graded Leibniz rule; its antisymmetrized
first-cochain products reproduce all affine second and third minors.  Positive
Whitney Hodge Gram matrices give one rotation-covariant action.  Contracting
the same ordered triple cup with the vertex field and the exact radial
Jacobian yields the normalized-affine degree cochain, so the energetic and
topological discretizations no longer use incompatible stencils.

Periodic null-Lagrangian identities and convexity prove the exact cell
density
`Q_(R,K)=|A|^2/2+R|M_2(A)|^2/2+K|M_3(A)|^2/2` with zero corrector on the
six-tet and both checkerboard five-tet cells.  A 13-background unanchored scan
passes the joint `(a,L)`, translation, and triangulation quotient thresholds,
with all three regular values giving `B=1`.

This closes the compatible-cell and relaxed-continuum gates, not the
classical one.  The proved Gamma limit is the fixed-degree lower-semicontinuous
relaxation.  Smooth fixed-degree density in the complete-minor graph norm, or
regularity and isolation of the selected relaxed minimizer, remains open.
Thus `C_B1_relaxed_continuum=true`, but the unrelaxed physical-background,
Hessian, determinant, and degree-one portal gates remain false.  The next
ordered task is that density/regularity theorem; Hessian work is still
prohibited until it closes.

## AP-E10 stencil/cell/translation update (2026-07-18; superseded above)

No Route-E dependency or portal is promoted.  AP-E10 instead proves that the
current AP-E7/AP-E8 mainline cannot supply the missing continuum background.
The exact three-periodic sequence
`n_a=(sqrt(1-eta^2 a^2), eta a cos(theta), eta a sin(theta), 0)` converges
uniformly to the vacuum while a one-corner forward Skyrme two-minor remains
`(3 sqrt(3)/2) eta^2`.  A separate two-periodic hemisphere-valued sequence
has geometric degree zero but forward topological current tending to
`-16 eta^3`.  An `L^(4/3)` estimate excludes singular concentration, but not
this diffuse oscillatory defect.  Hence the required compensated-compactness
and current-identification assumptions are false for the declared stencil.

The critical finite-`R` cell problem is fully evaluated: cycle balance and
Jensen force the periodic corrector to zero for both six- and five-tet graphs.
The resulting homogenized quartics are inequivalent, with six/five ratios
`4/3` and `3` on two gradients.  This proves finite-`R` mesh dependence rather
than removing it.

Eleven same-action quotient backgrounds extend the numerical range to
`a=0.208333`, `L=7.5`.  All are stationary/admissible `B=1`, but centered
size, profile, and cross-triangulation thresholds fail.  Thus
`C_B1_physical_continuum=false`, Hessian and determinant authorization remain
false, and the degree-one portal remains prohibited.  The required repair is
a compatible tetrahedral cochain/cup-product/Hodge-star Skyrme action whose
topological three-cochain agrees with normalized-affine degree.  This is a
new prerequisite before any further Route-E portal work.

## AP-E9 scaling/Gamma-limit update (2026-07-18; superseded above)

AP-E9 proves that the AP-E8 full action is equicoercive on every fixed box in
strong `L2 x L2`: the original coordinate-edge Dirichlet term bounds the
piecewise-affine `H1` norm.  This is ordinary compactness, not degree-sector
compactness.  With `d=1-epsilon`, `w=gamma*a`, and
`R=gamma*a^2/d^2`, the exact barrier calculus gives

`b_epsilon(n_x dot n_y) >= |n_x-n_y|^4/(8d^2)`

and a Lambert-`W` edge margin.  On every positive bounded sublevel, uniform
relative distance from the barrier floor holds iff `inf w>0`; smooth recovery
fields lose the barrier only when
`R->0`.  Thus for fixed epsilon and `gamma=c*a^(-p)` the only simultaneous
window is `1<=p<2`, with `p=1` preferred.  A shrinking degree-one bubble
shows that the vanishing barrier alone cannot preserve degree.  For
`inf R>0`, the barrier does give `W1,4` compactness and closes degree.
At critical `p=2`, the raw six-tet and five-tet smooth-sampling quartics are
not scalar multiples, but a checkerboard finite-`R` Gamma density must first
be computed by periodic cell homogenization; that identification remains open.

The complete same-action Gamma limit remains conditional on an unproved
Gamma-liminf/compactness theorem for the one-corner forward Skyrme minors and
energy-dense fixed-smooth recovery.  The finite-`R` barrier cell problem is a
separate missing theorem.  The 23-case production scan reinforces that
boundary: every field is stationary, admissible, and has `B=[1,1,1]`, but
the relaxed barrier slope (`0.2196`), large-volume radius spread (`8.521%`),
cross-triangulation energy spreads (`5.442%`, `5.355%`), and `p=0` versus
`p=1` energy spread (`6.446%`) fail their declared gates.

Therefore `C_B1_finite_grid=true` remains valid while
`C_B1_physical_continuum=false`,
`regulator_stable_continuum_background=false`, and both same-action Hessian
and determinant-variation authorizations remain false.  No Route-E dependency
or portal is promoted.  Ordered work is now the discrete Skyrme-minor theorem
and the periodic barrier cell formula, then a smaller-`a` adaptive-box scan
with translation tightness, and only after those succeed the full Riemann Hessian.
The finite GW/domain-wall kernel
and actual `SO(3)` mapping-torus mod-two index remain parallel lanes.

## AP-E8 topology-preserving finite-grid update (2026-07-17; superseded above)

AP-E8 deliberately chose the one route that changes the AP-E7 finite-site
blocker.  It adds a subtracted logarithmic barrier on every unique
Freudenthal co-cell pair.  If all pair dots exceed `epsilon=0.01`, the
normalized-affine denominator has the exact lower bound
`sqrt((1+3 epsilon)/4)=0.507444578255`; the finite-energy configuration space
therefore decomposes into degree components.  Compactness plus the coercive
breathing-field potential proves an interior unanchored minimizer in every
nonempty component, and a sampled compact-supported continuum degree-one map
proves such components exist on sufficiently fine grids.

The `13/13` numerical card independently re-relaxes seven no-pin backgrounds:
four spacings at fixed `L=6` and three volumes at fixed `a=0.4`.  Every direct
tangent-plus-scalar gradient density is below `7.8e-7`, every minimum pair dot
is above `0.0777`, and all three regular targets give `B=[1,1,1]`.  Four
barrier-parameter controls and one perturbed start remain stationary and
admissible.  The solver-continuation gradient and intrinsic barrier Hessian
are independently differentiated.

This changes one dependency statement only:
`C_B1_finite_grid=true`.  It does not import a Route-E portal or a physical
continuum soliton.  Equicoercivity, Gamma convergence, barrier/triangulation
independence, size/volume convergence, a physical QC2D action, the determinant
Hessian, interacting gauge/ghost blocks, BRST superdeterminant, four-dimensional
dynamics, and quantum continuum limit remain absent.  The fixed-`a` energy
spread is only `0.459%`, but the weighted radius changes by about `18%`, so
`C_B1_physical_continuum=false`.

The integrated authorization consequently remains
`any_preportal_route_closed=false`, `degree_one_portal_constructed=false`, and
`physics_promotion_allowed=false`.  Ordered work is: prove or refute a
regulator-stable continuum limit; test alternative triangulations and
regulator trajectories; only then construct the total same-action physical
Hessian.  The finite GW/domain-wall kernel and actual `SO(3)` mapping-torus
index remain parallel lanes, and the portal remains last.

## AP-E7 regulator/topology/family update (2026-07-17; superseded above)

No Route-E dependency is removed and no portal is authorized, but the next
three alternatives now have exact theorem boundaries.

- The common background-field quadratic APS/PV prescription and the pure
  two-flavour heavy gauge/gravity phase are explicit.  The latter is identically
  one because the relevant index is integral and doubled.  The physical APS
  sign uses a fixed zero cut with gapped endpoints; arbitrary-cut spectral flow
  is not used without its endpoint orientation.  The Higgsed `(+,+)` source is
  conditional and not yet a complete finite overlap Yukawa operator.  An exact
  separated rank-two selector cannot extend equivariantly through the restored
  `Sp(4)` point, while the full `5` has zero mixed charge trace.  Hence no
  original-action all-scale `k=+2` theorem is imported into Route E.
- At finite site number the unrestricted fixed-boundary `S3` configuration
  space is path-connected.  The five requested independent lattice
  relaxations unwind to the vacuum; guarded `B=1` descent hits the admissibility
  boundary rather than an interior critical point.  Positive fixed-core
  restricted spectra and `D_W dagger D_W` are controls, not a physical
  gauge--meson--ghost--fermion super-Hessian.
- On the physical collective orbit there are only the trivial line and one
  `Z2` FR torsion line.  Topology supplies a Real lift but not the microscopic
  Pfaffian class or CPT regulator.  Torsion cannot replace the free `O(2)`
  class.  A separate rank-three construction exactly realizes `c1=x+2y` with
  a uniform finite matrix gap and passes Berry-Chern numerics, but it still
  uses a nonnormalizable spectator `CP1` and is not the charged-two-colour
  same-soliton Yukawa family.

The integrated authorization remains
`any_preportal_route_closed=false`,
`degree_one_portal_constructed=false`, and
`physics_promotion_allowed=false`.  The three source cards pass `42/42`,
`16/16`, and `18/18`; the independent integrated gate passes `16/16`.

## AP-E6 same-background update (2026-07-16; superseded above)

No Route-E dependency is removed in this round, but three stronger facts now
replace the previous templates.

- The simply connected physics `Sp(4)` candidate has the exact gauge result
  `Omega_5^Spin(BSp(2))=Z2`; `4` detects its generator and `5,10` do not.
  This does not determine either target-space torsion sign.  The full
  Euclidean regulator and heavy/gravitational APS ratio remain open despite
  successful PV-moment and local-periodicity controls.
- A **20/20** canonical coupled continuum `B=1` card now exists with checksum
  `81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59`.
  Its analytically differentiated sparse cubic `n+s` Hessian is off shell;
  the diquark/gauge/ghost controls are not one same-grid aggregate.  The
  sampled cubic gradient is not small and the displayed lattice degree has
  not converged to one.  No physical stability or quantum-phase claim is
  imported into Route E.
- The **25/25** Callias card consumes that exact checksum.  It proves that the
  localized orbit is `SO(3)`, not `CP1`, and that the scalar-vacuum `CP1` is
  nonnormalizable.  The ordinary endpoint Callias index is zero.  A separate
  two-band theorem rules out `c1=x+2y` in a trivial rank-two gapped mass
  family because its square is nonzero.  The WZW expression is now a typed
  restricted ansatz, not full equivariant/gauge-basic descent.

The integrated card passes **15/15**.  Consequently the `Sp(4)`,
on-shell-Hessian, and same-soliton determinant-line
closures are all false; no degree-one portal is authorized or constructed.
The preferred continuation is to compute a complete APS threshold, re-relax
each cubic solution before its Hessian, and quantize the actual `SO(3)`/FR
collective coordinate (with a rank-three mass-family branch in parallel).

## AP-E completion-frontier update (2026-07-16; superseded above)

The next three pre-portal lanes have been executed and integrated without
promoting a Route-E family claim.

- The new four-dimensional note defines a Wilson/Symanzik-compatible charged
  two-colour action and the complete **formal**
  gauge--meson--ghost--fermion super-Hessian. Its **24/24** finite '2^4'
  diagnostic checks topology extrapolation, gauge/ghost cancellation,
  declared boson blocks, a Wilson--Dirac Schur factorization, one log-det
  curvature, and cutoff orders. It is not Monte Carlo, determinant
  positivity, the exact nonlinear Skyrme Jacobi operator, or a continuum
  quantum calculation.
- On explicit 'G3=S1_R x S3' and 'G2=T2_RR x S2' representatives, the
  ambient product spectra give only a reference '(+1,+1)' phase. Defect
  mod-two regulators realize all four torsion characters, so the present
  charged-QC2D fields do not select a unique APS/Dai--Freed pair. The
  **58/58** audit keeps both microscopic eta signs false/open.
- The strongest simple group in the scanned charge/global-form set is the
  simply connected 'Sp(4)=Spin(5)': its '5' contains
  '2_(+1)+2_(-1)+1_0' and its '4' contains
  '2_0+1_(+1)+1_(-1)'. Two complex '4' copies supply the separate global
  'SU(2)_phi' doublet, so 'qq(phi dagger)^2' is a neutral triplet at group
  level. The older exactly-two unit-cell Gauss law is not rederived.
  'b0=15/2>0', but radiative threshold protection, the heavy eta match,
  full gauge bordism, and the nonperturbative 'B=1' phase remain open.
- On a hypothetical same 'B=1' family, the conditional Callias formula is
  'rank=epsilon p', 'c1(det Ind D)=epsilon p q'; rank-one 'O(+2)' needs
  'epsilon p=1,q=2'. Fixed-polarization CPT requires the Serre-dual
  'K tensor O(2)^vee=O(-4)', not the one-state 'O(-2)' inverse. Spatial WZW
  fiber integration lands in differential degree two, but equivariant/basic
  descent, the physical Yukawa family, and the same-soliton pullback are not
  derived. This card passes **27/27** while all composition gates remain false.
- The integrated **20/20** gate computes all three complete-lane booleans as
  false. Therefore no degree-one portal has been constructed and
  'physics_promotion_allowed=false'.

Next: compute the complete 'Sp(4)' (direct product as backup) regulator and
mapping-torus eta phases; relax and scan the same charged 'B=1' solution with
the exact projected Hessian and dynamical lattice ensembles; then calculate
its Callias/CPT determinant and equivariant WZW descent. Only after one full
lane closes may a degree-one Route-E portal be written.

## Superseding AP-E3/AP-E4 global/SQM update (2026-07-16)

This section supersedes the active-continuation language in the 2026-07-15
checkpoint.  It does not close a Route-E dependency or authorize a family
interpretation.

- The charged-two-colour lane now has a reproducible nonlinear classical-EFT
  proxy: a meson/charge-two-diquark vacuum scan, complete homogeneous Hessian,
  nonlinear `B=1` Skyrme boundary-value solution, Derrick check, radial meson
  generalized Hessian, and charged-scalar `l=0,1` eigenvalue checks.  These are
  necessary low-energy stability tests.  They are neither a lattice solution
  of the four-dimensional gauge theory nor a full coupled 3D/quantum Hessian.
- The mixed WZW phase is globally defined on non-extendible maps by a
  degree-five differential character on `S3 x S2`.  The bosonic lift is unique
  at fixed integral curvature, but the spin theory has
  `Omega_4^Spin(S3 x S2)=Z + Z2 + Z2`; two reduced APS/Dai--Freed signs remain
  microscopic UV data.  Thus “non-extendible action defined” and “fermionic
  UV determinant fixed” are separate statuses.
- A strict `SU(3)` representation-parity theorem rules out the original
  charge-one colour-singlet scalar: restriction to
  `[SU(2) x U(1)]/Z2` requires `2j+q=0 mod 2`.  Consequently a pure-`SU(3)`
  completion cannot preserve the original two-scalar `O(2)` dressing.  The
  closest explicit adjoint-decoupling variant uses two `bar(3)` scalars with
  light `1_(+2)` components and vectorlike fundamentals with tunably heavy
  `1_(-2)` partners.  It is anomaly-safe at the perturbative vectorlike level
  but produces a one-scalar flavour doublet, not Route E's triplet.  Its
  covering-`U(1)` WZW normalization also remains conditional until the
  faithful-`U(2)` bundle and correlated centre-flux classes are matched.
- An independent supersymmetric non-Abelian-vortex construction now derives a
  physical `CP1` tangent Grassmann zero mode and its `N=2` worldline SQM.
  The declared canonical Spin-c re-quantization gives the untwisted
  Dolbeault module and one ground state; the source-selected half-form
  ordering has none.  Three states require a separate `O(2)` coefficient line.  The AP-E3
  WZW level can provide that line only if the charged-two-colour `B=1` soliton
  is proved to have the same supersymmetric collective-coordinate theory and
  its degree-five WZW character is spatially transgressed to that line;
  the current mother models differ.  Product compactification stays as the
  six-dimensional anomaly-gated backup.  The mother model also does not
  derive the canonical vacuum line or its CPT map: fixed-canonical `O(-2)`
  has one negative mode, whereas three conjugate modes require the
  anti-canonical polarization/effective `O(-4)`.
- The degree-one Route-E portal is deliberately deferred until the same-model
  SQM/WZW composition or the compactification alternative closes.  All new
  cards remain fail-closed with `physics_promotion_allowed=false`.

Canonical evidence is in the new `../route_f/tex/ap_e3_*` and
`../route_f/tex/ap_e4_moduli_space_sqm.tex` notes, their verifiers under
`../route_f/code/`, and generated JSON/Markdown cards under
`../route_f/output/`.

## Superseding AP-E3/AP-E4 update (2026-07-15)

This section supersedes the AP-E3/AP-E4 open-status statements in the
2026-07-14 historical checkpoint below; it does not supersede any Route-E
dynamics blocker.

- Bosonic Schwinger partons with canonical CCR and two independent compact
  cell constraints `G_r=N_r-1=0` prove an exact `(N_1,N_2)=(1,1)` physical
  Hilbert within the declared indivisible unit cell.  Odd, singleton, and
  `(2,0)/(0,2)` sectors are absent.  Rank-ten/rank-twenty negative controls
  refer to the declared `0<=N_r<=2` audit truncation; the physical rank four
  is cutoff-independent.  The proof covers complete-cell positive-parent
  boundaries, not arbitrary Haldane-like intercell phases.
- The physical electron Pauli sign selects the anti-aligned Hund-triplet
  line.  With the fixed AP-E1 connection `A_+`,
  `i hbar <Omega_-|d Omega_->=+2 hbar A_+`; the line is
  `Q tensor Q=O(2)` and the signed level is `k=+2`.  This closes the old
  orientation blocker without manually reversing a Berry term.  The
  exactly-two/orientation verifier passes `26/26`.
- An explicit `SU(2)_c x U(1)_g`, `N_f=2` mixed-WZW candidate cancels all
  dynamical local and Witten gauge anomalies and has
  `kappa_L=-kappa_R=2`.  Its five-dimensional coefficient is `n=2`, and a
  positive unit baryon gives `k=nB=+2`.  It is only an anomaly-consistent
  intermediate UV completion: strong-vacuum selection, Pauli--Guersey/
  diquark gaps, compact-`U(1)` and bordism effects, soliton stability, and an
  all-scale completion beyond the abelian Landau pole are open.  In addition,
  the five-dimensional action is proved extension-independent only on
  extendible sectors; a differential-cohomology definition for non-extendible
  sectors remains open.
- AP-E4 now independently proves `T^(1,0)CP1=O(2)` and solves the declared
  canonical Spin-c operator.  At `R=1/2`, it has three positive zero modes,
  no negative kernel, and
  `lambda_(n,+/-)=+/-sqrt(n(n+3))/R` with multiplicity `2n+3` per sign; the
  first gap is `4` and every massive level is paired.  The `22/22` audit also
  checks `lambda_(3)=6 sqrt(2)` and distinguishes eigenvalue pairing from
  chirality.  It proves the ordinary-spin negative control: twisting only by
  `T=O(2)` gives
  two zero modes and gap `2 sqrt(3)`.  The canonical Spin-c half-canonical
  shift, fermionic origin, higher-dimensional anomaly cancellation, and
  family interpretation remain underived.
- Route E still may not import a three-family conclusion.  A degree-one
  portal, a physical moduli/compactification origin for the AP-E4 operator,
  complete anomaly cancellation, and stability below all gaps remain
  mandatory.  `physics_promotion_allowed=false`.

Canonical new evidence:
`../route_f/tex/ap_e3_exact_two_mixed_wzw_uv.tex`,
`../route_f/tex/ap_e4_tangent_dirac_spectrum.tex`, and their JSON/Markdown
verifier outputs in `../route_f/output/`.

## Superseding blocker execution update (2026-07-14)

This section overrides every later historical `done`/`closed` label for the
dynamics lanes.  The authority is `code_dyn/dyn_claim_registry.json`.

- The canonical case-sensitive path is `route_E/code_dyn/`.  The 21-node
  isolated runner now passes a full required-card dry-run with all 19 recovered
  scripts plus DYN-5V/DYN-7F and no preflight errors.  Root DYN scripts are
  compatibility delegates so they cannot overwrite corrected ledgers with a
  second implementation.
- RE-SC3/4/5 exist and are tracked (`15/15`, `14/14`, `18/18` mechanical
  checks), closing only the artifact-absence blocker.  They remain unpromoted;
  RE-SC4's numerical selectivity gap inherits invalid DYN-5 input, and RE-SC5
  is a one-loop toy scan whose lower G_LR edge has `M_I~=101 GeV`.
- DYN-9b-2 is `preliminary` (`14/14` mechanical), not a flavor fit.  Its double-normalised
  hypercharge input was repaired: with `g1=0.462` and `b1=41/10`,
  `y_t(M_X)=0.44116481` rather than `0.41556342`.  The relations
  `Y_nu=h-3f`, `Y_u=h+f` do not force a top-like `Y_nu` (`h=3f` is a
  counterexample): the fixed archival-kernel suppressions are `19.5x/342.2x`,
  whereas `9.6x/169.2x` applies only after an extra top-like ansatz.  Exact
  zeta invariance is restricted to uniform positive-real rescaling; complex
  rescaling rotates its phase by `y^2/|y|^2`.
- DYN-9b-3 remains `blocked_missing_branch_thermal_inputs` (`11/11`
  mechanical).  The corrected SM
  Davidson--Ibarra bound is `2.307794855e-6`; `M1=2.38045e10 GeV` is in the
  tau-resolved regime, and removal of a gravitino-specific ceiling does not
  make reheating unconstrained.
- DYN-8 is `preliminary`: `30/30` mechanical/disclosure checks pass while
  `physics_promotion_allowed=false`.  Its 210/45 contradiction is repaired;
  rare sampled 210-only LR minima make the 45-adjoint route an alternative,
  not a necessity.
- The H3/Cartan logic blocker is resolved by theorem demotion plus an explicit
  H3+ axiom: the one-dimensional abelian branch has a non-degenerate
  ad-invariant form, so original H3 proves only `N_fam<=3`; `N=3` is now
  explicitly conditional on a non-degenerate adjoint-trace/Killing contact.
  The physical motivation or UV realization of H3+ remains open.
- The Route-F bridge audit sharpens the selected-branch mathematics:
  `B(x,x) = 2 Delta = 2 sqrt(3) x^T K_tr x` in Route-E spherical
  coordinates exactly identifies the two-center
  locus with regular semisimple/non-null Killing contact, while `Delta=0` is
  the nilpotent double-zero cone.  A complete heavy-adjoint current loop can
  conditionally produce the Killing tensor, but it does not prove the
  existence of that sector or its conversion into a Majorana contact; H3+
  therefore remains open.  Full derivation and rollback trail:
  `../route_f/tex/another_physics_route_e_derivation_ledger.tex`.
- AP-E1 now proves the local fixed-norm doublet quotient
  `S3/U(1)=CP1` and the induced Fubini--Study action, but also proves that a
  physical global Q-ball phase cannot be removed pointwise.  Fixed-charge
  reduction gives a monopole rotor with all Landau levels; the Hopf line has
  Chern number one, not Route E's tangent-bundle value two.  Identifying its
  level with the existing `Q=10^6` benchmark would give `1,000,001` states,
  not three.  Consequently Route E may import the CP1 geometry but may not
  import `O(2)` or a family count unless AP-E3 closes the exactly-two
  UV/singleton rule and AP-E4 derives the independent chiral operator and gap.
  The canonical proof is
  `route_f/tex/ap_e1_projective_doublet_action.tex` (`30/30`
  arithmetic/source regressions, non-promoting).  The corrected Branch-B
  first-order sign gives `+k F`; `k=0` is a point rather than `CP1`, and the
  signed rotor spectrum depends on `|k|`.
- AP-E2 is now an exact mandatory regression: `30/30` deterministic checks
  preserve `B_Kill=2 Delta=2 sqrt(3) x^T K_tr x`, `A_q^2=Delta I/4`, and the
  transvectant normalization across `SL(2)` transformations, complex and
  infinite roots, the nonzero nilpotent boundary, zero-section exclusion, and
  both basis conventions.  It still derives neither H3+ nor a Berry level and
  sets `physics_promotion_allowed=false`.
- AP-E3 supplies a conditional level-magnitude-two *candidate*, not a Route-E
  dependency closure.  For two declared Bose-Hubbard orbitals, `0<mu<U`
  selects unit filling only in the bare onsite problem.  In the complete
  declared `H_portal=0` model, writing `C=mu+h/2+J_H/4`, the sufficient
  interacting `(1,1)` plateau condition is `C<U-J_H/8`, equivalently
  `mu<U-h/2-3J_H/8`.  Ferromagnetic `-J_H S_1.S_2` then gives the diagonal
  symmetric pair `|s;2>=|s> tensor |s>` and doubled Berry curvature/metric.
  With `A=-i<s|ds>` and `+i hbar a^dagger dot(a)`, the declared aligned
  `-h n.S` action has `k=-2`: its ket eigenline is `O(-2)`, whereas the dual
  prequantum line is `O(2)`, has `c1=+2`, and supports the three-state
  spin-one carrier.  Reversing orientation/coupling gives `k=+2`; AP-E3 has
  therefore derived `|k|=2`, not the signed chirality.  The `27/27` audit also
  verifies benchmark plateau-condition margin `2.025`, full spinful ground
  sector `(1,1)`, interacting charge gap `1.85`, singlet gap `J_H`, and a
  coercive occupancy tail when `4U>J_H`.  However the theory has not derived
  exactly two orbitals, excluded the unwanted singleton `|k|=1` sector, built
  an anomaly-safe four-dimensional portal, selected the `k=+2` orientation,
  or identified the triplet with chiral families.  Thus
  `ap_e3_physics_closed=false` and Route E remains barred from importing the
  family count.
- The next active bridge gate is AP-E4: construct the actual fluctuation/Dirac
  operator, prove or refute tangent-valued chiral zero modes, enumerate all
  partners, and compute the spectral gap.  Filled-fermion, mixed-WZW, and
  tangent/Dirac level mechanisms remain alternatives until one closes its own
  UV and anomaly gates.

Reproducible promotion diagnostics are recorded in
`../route_f/output/blocker_promotion_gate.{json,md}` (`18/18`, deliberately
with `physics_promotion_allowed=false`).

## A. Removable dependencies (prove in-repo)

- **R1 [remove] CI-1, complex-representation classification.**  Which simple
  compact groups admit complex irreps, and the concrete self-conjugacy rules
  (A_n label reversal, D_odd spinor swap, E6 arm swap, all others
  self-conjugate).  Closed by computing the longest Weyl element w0 from the
  Cartan matrix for every simple algebra of rank 4-15 plus E6/E7/E8/F4 and
  reading off the involution -w0.
  Status: **done** (2026-07-04, `verify_route_e_dependency_closure.py` §2;
  tex Sec. III item (i)).
- **R2 [remove] Enumeration completeness.**  The dim<=16 scan now runs over
  ALL simple algebras of rank 4-15 (B, C, D-even, E-series, F4 included)
  with a generic Cartan-matrix Weyl dimension engine; the complex dim-15/16
  list is unchanged; the so(9) real 16 and so(15) real 15 near-misses are
  recorded (fail chirality, not dimension).  Rank>=16 corroborated by the
  rank-16 smallest dims 17/33/32/32 (machine) + Slansky minimality (cited).
  Status: **done** (2026-07-04, audit §3-4; tex Prop. III.4 +
  Rem. III.5 + Appendix B).
- **R3 [remove] CI-2, SU(N)-fundamental cubic anomaly nonzero.**  Exact
  integer traces: tr T^3 = 2730 (SU(15)), 3360 (SU(16)).
  Status: **done** (2026-07-04, audit §5; tex Sec. III item (ii)).
- **R4 [remove] CI-3, Spin(10)-16 anomaly freedom (perturbative).**  Chiral
  16 built from five fermionic modes (chirality = occupation parity);
  tr({T_a,T_b}T_c) = 0 exactly over all 45^3 triples; Cartan weights = the
  demicube; index(16)/index(10) = 2 sanity.  The Witten global anomaly
  (pi_4(Spin(10)) = 0) stays textbook-cited [see R7/R11].
  Status: **done** (2026-07-04, audit §6; tex Sec. III item (iii)).
- **R5 [remove] CI-4, exceptional minima.**  E6 = 27, E7 = 56, F4 = 26,
  E8 = 248 by scan.
  Status: **done** (2026-07-04, audit §4; tex Sec. III item (iv)).
- **R6 [remove] Route-D dependency in Theorem 5's witness family.**  The
  admissible model class is now the internal Route-B EFT family
  {M_lambda : lambda in C*} (P0 holds because X is a gauge singlet and the
  16 content is anomaly-free by the exact cubic-trace computation, tex Sec. III item (iii)); zeta ranges over C*.  Route D is
  cited only as an optional UV interpretation.
  Status: **done** (2026-07-04, tex Def. VII.2 + Sec. II.A closing note).
- **R15 [remove-by-attribution] Prior-art provenance.**  Relation-to-prior-
  work subsection added with verified citations: Georgi-Glashow PRD 6 (1972)
  429; Okubo PRD 16 (1977) 3528; Witten PLB 117 (1982) 324; Baez-Huerta
  Bull. AMS 47 (2010) 483 (arXiv:0904.1556); Candelas-Horowitz-Strominger-
  Witten NPB 258 (1985) 46; King-Luhn RPP 76 (2013) 056201
  (arXiv:1301.1340).  New-assembly claims stated without overclaiming.
  Status: **done** (2026-07-04, tex Sec. I.A; six entries appended to
  `paper/refs.bib`).

## B. Retained inputs (fully explained, tagged in the tex)

- **R7 [explain] The principle set itself.**  P0-P3 and E1-E2 are axioms
  (retained inputs (i), (ii)); P0 scope: 4D perturbative anomalies machine-checked (R3/R4),
  Witten global anomaly cited.
  Status: **done** (tagged, tex Sec. II.A provenance ledger).
- **R8 [explain] Minimality tie-break in Theorem III.7.**  Occam selection
  clause, tagged as retained input (iii); dropping it reopens e.g. the E6 27 with vectorlike matter.
  Status: **done** (tagged as retained input (iii)).
- **R9 [explain] Route-B EFT inputs.**  U(1)_R charge assignment (retained input (iv); the
  Majorana-only rule is derived from it in Route A), projector P_{nu^c}
  (v), mu != 0 and no-extra-unpaired-sector clause (vi).
  Status: **done** (tagged as retained inputs (iv)-(vi)).
- **R10 [explain] Route-A benchmark data.**  zeta digits, windows, Z_{N<=6}
  miss, N = 178, Audit-0.5 invariants: corroboration only, premises of
  nothing.
  Status: **done** (tagged, provenance ledger + Obs. VII.6).
- **R11 [explain] Textbook theorems.**  Weyl dimension formula, RR + Serre,
  Cartan criterion, CG/Schur, c_1 = zero divisor, Dolbeault-Hodge,
  pi_4(Spin(10)) = 0, rank>=17 minimal-dimension tables.
  Status: **done** (tagged, provenance ledger).

## C. Provably not closable

- **R12 [impossible] The value of zeta.**  Theorem VII.4: no zeta-blind
  principle set entails zeta = zeta_0; any closure is a new conditional
  axiom (moduli stabilization / hidden phase-radial data) that must be
  audited on its own.  Status: **permanently-open (by theorem)**.

## D. Out of Route-E scope (inherited, untouched)

- **R13 [scope]** Route-A companion audits: full flavor fit, d=5
  proton-decay Wilson tensors, thresholds, source-sector UV origin;
  dependency graph 0 -> (1 || 4a) -> 3 -> 2, 4b optional.
  R13 is now being executed as the dynamics program below (section E),
  which supersedes this placeholder.
- **R14 [scope]** Promotion decision for Route E into the main paper
  (author's call; criteria in README.md).

## E. Dynamics program (planned 2026-07-04; items DYN-0..DYN-8)

> **Historical execution log.**  The original `done`/`closed` labels and
> numerical interpretations below are retained for provenance only.  They are
> superseded by the 2026-07-14 blocker update and
> `code_dyn/dyn_claim_registry.json`; in particular DYN-5, DYN-7, and DYN-8
> must not be promoted from this section.

The reconstruction paper is deliberately kinematic (representation theory +
index theory only).  This program adds the action-level dynamics along the
audit graph fixed by Audit 0, `0 -> (1 || 4a) -> 3 -> 2`, extended by a
hidden-sector lane, an explicitly-new-axiom lane for zeta stabilization, a
cosmology lane, and a falsifiability section.

Standing discipline for every DYN item:
- plain `python3` + `numpy` (no sympy); each script writes a JSON + MD
  ledger with negative-boundary flags and prints a one-line summary;
- every physics result is conditional-on-benchmark; nothing enters the
  theorem core of the reconstruction paper;
- Theorem VII.4 caps the program: dynamics may *constrain* zeta
  observationally (DYN-7) or price stabilization axioms (DYN-6), but a
  zeta-blind principle set can never derive it — no lane may claim
  otherwise;
- history assets are recovered read-only via `git show HEAD:<path>`;
  the deleted working-tree files are never blanket-restored.

- **DYN-0 [scaffold + conventions] (1 session).**
  Deliverable: dynamics inventory ledger + Audit-0 conventions addendum.
  Steps: (1) inventory of reusable history assets (initial list recorded
  below per item); (2) conventions addendum: superpotential parameters
  (m, M, lambda, eta, gamma, gamma_bar, h), PS-singlet VEV names
  (p, a, omega, sigma, sigma_bar), reaffirm the fixed Aulakh-Girdhar
  convention x = -lambda*omega/m, xi = lambda*M/(eta*m),
  8x^3 - 15x^2 + 14x - 3 = -xi(1-x)^2; (3) rerun-from-history gates:
  the 4a1 vacuum-cubic validator, literature mass import, and triplet
  symbolic inverse must reproduce their HEAD ledgers.
  Decision points for the author: (a) file lanes — recommended: new
  scripts at root `code/` with outputs in `output/audit{1,2,3,4a,4a1}/`
  (extend the .gitignore whitelist for audit3 outputs deliberately;
  audit2 outputs likewise if new files are added there); alternative:
  a self-contained `dynamics/` sandbox; (b) SUSY-scale treatment —
  CMSGUT is supersymmetric, so fix M_S handling (single benchmark vs
  scanned 1-10 TeV).
  Assets: `code/audit0_conventions_card.py`,
  `code/audit4a1_cmsgut_vacuum_branches.py`,
  `code/audit4a1_cmsgut_literature_mass_import.py`,
  `code/audit4a1_triplet_symbolic_inverse.py`,
  `output/audit4a1/*` (all in HEAD).
  Status: **done** (2026-07-04; `code/audit0_dyn_conventions_and_inventory.py`,
  13/13 checks; ledger `output/audit0/dyn0_conventions_inventory.{json,md}`.
  Decisions adopted: (a) root `code/` + whitelisted `output/audit*/` lanes,
  `.gitignore` whitelist extended for `output/audit0/`; (b) M_S = single
  effective threshold, benchmark 3 TeV, DYN-2 fit window [1, 10] TeV, MSSM
  two-loop running below M_GUT.  Rerun-from-history gates: all four
  audit-4a1 scripts extracted from HEAD reproduce their ledgers up to
  volatile fields; 38 inventory assets verified present in HEAD; the four
  vacuum-cubic special points re-verified exactly.  DYN-1 starting
  checklist snapshotted: heavy spectrum placeholder + scalar-Hessian
  Goldstone eigenvectors pending.)

- **DYN-1 [Audit 4a completion: source-sector vacuum dynamics]
  (2-3 sessions).**
  Goal: the SM-face vacuum branch of the CMSGUT superpotential
  (210 + 126 + 126bar + 10) and a NON-placeholder
  `output/audit4a/heavy_spectrum.json`.
  Steps: (1) branch classification: solve the vacuum cubic over complex x
  on a benchmark xi grid; identify the branch preserving exactly G_SM;
  record all other branches (SU(5), flipped SU(5)xU(1), G_LR, PS) as
  negative-boundary entries; (2) scalar Hessian / Goldstone eigenvector
  audit — the item the companion note explicitly lists as pending:
  assemble the full PS-block scalar mass matrices (imported 4x4 doublet
  and 5x5 triplet blocks, mixed G/E/F/J/X blocks, plus the remaining
  blocks derived here), verify the 33 = 45 - 12 Goldstones as actual
  zero-eigenvectors of the Hessian, and no tachyons on the physical
  branch modulo flat directions; (3) heavy-spectrum export: diagonalize
  at benchmark couplings and emit mass / SM irrep / multiplicity /
  origin-block records conforming to the fixed audit-4a schema;
  (4) doublet-triplet tuning: impose det(M_doublet) = 0, export the
  light-doublet composition (alpha_i, alpha_bar_i) needed by DYN-3
  dressing and DYN-4 Yukawa matching.
  Gates: exact Goldstone count 33; eigen-residuals < 1e-10; the three
  literature special points reproduced; schema validation; flags:
  `no_unique_vacuum_claimed`, `tree_level_only`.
  Assets: the four audit4a1 scripts; `code/solve_spin10_vacuum_alignment.py`,
  `code/audit_orbit_superpotential_hessian.py`,
  `code/verify_spin10_component_hessian.py`,
  `code/verify_combined_superpotential_flatness.py` (HEAD).
  Unblocks: DYN-2, DYN-3.
  Status: **DYN-1a done** (2026-07-05;
  `code/audit4a_dyn1a_vacuum_goldstone_audit.py`, 23/23 checks; ledgers
  `output/audit4a/dyn1a_vacuum_goldstone.{json,md}` +
  `output/audit4a/heavy_spectrum.json`).
  Done in 1a: (i) the singlet superpotential is now DERIVED, not
  transcribed -- explicit W with c_sigma = 1 has all five F-terms vanishing
  identically on the AG branch (worst residual 1.6e-15 over 200 grid
  points) and its Hessian EQUALS the transcribed neutral G block under the
  field map (p, sqrt3 a, sqrt6 w, i sqrt2 sigma, -i sqrt2 bar_sigma),
  max deviation 4e-16; (ii) branch classification: SM window x in (0,1/3),
  and the three enhanced-symmetry boundaries read off the gaugino-column
  vector masses -- x->1/3: G+E massless (flipped SU(5)xU(1), 13 vectors),
  x->0: G+F massless (G_LR, 3 vectors), x->1/2 on the complex-sigma slice:
  X massless (SU(5), 12 vectors) -- which independently validates the
  sector assignments; (iii) the previously-pending Goldstone eigenvector
  audit is CLOSED: in all five sectors the chiral null vector IS the
  gauge-orbit direction (alignment residual exactly 0), J^dagger-J has one
  zero mode and no tachyons, super-Higgs full blocks non-singular,
  count 1+2+6+12+12 = 33 = 45-12; (iv) doublet-triplet: det-linear tuning
  M_H* = -0.6258+0.1426i, light pair composition |alpha| =
  (0.662, 0.119, 0.733, 0.098) exported, triplet block non-singular
  (singular values 4.80..0.086); (v) heavy_spectrum.json non-placeholder:
  29 benchmark-mass states across 7 sectors, schema-conform, partiality
  flagged.  Sector letters source-verified against the fetched AG arXiv
  source (items a-e, Y_AG = 2 Y_SM; sha256 4022df72... staged in
  scratchpad).
  **DYN-1b done** (2026-07-05; `code/audit4a_dyn1b_full_spectrum.py`,
  16/16 checks; ledgers `output/audit4a/dyn1b_full_spectrum.{json,md}`;
  `output/audit4a/heavy_spectrum.json` now COMPLETE: 57 states, all 26 AG
  sector types).  Method: two independent computations gated against each
  other -- (i) a PS->SM decomposition engine derived from first principles
  (SU(4)->SU(3)xU(1)_B-L and SU(2)_R->T3R branchings, Y_SM = T3R+(B-L)/2)
  applied to 210+126+126b+10 (472 states) and to the 45 (12+33); (ii) the
  transcribed AG census (Table I: 21 unmixed masses/19 letter types;
  R/h/t mixed pure chiral; G/E/F/J/X mixed chiral-gauge with m_lambda
  gauge-multiplet masses, items i-v).  Completeness gate: the two multisets
  agree exactly.  Extra gates: all 26 threshold-index rows {S3,S2,S1}
  recomputed from Dynkin indices match AG Table 2 incl. the printed
  S_W/S_X combinations; R closed-form eigenvalues match the matrix; the
  five m_lambda formulas equal the transcribed gaugino-column norms to
  1e-12; at the benchmark exactly one massless level exists (det-tuned
  Higgs pair); 9 accidental zero-mass loci in x in (0.18, 0.30) are
  scanned and DISCLOSED (benchmark x = 0.1 sits below all of them; DYN-2
  must respect them).  Vacuum-consistency insight recorded: on the branch
  M+eta(p+3a-6w) = 0 so m_A = 24 eta w exactly (= 2.4 at benchmark).
  The full derivation (source line ranges, Y-normalization, branching
  tables, conventions, formulas) is recorded in the ledger's
  `derivation_log` and the human-readable MD.  DYN-1 is CLOSED; remaining
  spectrum items (GeV normalization, vector-multiplet threshold
  bookkeeping, AG Table-2 b-vector usage) belong to DYN-2.

- **DYN-4 [Audit 1: flavor/seesaw dynamics] (2-3 sessions; runs in
  parallel with DYN-1 per the fixed graph).**
  Goal: fit the sl2-covariant Yukawa sector (Y_10, Y_126bar, optional
  Y_120) to charged-fermion masses + CKM + neutrino data with
  M_R = M_V + zeta K_tr; produce a data-driven posterior for zeta
  (modulus and phase) — the missing error budget.
  Steps: (1) target-table refresh (the audit-1 ledger marks targets as
  requiring refresh): current PDG/NuFIT central values + uncertainties at
  fixed conventions; (2) covariant parametrization in the
  (psi_0, psi_1, psi_2) basis, enforcing the Audit-1b rank contract
  (>= 3 aligned tensors, tensor-by-tensor residual reporting);
  (3) chi^2 fit + type-I seesaw replay; reproduce the benchmark contact
  decomposition (contact fraction 1.304275e-1) as a regression gate;
  (4) sensitivity refresh: recompute the Delta_s / Delta_delta windows
  from the new fit; publish the zeta posterior.
  Gates: benchmark replay bit-comparable; per-observable pull table;
  flags: `fit_not_prediction`, `conventions_fixed_by_audit0`.
  Assets: `code/audit1_flavor_target_conventions.py`,
  `code/audit1b_covariant_rank_contract.py`,
  `code/audit_flavor_fit_observables.py`,
  `code/scan_majorana_contact_sensitivity.py`,
  `code/audit_cp1_yukawa_covariant.py` (HEAD).
  Unblocks: DYN-3, DYN-5, DYN-7.
  Status: **DYN-4a done** (2026-07-05;
  `code/audit1_dyn4a_seesaw_replay_zeta_posterior.py`, 19/19 checks;
  ledgers `output/audit1/dyn4a_seesaw_zeta_posterior.{json,md}`).
  Provenance closed: the author restored the archival output set under
  `route_E/output/` (sha256 recorded); using ONLY the recipe preserved in
  HEAD code, the paper anchors REGENERATE bit-level -- zeta =
  0.1076472949+0.0736514853i, contact fraction 1.304275e-1, and the
  historical sensitivity windows to all printed digits (loose 3.208798e-5 /
  5.424119e-5 as HALF-WIDTHS of the two-sided intervals; tight 9.9379e-7 /
  1.7993e-6).  Conventions recorded: c_hat = -K_tr(paper) orientation,
  M_* = sigma_max = 3.9265e15 GeV, v_u = 100 GeV; heavy Majorana spectrum
  (3.93e15, 3.22e13, 2.38e10 GeV) exported for DYN-7.  Audit-1b rank
  contract satisfied with an orthonormal 6-tensor cascade (tail 1.8e-16).
  Target refresh: NuFit-6.0 NO (IC24-with-SK primary, IC19 variant); the
  sin^2 theta_23 octant flip (+6.9 sigma) dominates; all other shifts
  < 0.4 sigma.  FINDINGS: benchmark chi^2 vs NuFit-6.0 = 47.5, entirely
  theta_23-driven; the zeta-plane scan shows zeta_best = zeta_benchmark
  (chi^2 47.46) -- the zeta direction is ORTHOGONAL to the octant-flip
  tension and cannot absorb it, and the conditional posterior is
  razor-thin around the benchmark (Delta-chi^2 windows ~ 1e-6);
  delta_CP prediction (not fitted): branches {324, 216} deg, the 216 deg
  branch is 0.1 sigma from NuFit.  Consequence: the covariant Dirac-sector
  refit is demonstrably necessary to address theta_23.
  **DYN-4b done** (2026-07-05;
  `code/audit1_dyn4b_refreshed_card_unconditional_zeta.py`, 13/13 checks;
  ledgers `output/audit1/dyn4b_unconditional_zeta.{json,md}`).
  Design: DYN-4a proved zeta (M_V frozen) orthogonal to the theta_23 flip,
  so 4b frees the FULL Majorana sector -- the inverse seesaw absorbs
  the NuFit-6.0 central inputs exactly and propagates the covariant
  decomposition over declared nuisance priors.  No likelihood or draw
  reweighting was applied, so this is not a posterior.  Refreshed central
  card: zeta' = 0.105557+0.071595i
  (arg moved 75 old loose windows, |zeta| moved 689 -- the printed anchor
  digits are decisively superseded by DATA, exactly as the boundary
  theorem requires), contact fraction' = 0.1275, M_*' = 4.09e15 GeV,
  I' = 0.007790+0.017991i, J' = -9.143e-5+5.142e-4i; Audit-0.5 re-test
  still no-hit (miss 0.566 rad).  UNWEIGHTED PRIOR-DRAW REGRESSION
  (NuFit-6.0 Gaussian proposal factors x uniform Majorana phases x
  log-uniform m1 in [1e-4, 3e-2] eV, 4000 draws): |zeta| = 0.144
  [0.113, 0.184] (central 68% draw interval), arg zeta = 0.605
  [0.465, 0.749] rad (circular).  CONDITIONAL ROBUSTNESS RESULT: within
  this declared draw distribution, the sample fractions with contact
  fraction > 0.01 and > 0.05 are 1.000 and 0.988; this does not constitute
  a likelihood posterior or a global exclusion of the Veronese-only branch.
  Nuisance
  attribution: the alpha_21 Majorana phase dominates the width
  (corr +0.53), m1 +0.19.  DYN-7 feed: log10 M_1 = 10.44 [10.39, 10.52],
  the sample fraction with M_1 > 1e9 GeV is 1.000 (a necessary
  Davidson--Ibarra scale diagnostic, not a leptogenesis-viability claim).
  Old-card I/J regress to the Audit-0 paper digits.
  Remaining as **DYN-4c** (open): kernel-level Dirac refit -- rerun the
  covariant CP1-kernel fit (archival scan_cp1_o2_yukawa chain) against
  refreshed charged-fermion masses + CKM, then propagate through 4b;
  also the GUT-scale sum-rule variant using the DYN-1 light-doublet
  composition (alpha, bar_alpha).

- **DYN-2 [Audit 3: thresholds and unification] (1-2 sessions).**
  Goal: two-loop MSSM running + one-loop GUT-scale threshold matching
  driven by `heavy_spectrum.json`; output a (M_GUT, alpha_G, M_S)
  compatibility WINDOW with pulls, never a point claim.
  Steps: (1) b-coefficient generator summing one-loop contributions of
  spectrum entries active between scales; (2) two-loop runner (numpy RK
  integration; standard MSSM b_ij; SUSY threshold at M_S per DYN-0
  decision); (3) matching Delta_i = sum_R b_i^R ln(M_R/M_GUT)/(2 pi);
  chi^2 against alpha_em(M_Z), sin^2 theta_W, alpha_s with uncertainties;
  (4) perturbativity gate to M_Pl; benchmark-neighborhood scan of the
  DYN-1 couplings.
  Gates: no-threshold limit reproduces textbook MSSM unification;
  flags: `matching_one_loop_only`, `no_unique_scale_claimed`.
  Assets: `code/two_loop_spectrum_fit.py`,
  `code/scan_mediator_threshold_rge.py`, `code/scan_thresholds_item5.py`
  (HEAD).
  Depends: DYN-1.
  Status: **done** (2026-07-05; `code/audit3_dyn2_thresholds_unification.py`,
  15/15 checks; ledgers `output/audit3/dyn2_thresholds_unification.{json,md}`;
  `.gitignore` whitelist extended for `output/audit3/`).
  Built: (i) b-generator with derivation-level universality gates -- the
  chiral census S-sum is (127,127,127) = T(210)+2T(126)+T(10), the 45 sums
  to (8,8,8) = T(45), so b_SO(10) = 109 is i-independent; (ii) two-loop
  RK4 runner, self-tested, textbook no-threshold MSSM unification
  reproduced (M_12 = 1.04e16, alpha_G^-1 = 25.9, alpha_3 mismatch -0.13%);
  (iii) one-loop supermultiplet matching at mu* = M_X (AG's lightest
  p-decay vector convention; eaten pairs at m_lambda, vector superfields
  at weight -3*pair; sign fixed by the EFT-decoupling derivation recorded
  in the ledger); AG's coefficient algebra (120, 60, (4,-9.6,5.6),
  0.0167 = 2/120) machine-verified; (iv) damped multi-start Newton solves;
  MC window over PDG-vintage input errors (flagged for the DYN-4 refresh).
  FINDINGS (disclosed, fed to DYN-8): the alpha_3 pull landscape over the
  benchmark slice runs from -56 sigma (x = 0.05) through zero between
  x = 0.12-0.14 to +54 sigma (x = 0.18); a compatibility point exists at
  x* = 0.15 (pull +0.8 sigma) BUT with M_X ~ 1.9e13 GeV and
  alpha_G^-1 ~ 52.5, so d = 6 proton decay is catastrophically fast: the
  raw benchmark slice (lambda = eta = 1, benchmark gamma, M_S = 3 TeV) is
  EXCLUDED by unification + proton decay combined.  Viability requires the
  (lambda, eta, gamma, complex xi) cancellation-region scan -- opened as
  the DYN-2b extension below.  Secondary: exact 3-parameter solve at
  x = 0.1 wants M_S ~ 5.6e7 GeV (outside [1,10] TeV, disclosed).
- **DYN-2b [extension]: viability scan.**  Scan (lambda, eta,
  gamma, bar_gamma, complex xi/x) for regions where BOTH the alpha_3 pull
  is small AND M_X stays >= 1e16 GeV (proton-safe), using the DYN-2
  machinery as-is; joint with the DYN-3 tau_p floor once available.
  This is AG's "cancellation regions" question made reproducible.
  Status: **DYN-2b/4c-alpha done** (2026-07-05;
  `code/audit3_dyn2b_rescue_scan.py`, 6/6 checks; ledgers
  `output/audit3/dyn2b_rescue_scan.{json,md}`).
  Joint (x, eta, M_S) scan, 10 x 6 x 4 grid, three survival functions
  (|pull| < 3; tau_d5 > Super-K; tau_d6 > Super-K -- the d=6 trap):
  **ZERO living points**.  131 cells solved, 109 have no physical
  exact-matching solution (dead by construction); binding census:
  pull kills 101, d=6 kills all 30 unification-compatible cells; best
  margin -9.5 orders.  The gamma lever is BOUNDED OUT analytically
  (gamma never enters the E/X vector masses that set M_X; its h+t
  threshold weight (3.8, 3.0, 5.0) of 127 caps the indirect ln M_X
  shift at ~1.3 orders vs the 2-3 needed).  STRUCTURAL CONCLUSION: the
  spectrum SHAPE of the minimal 210+126+126b+10 slice family forces
  exact-matching unification down to M_X ~ 1e13-1e14 GeV, and d=6 then
  kills every unification-compatible cell independent of the mini-split
  d=5 rescue.  In-slice levers exhausted: x (DYN-2), zeta/Majorana
  (DYN-4), eta + M_S (this scan), gamma (bounded).  Surviving escapes
  are MODEL-LEVEL: (i) non-SUSY SO(10) with an intermediate Pati-Salam
  scale -- no higgsino d=5 at all, two-step unification raises M_X
  (opened below as DYN-9); (ii) extra Higgs content (120-plet; archival
  audits exist); (iii) proper split-spectrum running (refinement only;
  cannot lift M_X by 2-3 orders).
- **DYN-9 [new lane]: non-SUSY / intermediate-scale variant.**
  Rebuild DYN-2 for a non-supersymmetric SO(10) chain
  (SO(10) -> Pati-Salam or G_LR at M_I -> SM), where d=5 proton decay is
  absent (no higgsinos) and two-step running decouples M_I from M_X.
  Status: **DYN-9a done** (2026-07-05;
  `code/audit9_dyn9_nonsusy_intermediate.py`, 15/15 checks after the
  2026-07-05 seesaw-disclosure patch; ledgers
  `output/audit9/dyn9_nonsusy_intermediate.{json,md}`; whitelist extended
  for `output/audit9/`).
  Method: all one-loop b coefficients DERIVED from field content (ESH
  scalars, GUT-normalized abelian factors) and gated against hand values;
  SM segment two-loop; Newton per chain; anchored against the BDM
  two-loop results (PRD 80 015013): the G_LR solver lands at
  (9.43, 16.32, 45.6) vs BDM (9.5, 16.2, 45.5) -- one-loop fidelity is
  excellent.
  RESULTS (current bound tau(p -> e+ pi0) > 2.4e34 yr):
  | G_LR (via 45_H, 126b) | M_I = 1e9.4, M_X = 1e16.3 | tau ~ 1e36:
    **ALIVE** (above the 10-yr Hyper-K e+pi0 reach ~1e35; only partially
    probeable) | but seesaw/leptogenesis strained: archival M_R ceiling
    needs f ~ 1e6, and M_1 > M_I |
  | PS = 2L2R4C no D (the 210-COMPATIBLE chain) | M_I = 1e11.9,
    M_X = 1e15.7 | tau ~ 4.3e33: **MARGINAL** (5.6x below the bound =
    inside BDM's threshold-spread caveat) -- either already dead or
    sitting EXACTLY in the Hyper-K discovery window; M_1 < M_I holds
    (thermal N_1 production KINEMATICS only -- see the perturbativity
    gate below) |
  | PS + D (needs 54_H) | M_X = 1e15.0 | tau ~ 3e30: dead |
  | 2L1R4C | M_X = 1e14.4 | tau ~ 3e28: dead (BDM reproduced) |
  STRUCTURAL VERDICT: the rescue EXISTS -- d=5 (the DYN-3 killer) is
  SUSY-specific and absent here, and at least one chain clears the d=6
  bound.  The two live-ish options trade off cleanly: G_LR is safely
  alive but strains the seesaw and swaps the framework's 210_H for 45_H
  (conditional-input change; CORRECTED by DYN-9b-1: the (p, a, 0, 0)
  vev pattern of the 210 itself leaves exactly the left-right group
  unbroken with D-parity broken by the D-odd p, so the 45_H swap is
  OPTIONAL, not forced); the 210-compatible PS chain preserves the
  source choice but is pinned at the bound -- maximally falsifiable.
  The kinematic theorem core (16 with forced nu^c, N = 3,
  K_tr = Killing) is SUSY-agnostic and untouched; the 126bar (P_nu^c
  source) survives in every chain.
  PERTURBATIVITY GATE (added 2026-07-05, consistency sweep): full-chain
  f_needed = M_R,i^archival/M_I table -- G_LR [8.9, 1.2e4, 1.5e6],
  PS [0.029, 39, 4.8e3], PS+D [5.4e-4, 0.73, 89], 2L1R4C
  [0.066, 90, 1.1e4] -- NO chain that survives proton decay can host
  the archival M_R tower via renormalizable f v_R with f < 4pi (worst
  PS factor 379x above 4pi; PS+D comes closest but is proton-dead).
  Root cause is a SCALE CLASH: the archival contact card has
  M_* = 3.9e15 GeV while surviving chains have M_I <= 1e12.  Hence the
  archival zeta/M_* card is SUSY-slice-LOCAL and NOT transplantable
  (ledger flags archival_MR_tower_perturbatively_realizable_on_
  surviving_chains = false, archival_zeta_card_transplantable_to_
  nonsusy_chains = false); zeta's value being branch-local is CONSISTENT
  with the boundary theorem (zeta is an input, not a derived quantity),
  but the non-SUSY flavor refit is REQUIRED for the alive branch, not
  optional.  Escape routes for the tower itself (conditional, string
  lane, unpromoted): d=5 Majorana sources (16 16 16bar 16bar)/M_s or
  instanton-induced M_R decouple the ceiling from f v_R -- see the
  Route-D/DYN-6 pricing discipline.
  Remaining as **DYN-9b** (open, REQUIRED for the alive branch):
  non-SUSY scalar potential and spectrum of the 210/126 system (real
  thresholds replacing ESH); the non-SUSY flavor refit incl. zeta/M_*
  re-extraction at v_R ~ M_I (archival M_R was SUSY-convention and
  fails the perturbativity gate above); two-loop intermediate segment;
  Yukawa-sector viability (non-SUSY SO(10) fits need the 126bar + 10
  complex doublet structure); non-SUSY leptogenesis re-run (gravitino
  tension dissolves; non-SUSY loop function); non-SUSY messenger
  protection (holomorphy unavailable -- anomalous-U(1)/Stueckelberg
  selection rule is the string-motivated substitute, unpromoted).

- **DYN-3 [Audit 2: d=5 proton-decay pipeline] (2-3 sessions).**
  Goal: physical C_5L / C_5R Wilson tensors and channel widths.
  Steps: (1) triplet Yukawa maps Y_QQ, Y_QL, Y_UE, Y_UD from the DYN-4
  flavor basis and DYN-1 light-doublet composition; (2) contraction with
  the exported symbolic triplet-inverse entries S_i^j — the audit-2
  source-basis contract (already gated in HEAD) now fed with physical
  tensors; (3) running/dressing: C5 anomalous-dimension running
  M_GUT -> M_S, one-loop wino/higgsino dressing to C6, QCD running to
  2 GeV; (4) hadronic matrix elements as cited lattice inputs;
  (5) widths for p -> K+ nubar (and channel table), compared to
  Super-K bounds and Hyper-K reach; emit the kill-criterion entry for
  DYN-8.
  Gates: deterministic random-tensor contract regression; decoupling
  limit (M_T -> infinity gives vanishing widths); flags:
  `hadronic_inputs_cited_not_computed`, `dressing_one_loop_only`.
  Assets: `code/audit2_source_basis_wilson_builder.py`,
  `code/scan_dimension5_wilson_tensors.py`,
  `code/audit_dressed_dimension5_channels.py`,
  `code/audit_mssm_mixing_d5_dressing.py`, `code/audit_full_knu_width.py`,
  `code/scan_proton_channel_bounds.py` (HEAD).
  Depends: DYN-1, DYN-2, DYN-4 (graph order: after Audit 3).
  Status: **done** (2026-07-05;
  `code/audit2_dyn3_proton_d5_kill_criterion.py`, 13/13 checks; ledgers
  `output/audit2/dyn3_proton_d5_kill_criterion.{json,md}`).
  Three stacked layers: (1) archival Knu calibration replayed and
  internally exact (S_T anchors scale as sqrt(bound ratio) to machine
  precision: Gamma ~ S_T^2; worst channel LLLL up-up-down p -> K+ nubar);
  (2) the HEAD audit-2 Wilson contract regenerated (deterministic gate)
  and fed PHYSICAL tensors -- DYN-1 det-tuned triplet inverse entries
  (max |S_i^j| = 1.99/m at x = 0.1) and the archival Yukawa sectors;
  (3) the HARD KILL CRITERION via the robust scaling law
  tau = SK_bound (M_T_eff/1e17 GeV)^2 with the minimal-SUSY-GUT
  literature anchor: at the DYN-2 compatibility point x* = 0.15,
  M_T_eff = 1.45e13 GeV gives tau(p -> K+ nubar) ~ 1.2e26 yr --
  **7.7 orders below Super-K's 5.9e33 yr** (arXiv:1408.1195) -- and
  across the ENTIRE solved x-scan the unification-vs-proton gap in
  m_scale is 2.7e2 to 1.2e6.  COMBINED EXCLUSION of the raw benchmark
  slice (lambda = eta = 1, benchmark gamma, M_S = 3 TeV) recorded as the
  DYN-8 hard entry.  Escape levers quantified: mini-split sfermions
  (3-10 orders in tau), triplet retuning via the DYN-2b/4c joint scan.
  Refinement left open as **DYN-3b**: flavor rotations to the physical
  basis (DYN-4c), explicit one-loop dressing replacing the literature
  anchor, channel table beyond the worst channel.

- **DYN-5 [hidden-messenger dynamics at one loop] (1 session).**
  Goal: promote the tree-level matching-silence statement to audited
  one-loop order: delta Z_N = |lambda|^2/(16 pi^2) ln(M_*/M_X) at
  benchmark (the companion note's 8.259697e-4 anchor), apply the
  documented replay shifts to Y_nuD and M_R, verify Dirac/triplet-channel
  silence at one loop under the R-selection rule, and state the
  R-breaking spurion conditions to the audited order.
  Gates: caveat numbers reproduced; replay stability vs the DYN-4
  windows reported.
  Assets: `code/verify_hidden_zeta_origin.py`,
  `code/verify_routeB_mediator_grading.py`,
  `code/scan_mediator_r_window.py` (HEAD).
  Depends: DYN-4.
  Status: **done** (2026-07-05; `code/audit5_dyn5_messenger_one_loop.py`,
  21/21 checks; ledgers `output/audit5/dyn5_messenger_one_loop.{json,md}`;
  whitelist extended for `output/audit5/`; requires the archival closure
  card, same sha256 as DYN-4a).
  FINDINGS (the tree theorem upgrades cleanly, with two sharpenings and
  one disclosed non-silence):
  (1) delta Z_N is EXACTLY family-universal: the portal coupling is
  lambda * 1_3 and the messenger mass matrix M_* K_tr^-1 has all
  singular values sqrt(3) M_* (from K_tr^2 = I/3), so the leading log is
  deltaZ * identity; the caveat anchor 8.259697e-4 is the per-e-fold
  coefficient (reproduced to 1e-10), and the STRICT benchmark interval
  is ln(M_*/(sqrt(3) M_*)) = -ln(3)/2, giving |deltaZ| = 4.54e-4.
  (2) Light-sector silence is STRUCTURAL, not perturbative: the seesaw
  m_nu = -mD M_R^-1 mD^T is invariant under mD -> mD A,
  M_R -> A^T M_R A for ANY invertible A (machine-zero, universal AND
  non-universal delta Z_N) -- canonical-normalization corrections cancel
  identically in every light observable.  The paper's LINEARIZED replay
  formulas have a quantified validity range: the O(deltaZ^2) residual
  crosses the loose contact-sensitivity window at deltaZ ~ 1.1e-2
  (ln M ratio ~ 14) and the REFRESHED window at deltaZ ~ 2.0e-3
  (ln M ratio ~ 2.4 -- barely one decade of messenger interval); beyond
  that use the exact A-form.
  (3) The HEAVY sector is NOT silent (disclosed): M_R singular values,
  M_1 (leptogenesis input) and M_* rescale by (1 - deltaZ); worst
  scanned shift = 3.6% of the DYN-4b M_1 window.  The normalized zeta
  extraction is EXACTLY invariant (zeta anchor digits are deltaZ-safe).
  (4) R-selection at one loop: all X-decorated Dirac/triplet/Majorana
  operators have R = 2 + m; the portal forces R(X) = 1; silence fails
  only at the measure-zero assignment R(16) = 2; non-renormalization
  makes one-loop effects Kahler-only and the nu^c leg appears in no
  visible Yukawa or d=5 channel, so delta Y_{u,d,e} = delta C_5L =
  delta C_5R = 0 at one loop and the DYN-2 threshold vector is
  untouched.  Two-loop leakage estimated at 1.5e-4 in alpha^-1 units
  (negligible vs DYN-2 sensitivity; two-loop silence NOT claimed).
  (5) Spurion conditions (regeneration iff k r_s = -m): positive-R
  spurions safe at the holomorphic level; even-R spurions
  (constant-W/gravitino type) preserve a residual (-1)^R parity so all
  ODD decorations stay absent to all orders; odd-R spurions are the
  dangerous case.  Order-estimate ceilings vs the DYN-4 windows:
  eps_odd < 1.4e-2 (loose) / 4.4e-4 (refreshed).
  Boundary: zeta NOT derived; R symmetry NOT claimed anomaly free; no
  microscopic messenger completion constructed; ceilings are order
  estimates.

- **DYN-6 [zeta-stabilization axioms: new-axiom lane, optional].**
  Constraint from Theorem VII.4: any mechanism fixing zeta is a NEW
  conditional axiom; this lane prices axioms, it does not derive zeta.
  For each candidate — (a) axion-volume modulus potential (string
  sandbox continuation), (b) hidden radial/phase potential with
  radiative minimization (clockwork-Hessian machinery in HEAD),
  (c) discrete anchors (already disfavored by the Z_178 diagnostic) —
  produce: precise axiom statement, parameter count, threshold/UV audit
  cost, and at least one prediction beyond zeta itself (otherwise it is
  re-parametrization).  Promotion bar identical to the string sandbox's.
  Assets: `code/audit_clockwork_radial_driver_hessian.py`,
  `code/audit_majorana_monomial_clockwork_origin.py`,
  `code/audit_constrained_clockwork_source_hessian.py` (HEAD).

- **DYN-7 [cosmology lane: leptogenesis as an independent constraint on
  arg zeta] (1-2 sessions, after DYN-4).**
  The contact phase is the extra CP phase of M_R.  Compute the CP
  asymmetry epsilon_1 (Davidson-Ibarra-type estimate) from the DYN-4
  (Y_nuD, M_R) output, an analytic washout/efficiency estimate, and an
  order-of-magnitude eta_B; scan arg zeta over the DYN-4 posterior and
  report sign/magnitude compatibility with eta_B ~ 6e-10.  This does not
  derive zeta (it is a new empirical anchor, consistent with Theorem
  VII.4); it is the first observational constraint on the otherwise-free
  phase, and a kill criterion if incompatible across the whole window.
  Gates: order-of-magnitude flags explicit; sign check explicit.
  Status: **done** (2026-07-05; `code/audit7_dyn7_leptogenesis_argzeta.py`,
  13/13 checks; ledgers `output/audit7/dyn7_leptogenesis_argzeta.{json,md}`;
  whitelist extended for `output/audit7/`).
  Chain: Takagi (self-tested) -> h = Y_nu V -> epsilon_1 with the SUSY
  loop function (Davidson-Ibarra respected point-wise) -> BDP-type
  washout -> eta_B = -0.96e-2 eps1 kappa (sign convention explicit).
  FINDINGS (honest form): on the archival Dirac slice thermal unflavored
  leptogenesis is MARGINAL -- central chains give the WRONG SIGN
  (antimatter) and |eta_B| ~ 1e-11 (boost x50-68 needed); bottleneck =
  strong washout (m_tilde ~ 0.05 eV, kappa ~ 6e-3) times phase-suppressed
  epsilon_1 (0.4% of the DI ceiling); posterior P(success) = 0.4% (tails
  only); BUT the boost sensitivity shows the slice is only ONE order of
  magnitude from typical viability (x3 -> P = 34%, x5 -> 41%, x10 -> 45%),
  within reach of flavored-leptogenesis enhancements or a modified Dirac
  sector (DYN-4c).  arg zeta itself is NOT the lever on this slice
  (success nearly uncorrelated with the scanned variables; the
  14-sample conditioned window is noisy and flagged).  Soft kill-criterion
  entry recorded for DYN-8: thermal unflavored leptogenesis on this slice
  underproduces eta_B by ~1.5 orders and prefers the wrong sign at the
  central phase convention; escape routes = flavor effects (O(few)),
  Dirac-sector refit, or non-thermal production.  Gravitino/reheating
  tension disclosed, not resolved.

- **DYN-8 [falsifiability section] (1-2 sessions; collects all lanes;
  CLOSES the main line).**
  Add a Falsifiability section to the reconstruction paper + a ledger.
  UPDATED SCOPE (2026-07-05, after the consistency sweep): the section
  opens with a BRANCH MAP -- (i) kinematic theorem core (SUSY-agnostic,
  machine-verified), (ii) SUSY minimal slice EXCLUDED (DYN-2 x DYN-3
  combined + DYN-2b structural obstruction), (iii) non-SUSY alive branch
  (G_LR / 210-compatible PS; pending the REQUIRED DYN-9b refit),
  (iv) string-conditional lanes D3-D5 (priced, unpromoted).  Kill bank,
  each entry branch-tagged: fourth chiral family (kills N <= 3); Dirac
  neutrinos / absence of 0nubb (removes the canonical-contact target);
  the PS-chain Hyper-K window (tau ~ 4.3e33 vs bound 2.4e34: dead or
  discoverable); G_LR tau ~ 1e36 (partial Hyper-K reach); DYN-7
  leptogenesis soft kill (SUSY-scoped; gravitino tension dissolves on
  the non-SUSY branch); DYN-5 spurion ceilings; D4's instanton-action
  bound S' >= 7.7 (if D4 closes); D5's M_SS window (or its emptiness);
  Z_178 retirement on the refreshed card (data supersedes anchors,
  strengthening the continuous-phase reading).  Paper touch-ups: main
  paper status ledger deferred -> completed-with-findings; Z_178
  re-scope.  Deliverable tex stays SELF-CONTAINED (repo codenames only
  in the terminology remark + literal file paths).
  Status: **done** (2026-07-05;
  `code/audit8_dyn8_falsifiability_collection.py`, 22/22 checks;
  ledgers `output/audit8/dyn8_falsifiability_collection.{json,md}`;
  whitelist extended for `output/audit8/`).  THE MAIN LINE
  0 -> (1 || 4a) -> 2 -> 3 -> 5 -> 8 IS CLOSED.
  Delivered: (1) collector audit re-verifying every quoted number
  against its source ledger -- 12 ledgers, 210 upstream checks green,
  no ledger claims a derived zeta; (2) reconstruction paper gains the
  "Dynamics Audits and Falsifiability" section (branch map with four
  layers + kill bank K1-K11, all branch-tagged, self-contained
  language; terminology remark extended for the DYN/D file aliases;
  status ledger updated executed-audits + still-open; manifest extended
  with the collector); (3) main paper status ledger updated
  (deferred -> executed July 2026 with findings) and the Z_178
  observation carries the retirement note (misses the refreshed card by
  76x the loose window).  Both PDFs rebuilt clean.
  Kill bank: K1 fourth family (core); K2 Dirac neutrinos (conditional
  contact target; unweighted prior-draw fraction cf>0.01 = 1.000);
  K3 SUSY-slice exclusion (FIRED);
  K4 PS Hyper-K window (dead-or-discoverable, the sharpest entry);
  K5 channel hierarchy (e+ pi0 dominance, G_LR 1e36); K6 leptogenesis
  unflavored regression (branch-tagged and unpromoted); K7 historical
  messenger-ceiling diagnostic (invalid pending an interacting DYN-5
  action); K8 fixed-archival-tower N2,N3 > M_I ordering only (the
  instanton mechanism permits but does not require it); K9 historical-invalid
  Delta S = 4.66 / 3.18 arithmetic inherited from DYN-5; K10 preliminary
  bridge toy window / PS-grid obstruction (unpromoted); K11 Z_178 retirement.

Suggested sequencing: DYN-0 -> (DYN-1 || DYN-4) -> DYN-2 -> DYN-3 ->
DYN-5 -> DYN-8, with DYN-6/DYN-7 optional lanes after DYN-4.  Statuses
start `open`; update per session.

## F. String-conditional program (planned 2026-07-05; items D3-D5,
## then DYN-8 close-out and the staged DYN-9b)

> **Historical execution log.**  The DYN-9b `closed`, K8 `requirement`,
> unflavored posterior, reheating, and 210/45 statements below are superseded
> by the 2026-07-14 blocker update and claim registry.  RE-SC3/4/5 are present
> but unpromoted; the global fit, flavored kinetics, messenger rederivation,
> and uncertainty envelopes remain open.

Discipline: every item below is a CONDITIONAL string interpretation
under the Route-D promotion bar (precise statement, assumptions-vs-
derived list, reproducible script, non-overclaiming impact note) and
the DYN-6 pricing language (axioms are PRICED, not derived; each must
state at least one consequence beyond the quantity it sources, or
disclose "none").  Nothing here derives zeta (boundary theorem).
Scripts live in `route_d/code/verify_dN_*.py` with flat ledgers
`route_d/output/dN_*.{json,md}` (the string sandbox, physically
separated from the dynamics lanes).  Execution order:
D3 -> D4 -> D5 -> DYN-8 -> DYN-9b.

- **D3 [instanton / d=5 Majorana-source pricing card] (1 session).**
  Goal: price the two string escapes for the archival M_R tower that
  the DYN-9 perturbativity gate closed for renormalizable f v_R.
  Steps: (1) the d=5-operator source (16 16 16bar 16bar)/M_s with
  v_R ~ M_I gives f'_needed = M_R,i M_s / v_R^2 -- EXPECTED WORSE than
  the renormalizable case (extra suppression v_R/M_s; at M_s = 2e16,
  v_R = 8e11 the top entry needs f' ~ 1e8): compute and GATE the
  negative result, closing this escape at M_I-scale v_R and leaving
  instantons as the only string escape for the tower.
  (2) instanton source M_R,i ~ M_s e^{-S_i}: required actions
  S_i = ln(M_s/M_R,i) over an M_s grid (M_X of both surviving chains,
  2e16, 1e17, 1e18); at M_s = 2e16 the archival tower needs
  S = [13.6, 6.4, 1.6] and the D2 contact action is
  S_zeta = -ln|zeta| = 2.04 -- record the top-of-tower/contact
  same-ballpark coincidence explicitly (diagnostic, not evidence).
  (3) rank bookkeeping: a single instanton gives generically rank-1
  M_R (zero-mode counting), so the rank-3 tower needs >= 3 distinct
  cycles; the family texture becomes a NEW conditional input.
  (4) axiom price card: statement, parameter count (3 actions +
  phases + texture), what it buys (tower decoupled from the f < 4pi
  ceiling; B-L broken by the instanton itself, independent of M_I),
  what it does NOT buy (zeta's value; the texture).
  Gates: negative d=5-escape gate; S_i table; rank statement; price
  card with beyond-zeta consequence or explicit "none".
  Depends: DYN-9 (patched ledger).
  Status: **done** (2026-07-05;
  `route_d/code/verify_d3_instanton_majorana_pricing.py`, 15/15 checks;
  ledgers `route_d/output/d3_instanton_majorana_pricing.{json,md}`).
  FINDINGS: (1) d=5-operator escape CLOSED as expected -- f' = f *
  (M_s/v_R) strictly worse; PS chain at M_s = 2e16 needs f'_top ~
  1.2e8; renormalizable clash factor v_R_min/M_I = 379 can only grow.
  (2) Instanton escape priced: S = [1.63, 6.43, 13.64] at M_s = 2e16;
  D2 contact action S_zeta = 2.037 is same-ballpark as the tower top
  (recorded as diagnostic, NOT evidence); weak coupling (S_top >= 1)
  needs M_s >= 1.07e16 -- fine on G_LR (M_X = 2.1e16), but on the PS
  chain the string scale must sit ABOVE M_X = 5.4e15 (disclosed).
  (3) Takagi bookkeeping: the tower is exactly three rank-1 terms =>
  >= 3 distinct cycles; instanton data = 3 actions + 9 texture = 12 =
  full Majorana data -- the axiom buys SCALE (f < 4 pi evasion, M_R
  decoupled from M_I), ZERO texture compression.  (4) Beyond-zeta
  consequence: N_2, N_3 heavier than the intermediate gauge scale
  (PS: 39x, 4763x M_I) -- impossible renormalizably, falsifiable if
  the Z'/W_R scale and heavy-N masses are both measured.  Boundary:
  no zero modes / global divisor / moduli; texture is a NEW
  conditional input; unpromoted.

- **D4 [Stueckelberg / anomalous-U(1) protection card] (1 session).**
  Goal: the non-SUSY substitute for the DYN-5 R-selection rule -- a
  perturbatively exact anomalous U(1)_A (GS-cancelled, Stueckelberg-
  massive) charging the messenger sector; violations only from
  instantons e^{-S'}.
  Steps: (1) charge bookkeeping by brute-force enumeration (reuse the
  DYN-5 machinery): portal XN allowed, all Dirac/triplet decorations
  forbidden, AND the obstruction the SUSY version never faced -- the
  X mass term K_tr^-1(X,X) carries charge 2 q_A(X) != 0, so it must be
  sourced by a charged vev or by the SAME instanton class as the D3
  Majorana source (this couples D3+D4 into one bookkeeping; enumerate
  which assignments close, or record the no-go).
  (2) strength gate: perturbative exactness upgrades silence to ALL
  loop orders (stronger than the SUSY one-loop statement); instanton
  contamination must satisfy the DYN-5 ceilings:
  S' >= ln(1/eps) = 7.73 (refreshed window) / 4.27 (loose); compare
  with the D3 actions -- disclose whether the tower-generating and
  protection-violating instantons can be the same class consistently.
  (3) NOT claimed: anomaly inflow details, global embedding, moduli
  stabilization.
  Gates: a closing charge assignment or an explicit no-go; S' ceilings
  computed; D3-D4 shared-instanton consistency statement.
  Depends: D3, DYN-5.
  Status: **done** (2026-07-05;
  `route_d/code/verify_d4_stueckelberg_protection.py`, 14/14 checks;
  ledgers `route_d/output/d4_stueckelberg_protection.{json,md}`).
  FINDINGS (the enumeration SHARPENED the planned card):
  (1) GENERALIZED ABELIAN NO-GO: if the Dirac Yukawa and the portal are
  charge-allowed then q(X L H) = q(X X) IDENTICALLY, for any number of
  abelian factors (exhaustive k = 1, 2; B-L instance: the portal forces
  B-L(X) = -1, both operators sit at -2, so matter parity does not
  separate them).  The roadmap's "closing charge assignment" does NOT
  exist -- the no-go branch of the gate fired.
  (2) Non-SUSY Lorentz bookkeeping: all m = 1 SUSY-style decorations
  have ODD fermion count => forbidden outright in components; m = 2
  first arises at dimension 7.  The only renormalizable contamination
  is the direct X L Htilde -- exactly the operator the no-go ties to
  the X mass.
  (3) Historical-invalid numerical diagnostic (NOT a closure condition):
  under the conditional unit-prefactor XX benchmark at M_s = 2e16,
  eps_XX = sqrt(3) M_*/M_s = 0.34 and S_XX = 1.08.  The old DYN-5
  ceilings then replay Delta S = S'(XLH) - S(XX) = 4.66 / 3.18 and
  S'(XLH) = 5.74 / 4.26, but these are not physics bounds because the
  underlying DYN-5 messenger action is invalid.  Only the exact additive-
  charge no-go survives.  D3 consistency: the XX and NN operators carry
  opposite B-L charges (-2 and +2); treating their compensating sources as
  conjugate instanton orientations with independent actions is an extra
  ansatz, not a consequence of charge equality.  The K_tr texture of the
  X mass is likewise an explicit conditional input (the SUSY theorem
  postulated it).
  (4) Portability: DYN-5's exact light-sector silence is pure linear
  algebra and survives non-SUSY unchanged.
  (5) Beyond-window low-energy consequence: honest NONE (X at
  sqrt(3) M_* ~ 6.8e15 GeV); the card's testable content is internal
  zero-mode counting in a concrete embedding.
  Boundary: a valid numerical selectivity bound remains OPEN pending an
  interacting messenger action and explicit zero-mode calculation; no
  anomaly inflow / global embedding; zeta NOT derived; unpromoted.

- **D5 [high-scale SUSY-breaking bridge scan] (1-2 sessions).**
  Goal: the interpolating family between the EXCLUDED SUSY slice and
  the DYN-9 non-SUSY chains: SUSY broken at M_SS (non-SUSY below, SUSY
  above) -- a UV-SUSY toy interpolation.  The inequality M_SS < M_*
  is only a necessary matching-scale ordering; it does not validate
  the interacting Route-B/DYN-5 messenger rejected by DYN-5V.
  Steps: (1) two-segment running per surviving chain (G_LR and the
  210-compatible PS): DYN-9 non-SUSY b's below M_SS, SUSY versions
  above; solve (M_I, M_X, alpha_G) as functions of M_SS on a grid
  M_SS in [1e4, M_X] GeV.
  (2) proton ledger: tau_d6(M_X(M_SS)); d=5 REAPPEARS above M_SS --
  compute the mini-split scaling of the DYN-3 estimate with sfermions
  at M_SS (the "3-10 orders" lever becomes a computed curve).
  (3) window intersection: alive(tau_d6) AND tau_d5(M_SS) > Super-K
  AND M_SS < M_* (holomorphy at matching); disclose both orderings of
  M_SS vs M_I.  Output: the surviving [M_SS] window per chain, or
  EMPTY.
  (4) verdict: non-empty window => a falsifiability entry (tau_p
  pinned into a band); empty => the bridge is DEAD and the alive
  branch definitively loses Route-B holomorphy (making D4
  load-bearing for any messenger story there).
  Gates: DYN-9 limit reproduced as M_SS -> M_X; the SUSY-slice
  exclusion reproduced as M_SS -> TeV; explicit window verdict.
  Depends: DYN-9, DYN-3, DYN-5.
  Status: **done** (2026-07-05;
  `route_d/code/verify_d5_susy_breaking_bridge.py`, 18/18 checks;
  ledgers `route_d/output/d5_susy_breaking_bridge.{json,md}`).
  FINDINGS: SUSY intermediate b's DERIVED with the anomaly-forced
  content doubling (U(1)_{B-L} cubic anomaly of Delta_R = 24 forces
  Delta_R-bar; SUSY PS doubles Sigma_R): b2R jumps -7/3 -> +5 (LR) and
  11/3 -> +41 (PS).  Endpoint gates: DYN-9 ledger reproduced to 1e-4;
  DYN-3 calibration 1.24e26 reproduced.
  | **G_LR toy-grid window (unpromoted)** | log10 M_SS in [10.30, 15.55] |
    bound BELOW by no-physical-solution (solvability, near the DYN-9
    M_I = 1e9.4 -- NOT by d=5: tau_d5 ~ 1e46 at the lower edge, the
    mini-split lever never binds), bound ABOVE by M_SS < M_*
    (holomorphy) | M_X stays ~1e16.3, tau_d6 in [6.2e35, 9.3e35] |
    all Case A (M_I <= M_SS); M_SS < M_* is necessary but does not
    restore the invalid DYN-5 action |
  | **PS toy-grid window empty under stated assumptions** | only 18/263 points physical,
    all with a vanishing SUSY segment (longest case-A segment ~ 0
    dex); b2R = +41 is associated with this scan result but no causal
    counterfactual has been performed; no sampled M_SS gives
    tau_d6 > 2.4e34 |
  | TeV end: NO physical bridge solution for either chain --
    stronger than d=5-dead; the DYN-2/3 exclusion structure
    reproduced (the catastrophic DYN-3 magnitude was slice-specific,
    disclosed) |
  READING: this preliminary grid contains a G_LR arithmetic window under
  one-loop/ESH/degenerate-spectrum assumptions and no PS window.  It is not
  a physical bridge-survival theorem and supplies no valid messenger action.
  Boundary: one-loop intermediate segments; ESH-minimal content both
  phases; M_T = M_X and degenerate soft spectrum assumed; above-M_X
  perturbativity not audited; unpromoted.

- **DYN-8 close-out** as specified in section E (branch map + branch-
  tagged kill bank incl. the D3-D5 outputs).  Status: **done**
  (2026-07-05; see the section-E DYN-8 entry for the full record).

- **DYN-9b [non-SUSY core re-derivation -- REQUIRED for the alive
  branch] (staged, 3+ sessions).**
  9b-1: non-SUSY scalar potential and spectrum of the 210+126bar(+10)
  system (superpotential-level blocks carry over from DYN-1, the
  potential does NOT); redo the DYN-9 thresholds with the real
  spectrum replacing ESH.
  Status 9b-1: **first stage done** (2026-07-05;
  `code/audit9_dyn9b1_nonsusy_vacuum_thresholds.py`, 13/13; ledgers
  `output/audit9/dyn9b1_nonsusy_vacuum_thresholds.{json,md}`).
  FINDINGS: (1) little-group map on the (p, a, omega, sigma) slice via
  zero-counting of the SUSY-INDEPENDENT gauge masses (the DYN-1b
  m_lambda formulas are |D_mu vev|^2 objects): seven patterns verified
  -- (p,0,0,0) -> PS (21, D broken), (0,a,0,0) and (p,a,0,0) -> the
  left-right group (15; p breaks D), (0,0,w,0) -> 3c2L1R1BL (13),
  (t,t,t,0) -> flipped SU(5) (25), (-t,-t,t,0) -> SU(5) (25), AG
  benchmark -> SM (12).  KEY: **the left-right chain is 210-REALIZABLE;
  the 45_H swap is optional** (D-parity assignments cited, CMP 1984).
  (2) Eaten-Goldstone integers: PS step eats E+X = 24, LR step
  J+E+X = 30; second stage G+J+F = 9 (PS) / G+F = 3 (LR) from the
  126bar's 60 real dof.  (3) ESH SENSITIVITY (survivors at M_I, worst
  case, masses NOT derived): PS verdict threshold-fragile in BOTH
  directions among physical solutions -- (15,1,1) or (10,2,2)+(10b,2,2)
  lift tau to ~1e40 (alive), (15,2,2) lands at 3.6e34 (just above the
  bound, inside the Hyper-K window), (15,3,1)/(6,2,2) sink it to
  7.6e29/1.0e32 (dead); (15,1,3) has NO physical solution (inverted
  M_I > M_X).  G_LR NOT unconditionally robust: (1,3,1,0) pulls tau to
  1.2e33 (dead), (3,2,2,2/3) and a second bidoublet to ~1.3e34
  (marginal-dead); (8,1,1,0) is trans-Planckian (flagged unphysical);
  (1,1,3,0) harmless.  READING: the survivor SPECTRUM decides both
  verdicts -- the K4/K5 falsifiability entries now carry an
  in-framework threshold-spread quantification, and 9b-1b (quartic
  Clebsch descent -> survivor masses) is the decisive open item.
  (4) Tree-level pseudo-Goldstone landmine disclosed with the one-loop
  stabilization literature anchor (BDM PRD 81 035015), not re-derived.
  Remaining as **9b-1b**: the 210+126bar quartic invariant descent on
  and off the singlet slice, survivor masses, tree-vs-one-loop
  stability of the (p,a) and (p) vacua.
  Status 9b-1b: **done** (2026-07-05;
  `code/audit9_dyn9b1b_210_quartic_descent.py`, 14/14; ledgers
  `output/audit9/dyn9b1b_210_quartic_descent.{json,md}`).  The 210 is
  built as an explicit rank-4 antisymmetric tensor; every step gates
  against the SUSY-side anchors (I2 -> 24(p^2+3a^2+6w^2),
  I3 -> 48(a^3+3pw^2+6aw^2) = the AG superpotential structure
  bit-level; the 45x45 vector mass matrix reproduces the DYN-9b-1
  little-group map AND the DYN-1b J/F/E/X mass ratios with census
  multiplicities).  MACHINE GROUP THEORY: the five polynomial quartic
  contractions span only FOUR dimensions (the cyclic pattern Q5
  degenerates with the 2-2 pattern Q2 -- an identity of the
  antisymmetric 4-tensor found by the rank test), and the Hodge
  epsilon quartic is a genuine independent FIFTH invariant whose slice
  descent vanishes IDENTICALLY (it moves only off-slice masses; its
  coupling is NOT scanned -- flagged).  Exact Hessians (Richardson on
  cubic gradients); GAUGE GATE: exactly 24/30 Goldstone zeros at the
  PS/(p,a) vacua for random couplings.
  FINAL TREE-LEVEL VERDICTS (250 random O(1) couplings in full-sum
  units, 4 scanned quartics):
  | PS (p-only) vacuum | tree-level local minimum on 16% of coupling
    space (39/250); on EVERY positive sample the 210 survivors sit at
    or above the gauge scale, so the threshold corrections vanish and
    tau(p -> e+ pi0) = 4.25e33 at all percentiles = the ESH/DYN-9
    verdict EXACTLY.  **K4 (Hyper-K window, dead-or-discoverable)
    survives the real-spectrum test unchanged**; the light-survivor
    rescue/kill scenarios of 9b-1 require coupling tuning, not
    generic couplings. |
  | LR (p,a) vacuum | tree-level local minimum on 0/249 samples --
    the classic non-SUSY tree-level obstruction BITES: within the
    scanned coupling space the 210-only left-right vacuum is NEVER a
    tree minimum.  DYNAMICAL REVERSAL of the 9b-1 kinematic point:
    the 210 CAN host the LR little group, but its tree potential
    rejects it; the left-right chain therefore NEEDS the 45_H route,
    or one-loop stabilization (BDM PRD 81 035015 type), or the
    UNSCANNED epsilon quartic (which contributes exactly and only to
    off-slice masses -- the one direction that could rescue tree
    positivity; open). |
  Slice phase-diagram sample: deepest bounded slice minimum is
  PS-type 26 / LR-type 4 / SU(5)-or-SM-type 120 out of 150 -- tree
  level prefers the enhanced-symmetry vacua, consistent with the
  landmine literature.  Boundary: 126bar quartics deferred (9b-1c);
  epsilon coupling unscanned; tree level only; coupling sampling in
  full-sum-invariant units; zeta NOT derived.
  Status 9b-1c: **done** (2026-07-05;
  `code/audit9_dyn9b1c_eps_and_126_quartics.py`, 11/11; ledgers
  `output/audit9/dyn9b1c_eps_and_126_quartics.{json,md}`).
  (A) THE EPSILON LEVER IS EXHAUSTED: the full Q6 gradient (three
  direct placements + the dual-slot adjoint) is gated; the exact
  epsilon Hessian (cubic-homogeneity identity) passes the strongest
  available check -- exactly 24/30 Goldstone zeros persist at
  lambda_eps != 0.  Sweeping lambda_eps over [-1.5, 1.5] on 200
  random coupling points: **0 rescues** (median best min-eigenvalue
  -9.6e2, deeply tachyonic); PS positivity never flips (0/200).
  VERDICT: the 210-only left-right vacuum is TREE-DEAD with no
  remaining tree-level lever -- the left-right chain requires the
  45_H route or one-loop stabilization, full stop.
  (B) THE 126bar SECTOR: embedding gates extended -- the mixed cubic
  descends EXACTLY to the AG eta-term |sigma|^2 (p + 3a - 6w) (this
  also fixed the e_w sign convention, invisible to the omega-even 210
  invariants; 9b-1b's flipped/SU(5) locus labels swap, dims
  unaffected), and the full vector-mass gate with sigma passes with
  scale = 0.500000, sigma-norm c = 1.000000, J/F/E cross-checked.
  Duality projection: the sigma direction lies entirely in the 126bar
  half (norm 1 to 1e-12); the mixed cubic + four mixed quartics split
  the 126bar into exactly (6,1,1)+(10,1,3)+(10bar,3,1)+(15,2,2), the
  D-parity-odd CUBIC being what splits the two 30s; sigma is
  machine-located in the SU(2)_R-charged 30 (the SU(2) naming is
  ANCHORED by sigma: nu^c is right-handed).
  (C) K4 RETEST WITH THE 126bar REMNANTS: over random (eta + 4 mixed)
  couplings with the extended-survival tuning and positive remnants
  (45/150 viable): tau = 4.25e33 at ALL percentiles -- the ESH/DYN-9
  verdict again.  **K4 (the Hyper-K window) has now survived the real
  210 spectrum (9b-1b) AND the real 126bar spectrum (9b-1c) at tree
  level.**  Boundary: (Sigma Sigma*)^2 self-quartics enter only the
  M_I-scale stationarity, not the M_X masses (deferred with that
  argument); epsilon-mixed patterns not enumerated; generating set of
  mixed quartics only; tree level; zeta NOT derived.
  The 9b-1 series (kinematics -> 210 potential -> epsilon + 126bar)
  is CLOSED.
  9b-2: non-SUSY flavor refit at v_R ~ M_I: a perturbative M_R tower,
  zeta/M_* RE-EXTRACTION on the alive branch (branch-local zeta is
  consistent with the boundary theorem), contact-essentiality retest.
  Status 9b-2: **done** (2026-07-05;
  `code/audit9_dyn9b2_nonsusy_flavor_refit.py`, 14/14; ledgers
  `output/audit9/dyn9b2_nonsusy_flavor_refit.{json,md}`).  Refit-level
  structural results from the light-neutrino side (a FULL non-SUSY
  Yukawa fit is NOT performed -- the Dirac kernel shape is a retained
  conditional input, disclosed):
  (1) FIXED-KERNEL CEILING: the corrected one-loop replay gives
  y_t(M_X)=0.44116481.  The relations Y_nu=h-3f and Y_u=h+f do not force
  Y_nu~Y_u (h=3f gives Y_nu=0, Y_u=4f).  The actual archival-kernel
  suppressions are 19.5x/342.2x (PS/G_LR); 9.6x/169.2x applies only
  under an additional top-like ansatz.
  (2) TYPE-II ESCAPE FAILS (order estimate, generous lambda ~ 1 and
  f = 4 pi): deficit 3.7 (PS) / 7.3 (G_LR) orders with the
  Delta_L-type block at the gauge scale (9b-1c placement).
  (3) **POSITIVE-REAL RESCALING THEOREM**: for uniform y>0, zeta, the
  projective contact direction, and contact fraction are invariant and
  M_* -> y^2 M_*.  For complex y, M_* -> |y|^2 M_* and
  zeta -> (y^2/|y|^2) zeta, so its phase is not invariant.  Re-extracted windows at the
  ceiling: M_*' <= 1.0e13 (PS) / 3.4e10 (G_LR).  Contact essentiality
  (P(cf > 0.01) = 1.000) carries over verbatim (conditional on the
  kernel shape only).
  (4) CONDITIONAL SOURCE COMPARISON: the rescaled renormalizable
  source gives M_1' = 7e7 (PS) / 2e5 (G_LR) GeV -- below the
  Davidson-Ibarra floor (~5e8, order estimate); the scale-decoupled
  (instanton-type) benchmark keeps the archival M_1 = 2.8e10 above the
  adopted floor.  No source is selected without a global flavor fit and
  flavored thermal kinetics.  The N_2,N_3>M_I ordering belongs only to the
  fixed archival tower and is not a generic instanton requirement.
  9b-3: non-SUSY leptogenesis re-run on the 9b-2 tower (non-SUSY loop
  function; the gravitino constraint dissolves).
  Status 9b-3: **done** (2026-07-06;
  `code/audit9_dyn9b3_nonsusy_leptogenesis.py`, 11/11; ledgers
  `output/audit9/dyn9b3_nonsusy_leptogenesis.{json,md}`).  The DYN-7
  pipeline rerun with exactly three replacements: f_SUSY -> f_SM
  (verified: hierarchical ratio exactly 1/2; log1p ESSENTIAL -- at
  tower hierarchies x ~ 1e10 the SM function suffers catastrophic
  float cancellation with a naive log, which inflated eps1 by ~100x
  in the first pass and was caught by the Davidson-Ibarra gate);
  v_u = 100 -> v = 174 in the loop Yukawa (m_D fixed by data;
  m_tilde and kappa v-independent, match DYN-7 exactly); gravitino
  ceiling REMOVED.
  FINDINGS: net suppression eps1(non-SUSY)/eps1(DYN-7) = 0.165 (~6x
  harder); the unweighted prior-draw regression has 0/4000 target-band
  hits (not a posterior or success probability; DYN-7 diagnostic: 0.35%);
  median |eta_B| = 1.45e-11 vs observed 6.1e-10; boost table: x10 ->
  0.13, x30 -> 0.40, x60 -> 0.45, x100 -> 0.34.  These are regression
  target-band hit fractions, not viability probabilities.  The
  gravitino-specific ceiling is absent on this branch, but reheating
  history and initial abundance remain unspecified.
  K6 branch-tagged verdict for the 9b-4 refresh: a SOFT constraint on
  the source scenario, not a kill -- viability requires flavored +
  O(10) effects, non-thermal production, or a resonant D3-type tower
  (the instanton actions are free inputs, so near-degeneracy is
  available at texture cost).  Boundary: unflavored N_1-dominated
  thermal estimate; archival tower under the D3-type source is a
  conditional benchmark (D3 unpromoted); zeta NOT derived.
  9b-4: DYN-8 ledger refresh with the 9b numbers.
  Status 9b-4: **done** (2026-07-06; the DYN-8 collector
  `code/audit8_dyn8_falsifiability_collection.py` extended with the
  refresh section, now 30/30 checks against 21 ledgers / 304 upstream
  checks; ledgers `output/audit8/` regenerated).  Refresh entries:
  K4r (TRIPLE-TESTED: ESH one-loop -> computed 210 spectrum ->
  computed 126bar spectrum, tau = 4.25e33 at all percentiles in every
  treatment); K5r (45_H narrative corrected: LR 210-realizable
  kinematically, tree-dead dynamically, epsilon lever exhausted ->
  adjoint route or one-loop stabilization required); K6r (REPLACED by
  the branch-tagged 9b-3 verdict: ~6x harder, P = 0/4000 unflavored,
  gravitino ceiling dissolved -- soft constraint on the source
  scenario); K8r (fixed-tower conditional diagnostic; no source selection,
  D3 unpromoted); K11r (positive-real invariance plus complex phase
  covariance, restricted to the fixed-kernel family).
  Both papers synchronized: the reconstruction paper's branch map now
  records the executed vacuum/source/contact stages, K4/K5/K6/K8 and
  the remark carry the refresh, and the main paper's status note
  records the executed re-derivation.  Both PDFs rebuilt clean; the
  collector GATES on all of it.
  Depends: DYN-9 (done), DYN-8 (structure); D3/D4 as source
  dependencies (D3 now load-bearing via K8r).
  **Historical close-out label withdrawn:** the ledgers execute
  mechanically, but global flavor, flavored kinetics, messenger/source,
  and threshold uncertainties keep the physics lane open.

Backlog (optional, unchanged): DYN-6 (axiom-pricing generalization),
DYN-3b (physical-basis rotations, explicit dressing).
(DYN-4c is promoted out of the backlog into SUB-B below.)

## G. Submission program (planned 2026-07-06)

> **Historical submission log with current correction.**  The theorem Letter
> now proves only \(N_{\rm fam}\le3\) from original H3.  Equality,
> \(g=0\), the \(\mathcal O(2)\) carrier/two-center branch, and the Killing
> contact are conditional on the explicit H3+ Killing-contact axiom.  The
> corrected claim matrix is `SUBMISSION.md`.

Strategy: split into three submission units rather than one monolith;
two hard prerequisites gate everything.  Order:
SUB-A -> SUB-B -> SUB-1 -> (SUB-2 | SUB-3).

- **SUB-A [prior-art due diligence, round 2] (1-2 sessions).**
  Adversarial literature search BEFORE any submission, recorded in
  `route_E/SUBMISSION.md`:
  (i) the three-point assembly -- family symmetry GENERATED as the
  automorphism algebra of its own carrier (vs imposed); N_fam <= 3 a
  priori with equality by the Cartan criterion; Majorana contact
  direction = the Killing form; second-transvectant/Killing readings
  of neutrino textures;
  (ii) the non-SUSY 210 vacuum literature (He--Meljanac, Buccella et
  al., and successors): whether the tree-level status of the
  left-right little group for the 210, the quartic-invariant count,
  and the epsilon-invariant slice behavior are already known --
  DECIDES whether SUB-2 is submittable and how strongly SUB-1's
  claims may be worded;
  (iii) boundary-theorem genre precedents (theories proving their own
  inputs underivable).
  Gates: every candidate-novelty claim gets a verdict
  (clear / partial-overlap / covered) with citations.
  Status: **done** (2026-07-06; verdicts in `route_E/SUBMISSION.md`).
  HEADLINES: the three-point assembly survives (Killing-form contact
  and the transvectant reading CLEAR; automorphism-generated family
  and the index bound ADJACENT with a required
  distinguish-paragraph); BUT the 210 tree-vacuum claims hit a
  CONFLICT -- He-Meljanac PRD 33 2695 (1986) and the 1989 stability
  papers find left-right minima in finite parameter ranges, clashing
  with 9b-1b's fixed-vev-ratio scan.  Corrective audit **SUB-B0**
  inserted (mandatory): free the LR vev ratio, rescan positivity,
  reproduce-or-localize.  K5r downgraded to ratio-conditional pending
  SUB-B0.  Bonus candidate-novelty found: the literature's
  "most general" 210 potential uses 4 quartic invariants (the
  parity-even set); our machine-independent FIFTH (Hodge-epsilon,
  parity-odd, slice-vanishing, stability-relevant) may be absent
  there -- full-text verification required before wording any claim.
  SUB-1 is UNBLOCKED (no conflict touches it); SUB-3 must wait for
  the K5r correction.

- **SUB-B0 [corrective audit: the left-right vev ratio freed]
  (1 session; MANDATORY, inserted by SUB-A).**
  Rerun the 9b-1b/1c left-right analysis with the vev ratio a/p as a
  scanned variable (stationarity solved per ratio), epsilon coupling
  included, tree positivity over (ratio x couplings): reproduce the
  He-Meljanac stable left-right regions or localize the disagreement
  (their 4-invariant potential vs our 5-invariant one is a candidate
  source).  Propagate the outcome to K5r, the DYN-8 ledger, and both
  papers.
  Status: **done** (2026-07-06;
  `code/audit9_dyn9b1d_lr_ratio_scan.py`, 7/7; ledgers
  `output/audit9/dyn9b1d_lr_ratio_scan.{json,md}`).  Method:
  polarization decomposition of the Hessians (three builds per
  pattern give the exact Hessian at ANY vev ratio; assembly gated
  against direct builds; analytic stationarity solve, gradient
  residual 1.9e-12; Goldstone 30 at every ratio incl. random
  epsilon).  RESULTS:
  **Q1 = YES**: at lambda_eps = 0 (the classic four-invariant family)
  the left-right stationary family IS a tree-level local minimum at
  8/8800 sampled (ratio x coupling) points, ALL at vev ratios
  0.6-0.7 -- He-Meljanac (1986) reproduced at claim level; 9b-1b's
  negative verdict at fixed a/p = 0.8 CONFIRMED as a ratio artifact
  (0.8 vs 0.6-0.7: missed by that little).  **Q2 = NO**: the fifth
  (epsilon) invariant neither rescues any of the 8792 negative points
  nor kills any of the 8 positive ones -- the seductive hypothesis
  dies honestly; the epsilon invariant's stability relevance is nil
  in the sampled space (its candidate-novelty status shrinks to the
  invariant-count statement itself).  K5r REVERTED accordingly in the
  collector and both papers: the left-right chain is 210-realizable
  at tree level in rare coupling regions; the 45_H route is an
  alternative, not a necessity.  Prior-art CONFLICT of SUB-A(ii):
  RESOLVED as our artifact.  Boundary: local (not absolute)
  positivity; He-Meljanac full text still unverified (paywall);
  full-sum coupling units.

- **SUB-B [theta_23 handling = DYN-4c executed] (1-2 sessions).**
  The benchmark Dirac card carries a 6.9 sigma theta_23 tension vs
  NuFIT 6.0 (DYN-4a); a referee strikes here first.  Strong option:
  kernel-level covariant Dirac refit (free perturbations around the
  archival kernels, refit to NuFIT 6.0, zeta/contact re-extraction,
  essentiality retest).  Fallback if the refit fails: demote the
  benchmark card to reproducibility-anchor-only status throughout the
  tex, citing the DYN-4b orthogonality result.
  Gates: post-refit chi^2 published; contact-essentiality verdict
  under the refit; tex updated accordingly.
  Status: **done** (2026-07-06;
  `code/audit1_dyn4c_kernel_dirac_refit.py`, 7/7; ledgers
  `output/audit1/dyn4c_kernel_dirac_refit.{json,md}`).
  FINDINGS: (1) EXACT ABSORPTION -- the inverse seesaw refits the
  Majorana sector to the NuFIT 6.0 centrals identically at UNCHANGED
  Dirac kernels (residual 9.6e-12, machine level): the 6.9 sigma
  theta_23 pull is a frozen-anchor property; chi^2_osc = 0 by
  construction, |zeta_fit| = 0.1275 inside the DYN-4b band.
  (2) ESSENTIALITY LIFTED: over 1000 kernel-perturbed refits (eps up
  to 0.3) cf never falls below 0.01 (worst 0.0151) -- the DYN-4b
  conditionality on the exact kernel shape is removed at the
  perturbation level.  (3) CHARGED-LEPTON ROUTE: a 5% Y_e kernel
  perturbation alone absorbs theta_23 below 1 sigma (best 0.1 sigma)
  with the Majorana sector frozen.  (4) Leptogenesis feed: modest
  (kappa gain ~1.1x at kernel level) -- the escape is real but not a
  magic lever.  K2 strengthened in the tex and the collector (K2r;
  29/29 against 19 ledgers).  Boundary: the GUT-scale sum-rule refit
  (quark sector) remains the deferred full flavor fit; m_1/Majorana
  phases at benchmark values.

- **SUB-1 [the theorem letter -- submission unit (1)] (1-2 sessions).**
  Extract the kinematic core into a short, self-contained paper:
  enumerated Spin(10)-16 uniqueness with nu^c forced; the family
  space as the automorphism algebra of its own carrier; N_fam <= 3 a
  priori, = 3 by the Cartan criterion; the Majorana contact = Killing
  form (B = 2 sqrt3 K_tr).  No phenomenology claims; the audits cited
  as machine verification.  Venue: PLB-style letter or PRD.  Language
  discipline: fully self-contained, no repo aliases.
  Gates: builds clean; claims worded per SUB-A verdicts; every
  theorem statement traceable to a ledger.
  Status: **draft complete** (2026-07-06;
  `route_E/tex/three_families_killing_contact.tex`, builds clean with
  the shared bibliography, 8 pp preprint).  Content: hypotheses
  H1/H2/H3 stated as explicit physical choices; the enumerated
  one-family uniqueness with the forced sixteenth state and the
  so(9)/so(15) near-misses; the genus ladder with the a-priori bound;
  the Cartan selection; the Killing-contact identity B = 2 sqrt3 K_tr
  with full proof; the Poincare-Hopf two-center remark and the
  "why m = 1" self-carried remark; the one-principle-two-consequences
  remark; a What-is-not-claimed section (incl. the one-sentence
  pointer to the coefficient-underivability theorem of the companion);
  the machine-verification section citing the 68 checks and the
  repository.  The ADJACENT-literature paragraph distinguishes the
  imposed-family-symmetry, index-counting, and division-algebra
  traditions per the SUB-A verdicts.  Fully self-contained language
  (no repository aliases).  REMAINING (author-side): read-through,
  venue selection (PRD letter-style vs PLB reformat), arXiv category,
  cover letter.

- **SUB-2 [the 210 tree-vacuum technical paper -- unit (2);
  CONDITIONAL on SUB-A(ii)] (2 sessions).**
  The 9b-1b/1c content as a standalone technical paper: the quartic
  basis with the machine-found Q5 = Q2 identity, the independent
  epsilon invariant and its slice-vanishing, exact Hessians with
  Goldstone gates, tree-level positivity verdicts (PS viable / LR
  dead with the epsilon lever exhausted), and the 126bar mixed-sector
  analysis.  Natural sequel to the 45+126 literature.  Venue: PRD.

- **SUB-3 [the framework flagship -- unit (3)] (1 session).**
  The reconstruction paper (with the falsifiability section) to
  arXiv FIRST for priority and feedback; journal submission deferred
  until SUB-1/SUB-2 outcomes are known.

Statuses start `open`; update per session.
