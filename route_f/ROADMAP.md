# Route F Roadmap: One Action, One Evidence Chain

Created: 2026-07-13
Last updated: 2026-07-20

Status values: `open`, `in-progress`, `done`, `failed`, `permanently-open`.
All items start `open` unless marked otherwise.

## AP-E12 endpoint density/regularity checkpoint (2026-07-20; current authority)

Status: **the endpoint was audited but not closed.  Degree is continuous in
the complete-minor graph norm, the AP-E11 density is exactly strongly
quasiconvex with Legendre--Hadamard constant one, and the direct Malý dipole
counterexample mechanism is excluded by a sharp scaling incompatibility.
Neither full endpoint graph density nor the relaxed-minimizer fallback is
proved.  The Hessian, determinant, and portal gates remain false.**  The
production card passes `9/9` without promoting its finite-grid basin proxy.

### New exact results

- Strong graph convergence of `u`, `Du`, `M_2(Du)`, and `M_3(Du)` implies
  convergence of the degree integral.  Since smooth degrees are integers, a
  target-valued graph recovery sequence is automatically fixed-degree for all
  sufficiently large indices.  No separate topology patch is needed.
- For `W(F)=(|F|^2+R|M_2(F)|^2+K|M_3(F)|^2)/2` and rank-one `H`, every minor
  is affine along `F+tH`, giving
  `D2W(F)[H,H]=|H|^2+R|DM_2(F)H|^2+K|DM_3(F)H|^2 >= |H|^2`.
- Periodic null-Lagrangian means give an exact strong-quasiconvexity identity:
  the integral excess over an affine gradient is one half the sum of the
  squared first-, second-, and third-minor increments.
- A radial dipole `r^alpha gamma(theta)` has finite `L2` second-minor energy
  only for `alpha>1/2`, while a non-vanishing forced filling cost requires
  `alpha<=1/2`.  The logarithmic borderline has the same incompatibility.
  This rules out that counterexample family, not every possible endpoint
  concentration.

### Why the two requested routes remain open

- Malý's Cartesian approximation theorem loses exponent (`q<p`); Mucci's
  area approximation controls all minors in `L1`.  Ordinary `W1,2(S3)`
  density does not control the endpoint `L2` minors.  The fixed-degree
  complete-minor density theorem is therefore **not proved and not
  disproved**.
- The formal classical stress is
  `P=F+R[(tr G)F-FG]+K F cof(G)`, `G=F^T F`, but a minimizer of the relaxed
  fixed-degree functional is not yet authorized to satisfy its naive
  Euler--Lagrange equation without local fixed-trace recovery.
- The model has pointwise `(p,q)=(2,6)` growth in dimension three.  The 2024
  partial-regularity theorem for relaxed strongly quasiconvex integrands
  requires `q<min(np/(n-1),p+1)=3`, and does not contain the sphere constraint
  or this graph relaxation.  Thus weak EL, partial regularity, classicality,
  and continuum isolation remain false.

### Numerical evidence and ordered continuation

- On `(N,L)=(17,5.5),(21,6),(25,6.5)`, maximum one-tetrahedron energy
  fractions decrease `3.966e-3 -> 1.849e-3 -> 1.173e-3`; all fields are
  stationary, admissible, and have three-target `B=[1,1,1]`.
- Three amplitude-`0.24` tangent perturbations re-relax to relative energy
  spread `3.38e-11`, radial-CDF distance `5.70e-7`, and field RMS `1.31e-5`
  after quotienting lattice translations and proper target `SO(3)`.  This is
  finite-dimensional same-basin evidence only.
- Next prove an annular target-valued endpoint replacement/equiintegrability
  lemma.  The parallel analytic route is a reverse-Hölder estimate giving
  local `L^(2+delta)` control of the complete-minor vector, followed by a
  controlled `S3` projection.  Only after classicality and continuum
  isolation close may the same-action Riemann Hessian be assembled.

Canonical artifacts:

- `route_f/tex/ap_e12_graph_density_regular_minimizer.tex` and PDF;
- `route_f/code/verify_ap_e12_graph_density_regular_minimizer.py`;
- `route_f/output/ap_e12_graph_density_regular_minimizer.{json,md}`;
- AP-E12 section in the master Route-E derivation ledger.

## AP-E11 compatible cochain checkpoint (2026-07-19; superseded above)

Status: **the tetrahedral discrete-de-Rham replacement, its exact cell
formula, and the declared translation/triangulation quotient gate are
closed.  Gamma convergence is unconditional only to the fixed-degree
lower-semicontinuous relaxation.  Identification with the unrelaxed smooth
sector is still open, so the Hessian, determinant, and portal gates remain
false.**  The production card passes `11/11`; all 13 unanchored backgrounds
are stationary, admissible, and have three-target `B=[1,1,1]`.

### Compatible complex, Hodge star, and degree

- The Alexander--Whitney front/back cup product obeys the exact graded
  Leibniz identity and associativity.  For `a^A=dn^A`, full antisymmetrization
  of `a^A cup a^B` and `a^A cup a^B cup a^C` equals the affine Whitney
  two- and three-minor cochains, including the necessary `1/2!` and `1/3!`
  factors.
- The element Whitney mass matrices `M_k^T` for `k=1,2,3` are positive Gram
  matrices, rotation covariant, and exact on constant affine forms.  The
  smallest measured eigenvalue is `2.485e-2`; the maximum rotation residual
  is `3.553e-15`.
- The same ordered triple cup gives
  `C_h(T)=det[n_0,n_1,n_2,n_3]`.  Radial pullback of the normalized affine map
  gives the exact local degree cochain
  `Omega_h(T)=s_T C_h(T)/(2 pi^2) int_Delta |sum lambda_i n_i|^(-4) dlambda`.
  Independent Duffy--Gauss and three-regular-value calculations agree within
  `3.069e-7`.

### Cell formula and continuum statement

The compatible action is the positive Whitney-Hodge norm of the complete
first, second, and third affine minors, with coefficients `R=1`, `K=0.35`,
plus the mass potential.  Its exact periodic cell formula is

`Q_(R,K)(A)=|A|^2/2 + R|M_2(A)|^2/2 + K|M_3(A)|^2/2`.

First-, second-, and third-minor means are periodic null Lagrangians.  Jensen
in the complete-minor vector proves the zero corrector globally for the
uniform six-tet and both checkerboard five-tet cells.  Numerical cell
optimization agrees within `2.555e-14`.

Dirichlet coercivity, weak minor identification, the normalized-current
formula, polyconvex lower semicontinuity, and smooth nodal recovery prove
Gamma convergence to the fixed-degree lower-semicontinuous relaxation.  The
remaining theorem blocker is density of smooth fixed-degree maps in the
complete-minor graph norm, or an equivalent regularity/isolation theorem for
the selected relaxed minimizer.  Until one of those is proved,
`full_unrelaxed_density_theorem=false` and no classical Hessian may be built.

### Quotient gate and ordered continuation

- Six joint-limit backgrounds reach `(N,L,a)=(37,8,0.222222)`; three
  translated starts and four five-tet controls share `(29,7,0.25)`.
- Joint-tail energy, RMS-size, and radial-CDF spreads are respectively
  `0.437%`, `3.489%`, and `4.067%`; translation residuals are below
  `7.2e-7`; cross-triangulation spreads are `0.463%`, `0.752%`, and `3.556%`.
  All declared numerical thresholds pass.
- Therefore `compatible_cell_formula=true`, `quotient_gate=true`, and
  `relaxed_regulator_stable=true`, while classical `regulator_stable`,
  `hessian_gate_open`, `determinant_variation_gate_open`, and `portal` remain
  false.
- Next, prove graph-norm density or regularity/isolation of the relaxed
  minimizer.  Only then assemble the same-action Riemann Hessian and regulated
  determinant variation.  GW/domain-wall and the actual `SO(3)` mod-two index
  remain parallel lanes.

Canonical artifacts:

- `route_f/tex/ap_e11_compatible_cochain_action.tex` and PDF;
- `route_f/code/scan_ap_e11_compatible_cochain_action.py`;
- `route_f/output/ap_e11_compatible_cochain_action.{json,md}`;
- AP-E11 section in `route_f/tex/another_physics_route_e_derivation_ledger.tex`.

## AP-E10 stencil/cell/translation checkpoint (2026-07-18; superseded above)

Status: **the current one-corner forward Skyrme mainline is refuted as a
continuum regulator; the finite-`R` cell problem is closed with a negative
universality result; translation quotienting is implemented; the continuum,
Hessian, determinant, and portal gates remain false.**  The production card
passes `10/10`, and all 11 relaxed fields pass stationarity, admissibility,
and three-target degree.

### Exact compensated-compactness and current result

- A three-periodic exact `S3` sequence converges uniformly to the vacuum while
  one forward two-minor equals `(3 sqrt(3)/2) eta^2` at every site.  Uniform
  Dirichlet and Skyrme energy bounds therefore do not imply weak continuity
  of the declared one-corner minors.
- A two-periodic hemisphere-valued sequence has normalized-affine degree zero
  but forward current `J_a=-16 eta^3 sqrt(1-3 eta^2 a^2)`.  With the exact
  vacuum cutoff `chi=prod sin^2(pi x_i)`, the weak current defect is
  `-16 eta^3 chi^3 dx`.
- The singular-value inequality
  `integral |J_a|^(4/3) <= (1/3) integral |wedge^2 D^a n|^2` excludes singular
  concentration of the forward current.  It does not remove the diffuse
  oscillatory defect or control the geometric normalized-affine current.

### Exact finite-R cell formula

For periodic bonds `B_Y`,

`Q_hom(A)=|Y|^(-1) inf_phi sum_(p,r in B_Y) |A r + phi(p+r)-phi(p)|^4`.

Every bond-direction graph is cycle-balanced for the uniform six-tet and
checkerboard five-tet period-two cells.  Direction-wise Jensen proves
`phi=0` exactly.  Therefore the AP-E9 displayed `Q_6` and `Q_5` are the true
homogenized barrier densities.  They remain nonproportional: the six/five
ratios are `4/3` on `(v,0,0)` and `3` on `(v,v,0)`.  The two checkerboard
phases agree, but five- and six-tet regulators do not.

### Quotient scan and gate

- The joint sequence `(N,L)=(25,6),(29,6.5),(33,7),(37,7.5)` reaches
  `a=0.208333`.  Four translated starts and both five-tet phases are compared
  on the common `N=29,L=6.5` grid.
- All 11 fields have projected-gradient density below `1.71e-6`, edge margin
  at least `0.1064`, and `B=[1,1,1]`.
- Translation quotienting uses the `1-n0` barycentre, covariance, and aligned
  radial profile.  Dynamic shifts are accepted only if admissibility and
  degree persist, then the same action is re-relaxed.  Re-relaxation exhibits
  Peierls locking, so the quotient rather than raw center is authoritative.
- Joint-tail energy spread is `1.656%` and passes, but centered-RMS spread is
  `9.612%` and successive profile distance is `7.072%`.  Translated-start
  profile distance is `3.500%`.  Cross-mesh energy and radius spreads are
  `5.143%` and `11.378%`.  Those four gates fail.

### Ordered continuation

1. Replace the one-corner product by a complete tetrahedral discrete de Rham
   construction: primal one-cochain `dn`, graded shifted cup product with a
   discrete Leibniz identity, and positive Hodge star.
2. Construct the topological three-cochain from that same cup product and
   prove equality with normalized-affine degree on the admissible set.
3. Calibrate or symmetrize the Hodge star across body diagonals, solve the new
   cell problem, and demand triangulation-independent quartic response.
4. Re-run the translation-quotient `a,L` scan only for that replacement.
5. Build the same-action Riemann Hessian and regulated determinant variation
   only if the replacement closes both theorem and numerical gates.
6. Retain finite GW/domain-wall and actual `SO(3)` mod-two-index lanes in
   parallel; the degree-one Route-E portal remains last.

Canonical artifacts:

- `route_f/tex/ap_e10_compactness_homogenization_centering.tex` and PDF;
- `route_f/code/scan_ap_e10_compactness_homogenization_centering.py`;
- `route_f/output/ap_e10_compactness_homogenization_centering.{json,md}`;
- AP-E10 section in `route_f/tex/another_physics_route_e_derivation_ledger.tex`.

## AP-E9 scaling/Gamma-limit checkpoint (2026-07-18; superseded above)

Status: **fixed-box strong-`L2` equicoercivity, the barrier zero-Gamma-limit,
and positive-`R` degree compactness are closed; the finite-`R` homogenized
density, full-action Gamma limit, and a regulator-stable continuum background
are not.**  The production
card passes `8/8` implementation/provenance checks and all 23 relaxed cases
pass stationarity, admissibility, and three-target degree, but the predeclared
continuum gate returns false.  The Hessian/determinant embargo remains active.

### Exact scaling result

Write
`d(a)=1-epsilon(a)`, `w(a)=gamma(a)*a`, and
`R(a)=gamma(a)*a^2/d(a)^2`.  For
`x_e=(n_x dot n_y-epsilon)/d`, the edge function is
`phi(x)=-log(x)-1+x`, and

`b_epsilon(n_x dot n_y) >= |n_x-n_y|^4/(8 d(a)^2)`.

- Every full-action sublevel on a fixed box is strongly `L2 x L2`
  equicoercive for every positive `gamma(a)` and
  `-1/3<epsilon(a)<1`; the coordinate-edge Dirichlet term supplies a uniform
  piecewise-affine `H1` bound.  This compactness does **not** preserve degree.
- On every fixed `C>0` sublevel the exact Lambert bound is
  `x_e >= -W_0(-exp(-1-C/w)) >= exp(-1-C/w)`.  Uniform relative distance from
  the edge floor holds iff `inf w>0`; an absolute margin also needs
  `inf d>0`.  A one-vertex vacuum-collar sequence disproves the margin when
  `w->0`.
- For a smooth field, the barrier is
  `(R/8) integral Q_T(grad n) + o(R)`.  With fixed epsilon and
  `gamma=c*a^(-p)`, uniform edge interiority needs `p>=1`, while smooth
  disappearance needs `p<2`.  The unique window is therefore `1<=p<2`; the
  recommended minimal path is `p=1`, `gamma=bar_gamma/a`.
- With the barrier embedded by dual-cell piecewise-constant fields, its
  strong-`L2` Gamma limit for `R->0` is zero on `L2(S3)` and infinity outside;
  it forgets trace and degree.  For `inf R>0`, the quartic lower bound gives
  `W1,4`/Holder compactness and closes degree; `R->infinity` forces the vacuum.
  At finite nonzero `R`, however, the smooth-sampling `Q_T` is only a raw
  Cauchy--Born density.  Checkerboard meshes require a homogenized cell
  formula that has not been evaluated.  The six-tet and five-tet raw tensors
  have ratios `4/3` and `3` on two test gradients, so they are not scalar
  multiples before corrector minimization.
- A complete same-action Gamma theorem follows conditionally if the AP-E7
  base actions are equicoercive/Gamma-convergent, finite-energy limits admit
  energy-dense fixed smooth approximants, and nodal samples recover each fixed
  approximant, with `a^2/d->0` and `R->0`.  Those hypotheses remain unproved
  for the one-corner forward Skyrme minors.

### Production scan and negative continuum verdict

- The scan covers six spacings at fixed `L=6`, four volumes through `L=12` at
  fixed `a=0.4`, `p=0,1,2` controls, both checkerboard five-tet phases, the
  uniform six-tet mesh, and centered/noninteger-translated starts at
  `a=0.4` and `a=0.3`.
- All 23 fields are unanchored stationary points with projected-gradient
  density at most `1.723e-6`, pair margin at least `0.067786`, and
  `B=[1,1,1]`.  Translated starts return the same mesh-specific branches.
- The relaxed `p=1` fixed-box barrier has fitted slope `0.2196`, outside the
  required `1.0+-0.4`.  The largest-three-volume energy spread is only
  `0.1001%`, but the centered-radius spread is `8.521%`, above `5%`.
- Six-versus-five-tet total-energy spreads are `5.442%` at `a=0.4` and
  `5.355%` at `a=0.3`, both above `3%` and not visibly decreasing.  On the
  common `N=21,L=6` grid, the two vanishing-regulator controls `p=0,1` differ
  by `6.446%` in total energy, above `3%`.
- Therefore
  `numerical_regulator_stability=false`,
  `C_B1_physical_continuum=false`,
  `same_action_riemann_hessian_allowed=false`, and
  `regulated_determinant_variation_allowed=false`.

### Remaining theorem blockers and ordered continuation

1. Prove a discrete compensated-compactness/Gamma-liminf theorem for the
   forward Skyrme minors, smooth fixed-boundary fixed-degree density in the
   natural minor-bounded class, and exclusion/accounting of topological-current
   concentration defects.
2. Evaluate the finite-`R` periodic homogenized cell formula before treating
   either displayed `Q_T` as a Gamma density; this is mandatory if a critical
   quartic branch is retained.
3. Continue the `p=1`, fixed-epsilon branch to smaller `a` only with adaptive
   boxes and a translation quotient/physical centering condition; require the
   barrier, radius, triangulation, and regulator gates simultaneously.
4. If the base-minor theorem fails, declare the critical quartic as physical
   and replace the raw stencil by a rotationally calibrated,
   triangulation-independent discretization before restarting the continuum
   study.
5. Build the complete same-action Riemann Hessian and regulated determinant
   variation only after one continuum background passes all gates.
6. Retain the finite GW/domain-wall kernel and actual `SO(3)` Yukawa
   mapping-torus mod-two index as parallel lanes.  The degree-one Route-E
   portal remains last.

Canonical artifacts:

- `route_f/tex/ap_e9_gamma_limit_regulator.tex` and its PDF;
- `route_f/code/ap_e9_triangulation_tools.py`;
- `route_f/code/scan_ap_e9_gamma_triangulation.py`;
- `route_f/output/ap_e9_gamma_triangulation.{json,md}`;
- AP-E9 section in `route_f/tex/another_physics_route_e_derivation_ledger.tex`.

## AP-E8 topology-preserving finite-grid checkpoint (2026-07-17; superseded above)

Status: **the finite-grid topology/stationarity sub-blocker is closed; the
physical continuum Lane-B conjunction, the other two lanes, and the portal
remain open.**  The production card passes `13/13` and sets
`physics_promotion_allowed=false`, `portal_start_allowed=false`, and
`lane_closed=false`.

### Why this path was selected

- A finite GW/domain-wall target kernel would advance Lane A but still needs
  the actual target mass-family lift, determinant orientation, mixed heavy
  threshold, and all-scale `k=+2` mechanism.
- The actual `SO(3)` Yukawa mapping-torus mod-two index would select a torsion
  Pfaffian line in Lane C but would not change AP-E7's finite-site unwinding.
- The topology-preserving action places a finite-energy boundary strictly
  before the AP-E7 dislocation locus, creates exact finite-grid components, and therefore
  changes a proved blocker rather than adding another proxy calculation.

### What is now proved and computed

- **Action:** on the unique unordered vertex pairs of a fixed Freudenthal
  triangulation,
  `E_TP = E_APE7 + gamma*a*sum b_epsilon(n_x dot n_y)`, with
  `b_epsilon(t)=-log((t-epsilon)/(1-epsilon))+(t-1)/(1-epsilon)` for
  `t>epsilon` and `+infinity` otherwise.  Production values are
  `epsilon=0.01`, `gamma=1`.
- **Exact topology:** if every tetrahedron pair dot exceeds `epsilon`, then
  every normalized-affine denominator is at least
  `sqrt((1+3 epsilon)/4)=0.507444578255`.  Degree is locally constant and no
  continuous finite-energy path changes `B`.
- **Exact existence:** in every nonempty fixed-grid admissible component,
  barrier-bounded sublevels stay a positive distance from the exceptional
  set.  Compact `S3` variables plus the coercive `m_s>0` potential yield an
  interior minimizer with no core or centre pin.  Sampling a compact-supported
  continuum degree-one field proves nonempty `B=1` components on every
  sufficiently fine grid.
- **Seven direct stationary representatives:** fixed `L=6` at
  `N=15,17,19,21` and fixed `a=0.4` at `(N,L)=(16,6),(21,8),(26,10)` all have
  three-target `B=[1,1,1]`, positive pair margins, and direct tangent-plus-
  scalar gradient densities between `1.052e-7` and `7.704e-7`, below the
  declared `2e-6` tolerance.
- **Controls:** four nearby `(gamma,epsilon)` values and one deterministic
  tangent/scalar perturbation return admissible stationary `B=1` endpoints.
  The fixed-`a` energy spread is `0.459%`, but the weighted radius still changes
  by about `18%`; no volume convergence is claimed.  The strong gamma
  dependence is evidence that regulator independence remains open.
- **Calculus:** the full chart gradient, differentiable solver continuation,
  and intrinsic barrier Hessian have directional residuals `9.38e-11`,
  `2.81e-11`, and `8.16e-9`.  A boundary-compatible compact-supported profile
  has global barrier exponent `2.1591` and local slopes approaching two, in
  agreement with the analytic smooth-field `O(a^2)` expansion.

### Fail-closed boundary

The following remain false: numerical certification of the global minimizer,
complete line-search-segment admissibility, equicoercivity, Gamma convergence,
a continuum stationary limit, barrier/triangulation independence, a physical
charged-QC2D action, total same-action Riemann Hessian, regulated determinant
second variation, interacting gauge/meson/ghost blocks, BRST
superdeterminant, four-dimensional dynamics/HMC, and the quantum continuum
limit.  Accordingly
`C_B1_finite_grid=true` but `C_B1_physical_continuum=false`.

### Ordered continuation

1. Establish an equicoercive/Gamma-liminf estimate for a declared
   `gamma(a),epsilon(a)` scaling, or construct a barrier-core-collapse
   counterexample.  Track local gradient concentration and a topological
   radius in addition to energy.
2. Repeat the stationary family for at least two alternative cube
   triangulations, translated starts, wider `a,L`, and several regulator
   trajectories; extrapolate observables jointly rather than calling the
   barrier energy physical.
3. Only after a regulator-stable continuum background appears, assemble the
   full same-action tangent Hessian with the sphere-curvature shift, then add
   interacting gauge/ghost blocks and the bosonic second variation of one
   regulated fermion determinant.
4. In parallel, complete either the finite target GW/domain-wall operator or
   the actual `SO(3)` Yukawa mapping-torus mod-two index and microscopic CPT
   regulator.
5. Start a degree-one Route-E portal only after one complete physical lane
   closes; AP-E8 alone does not authorize it.

Canonical artifacts:

- `route_f/tex/ap_e8_topology_preserving_b1.tex` and its PDF;
- `route_f/code/scan_ap_e8_topology_preserving_b1.py`;
- `route_f/output/ap_e8_topology_preserving_b1.{json,md}`;
- AP-E8 section in `route_f/tex/another_physics_route_e_derivation_ledger.tex`.

## AP-E7 regulator/topology/family checkpoint (2026-07-17; superseded above)

Status: **five exact subresults, three physical lanes still open, portal not
started.**  The lane cards pass `42/42`, `16/16`, and `18/18`; the independent
execution-frontier gate passes `16/16`.

### Lane A -- common APS/PV and threshold

- Done: the exact background quadratic complex, common PV determinant
  prescription, finite local scheme, fixed-zero-cut APS domain, and pure
  two-flavour heavy gauge/gravity index.  The closed-spin phase is exactly
  `+1` for every bundle.
- Conditional only: the separated Higgs projectors define a GW-compatible
  target source and even mod-two counts `(+,+)`.  A complete finite
  overlap/domain-wall target operator and its determinant orientation are not
  yet present.
- Permanent design constraint: a continuous equivariant **exact separated
  rank-two projector** cannot survive the restored `Sp(4)` point.  New UV
  fields, a gap closing, or an explicitly emergent symmetry are required.
  The theorem does not rule out every non-idempotent coupling.
- Open: original target mass-family lift, actual mapping-torus pair, mixed
  heavy determinant, nonperturbative continuum measure, and all-scale `k=+2`.

### Lane B -- discrete topology before super-Hessian

- Done: the finite-site configuration space is path-connected, and the
  three-target Freudenthal degree is locally constant only on the admissible
  subset.  Five full-grid independent solves unwind to `B_geom=0`; guarded
  descent reaches each positive admissibility floor with nonzero gradient.
- Control only: independently re-relaxed fixed-core sectors retain
  `B_geom=1`, have free constrained gradients below `6.2e-8`, and positive
  exact restricted `n+s` Hessians.  Their full unanchored gradients are
  `0.24--0.66`.
- Open: an unanchored admissible `B=1` stationary sequence, interacting gauge
  and meson cross blocks, the second derivative of the regulated fermion
  determinant, BRST superdeterminant, dynamics/HMC, and continuum stability.

### Lane C -- physical `SO(3)` line or rank-three family

- Done: `H^2(SO(3);Z)=Z2`, the FR holonomy and determinant square, and a
  topological CPT Real lift.  The standard conditional `N_c=2,B=1` rule
  favours the trivial/bosonic line.
- Done mathematically: a basepoint-free `O(1,2)` triple defines the minimal
  three-band mass with exact gap and `c1(E+)=x+2y`; 11,664 samples and Berry
  meshes reproduce `(c1_x,c1_y)=(1,2)` to roundoff.
- Open physically: the actual mapping-torus mod-two index, microscopic CPT
  regulator, localized physical replacement for the spectator `CP1`,
  three-band charged-two-colour Yukawa embedding, uniform Fredholm gap, and
  gauge-basic descent.

### Ordered continuation

1. Complete one literal finite regulator/operator, including target Yukawa
   kernel and endpoint orientation, or change the UV field content explicitly.
2. Replace the naive lattice topology by a controlled admissible or
   spherical-volume action; re-relax at several `a,L`; only then compute the
   physical determinant Hessian and interacting gauge/ghost blocks.
3. Compute the actual `SO(3)` mapping-torus Pfaffian/CPT data, with the
   physical rank-three embedding as the parallel alternative.
4. Start a degree-one Route-E portal only after one full conjunction closes.

## AP-E6 same-background checkpoint (2026-07-16; superseded above)

Status: **three rigorous advances, three full lanes still open, and no
portal work started.**  This round strengthens the meaning of closure rather
than promoting a partial calculation.

### Lane A -- simply connected `Sp(4)`

- The revised **39/39** card computes, in physics/mathematics notation,
  `Omega_5^Spin(BSp(4)_phys)=Omega_5^Spin(BSp(2)_math)=Z2`.  A unit instanton
  on `S4` times the nonbounding spin circle generates it.  Restriction to the
  instanton `SU(2)` gives mod-two characters
  `(chi_4,chi_5,chi_10)=(-1,+1,+1)`, so the displayed two Dirac `5`s have
  trivial dynamical gauge-bordism character.
- This does **not** select the two independent target characters on
  `G3=S1_R x S3` and `G2=T2_Arf=1 x S2`.  Their common microscopic
  Dirac--Yukawa mass-family lift and APS problem are still absent.
- Fermionic gamma-five reality and three PV moment identities pass, but the
  full interpolation/kernel, regulator statistics, APS domain, finite local
  scheme, and gauge--scalar--ghost regulators are not specified.  Hence
  `sp4_euclidean_regulator_complete=false`.
- The light `B=+1` orientation is `k_IR=+2`.  The complete-`5` identity
  `+2-2=0` is only a conditional uniform-flavour obstruction.  A local
  pure-gauge theta-periodicity check is not an unbroken/gravitational APS
  determinant ratio; `heavy_threshold_eta_matched=false` and the lane stays
  open.

### Lane B -- one canonical relaxed `B=1` field

- The **20/20** card solves a coupled continuum `F(r),s(r)` hedgehog-sector
  BVP, not the old analytic profile.  It is
  solved at `R=8,10,12,16`.  At `R=16`,
  `E/(4 pi)=6.1282155`, `s(0)=0.9171289646`, and `B=+1`.  The canonical
  4097-point little-endian `(r,F,s)` checksum is
  `81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59`.
- The analytic CSR `n+s` Hessian is exactly differentiated from one declared
  cubic energy, includes the sphere multiplier, and projects translations and
  isorotations.  The diquark and free toron/ghost blocks use separately
  declared grids/boundaries and are audited rather than merged.  Consequently
  the aggregate same-grid sparse-Hessian gate is false.
- The crucial negative result is numerical, not semantic.  At fixed `L=8`,
  `(N,B_a,gradient_density)` is approximately
  `(9,.168,.348)`, `(13,.440,.430)`, `(17,.625,.422)`.  The sampled
  continuum solution is not stationary for the cubic action and its lattice
  degree is not converged.  Its negative projected curvatures are therefore
  off-shell diagnostics, not physical instability eigenvalues.  The sampled
  boundary also differs from the exact vacuum by about `.0493` at `L=8`.
- The `n_s_sparse_constrained_second_variation_complete` subgate is true, but
  `sparse_projected_hessian_complete_in_declared_sector`,
  `discrete_stationarity_achieved`, `lattice_topology_converged_to_B1`,
  `multigrid_multivolume_converged`, the full fermionic super-Hessian, 4d
  importance sampling, and the quantum-continuum gate all remain false.

### Lane C -- Yukawa/Callias/CPT/WZW on that checksum

- The revised **25/25** verifier imports Lane B's public solver and hard-matches the
  canonical `(r,F,s)` checksum.  Its declared Yukawa term couples to `F` and
  explicitly decouples `s`; it does not silently substitute another soliton.
- The localized internal stabilizer is exactly `Z2`, so the physical orbit is
  `SU(2)/Z2=SO(3)`, not `CP1`.  Combined space--isospin quotienting gives the
  same orbit.  The separate scalar-vacuum `CP1` has norm proportional to
  volume and is a bulk Goldstone direction.
- Constant asymptotic mass makes the ordinary static Callias boundary bundle
  trivial and its index zero.  Four-dimensional spectral flow can still be
  `B=1`; it does not force a static endpoint zero mode.
- A new two-band no-go is exact under its stated assumptions: a globally
  gapped `2x2` mass in a trivial rank-two bundle has `c1(E_+)^2=0`, whereas
  the desired `c1=x+2y` has square `4xy`.  A future model needs at least three
  bands, a nontrivial ambient bundle, or an infinite-dimensional family.
- The universal differential character is now typed on `SU(2) x CP1` before
  evaluation and fiber integration.  This is a local restricted ansatz, not
  the missing full equivariant cocycle/vertical-holonomy proof.  The physical
  CPT regulator, `O(2)` determinant line, gauge-basic microscopic descent,
  and same-soliton composition remain false.

### Portal decision and next ordered work

The independent **15/15** AP-E6 frontier card recomputes each lane as a conjunction.
All three are false, so `any_preportal_route_closed=false`, portal start is
unauthorized, `degree_one_portal.constructed=false`, and
`physics_promotion_allowed=false`.

1. **Regulator/threshold route:** write one explicit five-dimensional
   interpolation with APS domains and local scheme for every fermion, PV,
   gauge, scalar, and ghost field; compute both target generator phases plus
   the unbroken and gravitational heavy ratios.  A direct
   `SU(2)c x SU(2)H` completion remains the backup if complete-`5` threshold
   matching cannot preserve `+2`.
2. **On-shell Hessian route:** re-relax the field for each cubic action until
   the projected gradient meets a preregistered tolerance, use a geometric
   lattice-degree estimator, then repeat spacing/volume limits before adding
   Wilson/overlap determinant curvature and dynamical 4d ensembles.
   Before raising the Hessian-completeness gate, add pure-`n`, pure-`s`,
   `n-s` cross, and cell-local geodesic probes; replace raw IPR by a physical
   participation volume and enlarge the three-point diquark volume fit.
3. **Physical-moduli route:** preferentially quantize the actual `SO(3)` (or
   its `SU(2)` FR cover), calculate its torsion determinant line and CPT/FR
   holonomy, and ask whether that observable can replace the discarded free
   `O(2)` class.  In parallel, test a rank-three mass-family construction
   against the `c1^2` obstruction and a uniform Fredholm-gap gate.
4. Build a degree-one Route-E portal only after one entire route above closes.

Canonical evidence is under `route_f/tex/ap_e6_*`, `route_f/code/*ap_e6*`,
and `route_f/output/ap_e6_*`; the master derivation ledger records the full
formula chain and rollback boundary.

## AP-E completion-frontier checkpoint (2026-07-16; superseded above)

Status: **formal derivations and finite controls complete; all three
pre-portal physics lanes remain open.** The integrated promotion verifier is
**20/20** and records
'any_preportal_route_closed=false',
'degree_one_portal.constructed=false', and
'physics_promotion_allowed=false'.

### Lane 1 -- four-dimensional charged-QC2D/Hessian

- A Wilson/Symanzik-compatible four-dimensional action, background gauge
  fixing, FP operator, meson/charge-two-diquark sector, Wilson fermions, and
  the full formal graded Hessian are explicit.
- The **24/24** deterministic card tests a much smaller '2^4' trivial-link
  frozen-background block. Numerical anchors are
  'B(0)=0.999981468235',
  'lambda_min(K_mes,test)=0.621297848377',
  'lambda_min(K_Delta)=0.694069172619',
  the tachyonic control '-0.349308273811',
  'sigma_min(D)=0.702570765924', and Wilson/Symanzik orders
  '1.99861/3.99628'.
- This is not importance sampling, a dynamical charged-QC2D phase,
  determinant positivity, the complete nonlinear Skyrme Jacobi Hessian, a
  renormalized continuum limit, FR quantization, or global stability. All ten
  corresponding JSON gates remain false.

### Lane 2 -- APS generators and semisimple global form

- 'G3=S1_R x S3' and 'G2=T2_RR x S2' explicitly generate the two reduced
  spin-bordism factors. Their circle/Arf mod-two indices are one. The ambient
  chirality-paired product spectrum gives only a reference '(+1,+1)' phase;
  defect regulators realize all four torsion characters. This proves
  regulator dependence, not the unique charged-QC2D Dai--Freed determinant.
- The **58/58** search identifies simply connected
  'Sp(4)=USp(4)=Spin(5)' as the strongest simple candidate in the scanned
  set:
  '5 -> 2_(+1)+2_(-1)+1_0',
  '4 -> 2_0+1_(+1)+1_(-1)'.
  Two complex '4' copies preserve the separate 'SU(2)_phi' copy symmetry and
  give the neutral triplet 'qq(phi dagger)^2' at group level.
  Direct 'SU(2)c x SU(2)H' is the minimal backup; simply connected 'SU(4)'
  is secondary. Their diagonal/central quotients fail the odd singlet screen.
- The candidate is asymptotically free at one loop ('b0=15/2') and passes
  displayed vectorlike/Witten checks. It is not closed: the light threshold
  is tuned; full gauge bordism, heavy eta matching, radiative protection,
  monopoles, strong-vacuum selection, and the stable 'B=1' soliton are open.
  The exactly-two unit-cell Gauss rule is not rederived.

### Lane 3 -- same-soliton SQM/Callias/CPT/descent

- A low-energy tangent-SQM characterization and conservative sufficient
  same-model audit are explicit; microscopic bulk supercharges are treated
  as sufficient, not logically necessary, because emergent worldline 'N=2'
  is possible.
- For 'c1(E_+)=p x+q y', the families-Callias template gives
  '(rank,c1(det Ind D))=(epsilon p,epsilon p q)'. Rank-one 'O(+2)' requires
  'epsilon p=1,q=2'. The Berry/template control yields
  'c1=1.999999999999981' and a unit gap, but no physical same-model
  Dirac--Yukawa family or Fredholm gap has been supplied.
- Fixed-polarization CPT uses
  'E_-=K_CP1 tensor E_+^vee=O(-4)' and gives three opposite-chirality modes;
  raw 'O(-2)' gives one. Spatial differential-character integration
  'Hhat^5 -> Hhat^2' is defined, while equivariant refinement, basic
  curvature, large-gauge/stabilizer holonomy, and the same-soliton integral
  remain open. The card passes **27/27** with composition false.

### Portal gate and ordered continuation

The degree-one portal is authorized only if at least one complete lane above
closes. None has. Its future theorem must still prove a degree '+1' map,
pull back 'O(2)' with the physical orientation, give an anomaly-safe
operator, preserve CPT, add no zero modes, and lie below all retained gaps.

1. Complete the simply connected 'Sp(4)' Euclidean/regulator action, evaluate
   both target-generator mapping-torus eta phases and
   'Omega_5^Spin(BSp(4))', and match/protect every heavy threshold. Use the
   unquotiented product group as backup.
2. Relax the same charged 'B=1' solution for multiple spacings and volumes,
   assemble the exact sparse gauge-fixed nonlinear super-Hessian, project
   gauge/translation/isorotation directions, and run dynamical ensembles
   after the measure sign is proved.
3. Compute the same-soliton Yukawa--Callias family, determinant line, CPT
   regulator, equivariant WZW class, and vertical/large-gauge holonomies.
4. Construct the degree-one Route-E portal only after one full lane closes.

Canonical evidence:
'tex/ap_e5_4d_qc2d_lattice_hessian.tex',
'tex/ap_e3_aps_global_form_search.tex',
'tex/ap_e4_same_soliton_callias_descent.tex', their verifiers/cards/PDFs,
and 'code/verify_ap_e5_completion_frontier.py'. The master derivation ledger
has been synchronized.

## Blocker execution checkpoint (2026-07-14; supersedes older F0-A text)

- **Artifact recovery: done.** RE-SC3/4/5 scripts and ledgers exist, are
  tracked, and make the full `--string-cards required` dry-run preflight clean.
  They remain `unpromoted_pricing_only`; file presence is not a physics pass.
- **Single implementation: done.** `route_E/code_dyn/` is canonical and the
  19 root DYN paths delegate to it.  The runner/resolver use the exact
  case-sensitive `route_E` path in source lookup, isolated staging, and
  execution.
- **P0 numerical repairs: done.** DYN-9b-2 now uses already-normalised
  `g1=0.462` (`y_t(M_X)=0.44116481`); DYN-9b-3 uses one square root in the
  light spectrum (`epsilon_DI^SM=2.307794855e-6`) and no longer calls reheating
  unconstrained; DYN-8 resolves its 210/45 branch-map contradiction and reads
  DYN-5V/DYN-7F.
- **P0 inference-boundary repairs: done.** `Y_nu=h-3f`, `Y_u=h+f` no longer
  masquerades as an SO(10) top-like lock; `h=3f` is the explicit
  counterexample.  DYN-9b-2 separates the actual archival-kernel suppression
  (`19.5x/342.2x`) from the optional top-like tension (`9.6x/169.2x`) and
  restricts exact zeta invariance to uniform positive-real rescaling.  D3's
  `N_2,N_3>M_I` ordering is fixed-tower-specific, D4's numerical gap is
  historical-invalid, and D5's `M_SS<M_*` is only a necessary scale ordering.
- **Promotion guard: done.** `code/audit_blocker_promotion_gate.py` is `18/18`
  and deliberately returns `physics_promotion_allowed=false`.  DYN-8 is
  `30/30` mechanical/disclosure with the same non-promotion result.
- **H3 logic repair: done; physical origin open.**  The one-dimensional
  abelian counterexample is now machine-checked.  Original H3 gives only
  `N_fam<=3`; `N_fam=3` is explicitly conditional on the H3+
  nondegenerate adjoint-trace/Killing-contact axiom.  Motivating or realizing
  H3+ dynamically remains an open physics problem.
- **Still blocking F0-A/F5:** a global branch-local non-SUSY Spin(10) flavor
  fit, tau-resolved Boltzmann/density-matrix kinetics with spectator/reheating
  inputs, a valid interacting DYN-5 messenger action, recomputed RE-SC4
  pricing, and a threshold/experimental-bound envelope for RE-SC5.

## AP-E3 global/nonlinear and AP-E4 SQM checkpoint (2026-07-16; superseded above)

Status: **three requested lanes executed with theorem-level no-go and
fail-closed separation; no Route-E physics promotion.**

- **Charged two-colour nonlinear proxy:** the dedicated program scans a
  declared meson/charge-two-diquark linear-sigma vacuum over a two-parameter
  matching ansatz, records condensates and the six-field Hessian, and solves
  the full nonlinear `B=1` massive-Skyrme hedgehog.  Baryon number, Derrick
  virial, box/grid convergence, the radial generalized meson Hessian, charged
  scalar `l=0,1` ordering, and a tachyonic negative control are machine
  checked.  This closes a classical nonlinear-EFT necessary gate only:
  `full_3d_hessian_closed=false`, `lattice_qc2d_closed=false`, and the low-
  energy matching coefficients remain UV inputs.
  The `18/18` benchmark records phase onset
  `mu_lift^2=-m_pi^2=-0.25`, gaps `(m_pi,m_sigma,m_Delta)=(0.5,3.5,0.9)`,
  `|B-1|=4.06e-10`, relative Derrick residual `2.52e-8`, radial finite-box
  `omega_0^2=0.311400567`, charged `l=0` bound
  `omega_0^2=0.696506101<0.81`, and the tachyonic control `-1.430515810`.
- **Differential-cohomology definition:** normalize integral generators
  `omega_3` and `omega_2` on `S3` and `S2`.  The mixed curvature
  `n omega_3 wedge omega_2` defines a degree-five differential character,
  whose holonomy on the pushed-forward spacetime four-cycle is meaningful
  without choosing a five-dimensional spacetime extension.  Since
  `H^4(S3 x S2;R/Z)=0`, this bosonic lift is unique at fixed curvature/class.
- **APS/spin refinement:** stable splitting gives
  `Omega_4^Spin(S3 x S2)=Z + Z2 + Z2`.  The reduced generators are detected by
  a regular-point inverse-image spin one-manifold for the `S3` projection and
  the Arf invariant of the inverse-image spin surface for the `S2` projection.
  The globally defined spin action therefore also contains two signs
  `(epsilon_3,epsilon_2)`.  A UV Dai--Freed/APS determinant must determine
  them; the differential form alone cannot.  The global/UV card passes
  `38/38` algebraic and normalization checks while retaining those false
  gates.
- **`SU(3)` representation no-go:** adjoint breaking leaves the faithful
  group `[SU(2) x U(1)]/Z2`, so every irrep obeys `2j+q=0 mod 2`.  The original
  colour-singlet `1_(+1)` scalar and the resulting
  `qq(phi^dagger)^2 -> O(2)` dressing cannot descend from pure `SU(3)`.  A
  concrete nearest variant uses two `bar(3)` scalars,
  `bar(3)->2_(-1)+1_(+2)`, with adjoint mass splitting, plus vectorlike
  fundamentals `3->2_(+1)+1_(-2)` whose singlet partners are made heavy by a
  displayed mass/adjoint-Yukawa tuning.  It really decouples at tree level,
  but `qq phi^dagger` is a flavour doublet/O(1), not the Route-E triplet/O(2).
  Threshold, radiative, finite-monopole, and exact-emergent-2-group gates stay
  open.  The charge-two near-miss also inherits a global-form gate: its
  covering-`U(1)` level normalization is conditional until determinant lines
  and first Chern classes are matched on faithful `U(2)` bundles with
  correlated centre flux;
  `u2_quotient_global_bundle_normalization_proven=false`.
- **Moduli-space `N=2` SQM:** a half-BPS non-Abelian-vortex mother theory
  independently supplies `CP1` orientational moduli and a physical tangent
  fermion through BPS collective-coordinate supersymmetry.  Its chart
  covariance and finite `L2`
  metric are checked.  In the declared canonical Spin-c re-quantization,
  quantization gives `Omega^(0,*)` and one untwisted ground state; the
  source-selected half-form ordering has none.  A separate `E=O(2)`
  WZ/magnetic or Fermi-index line gives
  three positive-chirality states and the AP-E4 paired spectrum
  `lambda_(n,+/-)=+/-2 sqrt(n(n+3))/sqrt(C)`, `n>=1`, with
  `Delta_D=4/sqrt(C)` and `Delta_H=8/C`.  The SQM card passes `27/27`.
  This is not yet a
  completion of the charged-two-colour branch: the BPS vortex and its bulk
  supersymmetry have not been derived for the same `B=1` soliton.  The AP-E3
  level-two line is composable only after spatial transgression of the
  degree-five character and that same-moduli pullback theorem.
  The mother model does not select the canonical vacuum line used by the
  three-state theorem.  Its CPT map is also open: fixed-canonical `O(-2)`
  has one negative mode, not the conjugate three-state kernel; the physical
  antibaryon must derive an anti-canonical/effective `O(-4)` polarization.
- **Backup and handoff:** product compactification is retained only as an
  anomaly-polynomial-checked backup.  The degree-one Route-E portal is
  deliberately postponed until either same-model SQM/WZW composition or the
  compactification branch closes.  Every new card retains
  `physics_promotion_allowed=false`.

New canonical artifacts:

- `tex/ap_e3_charged_two_colour_soliton_proxy.tex`, its bibliography,
  `code/scan_ap_e3_charged_two_colour_proxy.py`, and generated proxy cards;
- `tex/ap_e3_nonextendible_wzw_su3_uv_audit.tex`, its bibliography,
  `code/verify_ap_e3_nonextendible_wzw_su3_uv.py`, and generated global-UV
  cards;
- `tex/ap_e4_moduli_space_sqm.tex`, its bibliography,
  `code/verify_ap_e4_moduli_space_sqm.py`, and generated SQM cards.

Ordered continuation:

1. Replace the nonlinear EFT matching ansatz by controlled charged-QC2D
   evidence and diagonalize the full coupled three-dimensional/quantum
   soliton Hessian.
2. Compute the two spin-torsion signs from a chosen microscopic regulator and
   replace pure `SU(3)` by a semisimple embedding compatible with the
   charge-one/exactly-two dressing, or accept the exact no-go.
3. Derive `N=2` collective-coordinate supersymmetry, its canonical vacuum
   line and CPT map, then transgress/gauge-descend the AP-E3 `O(2)` line for
   the same soliton; otherwise execute the anomaly-free product compactification backup.
4. Only then construct and orient the degree-one Route-E portal below every
   retained gap.

## AP-E3 UV and AP-E4 spectral checkpoint (2026-07-15; superseded where noted above)

Status: **exact-cell and canonical-operator mathematics done; mixed-WZW
intermediate completion anomaly-consistent; Route-E physics non-promoting.**

- **Exactly-two theorem:** use bosonic Schwinger partons with canonical CCR
  and impose two independent compact constraints `G_r=N_r-1=0`, `r=1,2`,
  per indivisible unit cell.  The physical Hilbert is exactly
  `C2 tensor C2`; odd, singleton, `(2,0)`, and `(0,2)` sectors do not exist.
  The rank-four result is cutoff-independent.  Rank ten for one total-number
  constraint and rank twenty for parity refer only to the verifier's declared
  bosonic audit truncation `0<=N_r<=2`; they are negative controls, not an
  additional continuum theorem.  Complete-cell
  positive-parent boundaries obey `Delta>=h`; arbitrary intercell topological
  phases remain an explicit boundary gate.
- **Physical signed orientation:** Coulomb exchange supplies `J_H>0`, and the
  negative electron magnetic moment gives
  `H_Z=+h n.(S_1+S_2)`.  The unique ground state is the anti-aligned triplet,
  with local gap `h`, quotient line `Q tensor Q=O(2)`, and
  `i hbar <Omega_-|d Omega_->=+2 hbar A_+`; hence `k=+2` without a manual sign
  flip.  `verify_ap_e3_exact_two_mixed_wzw.py` passes `26/26`.
- **Mixed-WZW intermediate UV candidate:**
  `SU(2)_c x U(1)_g` with `N_f=2` vectorlike charge-one Dirac flavours and a
  charge-one scalar doublet cancels every dynamical perturbative gauge anomaly
  and the color Witten anomaly.  Gauged `U(1)_g` removes off-diagonal
  Pauli--Guersey currents, while
  `kappa_L=-kappa_R=n_c X_q=2`.  The integral action
  `2 pi hbar 2 integral(omega_3 wedge omega_2)` reduces on `B=+1` to `k=+2`.
  The first local color singlet is `qq(phi^dagger)^2`, simultaneously an
  exactly-two dressing, an `SU(2)_phi` triplet, and `O(2)`.  This construction
  The anomaly ledger includes `SU(2)_c U(1)_g^2=0` and uses a Bardeen scheme
  preserving the dynamical gauge symmetries.  The five-dimensional proof is
  extension-independent on extendible sectors; a differential-character/
  Cech-bordism definition on non-extendible sectors is still open.  The model
  also needs a nonperturbative mesonic phase, positive diquark/PG gaps,
  compact-monopole and bordism audits, a stable unit soliton, and an all-scale
  completion beyond `b_U(1)=6`.
- **AP-E4 tangent theorem:** horizontal projective variation
  `eta=(1-zz^dagger)delta z` is a section of the pullback tangent bundle, and
  the two-chart Jacobian proves `T^(1,0)CP1=O(2)`.  This is a bosonic
  tangent-valued fluctuation theorem, not a fermion theorem.
- **AP-E4 canonical Spin-c spectrum:** for
  `D_T^c=sqrt(2)(dbar_T+dbar_T^dagger)` and AP-E1 `R=1/2`, Hodge/RR gives
  `(dim ker+,dim ker-,index)=(3,0,3)`.  The full nonzero tower is
  `lambda_(n,+/-)=+/-sqrt(n(n+3))/R`, multiplicity `2n+3` for each sign.  The
  first gap is `4`, with five states per sign; all massive modes are paired.
  As an explicit audit anchor, `n=3` is `6 sqrt(2)`, not `6`; the `+/-`
  eigenstates mix `W+` and `W-` and are not themselves chirality eigenstates.
  `verify_ap_e4_tangent_dirac_spectrum.py` passes `22/22` with an independent
  finite-`SU(2)` Casimir diagonalization.
- **Decisive Spin control:** the unique ordinary spin structure has
  `S^+=O(-1)`.  Ordinary spin Dirac twisted only by `T=O(2)` therefore acts on
  `O(1)`, has two zero modes, and gap `2 sqrt(3)`.  Canonical Spin-c on
  `O(2)` is equivalent to ordinary twist `O(3)`; the half-canonical `O(1)`
  shift must be physically derived.
- **Promotion gate:** no bosonic argument creates a target-space fermion.  A
  moduli-space SQM or anomaly-free higher-dimensional compactification, the
  canonical Spin-c determinant, the distinction between automorphism modes
  and matter, a degree-one Route-E portal, and all four-/six-dimensional
  anomaly checks remain open.  `ap_e3_full_uv_closed=false`,
  `ap_e4_physics_closed=false`, and `physics_promotion_allowed=false`.

Current artifacts:

- `tex/ap_e3_exact_two_mixed_wzw_uv.tex` and its bibliography: complete
  exactly-two, physical-sign, gauge-anomaly, two-group/WZW, dressing, and
  deep-UV boundary derivation;
- `code/verify_ap_e3_exact_two_mixed_wzw.py` and
  `output/ap_e3_exact_two_mixed_wzw.{json,md}`: `26/26` deterministic audit;
- `tex/ap_e4_tangent_dirac_spectrum.tex` and its bibliography: complete
  tangent projection, Spin/Spin-c distinction, chirality, full spectrum, gap,
  partner, and anomaly-gate derivation;
- `code/verify_ap_e4_tangent_dirac_spectrum.py` and
  `output/ap_e4_tangent_dirac_spectrum.{json,md}`: `22/22` deterministic
  exact/matrix audit.
- `output/pdf/ap_e3_exact_two_mixed_wzw_uv.pdf` and
  `output/pdf/ap_e4_tangent_dirac_spectrum.pdf`: clean-built, warning-free,
  page-by-page checked 12-page and 9-page notes.

Ordered continuation:

1. Run a charged-two-colour nonperturbative phase/soliton scan over
   `(e_g,v/Lambda_c)`, measuring meson and diquark condensates, `m_PG`, and
   `B=1` lifetime.
2. Complete the faithful-group, discrete-axial, compact-monopole, APS/bordism,
   non-extendible differential-character, and all-scale-embedding audit of the
   mixed-WZW branch.
3. Choose AP-E4's physical realization: derive either moduli-space `N=2` SQM
   or an anomaly-free product compactification; otherwise fall back to the
   ordinary-spin two-mode result.
4. Build and orient a degree-one Route-E portal below the AP-E3/AP-E4 gaps.
5. Continue to AP-E5 only with the selected action: solve the full Q-ball/
   Hopf boundary-value, Hessian, compactness, and emission problem.

## F0-D Another-Physics / Route-E bridge ledger

Status: `in-progress`, deliberately non-promoting.  This lane runs in
parallel with F0-A/B/C and does not block the F1 action freeze.

Current exact results (2026-07-14):

- On the H3+-selected carrier
  \(H^0(\mathbb{CP}^1,T_{\mathbb{CP}^1})\simeq
  H^0(\mathbb{CP}^1,\mathcal O(2))\simeq\mathfrak{sl}_2(\mathbb C)\),
  hence the complex section count is exactly three.  H3 alone still proves
  only `N_fam<=3`; this identity does not dynamically derive H3+.
- For \(q=(a+b\xi+c\xi^2)\partial_\xi\),
  \(\Delta=b^2-4ac\) obeys
  \(B(x,x)=2\Delta=2\sqrt3\,x^TK_{\rm tr}x\) in Route-E spherical
  coordinates.  Therefore two distinct
  centers, a regular semisimple \(\mathfrak{sl}_2\) element, and non-null
  Killing/contact norm are exactly equivalent on this branch.  The
  \(\Delta=0\) boundary of the nonzero theorem domain is the
  nilpotent/double-zero cone; the zero section is a separate excluded orbit.
- Fixed-norm complex doublets give the exact Hopf reduction
  \(S^3/U(1)=\mathbb{CP}^1\).  Its moment map supplies a bounded
  charge/orientation polarity \([-1,1]\); this safely replaces the literal
  positive/negative-energy interpretation, but its physical identification
  is a bridge axiom.
- Route E's rescaling parameter has weight two:
  \(\zeta\mapsto e^{2i\arg y}\zeta\).  Its absolute phase is removable.
  A physical phase requires an independent reference such as
  \(\mathcal R=M_V^{-1}M_C\) and basis-invariant quantities
  \(\operatorname{ImTr}(\mathcal R^n)\) or
  \(\arg\det(I+\mathcal R)\).  A field cannot serve as its own phase
  reference.
- A sextic Q-ball consistently realizes an ``energy bubble'' only
  conditionally.  For
  \(U=m^2f^2-\lambda f^4+\eta f^6/M^2\), existence requires
  \(m^2-\lambda^2M^2/(4\eta)<\omega^2<m^2\).  The benchmark
  \(m=M=\lambda=\eta=1\,{\rm GeV}\), \(Q=10^6\) gives
  \(f_0=1/\sqrt2\,{\rm GeV}\), \(R=12.8424\,{\rm fm}\), and
  \(E/Q=0.872679\,{\rm GeV}\).  Exact U(1) symmetry makes its global phase
  unobservable; a phase portal therefore reopens charge-leakage and lifetime
  gates.
- The candidate messenger action gives
  \(M_R^{\rm eff}=M_*[M_V+\lambda(\Phi)^2K_{\rm tr}]\) and, for
  \(\lambda(\Phi)=g\Phi/\Lambda\),
  \(\zeta_{\rm eff}=g^2\Phi^2/\Lambda^2\).  This is an explicit algebraic
  bridge, not a UV completion: charge assignments, a genuine interaction,
  the `XLH` selection rule, canonical normalization, the full six-by-six
  propagator, loop matching, Q-ball decay, and time-dependent flavor fits all
  remain mandatory.

### AP-E1 projective-doublet checkpoint (2026-07-14)

Status: **geometry closed; physical level selection and stability open.**

- **Local theorem:** one nonzero charge-one complex doublet with fixed norm
  and an auxiliary local common-phase redundancy has
  `S^3/U(1)=CP1`.  Eliminating the auxiliary connection gives the
  Fubini--Study action exactly.  The connection and metric follow directly
  from the vertical/horizontal decomposition of the flat `C^2`
  kinetic term.
- **Global-only no-go:** if the common phase is a physical global symmetry,
  it remains a local field.  Quotienting a constant phase does not give a
  pointwise `CP1`; formally integrating it out gives an additional
  nonlocal transverse term.
- **Fixed-charge theorem:** in the coherent-orientation Q-ball ansatz, Routh
  reduction at `Q=hbar*k` gives a Dirac-monopole rotor on `T*CP1`.  Its
  complete Hilbert space contains all Landau levels.  `H^0(CP1,O(k))`
  requires a separate
  first-order Kähler reduction or a controlled LLL projection with
  `Delta_LLL=2*hbar^2*(k+2)/I`.
- **Route-E obstruction:** the minimal Hopf bundle has Chern number one,
  whereas `T_CP1=O(2)`.  Thus CP1 alone does not derive the Route-E level two.
  If the existing `Q=10^6` Q-ball charge is identified with `k`, the
  holomorphic space has dimension `1,000,001`,
  not three.
- **Recommended branch:** keep macroscopic `U(1)_Q` charge separate from an
  independent level-two projective worldline sector.  This removes the charge
  contradiction, but why `k=2`, why the LLL is isolated, why the triplet is
  chiral-family space, and why the bubble is stable are still explicit gates.
- **Stability boundary:** the pure `3+1`-dimensional two-derivative CP1 model
  fails Derrick scaling.  A finite-`e` Skyrme/Faddeev term can balance the
  scale, but radial unwinding through `Z=0`, portal-induced charge loss, and
  the full fluctuation Hessian must be tested in AP-E5.

### AP-E2 exact regression and AP-E3 microscopic-level checkpoint (2026-07-14)

Status: **AP-E2 done and non-promoting; AP-E3 candidate level magnitude
derived, UV exactly-two and signed-chirality rules open.**

- **AP-E2:** `verify_ap_e2_discriminant_regression.py` passes `30/30`.  It
  freezes the exact projective/Killing/transvectant identities with rational
  arithmetic, 100-decimal complex tests, `SL(2)` covariance, polarized
  bilinears, finite/infinite root charts, the nonzero nilpotent boundary, the
  zero-section exception, the spherical/normalized-basis factor two, and four
  wrong-convention negative controls.  H3+, dynamics, and the Berry level
  remain underived; `physics_promotion_allowed=false`.
- **AP-E3 selected candidate:** declare two orbitals with
  \[
  H_{\rm Mott}=\sum_{r=1}^2\left[\frac U2N_r(N_r-1)-\mu N_r\right],
  \qquad 0<\mu<U,
  \]
  and add `-J_H S_1.S_2`, `J_H>0`.  The displayed Mott window guarantees unit
  occupancy only for the bare onsite Hamiltonian.  For the complete declared
  `H_portal=0` Mott--Hund--orientation model, define
  \[
  C\equiv\mu+\frac h2+\frac{J_H}{4};\qquad
  C<U-\frac{J_H}{8}
  \quad\Longleftrightarrow\quad
  \mu<U-\frac h2-\frac{3J_H}{8}.
  \]
  This is the concise sufficient condition for the interacting `(1,1)`
  plateau.  Unit occupancy on each orbital plus ferromagnetic locking then
  selects the symmetric triplet.  Along the diagonal coherent-state locus,
  \[
  |s;2\rangle=|s\rangle\otimes|s\rangle,
  \quad \mathcal A_2=2\mathcal A_1,
  \quad g_2=2g_1,
  \quad \mathcal L_{\rm ket}=\nu_2^*\mathcal O_{\mathbb{CP}^2}(-1)
       =\mathcal O_{\mathbb{CP}^1}(-2),
  \quad \mathcal L_{\rm pre}=\mathcal L_{\rm ket}^*
       =\mathcal O_{\mathbb{CP}^1}(2).
  \]
  With the fixed convention `A=-i<s|ds>` and the microscopic kinetic sign
  `+i hbar a^dagger dot(a)`, the aligned `-h n.S` model has signed action
  `k=-2`; reversing orientation/coupling gives `k=+2`.  Thus this declared
  dimer derives `|k|=2`, while its dual prequantum line has `c1=+2` and the
  three-state space `H0(CP1,O(2))`.  The `27/27` audit finds
  `spec(-J_H S_1.S_2)=(-1/4,-1/4,-1/4,3/4) J_H`, singlet gap `J_H`, Berry
  residual `2.24e-16`, metric residual `2.22e-16`, and final numerical
  curvature-magnitude estimate `2.000000501994128`.  At `U=4`, `mu=1.5`,
  `J_H=1`, `h=0.2`, the sufficient-condition margin
  `(U-J_H/8)-C=U-mu-h/2-3J_H/8` is `2.025`; the full spinful occupancy sector
  has unique ground `(1,1)` and interacting charge gap `1.85`.  The exact
  large-occupancy lower bound is coercive for `4U>J_H`.
- **Unclosed UV gates:** neither occupancy algebra nor large-gauge invariance
  explains why there are exactly two orbitals.  A permitted singleton retains
  an unwanted `|k|=1` doublet; without Hund locking the target is
  `CP1 x CP1`, and opposite orientations cancel to `k=0`.  The dimer must be
  embedded in an anomaly-consistent four-dimensional theory, odd/singleton
  sectors must be forbidden rather than assumed heavy, and every
  exchange/portal correction must keep the interacting charge gap open, as
  well as remain below the singlet and orientation gaps.  The Route-E portal
  and chiral-family interpretation are also open.  In addition, the declared
  aligned coupling gives `k=-2`; a microscopic orientation/coupling principle
  selecting `k=+2` has not been derived, so signed chirality is a separate
  blocker.  Hence
  `ap_e3_physics_closed=false` and `physics_promotion_allowed=false`.
- **Alternative branches:** keep (i) a filled-fermion determinant with fixed
  filling/signs, (ii) a mixed WZW completion whose quantized coefficient can
  reduce to `k=n_c B` but still needs the correct `n_c=2` UV coset, and
  (iii) the AP-E4 tangent/Dirac route.  None currently supplies an
  exactly-two UV theorem.

Reproducible artifacts:

- `tex/another_physics_route_e_derivation_ledger.tex` and its local BibTeX
  file contain the complete derivations, assumptions, no-go results, and
  rollback trail;
- `output/another_physics_route_e_derivation_ledger.pdf` is the visually
  checked compiled ledger;
- `code/verify_another_physics_route_e_bridge.py` writes
  `output/another_physics_route_e_bridge.{json,md}` and currently passes
  `27/27` algebraic/numerical checks while setting
  `physics_promotion_allowed=false`.
- `tex/ap_e1_projective_doublet_action.tex` and its local bibliography contain
  the complete local/global/fixed-charge derivations, the first-principles
  boundary, the `O(2)` obstruction, and three completion branches;
- `output/ap_e1_projective_doublet_action.pdf` is the clean-built, visually
  checked AP-E1 paper;
- `code/verify_ap_e1_projective_doublet.py` writes
  `output/ap_e1_projective_doublet.{json,md}` and passes `30/30`
  arithmetic/source regression checks while deliberately setting
  `physics_promotion_allowed=false`.  The audit now covers the corrected
  Branch-B symplectic sign, the `k=0` domain exception, flux reversal through
  `|k|`, fixed-charge orientation energy, and critical-source hashes/tokens.
- `code/verify_ap_e2_discriminant_regression.py` writes
  `output/ap_e2_discriminant_regression.{json,md}` and passes `30/30` exact,
  fail-closed regressions with no physics promotion;
- `tex/ap_e2_discriminant_regression.tex` contains the complete AP-E2
  two-chart, Killing, transvectant, basis, and boundary proof;
- `code/verify_ap_e3_level_two_microscopic.py` writes
  `output/ap_e3_level_two_microscopic.{json,md}` and passes `27/27`; its status
  is `ap_e3_hund_pair_derives_abs_k2_conditionally_uv_and_sign_open`, not a
  completed UV theory;
- `tex/ap_e3_level_two_microscopic_origin.tex` and its local bibliography
  contain the complete Mott--Hund, coercivity, Veronese, Berry/Chern,
  alternative-branch, and UV-blocker derivation;
- `output/pdf/ap_e2_discriminant_regression.pdf` and
  `output/pdf/ap_e3_level_two_microscopic_origin.pdf` are the clean-built,
  page-by-page checked AP-E2 (12-page) and AP-E3 (13-page) research notes.

Historical ordered follow-up (fail closed; superseded by the AP-E7 current
authority above):

The labels in this archived sequence record the plan as it stood before the
AP-E4--AP-E7 execution rounds.  They are not the current task numbering or
authorization state.

1. **AP-E1: done at geometry level, non-promoting.**  The local quotient is
   proved, the global-only shortcut is refuted, and the fixed-charge/LLL/O(2)
   distinction is explicit.  The selected continuation is the separated
   macroscopic-charge/level-two branch.  AP-E3 now realizes the magnitude of
   that level in a declared Mott/Hund pair, while the exactly-two UV and signed
   orientation rules remain open.
2. **AP-E2: done, exact and non-promoting.**  Retain the `30/30`
   discriminant/contact suite as a mandatory regression.
3. **AP-E3: candidate-level derivation done; physics open.**  The declared
   Mott/Hund pair derives `|k|=2`; now derive the exactly-two UV field content,
   exclude singleton/odd sectors, select the signed orientation/chirality, and
   construct the anomaly-safe Route-E portal.  Failure keeps the dimer as an
   illustrative EFT only.
4. **AP-E4: was the next active mathematical gate.**  Build the fluctuation/Dirac
   operator, prove whether the physical mode is tangent-valued, and audit all
   chiral zero modes, unwanted partners, and the spectral gap.  This is also
   the independent tangent/Dirac alternative to the AP-E3 dimer.
5. **AP-E5:** solve the radial Q-ball boundary-value problem and fixed-charge
   stability/compactness/emission thresholds.
6. **AP-E6:** construct an anomaly-consistent messenger/contact model and
   pass the complete operator, kinetic, and six-by-six matching gates.
7. **Legacy AP-E7 stage:** quotient all rephasings and isolate genuinely measurable CP
   invariants without assuming away bubble decay.
8. **AP-E8:** compute the proposed adjoint-current correlator and decide
   whether it yields a Majorana contact or only a kinetic Killing tensor.
9. **AP-E9:** rerun Route-E flavor, threshold, proton-decay, cosmology, and
   Floquet gates using one action; only here may promotion be requested.
10. **AP-E10:** keep the gravity/information-density tests independent until
    an explicit operator connects them to AP-E1--9.

## Definition of done

Route F is complete only when one active branch has:

- a fully normalized action and declared field content;
- a stable/metastable symmetry-breaking vacuum with a complete heavy spectrum;
- two-loop running plus threshold matching;
- a physical-basis flavor/seesaw fit with identifiable predictions;
- complete relevant proton-decay amplitudes and uncertainty budgets;
- amplitude-level Ward, crossing, high-energy, and positivity checks;
- a clean-clone command that regenerates all ledgers and papers.

The value of `zeta` need not be derived for the conditional theory to close.
If it remains a fit parameter, it must be labelled as such.  A UV completion is
optional and may remain permanently open.

## Dependency graph

```text
F0-A/B/C -> F1 branch/action freeze
F1 -> (F2 vacuum+spectrum || F3 generator package)
F2 -> F4 two-loop unification
F1 + F2 + F3 + F4 -> F5 flavor+seesaw
F2 + F3 + F4 + F5 -> F6 proton decay
F1 + F2 + F3 -> F7 full amplitudes
F4 + F5 + F6 + F7 -> F8 joint model comparison and predictions

F0-C runs in parallel; for the non-SUSY primary branch it must be translated
to ordinary Weyl mass/kinetic matching.  F9 family-carrier realization and F10
global string completion are optional after F1.  F11 synchronizes papers after
the relevant gates.

F0-D also runs in parallel.  Its exact identities may be imported immediately;
its Q-ball, phase, and H3+ mechanisms may enter F0-B/C or F5 only after their
declared action, symmetry, stability, matching, and phenomenology gates pass.
```

## P0: blockers

### F0-A Evidence recovery and single status registry

Status: `in-progress`.

Re-audit history (2026-07-13): all 19 DYN source scripts were recovered at
`route_E/code_dyn/`; 16 initially replayed with `219/219` internal checks.
The former missing-card/path/ledger statements are now closed by the checkpoint
above and are retained only as discovery history.  Scientific-status defects
remain catalogued in `route_f/CODE_DYN_REAUDIT.md`.

Execution update (2026-07-14):

- **F0-A1 canonical-path/cache subgate: done.**  All 19 recovered scripts use
  `route_E/code_dyn/route_e_paths.py`; the exact case-sensitive `route_E`
  layout works from both the repository root and `route_E/code_dyn/`.
  DYN-9b-1c/1d caches are isolated
  under `ROUTE_E_CACHE_DIR` and keyed by source plus NumPy version.
- **F0-A2 runner subgate: implemented, full closure still in-progress.**
  `run_route_e_dynamics.py` snapshots minimal inputs into an isolated
  workspace, executes a fail-closed DAG, records Git/Python/NumPy/BLAS/SciPy
  provenance and SHA-256 digests, and separates mechanical `all_pass` from
  physics status.  An isolated DYN-0 -> DYN-4a replay passes; a full optional
  dry-run blocks only declared missing dependencies/descendants.  The
  required-card full dry-run is clean; the expensive full numerical replay and
  clean-clone test remain open.
- **F0-A3 registry subgate: done.**  `dyn_claim_registry.json` is the canonical
  status source.  DYN-4a's interval defect is repaired but its fit remains
  conditional; DYN-5 is `invalid_pending_rederivation`; DYN-7 is
  `blocked_missing_branch_thermal_inputs`; DYN-9b-2 and DYN-8 are preliminary,
  while DYN-9b-3 remains blocked on branch-local flavored thermal inputs.
- **Scientific guard evidence:** DYN-4a now uses a two-optimizer candidate fit,
  a converged needle-basin stencil, and connected-local nuisance-minimized
  profiles (`23/23`; no global-basin completeness claim); DYN-5V verifies
  tree Kahler/full-6x6 matching and the selection-rule counterexample (`9/9`);
  DYN-7F classifies the tau-resolved regime and repairs the
  Davidson--Ibarra double-square-root bug (`7/7`).  Passing guard arithmetic
  intentionally does not close the DYN-5 or DYN-7 physics blockers.

Remaining tasks (supersedes the original discovery-stage task list):

- **F0-A1 path/provenance gate: done.**  The exact canonical path is
  `route_E/code_dyn/`; root paths are true delegate wrappers; no lowercase
  alias or manual symlink is part of the supported execution path.
- **F0-A2 clean-run gate: in-progress.**  The fail-closed 21-node DAG and
  required-card dry-run are implemented and clean.  The expensive full
  numerical clean-clone replay, content-addressed output publication, and
  end-to-end digest comparison remain open.
- **F0-A3 scientific-status gate: done.**  The registry is authoritative;
  arithmetic success and physics promotion are separate fields.
- **F0-A4 string-existence gate: done for existence.**  RE-SC3/4/5 are present
  and tracked but remain `unpromoted_pricing_only`; F10 is still required.
- **F0-A5 document synchronization: in-progress.**  Regenerate status tables
  from the registry and remove remaining stale historical claims.

Acceptance:

- direct execution from the canonical source location has zero path aliases or
  manual-symlink requirements;
- clean-clone rebuild has zero missing required cited paths, while optional
  cards are explicitly represented as optional rather than causing a crash;
- every quoted number resolves to one JSON field and one command;
- every numerical assertion fails with a nonzero exit status; "recorded" and
  "disclosed" entries are not counted as assertion passes;
- DYN-4 profile/marginal intervals pass adaptive-resolution and convergence
  tests; DYN-5 includes tree matching and a legitimate loop interaction;
  DYN-7/9b-3 pass a flavored calculation; random-scan hit fractions include
  confidence intervals and prior/measure sensitivity;
- no pair of documents labels the same item both `done` and `open`;
- current evidence audit changes from fail to pass.

### F0-B Repair the three-family theorem

Status: `done` for the logic repair; H3+ physical motivation remains `open`.

Tasks:

- [done] separate invariant contact (H3) from the nondegenerate
  adjoint-trace/Killing-contact strengthening (H3+);
- [done] handle the one-dimensional abelian `B=[1]` counterexample exactly;
- [done] demote the unconditional theorem to `N_fam<=3` and make
  `N_fam=3`, `g=0`, `O(2)`, and the two-center statement H3+-conditional;
- [open] provide a physical or UV derivation of H3+ if it is to be more than a
  transparent selection axiom;
- [retained] keep the minimal-dimension assumption in the `Spin(10):16` theorem.

Acceptance:

- the abelian form `B=[1]` is either admitted (the theorem becomes a bound) or
  excluded by a stated hypothesis with a valid proof;
- an independent symbolic/manual audit checks every genus/orbit case;
- theorem titles and abstracts state the same domain as the proof.

### F0-C Repair Route-B selection and normalization

Status: `in-progress`.

2026-07-14 progress: `audit5_dyn5_model_validity.py` proves the displayed
quadratic action has no cubic messenger tensor, verifies
`delta Z_tree=|zeta|/3`, checks the full six-by-six Schur/Weinberg identity to
better than `1e-12`, and shows both explicitly and for a single additive
Abelian charge that `X L H_u` is allowed.  This is a correct fail-closed
diagnosis, not the missing interacting completion; the latter remains open.

Tasks:

- construct a pre-`B-L` gauge-invariant messenger/source model;
- for the primary non-SUSY branch, replace the literal
  superpotential/Kahler/`U(1)_R` language by ordinary Weyl mass and kinetic
  matrices plus a genuine gauge or discrete symmetry; keep a superspace
  implementation comparison-only unless a SUSY branch is reactivated;
- enumerate all operators through at least dimension six and explicitly forbid
  `X L H_u` and triplet contamination;
- check continuous/discrete anomalies;
- compute the full `N-X` kinetic and mass matching, not only the
  superpotential Schur complement.

Acceptance:

\[
\frac{\left\|[\mathcal M_{NX}^{-1}]_{NN}
-(M_V+\lambda^2K_{\rm tr})^{-1}\right\|_F}
{\max\!\left(\epsilon_{M^{-1}},
\left\|[\mathcal M_{NX}^{-1}]_{NN}\right\|_F,
\left\|(M_V+\lambda^2K_{\rm tr})^{-1}\right\|_F\right)}
<10^{-10},
\]

and, independently,

\[
\frac{\left\|Y_{\nu D}[\mathcal M_{NX}^{-1}]_{NN}Y_{\nu D}^T
-Y'_{\nu D}(M'_R)^{-1}Y_{\nu D}^{\prime T}\right\|_F}
{\max\!\left(\epsilon_W,
\left\|Y_{\nu D}[\mathcal M_{NX}^{-1}]_{NN}Y_{\nu D}^T\right\|_F,
\left\|Y'_{\nu D}(M'_R)^{-1}Y_{\nu D}^{\prime T}\right\|_F\right)}
<10^{-10}.
\]

The floors \(\epsilon_{M^{-1}}\) and \(\epsilon_W\) must carry the same units
as their respective quantities and be fixed before the scan.  Export
\(Z_N\) and verify explicitly
\(M_c=Z_N^{-T/2}(M_V+\lambda^2K_{\rm tr})Z_N^{-1/2}\).

All heavy Takagi singular values must lie below the cutoff and above the EFT's
declared external-energy range; no unintended light sterile state may appear.
Any claimed
`|lambda|^2/(16 pi^2)` wavefunction term must be derived from explicit
trilinear/gauge interactions; the quadratic Gaussian model alone does not
qualify.

### F1 Freeze the active action

Status: `in-progress`.

Compare, using identical conventions:

- F-210: non-SUSY `210_H + 10_H + 120_H + overline{126}_H`;
- F-54: non-SUSY `54_H + 10_H,C + 126_H`;
- minimal `45_H + 126_H + 10_H,C` only as a stress-test branch because of the
  published perturbative/light-doublet problem.

Acceptance:

- write the complete renormalizable action/potential and all normalizations;
- state the breaking chain, light-doublet content, accidental/global
  symmetries, cutoff, and parameter count;
- choose one primary active branch and label all others comparison-only.

### F2 Vacuum and full spectrum

Status: `open`.

Acceptance:

- solve all stationarity equations, not a fixed-ratio slice;
- Hessian Goldstone count equals `dim Spin(10) - dim H`, with each residual
  divided by the declared scalar mass-squared scale and smaller than `1e-10`;
- the null vectors align with the gauge-orbit directions to the same relative
  tolerance;
- all physical scalar masses satisfy the declared stability/metastability
  criterion;
- export every SM irrep, multiplicity, mass, origin block, and uncertainty;
- an independent decomposition/census calculation agrees exactly.

A positive Hessian proves only a local minimum.  Any `metastable` label also
requires named competing vacua and a tunneling/bounce-action estimate; without
that calculation the status remains `local-minimum-only`.

### F3 Exact Spin(10) generator and Clebsch package

Status: `open`.

Acceptance:

\[
\max_{A,B}
\frac{\|[T_A,T_B]-if_{AB}{}^CT_C\|_F}
{\max(1,\|T_A\|_F\|T_B\|_F)}<10^{-12},
\]

\[
\max_A\|T_A-T_A^\dagger\|_F<10^{-12},
\qquad \operatorname{Tr}(T_AT_B)=I_{16}\delta_{AB},
\qquad I_{16}=2
\]

in the adopted long-root-length-squared-two convention (or an explicitly
translated equivalent convention).

P16 candidate rows must then be contracted with the actual root metric,
color/weak projectors, and Fierz identities to produce independent Wilson
coefficients.

The scalar Clebsch package is conditional on the branch selected at F1.  For
the F-210 continuity branch it must include normalized tensors and mixing for

\[
16\otimes16=10_s\oplus120_a\oplus126_s,
\]

with the correct family symmetries for the active
`10 + 120 + overline126` action.  For F-54 or any branch without a `120`,
include exactly the selected scalar representations and do not import the
F-210 flavor content by assumption.  In every case the package must supply
the actual Route-C Branch-S and scalar-mediated proton channels of F1.

## P1: physical closure

### F4 Two-loop running and threshold covariance

Status: `open`.

Acceptance:

- two-loop gauge and relevant Yukawa RGEs across every scale;
- one-loop matching from the actual F2 spectrum;
- perturbativity and Landau-pole checks;
- posterior/profile over threshold nuisance parameters;
- independent reproduction of `M_I`, `M_U`, and `alpha_U`.

### F5 Global flavor, seesaw, and identifiability

Status: `open`.

Acceptance:

- run cited quark/lepton/CKM/PMNS data with covariance to the matching scale;
- simultaneous fit of the chosen Yukawa sector, with `chi2/dof`, pulls, priors,
  parameter count, and Jacobian/Hessian rank;
- distinguish inverse reconstruction from prediction;
- blind the hidden-sector `zeta` calculation before comparison;
- report out-of-fit predictions such as `m_bb`, mass sum, CP phases, and heavy
  neutrino hierarchy with credible intervals.

### F6 Proton-decay closure

Status: `open`.

Acceptance:

- physical gauge/triplet eigenstates and flavor rotations;
- dimension-six gauge/scalar channels, plus dimension five only for any SUSY
  comparison branch;
- short- and long-distance RG, current lattice matrix elements with covariance,
  and current experimental likelihoods;
- channel table at least for `e+ pi0`, `K+ nubar`, and model-leading modes;
- uncertainty decomposition showing the `M_X^4` sensitivity.

### F7 Real amplitude/Ward/positivity audit

Status: `open`.

Start with at least two baryon-violating and two conserving channels, then cover
every symmetry-inequivalent channel class or prove that the remainder follows
by the group action and crossing before declaring amplitude closure.

Acceptance:

- complete tree helicity amplitudes with `s,t,u` crossing and identical-fermion
  signs;
- residues numerically extracted at poles and matched to the action;
- derive the convention-dependent Slavnov--Taylor/Goldstone identity from the
  frozen action, writing its phase/sign/tree normalization as \(c_\phi\), and
  require random on-shell points to satisfy

\[
\frac{\|p_\mu\mathcal M^\mu-c_\phi M_X\mathcal M(\phi)\|}
{\max(\epsilon_{\rm amp},\|p_\mu\mathcal M^\mu\|,
\|c_\phi M_X\mathcal M(\phi)\|)}<10^{-10},
\]

where \(\epsilon_{\rm amp}\) has the same units and is fixed before sampling;

- all unwanted `E^4` and `E^2` coefficients cancel or the branch fails;
- partial-wave eigenvalues and massive-vector crossing/positivity bounds pass
  over a declared energy range, after stating the dispersion assumptions,
  subtracting physical poles, and regulating/subtracting massless SM forward
  exchange.  The cited single-vector bounds are not automatically valid for a
  non-Abelian multiplet with massless SM exchange: derive the needed crossing
  matrices, restrict to a genuinely applicable gapped sector, or record
  `not_applicable` rather than forcing a false branch failure.

### F8 Joint predictions and model comparison

Status: `open`.

Acceptance:

- compare F-210, F-54, and null/contact-free variants with parameter penalties,
  not best-fit `chi2` alone;
- publish identifiable predictions and kill criteria with likelihoods;
- no datum used in a fit may be presented as a prediction.

## P2: optional UV and publication

### F9 Physical origin of the family carrier

Status: `open`.

Two creative options may be developed in parallel:

- a 6D `Spin(10) x U(1)_F` flux model on `P1`;
- a 4D defect/deconstructed index model producing three localized `16`s.

Acceptance: quantized flux/boundary data, anomaly cancellation or inflow,
exactly three complete chiral `16` zero modes, no additional massless chiral
exotics, and all vectorlike pairs lifted above the declared scale.  The
benchmark numerical Dirac spectrum must have three protected zeros and a
nonzero fourth-mode gap throughout a pre-registered bounded deformation
neighborhood; special loci with additional paired zero modes are not excluded
by the index and must be catalogued rather than treated as theorem failures.

### F10 Global string/instanton gate

Status: `open` and optional.

Structural acceptance:

- explicit base, GUT divisor, matter curve, resolved geometry;
- flux quantization, tadpole cancellation, chirality three, no additional
  massless chiral exotics (with vectorlike pairs lifted), and massless
  hypercharge;
- explicit rigid instanton divisor, universal/charged zero modes,
  Freed--Witten/GS consistency, and unwanted-operator veto;
- either match the non-SUSY F1 action, including a consistent high-scale
  SUSY-breaking bridge, or remain explicitly a comparison branch.

Separate precision-`zeta` gate (optional): calculate the Pfaffian and control
multi-instanton/moduli corrections well enough for the precision claimed.  A
structurally successful completion is not failed merely because it leaves
`zeta` as conditional data.

Failure leaves Route D as a checked local interpretation only; it does not
fail the conditional four-dimensional branch.

### F11 Paper synchronization

Status: `open`.

Acceptance:

- theorem Letter contains only repaired kinematic claims;
- dynamics paper contains only evidence rebuilt from the F0 registry;
- UV paper/appendix remains conditional until F10 passes;
- all PDFs build from one clean-clone command and the root `roadmap.md` points to
  the canonical Route-F statuses.

## Immediate next action

Run the 21-node full numerical DAG (19 recovered lanes plus DYN-5V/DYN-7F
guards) from a clean clone with the now-present RE-SC3/4/5 cards, and compare
all source/input/output digests with the current ledger set.  In parallel,
construct an explicit
interacting messenger/selection sector for DYN-5 and supply branch-local
thermal inputs for a two-flavor or density-matrix DYN-7/9b-3 calculation.
F0-B may proceed independently.  Do not spend more compute on new flavor or
proton scans until F1 fixes the action and F2 exports the actual spectrum.
