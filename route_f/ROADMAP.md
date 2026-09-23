# Route F Roadmap: Three-Layer Theory Program

Created: 2026-07-13
Last updated: 2026-09-23

Status values: `open`, `in-progress`, `done`, `failed`, `permanently-open`.
All items start `open` unless marked otherwise.

## Independent Route G, not a new Route-F prerequisite (2026-09-22)

The user selected Scheme A and requested a separate route. Its fixed
M4 x S1 common-spectrum/probe prototype, short derivation and 331/331
implementation checks are at [Route G](../route_g/README.md).
It has no measured-particle fit, chiral family derivation or P54
parameter map. A degenerate joint-dark state is an apparatus limitation,
not a proof of the energy-spectrum paper or a repair of P-PHYS1.

G2 now supplies a tree-level local collision with a localized dynamical
detector and explicit finite recoil (2026-09-23). The
[Route-G roadmap](../route_g/ROADMAP.md) now also records the G2-R
joint output/recoil prediction and mass-sign readout limitation,
plus G2-S finite-resolution forecasts and the user-activated G2-D
external-drive branch with an explicit work ledger.
G2-M fixes a model-unit readout contract, and G2-P audits the missing
physical scale and a candidate SM Higgs-density probe. The latter is
not activated: common mass shifts, added elastic amplitudes and a
coupling-versus-Higgs-width/stability check are recorded explicitly.
No hardware calibration or particle mass identification is inferred.
It does not become a blocker for layers P,
U or S, and does not overwrite their results or promote a Route-F fit.

## Governing decision (2026-08-30; current program authority)

Route F is no longer one serial chain in which a bespoke soliton regulator
must become classical before the four-dimensional theory may advance.  The
program is split into three logically independent layers:

1. **P -- four-dimensional testable mainline.**  Freeze one renormalizable
   non-supersymmetric `Spin(10)` action and compute its vacuum, full spectrum,
   running, flavor, proton decay, amplitudes, and out-of-fit predictions.
   This is the only layer that can close the conditional four-dimensional
   physical theory.
2. **U -- family-UV line.**  Seek a microscopic origin of three complete
   chiral `16`s.  This line is required only to upgrade the family carrier
   from a transparent geometric axiom to a first-principles result.  It does
   not block a correctly labelled conditional four-dimensional theory.
3. **S -- independent soliton mathematics.**  AP-E1--AP-E18, the AP-E11
   complete-minor action, GTA, the relaxation defect `mu`, classicality,
   Callias lines, and related Hessian questions remain valid standalone
   mathematical research.  They do not block P or U and cannot promote a
   four-dimensional or family-UV claim by themselves.

This section overrides every historical use below of `single mainline`,
`portal remains last`, or `Hessian/determinant/classicality blocks physics`.
Those phrases are retained only as records of the earlier AP-E program.  In
particular,

```text
GTA -> S-layer classicality only
GTA -/-> P action/vacuum/RG/flavor/proton/amplitudes
GTA -/-> U Spin^c index/anomaly/zero-mode construction
```

No new AP-E lattice scan is authorized.  A future soliton--family portal is
optional and starts a new hypothesis class; it is not part of the dependency
graph below.  Equality of `CP1`, `O(2)`, or `c1=2` is never sufficient to
identify two physical sectors.

The fail-closed program flags are now

```text
gta_blocks_layer_P = false
gta_blocks_layer_U = false
soliton_portal_required = false
soliton_results_promote_4d_physics = false
layer_P_conditional_closure_independent_of_U = true
```

## Layer P -- four-dimensional testable mainline

### P-PHYS1: physics-first low-order feasibility screen (`bounded screen done; no joint witness`; 2026-09-22)

**New governing priority, explicitly adopted by the user.** The frozen
P54PQ-v2 benchmark is now **regression-only**. Saving that point is not a
success criterion. Its kernels and exact identities remain valid tests,
not a physical vacuum selected by observations. This instruction supersedes
the next-step priorities in P-MOM1/P-COORD1 and the historical requirement
that all of P3 must close before any proton-decay screening can begin.

The next deliverable is one source-backed low-order feasibility figure
and explicit rejection/continuation conditions for the common two-matrix
flavor structure, gauge unification, seesaw scale and proton decay. No
new quartic scan, local-kernel rescue, coordinate-convergence certificate
or long derivation paper is required. Use one-loop gauge running with
the corrected actual PS parent census, tree matching, common two-matrix
Yukawas and declared neutrino/threshold approximations.

Three distinct outcomes must remain distinct: a tested leading-order
candidate passes specified screens; a tested candidate fails a specified
condition; or a condition has not been tested. A failed local optimizer
is not a whole-model no-go. An allowed effective Higgs overlap is not yet
an action-realized scalar vacuum. A proton amplitude left as a flavor
condition is not a passed lifetime test. Scenario spread is not a
confidence interval or a bound on omitted finite matching.

Full finite matching, Nielsen/BRST physical-pole certification and
higher-order control are **precision/promotion gates**, not prerequisites
for clearly labelled exploratory screening. Correct fields, shared
Yukawa structure, scale conventions, EFT applicability and explicit
uncertainty remain mandatory. Complete vacuum viability, PQ cosmology and
the full physical fit are not silently certified by this screen.

**Bounded deliverable completed.** Three literature-initialized flavor
starts (at most 1000 function evaluations each), one fixed-sigma refinement,
and nine correlated whole-parent threshold cards are saved. All older
source/benchmark files are preserved. The one figure is
`output/figures/p54_physics_first_feasibility.png` (also SVG); the concise
joint report `output/p54_physics_first_screen.md` includes the exact
overlap/scale relations and rejection conditions. No new TeX was added.

- Correct four-parent PS beta `(2/3,26/3,26/3)` gives LO `MI=5.08e13`,
  `MX=MU=1.90e15` GeV, `alphaU^-1=38.22` and
  `sigma=MI/g4(MI)=8.91e13` GeV. Nine whole-parent cards move MX through
  `(1.32--2.62)e15` GeV but preserve MI by D parity. Analytically,
  `log(MI/MZ)=2*pi*(alpha1^-1-.4*alpha3^-1-.6*alpha2^-1)/(44/5)`.
  This is a restricted LO identity, not a full-threshold/all-orders result.
- The same two complex symmetric matrices and four normalized overlaps
  obey `MR=i*2*sqrt(6)*sigma*f_D`. The shape-first optimum narrowly misses
  the down-Yukawa policy filter (-31.97% versus 30%); every other tested
  shape group passes. This small boundary miss is NOT the main conclusion.
  A 35% sensitivity is saved without changing the official gate.
- The substantive tension: the SAME shape allows
  `sigma=[7.04e11,1.80e13]` GeV with Dirac-unit norm caps <=1, whereas
  gauge unification requires 4.94 times the upper endpoint. Even without
  norm caps, overlap normalization bounds this profile at `2.06e13` GeV.
  This rejects that profile, not every texture. The single fixed-sigma
  retry meets the overlap/scale condition but misses charged flavor
  badly (muon Yukawa +80.28%), CKM angles, PMNS angles and the mass-squared
  ratio. Its 1000-evaluation budget was exhausted without more retries.
  **No joint LO witness found; not a global no-go.**
- Gauge-proton baseline condition:
  `3.09879^2 |F_L|^2 + 2.94180^2 |F_R|^2 < 1.92336^2`, using the report's
  Wilson convention and published Super-K/lattice inputs. Actual fitted
  flavor amplitudes, scalar exchange and other channels are NOT computed;
  ellipse interiors are conditional, not certified physical points.

The scale, flavor and joining scripts pass **52/52, 18/18 and 32/32**
algebra/numerical/provenance checks. These are implementation checks, not
physical passes. Full matrices, residuals, seeds, norm-cap sensitivity,
heavy-neutrino masses, omitted-log diagnostics and hashes are retained in
the `p54_scales_lo_screen`, `p54_flavor_lo_screen` and
`p54_physics_first_screen` JSON/Markdown reports.

**Stop/continue decision:** no full matching or pole work for the failed
LO candidates. Any further type-I exploration must target the shared
flavor--seesaw scale with a declared budget. If tension persists, compare
the SAME action's type-I+II option with the same f matrix, not a third
free flavor matrix. Actual nondegenerate lower-stage spectra shifting MI
are another sensitivity; repeating these D-parity cards cannot do that.
Neither alternative was executed or shown to solve the tension. After
a common LO witness exists, calculate its proton amplitudes and test a
small corresponding parent-mass/overlap target against a stable action.
Only still-promising cases earn full matching and physical-pole work.

Scope: tree sum rules at MI are a proxy, not completed PS Yukawa running.
Sequential sterile thresholds, neutrino RG, finite matching, fixed
intermediate-vector splitting logs, scalar-vacuum realization, type-I
dominance, collective perturbativity and PQ cosmology remain untested.
The fixed-scale retry has one `MN/MI=1.34`; its common sharp seesaw
threshold is approximate. The 10%/30% charged and 10%/0.12-rad CKM
tolerances are screening policies, not statistical errors or freely
adjustable finite shifts. See the report for exact rejection boundaries.

### P-COORD1: coordinate invariance and threshold-aware convergence (`bounded audit done; no physical promotion`; 2026-09-17)

This is an auxiliary test of the user's pi-coordinate suggestion, not a
replacement of P-MOM1's physical mixing/Nielsen priority. It separates
three different meanings of convergence: linear-solver conditioning,
momentum-series convergence, and quantum loop control. No quartic,
background, finite counterterm or matching input is changed.

- For any invertible local linear map `delta r=J q`, the pair transforms
  as `(Z,G)->(J^T Z J,J^T G J)`. Thus `det(Z-lambda G)` changes only by
  `det(J)^2`: the generalized loop/tree eigenvalues are invariant.
  Pi rescaling does not lower the previous **2.46659** insertion.
  Source responses and all six saved radial kernel spectra also agree.
  A nonlinear log chart requires its gradient/connection term off shell;
  it cannot change the inertia of a stationary covariant Hessian.
- Tree-metric whitening improves the condition number of the known
  positive local matrix `M=G+Z` from **6.90791 to 3.43569**. A controlled
  Richardson solve has 16-step energy error `1.671e-3 -> 1.340e-5`.
  Pure pi rescaling changes neither; relative-VEV normalization actually
  worsens the condition number to **50.58**. These are solver results,
  not physical perturbative-control results. Exact whitening of M gives
  condition one but leaves the relative loop/tree spectrum unchanged.
- A genuine momentum-series improvement is proved for the equal-mass
  master `F(z)=integral_0^1 log(1+z*x*(1-x)) dx`, `z=p_E^2/m^2`.
  Its Taylor radius is four. The physical threshold fixes
  `w=(sqrt(1+z/4)-1)/(sqrt(1+z/4)+1)`, with
  `F=sum_{n>=1} [8n/(4n^2-1)] w^n` and certified positive tail
  `R_N <= [8(N+1)/(4(N+1)^2-1)] w^(N+1)/(1-w)` for Euclidean `0<=w<1`.
  At the actual lightest-vector mass and `p_E^2=.05 omega^2`,
  `z=9.87154`, `w=.301243`: eight-term Taylor error is **34.63**, whereas
  the mapped error is **6.255e-6**, below the proven **6.517e-6** bound.
  Six original momentum probes and the unchanged P50 mass are tested.
- This is not an all-orders loop theorem, Borel reconstruction, or a new
  mass formula. P-MOM1 already evaluates exact masters; no old value is
  replaced. A resummed local metric can make `lambda/(1+lambda)=.7115`
  look smaller, but independent two-loop and vertex remainders remain
  unknown. Doubly massless logarithms must remain nonlocal.

Artifacts: `code/verify_p54_coordinate_convergence.py`, its JSON/Markdown
reports, `tex/p54_coordinate_spectral_convergence.tex` and PDF. The bounded
verifier passes **183/183** checks without new Hessians or parameter cards.
All three direct input/source hashes match. The six-page PDF compiles
without overfull or undefined-reference warnings and has been checked on
every page; Python compilation and whitespace checks pass.

**Recommended use:** adopt tree-metric whitening for numerical work and
continue reporting invariant spectra. If an analytic series is needed,
extend the threshold map to unequal-mass/mixed masters with an assembled
error bound, keeping massless logs and gauge cancellations exact. This
is optional numerical infrastructure, not a new blocker. The structural
alternative is a retained-light-block/nonlocal EFT with controlled power
counting. Fixed-bare Nielsen/BRST, the complete physical mixing block,
full matching, finite seesaw and physical fit gates remain unchanged.

### P-MOM1: common-prescription radial gauge kernel and exact action certificate (`two subgates done; physical pole/control open`; 2026-09-17)

This supersedes the **next-step instructions** of P-INV1. The action,
original broken background, quartics, common invariant finite mass CT and
matching inputs are unchanged. **No new coupling or vacuum scan** was
performed. The frozen point is still not an accepted physical benchmark.

**Completed and independently checked.**

1. `code/verify_p54_gauge_mixed_momentum.py` and its JSON/Markdown reports
   implement the radial one-loop bosonic background kernel in the explicit
   fixing `F=partial.a+xi*g*(T X).eta`, DR/MSbar. From the same kinetic
   action, the off-diagonal quadratic vertex is `W=-2g T partial X`.
   Vector seagull, vector bubble, scalar-vector bubble, all scalar pairs
   and finite-xi direct ghost terms are now evaluated. Minimal wave
   subtraction retains the dimensional rational terms `+2a`, `-2` and
   `-1/8`; no finite Z is fitted. **99/99** checks include independent
   convergent momentum integrals, the massless angular shell identity,
   dimensional gamma-function limits, tree gauge Ward identity, original
   CW/tadpole/common-CT recovery and all six prior scalar-plus-static-vector
   kernels. Ten cached action jets are reused; **zero new full Hessians**.
2. Finite-xi scalar gauge-fixing masses and their derivatives are included
   together with longitudinal vectors and ghosts. At `p_E^2=.05 omega^2`,
   `||K_xi-K_0||_F` falls from `.99037` at `xi=.1` to `9.7311e-5` at
   `xi=1e-6`. Direct scalar ghosts vanish in this limit. The original
   finite CT is held fixed and residual finite-xi tadpoles are saved.
   **This is a Landau-limit check, not a fixed-bare Nielsen certificate.**
3. `code/verify_p54_exact_transverse_certificate.py` and its JSON/Markdown
   reports give **109/109 exact-zero checks**, including all 58
   unit-invariant/projector identities on all 328 output covectors.
   Gaussian-integer five-form contractions at a rationalized background
   certify the actual action, including beta/chi4/eta1 mixed blocks;
   no floating Hessian/cache/tolerance enters the proof. Block
   multihomogeneity continues it to the fixed-orientation radial family.
   The maximum conservative integer-operation bound is 1,008,000, far
   below `2^62`. Thus the P50/P36 action identities are no longer merely
   numerically supported. The positive-hard-spectrum and scalar-one-loop
   qualifications of the inverse-gap theorem remain.

**What this changes physically, and what it does not.**

- The local gauge insertion has generalized eigenvalues
  `(-.1475756,0,.0003448)`. Adding it to the old hard scalar insertion
  changes the largest eigenvalue `2.61032 -> 2.46659`, not to a controlled
  small correction. This diagnostic excludes the nonlocal scalar
  soft-soft part and unfitted fermions; it is not a physical residue.
- At `p_E^2=.05 omega^2`, the bosonic radial Euclidean kernel has
  generalized eigenvalues `(-.589137,.142534,2.188925) omega^2`.
  The smallest was `-.584214` with only scalar momentum and static vectors.
  These are off-shell background eigenvalues, **not particle masses or
  a gauge-independent instability proof**.
- The exact action theorem now justifies `m50^2=12 tau sigma^2` and
  `m36^2=(58/5) tau sigma^2`, while the fixed-action fluctuation vertex
  retains `T_w=.36`. The previously derived scalar bound
  `rho_Z >= .01471609503/tau` therefore has an exact action foundation.
  It is not a bound on the signed full gauge/fermion kernel.

Derivation: `tex/p54_same_prescription_gauge_momentum.tex`, incorporating
`tex/p54_exact_transverse_certificate_fragment.tex`, and the PDF in
`output/pdf/`. The two new verifiers total **208/208** checks; numerical
regressions and exact integer certificates are explicitly distinguished.
All **15/15** recorded source hashes match the current files. The ten-page
PDF compiles without overfull/undefined-reference warnings and has been
visually checked on every page. Python compilation and whitespace checks
pass. No new commit/push was requested in this step.

**Next bounded choices; no blind quartic weakening.**

1. Prioritize the fixed-bare Nielsen/BRST insertions and the full neutral
   physical mixing block, including the background-to-quantum map and
   background displacement. Reapply the loop light-doublet condition in
   the same tadpole prescription. Do not identify a zero of this radial
   background block alone with a physical mass pole.
2. Alternatively construct a retained-light-block EFT with explicit
   momentum power counting from the certified vertex-to-gap criterion.
   A large finite field renormalization or integrating out a nonzero but
   small gap is not a proof of perturbative control. Any replacement
   requires mixed BFB, competing-vacuum and higher-order remainder tests.
3. The existing complete radial Majorana function can be added for three
   specified Takagi masses; no fitted masses are invented. Nonradial
   fermion kernels still need actual complex Yukawas. Full Wilson/box,
   lower gauge/ghost/Yukawa finite matching and the finite `C_HN` seesaw
   remain open. **Physical fit, determinant and portal are not promoted.**

### P-INV1: invariant attribution, inverse-gap obstruction and radial fermions (`scalar audit done; no controlled replacement`; 2026-09-14)

This supersedes the **next-step instructions** of P-IR1, not its results.
The unchanged P54PQ-v2 action has now been decomposed invariant by
invariant, on the original real-coupling radial slice. No observed mass
or flavor objective was fitted. No default parameter or old matching
trajectory was overwritten. The phrase original point refers to the
broken background `r=(1,.1265,.25)`, not the symmetric field origin.

Artifacts:

- `code/verify_p54_invariant_decomposition.py`, JSON/Markdown: **445/445**.
  All 29 coefficient directions, full ordered interference tensors,
  seagull/common-CT pieces and signed allocations are retained.
- `code/verify_p54_nonuniform_reselection.py`, JSON/Markdown: **183/183**.
  Six first-stage nonuniform cards, own spectra and tadpoles, necessary
  pure-H BFB gates, and independent final-action jet checks.
- `code/verify_p54_joint_portal_probes.py`, JSON/Markdown: **36/36**.
  Two explicitly declared second-stage joint mixed-coupling cards;
  zero new Hessians, not an adaptive optimizer.
- `code/verify_p54_transverse_gap_obstruction.py` and its reports:
  **112/112**, the two reducing scalar blocks and conditional inverse-gap
  bound. Both explicit sparse integer projectors and all 29 unit-action
  residuals are retained, with 27 cache-content hashes.
- `code/verify_p54_radial_fermion_momentum.py`, JSON/Markdown: **56/56**;
  `notes/p54_gauge_fermion_momentum_audit.md`.
- Full derivation: `tex/p54_invariant_attribution_nonuniform_control.tex`
  and PDF under `output/pdf/`.

The five new verifiers total **832/832** passing algebra/regression
checks; that includes correctly rejecting all eight engineering cards.
All **49/49** recorded source-hash entries match current files. The
eight-page PDF was compiled without overfull/undefined-reference warnings
and visually checked on every page. Python compilation and
`git diff --check` pass. No commit or push was requested in this step.

**What changed in our understanding.**

1. Multihomogeneity reconstructs every radial `H,T,U` block from a single
   unit-coefficient Hessian. The first build used 27 unit Hessians plus
   one new independent full-action radial point; `xi1,xi2` have algebraic
   zero radial jets. This is an exact polynomial method, not a fitted
   mass formula. It includes real/imaginary field coordinates, but new
   imaginary coupling directions still need their own basis.
2. At the original point, `lambda2,lambda4,lambda4p` carry **94.11%** of
   the leading scalar kinetic insertion and **96.13%** of the negative
   canonical-sigma scalar-plus-CT curvature in the declared equal-share
   interference ledger. The kinetic self terms sum to only **0.9071**,
   versus the full **2.6103**: interference cannot be dropped. These
   basis-dependent allocations are not separate observables or total
   parameter derivatives. Mass terms have zero vertex allocations but
   still determine propagators.
3. The three transverse quartics do not stiffen the tree radial
   potential, whereas `H_tree[sigma,sigma]=2 lambda0 sigma^2`. Nevertheless
   reducing them alone is not a repair: the small propagator gap exposes
   the still-large mixed-field vertices. From `tau=.1` to `.03`, the
   minimum hard gap falls `.00883111 -> .00267622 omega^2`, the leading
   canonical-Phi fraction rises `.667676 -> .824065`, and `rho_Z` rises
   `.713427 -> 1.099954`.
4. There is a conditional analytic explanation, not just a scan:
   pure-Sigma reducing blocks of ranks 50 and 36 have
   `m^2=(12,11.6) tau sigma^2`, while their fixed-parameter `T_w=.36`
   stays nonzero. The stationary mass counterterm relation cancels
   static portal contributions, **not their fluctuation derivatives**.
   Positive scalar-pair weights imply `rho_Z >= .01294651897/tau`;
   a fixed three-component witness strengthens it to
   `rho_Z >= .01471609503/tau`. Thus `tau<.04905365` fails the `.3`
   scalar engineering margin on this ray under those identities.
   Integer-over-eight projector algebra is exact; the action-to-block
   identities have strong numerical verification, **not yet a symbolic
   or interval certificate**. This is not an unconditional whole-model
   exclusion, nor a bound on the complete gauged pole kernel.
5. All **eight** nonuniform cards fail control. In the second stage,
   `lambda0=.2`, transverse quartics/chi4 multiplied by `.01`, and
   alpha/beta/chi2 by `.1` give, at `sigma=.5`, a positive lowest hard
   radial curvature `.0229910 omega^2`, but `rho_H=.538986` and
   `rho_Z=1.89903`. Positive curvature alone remains insufficient.
6. A structural blind spot is now explicit: changing `xi1` to a negative
   value leaves all these radial jets unchanged but makes the classical
   potential unbounded on a complex-null H direction. Necessary gates
   `xi1>=0`, `xi1+xi2>=0` are implemented. Complete mixed-field BFB,
   flat-direction lower-degree terms and competing vacua remain open.
7. The radial fermion subset no longer needs a global family fit:
   with three input Takagi masses `M_i=sigma f_i`, only the sigma-sigma
   entry is nonzero,
   `Pi_F^FV(s)=-sum[(M_i/sigma)^2(4 M_i^2+s)L(M_i^2,M_i^2;s)]/(16 pi^2)`.
   Its finite kinetic derivative and old Majorana cap are independently
   verified. Synthetic mass cards are not fitted physical families.
   The full nonradial fermion kernel still needs the complex Yukawas.
8. The explicit background gauge fixing produces a derivative mixed
   vertex `W=-2g T^a partial_mu X`. Its scalar-vector bubble cannot be
   recovered by differentiating the constant-background CW spectrum.
   Actual vector mass first/second vertices are now saved, but the
   finite vector/mixed integrals, dimensional rational terms and
   finite-xi/Nielsen consistency remain unfinished. Direct scalar ghost
   vertices may vanish in this strict Landau fixing; that does not close
   the general ghost/matching audit.

**Next priority -- bounded, not another blind scalar scan.**

1. Complete the same-prescription radial vector/scalar-vector momentum
   integrals with the already constructed vertices, common counterterms
   and finite-xi/BRST/Nielsen checks; combine with the existing scalar
   and exact radial fermion subsets. Determine the relevant physical
   block before interpreting poles or residues.
2. In parallel, obtain a symbolic/interval certificate for the two
   reducing action blocks. Use their vertex-to-gap criterion to design
   any further change; alternatively test an EFT retaining the light
   blocks explicitly with consistent momentum power counting. Merely
   weakening tau further, absorbing a large finite Z, or integrating out
   these blocks without checking the expansion does not prove control.
3. Any future candidate also needs mixed BFB/global-vacuum comparison,
   full nonradial stability, loop light-doublet retuning and its own
   running/thresholds. **No replacement benchmark or physical fit is
   promoted.** Full Wilson/box, lower gauge/ghost/Yukawa matching and the
   complete finite seesaw including `C_HN` remain open.

### P-IR1: actual Goldstone/IR audit, scalar momentum and bounded reselection (`scalar subgate done; no controlled replacement point selected`; 2026-09-13)

The requested first test was performed at the **unchanged original P54
action and scalar point**, with its common fixed-VEV DR/MSbar Landau
prescription. Only after a quantitative control failure was established
were four separate engineering parameter cards tested. The model's
operators/family content, default parameters and old matching trajectories
are not modified. This checkpoint does not close full gauged background
stability, physical pole masses or the finite matching chain.

Artifacts:

- `code/verify_p54_goldstone_ir.py` and its JSON/Markdown reports:
  **65/65** checks;
- `code/verify_p54_scalar_reselection.py` and its JSON/Markdown reports:
  **36/36** checks;
- detailed derivation: `tex/p54_goldstone_ir_background_control.tex`
  and its PDF under `output/pdf/`.

Both verifiers reuse the ten existing polynomial Hessian cache entries;
**zero new full Hessian evaluations** are needed. Scalar momentum
integrals use actual degenerate mass blocks, with independent integral
and quadrature checks. These are not fitted mass formulas. Test success
includes correctly rejecting physical promotion.

The six-page PDF was compiled and visually checked on all pages, with no
remaining overfull-box or unresolved-reference warnings. All 15 source-hash
entries in the two reports match current files; `git diff --check` passes.

**What the Goldstone test closes and what it does not:**

1. For an invariant functional, `H T X = T grad V`. Applying this to
   the same hard functional and common counterterms fixes the entire
   34-dimensional Goldstone shift: the two contributions cancel, with
   matrix residual `7.8e-15 omega^2`. The actual radial cubic tensors
   independently satisfy the differentiated Ward identity. A positive
   Goldstone mass cannot be introduced as a free repair parameter.
2. The soft one-loop tadpoles vanish as `epsilon log(epsilon)`, but the
   potential Hessian contains `B log(epsilon/mu0^2)` with
   `B_rs=Tr[(Q^T T_r Q)(Q^T T_s Q)]/(32 pi^2)`, a positive Gram matrix.
   Its projection on the old negative witness is `0.000313131 omega^2`;
   this does **not** say that every component of B is small. Restoring
   these modes does not give a finite positive mass correction. The
   zero-momentum full Hessian remains IR singular.
3. The complete one-loop **scalar** momentum dependence on the radial
   slice is now evaluated from cubic/quartic vertices, including all
   hard/soft pairs. At nonzero Euclidean momentum the singular log is
   replaced by `log(pE^2/mu0^2)-2`; regulator cancellation is explicit.
   This is not the complete gauge/ghost/fermion pole kernel. The old
   zero-momentum vector term is retained only in labelled diagnostics.
   Neither their low-momentum negative eigenvalues nor high-momentum
   positive eigenvalues decide the physical pole/background question.
4. The hard scalar kinetic correction relative to the tree radial metric
   has eigenvalues `(0.00900848,0.27686405,2.61031989)`. Thus the original
   scalar small-insertion expansion fails in this prescription even
   apart from the zero-momentum IR problem. Exact canonical normalization
   alone does not bound omitted orders or change a curvature's inertia.
   Full two-loop Goldstone resummation is not declared complete.

**A repair direction is now ruled out analytically, not by scanning.**
For `V0 -> z V0` with all existing scalar coefficients scaled, but the
old radii, gauge coupling and matching scale fixed,

```text
H_hard(z) = z H_tree + z^2 (S + log(z) L_S) + V
deltaZ_scalar,h(z) = z deltaZ_scalar,h(1).
```

The sigma coordinate obeys, for every `0<z<=1`,
`H_sigma,sigma(z)/omega^2 <= 0.02880405 z - 1.54130008 z^2 - 0.00402215`
`<= -0.00388757 < 0`. The omitted logarithmic term is nonpositive.
This excludes the **entire bosonic uniform-weakening ray**, not general
scalar choices or fermionic corrections at new points.

**Actual bounded reselection attempt -- no accepted point.**
Four coarse `(z,sigma/omega)` choices were tested: `(.03,.3)`,
`(.03,.5)`, `(.1,.3)`, `(.1,.5)`, with `w=1`, `vs=.25`. Each has its
own re-solved dependent mass parameters, tree light-doublet condition,
full ambient tree spectrum, actual Goldstone/SM Casimir classification,
new finite tadpoles and hard radial/kinetic matrices. No observation
enters the selection. Three have a negative hard radial eigenvalue.
The fourth has minimum `+0.00971366 omega^2`, but its maximum relative
radial insertion is `1.4323` and hard kinetic correction `0.30988`, so
it also fails control. The declared `0.3` engineering margin is not a
convergence theorem; the fourth point already fails the weaker mass
insertion criterion `<1`. Its positive spectrum is not promoted.
All parameter cards and rejections are saved; no replacement physical
benchmark is selected. Changing sigma would also require recomputing,
not reusing, the historical P2 unification/threshold solution.

**Revised next options:**

1. Before another parameter search, decompose the independent scalar
   invariants' contributions to the basis-invariant relative mass and
   kinetic corrections, particularly the independent 126 quartics.
   Test nonuniform coupling changes against both quantities, not just
   the signs of selected masses or the sizes of bare couplings.
2. To decide the physical fate of the original point, complete the
   same-prescription vector/Goldstone/ghost mixed and fermionic momentum
   maps with Nielsen consistency. The computed scalar subset alone
   is insufficient; the original point has not been excluded to all orders.
3. Only a controlled common point can feed a recomputed full nonradial
   response, loop light-doublet condition, upper PS saddle and running/
   thresholds. Full Wilson/box, lower gauge/ghost/Yukawa matching and
   finite CHN seesaw remain open. Their formal work need not wait for
   physical promotion, but physical flavor fit remains disabled.

### P-REN1: same-bare transport, radial obstruction and CHN operator closure (`scheme algebra done; frozen hard benchmark fails stability; full matching open`; 2026-09-12)

This is the current checkpoint for the requested order **common parameters
and tadpoles before complete Wilson/box, lower matching and finite seesaw**.
It changes the assessment of the existing scalar point without changing
the P54PQ-v2 action, scalar parameters or two UV family matrices. The
following are bounded completed calculations, not completion of that chain:

- `code/verify_p54_common_renormalization.py`: **46/46** checks;
- `code/verify_p54_chn_threshold_closure.py`: **20/20** checks;
- strengthened `code/verify_p54_finite_yukawa_interface.py`: **18/18**;
- existing upper-vector matching regression/source inventory: **38/38**;
- complete derivations: `tex/p54_common_renormalization_radial_bound_chn.tex`
  and its PDF in `output/pdf/`; numerical matrices, checks and source
  hashes are in the corresponding JSON/Markdown reports.

The seven-page PDF was compiled and visually checked on every page;
no overfull-box or unresolved-reference warnings remain. All 24 recorded
source-hash entries in these four reports match current files. The 122
checks include successful detection of a failed physical stability
condition; they are not 122 passed physical acceptance conditions.

**Common parameter/tadpole transport -- algebra closed, physical control open.**
The dependent mass parameters obey `p_MS = p_fixed + delta_p` at one loop
for the same bare action. The finite shifts, in units of omega squared,
are `(-0.05633071, -5.46974409, -0.03423771)` for `(mu2,nu2,mus2)`.
The complex-126 factor of one half is required and passes independent
radial/phase normalization checks. Removing its large counterterm at the
old numerical input is **not** a scheme conversion. Tree backgrounds,
mass kernels and Wilson responses must also transform: if
`q=f(p)+Delta`, `p'=p+a`, `q'=q+b`, then
`Delta'=Delta+b-Df*a`. The actual upper scalar-exchange response obeys
this identity, including the induced tree-valley displacement.
The broken MS tree displacement and explicit-loop displacement cancel
formally at first order. Their large individual coefficients do not
justify a small finite displacement or a resummed physical MS benchmark.
Tiny formal loop-parameter continuations test this identity only.

**New physical working-point obstruction -- not a decimal-level fit issue.**
The previously untested physical radial block has metric
`G=diag(12/5,2,1)`. The bosonic hard-plus-common-counterterm generalized
mass squares are `(-0.758448, 0.094119, 3.093105) omega^2`.
The negative direction is orthogonal to the actual gauge orbit and PQ
phase; doublet retuning has zero radial projection. Full spectral finite
differences at two steps verify the radial Fréchet Hessian. Large scalar
multiplicities dominate the effect; no spectrum formula is fitted.

For the actual radial fermion tensor, only the three masses
`M_i=sigma*f_i` are nonzero. Their fixed-VEV one-loop radial correction is

```text
F_sigma,sigma = -4 sum_i M_i^4 log(M_i^2/mu0^2)/(16 pi^2 sigma^2)
             <= 6 mu0^4/(e 16 pi^2 sigma^2)
             = 0.08751124 omega^2.
```

The bound is analytical for every three-family Yukawa matrix, not a
sampling bound. Even adding the maximal positive correction leaves a
negative generalized eigenvalue `-0.7146953 omega^2`. Thus a subsequent
family fit cannot rescue **this fixed-point one-loop hard-potential
truncation**. This is not an all-orders, soft-IR-complete or
gauge-independent pole-mass exclusion of the Spin(10) theory. The earlier
doublet-only positivity and the 6.76 upper insertion warning remain
historical subset results, not whole-background stability certificates.

**CHN finite subset and sequential operator closure.**
For `M(rho)=M-2 C_HN rho`, the direct linear-CHN Yukawa bubble gives
`deltaY = 2 Y diag[M_i(1+log(mu^2/M_i^2))] C_HN/(16 pi^2)` in a positive
Takagi chart. Its terminal C5 contribution and finite Higgs-mass tadpole
are calculated at leading `qH=0`; integral, RG, determinant and complex
family-covariance checks pass. These are not the complete finite seesaw
coefficients. At a partial threshold, exact field-dependent elimination
also gives `D6=2 Y_h M_h^-1 C_hr` multiplying
`-D6 (LH) N_r (Hdagger H)+h.c.`. Its mass down-mixing
`2(D6 M_r Y_r^T + Y_r M_r D6^T)/(16 pi^2)` closes a specific mixed-block
RG residual. The retained-line bubble belongs to the EFT subtraction,
not a second hard matching contribution.

Existing synthetic-family trajectories develop nonzero mixed CHN
Takagi entries despite input alignment with MR; their first two D6
norms are nonzero. No new terms were silently applied to these old
trajectories. A further controlled expansion in `M_r/M_h` can place this
feedback beyond leading dimension-five accuracy. Without that additional
hierarchy approximation, retained-mass finite matching needs the generated
operator or a remainder bound. No physical fitted-trajectory bound or
full dimension-six anomalous-dimension system is claimed.

**Remaining order and explicit acceptance gates:**

1. Keep the action and resolve background/control consistently: organize
   the scalar/Goldstone IR contribution and momentum/gauge treatment in
   the same parameter prescription. Do not delete the counterterm or
   delegate the failed radial block to flavor fitting. If this point
   cannot support a controlled background, a separately declared,
   consistently reselected scalar point is needed; none is chosen here.
2. Complete the source-dependent Wilson/box and lower gauge/ghost/Yukawa
   functional in that common convention. Formal diagram/operator work
   can proceed independently of physical background promotion; its
   algebraic completion must not be mistaken for a viable spectrum.
3. Close finite seesaw using either an enlarged basis with retained-mass
   feedback or an explicitly justified hierarchy expansion. Correlated
   finite MR, lambda, nonzero-qH, wave-function, pre-existing-C5 and mixed
   diagrams remain open. Only then return to the full light state/fit.

The shared fit consumer now additionally requires
`common_parameter_tadpole_consistency`, `full_background_stability`,
`mass_operator_control`, and `CHN_operator_closure`; a payload containing
all the older flags but missing these is rejected. Physical fit remains
disabled. No new vacuum/lattice scan, UV parameter, action change, commit
or push is made in this checkpoint. Only three new mixed polynomial
radial Hessians were required; subsequent verification reuses their cache.

### P-MATCH-SCALAR2: common tadpole anchor, scalar/Wilson subsets and finite seesaw (`subgates done; full matching and control open`; 2026-09-11)

This checkpoint advances the requested upper scalar/tadpole/Wilson,
lower finite matching and sequential seesaw chain. It **does not finish
the entire chain**. The P54 action, scalar parameter point and two UV
family matrices remain unchanged. No new lattice or vacuum scan was
performed. Missing cubic/quartic action vertices were obtained by exact
central differences of the polynomial Hessian, cached separately from
the original read-only Hessian cache; this is not fitting a spectrum.

Artifacts and bounded validation:

- `code/verify_p54_upper_scalar_matching.py`: **14/14**;
- `code/verify_p54_upper_wilson_exchange.py`: **16/16**;
- `code/verify_p54_lower_scalar_finite_subtraction.py`: **22/22**;
- `code/verify_p54_finite_seesaw_matching.py`: **24/24**;
- corresponding JSON/Markdown records in `output/`;
- `tex/p54_scalar_tadpole_wilson_seesaw_matching.tex` and its PDF;
- original default sequential-seesaw evolution rerun: **46/46**, refreshing
  its tree-threshold report/source hash (separate from the finite report).

The nine-page PDF was compiled and visually checked on every page;
no overfull-box or unresolved-reference warnings remain. The five
reported source-hash inventories match the current files (122/122
combined new/default-regression checks).

**New findings that change the next decision:**

1. Eight PS-invariant real bilinears reconstruct the full **248 by 248**
   scalar-cubic hard kinetic and seagull/bubble mass kernels. Both
   heavy-heavy and heavy-soft lines are included at leading retained-mass
   hard order, including the upper Landau Goldstones. Independent phase/
   parent probes, a second polynomial step and all PS Ward commutators
   validate the reconstruction. Scalar kinetic eigenvalues range from
   `0.0000938883` to `0.00450338`. This is not a retained-mass error bound.
2. The historical **broken-vacuum fixed-VEV counterterms are one common
   prescription**, not permission to set both backgrounds' tadpoles to
   zero independently. At their exact common scale
   `mu0/omega=0.5626028899853225`, the upper hard-plus-bosonic-CT radial
   tadpole is `(0.09810775,0.00660187) omega^3`. Its linear heavy-valley
   response is `(-0.01393190,-0.05559093) omega` in `(w,vs)`, about 22%
   of the upper `vs`. This is a factorizable hard matching insertion,
   **not the full quantum upper saddle**: retained EFT tadpoles,
   fermionic fixed-VEV counterterms and the final light retuning are absent.
3. All **24 Yukawa-coupled upper real six modes**, including ten real
   invariant mass mixings, now enter the factorizable scalar-exchange
   response `delta C = delta J^T D^-1 J + J^T D^-1 delta J
   - J^T D^-1 delta D D^-1 J`. Scalar/vector mass diagrams, the common
   bosonic CT and heavy-valley response accompany the previous proper
   Yukawa/vector vertices and fermion legs. Heavy scalar normalization
   cancels between D and J and must not be counted twice.
4. **A substantive control warning:** the computed subset has
   `rho(D^-1/2 deltaD D^-1/2)=6.755926`, with eigenvalue range
   `[-0.131468,6.755926]`. Its inverse-kernel Neumann series therefore
   fails at unit loop parameter, although the resummed subset mass is
   positive. The two synthetic source-projection ratios, `0.0441` and
   `0.3605`, do not certify control of the full propagator. This does not
   exclude the full model, since diagrams, parameter feedback and
   retained-mass corrections are missing. **Do not promote this point
   to a fit or present a selective resummation as a prediction.** Passing
   the algebra tests does not mean this perturbative-control test passes.
   On the dominant normalized mass direction the scalar, vector,
   common bosonic CT and valley contributions are respectively
   `(-0.538946,-0.016987,+7.346374,-0.034515)`: the transported finite
   counterterm dominates, not decimal-level loop-kernel error. These
   are contractions on one common direction, not separate eigenvalues.
5. At the actual lower broken background, the exact 303-dimensional
   scalar quotient determinant factorizes into the 55-dimensional
   heavy block and the nonlocal PS kernel. Its finite MSbar nonlocal-minus-
   local-two-derivative vacuum contribution is
   `-4.652429862e-5 omega^4`, with independently verified scale derivative
   and broken-PS covariance. This constant-background **scalar** result
   is not a full gauge/ghost/source functional. It demonstrates why a
   local kinetic truncation cannot simply replace the nonlocal kernel
   inside a loop integral; the corresponding Wilson subtraction is needed.
6. The canonical terminal type-I finite C5 formula is implemented and
   checked against its full beta difference, but **rejects actual P54
   input** with nonzero CHN/type-II C5. For that input only the universal
   gauge/quartic finite C5 term is now applied at each moving Takagi
   block, retaining CHN in the running/tree map. Two existing synthetic
   families have three live thresholds; relative C5 shifts versus tree
   thresholds are `0.00221120` and `0.00224041`. These are not physical
   fits, neutrino predictions or calibrated matching errors.

**Next order / exact remaining scope:**

1. Resolve the common tadpole/parameter and full mass-operator control
   issue, alongside upper source-dependent scalar three/four-point
   Wilson vertices and nonfactorizable boxes/direct-vector operators.
   The large insertion requires diagnosis of missing contributions and
   the common renormalized parameters. If it persists after completion,
   compare a consistently converted prescription or reject the benchmark;
   do not change the action or fit isolated mass formulas by default.
2. Complete the lower source-dependent hard-minus-EFT finite functional:
   derivatives of the scalar contribution, gauge/ghost/Goldstone package,
   Yukawa/kinetic/contact operators, tree/loop Wilson insertions and
   controlled retained-mass dependence. Do not relabel the constant
   scalar determinant subtraction as full covariant finite matching.
3. Complete finite seesaw with actual CHN and pre-existing C5 insertions,
   mixed removed/retained sterile graphs and correlated finite
   Yukawa/MR/Higgs/wave-function maps. The canonical terminal limit and
   block beta-difference tests are regressions, not substitutes.
4. Only then solve the common-prescription complete complex fermionic
   light eigenstate, single-light-doublet retuning and positive heavy
   block, followed by a constrained two-spurion global fit.

`full_upper_Wilson_matching=false`, `full_lower_matching=false`,
`full_P54_finite_seesaw_matching=false`, `physical_fit_enabled=false`.
The computed bosonic mass subset additionally has
`Neumann_series_converges_at_unit_loop_parameter=false`.
No commit or push was requested or performed in this checkpoint.

### P-MATCH-KJC1: one-loop source transport and upper vector contribution (`subgates done; full matching open`; 2026-09-11)

The next matching step preserves the frozen action and separates **complete
transport algebra** from **incomplete diagram content**. Artifacts:

- `code/verify_p54_upper_vector_matching.py`, `output/p54_upper_vector_matching.{json,md}`;
- `code/verify_p54_kjc_finite_matching.py`, `output/p54_kjc_finite_matching.{json,md}`;
- `tex/p54_source_preserving_finite_matching.tex` and its eight-page PDF.

The actual upper stationary PS background, not the lower broken vacuum,
is used throughout. Its 24 coset vectors have `M_V^2=0.317035178968` in
the historical reference units. The DR/MSbar background-field Landau
calculation closes the **upper heavy-vector active Yukawa vertex,
fermion kinetic metric and vector--scalar active kinetic contribution**
at leading retained-mass order; **38/38** checks pass. The dimensionally
finite numerator constants are retained before subtraction. Actual
Clifford/real-scalar state sums reproduce coset Casimirs `CF=3`,
`CS(H,F,L,R)=(3,7,6,6)` and the complete removed-vector beta difference.
This does not close all gauge-dependent potential or Wilson matching.

New results that must not be repeated or discarded:

1. The vector--scalar metric generates **mirror H/F Yukawa increments**,
   in addition to the earlier pure-Yukawa cubic-spurion mirrors. Their
   coefficients in the fixed mirror basis are
   `+0.0001280086893 i h_raw` and `-0.0003536675766 i f_raw`.
   Six PS invariants remain necessary; their coefficients are determined
   by two UV spurions, not new free family matrices.
2. The earlier same-colour diquark probes do not see the upper triplets
   because their antisymmetric colour contractions vanish. Adding the
   actual `uc(colour1)dc(colour2)` and `uL(colour1)dL(colour2)` vertices
   gives 36 real source slots and nonzero upper induced **delta C**.
   This extends a deliberately restricted probe set; it is not a claim
   that the earlier P-SPEC1 channel ledger was a complete operator basis.
3. Differentiating exact Gaussian elimination transports all of
   `(delta K,delta J,delta C)`. Actual upper `304 -> 249 -> 129` transport,
   direct and reversed order, canonical normalization and finite
   differences pass **108/108**. The extra elimination of whole L/R
   parents is an off-shell algebra regression, not a physical upper-scale
   decoupling. No source-free mass-only matching is promoted.
4. The unnormalized K carrier must not also contain scalar external-leg
   normalization in J. Canonicalization transports the curvature and J
   together. Double-counted scalar legs and dropped induced delta C are
   detected. The contribution composer rejects mismatched backgrounds,
   scales, source bases and prescriptions, plus duplicate diagrams.

`verify_p54_finite_yukawa_interface.py` now has a shared fail-closed fit
guard requiring full upper/lower KJC, all-six PS evolution, finite
sequential seesaw, complete fermionic light feedback, and retained-mass
error control, in addition to the existing requirements. Its canonical
algebra still passes **17/17**. These metadata guards do not prove a
diagram inventory complete. The existing upper pure-Yukawa verifier was
also rerun in memory and remains **44/44**, without overwriting its
historical numerical report. The geometry reader accepts an optional
read-only cache directory; no Hessian or lattice scan was recomputed.
The compiled eight-page PDF was rendered and visually checked on every
page; no unresolved references or overfull-box warnings remained.

**Remaining limits / next order:**

1. Complete upper scalar-cubic kinetic, source-dependent potential,
   fixed-VEV counterterms and Wilson/operator matching in this same PS
   chart. The vector bubble keeps all eliminated scalar masses, but
   retained masses are expanded: `max |m_retained^2|/M_V^2=0.328605`.
   Establish a remainder bound or add the required higher orders before
   assigning precision-fit accuracy.
2. Complete the lower covariant **hard minus EFT** finite functional,
   including contact/box operators, and the finite sequential seesaw map.
   Do not add a broken-background one-step scalar metric to the upper PS
   threshold. A limited scalar-current response does not replace a full
   gauge/Lorentz operator basis or complete scattering amplitudes.
3. Only after that, solve the same-prescription full complex fermionic
   light eigenstate with one retuning and positive heavy block, then the
   constrained two-spurion global fit. The precise generalized eigenstate
   and first-order retuning formulas are recorded in the new TeX.

`full_finite_matching_complete=false`, `physical_fit_enabled=false`.
The KJC algebra and this vector subset are closed, but the requested
complete matching -> full light state -> physical fit chain is **not**
finished. No missing finite term is set to zero for that chain. No commit
or push was requested or performed in this checkpoint.

### P-SPEC1: common spectrum and interaction-defined channels (`tree subgate done`; 2026-09-10)

The adopted energy-spectrum interpretation preserves P54PQ-v2, its
stationary benchmark and its two correlated UV family matrices. Description
layers denote retained field coordinates, **not** extra spatial dimensions.
No mechanical rotation, helix, arbitrary observation projector, new scalar
parameter or lattice scan is introduced. The next physical matching gates
below remain open; this entry is not an alternative flavor fit.

Artifacts:

- `code/verify_p54_spectral_channels.py`;
- `output/p54_spectral_channels.{json,md}` and generated numerical TeX table;
- `tex/p54_spectral_channels_nested_elimination.tex` and its PDF.

The verifier passes **153/153**, reading the existing content-addressed
Hessians without recomputation. It uses the **exact tree** light doublet,
not the mixed-order bosonic-improved ray. The 33 real current slots consist
of the actual `HdaggerH` cubic source, real/imaginary h/f coefficients for
seven fixed fermion-component pairs, and four complete-eigenspace-averaged
vector-pair mass vertices. The Yukawa coefficients retain their unknown
family spurions explicitly. All SM intertwiner Ward identities and the
component hypercharge selection rules are checked; coloured examples use
declared matched components, not every possible colour channel.

After excluding 33 eaten directions, all **295 physical real scalar modes**
are kept, including the four tree Higgs zeros and the physical PQ zero.
The artifact exports 24 mass-eigenvalue clusters and full cross-channel
residues `R_lambda=J^T Pi_lambda J`, with positivity and two spectral-moment
checks. The declared current set has source rank 19 and generates a
**21-dimensional H-invariant observable subspace**, seeing 12 pole clusters.
The full and minimal response agree. These counts belong to this current
set with unit-spurion tensors, not to a new 21-field physical theory.

Two structural findings supersede the earlier pictorial projection idea:

1. Exact elimination must propagate the triple `(K,J,C)`, not K alone:
   `K'=Krr-Krh Khh^-1 Khr`, `J'=Jr-Krh Khh^-1 Jh`,
   `C'=C+Jh^T Khh^-1 Jh`. The source-source term is nonlocal before a
   derivative expansion and may carry poles. Actual `295=9+231+55`
   factorization, direct vs staged vs opposite order, and passive metric/
   coordinate covariance agree at five spacelike/complex momenta. Omitting
   C loses **100% of the HH-to-HH scalar-exchange response** in the final
   chart. That chart is a factorization tool, not a new local EFT census.
2. The initial expectation of a nonzero NN coupling to the physical PQ
   mode was **rejected**. After gauge removal the PQ tangent lies wholly
   in S, whose direct Yukawa tensor vanishes. All declared linear probe
   residues on it vanish. This does not remove the physical axion or claim
   it is globally decoupled/massless after QCD. The massive spectral moment
   nevertheless reproduces `omega CHN=-0.0848692011507 i f_raw` exactly in
   the existing common-phase convention.

The singlet channels see common radial m2/omega2 values approximately
`0.01352047, 0.09971569, 2.82916586` with different residues. These are
tree benchmark poles, not loop-complete or experimentally fitted masses.

Next: carry this **source-dependent** functional, including its nonlocal
contact and induced vertices, into the common finite hard-minus-EFT
matching. Do not replace the missing gauge/Goldstone/ghost/scalar diagrams
or the fitted fermion-inclusive light state by pole-visibility arguments.
Do not repeat this unchanged tree audit as a new physical milestone.

**Extra-dimensional comparison-only contract:** a single fixed geometry,
boundary condition and action must predict an entire tower and its probe
overlaps. For the unadopted flat Neumann interval, two shared parameters
give `m_n^2=M5^2+(n*pi/L)^2` and the n=2/n=1 squared-gap ratio 4. A common
probe profile fixes all residues. No separate geometry, profile or freely
chosen boundary condition per particle is allowed. No extra-dimensional
model has been fitted or merged into P54, and it does not block P.

### P0 Claim and action freeze (`in-progress`; highest priority)

The baseline is the predictive two-Yukawa-matrix non-supersymmetric branch

```text
P54PQ-v2: Spin(10) x U(1)_PQ with
3 x 16_F + 54_H,R + 126_H,C + 10_H,C + 1_H,C.
```

The complex gauge singlet is part of the action, not optional decoration.  It
is required by the PQ branch.  The shorter no-PQ field list is now named

```text
P54N: Spin(10) with 3 x 16_F + 54_H,R + 126_H,C + 10_H,C,
      no PQ and three complex symmetric Yukawa matrices.
```

P54N is comparison-only.  It may not be described as the same two-matrix
theory with the singlet omitted.

The continuity branch

```text
P-210: 210_H + 10_H + 120_H + overline{126}_H
```

is comparison/rescue only until P-54 fails a predeclared vacuum, flavor, or
proton-decay gate.  The minimal `45_H + 126_H + complex 10_H` model remains a
stress test.  P0 must write the full normalized action, scalar potential,
Yukawa terms, symmetry-breaking chain, light-doublet assumptions, cutoff,
and physical parameter count.  Route-B messenger data are optional; an
invalid or absent messenger may not delay the baseline action.

Acceptance:

- one primary action and one convention card, with every other branch marked
  comparison-only;
- an observable map `F: Theta/G_basis -> O` and an initial Jacobian-rank and
  parameter-count ledger;
- compact physical family symmetry stated as `SU(2)_F` when used;
  `SL(2,C)` is its holomorphic complexification, not a positive-energy global
  symmetry;
- `M_V`, `zeta K_tr`, and every family-breaking spurion counted separately;
  an arbitrary symmetric `M_V` may not be advertised as a predictive
  `K_tr` texture.

#### P0-A action/parameter/convention card (`done`; 2026-08-30)

The frozen renormalizable candidate is `P54PQ-v2`.  Its card fixes Minkowski,
generator, self-duality, tensor, PQ-charge, vacuum, Clebsch, light-doublet,
cutoff, and basis-quotient conventions; writes gauge, kinetic, Yukawa, and all
declared scalar terms; and exposes the P1 spectrum interface.  Artifacts:

- `tex/p54_action_parameter_convention_card.tex` and
  `output/pdf/p54_action_parameter_convention_card.pdf`;
- `code/verify_p54_action_card.py`;
- `output/p54_action_parameter_convention_card.{json,md}`.

The verifier passes `20/20`.  It audits 29 scalar coefficient-monomial pairs,
operator dimensions, Hermiticity, continuous PQ charge, field content,
Pati--Salam dimensions, and exact rational fingerprints of the 54 vacuum.

Two definition-level corrections are authoritative.  With

```text
q_PQ(Psi,Phi,Sigma,phi,S)=(-1,0,+2,-2,-4),
```

the frequently printed `eta1 Sigma Sigma* Sigma* phi` term has charge `-4`.
The frozen term is the uniquely charge-neutral
`eta1 Sigma Sigma* Sigma phi + h.c.`.  Likewise the printed chi4 contraction
must use `Phi_mn`, the two free five-form indices; `Phi_ij` repeats indices
three times and leaves `m,n` free.

The independent basis audit found that version 1 omitted the allowed
renormalizable invariant `chi7 Phi_ij phi_i phi_j S* + h.c.`.  It is unique
because `54` occurs once in `Sym^2(10)`, and it is PQ neutral.  Version 1 is
therefore superseded, not retained as an alternative convention.

The corrected five complex-coupling phase vectors have exact rank two.
Consequently three invariant scalar CP phases survive:

```text
delta1 = arg(eta3) - 2 arg(eta1),
delta2 = arg(chi4) + arg(chi6) - 2 arg(eta1).
delta3 = arg(chi7) - arg(chi6).
```

The count is `60` raw real action coefficients, `49` after the classical
field-basis quotient, and `48` continuous zero-temperature observable
parameters after the anomalous PQ reparametrization.  This closes the P0-A
definition subgate.  The audit and P1 results below supersede the former open
items; the remaining P0 debt is the higher-dimension PQ-quality policy.

**Canonical complex Yukawa dictionary (2026-09-05, latest authority):**
the `20/20` absolute Clifford audit is now completed by the `28/28`
common-phase audit `code/verify_p54_common_yukawa_phase.py`. It conjugates
the actual scalar embedding and spinor chirality together and verifies
standard-matter hypercharges and all 45 covariance equations. With raw
couplings in that dictionary, `h_D=sqrt(2) h_raw`,
`f_D=2 f_raw/sqrt(3)`, `kappa_R=+i 4sqrt(2)` and
`kappa_LL=-i 4sqrt(2)`. Thus `f_M=i 2sqrt(6) f_D` at the Spin(10)
boundary; the old positive `2sqrt(6)` was an absolute-CG convention.
The `i` is transported component phase, not an extra CP parameter.
The action card records both conventions without adding a field/matrix.
In the actual retained copy order, up/nu use `(c1,c4)` and down/e use
`(conj(c2),conj(c3))`, with `-3` on lepton 126 vertices. Random full
copy/family transformations and independent neutrino rephasings preserve
the complete type-I+II contraction. The UV ratio remains protected on
the complete-parent one-loop parity submanifold proved below; do not
impose it on finite-matched off-locus data or after arbitrary decoupling.
Do not reuse old magnitude-only copy maps or historical external-matrix
CW numbers as physical predictions.

#### P0-B invariant basis and PQ anomaly/global form (`done`; 2026-08-30)

Artifacts:

- `tex/p54_spin10_pq_audit.tex` and
  `output/pdf/p54_spin10_pq_audit.pdf`;
- `code/verify_p54_spin10_pq_audit.py`;
- `output/p54_spin10_pq_audit.{json,md}`.

The independent verifier passes `17/17`.  It enumerates all PQ-neutral
multidegrees through degree four, proves the singlet multiplicities from
`Spin(10)` symmetric-product channels, and numerically checks that
`126 x 126bar` has no symmetric-traceless `54` bilinear.  The exact anomaly
ledger is

```text
A[Spin(10)^2-PQ] = -6,
A[SU(3)c^2-PQ] = -6,
Nhat_QCD = -12,
A[grav^2-PQ] = A[PQ^3] = -48.
```

The faithful global group is
`(Spin(10) x U(1)_PQ)/Z4_diag`: a PQ rotation by `pi/2` is undone by the
`Spin(10)` center on every field.  Hence the naive scalar-gcd value 6 is
reduced to the physical `N_DW=3`.  This closes the normalization audit but
opens an honest cosmology choice: pre-inflationary PQ breaking, extra
anomalous matter, or a bounded explicit bias.  `N_DW` may not be reset to one
by convention.

### P1 Vacuum, spectrum, and exact group data (`done` at tree-level algebraic gate; 2026-08-30)

Solve all stationarity equations of P0, verify Goldstone/gauge-orbit
alignment, and export every heavy and light SM irrep with its mass and
uncertainty.  Build the exact `Spin(10)` generator/Clebsch package for the
selected scalar content.  A local minimum is not a metastability theorem
without competing vacua and a bounce estimate.

Completed artifacts:

- `tex/p54_p1_stationary_hessian_spectrum.tex` and
  `output/pdf/p54_p1_stationary_hessian_spectrum.pdf`;
- `code/verify_p54_p1_hessian_spectrum.py` and
  `requirements-p54-p1.txt`;
- `output/p54_p1_stationary_hessian_spectrum.{json,md}`.

The JAX verifier passes `14/14`.  It differentiates the complete corrected
tensor action on all
`328 = 54 + 2*126 + 2*10 + 2*1` canonical real scalar coordinates.  The full
gradient residual is `3.80e-16`; the `328 x 328` Hessian has 38 zeros and no
negative modes.  Projection and Casimir tests identify them as exactly
`33 gauge + 1 PQ + 4 real light-doublet` modes.  The light zero subspace has

```text
(C3,C2,Y^2) = (0,3/4,1/4),
field weight = 0.9999920937 in 10_H + 7.9063e-6 in 126_H.
```

The 35-row SM-irrep ledger covers all 328 real coordinates with maximum
unbroken-algebra leakage `1.93e-10`.  The vector spectrum in `g10^2` units is

```text
0 x12, 0.1225 x8, 0.6125 x1, 1 x12, 1.1225 x12.
```

The radial solver finds 27 signed roots; the declared `omega=+1` orientation
is the lowest enumerated branch.  This is a dimensionless tree-level
existence benchmark, not a measured-scale phenomenological fit or a
nonperturbative lifetime theorem.  P2 must run/match this actual spectrum,
propagate covariance, and reimpose the doublet condition after loop
corrections.

### P2 Running, thresholds, and bosonic CW matching (`physical matching reopened 2026-09-05`)

**The following historical P2 numbers are reproducible algebraic outputs,
not a certified Wilsonian matching solution.** The new independent
`code/verify_p54_ps_active_census.py` / `output/p54_ps_active_census.json`
passes `14/14` audit checks and finds two failed historical physical gates:

- The four active parents `(1,2,2)10`, `(15,2,2)126`, `DeltaL`, `DeltaR`
  give `a_PS=(2/3,26/3,26/3)` and `b44=3551/6`. The stored table needs
  one extra complex `(6,1,1)` to obtain `a4=1,b44=1209/2`.
- The lower threshold includes only the SU(2)L-singlet DeltaR parent,
  so `lambda_I,L'=0`, but its own beta discontinuity requires `-71`.
  The matching residual derivative is exactly `71/(12 pi)`, not a small
  numerical tolerance. Old scale/covariance and PQ-extension scale replays
  therefore cannot be promoted to physical results.

A same-field-content four-parent repair closes all six logarithmic
identities: physical scalar indices are `(12,6,6)` at U and `(67,71,67)`
at I; `lambda_U'=(72,120,120)`, `lambda_I'=(-46,-71,-41/5)`.
Finite staged thresholds and a new scale solution remain open. At the
simultaneously broken background the parent/mass/Goldstone projectors need
not commute, so replacing a projector in the old formula is insufficient.
See `output/p54_ps_yukawa_matching_research.md` and the new full matching
TeX below. This correction supersedes every closure statement in the
historical P2 ledger that follows; no old outputs are erased.

Artifacts:

- `tex/p54_p2_running_thresholds_cosmology.tex` and
  `output/pdf/p54_p2_running_thresholds_cosmology.pdf`;
- `code/verify_p54_p2_running_thresholds.py`;
- `output/p54_p2_running_thresholds_cosmology.{json,md}`;
- `code/search_p54_hierarchical_p1.py` and
  `output/p54_hierarchical_p1_search.{json,md}`;
- `code/verify_p54_p2_two_site_matching.py` and
  `output/p54_p2_two_site_matching.{json,md}`;
- `tex/p54_p2_bosonic_cw.tex`,
  `output/pdf/p54_p2_bosonic_cw.pdf`,
  `code/verify_p54_p2_bosonic_cw.py`, and
  `output/p54_p2_bosonic_cw.{json,md}`.

The integrated verifier passes `21/21`; the parent-resolved fixed-point
verifier passes `11/11`.  The historical no-threshold regression is

```text
MI = 4.60943e13 GeV,
MU = 1.18126e15 GeV,
alphaU^-1 = 37.60198,
MI/MU = 0.0390213.
```

The original `sigma/omega=0.35` point remains a deliberately failed
regression: its no-threshold mismatch factor is `8.969`.  It is no longer the
active blocker.  Keeping every P1 cubic/quartic coupling fixed and deriving
only the three radial quadratic masses at each hierarchy gives an isolated
stationary branch with 38 zero modes, no tachyon, and exactly one light
doublet.  At the threshold fixed point,

```text
sigma/omega input = 0.1265,
MI/MU output       = 0.126524923,
MI                 = 1.07986e14 GeV,
MU                 = 8.53475e14 GeV,
alphaU^-1          = 39.70141.
```

The fixed-point mismatch is only `1.97e-4` fractionally.  Exact joint
Pati--Salam Casimir projectors decompose all 328 real coordinates and place
the 126 VEV uniquely in `Sigma126:(10-pair,1,3)`.  Thresholds use the spectral
operators `Tr(P_parent t_i^2 log M)`, so mixed broken-phase states are never
assigned by their nearest mass.  The adjoint ledger independently finds nine
intermediate and 24 GUT massive vectors.

The full 35-sector ledger and universal `(84,84,84)` index are retained.
Two-site covariance now propagates the experimental inputs and independent
10% log-mass errors on 27 exactly degenerate scalar/vector blocks.  In
`(log10 MI,log10 MU,alphaU^-1)`,

```text
sigma_exp   = (0.00794,0.02559,0.06130),
sigma_thr   = (0.03952,0.02848,0.24057),
sigma_total = (0.04031,0.03829,0.24826).
```

The loop-level one-doublet theorem is nevertheless fixed.  With

```text
kappa = d mh^2/d xi02 = -0.999997474290,
gap to the next doublet = 0.0305570303 omega^2,
delta xi02 = 1.00000252572 h^T Pi_D h/omega^2,
```

the projected self-energy norm must remain below the gap.

The requested bosonic calculation is now complete in a declared
background-field Landau-`MSbar` hard-matching scheme at `muU=gU*omega`.
The Fréchet second variation of the field-dependent `328 x 328` scalar
operator and exact vector orbit mass operator retains 290 hard scalars and
33 hard vectors, with the 38 scalar and 12 vector EFT modes removed.  It
passes `14/14` checks and gives

```text
Pi_scalar^hard/omega^2 = -0.0182525705951 I4,
Pi_vector^hard/omega^2 = -0.00198408049986 I4,
Pi_B^hard/omega^2      = -0.0202366510950 I4,
delta xi02,B(muU)      = -0.0202367022071.
```

This is a scheme-dependent zero-momentum matching curvature, not a pole
mass.  Its scale replay at `mu/muU=(0.5,1,2)` gives
`kappa_B/omega^2=(+0.0168186,-0.0202367,-0.0572919)`, which must be cancelled
by parameter running and the fermionic matching term rather than interpreted
as an error bar.  The heavy-Yukawa projection is therefore retained without
being set to zero:

```text
eta_Y(mu) = h^T Pi_heavy-Y(0;mu) h / omega^2,
delta xi02(mu) = [kappa_B(mu)/omega^2 + eta_Y(mu)] / w10.
```

The nearest-heavy-doublet norm inequality remains a P3 profiling check.

In parallel, the minimal domain-wall repair candidate `P54PQ-F10` adds two
left-handed `10_F` fields with `qPQ=+2` and mass operator `S 10_F 10_F`.
It shifts the mixed anomaly from `-6` to `-2`, so `Nhat=-4` and the same
diagonal `Z4` quotient gives physical `N_DW=1`.  Its one-loop beta shifts are
universal.  Its explicit two-loop benchmarks at `y_F=(0.5,1,2)` place
`M_F/MU=(0.1563,0.3132,0.6276)` inside the PS interval and retain
`MI/MU=(0.12548,0.12587,0.12626)`.  It remains a separate branch because
`y_F` is a new threshold parameter.

Running, two-site thresholds, hierarchy closure, covariance, the projected
bosonic CW coefficient, and first-order light-mass matching are computed.
The complete complex loop eigenpair and tadpole scheme remain open. A complete
gauge-independent pole mass is deliberately not claimed.  P3 may now open,
but it must fit/profile `eta_Y(mu)` and test the heavy-doublet norm bound;
neither quantity may be silently set to zero.

### P3 Global flavor, seesaw, and identifiability (`in-progress; corrected audit 2026-09-05`)

Authoritative artifacts:

- `tex/p54_theory_audit_and_repair.tex` and its PDF: complete derivations,
  corrections, EC assessment, and two constructive repair routes;
- `code/verify_p54_theory_audit.py` and
  `output/p54_theory_audit.{json,md}` (`14/14`);
- `code/verify_p54_relaxed_light_quartic.py` and
  `output/p54_relaxed_light_quartic.json` (`7/7`);
- corrected `code/extract_p54_p3_light_overlaps.py` and its ledgers
  (`5/5`);
- `tex/p54_p3_flavor_fermionic_cw_gate.tex` is now the corrected short
  P3 report. Its compatibility verifier/output delegates to the new audit.

**Latest phase correction, not a new fit:** the earlier magnitude-only
P3 report and its subsequent convention stress test are historical. The
explicit standard-matter Clifford dictionary now fixes the actual map:

```text
c = (phi_hol, phi_anti, Sigma_hol, Sigma_anti), with complex signs retained.
Yd=a h+d f; Ye=a h-3d f; Yu=b h+e f.
(a,b,d,e)=(conj(c2),c1,conj(c3),c4).
tree: r=-0.905014358012, s=-29.7554251165.
bosonic: (a,b,d,e)=(-.7631652014906,.6462011667549,
                    1.2183583687513e-5,-.001710897393866)
          up to the complete tiny imaginary parts recorded in JSON.
bosonic: r=-.846738249455, s=+165.843989136.
```

The exact two-matrix identity remains
`Yu=r[(3+s)Yd+(1-s)Ye]/4`, at one common PS matching scale before
finite corrections. But the earlier `yt<=0.00801414124`/`0.111193537`
bounds used obsolete maps or magnitude envelopes and are not current
constraints. The small coefficient `d` crosses zero under the bosonic
correction. Use `(a,b,d,e)` directly in the fit, not the near-singular
ratio `s=ae/(bd)`. Large `s` is not itself a large physical rotation.

The earlier factors `69.49`, repair target `|r|>=62.9`, and empty
`eta_Y` profile-domain claim are **withdrawn**. A large Gaussian
chi-square lower bound does not empty a parameter domain. The frozen
tree geometry has diagnostic tension, but complete physical PS
running/thresholds and the fitted fermion-inclusive scalar eigenvector
are still needed to assess it.
A user-chosen percentage envelope does not bound those effects.
No global flavor/seesaw fit, neutrino prediction interval, or full-model
exclusion is claimed.

#### P3-A: scalar/loop interface repair (`bosonic subgate done; full physical gate open`)

P2 computed only `P Pi_B P`. It did not compute `Q Pi_B P`,
`Q Pi_B Q`, or the tadpole/vacuum shift. For the full complex copy matrix
`M=[[A,b^dagger],[b,C]]`, with `C>0`, the exact zero condition is
`A=b^dagger C^-1 b`, with light vector proportional to
`c0-Q C^-1 b`. The missing mass is formally second order in loops, but
the missing eigenvector rotation is first order and matters for flavor.
A projected mass retune plus a Q-block norm bound is insufficient.
The new `tex/p54_fixed_vev_full_doublet_matching.tex` and PDF give the
complete derivation. `code/verify_p54_full_doublet_cw.py` and its JSON/MD
pass `18/18`. At the historical background, the SM-singlet tangent is
exactly three radial plus two gauge/PQ phase directions. The invariant
fixed-VEV prescription is

```text
delta_mu2=5 t_omega/(12 omega), delta_nu2=t_sigma/sigma,
delta_mus2=t_s/vs;
H_ct=-delta_mu2 P54-(delta_nu2/2)P126-delta_mus2 PS.
```

The factor `1/2` on the canonical 126 counterterm is essential. Two phase
Ward checks close below `1.3e-13 omega^2`. All 4 complex neutral doublet
copies are retained, including independently checked imaginary jets.
The retuned **bosonic one-loop-truncated** matrix gives

```text
delta_xi02=-0.0202627041913 omega^2,
eigenvalues/omega^2=(0,0.03025666325,0.19979763559,0.43263922302),
|c|=(0.6462011668,0.7631652015,0.00001218358,0.0017108974),
angle(c_tree,c_B)=0.0329814289 rad,
Schur residual=6.2e-18 omega^2.
```

The original coarse Q/gap bound **fails**, with ratio `4.79545`, but this
does not prove a tachyon. The sharper exact congruence test gives
`spec(C0^-1/2 DeltaC C0^-1/2)=(-0.0109064,0.4393631,0.5234555)`, hence
`C >= 0.9890936 C0 > 0` at this truncation. This is a positivity
certificate, not proof of small higher loops; the 52% relative correction
and large raw tadpole subtraction require a control audit. Other-irrep
loop stability and physical pole masses are not computed. Most of the
light rotation is in the 10 pair; `w126=2.92732e-6` remains small.

The missing heavy-relaxed *tree* quartic has now been computed from the
same 328-real-coordinate action:

```text
J_i=V3[e_i,q,q]; C=positive physical scalar Hessian.
lambda_EFT=V4[q,q,q,q]/6-J^T C^-1 J/2
          =0.429206300-0.0200527602=0.409153539>0.
```

All 290 massive scalars are included; two light directions, two derivative
steps, zero-mode source orthogonality and direct relaxed paths pass
`7/7`. This supports local tree stability modulo symmetries, not a
global boundedness theorem, loop vacuum, or physical Higgs mass.

#### P3-B: complex fermionic projection (`conditional diagnostic done; fit open`)

For fixed symmetric invertible `MR`, `D(z)=sum z_a Y_a`, the exact
complex hard-Majorana Hessian is
`Pi_ab=-Tr(Y_a^dagger Y_b W)/(8 pi^2)`,
where `W=MR^dagger MR [log(MR^dagger MR/mu^2)-1]`.
Its real form retains the imaginary off-diagonal blocks.
An independent four-real-coordinate finite difference test gives relative
error `3.7e-8`.

The published Mummidi--Patel matrix pair, transported with positive-real
corrected overlaps, gives conditional `eta_Y=1.28882418e-8` and
`||Q(Pi_F-delta_xi_F P10)Q||=9.68094794e-4 omega^2`, or
`0.03168157` of the tree gap. This is a complex-copy heavy-neutrino
diagnostic, not fitted matrices, an actual scalar-phase match, or the
total boson-plus-fermion Q certificate. New heavy fermions would add terms.

**Latest completion/correction:** use `f_M`, not the unit-Dirac `f_D`, in
`MR=sigma f_M`; the actual common-phase UV conversion is
`f_M=i 2sqrt(6) f_D` (the former positive ratio was a magnitude convention). The
fixed-VEV fermionic term is

```text
t_sigma,F=-Tr[X^2(log(X/mu^2)-1)]/(8 pi^2 sigma),
delta_nu2,F=t_sigma,F/sigma,
DeltaD_F=Pi_F-(delta_nu2,F/2)P126, X=MR^dagger MR.
```

The scale is held fixed under field differentiation.
`code/verify_p54_fermion_tadpole.py` and JSON/MD pass `11/11` on three
synthetic complex families, four copies/eight real directions, degenerate
Takagi masses and family/copy basis changes. A synthetic assembly checks
the Schur solve but is explicitly not fitted heavy Yukawa data.

`code/verify_p54_typeii_triplet_source.py` and JSON pass `10/10`.
The complete two-copy 54/126 neutral-triplet inverse gives
`omega*z126,canonical/v^2=-0.004700734233` in its exported phase convention.
The 54 source through mixing partly cancels the direct 126 source;
neglecting mixing would miss a 5.84% reduction relative to the unmixed
response. This uses tree K/V3 and the bosonic-improved light direction,
not a loop-accurate Wilson coefficient. Compute the canonical LL
intertwiner is now contracted in the common phase dictionary below; use
`C5II=-Y_LL K^-1 J`, together with
`C5I=-Ynu MR^-1 Ynu^T`; type II is not authorized to be set to zero.

#### P2/P3 covariant EFT, finite operators and moving thresholds (`latest authority`; 2026-09-06)

The preceding checkpoint was committed and pushed first as `4ad74e7` to
`origin/main` (`boypatrick/GUT_routeA`). The following work is subsequent
uncommitted research. Detailed derivations are in
`tex/p54_covariant_eft_finite_yukawa_seesaw.tex` and its PDF. This entry
supersedes the older four-holomorphic-matrix/zero-CHN matching assumptions,
without changing the frozen P54PQ-v2 action or adding UV family parameters.

**Structural results, not decimal fitting:**

1. **Local lower covariant tree action constructed.** The 55-heavy scalar
   valley and algebraic elimination of the 24 upper vectors give
   `g=T^T(1-O(O^TO)^-1O^T)T` and actual derivative vertices. The 24
   directions form a coset, not a subgroup. The full nine lower gauge
   masses/Goldstone norms and previous Schur metric are reproduced.
   The tree upper-relaxed quartic on the improved ray is `.407929900088`,
   with an accompanying `.024787549073 h²(∂h)²/omega²` operator; it is not
   an already matched SM quartic. The exact quadratic resolvent is
   exposed in the API. Only 25 retained eigenmodes couple linearly to the
   integrated block; 12 have a pole-expansion ratio `1.0276132` but small
   mixing. Handle those with the exact resolvent, not a new global
   scale-separation blocker or another parameter scan.
2. **Actual finite upper pure-Yukawa diagrams computed.** All 55 scalar
   eigenstates are included, with 24 nonzero six-parent couplings.
   Real/imaginary mass splitting generates two conjugate-bidoublet
   operators `Htilde,Ftilde`. Their coefficients are cubic polynomials
   in the original `h_raw,f_raw`, times exact split logarithms, not free
   nuisance matrices. The four-parent finite boundary is **not closed**.
   The six-invariant real-Weyl beta tensor **is closed**; the original
   four-parent routine remains a valid holomorphic zeroth-order subflow,
   not the full resummed finite-matched flow. The generated mirrors carry
   the required PQ/axion phase and vanish in the degenerate-mass limit.
3. **Actual scalar kinetic matching subset computed.** A same-action
   one-step scalar-cubic bubble gives the full four-complex-doublet
   kinetic matrix, including hard-soft pairs and charged partners.
   On `c_B`, `delta Z=.00243676715822`. All electroweak Ward tests pass.
   This is the background-field Landau scalar contribution, not a
   gauge-independent pole result or the complete two-site subtraction.
   Its 290-hard result must not be added to an upper subset a second time.
4. **The actual tree lower boundary has nonzero CHN.** Eliminating the
   three massive neutral radial scalars gives
   `omega C_HN/f_raw = -.0848692011507 i` on the exact tree light ray,
   for `L += C_HN NN H†H+h.c.` without a factor 1/2. Gauge/PQ zero modes
   are retained, never inverted. This is independently checked by
   solving the same radial potential at finite Higgs field. Thus
   `C_HN=0` is false for the frozen P54 tree matching, even though it
   defines a one-loop invariant sterile-EFT subspace. Improved-light
   frozen-tree projections are recorded separately as mixed-order data.
5. **Moving sterile thresholds are implemented.** Matrix running solves
   `mu=M_i(mu)` with level crossings and degenerate Takagi clusters,
   uses exact block-tree Schur matching, and continues below the final
   sterile state. The 2024 Weinberg feedback correction is included;
   old Antusch equations are negative controls only. Actual nonzero
   `C_HN` must accompany the lower boundary, including its feedback into
   `Ynu`, `MR`, the Higgs quadratic term and quartic. The current engine
   is a dipole-free dimension-five subsystem with the PQ axion held as
   a spectator; finite dipole matching and quantum axion effects are
   not proven absent. Neither these
   synthetic trajectories nor terminal mass proxies constitute a fit.
   Starting with diagnostic `qH=0` generates terminal values about
   `2.22e-6, 2.98e-6` in omega² units, above the terminal scale squared.
   Thus a light-Higgs endpoint is not consistent without the already
   required quadratic finite matching/retuning. The mass-independent
   MS-bar continuation is diagnostic, not a physical low-energy EFT.

Artifacts are `verify_p54_lower_covariant_eft.py`,
`verify_p54_upper_yukawa_thresholds.py`,
`verify_p54_scalar_kinetic_matching.py`, `verify_p54_scalar_chn.py`,
`verify_p54_sequential_seesaw.py`, with same-stem JSON/MD records.
The generated records retain exact matrices, source hashes, scope and
tests, so subsequent work does not reconstruct formulas from rounded text.
Final regression counts are respectively `24/24`, `44/44`, `24/24`,
`18/18`, `46/46` (156 passing checks). These count bounded mathematical
and numerical regressions, not completed physical matching conditions.

**Still open, and now precisely narrowed:** upper vector fermion
kinetic/vertex diagrams, vector/scalar/Goldstone/ghost scalar kinetic
matching, upper scalar/tadpole Wilson data, lower covariant hard-minus-EFT
matching, and finite sequential type-I/type-II matching in the enlarged
operator basis, including the appropriate PQ/axion Wilson couplings
before full physical promotion. A completed tree action or one-loop RGE is not a finite
matching coefficient. After these are assembled, rerun physical scales
and the constrained scalar/flavor fixed point. No physical fit, full-model
exclusion, Higgs pole, determinant or portal is promoted by this checkpoint.

#### P2/P3 complex phase, full PS flow and local-feedback checkpoint (`historical; read with the 2026-09-06 operator-basis correction`)

The prior fixed-VEV checkpoint was committed/pushed as `b8e7398` on
`origin/main` before these calculations. Current detailed derivations:
`tex/p54_complex_phase_ps_matching.tex` and its PDF. Reproduction/source
hashes are in the five same-named Python/JSON/MD artifact sets below.

| Artifact prefix | Tests | Precisely completed scope |
|---|---:|---|
| `p54_common_yukawa_phase` | 28/28 | Actual four-copy/standard-matter complex dictionary, including common type-I+II phases |
| `p54_ps_finite_thresholds` | 25/25 | Finite upper gauge threshold, actual lower tree heavy valley and one-step regression; not lower finite matching |
| `p54_ps_yukawa_flow` | 29/29 | All-active four-parent one-loop matrix flow and two-loop gauge Yukawa trace, actual tensors and nonzero ODE |
| `p54_finite_yukawa_interface` | 17/17 | Canonical finite matching algebra/kernels with synthetic finite arrays, not P54 diagram values |
| `p54_self_consistent_light` | 36/36 | Actual Clebsch-constrained local Schur feedback on two synthetic family choices; no data fit |

1. **Phase dictionary closed.** The common LL/singlet coefficients are
   `kappa_LL=-i4sqrt(2)`, `kappa_R=+i4sqrt(2)`. At the bosonic light
   point the mixed-order seesaw is
   `Mnu=(v^2/omega)[(.01628382104964 i) f_D +
   (.80681480328827 i) Ynu f_D^-1 Ynu^T]`.
   Both phases and the 54/126 scalar source are transported, not replaced
   by absolute values. This is not full loop Weinberg matching.
2. **Upper finite gauge matching closed locally.** Keep all scalar
   parameters fixed and solve the PS-symmetric heavy radial saddle:
   `(omega,vs)=(1.00081030592,.254277697084)`. The 55 integrated real
   modes have `m2_min=.102882040959`; the 124 tachyons lie exclusively in
   retained fields. At `mu=gU*omega_PS`,
   `lambda_U=(2.94405268936,7.60297883009,7.60297883009)`, including
   finite vector constants `(4,6,6)`. All six site-log identities still
   close. The retained/integrated masses have no uniform hierarchy;
   this does not certify the entire two-site renormalizable truncation.
3. **Lower tree geometry computed, finite functional open.** The same
   55-heavy block remains positive at the actual broken vacuum:
   `m2_min=.099450257411`. On `B=(I-PGupper)Bactive`, form
   `W=BH^T H B`, `E=-C^-1W`, `S=B^T H B-W^T C^-1W`,
   `G=B^T B+E^T E`. The metric spectrum is `[.984249788817,1.01800727050]`;
   `(S,G)` has 13 zeros and no negative modes. The API exports the full
   matrices. The finite lower coefficient must still include the
   momentum-dependent Schur operator and its induced covariant vertices;
   neither old weighted logs nor generalized `(S,G)` masses alone suffice.
   The one-step SM finite result is only an independent regression.
4. **Full all-active PS flow closed at the declared orders.** Canonical
   boundary units are `H=sqrt(2)h_raw`, `F=4f_raw`,
   `L=R=2sqrt(2)f_raw`; lower Dirac projection uses `F/(2sqrt(3))`,
   while `MR=+2i sigma_dimful R`, `ML=-2i Delta_L L`.
   Generic 248-real-scalar/48-Weyl beta tensors equal the closed four-matrix
   system to `2.28e-16`. SM-top, no-DeltaL, independent family covariance,
   and nonzero coupled ODE regressions pass. The two-loop gauge flow now
   includes its actual Yukawa trace, not only gauge/matter coefficients.
   **New exact insight:** with `gL=gR`, symmetric `H,F` and
   `F=sqrt(2)L=sqrt(2)R`, the full one-loop system preserves all these
   relations. Do not invent four independent UV family matrices to fit
   data. Off-locus finite matching and partial decoupling still require
   the general flow. Complete SM/sequential-seesaw evolution is not yet
   replaced by this all-active PS subroutine.
5. **Differentiable local scalar/flavor solve closed.** Two synthetic
   actual-Clebsch examples retune to positive gaps `.03026151664` and
   `.03026452622`. The analytic
   `xi'=c†DeltaF'c/(c†P10c)` and
   `c'=-D^+(DeltaF'-xi'P10)c` agree with re-solved finite differences.
   Rebuild `J_A(c)=T_Apq cR_p cR_q` for every new light vector using eight
   existing cached Hessians; do not freeze the old type-II source.
   Tree triplet K/J remains a mixed-order approximation. The finite
   kinetic interface transforms both D and Yukawa vertices together:
   a zero-mode's original-coordinate direction changes only by
   `1/sqrt(c†ZHc)`, not an independent kinetic-induced rotation.

**Remaining physical gate:** derive the upper-matched covariant PS scalar
and gauge functional, finish lower finite gauge matching and both-site
finite Yukawa vertex/kinetic diagrams, and decide the actual mass-ordered
decoupling intervals from this spectrum. Then propagate the matched
matrices through SM/sterile-neutrino/Weinberg EFTs and iterate the local
scalar solve inside the constrained likelihood. Missing diagram values
must not default to zero. No new physical scales, best fit, neutrino
prediction, pole Higgs mass or full-model exclusion is asserted.

#### P3-C: constructive inverse problem (`derived; not solved`)

Prefer an action-derived inverse scalar/flavor solve to a blind optimizer.
For effective `H,F,r,s` and declared norm caps, normalized overlap
magnitudes exist iff

```text
(1+|r|^2)||H||^2/h_*^2+(1+|rs|^2)||F||^2/f_*^2 <= 1.
```

For fixed VEVs and candidate light vector `c`, the action Hessian is affine
in real invariant couplings. Solve `M_D(p)c=0`,
`Q_c M_D(p) Q_c >= Delta Q_c`, stationarity and other-irrep positivity
as a semidefinite inner problem. Searching over `c` and recomputing loops
remain nonlinear. Every candidate must rerun the relaxed quartic, complete
spectrum, perturbativity, competing-vacuum tests and P2 thresholds.

The historical numerical overlap bound at the old diagnostic `r` is
superseded by the actual complex copy dictionary above. The abstract
norm/feasibility identity can be reused only with the correct evolved
normalizations and declared threshold terms; no old decimal bound is a
constraint on the repaired fit.

#### P3-D: explicit alternatives and no-go filters (`proposed; not adopted/fitted`)

- **Antisymmetric-only repair:** at fixed geometry, even complex
  `120` spurions obey `yt-yc<=2U`, since a complex 3-by-3
  antisymmetric matrix has singular values `(g,g,0)`. The diagnostic
  target fails. This does not exclude a full 120 scalar extension
  that also changes the light state.
- **PQ mirror-10:** the allowed dimension-five operator
  `S* (16_i 16_j)_10 phi / Lambda` adds a symmetric `g` and yields
  `Yu-r[(3+s)Yd+(1-s)Ye]/4=(|a|^2-|b|^2)g/a`.
  A rank-one third-family example reproduces its three diagonal targets
  with `(h33,f33,g33)=(3.41540,-0.55474,-3.08124)`.
  It is not a global fit. A vectorlike
  `16(+3)+anti16(-3)` pair realizes rank at most two.
  Exact VEV mixing rescales both original Yukawas and gives
  `h_UV=y_UV=4.11034` in the toy; the vertex loop factor is about
  `0.107` before multiplicities, so weak perturbativity is not certified.
  New thresholds, running, scalar corrections and fit are required.
  Its PQ anomalies cancel and `N_DW=3` is unchanged.
- **Two-doublet EFT:** remains an openly different branch with new
  spectrum, running and flavor-changing-neutral-current obligations.
- **Einstein--Cartan/central twist:** minimal Planck-scale algebraic
  torsion produces a suppressed dimension-six current interaction, not
  an automatic independent symmetric Yukawa. Its NDA estimate is not a
  calculated universal correction. A low-scale propagating torsion sector
  requires its own healthy action; a contact EFT cannot be used above its
  mediator mass. Discrete central quotients and Nieh--Yan terms do not
  automatically repair the PQ anomaly/domain-wall number.

Recent connections checked: Haba--Shimizu--Yamada (2023), PRD 108 095005
(three symmetric Yukawas without visible PQ); Erdmenger et al. (2024),
arXiv:2409.06766 (torsional anomaly/renormalization); Gao (2025), PRD 111
055013 (a distinct spontaneous-CP/two-doublet model). Vectorlike rank-one
flavor has older precedent, arXiv:0911.2242; no global novelty claim is made.

#### P3-E: retained global-fit/reporting contract (`open`)

Fit quarks, charged leptons, CKM, neutrinos and PMNS simultaneously with
the matched action. Publish likelihood, priors, pulls, parameter count,
`chi2/dof`, and observable-Jacobian rank/identifiability. The baseline
Majorana tensor is action-related to `f`; do not silently make it generic.
Retain the original comparison program as explicitly different hypotheses:

```text
P3-0: generic allowed Majorana matrix (null/comparison, not the frozen baseline),
P3-1: M_R = M_V + zeta K_tr (declared extra matching hypothesis),
P3-2: UV-restricted M_V and zeta only if Layer U derives them.
```

No fitted datum counts as a prediction. Only after a viable fit publish
out-of-fit intervals for `m_bb`, neutrino mass sum, CP phases, heavy-neutrino
hierarchy and flavor ratios relevant to proton decay.

### P4 Proton decay and amplitude consistency (`blocked by viable P3 matching/fit`)

Use physical mass eigenstates and P3 flavor rotations to calculate the
complete relevant dimension-six channels and any branch-valid scalar
channels.  Route C is reclassified as a post-action amplitude-consistency
audit: residues must be extracted from complete amplitudes, with Ward/
Goldstone identities, crossing, partial waves, and applicable positivity
bounds checked after P0--P1.

### P5 Model comparison and conditional closure (`blocked by P3/P4`)

Compare P-54, P-210, the generic-Majorana null, and the `K_tr`-restricted
variant with parameter penalties and uncertainties, not best-fit residuals
alone.  Layer P closes when one action has a stable or declared metastable
vacuum, complete spectrum, two-loop matched running, identifiable flavor fit,
proton amplitudes, and at least two predeclared out-of-fit predictions.

Layer-P closure establishes a **conditional four-dimensional theory**.  It
does not claim a first-principles origin of the family curve, `zeta`, or a
global compactification.

## Layer U -- family UV

### U0 Object separation and claim boundary (`in-progress`)

Use distinct symbols for the family curve `Sigma_F`, order-parameter target
`P(Z)`, soliton moduli `M_B`, Berry base, and string matter curve `C_m`.  A
map between them must specify bundle, gauge representation, chirality,
locality, scale, anomaly, and kinetic data.  A degree-one map or equal Chern
number alone is not a portal.

The Route-E theorem is retained in its exact conditional form:

```text
H3 gives N_fam <= 3;
H3+ selects the semisimple/Killing branch and hence N_fam = 3.
```

H3+ is a selection axiom until U derives it or replaces it with a more direct
flux/index mechanism.

### U1 Direct Spin^c family construction (`open`; primary U target)

Prioritize a six-dimensional or defect construction in which a field in the
`Spin(10)` half-spinor `16` produces three four-dimensional chiral zero
modes.  On `CP1`, the microscopic bookkeeping target is

```text
K_CP1 = O(-2),
K_CP1^(1/2) = O(-1),
L_F = O(3),
K_CP1^(1/2) tensor L_F = O(2),
h0(O(2)) = 3, h1(O(2)) = 0.
```

Thus `O(2)` is the effective positive-spinor bundle after including the spin
connection; it need not mean two colours, two orbitals, or two unit cells.
The target theorem is representation-valued,

```text
Ind D_UV = 3 x 16,
```

not merely `c1=2` for a quantum-mechanical determinant line.

### U2 Anomaly, flux selection, and exact spectrum (`open`)

Acceptance requires flux quantization or boundary data, dimension-appropriate
bulk/defect anomaly cancellation or inflow, exactly three complete chiral
`16`s, no other massless chiral SM-charged states, and a nonzero fourth-mode/
KK gap throughout a predeclared deformation neighborhood.  Vectorlike pairs
must be listed and lifted explicitly.  The reason for the `O(3)` gauge line/
flux sector must be stated as topology, tadpole, boundary data, or dynamics
rather than hidden in the final `O(2)` notation.

### U3 Yukawa overlaps and four-dimensional matching (`open`)

Derive normalized zero-mode wavefunctions, overlap-generated Yukawa
structures, family-symmetry breaking, and matching onto the P0 action.  U3
passes only if the construction restricts parameters or produces an
out-of-fit relation; reproducing a fitted matrix by adjustable overlaps is
inverse reconstruction.

### U4 Optional global string realization (`open`, permanently optional)

Route D remains an optional implementation after U1--U3: resolved geometry,
quantized flux, tadpole cancellation, massless hypercharge, exotic-free
spectrum, instanton zero modes, and unwanted-operator veto are required for
promotion.  Failure of U4 does not fail Layers P or a non-string U branch.

Layer-U closure upgrades the geometric family ansatz to a first-principles
family origin.  It is not required to publish or test the conditional P
theory.

## Layer S -- independent soliton mathematics

### S0 Scope and non-promotion rule (`done`)

AP-E1--AP-E18 and the AP-E11 action are preserved as a standalone research
program.  Its valid questions include Gamma relaxation, graph-norm density,
defect measures, GTA, regularity, isolated `B=1` solutions, Hessians,
Callias/APS determinant lines, and WZW descent.  Results may be published as
mathematics or as properties of a declared hidden-sector EFT.

They do not prove four-dimensional family replication unless a future common
microscopic action independently supplies the full U0 bundle/gauge/chirality
map.  In particular,

```text
mu = 0 or GTA  does not imply  Ind D_UV = 3 x 16,
c1(det Callias) = 2 does not imply three SM families,
classical B=1 isolation does not select the Spin(10) action.
```

### S1 Current mathematical frontier (`open`, optional priority)

AP-E18 proves that GTA is not a universal consequence of the available
additive tangent estimates.  A scale-polished recovery sequence or direct
energy-amplitude tightness theorem remains a legitimate S-layer goal.  It is
not an urgent Route-F blocker.  No P or U task waits on `mu=0`, classicality,
the same-action Hessian, determinant variation, or a soliton portal.

### S2 Portal retirement from the main dependency graph (`done`)

The former degree-one Route-E portal is retired as a required milestone.  A
future portal proposal must begin as a new optional model with one common
microscopic action and must pass U0--U3 independently.  Topological degree,
an `O(2)` isomorphism, or a matching Berry number is supporting evidence only.

## Authoritative dependency graph and immediate actions

```text
P0 action freeze -> P1 vacuum+spectrum -> P2 running+thresholds
P0 + P1 + projected P2 -> full complex loop eigenpair + PS Yukawa matching
corrected eigenpair + PS matching -> P3 inverse scalar/flavor feasibility
P3 failure -> declared alternative -> rerun P1 + P2 -> restart P3
P3 viable -> P4 proton decay+amplitudes
P2 + viable P3 + P4 -> P5 model comparison+conditional 4d closure

U0 -> U1 Spin^c construction -> U2 anomaly+exact spectrum
U1 + U2 -> U3 Yukawa matching -> first-principles family-origin claim
U4 string realization is optional after U2

S0 -> S1 mathematical regularity (optional)
S has no blocking arrow into P or U
```

Next execution order:

1. **P:** the P-SPEC1 common scalar spectrum and source-preserving nested
   elimination tree subgate above is complete. Reuse `(K,J,C)` and the
   interaction-derived residues; do not repeat the unchanged tree audit
   or infer that a massless physical mode must appear in every channel.
   Retain the completed fixed-VEV rule and complete bosonic complex
   doublet matrix; do not repeat the cached Hessian calculation unchanged.
   Reuse the common complex dictionary, six-invariant PS tensor flow,
   upper finite gauge and pure-Yukawa matching, lower covariant tree
   action, scalar kinetic subset, nonzero tree CHN and moving-threshold
   engine. Complete the missing vector/Goldstone/ghost and scalar
   Wilson diagrams, lower hard-minus-EFT matching, and finite sequential
   neutrino matching; carry the actual nonzero CHN through the enlarged
   dimension-five EFT. Use the exact heavy resolvent for the few coupled
   modes requiring it rather than a new global hierarchy condition.
   Only then rerun the physical scales and
   iterate the fitted Schur light eigenpair and constrained likelihood.
   Treat `eta_Y` as a derived correlated profile quantity, not a free scalar.
   The repaired six logarithmic identities are required regression tests;
   they do not replace finite matching or a physical fit.
2. **P alternatives:** only if a baseline repair fails, compare the explicit
   mirror-10/vectorlike completion, a two-doublet EFT or a full enlarged
   scalar model. Do not add these inside an unrestricted nuisance.
   Momentum-dependent pole and proton predictions require viable matched
   flavor matrices; the full Q-block computation itself is now necessary
   diagnostic work and is no longer incorrectly blocked by the tree no-go.
3. **U:** retain the `O(3) gauge line + O(-1) spin bundle -> O(2)`
   representation-valued index/anomaly ledger. Flux selection and physical
   carrier remain independent assumptions to derive.
4. **S:** freeze new scans. GTA/regularity is a separate mathematical
   task with no blocking arrow into P or U.
5. Keep clean-run/source hashes and document registry synchronized.
   Mechanical replay is not a substitute for a physical fit.

## S-layer archive: AP-E18 GTA/minimizing-diagonal checkpoint (2026-07-24; current S authority)

Status: **GTA is not a universal consequence of the AP-E17 additive
minimizing-diagonal estimates.  The standard reverse-Hölder imports fail,
and a boundary-straddling rank-two needle is invisible to those normalized
limits while violating every GTA thickness.  Existence of a scale-polished
recovery sequence remains open.**  The production card passes `14/14`.

### Reverse-Hölder verdict

- The AP-E17 comparison deficit is `o(r_l^3)` on each fixed tangent
  subball.  It is not controlled relative to an unresolved ball
  `s_l << r_l`, so it supplies no scale-uniform Gehring seed.
- The coarse derivative growth is `(p,q)=(2,6)`.  For `n=3,p=2`, the 2024
  relaxed strongly-quasiconvex theorem requires
  `q<min(np/(n-1),p+1)=3`; it cannot be imported.
- `Xi=(Du,sqrt(R)M2,sqrt(K)M3)` obeys curl/div identities and has quadratic
  energy, but admissible graph fields also obey nonlinear Plücker,
  `|u|=1`, and `u^T Du=0` constraints.  The map is not minimizing over the
  full linear `A`-free competitor class required by the 2025 theorem.

### Boundary-straddling rank-two needle

- Set
  `A_e=e^(-1/2)`, `b_e=h_e=e^2`.
  The tangent `L2` mass is `O(e^5)`, while the leading first- and
  second-minor energies are both `O(e)` and the pure third minor vanishes.
- The energy-amplitude correlations do not vanish:
  `int |z_e|^2|Dz_e|^2 ~ 1` and
  `int |z_e|^2|M2(Dz_e)|^2 ~ 1`.
- If the needle straddles the tangent-ball boundary, every shrinking
  thickness fails:
  `GTA(e,delta)>=c/delta^2` for `delta>=h_e`, and
  `GTA(e,delta)>=c/(e^2 delta)` for `delta<h_e`.
  The unweighted annular graph energy still tends to zero in both branches.
- At physical radius `r_e=e`, the target displacement is `sqrt(e)` and the
  physical added energy is `O(r_e^4)=o(r_e^3)`.  A compact chart homotopy
  preserves degree, and the perturbation is invisible to the tangent defect
  measure and normalized additive comparison error.
- Therefore GTA cannot be derived from the additive tangent estimates alone.
  This does not prove that the selected repository sequence contains the
  needle or disprove a specially selected, scale-polished, or exactly
  locally minimizing GTA sequence.

### Optional S-layer continuation

1. Construct a **scale-polished recovery sequence** with
   `epsilon_l(B_s(y)) <= omega_l s^3`, `omega_l->0`, uniformly down to a
   declared microscopic cutoff.
2. Alternatively prove energy-amplitude tightness
   `int |v_l-v_0|^2(1+R|Dv_l|^2+K|M2(Dv_l)|^2)->0`
   together with a one-hemisphere collar.  This condition plus no boundary
   mass yields a slow GTA annulus.
3. If both constructions fail, the relaxation must retain a
   representative-dependent boundary concentration coordinate; AP-E17
   Dirac rigidity remains only conditional.

Canonical artifacts:

- `route_f/tex/ap_e18_gta_minimizing_diagonal.tex` and PDF;
- `route_f/code/verify_ap_e18_gta_minimizing_diagonal.py`;
- `route_f/output/ap_e18_gta_minimizing_diagonal.{json,md}`;
- AP-E18 section in the master Route-E derivation ledger.

Within Layer S, full `mu=0`, local recovery, classicality, isolation, Hessian,
and determinant remain open.  They do not block Layers P or U.  The former
portal dependency is retired by S2.

## AP-E17 boundary-transfer/Dirac/capacity checkpoint (2026-07-24; superseded above)

Status: **normalized full local quasiminimality is proved on the AP-E16
diagonal.  A recovery-compatible sphere-valued boundary modification and
Dirac rigidity are proved under the explicit graph-tight annulus condition
(GTA).  GTA is not derived from the present endpoint graph bounds.  The
nonvolumetric WZ-capacity/topology alternative is proved unconditionally.**
The production card passes `12/12`.

### Volumetric transfer and rigidity

- The global minimizing excess already gives
  `epsilon_(j_l)(B_(rho r_l))/mu(B_(r_l))->0` for every fixed `rho<1`.
  This is full relative-homotopy quasiminimality, not phase stationarity.
- GTA requires a one-hemisphere transition layer, vanishing annular graph
  energy, and the exact weighted cutoff remainder tending to zero.
- Under GTA,
  `P[q+r_l(eta_l v_l+(1-eta_l)F_0 y)]` equals the original map near the
  boundary and preserves relative homotopy and degree.  Its normalized
  energy tends to `|B_1|W(F_0)`.
- Null-Lagrangian moments give the exact generalized Young-measure variance
  identity.  The affine comparison forces every ordinary graph variance and
  the complete concentration mass to vanish:
  `nu=delta_(F_0), lambda=0`.
- Hence the volumetric part of `mu` vanishes at every GTA point.

### Exact remaining endpoint blocker

- GTA does not follow from strong `L2` convergence plus bounded graph energy.
  A rank-two microball with radius `h` and amplitude `h^(1/4)` has
  `int|z_h|^2~h^(7/2)`, `int|Dz_h|^2~h^(3/2)`,
  `int|M2(Dz_h)|^2~1`, `M3=0`, but its thickness-`h` weighted gluing term is
  `h^(-3/2)`.
- Classical gradient decomposition/equiintegrability controls one
  `L^p` gradient scale.  It does not preserve the simultaneous endpoint
  `L2` graph of `Du,M2,M3`, the `S3` target, and relative degree.
- The next mainline theorem is therefore GTA for the actual minimizing
  diagonal, not another abstract Young-measure classification.

### Nonvolumetric alternative

- With `b(r)=(2pi^2)^(-1) int_(B_r) J_u`,
  `E_3(B_r)>=3K pi^3 |b(r)|^2/(2r^3)`.
- If local energy is `O(r^d)`, then
  `|b(r)|=O(r^((d+3)/2))->0`.  A nonzero integer point degree is impossible.
  The Hopf base also has zero monopole charge because `da` is exact on the
  ball.
- A normalized tangent is either WZ-capacity-active, with a positive sextic
  fraction but vanishing absolute charge, or WZ-neutral.  In the neutral
  branch any remaining defect is analytic, not protected by Route-E
  topology.

Canonical artifacts:

- `route_f/tex/ap_e17_boundary_transfer_dirac_capacity.tex` and PDF;
- `route_f/code/verify_ap_e17_boundary_transfer_dirac_capacity.py`;
- `route_f/output/ap_e17_boundary_transfer_dirac_capacity.{json,md}`;
- AP-E17 section in the master Route-E derivation ledger.

Full `mu=0`, local recovery, classicality, isolation, Hessian, determinant,
and portal remain closed.

## AP-E16 tangent/Caccioppoli/Young-measure checkpoint (2026-07-24; superseded above)

Status: **the common `mu`-a.e. tangent and a conditional purified Hopf-base
Caccioppoli inequality are proved.  Unconditional excess contraction from
phase purification plus weighted cutoff decay is disproved by an exact
homogeneous pure-base Young measure.  The counterexample is not
asymptotically minimizing, so the minimizer theorem remains open but its
missing hypothesis is now identified: normalized full local quasiminimality
must survive the blow-up.**  The production card passes `12/12`.

### Closed tangent statements

- The Hopf energy split defines five nonnegative defect measures
  `mu_a,mu_n,mu_an,mu_c,mu_cs`, whose sum is `mu`.  At `mu`-a.e. `x_0`,
  every component has the same normalized tangent `tau` with
  Radon--Nikodym weight `theta_beta(x_0)`.
- A diagonal recovery sequence realizes those component tangents and obeys
  `Delta_j(B_(Lr))/mu(B_r)->0` for each fixed `L`.  This proves normalized
  stationarity against compact exact fibre shifts, not strong convergence
  of `a_j`.
- If the target amplitude is `s=r^alpha` and `mu(B_r)~r^d`, simultaneous
  normalization of first-, second-, and third-minor energies requires
  `2alpha+1=4alpha-1=6alpha-3=d`.  The unique solution is
  `(alpha,d)=(1,3)`.  Nonvolumetric tangents must retain order-specific
  concentration coordinates.

### Conditional local inequality and exact no-go

- In one normal target ball, the geodesic cutoff
  `v=exp_q((1-eta)log_q u)` gives
  `Phi_j(s)<=theta Phi_j(t)+weighted cutoff+potential+C epsilon_j(B_t)`,
  where `theta<1` and `epsilon_j` is the full relative-homotopy local
  comparison deficit.  All rank-one second- and third-minor cutoff terms are
  displayed explicitly in the AP-E16 ledger.
- The exact maps
  `u_N=(sqrt(1-N^-2(sin^2 Nx1+sin^2 Nx2)),0,
  N^-1 sin Nx1,N^-1 sin Nx2)` lie in `S3`, have degree zero, and remain
  nontrivial after exact phase minimization.  They generate
  `F_31=cos(theta_1), F_42=cos(theta_2)`, with
  `<F>=<M2(F)>=0` but graph-energy gap `1/2+R/8`.
- Their weighted cutoff tends to zero while their Hopf-base curvature
  persists.  Thus the AP-E15 two-limit weighted estimate is necessary for
  the cutoff argument but not sufficient for `mu=0`.
- This is not a minimizing counterexample: the constant map removes the
  positive local energy.  Its role is to prove that an excess theorem must
  use normalized full local quasiminimality, not only phase stationarity.

### Ordered continuation

1. **Volumetric branch:** prove a recovery-compatible boundary-modification
   theorem transferring
   `epsilon_(j_l)(B_(r_l))/mu(B_(r_l))->0` to the tangent Young measure.
   Combine it with the exact strong-quasiconvex identity to force the tangent
   to be Dirac.
2. **Nonvolumetric branch:** retain the component tangents and prove a
   concentration-capacity or topological lower-bound alternative.
3. Only after both branches exclude nonzero defect may full `mu=0`, local
   recovery, regularity, isolation, and the same-action bosonic Hessian be
   promoted.
4. Parallel Dirac/Callias operator mathematics remains allowed.  Determinant
   and degree-one portal promotion remain false.

Canonical artifacts:

- `route_f/tex/ap_e16_tangent_caccioppoli_young_measure.tex` and PDF;
- `route_f/code/verify_ap_e16_tangent_caccioppoli_young_measure.py`;
- `route_f/output/ap_e16_tangent_caccioppoli_young_measure.{json,md}`;
- AP-E16 section in the master Route-E derivation ledger.

## AP-E15 relaxation-defect/fibre-purification checkpoint (2026-07-23; superseded above)

Status: **the AP-E11 action is frozen and no further lattice scans are
authorized.  The relaxation defect measure is defined, and every
compact-fibre-phase-reducible vertical defect is proved zero.  Weighted
cutoff decay is proved at Lebesgue-a.e. centres, but not yet on possible
singular defect support.  Full `mu=0`, classicality, isolation, Hessian,
determinant, and portal remain closed.**  The continuum card passes `11/11`.

### Closed continuum statements

- For `Xi=(Du,sqrt(R)M2,sqrt(K)M3)`,
  `|Xi_j|^2 dx/2 ->* |Xi_u|^2 dx/2+mu` defines a nonnegative Radon measure.
  Its mass is the relaxed-minus-naive gap, and `mu=0` is equivalent to strong
  complete-minor graph convergence of the recovery subsequence.
- For compact phase variations, `n` and `da` are unchanged and
  `a -> a+2dchi`.  The exact phase Hessian is bounded below by
  `(1-m^2 r^2/pi^2)||dphi||_2^2` on `B_r`.
- If `E(u_j)->I_B`, the local optimal phase drop obeys
  `Delta_j(B_r)<=E(u_j)-I_B->0`; strong convexity then gives
  `||dchi_j||_2->0`.  This proves that the phase-reducible vertical component
  is zero, without claiming that all mixed `a_j,Dn_j,da_j` defects vanish.
- At almost every Lebesgue point of a fixed graph map,
  `r^-2 int_Br |u-u(x)|^2|M2|^2=O(r^3)+o(r)->0`.

### Remaining theorem and strict ordering

1. Prove the sequence-uniform two-limit estimate
   `lim_(r->0) limsup_j r^-2 int_Br |u_j-q_(j,r)|^2|M2(Du_j)|^2=0`
   for `mu`-a.e. points, or a scale-invariant excess contraction which implies
   it.  Fixed-map a.e. differentiation neither exchanges the recovery limit
   nor excludes a diffuse base/mixed defect.
2. Use that estimate to prove full `mu=0`, then local graph recovery,
   regularity, and continuum isolation for the selected `B=1` minimizer.
3. Only after classicality and isolation close may the same-action bosonic
   Riemann Hessian be assembled.
4. Dirac/Callias operators on a prescribed smooth representative may be
   studied as parallel mathematics.  Determinant and portal promotion remain
   false.

Canonical artifacts:

- `route_f/tex/ap_e15_relaxation_defect_fibre_purification.tex` and PDF;
- `route_f/code/verify_ap_e15_relaxation_defect_fibre_purification.py`;
- `route_f/output/ap_e15_relaxation_defect_fibre_purification.{json,md}`;
- AP-E15 section in the master Route-E derivation ledger.

## AP-E14 WZ/Morrey and Hopf div--curl checkpoint (2026-07-20; superseded above)

Status: **WZ-flux decay is proved for every fixed complete-minor graph map,
and therefore for the selected relaxed `B=1` representative.  This decay
does not imply the required linear Morrey oscillation.  The exact Hopf
base/fibre split is closed, while a universal exact div--curl collapse of the
full positive cutoff remainder is disproved.  Hessian, determinant, and
portal remain closed.**  The production card passes `12/12`.

### Closed theorem: fixed-map WZ decay

- With `J_u=det[u,du]` and `|J_u|<=|M3(Du)|`, every fixed graph map obeys
  `sup_x dist(int_Br J_u,2 pi^2 Z)^2/r^3 <= (4 pi/3)
  sup_x int_Br |M3|^2 -> 0`.
- The convergence is uniform on compact interior sets by absolute continuity
  of the one fixed `L1` density `|M3|^2`.  Minimality is not used.
- AP-E13 is not a contradiction: its counterexample chooses a different,
  non-equiintegrable map at each scale.

### Morrey and Hopf decisions

- `u=(cos ell(r),sin ell(r),0,0)`, `ell=log log(e/r)`, has finite Dirichlet
  tail `2 pi/log(e/R)`, rank-one derivative, `M2=M3=WZ=0`, but oscillation
  diameter two at every axis scale.  It can be localized in a vacuum cylinder
  without changing an exterior smooth `B=1` degree.  Thus graph energy,
  degree, and WZ decay do not prove linear Morrey oscillation; actual relaxed
  minimality remains essential.
- For Hopf projection `n=h(u)` and contact form `a=2u*theta`, with the stated
  orientations, `da=-n*omega_S2` and
  `|Du|^2=(|a|^2+|Dn|^2)/4`,
  `|M2|^2=(|a wedge Dn|^2+|da|^2)/16`,
  `|M3|^2=|a wedge da|^2/64`, and
  `B=(16 pi^2)^(-1) int a wedge da`.
- Only the signed Chern--Simons integral changes by an exact form under fibre
  phase.  A periodic rank-one-base lift has `da=a wedge da=0` but
  `int |M2|^2=A^2 B^2 pi^3/2>0`, disproving universal full cancellation of
  the positive cutoff remainder.

### Numerical scope and ordered continuation

- Three newly relaxed compatible-action backgrounds (`N=17,21,25`) are
  stationary, admissible, and `B=[1,1,1]`.  Maximum edge-angle quotients
  decrease `3.7669 -> 3.4636 -> 3.3162`; the exact normalized-affine local WZ
  Cauchy ratio never exceeds `0.972327697`.  This is evidence only.
- Next derive an inner-variation/monotonicity identity for the lower-
  semicontinuous relaxation with explicit stress defect measure.  Then prove
  its rank-one vertical--horizontal Hopf defect vanishes, or remove it by a
  degree-preserving local comparison, before attempting epsilon regularity
  and an `L-infinity` Morrey upgrade.
- Same-action Hessian, determinant variation, and degree-one portal stay
  embargoed until local recovery, classicality, and continuum isolation close.

Canonical artifacts:

- `route_f/tex/ap_e14_wz_morrey_hopf_divcurl.tex` and PDF;
- `route_f/code/verify_ap_e14_wz_morrey_hopf_divcurl.py`;
- `route_f/output/ap_e14_wz_morrey_hopf_divcurl.{json,md}`;
- AP-E14 section in the master Route-E derivation ledger.

## AP-E13 annular replacement/reverse-Hölder checkpoint (2026-07-20; superseded above)

Status: **the unconditional good-sphere annular replacement is disproved by
an explicit Wess--Zumino filling obstruction.  A scale-compatible conditional
replacement is proved.  Exact strong quasiconvexity alone does not close the
parallel reverse-Hölder route.  Endpoint density itself remains undecided;
the Hessian, determinant, and portal gates remain false.**  The production
card passes `9/9`.

### Exact no-go theorem

- For `g_eps(omega)=(cos eps,sin eps omega)` and the radial-constant shell map
  on `A_(rho,2rho)`, the raw minor norms are
  `E1=8 pi rho sin^2(eps)`, `E2=2 pi sin^4(eps)/rho`, and `E3=0`.
- Every trace-preserving filling has Wess--Zumino volume
  `V(eps)+2 pi^2 k`, where
  `V(eps)=2 pi(eps-sin(2eps)/2)~4 pi eps^3/3`.  Hence
  `int_B |M3(DU)|^2 >= V(eps)^2/|B|`.
- At `eps=sqrt(rho)`, the shell graph energy tends to zero but the filling
  lower bound tends to `4 pi/3`.  Every sphere has the same trace, so choosing
  a good radius cannot help.  The filling densities keep positive mass on
  shrinking balls and therefore fail `L2` equiintegrability.
- More generally, for `eps=rho^alpha`, shell energy vanishes for
  `alpha>1/4`, while the WZ-forced cubic-minor filling cost vanishes only for
  `alpha>1/2`.  This is a genuine annular window, distinct from AP-E12's
  excluded full-core Malý dipole.

### Repaired theorem and reverse-Hölder boundary

- If the good-sphere trace lies in one normal ball and
  `||log_q g||_infinity <= Lambda rho`, geodesic contraction preserves the
  trace and gives
  `I1<=C(rho T1+Lambda^2 rho^3)`,
  `I2<=C(rho T2+Lambda^2 rho T1)`, and
  `I3<=C Lambda^2 rho T2`.  Coarea converts `rho Tk` to shell costs.  This
  closes the conditional target-valued endpoint replacement.
- A necessary weaker datum is
  `dist(V_WZ(g),2 pi^2 Z)^2/rho^3 -> 0`.  Shell graph energy alone does not
  imply it.
- The stress cutoff contains
  `rho^(-1)|u-q||M2||M3|`.  Young's inequality leaves
  `rho^(-2)|u-q|^2|M2|^2`, which requires the desired Morrey oscillation or
  prior higher integrability.  Therefore exact strong quasiconvexity is not
  an independent reverse-Hölder proof.

### Ordered continuation

1. Prove WZ-flux decay plus scale-compatible trace oscillation for the
   selected relaxed `B=1` minimizer, using a monotonicity/frequency formula.
2. In parallel, derive a genuinely compensated Noether-current/
   `S3=SU(2)` Maurer--Cartan reverse-Hölder estimate that removes the cutoff
   remainder without assuming its conclusion.  The 2026 Hopf-lifting theorem
   suggests a concrete variable split: `S3` field -> `S2` Hopf projection plus
   vertical gauge one-form; test whether the remaining cutoff product is an
   exact/div--curl pairing in those variables.
3. Only after local recovery, classicality, and continuum isolation close may
   the same-action Riemann Hessian and regulated determinant be assembled.

Canonical artifacts:

- `route_f/tex/ap_e13_annular_replacement_reverse_holder.tex` and PDF;
- `route_f/code/verify_ap_e13_annular_replacement_reverse_holder.py`;
- `route_f/output/ap_e13_annular_replacement_reverse_holder.{json,md}`;
- AP-E13 section in the master Route-E derivation ledger.

## AP-E12 endpoint density/regularity checkpoint (2026-07-20; superseded above)

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
   parallel; historically the degree-one Route-E portal remained last.  The
   required portal is now retired by S2.

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
   mapping-torus mod-two index as parallel lanes.  Historically the
   degree-one Route-E portal remained last; the required portal is now
   retired by S2.

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

## Legacy single-chain crosswalk (retained for provenance)

The F0--F11 plan below predates the three-layer decision.  Its calculations
and acceptance tests remain useful, but its ordering is superseded:

- `F1--F8` map to Layer P;
- `F9--F10` map to Layer U;
- the AP-E/GTA program maps to Layer S and has no blocking arrow into either;
- `F0-A` is parallel reproducibility work;
- `F0-C` Route-B messenger work is optional unless P0 selects a branch that
  actually uses it.

### Legacy definition of done

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

### Legacy dependency graph (superseded by the authoritative graph above)

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

## Legacy P0: historical blockers (remapped to P/U or parallel evidence)

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

Compare, using identical conventions, with F-54 now the primary baseline:

- F-54 primary: non-SUSY `54_H + 10_H,C + 126_H`;
- F-210 comparison/rescue: non-SUSY
  `210_H + 10_H + 120_H + overline{126}_H`;
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

## Legacy P1: historical physical closure (now Layer P)

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

## Legacy P2: historical optional UV and publication (now Layer U)

### F9 Physical origin of the family carrier

Status: `open`.

Two creative options may be developed in parallel:

- a 6D `Spin(10) x U(1)_F` flux model on `P1`;
- a 4D defect/deconstructed index model producing three localized `16`s.

For the 6D branch, distinguish the microscopic gauge line from the effective
positive-spinor bundle:

`K_P1^(1/2)=O(-1)`, `L_F=O(3)`, and
`K_P1^(1/2) tensor L_F=O(2)`.  The first target is therefore an anomaly-safe
degree-three gauge flux whose spin-twisted kernel is the Route-E triplet, not
an exactly-two colour/orbital rule inferred from the final `O(2)` label.

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

## Legacy immediate next action (superseded 2026-08-30)

Run the 21-node full numerical DAG (19 recovered lanes plus DYN-5V/DYN-7F
guards) from a clean clone with the now-present RE-SC3/4/5 cards, and compare
all source/input/output digests with the current ledger set.  In parallel,
construct an explicit
interacting messenger/selection sector for DYN-5 and supply branch-local
thermal inputs for a two-flavor or density-matrix DYN-7/9b-3 calculation.
F0-B may proceed independently.  Do not spend more compute on new flavor or
proton scans until F1 fixes the action and F2 exports the actual spectrum.
The current immediate actions are the P/U/S list in the authoritative section
at the top of this file; in particular the clean DAG is parallel evidence
work, not the next physics milestone.
