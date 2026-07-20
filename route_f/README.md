# Route F: Theory Closure and Evidence Integration

Route F is the integration route for Routes A--E.  It is not a sixth claim of
new fundamental physics.  Its job is to turn the useful representation,
geometry, messenger, amplitude, dynamics, and string ideas into one auditable
theory with one action and one evidence chain.

The review found three immediate blockers:

1. Route B's stated `U(1)_R` rule does not forbid the post-breaking operator
   `X L H_u`, and integrating out a canonical propagating messenger also gives
   a tree-level Kahler correction that is much larger than the quoted loop
   estimate.
2. Route E's original Cartan-selection proof did not exclude the genus-one
   branch: a one-dimensional abelian Lie algebra admits a non-degenerate
   symmetric ad-invariant bilinear form even though its Killing form vanishes.
   This logic blocker is repaired by proving only `N_fam<=3` from original
   H3 and making `N_fam=3` conditional on an explicit H3+
   Killing-contact axiom; a physical origin for H3+ remains open.
3. Route E's 19 DYN source scripts and Route-D pricing cards RE-SC3/4/5 are
   now tracked.  This closes the file-existence blocker, not the physics:
   RE-SC3/4/5 remain unpromoted; DYN-9b-2 is not a global non-SUSY flavor fit;
   DYN-9b-3 lacks tau-resolved kinetics; and DYN-8 is only a fail-closed
   collector.  Route C and D likewise still contain classification/arithmetic
   gates rather than complete amplitudes or a global compactification.  The
   evidence card checks a 29-item dynamics P0 subset plus 4 core artifacts,
   not every path mentioned anywhere in the roadmap.

## 2026-07-20 AP-E14 WZ/Morrey and Hopf div--curl audit

AP-E14 proves the Wess--Zumino decay requested after AP-E13 for every fixed
complete-minor graph map, hence also for the selected relaxed `B=1`
representative.  The proof is the uniform absolute continuity of
`|M3(Du)|^2`: local squared flux divided by `r^3` tends to zero without using
minimality.

The stronger linear Morrey estimate remains open.  A localized rank-one
circle-valued defect has finite Dirichlet energy, zero higher minors and zero
WZ flux, yet oscillates across the full circle at every scale.  Therefore the
missing estimate must exploit actual relaxed minimality, not only graph
bounds and topology.

The Hopf projection and fibre form give exact identities for all three minor
orders.  They show that signed Chern--Simons flux has exact phase
cancellation, but a positive vertical--horizontal shear survives even when
all curl/Chern--Simons terms vanish.  This disproves universal full div--curl
cancellation of the cutoff remainder.  The `12/12` card and three re-relaxed
grids provide algebraic and numerical checks without promoting continuum
regularity.  Hessian, determinant, and portal remain closed.  See
`tex/ap_e14_wz_morrey_hopf_divcurl.tex`,
`code/verify_ap_e14_wz_morrey_hopf_divcurl.py`, and
`output/ap_e14_wz_morrey_hopf_divcurl.{json,md}`.

## 2026-07-20 AP-E13 annular replacement and reverse-Hölder audit (superseded by AP-E14 for WZ/Morrey)

AP-E13 disproves the unconditional annular endpoint-replacement strategy with
an exact `S3` geodesic-cap family.  At cap radius `eps=sqrt(rho)`, the shell's
complete-minor energy tends to zero, while relative Wess--Zumino volume forces
every trace-preserving filling to carry asymptotic third-minor `L2` cost
`4 pi/3`.  Since the trace is identical on every shell sphere, good-radius
selection cannot repair it, and the filling family is not equiintegrable.

The corrected theorem assumes scale-compatible target oscillation
`||log_q g||_infinity<=Lambda rho`.  Its geodesic contraction controls all
three filling-minor norms by the good-sphere trace norms and hence by the
shell energy.  The reverse-Hölder alternative remains open: its highest
stress cutoff leaves a term requiring that same Morrey bound or already-known
higher integrability.  The `9/9` card therefore changes the blocker without
promoting endpoint density, classicality, Hessian, determinant, or portal.
See `tex/ap_e13_annular_replacement_reverse_holder.tex`,
`code/verify_ap_e13_annular_replacement_reverse_holder.py`, and
`output/ap_e13_annular_replacement_reverse_holder.{json,md}`.

## 2026-07-20 AP-E12 endpoint graph-density and regularity audit (superseded by AP-E13 for the annular gate)

AP-E12 attacks the last AP-E11 classical-continuum gate.  It proves three
useful facts: graph-norm convergence makes the degree continuous; the
complete-minor density has uniform Legendre--Hadamard constant one and an
exact strong-quasiconvexity identity; and the direct radial Malý dipole has no
endpoint scaling window.  For `r^alpha gamma(theta)`, finite `L2` second-minor
energy needs `alpha>1/2`, while non-vanishing forced filling cost needs
`alpha<=1/2`.

The endpoint density theorem is still neither proved nor disproved.  Known
Cartesian smooth approximation loses exponent, and the fallback regularity
theorem for relaxed quasiconvex minimizers does not reach this model's
dimension-three `(2,6)` growth or its sphere/degree constraints.  The `9/9`
card finds decreasing cell concentration and a tight finite-grid basin after
translation/target-`SO(3)` quotienting, but does not promote continuum
isolation.  Weak Euler--Lagrange, partial regularity, classicality, Hessian,
determinant, and portal gates remain false.  See
`tex/ap_e12_graph_density_regular_minimizer.tex`,
`code/verify_ap_e12_graph_density_regular_minimizer.py`, and
`output/ap_e12_graph_density_regular_minimizer.{json,md}`.

## 2026-07-19 AP-E11 compatible tetrahedral cochain action (superseded by AP-E12 for the classical gate)

AP-E11 replaces the AP-E10-refuted one-corner stencil by one compatible
discrete-de-Rham evidence chain.  Primal `dn` is combined with the
Alexander--Whitney front/back cup product, which obeys the exact graded
Leibniz rule.  Full antisymmetrization reproduces the affine two- and
three-minor cochains.  Positive Whitney mass matrices provide the calibrated
Hodge stars, and radial pullback of the same triple cup gives the exact
normalized-affine degree on every admissible tetrahedron.

The periodic cell problem is closed analytically: all complete-minor means
are periodic null Lagrangians, so convexity makes the zero corrector global
and gives the mesh-independent density
`|A|^2/2+R|M_2(A)|^2/2+K|M_3(A)|^2/2`.  The 13-case unanchored production scan
passes the declared joint-limit, translation, and five-/six-tet quotient
thresholds; all fields are stationary, admissible, and `B=1`.

The rigorous continuum conclusion is deliberately narrower.  Gamma
convergence is proved to the fixed-degree lower-semicontinuous relaxation.
Smooth fixed-degree density in the full minor graph norm (or regularity and
isolation of the selected relaxed minimizer) remains unproved.  Consequently
the relaxed regulator gate passes, but the classical background, Hessian,
determinant, and portal gates remain closed.  See
`tex/ap_e11_compatible_cochain_action.tex`,
`code/scan_ap_e11_compatible_cochain_action.py`, and
`output/ap_e11_compatible_cochain_action.{json,md}`.

## 2026-07-18 AP-E10 stencil no-go, cell formula, and translation quotient (superseded by AP-E11)

AP-E10 proves that the AP-E7 one-corner forward Skyrme stencil cannot support
the proposed continuum proof.  Exact three- and two-periodic sphere-valued
sequences converge uniformly to the vacuum while retaining respectively a
nonzero two-minor and a nonzero diffuse forward topological current.  The
Skyrme energy controls the current in `L^(4/3)`, so the failure is a diffuse
lattice alias rather than singular concentration.

The finite-`R` periodic cell problem is also solved.  Direction-wise cycle
balance and Jensen prove the zero corrector exactly for both six- and
five-tet graphs.  Their homogenized quartics remain nonproportional; finite-
`R` triangulation dependence is therefore physical to this regulator.

An 11-case production scan implements a barycentre-aligned translation
quotient, dynamic degree-preserving centering, and the joint range
`a=0.25 -> 0.208333`, `L=6 -> 7.5`.  All backgrounds are stationary,
admissible, and `B=1`, but centered-size, radial-profile, and cross-mesh gates
fail.  The same-action Hessian and regulated determinant variation are not
built.  The next mainline is a compatible tetrahedral cochain/cup-product/
Hodge-star Skyrme action, not another refinement of the one-corner stencil.
See `tex/ap_e10_compactness_homogenization_centering.tex`,
`code/scan_ap_e10_compactness_homogenization_centering.py`, and
`output/ap_e10_compactness_homogenization_centering.{json,md}`.

## 2026-07-18 AP-E9 scaling and regulator-stability result (superseded by AP-E10)

AP-E9 proves fixed-box strong-`L2` equicoercivity but disproves the stronger
claim that the AP-E8 barrier can both vanish in the continuum and by itself
close degree sectors.  With `d=1-epsilon`, `w=gamma*a`, and
`R=gamma*a^2/d^2`, uniform edge-floor separation requires `inf w>0`, while
the smooth barrier disappears only for `R->0`.  For fixed epsilon and
`gamma=c*a^(-p)`, the unique overlap is `1<=p<2`.  For `inf R>0` the
barrier gives `W1,4` compactness and closes degree.  At critical `p=2`, the
raw six-tet/five-tet Cauchy--Born quartics are not scalar multiples; their
finite-`R` Gamma densities require a still-unevaluated homogenized cell
problem, so no continuum-density identification is claimed.

The reproducible 23-case scan tests six spacings, four volumes, three scaling
exponents, two checkerboard five-tet phases, the original six-tet mesh, and
centered/translated starts at two spacings.  All cases are stationary,
admissible `B=1` backgrounds, but the barrier-slope, volume-radius,
cross-triangulation, and regulator-trajectory gates fail.  Consequently the
continuum background, full same-action Hessian, and regulated determinant
variation remain fail-closed.  See
`tex/ap_e9_gamma_limit_regulator.tex`,
`code/ap_e9_triangulation_tools.py`,
`code/scan_ap_e9_gamma_triangulation.py`, and
`output/ap_e9_gamma_triangulation.{json,md}`.

## 2026-07-17 AP-E8 topology-preserving finite-grid result (superseded by AP-E9 for continuum scaling)

AP-E8 selects the topology-preserving-action route and genuinely changes one
AP-E7 blocker.  A subtracted logarithmic barrier on all unique Freudenthal
co-cell pairs gives an exact normalized-affine no-zero bound, separates the
finite-energy configuration space into degree components, and proves an
interior unanchored minimizer in every nonempty component.  Sampling a smooth
compact-supported degree-one field makes `B=1` components nonempty on all
sufficiently fine grids.

The production card passes `13/13`.  Seven no-pin final fields at four
spacings and three volumes all have direct tangent-plus-scalar gradient below
`7.8e-7`, positive pair margins, and three-target `B=[1,1,1]`.  Four nearby
barrier choices and one perturbed start give the same finite-grid verdict.
The exact action calculus, differentiable solver continuation, and intrinsic
barrier Hessian pass independent directional checks.

This is not yet a continuum or physical-lattice result.  The weighted radius
has not volume-converged, gamma dependence is strong, and no equicoercivity,
Gamma convergence, regulator independence, physical QC2D action, determinant
Hessian, BRST superdeterminant, HMC, or quantum continuum limit is present.
Thus `C_B1_finite_grid=true` but `C_B1_physical_continuum=false`; every complete
pre-portal lane, portal authorization, and physics promotion remain false.
See `tex/ap_e8_topology_preserving_b1.tex`,
`code/scan_ap_e8_topology_preserving_b1.py`, and
`output/ap_e8_topology_preserving_b1.{json,md}`.

## 2026-07-17 AP-E7 result (superseded by AP-E8 for finite-grid Lane B)

The requested common-regulator, discrete re-relaxation, and actual-moduli
alternatives have been executed and integrated.  The outcome is intentionally
fail-closed:

- A common quadratic APS/PV gauge--scalar--ghost--fermion prescription and the
  pure two-flavour heavy gauge/gravity threshold are now explicit.  The pure
  phase is exactly `+1`.  The conditional Higgsed `(+,+)` source is not the
  still-missing original target-family mapping-torus computation.
- A theorem shows that the unrestricted finite-site `S3` configuration space
  is connected.  All five unconstrained grid/volume solves unwind to vacuum.
  The positive fixed-core spectra are exact only in that artificial restricted
  sector; the assembled fermion normal and free gauge/ghost blocks are not a
  physical super-Hessian.
- The actual `SO(3)` orbit has one nontrivial `Z2` FR line, but topology alone
  does not select the microscopic Pfaffian class.  A minimal rank-three mass
  family exactly realizes `c1=x+2y` and passes numerical Chern/gap checks, yet
  remains a spectator bundle rather than a same-soliton Yukawa derivation.

The three AP-E7 cards pass `42/42`, `16/16`, and `18/18`; the independent
execution-frontier card passes `16/16` while keeping all three complete lanes,
portal construction, and physics promotion false.  See
`tex/another_physics_route_e_derivation_ledger.tex` for the combined formulas
and `tex/ap_e7_*` plus `output/ap_e7_*` for the standalone proofs and cards.

## 2026-07-16 AP-E6 same-background result (superseded by AP-E7)

The next execution round is fail-closed and more restrictive than the prior
frontier templates.

- `ap_e6_sp4_eta_bordism_threshold` passes **39/39** and proves
  `Omega_5^Spin(BSp(2))=Z2` with characters `4:-1`, `5:+1`, `10:+1`.
  The full Euclidean regulator, both target mapping-torus signs, and the
  heavy/gravitational APS ratio remain open.  The checked `+2-2=0` is a
  conditional charge-trace obstruction, not a computed determinant.
- `ap_e6_relaxed_b1_multigrid_hessian` passes **20/20** and solves a genuine coupled continuum
  `B=1` BVP and publishes one canonical `(r,F,s)` checksum.  Its sparse
  constrained/projected cubic `n+s` Hessian is analytically exact, while the
  separate diquark/gauge/ghost blocks are not one aggregate matrix.  The
  field is evaluated off shell: the cubic gradient is still `O(.4)` and the
  lattice degree is not converged.  Its negative curvatures are diagnostics,
  not physical instability modes.
- `ap_e6_same_soliton_yukawa_callias` passes **25/25** and consumes that exact checksum.  It finds
  the physical localized orbit `SO(3)`, not `CP1`; the scalar `CP1` norm grows
  with volume, the static Callias index is zero, and a new `c1^2` theorem
  rules out the desired `x+2y` line in a trivial two-band mass family.  Its
  differential character is only a correctly typed restricted ansatz; the
  physical CPT regulator and gauge-basic microscopic descent remain open.

The AP-E6 frontier verifier passes **15/15** and independently recomputes all three conjunctions
as false.  No degree-one portal is started.  Full formulas, numerical
anchors, limitations, and suggested repairs are synchronized in
`ROADMAP.md` and `tex/another_physics_route_e_derivation_ledger.tex`.

## 2026-07-16 AP-E completion frontier (superseded above)

The latest research bundle adds three standalone derivations and one
integrated promotion gate:

- 'ap_e5_4d_qc2d_lattice_hessian' defines the four-dimensional lattice
  action and full formal graded Hessian, with a **24/24** finite
  frozen-background diagnostic;
- 'ap_e3_aps_global_form_search' computes the two product/defect mod-two
  character tables and finds simply connected 'Sp(4)=Spin(5)' as the best
  simple charge/global-form candidate in the scanned set (**58/58**);
- 'ap_e4_same_soliton_callias_descent' derives the conditional
  families-Callias determinant, fixed-polarization CPT/Serre map, and
  differential-character fiber integration/descent gate (**27/27**);
- 'verify_ap_e5_completion_frontier.py' integrates them (**20/20**).

The positive results are exact only within their declared levels. The lattice
calculation is not Monte Carlo or the full nonlinear/quantum Hessian; the
defect characters do not select a unique charged-QC2D Dai--Freed regulator;
the 'Sp(4)' split is tree-level and unprotected; and the Callias 'O(2)',
CPT 'O(-4)', and WZW degree-two line have not been derived from one physical
'B=1' mother model. Accordingly all three complete-lane booleans are false,
the degree-one Route-E portal has not been built, and
'physics_promotion_allowed=false'.

The canonical rollback trail is in
'tex/another_physics_route_e_derivation_ledger.tex'; current priorities are
the complete 'Sp(4)' eta/threshold/bordism calculation, a multi-spacing
dynamical charged-QC2D plus exact projected Hessian calculation, and the
same-soliton Yukawa--Callias/equivariant-WZW derivation. Portal construction
starts only after at least one of those lanes closes.

## 2026-07-14 Another-Physics / Route-E bridge status

The bridge audit now separates exact mathematics from added physics.  On the
H3+-selected CP1/O(2) branch, the binary-quadratic discriminant, the
sl2 Killing norm, the Route-E `K_tr` norm, and the existence of two distinct
centers are exactly equivalent.  Hopf reduction also gives a safe bounded
polarity without negative energy.  By contrast, the proposed sextic Q-ball
and the dynamical `zeta_eff = g^2 Phi^2/Lambda^2` messenger portal are
conditional constructions: exact U(1) symmetry hides the Q-ball phase, so an
observable phase needs a second spurion and must pass charge-leakage,
selection-rule, matching, stability, and phenomenology gates.

The expanded independent bridge verifier passes `27/27`
mathematical/numerical checks
and deliberately records `physics_promotion_allowed=false`.  The complete
formula ledger and the ordered AP-E1--AP-E10 program are canonical in
`tex/another_physics_route_e_derivation_ledger.tex` and `ROADMAP.md`.

## 2026-07-14 AP-E1 projective-doublet result

AP-E1 is closed for the quotient geometry and remains open for physical
selection.  A local fixed-norm charge-one doublet gives the exact Hopf quotient
`S3/U(1)=CP1` and the Fubini--Study kinetic term.  A physical global phase does
not: fixed-charge Routh reduction gives a monopole-coupled `T*CP1` rotor, whose
full Hilbert space contains higher Landau levels.  The finite-dimensional
`H0(CP1,O(k))` carrier therefore requires a first-order or controlled LLL gate.

The minimal Hopf bundle has Chern number one, while `T_CP1=O(2)`.  More
decisively, identifying the projective level with the existing `Q=10^6`
benchmark would give `dim H0(O(k))=1,000,001`, not three.  The recommended
continuation keeps the macroscopic Q-ball charge separate from an independent
level-two projective worldline sector.  Its level rule, LLL gap, chiral-family
interpretation, and full soliton stability remain fail-closed AP-E3--AP-E5
gates.  The AP-E1 arithmetic/source regression verifier passes `30/30` checks
and still records
`physics_promotion_allowed=false`.

## 2026-07-14 AP-E2/AP-E3 execution result

AP-E2 is complete as an exact conditional regression, not as a physical
derivation.  `verify_ap_e2_discriminant_regression.py` passes `30/30` and
preserves

```text
B_Kill(A_q,A_q)=2 Delta(p)=2 sqrt(3) x^T K_tr x
A_q^2=Delta(p) I/4
(p,p)_2=-2 Delta(p)
```

across exact rational and 100-decimal complex samples, `SL(2)` transformations,
finite/infinite roots and double roots, the zero-section exception, both basis
normalizations, and wrong-convention negative controls.  It does not derive
H3+, a microscopic carrier, dynamics, or `|k|=2`, and sets
`physics_promotion_allowed=false`.

AP-E3 now has a conditional microscopic *candidate*.  Two declared
Bose-Hubbard orbitals with
`H_Mott=sum_r[U N_r(N_r-1)/2-mu N_r]` carry one doublet each only in the bare
onsite window `0<mu<U`.  For the complete declared `H_portal=0` model, put
`C=mu+h/2+J_H/4`; a concise sufficient interacting `(1,1)` plateau condition
is `C<U-J_H/8`, equivalently `mu<U-h/2-3J_H/8`.  Ferromagnetic
`-J_H S_1.S_2` then locks the pair to the symmetric spin-one triplet.  On the
diagonal coherent locus, `|s;2>=|s> tensor |s>`, the Berry connection and
quantum metric double.  With `A=-i<s|ds>` and the microscopic kinetic sign
`+i hbar a^dagger dot(a)`, the declared aligned `-h n.S` action is `k=-2`:
the ket eigenline is `O(-2)`.  Its dual prequantum line is `O(2)`, has
`c1=+2`, and `dim H0(CP1,O(2))=3`; reversing orientation/coupling gives
`k=+2`.  Hence the derived invariant statement is `|k|=2`, while signed
chirality remains open.  The `27/27` audit obtains singlet gap `J_H`, Berry
residual `2.24e-16`, and final curvature-magnitude estimate
`2.000000501994128`.  At `U=4`, `mu=1.5`, `J_H=1`, `h=0.2`, the sufficient
condition margin is `2.025`; direct spinful occupancy enumeration retains the
unique `(1,1)` ground sector with interacting charge gap `1.85`, and the
analytic tail is coercive when `4U>J_H`.

This is not yet the UV origin of Route E.  The exactly-two orbital rule and
singleton/odd-sector exclusion are not derived; a singleton would leave a
`|k|=1` doublet.  The dimer also lacks an anomaly-consistent four-dimensional
embedding, a derived Route-E portal below all charge/singlet/orientation gaps,
an orientation/coupling principle selecting `k=+2` rather than the declared
aligned action's `k=-2`, and a proof that the triplet is chiral-family data.
The Mott number quoted by the audit starts from a bare onsite gap; exchange and
portal corrections must remain smaller than the retained interacting gap.
Fermion-determinant, mixed-WZW (`k=n_c B`), and AP-E4 tangent/Dirac
constructions remain alternatives.  AP-E4 is next: construct the physical
Dirac/fluctuation operator, its chirality, partner spectrum, and gap.

## 2026-07-15 AP-E3 UV/AP-E4 execution result

The old AP-E3 exactly-two and orientation blockers are now closed inside a
strictly declared constrained-cell model.  Bosonic Schwinger partons obey
canonical CCR, and two independent compact Gauss laws `G_r=N_r-1=0` make
`(N_1,N_2)=(1,1)` an exact physical-state rule, not a finite Mott penalty.
The negative electron magnetic moment fixes
`H_Z=+h n.(S_1+S_2)`, selects the anti-aligned triplet, and gives
`i hbar <Omega_-|d Omega_->=+2 hbar A_+`.  Its positive line is
`Q tensor Q=O(2)`.  The `26/26` verifier checks projector rank, negative
controls, the exact Hund-Zeeman spectrum, random-direction ground projectors,
Berry sign, quantum metric, Chern number, gauge anomalies, and WZW reduction.

An independent anomaly-consistent mixed-WZW candidate is now explicit:
`SU(2)_c x U(1)_g` with two vectorlike charge-one Dirac flavours and a
charge-one scalar doublet.  It has `kappa_L=-kappa_R=2` and a `B=+1` soliton
reduces to signed `k=+2`; the first local baryon dressing is
`qq(phi^dagger)^2`, an `O(2)` triplet.  This remains an intermediate UV
completion because strong-vacuum selection, diquark/Pauli--Guersey gaps,
compact-monopole/bordism consistency, soliton stability, and the abelian
Landau pole are unresolved.  The five-dimensional construction is
extension-independent on extendible sectors; non-extendible sectors still
need a differential-character/Cech-bordism definition.

AP-E4 is mathematically solved for the declared canonical Spin-c operator.
Horizontal linearization proves `T^(1,0)CP1=O(2)`.  At AP-E1 radius `R=1/2`,
`D_T^c=sqrt(2)(dbar_T+dbar_T^dagger)` has three positive zero modes, no
negative kernel, and
`lambda_(n,+/-)=+/-sqrt(n(n+3))/R`, multiplicity `2n+3` per sign.  The first
gap is `4` and every massive level has a partner.  The `22/22` verifier checks
this with exact cohomology and an independent finite-`SU(2)` Casimir
diagonalization, including the corrected `n=3` value `6 sqrt(2)` and the
distinction between eigenvalue sign and chirality.  Ordinary spin Dirac
twisted only by `T=O(2)` instead has two
zero modes and gap `2 sqrt(3)`: the canonical Spin-c half-canonical shift must
still be derived from a moduli-space or compactification UV theory.

Neither result identifies the triplet/zero modes with Route-E chiral
families.  A degree-one portal, physical fermionic origin, four-/six-
dimensional anomaly cancellation, and stability below all gaps remain open;
`physics_promotion_allowed=false`.

## 2026-07-16 AP-E3 global/nonlinear and AP-E4 SQM result

The requested charged-two-colour computation has been implemented as a
strictly labelled nonlinear low-energy proxy.  It minimizes a six-real-field
meson/charge-two-diquark potential on a declared `(e_g,v_g/Lambda_c)` matching
grid and solves the full nonlinear massive-Skyrme `B=1` hedgehog.  Its
deterministic card includes the vacuum Hessian and mass gaps, baryon-number
and Derrick residuals, grid/volume convergence, a generalized radial meson
Hessian, charged-diquark `l=0,1` blocks, and an unstable negative control.
This is useful necessary evidence, but not a four-dimensional lattice gauge-
theory calculation or a complete coupled 3D/quantum soliton Hessian.

The non-extendible mixed-WZW sector now has an intrinsic global definition.
For `X=S3 x S2`, the integral five-form is the curvature of a degree-five
differential character, evaluated directly on every mapped spacetime
four-cycle.  At fixed curvature/class the bosonic lift is unique because
`H^4(X;R/Z)=0`.  Spin refinement is not unique:
`Omega_4^Spin(X)=Z + Z2 + Z2`, and the two reduced signs require a microscopic
APS/Dai--Freed determinant.  The note gives both inverse-image bordism
invariants explicitly instead of hiding them in an extension choice.

The proposed pure-`SU(3)` all-scale embedding has also been decided.  The
faithful subgroup after adjoint breaking is `[SU(2) x U(1)]/Z2`, whose
representations satisfy `2j+q=0 mod 2`.  Therefore the original charge-one
colour singlet and its two-scalar `O(2)` baryon dressing cannot descend from
`SU(3)`.  A concrete charge-two variant with
`bar(3)->2_(-1)+1_(+2)` and explicit adjoint mass splitting does decouple its
coloured scalars and the `1_(-2)` partners of vectorlike fundamentals at tree
level.  Its baryon is dressed by one scalar and is an `O(1)` doublet, so it is
a controlled near-miss rather than Route-E closure.  Moreover, its covering-
`U(1)` mixed-level normalization is not yet a faithful-`U(2)` all-bundle
theorem; the correlated centre-flux/determinant-line match is an explicit
false gate.

Finally, the preferred AP-E4 route now independently derives a physical
tangent fermion: a half-BPS non-Abelian vortex has `CP1` orientational moduli
and BPS collective-coordinate supersymmetry, giving `N=2` moduli-space SQM.
The declared canonical Spin-c re-quantization has one untwisted Dolbeault
ground state, whereas the source-selected half-form ordering has none.  A separate `O(2)`
coefficient line reproduces the three positive-chirality modes and paired
gap.  AP-E3's WZW level can supply that line only after a same-mother-model
theorem identifies the charged-two-colour `B=1` moduli, spatially transgresses
the degree-five WZW character, and fixes the signed line pullback.
The independent BPS vortex does not prove this composition.  Product
compactification remains the anomaly-checked backup, and the degree-one
Route-E portal remains deliberately deferred.  All three cards set
`physics_promotion_allowed=false`.

The vacuum-line and CPT choices are additional gates.  At fixed canonical
polarization, the AP-E3 antibaryon sign `k=-2` gives one negative mode, not
the conjugate of the three positive modes.  A physical completion must
derive the anti-canonical polarization map (or its effective `O(-4)` form)
rather than infer it from the sign of the WZW integer.

## 2026-07-14 execution status

- `route_E/code_dyn/run_route_e_dynamics.py` executes the 19 recovered lanes
  plus DYN-5V/DYN-7F guards in an isolated snapshot and records a digest-rich
  manifest.  Mechanical `all_pass` and physics status are separate fields.
- An isolated DYN-0 -> DYN-4a run succeeds.  DYN-4a now reports a converged
  connected-local profile (`23/23`) with
  `chi2: 47.455523 -> 39.193276`; it explicitly does not claim global-basin
  or disconnected-component completeness.
- DYN-5V passes `9/9` matching checks yet deliberately classifies DYN-5 as
  `invalid_pending_rederivation`: `delta Z_tree=|zeta|/3=0.0434773011`, the
  full six-by-six/Schur Weinberg tensors agree below `1e-12`, the displayed
  quadratic action has no cubic loop interaction, and `X L H_u` is allowed.
- DYN-7F passes `7/7` diagnostics yet classifies DYN-7 as
  `blocked_missing_branch_thermal_inputs`.  At
  `M1=2.38045e10 GeV` the non-SUSY lane is tau-resolved/two-flavor; the old
  Davidson--Ibarra bound also contained an extra square root (factor `3.91`).
- The required-card 21-node dry-run has no preflight errors.  DYN-9b-2's
  hypercharge normalization is repaired (`y_t(M_X)=0.44116481`), DYN-9b-3
  reproduces the corrected SM DI ceiling `2.307794855e-6`, and DYN-8 passes
  `30/30` disclosure checks while explicitly setting
  `physics_promotion_allowed=false`.
- `audit_blocker_promotion_gate.py` passes `18/18`: it verifies the
  one-dimensional Abelian H3 counterexample and exact additive-charge no-go,
  blocks RE-SC4's DYN-5-derived numerical gap, separates fixed-kernel from
  optional top-like DYN-9b-2 scaling, checks complex-rescaling covariance,
  exposes RE-SC5's `M_I~=101 GeV` toy lower edge, and preserves the remaining
  fit and flavored-kinetics blockers.

Artifacts:

- `THEORY_AUDIT.md`: Route-by-route findings, derivations, and recommendations.
- `CODE_DYN_REAUDIT.md`: restored-source replay, scientific QA, and revised
  evidence classification for Route E.
- `ROADMAP.md`: ordered work packages with pass/fail gates.
- `LITERATURE_2023_2026.md`: recent primary-literature matrix.
- `code/verify_route_f_diagnostics.py`: a small diagnostic arithmetic card.
- `code/audit_evidence_manifest.py`: checks whether the evidence cited for the
  Route-E dynamics claims exists in this checkout.
- `code/audit_blocker_promotion_gate.py`: recomputes the P0 numerical repairs
  and proves that mechanically green incomplete ledgers remain non-promotable.
- `code/verify_another_physics_route_e_bridge.py`: checks the exact
  discriminant/Killing identity, Hopf moment map, zeta covariance, and the
  Q-ball benchmark without promoting the candidate bridge.
- `code/verify_ap_e1_projective_doublet.py`: checks the Hopf projector and
  metric, Chern number, Routh/Legendre transforms, Schwinger triplet, monopole
  spectrum for both flux signs, Branch-B symplectic sign, the exceptional
  zero-level orbit, fixed-charge orientation energy, Q-ball charge mismatch,
  Derrick scaling, and critical TeX/source hashes (`30/30`).
- `code/verify_ap_e2_discriminant_regression.py`: exact discriminant, Killing,
  transvectant, basis, orbit-boundary, covariance, and negative-control
  regression (`30/30`, non-promoting).
- `code/verify_ap_e3_level_two_microscopic.py`: full spinful Mott occupancy,
  Hund triplet, Veronese/Berry/metric doubling, `|k|=2`, ket/prequantum-line
  sign separation, representation, scale hierarchy, coercive-tail, and failure
  controls (`27/27`; candidate-level only).
- `code/verify_ap_e3_exact_two_mixed_wzw.py`: exact separate-Gauss projector,
  physical `k=+2` orientation, anomaly ledger, mixed-WZW quantization and
  worldline reduction, dressing/evenness, running, and negative controls
  (`26/26`; all-scale UV and Route-E portal remain open).
- `code/verify_ap_e4_tangent_dirac_spectrum.py`: tangent transition/flux,
  canonical Spin-c cohomology, exact and matrix-derived complete spectrum,
  chirality/partner tower, gap, ordinary-spin control, and fail-closed
  physical-origin gates (`22/22`).
- `code/verify_ap_e6_sp4_eta_bordism_threshold.py`: exact low-degree
  `BSp(2)` AHSS and representation character, partial regulator controls,
  target-sign nonselection, and fail-closed heavy-threshold audit (`39/39`).
- `code/scan_ap_e6_relaxed_b1_multigrid_hessian.py`: canonical coupled
  hedgehog BVP, exact sparse `n+s` second variation, off-shell multigrid/
  multivolume scan, and aggregate-Hessian negative gate (`20/20`).
- `code/verify_ap_e6_same_soliton_yukawa_callias.py`: frozen-profile
  provenance, actual `SO(3)` orbit, Callias/two-band no-go, algebraic CPT,
  and typed-but-undescended WZW character (`25/25`).
- `code/verify_ap_e6_execution_frontier.py`: independently checks fresh
  manifests, the canonical same-background checksum, every lane conjunction,
  and the no-portal rule (`15/15`).
- `code/verify_ap_e7_common_sp4_aps_pv_regulator.py`,
  `code/scan_ap_e7_discrete_rerelax_superhessian.py`, and
  `code/verify_ap_e7_so3_fr_rank3_family.py`: common APS/PV, finite-site
  topology/re-relaxation, and actual-`SO(3)`/rank-three audits (`42/42`,
  `16/16`, and `18/18`; all physical lane flags false).
- `code/verify_ap_e7_execution_frontier.py`: independent fresh-manifest,
  conjunctive-lane and alternative-family gate (`16/16`; no portal).
- `tex/another_physics_route_e_derivation_ledger.tex`: complete reversible
  derivation ledger for the Another-Physics / Route-E relationship.
- `tex/another_physics_route_e_bridge.bib`: local primary-source bibliography
  for the bridge ledger.
- `tex/ap_e1_projective_doublet_action.tex` and
  `tex/ap_e1_projective_doublet_action.bib`: complete AP-E1 derivation and
  primary-source bibliography.
- `tex/ap_e2_discriminant_regression.tex`: complete 12-page AP-E2
  first-principles theorem, convention ledger, and non-promotion handoff.
- `tex/ap_e3_level_two_microscopic_origin.tex` and
  `tex/ap_e3_level_two_microscopic.bib`: complete AP-E3 Mott--Hund,
  Veronese/Berry/Chern derivation, alternatives, failure gates, and
  primary-source bibliography.
- `tex/ap_e3_exact_two_mixed_wzw_uv.tex` and its bibliography: complete
  constrained-cell, Pauli-orientation, anomaly-consistent mixed-WZW, local
  dressing, and deep-UV boundary derivation.
- `tex/ap_e4_tangent_dirac_spectrum.tex` and its bibliography: complete
  tangent-valued geometry, canonical/ordinary Spin distinction, index,
  chirality, partner spectrum, gap, and anomaly-gate derivation.
- `output/pdf/ap_e3_exact_two_mixed_wzw_uv.pdf` and
  `output/pdf/ap_e4_tangent_dirac_spectrum.pdf`: clean-built, warning-free,
  page-by-page checked AP-E3 UV (12-page) and AP-E4 (9-page) notes.
- `output/another_physics_route_e_derivation_ledger.pdf`: clean-built,
  visually checked 41-page rendering of the derivation ledger; the same
  current build is mirrored under `output/pdf/`.
- `output/pdf/ap_e12_graph_density_regular_minimizer.pdf`: warning-free,
  page-by-page checked 6-page endpoint-density/regularity note.
- `tex/ap_e6_sp4_eta_bordism_threshold.tex`,
  `tex/ap_e6_relaxed_b1_multigrid_hessian.tex`, and
  `tex/ap_e6_same_soliton_yukawa_callias.tex`: complete AP-E6 formula,
  numerical, and rollback records with their local bibliographies.
- `output/pdf/ap_e6_sp4_eta_bordism_threshold.pdf`,
  `output/pdf/ap_e6_relaxed_b1_multigrid_hessian.pdf`, and
  `output/pdf/ap_e6_same_soliton_yukawa_callias.pdf`: warning-free,
  page-by-page checked 12-, 7-, and 14-page AP-E6 notes.
- `output/ap_e1_projective_doublet_action.pdf`: clean-built, visually checked
  15-page AP-E1 note.
- `output/pdf/ap_e2_discriminant_regression.pdf` and
  `output/pdf/ap_e3_level_two_microscopic_origin.pdf`: clean-built,
  page-by-page checked AP-E2 (12-page) and AP-E3 (13-page) notes.
- `route_E/code_dyn/dyn_claim_registry.json`: single physics-status registry.
- `route_E/code_dyn/run_route_e_dynamics.py`: fail-closed isolated DAG runner.
- `output/another_physics_route_e_bridge.{json,md}`: `27/27` bridge audit
  with explicit non-promotion status.
- `code/scan_ap_e3_charged_two_colour_proxy.py` and
  `output/ap_e3_charged_two_colour_proxy.{json,md}`: nonlinear EFT vacuum,
  soliton, Hessian, and convergence card (`18/18`; microscopic QC2D and full
  coupled 3D/quantum Hessian remain open).
- `code/verify_ap_e3_nonextendible_wzw_su3_uv.py` and
  `output/ap_e3_nonextendible_wzw_su3_uv.{json,md}`: differential-character,
  spin-bordism, representation no-go, and conditional `SU(3)` decoupling
  audit (`38/38`).
- `code/verify_ap_e4_moduli_space_sqm.py` and
  `output/ap_e4_moduli_space_sqm.{json,md}`: intrinsic tangent-fermion SQM,
  declared-vacuum-line spectrum, WZW-transgression/CPT gates, and backup
  compactification audit (`27/27`).
- `output/pdf/ap_e3_charged_two_colour_soliton_proxy.pdf`,
  `output/pdf/ap_e3_nonextendible_wzw_su3_uv_audit.pdf`, and
  `output/pdf/ap_e4_moduli_space_sqm.pdf`: 2026-07-16 derivation notes,
  clean-built and page-by-page checked at 9, 12, and 11 pages respectively.
- `output/ap_e1_projective_doublet.{json,md}`: `30/30` AP-E1 audit with the
  same explicit non-promotion status.
- `output/ap_e2_discriminant_regression.{json,md}`: `30/30` AP-E2 exact
  regression with `physics_promotion_allowed=false`.
- `output/ap_e3_level_two_microscopic.{json,md}`: `27/27` AP-E3 candidate
  audit with `ap_e3_physics_closed=false`.
- `output/ap_e3_exact_two_mixed_wzw.{json,md}`: `26/26` exact-cell and
  anomaly-consistent intermediate-UV audit with full physics promotion false.
- `output/ap_e4_tangent_dirac_spectrum.{json,md}`: `22/22` AP-E4 canonical
  operator audit with mathematical closure and physical closure false.
- `output/ap_e6_sp4_eta_bordism_threshold.{json,md}`,
  `output/ap_e6_relaxed_b1_multigrid_hessian.{json,md}`, and
  `output/ap_e6_same_soliton_yukawa_callias.{json,md}`: `39/39`, `20/20`,
  and `25/25` mechanical cards with every full-lane flag false.
- `output/ap_e6_execution_frontier.{json,md}`: `15/15` fresh-manifest and
  same-background integration card with portal authorization false.
- `tex/ap_e7_common_sp4_aps_pv_regulator.tex`,
  `tex/ap_e7_discrete_rerelax_superhessian.tex`, and
  `tex/ap_e7_so3_fr_rank3_family.tex`, with their local bibliographies:
  complete AP-E7 derivations and explicit proof/assumption boundaries.
- `output/pdf/ap_e7_common_sp4_aps_pv_regulator.pdf`,
  `output/pdf/ap_e7_discrete_rerelax_superhessian.pdf`, and
  `output/pdf/ap_e7_so3_fr_rank3_family.pdf`: warning-free, page-by-page
  checked 15-, 6-, and 8-page AP-E7 notes.
- `output/ap_e7_common_sp4_aps_pv_regulator.{json,md}`,
  `output/ap_e7_discrete_rerelax_superhessian.{json,md}`, and
  `output/ap_e7_so3_fr_rank3_family.{json,md}`: `42/42`, `16/16`, and
  `18/18` fail-closed evidence cards.
- `output/ap_e7_execution_frontier.{json,md}`: `16/16` integrated gate with
  every complete-lane, portal, and promotion boolean false.
- `output/`: generated JSON, Markdown, and compiled-paper artifacts.

The recommended active physics branch is a non-supersymmetric intermediate-
scale `Spin(10)` model, with the kinematic Route-A/E core treated as a
conditional organizing principle.  The exact Higgs content remains a decision
gate: Route F compares a `210_H` continuity branch with the better-established
`54_H` branch before freezing the action.  Route B and Route D remain optional
UV mechanisms until their selection-rule and global-consistency gates pass.

Workspace provenance note (updated 2026-07-17): the initial Route-E/Route-F evidence
set was committed and pushed at `6662dd6`; the blocker corrections were
committed and pushed at `e7ba020`; the reversible Another-Physics bridge ledger
was committed and pushed at `6c61d56`.  AP-E1 was then executed as a separate
non-promoting derivation and committed/pushed at `cce4913`.  AP-E2 is now
frozen as an exact non-promoting regression.  AP-E3's constrained-cell and
signed-orientation subgates and AP-E4's canonical spectral mathematics are now
complete.  AP-E5, AP-E6, and the present AP-E7 alternatives have since been
executed at their declared mechanical scopes, without closing a physical
pre-portal lane.  The mixed-WZW strong/all-scale gates, AP-E4 fermionic
Spin-c origin/anomalies, the degree-one Route-E portal, and the remaining
downstream phenomenology gates must still pass before this bridge can alter
any Route-E physics claim.
