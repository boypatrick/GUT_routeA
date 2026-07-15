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

The independent bridge verifier passes `22/22` mathematical/numerical checks
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
  visually checked 17-page rendering of the derivation ledger.
- `output/ap_e1_projective_doublet_action.pdf`: clean-built, visually checked
  15-page AP-E1 note.
- `output/pdf/ap_e2_discriminant_regression.pdf` and
  `output/pdf/ap_e3_level_two_microscopic_origin.pdf`: clean-built,
  page-by-page checked AP-E2 (12-page) and AP-E3 (13-page) notes.
- `route_E/code_dyn/dyn_claim_registry.json`: single physics-status registry.
- `route_E/code_dyn/run_route_e_dynamics.py`: fail-closed isolated DAG runner.
- `output/another_physics_route_e_bridge.{json,md}`: `24/24` bridge audit
  with explicit non-promotion status.
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
- `output/`: generated JSON, Markdown, and compiled-paper artifacts.

The recommended active physics branch is a non-supersymmetric intermediate-
scale `Spin(10)` model, with the kinematic Route-A/E core treated as a
conditional organizing principle.  The exact Higgs content remains a decision
gate: Route F compares a `210_H` continuity branch with the better-established
`54_H` branch before freezing the action.  Route B and Route D remain optional
UV mechanisms until their selection-rule and global-consistency gates pass.

Workspace provenance note (updated 2026-07-15): the initial Route-E/Route-F evidence
set was committed and pushed at `6662dd6`; the blocker corrections were
committed and pushed at `e7ba020`; the reversible Another-Physics bridge ledger
was committed and pushed at `6c61d56`.  AP-E1 was then executed as a separate
non-promoting derivation and committed/pushed at `cce4913`.  AP-E2 is now
frozen as an exact non-promoting regression.  AP-E3's constrained-cell and
signed-orientation subgates and AP-E4's canonical spectral mathematics are now
complete.  The mixed-WZW strong/all-scale gates, AP-E4 fermionic Spin-c
origin/anomalies, the degree-one Route-E portal, and AP-E5--AP-E9 must still
pass before this bridge can alter any Route-E physics claim.
