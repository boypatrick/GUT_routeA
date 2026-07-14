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
- `tex/another_physics_route_e_derivation_ledger.tex`: complete reversible
  derivation ledger for the Another-Physics / Route-E relationship.
- `tex/another_physics_route_e_bridge.bib`: local primary-source bibliography
  for the bridge ledger.
- `output/another_physics_route_e_derivation_ledger.pdf`: clean-built,
  visually checked 17-page rendering of the derivation ledger.
- `route_E/code_dyn/dyn_claim_registry.json`: single physics-status registry.
- `route_E/code_dyn/run_route_e_dynamics.py`: fail-closed isolated DAG runner.
- `output/another_physics_route_e_bridge.{json,md}`: `22/22` bridge audit
  with explicit non-promotion status.
- `output/`: generated JSON, Markdown, and compiled-paper artifacts.

The recommended active physics branch is a non-supersymmetric intermediate-
scale `Spin(10)` model, with the kinematic Route-A/E core treated as a
conditional organizing principle.  The exact Higgs content remains a decision
gate: Route F compares a `210_H` continuity branch with the better-established
`54_H` branch before freezing the action.  Route B and Route D remain optional
UV mechanisms until their selection-rule and global-consistency gates pass.

Workspace provenance note (2026-07-14): the initial Route-E/Route-F evidence
set was committed and pushed at `6662dd6`; the blocker corrections were
committed and pushed at `e7ba020`.  The subsequent Another-Physics bridge
ledger is explicitly non-promoting and must pass AP-E1--AP-E9 before it can
alter any Route-E physics claim.
