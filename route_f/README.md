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
2. Route E's Cartan-selection proof does not exclude the genus-one branch as
   written: a one-dimensional abelian Lie algebra admits a non-degenerate
   symmetric ad-invariant bilinear form even though its Killing form vanishes.
3. Route E's 19 DYN source scripts have been recovered under
   `route_E/code_dyn/`.  Canonical paths and a fail-closed isolated DAG now
   exist, but Route-D string-pricing-card `RE-SC3/4/5` evidence is still
   missing and prevents a complete run.  Several green arithmetic gates also
   fail their physical interpretation.  Route C and D likewise still contain
   classification/arithmetic gates rather than complete amplitudes or a global
   compactification.  The evidence card checks a 29-item dynamics P0 subset
   plus 4 core artifacts, not every path mentioned anywhere in the roadmap.

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

Artifacts:

- `THEORY_AUDIT.md`: Route-by-route findings, derivations, and recommendations.
- `CODE_DYN_REAUDIT.md`: restored-source replay, scientific QA, and revised
  evidence classification for Route E.
- `ROADMAP.md`: ordered work packages with pass/fail gates.
- `LITERATURE_2023_2026.md`: recent primary-literature matrix.
- `code/verify_route_f_diagnostics.py`: a small diagnostic arithmetic card.
- `code/audit_evidence_manifest.py`: checks whether the evidence cited for the
  Route-E dynamics claims exists in this checkout.
- `route_E/code_dyn/dyn_claim_registry.json`: single physics-status registry.
- `route_E/code_dyn/run_route_e_dynamics.py`: fail-closed isolated DAG runner.
- `output/`: generated JSON and Markdown ledgers.

The recommended active physics branch is a non-supersymmetric intermediate-
scale `Spin(10)` model, with the kinematic Route-A/E core treated as a
conditional organizing principle.  The exact Higgs content remains a decision
gate: Route F compares a `210_H` continuity branch with the better-established
`54_H` branch before freezing the action.  Route B and Route D remain optional
UV mechanisms until their selection-rule and global-consistency gates pass.

Workspace provenance note (2026-07-14): Route-E/Route-F sources and the compact
evidence set are being prepared for intentional staging.  Until staging is
completed, a clean clone does not contain this integration lane; consult Git
status rather than assuming the local workspace equals `HEAD`.
