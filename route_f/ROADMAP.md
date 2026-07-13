# Route F Roadmap: One Action, One Evidence Chain

Created: 2026-07-13
Last updated: 2026-07-14

Status values: `open`, `in-progress`, `done`, `failed`, `permanently-open`.
All items start `open` unless marked otherwise.

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
```

## P0: blockers

### F0-A Evidence recovery and single status registry

Status: `in-progress`.

Re-audit update (2026-07-13): all 19 DYN source scripts are present at
`route_E/code_dyn/`.  In a temporary historical-layout replay, 16 scripts
complete with `219/219` internal checks.  DYN-9b-2, DYN-9b-3, and DYN-8 remain
blocked by the missing Route-E string-card `RE-SC3/4/5` evidence chain
(historical D3--D5 filenames under `route_d/`).  The current source
location itself is not runnable as one DAG, no DYN ledgers are present in the
workspace, and scientific-status defects are catalogued in
`route_f/CODE_DYN_REAUDIT.md`.

Execution update (2026-07-14):

- **F0-A1 canonical-path/cache subgate: done.**  All 19 recovered scripts use
  `route_E/code_dyn/route_e_paths.py`; lowercase paths work from both the
  repository root and `route_E/code_dyn/`.  DYN-9b-1c/1d caches are isolated
  under `ROUTE_E_CACHE_DIR` and keyed by source plus NumPy version.
- **F0-A2 runner subgate: implemented, full closure still in-progress.**
  `run_route_e_dynamics.py` snapshots minimal inputs into an isolated
  workspace, executes a fail-closed DAG, records Git/Python/NumPy/BLAS/SciPy
  provenance and SHA-256 digests, and separates mechanical `all_pass` from
  physics status.  An isolated DYN-0 -> DYN-4a replay passes; a full optional
  dry-run blocks only declared missing dependencies/descendants.  The
  expensive full replay and clean-clone test remain open.
- **F0-A3 registry subgate: done.**  `dyn_claim_registry.json` is the canonical
  status source.  DYN-4a's interval defect is repaired but its fit remains
  conditional; DYN-5 is `invalid_pending_rederivation`; DYN-7 is
  `blocked_missing_branch_thermal_inputs`; DYN-9b-3 and DYN-8 remain blocked.
- **Scientific guard evidence:** DYN-4a now uses a two-optimizer candidate fit,
  a converged needle-basin stencil, and connected-local nuisance-minimized
  profiles (`23/23`; no global-basin completeness claim); DYN-5V verifies
  tree Kahler/full-6x6 matching and the selection-rule counterexample (`9/9`);
  DYN-7F classifies the tau-resolved regime and repairs the
  Davidson--Ibarra double-square-root bug (`7/7`).  Passing guard arithmetic
  intentionally does not close the DYN-5 or DYN-7 physics blockers.

Tasks:

- **F0-A1 path/provenance gate:** move or wrap the restored scripts behind one
  repository-root resolver; use lowercase `route_e`; add `--repo-root`,
  `--output-dir`, seed, and cache arguments; remove dependence on manual
  symlinks; track the complete input/source tree in Git; and record commit plus
  dirty-tree state.
- **F0-A2 clean-run gate:** build a fail-closed DAG runner for all required
  recovered lanes (`DYN-0--5`, `DYN-7--9b`, explicitly excluding optional
  `DYN-6`), generate every required JSON/MD ledger from a clean clone, use
  content-addressed and
  atomically written caches, and store Python/NumPy/BLAS plus full
  input/source/output digests.  Decide whether `RE-SC3/4/5` are required
  dependencies or optional string cards; DYN-8 must not load an optional card
  unconditionally.  DYN-6 remains optional/open and is not missing evidence.
- **F0-A3 scientific-status gate:** create one claim registry using at least
  `proved`, `conditional`, `numerical`, `preliminary`,
  `invalid_pending_rederivation`, `blocked_missing_evidence`, `open`, and
  `no-go`.  At present DYN-1/2 are conditional SUSY-benchmark evidence; DYN-3
  is preliminary until physical flavor rotations/dressing; DYN-4 is inverse
  reconstruction rather than a global fit; DYN-5 is invalid pending an
  interacting messenger rederivation; DYN-7/9b-3 require flavored
  leptogenesis; and DYN-9/9b remain preliminary.
- **F0-A4 string-existence gate:** recover or reimplement Route-E string cards
  `RE-SC3/4/5` without confusing them with Route-D `RD-D3H/RD-D4M/RD-D5P`.
  Until their scripts and ledgers exist, K8--K10 are
  `blocked_missing_evidence`; after recovery they remain
  `unpromoted_pricing_only` until F10 passes.
- generate README/ROADMAP/paper status tables from the registry and remove
  stale contradictions, including the DYN-8 210-vs-45 branch-map conflict.

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

Status: `open`.

Tasks:

- define `canonical` precisely;
- either require the Killing form or prove functorial `Aut(g)` naturality;
- handle the one-dimensional abelian counterexample explicitly;
- weaken "two distinct centers" to a degree-two divisor unless a semisimple
  orbit-selection principle is supplied;
- keep the minimal-dimension assumption in the `Spin(10):16` theorem.

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

Status: `open`.

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

Run the 21-node full DAG (19 recovered lanes plus DYN-5V/DYN-7F guards) from a
clean clone after either recovering `RE-SC3/4/5` or formally excluding their
consumers from the required target.  In parallel, construct an explicit
interacting messenger/selection sector for DYN-5 and supply branch-local
thermal inputs for a two-flavor or density-matrix DYN-7/9b-3 calculation.
F0-B may proceed independently.  Do not spend more compute on new flavor or
proton scans until F1 fixes the action and F2 exports the actual spectrum.
