# Route E Roadmap: Dependency Closure

Goal (2026-07-04): remove or fully explain every dependency of
`tex/route_e_first_principles.tex` on Route A/B/D material and on classical
citations; keep the proved-vs-retained boundary explicit.  Items are ordered;
each is `open` / `done` / `permanently-open (by theorem)` / `out-of-scope`.

Categories: **[remove]** = prove in-repo, dependency eliminated;
**[explain]** = retained input, justified and tagged, not eliminable by
principles of the P0-P3 kind; **[impossible]** = provably not closable;
**[scope]** = inherited Route-A obligations, untouched by Route E.

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
  NuFit-6.0 exactly and the posterior of the covariant decomposition is
  computed.  Refreshed central card: zeta' = 0.105557+0.071595i
  (arg moved 75 old loose windows, |zeta| moved 689 -- the printed anchor
  digits are decisively superseded by DATA, exactly as the boundary
  theorem requires), contact fraction' = 0.1275, M_*' = 4.09e15 GeV,
  I' = 0.007790+0.017991i, J' = -9.143e-5+5.142e-4i; Audit-0.5 re-test
  still no-hit (miss 0.566 rad).  UNCONDITIONAL posterior (NuFit-6.0
  gaussians x uniform Majorana phases x log-uniform m1 in [1e-4, 3e-2] eV,
  4000 draws): |zeta| = 0.144 [0.113, 0.184] (68%), arg zeta = 0.605
  [0.465, 0.749] rad (circular).  KEY ROBUSTNESS RESULT: contact
  essentiality survives the refresh and full marginalization --
  P(contact fraction > 0.01) = 1.000, P(> 0.05) = 0.988; the
  Veronese-only branch stays excluded posterior-wide.  Nuisance
  attribution: the alpha_21 Majorana phase dominates the width
  (corr +0.53), m1 +0.19.  DYN-7 feed: log10 M_1 = 10.44 [10.39, 10.52],
  P(M_1 > 1e9 GeV) = 1.000 (Davidson-Ibarra viable).  Old-card I/J
  regress to the Audit-0 paper digits.
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
  (conditional-input change); the 210-compatible PS chain preserves the
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

Suggested sequencing: DYN-0 -> (DYN-1 || DYN-4) -> DYN-2 -> DYN-3 ->
DYN-5 -> DYN-8, with DYN-6/DYN-7 optional lanes after DYN-4.  Statuses
start `open`; update per session.

## F. String-conditional program (planned 2026-07-05; items D3-D5,
## then DYN-8 close-out and the staged DYN-9b)

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
  Depends: DYN-9 (patched ledger).  Status: open.

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
  Depends: D3, DYN-5.  Status: open.

- **D5 [high-scale SUSY-breaking bridge scan] (1-2 sessions).**
  Goal: the interpolating family between the EXCLUDED SUSY slice and
  the DYN-9 non-SUSY chains: SUSY broken at M_SS (non-SUSY below, SUSY
  above) -- string-natural (UV SUSY), and it preserves the
  Route-B/DYN-5 holomorphy machinery at the matching scale iff
  M_SS < M_*.
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
  Depends: DYN-9, DYN-3, DYN-5.  Status: open.

- **DYN-8 close-out** as specified in section E (branch map + branch-
  tagged kill bank incl. the D3-D5 outputs).  Status: open.

- **DYN-9b [non-SUSY core re-derivation -- REQUIRED for the alive
  branch] (staged, 3+ sessions).**
  9b-1: non-SUSY scalar potential and spectrum of the 210+126bar(+10)
  system (superpotential-level blocks carry over from DYN-1, the
  potential does NOT); redo the DYN-9 thresholds with the real
  spectrum replacing ESH.
  9b-2: non-SUSY flavor refit at v_R ~ M_I: a perturbative M_R tower,
  zeta/M_* RE-EXTRACTION on the alive branch (branch-local zeta is
  consistent with the boundary theorem), contact-essentiality retest.
  9b-3: non-SUSY leptogenesis re-run on the 9b-2 tower (non-SUSY loop
  function; the gravitino constraint dissolves).
  9b-4: DYN-8 ledger refresh with the 9b numbers.
  Depends: DYN-9 (done), DYN-8 (structure); D3/D4 optionally as
  source alternatives.  Status: open.

Backlog (optional, unchanged): DYN-6 (axiom-pricing generalization),
DYN-4c (kernel-level Dirac refit -- the DYN-7 escape), DYN-3b
(physical-basis rotations, explicit dressing).
