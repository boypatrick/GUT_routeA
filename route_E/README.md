# Route E Workspace: First-Principles Reconstruction

Route E is a separate sandbox that rewrites the Route-A note as a derivation
from an explicit four-principle set (P0 consistency, P1 face projection,
P2 index protection with a self-carried tangent family bundle, P3 canonical
contact), plus a separately disclosed H3+ Killing-contact selection axiom.
It asks how much of the Route-A skeleton and of the Route-A
conditional-input list can be turned into theorems of that principle set, and
it proves that the residue (notably the value of `zeta`) cannot.

Route E remains a candidate restatement.  The main paper's status ledger is
synchronized with its blocker findings, but no conditional Route-E claim is
thereby promoted into the Route-A theorem core.

## Current Status

### Superseding promotion correction (2026-07-14)

The theorem-core audits and the dynamics program have different status.  The
machine-readable authority for dynamics is `code_dyn/dyn_claim_registry.json`;
a script exit code or ledger `all_pass=true` is mechanical evidence only.

- RE-SC3/4/5 have been restored, so their file-existence blocker is closed,
  but all three remain `unpromoted_pricing_only`.
- DYN-5 is `invalid_pending_rederivation`: the displayed `XN` term is
  bilinear, its asserted Yukawa loop is absent, tree matching gives
  `delta Z=|zeta|/3`, and the displayed charges allow `X L H_u`.
- DYN-9b-2 is a preliminary fixed-kernel scaling study, not the missing
  global non-SUSY flavor fit.  The displayed SO(10) matrix relations do not
  force a top-like neutrino Dirac coupling; exact zeta invariance is limited
  to uniform positive-real rescaling, while complex rescaling rotates its
  phase.  DYN-7 and DYN-9b-3 remain blocked on
  tau-resolved flavored kinetics and branch-local thermal inputs.  DYN-8 is a
  fail-closed status collector with `physics_promotion_allowed=false`, not a
  closure theorem.
- The genus-selection step is repaired rather than promoted as originally written:
  a one-dimensional
  abelian Lie algebra admits a non-degenerate symmetric ad-invariant form even
  though its Killing form vanishes.  The current argument proves the
  Riemann--Roch upper bound `N_fam<=3`; selecting `g=0` and `N_fam=3`
  is now explicitly conditional on the H3+ semisimplicity/non-degenerate-
  Killing axiom.  A physical origin for H3+ remains open.

The numerical blocker gate is `route_f/code/audit_blocker_promotion_gate.py`.

- Status: complete draft with two machine audits (first-principles 44/44;
  dependency-closure 24/24) and a dependency roadmap (`ROADMAP.md`, items
  R1-R15).
- Dependency closure (2026-07-04): the four former classical inputs are now
  machine-verified (w0 involutions for every simple algebra of rank 4-15
  plus exceptional; exact SU(15)/SU(16) cubic traces; the exactly vanishing
  cubic tensor of the Fock-constructed chiral 16; exceptional minima
  27/56/26/248).  Theorem 5's witness family is internal (the Route-B EFT at
  variable coupling), so Route E has NO dependency on Route-D material;
  Route D is cited only as optional UV interpretation.  Every remaining
  external ingredient is tagged in the paper's provenance ledger as
  textbook-cited, retained inputs (including H3+), or benchmark corroboration.
- Claim boundary: Route E strengthens the *derivation* of the Route-A
  skeleton; it does not add phenomenology.  The separate companion dynamics
  lanes have been mechanically replayed but remain unpromoted/open, and Route E
  does not claim a PSLT-only unconditional GUT proof.
  Its Theorem 5 proves that last claim impossible within the stated model
  class, which is the Route-A boxed corollary upgraded to a theorem.

Checked so far (see `output/route_e_first_principles.{json,md}`):

- Enumerated six-face uniqueness: among rank >= 4 simple compact groups with
  complex irreps, the only complex irreps of dimension 16 are the SU(16)
  fundamental (gauge-anomalous, excluded) and the Spin(10) half-spinors;
  the only complex irreps of dimension 15 (SU(5) Sym^2 5, SU(6) Lambda^2 6,
  SU(15) fundamental) all fail content or anomaly, so the sixteenth state
  (`nu^c`) is forced, not assumed.
- Genus ladder `h^0(Sigma_g, T) = 3, 1, 0`: three families is an a-priori
  maximum.  The former claim that the Cartan criterion alone selects `g=0`
  is blocked by the one-dimensional abelian invariant-form counterexample;
  `N=3`, `T~O(2)`, and the two-center divisor are conditional on the added
  genus-zero/semisimplicity gate described above.
- On the H3+-selected `g=0` branch, the invariant contact is unique
  (multiplicity one) and equals the Killing form of the family algebra in the
  spherical basis: `B = 2 sqrt(3) K_tr`.
- Route-B form-uniqueness: every term of the messenger superpotential is
  fixed by invariance; integrating out gives `Delta M_R = lambda^2 K_tr`.
- Boundary theorem: a two-model argument shows no zeta value is a consequence
  of the principles; the Route-A negative scans (Z_{N<=6} miss, first window
  hit at N = 178, Audit 0.5 no-hit) are finite-class corroborations, and
  17/178, 70/733, 87/911 are ordinary continued-fraction convergents.

## Directory Layout

- `code/`: the Route-E verification audit (plain python3 + numpy, no sympy).
- `code_dyn/`: the single canonical DYN implementations, validity guards,
  claim registry, and isolated fail-closed runner.  Historical root `code/`
  DYN paths are compatibility delegates, not independent implementations.
- `output/`: generated JSON + Markdown ledgers with negative-boundary flags.
- `tex/`: the rewritten paper draft `route_e_first_principles.tex`.

## Promotion Criteria

Identical to Route D:

1. A precise theorem, proposition, or falsifiable audit question.
2. A list of assumptions separated from derived statements.
3. A reproducible script or hand-checkable derivation.
4. A statement of how it affects Route A/B/C/D without overclaiming.

## Default Paper Use

Until promoted, Route E enters the main paper at most as a cited companion
draft.  In particular the main paper's Lean-Note Status Ledger stays as is;
the Route-E ledger relocations (conditional inputs -> theorems, "not claimed"
-> "provably impossible in the stated class") take effect only if promotion
is accepted.
