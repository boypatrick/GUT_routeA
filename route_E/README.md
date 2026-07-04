# Route E Workspace: First-Principles Reconstruction

Route E is a separate sandbox that rewrites the Route-A note as a derivation
from an explicit four-principle set (P0 consistency, P1 face projection,
P2 index protection with a self-carried tangent family bundle, P3 canonical
contact).  It asks how much of the Route-A skeleton and of the Route-A
conditional-input list can be turned into theorems of that principle set, and
it proves that the residue (notably the value of `zeta`) cannot.

Route E does not modify `paper/gut_framework.tex`.  It is a candidate
restatement; promotion would move items between the Route-A ledger
categories, so it must clear the same bar as Route D before any of it is
cited by the main paper.

## Current Status

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
  textbook-cited, retained inputs (i)-(vi), or benchmark corroboration.
- Claim boundary: Route E strengthens the *derivation* of the Route-A
  skeleton; it does not add phenomenology, does not replay the deferred
  companion audits, and does not claim a PSLT-only unconditional GUT proof.
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
  maximum; the Cartan criterion (Killing form nondegenerate iff semisimple)
  selects g = 0, hence N = 3, `T ~ O(2)`, and the two-center divisor as the
  Poincare-Hopf zero divisor of a Cartan flow (deg = chi(S^2) = 2).
- The invariant contact is unique (multiplicity one) and equals the Killing
  form of the family algebra in the spherical basis: `B = 2 sqrt(3) K_tr`.
- Route-B form-uniqueness: every term of the messenger superpotential is
  fixed by invariance; integrating out gives `Delta M_R = lambda^2 K_tr`.
- Boundary theorem: a two-model argument shows no zeta value is a consequence
  of the principles; the Route-A negative scans (Z_{N<=6} miss, first window
  hit at N = 178, Audit 0.5 no-hit) are finite-class corroborations, and
  17/178, 70/733, 87/911 are ordinary continued-fraction convergents.

## Directory Layout

- `code/`: the Route-E verification audit (plain python3 + numpy, no sympy).
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
