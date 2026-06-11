# Route D Workspace: String-Constrained Lift

Route D is an intentionally separate sandbox for a string-derived or
string-constrained lift of the Route-A framework.  It asks whether the current
conditional inputs can be replaced by standard string compactification data:
singular fibers, matter curves, flux indices, Euclidean-brane instantons,
Green-Schwarz selection rules, Wilson lines, and modular/geometric functions.

Route D is not part of the Route-A theorem core, the optional Route-B
hidden-zeta mechanism, or the active Route-C amplitude-bootstrap audit until a
precise statement is verified.

## Current Status

- Status: scoped as a string-constrained lift, with D1-F local placement
  checked.
- Claim boundary: no Route-D result should be cited as evidence in the paper
  until it has a mathematical statement, a verification plan, and a roadmap
  entry.
- Default rule: keep Route D independent from `paper/gut_framework.tex` unless
  the route produces a checked result worth promoting.

Checked so far:

- D1-F local F-theory placement audit:
  `E6 -> SO(10) x U(1)` gives
  `78 = 45_0 + 1_0 + 16_-3 + bar16_+3`, and
  `K_{P1}^{1/2} tensor O(3) = O(2)` gives `h0=3`, `h1=0`.
  The audit is local only; it does not construct a global compactification or
  verify the massless-hypercharge flux condition.
- D2-E3 E3-instanton/Green-Schwarz audit:
  the benchmark coefficient can be interpreted as
  `zeta = A exp(2 pi i T_Gamma)` with hidden phase/radial data supplied by an
  axion-volume modulus.  A toy GS charge check shows how an instanton factor
  can neutralize a perturbatively forbidden Majorana insertion.  The audit does
  not construct an E3 divisor, count charged zero modes, or claim a ten-digit
  prediction of `zeta`.

## Candidate String Grafts

1. F-theory local GUT placement:
   `D5` singularity on a divisor gives `Spin(10)`, and codimension-two
   enhancement `D5 -> E6` gives the `16` matter curve.
2. Matter-curve index:
   on `Sigma = P1`, zero modes of `K_Sigma^(1/2) tensor L_flux` reproduce
   `O(2)` when `deg L_flux = 3`.
3. E3-instanton origin for the Majorana contact:
   `zeta ~ A exp(2 pi i tau_Gamma)` with phase/radial data supplied by the
   axion-volume pair of a Euclidean brane cycle.
4. Heterotic comparison branch:
   standard embedding and `SU(4)` bundle deformations can explain why tangent
   or tangent-related bundles appear, and why three generations are topological.
5. Modular/flavor check:
   compare the benchmark elliptic data `(I,J)` or `j(E_M)` to modular functions
   of a compactification modulus.

## Directory Layout

- `code/`: scripts or symbolic checks for Route D.
- `output/`: generated data, tables, and machine-readable ledgers.
- `tex/`: Route-D notes, derivations, or candidate appendix drafts.

## Promotion Criteria

A Route-D result can be considered for the main paper only after it has:

1. A precise theorem, proposition, or falsifiable audit question.
2. A list of assumptions separated from derived statements.
3. A reproducible script or hand-checkable derivation.
4. A statement of how it affects Route A, B, or C without overclaiming.

## Default Paper Use

Until the grafts are verified, Route D should enter at most as an Outlook or
candidate appendix.  It should not be used to claim a unique string vacuum or a
PSLT-only derivation of the GUT.
