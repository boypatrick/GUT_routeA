# Route D String-Graft Assessment

This note evaluates which string-theory grafts are worth developing and which
claims must be weakened before they can enter the paper.

## Executive View

The string comments are valuable, but they should not trigger a full rewrite of
the Route-A paper yet.  The strongest strategy is a targeted Route-D appendix or
companion note:

- keep Route A as the representation/index/contact skeleton,
- use Route D to explain possible UV origins of the conditional inputs,
- preserve the theorem boundary that no unique string vacuum is derived.

## Graft Ranking

### A. F-theory local placement: high value

Useful claims:

- `D5` singularity on a divisor realizes `Spin(10)`.
- Codimension-two enhancement `D5 -> E6` realizes a `16` matter curve through
  the adjoint branching
  `78 = 45_0 + 1_0 + 16_-3 + bar16_+3`.
- On a matter curve `Sigma = P1`, a flux bundle with degree 3 gives
  `K_Sigma^(1/2) tensor L_flux = O(-1) tensor O(3) = O(2)`.

Why it matters:

- This gives a string-local origin for the Route-A six-face object and for the
  `O(2)` family carrier.

Required caveats:

- chirality and multiplicity depend on flux data,
- hypercharge flux/Wilson-line breaking must avoid exotics and keep hypercharge
  massless,
- this is not yet a global compactification.

Recommended paper role:

- strong candidate for a short appendix if written conservatively.

### B. E3 instanton for `zeta`: highest value, but conditional

Useful claims:

- An E3 instanton naturally gives
  `zeta ~ A exp(2 pi i tau_Gamma)`.
- The hidden phase/radial sector becomes an axion-volume pair.
- Charged zero-mode counting plus Green-Schwarz selection rules can generate a
  Majorana-only operator while leaving perturbative Dirac/triplet operators
  untouched.

Why it matters:

- This directly targets Route B's biggest open datum: the origin of the complex
  coefficient multiplying `K_tr`.

Required weakening:

- Replace "zeta != 0 is inevitable" by "zeta != 0 is generic/expected if no
  exact gauge/discrete selection rule forbids the instanton and the required
  zero-mode structure exists."
- Do not claim ten-digit prediction.  Instanton sums and moduli stabilization
  can easily change the last many digits.

Recommended paper role:

- best candidate for the core of a Route-D string appendix or companion note.

### C. Heterotic standard embedding: medium value as comparison

Useful claims:

- The standard embedding `V=TX` is a textbook solution of the heterotic Bianchi
  identity and explains tangent-bundle language.
- An `SU(4)` bundle has `SO(10)` commutant in `E8`.
- Generation number is topological:
  `N_gen = 1/2 |int_X c3(V)|`.

Risks:

- The `TX + O` to stable `SU(4)` branch needs careful bundle-stability and
  Chern-class bookkeeping.
- This is a different global construction from the F-theory local picture.

Recommended paper role:

- comparison/outlook, not the first Route-D appendix unless F-theory fails.

### D. Modular flavor and elliptic `I,J`: high upside, speculative

Useful claim:

- Benchmark quartic data can define an elliptic curve; its `j` invariant may
  be compared to modular functions of a compactification modulus.

Why it matters:

- This is falsifiable with a short calculation.

Required caveats:

- No claim until a script computes `j(E_M)`, inverts to `tau`, and compares to
  a small library of modular forms/functions.

Recommended paper role:

- Route-D numerical audit.  Add only if it passes or as an explicit negative
  result in a companion note.

### E. Rubik covariance as string moduli/Wilson-line geometry: useful language

Useful claims:

- Wilson lines, fluxes, and bundle choices make the "face" idea geometric.
- Legal face changes are constrained by topology, anomaly matching, and
  Green-Schwarz mechanisms.

Risk:

- Too broad if written as a theorem.

Recommended paper role:

- one explanatory paragraph, not a theorem.

## Minimal Upgrade Path

1. Write a Route-D appendix draft with only:
   - F-theory `E6` enhancement branching,
   - `P1` matter-curve index giving `O(2)`,
   - E3 instanton form for `zeta`,
   - Green-Schwarz/zero-mode interpretation of Majorana-only selection.
2. Keep all global compactification and vacuum-selection claims conditional.
3. Do not insert modular flavor or heterotic claims into the main paper until
   their audits are complete.

## Bottom Line

Worth developing:

- F-theory local placement,
- E3 instanton origin of `zeta`,
- matter-curve index origin of `O(2)`.

Worth testing, but not yet promoting:

- modular `I,J`/`j(E_M)` identification,
- heterotic `SU(4)` comparison branch.

Worth avoiding as a strong claim:

- "string theory makes the whole GUT unconditional",
- "zeta is predicted to ten digits",
- "swampland forces exactly this contact term",
- "SL(2)_F is an exact fundamental global symmetry".
