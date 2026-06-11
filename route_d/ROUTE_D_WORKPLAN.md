# Route D Workplan: String-Constrained Lift

## Purpose

Route D tests whether the conditional Route-A/Route-B data can be embedded in
standard string compactification structures.  The goal is not to claim a unique
string vacuum.  The goal is to replace some arbitrary EFT choices by discrete
topological data, index formulae, and computable moduli-dependent functions.

## Core Question

Can the Route-A skeleton

- `Spin(10):16` as one six-face family object,
- `H^0(P1,O(2))` as three protected family sections,
- `K_tr` as the unique second-transvectant contact direction,

be realized as a string-local or string-constrained structure without
overstating vacuum uniqueness?

## D0 Scope and Boundary

Route D may explain or constrain the following open data:

- `P_{nu^c}` and the source-sector UV origin,
- the origin of the `O(2)` family carrier,
- the UV meaning of the Route-B Majorana-only selection rule,
- the hidden phase/radial sector supplying `zeta`,
- whether the benchmark elliptic/modular data match a compactification modulus.

Route D must not claim:

- a unique string vacuum,
- an exact prediction of the ten-digit benchmark `zeta`,
- full flavor/proton/threshold closure,
- or a PSLT-only derivation of the GUT.

## Milestones

### D1-F: F-theory local placement audit

Status: checked locally on 2026-06-11.

Target statements:

- `E6 -> SO(10) x U(1)` adjoint branching:
  `78 = 45_0 + 1_0 + 16_-3 + bar16_+3`.
- Matter curve zero-mode carrier:
  `K_{P1}^{1/2} tensor O(3) = O(2)`, hence `h0=3`, `h1=0`.
- Hypercharge flux/Wilson-line breaking remains conditional and must satisfy
  the standard massless-hypercharge topological condition.

Deliverable:

- derivation note:
  `tex/d1_f_theory_local_placement.tex`,
- compiled PDF:
  `tex/d1_f_theory_local_placement.pdf`,
- verification script:
  `code/verify_d1_f_theory_local_placement.py`,
- output:
  `output/d1_f_theory_local_placement.{json,md}`.

Result:

- The local branching and `P1` cohomology checks pass:
  `45+1+16+16=78`, and
  `K_{P1}^{1/2} tensor O(3)=O(2)` has `h0=3`, `h1=0`.
- This supports a conservative Route-D appendix paragraph about the local
  origin of the `Spin(10):16` matter curve and the `O(2)` carrier.
- It does not verify a global compactification, net chirality beyond the
  selected flux datum, absence of exotics, or the massless-hypercharge flux
  topological condition.

### D2-E3: E3-instanton origin for `zeta K_tr`

Status: checked as a conditional instanton-origin audit on 2026-06-11.

Target statements:

- `zeta = A exp(2 pi i tau_Gamma)` supplies phase/radial data from
  `(int_Gamma C4, Vol(Gamma))`.
- Majorana-only generation follows from charged zero-mode counting and
  Green-Schwarz selection rules, not from an anomaly-free global symmetry.
- Family-blind instanton data can only contribute along the invariant
  transvectant direction, if the relevant equivariance/locality assumptions
  hold.

Required caveats:

- existence of the instanton and its zero-mode structure is an assumption until
  a concrete compactification is given,
- multi-instanton/multi-cover corrections are too large to justify ten-digit
  prediction of `zeta`,
- `zeta != 0` is generic/expected under quantum gravity assumptions, not a
  universal theorem without specifying gauge/discrete symmetries.

Deliverable:

- derivation note:
  `tex/d2_e3_instanton_zeta_audit.tex`,
- compiled PDF:
  `tex/d2_e3_instanton_zeta_audit.pdf`,
- verification script:
  `code/verify_d2_e3_instanton_zeta.py`,
- output:
  `output/d2_e3_instanton_zeta.{json,md}`.

Result:

- The benchmark coefficient has
  `|zeta|=0.13043190325293763`,
  `arg(zeta)=0.600038020318215`,
  and effective action
  `S_eff=-log|zeta|=2.0369040025655094`
  for unit prefactor.
- In the convention `zeta=A exp(2 pi i T_Gamma)`, the hidden phase/radial
  sector is naturally an axion-volume pair.
- A minimal toy GS charge ledger verifies the algebraic target:
  if `q(N)=+1`, then `q(NN)=+2`, so an instanton factor with charge `-2`
  makes the Majorana insertion neutral.
- The two-instanton scale `|zeta|^2=0.01701248138618368` is about
  `530.182` times the recorded `Delta_s=3.208798e-5` sensitivity window,
  so this route cannot justify a ten-digit prediction of `zeta` without a
  concrete compactification and correction audit.
- The audit does not construct an E3 divisor or count charged zero modes; it
  supports only a conditional string-origin interpretation of the Route-B
  hidden coefficient and selection rule.

### D3-H: heterotic comparison branch

Target statements:

- standard embedding `V=TX` solves the heterotic Bianchi identity at the
  textbook level and explains tangent-bundle language,
- an `SU(4)` bundle branch can have `SO(10)` commutant in `E8`,
- generation number is topological:
  `N_gen = 1/2 |int_X c3(V)|`,
- `chi(X)=+-6` gives three generations in known constructions.

Boundary:

- this is a comparison branch, not automatically the same compactification as
  the F-theory local picture.

### D4-M: modular/flavor benchmark audit

Target statement:

- compute the benchmark elliptic invariant from `(I,J)`,
  `j(E_M) proportional to I^3/(I^3 - 27 J^2)`,
  invert to a candidate `tau_M`, and test whether `zeta` or related contact
  data match low-weight modular functions.

Deliverable:

- a script in `code/`,
- output tables under `output/`,
- a pass/fail note.  If no simple modular match appears, record the negative
  result.

### D5-P: paper promotion decision

Status: decided on 2026-06-11.

Possible outcomes:

- Outlook paragraph only,
- short string-embedding appendix,
- separate Route-D companion note,
- or retire the string lift if the audits fail.

Decision:

- Promote D1-F and D2-E3 to a short optional appendix in the Route-A
  manuscript.
- Do not promote them into the theorem core.
- Keep all global compactification, hypercharge-flux, charged-zero-mode, and
  high-precision `zeta` caveats in the appendix text.

Implemented in:

- `paper/gut_framework.tex`, Appendix B,
  "Optional Route-D String-Constrained Placement".
- Supporting references added to `paper/refs.bib`.
- Main-paper reproducibility manifest updated to include the two Route-D
  scripts and outputs.

## Boundary

D5-P has promoted D1-F and D2-E3 to a conservative optional appendix.  Route D
still does not support a claim of global F-theory completion, exact string
prediction of the benchmark, or PSLT-only GUT derivation.
