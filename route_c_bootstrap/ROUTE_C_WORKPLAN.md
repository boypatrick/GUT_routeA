# Route C Workplan: GUT Amplitude Bootstrap

This is the working plan for Route C.  The goal is not to declare that the GUT
has already been derived from bootstrap axioms.  The goal is to build a
sequence of increasingly sharp consistency filters that can narrow, support, or
refute conditional action-level GUT branches.

## Core Question

Can the Route-A representation/index skeleton be extended to an action-level
GUT branch by on-shell consistency?

In symbolic form, Route C asks whether

```text
Lorentz + unitarity + factorization + Ward identities
+ anomaly cancellation + controlled high-energy behavior
```

is strong enough to select or strongly prefer a small set of GUT branches from

```text
SU(5), Pati-Salam, Spin(10), E6, superstring-derived completions.
```

The expected honest outcome is one of:

1. `Spin(10)` is preferred only after explicit minimality/no-exotics
   assumptions.
2. several branches survive, with a precise list of missing assumptions.
3. the current conditional action assumptions fail a concrete consistency test.
4. a string-like UV completion becomes necessary because finite-pole field
   theory cannot satisfy the required high-energy softness.

## Non-Negotiable Boundary

Route C must not be written as:

```text
four axioms => our full GUT.
```

It should be written as:

```text
on-shell consistency => filters on allowed GUT action branches.
```

Route A remains the representation/index/contact skeleton.  Route C is an
action-level audit program.

## Stage C0: Representation Pre-Ledger

### Purpose

Fix the external-state and candidate-group bookkeeping before writing any
amplitude.

### Current artifact

```text
code/bootstrap_candidate_ledger.py
output/candidate_ledger.json
output/candidate_ledger.md
```

### Verification conditions

- The visible family content is

```text
Q + L + u^c + d^c + nu^c + e^c.
```

- The filter is explicitly stated:

```text
simple group
+ one irreducible family object
+ contains nu^c
+ no chiral SM exotics
+ minimal single-object realization.
```

- The pre-ledger does not claim to be an amplitude bootstrap.

### Current result

Under the stated filter, the standard minimal survivor is

```text
Spin(10):16.
```

`SU(5)`, Pati-Salam, and `E6` are not refuted as theories; they fail this
specific pre-amplitude filter for stated reasons.

## Stage C1: External State Algebra and Charges

### Status

Implemented as the first exact charge/anomaly ledger.  This remains
pre-amplitude bookkeeping; it fixes conventions before pole residues or Ward
identities are introduced.

### Current artifacts

```text
code/external_state_charge_table.py
output/external_state_charges.json
output/external_state_charges.md
```

### Work items

1. Build a machine-readable external-state table:

```text
state, SU(3)_C, SU(2)_L, Y, B-L, chirality, family index.
```

2. Add a table of possible bilinear currents:

```text
J^a_mu = psi_bar gamma_mu T^a psi.
```

3. Mark which currents are SM-preserving and which are broken-GUT currents.

### Verification conditions

- Hypercharge normalization reproduces

```text
Tr Y^2 / Tr T_3L^2 = 5/3
```

for one full family plus `nu^c`.

- All external states have consistent conjugation conventions.

- Gauge-invariant external amplitudes have total SM charge zero.

### Current result

For one left-handed family

```text
Q + L + u^c + d^c + nu^c + e^c,
```

the generated ledger verifies exactly

```text
Tr Y^2 = 10/3,
Tr T_3L^2 = 2,
Tr Y^2 / Tr T_3L^2 = 5/3.
```

The following anomaly checks vanish in exact rational arithmetic:

```text
SU(3)^3,
SU(3)^2 U(1)_Y,
SU(2)^2 U(1)_Y,
grav^2 U(1)_Y,
U(1)_Y^3,
SU(3)^2(B-L),
SU(2)^2(B-L),
grav^2(B-L),
(B-L)^3,
U(1)_Y^2(B-L),
U(1)_Y(B-L)^2.
```

The Witten `SU(2)` check also passes with four left-handed doublets.

### Failure modes

- A candidate requires external states outside the declared visible family.
- A candidate silently changes charge normalization.
- A channel uses inconsistent left-handed/right-handed conjugation.

## Stage C2: Minimal Four-Point Pole Ansatz

### Status

Implemented as the first charge-filtered symbolic pole ledger.  It does not
assert mediator existence or residue positivity; it only lists pole components
not forbidden by the P1 SM-face quantum numbers.

### Current artifacts

```text
code/four_point_pole_ansatz.py
output/four_point_pole_ansatz.json
output/four_point_pole_ansatz.md
```

### Work items

Construct the first symbolic tree-level four-point ansatz:

```text
A_4(1,2,3,4)
  = sum_X R_s(X)/(s - M_X^2)
  + sum_Y R_t(Y)/(t - M_Y^2)
  + sum_Z R_u(Z)/(u - M_Z^2)
  + contact.
```

At first, keep the residues symbolic:

```text
R_s(X) = g_12X g_34X.
```

The implemented all-incoming chiral-bilinear seed is

```text
A_4(12|34) = sum_X R_s(X)/(s-M_X^2) + C6(12|34),
R_s(X) = g_12X g_34Xbar.
```

### Verification conditions

- Each pole has a declared intermediate representation.

- Each residue factorizes into two three-point couplings.

- No pole appears in a channel where charge conservation forbids it.

- Contact terms are separated from pole terms and assigned a mass dimension.

### Current result

Using the P1 external states, the generated ledger finds:

```text
pair types = 21,
pair-product components = 31,
tested pair-component pairings = 496,
allowed symbolic pole components = 18,
forbidden by Abelian charge = 455,
forbidden by non-Abelian conjugacy after Abelian pass = 23.
```

The diagnostic proton-operator seeds are:

```text
Q Q Q L: 2 allowed symbolic pole components,
u^c u^c d^c e^c: 2 allowed symbolic pole components,
Q Q u^c e^c: 0,
Q L u^c d^c: 0.
```

These counts only mean that charge conservation and pair-product conjugacy do
not forbid the symbolic pole.  P3--P7 must still check residue positivity,
Ward identities, high-energy behavior, and low-energy bounds.

### Failure modes

- A residue cannot be factorized into legal three-point data.

- A pole requires a mediator outside the candidate action branch.

- A contact term is used to hide a Ward-identity failure.

## Stage C3: Factorization and Residue Positivity

### Status

Implemented as a conditional Gram-matrix residue ledger.  This stage verifies
that the P2 symbolic poles factorize and that the resulting residue blocks are
positive semidefinite if the declared mediator has positive Hilbert-space norm.
It does not yet prove mediator existence, generator algebra, spin-statistics,
Ward identities, high-energy softness, or proton safety.

### Current artifacts

```text
code/residue_positivity_checker.py
output/residue_positivity.json
output/residue_positivity.md
```

### Work items

For each allowed pole,

```text
Res_{s=M_X^2} A_4 = g_12X g_34X.
```

Build a residue matrix in the space of two-particle channels.

### Verification conditions

- Physical residues have positive norm in unitary channels.

- Negative norm or wrong-sign residues are flagged.

- Degenerate mediator multiplets have consistent residue patterns.

- Residues respect the candidate group's generator algebra.

### Mathematical test

For a channel basis `alpha,beta`,

```text
R_{alpha beta} = sum_X g_{alpha X} g^*_{beta X}
```

must be positive semidefinite for physical intermediate states.

### Current result

The generated P3 ledger verifies:

```text
mediator sectors = 26,
total vertex channels = 31,
P2 allowed poles checked = 18,
factorization checks passed = 18,
PSD blocks under positive metric = 26,
largest channel-block size = 3.
```

For each mediator block, the unit-coupling positive-metric spectrum is of the
form

```text
[n,0,...,0],
```

while the wrong-sign mediator stress test gives

```text
[-n,0,...,0],
```

which flags the required positive-norm assumption.

## Stage C4: Ward Identities

### Status

Implemented as a Ward-identity readiness ledger.  P4 separates unbroken SM
current sectors that are bookkeeping-ready from broken/off-face sectors that
require explicit massive-vector, Goldstone, Higgs, or source data before a real
Ward-identity cancellation can be claimed.

### Current artifacts

```text
code/ward_identity_ledger.py
output/ward_identity_ledger.json
output/ward_identity_ledger.md
```

### Work items

Implement symbolic Ward checks for amplitudes with external gauge bosons or
internal mediator currents.

For gauge-boson polarization,

```text
epsilon_mu -> p_mu
```

the full amplitude must vanish after summing all required diagrams.

### Verification conditions

- Ward identities hold for unbroken SM gauge bosons.

- Broken-generator amplitudes include the required Goldstone/Higgs-sector
  contributions.

- Missing diagrams are reported as missing action-level data, not ignored.

### Current result

The generated P4 ledger verifies the following bookkeeping classification:

```text
unbroken SM current sectors = 3,
current transitions checked = 36,
massless diagonal transitions ready = 6,
broken/off-face transitions needing completion = 30,
mediator sectors inspected = 26,
mediator sectors needing completion if vector = 26.
```

For unbroken currents the target remains

```text
epsilon_mu -> p_mu,    p_mu M^mu = 0.
```

For broken massive-vector currents the required target is the generalized
Ward/Slavnov-Taylor relation

```text
p_mu M^mu(V_X) = M_X M(phi_X).
```

Thus P4 currently reports the missing action-level data rather than hiding it:

```text
component-level SU(3)_C and SU(2)_L generator matrices,
explicit broken Spin(10)/Pati-Salam generator matrices,
mass matrix for broken vector mediators,
Goldstone/Higgs/source sector defining p_mu M^mu(V)=M M(phi),
spin-statistics and numerator structures for the actual amplitudes.
```

### Failure modes

- An isolated broken-gauge-boson exchange grows badly at high energy.

- A candidate requires a Higgs/source sector not present in the branch.

- The current is not conserved in the massless limit.

## Stage C5: Anomaly Cancellation and Chiral Consistency

### Status

Implemented as an anomaly and chiral-consistency comparison ledger.  P5 uses
the P1 one-family `SM + nu^c` anomaly cancellation as the baseline and compares
`SU(5)`, Pati-Salam, `Spin(10)`, `E6`, and superstring-derived completions.

### Current artifacts

```text
code/anomaly_chiral_consistency_ledger.py
output/anomaly_chiral_consistency.json
output/anomaly_chiral_consistency.md
```

The P1--P5 paper-facing derivations are consolidated in:

```text
tex/route_c_derivation_ledger.tex
```

### Work items

Build an anomaly ledger for each candidate representation:

```text
G^3 anomaly
G^2 U(1) anomaly
U(1)^3 anomaly
mixed gravitational anomaly
Witten SU(2) anomaly
```

### Verification conditions

- One family plus `nu^c` is SM-anomaly free.

- Candidate GUT representation is chiral but anomaly consistent.

- Vectorlike extra states are explicitly paired and liftable.

### Candidate expectations

- `SU(5)`: anomaly cancellation works for `10 + bar5`, but not as one
  irreducible family object.
- Pati-Salam: natural family organization, but product group rather than
  simple group.
- `Spin(10)`: `16` is the standard minimal anomaly-safe single object.
- `E6`: `27` is consistent but includes extra states needing a lifting audit.
- Superstring-derived models: anomaly cancellation may follow from higher
  dimensional Green-Schwarz or modular consistency, but compactification data
  must still produce the visible spectrum.

### Current result

The generated P5 ledger reports:

```text
candidates inspected = 5,
field-theory candidates = 4,
UV-completion classes = 1,
Route-C minimal single-object survivors = Spin(10),
baseline all anomalies pass = true.
```

Candidate status:

```text
SU(5): anomaly-consistent as 10 + bar5 + 1, but not one irreducible object.
Pati-Salam: anomaly-consistent and natural, but product group and not one object.
Spin(10): anomaly-safe standard minimal single-object survivor via 16.
E6: anomaly-consistent but requires 10 + 1 extra-state lifting audit.
Superstring-derived completions: UV-completion class; compactification-dependent.
```

P5 therefore does not claim that anomaly cancellation alone selects the GUT.
It separates anomaly consistency from Route-C structural filters and from
compactification/exotics audits.

The TeX derivation ledger records the P1 charge/anomaly formulas, P2 symbolic
pole ansatz, P3 Gram-matrix positivity argument, P4 Ward-completion boundary,
and P5 candidate anomaly/chiral-consistency comparison so these derivations do
not live only in generated Markdown and JSON outputs.

## Stage C6: Controlled High-Energy Behavior

### Status

Implemented as a controlled high-energy-growth audit.  P6 consumes the P3
residue-factorization ledger and the P4 Ward-completion ledger.  It classifies
which symbolic branches would have uncontrolled longitudinal-vector growth if
treated as isolated broken-vector exchange, which scalar/auxiliary
interpretations avoid that particular Ward problem, and when a tower-like UV
completion should become a comparison point.

### Current artifacts

```text
code/high_energy_growth_audit.py
output/high_energy_growth_audit.json
output/high_energy_growth_audit.md
```

### Work items

Expand candidate amplitudes at

```text
s >> M_X^2
```

and classify high-energy growth.

### Verification conditions

- No uncancelled growth of the schematic form

```text
A(s) ~ s^2/M_X^4
```

in channels expected to remain perturbative.

- Longitudinal amplitudes obey equivalence-theorem checks.

- Any required cancellation is traced to a specific gauge/Higgs/source
  relation.

### Failure modes

- A branch needs unaudited states to cancel high-energy growth.

- A contact term is tuned without a symmetry reason.

- A finite-pole field-theory branch fails, suggesting a tower or UV completion
  may be needed.

### Current result

The generated P6 ledger reports:

```text
mediator sectors inspected = 26,
finite-vector branches needing completion = 26,
unbroken-vector branches bookkeeping-ready = 0,
scalar/auxiliary interpretations without longitudinal-vector risk = 26,
contact-only components needing UV completion = 18,
broken transitions needing completion = 30,
tower-like UV forced now = false.
```

Thus every P3 mediator sector, if interpreted as a finite broken vector,
requires explicit Goldstone/Higgs/source data satisfying

```text
p_mu M^mu(V_X)=M_X M(phi_X).
```

The contact-only EFT interpretation is not a UV completion because a
four-fermion contact normalization generically scales as

```text
A_contact ~ s/Lambda^2.
```

P6 does not force a superstring/tower completion.  It says a tower-like UV
completion becomes a serious comparison class only if finite-pole
Goldstone/Higgs/source completions cannot cancel the high-energy growth.

## Stage C7: Low-Energy Matching and Proton Bounds

### Status

Implemented as a low-energy matching and proton-bound readiness report.  P7
uses P2/P3 symbolic pole components only after attaching the P6 completion
status.  It classifies operator `B` and `L` content, writes symbolic
dimension-six matching coefficients, and marks physical proton bounds as
conditional until the missing completion/flavor/hadronic inputs are supplied.

### Current artifacts

```text
code/low_energy_matching_proton_report.py
output/low_energy_matching_proton_report.json
output/low_energy_matching_proton_report.md
```

### Work items

Integrate out heavy GUT mediators to obtain effective operators:

```text
L_eff ~ (g_GUT^2/M_X^2) qqql
```

and, where applicable,

```text
C_5 ~ Y_T M_T^{-1} Y_T.
```

### Verification conditions

- Operator basis is gauge invariant under the SM.

- Wilson coefficients are expressed in physical flavor bases or marked as
  pre-flavor-audit quantities.

- Proton-decay limits are applied to the correct channels:

```text
p -> e^+ pi^0
p -> K^+ anti-nu
```

### Failure modes

- A branch produces proton decay above experimental bounds.

- Flavor rotations are missing but the result is presented as physical.

- A cancellation is asserted before Wilson tensors are computed.

### Current result

The generated P7 ledger reports:

```text
allowed symbolic poles processed = 18,
B and L conserving poles = 12,
standard BNV seed poles = 4,
BNV-with-sterile-neutrino poles = 2,
baryon-violating symbolic poles = 6,
BNV vector poles blocked by P6 completion = 6,
physical proton bounds evaluable now = 0,
all physical bounds conditional = true.
```

The symbolic matching template is

```text
C6[operator; X] = g_leftX g_rightX^*/M_X^2.
```

The physical proton-bound interface is only a template at this stage:

```text
Gamma(p -> channel a) = sum_ij C_i^phys C_j^{phys *} H_ij^(a).
```

Numerical bounds require a P6-completed high-energy branch, mediator masses and
couplings, physical flavor rotations, RG factors, lattice/chiral matrix
elements, and the experimental lifetime limit for the selected channel.  P7
therefore records bound readiness rather than quoting lifetimes.

## Stage C8: Candidate-Theory Comparison

### Status

Implemented as a comparative theory scorecard.  P8 compares `SU(5)`,
Pati-Salam, `Spin(10)`, `E6`, and superstring-derived completions using the
explicit P1--P7 gates rather than anomaly cancellation or minimality alone.
It is a hard-gate ledger, not a weighted score.

### Current artifacts

```text
code/theory_comparison_scorecard.py
output/theory_comparison_scorecard.json
output/theory_comparison_scorecard.md
```

The final comparison should happen only after Stages C1-C7 produce concrete
outputs.

| theory | what it tests well | expected weakness | Route-C decision criterion |
| --- | --- | --- | --- |
| `SU(5)` | simple GUT economy; classic proton channels | family split as `10+bar5+1`; no single six-face object | survives only if single-object minimality is relaxed |
| Pati-Salam | quark-lepton organization; `B-L` clarity | not simple; family split into two reps | survives as an intermediate face, not as final simple object |
| `Spin(10)` | one `16` contains all six faces plus `nu^c` | action-level Higgs/source/proton data still conditional | preferred if amplitudes close without unaudited exotics |
| `E6` | larger simple unification; contains `Spin(10)` | `27` brings extra states | survives only with explicit lifting and threshold audit |
| superstring theory | UV-soft amplitudes, infinite tower, anomaly/modular consistency | compactification choices, moduli, exotics, landscape non-uniqueness | relevant if finite-pole GUT bootstrap needs a tower or if string compactification yields the Route-A skeleton |

### Current result

The generated P8 scorecard reports:

```text
candidates compared = 5,
leading conditional field-theory branch = Spin(10),
tower-like UV forced now = false.
```

Final status counts:

```text
conditional_branch_if_single_object_filter_is_relaxed = 1,
conditional_intermediate_face_not_final_simple_object = 1,
leading_conditional_route_c_branch = 1,
conditional_branch_with_exotics_lifting_debt = 1,
uv_completion_class_not_forced_by_p1_p7 = 1.
```

Candidate verdicts:

```text
SU(5): conditional branch if the single-object filter is relaxed.
Pati-Salam: conditional intermediate face, not a final simple object.
Spin(10): leading conditional Route-C finite field-theory branch.
E6: conditional branch with exotics-lifting debt.
Superstring-derived completions: UV-completion class, not forced by P1-P7.
```

P8 therefore closes the first Route-C pass as a comparative filter.  It does
not convert `Spin(10)` into a completed action-level GUT; it identifies it as
the leading conditional finite field-theory branch under the current gates.

## Superstring-Theory Comparison

Superstring theory should be treated as a UV-completion class, not simply as
one more four-dimensional finite GUT group.

### What to compare

1. **Amplitude softness**

String amplitudes are not finite-pole rational functions.  They contain an
infinite tower and often show Regge-soft behavior.  If Route C finds that a
finite set of GUT poles cannot satisfy high-energy softness, a string-like
tower becomes a serious comparison point.

2. **Factorization**

String amplitudes factorize on infinitely many massive poles.  Route C should
ask whether the GUT mediator pole is the first pole of a tower or an isolated
field-theory state.

3. **Anomaly cancellation**

String constructions may cancel anomalies through higher-dimensional
mechanisms such as Green-Schwarz terms.  The four-dimensional compactified
spectrum must still satisfy the low-energy anomaly and exotics checks.

4. **Compactification output**

For this project, a useful string comparison must show whether a compactified
model can produce:

```text
Spin(10) or a Spin(10)-like face,
three protected families,
no light chiral exotics,
acceptable proton bounds,
and a controlled Majorana/contact sector.
```

### What not to claim

Do not claim:

```text
Route C proves superstring theory.
```

The correct statement is:

```text
If finite-field-theory GUT branches fail high-energy bootstrap constraints,
then string-like towers become a candidate UV-completion class to compare.
```

## Deliverables and Order

### P0. Candidate ledger

Status: done.

Deliverables:

```text
output/candidate_ledger.json
output/candidate_ledger.md
```

### P1. External state and charge table

Deliverables:

```text
output/external_state_table.json
output/external_state_table.md
```

### P2. Symbolic pole ansatz

Deliverables:

```text
code/bootstrap_pole_ansatz.py
output/pole_ansatz.json
```

### P3. Factorization and positivity checker

Deliverables:

```text
code/residue_positivity_checker.py
output/residue_positivity.json
output/residue_positivity.md
```

### P4. Ward-identity checker

Deliverables:

```text
code/ward_identity_ledger.py
output/ward_identity_ledger.json
output/ward_identity_ledger.md
```

### P5. Anomaly ledger

Deliverables:

```text
code/anomaly_chiral_consistency_ledger.py
output/anomaly_chiral_consistency.json
output/anomaly_chiral_consistency.md
```

### P6. High-energy behavior audit

Deliverables:

```text
code/high_energy_growth_audit.py
output/high_energy_growth_audit.json
output/high_energy_growth_audit.md
```

### P7. Low-energy matching and proton-bound report

Deliverables:

```text
code/low_energy_matching_proton_report.py
output/low_energy_matching_proton_report.json
output/low_energy_matching_proton_report.md
```

### P8. Comparative theory scorecard

Deliverables:

```text
code/theory_comparison_scorecard.py
output/theory_comparison_scorecard.json
output/theory_comparison_scorecard.md
```

The comparison must include at least:

```text
SU(5), Pati-Salam, Spin(10), E6, superstring-derived completions.
```

### P9. Spin(10)-locked action-completion ledger

Status: done.

Deliverables:

```text
code/spin10_action_completion_ledger.py
output/spin10_action_completion_ledger.json
output/spin10_action_completion_ledger.md
```

P9 stops broadening the comparison and locks the finite field-theory branch to
the P8-leading `Spin(10)` option.  It supplies the first action-completion data
requested by P4--P7:

```text
16 -> (4,2,1) + (bar4,1,2),
```

with a 16-state Pati-Salam/SM-face basis and sparse broken-generator maps.

Current generated data:

```text
basis dimension = 16,
broken generator sparse maps = 32,
SU(4)_C leptoquark root maps = 6,
SU(2)_R charged root maps = 2,
Spin(10)/Pati-Salam (6,2,2) Clebsch maps = 24.
```

Verification:

```text
P4 broken transition pairs = 30,
Spin(10) adjoint vector transition pairs generated = 22,
non-adjoint P4 completion pairs = 8,
generated pairs are subset of P4 broken pairs = true,
P9 generator layer passes = true.
```

Boundary:

P9 does not turn every P4 charge-bookkeeping transition into a `D5` gauge
generator.  The 8 non-adjoint pairs remain possible scalar/source/auxiliary
completion channels.  P9 also supplies Goldstone/Higgs completion formulae,
mediator mass/coupling symbols, physical flavor rotation templates, and the
proton-bound input interface, but it does not choose a Higgs potential,
mediator spectrum, flavor fit, or numerical proton lifetime.

Next:

```text
P10_choose_spin10_breaking_branch_and_mass_spectrum
```

Choose a concrete `Spin(10)` breaking order-parameter/source branch, diagonalize
the broken-vector mass matrix, and replay P6/P7 with explicit masses,
couplings, and flavor rotations.

### P10. Spin(10) staged breaking mass-matrix ledger

Status: done.

Deliverables:

```text
code/spin10_breaking_mass_matrix.py
output/spin10_breaking_mass_matrix.json
output/spin10_breaking_mass_matrix.md
```

P10 chooses the conservative `staged_orthogonal_source_branch`.  It is not a
full Higgs potential.  It is the first concrete source branch that promotes the
P9 symbolic mass inputs to a diagonalizable broken-vector mass matrix while
preserving the hypercharge zero mode.

Chosen breaking chain:

```text
Spin(10) -> SU(4)_C x SU(2)_L x SU(2)_R
SU(4)_C -> SU(3)_C x U(1)_{B-L}
SU(2)_R -> U(1)_{T3R}
U(1)_{T3R} x U(1)_{B-L} -> U(1)_Y
```

Mass blocks:

```text
M_(6,2,2)^2 = kappa_PS g_10^2 v_PS^2        dimension 24
M_LQ^2      = (2/3) g_4^2 v_4^2             dimension 6
M_WR^2      = g_R^2 v_R^2                   dimension 2
M_neutral^2 = v_Y^2/4 [[g_R^2, -g_R g_BL],
                       [-g_R g_BL, g_BL^2]] dimension 2
```

The `2/3` leptoquark coefficient follows from

```text
T_BL = sqrt(3/8) diag(1/3,1/3,1/3,-1),
Delta T_BL^2 = (4/3)^2 (3/8) = 2/3.
```

Verification:

```text
active matrix dimension including hypercharge zero mode = 34,
positive broken eigenvalue count = 33,
zero eigenvalue count = 1,
all P9 generator maps assigned to mass blocks = true,
P10 mass matrix layer passes = true.
```

Boundary:

P10 does not choose numerical scales, threshold spectrum, physical flavor
rotations, RG factors, hadronic matrix elements, current proton lifetime
limits, or the scalar/source branch for the 8 non-adjoint P4 pairs.

Next:

```text
P11_replay_P6_P7_with_staged_spin10_masses
```

Choose symbolic or numerical ranges for `v_PS`, `v_4`, `v_R`, `v_Y`, and the
couplings, then replay the P6/P7 matching gates with these explicit mass
blocks.

### P11. P6/P7 replay with staged Spin(10) mass blocks

Status: done.

Deliverables:

```text
code/spin10_mass_replay_matching_gate.py
output/spin10_mass_replay_matching_gate.json
output/spin10_mass_replay_matching_gate.md
```

P11 attaches the P10 staged mass blocks to the matching gates:

```text
M_(6,2,2) -> sqrt(kappa_PS) g_10 v_PS
M_LQ      -> sqrt(2/3) g_4 v_4
M_WR      -> g_R v_R
M_Zprime  -> 0.5 v_Y sqrt(g_R^2+g_BL^2)
```

It then checks whether the 18 P7 chiral-pair pole rows can legitimately use
those adjoint-current vector mass blocks.

Verification:

```text
current vector mass blocks available = 4,
P7 rows replayed = 18,
P7 rows with candidate adjoint mass blocks = 0,
P7 rows remaining scalar/source mass debt = 18,
baryon-violating or sterile-BNV rows = 6,
physical proton bounds evaluable now = 0.
```

Boundary:

P11 improves the adjoint current-vector side of the matching problem, but it
does not turn the P7 chiral-pair pole ledger into adjoint-vector exchange.  The
P7 rows remain scalar/source mass debt unless a later branch supplies the
appropriate scalar/source masses and Wilson-basis map.

Next:

```text
P12_choose_vector_current_or_scalar_source_proton_branch
```

Decide whether proton matching proceeds through completed `Spin(10)`
current-vector exchange or through scalar/source chiral-pair poles.  The two
branches require different Wilson-basis maps.

### P12-S. Scalar/source chiral-pair branch

Status: done for the first scalar/source bookkeeping layer.

Deliverables:

```text
code/scalar_source_chiral_pair_branch.py
output/scalar_source_chiral_pair_branch.json
output/scalar_source_chiral_pair_branch.md
```

P12 records both required post-P11 directions:

```text
Branch V: completed Spin(10) current-current vector exchange using P10
          adjoint masses and an explicit current Wilson basis.
Branch S: scalar/source chiral-pair poles using the P7 all-incoming
          chiral-pair ledger and scalar/source mass denominators.
```

The branch activated first is Branch S, because P11 showed that all 18 P7
chiral-pair rows remain scalar/source mass debt rather than adjoint-vector
matches.

Scalar/source ansatz:

```text
L ⊃ lambda_{ab,R} psi_a psi_b S_R + h.c. + M_{S_R}^2 |S_R|^2
C6(12|34;R)=lambda_{12,R} lambda^*_{34,R}/M_{S_R}^2
M_{S_R}^2 > 0
```

Generated summary:

```text
P7 rows converted to scalar/source rows = 18,
unique scalar/source sectors = 12,
B and L conserving rows = 12,
standard BNV seed rows = 4,
BNV-with-sterile-neutrino rows = 2,
baryon-violating or sterile-BNV rows = 6,
physical proton bounds evaluable now = 0.
```

Verification:

```text
all P11 scalar/source debt rows accounted for = true,
every scalar row has a mass symbol = true,
every source sector has a positive mass condition = true,
P12 scalar/source layer passes = true.
```

Boundary:

P12-S removes the ambiguity that the P7 chiral-pair rows were secretly adjoint
vector rows.  It does not compute proton lifetimes.  Numerical proton limits
still require scalar/source flavor tensors, a chiral/Fierz Wilson basis,
physical flavor rotations, RG evolution, hadronic matrix elements, and
experimental channel limits.  Branch V remains a required future attempt and
must not be discarded.

Next:

```text
P13_scalar_source_flavor_tensor_and_operator_basis
```

Choose scalar/source flavor tensors and an operator basis for the six BNV or
sterile-BNV rows, or return to Branch V before inserting proton limits.

### P13-S. Scalar/source flavor tensors and Wilson basis

Status: done for the first symbolic flavor/Wilson-basis layer.

Deliverables:

```text
code/scalar_source_flavor_wilson_basis.py
output/scalar_source_flavor_wilson_basis.json
output/scalar_source_flavor_wilson_basis.md
```

P13-S chooses a conservative symbolic flavor tensor basis rather than a numeric
flavor fit.  It maps the six P12-S BNV or sterile-BNV rows to four all-left
chiral operator classes:

```text
O_QQQL_S,
O_QQQL_T,
O_UUDE,
O_UDDN.
```

The chosen flavor tensors are generic complex family tensors subject only to
the pair-exchange symmetries forced by the source channel:

```text
QQ singlet source: symmetric family pair,
QQ triplet source: antisymmetric family pair,
u^c u^c source: antisymmetric family pair,
d^c d^c source: antisymmetric family pair,
distinct-field pairs: general 3x3 family matrix.
```

Generated summary:

```text
BNV or sterile-BNV rows processed = 6,
canonical Wilson basis count = 4,
direct source-channel rows = 4,
Fierz/color recoupling required rows = 2,
canonical family projection required rows = 1,
flavor tensors introduced = 11,
independent complex vertex components = 78,
physical proton bounds evaluable now = 0.
```

Verification:

```text
all P12 BNV rows mapped = true,
every row has a basis id = true,
every row has left and right flavor tensors = true,
Branch V retained = true,
P13 scalar/source basis layer passes = true.
```

Boundary:

P13-S chooses the symbolic Wilson basis and tensor shapes only.  It does not
compute Fierz/color recoupling matrices, physical flavor rotations, RG factors,
hadronic matrix elements, or proton lifetimes.  Branch V remains open and must
still be tried before any final proton-matching conclusion.

Next:

```text
P14_vector_branch_matching_gate_or_scalar_source_recoupling
```

Either return to the completed current-current vector Branch V before inserting
proton limits, or compute the two non-direct scalar/source chiral/color Fierz
recoupling matrices and the direct `O_UDDN` antisymmetric family projection.

### P14-V. Vector-branch matching gate before proton limits

Status: done for the first Branch V current-current matching gate.

Deliverables:

```text
code/spin10_vector_branch_matching_gate.py
output/spin10_vector_branch_matching_gate.json
output/spin10_vector_branch_matching_gate.md
```

P14-V returns to the Spin(10) current-current vector branch instead of
continuing directly to proton limits.  It consumes the P9 broken-generator maps,
the P10 staged vector mass blocks, and the P11 replay gate.  The matching
template is

```text
J_X^mu = sum_a c_a psi_target(a)^dagger bar_sigma^mu psi_source(a),
L_eff[X] = - g_X^2 J_X^mu J_X,mu^dagger / M_X^2.
```

Generated summary:

```text
broken vector maps = 32,
current transition entries = 80,
current pair rows = 72,
Pati-Salam SU(4)_C leptoquark pair rows = 36,
Spin(10)/Pati-Salam off-face pair rows = 24,
broken SU(2)_R charged-current pair rows = 12,
physical proton bounds evaluable now = 0,
Branch V matching gate passes = true.
```

Gate counts:

```text
Pati-Salam leptoquark mixed pairs requiring crossing/Fierz = 24,
Pati-Salam leptoquark same-sector noncanonical pairs = 12,
off-face diquark-like pairs requiring crossing/Fierz = 12,
off-face mixed pairs requiring crossing/Fierz = 12,
SU(2)_R charged-current non-proton-seed rows = 12.
```

Boundary:

P14-V does not reuse the P7 scalar/source chiral-pair rows as adjoint vectors.
It explicitly keeps the P11 result:

```text
P7 rows with candidate adjoint mass blocks = 0,
P7 rows remaining scalar/source mass debt = 18.
```

Thus Branch V is now an explicit current-current matching gate, not a physical
proton lifetime calculation.  Numerical proton limits remain blocked by:

```text
current-current to all-left chiral Wilson basis recoupling,
physical flavor rotations,
numerical vector masses and couplings,
RG running,
hadronic matrix elements,
experimental channel limits.
```

Next:

```text
P15_V_current_to_all_left_fierz_and_flavor_gate
```

Choose the convention-fixed current-current to all-left operator map for the
vector branch, or keep P13-S scalar/source recoupling as the active
proton-matching path.  Do not insert proton limits until one branch has flavor,
RG, and hadronic inputs.

### P15-V. Vector current to all-left Fierz/flavor gate

Status: done for the first fixed Fierz/crossing convention gate.

Deliverables:

```text
code/spin10_vector_all_left_fierz_gate.py
output/spin10_vector_all_left_fierz_gate.json
output/spin10_vector_all_left_fierz_gate.md
```

P15-V fixes the two-component Weyl identity

```text
(chi^dagger bar_sigma^mu psi)(eta^dagger bar_sigma_mu xi)
  = 2 (chi^dagger eta^dagger)(psi xi)
```

and applies it to the P14-V vector-current product

```text
(psi_t1^dagger bar_sigma^mu psi_s1)
(psi_s2^dagger bar_sigma_mu psi_t2)
```

using the crossed all-left skeleton

```text
(psi_s1 psi_t2)(bar(psi_t1) bar(psi_s2))
```

with symbolic coefficient prefactor

```text
-2 g_X^2 c1 c2^*/M_X^2.
```

Generated summary:

```text
P14 current-pair rows consumed = 72,
P15 rows generated = 72,
all-left status counts = {B_and_L_conserving_crossed_skeleton: 72},
canonical BNV basis counts = {none: 72},
canonical BNV rows evaluable now = 0,
physical proton bounds evaluable now = 0,
all rows have crossed fields = true,
all rows have charge audit = true,
Branch V all-left gate passes = true.
```

Boundary:

P15-V does not produce a canonical BNV proton operator from the current sparse
P9/P14 vector-map pairing.  It means Branch V is now in a fixed crossed
all-left bookkeeping basis, but physical proton limits are still blocked.  The
missing vector-side data are:

```text
physical hermitian vector multiplet assembly,
cross-root Clebsch phases,
external-state crossing convention,
projection to canonical BNV Wilson basis if applicable,
physical flavor rotations,
RG running,
hadronic matrix elements,
experimental channel limits.
```

P13-S / Branch S remains a separate scalar/source Wilson-basis branch and must
not be discarded.

Next:

```text
P16_choose_vector_physical_multiplet_or_return_to_scalar_source
```

Either assemble physical hermitian vector multiplets and their cross-root
Clebsch phases to search for canonical BNV vector operators, or return to the
P13-S scalar/source recoupling path.  Do not insert proton limits yet.

### P16-V. Physical Hermitian vector multiplet assembly

Status: done for the first cross-root physical vector assembly gate.

Deliverables:

```text
code/spin10_vector_physical_multiplet_assembly.py
output/spin10_vector_physical_multiplet_assembly.json
output/spin10_vector_physical_multiplet_assembly.md
```

P16-V tests the missing vector-side datum identified by P15-V.  Instead of
using the same-map product `J_X J_X^\dagger`, it assembles charge-conjugate
root maps into symbolic physical Hermitian vector multiplets.  The physical
product is

```text
J_q^mu J_-q,mu
```

with an unfixed cross-root Clebsch/phase convention `xi`.  Applying the same
two-component all-left Fierz gate gives the symbolic coefficient template

```text
C_phys = -2 g_X^2 xi c_q c_-q / M_X^2.
```

Generated summary:

```text
P15 same-map canonical BNV rows = 0,
physical cross-root rows = 448,
Pati-Salam SU(4)_C leptoquark rows = 144,
Spin(10)/Pati-Salam off-face rows = 288,
broken SU(2)_R charged-current rows = 16,
B,L conserving physical pairs = 160,
B-L preserving BNV noncanonical pairs = 144,
canonical vector BNV basis candidates = 144,
O_V_QQ_DbarNbar candidates = 72,
O_V_QQ_UbarEbar candidates = 72,
physical proton bounds evaluable now = 0,
Branch V physical multiplet gate passes = true.
```

Interpretation:

P16-V changes the Branch-V status.  The same-map sparse-current gate in P15
found no canonical BNV operator.  The physical Hermitian cross-root gate does
generate canonical vector BNV skeletons, but only in the
`Spin(10)/Pati-Salam` off-face `M_(6,2,2)` block.  The Pati-Salam leptoquark
and broken `SU(2)_R` physical pairs remain `B,L` conserving in this gate.

Boundary:

P16-V is still not a proton-limit calculation.  The canonical vector BNV rows
are pre-flavor, pre-RG, and pre-hadronic.  The following data remain missing:

```text
physical Hermitian generator normalization,
cross-root Clebsch phase convention,
color and weak projection to the canonical vector BNV basis,
physical flavor rotations,
numerical masses and couplings,
RG running,
hadronic matrix elements,
experimental channel limits.
```

P13-S / Branch S remains a separate scalar/source Wilson-basis branch and must
not be discarded.

Next:

```text
P17_vector_clebsch_projection_and_flavor_interface
```

Fix the Hermitian generator normalization and cross-root Clebsch phase
convention for the 144 canonical vector BNV candidates, then export a
flavor-rotation interface.  Do not insert proton limits until RG and hadronic
inputs are supplied.

## First Decision Point

After P2-P4, decide whether the first mediator is modeled as:

1. a generic broken GUT gauge boson with symbolic residues; or
2. an explicitly `Spin(10)` broken generator.

The safer path is option 1 first.  It avoids assuming the answer before the
bootstrap filter has done any work.

Status after P16-V: the first broad comparison pass is complete, P9 locked the
finite field-theory branch to `Spin(10)`, P10/P11 supplied staged adjoint-vector
mass gates, P12-S activated the scalar/source chiral-pair interpretation for
the P7 rows, P13-S selected the symbolic scalar/source flavor/Wilson basis, and
P14-V returned to Branch V with an explicit current-current matching gate.
P15-V then fixed the first all-left Fierz/crossing convention and found that the
current sparse-map pairing gives only B/L-conserving crossed skeletons.  This is
not the end of Branch V: P16-V assembled physical Hermitian cross-root vector
pairs and found 144 canonical vector BNV skeletons in the `M_(6,2,2)` off-face
block.  Physical proton bounds remain unavailable until Clebsch phases, flavor,
RG, hadronic, and channel-limit data are supplied.  This is still a conditional
branch, not a completed GUT proof or a proton-lifetime calculation.
