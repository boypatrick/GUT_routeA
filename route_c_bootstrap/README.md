# Route C: GUT Amplitude Bootstrap

This directory is a separate workspace for testing a possible Route-C program:
use on-shell consistency conditions to narrow or refute conditional GUT
branches.  It is not part of the Route-A theorem core and should not be cited
as a completed derivation.

## Boundary

Route A currently establishes the representation/index/contact skeleton:

- `Spin(10):16` as the standard minimal six-face object.
- `H^0(CP^1,O(2))` as the protected three-family carrier.
- `K_tr` as the unique `SL(2)` second-transvectant contact direction.

Route C is different.  It asks whether action-level data can be constrained by
amplitude consistency:

- Lorentz invariance.
- Unitarity and positive residues.
- Pole factorization.
- Crossing symmetry where applicable.
- Gauge Ward identities.
- Anomaly cancellation.
- Controlled high-energy behavior.

## First target

The first nontrivial target is a four-point amplitude ledger with external
SM-plus-`nu^c` states embedded in a candidate family object,

```text
A(16_i 16_j -> 16_k 16_l).
```

A heavy mediator pole must factorize schematically as

```text
A_4 ~ g_ijX g_klX / (s - M_X^2) + regular.
```

The low-energy limit should match baryon-violating operators such as

```text
(g_GUT^2/M_X^2) qqql.
```

Before the pole ansatz, P1 fixes the external-state conventions and anomaly
bookkeeping:

```text
code/external_state_charge_table.py
output/external_state_charges.json
output/external_state_charges.md
```

The P1 ledger verifies `Tr Y^2 / Tr T3L^2 = 5/3`, one-family SM anomaly
cancellation with `nu^c`, the `B-L` anomaly checks, and the Witten `SU(2)`
even-doublet condition.

P2 then builds the first charge-filtered symbolic pole ledger:

```text
code/four_point_pole_ansatz.py
output/four_point_pole_ansatz.json
output/four_point_pole_ansatz.md
```

It declares pair-product mediator quantum numbers for

```text
A_4(12|34) = sum_X R_s(X)/(s-M_X^2) + C6(12|34),
```

and keeps every residue symbolic.  A P2-allowed pole is only allowed by
SM-face charge conservation and non-Abelian conjugacy; mediator existence,
positive residues, and Ward identities remain later checks.

P3 groups the P2 pole components into conditional residue blocks:

```text
code/residue_positivity_checker.py
output/residue_positivity.json
output/residue_positivity.md
```

For each declared mediator sector it records the formal Gram matrix

```text
R_{alpha beta}=g_{alpha X} g^*_{beta X}.
```

The P3 positivity statement is conditional on a positive-norm mediator.  A
wrong-sign mediator metric is explicitly recorded as a failing stress test.

P4 adds Ward-identity bookkeeping:

```text
code/ward_identity_ledger.py
output/ward_identity_ledger.json
output/ward_identity_ledger.md
```

It separates unbroken SM current sectors from broken/off-face current sectors.
For unbroken vectors the target is `p_mu M^mu=0`; for broken massive vectors
the target is the generalized relation
`p_mu M^mu(V_X)=M_X M(phi_X)`, so the missing Goldstone/Higgs/source data are
reported explicitly.

P5 adds the anomaly and chiral-consistency comparison ledger:

```text
code/anomaly_chiral_consistency_ledger.py
output/anomaly_chiral_consistency.json
output/anomaly_chiral_consistency.md
```

It compares `SU(5)`, Pati-Salam, `Spin(10)`, `E6`, and superstring-derived
completions against the P1 `SM + nu^c` anomaly baseline.  Passing anomaly
checks is treated as necessary but not sufficient; exotics, single-object
minimality, compactification data, and later amplitude checks remain separate.

P6 adds the controlled high-energy-growth audit:

```text
code/high_energy_growth_audit.py
output/high_energy_growth_audit.json
output/high_energy_growth_audit.md
```

It consumes the P3 residue blocks and P4 Ward-completion ledger.  A finite
broken-vector interpretation is marked incomplete until the corresponding
Goldstone/Higgs/source data enforce

```text
p_mu M^mu(V_X)=M_X M(phi_X).
```

Contact-only components are treated as low-energy EFT data, not UV softness
proofs.  Tower-like UV completion is not forced by P6; it is retained as a
comparison class if finite-pole completions fail.

P7 adds the low-energy matching and proton-bound readiness report:

```text
code/low_energy_matching_proton_report.py
output/low_energy_matching_proton_report.json
output/low_energy_matching_proton_report.md
```

It converts P2/P3 symbolic poles into conditional dimension-six matching rows,
classifies their `B` and `L` content, and keeps every physical proton-bound
claim conditional until P6 completion, flavor rotations, mediator masses and
couplings, RG factors, hadronic matrix elements, and experimental channel
limits are supplied.

P8 adds the comparative theory scorecard:

```text
code/theory_comparison_scorecard.py
output/theory_comparison_scorecard.json
output/theory_comparison_scorecard.md
```

It compares `SU(5)`, Pati-Salam, `Spin(10)`, `E6`, and superstring-derived
completions using explicit P1--P7 gates.  It reports `Spin(10)` as the leading
conditional finite field-theory branch under the current minimal single-object
filter, while keeping all action-level and phenomenological completion debts
visible.

P9 locks the branch to the P8-leading `Spin(10)` option and starts the
action-level completion ledger:

```text
code/spin10_action_completion_ledger.py
output/spin10_action_completion_ledger.json
output/spin10_action_completion_ledger.md
```

It builds a 16-state Pati-Salam/SM-face basis and sparse broken-generator maps
for the `Spin(10)` half-spinor.  The ledger supplies:

- 6 `SU(4)_C` leptoquark root maps;
- 2 broken `SU(2)_R` charged root maps;
- 24 `Spin(10)`/Pati-Salam `(6,2,2)` Clebsch maps.

The P9 check finds that these maps generate 22 `Spin(10)` adjoint-vector
transition pairs.  P4 has 30 charge-bookkeeping broken pairs, so the remaining
8 are recorded as non-adjoint scalar/source/auxiliary completion candidates
rather than being misidentified as `D5` gauge generators.

P10 chooses the first concrete staged/source breaking branch and turns the P9
symbolic mass inputs into a diagonalizable mass-matrix ledger:

```text
code/spin10_breaking_mass_matrix.py
output/spin10_breaking_mass_matrix.json
output/spin10_breaking_mass_matrix.md
```

The branch is:

```text
Spin(10) -> SU(4)_C x SU(2)_L x SU(2)_R
SU(4)_C -> SU(3)_C x U(1)_{B-L}
SU(2)_R -> U(1)_{T3R}
U(1)_{T3R} x U(1)_{B-L} -> U(1)_Y
```

It produces the block form

```text
diag(M_(6,2,2)^2 I_24, M_LQ^2 I_6, M_WR^2 I_2, M_neutral^2),
```

with

```text
M_(6,2,2)^2 = kappa_PS g_10^2 v_PS^2,
M_LQ^2 = (2/3) g_4^2 v_4^2,
M_WR^2 = g_R^2 v_R^2,
M_neutral^2 = v_Y^2/4 [[g_R^2, -g_R g_BL],[-g_R g_BL, g_BL^2]].
```

The demo diagonalization has 33 positive broken eigenvalues and one
hypercharge zero mode.

P11 replays the P6/P7 matching gates with the staged `Spin(10)` mass blocks:

```text
code/spin10_mass_replay_matching_gate.py
output/spin10_mass_replay_matching_gate.json
output/spin10_mass_replay_matching_gate.md
```

It supplies the mass replacement dictionary

```text
M_(6,2,2) -> sqrt(kappa_PS) g_10 v_PS
M_LQ      -> sqrt(2/3) g_4 v_4
M_WR      -> g_R v_R
M_Zprime  -> 0.5 v_Y sqrt(g_R^2+g_BL^2)
```

and then applies an adjoint-current gate to the 18 P7 chiral-pair rows.  The
result is conservative:

```text
current-vector mass blocks available = 4,
P7 rows replayed = 18,
P7 rows with candidate adjoint mass blocks = 0,
P7 rows remaining scalar/source mass debt = 18,
physical proton bounds evaluable now = 0.
```

Thus P11 improves the vector-current side but does not unlock proton bounds for
the P7 chiral-pair pole ledger.  The next decision is whether proton matching
continues in a completed current-current vector basis or in a scalar/source
chiral-pair basis.

P12 records both post-P11 branch directions and activates the scalar/source
chiral-pair branch first:

```text
code/scalar_source_chiral_pair_branch.py
output/scalar_source_chiral_pair_branch.json
output/scalar_source_chiral_pair_branch.md
```

The two required future attempts are:

```text
Branch V: completed Spin(10) current-current vector exchange
Branch S: scalar/source chiral-pair poles
```

P12-S assigns scalar/source denominators to all 18 P7 chiral-pair rows:

```text
L ⊃ lambda_{ab,R} psi_a psi_b S_R + h.c. + M_{S_R}^2 |S_R|^2
C6(12|34;R)=lambda_{12,R} lambda^*_{34,R}/M_{S_R}^2
```

The generated scalar/source ledger has 12 source sectors, converts all 18 P7
rows, and keeps 6 baryon-violating or sterile-BNV rows conditional.  Physical
proton bounds remain unevaluable because scalar masses, flavor tensors, Fierz
and chiral contractions, RG factors, hadronic matrix elements, and experimental
channel limits are still missing.

The P1--P12 derivations are also consolidated in a standalone TeX ledger:

```text
tex/route_c_derivation_ledger.tex
```

This TeX file is the paper-facing derivation record for the current Route-C
bootstrap chain.  The generated `.md`/`.json` files remain the machine-readable
audit surface.

## Workplan

The full staged plan is maintained in
[`ROUTE_C_WORKPLAN.md`](ROUTE_C_WORKPLAN.md).  It orders the work from the
representation pre-ledger through pole factorization, Ward identities,
anomaly checks, high-energy behavior, low-energy matching, and final comparison
against `SU(5)`, Pati-Salam, `Spin(10)`, `E6`, and superstring-derived
completions.

## What counts as progress

Progress is not a slogan that bootstrap "derives the GUT."  Progress means one
of the following concrete outcomes:

1. a candidate group or action branch is ruled out by consistency;
2. several branches remain, with a precise list of missing assumptions;
3. `Spin(10)` is preferred only after clearly stated minimality/no-exotics
   assumptions;
4. a measurable low-energy operator or bound is derived from the surviving
   amplitude data.
