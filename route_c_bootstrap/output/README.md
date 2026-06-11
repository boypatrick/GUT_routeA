# Route-C Output Directory

Generated Route-C bootstrap ledgers and scratch outputs go here.

The first generated files are expected to be:

- `candidate_ledger.json`
- `candidate_ledger.md`
- `external_state_charges.json`
- `external_state_charges.md`
- `four_point_pole_ansatz.json`
- `four_point_pole_ansatz.md`
- `residue_positivity.json`
- `residue_positivity.md`
- `ward_identity_ledger.json`
- `ward_identity_ledger.md`
- `anomaly_chiral_consistency.json`
- `anomaly_chiral_consistency.md`
- `high_energy_growth_audit.json`
- `high_energy_growth_audit.md`
- `low_energy_matching_proton_report.json`
- `low_energy_matching_proton_report.md`
- `theory_comparison_scorecard.json`
- `theory_comparison_scorecard.md`
- `spin10_action_completion_ledger.json`
- `spin10_action_completion_ledger.md`
- `spin10_breaking_mass_matrix.json`
- `spin10_breaking_mass_matrix.md`
- `spin10_mass_replay_matching_gate.json`
- `spin10_mass_replay_matching_gate.md`
- `scalar_source_chiral_pair_branch.json`
- `scalar_source_chiral_pair_branch.md`
- `scalar_source_flavor_wilson_basis.json`
- `scalar_source_flavor_wilson_basis.md`
- `spin10_vector_branch_matching_gate.json`
- `spin10_vector_branch_matching_gate.md`
- `spin10_vector_all_left_fierz_gate.json`
- `spin10_vector_all_left_fierz_gate.md`
- `spin10_vector_physical_multiplet_assembly.json`
- `spin10_vector_physical_multiplet_assembly.md`

The paper-facing derivation consolidation for P1--P16 is not generated in this
directory.  It lives at:

```text
../tex/route_c_derivation_ledger.tex
```

P12-S is the first post-P11 branch activation.  It records both future
directions, activates the scalar/source chiral-pair branch first, converts all
18 P7 chiral-pair rows into symbolic scalar/source matching coefficients, and
leaves the completed current-current vector branch as a required future
attempt.

P13-S maps the six BNV or sterile-BNV scalar/source rows to a symbolic
all-left chiral Wilson basis.  It introduces flavor tensors, marks direct rows,
marks rows needing Fierz/color recoupling, and still keeps physical proton
bounds unevaluable.

P14-V returns to the current-current vector branch.  It builds 72 vector
current-pair matching rows from the P9 broken-generator maps and P10 staged
mass denominators, while preserving the P11 result that the P7 chiral-pair rows
remain scalar/source mass debt.  It is a matching gate, not a proton-lifetime
calculation.

P15-V fixes a two-component Fierz/crossing convention for Branch V and maps all
72 current-pair rows to crossed all-left skeletons.  The current sparse-map
pairing yields 72 `B_and_L_conserving_crossed_skeleton` rows and zero canonical
BNV proton rows, so physical proton bounds remain blocked.

P16-V assembles physical Hermitian vector multiplets from charge-conjugate
root-map pairs.  It generates 448 cross-root rows and finds 144 canonical vector
BNV candidates before flavor, all in the `M_(6,2,2)` off-face block:
72 `O_V_QQ_DbarNbar` and 72 `O_V_QQ_UbarEbar` rows.  Cross-root Clebsch phases,
Hermitian generator normalization, flavor rotations, RG factors, hadronic
matrix elements, and proton-limit inputs remain unresolved, so physical proton
bounds are still not evaluable.
