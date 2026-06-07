# Route-C P8 Comparative Theory Scorecard

P8 compares candidate branches using the explicit P1--P7 gates.  It is a
hard-gate ledger, not a weighted score and not a completed GUT proof.

## Summary

| quantity | value |
| --- | ---: |
| candidates compared | 5 |
| leading conditional field-theory branch | Spin(10) |
| tower-like UV forced now | False |

Final status counts:

| final status | count |
| --- | ---: |
| conditional_branch_if_single_object_filter_is_relaxed | 1 |
| conditional_intermediate_face_not_final_simple_object | 1 |
| leading_conditional_route_c_branch | 1 |
| conditional_branch_with_exotics_lifting_debt | 1 |
| uv_completion_class_not_forced_by_p1_p7 | 1 |

## Candidate Gate Matrix

| candidate | C0 | P5 anomaly/chiral | P6 high energy | P7 proton readiness | final status |
| --- | --- | --- | --- | --- | --- |
| SU(5) | fail_or_relax | pass | completion_required | not_evaluable_yet | conditional_branch_if_single_object_filter_is_relaxed |
| Pati-Salam | fail_or_relax | pass | completion_required | not_evaluable_yet | conditional_intermediate_face_not_final_simple_object |
| Spin(10) | pass | pass | completion_required | not_evaluable_yet | leading_conditional_route_c_branch |
| E6 | fail_or_relax | conditional | completion_required | not_evaluable_yet | conditional_branch_with_exotics_lifting_debt |
| Superstring-derived completions | not_applicable | conditional | completion_required | not_evaluable_yet | uv_completion_class_not_forced_by_p1_p7 |

## Candidate Notes

### SU(5)

- family realization: `10 + bar5 + 1`
- final status: `conditional_branch_if_single_object_filter_is_relaxed`
- note: Anomaly-consistent and simple, but the family is split as 10 + bar5 + 1 rather than one irreducible six-face object.
- next required audits:
  - explicit statement that Route-A single-object criterion is relaxed
  - SU(5)-specific broken-vector completion
  - classic X/Y proton matching and thresholds

### Pati-Salam

- family realization: `(4,2,1) + (bar4,1,2)`
- final status: `conditional_intermediate_face_not_final_simple_object`
- note: Anomaly-consistent and physically natural, but product-group rather than simple single-object unification.
- next required audits:
  - explicit product-group bootstrap branch
  - Pati-Salam leptoquark and SU(2)_R broken Ward completion
  - proton and threshold audit if used as more than an intermediate face

### Spin(10)

- family realization: `16 half-spinor`
- final status: `leading_conditional_route_c_branch`
- note: Only field-theory candidate passing the current minimal single-object filter, but still blocked from physical high-energy/proton claims by P6/P7 completion debt.
- next required audits:
  - explicit broken D5 generator matrices
  - Spin(10)-breaking Higgs/Goldstone/source sector
  - mediator masses and couplings
  - flavor-basis Wilson tensors and proton bounds

### E6

- family realization: `27 -> 16 + 10 + 1 under Spin(10)`
- final status: `conditional_branch_with_exotics_lifting_debt`
- note: Simple and contains a Spin(10) 16 inside 27, but the extra 10 + 1 states must be lifted or audited.
- next required audits:
  - mass/lifting mechanism for 10 + 1
  - check no light chiral exotics remain
  - E6 threshold/proton audit including lifted states

### Superstring-derived completions

- family realization: `model-dependent; may realize SU(5), Pati-Salam, Spin(10), E6, or other quiver/brane/heterotic branches`
- final status: `uv_completion_class_not_forced_by_p1_p7`
- note: Not a single finite four-dimensional GUT candidate.  Relevant as a comparison class if finite-pole GUT completions fail high-energy softness or require a tower.
- next required audits:
  - explicit compactification and massless spectrum
  - Green-Schwarz/tadpole/modular consistency
  - chiral index, exotics lifting, and matching to the Route-A skeleton

## P8 Decision Boundary

Spin(10) is the leading conditional finite field-theory branch under the current minimal single-object filter.  This is not yet a full action-level or phenomenological proof.  Superstring-derived completions remain a UV comparison class, not a P8 conclusion.

P8 is a hard-gate scorecard over existing Route-C ledgers.  It does not upgrade conditional branches into completed GUTs.
