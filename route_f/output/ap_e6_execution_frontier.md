# AP-E6 execution-frontier gate

- Mechanical checks: **15/15**
- Sp(4) eta/threshold lane closed: **False**
- Relaxed B=1/full-Hessian lane closed: **False**
- Same-soliton Callias/CPT/descent lane closed: **False**
- Canonical same-background SHA-256: **81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59**
- Any pre-portal route closed: **False**
- Portal start authorized: **False**
- Portal constructed: **False**
- Physics promotion allowed: **False**

## Decision

No complete pre-portal lane is closed.  In accordance with the ordered research rule, the degree-one Route-E portal is not started.

## Recomputed closure inputs

### sp4_eta_threshold

- `sp4_euclidean_regulator_complete`: `false`
- `omega5_spin_bsp4_computed`: `true`
- `target_mapping_torus_eta_pair_selected`: `false`
- `heavy_threshold_eta_matched`: `false`
- `threshold_radiatively_protected`: `false`
- `strong_b1_phase_proven`: `false`

### relaxed_b1_full_hessian

- `same_relaxed_b1_solution_computed`: `true`
- `multigrid_multivolume_converged`: `false`
- `sparse_projected_hessian_complete_in_declared_sector`: `false`
- `all_blocks_assembled_on_one_same_grid_and_boundary_condition`: `false`
- `discrete_stationarity_achieved`: `false`
- `lattice_topology_converged_to_B1`: `false`
- `physical_projected_hessian_at_stationary_discrete_solution`: `false`
- `projected_hessian_stability_gate_pass`: `false`
- `full_gauge_meson_ghost_fermion_hessian_complete`: `false`
- `dynamical_4d_importance_sampling_performed`: `false`
- `continuum_quantum_phase_proven`: `false`

### same_soliton_callias_descent

- `same_physical_moduli_space_derived`: `false`
- `actual_yukawa_callias_operator_derived`: `false`
- `uniform_fredholm_gap_proven`: `false`
- `determinant_line_o2_derived`: `false`
- `physical_cpt_regulator_derived`: `false`
- `gauge_basic_wzw_descent_proven`: `false`
- `same_soliton_composition_closed`: `false`

## Checks

- [PASS] `provenance` - all AP-E6 source sets exist and are hashed: hashed=10/10
- [PASS] `source_cards` - sp4 mechanical card is green: checks=39/39
- [PASS] `source_cards` - sp4 card forbids premature physics promotion: physics_promotion_allowed=False
- [PASS] `provenance` - sp4 source-card manifest is fresh: fresh=3/3; stale=[]
- [PASS] `source_cards` - b1 mechanical card is green: checks=20/20
- [PASS] `source_cards` - b1 card forbids premature physics promotion: physics_promotion_allowed=False
- [PASS] `provenance` - b1 source-card manifest is fresh: fresh=3/3; stale=[]
- [PASS] `source_cards` - callias mechanical card is green: checks=25/25
- [PASS] `source_cards` - callias card forbids premature physics promotion: physics_promotion_allowed=False
- [PASS] `provenance` - callias source-card manifest is fresh: fresh=5/5; stale=[]
- [PASS] `same_background` - Callias/CPT/WZW card consumes the canonical relaxed B=1 profile: B1=81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59; Callias=81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59
- [PASS] `closure_logic` - sp4 source lane flag equals the independent conjunction: source=False; recomputed=False; inputs={'sp4_euclidean_regulator_complete': False, 'omega5_spin_bsp4_computed': True, 'target_mapping_torus_eta_pair_selected': False, 'heavy_threshold_eta_matched': False, 'threshold_radiatively_protected': False, 'strong_b1_phase_proven': False}
- [PASS] `closure_logic` - b1 source lane flag equals the independent conjunction: source=False; recomputed=False; inputs={'same_relaxed_b1_solution_computed': True, 'multigrid_multivolume_converged': False, 'sparse_projected_hessian_complete_in_declared_sector': False, 'all_blocks_assembled_on_one_same_grid_and_boundary_condition': False, 'discrete_stationarity_achieved': False, 'lattice_topology_converged_to_B1': False, 'physical_projected_hessian_at_stationary_discrete_solution': False, 'projected_hessian_stability_gate_pass': False, 'full_gauge_meson_ghost_fermion_hessian_complete': False, 'dynamical_4d_importance_sampling_performed': False, 'continuum_quantum_phase_proven': False}
- [PASS] `closure_logic` - callias source lane flag equals the independent conjunction: source=False; recomputed=False; inputs={'same_physical_moduli_space_derived': False, 'actual_yukawa_callias_operator_derived': False, 'uniform_fredholm_gap_proven': False, 'determinant_line_o2_derived': False, 'physical_cpt_regulator_derived': False, 'gauge_basic_wzw_descent_proven': False, 'same_soliton_composition_closed': False}
- [PASS] `portal` - portal construction obeys the one-complete-lane prerequisite: any_closed=False; constructed=False
