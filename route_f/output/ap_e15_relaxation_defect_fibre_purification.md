# AP-E15 relaxation-defect and fibre-purification audit

- Status: **PASS**
- Checks: **11/11**
- Phase Hessian maximum residual: `8.423e-07`
- Strong-convexity critical radius: `9.068996821`
- Weighted-cutoff affine slope: `3.000000000000`

## Theorem ledger

- `mu` is a nonnegative Radon relaxation-defect measure.
- `mu=0` iff the selected recovery subsequence converges strongly in the complete-minor graph norm.
- Every compact-fibre-phase-reducible vertical defect vanishes for an asymptotically minimizing smooth fixed-degree sequence.
- Weighted cutoff decay holds at Lebesgue-a.e. centres for each fixed map, but the sequence-uniform `lim_r limsup_j` estimate is open.

## Gates

- `ap_e11_action_frozen`: `True`
- `new_lattice_scan`: `False`
- `relaxation_defect_measure_defined`: `True`
- `compact_fibre_phase_strong_convexity`: `True`
- `phase_reducible_vertical_defect_zero`: `True`
- `mixed_hopf_base_defect_zero`: `False`
- `weighted_cutoff_decay_lebesgue_a.e.`: `True`
- `weighted_cutoff_decay_on_defect_support`: `False`
- `sequence_uniform_weighted_cutoff_mu_a.e.`: `False`
- `full_mu_zero`: `False`
- `local_recovery`: `False`
- `regularity`: `False`
- `continuum_isolation`: `False`
- `bosonic_hessian_authorized`: `False`
- `parallel_dirac_callias_mathematics_allowed`: `True`
- `determinant_promotion`: `False`
- `degree_one_portal_promotion`: `False`
- `physics_promotion_allowed`: `False`

## Check details

- **algebra / phase Hessian coefficients**: `PASS` — maximum five-point finite-difference residual 8.423e-07
- **analysis / production convexity radius**: `PASS` — lambda_r=1-m^2 r^2/pi^2 and m^2=0.12
- **analysis / small-ball strong convexity**: `PASS` — r=0.25 and r=9 lie below the exact critical radius
- **defect / Hilbert polarization identity**: `PASS` — residual 0.000e+00
- **defect / nonnegative toy defect masses**: `PASS` — quadratic norm gaps are nonnegative
- **purification / degree-preserving local drop bound**: `PASS` — Delta_j(Q)<=E(u_j)-I_B -> 0 for every fixed admissible ball
- **cutoff / smooth weighted-cutoff exponent**: `PASS` — fitted affine exponent 3.000000000000
- **cutoff / a.e. versus singular-support scope**: `PASS` — fixed-map Lebesgue differentiation does not exchange r and recovery limits
- **policy / no new lattice scan**: `PASS` — forbidden scan calls found: []
- **gates / classicality chain remains closed**: `PASS` — full mu=0, recovery, regularity, isolation, and Hessian remain false
- **gates / determinant and portal embargo**: `PASS` — parallel operator mathematics is allowed without downstream promotion
