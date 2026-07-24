# AP-E16 tangent/Caccioppoli/Young-measure audit

- Status: **PASS**
- Checks: **12/12**
- Simultaneous field scaling: `alpha=1.000000000000, d=3.000000000000`
- Homogeneous Young-measure gap: `0.625000000000`
- Weighted-cutoff slope: `2.006692368`

## Main conclusions

- At `mu`-a.e. points, normalized tangent measures and their Hopf Radon-Nikodym components exist; a diagonal recovery sequence can be chosen with vanishing normalized compact-phase drop.
- One ordinary field blow-up retains all three minor orders only in the volumetric branch `mu(B_r)~r^3`, with target amplitude `s~r`.
- The normal-ball Hopf-base Caccioppoli inequality is proved with explicit weighted and local-comparison-deficit terms.
- A phase-purified homogeneous pure-base Young measure has positive defect while its weighted cutoff tends to zero.  It refutes unconditional excess contraction, but is not an asymptotically minimizing sequence.

## Gates

- `tangent_measure_blowup`: `True`
- `hopf_defect_radon_nikodym_split`: `True`
- `normalized_phase_stationarity`: `True`
- `single_field_all_order_tangent_general`: `False`
- `single_field_volumetric_branch`: `True`
- `conditional_hopf_base_caccioppoli`: `True`
- `unconditional_excess_contraction`: `False`
- `homogeneous_base_young_measure_counterexample`: `True`
- `counterexample_is_asymptotically_minimizing`: `False`
- `two_limit_weighted_decay_sufficient_alone`: `False`
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

- **scaling / unique simultaneous complete-minor scaling**: `PASS` — the three equations give target amplitude s~r and mu(B_r)~r^3
- **scaling / nonvolumetric branches require order-specific tangents**: `PASS` — a first-order field normalization does not retain all higher orders for d!=3
- **cutoff / rank-one second-minor expansion**: `PASS` — maximum residual 1.421e-14
- **cutoff / rank-one third-minor expansion**: `PASS` — maximum residual 2.132e-14
- **young measure / zero barycentric gradient**: `PASS` — means [-3.0839528461809905e-17, -7.37257477290143e-17]
- **young measure / zero barycentric second minor**: `PASS` — mean minor 6.264e-19
- **young measure / positive exact Jensen gap**: `PASS` — gap 0.625000000000
- **sphere generator / exact target constraint**: `PASS` — all sampled generators lie on S3 to machine precision
- **sphere generator / pure-base leading limits**: `PASS` — first/second energy persist while the Hopf connection vanishes
- **sphere generator / weighted cutoff vanishes with positive defect**: `PASS` — weighted slope 2.006692368
- **analysis / conditional Caccioppoli versus unconditional contraction**: `PASS` — Caccioppoli needs normalized full quasiminimality for excess contraction
- **policy / no new lattice scan and downstream embargo**: `PASS` — forbidden scan calls found: []
