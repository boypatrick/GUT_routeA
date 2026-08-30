# AP-E18 GTA minimizing-diagonal audit

- Status: **PASS**
- Checks: **14/14**
- Maximum exact-scaling residual: `1.088e-12`
- Universal GTA from the AP-E17 additive diagonal: **false**
- Scale-polished/exact-local-minimizer GTA: **open**

## Main conclusions

- Standard reverse-Hölder import fails both the `(p,q)` exponent and nonlinear-graph hypotheses.
- A boundary-straddling rank-two needle has strong tangent `L2` convergence and only `O(epsilon)` normalized graph energy.
- Its first- and second-minor energy-amplitude correlations remain order one, so every shrinking GTA thickness fails.
- The physical perturbation costs `O(rho^4)=o(rho^3)`, preserves degree, and is invisible to the tangent defect measure.
- Additive global near-minimality is therefore too weak at unresolved microscales; a uniformly scale-aware local gauge is the next gate.

## Gates

- `GTA_implied_by_AP_E17_additive_diagonal`: `False`
- `universal_reverse_holder_from_current_hypotheses`: `False`
- `rank_two_boundary_needle_counterexample`: `True`
- `counterexample_changes_tangent_defect_measure`: `False`
- `counterexample_changes_degree`: `False`
- `existence_of_a_scale_polished_GTA_diagonal`: `False`
- `required_new_gate_uniform_microscale_quasiminimality`: `True`
- `full_mu_zero`: `False`
- `local_recovery`: `False`
- `classicality`: `False`
- `continuum_isolation`: `False`
- `bosonic_hessian_authorized`: `False`
- `parallel_dirac_callias_mathematics_allowed`: `True`
- `determinant_promotion`: `False`
- `degree_one_portal_promotion`: `False`
- `physics_promotion_allowed`: `False`

## Check details

- **needle / exact scaling exponents**: `PASS` — maximum log-slope residual 1.088e-12
- **needle / strong tangent L2 convergence**: `PASS` — the primary L2 mass is epsilon^5
- **needle / normalized graph-energy perturbation is o(1)**: `PASS` — both leading first- and second-minor energies are O(epsilon)
- **needle / energy-amplitude correlation does not vanish**: `PASS` — both leading weighted correlations are order one
- **geometry / the perturbation is genuinely rank two**: `PASS` — rank=2, one second minor is nonzero, and all pure third minors vanish
- **geometry / sphere-chart amplitude remains in one hemisphere**: `PASS` — rho_epsilon A_epsilon=epsilon^(1/2)->0
- **minimizing diagonal / physical added energy is o(rho^3)**: `PASS` — added energy is O(rho^4), hence normalized additive deficit is O(rho)
- **GTA / annular graph energy still vanishes**: `PASS` — both delta>=h and delta<h branches tend to zero
- **GTA / every boundary-straddling thickness has divergent weighted quotient**: `PASS` — delta>=h: Q>=c/delta^2; delta<h: Q>=c/(epsilon^2 delta)
- **reverse Holder / no positive fixed-ball integrability gain**: `PASS` — the ratio diverges for every tested sigma>0 with the exact negative exponent
- **literature / standard (p,q) theorem is outside range**: `PASS` — n=3,p=2 gives q<3, whereas the ambient upper growth is q=6
- **literature / A-free theorem is not directly transferable**: `PASS` — curl/div constraints do not remove nonlinear Plucker and S3 restrictions
- **gates / fail-closed decision**: `PASS` — GTA is refuted for arbitrary additive diagonals; scale-polished recovery is open
- **policy / no lattice relaxation or mesh scan**: `PASS` — forbidden execution calls found: []
