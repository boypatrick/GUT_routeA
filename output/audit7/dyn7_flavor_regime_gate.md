# DYN-7 / DYN-9b-3 flavor-regime gate

- arithmetic checks: **7/7 pass**
- benchmark `M1 = 2.380448e+10 GeV`
- non-SUSY classifier: **two_flavor_tau_resolved**
- DYN-7 status: **blocked_missing_branch_thermal_inputs**
- DYN-9b-3 status: **blocked_requires_two_flavor**
- historical success probabilities and conditioned `arg(zeta)` windows are diagnostics only

## Davidson-Ibarra correction

- `m3 = 5.017967716118e-02 eV` (one square root)
- SUSY, `v_u=100 GeV`: `1.397415940654e-05`
- SM, `v=174 GeV`: `2.307794855090e-06`
- old double-square-root overestimate: `3.911889`

## Checks

- [PASS] DYN-4a exports a positive hierarchical heavy-neutrino spectrum: M1=2.380448e+10 GeV
- [PASS] the non-SUSY benchmark lies in the tau-resolved two-flavor regime: classifier=two_flavor_tau_resolved
- [PASS] MSSM flavor certification fails closed when branch-local thermal inputs are absent: missing=tan_beta,charged_lepton_thermal_rates,spectator_matrix
- [PASS] light masses are extracted with one square root, not two: m3=5.017967716118e-02 eV
- [PASS] corrected SUSY Davidson-Ibarra central bound is reproduced: epsilon_DI_SUSY=1.397415940654e-05
- [PASS] historical double-square-root bound was overestimated by about 3.91: factor=3.91188886
- [PASS] neither historical unflavored scan is promoted to a physical posterior
