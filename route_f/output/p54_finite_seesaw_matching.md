# Finite seesaw matching: canonical limit and actual sequential subset

Checks: 24/24.

**Full P54 finite seesaw matching remains incomplete.** The canonical all-N-decoupled C5 formula is implemented only for CHN=C5_UV=0 and leading qH/M^2=0. Actual P54 input has nonzero CHN and type-II C5, so that restricted formula explicitly rejects it.

The universal gauge/quartic finite C5 term is now applied at each moving Takagi mass block while retaining the existing nonzero CHN in the trajectory. Finite CHN/pre-existing-C5 insertions, mixed removed/retained-N graphs and correlated finite parameter matching are still absent.

| Diagnostic case | Live finite thresholds | Relative C5 shift from tree thresholds |
|---|---:|---:|
| 0 | 3 | 0.0022111964 |
| 1 | 3 | 0.0022404106 |

These shifts use the old synthetic families in the actual P54 group/light dictionary. They are not fitted neutrino predictions. The matching-scale variations are diagnostics of an incomplete map, not calibrated theory errors.

The canonical full-C5 scale derivative and the partial-block universal scale derivatives reproduce the corresponding beta differences. Degenerate blocks and generic family transformations are checked.

Sources: [Zhang-Zhou, Eqs. 49/51/91](https://arxiv.org/html/2107.12133v3); [Ohlsson-Pernow, Eqs. 42/45](https://arxiv.org/html/2201.00840v2).

| Check | Residual | Pass |
|---|---:|:---:|
| actual_CHN_and_typeII_point_rejected_by_canonical_formula | 0.000e+00 | True |
| full_canonical_C5_matching_scale_derivative | 8.607e-13 | True |
| canonical_finite_C5_complex_symmetric | 0.000e+00 | True |
| canonical_finite_C5_full_family_covariance | 1.464e-15 | True |
| omitted_finite_legs_detected | 0.000e+00 | True |
| degenerate_Higgs_wave_constant | 0.000e+00 | True |
| degenerate_lepton_wave_constant | 0.000e+00 | True |
| degenerate_Majorana_block_O3_invariance | 3.418e-16 | True |
| nondegenerate_universal_block_additivity | 3.657e-18 | True |
| partial_block_gauge_quartic_RG_identity_1 | 9.392e-15 | True |
| partial_block_gauge_quartic_RG_identity_2 | 1.440e-14 | True |
| partial_block_gauge_quartic_RG_identity_3 | 2.299e-14 | True |
| case0_three_live_finite_events | 0.000e+00 | True |
| case0_all_event_roots | 8.882e-16 | True |
| case0_CHN_kept_nonzero_at_all_thresholds | 0.000e+00 | True |
| case0_no_event_claims_complete_matching | 0.000e+00 | True |
| case0_finite_result_differs_from_tree | 0.000e+00 | True |
| case0_finite_event_trajectory_covariance | 6.835e-16 | True |
| case1_three_live_finite_events | 0.000e+00 | True |
| case1_all_event_roots | 4.441e-16 | True |
| case1_CHN_kept_nonzero_at_all_thresholds | 0.000e+00 | True |
| case1_no_event_claims_complete_matching | 0.000e+00 | True |
| case1_finite_result_differs_from_tree | 0.000e+00 | True |
| case1_finite_event_trajectory_covariance | 3.651e-16 | True |
