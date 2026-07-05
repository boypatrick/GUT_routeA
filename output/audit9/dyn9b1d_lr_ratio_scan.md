# DYN-9b-1d (SUB-B0): the left-right vev ratio freed

7/7 checks pass.

- Q1 (He-Meljanac reproduction): positive points 8/8800 at lambda_eps = 0.
- Q2 (the epsilon question): rescued 0, killed 0, positive-at-some-eps 8.
- Conclusion: 9b-1b's negative left-right verdict WAS a fixed-ratio artifact: freeing the ratio reproduces stable left-right tree minima in finite coupling regions, consistent with He-Meljanac (1986) at claim level.  K5r REVERTS: the left-right chain is 210-realizable at tree level in those regions; the 45_H route is an alternative, not a necessity.

## Checks

- [PASS] POLARIZATION GATE: the assembled Hessian equals a direct build at r = 0.8 for every polynomial pattern (the epsilon pattern is certified by the Goldstone gate below)
- [PASS] FULL-STATIONARITY GATE: with (m2, kappa) solved on the slice, the total gradient vanishes in ALL 210 directions (little-group invariance argument verified numerically; epsilon contributes no gradient at any slice point)
- [PASS] GOLDSTONE GATE at four ratios with random couplings AND random epsilon: exactly 30 zeros every time -- certifies the polarization assembly including the epsilon block
- [PASS] positivity scan executed over 16 ratios x 300 coupling points x 13 epsilon values
- [PASS] Q1 VERDICT (the He-Meljanac reproduction question): does the left-right stationary family admit tree-level LOCAL MINIMA somewhere in (ratio x couplings) at lambda_eps = 0 (the classic four-invariant potential family)?
- [PASS] Q2 VERDICT (the epsilon question): rescue and kill counts under the fifth invariant recorded
- [PASS] K5r PROPAGATION statement recorded for the ledger refresh and both papers
