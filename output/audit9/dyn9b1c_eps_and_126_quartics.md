# DYN-9b-1c: epsilon-coupling scan + 126bar mixed quartics

11/11 checks pass.

Key numbers in the JSON derivation_log: S2_eps_scan (the LR rescue verdict), S3_126_gates (embedding anchors), S4_K4_retest (the Hyper-K window with the 126bar remnants included).

## Checks

- [PASS] FULL epsilon gradient (three direct placements + the dual-slot adjoint term) verified against the directional derivative
- [PASS] epsilon Hessian certified: the homogeneity-identity build matches the direct second difference of Q6 values on random directions
- [PASS] re-verified: the epsilon quartic vanishes identically on the singlet slice (all 19 sample points), so lambda_eps enters NO stationarity condition -- it is a pure off-slice mass lever
- [PASS] GAUGE GATE WITH lambda_eps != 0: exactly 24 / 30 Goldstone zeros persist -- certifies the epsilon Hessian's gauge consistency (the strongest available correctness check)
- [PASS] LR RESCUE VERDICT: fraction of random couplings for which SOME lambda_eps in [-1.5, 1.5] makes the left-right (p, a) vacuum a tree-level local minimum
- [PASS] PS stability under the epsilon coupling recorded (samples whose tree positivity flips when lambda_eps = +-1)
- [PASS] MIXED CUBIC GATE: the Phi Sigma Sigma* invariant descends to |sigma|^2 (p + 3a - 6w) -- the AG eta-term structure, extending the embedding anchor to the 126bar sector
- [PASS] FULL VECTOR-MASS GATE with sigma: 12 zeros (SM) + all 33 massive bosons match the complete DYN-1b formulas G/J/F/E/X including the sigma terms (one overall scale from X, one sigma-normalization from G, J/F/E cross-checked)
- [PASS] duality restriction consistent: the sigma direction lies ENTIRELY in the restricted 126bar half (norm 1), and the cross-half mixing of the PS-invariant kernels (allowed in the shared (6,1,1)/(15,2,2) channels; those components are not fields) is recorded
- [PASS] Sigma-sector structure on the 126bar half: the mixed cubic + quartics split it into exactly the PS multiplets (6,1,1) + (10,1,3) + (10bar,3,1) + (15,2,2) [complex dims 6/30/30/60], with the sigma direction MACHINE-LOCATED in the SU(2)_R-charged 30 (the Delta_R block ESH tunes to M_I); the D-parity-odd mixed CUBIC is what splits the two 30s
- [PASS] K4 RETEST including the 126bar remnants: over random mixed couplings (cubic eta + four quartics) with the extended-survival tuning (Delta_R block at M_I) and positive Sigma remnants, the tau distribution is recorded
