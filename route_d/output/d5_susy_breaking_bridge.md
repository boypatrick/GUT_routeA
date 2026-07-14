# Route-D D5: high-scale SUSY-breaking bridge scan

18/18 checks pass.  CONDITIONAL scenario, NOT promoted; zeta NOT derived.

## Verdicts

- G_LR bridge window: log10 M_SS in [10.30, 15.55], binding: ['no_physical_solution'] / ['M_SS >= M_*'], cases ['A(MI<=MSS)']
- PS bridge window: EMPTY
- PS d=6 rescue anywhere on the bridge: False
- SUSY-segment content doubling (Delta-bar, Sigma-bar) forced by the U(1)_{B-L} cubic anomaly: b2R jumps -7/3 -> +5 (LR) and 11/3 -> +41 (PS).
- Endpoint gates: DYN-9 ledger reproduced at the high end; DYN-3 kill number 1.2e26 reproduced at (1.45e13, 3 TeV); TeV end d=5 dead on both chains.
- Every alive point has M_SS < M_*.  This is only a necessary matching-scale ordering and does not restore the invalid DYN-5 messenger action.

## Boundary (NOT claimed)

- One-loop intermediate segments; ESH-minimal content; M_T = M_X; degenerate soft spectrum; above-M_X perturbativity not audited.
- Route-D promotion bar NOT passed; must not be cited as evidence.

## Checks

- [PASS] non-SUSY intermediate b's reused as gated in DYN-9: G_LR (-7,-3,-7/3,11/2), PS (-23/3,-3,11/3)
- [PASS] SUSY G_LR b's DERIVED: (b3, b2L, b2R, bBL) = (-3, 1, 5, 15)
- [PASS] SUSY PS b's DERIVED: (b4, b2L, b2R) = (12, 1, 41)
- [PASS] MSSM segment b's are the textbook (33/5, 1, -3)
- [PASS] anomaly gate: gauged U(1)_{B-L} cubic anomaly of the DeltaR superfield fermions is (+2)^3 x 3 = 24 != 0, so the SUSY chain REQUIRES the conjugate DeltaR-bar (and SUSY PS likewise doubles SigmaR) -- a content doubling relative to the non-SUSY ESH set, which is what makes the SUSY-segment b's so much larger
- [PASS] dense two-loop SM trajectory reproduces the independent 160-step RK4 oracle at 1e5, 1e10, and 1e15 GeV
- [PASS] ENDPOINT GATE (DYN-9 limit), G_LR: M_SS above M_X reproduces the DYN-9 ledger (log10 M_I, log10 M_X, alpha_G^-1)
- [PASS] ENDPOINT GATE (DYN-9 limit), PS: M_SS above M_X reproduces the DYN-9 ledger (log10 M_I, log10 M_X, alpha_G^-1)
- [PASS] DYN-3 CALIBRATION GATE: at the DYN-2 compatibility point (M_T_eff = 1.45e13, M_SS = 3 TeV) the mini-split curve returns the DYN-3 kill number tau ~ 1.2e26 yr
- [PASS] bridge scan bookkeeping: every grid point either yields a physical solution or is DISCLOSED as a no-physical-solution obstruction region (multi-start Newton; obstruction = part of the map, as in DYN-2b)
- [PASS] TeV-END GATE: at M_SS ~ 3 TeV every bridge is EXCLUDED -- physical-and-d=5-dead or no physical solution (the DYN-2/3 exclusion structure; the catastrophic DYN-3 magnitude was slice-specific, disclosed)
- [PASS] G_LR bridge window verdict computed
- [PASS] PS bridge window verdict computed (incl. whether the SUSY segment RESCUES the marginal chain past the 2.4e34 bound)
- [PASS] PS d=6 rescue question answered explicitly across the whole bridge
- [PASS] PS sampled-segment diagnostic: with the declared doubled-Sigma content (b2R=41), no physical grid point contains more than 0.05 dex of SUSY PS running below M_X.  This records the result but does not by itself prove b2R=41 is the unique cause
- [PASS] necessary matching-scale ordering intersected: every alive point has M_SS < M_* = 3.93e15 GeV; this is not sufficient to validate the Route-B/DYN-5 messenger action rejected by DYN-5V
- [PASS] M_T = M_X assumption FLAGGED with product scaling: tau_d5 scales as (M_T M_SS)^2, so each decade of triplet lightness shifts the d=5 window edge up one decade in M_SS
- [PASS] BRIDGE VERDICT: the high-scale SUSY-breaking bridge has a non-empty survival window on at least one chain -- the string-natural interpolation benchmark exists on the sampled G_LR grid
