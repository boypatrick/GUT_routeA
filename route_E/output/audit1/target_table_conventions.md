# Audit 1 Flavor / Seesaw Target-Table Convention

Audit 1 has been opened as a convention scaffold, not a global fit.

## Digest

- card sha256: `fa582e4c56abcc5f04e167fd28cac75f684d9d5417aab1315d5284b4a7fd19f9`
- audit0 sha256: `a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc`
- source artifacts: 5/5 present

## Fit Observable Contract

- 19 fit observables: 6 quark Yukawas, 3 charged-lepton Yukawas, 4 CKM quantities, 6 PMNS/neutrino quantities.
- Pure predictions: lightest mass, two Majorana phases, m_beta_beta, and heavy-neutrino spectrum.
- Publication targets are intentionally marked `REQUIRES_PUBLICATION_REFRESH`.

## Target Table Version Policy

- `target_table_v1_regression`: frozen local regression anchor; `zeta_role = fixed_benchmark`.
- `target_table_v2_publication`: refreshed publication target table; `zeta_role = free_parameter`.
- Every future fit output must state both `target_table_version` and `zeta_role`.

## Zeta Phase Policy

- Audit 0.5 found no hit in the finite visible-spurion phase-transfer table.
- Publication fits should treat `arg(zeta)` as an independent hidden CP parameter unless a later UV mechanism predicts it before seeing the fit.
- Tree-level locking to CKM/PMNS phases or to the visible invariants `I,J` is forbidden by default.

## Local Anchors

- legacy `V_us = 0.22499845986972888`
- legacy `V_cb = 0.04099971935403949`
- legacy `V_ub = 0.0037`
- legacy `J_CKM = 3.097061593767482e-05`
- local seesaw `sin2_theta12 = 0.304`
- local seesaw `sin2_theta23 = 0.573`
- local seesaw `sin2_theta13 = 0.0222`

## Boundary

This file does not claim a completed flavor fit. It only fixes the target-table contract, legacy local anchors, basis conventions, and downstream outputs required by Audit 2 and Audit 3.
