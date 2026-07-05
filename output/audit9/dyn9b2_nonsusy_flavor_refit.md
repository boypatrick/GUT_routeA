# DYN-9b-2: non-SUSY flavor refit + zeta re-extraction

11/11 checks pass.

## Verdicts

1. **SO(10) lock tension**: y_t(M_X) = 0.416 (machine-derived); the perturbative ceiling caps the nu-Dirac scale at y_max = 0.0458 (PS) / 0.0026 (G_LR) -- tension factors 9x / 159x.
2. **Type-II fails** by 3.7 / 7.3 orders (PS / G_LR).
3. **Zeta-invariance theorem**: zeta, contact direction and contact fraction are EXACTLY invariant under the Dirac rescale; only M_* moves (M_*' = 4 pi M_I at the ceiling: 1.04e+13 PS, 3.35e+10 G_LR).
4. **Leptogenesis selects the scale-decoupled source**: rescaled M_1 falls below the Davidson-Ibarra floor on both chains; the archival M_1 survives only under an instanton-type tower -- the D3 coexistence prediction (K8) becomes a REQUIREMENT of the alive branch with viable leptogenesis (D3 stays unpromoted).

## Boundary (NOT claimed)

- No full non-SUSY Yukawa fit; the Dirac kernel shape is a retained conditional input.
- The SO(10) lock is quantified from the third-generation relation, not re-derived from a global fit.
- Type-II and Davidson-Ibarra numbers are order estimates.
- zeta's value is NOT derived (boundary theorem intact).

## Checks

- [PASS] archival closure card sha256 matches the DYN-4a provenance; the DYN-9 ledger declares the card slice-local (the premise this audit sharpens)
- [PASS] archival replay reproduces the paper anchors (zeta, M_*, tower)
- [PASS] y_top(M_X) machine-derived from the one-loop SM running: the SO(10) lock scale for the third-generation nu-Dirac coupling
- [PASS] archival Dirac kernel scale recorded: sigma_max(Y_nu) is O(y_t) -- the archival card RESPECTS the SO(10) lock
- [PASS] SO(10) LOCK TENSION quantified: the perturbative ceiling 4 pi M_I caps the nu-Dirac scale at y_max; the renormalizable minimal source needs the third-generation coupling ~9x (PS) / ~160x (G_LR) BELOW its locked scale y_t(M_X) -- i.e. 2 / 4.4 orders of magnitude in the tower, a structural tension only removable by that much Yukawa tuning
- [PASS] TYPE-II ESCAPE FAILS (order estimate, flagged; lambda ~ 1 and f at 4 pi, i.e. GENEROUS): with the Delta_L-type block at the gauge scale (DYN-9b-1c) the induced triplet vev falls short of the atmospheric scale by 3.7 (PS) / 7.3 (G_LR) orders
- [PASS] ZETA-INVARIANCE THEOREM (machine-verified): under Y_nu -> y Y_nu the tower rescales uniformly, so zeta, the contact direction and the contact fraction are EXACTLY invariant; only M_* = y^2 M_*_arch moves.  The contact card's normalized content is scale-covariant; M_* is the only branch-local number
- [PASS] re-extracted M_* windows at v_R ~ M_I recorded per chain (renormalizable-source scenario at the perturbative ceiling): M_*' = 4 pi M_I by construction
- [PASS] contact essentiality carries over VERBATIM to the rescaled card (cf invariant; the DYN-4b posterior statement P(cf > 0.01) = 1.000 remains conditional on the Dirac kernel SHAPE only, which the rescale does not touch)
- [PASS] LEPTOGENESIS SOURCE SELECTION: under the rescaled renormalizable source M_1 falls BELOW the Davidson-Ibarra floor on both chains (PS ~ 7e7, G_LR ~ 2e5 GeV); under a scale-decoupled (instanton-type) source the archival M_1 ~ 2.8e10 survives -- thermal leptogenesis SELECTS the scale-decoupled source on the alive branch (DI floor is an order estimate, flagged)
- [PASS] K8 UPGRADE: the only quantified source scenario that survives S1 (lock tension) + S2 (type-II fails) + S4 (leptogenesis) is the scale-decoupled instanton-type tower of the D3 card, whose coexistence prediction (N_2, N_3 above the intermediate gauge scale) therefore upgrades from an optional pricing card to a REQUIREMENT of the alive branch with viable leptogenesis
