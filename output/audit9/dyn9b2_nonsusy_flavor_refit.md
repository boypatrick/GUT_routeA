# DYN-9b-2: fixed-kernel non-SUSY rescaling audit

14/14 checks pass.

## Verdicts

1. **Fixed-kernel ceiling**: y_t(M_X) = 0.441 is an optional top-like reference, not an SO(10) theorem.  The perturbative ceiling caps the nu-Dirac scale at y_max = 0.0458 (PS) / 0.0026 (G_LR).  The actual archival-kernel suppressions are 19.5x / 342.2x; the additional top-like ansatz gives 9.6x / 169.2x.
2. **Type-II fails** by 3.7 / 7.3 orders (PS / G_LR).
3. **Positive-real rescaling theorem**: for uniform y>0, zeta, projective contact direction, and contact fraction are exactly invariant while M_* moves (M_*' = 4 pi M_I at the ceiling: 1.04e+13 PS, 3.35e+10 G_LR).  For complex y, zeta acquires the phase y^2/|y|^2.
4. **Conditional source comparison (not a selection)**: rescaled M_1 falls below the Davidson-Ibarra floor on both chains, while the archival M_1 lies above the adopted floor in the instanton-type benchmark.  K8 remains unpromoted because the comparison retains a fixed Dirac kernel and omits flavored thermal kinetics.

## Boundary (NOT claimed)

- No full non-SUSY Yukawa fit; the Dirac kernel shape is a retained conditional input.
- Y_nu=h-3f and Y_u=h+f do not force a top-like Y_nu; the top-like comparison is an additional benchmark ansatz.
- Exact complex-zeta invariance is restricted to uniform positive-real rescaling; complex rescaling is phase-covariant.
- Type-II and Davidson-Ibarra numbers are order estimates.
- zeta's value is NOT derived (boundary theorem intact).
- Physics promotion is blocked pending the global fit, threshold matching, flavored kinetics, and a valid UV source construction.

## Checks

- [PASS] archival closure card sha256 matches the DYN-4a provenance; the DYN-9 ledger declares the card slice-local (the premise this audit sharpens)
- [PASS] archival replay reproduces the paper anchors (zeta, M_*, tower)
- [PASS] y_top(M_X) machine-derived from the one-loop SM running for the optional top-like third-generation benchmark
- [PASS] algebraic counterexample: Y_nu=h-3f and Y_u=h+f do not force Y_nu ~ Y_u; h=3f gives Y_nu=0 and Y_u=4f
- [PASS] archival Dirac kernel scale recorded: sigma_max(Y_nu) happens to be top-like within an order-one factor in this card; this is a benchmark property, not a consequence of the SO(10) relations
- [PASS] FIXED-ARCHIVAL-KERNEL diagnostic: the perturbative 4 pi M_I ceiling requires a uniform norm suppression y_arch/y_max of about 19 (PS) / 342 (G_LR)
- [PASS] OPTIONAL TOP-LIKE ANSATZ diagnostic: if one additionally imposes y_nu,3 ~ y_t(M_X), the same ceilings correspond to tensions of about 10 (PS) / 169 (G_LR); this is not an SO(10) theorem
- [PASS] TYPE-II ESCAPE FAILS (order estimate, flagged; lambda ~ 1 and f at 4 pi, i.e. GENEROUS): with the Delta_L-type block at the gauge scale (DYN-9b-1c) the induced triplet vev falls short of the atmospheric scale by 3.7 (PS) / 7.3 (G_LR) orders
- [PASS] POSITIVE-REAL RESCALING THEOREM (machine-verified): for y>0, Y_nu -> y Y_nu gives M_R -> y^2 M_R and M_* -> y^2 M_*, so zeta, the projective contact direction, and contact fraction are exactly invariant within this fixed-kernel one-parameter family
- [PASS] COMPLEX-RESCALING BOUNDARY (machine-verified): for complex y, M_* -> |y|^2 M_* and zeta -> (y^2/|y|^2) zeta; the contact fraction is invariant but the complex zeta phase is not
- [PASS] re-extracted M_* windows at v_R ~ M_I recorded per chain (renormalizable-source scenario at the perturbative ceiling): M_*' = 4 pi M_I by construction
- [PASS] contact essentiality carries over to the positive-real rescaled card (cf invariant; the DYN-4b posterior statement P(cf > 0.01) = 1.000 remains conditional on the Dirac kernel SHAPE only, which the rescale does not touch)
- [PASS] CONDITIONAL SOURCE COMPARISON: under the fixed-kernel rescaled source M_1 falls BELOW the Davidson-Ibarra floor on both chains (PS ~ 7e7, G_LR ~ 2e5 GeV); under a scale-decoupled (instanton-type) benchmark the archival M_1 ~ 2.8e10 stays above the order-estimate floor.  This does not select a source without a global flavor fit and flavored thermal kinetics
- [PASS] K8 NON-PROMOTION: the fixed archival tower used by the D3 benchmark has N_2 and N_3 above the intermediate gauge scale; the instanton mechanism permits but does not generically require this ordering, and the D3 card is explicitly unpromoted
