# DYN-9b-3: non-SUSY leptogenesis rerun

10/10 checks pass.

## Verdicts

1. **Net suppression vs DYN-7**: eps1 scales by 0.165 (SM loop function x SM Higgs normalization) -- the non-SUSY branch is ~6x HARDER for thermal unflavored leptogenesis.
2. **Posterior**: P(success) = 0.0000 (DYN-7: 0.0035); median |eta_B| = 1.45e-11; boost table {'x3': 0.0, 'x5': 0.0, 'x10': 0.1265, 'x30': 0.39525, 'x60': 0.44675, 'x100': 0.3415, 'x300': 0.09425, 'x1000': 0.02675}.
3. **The gravitino gain**: the reheating ceiling is gone -- thermal production is unconstrained and NON-THERMAL production opens as a qualitatively new escape.
4. **K6 update**: a soft constraint on the source scenario (flavored + O(10) effects, non-thermal production, or a resonant D3-type tower), not a kill.

## Boundary (NOT claimed)

- Unflavored, N_1-dominated, thermal-fit order-of-magnitude estimate; flavored effects not computed.
- The archival tower under the D3-type source is a conditional benchmark; the source card remains unpromoted.
- zeta's value is NOT derived.

## Checks

- [PASS] premise chain-of-custody: DYN-9b-2 selected the scale-decoupled (D3-type) source, so the ARCHIVAL tower is the benchmark input and the K8 coexistence prediction is a standing requirement; DYN-7's gravitino tension was disclosed, not resolved
- [PASS] SM loop function verified: hierarchical limit f_SM -> -3/(2 sqrt x), exactly HALF the SUSY function used in DYN-7
- [PASS] central chain (benchmark card) computed with the SM function and SM Higgs normalization; m_tilde and kappa are v-independent and match DYN-7 exactly
- [PASS] NET SUPPRESSION vs the SUSY-slice pipeline quantified: eps1(non-SUSY)/eps1(DYN-7) = (f_SM/f_SUSY) x (v_arch/v_SM)^2 ~ 0.17 -- the non-SUSY branch makes thermal unflavored leptogenesis ~6x HARDER, before the gravitino gain
- [PASS] posterior scan completed on >= 3800/4000 draws and the SM Davidson-Ibarra ceiling is respected point-wise
- [PASS] non-SUSY posterior statistics published and compared to DYN-7: P(success) drops with the ~6x suppression
- [PASS] boost sensitivity recomputed on the non-SUSY pipeline: typical viability is reachable at SOME enhancement, and the sweet spot sits ~6x above DYN-7's (x3-10 there)
- [PASS] GRAVITINO CONSTRAINT DISSOLVED (the one unambiguous non-SUSY gain): without gravitinos there is no reheating ceiling, so thermal N_1 production at M_1 ~ 2.8e10 GeV is UNCONSTRAINED and non-thermal production channels are also freed -- DYN-7's disclosed tension is resolved BY THE BRANCH, not by tuning
- [PASS] escape hierarchy UPDATED for the K6 falsifiability entry: the needed enhancement grows to ~x30-100 (unflavored thermal), but TWO qualitatively new escapes open on the non-SUSY branch (non-thermal production unconstrained; a D3-type tower with engineered near-degeneracy for resonance)
- [PASS] K6 branch-tagged verdict recorded for the DYN-8 refresh (9b-4)
