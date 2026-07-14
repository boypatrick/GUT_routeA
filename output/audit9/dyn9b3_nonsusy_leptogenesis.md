# DYN-9b-3: historical non-SUSY unflavored diagnostic

11/11 checks pass.

## Verdicts

1. **Net suppression vs DYN-7**: eps1 scales by 0.165 (SM loop function x SM Higgs normalization) -- the non-SUSY branch is ~6x HARDER for thermal unflavored leptogenesis.
2. **Unweighted prior-draw regression (no likelihood)**: target-band hit fraction = 0.0000 (DYN-7 legacy regression reference: 0.0035); median |eta_B| = 1.45e-11; boost table {'x3': 0.0, 'x5': 0.0, 'x10': 0.1265, 'x30': 0.39525, 'x60': 0.44675, 'x100': 0.3415, 'x300': 0.09425, 'x1000': 0.02675}.
3. **Reheating boundary**: only the gravitino-specific ceiling is absent.  Thermal and non-thermal production remain conditional on an explicit reheating history and initial abundance.
4. **Physics status**: blocked_missing_branch_thermal_inputs; K6 cannot be updated from these unflavored hit fractions; physics promotion is forbidden.

## Boundary (NOT claimed)

- Unflavored, N_1-dominated regression only; flavored effects, spectator rates, reheating history, and initial abundance not computed.
- No likelihood or importance weights are applied; target-band hit fractions are regression diagnostics, not posterior or viability probabilities.
- The archival tower under the D3-type source is a conditional benchmark; the source card remains unpromoted.
- zeta's value is NOT derived.

## Checks

- [PASS] premise chain-of-custody: the archival D3-type tower is an explicit conditional benchmark, not a selected source; K8 and D3 remain unpromoted while DYN-7's gravitino tension was disclosed
- [PASS] SM loop function verified: hierarchical limit f_SM -> -3/(2 sqrt x), exactly HALF the SUSY function used in DYN-7
- [PASS] central chain (benchmark card) computed with the SM function and SM Higgs normalization; m_tilde and kappa are v-independent and match DYN-7 exactly
- [PASS] NET SUPPRESSION vs the SUSY-slice pipeline quantified: eps1(non-SUSY)/eps1(DYN-7) = (f_SM/f_SUSY) x (v_arch/v_SM)^2 ~ 0.17 -- the non-SUSY branch makes thermal unflavored leptogenesis ~6x HARDER, before the gravitino gain
- [PASS] prior-draw regression completed on >= 3800/4000 draws and the SM Davidson-Ibarra ceiling is respected point-wise
- [PASS] unweighted non-SUSY prior-draw target-band hit fraction recomputed and compared with the like-for-like DYN-7 regression reference; no likelihood or draw reweighting is applied
- [PASS] fixed multiplicative boost stress test recomputed as unweighted target-band hit fractions; its maximum is a regression target, not evidence that a physical enhancement or viability is realized
- [PASS] the gravitino-specific SUSY reheating ceiling is absent on the non-SUSY branch, but no branch-local reheating history is supplied; thermal or non-thermal production is therefore OPEN, not unconstrained
- [PASS] escape hierarchy recorded for the K6 diagnostic: the needed unflavored enhancement is ~x30-100; flavored, resonant and non-thermal mechanisms remain conditional repair hypotheses rather than computed escapes
- [PASS] K6 branch-tagged non-promotion verdict recorded for DYN-8
- [PASS] promotion predicate fails closed: the source is unselected, the unweighted hit fraction is not publishable as a success probability, and physics promotion remains forbidden
