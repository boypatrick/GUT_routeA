# P54PQ-v2 action-card verification

Status: **20/20 checks passed**.

The primary P-layer action is `Spin(10) x U(1)_PQ` with `3x16F + 54R + 126C + 10C + 1C`. The singlet is required; the no-PQ short-field-list theory is the comparison branch `P54N`.

## Definition-level correction

The printed `Sigma Sigma* Sigma* phi` eta1 monomial has PQ charge `-4`. The frozen card uses `Sigma Sigma* Sigma phi + h.c.`, which has charge zero.

## Parameter count

- Scalar sector: `34` raw reals, `32` after the rank-`2` scalar rephasing orbit.
- Yukawa sector: `24` raw reals and `15` after the generic `U(3)` family quotient.
- Full action: `60` raw, `49` classically basis-quotiented, `48` continuous zero-temperature observable parameters after the anomalous PQ reparametrization.
- Three scalar CP phases survive: `arg(eta3)-2arg(eta1)`, `arg(chi4)+arg(chi6)-2arg(eta1)`, and `arg(chi7)-arg(chi6)`.

## Checks

| Group | Check | Result |
|---|---|---|
| action | all scalar operators have dimension four | PASS |
| action | all scalar operators are PQ neutral | PASS |
| action | complex scalar operators carry an explicit Hermitian partner | PASS |
| correction | printed eta1 monomial fails continuous PQ | PASS |
| correction | corrected eta1 monomial is PQ neutral | PASS |
| correction | corrected eta1 and chi4 contractions are index balanced | PASS |
| parameters | scalar raw count | PASS |
| parameters | scalar rephasing rank | PASS |
| parameters | three scalar CP phases survive | PASS |
| parameters | two symmetric Yukawa matrices give 15 physical reals | PASS |
| parameters | joint classical and quantum counts | PASS |
| vacuum | 54 Pati--Salam direction is traceless | PASS |
| vacuum | 54 radial potential fingerprints | PASS |
| vacuum | 54 commutant has Pati--Salam dimension | PASS |
| representations | Pati--Salam decompositions preserve dimensions | PASS |
| representations | Yukawa symmetry follows from representation channel | PASS |
| branch | primary field contract is exact | PASS |
| branch | comparison representations are absent from primary | PASS |
| branch | exactly two PQ-allowed Yukawa channels | PASS |
| artifact | TeX card contains all required interfaces | PASS |

## Remaining downstream debts

- PQ-quality policy for higher-dimension operators
- P2 running and threshold matching on the completed P1 spectrum
- scaled-SVD Jacobian rank after the spectrum and matching maps exist
