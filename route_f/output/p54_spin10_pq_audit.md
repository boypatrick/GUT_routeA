# P54PQ-v2 Spin(10) invariant-basis and PQ audit

Status: **17/17 checks passed**.

## Action-level result

The independent basis audit found one omitted operator: `chi7 Phi_ij phi_i phi_j S* + h.c.`. Therefore `P54PQ-v1` is superseded by `P54PQ-v2`.

The corrected counts are `60` raw, `49` classical-quotient, and `48` continuous quantum-PQ observable parameters. Three scalar CP phases survive.

## PQ anomaly and domain walls

- `A[Spin(10)^2-PQ] = -6`.
- `A[SU(3)c^2-PQ] = -6` and `Nhat = -12`.
- The naive scalar-gcd count is `6`; quotienting by the diagonal Spin(10) center gives physical `N_DW = 3`.
- `N_DW=3` is not cosmologically harmless: a post-inflation PQ transition needs an additional repair.

## Verification

| Group | Check | Result |
|---|---|---|
| basis | every admitted multidegree is PQ neutral | PASS |
| basis | every neutral multidegree is admitted or explicitly excluded by Spin(10) channels | PASS |
| basis | coefficient multiplicity matches invariant multiplicity | PASS |
| parameters | corrected scalar action has 34 raw real coefficients | PASS |
| parameters | five coupling phases have rank-two rephasing orbit | PASS |
| parameters | three scalar CP invariants remain | PASS |
| parameters | full corrected counts are 60/49/48 | PASS |
| tensor | 126 x 126bar has no symmetric-traceless 54 bilinear | PASS |
| tensor | the absent 54 test is non-vacuous because the 45 bilinear survives | PASS |
| alignment | 54 breaks Spin(10) to the 21-generator Pati-Salam algebra | PASS |
| alignment | 54 plus 126 leaves exactly the 12-generator SM algebra | PASS |
| anomaly | Spin(10)^2-PQ coefficient is -6 | PASS |
| anomaly | QCD instanton coefficient is Nhat=-12 | PASS |
| global-form | a pi/2 PQ rotation is the inverse Spin(10) Z4 center action | PASS |
| domain-wall | naive N_DW=6 is reduced to physical N_DW=3 | PASS |
| domain-wall | minimal model retains a post-inflation domain-wall problem | PASS |
| artifact | companion TeX contains basis, anomaly, and global-form interfaces | PASS |
