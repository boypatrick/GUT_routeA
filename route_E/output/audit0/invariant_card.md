# Audit 0 Convention / Invariant Card

Audit 0 fixes the common convention layer for the deferred companion audits. It is not a flavor fit, threshold closure, proton lifetime calculation, or UV completion.

## Digest

- card sha256: `a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc`
- source artifacts: 7/7 present

## Core Anchors

- `zeta = 0.1076472949 + 0.0736514853 i`
- `|zeta| = 0.13043190325293763`
- `arg(zeta) = 0.600038020318215 rad`
- `|sqrt(zeta)| = 0.3611535729477664`
- `|lambda|^2/(16 pi^2) = 8.259696764e-04`
- `K_tr inverse check = 3.846e-16`
- `M_R contact fraction = 1.304275167e-01`
- `Veronese + contact residual = 2.687e-17`

## Dependency Graph

`0 -> (1 || 4a) -> 3 -> 2`, with `4b` optional and parallel.

| audit | role |
|---|---|
| 0 | convention/invariant card |
| 1 | full flavor and seesaw fit |
| 4a | field-theory source-sector UV/heavy-spectrum interface |
| 4b | optional Route-D global/string geometry track |
| 3 | threshold and unification audit using source spectrum |
| 2 | d=5 proton Wilson/lifetime audit consuming all previous outputs |

## Precision Policy

Long decimals are retained as reproducibility anchors. They are not ten-digit physical predictions; follow-up audits should treat only the first few significant figures as physically meaningful unless a concrete UV replay justifies more.

## Next Required Inputs

- Audit 1 must set the final CKM/PMNS/mass target table and fit convention.
- Audit 4a must provide a heavy-spectrum schema before Audit 3 can be paper-grade.
- Audit 2 must not be run as a final proton-lifetime audit until flavor rotations, source spectrum, and threshold-fixed scales are supplied.
