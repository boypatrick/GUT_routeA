# Audit 0.5 Phase-Transfer Diagnostic

This diagnostic compares `arg(zeta)` against the finite visible-spurion table `{arg I, 2 arg I, arg J, arg I + arg J}`.

Real single-coefficient monomials are tested modulo pi because the real coefficient may have either sign. The `2 arg I` row is tested modulo 2 pi because it is interpreted as `lambda = c I`, followed by `zeta = lambda^2` with `c^2 > 0`.

## Digest

- card sha256: `cb8f3b3dbec52185b7facffcfa4cd1567fc3a40c983bbf74fe62fe1f8cf97506`
- audit0 sha256: `a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc`
- source artifacts: 2/2 present

## Result

- loose hit count: 0
- tight hit count: 0
- closest candidate: `arg_I`
- closest absolute phase error: `5.372532875265e-01` rad
- loose window: `5.424119231661e-05` rad
- tight window: `1.799256517287e-06` rad

## Two-Term Diagnostic

The real two-term ansatz `zeta = c1 I + c2 J` is solved exactly when the real 2x2 system is nonsingular. This is a fit-cost diagnostic, not a prediction.

Raw `c1,c2` sizes are convention-dependent because they scale with the chosen normalization of `I,J`; interpret them together with the recorded basis, invariant convention, and term magnitudes.

- status: `solved`
- c1: `1.037751489216e+01`
- c2: `-2.163199491664e+02`
- residual abs: `1.387778780781e-17`
- naturalness index max(|c1 I|, |c2 J|)/|zeta|: `1.646438908603e+00`
- cancellation index (|c1 I|+|c2 J|)/|zeta|: `2.585476554589e+00`

## Candidate Table

| candidate | distance | phase rad | signed delta rad | loose | tight |
|---|---|---:|---:|---|---|
| `arg_I` | `mod_pi` | 1.137291307845e+00 | 5.372532875265e-01 | False | False |
| `two_arg_I` | `mod_2pi` | 2.274582615689e+00 | 1.674544595371e+00 | False | False |
| `arg_J` | `mod_pi` | 1.713683799838e+00 | 1.113645779520e+00 | False | False |
| `arg_I_plus_arg_J` | `mod_pi` | 2.850975107683e+00 | -8.906555662251e-01 | False | False |

## Boundary

Audit 0.5 only tests a finite visible-spurion phase-transfer table against the local benchmark.
