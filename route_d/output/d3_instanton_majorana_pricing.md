# Route-D D3: instanton / d=5 Majorana-source pricing card

15/15 checks pass.  CONDITIONAL string interpretation, NOT promoted; zeta NOT derived.

## Verdicts

1. **d=5-operator escape CLOSED** (negative gate): `f' = f_renorm * (M_s/v_R)` is strictly worse; PS chain at `M_s = 2e16` needs `f'_top = 1.2e+08`.  The renormalizable clash factor is `v_R_min/M_I ~ 379`.
2. **Instanton escape priced**: required actions at `M_s = 2e16` are `S = [1.63, 6.43, 13.64]` (top/mid/bottom); D2 contact action `S_zeta = 2.0369` is same-ballpark as the tower top (diagnostic, not evidence).  The diagnostic `S_top >= 1` needs `M_s >= 1.07e+16` GeV: OK on G_LR and above `M_X` on PS; it does not prove semiclassical control.
3. **Rank/texture (conditional source ansatz)**: the tower is exactly three rank-1 Takagi terms, so the rank-1-per-cycle ansatz needs at least three cycles; instanton data = 3 + 9 = 12 real parameters = the full Majorana data.  **The axiom buys scale, not texture.**
4. **Fixed archival-tower consequence**: N_2, N_3 are heavier than the intermediate gauge scale (PS: 39x, 4763x M_I) -- impossible for perturbative f v_R, but permitted rather than generically required by the instanton benchmark.

## Boundary (NOT claimed)

- No global compactification, divisor, zero-mode spectrum, or moduli stabilization is constructed.
- The family texture u_k is a NEW conditional input.
- zeta's value is NOT derived (boundary theorem intact).
- Route-D promotion bar NOT passed; this card must not be cited as evidence.

## Checks

- [PASS] premise: the patched DYN-9 ledger declares the archival zeta card NOT transplantable and the tower NOT perturbatively realizable on surviving chains
- [PASS] archival closure card sha256 matches the DYN-4a provenance
- [PASS] f'_needed = f_renorm * (M_s / v_R) is STRICTLY worse than the renormalizable coupling for every surviving chain and every M_s in the grid, and always violates 4 pi at the tower top
- [PASS] anchor: on the 210-compatible PS chain at M_s = 2e16 the d=5 source needs f' ~ 1.2e8 at the tower top
- [PASS] structural clash quantified: a renormalizable tower top needs v_R >= sigma_top/(4 pi) = 3.1e14 GeV, while unification fixes M_I(PS) = 8.2e11 -- a factor ~380 that the d=5 operator can only worsen; the escape is CLOSED at intermediate-scale v_R
- [PASS] D2 link: S_zeta = -ln|zeta| reproduces the D2 ledger value effective_action_if_prefactor_one
- [PASS] required actions S_i = ln(M_s/M_R,i) at M_s = 2e16: [top, mid, bottom] = [1.63, 6.43, 13.64]
- [PASS] minimal exponential-suppression diagnostic: S_top >= 1 requires M_s >= e * M_R,top = 1.07e16 GeV -- satisfied by the G_LR chain M_X = 2.1e16 but NOT by the PS chain M_X = 5.4e15; this is not a proof of string-loop or alpha-prime control
- [PASS] same-ballpark DIAGNOSTIC (not evidence): at M_s = 2e16 the top-of-tower action and the D2 contact action differ by < 1
- [PASS] replayed archival M_R singular values match the DYN-4a ledger
- [PASS] Takagi factorization M_R = U diag(sigma) U^T verified (self-implemented, residual machine-zero)
- [PASS] under the explicit rank-1-per-cycle source ansatz, the Takagi tower uses exactly three positive rank-1 terms; this is a conditional UV cost, not a theorem about general instanton sectors
- [PASS] parameter count: instanton data (3 actions + 9 texture) = 12 = the FULL complex-symmetric Majorana data -- the axiom buys SCALE (f < 4 pi evasion), ZERO texture compression
- [PASS] fixed-tower consequence: on BOTH surviving chains the replayed archival N_2 and N_3 are heavier than the intermediate gauge scale and in fact satisfy sigma/M_I > 4 pi; the instanton benchmark permits this coexistence without a perturbative f v_R source, but does not impose it on a generic independently fitted tower
- [PASS] price card assembled: statement, 12-parameter cost, buys / does-not-buy lists, beyond-zeta consequence
