# DYN-5: hidden-messenger dynamics at one loop

21/21 checks pass.

## Result

The tree-level matching-silence theorem UPGRADES cleanly to one loop,
with one structural sharpening and one disclosed non-silence:

1. **delta Z_N is exactly family-universal.**  The portal coupling is
   `lambda * 1_3` and the messenger mass matrix `M_* K_tr^-1` has all
   singular values `sqrt(3) M_*` (from `K_tr^2 = I/3`), so the
   leading-log correction is `deltaZ * identity`.  Anchor:
   `|lambda|^2/(16 pi^2) = 8.259697e-04` per e-fold (tex caveat
   digits reproduced).  Strict-benchmark interval `-ln(3)/2` gives
   `|deltaZ| = 4.537e-04`.
2. **Light-sector silence is STRUCTURAL, not perturbative.**  The
   seesaw `m_nu = -mD M_R^-1 mD^T` is invariant under `mD -> mD A`,
   `M_R -> A^T M_R A` for ANY invertible `A` (machine-zero check,
   universal and non-universal): canonical-normalization corrections
   cancel identically in every light observable, to all orders in
   `deltaZ` at fixed matching order.
3. **The paper's linearized replay formulas have a quantified
   validity range.**  Residual is O(deltaZ^2); it crosses the loose
   contact-sensitivity window at `deltaZ = 1.13e-02`
   (ln M ratio ~ 13.7) and the refreshed
   window at `deltaZ = 2.00e-03`
   (ln M ratio ~ 2.4).  Beyond that the
   exact A-form must be used.
4. **The heavy sector is NOT silent** (disclosed): `M_R` singular
   values, `M_1` (leptogenesis input) and `M_*` rescale by
   `(1 - deltaZ)`; worst scanned shift = 0.0025 dex =
   3.6% of the DYN-4b M_1 window.
   The normalized zeta extraction is exactly invariant.
5. **R-selection silence at one loop.**  All X-decorated
   Dirac/triplet/Majorana operators have `R = 2 + m`; the portal
   forces `R(X) = 1`; silence fails only at the measure-zero
   assignment `R(16) = 2`.  Non-renormalization makes one-loop
   effects Kahler-only, and the `nu^c` leg appears in no visible
   Yukawa or d=5 channel, so `delta Y_{u,d,e} = delta C_5L =
   delta C_5R = 0` at one loop and the DYN-2 threshold vector is
   untouched (`X`, `N` visible singlets).
6. **R-breaking spurion conditions** (regeneration iff `k r_s = -m`):
   positive-R spurions are safe at the holomorphic level; even-R
   spurions (constant-W/gravitino type, `r_s = -2`) preserve a
   residual `(-1)^R` parity so all ODD decorations stay absent to all
   orders; odd-R spurions are the dangerous case.  Order-estimate
   ceilings: `eps_odd_loose = 1.4e-02`, `eps_odd_tight = 4.4e-04`, `eps_even_loose = 5.1e-03`, `eps_even_tight = 1.6e-04`

## Boundary (NOT claimed)

- zeta is NOT derived (Theorem VII.4 boundary intact).
- The R symmetry is NOT claimed anomaly free; no UV implementation.
- NO microscopic messenger completion is constructed; the interval
  `ln(M_*/M_X)` is scanned, not derived.
- Two-loop silence NOT claimed (leakage estimate 1.5e-04 in alpha^-1 units, via the
  `Tr(Ynu+ Ynu)` trace term; flagged order estimate).
- Spurion ceilings are order-of-magnitude estimates.

## Checks

- [PASS] lambda = sqrt(zeta): lambda^2 - zeta machine zero and tex digits match
- [PASS] |lambda|, |zeta|, arg zeta reproduce the tex anchors
- [PASS] companion caveat anchor |lambda|^2/(16 pi^2) = 8.259697e-4
- [PASS] higher-order remark |zeta|^2 ~ 1.7e-2 reproduced
- [PASS] K_tr^2 = I/3 and K_tr^-1 = 3 K_tr (fixed convention)
- [PASS] messenger mass matrix M_* K_tr^-1 has ALL singular values sqrt(3) M_* => exact family degeneracy => delta Z_N proportional to identity at leading log
- [PASS] portal coupling matrix is lambda * 1_3: (Y+ Y)/(16 pi^2) = per-efold coefficient times identity
- [PASS] strict Route-B benchmark interval ln(M_*/M_X) = -ln(3)/2 (X sits AT sqrt(3) M_*): |delta Z_N| = 4.537102e-4
- [PASS] archival closure card present with the DYN-4a provenance sha256
- [PASS] benchmark replay reproduces the zeta anchor (links S3 to DYN-4a)
- [PASS] EXACT light-sector silence: m_nu invariant under the canonical redefinition (universal AND non-universal delta Z_N) at machine level given cond(M_R) ~ 1.6e5
- [PASS] linearized replay residual scales as O(deltaZ^2) (halving deltaZ quarters the residual)
- [PASS] at the strict benchmark |deltaZ| = 4.53e-4 the linearized-formula residual sits BELOW even the refreshed window
- [PASS] heavy sector NOT silent: all M_R singular values shift by exactly (1 - deltaZ); worst scanned shift stays inside the DYN-4b M_1 window
- [PASS] normalized zeta extraction EXACTLY invariant under delta Z_N (only M_* rescales by 1 - deltaZ)
- [PASS] R enumeration: all base channels have R = 2; every X-decorated Dirac/triplet/Majorana operator has R = 2 + m (zero violations)
- [PASS] portal forces R(X) = 2 - R(16); silence fails ONLY at the measure-zero assignment R(16) = 2 (benchmark R(16) = 1 is safe)
- [PASS] SUSY non-renormalization => one-loop effects are Kahler-only; the nu^c leg appears ONLY in Y_nuD => delta Y_u = delta Y_d = delta Y_e = delta C_5L = delta C_5R = 0 at one loop; X, N visible singlets => Delta b = (0,0,0), DYN-2 threshold vector untouched
- [PASS] two-loop leakage (ORDER ESTIMATE): shift in alpha^-1 over the full run < 1e-3, i.e. far below the DYN-2 threshold-sensitivity scale O(0.1); two-loop silence itself is NOT claimed
- [PASS] spurion regeneration condition k r_s = -m: r_s = -1 regenerates m = 1 at first order; r_s = -2 (constant-W / gravitino type) regenerates ONLY even m (m = 2 at k = 1), never the m = 1 portal decoration; positive-R spurions never regenerate at the holomorphic level
- [PASS] spurion ceilings vs DYN-4 windows computed (ORDER estimates): loose window tolerates eps ~ 1e-2, refreshed window forces eps below ~ 1e-3
