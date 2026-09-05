# PS finite gauge-threshold audit

The four-parent active EFT is unchanged. The upper PS-symmetric coefficient and a full one-step SM regression are computed; the finite lower staged matching is still open.

Checks: **25/25**.

## Upper same-action PS background

- VEVs: `{'omega': 1.0008103059215108, 'sigma': 0.0, 'vs': 0.2542776970837077}`; full gradient norm `7.729e-16`.
- Retained active fields: `248` real, `124` tachyonic modes.
- Integrated physical fields: `55` real; smallest mass squared `0.102882040959`.
- Active tachyons describe the lower symmetry breaking and are not included in upper hard logarithms.
- This is a tree stationary saddle of the frozen parameters. It is not claimed to be stationary for the one-loop fixed-VEV retuned action.
- Largest retained / smallest integrated squared-mass ratio: `1.01261`. There is no uniform hierarchy between every retained and integrated scalar.
- Finite upper lambda in `(4,L,R)`: `[2.944052689359527, 7.602978830094123, 7.602978830094121]`.
- Upper scale slopes: `[71.99999999999996, 119.99999999999994, 119.99999999999993]`.

## Lower matching remains a functional calculation

At simultaneous nonzero VEVs, `||[Pactive,H]||=0.265322` and `||[Pactive,P_G,U]||=0.60996`.
The same site-resolved gauge orbit is used for vector and Goldstone projectors. Its upper Goldstone image has active-parent support, so a naive complementary scalar projector is not a staged EFT.
The exact scalar Schur operator contains momentum dependence and induced gauge vertices. A weighted mass logarithm, or only the zero-momentum Schur mass, does not compute that finite matching.

## Actual-vacuum tree heavy valley

- The same 55 heavy scalar directions have positive `C`, with minimum mass squared `0.0994502574109`.
- The upper gauge slice uses `B=(I-PGupper)Bactive`, with source `W=BH^T H B` and response `D=-C^-1 W`.
- Its tree Schur mass and pullback metric are `S=B^T H B-W^T C^-1 W` and `G=B^T B+D^T D`.
- Metric eigenvalues range from `0.984249788817` to `1.0180072705`; `||D||2=0.134191171461`.
- The generalized pair `(S,G)` has `13` zeros (nine lower Goldstones and four Higgs coordinates) and no negative eigenvalues. Its lower-orbit Ward residual is `3.207e-16`.
- `heavy_valley_api(...)` exports all quadratic source, metric and embedding matrices. These derivative-expansion masses do not replace the finite lower functional or its gauge vertices.

## One-step unbroken-SM regression

- Complete finite lambda `(3,2,1)`: `[36.63010815544616, -6.547522795837455, 119.73803707594224]`.
- Finite vector constants: `[4.999999999999997, 5.9999999999999964, 7.999999999999995]`.
- Measured scale slopes: `[25.999999999999954, 48.99999999999997, 92.5999999999999]`.
- All six staged logarithmic identities close and sum to this one-step slope. No physical two-stage scale refit is claimed.

The fixed-VEV mass counterterms retain the original invariant normalization, including `deltaH_nu=-delta_nu2 P126/2`. One-loop counterterm insertion into these one-loop threshold mass logs is beyond the stated order.

Finite-threshold convention: [Hall](https://doi.org/10.1016/0550-3213(81)90498-3), [SO(10) threshold calculation](https://doi.org/10.1140/epjc/s10052-020-8308-9).
