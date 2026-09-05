# P54 fixed-VEV and matching checkpoint

Date: 2026-09-05. This is a completed bounded calculation, not a completed physical global fit.

Before this work the prior theory audit was committed and pushed as `cb60bdf` on `origin/main`. These are subsequent working-tree changes.

## Results and boundaries

| Calculation | Regressions | Result |
|---|---:|---|
| Fixed-VEV bosonic full complex doublet | 18/18 | All 4 complex copies, invariant tadpoles and Schur-retuned bosonic light state computed |
| PS active census and matching logs | 14/14 | Two historical physical gates fail; all 6 repaired logarithmic identities close |
| Spin(10) Clifford/CG audit | 20/20 | Absolute normalization, relative -3 and triplet alignment checked; common copy phases still open |
| Fermionic fixed-VEV interface | 11/11 | Synthetic complex-family tests only, not fitted matrices |
| Complete mixed triplet scalar source | 10/10 | Tree K/V3 with bosonic-improved light direction; not a loop-accurate neutrino coefficient |

The invariant finite tadpole rule is

\[
\delta\mu^2=\frac{5t_\omega}{12\omega},\quad
\delta\nu^2=\frac{t_\sigma}{\sigma},\quad
\delta\mu_s^2=\frac{t_s}{v_s},\qquad
H_{\rm ct}=-\delta\mu^2P_{54}-\frac{\delta\nu^2}{2}P_{126}-\delta\mu_s^2P_S.
\]

The bosonic one-loop-truncated heavy doublet eigenvalues are
`(0.03025666325,0.19979763559,0.43263922302) omega^2`.
The light-state magnitudes are
`(0.6462011668,0.7631652015,0.00001218358,0.0017108974)`.
The original coarse Q/gap criterion fails; the exact relative matrix test instead gives `C >= 0.9890936 C0 > 0`. This is positivity at the truncation, not a proof of perturbative convergence or all-sector loop stability.

The two structural repairs are:

1. The PS beta table and threshold census must describe the same active fields. The old lower SU(2)L matching-log slope is zero but must be `-71`. The four-parent repair fixes all six log identities; finite staged matching and replacement scales are still needed.
2. With canonical P1 scalars, the positive CG normalizations obey `h_D=sqrt(2)Y10`, `f_D=2Y126/sqrt(3)` and `f_M=2sqrt(6)f_D`. Thus `MR=sigma f_M`, not `sigma f_D`. This changes conventions, not field content. All scalar/spinor conjugations and phases must be transported together.

The action card, Route-F roadmap and root roadmap are synchronized. The detailed derivations are in `route_f/tex/p54_fixed_vev_full_doublet_matching.tex` and `route_f/output/pdf/p54_fixed_vev_full_doublet_matching.pdf`.

## Reproduction

From the repository root, use an isolated Python 3.9 environment with `route_f/code/requirements_p54_matching.txt`. The recorded run used NumPy 2.0.2, SciPy 1.13.1 and JAX/JAXLIB 0.4.30 with x64 enabled by P1. Then execute in this order:

```sh
python route_f/code/verify_p54_ps_active_census.py
python route_f/code/verify_p54_full_doublet_cw.py
python route_f/code/verify_p54_typeii_triplet_source.py
python route_f/code/verify_p54_spinor_intertwiners.py
python route_f/code/verify_p54_fermion_tadpole.py
python route_f/code/verify_p54_action_card.py
```

The full-doublet run evaluates the actual 328-real scalar Hessian at polynomial jet points; a cold run is substantially more expensive than the other tests. Its ignored `tmp/p54_full_doublet_cw` cache is keyed by action source, parameters and exact field coordinates. An unchanged replay does not perform a new lattice or parameter scan. JSON ledgers record source hashes, complete complex bases/matrices and residuals. Runtime/cache-count metadata can differ across replays without changing the physics tensors.

## Next physical gate

Export the common phase-resolved scalar/spinor projections; construct the finite, site-resolved four-parent PS gauge and Yukawa thresholds and full left/right beta system; rerun the matching scales; then iterate actual fermionic tadpoles/CW, the Schur light eigenpair and both seesaw contributions in the constrained likelihood. Do not treat an undefined EFT map as an empty fit domain. No physical best fit, full-model exclusion, pole mass or soliton-to-family promotion is asserted here.
