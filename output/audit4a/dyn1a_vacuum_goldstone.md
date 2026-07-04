# DYN-1a: Vacuum Dynamics and Goldstone Eigenvector Audit

`23/23` checks passed.

- The singlet superpotential is now DERIVED and literature-matched: the
  explicit W reproduces the transcribed neutral G block exactly
  (G_chiral = 2 Hess W under a fixed field map) and its five F-terms
  vanish identically on the AG branch over the whole SM window.
- Branch classification: SM window x in (0, 1/3); the four exact
  special points; boundary behavior read off the full mixed blocks
  (x->1/3: only G degenerates = flipped SU(5)xU(1); x->0: G and F
  degenerate = G_LR).
- Goldstone eigenvector audit (previously pending): in every mixed
  sector the chiral null vector IS the gauge-orbit direction
  (alignment residual < 1e-10), J^dagger J has exactly one zero mode
  and no tachyons, all full blocks are non-singular (super-Higgs),
  and the count closes: 1+2+6+12+12 = 33 = 45-12.
- Doublet-triplet: det-linear fine-tuning gives M_H*, one light
  doublet pair with exported composition (alpha, bar_alpha); the
  triplet block stays non-singular.
- heavy_spectrum.json is now non-placeholder: 29 benchmark-mass states
  across 7 sectors, schema-conform, with partiality explicitly
  flagged (remaining AG sectors = DYN-1b).

Boundary: benchmark units of m (GeV normalization = DYN-2); sector
letter -> SM irrep assignment provisional until the AG source
cross-check (DYN-1b); tree level only; no unique vacuum claimed.
