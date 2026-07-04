# DYN-2: Thresholds and Unification -- Derivation Record

`15/15` checks passed.

## Pipeline

1. **b-coefficient generator** from the DYN-1b spectrum: every heavy level
   carries its (S1,S2,S3).  Universality gates passed: chiral census sums
   to (127,127,127) = T(210)+2T(126)+T(10) and the 45 sums to (8,8,8), so
   the full-theory beta is i-independent, b_SO(10) = 109.
2. **Two-loop gauge running** (RK4 on alpha^-1): SM below M_S, MSSM above;
   self-tested against the one-loop analytic solution; the no-threshold
   limit reproduces textbook MSSM unification
   (M_12 = 1.04e+16 GeV, alpha_G^-1 = 25.87,
   alpha_3 mismatch -0.13%).
3. **One-loop supermultiplet matching** at mu* = M_X (lightest proton-decay
   vector, AG convention): Delta_i = (1/2pi) sum w S_i ln(M/M_X), eaten
   pairs at m_lambda, vector superfields at weight -3*pair.  AG's
   coefficient algebra (120, 60, (4,-9.6,5.6), 0.0167 = 2/120) verified.
4. **Viability landscape (M_S = 3 TeV benchmark).**
   The raw benchmark x = 0.1 point is solved;
   the fine x-scan records the alpha_3 pull landscape, and a
   COMPATIBILITY POINT exists at x* = 0.15:
   alpha_G^-1 = 52.515, M_X = 1.907e+13 GeV,
   m_scale = 2.662e+13 GeV,
   alpha_3 pull = +0.81 sigma.
   MC window at x* (16/50/84): log10 M_X = 13.281
   [13.276, 13.284],
   alpha_G^-1 = 52.514
   [52.507, 52.526],
   pull = +0.8
   [-0.0, +1.8] sigma.
   Away from x*, pulls reach O(100) sigma: one-loop GUT thresholds in this
   model are LARGE (consistent with AG's own conclusion that only specific
   xi regions keep corrections small); the window statement, not a point
   claim, is the deliverable.
   Secondary 3-parameter exact solve at x = 0.1 (M_S free):
   M_S_exact = 5.576e+07 GeV
   (OUTSIDE the [1,10] TeV window),
   alpha_G^-1 = 60.127.
5. **Perturbativity disclosed**: b = 109 above M_X puts the Landau pole at
   mu/M_X ~ 20.64 (famous MSGUT feature).
6. **x-scan** near the benchmark: see ledger (x = 0.20, 0.25 probe the
   disclosed accidental-zero region).

## Boundary

One-loop thresholds, gauge-only two-loop running; single effective M_S;
scheme constants dropped (supermultiplet-level DR-bar-style); PDG-vintage
electroweak inputs flagged for the DYN-4 refresh; the GeV normalization is
fixed only up to these approximations; no unique scale is claimed.
