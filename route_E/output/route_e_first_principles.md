# Route-E First-Principles Reconstruction Audit

`44/44` checks passed.

- Enumerated six-face uniqueness: dim-16 complex irreps are only the
  SU(16) fundamental (anomalous) and the Spin(10) half-spinors; all
  dim-15 candidates fail content or anomaly, so `nu^c` is forced.
- Genus ladder `h^0(Sigma_g, T) = 3, 1, 0` proves only `N <= 3`.
  The original invariant-contact hypothesis does not select a genus:
  for the one-dimensional abelian branch, `B(e,e)=1` is invariant and
  nondegenerate although the Killing form vanishes.  The added
  Killing-contact axiom conditionally selects `g = 0`, `N = 3`.
- The unique invariant contact equals the Killing form:
  `B = 2 sqrt(3) K_tr` in the spherical basis; `K_tr^2 = I/3`.
- Route-B Schur complement replays exactly: `Delta M_R = lambda^2 K_tr`,
  `lambda = sqrt(zeta)` at benchmark digits.
- All Route-A Route-B/Route-D table numbers replay, including the
  `Z_{N<=6}` miss (0.4471595 rad), the first window hit at N = 178, and
  the convergents 17/178, 70/733, 87/911.
- Two-model witness: the structural checks are zeta-independent, so no
  zeta value follows from the principle set.

Boundary: `zeta` is not derived; the no-extra-chiral-sector assumption
and the Route-B R-selection rule are not discharged; this core audit
does not replay the separate dynamics lanes; no PSLT-only unconditional GUT proof
is claimed (the audit's Theorem-5 witness proves the opposite).
