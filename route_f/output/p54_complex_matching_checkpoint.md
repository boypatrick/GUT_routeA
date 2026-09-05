# P54 complex matching and PS evolution checkpoint

Date: 2026-09-05. Previous checkpoint committed/pushed first as `b8e7398`
on `origin/main` at `https://github.com/boypatrick/GUT_routeA`.
These are subsequent working-tree changes, not an additional push.

## Completed bounded calculations

- Actual common complex scalar/spinor dictionary, with 28/28 checks.
  Up/nu uses `(c1,c4)`, down/e uses `(c2*,c3*)`; the Majorana and LL
  coefficients are `+i4sqrt(2)` and `-i4sqrt(2)`. This is a transported
  convention, not a new CP parameter. Old magnitude-only maps are obsolete.
- Same-action upper PS stationary saddle and finite gauge threshold,
  25/25 checks. `lambda_U=(2.94405268936,7.60297883009,7.60297883009)`
  includes vector constants. The retained tachyons are not integrated.
  Actual-vacuum heavy-valley Schur mass and kinetic metric are exported,
  with 13 zero modes and no negative generalized eigenvalues.
- Complete four-parent all-active one-loop PS Yukawa contraction and
  coupled ODE, plus the Yukawa contribution to the two-loop gauge flow.
  The exact generic/closed tensor residual is `2.28e-16`.
  Parity protects `F=sqrt(2)L=sqrt(2)R` at one loop on the Spin(10)
  boundary locus. Finite-matched off-locus data require the general flow.
- Finite canonical Yukawa matching/kernels, 17/17 tests on synthetic
  finite inputs with actual scalar geometry. Missing physical diagrams
  are rejected by the fit-readiness interface, not set to zero.
- Actual-Clebsch local scalar/flavor feedback, 36/36 tests on two
  synthetic complex family choices. Analytic light-state derivatives
  agree with re-solved finite differences. The full cubic triplet source
  is recomputed for each new light vector using eight existing cached
  Hessians, rather than freezing the previous type-II coefficient.

## What remains open

The positive heavy block gives a local EFT branch; it does not establish
a uniform hierarchy of all retained/integrated masses. The lower finite
threshold needs the upper-matched covariant action, including the
momentum-dependent Schur operator and induced gauge vertices. Neither
parent-weighted logarithms nor generalized zero-momentum masses are that
functional. The one-step SM result is a regression, not a silently changed
physical branch.

Both-site P54 finite Yukawa vertex/kinetic state sums, actual mass-ordered
decoupling, SM/sterile-neutrino/Weinberg evolution, and complete loop scalar
feedback are still required for a physical self-consistent global fit.
No new physical scales, best fit, exclusion, Higgs pole mass, determinant
or portal are claimed. The P/U/S separation is unchanged.

## Reproduce

Use Python 3.9 with `route_f/code/requirements_p54_matching.txt`. The
recorded isolated dependency directory is `/private/tmp/p54_matching_deps`.
From the repository root run, in order:

```sh
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_common_yukawa_phase.py
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_ps_finite_thresholds.py
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_ps_yukawa_flow.py
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_finite_yukawa_interface.py
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_self_consistent_light.py
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/p54_matching_deps /usr/bin/python3 route_f/code/verify_p54_action_card.py
```

Each verifier exports a same-named JSON/MD under `route_f/output`, including
input/source SHA256 values and individual regression results. Runtime/cache
metadata can vary between equivalent replays. The PS code exports both
actual scalar embeddings and the compact beta API for finite matching.

Derivations and claim boundaries:
`route_f/tex/p54_complex_phase_ps_matching.tex` and
`route_f/output/pdf/p54_complex_phase_ps_matching.pdf`.
The action card, Route-F roadmap and root roadmap are synchronized.
