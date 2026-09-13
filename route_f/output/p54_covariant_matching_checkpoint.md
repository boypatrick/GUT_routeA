# P54 covariant matching checkpoint

Date: 2026-09-06, Asia/Taipei.

The previous checkpoint was committed and pushed before this work:
`4ad74e7` (`Complete common phases, PS Yukawa flow, and local scalar feedback`),
`main`, remote `https://github.com/boypatrick/GUT_routeA`.
This document records the subsequent working-tree research, not a second push.

## What actually changed

- The same-action lower scalar valley and heavy-vector response now define
  a local covariant two-derivative action with actual vertex tensors.
  Near-crossing quadratic modes are handled by an exact Schur resolvent;
  they do not justify a new global hierarchy requirement or parameter scan.
- Finite upper pure-Yukawa graphs generate two conjugate-bidoublet
  invariant directions. Their coefficients are fixed cubic polynomials in
  the original two UV matrices times split mass logarithms. The six-invariant
  real-Weyl beta system closes; a four-matrix finite boundary does not.
- The broken-background scalar-cubic hard bubble yields a complete
  four-complex-doublet kinetic matrix. Both mixed hard-soft orientations,
  all charged partners and actual imaginary components are retained.
- The exact tree light ray has nonzero
  `omega C_HN = -0.0848692011507 i f_raw`, independently confirmed by
  nonlinear stationary heavy-scalar relaxation. A zero `NN H†H` boundary
  is not the actual P54 matching condition.
- A moving-threshold matrix engine replaces frozen-mass decoupling.
  Degenerate Takagi clusters and running level crossings are supported.
  The 2024 dimension-five feedback correction is used, not the historical
  Weinberg-only anomalous dimension. Nonzero `C_HN` and the Higgs quadratic
  coupling must travel with the enlarged lower EFT.

No extra UV family matrix, field multiplet, torsion term, lattice scan,
physical scale or experimental fit was introduced.

## Reproduce

Final regression ledger (156/156 passing checks):

| Verifier | Checks | Scope |
|---|---:|---|
| `p54_lower_covariant_eft` | 24/24 | local covariant tree action and exact quadratic resolvent |
| `p54_upper_yukawa_thresholds` | 44/44 | actual finite pure-Yukawa diagrams and six-invariant RG closure |
| `p54_scalar_kinetic_matching` | 24/24 | actual one-step scalar-cubic kinetic subset |
| `p54_scalar_chn` | 18/18 | nonzero exact-tree sterile-Higgs coefficient |
| `p54_sequential_seesaw` | 46/46 | moving thresholds and the single-insertion, dipole-free dimension-five flow |

Passing these tests does not promote any of the omitted finite matching
or physical-endpoint conditions below to a completed result.

Run from `/Users/boypatrick/codex/another_physics` with the existing runtime:

```sh
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/private/tmp/p54_matching_deps
/usr/bin/python3 route_f/code/verify_p54_lower_covariant_eft.py
/usr/bin/python3 route_f/code/verify_p54_upper_yukawa_thresholds.py
/usr/bin/python3 route_f/code/verify_p54_scalar_kinetic_matching.py
/usr/bin/python3 route_f/code/verify_p54_scalar_chn.py
/usr/bin/python3 route_f/code/verify_p54_sequential_seesaw.py
```

The first, third and fourth programs reuse the content-addressed Hessian
cache under `tmp/p54_full_doublet_cw`. They reject absent required cache
entries rather than silently recompute or substitute a different action.
The neutral-CHN verifier also independently differentiates and solves the
three-radial-variable restriction of the original potential; this is not
another full 328-field Hessian evaluation or new vacuum scan.

Every same-stem JSON records the complete numerical matrices, test results
and source hashes. Every same-stem Markdown records its derivation and
scope. The combined derivation is
`route_f/tex/p54_covariant_eft_finite_yukawa_seesaw.tex`; its PDF is under
`route_f/output/pdf/`. Build with the existing `pdflatex`, repeating until
cross-references settle, then render all pages for layout review.

## Accuracy and claim boundaries

The upper finite pure-Yukawa result and the one-step broken scalar-cubic
kinetic result are distinct subsets. They must not be added as though
they were independently complete two-site coefficients. Goldstone
propagators in the latter inherit background-field Landau gauge. The
exact square-root canonical normalization tests tensor algebra, not a
two-loop prediction.

The exact-tree CHN coefficient is the primary result. Its projections on
bosonic/local-loop-improved light rays use frozen tree vertices and are
explicitly mixed-order. Synthetic gauge/quartic/reference-scale inputs
in the evolution engine are not matching solutions; terminal CKM, PMNS
and mass proxies are diagnostics, not predictions. In particular,
`qH=0` is regenerated to values above the final scale squared; the
positive Higgs mass cannot consistently be treated as light there
without quadratic matching/retuning. The MS-bar continuation is a
mathematical diagnostic, not a physically valid light-Higgs endpoint.

Still uncomputed: the heavy-vector fermion kinetic/vertex package,
vector/scalar/Goldstone/ghost scalar kinetic package, upper scalar/tadpole
Wilson matching, lower covariant finite hard-minus-EFT subtraction, and
finite sequential type-I/type-II matching with the surviving operators.
Dimension-six effects, finite dipole matching, quantum axion effects
(the axion is presently a background spectator), Higgs pole self-energies
and Nielsen consistency are not supplied by these tests.

The next work is the remaining common-scheme diagram/Wilson calculation,
followed by physical scale solving and the constrained scalar/flavor
fixed point. The latest sections of `route_f/ROADMAP.md`, `route_f/README.md`
and the root `roadmap.md` are synchronized with these operator-basis changes.
