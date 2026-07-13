# Route-E `code_dyn` Re-audit

Date: 2026-07-13
Execution addendum: 2026-07-14

## Superseding execution addendum

The 2026-07-13 path/deployment findings below are retained as the discovery
record but are no longer current:

- all recovered scripts now import `route_e_paths.py`; canonical direct
  execution and isolated output/cache work without aliases;
- `run_route_e_dynamics.py` is a 21-node fail-closed DAG (19 recovered + two
  validity guards) with a digest/provenance manifest and a separate
  `dyn_claim_registry.json` physics status;
- an isolated DYN-0 -> DYN-4a replay succeeds; DYN-5/DYN-7 guard lanes also
  pass mechanically while correctly prohibiting physics promotion;
- DYN-4a is now a converged connected-local profile (`23/23`), not the former
  fixed-coordinate pseudo-profile; no global-basin completeness is claimed;
- DYN-5V (`9/9`) establishes the tree matching and invalidates the historical
  bilinear-as-Yukawa loop; DYN-7F (`7/7`) enforces the flavored-regime blocker
  and repairs the Davidson--Ibarra double square root.

Still open: the heavy full replay/clean-clone gate, `RE-SC3/4/5`, an explicit
interacting DYN-5 completion, flavored DYN-7/9b-3 kinetics, and DYN-8's stale
branch map.  When a statement below conflicts with this addendum or the claim
registry, the addendum/registry wins.

## Verdict

The previous existence finding must be revised: all 19 DYN source scripts have
now been recovered under `route_E/code_dyn/`.  They are substantial programs,
not empty placeholders.  In a compatibility replay tree, 16 scripts completed
with `219/219` internal checks marked pass.

That does **not** close the Route-E evidence chain:

- the papers cite `code/...`, while the recovered sources are at
  `route_E/code_dyn/...`;
- the moved scripts resolve their root incorrectly and the direct DAG does not
  run from the new location;
- no generated `dyn*.json`/`dyn*.md` ledger is present in the workspace;
- the entire untracked `route_E/` tree is absent from a clean clone;
- Route-E string-card `RE-SC3/4/5` scripts and ledgers are still missing,
  which blocks DYN-9b-2,
  DYN-9b-3, and the DYN-8 collector;
- several green internal checks certify arithmetic or bookkeeping while the
  associated physical inference remains invalid or preliminary.

The strongest newly recoverable result is therefore narrower than the paper's
claim: a **conditional SUSY CMSGUT benchmark** can regenerate its transcribed
spectrum and a script-level obstruction on the declared finite slice.  The
active non-SUSY branch, messenger loop,
flavor posterior, leptogenesis, and string completion are not closed.

The naming is fixed here to avoid a collision: `RE-SC3` is the instanton
Majorana-pricing file named `d3_instanton_majorana_pricing`, `RE-SC4` is the
Stueckelberg-protection file named `d4_stueckelberg_protection`, and `RE-SC5`
is the SUSY-breaking-bridge file named `d5_susy_breaking_bridge`.  Their
historical paths are under `route_d/`, but they are not the original Route-D
workplan items `RD-D3H/RD-D4M/RD-D5P`.

## Replay protocol and numerical result

All 19 files pass `py_compile` with Python `3.9.6`; the numerical replay used
NumPy `2.0.2`.
Direct execution was first tested at the restored location.  The common line

```python
ROOT = Path(__file__).resolve().parents[1]
```

now makes `ROOT=route_e`.  A script that subsequently asks for
`ROOT/route_E/output` therefore searches `route_E/route_E/output`, which does
not exist.  DYN-4a fails this way before producing a ledger.  DYN-8 also looks
for `ROOT/route_d` and `ROOT/paper`, which are wrong after the move.

To separate a deployment failure from an algorithm failure, the sources were
copied read-only to `/tmp/route_e_code_dyn_replay_20260713/code/`, with the
archival `route_E`, `route_d`, and `paper` inputs exposed at the historical
relative locations.  Outputs were written only under that temporary tree.

| Lane | Replay | Internal checks | Evidence class after review |
|---|---:|---:|---|
| DYN-0 | pass | 13/13 | provenance/history regression |
| DYN-1a | pass | 23/23 | conditional SUSY vacuum/Goldstone benchmark |
| DYN-1b | pass | 16/16 | conditional complete CMSGUT spectrum census |
| DYN-2 | pass | 15/15 | conditional SUSY threshold scan |
| DYN-2b | pass | 6/6 | finite benchmark-slice obstruction scan |
| DYN-3 | pass | 13/13 | preliminary proton estimate; archival flavor/dressing |
| DYN-4a | pass | 19/19 | inverse reconstruction; interval gate defective |
| DYN-4b | pass | 13/13 | conditional nuisance-prior distribution, not unconditional |
| DYN-4c | pass | 7/7 | algebraic/kernel perturbation study, not a global fit |
| DYN-5 | arithmetic pass | 21/21 | physical loop claim invalid pending rederivation |
| DYN-7 | arithmetic pass | 13/13 | wrong flavor regime pending flavored replay |
| DYN-9 | pass | 15/15 | one-loop ESH preliminary branch scan |
| DYN-9b-1 | pass | 13/13 | little-group and threshold sensitivity |
| DYN-9b-1b | pass | 14/14 | explicit-tensor exploratory scan |
| DYN-9b-1c | pass | 11/11 | incomplete tree-level quartic scan |
| DYN-9b-1d | pass | 7/7 | rare sampled LR points, not a global vacuum theorem |
| DYN-9b-2 | blocked | exit 1 | missing `RE-SC3` ledger |
| DYN-9b-3 | blocked | exit 1 | missing DYN-9b-2 ledger |
| DYN-8 | blocked | exit 1 | missing `RE-SC3` first; `RE-SC4/5` and downstream ledgers also absent |

The successful total is

\[
13+23+16+15+6+13+19+13+7+21+13+15+13+14+11+7=219.
\]

This count measures the programs' own gates.  Static AST inspection finds 246
`check(...)` calls, of which 39 (`15.85%`) pass a literal `True`; many such
entries mean "recorded" or "disclosed", not an independently tested
proposition.  Internal pass counts must not be read as 219 independent physics
tests.

## What is now numerically supported

### Conditional SUSY benchmark

DYN-1a regenerates the explicit singlet F-flat branch with worst residual
`1.61e-15`, finds 33 eaten Goldstone directions, and exports a non-placeholder
spectrum.  Within the transcribed Aulakh--Girdhar sector formulas and the
specified benchmark, DYN-1b closes the 472-state Higgs census and covers all 26
listed sector types.  This is not an independent full-action Hessian
derivation of the complete spectrum.

DYN-2 finds, on its declared slice, a compatibility point near

\[
x_*=0.15,\qquad
M_X=1.907\times10^{13}\ \mathrm{GeV},\qquad
\Delta\alpha_3\simeq0.81\sigma .
\]

DYN-2b scans 240 `(x, eta, M_S)` combinations, obtains no living point, and
reports best log-margin `-9.46`.  DYN-3 then obtains at the compatibility point

\[
\tau(p\to K^+\bar\nu)\simeq1.24\times10^{26}\ \mathrm{yr}
\]

against the script's bound `5.9e33 yr`, a gap

\[
\log_{10}\!\left(\frac{5.9\times10^{33}}
{1.24\times10^{26}}\right)\simeq7.68.
\]

The finite DYN-2b scan supports a no-go for the **specified minimal SUSY
slice**, with its best point killed primarily by the DYN-2b dimension-six
constraint.  DYN-3's dimension-five estimate is corroborating evidence, not
the primary no-go: it still uses archival Yukawa maps and does not complete the
physical-basis rotations and explicit dressing.  This remains a finite,
branch- and threshold-dependent scan, not a model-class theorem.

### Inverse-seesaw reconstruction

For invertible `m_D` and `m_nu`, the type-I relation

\[
m_\nu=-m_D M_R^{-1}m_D^T
\]

can be inverted algebraically:

\[
M_R=-m_D^T m_\nu^{-1}m_D.
\]

Consequently DYN-4c's exact absorption of refreshed oscillation central values
is a consistency identity once `m_D` and `m_nu` are supplied; it is not an
independent flavor prediction.  The contact decomposition and its numerical
stability under the stated random kernel perturbations are useful
characterizations, but a covariant `Spin(10)` Yukawa sum-rule fit remains open.
If the light-neutrino matrix has a zero mode, this inverse formula must be
replaced by the appropriate rank conditions and generalized-inverse analysis.

The name "unconditional posterior" in DYN-4b is too strong.  It remains
conditional on the archival Dirac and charged-lepton kernels, normal ordering,
`v_u=100 GeV`, and the chosen priors for `m_1` and Majorana phases.  A precise
label is "conditional inverse-seesaw zeta distribution under the stated
nuisance priors".

## Physics failures hidden by green arithmetic gates

### DYN-5 treats a bilinear mass as a Yukawa interaction

The displayed messenger system contains only quadratic terms,

\[
\frac{W}{M_*}=\frac12N^TM_VN+\lambda X^TN
-\frac12X^TK_{\rm tr}^{-1}X.
\]

DYN-5 sets `Y_portal=lambda I` and applies a Yukawa-like anomalous dimension
`Y^dagger Y/(16 pi^2)`.  In a renormalizable supersymmetric theory the standard
one-loop wavefunction term is generated by dimensionless cubic couplings
`y_ijk Phi_i Phi_j Phi_k`; a dimensionful quadratic mass mixing does not become
such a Yukawa tensor merely by calling its coefficient `Y_portal`.  A genuine
interaction such as `S X N`, or the corresponding non-SUSY interaction in the
active branch, must be specified before a loop calculation exists.

Moreover, with canonical Kähler terms, eliminating `X` at tree level in the
low-momentum derivative expansion `p^2/M_X^2 << 1`, while retaining `N` in the
low-energy theory, gives

\[
X=\lambda K_{\rm tr}N,
\qquad
K_{\rm eff}=N^\dagger
\left(I+|\lambda|^2K_{\rm tr}^\dagger K_{\rm tr}\right)N.
\]

Since the explicit `K_tr` is real symmetric with
`K_tr^dagger K_tr=K_tr^2=I/3`, and `lambda^2=zeta`,

\[
\delta Z_{\rm tree}=\frac{|\zeta|}{3}=0.0434773011,
\qquad
\delta Z_{\rm putative\ loop/efold}=\frac{|\zeta|}{16\pi^2}
=8.2597\times10^{-4},
\]

and

\[
\frac{\delta Z_{\rm tree}}{\delta Z_{\rm putative\ loop/efold}}
=\frac{16\pi^2}{3}=52.6379.
\]

The second number is only the coefficient asserted by DYN-5; the displayed
quadratic action does not contain an interaction that generates this loop.

If `N` and `X` are integrated at comparable scales, this derivative expansion
is insufficient; the physical matching object is the fully normalized
six-by-six mass/propagator system and its Weinberg coefficient.

The R-charge enumeration also omits the directly dangerous post-breaking
operator `X L H_u`, whose displayed charges give `1+1+0=2`.  DYN-5 therefore
does not repair Route B and must be classified
`INVALID_PENDING_REDERIVATION`.

### DYN-7 and DYN-9b-3 use the unflavored regime at the wrong scale

The replayed lightest heavy-neutrino scale is about

\[
M_1\simeq2.4\times10^{10}\ \mathrm{GeV}.
\]

For the standard thermal history this lies in the tau-resolved, two-flavor
regime (the exact transition depends on the branch and, in SUSY, on
`tan(beta)`).  The scripts explicitly use an unflavored approximation.
Therefore `14/4000`, `0/4000`, boost probabilities, and a conditioned
`arg(zeta)` interval are order-estimate diagnostics, not physical posterior
probabilities.  The repair is at least a two-flavor Boltzmann calculation and
preferably a density-matrix treatment with spectator effects and branch-local
thermal rates.

## Numerical and software defects

### DYN-4a interval is unresolved and is not a profile

The interval scan fixes the phase while scanning the modulus and fixes the
modulus while scanning the phase.  It does not compute

\[
\chi^2_{\rm prof}(s)=\min_\delta\chi^2(s,\delta),
\qquad
\chi^2_{\rm prof}(\delta)=\min_s\chi^2(s,\delta),
\]

and it does not marginalize a normalized posterior.  The code comment mentions
joint two-parameter thresholds `2.30/5.99`, while the implementation uses
`Delta chi2 <= 1`, appropriate only for a regular one-parameter profile under
additional assumptions.

The modulus-scale grid spacing is `0.005`, compared with a later bisection
width `Delta_s=1.008e-6`; the phase spacing is `0.0025 rad`, compared with
`Delta_delta=1.506e-6 rad`.  Thus the displayed interval grids are
under-resolved by factors

\[
\frac{0.005}{1.008\times10^{-6}}=4.96\times10^3,
\qquad
\frac{0.0025}{1.506\times10^{-6}}=1.66\times10^3.
\]

Both reported 68% intervals collapse to a singleton, but the acceptance gate
only checks that the lists are nonempty.  Replace this with adaptive
minimization, genuine profiling or marginalization, grid-refinement
convergence, and a Hessian/bisection cross-check.

### Rare LR points are an existence hint, not a probability theorem

DYN-9b-1d finds 8 positive points among 8800 sampled `(ratio,coupling)` pairs,
all near ratios `0.6--0.7`.  The point estimate is

\[
\hat p=8/8800=9.09\times10^{-4}.
\]

Even treating the samples as IID Bernoulli trials, a 95% Wilson interval is

\[
p\in[4.61\times10^{-4},1.79\times10^{-3}].
\]

More importantly, this frequency depends on the arbitrary coupling measure.
The scan demonstrates candidate existence within the sampled box, not a
measure-independent viable fraction or a classification of the full potential.
The scripts also retain tree-level-only and omitted-invariant/stationarity
boundaries, so `DYN-9b IS CLOSED` is not justified.

### Evidence deployment and collector integrity

- The historical `code/...` manifest and the restored `route_E/code_dyn/...`
  location disagree.  Case-sensitive Linux also distinguishes `route_E` from
  `route_e`.
- DYN-9b-1c and DYN-9b-1d use fixed global caches
  `/tmp/dyn9b1c_heps_cache.npz` and `/tmp/dyn9b1d_polar_cache.npz` without a
  source, parameter, input, NumPy, or BLAS digest.  Stale, cross-repository,
  and concurrent-cache contamination are possible.
- The scripts have no common CLI/root/output contract and execute work on
  import.  Several do not fail with a nonzero exit status when a check is
  false.
- DYN-8's K5r check correctly says the freed ratio restores rare 210-only LR
  minima, but its emitted `branch_map` still says the branch requires a
  `45_H` route because the 210-only vacuum is dead.  Even a green collector
  would therefore publish a stale conclusion.
- `RE-SC3`, `RE-SC4`, and `RE-SC5` are absent.  They are pricing cards even
  if later restored;
  they are not substitutes for a global divisor, flux, zero-mode, Pfaffian,
  and SUSY-breaking construction.

## Required repair gates

Three implementation strategies can proceed in parallel:

1. **Evidence recovery.**  Introduce a single repository-root resolver and
   lowercase path contract; add `--repo-root`, `--output-dir`, seed, and cache
   arguments; build a dependency-aware runner; make every numerical assertion
   fail closed; store source/input/output/environment digests; and generate all
   ledgers from a clean clone.  Conditional `RE-SC3/4/5` cards must either be
   restored or made optional dependencies in DYN-8.
2. **Physics repair.**  Rebuild DYN-5 from an explicit interacting messenger
   action including tree matching and the complete operator basis; replace
   DYN-7/9b-3 with flavored density-matrix leptogenesis; replace DYN-4 with a
   genuine global `Spin(10)` flavor fit; and complete the active non-SUSY
   scalar potential, vacuum, spectrum, and two-loop thresholds.
3. **Statistical robustness.**  Use adaptive profile/marginal inference,
   multiple seed ensembles, convergence diagnostics, confidence intervals,
   and sensitivity to the coupling prior/measure.  A precise random hit count
   from one seed must not be a physics acceptance gate.

Until these gates pass, Route E has recoverable source code and useful
conditional benchmark diagnostics, but the full theory is **not** closed.
