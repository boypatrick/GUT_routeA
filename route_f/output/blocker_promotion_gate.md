# Route-F blocker promotion gate

18/18 fail-closed checks pass.  This is an expected mechanical pass with `physics_promotion_allowed=false`.

## Resolved artifact blockers

- RE-SC3, RE-SC4, and RE-SC5 now exist and are hash-tracked.
- They remain `unpromoted_pricing_only`; existence is not theory closure.

## Resolved logical blocker

- Original H3 proves no genus selection: `B(e,e)=1` is a nondegenerate invariant contact on the one-dimensional abelian branch.
- The theorem now states `N_fam<=3` unconditionally and `N_fam=3` only under the explicit H3+ Killing-contact axiom.
- `Y_nu=h-3f`, `Y_u=h+f` does not force a top-like neutrino Dirac coupling (`h=3f` is a counterexample).
- Exact zeta invariance holds for uniform positive-real rescaling; complex rescaling rotates its phase by `y^2/|y|^2`.

## Open physics blockers

- H3+: physical motivation or UV realization of the Killing-contact axiom
- DYN-5: explicit interacting messenger action and operator-complete symmetry
- RE-SC4: recompute numerical selectivity pricing after DYN-5 rederivation
- RE-SC5: branch-specific experimental MI floor and threshold/proton envelopes
- DYN-9b-2: global branch-local non-SUSY Spin(10) flavor fit
- DYN-9b-3: tau-resolved kinetic or density-matrix evolution with thermal inputs

## Numerical repairs

- Correct `y_t(M_X)` = 0.44116481; historical double-normalised result = 0.41556342.
- Correct SM Davidson--Ibarra bound = 2.307794855090e-06.
- Fixed archival-kernel suppression = 19.5x (PS) / 342.2x (G_LR); optional top-like tension = 9.6x / 169.2x.
- RE-SC5 toy lower edge has `M_I=101.31` GeV; the illustrative `M_I>=10 TeV` grid edge begins at `log10(M_SS)=11.85`.

## Checks

- [PASS] one-dimensional abelian algebra defeats original H3 but fails H3+ — B(e,e)=1 is invariant/nondegenerate; Killing=0; only N_fam<=3 is unconditional
- [PASS] RE-SC3/4/5 exist and are mechanically green — re_sc3=15/15, re_sc4=14/14, re_sc5=18/18
- [PASS] the recovered string cards remain unpromoted pricing cards
- [PASS] string-card boundary flags expose their non-derived inputs
- [PASS] allowed NLH and XN force q(XLH)=q(XX) for every sampled additive charge — 121 exact integer assignments
- [PASS] RE-SC4 numerical selectivity gaps are blocked by invalid DYN-5 input — the Abelian no-go survives; Delta-S pricing does not
- [PASS] DYN-9b-2 uses already-normalised g1=0.462 with b1=41/10 — correct=0.44116481, double-normalised=0.41556342
- [PASS] SO(10) matrix relations do not promote the optional top-like ansatz — h=3f gives Y_nu=0 while Y_u=4f
- [PASS] fixed archival-kernel suppression is separated from top-like tension — archival=19.5/342.2; top-like=9.6/169.2
- [PASS] zeta invariance is domain-limited and complex rescaling is phase-covariant — complex covariance error=1.20e-16
- [PASS] SM Davidson--Ibarra bound uses physical m3-m1 after one square root — epsilon_DI_SM=2.307794855090e-06
- [PASS] DYN-9b-3 remains blocked in the tau-resolved regime
- [PASS] published toy-window lower edge is only about 101 GeV — log10(MSS)=10.30, MI=101.31 GeV
- [PASS] illustrative MI>=10 TeV diagnostic moves the grid edge to log10(MSS)=11.85 — diagnostic only; an experimental branch-specific bound is still required
- [PASS] RE-SC5 toy window remains non-promotable despite its mechanical scan
- [PASS] DYN-9b-2 is preliminary rather than a claimed global flavor fit
- [PASS] DYN-8 passes disclosure checks while forbidding physics promotion
- [PASS] DYN-8 contains no stale promotion or exhaustive-source wording
