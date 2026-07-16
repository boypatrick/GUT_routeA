# AP-E2 Discriminant / Contact Regression

Status: **ap_e2_exact_regression_pass_no_physics_promotion** — `30/30` checks pass.

This is a deterministic, fail-closed algebraic regression.  A green result
preserves the convention

```text
B_Kill(A_q,A_q) = 2 Delta(p) = 2 sqrt(3) x^T K_tr x,
A_q^2 = Delta(p) I / 4,
(p,p)_2 = -2 Delta(p).
```

It does not derive H3+, a microscopic carrier, or any dynamics, and therefore
sets `physics_promotion_allowed=false`.

## Coverage

- exact rational anchor and `256` seeded exact-rational samples;
- `256` seeded complex samples at `100` decimal places;
- exact `SL(2)` covariance and matrix conjugacy, polarized bilinear identities,
  scalar covariance, and uniqueness of the invariant symmetric line;
- real, complex, and infinite regular roots, the finite and infinite
  double-root charts, a `1e-70` near-null perturbation, and the theorem's
  explicit exclusion of the zero element;
- Route-E `K_tr` source-card normalization and full `H,E,F` invariance;
- explicit spherical-to-normalized-polynomial basis conversion;
- four required negative controls: missing `sqrt(2)`, missing factor two,
  using the spherical coefficient on the normalized basis, and replacing the
  complex bilinear contact by a Hermitian form.

Maximum intrinsic 100-digit complex-sample residual:
`1.67930436345868921e-101`.
The separately reported Route-E source-card residual is
`6.14719575160384581e-17` because that JSON
anchor stores ordinary double-precision decimals.

## Negative controls

- missing `sqrt(2)` detected: `true`
- missing Killing factor two detected: `true`
- normalized-basis factor mismatch detected: `true`
- Hermitian-for-bilinear substitution detected: `true`

These are tests that must disagree with the correct theorem.  Their detection
is a PASS; accidental agreement would fail the run.

## Mechanical checks

- [PASS] `A0_provenance` — all critical AP-E2 sources are present and hashed: hashed=6/6
- [PASS] `A0_provenance` — Route-E invariant card verifies its internal stable digest: stored=a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc; recomputed=a68f7c93e3eca0b63c95d7f7ac8a281fc58c429040fc98f1eca1846efb1d18bc
- [PASS] `A0_provenance` — Route-E card binds the current Audit-0 builder hash: card=0e4215ba4f35c00623045da89cb5b4240e3bd648b1ba3205c0772de00a11ec81; current=0e4215ba4f35c00623045da89cb5b4240e3bd648b1ba3205c0772de00a11ec81
- [PASS] `A0_provenance` — the AP-E2 theorem/proof block matches the frozen semantic source: block_sha256=3999b88f3fe7e396948c9695939ec5826d74cfc6d06f6b26ed05d3b0e0c4d094
- [PASS] `A1_exact` — exact rational anchor obeys A_q^2=Delta I/4: (a,b,c)=(Fraction(3, 2), Fraction(-7, 3), Fraction(5, 4)); Delta=-37/18
- [PASS] `A1_exact` — exact rational anchor obeys B=2 Delta and (p,p)_2=-2 Delta: Delta=-37/18; B=-37/9; transvectant=37/9
- [PASS] `A1_exact` — polarized anchor obeys B(A_p,A_r)=2D(p,r) and (p,r)_2=-B(A_p,A_r): p=(Fraction(3, 2), Fraction(-7, 3), Fraction(5, 4)); r=(Fraction(-2, 5), Fraction(11, 7), Fraction(4, 3))
- [PASS] `A1_exact` — the symmetric sl2-invariant contact line is exactly one-dimensional: constraint rank=5/6; kernel dimension=1
- [PASS] `A1_exact` — seeded rational samples obey all projective/Killing/transvectant identities exactly: seed=20260714; samples=256; failures=0; regular=256
- [PASS] `A1_exact` — seeded rational pairs obey the polarized contact identity exactly: samples=256; failures=0
- [PASS] `A1_exact` — det(z,A_p z)=p(z) fixes the homogeneous and affine-chart convention: samples=256; failures=0; affine xi=v/u
- [PASS] `A1_exact` — binary discriminant is exactly invariant under seeded SL(2,Z) substitutions: samples=256; failures=0
- [PASS] `A1_exact` — the projective matrix transforms by exact conjugacy under the declared SL(2) substitution: A[p(gz)]=g^-1 A[p] g; samples=256; failures=0
- [PASS] `A1_exact` — discriminant and contact have exact weight two under section rescaling: samples=256; failures=0
- [PASS] `A2_boundary` — finite and infinite double-root charts lie on the nonzero nilpotent/contact-null cone: cases=4; failures=0; includes u^2 and v^2
- [PASS] `A2_boundary` — real, complex, asymmetric, and infinite roots match eigenlines in the xi=v/u chart: cases=4; max residual=0.00e+00; min projective separation=0.316228
- [PASS] `A2_boundary` — finite zero, finite nonzero, and infinite double roots have vanishing first jet: cases=3; max section/eigenline/gradient residual=0.00e+00
- [PASS] `A2_boundary` — the zero section is detected and excluded from the theorem's nonzero nilpotent clause: Delta(0)=0 and A_0=0; theorem domain is explicitly X_q != 0
- [PASS] `A2_boundary` — 100-digit arithmetic resolves a 1e-70 perturbation away from the null cone: Delta=-4.0e-70; identity residual=0.0
- [PASS] `A3_route_e` — Route-E card uses the canonical spherical K_tr normalization: max|K_source-K_canonical|=3.54908512195e-17; max|K^2-I/3|=4.09813050107e-17
- [PASS] `A3_route_e` — K_tr is symmetric, invertible, and may have nonzero isotropic vectors without being degenerate: sym=0.0; inverse=1.22943915032e-16; det=(-1.9245008973e-1 + 0.0j); isotropic norm=(0.0 + 0.0j); |Kx|_max=5.7735026919e-1
- [PASS] `A3_route_e` — spherical and normalized-polynomial bases carry the required factor-two conversion: y=R*x/sqrt(2); B=2sqrt(3)x^TKx=4sqrt(3)y^TKy; residual=3.31036622241e-15
- [PASS] `A3_route_e` — loaded K_tr is invariant under the complete H,E,F generator triple: residuals={"E": "0.0", "F": "0.0", "H": "0.0"}
- [PASS] `A4_complex` — seeded 100-digit complex samples reproduce the discriminant/contact theorem: seed=20260714; samples=256; regular=256; failures=0; intrinsic max=1.67930436346e-101; source-card max=6.1471957516e-17
- [PASS] `A5_negative_control` — missing sqrt(2) in the spherical section map is rejected: intentional wrong-convention residual=1.0e+1
- [PASS] `A5_negative_control` — missing factor two in B_Kill=2 Delta is rejected: intentional wrong-convention residual=1.2e+1
- [PASS] `A5_negative_control` — using the spherical 2sqrt(3) factor on normalized polynomial coefficients is rejected: intentional wrong-basis residual=1.2e+1
- [PASS] `A5_negative_control` — Hermitian substitution for the complex symmetric contact is rejected: intentional wrong-convention residual=1.56204993518e+1
- [PASS] `A6_boundary` — the one-dimensional abelian counterexample passes invariant-pairing H3 but fails Killing-contact H3+: ad_e=0 gives Killing=0 while B(e,e)=1 is invariant and nondegenerate
- [PASS] `A6_boundary` — a green AP-E2 algebraic regression cannot promote a physical bridge: physics_promotion_allowed=false by construction

## Source hashes

- `route_E/output/audit0/invariant_card.json` — `891df84534327029b3eed227381b4cc0d24beb41ec1988a610051d0e53c6f131` (10028 bytes)
- `route_E/code/audit0_conventions_card.py` — `0e4215ba4f35c00623045da89cb5b4240e3bd648b1ba3205c0772de00a11ec81` (15657 bytes)
- `route_f/tex/another_physics_route_e_derivation_ledger.tex` — `e514c4b2bde0f7b9c1c2e07e5f86f46d76b1d7797262bbdd449abe8b095997f6` (88057 bytes)
- `route_f/tex/ap_e2_discriminant_regression.tex` — `0c54f461b6e4f1d0f5cebd88dc384b14295add78fd64aeeddcf2c1d5c4b59092` (31296 bytes)
- `route_f/code/verify_another_physics_route_e_bridge.py` — `a5b070b09a1b82deb8f40cd93ff032df45368eb0085b2b4b2a8212442f4c6ada` (42757 bytes)
- `route_f/code/verify_ap_e2_discriminant_regression.py` — `a658e2804509502a66c37cb4f8e2318329a7fd285266bcfd2da70ed92e489276` (40639 bytes)

The theorem-block hash lock is
`3999b88f3fe7e396948c9695939ec5826d74cfc6d06f6b26ed05d3b0e0c4d094`.

## Promotion boundary

- `physics_promotion_allowed=false`
- AP-E2 freezes an exact conditional invariant theorem only.  AP-E3 must
  independently derive any Berry level or `O(2)` prequantum bundle.
